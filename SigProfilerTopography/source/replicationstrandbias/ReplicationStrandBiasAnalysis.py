# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018-2020 Burcak Otlu


# Version2
# This version use np.arrays
# Right now replication strand bias analysis works for single point mutations and signatures.
# This python code analyses the Replication Strand Bias

import multiprocessing
import numpy as np
import pandas as pd

from SigProfilerTopography.source.commons.TopographyCommons import CHROM
from SigProfilerTopography.source.commons.TopographyCommons import START
from SigProfilerTopography.source.commons.TopographyCommons import END
from SigProfilerTopography.source.commons.TopographyCommons import SIGNAL

from SigProfilerTopography.source.commons.TopographyCommons import PYRAMIDINESTRAND
from SigProfilerTopography.source.commons.TopographyCommons import SAMPLE

from SigProfilerTopography.source.commons.TopographyCommons import SUBS
from SigProfilerTopography.source.commons.TopographyCommons import INDELS
from SigProfilerTopography.source.commons.TopographyCommons import DINUCS

from SigProfilerTopography.source.commons.TopographyCommons import MUTATION
from SigProfilerTopography.source.commons.TopographyCommons import LENGTH

from SigProfilerTopography.source.commons.TopographyCommons import LEADING
from SigProfilerTopography.source.commons.TopographyCommons import LAGGING
from SigProfilerTopography.source.commons.TopographyCommons import REPLICATIONSTRANDBIAS


from SigProfilerTopography.source.commons.TopographyCommons import updateDictionaries_simulations_integrated

from SigProfilerTopography.source.commons.TopographyCommons import readWig_with_fixedStep_variableStep

from SigProfilerTopography.source.commons.TopographyCommons import readFileInBEDFormat


from SigProfilerTopography.source.commons.TopographyCommons import getDictionary
from SigProfilerTopography.source.commons.TopographyCommons import readChrBasedMutationsDF
from SigProfilerTopography.source.commons.TopographyCommons import accumulate_simulations_integrated
from SigProfilerTopography.source.commons.TopographyCommons import accumulate_simulations_integrated_for_each_tuple

from SigProfilerTopography.source.commons.TopographyCommons import writeDictionary

from SigProfilerTopography.source.commons.TopographyCommons import SubsSignature2PropertiesListDictFilename
from SigProfilerTopography.source.commons.TopographyCommons import IndelsSignature2PropertiesListDictFilename
from SigProfilerTopography.source.commons.TopographyCommons import DinucsSignature2PropertiesListDictFilename

from SigProfilerTopography.source.commons.TopographyCommons import Type2ReplicationStrand2CountDict_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Signature2MutationType2ReplicationStrand2CountDict_Filename

from SigProfilerTopography.source.commons.TopographyCommons import Type2Sample2ReplicationStrand2CountDict_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Sample2Type2ReplicationStrand2CountDict_Filename

from SigProfilerTopography.source.commons.TopographyCommons import COMPUTATION_CHROMOSOMES_SEQUENTIAL_ALL_SIMULATIONS_PARALLEL
from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC

#For Supp Fig2B
CHR10_THRESHOLD_START = 16400000
CHR10_THRESHOLD_END = 26400000

#For Supp Fig2A
CHR20_START = 36260000
CHR20_END = 36830000

# FOR FINDING TRANSITION ZONES (LEADING or LAGGING)
# THRESHOLD_CONSECUTIVE_LONG_STRETCH_LENGTH= 250000 #used in Supp Fig2B
# THRESHOLD_CONSECUTIVE_LONG_STRETCH_LENGTH= 150000
THRESHOLD_CONSECUTIVE_LONG_STRETCH_LENGTH= 10000

# THRESHOLD_DISCARD_LATEST_TRANSITION_ZONE = 100000 #used in Supp Fig2B
THRESHOLD_DISCARD_LATEST_TRANSITION_ZONE = 25000
# THRESHOLD_LATEST_TRANSITION_ZONE = 0


########################################################################
def checkForSameSignedSlopeBetweenConsecutivePeakandValley(chrLong,peakorValleyStart, peakorValleyEnd, chrBasedSmoothedWaveletReplicationTimeSignalDF):
    transitionZoneList =[]

    # print('################ checkForConsecutive starts ############ fromStart: %s toEnd: %s' %(peakorValleyStart,peakorValleyEnd))
    subset_df = chrBasedSmoothedWaveletReplicationTimeSignalDF[(chrBasedSmoothedWaveletReplicationTimeSignalDF[START]>=peakorValleyStart) & (chrBasedSmoothedWaveletReplicationTimeSignalDF[END]<=peakorValleyEnd)]

    consecutiveLength = 0
    formerRow= None
    formerSlopeDirection = None

    start = peakorValleyStart

    for index,row in subset_df.iterrows():
        if formerRow is None:
            #We read the row for the first time
            formerRow = row
            consecutiveLength += 1000
        else:
            slope = (row.get(SIGNAL) - formerRow.get(SIGNAL)) / 1000
            formerRow = row

            if (formerSlopeDirection is None):
                formerSlopeDirection = np.sign(slope)
                consecutiveLength += 1000
            elif (formerSlopeDirection==np.sign(slope)):
                consecutiveLength += 1000
            else:
                #They have different signs
                if (consecutiveLength>=THRESHOLD_CONSECUTIVE_LONG_STRETCH_LENGTH):
                    # print('Slope sign changed -- Found one: from %d to %d with %d bases with slope sign %s' %(start,((row.get('start') + row.get('end'))//2), consecutiveLength, formerSlopeDirection))
                    transitionZoneList.append((chrLong,start,(row.get(START) + row.get(END))//2,formerSlopeDirection,consecutiveLength))
                #initialize and start again
                consecutiveLength = 1000
                start = (row.get(START) + row.get(END))//2
                formerRow= row
                formerSlopeDirection= np.sign(slope)
                continue

            # print('slope: %f - np.sign(slope): %f -  consecutiveLength: %d ' %(slope,np.sign(slope),consecutiveLength))
            formerSlopeDirection = np.sign(slope)

    #This is for the last probable transition zone.
    if (consecutiveLength >= THRESHOLD_CONSECUTIVE_LONG_STRETCH_LENGTH):
        # print('After for loop ends, found one: from %d to %s with %d bases with slope sign %s' % (start, (row.get('start') + row.get('end'))//2, consecutiveLength, formerSlopeDirection))
        transitionZoneList.append((chrLong,start,(row.get(START) + row.get(END))//2,formerSlopeDirection,consecutiveLength))

    # print('################ checkForConsecutive ends ############ fromStart: %s toEnd: %s' % (peakorValleyStart,peakorValleyEnd))
    return transitionZoneList
########################################################################


########################################################################
# chr10_subset_wavelet_processed_df
#           chr     start       end   signal
# 265577  chr10  16400500  16401499  24.9438
# 265578  chr10  16401500  16402499  24.9585

# valleys_peaks_df
#         chr     start       end    type
# 415     chr10  16454500  16455500    Peak
# 415  chr10  16528500  16529500  Valley

def findLongStretchesofConsistentTransitionZones(chrLong,fromStart,toEnd,chrBasedSmoothedWaveletReplicationTimeSignalDF,valleys_peaks_df):
    transitionZonesList =[]
    for index,row in  valleys_peaks_df.iterrows():
        peakorValleyStart = row[START]
        peakorValleyEnd = row[END]
        peakorValleyMidpoint = (peakorValleyStart+peakorValleyEnd)//2

        type = row['type']
        if (type =='Peak'):
            if (peakorValleyMidpoint>fromStart):
                # print('from: %d - to: %d - difference: %d'  %(fromStart,peakorValleyMidpoint, (peakorValleyMidpoint-fromStart)))
                found = checkForSameSignedSlopeBetweenConsecutivePeakandValley(chrLong,fromStart, peakorValleyMidpoint, chrBasedSmoothedWaveletReplicationTimeSignalDF)
                transitionZonesList.extend(found)
                # print('found %s' %found)
            fromStart=peakorValleyMidpoint
        elif (type=='Valley'):
            valleyStart =row[START]
            valleyEnd = row[END]
            valleyMidpoint = (valleyStart+valleyEnd)//2
            # This is something special to valley
            newValleyStart1 = valleyMidpoint - THRESHOLD_DISCARD_LATEST_TRANSITION_ZONE
            newValleyStart2 = valleyMidpoint + THRESHOLD_DISCARD_LATEST_TRANSITION_ZONE
            if (newValleyStart1>fromStart):
                # print('from: %d - to: %d - difference: %d' % (fromStart, newValleyStart1, (newValleyStart1 - fromStart)))
                found = checkForSameSignedSlopeBetweenConsecutivePeakandValley(chrLong,fromStart, newValleyStart1,chrBasedSmoothedWaveletReplicationTimeSignalDF)
                transitionZonesList.extend(found)
                # print('found %s' % found)
            # bypass the genome region between newValleyStart1 and newValleyStart2
            fromStart = newValleyStart2
    #
    #For the last interval
    if (toEnd>fromStart):
        # print('last one from: %d - to: %d -difference: %d' %(fromStart,toEnd,(toEnd-fromStart)))
        found = checkForSameSignedSlopeBetweenConsecutivePeakandValley(chrLong,fromStart, toEnd, chrBasedSmoothedWaveletReplicationTimeSignalDF)
        transitionZonesList.extend(found)
        # print('found %s' %found)

    return transitionZonesList
########################################################################


########################################################################
#TODO Is (replicationStrand_row['end']+1) okey?
#We assume that there are no overlapping intervals with positive and negative slopes.
#To test it have one array for positive slope fill with 1s
#                one array for negative slope fill with -2a
#                add them if you habe any -1 that means that you contradict this assumption.
def fillReplicationStrandArray(replicationStrand_row,chrBased_replication_array):
    # e.g.: replicationStrand_row
    # chr chrX
    # start   154861998
    # end 155096999
    # slopeDirection  1 (1 means leading strand -1 means lagging strand on positive strand)
    # length  235000

    # labels = ['chr', 'start', 'end', 'slopeDirection', 'length']
    chrBased_replication_array[replicationStrand_row['start']:replicationStrand_row['end']+1] = replicationStrand_row['slopeDirection']
########################################################################

########################################################################
# Summary:
#   if mutationPyramidineStrand and slope have the same sign increase LEADING STRAND count
#   else mutationPyramidineStrand and slope have the opposite sign increase LAGGING STRAND count
def searchMutationOnReplicationStrandArray_simulations_integrated(
        mutation_row,
        chrBasedReplicationArray,
        simNum2Type2ReplicationStrand2CountDict,
        simNum2Sample2Type2ReplicationStrand2CountDict,
        simNum2Type2Sample2ReplicationStrand2CountDict,
        simNum2Signature2MutationType2ReplicationStrand2CountDict,
        signature2PropertiesListDict,
        type,
        sample_based):

    start = mutation_row[START]
    mutationType = None

    pyramidineStrand = mutation_row[PYRAMIDINESTRAND]
    sample = mutation_row[SAMPLE]

    if(type==SUBS):
        end = start+1
        #e.g.: C>A
        mutationType = mutation_row[MUTATION]
    if (type==INDELS):
        end = start+mutation_row[LENGTH]
    elif (type==DINUCS):
        end = start+2

    #############################################################################################################
    #if there is overlap with chrBasedReplicationArray
    slicedArray = chrBasedReplicationArray[start:end]

    if (np.any(slicedArray)):

        #It must be full with at most -1 and +1
        uniqueValueArray = np.unique(slicedArray[np.nonzero(slicedArray)])

        if (uniqueValueArray.size>2):
            print('There is a situation!!!')

        elif ((uniqueValueArray.size==2) and (pyramidineStrand!=0)):
            #Increment both LEADING and LAGGING
            updateDictionaries_simulations_integrated(mutation_row,
                                        mutationType,
                                        sample,
                                        sample_based,
                                        simNum2Type2ReplicationStrand2CountDict,
                                        simNum2Sample2Type2ReplicationStrand2CountDict,
                                        simNum2Type2Sample2ReplicationStrand2CountDict,
                                        simNum2Signature2MutationType2ReplicationStrand2CountDict,
                                        LAGGING,
                                        signature2PropertiesListDict)

            updateDictionaries_simulations_integrated(mutation_row,
                                        mutationType,
                                        sample,
                                        sample_based,
                                        simNum2Type2ReplicationStrand2CountDict,
                                        simNum2Sample2Type2ReplicationStrand2CountDict,
                                        simNum2Type2Sample2ReplicationStrand2CountDict,
                                        simNum2Signature2MutationType2ReplicationStrand2CountDict,
                                        LEADING,
                                        signature2PropertiesListDict)


        # I expect the value of 1 (LEADING on the positive strand) or -1 (LAGGING on the positive strand) so size must be one.
        elif (uniqueValueArray.size == 1):
            for uniqueValue in np.nditer(uniqueValueArray):
                # type(decileIndex) is numpy.ndarray
                slope = int(uniqueValue)

                #They have the same sign, multiplication (1,1) (-1,-1) must be 1
                if (slope*pyramidineStrand > 0):
                    updateDictionaries_simulations_integrated(mutation_row,
                                            mutationType,
                                            sample,
                                            sample_based,
                                            simNum2Type2ReplicationStrand2CountDict,
                                            simNum2Sample2Type2ReplicationStrand2CountDict,
                                            simNum2Type2Sample2ReplicationStrand2CountDict,
                                            simNum2Signature2MutationType2ReplicationStrand2CountDict,
                                            LEADING,
                                            signature2PropertiesListDict)

                # They have the opposite sign, multiplication(1,-1) (-1,-)  must be -1
                elif (slope*pyramidineStrand < 0):
                    updateDictionaries_simulations_integrated(mutation_row,
                                            mutationType,
                                            sample,
                                            sample_based,
                                            simNum2Type2ReplicationStrand2CountDict,
                                            simNum2Sample2Type2ReplicationStrand2CountDict,
                                            simNum2Type2Sample2ReplicationStrand2CountDict,
                                            simNum2Signature2MutationType2ReplicationStrand2CountDict,
                                            LAGGING,
                                            signature2PropertiesListDict)
        else:
            print('There is a situation!!!')
    #############################################################################################################
########################################################################


########################################################################
def searchMutationsOnReplicationStrandArrayForApplyAsync(
    chrBased_replication_array,
    chrBased_simBased_subs_df,
    chrBased_simBased_indels_df,
    chrBased_simBased_dinucs_df,
    numofSimulations,
    sample_based,
    subsSignature2PropertiesListDict,
    indelsSignature2PropertiesListDict,
    dinucsSignature2PropertiesListDict):

    simNum2Type2ReplicationStrand2CountDict= {}
    simNum2Sample2Type2ReplicationStrand2CountDict= {}
    simNum2Type2Sample2ReplicationStrand2CountDict = {}
    simNum2Signature2MutationType2ReplicationStrand2CountDict = {}

    for simNum in range(0,numofSimulations+1):
        simNum2Type2ReplicationStrand2CountDict[simNum]={}
        simNum2Sample2Type2ReplicationStrand2CountDict[simNum]={}
        simNum2Type2Sample2ReplicationStrand2CountDict[simNum]={}
        simNum2Signature2MutationType2ReplicationStrand2CountDict[simNum]={}

    ##############################  Fill dictionaries for subs  starts ####################
    if ((chrBased_simBased_subs_df is not None) and (not chrBased_simBased_subs_df.empty)):
        chrBased_simBased_subs_df.apply(searchMutationOnReplicationStrandArray_simulations_integrated,
                                chrBasedReplicationArray=chrBased_replication_array,
                                simNum2Type2ReplicationStrand2CountDict=simNum2Type2ReplicationStrand2CountDict,
                                simNum2Sample2Type2ReplicationStrand2CountDict=simNum2Sample2Type2ReplicationStrand2CountDict,
                                simNum2Type2Sample2ReplicationStrand2CountDict=simNum2Type2Sample2ReplicationStrand2CountDict,
                                simNum2Signature2MutationType2ReplicationStrand2CountDict=simNum2Signature2MutationType2ReplicationStrand2CountDict,
                                signature2PropertiesListDict=subsSignature2PropertiesListDict,
                                type=SUBS,
                                sample_based=sample_based,
                                axis=1)
    ##############################  Fill dictionaries for subs  ends ######################


    ##############################  Fill dictionaries for indels  starts ####################
    if ((chrBased_simBased_indels_df is not None) and (not chrBased_simBased_indels_df.empty)):
        chrBased_simBased_indels_df.apply(searchMutationOnReplicationStrandArray_simulations_integrated,
                                chrBasedReplicationArray=chrBased_replication_array,
                                simNum2Type2ReplicationStrand2CountDict=simNum2Type2ReplicationStrand2CountDict,
                                simNum2Sample2Type2ReplicationStrand2CountDict=simNum2Sample2Type2ReplicationStrand2CountDict,
                                simNum2Type2Sample2ReplicationStrand2CountDict=simNum2Type2Sample2ReplicationStrand2CountDict,
                                simNum2Signature2MutationType2ReplicationStrand2CountDict=simNum2Signature2MutationType2ReplicationStrand2CountDict,
                                signature2PropertiesListDict=indelsSignature2PropertiesListDict,
                                type=INDELS,
                                sample_based=sample_based,
                                axis=1)
    ##############################  Fill dictionaries for indels  ends ######################

    ##############################  Fill dictionaries for indels  starts ####################
    if ((chrBased_simBased_dinucs_df is not None) and (not chrBased_simBased_dinucs_df.empty)):
        chrBased_simBased_dinucs_df.apply(searchMutationOnReplicationStrandArray_simulations_integrated,
                                chrBasedReplicationArray=chrBased_replication_array,
                                simNum2Type2ReplicationStrand2CountDict=simNum2Type2ReplicationStrand2CountDict,
                                simNum2Sample2Type2ReplicationStrand2CountDict=simNum2Sample2Type2ReplicationStrand2CountDict,
                                simNum2Type2Sample2ReplicationStrand2CountDict=simNum2Type2Sample2ReplicationStrand2CountDict,
                                simNum2Signature2MutationType2ReplicationStrand2CountDict=simNum2Signature2MutationType2ReplicationStrand2CountDict,
                                signature2PropertiesListDict=dinucsSignature2PropertiesListDict,
                                type=DINUCS,
                                sample_based=sample_based,
                                axis=1)
    ##############################  Fill dictionaries for indels  ends ######################


    #Fill the type2replicationStranCount dictionaries and return them
    return (simNum2Type2ReplicationStrand2CountDict,
            simNum2Sample2Type2ReplicationStrand2CountDict,
            simNum2Type2Sample2ReplicationStrand2CountDict,
            simNum2Signature2MutationType2ReplicationStrand2CountDict)
########################################################################


########################################################################
def  searchMutationsOnReplicationStrandArray(inputList):
    chrBased_replication_array = inputList[0]
    chrBased_subs_split_df = inputList[1]
    chrBased_indels_split_df = inputList[2]
    chrBased_dinucs_split_df = inputList[3]
    numofSimulations = inputList[4]
    sample_based = inputList[5]
    subsSignature2PropertiesListDict = inputList[6]
    indelsSignature2PropertiesListDict = inputList[7]
    dinucsSignature2PropertiesListDict = inputList[8]

    simNum2Type2ReplicationStrand2CountDict= {}
    simNum2Sample2Type2ReplicationStrand2CountDict= {}
    simNum2Type2Sample2ReplicationStrand2CountDict = {}
    simNum2Signature2MutationType2ReplicationStrand2CountDict = {}

    for simNum in range(0,numofSimulations+1):
        simNum2Type2ReplicationStrand2CountDict[simNum]={}
        simNum2Sample2Type2ReplicationStrand2CountDict[simNum]={}
        simNum2Type2Sample2ReplicationStrand2CountDict[simNum]={}
        simNum2Signature2MutationType2ReplicationStrand2CountDict[simNum]={}

    ##############################  Fill dictionaries for subs  starts ####################
    if ((chrBased_subs_split_df is not None) and (not chrBased_subs_split_df.empty)):
        chrBased_subs_split_df.apply(searchMutationOnReplicationStrandArray_simulations_integrated,
                                chrBasedReplicationArray=chrBased_replication_array,
                                simNum2Type2ReplicationStrand2CountDict=simNum2Type2ReplicationStrand2CountDict,
                                simNum2Sample2Type2ReplicationStrand2CountDict=simNum2Sample2Type2ReplicationStrand2CountDict,
                                simNum2Type2Sample2ReplicationStrand2CountDict=simNum2Type2Sample2ReplicationStrand2CountDict,
                                simNum2Signature2MutationType2ReplicationStrand2CountDict=simNum2Signature2MutationType2ReplicationStrand2CountDict,
                                signature2PropertiesListDict=subsSignature2PropertiesListDict,
                                type=SUBS,
                                sample_based=sample_based,
                                axis=1)
    ##############################  Fill dictionaries for subs  ends ######################


    ##############################  Fill dictionaries for indels  starts ####################
    if ((chrBased_indels_split_df is not None) and (not chrBased_indels_split_df.empty)):
        chrBased_indels_split_df.apply(searchMutationOnReplicationStrandArray_simulations_integrated,
                                chrBasedReplicationArray=chrBased_replication_array,
                                simNum2Type2ReplicationStrand2CountDict=simNum2Type2ReplicationStrand2CountDict,
                                simNum2Sample2Type2ReplicationStrand2CountDict=simNum2Sample2Type2ReplicationStrand2CountDict,
                                simNum2Type2Sample2ReplicationStrand2CountDict=simNum2Type2Sample2ReplicationStrand2CountDict,
                                simNum2Signature2MutationType2ReplicationStrand2CountDict=simNum2Signature2MutationType2ReplicationStrand2CountDict,
                                signature2PropertiesListDict=indelsSignature2PropertiesListDict,
                                type=INDELS,
                                sample_based=sample_based,
                                axis=1)
    ##############################  Fill dictionaries for indels  ends ######################

    ##############################  Fill dictionaries for indels  starts ####################
    if ((chrBased_dinucs_split_df is not None) and (not chrBased_dinucs_split_df.empty)):
        chrBased_dinucs_split_df.apply(searchMutationOnReplicationStrandArray_simulations_integrated,
                                chrBasedReplicationArray=chrBased_replication_array,
                                simNum2Type2ReplicationStrand2CountDict=simNum2Type2ReplicationStrand2CountDict,
                                simNum2Sample2Type2ReplicationStrand2CountDict=simNum2Sample2Type2ReplicationStrand2CountDict,
                                simNum2Type2Sample2ReplicationStrand2CountDict=simNum2Type2Sample2ReplicationStrand2CountDict,
                                simNum2Signature2MutationType2ReplicationStrand2CountDict=simNum2Signature2MutationType2ReplicationStrand2CountDict,
                                signature2PropertiesListDict=dinucsSignature2PropertiesListDict,
                                type=DINUCS,
                                sample_based=sample_based,
                                axis=1)
    ##############################  Fill dictionaries for indels  ends ######################


    #Fill the type2replicationStranCount dictionaries and return them
    return (simNum2Type2ReplicationStrand2CountDict,
            simNum2Sample2Type2ReplicationStrand2CountDict,
            simNum2Type2Sample2ReplicationStrand2CountDict,
            simNum2Signature2MutationType2ReplicationStrand2CountDict)
########################################################################



########################################################################
#This code checks whether valleys and peaks are one after another, not two consecutive elements are both valley and peak.
def checkforValidness(chrBased_valleys_peaks_df):
    formerRowType = None

    for index, row in chrBased_valleys_peaks_df.iterrows():
        if formerRowType is None:
            formerRowType = row['type']
        elif (row['type']== formerRowType):
            return False
        else:
            formerRowType = row['type']

    return True
########################################################################



########################################################################
def fill_chr_based_replication_strand_array(chrLong,
                                            chromSize,
                                            chrBasedSmoothedWaveletReplicationTimeSignalDF,
                                            chrBased_valleys_peaks_df):
    # +1 means leading strand, -1 means lagging strand
    # we will fill this array using smoothedSignal, peaks and valleys for each chromosome
    chrBased_replication_array = np.zeros(chromSize, dtype=np.int8)

    firstIndex = chrBasedSmoothedWaveletReplicationTimeSignalDF.index[0]
    lastIndex = chrBasedSmoothedWaveletReplicationTimeSignalDF.index[-1]

    startColumnIndex = chrBasedSmoothedWaveletReplicationTimeSignalDF.columns.get_loc(START)
    endColumnIndex = chrBasedSmoothedWaveletReplicationTimeSignalDF.columns.get_loc(END)

    start = chrBasedSmoothedWaveletReplicationTimeSignalDF.iloc[0, startColumnIndex]  # get the first row start
    end = chrBasedSmoothedWaveletReplicationTimeSignalDF.iloc[-1, endColumnIndex]  # get the last row end

    # Step1 Find the transition zones
    chrBasedTransitionZonesList = findLongStretchesofConsistentTransitionZones(chrLong,start,end,
                                                                               chrBasedSmoothedWaveletReplicationTimeSignalDF,
                                                                               chrBased_valleys_peaks_df)

    labels = ['chr', 'start', 'end', 'slopeDirection', 'length']
    chrBasedTransitionZonesDF = pd.DataFrame.from_records(chrBasedTransitionZonesList, columns=labels)

    # Step2 Fill the replication array using transition zones
    chrBasedTransitionZonesDF.apply(fillReplicationStrandArray, chrBased_replication_array=chrBased_replication_array,axis=1)

    return chrBased_replication_array
########################################################################


########################################################################
def read_repliseq_dataframes(smoothedWaveletRepliseqDataFilename,valleysBEDFilename,peaksBEDFilename):

    ################### Read the Smoothed Wavelet Replication Time Signal starts ###########################
    #new way, JAN 7, 2020
    repliseq_wavelet_signal_df =readWig_with_fixedStep_variableStep(smoothedWaveletRepliseqDataFilename)

    # #old way starts
    # # Do not use sum, GSM923442_hg19_wgEncodeUwRepliSeqMcf7SumSignalRep1.wig contains values greater than 600
    # # Use Smoothed Wavelet Signal, GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.wig
    # unprocessed_df = readRepliSeqSignal(smoothedWaveletRepliseqDataFilename)
    #
    # #Process the signal, convert into interval version
    # repliseq_wavelet_signal_df = processSmoothedWaveletSignal(unprocessed_df)
    # #old way ends

    print('Chromosome names in replication time signal data: %s' % (repliseq_wavelet_signal_df[CHROM].unique()))
    # print('repliseq_wavelet_signal_df[chr].unique')
    # print(repliseq_wavelet_signal_df['chr'].unique())
    ################### Read the Smoothed Wavelet Replication Time Signal ends #############################


    ############## Read the Valleys and Peaks starts #######################################
    #read Valleys (local minima) bed file and read Peaks (local maxima) bed file
    # valleysBEDFilename = 'GSM923442_hg19_wgEncodeUwRepliSeqMcf7ValleysRep1.bed'
    # peaksBEDFilename = 'GSM923442_hg19_wgEncodeUwRepliSeqMcf7PkRep1.bed'

    # #old way starts
    # valleys_df= readBED(valleysBEDFilename)
    # print('Chromosome names in replication time valleys data: %s' % (valleys_df['chr'].unique()))
    # # print('valleys_df[chr].unique()')
    # # print(valleys_df['chr'].unique())
    #
    # peaks_df = readBED(peaksBEDFilename)
    # print('Chromosome names in replication time peaks data: %s' % (peaks_df['chr'].unique()))
    # # print('peaks_df[chr].unique()')
    # # print(peaks_df['chr'].unique())
    #
    # valleys_df.drop(valleys_df.columns[[3,4,5,6,7,8]], axis=1, inplace=True)
    # peaks_df.drop(peaks_df.columns[[3,4,5,6,7,8]], axis=1, inplace=True)
    # #old way ends

    #new way starts JAN 7, 2020
    discard_signal=True
    valleys_df= readFileInBEDFormat(valleysBEDFilename,discard_signal)
    valleys_df[END] = valleys_df[END] - 1
    print('Chromosome names in replication time valleys data: %s' % (valleys_df[CHROM].unique()))

    peaks_df = readFileInBEDFormat(peaksBEDFilename,discard_signal)
    peaks_df[END] = peaks_df[END] - 1
    print('Chromosome names in replication time peaks data: %s' % (peaks_df[CHROM].unique()))
    #new way ends JAN 7, 2020


    ############## Read the Valleys and Peaks ends ########################################

    return repliseq_wavelet_signal_df, valleys_df, peaks_df
########################################################################

########################################################################
def replicationStrandBiasAnalysis(computationType,sample_based,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,smoothedWaveletRepliseqDataFilename,valleysBEDFilename, peaksBEDFilename,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,verbose):

    print('\n#################################################################################')
    print('--- ReplicationStrandBias Analysis starts')

    ###############################################
    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)
    ###############################################

    ###############################################
    repliseq_signal_df, valleys_df, peaks_df = read_repliseq_dataframes(smoothedWaveletRepliseqDataFilename,valleysBEDFilename,peaksBEDFilename)
    ###############################################

    ############################Chr based parallel code starts ################################################
    #prepare the input for parallel lines starts
    strandBias = REPLICATIONSTRANDBIAS

    #Accumulate chrBased Results
    accumulatedAllChromosomesType2ReplicationStrand2CountDict = {}
    accumulatedAllChromosomesSample2Type2ReplicationStrand2CountDict = {}
    accumulatedAllChromosomesType2Sample2ReplicationStrand2CountDict = {}
    accumulatedAllChromosomesSignature2MutationType2ReplicationStrand2CountDict = {}

    if (computationType == USING_APPLY_ASYNC):
        ################################################################################
        for chrLong in chromNamesList:
            chromSize = chromSizesDict[chrLong]

            #Read chrBasedSmoothedWaveletReplicationTimeSignalDF
            chrBased_SmoothedWaveletReplicationTimeSignal_df = repliseq_signal_df[repliseq_signal_df[CHROM] == chrLong]

            chrBasedValleysDF = valleys_df[valleys_df[CHROM] == chrLong].copy()
            chrBasedValleysDF['type'] = 'Valley'
            chrBasedValleysDF.astype(dtype={START: int, END: int})

            chrBasedPeaksDF = peaks_df[peaks_df[CHROM] == chrLong].copy()
            chrBasedPeaksDF['type'] = 'Peak'
            chrBasedPeaksDF.astype(dtype={START: int, END: int})

            # Concat Peaks and Valleys
            chrBased_valleys_peaks_df = pd.concat([chrBasedValleysDF, chrBasedPeaksDF], axis=0)

            # Sort Valleys and peaks
            chrBased_valleys_peaks_df.sort_values(START, inplace=True)

            ################################################################################
            if ((chrBased_SmoothedWaveletReplicationTimeSignal_df is not None) and (not chrBased_SmoothedWaveletReplicationTimeSignal_df.empty) and (checkforValidness(chrBased_valleys_peaks_df))):
                chrBased_replication_array = fill_chr_based_replication_strand_array(chrLong,chromSize,chrBased_SmoothedWaveletReplicationTimeSignal_df,chrBased_valleys_peaks_df)

                ####################################################################
                def accumulate_apply_async_result(result_tuple):
                    chrBased_SimNum2Type2Strand2CountDict = result_tuple[0]
                    chrBased_SimNum2Sample2Type2Strand2CountDict = result_tuple[1]
                    chrBased_SimNum2Type2Sample2Strand2CountDict = result_tuple[2]
                    chrBased_SimNum2Signature2MutationType2Strand2CountDict = result_tuple[3]

                    accumulate_simulations_integrated_for_each_tuple(
                        chrBased_SimNum2Type2Strand2CountDict,
                        chrBased_SimNum2Sample2Type2Strand2CountDict,
                        chrBased_SimNum2Type2Sample2Strand2CountDict,
                        chrBased_SimNum2Signature2MutationType2Strand2CountDict,
                        accumulatedAllChromosomesType2ReplicationStrand2CountDict,
                        accumulatedAllChromosomesSample2Type2ReplicationStrand2CountDict,
                        accumulatedAllChromosomesType2Sample2ReplicationStrand2CountDict,
                        accumulatedAllChromosomesSignature2MutationType2ReplicationStrand2CountDict)
                ####################################################################

                ################################################################################
                for simNum in range(0,numofSimulations+1):
                    chrBased_simBased_subs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, SUBS, simNum)
                    chrBased_simBased_indels_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, INDELS, simNum)
                    chrBased_simBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, DINUCS, simNum)
                    pool.apply_async(searchMutationsOnReplicationStrandArrayForApplyAsync,(chrBased_replication_array,
                                                                              chrBased_simBased_subs_df,
                                                                              chrBased_simBased_indels_df,
                                                                              chrBased_simBased_dinucs_df,
                                                                              numofSimulations,
                                                                              sample_based,
                                                                              subsSignature2PropertiesListDict,
                                                                              indelsSignature2PropertiesListDict,
                                                                              dinucsSignature2PropertiesListDict),callback=accumulate_apply_async_result)
                ################################################################################

        ################################################################################


    elif (computationType == COMPUTATION_CHROMOSOMES_SEQUENTIAL_ALL_SIMULATIONS_PARALLEL):
        ################################################################################
        for chrLong in chromNamesList:
            chromSize = chromSizesDict[chrLong]

            #Read chrBasedSmoothedWaveletReplicationTimeSignalDF
            chrBased_SmoothedWaveletReplicationTimeSignal_df = repliseq_signal_df[repliseq_signal_df[CHROM] == chrLong]

            chrBasedValleysDF = valleys_df[valleys_df[CHROM] == chrLong].copy()
            chrBasedValleysDF['type'] = 'Valley'
            chrBasedValleysDF.astype(dtype={START: int, END: int})

            chrBasedPeaksDF = peaks_df[peaks_df[CHROM] == chrLong].copy()
            chrBasedPeaksDF['type'] = 'Peak'
            chrBasedPeaksDF.astype(dtype={START: int, END: int})

            # Concat Peaks and Valleys
            chrBased_valleys_peaks_df = pd.concat([chrBasedValleysDF, chrBasedPeaksDF], axis=0)

            # Sort Valleys and peaks
            chrBased_valleys_peaks_df.sort_values(START, inplace=True)

            ################################################################################
            if ((chrBased_SmoothedWaveletReplicationTimeSignal_df is not None) and (not chrBased_SmoothedWaveletReplicationTimeSignal_df.empty) and (checkforValidness(chrBased_valleys_peaks_df))):
                chrBased_replication_array = fill_chr_based_replication_strand_array(chrLong,chromSize,chrBased_SmoothedWaveletReplicationTimeSignal_df,chrBased_valleys_peaks_df)

                #For each chrom
                poolInputList = []

                ################################################################################
                for simNum in range(0,numofSimulations+1):
                    chrBased_subs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, SUBS, simNum)
                    chrBased_indels_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, INDELS, simNum)
                    chrBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, DINUCS, simNum)

                    inputList = []
                    inputList.append(chrBased_replication_array) #same for all
                    inputList.append(chrBased_subs_df) # different split each time
                    inputList.append(chrBased_indels_df)  # different split each time
                    inputList.append(chrBased_dinucs_df)  # different split each time
                    inputList.append(numofSimulations)
                    inputList.append(sample_based)
                    inputList.append(subsSignature2PropertiesListDict)
                    inputList.append(indelsSignature2PropertiesListDict)
                    inputList.append(dinucsSignature2PropertiesListDict)

                    poolInputList.append(inputList)
                ################################################################################

                listofTuples = pool.map(searchMutationsOnReplicationStrandArray, poolInputList)

                accumulate_simulations_integrated(listofTuples,
                            accumulatedAllChromosomesType2ReplicationStrand2CountDict,
                            accumulatedAllChromosomesSample2Type2ReplicationStrand2CountDict,
                            accumulatedAllChromosomesType2Sample2ReplicationStrand2CountDict,
                            accumulatedAllChromosomesSignature2MutationType2ReplicationStrand2CountDict)
            ################################################################################

    ################################
    pool.close()
    pool.join()
    ################################

    if verbose: print('ReplicationStrandBiasAnalysis Results %s starts' %(computationType))
    if verbose: print('accumulatedAllChromosomesType2ReplicationStrand2CountDict[0]')
    if verbose: print(accumulatedAllChromosomesType2ReplicationStrand2CountDict[0])
    if verbose: print('accumulatedAllChromosomesSample2Type2ReplicationStrand2CountDict[0]')
    if verbose: print(accumulatedAllChromosomesSample2Type2ReplicationStrand2CountDict[0])
    if verbose: print('accumulatedAllChromosomesType2Sample2ReplicationStrand2CountDict[0]')
    if verbose: print(accumulatedAllChromosomesType2Sample2ReplicationStrand2CountDict[0])
    if verbose: print('accumulatedAllChromosomesSignature2MutationType2ReplicationStrand2CountDict[0]')
    if verbose: print(accumulatedAllChromosomesSignature2MutationType2ReplicationStrand2CountDict[0])
    if verbose: print('ReplicationStrandBiasAnalysis Results %s ends' %(computationType))

    ############################################################################################################
    #####################################       Output starts      #############################################
    ############################################################################################################
    writeDictionary(accumulatedAllChromosomesType2ReplicationStrand2CountDict,outputDir,jobname,Type2ReplicationStrand2CountDict_Filename,strandBias,None)
    writeDictionary(accumulatedAllChromosomesSignature2MutationType2ReplicationStrand2CountDict, outputDir, jobname,Signature2MutationType2ReplicationStrand2CountDict_Filename, strandBias, None)

    if sample_based:
        writeDictionary(accumulatedAllChromosomesSample2Type2ReplicationStrand2CountDict,outputDir,jobname,Sample2Type2ReplicationStrand2CountDict_Filename,strandBias,None)
        writeDictionary(accumulatedAllChromosomesType2Sample2ReplicationStrand2CountDict,outputDir,jobname,Type2Sample2ReplicationStrand2CountDict_Filename,strandBias,None)
    ############################################################################################################
    #####################################       Output ends      ###############################################
    ############################################################################################################

    print('--- ReplicationStrandBias Analysis ends')
    print('#################################################################################\n')
########################################################################