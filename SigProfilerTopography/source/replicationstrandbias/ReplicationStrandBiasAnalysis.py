# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018 Burcak Otlu


# Version2
# This version use np.arrays
# Right now replication strand bias analysis works for single point mutations and signatures.
# This python code analyses the Replication Strand Bias


import sys
import os

#############################################################
current_abs_path = os.path.abspath(os.path.dirname(__file__))
print('ReplicationStrandBiasAnalysis.py current_abs_path:%s' %(current_abs_path))
#############################################################

commonsPath = os.path.join(current_abs_path, '..','commons')
sys.path.append(commonsPath)

from SigProfilerTopography.source.commons.TopographyCommons import *

#For Supp Fig2B
CHR10_THRESHOLD_START = 16400000
CHR10_THRESHOLD_END = 26400000

#For Supp Fig2A
CHR20_START = 36260000
CHR20_END = 36830000

# Depreceated May 1, 2018
# THRESHOLD_SIGNATURE_MUTATION_PROBABILITY = 0.6

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
    subset_df = chrBasedSmoothedWaveletReplicationTimeSignalDF[(chrBasedSmoothedWaveletReplicationTimeSignalDF['start']>=peakorValleyStart) & (chrBasedSmoothedWaveletReplicationTimeSignalDF['end']<=peakorValleyEnd)]

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
            slope = (row.get('signal') - formerRow.get('signal')) / 1000
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
                    transitionZoneList.append((chrLong,start,(row.get('start') + row.get('end'))//2,formerSlopeDirection,consecutiveLength))
                #initialize and start again
                consecutiveLength = 1000
                start = (row.get('start') + row.get('end'))//2
                formerRow= row
                formerSlopeDirection= np.sign(slope)
                continue

            # print('slope: %f - np.sign(slope): %f -  consecutiveLength: %d ' %(slope,np.sign(slope),consecutiveLength))
            formerSlopeDirection = np.sign(slope)

    #This is for the last probable transition zone.
    if (consecutiveLength >= THRESHOLD_CONSECUTIVE_LONG_STRETCH_LENGTH):
        # print('After for loop ends, found one: from %d to %s with %d bases with slope sign %s' % (start, (row.get('start') + row.get('end'))//2, consecutiveLength, formerSlopeDirection))
        transitionZoneList.append((chrLong,start,(row.get('start') + row.get('end'))//2,formerSlopeDirection,consecutiveLength))

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
        peakorValleyStart = row['start']
        peakorValleyEnd = row['end']
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
            valleyStart =row['start']
            valleyEnd = row['end']
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
def searchMutationOnReplicationStrandArray(mutation_row,
                                            chrBasedReplicationArray,
                                            mutationType2ReplicationStrand2CountDict,
                                            mutationType2Sample2ReplicationStrand2CountDict,
                                            mutationProbability2Signature2ReplicationStrand2CountDict,
                                            mutationProbability2Signature2Sample2ReplicationStrand2CountDict,
                                            signatureList,
                                            mutationProbabilityList):


    mutationStart = mutation_row[START]
    mutationEnd = mutation_row[END]
    mutationPyramidineStrand = mutation_row[PYRAMIDINESTRAND]
    mutationType = mutation_row[MUTATION]
    mutationSample = mutation_row[SAMPLE]

    #############################################################################################################
    #if there is overlap with chrBasedReplicationArray
    slicedArray = chrBasedReplicationArray[mutationStart:(mutationEnd+1)]

    if (np.any(slicedArray)):
        uniqueValueArray = np.unique(slicedArray[np.nonzero(slicedArray)])

        # I expect the value of 1 (LEADING on the positive strand) or -1 (LAGGING on the positive strand) so size must be one.
        if (uniqueValueArray.size == 1):
            for uniqueValue in np.nditer(uniqueValueArray):
                # type(decileIndex) is numpy.ndarray
                slope = int(uniqueValue)

                #They have the same sign, multiplication (1,1) (-1,-1) must be 1
                if (slope*mutationPyramidineStrand > 0):
                    updateDictionaries(mutation_row,
                                       mutationType,
                                       mutationSample,
                                       mutationType2ReplicationStrand2CountDict,
                                       mutationType2Sample2ReplicationStrand2CountDict,
                                       mutationProbability2Signature2ReplicationStrand2CountDict,
                                       mutationProbability2Signature2Sample2ReplicationStrand2CountDict,
                                       LEADING,
                                       signatureList,
                                       mutationProbabilityList)

                # They have the opposite sign, multiplication(1,-1) (-1,-)  must be -1
                elif (slope*mutationPyramidineStrand < 0):
                    updateDictionaries(mutation_row,
                                       mutationType,
                                       mutationSample,
                                       mutationType2ReplicationStrand2CountDict,
                                       mutationType2Sample2ReplicationStrand2CountDict,
                                       mutationProbability2Signature2ReplicationStrand2CountDict,
                                       mutationProbability2Signature2Sample2ReplicationStrand2CountDict,
                                       LAGGING,
                                       signatureList,
                                       mutationProbabilityList)
        else:
            print('There is a situation!!!')
    #############################################################################################################

########################################################################


########################################################################
def  searchMutationsOnReplicationArray(inputList):
    chrBased_replication_array = inputList[0]
    chrBasedSPMsSplitDF = inputList[1]
    signatureList = inputList[2]
    mutationProbabilityList = inputList[3]

    chrBasedMutationType2ReplicationStrand2CountDict = {}
    chrBasedMutationType2Sample2ReplicationStrand2CountDict = {}
    mutationProbability2Signature2ReplicationStrand2CountDict = {}
    mutationProbability2Signature2Sample2ReplicationStrand2CountDict = {}

    ##############################  Fill the type2replicationStranCount dictionaries  starts ########
    if ((chrBasedSPMsSplitDF is not None) and (not chrBasedSPMsSplitDF.empty)):
        chrBasedSPMsSplitDF.apply(searchMutationOnReplicationStrandArray,
                                chrBasedReplicationArray=chrBased_replication_array,
                                mutationType2ReplicationStrand2CountDict=chrBasedMutationType2ReplicationStrand2CountDict,
                                mutationType2Sample2ReplicationStrand2CountDict=chrBasedMutationType2Sample2ReplicationStrand2CountDict,
                                mutationProbability2Signature2ReplicationStrand2CountDict=mutationProbability2Signature2ReplicationStrand2CountDict,
                                mutationProbability2Signature2Sample2ReplicationStrand2CountDict=mutationProbability2Signature2Sample2ReplicationStrand2CountDict,
                                signatureList=signatureList,
                                mutationProbabilityList=mutationProbabilityList,
                                axis=1)
    ##############################  Fill the type2replicationStranCount dictionaries  ends ##########

    #Fill the type2replicationStranCount dictionaries and return them
    return (chrBasedMutationType2ReplicationStrand2CountDict,
            chrBasedMutationType2Sample2ReplicationStrand2CountDict,
            mutationProbability2Signature2ReplicationStrand2CountDict,
            mutationProbability2Signature2Sample2ReplicationStrand2CountDict)
########################################################################

########################################################################
def fillReplicationStrandArrayAndSearchMutationsOnThisArray(inputList):

    chrLong = inputList[0]
    chrBasedSmoothedWaveletReplicationTimeSignalDF = inputList[1]
    chrBasedValleysPeaksDF = inputList[2]
    chrBasedSPMsDF = inputList[3]
    signatureList = inputList[4]
    mutationProbabilityList = inputList[5]

    chrBasedMutationType2ReplicationStrand2CountDict = {}
    chrBasedMutationType2Sample2ReplicationStrand2CountDict = {}
    mutationProbability2Signature2ReplicationStrand2CountDict = {}
    mutationProbability2Signature2Sample2ReplicationStrand2CountDict = {}

    # +1 means leading strand, -1 means lagging strand
    chrBased_replication_array = np.zeros(MAXIMUM_CHROMOSOME_LENGTH, dtype=np.int8)

    ##############################  Fill the replication array starts ##############################
    if (chrBasedSmoothedWaveletReplicationTimeSignalDF is not None and not chrBasedSmoothedWaveletReplicationTimeSignalDF.empty):

        firstIndex = chrBasedSmoothedWaveletReplicationTimeSignalDF.index[0]
        lastIndex = chrBasedSmoothedWaveletReplicationTimeSignalDF.index[-1]

        start = chrBasedSmoothedWaveletReplicationTimeSignalDF.ix[firstIndex,'start']  # get the first row start
        end = chrBasedSmoothedWaveletReplicationTimeSignalDF.ix[lastIndex,'end']  # get the last row end

        #Step1 Find the transition zones
        chrBasedTransitionZonesList = findLongStretchesofConsistentTransitionZones(chrLong,
                                                                                   start,
                                                                                   end,
                                                                                   chrBasedSmoothedWaveletReplicationTimeSignalDF,
                                                                                   chrBasedValleysPeaksDF)

        labels = ['chr', 'start', 'end', 'slopeDirection', 'length']
        chrBasedTransitionZonesDF = pd.DataFrame.from_records(chrBasedTransitionZonesList, columns=labels)


        #Step2 Fill the replication array using transition zones
        chrBasedTransitionZonesDF.apply(fillReplicationStrandArray,chrBased_replication_array=chrBased_replication_array,axis=1)
    ##############################  Fill the replication array ends ################################

    ##############################  Fill the type2replicationStranCount dictionaries  starts ########
    if (chrBasedSPMsDF is not None and (not chrBasedSPMsDF.empty)):
        chrBasedSPMsDF.apply(searchMutationOnReplicationStrandArray,
                                chrBasedReplicationArray=chrBased_replication_array,
                                mutationType2ReplicationStrand2CountDict=chrBasedMutationType2ReplicationStrand2CountDict,
                                mutationType2Sample2ReplicationStrand2CountDict=chrBasedMutationType2Sample2ReplicationStrand2CountDict,
                                mutationProbability2Signature2ReplicationStrand2CountDict=mutationProbability2Signature2ReplicationStrand2CountDict,
                                mutationProbability2Signature2Sample2ReplicationStrand2CountDict=mutationProbability2Signature2Sample2ReplicationStrand2CountDict,
                                signatureList=signatureList,
                                mutationProbabilityList=mutationProbabilityList,
                                axis=1)
    ##############################  Fill the type2replicationStranCount dictionaries  ends ##########

    #Fill the type2replicationStranCount dictionaries and return them
    return (chrBasedMutationType2ReplicationStrand2CountDict,
            chrBasedMutationType2Sample2ReplicationStrand2CountDict,
            mutationProbability2Signature2ReplicationStrand2CountDict,
            mutationProbability2Signature2Sample2ReplicationStrand2CountDict)

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
def replicationStrandBiasAnalysis(outputDir,jobname,singlePointMutationsFilename,smoothedWaveletRepliseqDataFilename,valleysBEDFilename, peaksBEDFilename,startMutationProbability,endMutationProbability,step):

    print('########################## ReplicationStrandBias Analysis starts ##########################')
    print('#################### ReplicationStrandBias Analysis system arguments: #####################')
    print(sys.argv)

    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)

    ############################################
    # jobname = sys.argv[1]
    # singlePointMutationsFilename = sys.argv[2]

    # smoothedWaveletRepliseqDataFilename = sys.argv[3]
    # valleysBEDFilename = sys.argv[4]
    # peaksBEDFilename = sys.argv[5]

    # startMutationProbability = float(sys.argv[6])
    # endMutationProbability = float(sys.argv[7])
    # step = float(sys.argv[8])
    ############################################

    #Load signatureList
    signatureList = []
    SignaturesFilePath = os.path.join(outputDir,jobname,DATA,SignatureFilename)
    if (os.path.exists(SignaturesFilePath)):
        signaturesArray = np.loadtxt(SignaturesFilePath,dtype=str, delimiter='\t')
        signatureList = list(signaturesArray)

    #Load the chrnames in single point mutations data
    ChrNamesFile = os.path.join(outputDir,jobname,DATA,ChrNamesInSPMsFilename)
    if (os.path.exists(ChrNamesFile)):
        chrNamesArray = np.loadtxt(ChrNamesFile, dtype=str, delimiter='\t')
        chrNamesInSPMs = chrNamesArray.tolist()

    ################### Read the Smoothed Wavelet Replication Time Signal starts ###########################
    # Do not use sum, GSM923442_hg19_wgEncodeUwRepliSeqMcf7SumSignalRep1.wig contains values greater than 600
    # Use Smoothed Wavelet Signal, GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.wig
    unprocessed_df = readRepliSeqSignal(smoothedWaveletRepliseqDataFilename)

    #Process the signal, convert into interval version
    processed_df = processSmoothedWaveletSignal(unprocessed_df)
    print('repliseq_wavelet_signal_df[chr].unique')
    print(processed_df['chr'].unique())
    # print('processed_df.shape interval version')
    # print(processed_df.shape)

    # print('processed_df.info()')
    # print(processed_df.info())
    ################### Read the Smoothed Wavelet Replication Time Signal ends #############################


    ############## Read the Valleys and Peaks starts #######################################
    #read Valleys (local minima) bed file and read Peaks (local maxima) bed file
    # valleysBEDFilename = 'GSM923442_hg19_wgEncodeUwRepliSeqMcf7ValleysRep1.bed'
    # peaksBEDFilename = 'GSM923442_hg19_wgEncodeUwRepliSeqMcf7PkRep1.bed'

    valleys_df= readBED(valleysBEDFilename)
    print('valleys_df[chr].unique()')
    print(valleys_df['chr'].unique())

    peaks_df = readBED(peaksBEDFilename)
    print('peaks_df[chr].unique()')
    print(peaks_df['chr'].unique())

    # print('debug starts')
    # print('-------------------------valleys_df.dtypes------------------------')
    # print(valleys_df.dtypes)
    # print('-------------------------peaks_df.dtypes------------------------')
    # print(peaks_df.dtypes)
    # print('debug ends')

    valleys_df.drop(valleys_df.columns[[3,4,5,6,7,8]], axis=1, inplace=True)
    peaks_df.drop(peaks_df.columns[[3,4,5,6,7,8]], axis=1, inplace=True)
    ############## Read the Valleys and Peaks ends ########################################


    #################### Prepare mutation Probability List starts ################################################
    mutationProbabilityList = prepareMutationProbabilityList(startMutationProbability, endMutationProbability, step)
    #################### Prepare mutation Probability List ends ##################################################

    ################################### Part3 starts ##################################
    ################## ncomms11383 Fig2b starts #######################################
    ######################### MutationType Based starts ###############################
    ######################### Signature Based starts ##################################
    ### In this plot we will analyze signature based the replication strand bias
    # Get the sample based number of mutations that has probability greater than equal to the 0.5

    # get Leading and Lagging Intervals
    # create an interval tree from these intervals
    # get for each mutation whether if it overlaps with leading or lagging replication strand.
    # update the mutation type and signature accordingly

    ############################Chr based parallel code starts ################################################
    #prepare the input for parallel lines starts
    replicationStrands = [LAGGING, LEADING]
    strandBias = REPLICATIONSTRANDBIAS

    #Accumulate chrBased Results
    accumulatedAllChromosomesMutationType2LeadingLaggingStrand2CountDict = {}
    accumulatedAllChromosomesMutationType2Sample2LeadingLaggingStrand2CountDict = {}
    accumulatedAllChromosomesMutationProbability2Signature2LeadingLaggingStrand2CountDict = {}
    accumulatedAllChromosomesMutationProbability2Signature2Sample2LeadingLaggingStrand2CountDict = {}

    #For more information
    accumulatedAllChromosomesMutationProbability2Signature2RatioDict = {}

    # ############################################################################################################
    # #####################################      Version 1 starts      ###########################################
    # ###############################      All Chromosomes in parallel      ######################################
    # ############################################################################################################
    # poolInputList = []
    # for chrShort in chrNamesInSPMs:
    #     chrLong = 'chr' + chrShort
    #
    #     #Read chrBased spms dataframe
    #     chrBased_spms_df = readChrBasedMutationDF(jobname, chrLong, singlePointMutationsFilename)
    #
    #     chrBasedSmoothedWaveletReplicationTimeSignalDF = processed_df[processed_df['chr'] == chrLong]
    #
    #     chrBasedValleysDF = valleys_df[valleys_df['chr'] == chrLong].copy()
    #     chrBasedValleysDF['type'] = 'Valley'
    #     chrBasedValleysDF.astype(dtype={'start': int,'end': int})
    #
    #     chrBasedPeaksDF = peaks_df[peaks_df['chr'] == chrLong].copy()
    #     chrBasedPeaksDF['type'] = 'Peak'
    #     chrBasedPeaksDF.astype(dtype={'start': int, 'end': int})
    #
    #     # Concat Peaks and Valleys
    #     chrBased_valleys_peaks_df = pd.concat([chrBasedValleysDF, chrBasedPeaksDF], axis=0)
    #
    #     # Sort Valleys and peaks
    #     chrBased_valleys_peaks_df.sort_values('start', inplace=True)
    #
    #     inputList = []
    #
    #     #After sorting check for consecutive rows always differ in terms of their types (peak is followed by valley or vice versa)
    #     if checkforValidness(chrBased_valleys_peaks_df):
    #         inputList.append(chrLong)
    #         inputList.append(chrBasedSmoothedWaveletReplicationTimeSignalDF)
    #         inputList.append(chrBased_valleys_peaks_df)
    #         inputList.append(chrBased_spms_df)
    #         inputList.append(signatureList) # same for all
    #         inputList.append(mutationProbabilityList) # same for all
    #         poolInputList.append(inputList)
    # #prepare the input for parallel lines ends
    #
    # listofTuples = pool.map(fillReplicationStrandArrayAndSearchMutationsOnThisArray, poolInputList)
    #
    # accumulate(listofTuples,
    #            accumulatedAllChromosomesMutationType2LeadingLaggingStrand2CountDict,
    #            accumulatedAllChromosomesMutationType2Sample2LeadingLaggingStrand2CountDict,
    #            accumulatedAllChromosomesMutationProbability2Signature2LeadingLaggingStrand2CountDict,
    #            accumulatedAllChromosomesMutationProbability2Signature2Sample2LeadingLaggingStrand2CountDict)
    # ############################################################################################################
    # #####################################      Version 1 ends      #############################################
    # ###############################      All Chromosomes in parallel      ######################################
    # ############################################################################################################


    ############################################################################################################
    #####################################      Version2  starts      ###########################################
    ###############################       Chromosomes sequentially      ########################################
    ###############################      All ChrBased Splits in parallel     ###################################
    ############################################################################################################
    for chrShort in chrNamesInSPMs:
        chrLong = 'chr' + chrShort

        # Read chrBased spms dataframe
        chrBased_spms_df = readChrBasedMutationDF(outputDir,jobname, chrLong, singlePointMutationsFilename)

        #Read chrBasedSmoothedWaveletReplicationTimeSignalDF
        chrBasedSmoothedWaveletReplicationTimeSignalDF = processed_df[processed_df['chr'] == chrLong]

        # +1 means leading strand, -1 means lagging strand
        # we will fill this array using smoothedSignal, peaks and valleys for each chromosome
        chrBased_replication_array = np.zeros(MAXIMUM_CHROMOSOME_LENGTH, dtype=np.int8)

        chrBasedValleysDF = valleys_df[valleys_df['chr'] == chrLong].copy()
        chrBasedValleysDF['type'] = 'Valley'
        chrBasedValleysDF.astype(dtype={'start': int, 'end': int})

        chrBasedPeaksDF = peaks_df[peaks_df['chr'] == chrLong].copy()
        chrBasedPeaksDF['type'] = 'Peak'
        chrBasedPeaksDF.astype(dtype={'start': int, 'end': int})

        # Concat Peaks and Valleys
        chrBased_valleys_peaks_df = pd.concat([chrBasedValleysDF, chrBasedPeaksDF], axis=0)

        # Sort Valleys and peaks
        chrBased_valleys_peaks_df.sort_values('start', inplace=True)

        if ((chrBasedSmoothedWaveletReplicationTimeSignalDF is not None) and (not chrBasedSmoothedWaveletReplicationTimeSignalDF.empty) and (checkforValidness(chrBased_valleys_peaks_df))):
            firstIndex = chrBasedSmoothedWaveletReplicationTimeSignalDF.index[0]
            lastIndex = chrBasedSmoothedWaveletReplicationTimeSignalDF.index[-1]

            startColumnIndex = chrBasedSmoothedWaveletReplicationTimeSignalDF.columns.get_loc('start')
            endColumnIndex = chrBasedSmoothedWaveletReplicationTimeSignalDF.columns.get_loc('end')

            start = chrBasedSmoothedWaveletReplicationTimeSignalDF.iloc[0, startColumnIndex]  # get the first row start
            end = chrBasedSmoothedWaveletReplicationTimeSignalDF.iloc[-1, endColumnIndex]  # get the last row end

            # Step1 Find the transition zones
            chrBasedTransitionZonesList = findLongStretchesofConsistentTransitionZones(chrLong,
                                                                                       start,
                                                                                       end,
                                                                                       chrBasedSmoothedWaveletReplicationTimeSignalDF,
                                                                                       chrBased_valleys_peaks_df)

            labels = ['chr', 'start', 'end', 'slopeDirection', 'length']
            chrBasedTransitionZonesDF = pd.DataFrame.from_records(chrBasedTransitionZonesList, columns=labels)

            # Step2 Fill the replication array using transition zones
            chrBasedTransitionZonesDF.apply(fillReplicationStrandArray,chrBased_replication_array=chrBased_replication_array, axis=1)


            if ((chrBased_spms_df is not None) and (not chrBased_spms_df.empty)):
                chrBasedSPMsDFSplits = np.array_split(chrBased_spms_df, numofProcesses)

                poolInputList = []

                for chrBasedSPMsDFSplit in chrBasedSPMsDFSplits:
                    inputList = []
                    inputList.append(chrBased_replication_array) #same for all
                    inputList.append(chrBasedSPMsDFSplit) # different split each time
                    inputList.append(signatureList)  # same for all
                    inputList.append(mutationProbabilityList)  # same for all
                    poolInputList.append(inputList)

                listofTuples = pool.map(searchMutationsOnReplicationArray, poolInputList)

                accumulate(listofTuples,
                           accumulatedAllChromosomesMutationType2LeadingLaggingStrand2CountDict,
                           accumulatedAllChromosomesMutationType2Sample2LeadingLaggingStrand2CountDict,
                           accumulatedAllChromosomesMutationProbability2Signature2LeadingLaggingStrand2CountDict,
                           accumulatedAllChromosomesMutationProbability2Signature2Sample2LeadingLaggingStrand2CountDict)
    ############################################################################################################
    #####################################      Version2  ends      #############################################
    ###############################       Chromosomes sequentially      ########################################
    ###############################      All ChrBased Splits in parallel     ###################################
    ############################################################################################################

    ######################### MutationType Based ends ###############################
    ######################### Signature Based ends ##################################
    ################## ncomms11383 Fig2A and Fig2B ends #############################
    ################################### Part3 ends ##################################




    ############################################################################################################
    #####################################       Output starts      #############################################
    ############################################################################################################
    calculateRatio(accumulatedAllChromosomesMutationProbability2Signature2LeadingLaggingStrand2CountDict,
                   accumulatedAllChromosomesMutationProbability2Signature2RatioDict,
                   replicationStrands)

    #There is no order in keys
    # print('Keys are not sorted.')
    print('accumulatedAllChromosomesMutationType2LeadingLaggingStrandCountDict')
    print(accumulatedAllChromosomesMutationType2LeadingLaggingStrand2CountDict)

    #############################################################################
    # Highligh the results for mutation probability threshold
    print('###############################################################')


    print('For mutation probability threshold: %f' %MUTATION_SIGNATURE_PROBABILITY_THRESHOLD)
    print('Signature\t%s\t%s' % (LAGGING,LEADING))

    print('For debug March 13, 2019')
    print('accumulatedAllChromosomesMutationProbability2Signature2LeadingLaggingStrand2CountDict')
    print(accumulatedAllChromosomesMutationProbability2Signature2LeadingLaggingStrand2CountDict)

    signature2ReplicationStrand2CountDict = accumulatedAllChromosomesMutationProbability2Signature2LeadingLaggingStrand2CountDict[MUTATION_SIGNATURE_PROBABILITY_THRESHOLD]

    for signature in signatureList:
        if (signature in signature2ReplicationStrand2CountDict.keys()):
            replicationStrand2CountDict = signature2ReplicationStrand2CountDict[signature]
            if ((LAGGING in replicationStrand2CountDict.keys()) and (LEADING in replicationStrand2CountDict.keys())):
                print('%s\t%d\t%d' %(signature, replicationStrand2CountDict[LAGGING], replicationStrand2CountDict[LEADING]))
    print('ReplicationStrandBias ratio: number of mutations on lagging/(lagging + leading)')
    print('###############################################################')
    #############################################################################

    #Convert
    signature2WeightedAverageRatioDict, signature2StdErrorDict, signature2SumofMutationProbabilitiesDict = convert(accumulatedAllChromosomesMutationProbability2Signature2LeadingLaggingStrand2CountDict, replicationStrands)

    ##############################################################################
    #To be used for plotting starts
    writeDictionary(accumulatedAllChromosomesMutationType2LeadingLaggingStrand2CountDict,outputDir,jobname,MutationType2ReplicationStrand2CountDict_Filename,strandBias,None)
    writeDictionary(accumulatedAllChromosomesMutationType2Sample2LeadingLaggingStrand2CountDict,outputDir,jobname,MutationType2Sample2ReplicationStrand2CountDict_Filename,strandBias,None)
    writeDictionary(accumulatedAllChromosomesMutationProbability2Signature2LeadingLaggingStrand2CountDict,outputDir,jobname,MutationProbability2Signature2ReplicationStrand2CountDict_Filename,strandBias,None)
    writeDictionary(accumulatedAllChromosomesMutationProbability2Signature2Sample2LeadingLaggingStrand2CountDict,outputDir,jobname,MutationProbability2Signature2Sample2ReplicationStrand2CountDict_Filename,strandBias,None)

    writeDictionary(signature2WeightedAverageRatioDict,outputDir,jobname,Signature2ReplicationWeightedAverageRatioDict_Filename,strandBias,None)
    writeDictionary(signature2StdErrorDict,outputDir,jobname,Signature2ReplicationStdErrorDict_Filename,strandBias,None)
    writeDictionary(signature2SumofMutationProbabilitiesDict,outputDir,jobname, Signature2ReplicationSumofMutationProbabilitiesDict_Filename,strandBias,None)
    #To be used for plotting ends
    ##############################################################################

    ############################################################################################################
    #####################################       Output ends      ###############################################
    ############################################################################################################

    print('########################## ReplicationStrandBias Analysis ends ############################')

########################################################################


