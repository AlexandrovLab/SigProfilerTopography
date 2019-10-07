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

from SigProfilerTopography.source.commons.TopographyCommons import *

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
#legacy code
# Summary:
#   if mutationPyramidineStrand and slope have the same sign increase LEADING STRAND count
#   else mutationPyramidineStrand and slope have the opposite sign increase LAGGING STRAND count
def searchMutationOnReplicationStrandArray(
        mutation_row,
        chrBasedReplicationArray,
        type2ReplicationStrand2CountDict,
        sample2Type2ReplicationStrand2CountDict,
        type2Sample2ReplicationStrand2CountDict,
        signature2MutationType2ReplicationStrand2CountDict,
        signature2NumberofMutationsDict,
        mutationProbabilityThreshold,
        type):

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
            updateDictionaries(mutation_row,
                                        mutationType,
                                        sample,
                                        type2ReplicationStrand2CountDict,
                                        sample2Type2ReplicationStrand2CountDict,
                                        type2Sample2ReplicationStrand2CountDict,
                                        signature2MutationType2ReplicationStrand2CountDict,
                                        LAGGING,
                                        signature2NumberofMutationsDict,
                                        mutationProbabilityThreshold)

            updateDictionaries(mutation_row,
                                        mutationType,
                                        sample,
                                        type2ReplicationStrand2CountDict,
                                        sample2Type2ReplicationStrand2CountDict,
                                        type2Sample2ReplicationStrand2CountDict,
                                        signature2MutationType2ReplicationStrand2CountDict,
                                        LEADING,
                                        signature2NumberofMutationsDict,
                                        mutationProbabilityThreshold)


        # I expect the value of 1 (LEADING on the positive strand) or -1 (LAGGING on the positive strand) so size must be one.
        elif (uniqueValueArray.size == 1):
            for uniqueValue in np.nditer(uniqueValueArray):
                # type(decileIndex) is numpy.ndarray
                slope = int(uniqueValue)

                #They have the same sign, multiplication (1,1) (-1,-1) must be 1
                if (slope*pyramidineStrand > 0):
                    updateDictionaries(mutation_row,
                                            mutationType,
                                            sample,
                                            type2ReplicationStrand2CountDict,
                                            sample2Type2ReplicationStrand2CountDict,
                                            type2Sample2ReplicationStrand2CountDict,
                                            signature2MutationType2ReplicationStrand2CountDict,
                                            LEADING,
                                            signature2NumberofMutationsDict,
                                            mutationProbabilityThreshold)

                # They have the opposite sign, multiplication(1,-1) (-1,-)  must be -1
                elif (slope*pyramidineStrand < 0):
                    updateDictionaries(mutation_row,
                                            mutationType,
                                            sample,
                                            type2ReplicationStrand2CountDict,
                                            sample2Type2ReplicationStrand2CountDict,
                                            type2Sample2ReplicationStrand2CountDict,
                                            signature2MutationType2ReplicationStrand2CountDict,
                                            LAGGING,
                                            signature2NumberofMutationsDict,
                                            mutationProbabilityThreshold)
        else:
            print('There is a situation!!!')
    #############################################################################################################
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

    startColumnIndex = chrBasedSmoothedWaveletReplicationTimeSignalDF.columns.get_loc('start')
    endColumnIndex = chrBasedSmoothedWaveletReplicationTimeSignalDF.columns.get_loc('end')

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
    # Do not use sum, GSM923442_hg19_wgEncodeUwRepliSeqMcf7SumSignalRep1.wig contains values greater than 600
    # Use Smoothed Wavelet Signal, GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.wig
    unprocessed_df = readRepliSeqSignal(smoothedWaveletRepliseqDataFilename)

    #Process the signal, convert into interval version
    repliseq_wavelet_signal_df = processSmoothedWaveletSignal(unprocessed_df)
    print('Chromosome names in replication time signal data: %s' % (repliseq_wavelet_signal_df['chr'].unique()))
    # print('repliseq_wavelet_signal_df[chr].unique')
    # print(repliseq_wavelet_signal_df['chr'].unique())
    ################### Read the Smoothed Wavelet Replication Time Signal ends #############################


    ############## Read the Valleys and Peaks starts #######################################
    #read Valleys (local minima) bed file and read Peaks (local maxima) bed file
    # valleysBEDFilename = 'GSM923442_hg19_wgEncodeUwRepliSeqMcf7ValleysRep1.bed'
    # peaksBEDFilename = 'GSM923442_hg19_wgEncodeUwRepliSeqMcf7PkRep1.bed'
    valleys_df= readBED(valleysBEDFilename)
    print('Chromosome names in replication time valleys data: %s' % (valleys_df['chr'].unique()))
    # print('valleys_df[chr].unique()')
    # print(valleys_df['chr'].unique())

    peaks_df = readBED(peaksBEDFilename)
    print('Chromosome names in replication time peaks data: %s' % (peaks_df['chr'].unique()))
    # print('peaks_df[chr].unique()')
    # print(peaks_df['chr'].unique())

    valleys_df.drop(valleys_df.columns[[3,4,5,6,7,8]], axis=1, inplace=True)
    peaks_df.drop(peaks_df.columns[[3,4,5,6,7,8]], axis=1, inplace=True)
    ############## Read the Valleys and Peaks ends ########################################

    return repliseq_wavelet_signal_df, valleys_df, peaks_df
########################################################################

########################################################################
def replicationStrandBiasAnalysis(computationType,sample_based,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,smoothedWaveletRepliseqDataFilename,valleysBEDFilename, peaksBEDFilename,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict):

    print('\n#################################################################################')
    print('--- ReplicationStrandBias Analysis starts')
    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)

    subsSignature2NumberofMutationsDict = getSubsSignature2NumberofMutationsDict(outputDir,jobname)
    indelsSignature2NumberofMutationsDict = getIndelsSignature2NumberofMutationsDict(outputDir,jobname)
    dinucsSignature2NumberofMutationsDict = getDictionary(outputDir,jobname,DinucsSignature2NumberofMutationsDictFilename)

    repliseq_signal_df, valleys_df, peaks_df = read_repliseq_dataframes(smoothedWaveletRepliseqDataFilename,valleysBEDFilename,peaksBEDFilename)

    ############################Chr based parallel code starts ################################################
    #prepare the input for parallel lines starts
    replicationStrands = [LAGGING, LEADING]
    strandBias = REPLICATIONSTRANDBIAS

    #Accumulate chrBased Results
    accumulatedAllChromosomesType2ReplicationStrand2CountDict = {}
    accumulatedAllChromosomesSample2Type2ReplicationStrand2CountDict = {}
    accumulatedAllChromosomesType2Sample2ReplicationStrand2CountDict = {}
    accumulatedAllChromosomesSignature2MutationType2ReplicationStrand2CountDict = {}

    if (computationType == COMPUTATION_ALL_CHROMOSOMES_PARALLEL):
        ############################################################################################################
        ###############################    All Chromosomes in parallel starts    ###################################
        ############################################################################################################
        poolInputList = []
        for chrLong in chromNamesList:
            chromSize = chromSizesDict[chrLong]

            original_chrBased_subs_df = readChrBasedSubsDF(outputDir, jobname, chrLong, SUBS, 0)
            original_chrBased_indels_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, INDELS, 0)
            original_chrBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, DINUCS, 0)

            combined_chrBased_subs_df = getCombinedChrBasedDF(outputDir, jobname, chrLong, original_chrBased_subs_df,SUBS, numofSimulations)
            combined_chrBased_indels_df = getCombinedChrBasedDF(outputDir, jobname, chrLong,original_chrBased_indels_df, INDELS, numofSimulations)
            combined_chrBased_dinucs_df = getCombinedChrBasedDF(outputDir, jobname, chrLong,original_chrBased_dinucs_df, DINUCS, numofSimulations)

            chrBased_SmoothedWaveletReplicationTimeSignal_df = repliseq_signal_df[repliseq_signal_df['chr'] == chrLong]

            chrBasedValleysDF = valleys_df[valleys_df['chr'] == chrLong].copy()
            chrBasedValleysDF['type'] = 'Valley'
            chrBasedValleysDF.astype(dtype={'start': int,'end': int})

            chrBasedPeaksDF = peaks_df[peaks_df['chr'] == chrLong].copy()
            chrBasedPeaksDF['type'] = 'Peak'
            chrBasedPeaksDF.astype(dtype={'start': int, 'end': int})

            # Concat Peaks and Valleys
            chrBased_valleys_peaks_df = pd.concat([chrBasedValleysDF, chrBasedPeaksDF], axis=0)

            # Sort Valleys and peaks
            chrBased_valleys_peaks_df.sort_values('start', inplace=True)

            if ((chrBased_SmoothedWaveletReplicationTimeSignal_df is not None) and (not chrBased_SmoothedWaveletReplicationTimeSignal_df.empty) and (checkforValidness(chrBased_valleys_peaks_df))):
                chrBased_replication_array = fill_chr_based_replication_strand_array(chrLong,chromSize,chrBased_SmoothedWaveletReplicationTimeSignal_df,chrBased_valleys_peaks_df)

                inputList = []
                inputList.append(chrBased_replication_array)  # same for all
                inputList.append(combined_chrBased_subs_df)  # different split each time
                inputList.append(combined_chrBased_indels_df)  # different split each time
                inputList.append(combined_chrBased_dinucs_df)  # different split each time
                inputList.append(subsSignature2NumberofMutationsDict)  # same for all
                inputList.append(indelsSignature2NumberofMutationsDict)
                inputList.append(dinucsSignature2NumberofMutationsDict)
                inputList.append(numofSimulations)
                poolInputList.append(inputList)

        listofTuples = pool.map(searchMutationsOnReplicationStrandArray, poolInputList)

        accumulate_simulations_integrated(listofTuples,
                        accumulatedAllChromosomesType2ReplicationStrand2CountDict,
                        accumulatedAllChromosomesSample2Type2ReplicationStrand2CountDict,
                        accumulatedAllChromosomesType2Sample2ReplicationStrand2CountDict,
                        accumulatedAllChromosomesSignature2MutationType2ReplicationStrand2CountDict)


        ############################################################################################################
        ###############################      All Chromosomes in parallel  ends    ##################################
        ############################################################################################################

    elif (computationType == COMPUTATION_CHROMOSOMES_SEQUENTIAL_SIMULATIONS_SEQUENTIAL):

        ################################################################################
        for chrLong in chromNamesList:
            chromSize = chromSizesDict[chrLong]

            # Read chrBasedSmoothedWaveletReplicationTimeSignalDF
            chrBased_SmoothedWaveletReplicationTimeSignal_df = repliseq_signal_df[repliseq_signal_df['chr'] == chrLong]

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

            ################################################################################
            if ((chrBased_SmoothedWaveletReplicationTimeSignal_df is not None) and (not chrBased_SmoothedWaveletReplicationTimeSignal_df.empty) and (checkforValidness(chrBased_valleys_peaks_df))):
                chrBased_replication_array = fill_chr_based_replication_strand_array(chrLong, chromSize,
                                                                                     chrBased_SmoothedWaveletReplicationTimeSignal_df,
                                                                                     chrBased_valleys_peaks_df)


                ################################################################################
                for simNum in range(0, numofSimulations + 1):
                    chrBased_subs_df = readChrBasedSubsDF(outputDir, jobname, chrLong, SUBS, simNum)
                    chrBased_indels_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, INDELS, simNum)
                    chrBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, DINUCS, simNum)

                    inputList = []
                    inputList.append(chrBased_replication_array)  # same for all
                    inputList.append(chrBased_subs_df)  # different split each time
                    inputList.append(chrBased_indels_df)  # different split each time
                    inputList.append(chrBased_dinucs_df)  # different split each time
                    inputList.append(subsSignature2NumberofMutationsDict)  # same for all
                    inputList.append(indelsSignature2NumberofMutationsDict)
                    inputList.append(dinucsSignature2NumberofMutationsDict)
                    inputList.append(numofSimulations)
                    simBased_tuple = searchMutationsOnReplicationStrandArray(inputList)
                    tupleList = []
                    tupleList.append(simBased_tuple)
                    accumulate_simulations_integrated(tupleList,
                                                      accumulatedAllChromosomesType2ReplicationStrand2CountDict,
                                                      accumulatedAllChromosomesSample2Type2ReplicationStrand2CountDict,
                                                      accumulatedAllChromosomesType2Sample2ReplicationStrand2CountDict,
                                                      accumulatedAllChromosomesSignature2MutationType2ReplicationStrand2CountDict)
                ################################################################################
            ################################################################################

    elif (computationType == COMPUTATION_CHROMOSOMES_SEQUENTIAL_ALL_SIMULATIONS_PARALLEL):
        ################################################################################
        for chrLong in chromNamesList:
            chromSize = chromSizesDict[chrLong]

            #Read chrBasedSmoothedWaveletReplicationTimeSignalDF
            chrBased_SmoothedWaveletReplicationTimeSignal_df = repliseq_signal_df[repliseq_signal_df['chr'] == chrLong]

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

            ################################################################################
            if ((chrBased_SmoothedWaveletReplicationTimeSignal_df is not None) and (not chrBased_SmoothedWaveletReplicationTimeSignal_df.empty) and (checkforValidness(chrBased_valleys_peaks_df))):
                chrBased_replication_array = fill_chr_based_replication_strand_array(chrLong,chromSize,chrBased_SmoothedWaveletReplicationTimeSignal_df,chrBased_valleys_peaks_df)

                #For each chrom
                poolInputList = []

                ################################################################################
                for simNum in range(0,numofSimulations+1):
                    chrBased_subs_df = readChrBasedSubsDF(outputDir, jobname, chrLong, SUBS, simNum)
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


    elif (computationType == COMPUTATION_CHROMOSOMES_SEQUENTIAL_CHROMOSOME_SPLITS_PARALLEL):

        ############################################################################################################
        ###############################       Chromosomes sequentially      ########################################
        ###############################      All ChrBased Splits in parallel starts    #############################
        ############################################################################################################
        for chrLong in chromNamesList:
            chromSize = chromSizesDict[chrLong]

            original_chrBased_subs_df = readChrBasedSubsDF(outputDir, jobname, chrLong, SUBS, 0)
            original_chrBased_indels_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, INDELS, 0)
            original_chrBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, DINUCS, 0)

            combined_chrBased_subs_df = getCombinedChrBasedDF(outputDir, jobname, chrLong, original_chrBased_subs_df,SUBS, numofSimulations)
            combined_chrBased_indels_df = getCombinedChrBasedDF(outputDir, jobname, chrLong,original_chrBased_indels_df, INDELS, numofSimulations)
            combined_chrBased_dinucs_df = getCombinedChrBasedDF(outputDir, jobname, chrLong,original_chrBased_dinucs_df, DINUCS, numofSimulations)

            #Read chrBasedSmoothedWaveletReplicationTimeSignalDF
            chrBased_SmoothedWaveletReplicationTimeSignal_df = repliseq_signal_df[repliseq_signal_df['chr'] == chrLong]

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

            if ((chrBased_SmoothedWaveletReplicationTimeSignal_df is not None) and (not chrBased_SmoothedWaveletReplicationTimeSignal_df.empty) and (checkforValidness(chrBased_valleys_peaks_df))):

                chrBased_replication_array = fill_chr_based_replication_strand_array(chrLong,chromSize,chrBased_SmoothedWaveletReplicationTimeSignal_df,chrBased_valleys_peaks_df)

                chrBased_subs_df_splits_list = None
                chrBased_indels_df_splits_list = None
                chrBased_dinucs_df_splits_list = None

                if ((combined_chrBased_subs_df is not None) and (not combined_chrBased_subs_df.empty)):
                    chrBased_subs_df_splits_list = np.array_split(combined_chrBased_subs_df, numofProcesses)

                if ((combined_chrBased_indels_df is not None) and (not combined_chrBased_indels_df.empty)):
                    chrBased_indels_df_splits_list = np.array_split(combined_chrBased_indels_df, numofProcesses)

                if ((combined_chrBased_dinucs_df is not None) and (not combined_chrBased_dinucs_df.empty)):
                    chrBased_dinucs_df_splits_list = np.array_split(combined_chrBased_dinucs_df, numofProcesses)

                poolInputList = []

                ##########################################################################
                for split_index in range(numofProcesses):
                    chrBased_subs_split_df = None
                    chrBased_indels_split_df = None
                    chrBased_dinucs_split_df = None

                    if ((chrBased_subs_df_splits_list is not None) and (len(chrBased_subs_df_splits_list))):
                        chrBased_subs_split_df = chrBased_subs_df_splits_list[split_index]

                    if ((chrBased_indels_df_splits_list is not None) and (len(chrBased_indels_df_splits_list))):
                        chrBased_indels_split_df = chrBased_indels_df_splits_list[split_index]

                    if ((chrBased_dinucs_df_splits_list is not None) and len(chrBased_dinucs_df_splits_list)):
                        chrBased_dinucs_split_df = chrBased_dinucs_df_splits_list[split_index]

                    inputList = []
                    inputList.append(chrBased_replication_array) #same for all
                    inputList.append(chrBased_subs_split_df) # different split each time
                    inputList.append(chrBased_indels_split_df)  # different split each time
                    inputList.append(chrBased_dinucs_split_df)  # different split each time
                    inputList.append(subsSignature2NumberofMutationsDict)  # same for all
                    inputList.append(indelsSignature2NumberofMutationsDict)
                    inputList.append(dinucsSignature2NumberofMutationsDict)
                    inputList.append(numofSimulations)
                    poolInputList.append(inputList)
                ##########################################################################

                listofTuples = pool.map(searchMutationsOnReplicationStrandArray, poolInputList)

                accumulate_simulations_integrated(listofTuples,
                            accumulatedAllChromosomesType2ReplicationStrand2CountDict,
                            accumulatedAllChromosomesSample2Type2ReplicationStrand2CountDict,
                            accumulatedAllChromosomesType2Sample2ReplicationStrand2CountDict,
                            accumulatedAllChromosomesSignature2MutationType2ReplicationStrand2CountDict)
        ############################################################################################################
        ###############################       Chromosomes sequentially      ########################################
        ###############################      All ChrBased Splits in parallel ends    ###############################
        ############################################################################################################

    print('ReplicationStrandBiasAnalysis Results %s starts' %(computationType))
    print('accumulatedAllChromosomesType2ReplicationStrand2CountDict[0]')
    print(accumulatedAllChromosomesType2ReplicationStrand2CountDict[0])
    print('accumulatedAllChromosomesSample2Type2ReplicationStrand2CountDict[0]')
    print(accumulatedAllChromosomesSample2Type2ReplicationStrand2CountDict[0])
    print('accumulatedAllChromosomesType2Sample2ReplicationStrand2CountDict[0]')
    print(accumulatedAllChromosomesType2Sample2ReplicationStrand2CountDict[0])
    print('accumulatedAllChromosomesSignature2MutationType2ReplicationStrand2CountDict[0]')
    print(accumulatedAllChromosomesSignature2MutationType2ReplicationStrand2CountDict[0])
    print('ReplicationStrandBiasAnalysis Results %s ends' %(computationType))

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

    ################################
    pool.close()
    pool.join()
    ################################

    print('--- ReplicationStrandBias Analysis ends')
    print('#################################################################################\n')
########################################################################
