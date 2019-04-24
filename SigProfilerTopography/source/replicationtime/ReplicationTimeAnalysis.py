# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018 Burcak Otlu

#############################################################
# This python code generates normalized mutation density data for different analyses types which are
# Using all mutations (all chromosomes, all samples, signatures) which leads to Aggregated Substitutions
# Using all mutations but in signature based manner leads to Signature Based
# Using all mutations in signature based and sample based manner leads to Signature based and Sample Based
# Using all indels (all chromosomes, all samples) which leads to Aggregated Indels
# Using all indels (diving indels of length >=3 as microhomology-mediated indels and indels of length <3 as repeat-mediated  indels)

## In this python code
## I read the wavelet-transformed signal wig across whole genome
## Sort the signal in descending order
## Divide the data into 10 equal deciles
## Return 10 deciles as data frames
## Read the mutation file
## Group by mutation file by chromosome
## Group by decile file by chromosome
## For each chromosome (Done in the loop sequentially)
## Make a list of input
## chromosome based mutation file mutation file
## chromosome based decile files coming from 10 different deciles (Done in parallel for all of the 10 deciles)
## In the function I will generate interval tree from decile file
## Overlap the mutations with the interval tree
## Return the total number of mutations that overlap between mutation file and decile file
## Get chromosome based mutation file
## Get chromosome based decile file
## Accumulate the total number of mutations that overlap for each decile file
#############################################################

#############################################################
# Constraints, Thresholds
# Please note that for sample based and signature based replication time analysis
# We consider samples with at least 3000 mutations at total in all deciles.
#############################################################

import twobitreader
from SigProfilerTopography.source.commons.TopographyCommons import *


##################################################################
def fillReplicationTimeNPArrays(inputList):
    chrLong = inputList[0]
    chromSize = inputList[1]
    decile_df_list = inputList[2]
    replicationTimeFilename_wo_extension = inputList[3]

    chrBasedReplicationTimeDataNPArray = fillChrBasedReplicationTimeNPArrayNewVersion(chrLong, chromSize,decile_df_list)
    replicationTimeArrayFilename = '%s_%s' % (chrLong, replicationTimeFilename_wo_extension)
    chrBasedReplicationTimeFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, REPLICATION,replicationTimeFilename_wo_extension, replicationTimeArrayFilename)
    np.save(chrBasedReplicationTimeFile, chrBasedReplicationTimeDataNPArray)
##################################################################


##################################################################
def writeNumberofAttributableBases(decile_df_list,replicationTimeFilename_wo_extension):
    decileIndex2NumberofAttributableBasesDict = {}
    for i,decile_df in enumerate(decile_df_list,1):
        numofAttBases = decile_df[NUMOFBASES].sum()
        decileIndex2NumberofAttributableBasesDict[i] = numofAttBases

    filepath = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,REPLICATION,replicationTimeFilename_wo_extension,DecileIndex2NumfAttributableBasesDictFilename)
    writeDictionaryUsingPickle(decileIndex2NumberofAttributableBasesDict,filepath)
##################################################################

##################################################################
#March 21, 2019 starts
def readReplicationTimeDataAndWriteChrBasedReplicationTimeNPArrays(genome,chromNamesList,chromSizesDict,replicationTimeFilename):

    replicationTimeFilename_wo_extension = os.path.basename(replicationTimeFilename)[0:-4]
    os.makedirs(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, REPLICATION, replicationTimeFilename_wo_extension),exist_ok=True)

    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)

    #Fist decile_df in decile_df_list contains the intervals that are replicated the earliest.
    #Last decile_df in decile_df_list contains the intervals that are replicated the latest.
    #Please note that each decile_df contains intervasl from all chroms (mixed chroms)
    chrNamesInReplicationTimeDataArray, decile_df_list = readRepliSeqTimeData(genome,pool,replicationTimeFilename)
    writeNumberofAttributableBases(decile_df_list,replicationTimeFilename_wo_extension)

    #Get this chrLong of each decile
    #The first decile is the earliest one with index 1
    #The last decile is the latest one with index 10
    #We provide these indexes later in the fillChrBasedDecileIndexArray function by enumerate

    poolInputList = []

    for chrLong in chromNamesList:
        inputList = []
        chromSize = chromSizesDict[chrLong]
        inputList.append(chrLong)
        inputList.append(chromSize)
        inputList.append(decile_df_list)
        inputList.append(replicationTimeFilename_wo_extension)
        poolInputList.append(inputList)

    pool.map(fillReplicationTimeNPArrays,poolInputList)

    ################################
    pool.close()
    pool.join()
    ################################

#March 21, 2019 ends
##################################################################

##################################################################
def generateIntervalVersion(replication_time_wavelet_signal_unprocessed_df):
    #Read the file and generate chr start end signal wavelet_smoothed_signal
    columns = [CHROM, START, END,SIGNAL]

    #Create an empty dataframe
    replication_time_wavelet_signal_interval_version_df = pd.DataFrame(columns=columns)

    #dummy initialization
    i = 0
    chrom = 'chr1'
    start = 0
    step = 1000

    rows_list = []

    # e.g. rows
    # fixedStep chrom=chr1 start=24500 step=1000 span=1000
    # 57.4679
    # 57.467
    # 57.4651
    # 57.4623
    # 57.4586
    # 57.454
    # 57.4484
    # 57.442
    # 57.4347
    # 57.4266
    # 57.4176

    for row in replication_time_wavelet_signal_unprocessed_df.itertuples(index=True, name='Pandas'):
        # row's type is <class 'pandas.core.frame.Pandas'>
        # row[0] is the index
        # row[1] is fixedStep chrom=chr1 start=24500 step=1000 span=1000 or 57.4679
        if (row[1].startswith('fixedStep')):
            #This is the information line
            chrom = row[1].split()[1].split('=')[1]
            start = int(row[1].split()[2].split('=')[1])
            step = int(row[1].split()[3].split('=')[1])
        else:
            signal = float(row[1])
            chr = chrom
            start = start
            end = start + step-1
            dict = {CHROM:chr, START:start, END:end, SIGNAL:signal}
            rows_list.append(dict)
            start += step

    # print('Number of intervals to be inserted in wavelet_processed_df: %d' %len(rows_list))
    #rows_list contain the list of row where each row is a dictionary

    replication_time_wavelet_signal_interval_version_df = pd.DataFrame(rows_list, columns=[CHROM,START,END,SIGNAL])

    # print('replication_time_wavelet_signal_interval_version_df.dtypes')
    # print(replication_time_wavelet_signal_interval_version_df.dtypes)

    replication_time_wavelet_signal_interval_version_df[CHROM] = replication_time_wavelet_signal_interval_version_df[CHROM].astype(str)
    replication_time_wavelet_signal_interval_version_df[START] = replication_time_wavelet_signal_interval_version_df[START].astype(np.int32)
    replication_time_wavelet_signal_interval_version_df[END] = replication_time_wavelet_signal_interval_version_df[END].astype(np.int32)
    replication_time_wavelet_signal_interval_version_df[SIGNAL] = replication_time_wavelet_signal_interval_version_df[SIGNAL].astype(np.float32)

    # print('replication_time_wavelet_signal_interval_version_df.dtypes')
    # print(replication_time_wavelet_signal_interval_version_df.dtypes)

    return replication_time_wavelet_signal_interval_version_df
##################################################################


##################################################################
def readRepliSeqTimeData(genome,pool,repliseqDataFilename):
    ###################################################################
    ############### Read RepliSeq Time data starts ####################
    ###################################################################
    #Read the wavelet signal
    replication_time_wavelet_signal_unprocessed_df = readWaveletSmoothedSignalReplicationTime(repliseqDataFilename)

    #Process the wavelet signal, convert into interval version
    # here column names are added
    replication_time_interval_version_df = generateIntervalVersion(replication_time_wavelet_signal_unprocessed_df)

    chrNamesInReplicationTimeDataArray = replication_time_interval_version_df[CHROM].unique()
    print('Chromosome names in replication time signal data: %s' %(chrNamesInReplicationTimeDataArray))

    #Augment wavelet_processed_df with numberofAttributableBases
    wavelet_processed_augmented_df = augment(genome,pool,replication_time_interval_version_df)

    #Return 10 deciles: the first decile is the earliest one and the tenth decile is the latest one

    #Sort in descending order
    #Higher the replication time signal earlier the replication is
    wavelet_processed_augmented_df.sort_values(SIGNAL, ascending=False, inplace=True)

    # print('############ after sort wavelet_processed_augmented_df ###################')
    # print(wavelet_processed_augmented_df.head())
    # print('############ after sort wavelet_processed_augmented_df ###################')

    #Split wavelet_processed_augmented_df into 10 deciles
    deciles_df_list = np.array_split(wavelet_processed_augmented_df,10)
    # print('Number of decile:%d' %len(deciles))
    #deciles is a list and each decile is a dataframe <class 'pandas.core.frame.DataFrame'>
    #The first decile is the earliest one
    #The last decile is the latest one
    # print('type(deciles_df_list):%s' %type(deciles_df_list))
    # type(deciles_df_list) --> <class 'list'>

    # ############################################################
    # #For information
    # totalNumberofIntervals = 0
    # for decileIndex, decile_df in enumerate(deciles_df_list,1):
    #     # print('decileIndex: %d' %decileIndex)
    #     # print('type(decile_df)')
    #     # print(type(decile_df))
    #     # type(decile_df) --> <class 'pandas.core.frame.DataFrame'>
    #     totalNumberofIntervals += len(decile_df)
    # print('totalNumberofIntervals in all deciles : %d' %totalNumberofIntervals)
    # ############################################################

    return chrNamesInReplicationTimeDataArray, deciles_df_list
    ###################################################################
    ############### Read RepliSeq Time data ends ######################
    ###################################################################

##################################################################



##################################################################
def getNumberofAttributableBases(wavelet_row, chrBasedGenome):
    start =wavelet_row[1]
    end = wavelet_row[2]
    #In my code ends are inclusive
    #twobitreader uses ends exlusive
    seq = chrBasedGenome.get_slice(start, end+1)
    numofAttributableBases = seq.count('A') + seq.count('T') + seq.count('G') + seq.count('C') + seq.count('a') + seq.count('t') + seq.count('g') + seq.count('c')
    # print('######### debug starts ##############')
    # print(wavelet_row)
    # print('len(seq):%d' %len(seq))
    # print('numofAttributableBases:%d' %numofAttributableBases)
    # print('##########  debug ends #############')
    return numofAttributableBases
##################################################################

##################################################################
def addNumofAttributableBasesColumn(inputList):
    chrLong = inputList[0]
    chrBased_wavelet_processed_df_group =inputList[1]
    wholeGenome = inputList[2]
    resulting_df = chrBased_wavelet_processed_df_group.apply(getNumberofAttributableBases, chrBasedGenome = wholeGenome[chrLong], axis= 1)

    # print('######## debug numberofAttributableBases starts ########')
    # print('for %s starts' %chrLong)
    # print('len(chrBased_wavelet_processed_df_group): %d ' %len(chrBased_wavelet_processed_df_group))
    # print('len(resulting_df): %d' %len(resulting_df))
    # print('type(chrBased_wavelet_processed_df_group): %s' % type(chrBased_wavelet_processed_df_group))
    # print('type(resulting_df): %s' %type(resulting_df))

    if (len(chrBased_wavelet_processed_df_group)!=len(resulting_df)):
        print('There is a situation: len(chrBased_wavelet_processed_df_group) is not equal  to len(resulting_df)')

    # print('debug starts')
    # print(resulting_df)
    # print(chrBased_wavelet_processed_df_group.columns.values.tolist())
    # print('debug ends')

    #Add the resulting_df as a new column to chrBased_wavelet_processed_df_group
    chrBased_wavelet_processed_df_group[NUMOFBASES] = resulting_df

    # print('debug starts')
    # print(chrBased_wavelet_processed_df_group.columns.values.tolist())
    # print('debug ends')

    # print('for %s ends' % chrLong)
    # print('######## debug numberofAttributableBases ends ########')

    return (chrLong,chrBased_wavelet_processed_df_group)
##################################################################


##################################################################
def searchMutation(mutation_row,
    sample2NumberofMutationsDict,
    signature2NumberofMutationsDict,
    sample2Signature2NumberofMutationsDict,
    chrBasedReplicationTimeDataArrayWithDecileIndex,
    type2DecileIndex2NumberofMutationsDict,
    sample2Type2DecileIndex2NumberofMutationsDict,
    MUTATION_SIGNATURE_PROBABILITY_THRESHOLD,
    type):

    # For single point mutations start and end are the same, therefore we need to add 1 to mutation_row[END]

    #end is exclusive for subs, indels and dincuc provided by readChrBased methods
    start= mutation_row[START]

    if (type==AGGREGATEDINDELS):
        end = start + mutation_row[LENGTH]
    elif (type==AGGREGATEDSUBSTITUTIONS):
        end = start + 1
    elif (type==AGGREGATEDDINUCS):
        end = start + 2

    slicedArray = chrBasedReplicationTimeDataArrayWithDecileIndex[start:end]

    #get the samplename
    sample = mutation_row[SAMPLE]

    # np.nonzero returns the indices of the elements that are non-zero.
    # np.unique finds the unique elements of an array returns ndarray the sorted unique values.
    # np.nditer efficient multi-dimensional iterator object to iterate over arrays.
    uniqueIndexesArray = np.unique(slicedArray[np.nonzero(slicedArray)])

    if (uniqueIndexesArray.size>0):
        for decileIndex in np.nditer(uniqueIndexesArray):
            # type(decileIndex) is numpy.ndarray
            decileIndexNum = int(decileIndex)

            ################### AGGREGATED starts #####################
            if decileIndexNum in type2DecileIndex2NumberofMutationsDict[type]:
                type2DecileIndex2NumberofMutationsDict[type][decileIndexNum] += 1
            else:
                type2DecileIndex2NumberofMutationsDict[type][decileIndexNum] = 1

            if (type==AGGREGATEDINDELS):
                if mutation_row[LENGTH] >= 3:
                    if decileIndexNum in type2DecileIndex2NumberofMutationsDict[MICROHOMOLOGY]:
                        type2DecileIndex2NumberofMutationsDict[MICROHOMOLOGY][decileIndexNum] += 1
                    else:
                        type2DecileIndex2NumberofMutationsDict[MICROHOMOLOGY][decileIndexNum] = 1
                else:
                    if decileIndexNum in type2DecileIndex2NumberofMutationsDict[REPEAT]:
                        type2DecileIndex2NumberofMutationsDict[REPEAT][decileIndexNum] += 1
                    else:
                        type2DecileIndex2NumberofMutationsDict[REPEAT][decileIndexNum] = 1
            ################### AGGREGATED ends #######################

            ########################### Signatures start ###########################
            for signature in signature2NumberofMutationsDict:
                if mutation_row[signature] >= MUTATION_SIGNATURE_PROBABILITY_THRESHOLD:
                    if decileIndexNum in type2DecileIndex2NumberofMutationsDict[signature]:
                        type2DecileIndex2NumberofMutationsDict[signature][decileIndexNum] += 1
                    else:
                        type2DecileIndex2NumberofMutationsDict[signature][decileIndexNum] = 1
            ########################### Signatures end #############################

            ############################# Sample Based starts ###########################
            if sample in sample2NumberofMutationsDict:
                ############## Sample Based Aggregated Substitutions starts ############
                if decileIndexNum in sample2Type2DecileIndex2NumberofMutationsDict[sample][type]:
                    sample2Type2DecileIndex2NumberofMutationsDict[sample][type][decileIndexNum] += 1
                else:
                    sample2Type2DecileIndex2NumberofMutationsDict[sample][type][decileIndexNum] = 1

                if (type==AGGREGATEDINDELS):
                    if mutation_row[LENGTH] >= 3:
                        if decileIndexNum in sample2Type2DecileIndex2NumberofMutationsDict[sample][MICROHOMOLOGY]:
                            sample2Type2DecileIndex2NumberofMutationsDict[sample][MICROHOMOLOGY][decileIndexNum] += 1
                        else:
                            sample2Type2DecileIndex2NumberofMutationsDict[sample][MICROHOMOLOGY][decileIndexNum] = 1

                    else:
                        if decileIndexNum in sample2Type2DecileIndex2NumberofMutationsDict[sample][REPEAT]:
                            sample2Type2DecileIndex2NumberofMutationsDict[sample][REPEAT][decileIndexNum] += 1
                        else:
                            sample2Type2DecileIndex2NumberofMutationsDict[sample][REPEAT][decileIndexNum] = 1
                ############## Sample Based Aggregated Substitutions ends ##############

            if sample in sample2Signature2NumberofMutationsDict:
                for signature in sample2Signature2NumberofMutationsDict[sample]:
                    if mutation_row[signature] >= MUTATION_SIGNATURE_PROBABILITY_THRESHOLD:
                        if decileIndexNum in sample2Type2DecileIndex2NumberofMutationsDict[sample][signature]:
                            sample2Type2DecileIndex2NumberofMutationsDict[sample][signature][decileIndexNum] += 1
                        else:
                            sample2Type2DecileIndex2NumberofMutationsDict[sample][signature][decileIndexNum] = 1
                ############## Sample Based Signatures ends ############################

            ############################# Sample Based ends ############################

##################################################################


##################################################################
def searchforMutations(sample2NumberofMutationsDict,
        signature2NumberofMutationsDict,
        sample2Signature2NumberofMutationsDict,
        chrBased_mutations_df_split,
        chrBasedReplicationTimeDataArrayWithDecileIndex,
        MUTATION_SIGNATURE_PROBABILITY_THRESHOLD,
        type):

    ############################################################
    # We will fill type2DecileIndex2NumberofMutationsDict for AGGREGATEDSUBSTITUTIONS and for each signature and return it
    type2DecileIndex2NumberofMutationsDict= {}

    # Initialize for AGGREGATEDSUBSTITUTIONS and subs signatures
    # type2DecileIndex2NumberofMutationsDict[AGGREGATEDSUBSTITUTIONS] = {}
    type2DecileIndex2NumberofMutationsDict[type] = {}
    if (type==AGGREGATEDINDELS):
        type2DecileIndex2NumberofMutationsDict[MICROHOMOLOGY] = {}
        type2DecileIndex2NumberofMutationsDict[REPEAT] = {}

    for signature in signature2NumberofMutationsDict:
        type2DecileIndex2NumberofMutationsDict[signature] = {}
    ############################################################

    ############################################################
    #Initialize sample2Type2DecileIndex2NumberofMutationsDict
    sample2Type2DecileIndex2NumberofMutationsDict = {}

    for sample in  sample2NumberofMutationsDict:
        sample2Type2DecileIndex2NumberofMutationsDict[sample] = {}
        sample2Type2DecileIndex2NumberofMutationsDict[sample][type] = {}
        if (type == AGGREGATEDINDELS):
            sample2Type2DecileIndex2NumberofMutationsDict[sample][MICROHOMOLOGY] = {}
            sample2Type2DecileIndex2NumberofMutationsDict[sample][REPEAT] = {}

    for sample in sample2Signature2NumberofMutationsDict:
        for signature in sample2Signature2NumberofMutationsDict[sample]:
            sample2Type2DecileIndex2NumberofMutationsDict[sample][signature] = {}
    ############################################################

    chrBased_mutations_df_split.apply(searchMutation,
                                sample2NumberofMutationsDict = sample2NumberofMutationsDict,
                                signature2NumberofMutationsDict=signature2NumberofMutationsDict,
                                sample2Signature2NumberofMutationsDict= sample2Signature2NumberofMutationsDict,
                                chrBasedReplicationTimeDataArrayWithDecileIndex= chrBasedReplicationTimeDataArrayWithDecileIndex,
                                type2DecileIndex2NumberofMutationsDict=type2DecileIndex2NumberofMutationsDict,
                                sample2Type2DecileIndex2NumberofMutationsDict=sample2Type2DecileIndex2NumberofMutationsDict,
                                MUTATION_SIGNATURE_PROBABILITY_THRESHOLD = MUTATION_SIGNATURE_PROBABILITY_THRESHOLD,
                                type= type,
                                axis=1)

    return type2DecileIndex2NumberofMutationsDict, sample2Type2DecileIndex2NumberofMutationsDict
##################################################################

##################################################################
#Please notice that replication time data are not overlappig data therefore only setting one decileIndex will be correct.
# e.g.:
# start end
#   10 1009
# 1010 2009
def fillArray(chrBased_replicationtimedata_row,chrBasedDecileIndexArray,decileIndex):
    start= chrBased_replicationtimedata_row[START]
    end = chrBased_replicationtimedata_row[END] + 1
    chrBasedDecileIndexArray[start:end] = decileIndex
##################################################################


##################################################################
def  fillChrBasedReplicationTimeNPArrayNewVersion(chrLong,chromSize,decile_df_list):
    # int8	Byte (-128 to 127)
    chrBasedReplicationTimeDataArrayWithDecileIndex = np.zeros(chromSize, dtype=np.int8)

    #First decileIndex is 1, last decile index is 10.
    for decileIndex, decile_df in enumerate(decile_df_list,1):
        chrBased_grouped_decile_df = decile_df.groupby(CHROM)
        # print('debug: len(chrBased_grouped_decile_df): %d' %len(chrBased_grouped_decile_df))
        if chrLong in chrBased_grouped_decile_df.groups.keys():
            chrBased_replicationtimedata_df = chrBased_grouped_decile_df.get_group(chrLong)
            #what is chrBased_decile's type? DataFrame
            if ((chrBased_replicationtimedata_df is not None) and (not chrBased_replicationtimedata_df.empty)):
                chrBased_replicationtimedata_df.apply(fillArray,chrBasedDecileIndexArray=chrBasedReplicationTimeDataArrayWithDecileIndex,decileIndex=decileIndex, axis=1)

    return chrBasedReplicationTimeDataArrayWithDecileIndex
##################################################################


##################################################################
def  fillChrBasedReplicationTimeNPArray(chrLong,chromSize,chrBased_grouped_decile_df_list):
    #We can set the starring index as 1 in builtin function enumerate
    #First chrBased_grouped_decile has index of 1
    #Last chrBased_grouped_decile has index of 10

    # int8	Byte (-128 to 127)
    chrBasedReplicationTimeDataArrayWithDecileIndex = np.zeros(chromSize, dtype=np.int8)

    #First decileIndex is 1, last decile index is 10.
    for decileIndex, chrBased_grouped_decile_df in enumerate(chrBased_grouped_decile_df_list,1):
        if chrLong in chrBased_grouped_decile_df.groups.keys():
            chrBased_replicationtimedata_df = chrBased_grouped_decile_df.get_group(chrLong)
            #what is chrBased_decile's type? DataFrame
            if ((chrBased_replicationtimedata_df is not None) and (not chrBased_replicationtimedata_df.empty)):
                chrBased_replicationtimedata_df.apply(fillArray,chrBasedDecileIndexArray=chrBasedReplicationTimeDataArrayWithDecileIndex,decileIndex=decileIndex, axis=1)

    return chrBasedReplicationTimeDataArrayWithDecileIndex
##################################################################


##################################################################
def fillDictionaries(sample2NumberofSubsDict,
                    sample2NumberofIndelsDict,
                    sample2NumberofDinucsDict,
                    subsSignature2NumberofMutationsDict,
                    indelsSignature2NumberofMutationsDict,
                    dinucsSignature2NumberofMutationsDict,
                    sample2SubsSignature2NumberofMutationsDict,
                    sample2IndelsSignature2NumberofMutationsDict,
                    sample2DinucsSignature2NumberofMutationsDict,
                    chrBasedReplicationTimeDataArrayWithDecileIndex,
                    chrBased_subs_df_split,
                    chrBased_indels_df_split,
                    chrBased_dinucs_df_split):
    type2DecileIndex2NumberofMutationsDict = {}
    sample2Type2DecileIndex2NumberofMutationsDict = {}

    ##################################################################################################################
    if ((chrBased_subs_df_split is not None) and (not chrBased_subs_df_split.empty)):
        subsType2DecileIndex2NumberofMutationsDict, sample2SubsType2DecileIndex2NumberofMutationsDict = searchforMutations(
            sample2NumberofSubsDict,
            subsSignature2NumberofMutationsDict,
            sample2SubsSignature2NumberofMutationsDict,
            chrBased_subs_df_split,
            chrBasedReplicationTimeDataArrayWithDecileIndex,
            SUBSTITUTION_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD,
            AGGREGATEDSUBSTITUTIONS)

        type2DecileIndex2NumberofMutationsDict.update(subsType2DecileIndex2NumberofMutationsDict)

        for sample in sample2SubsType2DecileIndex2NumberofMutationsDict:
            if sample in sample2Type2DecileIndex2NumberofMutationsDict:
                sample2Type2DecileIndex2NumberofMutationsDict[sample].update(sample2SubsType2DecileIndex2NumberofMutationsDict[sample])
            else:
                sample2Type2DecileIndex2NumberofMutationsDict[sample]  = sample2SubsType2DecileIndex2NumberofMutationsDict[sample]
    ##################################################################################################################

    ##################################################################################################################
    if ((chrBased_indels_df_split is not None) and (not chrBased_indels_df_split.empty)):
        indelType2DecileIndex2NumberofMutationsDict, sample2IndelsType2DecileIndex2NumberofMutationsDict = searchforMutations(
            sample2NumberofIndelsDict,
            indelsSignature2NumberofMutationsDict,
            sample2IndelsSignature2NumberofMutationsDict,
            chrBased_indels_df_split,
            chrBasedReplicationTimeDataArrayWithDecileIndex,
            INDEL_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD,
            AGGREGATEDINDELS)

        type2DecileIndex2NumberofMutationsDict.update(indelType2DecileIndex2NumberofMutationsDict)

        for sample in sample2IndelsType2DecileIndex2NumberofMutationsDict:
            if sample in sample2Type2DecileIndex2NumberofMutationsDict:
                sample2Type2DecileIndex2NumberofMutationsDict[sample].update(sample2IndelsType2DecileIndex2NumberofMutationsDict[sample])
            else:
                sample2Type2DecileIndex2NumberofMutationsDict[sample]  = sample2IndelsType2DecileIndex2NumberofMutationsDict[sample]
    ##################################################################################################################

    ##################################################################################################################
    if ((chrBased_dinucs_df_split is not None) and (not chrBased_dinucs_df_split.empty)):
        dinucType2DecileIndex2NumberofMutationsDict, sample2DinucsType2DecileIndex2NumberofMutationsDict = searchforMutations(
            sample2NumberofDinucsDict,
            dinucsSignature2NumberofMutationsDict,
            sample2DinucsSignature2NumberofMutationsDict,
            chrBased_dinucs_df_split,
            chrBasedReplicationTimeDataArrayWithDecileIndex,
            DINUC_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD,
            AGGREGATEDDINUCS)

        type2DecileIndex2NumberofMutationsDict.update(dinucType2DecileIndex2NumberofMutationsDict)

        for sample in sample2DinucsType2DecileIndex2NumberofMutationsDict:
            if sample in sample2Type2DecileIndex2NumberofMutationsDict:
                sample2Type2DecileIndex2NumberofMutationsDict[sample].update(sample2DinucsType2DecileIndex2NumberofMutationsDict[sample])
            else:
                sample2Type2DecileIndex2NumberofMutationsDict[sample]  = sample2DinucsType2DecileIndex2NumberofMutationsDict[sample]
    ##################################################################################################################

    #Question: can we update dictionaries using built in update method?
    # Answer: Be careful. It can used when keys are not overlapping, otherwise it overwrites

    return  type2DecileIndex2NumberofMutationsDict, sample2Type2DecileIndex2NumberofMutationsDict
##################################################################

##################################################################
def searchMutationsOnReplicationTimeNPArray(inputList):
    sample2NumberofSubsDict = inputList[0]
    sample2NumberofIndelsDict = inputList[1]
    sample2NumberofDinucsDict = inputList[2]
    subsSignature2NumberofMutationsDict = inputList[3]
    indelsSignature2NumberofMutationsDict = inputList[4]
    dinucsSignature2NumberofMutationsDict = inputList[5]
    sample2SubsSignature2NumberofMutationsDict = inputList[6]
    sample2IndelsSignature2NumberofMutationsDict = inputList[7]
    sample2DinucsSignature2NumberofMutationsDict = inputList[8]
    chrBased_subs_df_split = inputList[9]
    chrBased_indels_df_split = inputList[10]
    chrBased_dinucs_df_split = inputList[11]
    chrBased_replication_time_np_array = inputList[12]


    type2DecileIndex2NumberofMutationsDict, sample2Type2DecileIndex2NumberofMutationsDict = fillDictionaries(sample2NumberofSubsDict,sample2NumberofIndelsDict,sample2NumberofDinucsDict,
                    subsSignature2NumberofMutationsDict,indelsSignature2NumberofMutationsDict,dinucsSignature2NumberofMutationsDict,
                    sample2SubsSignature2NumberofMutationsDict,sample2IndelsSignature2NumberofMutationsDict,sample2DinucsSignature2NumberofMutationsDict,
                    chrBased_replication_time_np_array,
                    chrBased_subs_df_split,chrBased_indels_df_split,chrBased_dinucs_df_split)

    return type2DecileIndex2NumberofMutationsDict, sample2Type2DecileIndex2NumberofMutationsDict
##################################################################


##################################################################
#Subs, Indels, Dinucs
def generateReplicationTimeNPArrayAndSearchMutationsOnNPArray(inputList):
    chrLong = inputList[0]
    chromSize = inputList[1]

    sample2NumberofSubsDict = inputList[2]
    sample2NumberofIndelsDict = inputList[3]
    sample2NumberofDinucsDict = inputList[4]

    subsSignature2NumberofMutationsDict = inputList[5]
    indelsSignature2NumberofMutationsDict = inputList[6]
    dinucsSignature2NumberofMutationsDict = inputList[7]

    sample2SubsSignature2NumberofMutationsDict = inputList[8]
    sample2IndelsSignature2NumberofMutationsDict = inputList[9]
    sample2DinucsSignature2NumberofMutationsDict = inputList[10]

    chrBased_subs_df_split = inputList[11]
    chrBased_indels_df_split = inputList[12]
    chrBased_dinucs_df_split = inputList[13]

    chrBased_grouped_decile_df_list = inputList[14]

    #fill nparray slices with decileIndex and return it in the function fillChrBasedDecileIndexArray
    chrBasedReplicationTimeDataArrayWithDecileIndex = fillChrBasedReplicationTimeNPArray(chrLong,chromSize,chrBased_grouped_decile_df_list)

    type2DecileIndex2NumberofMutationsDict, sample2Type2DecileIndex2NumberofMutationsDict = fillDictionaries(sample2NumberofSubsDict,sample2NumberofIndelsDict,sample2NumberofDinucsDict,
                    subsSignature2NumberofMutationsDict,indelsSignature2NumberofMutationsDict,dinucsSignature2NumberofMutationsDict,
                    sample2SubsSignature2NumberofMutationsDict,sample2IndelsSignature2NumberofMutationsDict,sample2DinucsSignature2NumberofMutationsDict,
                    chrBasedReplicationTimeDataArrayWithDecileIndex,
                    chrBased_subs_df_split,chrBased_indels_df_split,chrBased_dinucs_df_split)

    return type2DecileIndex2NumberofMutationsDict, sample2Type2DecileIndex2NumberofMutationsDict
##################################################################



##################################################################
#Accumulate the counts coming from each chromosome split in type2DecileBasedAllChrAccumulatedCountDict and sample2Type2DecileBasedAllChrAccumulatedCountDict
def accumulateTypeBasedDictionaries(listofDictionaries,type2DecileBasedAllChrAccumulatedCountDict,sample2Type2DecileBasedAllChrAccumulatedCountDict):
    for dictList in listofDictionaries:

        typeBasedDict = dictList[0]
        sampleBasedTypeDict = dictList[1]

        ###########################################################################
        for type, decileBasedAllChrAccumulatedCountDict in typeBasedDict.items():

            if type not in type2DecileBasedAllChrAccumulatedCountDict:
                type2DecileBasedAllChrAccumulatedCountDict[type] = {}

            for decileIndex, decileCount in decileBasedAllChrAccumulatedCountDict.items():
                #if key in dictionary
                if decileIndex in type2DecileBasedAllChrAccumulatedCountDict[type]:
                    type2DecileBasedAllChrAccumulatedCountDict[type][decileIndex] += decileCount
                else:
                    type2DecileBasedAllChrAccumulatedCountDict[type][decileIndex] = decileCount
        ###########################################################################

        ######################## Sample Based starts #######################################
        for sample in sampleBasedTypeDict:
            if sample not in sample2Type2DecileBasedAllChrAccumulatedCountDict:
                sample2Type2DecileBasedAllChrAccumulatedCountDict[sample] = {}

            for type, decileBasedAllChrAccumulatedCountDict in sampleBasedTypeDict[sample].items():
                if type not in sample2Type2DecileBasedAllChrAccumulatedCountDict[sample]:
                    sample2Type2DecileBasedAllChrAccumulatedCountDict[sample][type] = {}

                for decileIndex, decileCount in decileBasedAllChrAccumulatedCountDict.items():
                    # if key in dictionary
                    if decileIndex in sample2Type2DecileBasedAllChrAccumulatedCountDict[sample][type]:
                        sample2Type2DecileBasedAllChrAccumulatedCountDict[sample][type][decileIndex] += decileCount
                    else:
                        sample2Type2DecileBasedAllChrAccumulatedCountDict[sample][type][decileIndex] = decileCount
        ######################## Sample Based ends #########################################

##################################################################


##################################################################
def getNormalizedMutationDensityList(mutationDensityDict):
    densities = mutationDensityDict.values()
    maxDensity = max(densities)
    if maxDensity>0:
        normalizedMutationDensities = [x/maxDensity for x in densities]

    else:
        normalizedMutationDensities = densities

    return normalizedMutationDensities
##################################################################



##################################################################
#March 22, 2019 starts
def getMutationDensityDictNewVersion(decileIndex2NumberofAttributableBasesDict,decileBasedAllChrAccumulatedCountDict):
    decileBasedMutationDensityDict = {}
    numberofMutations = 0

    for decileIndex in decileIndex2NumberofAttributableBasesDict:
        if (decileIndex in decileBasedAllChrAccumulatedCountDict):
            count = decileBasedAllChrAccumulatedCountDict[decileIndex]
            numberofMutations += count
            numofAttBases = decileIndex2NumberofAttributableBasesDict[decileIndex]
            mutationDensity = float(count) / numofAttBases
            decileBasedMutationDensityDict[decileIndex] = mutationDensity
        else:
            decileBasedMutationDensityDict[decileIndex] = 0

    return numberofMutations, decileBasedMutationDensityDict
#March 22, 2019 starts
##################################################################

##################################################################
def getMutationDensityDict(decile_df_list,decileBasedAllChrAccumulatedCountDict):
    decileBasedMutationDensityDict = {}
    numberofMutations = 0

    #Modifiled as enumerate(deciles,1) formerly it was enumerate(deciles,0)
    for i,decile_df in enumerate(decile_df_list,1):
        if (i in decileBasedAllChrAccumulatedCountDict):
            count = decileBasedAllChrAccumulatedCountDict[i]
            numberofMutations += count
            numofAttBases = decile_df[NUMOFBASES].sum()
            mutationDensity=float(count)/numofAttBases
            decileBasedMutationDensityDict[i] = mutationDensity
        else:
            decileBasedMutationDensityDict[i] = 0

        # print('decile: %d numofAttBases: %d' %(i,numofAttBases))

    return numberofMutations, decileBasedMutationDensityDict
##################################################################

##################################################################
def loadChrBasedReplicationTimeNPArray(chrLong,replicationTimeDataFilename):
    replicationTimeFilename_wo_extension = os.path.basename(replicationTimeDataFilename)[0:-4]

    chrBasedReplicationTimeArrayFilename = '%s_%s.npy' % (chrLong, replicationTimeFilename_wo_extension)

    chrBasedReplicationTimeArrayFilenamePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, REPLICATION,
                                                            replicationTimeFilename_wo_extension,chrBasedReplicationTimeArrayFilename )
    if (os.path.exists(chrBasedReplicationTimeArrayFilenamePath)):
        chrBased_replication_time_np_array = np.load(chrBasedReplicationTimeArrayFilenamePath)
    return chrBased_replication_time_np_array
##################################################################

##################################################################
def calculateCountsForMutationsUsingReplicationTimeNPArray(mutationTypes,computationType,outputDir,jobname,numofProcesses,pool,
                                    sample2NumberofSubsDict,
                                    sample2NumberofIndelsDict,
                                    sample2NumberofDinucsDict,
                                    subsSignature2NumberofMutationsDict,
                                    indelsSignature2NumberofMutationsDict,
                                    dinucsSignature2NumberofMutationsDict,
                                    sample2SubsSignature2NumberofMutationsDict,
                                    sample2IndelsSignature2NumberofMutationsDict,
                                    sample2DinucsSignature2NumberofMutationsDict,
                                    replicationTimeDataFilename,
                                    chromNamesList):

    type2DecileBasedAllChrAccumulatedCountDict = {}
    sample2Type2DecileBasedAllChrAccumulatedCountDict ={}

    if (computationType == COMPUTATION_ALL_CHROMOSOMES_PARALLEL):
        #It seems that chrNames in replicationTimeData are long chr names such as chr1, chrX

        ###################################################################################################################################################
        poolInputList = []
        for chrLong in chromNamesList:
            if (SUBS in mutationTypes):
                chrBased_subs_df = readChrBasedSubsDF(outputDir, jobname, chrLong, SUBS)

            if (INDELS in mutationTypes):
                chrBased_indels_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, INDELS)

            if (DINUCS in mutationTypes):
                chrBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, DINUCS)

            chrBased_replication_time_np_array = loadChrBasedReplicationTimeNPArray(chrLong,replicationTimeDataFilename)

            inputList = []
            inputList.append(sample2NumberofSubsDict)
            inputList.append(sample2NumberofIndelsDict)
            inputList.append(sample2NumberofDinucsDict)
            inputList.append(subsSignature2NumberofMutationsDict)
            inputList.append(indelsSignature2NumberofMutationsDict)
            inputList.append(dinucsSignature2NumberofMutationsDict)
            inputList.append(sample2SubsSignature2NumberofMutationsDict)
            inputList.append(sample2IndelsSignature2NumberofMutationsDict)
            inputList.append(sample2DinucsSignature2NumberofMutationsDict)
            inputList.append(chrBased_subs_df) #Different
            inputList.append(chrBased_indels_df) #Different
            inputList.append(chrBased_dinucs_df)  # Different

            inputList.append(chrBased_replication_time_np_array)  #Same
            poolInputList.append(inputList)
        ###################################################################################################################################################

        # Call the parallel code for poolInputList
        # It returns a list of whatever generateNPArrayAndSearchMutationsOnNPArray will return
        # Since generateNPArrayAndSearchMutationsOnNPArray return tuple
        # It must return list of tuples
        # Fill np.array with 0 if there is no signal, otherwise fill with the decileIndex if there is a signal
        listofDictionaries = pool.map(searchMutationsOnReplicationTimeNPArray,poolInputList)
        accumulateTypeBasedDictionaries(listofDictionaries,type2DecileBasedAllChrAccumulatedCountDict,sample2Type2DecileBasedAllChrAccumulatedCountDict)
        ###################################################################################################################################################

    elif (computationType == COMPUTATION_CHROMOSOMES_SEQUENTIAL_CHROMOSOME_SPLITS_PARALLEL):
        #It seems that chrNames in replicationTimeData are long chr names such as chr1, chrX
        for chrLong in chromNamesList:
            if (SUBS in mutationTypes):
                chrBased_subs_df = readChrBasedSubsDF(outputDir, jobname, chrLong, SUBS)

            if (INDELS in mutationTypes):
                chrBased_indels_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, INDELS)

            if (DINUCS in mutationTypes):
                chrBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, DINUCS)

            chrBased_replication_time_np_array = loadChrBasedReplicationTimeNPArray(chrLong,replicationTimeDataFilename)

            chrBased_subs_df_splits_list = []
            chrBased_indels_df_splits_list = []

            if ((chrBased_subs_df is not None) and (not chrBased_subs_df.empty)):
                chrBased_subs_df_splits_list = np.array_split(chrBased_subs_df, numofProcesses)

            if ((chrBased_indels_df is not None) and (not chrBased_indels_df.empty)):
                chrBased_indels_df_splits_list = np.array_split(chrBased_indels_df, numofProcesses)

            if ((chrBased_dinucs_df is not None) and (not chrBased_dinucs_df.empty)):
                chrBased_dinucs_df_splits_list = np.array_split(chrBased_dinucs_df, numofProcesses)

            poolInputList = []
            ###################################################################################################################################################
            for split_index in range(numofProcesses):
                chrBased_subs_df_split_array = None
                chrBased_indels_df_split_array = None
                chrBased_dinucs_df_split_array = None

                if (len(chrBased_subs_df_splits_list)):
                    chrBased_subs_df_split_array = chrBased_subs_df_splits_list[split_index]

                if (len(chrBased_indels_df_splits_list)):
                    chrBased_indels_df_split_array = chrBased_indels_df_splits_list[split_index]

                if (len(chrBased_dinucs_df_splits_list)):
                    chrBased_dinucs_df_split_array = chrBased_dinucs_df_splits_list[split_index]

                inputList = []
                inputList.append(sample2NumberofSubsDict)
                inputList.append(sample2NumberofIndelsDict)
                inputList.append(sample2NumberofDinucsDict)
                inputList.append(subsSignature2NumberofMutationsDict)
                inputList.append(indelsSignature2NumberofMutationsDict)
                inputList.append(dinucsSignature2NumberofMutationsDict)
                inputList.append(sample2SubsSignature2NumberofMutationsDict)
                inputList.append(sample2IndelsSignature2NumberofMutationsDict)
                inputList.append(sample2DinucsSignature2NumberofMutationsDict)
                inputList.append(chrBased_subs_df_split_array) #Different
                inputList.append(chrBased_indels_df_split_array) #Different
                inputList.append(chrBased_dinucs_df_split_array)  # Different
                inputList.append(chrBased_replication_time_np_array)  #Same
                poolInputList.append(inputList)
            ###################################################################################################################################################

            # Call the parallel code for poolInputList
            # It returns a list of whatever generateNPArrayAndSearchMutationsOnNPArray will return
            # Since generateNPArrayAndSearchMutationsOnNPArray return tuple
            # It must return list of tuples
            # Fill np.array with 0 if there is no signal, otherwise fill with the decileIndex if there is a signal
            listofDictionaries = pool.map(searchMutationsOnReplicationTimeNPArray,poolInputList)
            accumulateTypeBasedDictionaries(listofDictionaries,type2DecileBasedAllChrAccumulatedCountDict,sample2Type2DecileBasedAllChrAccumulatedCountDict)
            ###################################################################################################################################################

    elif (computationType == COMPUTATION_CHROMOSOMES_SEQUENTIAL):
        #It seems that chrNames in replicationTimeData are long chr names such as chr1, chrX
        for chrLong in chromNamesList:
            if (SUBS in mutationTypes):
                chrBased_subs_df = readChrBasedSubsDF(outputDir, jobname, chrLong, SUBS)

            if (INDELS in mutationTypes):
                chrBased_indels_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, INDELS)

            if (DINUCS in mutationTypes):
                chrBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, DINUCS)

            chrBased_replication_time_np_array = loadChrBasedReplicationTimeNPArray(chrLong,replicationTimeDataFilename)

            inputList = []
            inputList.append(sample2NumberofSubsDict)
            inputList.append(sample2NumberofIndelsDict)
            inputList.append(sample2NumberofDinucsDict)
            inputList.append(subsSignature2NumberofMutationsDict)
            inputList.append(indelsSignature2NumberofMutationsDict)
            inputList.append(dinucsSignature2NumberofMutationsDict)
            inputList.append(sample2SubsSignature2NumberofMutationsDict)
            inputList.append(sample2IndelsSignature2NumberofMutationsDict)
            inputList.append(sample2DinucsSignature2NumberofMutationsDict)
            inputList.append(chrBased_subs_df) #Different
            inputList.append(chrBased_indels_df) #Different
            inputList.append(chrBased_dinucs_df)  # Different
            inputList.append(chrBased_replication_time_np_array)  #Same
            ###################################################################################################################################################

            type2DecileIndex2NumberofMutationsDict, sample2Type2DecileIndex2NumberofMutationsDict = searchMutationsOnReplicationTimeNPArray(inputList)
            listofDictionaries = []
            listofDictionaries.append((type2DecileIndex2NumberofMutationsDict, sample2Type2DecileIndex2NumberofMutationsDict))
            accumulateTypeBasedDictionaries(listofDictionaries,type2DecileBasedAllChrAccumulatedCountDict,sample2Type2DecileBasedAllChrAccumulatedCountDict)
            ###################################################################################################################################################

    return type2DecileBasedAllChrAccumulatedCountDict, sample2Type2DecileBasedAllChrAccumulatedCountDict
##################################################################


##################################################################
def calculateCountsForMutationsFillingReplicationTimeNPArrayRuntime(mutationTypes,computationType,outputDir,jobname,numofProcesses,pool,
                                    sample2NumberofSubsDict,
                                    sample2NumberofIndelsDict,
                                    sample2NumberofDinucsDict,
                                    subsSignature2NumberofMutationsDict,
                                    indelsSignature2NumberofMutationsDict,
                                    dinucsSignature2NumberofMutationsDict,
                                    sample2SubsSignature2NumberofMutationsDict,
                                    sample2IndelsSignature2NumberofMutationsDict,
                                    sample2DinucsSignature2NumberofMutationsDict,
                                    chromSizesDict,
                                    decile_df_list,
                                    chrNamesList):

    type2DecileBasedAllChrAccumulatedCountDict = {}
    sample2Type2DecileBasedAllChrAccumulatedCountDict ={}

    #Get chrBased grouped deciles
    chrBased_grouped_decile_df_list = []
    #Get this chrLong of each decile
    #The first decile is the earliest one with index 1
    #The last decile is the latest one with index 10
    #We provide these indexes later in the fillChrBasedDecileIndexArray function by enumerate
    for decile_df in decile_df_list:
        chrBased_grouped_decile_df = decile_df.groupby(CHROM)
        chrBased_grouped_decile_df_list.append(chrBased_grouped_decile_df)

    if (computationType == COMPUTATION_ALL_CHROMOSOMES_PARALLEL):
        # It seems that chrNames in replicationTimeData are long chr names such as chr1, chrX
        poolInputList = []
        for chrLong in chrNamesList:
            chromSize = chromSizesDict[chrLong]

            if (SUBS in mutationTypes):
                chrBased_subs_df = readChrBasedSubsDF(outputDir, jobname, chrLong, SUBS)

            if (INDELS in mutationTypes):
                chrBased_indels_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, INDELS)

            if (DINUCS in mutationTypes):
                chrBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong,DINUCS)

            ###################################################################################################################################################
            inputList = []
            inputList.append(chrLong)
            inputList.append(chromSize)
            inputList.append(sample2NumberofSubsDict)
            inputList.append(sample2NumberofIndelsDict)
            inputList.append(sample2NumberofDinucsDict)
            inputList.append(subsSignature2NumberofMutationsDict)
            inputList.append(indelsSignature2NumberofMutationsDict)
            inputList.append(dinucsSignature2NumberofMutationsDict)
            inputList.append(sample2SubsSignature2NumberofMutationsDict)
            inputList.append(sample2IndelsSignature2NumberofMutationsDict)
            inputList.append(sample2DinucsSignature2NumberofMutationsDict)
            inputList.append(chrBased_subs_df)  # Different
            inputList.append(chrBased_indels_df)  # Different
            inputList.append(chrBased_dinucs_df)  # Different
            inputList.append(chrBased_grouped_decile_df_list)  # Same
            poolInputList.append(inputList)

        # Call the parallel code for poolInputList
        # It returns a list of whatever generateNPArrayAndSearchMutationsOnNPArray will return
        # It must return list of tuples
        # Since generateNPArrayAndSearchMutationsOnNPArray return tuple
        # Fill np.array with 0 if there is no signal, otherwise fill with the decileIndex if there is a signal
        listofDictionaries = pool.map(generateReplicationTimeNPArrayAndSearchMutationsOnNPArray, poolInputList)
        accumulateTypeBasedDictionaries(listofDictionaries, type2DecileBasedAllChrAccumulatedCountDict,
                                        sample2Type2DecileBasedAllChrAccumulatedCountDict)
        ###################################################################################################################################################

    elif (computationType == COMPUTATION_CHROMOSOMES_SEQUENTIAL_CHROMOSOME_SPLITS_PARALLEL):
        #It seems that chrNames in replicationTimeData are long chr names such as chr1, chrX
        for chrLong in chrNamesList:
            # chrLong = 'chr%s' %(chr)
            chromSize = chromSizesDict[chrLong]

            if (SUBS in mutationTypes):
                chrBased_subs_df = readChrBasedSubsDF(outputDir, jobname, chrLong, SUBS)

            if (INDELS in mutationTypes):
                chrBased_indels_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, INDELS)

            if (DINUCS in mutationTypes):
                chrBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, DINUCS)

            chrBased_subs_df_splits_list = []
            chrBased_indels_df_splits_list = []
            chrBased_dinucs_df_splits_list = []

            if ((chrBased_subs_df is not None) and (not chrBased_subs_df.empty)):
                chrBased_subs_df_splits_list = np.array_split(chrBased_subs_df, numofProcesses)

            if ((chrBased_indels_df is not None) and (not chrBased_indels_df.empty)):
                chrBased_indels_df_splits_list = np.array_split(chrBased_indels_df, numofProcesses)

            if ((chrBased_dinucs_df is not None) and (not chrBased_dinucs_df.empty)):
                chrBased_dinucs_df_splits_list = np.array_split(chrBased_dinucs_df, numofProcesses)

            ###################################################################################################################################################
            poolInputList = []
            for split_index in range(numofProcesses):
                if (len(chrBased_subs_df_splits_list)):
                    chrBased_subs_df_split_array = chrBased_subs_df_splits_list[split_index]
                else:
                    chrBased_subs_df_split_array = None

                if (len(chrBased_indels_df_splits_list)):
                    chrBased_indels_df_split_array = chrBased_indels_df_splits_list[split_index]
                else:
                    chrBased_indels_df_split_array = None

                if (len(chrBased_dinucs_df_splits_list)):
                    chrBased_dinucs_df_split_array = chrBased_dinucs_df_splits_list[split_index]
                else:
                    chrBased_dinucs_df_split_array = None

                inputList = []
                inputList.append(chrLong) #Same
                inputList.append(chromSize)
                inputList.append(sample2NumberofSubsDict)
                inputList.append(sample2NumberofIndelsDict)
                inputList.append(sample2NumberofDinucsDict)
                inputList.append(subsSignature2NumberofMutationsDict)
                inputList.append(indelsSignature2NumberofMutationsDict)
                inputList.append(dinucsSignature2NumberofMutationsDict)
                inputList.append(sample2SubsSignature2NumberofMutationsDict)
                inputList.append(sample2IndelsSignature2NumberofMutationsDict)
                inputList.append(sample2DinucsSignature2NumberofMutationsDict)
                inputList.append(chrBased_subs_df_split_array) #Different
                inputList.append(chrBased_indels_df_split_array) #Different
                inputList.append(chrBased_dinucs_df_split_array)  # Different
                inputList.append(chrBased_grouped_decile_df_list)  #Same
                poolInputList.append(inputList)

            # Call the parallel code for poolInputList
            # It returns a list of whatever generateNPArrayAndSearchMutationsOnNPArray will return
            # It must return list of tuples
            # Since generateNPArrayAndSearchMutationsOnNPArray return tuple
            # Fill np.array with 0 if there is no signal, otherwise fill with the decileIndex if there is a signal
            listofDictionaries = pool.map(generateReplicationTimeNPArrayAndSearchMutationsOnNPArray,poolInputList)
            accumulateTypeBasedDictionaries(listofDictionaries,type2DecileBasedAllChrAccumulatedCountDict,sample2Type2DecileBasedAllChrAccumulatedCountDict)
            ###################################################################################################################################################

    return type2DecileBasedAllChrAccumulatedCountDict, sample2Type2DecileBasedAllChrAccumulatedCountDict
##################################################################



##################################################################
def augment(genome,pool,wavelet_processed_df):
    if (genome==GRCh37):
        wholeGenome = twobitreader.TwoBitFile(os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,UCSCGENOME,HG19_2BIT))
    elif (genome==GRCh38):
        wholeGenome = twobitreader.TwoBitFile(os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,UCSCGENOME,HG38_2BIT))

    #Augment in parallel for each chromosome
    poolInputList = []
    chrBased_wavelet_processed_df_groups = wavelet_processed_df.groupby(CHROM)
    for chr, chrBased_wavelet_processed_df_group in chrBased_wavelet_processed_df_groups:
        inputList = []
        inputList.append(chr)
        inputList.append(chrBased_wavelet_processed_df_group)
        #Please note that when you provide the chr based hg19_genome it gives error
        inputList.append(wholeGenome)
        poolInputList.append(inputList)

    # print('Augmentation starts')
    #Each tuple contains chrLong and the dataframe with augmented column with number of attributable bases
    listofTuples = pool.map(addNumofAttributableBasesColumn,poolInputList)
    # print('len(poolInputList): %d' %len(poolInputList))
    # print('len(listofDFs): %d' %len(listofTuples))

    #Define frames which will be a list of dataframes
    frames = []

    for tuple in listofTuples:
        #chrLong = tuple[0]
        chrBased_augmented_df = tuple[1]
        frames.append(chrBased_augmented_df)

    augment_df = pd.concat(frames,ignore_index=True)

    return augment_df
##################################################################


##################################################################
def writeReplicationTimeData(outputDir,jobname,decile_df_list,decileIndex2NumberofAttributableBasesDict,type2DecileBasedAllChrAccumulatedCountDict,sample2Type2DecileBasedAllChrAccumulatedCountDict):
    os.makedirs(os.path.join(outputDir, jobname, DATA, REPLICATIONTIME), exist_ok=True)

    ##############################################################
    for type in type2DecileBasedAllChrAccumulatedCountDict:
        decileBasedAllChrAccumulatedCountDict = type2DecileBasedAllChrAccumulatedCountDict[type]

        normalizedMutationDensityFilename = type + '_NormalizedMutationDensity.txt'
        normalizedMutationDensityFilePath = os.path.join(outputDir,jobname, DATA, REPLICATIONTIME,normalizedMutationDensityFilename)

        # If decileBasedAllChrAccumulatedCountDict is not empty
        if (decileBasedAllChrAccumulatedCountDict):

            if (decile_df_list is not None):
                numberofMutations, mutationDensityDict = getMutationDensityDict(decile_df_list,decileBasedAllChrAccumulatedCountDict)
            elif (decileIndex2NumberofAttributableBasesDict is not None):
                numberofMutations, mutationDensityDict = getMutationDensityDictNewVersion(decileIndex2NumberofAttributableBasesDict, decileBasedAllChrAccumulatedCountDict)


            normalizedMutationDensityList = getNormalizedMutationDensityList(mutationDensityDict)

            with open(normalizedMutationDensityFilePath, 'w') as file:
                for normalizedMutationDensity in normalizedMutationDensityList:
                    file.write(str(normalizedMutationDensity) + ' ')
                file.write('\n')
    ##############################################################


    ######################### Sample Based starts #####################
    for sample in sample2Type2DecileBasedAllChrAccumulatedCountDict:
        for type in sample2Type2DecileBasedAllChrAccumulatedCountDict[sample]:
            decileBasedAllChrAccumulatedCountDict = sample2Type2DecileBasedAllChrAccumulatedCountDict[sample][type]

            normalizedMutationDensityFilename = '%s_%s_NormalizedMutationDensity.txt' %(sample,type)
            normalizedMutationDensityFilePath = os.path.join(outputDir,jobname, DATA, REPLICATIONTIME,normalizedMutationDensityFilename)

            #If decileBasedAllChrAccumulatedCountDict is not empty
            if (decileBasedAllChrAccumulatedCountDict):
                if (decile_df_list is not None):
                    numberofMutations, mutationDensityDict = getMutationDensityDict(decile_df_list, decileBasedAllChrAccumulatedCountDict)
                elif (decileIndex2NumberofAttributableBasesDict is not None):
                    numberofMutations, mutationDensityDict = getMutationDensityDictNewVersion(decileIndex2NumberofAttributableBasesDict, decileBasedAllChrAccumulatedCountDict)
                normalizedMutationDensityList = getNormalizedMutationDensityList(mutationDensityDict)

                with open(normalizedMutationDensityFilePath, 'w') as file:
                    for normalizedMutationDensity in normalizedMutationDensityList:
                        file.write(str(normalizedMutationDensity) + ' ')
                    file.write('\n')
    ######################### Sample Based ends #######################

##################################################################


##################################################################
#main function
def replicationTimeAnalysis(mutationTypes,computationType,replication_time_np_arrays_fill_runtime,genome,chromSizesDict,chromNamesList,outputDir,jobname,repliseqDataFilename):

    print('########################## ReplicationTimeAnalysis starts #########################')
    #########################################################################
    # Analysis Type can be
    # AggregatedSubstitutions: All in one
    # AggregatedIndels : All in one
    # AggregatedDinucs : All in one
    # IndelsBased : Microhomology, Repeat
    # SignatureBased: Subs Signatures Sig1, Sig2, ... and Indels Signatures  ID1, ID2, ..., ... and Dinucs Signatures  DBS1, DBS2, ...

    # We know the indels type we are interested in.
    # Microhomology indels --- len(indels) >= 3
    # Repeat indels --- len(indels) < 3
    #########################################################################

    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)

    ##########################################################################################
    sample2NumberofSubsDict = getSample2NumberofSubsDict(outputDir,jobname)
    sample2NumberofIndelsDict = getSample2NumberofIndelsDict(outputDir,jobname)
    sample2NumberofDinucsDict = getDictionary(outputDir,jobname,Sample2NumberofDinucsDictFilename)

    subsSignature2NumberofMutationsDict = getSubsSignature2NumberofMutationsDict(outputDir,jobname)
    indelsSignature2NumberofMutationsDict = getIndelsSignature2NumberofMutationsDict(outputDir,jobname)
    dinucsSignature2NumberofMutationsDict = getDictionary(outputDir,jobname,DinucsSignature2NumberofMutationsDictFilename)

    sample2SubsSignature2NumberofMutationsDict = getSample2SubsSignature2NumberofMutationsDict(outputDir,jobname)
    sample2IndelsSignature2NumberofMutationsDict = getSample2IndelsSignature2NumberofMutationsDict(outputDir,jobname)
    sample2DinucsSignature2NumberofMutationsDict = getDictionary(outputDir,jobname,Sample2DinucsSignature2NumberofMutationsDictFilename)
    ##########################################################################################

    # Fill replication np arrays runtime
    if (replication_time_np_arrays_fill_runtime):

        ###################################################################
        ############### Read MCF-7 RepliSeq Time data starts ##############
        ###################################################################
        # Fist decile_df in decile_df_list contains the intervals that are replicated the earliest.
        # Last decile_df in decile_df_list contains the intervals that are replicated the latest.
        # Please note that each decile_df contains intervasl from all chroms (mixed chroms)
        #What is the type of deciles? Deciles is a list of dataframes.
        chrNamesInReplicationTimeDataArray, decile_df_list = readRepliSeqTimeData(genome,pool,repliseqDataFilename)

        type2DecileBasedAllChrAccumulatedCountDict , sample2Type2DecileBasedAllChrAccumulatedCountDict = calculateCountsForMutationsFillingReplicationTimeNPArrayRuntime(
                                                                                                                                mutationTypes,
                                                                                                                                computationType,
                                                                                                                                outputDir,jobname,numofProcesses,pool,
                                                                                                                                sample2NumberofSubsDict,
                                                                                                                                sample2NumberofIndelsDict,
                                                                                                                                sample2NumberofDinucsDict,
                                                                                                                                subsSignature2NumberofMutationsDict,
                                                                                                                                indelsSignature2NumberofMutationsDict,
                                                                                                                                dinucsSignature2NumberofMutationsDict,
                                                                                                                                sample2SubsSignature2NumberofMutationsDict,
                                                                                                                                sample2IndelsSignature2NumberofMutationsDict,
                                                                                                                                sample2DinucsSignature2NumberofMutationsDict,
                                                                                                                                chromSizesDict,
                                                                                                                                decile_df_list,
                                                                                                                                chromNamesList)
        writeReplicationTimeData(outputDir,jobname,decile_df_list,None,type2DecileBasedAllChrAccumulatedCountDict,sample2Type2DecileBasedAllChrAccumulatedCountDict)
        ###################################################################
        ############### Read MCF-7 RepliSeq Time data ends ################
        ###################################################################

    # Use offline or online prepared replication np arrays
    else:
        replicationTimeFilename_wo_extension = os.path.basename(repliseqDataFilename)[0:-4]
        decileIndex2NumberofAttributableBasesDict = getDecileIndex2NumberofAttributableBasesDict(replicationTimeFilename_wo_extension)

        #loading offline prepared numpy array
        type2DecileBasedAllChrAccumulatedCountDict , sample2Type2DecileBasedAllChrAccumulatedCountDict = calculateCountsForMutationsUsingReplicationTimeNPArray(mutationTypes,
                                                                                                                                computationType,
                                                                                                                                outputDir,jobname,numofProcesses,pool,
                                                                                                                                sample2NumberofSubsDict,
                                                                                                                                sample2NumberofIndelsDict,
                                                                                                                                sample2NumberofDinucsDict,
                                                                                                                                subsSignature2NumberofMutationsDict,
                                                                                                                                indelsSignature2NumberofMutationsDict,
                                                                                                                                dinucsSignature2NumberofMutationsDict,
                                                                                                                                sample2SubsSignature2NumberofMutationsDict,
                                                                                                                                sample2IndelsSignature2NumberofMutationsDict,
                                                                                                                                sample2DinucsSignature2NumberofMutationsDict,
                                                                                                                                repliseqDataFilename,
                                                                                                                                chromNamesList)
        writeReplicationTimeData(outputDir,jobname,None,decileIndex2NumberofAttributableBasesDict,type2DecileBasedAllChrAccumulatedCountDict,sample2Type2DecileBasedAllChrAccumulatedCountDict)

    #######################################################################################################
    ############### Carry out Replication Time Data Analysis for each analysis type starts ################
    #######################################################################################################

    ################################
    pool.close()
    pool.join()
    ################################

    print('########################## ReplicationTimeAnalysis ends ############################')

##################################################################