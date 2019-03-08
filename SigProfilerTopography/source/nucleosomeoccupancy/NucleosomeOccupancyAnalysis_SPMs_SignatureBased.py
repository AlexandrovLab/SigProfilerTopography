# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018 Burcak Otlu

###############################################################################################################
# In this python code, nucleosome occupancy analysis is carried out
#   for all single point mutations
#   for all signatures with all single point mutations having probability >= 0.5 for that signature
###############################################################################################################

import os
import sys

#############################################################
current_abs_path = os.path.abspath(os.path.dirname(__file__))
print('NucleosomeOccupancyAnalysis_SPMs_SignatureBased.py current_abs_path:%s' %(current_abs_path))
#############################################################


commonsPath = os.path.join(current_abs_path, '..','commons')
sys.path.append(commonsPath)

from SigProfilerTopography.source.commons.TopographyCommons import *

########################################################################################
#This method is specific to this python file.
def fillSignalArrayAndCountArray(combinedList):
    chrBased_mutation_df = combinedList[0]
    chrBased_nucleosome_df_split = combinedList[1]
    chrLong = combinedList[2]
    signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = combinedList[3]

    ###########################################################################
    # Initialize chrBased_nucleosome_split_array with np array of zeros
    # Fill the nucleosome split array using chrBased_nucleosome_df_split
    chrBased_nucleosome_split_array = np.zeros(MAXIMUM_CHROMOSOME_LENGTH, dtype=np.float32)
    chrBased_nucleosome_df_split.apply(fillNucleosomeSignalArray, nucleosome_array=chrBased_nucleosome_split_array, axis=1)
    ###########################################################################

    ###########################################################################
    print('#### debug starts ####')
    print('chr:%s -- '
            'size(chrBased_mutation_df):%d size(chrBased_mutation_df):%f GB --- '
            'size(chrBased_nucleosome_df_split):%d size(chrBased_nucleosome_df_split):%f GB  --- '
            'size(chrBased_nucleosome_split_array):%d size(chrBased_nucleosome_split_array):%f GB' %(chrLong,
            sys.getsizeof(chrBased_mutation_df),sys.getsizeof(chrBased_mutation_df)/GIGABYTE_IN_BYTES,
            sys.getsizeof(chrBased_nucleosome_df_split),sys.getsizeof(chrBased_nucleosome_df_split) / GIGABYTE_IN_BYTES,
            sys.getsizeof(chrBased_nucleosome_split_array), sys.getsizeof(chrBased_nucleosome_split_array)/ GIGABYTE_IN_BYTES ))
    print('#### debug ends ####')
    ###########################################################################

    ##############################################
    # Initialization for each signature
    signature2SignalArrayDict = {}
    signature2CountArrayDict = {}

    for signature in signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
        signature2SignalArrayDict[signature] = np.zeros(windowSize)
        signature2CountArrayDict[signature] = np.zeros(windowSize, dtype=int)
    ##############################################

    ##############################################
    # Initialization for all single point mutations for split
    allSinglePointMutationsSignalArray =  np.zeros(windowSize)
    allSinglePointMutationsCountArray = np.zeros(windowSize, dtype=int)
    ##############################################

    #Initialzie the list, you will return this list
    split_chrBased_SignalArrayAndCountArray_DictionaryList = []

    # Fill signature2SignalArrayDict and signature2CountArrayDict using chrBased_mutation_df and chrBased_nucleosome_split_array
    chrBased_mutation_df.apply(fillSplitBasedSignalArrayAndCountArrayForSPMs, nucleosome_array=chrBased_nucleosome_split_array, signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict= signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict, signature2SignalArrayDict=signature2SignalArrayDict,signature2CountArrayDict=signature2CountArrayDict, allSinglePointMutationsSignalArray=allSinglePointMutationsSignalArray,allSinglePointMutationsCountArray=allSinglePointMutationsCountArray, axis=1)

    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(signature2SignalArrayDict)
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(signature2CountArrayDict)

    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(allSinglePointMutationsSignalArray)
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(allSinglePointMutationsCountArray)

    return split_chrBased_SignalArrayAndCountArray_DictionaryList
########################################################################################


########################################################################################
#This method is specific to this python file.
def fillSignalArrayAndCountArrayWithExtraSampleBased(combinedList):
    chrBased_mutation_df = combinedList[0]
    chrBased_nucleosome_df_split = combinedList[1]
    chrLong = combinedList[2]
    maximum_chrom_size = combinedList[3]
    signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = combinedList[4]
    samplesWithAtLeast10KMutations2NumberofMutationsDict = combinedList[5]
    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = combinedList[6]

    ###########################################################################
    # Initialize chrBased_nucleosome_split_array with np array of zeros
    # Fill the nucleosome split array using chrBased_nucleosome_df_split
    chrBased_nucleosome_split_array = np.zeros(MAXIMUM_CHROMOSOME_LENGTH, dtype=np.float32)
    #In fact here, we are increasing our memory usage
    chrBased_nucleosome_df_split.apply(fillNucleosomeSignalArray, nucleosome_array=chrBased_nucleosome_split_array, axis=1)
    ###########################################################################

    ###########################################################################
    print('#### debug starts ####')
    print('chr:%s -- '
            'size(chrBased_mutation_df):%d size(chrBased_mutation_df):%f GB --- '
            'size(chrBased_nucleosome_df_split):%d size(chrBased_nucleosome_df_split):%f GB  --- '
            'size(chrBased_nucleosome_split_array):%d size(chrBased_nucleosome_split_array):%f GB' %(chrLong,
            sys.getsizeof(chrBased_mutation_df),sys.getsizeof(chrBased_mutation_df)/GIGABYTE_IN_BYTES,
            sys.getsizeof(chrBased_nucleosome_df_split),sys.getsizeof(chrBased_nucleosome_df_split) / GIGABYTE_IN_BYTES,
            sys.getsizeof(chrBased_nucleosome_split_array), sys.getsizeof(chrBased_nucleosome_split_array)/GIGABYTE_IN_BYTES ))
    print('#### debug ends ####')
    ###########################################################################

    ##############################################
    # Initialization for each signature
    signature2SignalArrayDict = {}
    signature2CountArrayDict = {}

    # for signature in signatures:
    for signature in signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
        signature2SignalArrayDict[signature] = np.zeros(windowSize)
        signature2CountArrayDict[signature] = np.zeros(windowSize, dtype=int)
    ##############################################

    ##############################################
    # Initialization for all single point mutations for split
    allSinglePointMutationsSignalArray =  np.zeros(windowSize)
    allSinglePointMutationsCountArray = np.zeros(windowSize, dtype=int)
    ##############################################

    ########################################################
    # Initialization2 for sample based for each split
    sample2Signature2SignalArrayDict = {}
    sample2Signature2CountArrayDict = {}

    sample2AllSinglePointMutationsSignalArrayDict = {}
    sample2AllSinglePointMutationsCountArrayDict = {}

    for sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
        sample2AllSinglePointMutationsSignalArrayDict[sample] = np.zeros(windowSize)
        sample2AllSinglePointMutationsCountArrayDict[sample] = np.zeros(windowSize, dtype=int)

        sample2Signature2SignalArrayDict[sample] = {}
        sample2Signature2CountArrayDict[sample] ={}

        for signature in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict[sample]:
            sample2Signature2SignalArrayDict[sample][signature] = np.zeros(windowSize)
            sample2Signature2CountArrayDict[sample][signature] = np.zeros(windowSize, dtype=int)

    # Sep 18, 2018 ends
    ########################################################

    #Initialzie the list, you will return this list
    split_chrBased_SignalArrayAndCountArray_DictionaryList = []

    # Fill signature2SignalArrayDict and signature2CountArrayDict using chrBased_mutation_df and chrBased_nucleosome_split_array
    chrBased_mutation_df.apply(fillSplitBasedSignalArrayAndCountArrayForSPMsWithExtraSampleBased,
                               nucleosome_array=chrBased_nucleosome_split_array,
                               maximum_chrom_size=maximum_chrom_size,
                               signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict=signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
                               samplesWithAtLeast10KMutations2NumberofMutationsDict=samplesWithAtLeast10KMutations2NumberofMutationsDict,
                               sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict= sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
                               signature2SignalArrayDict=signature2SignalArrayDict,
                               signature2CountArrayDict=signature2CountArrayDict,
                               allSinglePointMutationsSignalArray=allSinglePointMutationsSignalArray,
                               allSinglePointMutationsCountArray=allSinglePointMutationsCountArray,
                               sample2Signature2SignalArrayDict= sample2Signature2SignalArrayDict,
                               sample2Signature2CountArrayDict=sample2Signature2CountArrayDict,
                               sample2AllSinglePointMutationsSignalArrayDict=sample2AllSinglePointMutationsSignalArrayDict,
                               sample2AllSinglePointMutationsCountArrayDict=sample2AllSinglePointMutationsCountArrayDict,
                               axis=1)

    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(signature2SignalArrayDict)
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(signature2CountArrayDict)

    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(allSinglePointMutationsSignalArray)
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(allSinglePointMutationsCountArray)

    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(sample2Signature2SignalArrayDict)
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(sample2Signature2CountArrayDict)

    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(sample2AllSinglePointMutationsSignalArrayDict)
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(sample2AllSinglePointMutationsCountArrayDict)

    return split_chrBased_SignalArrayAndCountArray_DictionaryList
########################################################################################



########################################################################################
#TODO let's consider all the signaturewWithAtLeast10KEligibleMutations
#This function is specific to this python file
def accumulateSplitArrays(signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,allSplits_chrBased_SignalArrayAndCountArray_DictionaryList):

    ##############################################
    # Fill these dictionaries
    signature2AccumulatedSplitsChrBasedSignalArrayDict = {}
    signature2AccumulatedSplitsChrBasedCountArrayDict = {}

    for signature in signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
        signature2AccumulatedSplitsChrBasedSignalArrayDict[signature] = np.zeros(windowSize)
        signature2AccumulatedSplitsChrBasedCountArrayDict[signature] = np.zeros(windowSize,dtype=int)
    ##############################################

    ##############################################
    # Fill these arrays
    allSingleMutationsAccumulatedSplitsSignalArray = np.zeros(windowSize)
    allSingleMutationsAccumulatedSplitsCountArray = np.zeros(windowSize, dtype=int)
    ##############################################

    ##############################################
    for split_chrBased_SignalArrayAndCountArray_DictionaryList in allSplits_chrBased_SignalArrayAndCountArray_DictionaryList:
        signature2SplitSignalArrayDict = split_chrBased_SignalArrayAndCountArray_DictionaryList[0]
        signature2SplitCountArrayDict = split_chrBased_SignalArrayAndCountArray_DictionaryList[1]

        allSingePointMutationsSplitSignalArray = split_chrBased_SignalArrayAndCountArray_DictionaryList[2]
        allSingePointMutationsSplitCountArray = split_chrBased_SignalArrayAndCountArray_DictionaryList[3]

        accumulateSignatureBasedArrays(signature2AccumulatedSplitsChrBasedSignalArrayDict,signature2SplitSignalArrayDict)
        accumulateSignatureBasedArrays(signature2AccumulatedSplitsChrBasedCountArrayDict,signature2SplitCountArrayDict)

        allSingleMutationsAccumulatedSplitsSignalArray += allSingePointMutationsSplitSignalArray
        allSingleMutationsAccumulatedSplitsCountArray += allSingePointMutationsSplitCountArray
    ##############################################

    return signature2AccumulatedSplitsChrBasedSignalArrayDict, signature2AccumulatedSplitsChrBasedCountArrayDict, allSingleMutationsAccumulatedSplitsSignalArray, allSingleMutationsAccumulatedSplitsCountArray
########################################################################################


########################################################################################
#This function is specific to this python file
def accumulateSplitArraysWithExtraSampleBased(
        signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
        sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
        allSplits_chrBased_SignalArrayAndCountArray_DictionaryList):

    ##############################################
    #Initialization
    # Fill these dictionaries
    signature2AccumulatedSplitsChrBasedSignalArrayDict = {}
    signature2AccumulatedSplitsChrBasedCountArrayDict = {}

    for signature in signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
        signature2AccumulatedSplitsChrBasedSignalArrayDict[signature] = np.zeros(windowSize)
        signature2AccumulatedSplitsChrBasedCountArrayDict[signature] = np.zeros(windowSize,dtype=int)
    ##############################################

    ##############################################
    #Initialization
    # Fill these arrays

    allSingleMutationsAccumulatedSplitsSignalArray = np.zeros(windowSize)
    allSingleMutationsAccumulatedSplitsCountArray = np.zeros(windowSize, dtype=int)
    ##############################################

    ##############################################
    # Initialization3 for sample based for all splits coming from the same chromosome
    #Fill sampleBased signatureBased dictionaries
    sample2Signature2AccumulatedSplitsChrBasedSignalArrayDict = {}
    sample2Signature2AccumulatedSplitsChrBasedCountArrayDict = {}

    sample2AllSinglePointMutationsAccumulatedSplitsChrBasedSignalArrayDict= {}
    sample2AllSinglePointMutationsAccumulatedSplitsChrBasedCountArrayDict= {}

    for sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
        sample2Signature2AccumulatedSplitsChrBasedSignalArrayDict[sample] ={}
        sample2Signature2AccumulatedSplitsChrBasedCountArrayDict[sample] ={}

        sample2AllSinglePointMutationsAccumulatedSplitsChrBasedSignalArrayDict[sample] = np.zeros(windowSize)
        sample2AllSinglePointMutationsAccumulatedSplitsChrBasedCountArrayDict[sample] = np.zeros(windowSize, dtype=int)

        for signature in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict[sample]:
            sample2Signature2AccumulatedSplitsChrBasedSignalArrayDict[sample][signature] = np.zeros(windowSize)
            sample2Signature2AccumulatedSplitsChrBasedCountArrayDict[sample][signature] = np.zeros(windowSize, dtype=int)
    ##############################################


    ##############################################
    for split_chrBased_SignalArrayAndCountArray_DictionaryList in allSplits_chrBased_SignalArrayAndCountArray_DictionaryList:
        signature2SplitSignalArrayDict = split_chrBased_SignalArrayAndCountArray_DictionaryList[0]
        signature2SplitCountArrayDict = split_chrBased_SignalArrayAndCountArray_DictionaryList[1]

        allSingePointMutationsSplitSignalArray = split_chrBased_SignalArrayAndCountArray_DictionaryList[2]
        allSingePointMutationsSplitCountArray = split_chrBased_SignalArrayAndCountArray_DictionaryList[3]

        sample2Signature2SplitSignalArrayDict = split_chrBased_SignalArrayAndCountArray_DictionaryList[4]
        sample2Signature2SplitCountArrayDict = split_chrBased_SignalArrayAndCountArray_DictionaryList[5]

        sample2AllSinglePointMutationsSplitSignalArrayDict = split_chrBased_SignalArrayAndCountArray_DictionaryList[6]
        sample2AllSinglePointMutationsSplitCountArrayDict = split_chrBased_SignalArrayAndCountArray_DictionaryList[7]
        #################################################################


        ###############################################################
        accumulateSignatureBasedArrays(signature2AccumulatedSplitsChrBasedSignalArrayDict,signature2SplitSignalArrayDict)
        accumulateSignatureBasedArrays(signature2AccumulatedSplitsChrBasedCountArrayDict,signature2SplitCountArrayDict)

        allSingleMutationsAccumulatedSplitsSignalArray += allSingePointMutationsSplitSignalArray
        allSingleMutationsAccumulatedSplitsCountArray += allSingePointMutationsSplitCountArray
        ###############################################################

        ###############################################################
        accumulateSampleBasedSignatureBasedArrays(sample2Signature2AccumulatedSplitsChrBasedSignalArrayDict,sample2Signature2SplitSignalArrayDict)
        accumulateSampleBasedSignatureBasedArrays(sample2Signature2AccumulatedSplitsChrBasedCountArrayDict,sample2Signature2SplitCountArrayDict)

        accumulateSampleBasedArrays(sample2AllSinglePointMutationsAccumulatedSplitsChrBasedSignalArrayDict,sample2AllSinglePointMutationsSplitSignalArrayDict)
        accumulateSampleBasedArrays(sample2AllSinglePointMutationsAccumulatedSplitsChrBasedCountArrayDict,sample2AllSinglePointMutationsSplitCountArrayDict)
        ###############################################################


    ##############################################

    return  signature2AccumulatedSplitsChrBasedSignalArrayDict, \
            signature2AccumulatedSplitsChrBasedCountArrayDict, \
            allSingleMutationsAccumulatedSplitsSignalArray, \
            allSingleMutationsAccumulatedSplitsCountArray,\
            sample2Signature2AccumulatedSplitsChrBasedSignalArrayDict,\
            sample2Signature2AccumulatedSplitsChrBasedCountArrayDict,\
            sample2AllSinglePointMutationsAccumulatedSplitsChrBasedSignalArrayDict,\
            sample2AllSinglePointMutationsAccumulatedSplitsChrBasedCountArrayDict
########################################################################################


########################################################################################
#main function
def nucleosomeOccupancyAnalysis_SPMs_SignatureBased(jobname, singlePointMutationsFilename, nucleosomeFilename):

    withExtraSampleBasedAnalysis = True

    ##########################################################################################
    GRCh37ChromSizesDict = {}

    GRCh37ChromSizesDictPath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB,UCSCGENOME,GRCh37ChromSizesDictFilename)

    if (os.path.exists(GRCh37ChromSizesDictPath)):
        GRCh37ChromSizesDict = readDictionary(GRCh37ChromSizesDictPath)
    ##########################################################################################


    ##########################################################################################
    #These lists are filled in py in readMutations function
    #Load Chrnames in nucleosome file
    chrNamesInNucleosomeList = []
    chrNamesInNucleosomeFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME, ChrNamesInNucleosomesFilename)
    if (os.path.exists(chrNamesInNucleosomeFilePath)):
        chrNamesInNucleosomeArray = np.loadtxt(chrNamesInNucleosomeFilePath,dtype=str, delimiter='\t')
        chrNamesInNucleosomeList = list(chrNamesInNucleosomeArray)

    #Load samplesWithAtLeast10KMutations2NumberofMutationsDict
    samplesWithAtLeast10KMutations2NumberofMutationsDict = {}

    samplesWithAtLeast10KMutations2NumberofMutationsDictFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT,jobname,DATA,SamplesWithAtLeast10KMutations2NumberofMutationsDictFilename)

    if (os.path.exists(samplesWithAtLeast10KMutations2NumberofMutationsDictFilePath)):
        samplesWithAtLeast10KMutations2NumberofMutationsDict = readDictionary(samplesWithAtLeast10KMutations2NumberofMutationsDictFilePath)

    #Load signaturesWithAtLeast10KEligibleMutations
    signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = {}

    SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT,jobname,DATA,SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilename)

    if (os.path.exists(SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath)):
        signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = readDictionary(SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath)

    #Load sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict
    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict= {}

    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath = os.path.join(current_abs_path,ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT,jobname,DATA,Sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilename)

    if (os.path.exists(sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath)):
        sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = readDictionary(sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath)
    ##########################################################################################

    ########################################################################################
    # Initialization for each signature
    signature2AccumulatedAllChromsSignalArrayDict = {}
    signature2AccumulatedAllChromsCountArrayDict = {}

    # for signature in signatureList:
    for signature in signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
        signature2AccumulatedAllChromsSignalArrayDict[signature] = np.zeros(windowSize)
        signature2AccumulatedAllChromsCountArrayDict[signature] = np.zeros(windowSize, dtype=int)
    ########################################################################################


    ########################################################################################
    # Initialization for single point mutations
    allSinglePointMutationsAccumulatedAllChromsSignalArray = np.zeros(windowSize)
    allSinglePointMutationsAccumulatedAllChromsCountArray = np.zeros(windowSize, dtype=int)
    ########################################################################################

    ########################################################################################
    # Initialization1 for sample based for all chromosomes
    sample2Signature2AccumulatedAllChromsSignalArrayDict = {}
    sample2Signature2AccumulatedAllChromsCountArrayDict = {}

    sample2AllSinglePointMutationsAccumulatedAllChromsSignalArrayDict = {}
    sample2AllSinglePointMutationsAccumulatedAllChromsCountArrayDict = {}

    for sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
        sample2Signature2AccumulatedAllChromsSignalArrayDict[sample] = {}
        sample2Signature2AccumulatedAllChromsCountArrayDict[sample]= {}

        for signature in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict[sample]:
            sample2Signature2AccumulatedAllChromsSignalArrayDict[sample][signature] =  np.zeros(windowSize)
            sample2Signature2AccumulatedAllChromsCountArrayDict[sample][signature] = np.zeros(windowSize, dtype=int)

        sample2AllSinglePointMutationsAccumulatedAllChromsSignalArrayDict[sample] = np.zeros(windowSize)
        sample2AllSinglePointMutationsAccumulatedAllChromsCountArrayDict[sample] = np.zeros(windowSize, dtype=int)
    ########################################################################################

    ######################################################
    print('chrNamesInNucleosomeList')
    print(chrNamesInNucleosomeList)
    # short chromosomes names e.g.: chrNamesList = ['1','10','11','12','13', '14', '15', '16', '17', '18', '19', '2', '20', '21', '22', '3', '4', '5', '6', '7', '8', '9', 'MT', 'X', 'Y']
    ######################################################

    ######################################################
    #Better practice
    numofProcesses = multiprocessing.cpu_count()
    print('multiprocessing.cpu_count() --> number of processes: %d' % numofProcesses)
    pool = multiprocessing.Pool(numofProcesses)
    ######################################################

    ###################################################################################
    ##################  For each chromsome sequential starts ##########################
    ###################################################################################
    for chrLong in chrNamesInNucleosomeList:

        maximum_chrom_size = GRCh37ChromSizesDict[chrLong]

        #FIRST READ CHRBASED NUCLEOSOME OCCUPANCY
        #TODO read chrbased nucleosome data as a numpy array
        chrBased_nucleosome_df = readChrBasedNuclesomeDF(chrLong,nucleosomeFilename)
        print('chromosome %s  -- chrBased_nucleosome_df: %d in bytes %f in GB' % (chrLong,sys.getsizeof(chrBased_nucleosome_df), sys.getsizeof(chrBased_nucleosome_df)/GIGABYTE_IN_BYTES))

        #TODO: This is specific to out data right now
        #Nucleosomes have chrM
        #SinglePointMutations and Indels have chrMT
        if (chrLong=='chrM'):
            chrLong='chrMT'

        #THEN READ CHRBASED MUTATION
        chrBased_mutation_df=readChrBasedMutationDF(jobname,chrLong,singlePointMutationsFilename)
        print('chromosome %s  -- chrBased_mutation_df: %d in Bytes %f in GigaBytes' %(chrLong,sys.getsizeof(chrBased_mutation_df),sys.getsizeof(chrBased_mutation_df)/GIGABYTE_IN_BYTES))

        if ((chrBased_mutation_df is not None) and (chrBased_nucleosome_df is not None)):

            #################################################
            #Question: Shall I split chrBased_mutation_df or chrBased_nucleosome_df?
            #Answer: Since chrBased_nucleosome_df is bigger than chrBased_mutation_df, therefore split chrBased_nucleosome_df

            # Question: Shall I do fill nucleosome array here or in each core separately?
            #for all rows of nucleosome_df update the np.array of zero using slicing
            #Answer: Let's do it in each core separately
            # chrBased_nucleosome_df.apply(fillNucleosomeArray, nucleosome_array=chrBased_nucleosome_array, axis=1)
            #################################################

            #################################################
            # split chrBased_mutation_df into number of available cores
            # for each split_chrBased_nucleosome_df first fill nucleosome_array
            # then fill split_signal and split_count array
            # accumulate all the split_signal and split_count arrays in chrBased signal array and count arrays

            #TODO: later on split the one that is bigger whether it is chrBasedMutationDF or chrBasedNucleosomeDF
            #whether you  split nucleosome or mutations does not matter
            # you have to accumulate in the same manner.
            # chrBased_mutation_df_splits = np.array_split(chrBased_skin_df, numofProcesses)
            chrBased_nucleosome_df_splits = np.array_split(chrBased_nucleosome_df,numofProcesses)

            #Prepare the poolInputList
            poolInputList = []

            if (withExtraSampleBasedAnalysis):

                ###########################################################################
                ##############   SPMs with extra sampleBased  analysis starts #############
                ###########################################################################

                # For each split chrBased mutation df do it in parallel
                for chrBased_nucleosome_df_split in chrBased_nucleosome_df_splits:
                    inputList = []

                    # Question: Does each split has header?
                    # Answer: Yes
                    inputList.append(chrBased_mutation_df) #Same chrBased_mutation_df
                    inputList.append(chrBased_nucleosome_df_split)  #Different chrBased_nucleosome_df split
                    inputList.append(chrLong) #same chrLong
                    inputList.append(maximum_chrom_size)
                    inputList.append(signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict) #same signature
                    inputList.append(samplesWithAtLeast10KMutations2NumberofMutationsDict) #Sep 18, 2018
                    inputList.append(sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict) #Sep 19, 2018
                    poolInputList.append(inputList)
                #########################################################################

                #########################################################################
                allSplits_chrBased_SignalArrayAndCountArray_DictionaryList = pool.map(fillSignalArrayAndCountArrayWithExtraSampleBased,poolInputList)

                #Accumulate the results coming from each split
                #allSplits_chrBased_SignalArrayAndCountArray_DictionaryList is a list of lists.
                #Each list consist of two dictionaries
                signature2AccumulatedSplitsChrBasedSignalArrayDict, \
                signature2AccumulatedSplitsChrBasedCountArrayDict, \
                allSingleMutationsAccumulatedSplitsSignalArray, \
                allSingleMutationsAccumulatedSplitsCountArray, \
                sample2Signature2AccumulatedSplitsChrBasedSignalArrayDict, \
                sample2Signature2AccumulatedSplitsChrBasedCountArrayDict, \
                sample2AllSinglePointMutationsAccumulatedSplitsChrBasedSignalArrayDict, \
                sample2AllSinglePointMutationsAccumulatedSplitsChrBasedCountArrayDict = accumulateSplitArraysWithExtraSampleBased(signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,allSplits_chrBased_SignalArrayAndCountArray_DictionaryList)

                #################################################

                #################################################
                ##Accumulate the results coming from each chromosome
                accumulateSignatureBasedArrays(signature2AccumulatedAllChromsSignalArrayDict,signature2AccumulatedSplitsChrBasedSignalArrayDict)
                accumulateSignatureBasedArrays(signature2AccumulatedAllChromsCountArrayDict,signature2AccumulatedSplitsChrBasedCountArrayDict)

                ##Accumulate the results coming from each chromosome
                allSinglePointMutationsAccumulatedAllChromsSignalArray += allSingleMutationsAccumulatedSplitsSignalArray
                allSinglePointMutationsAccumulatedAllChromsCountArray += allSingleMutationsAccumulatedSplitsCountArray
                #################################################

                ###############################################################
                ##Accumulate the results coming from each chromosome
                accumulateSampleBasedSignatureBasedArrays(sample2Signature2AccumulatedAllChromsSignalArrayDict,sample2Signature2AccumulatedSplitsChrBasedSignalArrayDict)
                accumulateSampleBasedSignatureBasedArrays(sample2Signature2AccumulatedAllChromsCountArrayDict,sample2Signature2AccumulatedSplitsChrBasedCountArrayDict)
                accumulateSampleBasedArrays(sample2AllSinglePointMutationsAccumulatedAllChromsSignalArrayDict,sample2AllSinglePointMutationsAccumulatedSplitsChrBasedSignalArrayDict)
                accumulateSampleBasedArrays(sample2AllSinglePointMutationsAccumulatedAllChromsCountArrayDict,sample2AllSinglePointMutationsAccumulatedSplitsChrBasedCountArrayDict)
                ###############################################################

                ###########################################################################
                ##############   SPMs with extra sampleBased  analysis ends ###############
                ###########################################################################

            else:

                ###########################################################################
                ##############   SPMs analysis starts #####################################
                ###########################################################################

                #########################################################################
                # For each split chrBased mutation df do it in parallel
                for chrBased_nucleosome_df_split in chrBased_nucleosome_df_splits:
                    inputList = []

                    # Question: Does each split has header?
                    # Answer: Yes
                    inputList.append(chrBased_mutation_df)  # Same chrBased_mutation_df
                    inputList.append(chrBased_nucleosome_df_split)  # Different chrBased_nucleosome_df split
                    inputList.append(chrLong)  # same chrLong
                    inputList.append(signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict)  # same signature
                    poolInputList.append(inputList)
                #########################################################################

                #########################################################################
                allSplits_chrBased_SignalArrayAndCountArray_DictionaryList = pool.map(fillSignalArrayAndCountArray,poolInputList)

                #################################################
                # Accumulate the results coming from each split
                # allSplits_chrBased_SignalArrayAndCountArray_DictionaryList is a list of lists.
                # Each list consist of two dictionaries
                signature2AccumulatedSplitsChrBasedSignalArrayDict, \
                signature2AccumulatedSplitsChrBasedCountArrayDict, \
                allSingleMutationsAccumulatedSplitsSignalArray, \
                allSingleMutationsAccumulatedSplitsCountArray = accumulateSplitArrays(signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict, allSplits_chrBased_SignalArrayAndCountArray_DictionaryList)
                #################################################

                #################################################
                ##Accumulate the results coming from each chromosome
                accumulateSignatureBasedArrays(signature2AccumulatedAllChromsSignalArrayDict,signature2AccumulatedSplitsChrBasedSignalArrayDict)
                accumulateSignatureBasedArrays(signature2AccumulatedAllChromsCountArrayDict,signature2AccumulatedSplitsChrBasedCountArrayDict)

                ##Accumulate the results coming from each chromosome
                allSinglePointMutationsAccumulatedAllChromsSignalArray += allSingleMutationsAccumulatedSplitsSignalArray
                allSinglePointMutationsAccumulatedAllChromsCountArray += allSingleMutationsAccumulatedSplitsCountArray
                #################################################

                ###########################################################################
                ##############   SPMs analysis ends #######################################
                ###########################################################################


    ###################################################################################
    ##################  For each chromsome sequential ends ############################
    ###################################################################################

    ################################
    pool.close()
    pool.join()
    ################################


    ####################################################################################
    ########### Write All Single Mutations Average Nucleosome Occupancy starts #########
    ####################################################################################
    writeAverageNucleosomeOccupancyFiles(allSinglePointMutationsAccumulatedAllChromsSignalArray,allSinglePointMutationsAccumulatedAllChromsCountArray, jobname,AGGREGATEDSUBSTITUTIONS)
    ####################################################################################
    ########### Write All Single Mutations Average Nucleosome Occupancy ends ###########
    ####################################################################################


    ####################################################################################
    ############## Write Signature Based Average Nucleosome Occupancy starts ###########
    ####################################################################################
    writeSignatureBasedAverageNucleosomeOccupancyFiles(signature2AccumulatedAllChromsSignalArrayDict,signature2AccumulatedAllChromsCountArrayDict,jobname)
    ####################################################################################
    ############## Write Signature Based Average Nucleosome Occupancy ends #############
    ####################################################################################


    ####################################################################################
    ############################# Write Sample Based starts ############################
    ####################################################################################

    if (withExtraSampleBasedAnalysis):
        ####################################################################################
        writeSampleBasedSignatureBasedAverageNucleosomeOccupancyFiles(sample2Signature2AccumulatedAllChromsSignalArrayDict,sample2Signature2AccumulatedAllChromsCountArrayDict,jobname)
        ####################################################################################

        ####################################################################################
        writeSampleBasedAverageNucleosomeOccupancyFiles(sample2AllSinglePointMutationsAccumulatedAllChromsSignalArrayDict,sample2AllSinglePointMutationsAccumulatedAllChromsCountArrayDict,jobname,AGGREGATEDSUBSTITUTIONS)
        ####################################################################################

    ####################################################################################
    ############################# Write Sample Based ends ##############################
    ####################################################################################

########################################################################################
