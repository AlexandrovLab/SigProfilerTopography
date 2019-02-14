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
#   for all indels
###############################################################################################################


import os
import sys


#############################################################
current_abs_path = os.path.abspath(os.path.dirname(__file__))
print('NucleosomeOccupancyAnalysis_SPMs_SignatureBased_Indels.py current_abs_path:%s' %(current_abs_path))
#############################################################


commonsPath = os.path.join(current_abs_path, '..','commons')
sys.path.append(commonsPath)

from SigProfilerTopography.source.commons.TopographyCommons import *


########################################################################################
def fillSignalArrayAndCountArrayWithExtraSampleBased(inputList):
    chrBased_spms_df = inputList[0]
    chrBased_indels_df = inputList[1]
    chrBased_nucleosome_df_split = inputList[2]
    chrLong = inputList[3]
    signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = inputList[4]
    samplesWithAtLeast10KMutations2NumberofMutationsDict = inputList[5]
    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = inputList[6]


    ###########################################################################
    # Initialize chrBased_nucleosome_split_array with np array of zeros
    # Fill the nucleosome split array using chrBased_nucleosome_df_split
    chrBased_nucleosome_split_array = np.zeros(MAXIMUM_CHROMOSOME_LENGTH, dtype=np.float32)
    chrBased_nucleosome_df_split.apply(fillNucleosomeSignalArray, nucleosome_array=chrBased_nucleosome_split_array, axis=1)
    ###########################################################################


    ###############################################################################
    ################### Initialization2 for each Split starts #####################
    ###############################################################################

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

    ##############################################
    # Initialization for all indels for split
    allIndelsSignalArray = np.zeros(windowSize)
    allIndelsCountArray = np.zeros(windowSize, dtype=int)
    ##############################################

    ##############################################
    # Initialization sample based

    # Initialization2 for sample based for each split
    sample2Signature2SignalArrayDict = {}
    sample2Signature2CountArrayDict = {}

    sample2AllSinglePointMutationsSignalArrayDict = {}
    sample2AllSinglePointMutationsCountArrayDict = {}

    sample2AllIndelsSignalArrayDict = {}
    sample2AllIndelsCountArrayDict = {}

    for sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
        sample2Signature2SignalArrayDict[sample] = {}
        sample2Signature2CountArrayDict[sample] = {}

        sample2AllSinglePointMutationsSignalArrayDict[sample] = np.zeros(windowSize)
        sample2AllSinglePointMutationsCountArrayDict[sample] = np.zeros(windowSize, dtype=int)

        sample2AllIndelsSignalArrayDict[sample] = np.zeros(windowSize)
        sample2AllIndelsCountArrayDict[sample] = np.zeros(windowSize, dtype=int)

        for signature in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict[sample]:
            sample2Signature2SignalArrayDict[sample][signature] = np.zeros(windowSize)
            sample2Signature2CountArrayDict[sample][signature] = np.zeros(windowSize, dtype=int)
    ##############################################

    ###############################################################################
    ################### Initialization2 for each split ends #######################
    ###############################################################################


    ###############################################################################
    ################### Fill signal and count array starts ########################
    ###############################################################################
    #Fill for single point mutations
    chrBased_spms_df.apply(fillSplitBasedSignalArrayAndCountArrayForSPMsWithExtraSampleBased,
                           nucleosome_array=chrBased_nucleosome_split_array,
                           signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
                           samplesWithAtLeast10KMutations2NumberofMutationsDict=samplesWithAtLeast10KMutations2NumberofMutationsDict,
                           sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict=sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
                           signature2SignalArrayDict=signature2SignalArrayDict,
                           signature2CountArrayDict=signature2CountArrayDict,
                           allSinglePointMutationsSignalArray=allSinglePointMutationsSignalArray,
                           allSinglePointMutationsCountArray=allSinglePointMutationsCountArray,
                           sample2Signature2SignalArrayDict=sample2Signature2SignalArrayDict,
                           sample2Signature2CountArrayDict=sample2Signature2CountArrayDict,
                           sample2AllSinglePointMutationsSignalArrayDict=sample2AllSinglePointMutationsSignalArrayDict,
                           sample2AllSinglePointMutationsCountArrayDict=sample2AllSinglePointMutationsCountArrayDict,
                           axis=1)

    #Fill for indels
    if (chrBased_indels_df is not None):
        chrBased_indels_df.apply(fillSplitBasedSignalArrayAndCountArrayForIndelsWithExtraSampleBased,
                                 nucleosome_array=chrBased_nucleosome_split_array,
                                 sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict=sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
                                 allIndelsSignalArray=allIndelsSignalArray,
                                 allIndelsCountArray=allIndelsCountArray,
                                 sample2AllIndelsSignalArrayDict=sample2AllIndelsSignalArrayDict,
                                 sample2AllIndelsCountArrayDict =sample2AllIndelsCountArrayDict,
                                 axis=1)

    ###############################################################################
    ################### Fill signal and count array ends ##########################
    ###############################################################################

    # Initialzie the list, you will return this list
    split_chrBased_SignalArrayAndCountArray_DictionaryList = []

    # Append signature arrays
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(signature2SignalArrayDict)
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(signature2CountArrayDict)

    # Append single point mutations arrays
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(allSinglePointMutationsSignalArray)
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(allSinglePointMutationsCountArray)

    # Append indels arrays
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(allIndelsSignalArray)
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(allIndelsCountArray)

    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(sample2Signature2SignalArrayDict)
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(sample2Signature2CountArrayDict)

    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(sample2AllSinglePointMutationsSignalArrayDict)
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(sample2AllSinglePointMutationsCountArrayDict)

    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(sample2AllIndelsSignalArrayDict)
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(sample2AllIndelsCountArrayDict)


    return split_chrBased_SignalArrayAndCountArray_DictionaryList
########################################################################################

########################################################################################
#This method is specific to this python file.
def fillNucleosomeSignalArrayAndCountArray(inputList):
    chrBased_spms_df = inputList[0]
    chrBased_indels_df = inputList[1]
    chrBased_nucleosome_df_split = inputList[2]
    chrLong = inputList[3]
    signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = inputList[4]

    ###########################################################################
    # Initialize chrBased_nucleosome_split_array with np array of zeros
    # Fill the nucleosome split array using chrBased_nucleosome_df_split
    chrBased_nucleosome_split_array = np.zeros(MAXIMUM_CHROMOSOME_LENGTH, dtype=np.float32)
    chrBased_nucleosome_df_split.apply(fillNucleosomeSignalArray, nucleosome_array=chrBased_nucleosome_split_array, axis=1)
    ###########################################################################

    ###########################################################################
    print('#### debug starts ####')
    print('chr:%s -- '
            'size(chrBased_spms_df):%d size(chrBased_spms_df):%f GB --- '
            'size(chrBased_indels_df):%d size(chrBased_indels_df):%f GB --- '
            'size(chrBased_nucleosome_df_split):%d size(chrBased_nucleosome_df_split):%f GB  --- '
            'size(chrBased_nucleosome_split_array):%d size(chrBased_nucleosome_split_array):%f GB' %(chrLong,
            sys.getsizeof(chrBased_spms_df),sys.getsizeof(chrBased_spms_df)/GIGABYTE_IN_BYTES,
            sys.getsizeof(chrBased_indels_df),sys.getsizeof(chrBased_indels_df)/GIGABYTE_IN_BYTES,
            sys.getsizeof(chrBased_nucleosome_df_split),sys.getsizeof(chrBased_nucleosome_df_split) / GIGABYTE_IN_BYTES,
            sys.getsizeof(chrBased_nucleosome_split_array), sys.getsizeof(chrBased_nucleosome_split_array)/GIGABYTE_IN_BYTES ))
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

    ##############################################
    # Initialization for all indels for split
    allIndelsSignalArray = np.zeros(windowSize)
    allIndelsCountArray = np.zeros(windowSize, dtype=int)
    ##############################################

    #Initialzie the list, you will return this list
    split_chrBased_SignalArrayAndCountArray_DictionaryList = []

    #Fill for single point mutations
    # Fill signature2SignalArrayDict and signature2CountArrayDict using chrBased_spms_df and chrBased_nucleosome_split_array
    chrBased_spms_df.apply(fillSplitBasedSignalArrayAndCountArrayForSPMs, nucleosome_array=chrBased_nucleosome_split_array, signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict, signature2SignalArrayDict=signature2SignalArrayDict,signature2CountArrayDict=signature2CountArrayDict, allSinglePointMutationsSignalArray=allSinglePointMutationsSignalArray,allSinglePointMutationsCountArray=allSinglePointMutationsCountArray, axis=1)

    #Fill for indels
    if (chrBased_indels_df is not None):
        chrBased_indels_df.apply(fillSplitBasedSignalArrayAndCountArrayForIndels, nucleosome_array=chrBased_nucleosome_split_array, allIndelsSignalArray=allIndelsSignalArray,allIndelsCountArray=allIndelsCountArray, axis=1)

    #Append signature arrays
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(signature2SignalArrayDict)
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(signature2CountArrayDict)

    #Append single point mutations arrays
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(allSinglePointMutationsSignalArray)
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(allSinglePointMutationsCountArray)

    #Append indels arrays
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(allIndelsSignalArray)
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(allIndelsCountArray)

    return split_chrBased_SignalArrayAndCountArray_DictionaryList
########################################################################################


########################################################################################
def accumulateSplitArraysWithExtraSampleBased(signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,allSplits_chrBased_SignalArrayAndCountArray_DictionaryList):

    #Initialize them
    #Fill them using allSplits_chrBased_SignalArrayAndCountArray_DictionaryList
    #Return them

    ########################################################################
    ############ Initialization3 for all splits starts #####################
    ########################################################################
    signature2AccumulatedSplitsChrBasedSignalArrayDict = {}
    signature2AccumulatedSplitsChrBasedCountArrayDict = {}

    allSingleMutationsAccumulatedSplitsSignalArray = np.zeros(windowSize)
    allSingleMutationsAccumulatedSplitsCountArray = np.zeros(windowSize, dtype=int)

    allIndelsAccumulatedSplitsSignalArray = np.zeros(windowSize)
    allIndelsAccumulatedSplitsCountArray = np.zeros(windowSize, dtype=int)

    for signature in signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
        signature2AccumulatedSplitsChrBasedSignalArrayDict[signature] = np.zeros(windowSize)
        signature2AccumulatedSplitsChrBasedCountArrayDict[signature] = np.zeros(windowSize, dtype=int)

    sample2Signature2AccumulatedSplitsChrBasedSignalArrayDict = {}
    sample2Signature2AccumulatedSplitsChrBasedCountArrayDict = {}

    sample2AllSinglePointMutationsAccumulatedSplitsChrBasedSignalArrayDict = {}
    sample2AllSinglePointMutationsAccumulatedSplitsChrBasedCountArrayDict = {}

    sample2AllIndelsAccumulatedSplitsChrBasedSignalArrayDict = {}
    sample2AllIndelsAccumulatedSplitsChrBasedCountArrayDict = {}


    for sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
        sample2Signature2AccumulatedSplitsChrBasedSignalArrayDict[sample] ={}
        sample2Signature2AccumulatedSplitsChrBasedCountArrayDict[sample] ={}

        sample2AllSinglePointMutationsAccumulatedSplitsChrBasedSignalArrayDict[sample] = np.zeros(windowSize)
        sample2AllSinglePointMutationsAccumulatedSplitsChrBasedCountArrayDict[sample] = np.zeros(windowSize, dtype=int)

        sample2AllIndelsAccumulatedSplitsChrBasedSignalArrayDict[sample] = np.zeros(windowSize)
        sample2AllIndelsAccumulatedSplitsChrBasedCountArrayDict[sample] = np.zeros(windowSize, dtype=int)

        for signature in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict[sample]:
            sample2Signature2AccumulatedSplitsChrBasedSignalArrayDict[sample][signature] = np.zeros(windowSize)
            sample2Signature2AccumulatedSplitsChrBasedCountArrayDict[sample][signature] = np.zeros(windowSize, dtype=int)
    ########################################################################
    ############ Initialization3 for all splits ends #######################
    ########################################################################

    ##############################################
    for split_chrBased_SignalArrayAndCountArray_DictionaryList in allSplits_chrBased_SignalArrayAndCountArray_DictionaryList:

        signature2SplitSignalArrayDict = split_chrBased_SignalArrayAndCountArray_DictionaryList[0]
        signature2SplitCountArrayDict = split_chrBased_SignalArrayAndCountArray_DictionaryList[1]

        allSingePointMutationsSplitSignalArray = split_chrBased_SignalArrayAndCountArray_DictionaryList[2]
        allSingePointMutationsSplitCountArray = split_chrBased_SignalArrayAndCountArray_DictionaryList[3]

        allIndelsSplitSignalArray = split_chrBased_SignalArrayAndCountArray_DictionaryList[4]
        allIndelsSplitCountArray = split_chrBased_SignalArrayAndCountArray_DictionaryList[5]

        sample2Signature2SplitSignalArrayDict = split_chrBased_SignalArrayAndCountArray_DictionaryList[6]
        sample2Signature2SplitCountArrayDict = split_chrBased_SignalArrayAndCountArray_DictionaryList[7]

        sample2AllSinglePointMutationsSplitSignalArrayDict = split_chrBased_SignalArrayAndCountArray_DictionaryList[8]
        sample2AllSinglePointMutationsSplitCountArrayDict = split_chrBased_SignalArrayAndCountArray_DictionaryList[9]

        sample2AllIndelsSplitSignalArrayDict = split_chrBased_SignalArrayAndCountArray_DictionaryList[10]
        sample2AllIndelsSplitCountArrayDict = split_chrBased_SignalArrayAndCountArray_DictionaryList[11]

        #################################################################

        #####################################################################################################
        ######################### Accumulate right in the left starts  ######################################
        #####################################################################################################

        ######################### Accumulate starts ######################################
        accumulateSignatureBasedArrays(signature2AccumulatedSplitsChrBasedSignalArrayDict,signature2SplitSignalArrayDict)
        accumulateSignatureBasedArrays(signature2AccumulatedSplitsChrBasedCountArrayDict,signature2SplitCountArrayDict)

        allSingleMutationsAccumulatedSplitsSignalArray += allSingePointMutationsSplitSignalArray
        allSingleMutationsAccumulatedSplitsCountArray += allSingePointMutationsSplitCountArray

        allIndelsAccumulatedSplitsSignalArray += allIndelsSplitSignalArray
        allIndelsAccumulatedSplitsCountArray += allIndelsSplitCountArray
        ######################### Accumulate starts ######################################

        ############################Accumulate Sample Based starts ###################################
        accumulateSampleBasedSignatureBasedArrays(sample2Signature2AccumulatedSplitsChrBasedSignalArrayDict,sample2Signature2SplitSignalArrayDict)
        accumulateSampleBasedSignatureBasedArrays(sample2Signature2AccumulatedSplitsChrBasedCountArrayDict,sample2Signature2SplitCountArrayDict)

        accumulateSampleBasedArrays(sample2AllSinglePointMutationsAccumulatedSplitsChrBasedSignalArrayDict,sample2AllSinglePointMutationsSplitSignalArrayDict)
        accumulateSampleBasedArrays(sample2AllSinglePointMutationsAccumulatedSplitsChrBasedCountArrayDict,sample2AllSinglePointMutationsSplitCountArrayDict)

        accumulateSampleBasedArrays(sample2AllIndelsAccumulatedSplitsChrBasedSignalArrayDict,sample2AllIndelsSplitSignalArrayDict)
        accumulateSampleBasedArrays(sample2AllIndelsAccumulatedSplitsChrBasedCountArrayDict,sample2AllIndelsSplitCountArrayDict)
        ############################Accumulate Sample Based ends #####################################

        #####################################################################################################
        ######################### Accumulate right in the left ends  ########################################
        #####################################################################################################


    ##############################################

    return  signature2AccumulatedSplitsChrBasedSignalArrayDict, \
            signature2AccumulatedSplitsChrBasedCountArrayDict, \
            allSingleMutationsAccumulatedSplitsSignalArray, \
            allSingleMutationsAccumulatedSplitsCountArray, \
            allIndelsAccumulatedSplitsSignalArray,\
            allIndelsAccumulatedSplitsCountArray,\
            sample2Signature2AccumulatedSplitsChrBasedSignalArrayDict,\
            sample2Signature2AccumulatedSplitsChrBasedCountArrayDict,\
            sample2AllSinglePointMutationsAccumulatedSplitsChrBasedSignalArrayDict,\
            sample2AllSinglePointMutationsAccumulatedSplitsChrBasedCountArrayDict,\
            sample2AllIndelsAccumulatedSplitsChrBasedSignalArrayDict,\
            sample2AllIndelsAccumulatedSplitsChrBasedCountArrayDict
########################################################################################


########################################################################################
#This function is specific to this python file
def accumulateSplitArrays(signatures,allSplits_chrBased_SignalArrayAndCountArray_DictionaryList):

    ##############################################
    # Fill these signature dictionaries
    signature2AccumulatedSplitsChrBasedSignalArrayDict = {}
    signature2AccumulatedSplitsChrBasedCountArrayDict = {}

    for signature in signatures:
        signature2AccumulatedSplitsChrBasedSignalArrayDict[signature] = np.zeros(windowSize)
        signature2AccumulatedSplitsChrBasedCountArrayDict[signature] = np.zeros(windowSize,dtype=int)
    ##############################################

    ##############################################
    # Fill these single point mutations arrays
    allSingleMutationsAccumulatedSplitsSignalArray = np.zeros(windowSize)
    allSingleMutationsAccumulatedSplitsCountArray = np.zeros(windowSize, dtype=int)
    ##############################################

    ##############################################
    # Fill these indels arrays
    allIndelsAccumulatedSplitsSignalArray = np.zeros(windowSize)
    allIndelsAccumulatedSplitsCountArray = np.zeros(windowSize, dtype=int)
    ##############################################

    ##############################################
    for split_chrBased_SignalArrayAndCountArray_DictionaryList in allSplits_chrBased_SignalArrayAndCountArray_DictionaryList:

        ############################################################################################
        signature2SplitSignalArrayDict = split_chrBased_SignalArrayAndCountArray_DictionaryList[0]
        signature2SplitCountArrayDict = split_chrBased_SignalArrayAndCountArray_DictionaryList[1]

        allSingePointMutationsSplitSignalArray = split_chrBased_SignalArrayAndCountArray_DictionaryList[2]
        allSingePointMutationsSplitCountArray = split_chrBased_SignalArrayAndCountArray_DictionaryList[3]

        allIndelsSplitSignalArray = split_chrBased_SignalArrayAndCountArray_DictionaryList[4]
        allIndelsSplitCountArray = split_chrBased_SignalArrayAndCountArray_DictionaryList[5]
        ############################################################################################

        ############################################################################################
        #Accumulate Signatures
        accumulateSignatureBasedArrays(signature2AccumulatedSplitsChrBasedSignalArrayDict,signature2SplitSignalArrayDict)
        accumulateSignatureBasedArrays(signature2AccumulatedSplitsChrBasedCountArrayDict,signature2SplitCountArrayDict)

        #Accumulate Single Point Mutations
        allSingleMutationsAccumulatedSplitsSignalArray += allSingePointMutationsSplitSignalArray
        allSingleMutationsAccumulatedSplitsCountArray += allSingePointMutationsSplitCountArray

        #Accumulate Indels
        allIndelsAccumulatedSplitsSignalArray += allIndelsSplitSignalArray
        allIndelsAccumulatedSplitsCountArray += allIndelsSplitCountArray
        ############################################################################################

    return signature2AccumulatedSplitsChrBasedSignalArrayDict, signature2AccumulatedSplitsChrBasedCountArrayDict, allSingleMutationsAccumulatedSplitsSignalArray, allSingleMutationsAccumulatedSplitsCountArray, allIndelsAccumulatedSplitsSignalArray, allIndelsAccumulatedSplitsCountArray
########################################################################################


########################################################################################
def nucleosomeOccupancyAnalysis_SPMs_SignatureBased_Indels(jobname,singlePointMutationsFilename,indelsFilename,nucleosomeFilename):

    withExtraSampleBasedAnalysis = True

    ################# input starts #############################
    # print(sys.argv)

    # sys.argv[0] contains the python program name
    # jobname
    # jobname = sys.argv[1]

    # breast_cancer_mutation_probabilities_final.txt
    # singlePointMutationsFilename = sys.argv[2]
    # singlePointMutationsFilename_wo_extension = os.path.splitext(singlePointMutationsFilename)[0]

    # insertions and deletions file
    # indelsFilename = sys.argv[3]

    # wgEncodeSydhNsomeK562Sig.wig
    # nucleosomeFilename = sys.argv[4]
    #################### input ends ############################

    ##########################################################################################
    #Load Chrnames from nucleosome file
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
    ########### Initialization1 All Chroms starts ##########################################
    ########################################################################################

    # Initialization for each signature
    signature2AccumulatedAllChromsSignalArrayDict = {}
    signature2AccumulatedAllChromsCountArrayDict = {}

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
    # Initialization for indels
    allIndelsAccumulatedAllChromsSignalArray = np.zeros(windowSize)
    allIndelsAccumulatedAllChromsCountArray = np.zeros(windowSize, dtype=int)
    ########################################################################################

    ########################################################################################
    # Sample Based
    sample2Signature2AccumulatedAllChromsSignalArrayDict = {}
    sample2Signature2AccumulatedAllChromsCountArrayDict = {}

    sample2AllSinglePointMutationsAccumulatedAllChromsSignalArrayDict = {}
    sample2AllSinglePointMutationsAccumulatedAllChromsCountArrayDict = {}

    sample2AllIndelsAccumulatedAllChromsSignalArrayDict = {}
    sample2AllIndelsAccumulatedAllChromsCountArrayDict = {}

    for sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
        sample2Signature2AccumulatedAllChromsSignalArrayDict[sample] = {}
        sample2Signature2AccumulatedAllChromsCountArrayDict[sample]= {}

        sample2AllSinglePointMutationsAccumulatedAllChromsSignalArrayDict[sample] = np.zeros(windowSize)
        sample2AllSinglePointMutationsAccumulatedAllChromsCountArrayDict[sample] = np.zeros(windowSize, dtype=int)

        sample2AllIndelsAccumulatedAllChromsSignalArrayDict[sample] = np.zeros(windowSize)
        sample2AllIndelsAccumulatedAllChromsCountArrayDict[sample] = np.zeros(windowSize, dtype=int)

        for signature in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict[sample]:
            sample2Signature2AccumulatedAllChromsSignalArrayDict[sample][signature] =  np.zeros(windowSize)
            sample2Signature2AccumulatedAllChromsCountArrayDict[sample][signature] = np.zeros(windowSize, dtype=int)
    ########################################################################################

    ########################################################################################
    ########### Initialization1 All Chroms ends ############################################
    ########################################################################################


    ######################################################
    print('chrNamesInNucleosomeList')
    print(chrNamesInNucleosomeList)
    # e.g.: chrNamesList = ['1','10','11','12','13', '14', '15', '16', '17', '18', '19', '2', '20', '21', '22', '3', '4', '5', '6', '7', '8', '9', 'MT', 'X', 'Y']
    ######################################################

    ######################################################
    #Better practice
    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)
    print('multiprocessing.cpu_count() --> number of processes: %d' % numofProcesses)
    ######################################################

    ###################################################################################
    ##################  For each chromsome sequential starts ##########################
    ###################################################################################
    for chrLong in chrNamesInNucleosomeList:

        #FIRST READ CHRBASED NUCLEOSOME OCCUPANCY
        chrBased_nucleosome_df = readChrBasedNuclesomeDF(chrLong,nucleosomeFilename)
        print('chromosome %s  -- chrBased_nucleosome_df: %d in bytes %f in GB' % (chrLong,sys.getsizeof(chrBased_nucleosome_df), sys.getsizeof(chrBased_nucleosome_df)/GIGABYTE_IN_BYTES))

        #TODO: This is specific to our data right now
        #Nucleosomes have chrM
        #SinglePointMutations and Indels have chrMT
        if (chrLong=='chrM'):
            chrLong='chrMT'

        #THEN READ CHRBASED SINGLE POINT MUTATIONS
        chrBased_spms_df=readChrBasedMutationDF(jobname,chrLong,singlePointMutationsFilename)
        print('chromosome %s  -- chrBased_spms_df: %d in Bytes %f in GigaBytes' %(chrLong,sys.getsizeof(chrBased_spms_df),sys.getsizeof(chrBased_spms_df)/GIGABYTE_IN_BYTES))

        #THEN READ CHRBASED INDELS
        chrBased_indels_df = readChrBasedIndelsDF(jobname,chrLong,indelsFilename)

        if (chrBased_spms_df is not None):

            #################################################
            #Question: Shall I split chrBased_spms_df or chrBased_nucleosome_df?
            #Answer: Since chrBased_nucleosome_df is bigger than chrBased_spms_df, therefore split chrBased_nucleosome_df

            # Question: Shall I do fill nucleosome array here or in each core separately?
            #for all rows of nucleosome_df update the np.array of zero using slicing
            #Answer: Let's do it in each core separately
            # chrBased_nucleosome_df.apply(fillNucleosomeArray, nucleosome_array=chrBased_nucleosome_array, axis=1)
            #################################################

            #################################################
            # split chrBased_spms_df into number of available cores
            # for each split_chrBased_nucleosome_df first fill nucleosome_array
            # then fill split_signal and split_count array
            # accumulate all the split_signal and split_count arrays in chrBased signal array and count arrays

            #TODO: later on split the one that is bigger whether it is chrBasedMutationDF or chrBasedNucleosomeDF
            # chrBased_mutation_df_splits = np.array_split(chrBased_skin_df, numofProcesses)


            ###############################################################################
            if (chrBased_nucleosome_df is not None):

                chrBased_nucleosome_df_splits = np.array_split(chrBased_nucleosome_df,numofProcesses)

                #Prepare the poolInputList
                poolInputList = []

                #########################################################################
                if (withExtraSampleBasedAnalysis):

                    ######################################################################################
                    ##############   SPMs and Indels with extra sampleBased  analysis starts #############
                    ######################################################################################


                    # For each split chrBased mutation df do it in parallel
                    for chrBased_nucleosome_df_split in chrBased_nucleosome_df_splits:
                        inputList = []

                        # Question: Does each split has header?
                        # Answer: Yes
                        inputList.append(chrBased_spms_df) #Same chrBased_mutation_df
                        inputList.append(chrBased_indels_df)  # Same chrBased_indels_df
                        inputList.append(chrBased_nucleosome_df_split)  #Different chrBased_nucleosome_df split
                        inputList.append(chrLong) #same chrLong
                        inputList.append(signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict) #same signature
                        inputList.append(samplesWithAtLeast10KMutations2NumberofMutationsDict) #Sep 18, 2018
                        inputList.append(sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict) #Sep 19, 2018
                        poolInputList.append(inputList)

                    #########################################################################
                    allSplits_chrBased_SignalArrayAndCountArray_DictionaryList = pool.map(fillSignalArrayAndCountArrayWithExtraSampleBased,poolInputList)


                    #Accumulate the results coming from each split
                    #allSplits_chrBased_SignalArrayAndCountArray_DictionaryList is a list of lists.
                    #Each list consist of two dictionaries
                    signature2AccumulatedSplitsChrBasedSignalArrayDict, \
                    signature2AccumulatedSplitsChrBasedCountArrayDict, \
                    allSingleMutationsAccumulatedSplitsSignalArray, \
                    allSingleMutationsAccumulatedSplitsCountArray, \
                    allIndelsAccumulatedSplitsSignalArray, \
                    allIndelsAccumulatedSplitsCountArray, \
                    sample2Signature2AccumulatedSplitsChrBasedSignalArrayDict, \
                    sample2Signature2AccumulatedSplitsChrBasedCountArrayDict, \
                    sample2AllSinglePointMutationsAccumulatedSplitsChrBasedSignalArrayDict, \
                    sample2AllSinglePointMutationsAccumulatedSplitsChrBasedCountArrayDict, \
                    sample2AllIndelsAccumulatedSplitsSignalArrayDict, \
                    sample2AllIndelsAccumulatedSplitsCountArrayDict = accumulateSplitArraysWithExtraSampleBased(
                        signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
                        sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
                        allSplits_chrBased_SignalArrayAndCountArray_DictionaryList)

                    #################################################


                    ###############################################################################################################
                    ###################  Accumulate the results coming from each chromosome starts  ###############################
                    ###############################################################################################################

                    #Accumulate right in the left
                    #################################################
                    accumulateSignatureBasedArrays(signature2AccumulatedAllChromsSignalArrayDict,signature2AccumulatedSplitsChrBasedSignalArrayDict)
                    accumulateSignatureBasedArrays(signature2AccumulatedAllChromsCountArrayDict,signature2AccumulatedSplitsChrBasedCountArrayDict)

                    ##Accumulate the results coming from each chromosome
                    allSinglePointMutationsAccumulatedAllChromsSignalArray += allSingleMutationsAccumulatedSplitsSignalArray
                    allSinglePointMutationsAccumulatedAllChromsCountArray += allSingleMutationsAccumulatedSplitsCountArray

                    allIndelsAccumulatedAllChromsSignalArray += allIndelsAccumulatedSplitsSignalArray
                    allIndelsAccumulatedAllChromsCountArray += allIndelsAccumulatedSplitsCountArray
                    #################################################

                    ###########################Accumulate sample based starts ####################################
                    ##Accumulate the results coming from each chromosome
                    accumulateSampleBasedSignatureBasedArrays(sample2Signature2AccumulatedAllChromsSignalArrayDict,sample2Signature2AccumulatedSplitsChrBasedSignalArrayDict)
                    accumulateSampleBasedSignatureBasedArrays(sample2Signature2AccumulatedAllChromsCountArrayDict,sample2Signature2AccumulatedSplitsChrBasedCountArrayDict)

                    accumulateSampleBasedArrays(sample2AllSinglePointMutationsAccumulatedAllChromsSignalArrayDict,sample2AllSinglePointMutationsAccumulatedSplitsChrBasedSignalArrayDict)
                    accumulateSampleBasedArrays(sample2AllSinglePointMutationsAccumulatedAllChromsCountArrayDict,sample2AllSinglePointMutationsAccumulatedSplitsChrBasedCountArrayDict)

                    accumulateSampleBasedArrays(sample2AllIndelsAccumulatedAllChromsSignalArrayDict,sample2AllIndelsAccumulatedSplitsSignalArrayDict)
                    accumulateSampleBasedArrays(sample2AllIndelsAccumulatedAllChromsCountArrayDict,sample2AllIndelsAccumulatedSplitsCountArrayDict)
                    ###########################Accumulate sample based ends ######################################

                    ###############################################################################################################
                    ###################  Accumulate the results coming from each chromosome ends  #################################
                    ###############################################################################################################



                    ######################################################################################
                    ##############   SPMs and Indels with extra sampleBased  analysis ends ###############
                    ######################################################################################


                #########################################################################

                #########################################################################
                else:

                    ###########################################################################
                    ##############   SPMs and INDELs analysis starts ##########################
                    ###########################################################################

                    #########################################################################
                    # For each split chrBased mutation df do it in parallel
                    for chrBased_nucleosome_df_split in chrBased_nucleosome_df_splits:
                        inputList = []

                        # Question: Does each split has header?
                        # Answer: Yes
                        inputList.append(chrBased_spms_df)  # Same chrBased_spms_df
                        inputList.append(chrBased_indels_df)  # Same chrBased_indels_df
                        inputList.append(chrBased_nucleosome_df_split)  # Different chrBased_nucleosome_df split
                        inputList.append(chrLong)  # same chrLong
                        inputList.append(signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict)  # same signature
                        poolInputList.append(inputList)
                    #########################################################################

                    #########################################################################
                    allSplits_chrBased_SignalArrayAndCountArray_DictionaryList = pool.map(fillNucleosomeSignalArrayAndCountArray, poolInputList)

                    # Accumulate the results coming from each split
                    # allSplits_chrBased_SignalArrayAndCountArray_DictionaryList is a list of lists.
                    # Each list consist of two dictionaries
                    signature2AccumulatedSplitsChrBasedSignalArrayDict, \
                    signature2AccumulatedSplitsChrBasedCountArrayDict, \
                    allSingleMutationsAccumulatedSplitsSignalArray, \
                    allSingleMutationsAccumulatedSplitsCountArray, \
                    allIndelsAccumulatedSplitsSignalArray, \
                    allIndelsAccumulatedSplitsCountArray = accumulateSplitArrays(signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict, allSplits_chrBased_SignalArrayAndCountArray_DictionaryList)
                    #################################################

                    #################################################
                    ##Accumulate the signature based results coming from each chromosome
                    accumulateSignatureBasedArrays(signature2AccumulatedAllChromsSignalArrayDict,signature2AccumulatedSplitsChrBasedSignalArrayDict)
                    accumulateSignatureBasedArrays(signature2AccumulatedAllChromsCountArrayDict,signature2AccumulatedSplitsChrBasedCountArrayDict)

                    ##Accumulate the single point mutations results coming from each chromosome
                    allSinglePointMutationsAccumulatedAllChromsSignalArray += allSingleMutationsAccumulatedSplitsSignalArray
                    allSinglePointMutationsAccumulatedAllChromsCountArray += allSingleMutationsAccumulatedSplitsCountArray

                    ##Accumulate the indels results coming from each chromosome
                    allIndelsAccumulatedAllChromsSignalArray += allIndelsAccumulatedSplitsSignalArray
                    allIndelsAccumulatedAllChromsCountArray += allIndelsAccumulatedSplitsCountArray
                    #################################################

                    ###########################################################################
                    ##############   SPMs and INDELs analysis ends ############################
                    ###########################################################################

                #########################################################################

            ##############################################################################
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
    writeAverageNucleosomeOccupancyFiles(allSinglePointMutationsAccumulatedAllChromsSignalArray,allSinglePointMutationsAccumulatedAllChromsCountArray, jobname, AGGREGATEDSUBSTITUTIONS)
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
    ################ Write All Indels Average Nucleosome Occupancy starts ##############
    ####################################################################################
    writeAverageNucleosomeOccupancyFiles(allIndelsAccumulatedAllChromsSignalArray,allIndelsAccumulatedAllChromsCountArray, jobname, AGGREGATEDINDELS)
    ####################################################################################
    ################ Write All Indels Average Nucleosome Occupancy ends ################
    ####################################################################################


    if (withExtraSampleBasedAnalysis):

        ####################################################################################
        writeSampleBasedSignatureBasedAverageNucleosomeOccupancyFiles(sample2Signature2AccumulatedAllChromsSignalArrayDict,sample2Signature2AccumulatedAllChromsCountArrayDict,jobname)
        ####################################################################################

        ####################################################################################
        writeSampleBasedAverageNucleosomeOccupancyFiles(sample2AllSinglePointMutationsAccumulatedAllChromsSignalArrayDict,sample2AllSinglePointMutationsAccumulatedAllChromsCountArrayDict,jobname,AGGREGATEDSUBSTITUTIONS)
        ####################################################################################

        ####################################################################################
        writeSampleBasedAverageNucleosomeOccupancyFiles(sample2AllIndelsAccumulatedAllChromsSignalArrayDict,sample2AllIndelsAccumulatedAllChromsCountArrayDict,jobname,AGGREGATEDINDELS)
        ####################################################################################
########################################################################################
