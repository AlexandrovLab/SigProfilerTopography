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
#   for all signatures with all single point mutations with a certain probability for that signature
#   for all indels
#   TODO for all signatures with all indels with a certain probability
#   TODO for all dinucs
###############################################################################################################

import os
import sys

#############################################################
current_abs_path = os.path.abspath(os.path.dirname(__file__))
print('NucleosomeOccupancyAnalysis_Subs_Indels.py current_abs_path:%s' %(current_abs_path))
#############################################################

commonsPath = os.path.join(current_abs_path, '..','commons')
sys.path.append(commonsPath)

from SigProfilerTopography.source.commons.TopographyCommons import *


########################################################################################
#March 12, 2019 starts
def fillSignalArrayAndCountArraysIndelsSignatureBasedAdded(inputList):
    pass
#Tomorrow fill signature also for indels
# subsSignature
# indelsSignature
# when done we can use this function
# left here
#March 12, 2019 ends
########################################################################################


########################################################################################
#March 7, 2019 starts
def fillSignalArrayAndCountArrays(inputList):

    chrbased_nucleosome_signal_array = inputList[0]
    chrBased_spms_df = inputList[1]
    chrBased_indels_df =  inputList[2]
    maximum_chrom_size =  inputList[3]
    signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict =  inputList[4]
    samplesWithAtLeast10KMutations2NumberofMutationsDict =  inputList[5]
    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict =  inputList[6]

    ##############################################################
    signature2SignalArrayDict, \
    signature2CountArrayDict, \
    allSinglePointMutationsSignalArray, \
    allSinglePointMutationsCountArray, \
    allIndelsSignalArray, \
    allIndelsCountArray, \
    sample2Signature2SignalArrayDict, \
    sample2Signature2CountArrayDict, \
    sample2AllSinglePointMutationsSignalArrayDict, \
    sample2AllSinglePointMutationsCountArrayDict, \
    sample2AllIndelsSignalArrayDict, \
    sample2AllIndelsCountArrayDict = initializationForEachChrom(
        signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
        sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict)
    ##############################################################


    ###############################################################################
    ################### Fill signal and count array starts ########################
    ###############################################################################
    #Fill for single point mutations
    if ((chrBased_spms_df is not None) and (not chrBased_spms_df.empty)):
        chrBased_spms_df.apply(fillSplitBasedSignalArrayAndCountArrayForSPMsWithExtraSampleBased,
                           nucleosome_array=chrbased_nucleosome_signal_array,
                           maximum_chrom_size =maximum_chrom_size,
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
    if ((chrBased_indels_df is not None) and (not chrBased_indels_df.empty)):
        chrBased_indels_df.apply(fillSplitBasedSignalArrayAndCountArrayForIndelsWithExtraSampleBased,
                                 nucleosome_array=chrbased_nucleosome_signal_array,
                                 maximum_chrom_size=maximum_chrom_size,
                                 sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict=sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
                                 allIndelsSignalArray=allIndelsSignalArray,
                                 allIndelsCountArray=allIndelsCountArray,
                                 sample2AllIndelsSignalArrayDict=sample2AllIndelsSignalArrayDict,
                                 sample2AllIndelsCountArrayDict =sample2AllIndelsCountArrayDict,
                                 axis=1)
    ###############################################################################
    ################### Fill signal and count array ends ##########################
    ###############################################################################


    ###############################################################################
    ################### Return  starts ############################################
    ###############################################################################
    # Initialzie the list, you will return this list
    chrBased_SignalArrayAndCountArray_List = []

    # Append signature arrays
    chrBased_SignalArrayAndCountArray_List.append(signature2SignalArrayDict)
    chrBased_SignalArrayAndCountArray_List.append(signature2CountArrayDict)

    # Append single point mutations arrays
    chrBased_SignalArrayAndCountArray_List.append(allSinglePointMutationsSignalArray)
    chrBased_SignalArrayAndCountArray_List.append(allSinglePointMutationsCountArray)

    # Append indels arrays
    chrBased_SignalArrayAndCountArray_List.append(allIndelsSignalArray)
    chrBased_SignalArrayAndCountArray_List.append(allIndelsCountArray)

    chrBased_SignalArrayAndCountArray_List.append(sample2Signature2SignalArrayDict)
    chrBased_SignalArrayAndCountArray_List.append(sample2Signature2CountArrayDict)

    chrBased_SignalArrayAndCountArray_List.append(sample2AllSinglePointMutationsSignalArrayDict)
    chrBased_SignalArrayAndCountArray_List.append(sample2AllSinglePointMutationsCountArrayDict)

    chrBased_SignalArrayAndCountArray_List.append(sample2AllIndelsSignalArrayDict)
    chrBased_SignalArrayAndCountArray_List.append(sample2AllIndelsCountArrayDict)

    return chrBased_SignalArrayAndCountArray_List
    ###############################################################################
    ################### Return  ends ##############################################
    ###############################################################################

# March 7, 2019 starts
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
#March 7, 2019 starts
def fillNecessaryDataStructures(genome,outputDir,jobname):

    ##########################################################################################
    chromSizesDict = getChromSizesDict(genome)
    ##########################################################################################

    ##########################################################################################
    #Load Chrnames from nucleosome file
    chrNamesInNucleosomeList = []

    chrNamesInNucleosomeFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME, ChrNamesInNucleosomesFilename)

    if (os.path.exists(chrNamesInNucleosomeFilePath)):
        chrNamesInNucleosomeArray = np.loadtxt(chrNamesInNucleosomeFilePath,dtype=str, delimiter='\t')
        chrNamesInNucleosomeList = list(chrNamesInNucleosomeArray)
    ##########################################################################################

    ##########################################################################################
    #Load signaturesWithAtLeast10KEligibleMutations
    signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = {}

    SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath = os.path.join(outputDir,jobname,DATA,SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilename)

    if (os.path.exists(SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath)):
        signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = readDictionary(SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath)
    ##########################################################################################

    ##########################################################################################
    #Load sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict
    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict= {}

    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath = os.path.join(outputDir,jobname,DATA,Sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilename)

    if (os.path.exists(sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath)):
        sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = readDictionary(sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath)
    ##########################################################################################

    ##########################################################################################
    #Load samplesWithAtLeast10KMutations2NumberofMutationsDict
    samplesWithAtLeast10KMutations2NumberofMutationsDict = {}

    samplesWithAtLeast10KMutations2NumberofMutationsDictFilePath = os.path.join(outputDir,jobname,DATA,SamplesWithAtLeast10KMutations2NumberofMutationsDictFilename)

    if (os.path.exists(samplesWithAtLeast10KMutations2NumberofMutationsDictFilePath)):
        samplesWithAtLeast10KMutations2NumberofMutationsDict = readDictionary(samplesWithAtLeast10KMutations2NumberofMutationsDictFilePath)
    ##########################################################################################

    return  chromSizesDict, \
            chrNamesInNucleosomeList,\
            signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,\
            sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,\
            samplesWithAtLeast10KMutations2NumberofMutationsDict

#March 7, 2019 ends
########################################################################################

########################################################################################
#March 7, 2019 starts
def initializationForAllChroms(signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict):

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

    return  signature2AccumulatedAllChromsSignalArrayDict,\
            signature2AccumulatedAllChromsCountArrayDict, \
            allSinglePointMutationsAccumulatedAllChromsSignalArray,\
            allSinglePointMutationsAccumulatedAllChromsCountArray,\
            allIndelsAccumulatedAllChromsSignalArray,\
            allIndelsAccumulatedAllChromsCountArray,\
            sample2Signature2AccumulatedAllChromsSignalArrayDict,\
            sample2Signature2AccumulatedAllChromsCountArrayDict,\
            sample2AllSinglePointMutationsAccumulatedAllChromsSignalArrayDict,\
            sample2AllSinglePointMutationsAccumulatedAllChromsCountArrayDict,\
            sample2AllIndelsAccumulatedAllChromsSignalArrayDict,\
            sample2AllIndelsAccumulatedAllChromsCountArrayDict

#March 7, 2019 ends
########################################################################################



########################################################################################
#March 7, 2019 starts
def initializationForEachChrom(signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict):
    ########################################################################################
    ########### Initialization2 For each chrom starts ######################################
    ########################################################################################

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
    allSinglePointMutationsSignalArray = np.zeros(windowSize)
    allSinglePointMutationsCountArray = np.zeros(windowSize, dtype=int)
    ##############################################

    ##############################################
    # Initialization for all indels for split
    allIndelsSignalArray = np.zeros(windowSize)
    allIndelsCountArray = np.zeros(windowSize, dtype=int)
    ##############################################

    # Initialization sample based
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
    ########################################################################################
    ########### Initialization2 For each chrom ends ########################################
    ########################################################################################

    return  signature2SignalArrayDict,\
            signature2CountArrayDict, \
            allSinglePointMutationsSignalArray,\
            allSinglePointMutationsCountArray, \
            allIndelsSignalArray, \
            allIndelsCountArray, \
            sample2Signature2SignalArrayDict, \
            sample2Signature2CountArrayDict, \
            sample2AllSinglePointMutationsSignalArrayDict, \
            sample2AllSinglePointMutationsCountArrayDict, \
            sample2AllIndelsSignalArrayDict, \
            sample2AllIndelsCountArrayDict


#March 7, 2019 ends
########################################################################################


########################################################################################
#March 7, 2019 starts
#For all chromosome parallel starts
def nucleosome_occupancy_analysis_subs_indels_all_chroms_parallel(genome,outputDir, jobname,singlePointMutationsFilename,indelsFilename,nucleosomeFilename):

    ##########################################################################
    chromSizesDict, \
    chrNamesInNucleosomeList, \
    signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict, \
    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict, \
    samplesWithAtLeast10KMutations2NumberofMutationsDict = fillNecessaryDataStructures(genome,outputDir,jobname)
    ##########################################################################

    ##########################################################################
    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)
    poolInputList = []
    ##########################################################################

    ###################################################################################
    ##################  For all chromsomes parallel starts ############################
    ###################################################################################
    for chrLong in chrNamesInNucleosomeList:

        inputListForSubs = []
        inputListForIndels = []

        maximum_chrom_size = chromSizesDict[chrLong]
        print('For debug %s chromSize:%d' %(chrLong,maximum_chrom_size))

        #FIRST READ CHRBASED NUCLEOSOME OCCUPANCY
        nucleosomeFilenameWoExtension = os.path.basename(nucleosomeFilename)[0:-4]

        ##############################################################
        signalArrayFilename = '%s_signal_%s.npy' % (chrLong, nucleosomeFilenameWoExtension)
        chrBasedSignalNucleosmeFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB,NUCLEOSOME, CHRBASED, signalArrayFilename)

        #################################################################################################################
        # if (os.path.exists(chrBasedSignalNucleosmeFile) and os.path.exists(chrBasedCountNucleosmeFile)):
        if (os.path.exists(chrBasedSignalNucleosmeFile)):
            chrbased_nucleosome_signal_array = np.load(chrBasedSignalNucleosmeFile)
            print('chromosome %s  -- signal_array_npy: %d in bytes %f in GB' % (chrLong,sys.getsizeof(chrbased_nucleosome_signal_array), sys.getsizeof(chrbased_nucleosome_signal_array)/GIGABYTE_IN_BYTES))

            #TODO: This is specific to our data right now
            #Nucleosomes have chrM
            #SinglePointMutations and Indels have chrMT
            if (chrLong=='chrM'):
                chrLong='chrMT'

            #THEN READ CHRBASED SINGLE POINT MUTATIONS
            if  (singlePointMutationsFilename!= NOTSET):
                chrBased_spms_df=readChrBasedMutationDF(outputDir,jobname,chrLong,singlePointMutationsFilename)
                print('chromosome %s  -- chrBased_spms_df: %d in Bytes %f in GigaBytes' %(chrLong,sys.getsizeof(chrBased_spms_df),sys.getsizeof(chrBased_spms_df)/GIGABYTE_IN_BYTES))
            else:
                chrBased_spms_df = None


            #THEN READ CHRBASED INDELS
            if (indelsFilename!=NOTSET):
                chrBased_indels_df = readChrBasedIndelsDF(outputDir,jobname,chrLong,indelsFilename)
                print('chromosome %s  -- chrBased_indels_df: %d in Bytes %f in GigaBytes' %(chrLong,sys.getsizeof(chrBased_indels_df),sys.getsizeof(chrBased_indels_df)/GIGABYTE_IN_BYTES))
            else:
                chrBased_indels_df= None

            inputListForSubs.append(chrbased_nucleosome_signal_array)
            inputListForSubs.append(chrBased_spms_df)
            inputListForSubs.append(chrBased_indels_df)
            inputListForSubs.append(maximum_chrom_size)
            inputListForSubs.append(signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict)
            inputListForSubs.append(samplesWithAtLeast10KMutations2NumberofMutationsDict)
            inputListForSubs.append(sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict)
            poolInputList.append(inputListForSubs)
    ###################################################################################
    ##################  For all chromsomes parallel ends ##############################
    ###################################################################################

    ########################################################################
    allChroms_SignalArrayAndCountArrayList_List = pool.map(fillSignalArrayAndCountArrays,poolInputList)
    ########################################################################

    ################################
    pool.close()
    pool.join()
    ################################

    ###################################################################################
    ##############################  Accumulation starts  ##############################
    ###################################################################################
    # Accumulate the results coming from each chrom
    # allChroms_SignalArrayAndCountArrayList_List is a list of lists.

    signature2AccumulatedAllChromsSignalArrayDict, \
    signature2AccumulatedAllChromsCountArrayDict, \
    allSinglePointMutationsAccumulatedAllChromsSignalArray, \
    allSinglePointMutationsAccumulatedAllChromsCountArray, \
    allIndelsAccumulatedAllChromsSignalArray, \
    allIndelsAccumulatedAllChromsCountArray, \
    sample2Signature2AccumulatedAllChromsSignalArrayDict, \
    sample2Signature2AccumulatedAllChromsCountArrayDict, \
    sample2AllSinglePointMutationsAccumulatedAllChromsSignalArrayDict, \
    sample2AllSinglePointMutationsAccumulatedAllChromsCountArrayDict, \
    sample2AllIndelsAccumulatedAllChromsSignalArrayDict, \
    sample2AllIndelsAccumulatedAllChromsCountArrayDict = accumulateSplitArraysWithExtraSampleBased(
        signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
        sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
        allChroms_SignalArrayAndCountArrayList_List)
    ###################################################################################
    ##############################  Accumulation ends  ################################
    ###################################################################################


    ####################################################################################
    ########### Write All Single Mutations Average Nucleosome Occupancy starts #########
    ####################################################################################
    writeAverageNucleosomeOccupancyFiles(allSinglePointMutationsAccumulatedAllChromsSignalArray,allSinglePointMutationsAccumulatedAllChromsCountArray,outputDir, jobname, AGGREGATEDSUBSTITUTIONS)
    ####################################################################################
    ########### Write All Single Mutations Average Nucleosome Occupancy ends ###########
    ####################################################################################

    ####################################################################################
    ############## Write Signature Based Average Nucleosome Occupancy starts ###########
    ####################################################################################
    writeSignatureBasedAverageNucleosomeOccupancyFiles(signature2AccumulatedAllChromsSignalArrayDict,signature2AccumulatedAllChromsCountArrayDict,outputDir,jobname)
    ####################################################################################
    ############## Write Signature Based Average Nucleosome Occupancy ends #############
    ####################################################################################

    ####################################################################################
    ################ Write All Indels Average Nucleosome Occupancy starts ##############
    ####################################################################################
    writeAverageNucleosomeOccupancyFiles(allIndelsAccumulatedAllChromsSignalArray,allIndelsAccumulatedAllChromsCountArray, outputDir,jobname, AGGREGATEDINDELS)
    ####################################################################################
    ################ Write All Indels Average Nucleosome Occupancy ends ################
    ####################################################################################

    ####################################################################################
    writeSampleBasedSignatureBasedAverageNucleosomeOccupancyFiles(sample2Signature2AccumulatedAllChromsSignalArrayDict,sample2Signature2AccumulatedAllChromsCountArrayDict,outputDir,jobname)
    ####################################################################################

    ####################################################################################
    writeSampleBasedAverageNucleosomeOccupancyFiles(sample2AllSinglePointMutationsAccumulatedAllChromsSignalArrayDict,sample2AllSinglePointMutationsAccumulatedAllChromsCountArrayDict,outputDir,jobname,AGGREGATEDSUBSTITUTIONS)
    ####################################################################################

    ####################################################################################
    writeSampleBasedAverageNucleosomeOccupancyFiles(sample2AllIndelsAccumulatedAllChromsSignalArrayDict,sample2AllIndelsAccumulatedAllChromsCountArrayDict,outputDir,jobname,AGGREGATEDINDELS)
    ####################################################################################

#For all chromosome parallel starts
#March 7, 2019 ends
########################################################################################




########################################################################################
#If chr based subs or indels dataframes are too big we can use this version
#March 6, 2019 starts
#For each chromosome sequential starts
def nucleosome_occupancy_analysis_subs_indels_each_chrom_sequential(genome,outputDir,jobname,singlePointMutationsFilename,indelsFilename,nucleosomeFilename):

    ##########################################################################
    chromSizesDict, \
    chrNamesInNucleosomeList, \
    signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict, \
    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict, \
    samplesWithAtLeast10KMutations2NumberofMutationsDict = fillNecessaryDataStructures(genome,outputDir,jobname)
    ##########################################################################

    ##########################################################################
    signature2AccumulatedAllChromsSignalArrayDict, \
    signature2AccumulatedAllChromsCountArrayDict, \
    allSinglePointMutationsAccumulatedAllChromsSignalArray, \
    allSinglePointMutationsAccumulatedAllChromsCountArray, \
    allIndelsAccumulatedAllChromsSignalArray, \
    allIndelsAccumulatedAllChromsCountArray, \
    sample2Signature2AccumulatedAllChromsSignalArrayDict, \
    sample2Signature2AccumulatedAllChromsCountArrayDict, \
    sample2AllSinglePointMutationsAccumulatedAllChromsSignalArrayDict, \
    sample2AllSinglePointMutationsAccumulatedAllChromsCountArrayDict, \
    sample2AllIndelsAccumulatedAllChromsSignalArrayDict, \
    sample2AllIndelsAccumulatedAllChromsCountArrayDict = initializationForAllChroms(signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict)
    ##########################################################################


    ###################################################################################
    ##################  For each chromsome sequential starts ##########################
    ###################################################################################
    for chrLong in chrNamesInNucleosomeList:

        maximum_chrom_size = chromSizesDict[chrLong]
        print('For debug %s chromSize:%d' %(chrLong,maximum_chrom_size))

        #FIRST READ CHRBASED NUCLEOSOME OCCUPANCY
        nucleosomeFilenameWoExtension = nucleosomeFilename[0:-4]

        ##############################################################
        signalArrayFilename = '%s_signal_%s.npy' % (chrLong, nucleosomeFilenameWoExtension)
        chrBasedSignalNucleosmeFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB,NUCLEOSOME, CHRBASED, signalArrayFilename)

        # countArrayFilename = '%s_count_%s.npy' % (chrLong, nucleosomeFilenameWoExtension)
        # chrBasedCountNucleosmeFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME,CHRBASED, countArrayFilename)

        #################################################################################################################
        # if (os.path.exists(chrBasedSignalNucleosmeFile) and os.path.exists(chrBasedCountNucleosmeFile)):
        if (os.path.exists(chrBasedSignalNucleosmeFile)):
            chrbased_nucleosome_signal_array = np.load(chrBasedSignalNucleosmeFile)
            print('chromosome %s  -- signal_array_npy: %d in bytes %f in GB' % (chrLong,sys.getsizeof(chrbased_nucleosome_signal_array), sys.getsizeof(chrbased_nucleosome_signal_array)/GIGABYTE_IN_BYTES))

            # Do we need count array? Actually no.
            # chrbased_nucleosome_count_array = np.load(chrBasedCountNucleosmeFile)
            # print('chromosome %s  -- count_array_npy: %d in bytes %f in GB' % (chrLong, sys.getsizeof(chrbased_nucleosome_count_array), sys.getsizeof(chrbased_nucleosome_count_array) / GIGABYTE_IN_BYTES))

            #TODO: This is specific to our data right now
            #Nucleosomes have chrM
            #SinglePointMutations and Indels have chrMT
            if (chrLong=='chrM'):
                chrLong='chrMT'

            #THEN READ CHRBASED SINGLE POINT MUTATIONS
            chrBased_spms_df=readChrBasedMutationDF(outputDir,jobname,chrLong,singlePointMutationsFilename)
            print('chromosome %s  -- chrBased_spms_df: %d in Bytes %f in GigaBytes' %(chrLong,sys.getsizeof(chrBased_spms_df),sys.getsizeof(chrBased_spms_df)/GIGABYTE_IN_BYTES))

            #THEN READ CHRBASED INDELS
            chrBased_indels_df = readChrBasedIndelsDF(outputDir,jobname,chrLong,indelsFilename)
            print('chromosome %s  -- chrBased_indels_df: %d in Bytes %f in GigaBytes' %(chrLong,sys.getsizeof(chrBased_indels_df),sys.getsizeof(chrBased_indels_df)/GIGABYTE_IN_BYTES))

            ##############################################################
            signature2SignalArrayDict, \
            signature2CountArrayDict, \
            allSinglePointMutationsSignalArray, \
            allSinglePointMutationsCountArray, \
            allIndelsSignalArray, \
            allIndelsCountArray, \
            sample2Signature2SignalArrayDict, \
            sample2Signature2CountArrayDict, \
            sample2AllSinglePointMutationsSignalArrayDict, \
            sample2AllSinglePointMutationsCountArrayDict, \
            sample2AllIndelsSignalArrayDict, \
            sample2AllIndelsCountArrayDict= initializationForEachChrom(signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict)
            ##############################################################

            ###############################################################################
            ################### Fill signal and count array starts ########################
            ###############################################################################

            if ((chrBased_spms_df is not None) and (not chrBased_spms_df.empty)):
                # Fill for single point mutations
                chrBased_spms_df.apply(fillSplitBasedSignalArrayAndCountArrayForSPMsWithExtraSampleBased,
                                       nucleosome_array=chrbased_nucleosome_signal_array,
                                       maximum_chrom_size = maximum_chrom_size,
                                       signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict=signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
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

            if ((chrBased_indels_df is not None) and (not chrBased_indels_df.empty)):
                # Fill for indels
                chrBased_indels_df.apply(fillSplitBasedSignalArrayAndCountArrayForIndelsWithExtraSampleBased,
                                         nucleosome_array=chrbased_nucleosome_signal_array,
                                         maximum_chrom_size=maximum_chrom_size,
                                         sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict=sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
                                         allIndelsSignalArray=allIndelsSignalArray,
                                         allIndelsCountArray=allIndelsCountArray,
                                         sample2AllIndelsSignalArrayDict=sample2AllIndelsSignalArrayDict,
                                         sample2AllIndelsCountArrayDict=sample2AllIndelsCountArrayDict,
                                         axis=1)
            ###############################################################################
            ################### Fill signal and count array ends ##########################
            ###############################################################################

            #Accumulate chr based dictionaries
            ###############################################################################################################
            ###################  Accumulate the results coming from each chromosome starts  ###############################
            ###############################################################################################################

            # Accumulate right in the left
            #################################################
            accumulateSignatureBasedArrays(signature2AccumulatedAllChromsSignalArrayDict,signature2SignalArrayDict)
            accumulateSignatureBasedArrays(signature2AccumulatedAllChromsCountArrayDict,signature2CountArrayDict)

            ##Accumulate the results coming from each chromosome
            allSinglePointMutationsAccumulatedAllChromsSignalArray += allSinglePointMutationsSignalArray
            allSinglePointMutationsAccumulatedAllChromsCountArray += allSinglePointMutationsCountArray

            allIndelsAccumulatedAllChromsSignalArray += allIndelsSignalArray
            allIndelsAccumulatedAllChromsCountArray += allIndelsCountArray
            #################################################

            ###########################Accumulate sample based starts ####################################
            ##Accumulate the results coming from each chromosome
            accumulateSampleBasedSignatureBasedArrays(sample2Signature2AccumulatedAllChromsSignalArrayDict,sample2Signature2SignalArrayDict)
            accumulateSampleBasedSignatureBasedArrays(sample2Signature2AccumulatedAllChromsCountArrayDict,sample2Signature2CountArrayDict)

            accumulateSampleBasedArrays(sample2AllSinglePointMutationsAccumulatedAllChromsSignalArrayDict,sample2AllSinglePointMutationsSignalArrayDict)
            accumulateSampleBasedArrays(sample2AllSinglePointMutationsAccumulatedAllChromsCountArrayDict,sample2AllSinglePointMutationsCountArrayDict)

            accumulateSampleBasedArrays(sample2AllIndelsAccumulatedAllChromsSignalArrayDict,sample2AllIndelsSignalArrayDict)
            accumulateSampleBasedArrays(sample2AllIndelsAccumulatedAllChromsCountArrayDict,sample2AllIndelsCountArrayDict)
            ###########################Accumulate sample based ends ######################################

            ###############################################################################################################
            ###################  Accumulate the results coming from each chromosome ends  #################################
            ###############################################################################################################

            ######################################################################################
            ##############   SPMs and Indels with extra sampleBased  analysis ends ###############
            ######################################################################################

    ###################################################################################
    ##################  For each chromsome sequential ends ############################
    ###################################################################################



    ####################################################################################
    ########### Write All Single Mutations Average Nucleosome Occupancy starts #########
    ####################################################################################
    writeAverageNucleosomeOccupancyFiles(allSinglePointMutationsAccumulatedAllChromsSignalArray,allSinglePointMutationsAccumulatedAllChromsCountArray,outputDir, jobname, AGGREGATEDSUBSTITUTIONS)
    ####################################################################################
    ########### Write All Single Mutations Average Nucleosome Occupancy ends ###########
    ####################################################################################

    ####################################################################################
    ############## Write Signature Based Average Nucleosome Occupancy starts ###########
    ####################################################################################
    writeSignatureBasedAverageNucleosomeOccupancyFiles(signature2AccumulatedAllChromsSignalArrayDict,signature2AccumulatedAllChromsCountArrayDict,outputDir,jobname)
    ####################################################################################
    ############## Write Signature Based Average Nucleosome Occupancy ends #############
    ####################################################################################

    ####################################################################################
    ################ Write All Indels Average Nucleosome Occupancy starts ##############
    ####################################################################################
    writeAverageNucleosomeOccupancyFiles(allIndelsAccumulatedAllChromsSignalArray,allIndelsAccumulatedAllChromsCountArray, outputDir,jobname, AGGREGATEDINDELS)
    ####################################################################################
    ################ Write All Indels Average Nucleosome Occupancy ends ################
    ####################################################################################

    ####################################################################################
    writeSampleBasedSignatureBasedAverageNucleosomeOccupancyFiles(sample2Signature2AccumulatedAllChromsSignalArrayDict,sample2Signature2AccumulatedAllChromsCountArrayDict,outputDir,jobname)
    ####################################################################################

    ####################################################################################
    writeSampleBasedAverageNucleosomeOccupancyFiles(sample2AllSinglePointMutationsAccumulatedAllChromsSignalArrayDict,sample2AllSinglePointMutationsAccumulatedAllChromsCountArrayDict,outputDir,jobname,AGGREGATEDSUBSTITUTIONS)
    ####################################################################################

    ####################################################################################
    writeSampleBasedAverageNucleosomeOccupancyFiles(sample2AllIndelsAccumulatedAllChromsSignalArrayDict,sample2AllIndelsAccumulatedAllChromsCountArrayDict,outputDir,jobname,AGGREGATEDINDELS)
    ####################################################################################


#For each chromosome sequential ends
#March 6, 2019 ends
########################################################################################


