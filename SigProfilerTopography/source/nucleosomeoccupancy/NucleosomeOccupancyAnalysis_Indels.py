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
#   for all indels
###############################################################################################################

import os
import sys

#############################################################
current_abs_path = os.path.abspath(os.path.dirname(__file__))
print('NucleosomeOccupancyAnalysis_Indels.py current_abs_path:%s' %(current_abs_path))
#############################################################

commonsPath = os.path.join(current_abs_path, '..','commons')
sys.path.append(commonsPath)

from SigProfilerTopography.source.commons.TopographyCommons import *

########################################################################################
#This method is specific to this python file.
def fillSignalArrayAndCountArrayWithExtraSampleBased(inputList):
    chrBased_indels_df = inputList[0]
    chrBased_nucleosome_df_split = inputList[1]
    chrLong = inputList[2]
    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = inputList[3]


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
    # Initialization for all indels for split
    allIndelsSignalArray = np.zeros(windowSize)
    allIndelsCountArray = np.zeros(windowSize, dtype=int)
    ##############################################

    ##############################################
    # Initialization sample based
    sample2AllIndelsSignalArrayDict = {}
    sample2AllIndelsCountArrayDict = {}

    for sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
        sample2AllIndelsSignalArrayDict[sample] = np.zeros(windowSize)
        sample2AllIndelsCountArrayDict[sample] = np.zeros(windowSize, dtype=int)
    ##############################################

    ###############################################################################
    ################### Initialization2 for each split ends #######################
    ###############################################################################


    ###############################################################################
    ################### Fill signal and count array starts ########################
    ###############################################################################
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

    # Append indels arrays
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(allIndelsSignalArray)
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(allIndelsCountArray)

    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(sample2AllIndelsSignalArrayDict)
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(sample2AllIndelsCountArrayDict)

    return split_chrBased_SignalArrayAndCountArray_DictionaryList
########################################################################################

########################################################################################
#This method is specific to this python file.
def fillSignalArrayAndCountArray(inputList):
    chrBased_indels_df = inputList[0]
    chrBased_nucleosome_df_split = inputList[1]
    chrLong = inputList[2]

    ###########################################################################
    # Initialize chrBased_nucleosome_split_array with np array of zeros
    # Fill the nucleosome split array using chrBased_nucleosome_df_split
    chrBased_nucleosome_split_array = np.zeros(MAXIMUM_CHROMOSOME_LENGTH, dtype=np.float32)
    chrBased_nucleosome_df_split.apply(fillNucleosomeSignalArray, nucleosome_array=chrBased_nucleosome_split_array, axis=1)
    ###########################################################################

    ###########################################################################
    print('#### debug starts ####')
    print('chr:%s -- '
            'size(chrBased_indels_df):%d size(chrBased_indels_df):%f GB --- '
            'size(chrBased_nucleosome_df_split):%d size(chrBased_nucleosome_df_split):%f GB  --- '
            'size(chrBased_nucleosome_split_array):%d size(chrBased_nucleosome_split_array):%f GB' %(chrLong,
            sys.getsizeof(chrBased_indels_df),sys.getsizeof(chrBased_indels_df)/GIGABYTE_IN_BYTES,
            sys.getsizeof(chrBased_nucleosome_df_split),sys.getsizeof(chrBased_nucleosome_df_split) / GIGABYTE_IN_BYTES,
            sys.getsizeof(chrBased_nucleosome_split_array), sys.getsizeof(chrBased_nucleosome_split_array)/ GIGABYTE_IN_BYTES ))
    print('#### debug ends ####')
    ###########################################################################

    ##############################################
    # Initialization for all indels for split
    allIndelsSignalArray = np.zeros(windowSize)
    allIndelsCountArray = np.zeros(windowSize, dtype=int)
    ##############################################

    #Initialzie the list, you will return this list
    split_chrBased_SignalArrayAndCountArray_DictionaryList = []

    #Fill for indels
    chrBased_indels_df.apply(fillSplitBasedSignalArrayAndCountArrayForIndels, nucleosome_array=chrBased_nucleosome_split_array, allIndelsSignalArray=allIndelsSignalArray,allIndelsCountArray=allIndelsCountArray, axis=1)

    #Append indels arrays
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(allIndelsSignalArray)
    split_chrBased_SignalArrayAndCountArray_DictionaryList.append(allIndelsCountArray)

    return split_chrBased_SignalArrayAndCountArray_DictionaryList
########################################################################################

########################################################################################
def accumulateSplitArraysWithExtraSampleBased(sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,allSplits_chrBased_SignalArrayAndCountArray_DictionaryList):

    #Initialize them
    #Fill them using allSplits_chrBased_SignalArrayAndCountArray_DictionaryList
    #Return them

    ########################################################################
    ############ Initialization3 for all splits starts #####################
    ########################################################################
    allIndelsAccumulatedSplitsSignalArray = np.zeros(windowSize)
    allIndelsAccumulatedSplitsCountArray = np.zeros(windowSize, dtype=int)

    sample2AllIndelsAccumulatedSplitsChrBasedSignalArrayDict = {}
    sample2AllIndelsAccumulatedSplitsChrBasedCountArrayDict = {}

    for sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
        sample2AllIndelsAccumulatedSplitsChrBasedSignalArrayDict[sample] = np.zeros(windowSize)
        sample2AllIndelsAccumulatedSplitsChrBasedCountArrayDict[sample] = np.zeros(windowSize, dtype=int)
    ########################################################################
    ############ Initialization3 for all splits ends #######################
    ########################################################################

    ##############################################
    for split_chrBased_SignalArrayAndCountArray_DictionaryList in allSplits_chrBased_SignalArrayAndCountArray_DictionaryList:

        print('For Debug Feb 22, 2019 len(split_chrBased_SignalArrayAndCountArray_DictionaryList)')
        print(len(split_chrBased_SignalArrayAndCountArray_DictionaryList))

        allIndelsSplitSignalArray = split_chrBased_SignalArrayAndCountArray_DictionaryList[0]
        allIndelsSplitCountArray = split_chrBased_SignalArrayAndCountArray_DictionaryList[1]

        sample2AllIndelsSplitSignalArrayDict = split_chrBased_SignalArrayAndCountArray_DictionaryList[2]
        sample2AllIndelsSplitCountArrayDict = split_chrBased_SignalArrayAndCountArray_DictionaryList[3]

        #####################################################################################################
        ######################### Accumulate right in the left starts  ######################################
        #####################################################################################################

        ######################### Accumulate starts ######################################
        allIndelsAccumulatedSplitsSignalArray += allIndelsSplitSignalArray
        allIndelsAccumulatedSplitsCountArray += allIndelsSplitCountArray
        ######################### Accumulate starts ######################################

        ############################Accumulate Sample Based starts ###################################
        accumulateSampleBasedArrays(sample2AllIndelsAccumulatedSplitsChrBasedSignalArrayDict,sample2AllIndelsSplitSignalArrayDict)
        accumulateSampleBasedArrays(sample2AllIndelsAccumulatedSplitsChrBasedCountArrayDict,sample2AllIndelsSplitCountArrayDict)
        ############################Accumulate Sample Based ends #####################################

        #####################################################################################################
        ######################### Accumulate right in the left ends  ########################################
        #####################################################################################################
    ##############################################

    return  allIndelsAccumulatedSplitsSignalArray,\
            allIndelsAccumulatedSplitsCountArray,\
            sample2AllIndelsAccumulatedSplitsChrBasedSignalArrayDict,\
            sample2AllIndelsAccumulatedSplitsChrBasedCountArrayDict
########################################################################################


########################################################################################
#This function is specific to this python file
def accumulateSplitArrays(allSplits_chrBased_SignalArrayAndCountArray_DictionaryList):

    ##############################################
    # Fill these indels arrays
    allIndelsAccumulatedSplitsSignalArray = np.zeros(windowSize)
    allIndelsAccumulatedSplitsCountArray = np.zeros(windowSize, dtype=int)
    ##############################################

    ##############################################
    for split_chrBased_SignalArrayAndCountArray_DictionaryList in allSplits_chrBased_SignalArrayAndCountArray_DictionaryList:

        ############################################################################################
        allIndelsSplitSignalArray = split_chrBased_SignalArrayAndCountArray_DictionaryList[0]
        allIndelsSplitCountArray = split_chrBased_SignalArrayAndCountArray_DictionaryList[1]
        ############################################################################################

        ############################################################################################
        #Accumulate Indels
        allIndelsAccumulatedSplitsSignalArray += allIndelsSplitSignalArray
        allIndelsAccumulatedSplitsCountArray += allIndelsSplitCountArray
        ############################################################################################

    return allIndelsAccumulatedSplitsSignalArray, allIndelsAccumulatedSplitsCountArray
########################################################################################


########################################################################################
#main funcyion
def nucleosomeOccupancyAnalysis_Indels(jobname,indelsFilename,nucleosomeFilename):

    withExtraSampleBasedAnalysis = True

    ##########################################################################################
    #Load Chrnames from nucleosome file
    chrNamesInNucleosomeList = []
    chrNamesInNucleosomeFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME, ChrNamesInNucleosomesFilename)
    if (os.path.exists(chrNamesInNucleosomeFilePath)):
        chrNamesInNucleosomeArray = np.loadtxt(chrNamesInNucleosomeFilePath,dtype=str, delimiter='\t')
        chrNamesInNucleosomeList = list(chrNamesInNucleosomeArray)

    #Load sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict
    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict= {}

    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath = os.path.join(current_abs_path,ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT,jobname,DATA,Sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilename)

    if (os.path.exists(sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath)):
        sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = readDictionary(sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath)
    ##########################################################################################

    ########################################################################################
    # Initialization1 for indels
    allIndelsAccumulatedAllChromsSignalArray = np.zeros(windowSize)
    allIndelsAccumulatedAllChromsCountArray = np.zeros(windowSize, dtype=int)

    sample2AllIndelsAccumulatedAllChromsSignalArrayDict = {}
    sample2AllIndelsAccumulatedAllChromsCountArrayDict = {}

    for sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
        sample2AllIndelsAccumulatedAllChromsSignalArrayDict[sample] = np.zeros(windowSize)
        sample2AllIndelsAccumulatedAllChromsCountArrayDict[sample] = np.zeros(windowSize, dtype=int)
    ########################################################################################

    ######################################################
    print('chrNamesInNucleosomeList')
    print(chrNamesInNucleosomeList)
    # short ones e.g.: chrNamesList = ['1','10','11','12','13', '14', '15', '16', '17', '18', '19', '2', '20', '21', '22', '3', '4', '5', '6', '7', '8', '9', 'MT', 'X', 'Y']
    ######################################################

    ######################################################
    #Better practice
    numofProcesses = multiprocessing.cpu_count()
    print('multiprocessing.cpu_count() --> number of processes: %d' %numofProcesses)
    pool = multiprocessing.Pool(numofProcesses)
    ######################################################

    ###################################################################################
    ##################  For each chromsome sequential starts ##########################
    ###################################################################################
    for chrLong in chrNamesInNucleosomeList:

        #FIRST READ CHRBASED NUCLEOSOME OCCUPANCY
        chrBased_nucleosome_df = readChrBasedNuclesomeDF(chrLong,nucleosomeFilename)

        #TODO: This is specific to out data right now
        #Nucleosome has chrM
        #SinglePointMutations and Indels have chrMT
        if (chrLong=='chrM'):
            chrLong='chrMT'

        #THEN READ CHRBASED INDELS
        chrBased_indels_df = readChrBasedIndelsDF(jobname,chrLong,indelsFilename)

        print('chromosome %s  -- chrBased_nucleosome_df: %d in bytes %f in GB' % (chrLong,sys.getsizeof(chrBased_nucleosome_df), sys.getsizeof(chrBased_nucleosome_df)/GIGABYTE_IN_BYTES))

        if (chrBased_indels_df is not None):

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
            chrBased_nucleosome_df_splits = np.array_split(chrBased_nucleosome_df,numofProcesses)

            #Prepare the poolInputList
            poolInputList = []

            #########################################################################
            if (withExtraSampleBasedAnalysis):

                ######################################################################################
                ####################     Indels with extra sampleBased  analysis starts ##############
                ######################################################################################

                # For each split chrBased mutation df do it in parallel
                for chrBased_nucleosome_df_split in chrBased_nucleosome_df_splits:
                    inputList = []

                    # Question: Does each split has header?
                    # Answer: Yes
                    inputList.append(chrBased_indels_df)  # Same chrBased_indels_df
                    inputList.append(chrBased_nucleosome_df_split)  #Different chrBased_nucleosome_df split
                    inputList.append(chrLong) #same chrLong
                    inputList.append(sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict) #Sep 19, 2018
                    poolInputList.append(inputList)

                #########################################################################
                allSplits_chrBased_SignalArrayAndCountArray_DictionaryList = pool.map(fillSignalArrayAndCountArrayWithExtraSampleBased,poolInputList)


                #Accumulate the results coming from each split
                #allSplits_chrBased_SignalArrayAndCountArray_DictionaryList is a list of lists.
                #Each list consist of two dictionaries
                allIndelsAccumulatedSplitsSignalArray, \
                allIndelsAccumulatedSplitsCountArray, \
                sample2AllIndelsAccumulatedSplitsSignalArrayDict, \
                sample2AllIndelsAccumulatedSplitsCountArrayDict = accumulateSplitArraysWithExtraSampleBased(
                    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
                    allSplits_chrBased_SignalArrayAndCountArray_DictionaryList)
                #################################################


                ###############################################################################################################
                ###################  Accumulate the results coming from each chromosome starts  ###############################
                ###############################################################################################################

                #Accumulate right in the left
                #################################################

                ##Accumulate the results coming from each chromosome
                allIndelsAccumulatedAllChromsSignalArray += allIndelsAccumulatedSplitsSignalArray
                allIndelsAccumulatedAllChromsCountArray += allIndelsAccumulatedSplitsCountArray
                #################################################

                ###########################Accumulate sample based starts ####################################
                ##Accumulate the results coming from each chromosome
                accumulateSampleBasedArrays(sample2AllIndelsAccumulatedAllChromsSignalArrayDict,sample2AllIndelsAccumulatedSplitsSignalArrayDict)
                accumulateSampleBasedArrays(sample2AllIndelsAccumulatedAllChromsCountArrayDict,sample2AllIndelsAccumulatedSplitsCountArrayDict)
                ###########################Accumulate sample based ends ######################################

                ###############################################################################################################
                ###################  Accumulate the results coming from each chromosome ends  #################################
                ###############################################################################################################


                ######################################################################################
                ######################## Indels with extra sampleBased  analysis ends ################
                ######################################################################################

            #########################################################################

            #########################################################################
            else:

                ######################################################################################
                ######################## Indels with extra sampleBased  analysis starts ##############
                ######################################################################################

                #########################################################################
                # For each split chrBased mutation df do it in parallel
                for chrBased_nucleosome_df_split in chrBased_nucleosome_df_splits:
                    inputList = []

                    # Question: Does each split has header?
                    # Answer: Yes
                    inputList.append(chrBased_indels_df)  # Same chrBased_indels_df
                    inputList.append(chrBased_nucleosome_df_split)  # Different chrBased_nucleosome_df split
                    inputList.append(chrLong)  # same chrLong
                    poolInputList.append(inputList)
                #########################################################################

                #########################################################################
                allSplits_chrBased_SignalArrayAndCountArray_DictionaryList = pool.map(fillSignalArrayAndCountArray,poolInputList)

                # Accumulate the results coming from each split
                # allSplits_chrBased_SignalArrayAndCountArray_DictionaryList is a list of lists.
                # Each list consist of two dictionaries
                allIndelsAccumulatedSplitsSignalArray, allIndelsAccumulatedSplitsCountArray = accumulateSplitArrays(allSplits_chrBased_SignalArrayAndCountArray_DictionaryList)
                #################################################

                #################################################
                ##Accumulate the indels results coming from each chromosome
                allIndelsAccumulatedAllChromsSignalArray += allIndelsAccumulatedSplitsSignalArray
                allIndelsAccumulatedAllChromsCountArray += allIndelsAccumulatedSplitsCountArray
                #################################################

                ######################################################################################
                ######################## Indels with extra sampleBased  analysis ends ################
                ######################################################################################

            #########################################################################

    ###################################################################################
    ##################  For each chromsome sequential ends ############################
    ###################################################################################

    ################################
    pool.close()
    pool.join()
    ################################


    ####################################################################################
    ################ Write All Indels Average Nucleosome Occupancy starts ##############
    ####################################################################################


    writeAverageNucleosomeOccupancyFiles(allIndelsAccumulatedAllChromsSignalArray,allIndelsAccumulatedAllChromsCountArray, jobname, AGGREGATEDINDELS)
    ####################################################################################
    ################ Write All Indels Average Nucleosome Occupancy ends ################
    ####################################################################################

    if (withExtraSampleBasedAnalysis):

        ####################################################################################
        writeSampleBasedAverageNucleosomeOccupancyFiles(sample2AllIndelsAccumulatedAllChromsSignalArrayDict,sample2AllIndelsAccumulatedAllChromsCountArrayDict,jobname,AGGREGATEDINDELS)
        ####################################################################################

########################################################################################
