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
#   for subs, indels and dinucs sample based and all samples pooled
#   for all subs signatures with all single point mutations with a certain probability for that signature
#   for all indels signatures with all indels with a certain probability for that signature
#   for all dinucs signatures with all dinucs with a certain probability for that signature
###############################################################################################################

# #############################################################
# current_abs_path = os.path.abspath(os.path.dirname(__file__))
# commonsPath = os.path.join(current_abs_path, '..','commons')
# sys.path.append(commonsPath)
# #############################################################

from SigProfilerTopography.source.commons.TopographyCommons import *

##############################################################################################################
#main function
def nucleosomeOccupancyAnalysis(mutationTypes,computationType,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,nucleosomeFilename):
    print('########################## NucleosomeOccupancyAnalysis starts ##########################')
    if  (computationType == COMPUTATION_ALL_CHROMOSOMES_PARALLEL):
        nucleosome_occupancy_analysis_all_chroms_parallel(mutationTypes,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,nucleosomeFilename)
    elif (computationType == COMPUTATION_CHROMOSOMES_SEQUENTIAL):
        nucleosome_occupancy_analysis_each_chrom_sequential(mutationTypes,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,nucleosomeFilename)
    print('########################## NucleosomeOccupancyAnalysis ends ############################')
##############################################################################################################


########################################################################################
def fillSignalArrayAndCountArrays(inputList):
    chrbased_nucleosome_signal_array = inputList[0]
    combined_chrBased_subs_df = inputList[1]
    combined_chrBased_indels_df =  inputList[2]
    combined_chrBased_dinucs_df = inputList[3]
    maximum_chrom_size =  inputList[4]
    sample2NumberofSubsDict = inputList[5]
    sample2NumberofIndelsDict = inputList[6]
    sample2NumberofDinucsDict = inputList[7]
    subsSignature2NumberofMutationsDict =  inputList[8]
    indelsSignature2NumberofMutationsDict = inputList[9]
    dinucsSignature2NumberofMutationsDict = inputList[10]
    sample2SubsSignature2NumberofMutationsDict =  inputList[11]
    sample2IndelsSignature2NumberofMutationsDict = inputList[12]
    sample2DinucsSignature2NumberofMutationsDict = inputList[13]

    ##############################################################
    # old way
    # type2SignalArrayDict = {}
    # type2CountArrayDict = {}
    # sample2Type2SignalArrayDict = {}
    # sample2Type2CountArrayDict = {}
    ##############################################################

    ##############################################################
    simNum2Type2SignalArrayDict = {}
    simNum2Type2CountArrayDict = {}
    simNum2Sample2Type2SignalArrayDict = {}
    simNum2Sample2Type2CountArrayDict = {}
    ##############################################################

    ###############################################################################
    ################### Fill signal and count array starts ########################
    ###############################################################################
    #Fill for single point mutations
    if ((combined_chrBased_subs_df is not None) and (not combined_chrBased_subs_df.empty)):
        combined_chrBased_subs_df.apply(fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated(),
                           nucleosome_array=chrbased_nucleosome_signal_array,
                           maximum_chrom_size =maximum_chrom_size,
                           sample2NumberofMutationsDict = sample2NumberofSubsDict,
                           signature2NumberofMutationsDict = subsSignature2NumberofMutationsDict,
                           sample2Signature2NumberofMutationsDict=sample2SubsSignature2NumberofMutationsDict,
                            simNum2Type2SignalArrayDict=simNum2Type2SignalArrayDict,
                            simNum2Type2CountArrayDict = simNum2Type2CountArrayDict,
                            simNum2Sample2Type2SignalArrayDict = simNum2Sample2Type2SignalArrayDict,
                            simNum2Sample2Type2CountArrayDict = simNum2Sample2Type2CountArrayDict,
                            MUTATION_SIGNATURE_PROBABILITY_THRESHOLD = SUBSTITUTION_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD,
                           type = AGGREGATEDSUBSTITUTIONS,
                           axis=1)

    #Fill for indels
    if ((combined_chrBased_indels_df is not None) and (not combined_chrBased_indels_df.empty)):
        combined_chrBased_indels_df.apply(fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated,
                                 nucleosome_array=chrbased_nucleosome_signal_array,
                                 maximum_chrom_size=maximum_chrom_size,
                                 sample2NumberofMutationsDict = sample2NumberofIndelsDict,
                                 signature2NumberofMutationsDict=indelsSignature2NumberofMutationsDict,
                                 sample2Signature2NumberofMutationsDict = sample2IndelsSignature2NumberofMutationsDict,
                                  simNum2Type2SignalArrayDict=simNum2Type2SignalArrayDict,
                                  simNum2Type2CountArrayDict=simNum2Type2CountArrayDict,
                                  simNum2Sample2Type2SignalArrayDict=simNum2Sample2Type2SignalArrayDict,
                                  simNum2Sample2Type2CountArrayDict=simNum2Sample2Type2CountArrayDict,
                                  MUTATION_SIGNATURE_PROBABILITY_THRESHOLD=INDEL_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD,
                                 type=AGGREGATEDINDELS,
                                 axis=1)

    #Fill for dinucs
    if ((combined_chrBased_dinucs_df is not None) and (not combined_chrBased_dinucs_df.empty)):
        combined_chrBased_dinucs_df.apply(fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated,
            nucleosome_array=chrbased_nucleosome_signal_array,
            maximum_chrom_size=maximum_chrom_size,
            sample2NumberofMutationsDict=sample2NumberofDinucsDict,
            signature2NumberofMutationsDict=dinucsSignature2NumberofMutationsDict,
            sample2Signature2NumberofMutationsDict=sample2DinucsSignature2NumberofMutationsDict,
          simNum2Type2SignalArrayDict=simNum2Type2SignalArrayDict,
          simNum2Type2CountArrayDict=simNum2Type2CountArrayDict,
          simNum2Sample2Type2SignalArrayDict=simNum2Sample2Type2SignalArrayDict,
          simNum2Sample2Type2CountArrayDict=simNum2Sample2Type2CountArrayDict,
          MUTATION_SIGNATURE_PROBABILITY_THRESHOLD=DINUC_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD,
            type=AGGREGATEDDINUCS,
            axis=1)

    ###############################################################################
    ################### Fill signal and count array ends ##########################
    ###############################################################################


    ###############################################################################
    ################### Return  starts ############################################
    ###############################################################################
    # Initialzie the list, you will return this list
    combined_chrBased_SignalArrayAndCountArray_List = []
    combined_chrBased_SignalArrayAndCountArray_List.append(simNum2Type2SignalArrayDict)
    combined_chrBased_SignalArrayAndCountArray_List.append(simNum2Type2CountArrayDict)
    combined_chrBased_SignalArrayAndCountArray_List.append(simNum2Sample2Type2SignalArrayDict)
    combined_chrBased_SignalArrayAndCountArray_List.append(simNum2Sample2Type2CountArrayDict)

    return combined_chrBased_SignalArrayAndCountArray_List
    ###############################################################################
    ################### Return  ends ##############################################
    ###############################################################################
########################################################################################


########################################################################################
def accumulateSignalCountArrays(all_partials_chrBased_SignalArrayAndCountArray_DictionaryList):

    #Initialize them
    type2AccumulatedSignalArrayDict = {}
    type2AccumulatedCountArrayDict = {}
    sample2Type2AccumulatedSignalArrayDict = {}
    sample2Type2AccumulatedCountArrayDict = {}

    #Fill them
    for partial_chrBased_SignalArrayAndCountArray_DictionaryList in all_partials_chrBased_SignalArrayAndCountArray_DictionaryList:

        ######################### Accumulate right in the left starts  ######################################
        type2PartialSignalArrayDict = partial_chrBased_SignalArrayAndCountArray_DictionaryList[0]
        type2PartialCountArrayDict = partial_chrBased_SignalArrayAndCountArray_DictionaryList[1]

        sample2Type2PartialSignalArrayDict = partial_chrBased_SignalArrayAndCountArray_DictionaryList[2]
        sample2Type2PartialCountArrayDict = partial_chrBased_SignalArrayAndCountArray_DictionaryList[3]
        ######################### Accumulate right in the left ends  ########################################

        ######################### Accumulate right in the left starts  ######################################
        accumulateTypeBasedArrays(type2AccumulatedSignalArrayDict,type2PartialSignalArrayDict)
        accumulateTypeBasedArrays(type2AccumulatedCountArrayDict,type2PartialCountArrayDict)

        accumulateSampleBasedTypeBasedArrays(sample2Type2AccumulatedSignalArrayDict,sample2Type2PartialSignalArrayDict)
        accumulateSampleBasedTypeBasedArrays(sample2Type2AccumulatedCountArrayDict,sample2Type2PartialCountArrayDict)
        ######################### Accumulate right in the left ends  ########################################

    return  type2AccumulatedSignalArrayDict, \
            type2AccumulatedCountArrayDict, \
            sample2Type2AccumulatedSignalArrayDict, \
            sample2Type2AccumulatedCountArrayDict
########################################################################################

########################################################################################
#For all chromosome parallel starts
def nucleosome_occupancy_analysis_all_chroms_parallel(mutationTypes,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,nucleosomeFilename):

    ##########################################################################
    sample2NumberofSubsDict = getSample2NumberofSubsDict(outputDir,jobname)
    sample2NumberofIndelsDict = getSample2NumberofIndelsDict(outputDir,jobname)
    sample2NumberofDinucsDict = getDictionary(outputDir,jobname, Sample2NumberofDinucsDictFilename)

    subsSignature2NumberofMutationsDict = getSubsSignature2NumberofMutationsDict(outputDir,jobname)
    indelsSignature2NumberofMutationsDict = getIndelsSignature2NumberofMutationsDict(outputDir,jobname)
    dinucsSignature2NumberofMutationsDict = getDictionary(outputDir, jobname, DinucsSignature2NumberofMutationsDictFilename)

    sample2SubsSignature2NumberofMutationsDict = getSample2SubsSignature2NumberofMutationsDict(outputDir,jobname)
    sample2IndelsSignature2NumberofMutationsDict = getSample2IndelsSignature2NumberofMutationsDict(outputDir,jobname)
    sample2DinucsSignature2NumberofMutationsDict = getDictionary(outputDir, jobname,Sample2DinucsSignature2NumberofMutationsDictFilename)
    ##########################################################################

    ##########################################################################
    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)
    poolInputList = []
    ##########################################################################

    ###################################################################################
    ##################  For all chromsomes parallel starts ############################
    ###################################################################################
    for chrLong in chromNamesList:
        inputList = []

        chromSize = chromSizesDict[chrLong]

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

            #Read original data
            #READ CHRBASED SINGLE POINT MUTATIONS
            original_chrBased_spms_df=readChrBasedSubsDF(outputDir,jobname,chrLong,SUBS,0)
            if (original_chrBased_spms_df is not None):
                print('chromosome %s  -- chrBased_spms_df: %d in Bytes %f in GigaBytes' %(chrLong,sys.getsizeof(original_chrBased_spms_df),sys.getsizeof(original_chrBased_spms_df)/GIGABYTE_IN_BYTES))

            #READ CHRBASED INDELS
            original_chrBased_indels_df = readChrBasedMutationsDF(outputDir,jobname,chrLong,INDELS,0)
            if (original_chrBased_indels_df is not None):
                print('chromosome %s  -- chrBased_indels_df: %d in Bytes %f in GigaBytes' %(chrLong,sys.getsizeof(original_chrBased_indels_df),sys.getsizeof(original_chrBased_indels_df)/GIGABYTE_IN_BYTES))

            #READ CHRBASED_DINUCS
            original_chrBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong,DINUCS,0)
            if (original_chrBased_dinucs_df is not None):
                print('chromosome %s  -- chrBased_dinucss_df: %d in Bytes %f in GigaBytes' %(chrLong, sys.getsizeof(original_chrBased_dinucs_df), sys.getsizeof(original_chrBased_dinucs_df) / GIGABYTE_IN_BYTES))

            combined_chrBased_subs_df= getCombinedChrBasedDF(outputDir, jobname, chrLong,original_chrBased_spms_df,SUBS,numofSimulations)
            combined_chrBased_indels_df = getCombinedChrBasedDF(outputDir, jobname, chrLong,original_chrBased_indels_df,INDELS,numofSimulations)
            combined_chrBased_dinucs_df= getCombinedChrBasedDF(outputDir, jobname, chrLong,original_chrBased_dinucs_df,DINUCS,numofSimulations)

            inputList.append(chrbased_nucleosome_signal_array)
            inputList.append(combined_chrBased_subs_df)
            inputList.append(combined_chrBased_indels_df)
            inputList.append(combined_chrBased_dinucs_df)
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
            poolInputList.append(inputList)
    ###################################################################################
    ##################  For all chromsomes parallel ends ##############################
    ###################################################################################

    ########################################################################
    allChroms_SignalArrayAndCountArrayList_List = pool.map(fillSignalArrayAndCountArrays, poolInputList)
    ########################################################################

    ################################
    pool.close()
    pool.join()
    ################################

    ###################################################################################
    ##############################  Accumulation starts  ##############################
    ###################################################################################
    #TODO update for combines
    # Accumulate the results coming from each chrom
    # allChroms_SignalArrayAndCountArrayList_List is a list of lists.

    type2AccumulatedSignalArrayDict, \
    type2AccumulatedCountArrayDict, \
    sample2Type2AccumulatedSignalArrayDict, \
    sample2Type2AccumulatedCountArrayDict = accumulateSignalCountArrays(allChroms_SignalArrayAndCountArrayList_List)
    ###################################################################################
    ##############################  Accumulation ends  ################################
    ###################################################################################

    #TODO update for combines
    writeAverageNucleosomeOccupancy(type2AccumulatedSignalArrayDict,type2AccumulatedCountArrayDict,sample2Type2AccumulatedSignalArrayDict,sample2Type2AccumulatedCountArrayDict,outputDir,jobname)

#For all chromosome parallel ends
########################################################################################



########################################################################################
#If chr based subs or indels dataframes are too big we can use this version
#For each chromosome sequential starts
def nucleosome_occupancy_analysis_each_chrom_sequential(mutationTypes,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,nucleosomeFilename):

    ##########################################################################
    sample2NumberofSubsDict = getSample2NumberofSubsDict(outputDir,jobname)
    sample2NumberofIndelsDict = getSample2NumberofIndelsDict(outputDir,jobname)
    sample2NumberofDinucsDict = getDictionary(outputDir,jobname,Sample2NumberofDinucsDictFilename)

    subsSignature2NumberofMutationsDict = getSubsSignature2NumberofMutationsDict(outputDir,jobname)
    indelsSignature2NumberofMutationsDict = getIndelsSignature2NumberofMutationsDict(outputDir,jobname)
    dinucsSignature2NumberofMutationsDict = getDictionary(outputDir,jobname,DinucsSignature2NumberofMutationsDictFilename)

    sample2SubsSignature2NumberofMutationsDict = getSample2SubsSignature2NumberofMutationsDict(outputDir,jobname)
    sample2IndelsSignature2NumberofMutationsDict = getSample2IndelsSignature2NumberofMutationsDict(outputDir,jobname)
    sample2DinucsSignature2NumberofMutationsDict =getDictionary(outputDir,jobname,Sample2DinucsSignature2NumberofMutationsDictFilename)
    ##########################################################################

    ##############################################################
    # old way
    # type2AccumulatedSignalArrayDict = {}
    # type2AccumulatedCountArrayDict = {}
    # sample2Type2AccumulatedSignalArrayDict = {}
    # sample2Type2AccumulatedCountArrayDict = {}
    ##############################################################

    ##############################################################
    simNum2Type2AccumulatedSignalArrayDict = {}
    simNum2Type2AccumulatedCountArrayDict = {}
    simNum2Sample2Type2AccumulatedSignalArrayDict = {}
    simNum2Sample2Type2AccumulatedCountArrayDict = {}
    ##############################################################

    ###################################################################################
    ##################  For each chromsome sequential starts ##########################
    ###################################################################################
    for chrLong in chromNamesList:
        chromSize = chromSizesDict[chrLong]

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

            #Read original data
            #READ CHRBASED SINGLE POINT MUTATIONS
            original_chrBased_spms_df=readChrBasedSubsDF(outputDir,jobname,chrLong,SUBS,0)
            if (original_chrBased_spms_df is not None):
                print('chromosome %s  -- chrBased_spms_df: %d in Bytes %f in GigaBytes' %(chrLong,sys.getsizeof(original_chrBased_spms_df),sys.getsizeof(original_chrBased_spms_df)/GIGABYTE_IN_BYTES))

            #READ CHRBASED INDELS
            original_chrBased_indels_df = readChrBasedMutationsDF(outputDir,jobname,chrLong,INDELS,0)
            if (original_chrBased_indels_df is not None):
                print('chromosome %s  -- chrBased_indels_df: %d in Bytes %f in GigaBytes' %(chrLong,sys.getsizeof(original_chrBased_indels_df),sys.getsizeof(original_chrBased_indels_df)/GIGABYTE_IN_BYTES))

            #READ CHRBASED_DINUCS
            original_chrBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong,DINUCS,0)
            if (original_chrBased_dinucs_df is not None):
                print('chromosome %s  -- chrBased_dinucss_df: %d in Bytes %f in GigaBytes' %(chrLong, sys.getsizeof(original_chrBased_dinucs_df), sys.getsizeof(original_chrBased_dinucs_df) / GIGABYTE_IN_BYTES))

            combined_chrBased_subs_df= getCombinedChrBasedDF(outputDir, jobname, chrLong,original_chrBased_spms_df,SUBS,numofSimulations)
            combined_chrBased_indels_df = getCombinedChrBasedDF(outputDir, jobname, chrLong,original_chrBased_indels_df,INDELS,numofSimulations)
            combined_chrBased_dinucs_df= getCombinedChrBasedDF(outputDir, jobname, chrLong,original_chrBased_dinucs_df,DINUCS,numofSimulations)

            ##############################################################
            #old way
            # type2SignalArrayDict = {}
            # type2CountArrayDict = {}
            # sample2Type2SignalArrayDict = {}
            # sample2Type2CountArrayDict = {}
            ##############################################################

            ##############################################################
            simNum2Type2SignalArrayDict = {}
            simNum2Type2CountArrayDict = {}
            simNum2Sample2Type2SignalArrayDict = {}
            simNum2Sample2Type2CountArrayDict = {}
            ##############################################################

            ###############################################################################
            ################### Fill signal and count array starts ########################
            ###############################################################################
            # Fill for single point mutations
            if ((combined_chrBased_subs_df is not None) and (not combined_chrBased_subs_df.empty)):
                combined_chrBased_subs_df.apply(fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated,
                                    nucleosome_array=chrbased_nucleosome_signal_array,
                                    maximum_chrom_size=chromSize,
                                    sample2NumberofMutationsDict=sample2NumberofSubsDict,
                                    signature2NumberofMutationsDict=subsSignature2NumberofMutationsDict,
                                    sample2Signature2NumberofMutationsDict=sample2SubsSignature2NumberofMutationsDict,
                                    simNum2Type2SignalArrayDict=simNum2Type2SignalArrayDict,
                                    simNum2Type2CountArrayDict=simNum2Type2CountArrayDict,
                                    simNum2Sample2Type2SignalArrayDict=simNum2Sample2Type2SignalArrayDict,
                                    simNum2Sample2Type2CountArrayDict=simNum2Sample2Type2CountArrayDict,
                                    MUTATION_SIGNATURE_PROBABILITY_THRESHOLD = SUBSTITUTION_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD,
                                    type = AGGREGATEDSUBSTITUTIONS,
                                    axis=1)

                # former version
                # chrBased_subs_df.apply(fillSignalArrayAndCountArrayForMutations,
                #                     nucleosome_array=chrbased_nucleosome_signal_array,
                #                     maximum_chrom_size=chromSize,
                #                     sample2NumberofMutationsDict=sample2NumberofSubsDict,
                #                     signature2NumberofMutationsDict=subsSignature2NumberofMutationsDict,
                #                     sample2Signature2NumberofMutationsDict=sample2SubsSignature2NumberofMutationsDict,
                #                     type2SignalArrayDict=type2SignalArrayDict,
                #                     type2CountArrayDict=type2CountArrayDict,
                #                     sample2Type2SignalArrayDict=sample2Type2SignalArrayDict,
                #                     sample2Type2CountArrayDict=sample2Type2CountArrayDict,
                #                     MUTATION_SIGNATURE_PROBABILITY_THRESHOLD = SUBSTITUTION_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD,
                #                     type = AGGREGATEDSUBSTITUTIONS,
                #                     axis=1)

            # Fill for indels
            if ((combined_chrBased_indels_df is not None) and (not combined_chrBased_indels_df.empty)):
                combined_chrBased_indels_df.apply(
                    fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated,
                    nucleosome_array=chrbased_nucleosome_signal_array,
                    maximum_chrom_size=chromSize,
                    sample2NumberofMutationsDict=sample2NumberofIndelsDict,
                    signature2NumberofMutationsDict=indelsSignature2NumberofMutationsDict,
                    sample2Signature2NumberofMutationsDict=sample2IndelsSignature2NumberofMutationsDict,
                    simNum2Type2SignalArrayDict=simNum2Type2SignalArrayDict,
                    simNum2Type2CountArrayDict=simNum2Type2CountArrayDict,
                    simNum2Sample2Type2SignalArrayDict=simNum2Sample2Type2SignalArrayDict,
                    simNum2Sample2Type2CountArrayDict=simNum2Sample2Type2CountArrayDict,
                    MUTATION_SIGNATURE_PROBABILITY_THRESHOLD=INDEL_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD,
                    type=AGGREGATEDINDELS,
                    axis=1)

            # Fill for Dinucs
            if ((combined_chrBased_dinucs_df is not None) and (not combined_chrBased_dinucs_df.empty)):
                combined_chrBased_dinucs_df.apply(
                    fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated,
                    nucleosome_array=chrbased_nucleosome_signal_array,
                    maximum_chrom_size=chromSize,
                    sample2NumberofMutationsDict=sample2NumberofDinucsDict,
                    signature2NumberofMutationsDict=dinucsSignature2NumberofMutationsDict,
                    sample2Signature2NumberofMutationsDict=sample2DinucsSignature2NumberofMutationsDict,
                    simNum2Type2SignalArrayDict=simNum2Type2SignalArrayDict,
                    simNum2Type2CountArrayDict=simNum2Type2CountArrayDict,
                    simNum2Sample2Type2SignalArrayDict=simNum2Sample2Type2SignalArrayDict,
                    simNum2Sample2Type2CountArrayDict=simNum2Sample2Type2CountArrayDict,
                    MUTATION_SIGNATURE_PROBABILITY_THRESHOLD=DINUC_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD,
                    type=AGGREGATEDDINUCS,
                    axis=1)
            ###############################################################################
            ################### Fill signal and count array ends ##########################
            ###############################################################################

            ###############################################################################################################
            ###################  Accumulate the results coming from each chromosome starts  ###############################
            ###############################################################################################################

            #####################################################################################################
            ######################### Accumulate right in the left starts  ######################################
            #####################################################################################################
            accumulateSimulationBasedTypeBasedArrays(simNum2Type2AccumulatedSignalArrayDict, simNum2Type2SignalArrayDict)
            accumulateSimulationBasedTypeBasedArrays(simNum2Type2AccumulatedCountArrayDict, simNum2Type2CountArrayDict)

            accumulateSimulationBasedSampleBasedTypeBasedArrays(simNum2Sample2Type2AccumulatedSignalArrayDict,simNum2Sample2Type2SignalArrayDict)
            accumulateSimulationBasedSampleBasedTypeBasedArrays(simNum2Sample2Type2AccumulatedCountArrayDict,simNum2Sample2Type2CountArrayDict)
            #####################################################################################################
            ######################### Accumulate right in the left ends  ########################################
            #####################################################################################################

            ###############################################################################################################
            ###################  Accumulate the results coming from each chromosome ends  #################################
            ###############################################################################################################

    ###################################################################################
    ##################  For each chromsome sequential ends ############################
    ###################################################################################

    writeSimulationBasedAverageNucleosomeOccupancy(simNum2Type2AccumulatedSignalArrayDict,simNum2Type2AccumulatedCountArrayDict,simNum2Sample2Type2AccumulatedSignalArrayDict,simNum2Sample2Type2AccumulatedCountArrayDict,outputDir,jobname)

#For each chromosome sequential ends
########################################################################################