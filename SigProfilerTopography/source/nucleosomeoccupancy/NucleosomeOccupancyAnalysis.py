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
def nucleosomeOccupancyAnalysis(computationType,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,nucleosomeFilename,subs_sig_prob,indels_sig_prob,dinuc_sig_prob):
    print('########################## NucleosomeOccupancyAnalysis starts ##########################')
    if (computationType == COMPUTATION_ALL_CHROMOSOMES_PARALLEL):
        nucleosome_occupancy_analysis_all_chroms_parallel(chromSizesDict, chromNamesList, outputDir,jobname,numofSimulations,nucleosomeFilename)
    elif  (computationType == COMPUTATION_CHROMOSOMES_SEQUENTIAL_ALL_SIMULATIONS_PARALLEL):
        nucleosome_occupancy_analysis_eachChrSequential_allSimsParallel(chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,nucleosomeFilename,subs_sig_prob,indels_sig_prob,dinuc_sig_prob)
    elif (computationType == COMPUTATION_CHROMOSOMES_SEQUENTIAL_SIMULATIONS_SEQUENTIAL):
        nucleosome_occupancy_analysis_eachChrSequential_eachSimSequential(chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,nucleosomeFilename)
    print('########################## NucleosomeOccupancyAnalysis ends ############################')
##############################################################################################################


########################################################################################
def fillSignalArrayAndCountArrays(inputList):
    chrbased_nucleosome_signal_array = inputList[0]
    subs_df = inputList[1]
    indels_df =  inputList[2]
    dinucs_df = inputList[3]
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
    subs_sig_prob = inputList[14]
    indels_sig_prob = inputList[15]
    dinuc_sig_prob = inputList[16]

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
    if ((subs_df is not None) and (not subs_df.empty)):
        subs_df.apply(fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated,
                           nucleosome_array=chrbased_nucleosome_signal_array,
                           maximum_chrom_size =maximum_chrom_size,
                           sample2NumberofMutationsDict = sample2NumberofSubsDict,
                           signature2NumberofMutationsDict = subsSignature2NumberofMutationsDict,
                           sample2Signature2NumberofMutationsDict=sample2SubsSignature2NumberofMutationsDict,
                            simNum2Type2SignalArrayDict=simNum2Type2SignalArrayDict,
                            simNum2Type2CountArrayDict = simNum2Type2CountArrayDict,
                            simNum2Sample2Type2SignalArrayDict = simNum2Sample2Type2SignalArrayDict,
                            simNum2Sample2Type2CountArrayDict = simNum2Sample2Type2CountArrayDict,
                            MUTATION_SIGNATURE_PROBABILITY_THRESHOLD = subs_sig_prob,
                           type = AGGREGATEDSUBSTITUTIONS,
                           axis=1)

    #Fill for indels
    if ((indels_df is not None) and (not indels_df.empty)):
        indels_df.apply(fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated,
                                 nucleosome_array=chrbased_nucleosome_signal_array,
                                 maximum_chrom_size=maximum_chrom_size,
                                 sample2NumberofMutationsDict = sample2NumberofIndelsDict,
                                 signature2NumberofMutationsDict=indelsSignature2NumberofMutationsDict,
                                 sample2Signature2NumberofMutationsDict = sample2IndelsSignature2NumberofMutationsDict,
                                  simNum2Type2SignalArrayDict=simNum2Type2SignalArrayDict,
                                  simNum2Type2CountArrayDict=simNum2Type2CountArrayDict,
                                  simNum2Sample2Type2SignalArrayDict=simNum2Sample2Type2SignalArrayDict,
                                  simNum2Sample2Type2CountArrayDict=simNum2Sample2Type2CountArrayDict,
                                  MUTATION_SIGNATURE_PROBABILITY_THRESHOLD=indels_sig_prob,
                                 type=AGGREGATEDINDELS,
                                 axis=1)

    #Fill for dinucs
    if ((dinucs_df is not None) and (not dinucs_df.empty)):
        dinucs_df.apply(fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated,
            nucleosome_array=chrbased_nucleosome_signal_array,
            maximum_chrom_size=maximum_chrom_size,
            sample2NumberofMutationsDict=sample2NumberofDinucsDict,
            signature2NumberofMutationsDict=dinucsSignature2NumberofMutationsDict,
            sample2Signature2NumberofMutationsDict=sample2DinucsSignature2NumberofMutationsDict,
          simNum2Type2SignalArrayDict=simNum2Type2SignalArrayDict,
          simNum2Type2CountArrayDict=simNum2Type2CountArrayDict,
          simNum2Sample2Type2SignalArrayDict=simNum2Sample2Type2SignalArrayDict,
          simNum2Sample2Type2CountArrayDict=simNum2Sample2Type2CountArrayDict,
          MUTATION_SIGNATURE_PROBABILITY_THRESHOLD=dinuc_sig_prob,
            type=AGGREGATEDDINUCS,
            axis=1)
    ###############################################################################
    ################### Fill signal and count array ends ##########################
    ###############################################################################


    ###############################################################################
    ################### Return  starts ############################################
    ###############################################################################
    # Initialzie the list, you will return this list
    SignalArrayAndCountArray_List = []
    SignalArrayAndCountArray_List.append(simNum2Type2SignalArrayDict)
    SignalArrayAndCountArray_List.append(simNum2Type2CountArrayDict)
    SignalArrayAndCountArray_List.append(simNum2Sample2Type2SignalArrayDict)
    SignalArrayAndCountArray_List.append(simNum2Sample2Type2CountArrayDict)

    return SignalArrayAndCountArray_List
    ###############################################################################
    ################### Return  ends ##############################################
    ###############################################################################
########################################################################################

########################################################################################
#For all chromosome parallel starts
def nucleosome_occupancy_analysis_all_chroms_parallel(chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,nucleosomeFilename):

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

    simNum2Type2AccumulatedSignalArrayDict = {}
    simNum2Type2AccumulatedCountArrayDict = {}
    simNum2Sample2Type2AccumulatedSignalArrayDict = {}
    simNum2Sample2Type2AccumulatedCountArrayDict = {}

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
    # Accumulate the results coming from each chrom
    # allChroms_SignalArrayAndCountArrayList_List is a list of lists.

    for chromosomeBased_SignalArrayAndCountArrayList in allChroms_SignalArrayAndCountArrayList_List:
        simNum2Type2SignalArrayDict = chromosomeBased_SignalArrayAndCountArrayList[0]
        simNum2Type2CountArrayDict = chromosomeBased_SignalArrayAndCountArrayList[1]
        simNum2Sample2Type2SignalArrayDict = chromosomeBased_SignalArrayAndCountArrayList[2]
        simNum2Sample2Type2CountArrayDict = chromosomeBased_SignalArrayAndCountArrayList[3]

        accumulateSimulationBasedTypeBasedArrays(simNum2Type2AccumulatedSignalArrayDict, simNum2Type2SignalArrayDict)
        accumulateSimulationBasedTypeBasedArrays(simNum2Type2AccumulatedCountArrayDict, simNum2Type2CountArrayDict)
        accumulateSimulationBasedSampleBasedTypeBasedArrays(simNum2Sample2Type2AccumulatedSignalArrayDict,simNum2Sample2Type2SignalArrayDict)
        accumulateSimulationBasedSampleBasedTypeBasedArrays(simNum2Sample2Type2AccumulatedCountArrayDict,simNum2Sample2Type2CountArrayDict)

    ###################################################################################
    ##############################  Accumulation ends  ################################
    ###################################################################################

    print('NucleosomeOccupancyAnalysis Results AllChromsParallel starts')
    print('simNum2Type2AccumulatedSignalArrayDict[0]')
    print(simNum2Type2AccumulatedSignalArrayDict[0])
    print('simNum2Type2AccumulatedCountArrayDict[0]')
    print(simNum2Type2AccumulatedCountArrayDict[0])
    print('simNum2Sample2Type2AccumulatedSignalArrayDict[0]')
    print(simNum2Sample2Type2AccumulatedSignalArrayDict[0])
    print('simNumSample2Type2AccumulatedCountArrayDict[0]')
    print(simNum2Sample2Type2AccumulatedCountArrayDict[0])
    print('NucleosomeOccupancyAnalysis Results AllChromsParallel ends')

    writeSimulationBasedAverageNucleosomeOccupancy(simNum2Type2AccumulatedSignalArrayDict,
                                                   simNum2Type2AccumulatedCountArrayDict,
                                                   simNum2Sample2Type2AccumulatedSignalArrayDict,
                                                   simNum2Sample2Type2AccumulatedCountArrayDict, outputDir, jobname)



#For all chromosome parallel ends
########################################################################################


########################################################################################
#For all chromosome parallel starts
def nucleosome_occupancy_analysis_eachChrSequential_allSimsParallel(chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,nucleosomeFilename,subs_sig_prob,indels_sig_prob,dinuc_sig_prob):

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
    ##########################################################################

    simNum2Type2AccumulatedSignalArrayDict = {}
    simNum2Type2AccumulatedCountArrayDict = {}
    simNum2Sample2Type2AccumulatedSignalArrayDict = {}
    simNum2Sample2Type2AccumulatedCountArrayDict = {}

    ###################################################################################
    ##################  For all chromsomes sequential starts ##########################
    ###################################################################################
    for chrLong in chromNamesList:
        poolInputList = []
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

            ###################################################################################
            ##################  For all simulations in parallel starts ########################
            ###################################################################################
            for simNum in range(0,numofSimulations+1):
                inputList = []

                #READ CHRBASED SINGLE POINT MUTATIONS
                chrBased_spms_df=readChrBasedSubsDF(outputDir,jobname,chrLong,SUBS,simNum)
                if (chrBased_spms_df is not None):
                    print('For sim%d chromosome %s  -- chrBased_spms_df: %d in Bytes %f in GigaBytes' %(simNum,chrLong,sys.getsizeof(chrBased_spms_df),sys.getsizeof(chrBased_spms_df)/GIGABYTE_IN_BYTES))

                #READ CHRBASED INDELS
                chrBased_indels_df = readChrBasedMutationsDF(outputDir,jobname,chrLong,INDELS,simNum)
                if (chrBased_indels_df is not None):
                    print('For sim%d chromosome %s  -- chrBased_indels_df: %d in Bytes %f in GigaBytes' %(simNum,chrLong,sys.getsizeof(chrBased_indels_df),sys.getsizeof(chrBased_indels_df)/GIGABYTE_IN_BYTES))

                #READ CHRBASED_DINUCS
                chrBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong,DINUCS,simNum)
                if (chrBased_dinucs_df is not None):
                    print('For sim%d chromosome %s  -- chrBased_dinucss_df: %d in Bytes %f in GigaBytes' %(simNum,chrLong,sys.getsizeof(chrBased_dinucs_df), sys.getsizeof(chrBased_dinucs_df)/GIGABYTE_IN_BYTES))

                inputList.append(chrbased_nucleosome_signal_array)
                inputList.append(chrBased_spms_df)
                inputList.append(chrBased_indels_df)
                inputList.append(chrBased_dinucs_df)
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
                inputList.append(subs_sig_prob)
                inputList.append(indels_sig_prob)
                inputList.append(dinuc_sig_prob)
                poolInputList.append(inputList)
            ###################################################################################
            ##################  For all simulations in parallel ends ##########################
            ###################################################################################

            ########################################################################
            #Provides list of  lists (each sub list comes from a simulation)
            allSimulations_SignalArrayAndCountArrayList_List = pool.map(fillSignalArrayAndCountArrays,poolInputList)
            ########################################################################

            #####################################################################################################
            ######################### Accumulate right in the left starts  ######################################
            #####################################################################################################
            for simulatonBased_SignalArrayAndCountArrayList in  allSimulations_SignalArrayAndCountArrayList_List:
                simNum2Type2SignalArrayDict = simulatonBased_SignalArrayAndCountArrayList[0]
                simNum2Type2CountArrayDict = simulatonBased_SignalArrayAndCountArrayList[1]
                simNum2Sample2Type2SignalArrayDict = simulatonBased_SignalArrayAndCountArrayList[2]
                simNum2Sample2Type2CountArrayDict = simulatonBased_SignalArrayAndCountArrayList[3]

                accumulateSimulationBasedTypeBasedArrays(simNum2Type2AccumulatedSignalArrayDict,simNum2Type2SignalArrayDict)
                accumulateSimulationBasedTypeBasedArrays(simNum2Type2AccumulatedCountArrayDict, simNum2Type2CountArrayDict)
                accumulateSimulationBasedSampleBasedTypeBasedArrays(simNum2Sample2Type2AccumulatedSignalArrayDict,simNum2Sample2Type2SignalArrayDict)
                accumulateSimulationBasedSampleBasedTypeBasedArrays(simNum2Sample2Type2AccumulatedCountArrayDict,simNum2Sample2Type2CountArrayDict)

            #####################################################################################################
            ######################### Accumulate right in the left ends  ########################################
            #####################################################################################################

        ###################################################################################
        ##################  For all chromsomes sequential ends ############################
        ###################################################################################

    ################################
    pool.close()
    pool.join()
    ################################

    print('NucleosomeOccupancyAnalysis Results EachChrSequential_AllSimsParallel starts')
    print('simNum2Type2AccumulatedSignalArrayDict[0]')
    print(simNum2Type2AccumulatedSignalArrayDict[0])
    print('simNum2Type2AccumulatedCountArrayDict[0]')
    print(simNum2Type2AccumulatedCountArrayDict[0])
    print('simNum2Sample2Type2AccumulatedSignalArrayDict[0]')
    print(simNum2Sample2Type2AccumulatedSignalArrayDict[0])
    print('simNumSample2Type2AccumulatedCountArrayDict[0]')
    print(simNum2Sample2Type2AccumulatedCountArrayDict[0])
    print('NucleosomeOccupancyAnalysis Results EachChrSequential_AllSimsParallel ends')

    writeSimulationBasedAverageNucleosomeOccupancy(simNum2Type2AccumulatedSignalArrayDict,
                                                   simNum2Type2AccumulatedCountArrayDict,
                                                   simNum2Sample2Type2AccumulatedSignalArrayDict,
                                                   simNum2Sample2Type2AccumulatedCountArrayDict, outputDir, jobname)


#For all chromosome parallel ends
########################################################################################



########################################################################################
#If chr based subs or indels dataframes are too big we can use this version
#For each chromosome sequential starts
def nucleosome_occupancy_analysis_eachChrSequential_eachSimSequential(chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,nucleosomeFilename):

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

            ###################################################################################
            ##################  For each simulation sequential starts #########################
            ###################################################################################
            for simNum in range(0,numofSimulations+1):
                # READ CHRBASED SINGLE POINT MUTATIONS
                chrBased_spms_df = readChrBasedSubsDF(outputDir, jobname, chrLong, SUBS,simNum)
                if (chrBased_spms_df is not None):
                    print('for sim%d chromosome %s  -- chrBased_spms_df: %d in Bytes %f in GigaBytes' %(simNum,chrLong,sys.getsizeof(chrBased_spms_df),sys.getsizeof(chrBased_spms_df)/GIGABYTE_IN_BYTES))

                #READ CHRBASED INDELS
                chrBased_indels_df = readChrBasedMutationsDF(outputDir,jobname,chrLong,INDELS,simNum)
                if (chrBased_indels_df is not None):
                    print('for sim%d chromosome %s  -- chrBased_indels_df: %d in Bytes %f in GigaBytes' %(simNum,chrLong,sys.getsizeof(chrBased_indels_df),sys.getsizeof(chrBased_indels_df)/GIGABYTE_IN_BYTES))

                #READ CHRBASED_DINUCS
                chrBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong,DINUCS,simNum)
                if (chrBased_dinucs_df is not None):
                    print('for sim%d chromosome %s  -- chrBased_dinucss_df: %d in Bytes %f in GigaBytes' %(simNum,chrLong,sys.getsizeof(chrBased_dinucs_df), sys.getsizeof(chrBased_dinucs_df) / GIGABYTE_IN_BYTES))

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
                if ((chrBased_spms_df is not None) and (not chrBased_spms_df.empty)):
                    chrBased_spms_df.apply(fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated,
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


                # Fill for indels
                if ((chrBased_indels_df is not None) and (not chrBased_indels_df.empty)):
                    chrBased_indels_df.apply(
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
                if ((chrBased_dinucs_df is not None) and (not chrBased_dinucs_df.empty)):
                    chrBased_dinucs_df.apply(
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
                ################  Accumulate the results coming from each simulation and  each chromosome starts  #############
                ###############################################################################################################
                ######################### Accumulate right in the left starts  ######################################
                accumulateSimulationBasedTypeBasedArrays(simNum2Type2AccumulatedSignalArrayDict, simNum2Type2SignalArrayDict)
                accumulateSimulationBasedTypeBasedArrays(simNum2Type2AccumulatedCountArrayDict, simNum2Type2CountArrayDict)
                accumulateSimulationBasedSampleBasedTypeBasedArrays(simNum2Sample2Type2AccumulatedSignalArrayDict,simNum2Sample2Type2SignalArrayDict)
                accumulateSimulationBasedSampleBasedTypeBasedArrays(simNum2Sample2Type2AccumulatedCountArrayDict,simNum2Sample2Type2CountArrayDict)
                ###############################################################################################################
                ################  Accumulate the results coming from each simulation and  each chromosome ends  ###############
                ###############################################################################################################

            ###################################################################################
            ##################  For each simulation sequential starts #########################
            ###################################################################################

    ###################################################################################
    ##################  For each chromsome sequential ends ############################
    ###################################################################################

    print('NucleosomeOccupancyAnalysis Results EachChrSequential_EachSimSequential starts')
    print('simNum2Type2AccumulatedSignalArrayDict[0]')
    print(simNum2Type2AccumulatedSignalArrayDict[0])
    print('simNum2Type2AccumulatedCountArrayDict[0]')
    print(simNum2Type2AccumulatedCountArrayDict[0])
    print('simNum2Sample2Type2AccumulatedSignalArrayDict[0]')
    print(simNum2Sample2Type2AccumulatedSignalArrayDict[0])
    print('simNumSample2Type2AccumulatedCountArrayDict[0]')
    print(simNum2Sample2Type2AccumulatedCountArrayDict[0])
    print('NucleosomeOccupancyAnalysis Results EachChrSequential_EachSimSequential ends')

    writeSimulationBasedAverageNucleosomeOccupancy(simNum2Type2AccumulatedSignalArrayDict,simNum2Type2AccumulatedCountArrayDict,simNum2Sample2Type2AccumulatedSignalArrayDict,simNum2Sample2Type2AccumulatedCountArrayDict,outputDir,jobname)

#For each chromosome sequential ends
########################################################################################