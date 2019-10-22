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
import pyBigWig


##############################################################################################################
#main function
def occupancyAnalysis(computationType,occupancy_type,sample_based,plusorMinus,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,library_file_with_path,library_file_memo,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict):
    print('\n#################################################################################')
    print('--- %s Analysis starts' %(occupancy_type))

    print('--- Computation Type:%s' % (computationType))
    print('--- Occupancy Type:%s' % (occupancy_type))
    print('--- Library file with path: %s\n' %library_file_with_path)

    if (computationType == COMPUTATION_ALL_CHROMOSOMES_PARALLEL):
        #depreceated, not maintained
        nucleosome_occupancy_analysis_all_chroms_parallel(chromSizesDict, chromNamesList, outputDir,jobname,numofSimulations,library_file_with_path,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict)
    elif  (computationType == COMPUTATION_CHROMOSOMES_SEQUENTIAL_ALL_SIMULATIONS_PARALLEL):

        #Using offline prepared npy files
        # occupancy_analysis_eachChrSequential_allSimsParallel(chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,libraryFilename_with_dir,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,occupancy_type)

        #Using pyBigWig for BigWig and BigBed files
        #Using Bed files preparing chr based signal array online
        occupancy_analysis_eachChrSequential_allSimsParallel_using_pyBigWig(occupancy_type,sample_based,plusorMinus,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,library_file_with_path,library_file_memo,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict)

    elif (computationType == COMPUTATION_CHROMOSOMES_SEQUENTIAL_SIMULATIONS_SEQUENTIAL):
        #depreceated, not maintained
        nucleosome_occupancy_analysis_eachChrSequential_eachSimSequential(chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,library_file_with_path)
    print('--- %s Analysis ends' %(occupancy_type))
    print('#################################################################################\n')
##############################################################################################################


########################################################################################
#Using pyBigWig for bigBed and bigWig files starts
#Using bed files prepared on the fly starts
def fillSignalArrayAndCountArrays_using_pyBigWig(inputList):
    chrLong = inputList[0]
    signal_file_with_path=inputList[1]
    library_file_type=inputList[2]
    chrBasedSignalArray=inputList[3]
    subs_df = inputList[4]
    indels_df =  inputList[5]
    dinucs_df = inputList[6]
    maximum_chrom_size =  inputList[7]
    sample2NumberofSubsDict = inputList[8]
    sample2NumberofIndelsDict = inputList[9]
    sample2NumberofDinucsDict = inputList[10]
    sample2SubsSignature2NumberofMutationsDict =  inputList[11]
    sample2IndelsSignature2NumberofMutationsDict = inputList[12]
    sample2DinucsSignature2NumberofMutationsDict = inputList[13]
    subsSignature2PropertiesListDict = inputList[14]
    indelsSignature2PropertiesListDict = inputList[15]
    dinucsSignature2PropertiesListDict = inputList[16]
    plusorMinus=inputList[17]

    # print('Debug We are in fillSignalArrayAndCountArrays_using_pyBigWig for chrLong:%s' %(chrLong))

    ##############################################################
    simNum2Type2SignalArrayDict = {}
    simNum2Type2CountArrayDict = {}
    simNum2Sample2Type2SignalArrayDict = {}
    simNum2Sample2Type2CountArrayDict = {}
    ##############################################################

    library_file=None
    my_upperBound=None
    signal_index=None
    if (library_file_type==BIGWIG):
        try:
            library_file = pyBigWig.open(signal_file_with_path)
            if chrLong in library_file.chroms():
                maximum_chrom_size = library_file.chroms()[chrLong]
            # For BigWig Files information in header is correct
            if ('sumData' in library_file.header()) and ('nBasesCovered' in library_file.header()):
                my_mean = library_file.header()['sumData'] / library_file.header()['nBasesCovered']
                std_dev=(library_file.header()['sumSquared'] - 2 * my_mean * library_file.header()['sumData'] + library_file.header()['nBasesCovered'] * my_mean * my_mean) ** (0.5) / (library_file.header()['nBasesCovered'] ** (0.5))
                # Scientific definition of outlier
                my_upperBound = std_dev * 3
            else:
                # Undefined
                my_upperBound = np.iinfo(np.int16).max
        except:
            print(signal_file_with_path)

    elif (library_file_type==BIGBED):
        try:
            library_file = pyBigWig.open(signal_file_with_path)
            if BED_6PLUS4 in str(library_file.SQL()):
                signal_index=3
            elif BED_9PLUS2 in str(library_file.SQL()):
                signal_index=7
            if chrLong in library_file.chroms():
                # For BigBed Files information in header is not meaningful
                maximum_chrom_size = library_file.chroms()[chrLong]
                my_mean=np.mean([float(entry[2].split('\t')[signal_index]) for entry in library_file.entries(chrLong, 0, maximum_chrom_size)])
                #Not scientific definition of outlier
                my_upperBound = my_mean * 10
            else:
                # Undefined
                my_upperBound = np.iinfo(np.int16).max
        except:
            print(signal_file_with_path)

    # For information
    # print('signal_file_with_path: %s' %(signal_file_with_path))
    # print('library_file_type: %s' %(library_file_type))
    # print('For %s maximum_chrom_size' %(chrLong))
    # print(maximum_chrom_size)
    # print('my_upperBound')
    # print(my_upperBound)

    ###################################################################################
    if ((((library_file_type==BIGWIG) or (library_file_type==BIGBED)) and (library_file is not None) and (chrLong in library_file.chroms()))
       or
        ((library_file_type==BED) and (chrBasedSignalArray is not None)) ):
        ######################################################## #######################
        ################### Fill signal and count array starts ########################
        ###############################################################################
        #Fill for single point mutations
        if ((subs_df is not None) and (not subs_df.empty)):
            # subs_df.apply(fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated_using_pyBigWig,
            # chrLong=chrLong,
            # library_file=library_file,
            # chrBasedSignalArray=chrBasedSignalArray,
            # library_file_type=library_file_type,
            # my_upperBound=my_upperBound,
            # maximum_chrom_size =maximum_chrom_size,
            # sample2NumberofMutationsDict = sample2NumberofSubsDict,
            # sample2Signature2NumberofMutationsDict=sample2SubsSignature2NumberofMutationsDict,
            # simNum2Type2SignalArrayDict=simNum2Type2SignalArrayDict,
            # simNum2Type2CountArrayDict = simNum2Type2CountArrayDict,
            # simNum2Sample2Type2SignalArrayDict = simNum2Sample2Type2SignalArrayDict,
            # simNum2Sample2Type2CountArrayDict = simNum2Sample2Type2CountArrayDict,
            # signature2PropertiesListDict=subsSignature2PropertiesListDict,
            # my_type = AGGREGATEDSUBSTITUTIONS,
            # plusorMinus=plusorMinus,
            # axis=1)

            columnNamesList = list(subs_df.columns.values)
            # We assume that after the column named 'Mutation' there are the signature columns in tab separated way.
            # Last column is the simulation number
            mutationtIndex = columnNamesList.index(MUTATION)
            signatures = columnNamesList[(mutationtIndex + 1):-1]
            mycolumns=[SAMPLE,START,SIMULATION_NUMBER]
            mycolumns.extend(signatures)
            [fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated_using_pyBigWig_using_list_comp(
                row,
                chrLong,
                library_file,
                chrBasedSignalArray,
                library_file_type,
                signal_index,
                my_upperBound,
                maximum_chrom_size,
                sample2NumberofSubsDict,
                sample2SubsSignature2NumberofMutationsDict,
                simNum2Type2SignalArrayDict,
                simNum2Type2CountArrayDict,
                simNum2Sample2Type2SignalArrayDict,
                simNum2Sample2Type2CountArrayDict,
                subsSignature2PropertiesListDict,
                AGGREGATEDSUBSTITUTIONS,
                plusorMinus,
                mycolumns) for row in subs_df[mycolumns].values]

        #Fill for indels
        if ((indels_df is not None) and (not indels_df.empty)):
            # indels_df.apply(fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated_using_pyBigWig,
            # chrLong=chrLong,
            # library_file=library_file,
            # chrBasedSignalArray=chrBasedSignalArray,
            # library_file_type=library_file_type,
            # my_upperBound=my_upperBound,
            # maximum_chrom_size=maximum_chrom_size,
            # sample2NumberofMutationsDict = sample2NumberofIndelsDict,
            # sample2Signature2NumberofMutationsDict = sample2IndelsSignature2NumberofMutationsDict,
            # simNum2Type2SignalArrayDict=simNum2Type2SignalArrayDict,
            # simNum2Type2CountArrayDict=simNum2Type2CountArrayDict,
            # simNum2Sample2Type2SignalArrayDict=simNum2Sample2Type2SignalArrayDict,
            # simNum2Sample2Type2CountArrayDict=simNum2Sample2Type2CountArrayDict,
            # signature2PropertiesListDict=indelsSignature2PropertiesListDict,
            # my_type=AGGREGATEDINDELS,
            # plusorMinus=plusorMinus,
            # axis=1)

            columnNamesList = list(indels_df.columns.values)
            # We assume that after the column named 'Mutation' there are the signature columns in tab separated way.
            # Last column is the simulation number
            mutationtIndex = columnNamesList.index(MUTATION)
            signatures = columnNamesList[(mutationtIndex + 1):-1]
            mycolumns=[SAMPLE,START,SIMULATION_NUMBER]
            mycolumns.extend(signatures)
            [fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated_using_pyBigWig_using_list_comp(
                row,
                chrLong,
                library_file,
                chrBasedSignalArray,
                library_file_type,
                signal_index,
                my_upperBound,
                maximum_chrom_size,
                sample2NumberofIndelsDict,
                sample2IndelsSignature2NumberofMutationsDict,
                simNum2Type2SignalArrayDict,
                simNum2Type2CountArrayDict,
                simNum2Sample2Type2SignalArrayDict,
                simNum2Sample2Type2CountArrayDict,
                indelsSignature2PropertiesListDict,
                AGGREGATEDINDELS,
                plusorMinus,
                mycolumns) for row in indels_df[mycolumns].values]

        #Fill for dinucs
        if ((dinucs_df is not None) and (not dinucs_df.empty)):
            # dinucs_df.apply(fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated_using_pyBigWig,
            # chrLong=chrLong,
            # library_file=library_file,
            # chrBasedSignalArray=chrBasedSignalArray,
            # library_file_type=library_file_type,
            # my_upperBound=my_upperBound,
            # maximum_chrom_size=maximum_chrom_size,
            # sample2NumberofMutationsDict=sample2NumberofDinucsDict,
            # sample2Signature2NumberofMutationsDict=sample2DinucsSignature2NumberofMutationsDict,
            # simNum2Type2SignalArrayDict=simNum2Type2SignalArrayDict,
            # simNum2Type2CountArrayDict=simNum2Type2CountArrayDict,
            # simNum2Sample2Type2SignalArrayDict=simNum2Sample2Type2SignalArrayDict,
            # simNum2Sample2Type2CountArrayDict=simNum2Sample2Type2CountArrayDict,
            # signature2PropertiesListDict=dinucsSignature2PropertiesListDict,
            # my_type=AGGREGATEDDINUCS,
            # plusorMinus=plusorMinus,
            # axis=1)

            columnNamesList = list(dinucs_df.columns.values)
            # We assume that after the column named 'Mutation' there are the signature columns in tab separated way.
            # Last column is the simulation number
            mutationtIndex = columnNamesList.index(MUTATION)
            signatures = columnNamesList[(mutationtIndex + 1):-1]
            mycolumns=[SAMPLE,START,SIMULATION_NUMBER]
            mycolumns.extend(signatures)
            [fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated_using_pyBigWig_using_list_comp(
                row,
                chrLong,
                library_file,
                chrBasedSignalArray,
                library_file_type,
                signal_index,
                my_upperBound,
                maximum_chrom_size,
                sample2NumberofDinucsDict,
                sample2DinucsSignature2NumberofMutationsDict,
                simNum2Type2SignalArrayDict,
                simNum2Type2CountArrayDict,
                simNum2Sample2Type2SignalArrayDict,
                simNum2Sample2Type2CountArrayDict,
                dinucsSignature2PropertiesListDict,
                AGGREGATEDDINUCS,
                plusorMinus,
                mycolumns) for row in dinucs_df[mycolumns].values]
            ###############################################################################
            ################### Fill signal and count array ends ##########################
            ###############################################################################

    if (library_file is not None):
        library_file.close()

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

#Using pyBigWig for bigBed and bigWig files ends
#Using bed files prepared on the fly ends
########################################################################################


########################################################################################
def fillSignalArrayAndCountArrays_using_pyBigWig_using_apply_async(chrLong,signal_file_with_path,library_file_type,chrBasedSignalArray,
                                                                   subs_df,indels_df,dinucs_df,maximum_chrom_size,
                                                                   sample2NumberofSubsDict,sample2NumberofIndelsDict,sample2NumberofDinucsDict,
                                                                   sample2SubsSignature2NumberofMutationsDict,
                                                                   sample2IndelsSignature2NumberofMutationsDict,
                                                                   sample2DinucsSignature2NumberofMutationsDict,
                                                                   subsSignature2PropertiesListDict,
                                                                   indelsSignature2PropertiesListDict,
                                                                   dinucsSignature2PropertiesListDict,
                                                                   plusorMinus):
    ##############################################################
    simNum2Type2SignalArrayDict = {}
    simNum2Type2CountArrayDict = {}
    simNum2Sample2Type2SignalArrayDict = {}
    simNum2Sample2Type2CountArrayDict = {}
    ##############################################################

    library_file = None
    my_upperBound = None
    signal_index = None
    if (library_file_type == BIGWIG):
        try:
            library_file = pyBigWig.open(signal_file_with_path)
            if chrLong in library_file.chroms():
                maximum_chrom_size = library_file.chroms()[chrLong]
            # For BigWig Files information in header is correct
            if ('sumData' in library_file.header()) and ('nBasesCovered' in library_file.header()):
                my_mean = library_file.header()['sumData'] / library_file.header()['nBasesCovered']
                std_dev = (library_file.header()['sumSquared'] - 2 * my_mean * library_file.header()['sumData'] +
                           library_file.header()['nBasesCovered'] * my_mean * my_mean) ** (0.5) / (
                                      library_file.header()['nBasesCovered'] ** (0.5))
                # Scientific definition of outlier
                my_upperBound = std_dev * 3
            else:
                # Undefined
                my_upperBound = np.iinfo(np.int16).max
        except:
            print(signal_file_with_path)

    elif (library_file_type == BIGBED):
        try:
            library_file = pyBigWig.open(signal_file_with_path)
            if BED_6PLUS4 in str(library_file.SQL()):
                signal_index = 3
            elif BED_9PLUS2 in str(library_file.SQL()):
                signal_index = 7
            if chrLong in library_file.chroms():
                # For BigBed Files information in header is not meaningful
                maximum_chrom_size = library_file.chroms()[chrLong]
                my_mean = np.mean([float(entry[2].split('\t')[signal_index]) for entry in
                                   library_file.entries(chrLong, 0, maximum_chrom_size)])
                # Not scientific definition of outlier
                my_upperBound = my_mean * 10
            else:
                # Undefined
                my_upperBound = np.iinfo(np.int16).max
        except:
            print(signal_file_with_path)

    ###################################################################################
    if ((((library_file_type == BIGWIG) or (library_file_type == BIGBED)) and (library_file is not None) and (
            chrLong in library_file.chroms()))
            or
            ((library_file_type == BED) and (chrBasedSignalArray is not None))):
        ######################################################## #######################
        ################### Fill signal and count array starts ########################
        ###############################################################################
        # Fill for single point mutations
        if ((subs_df is not None) and (not subs_df.empty)):
            columnNamesList = list(subs_df.columns.values)
            # We assume that after the column named 'Mutation' there are the signature columns in tab separated way.
            # Last column is the simulation number
            mutationtIndex = columnNamesList.index(MUTATION)
            signatures = columnNamesList[(mutationtIndex + 1):-1]
            mycolumns = [SAMPLE, START, SIMULATION_NUMBER]
            mycolumns.extend(signatures)
            [fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated_using_pyBigWig_using_list_comp(
                row,
                chrLong,
                library_file,
                chrBasedSignalArray,
                library_file_type,
                signal_index,
                my_upperBound,
                maximum_chrom_size,
                sample2NumberofSubsDict,
                sample2SubsSignature2NumberofMutationsDict,
                simNum2Type2SignalArrayDict,
                simNum2Type2CountArrayDict,
                simNum2Sample2Type2SignalArrayDict,
                simNum2Sample2Type2CountArrayDict,
                subsSignature2PropertiesListDict,
                AGGREGATEDSUBSTITUTIONS,
                plusorMinus,
                mycolumns) for row in subs_df[mycolumns].values]

        # Fill for indels
        if ((indels_df is not None) and (not indels_df.empty)):
            columnNamesList = list(indels_df.columns.values)
            # We assume that after the column named 'Mutation' there are the signature columns in tab separated way.
            # Last column is the simulation number
            mutationtIndex = columnNamesList.index(MUTATION)
            signatures = columnNamesList[(mutationtIndex + 1):-1]
            mycolumns = [SAMPLE, START, SIMULATION_NUMBER]
            mycolumns.extend(signatures)
            [fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated_using_pyBigWig_using_list_comp(
                row,
                chrLong,
                library_file,
                chrBasedSignalArray,
                library_file_type,
                signal_index,
                my_upperBound,
                maximum_chrom_size,
                sample2NumberofIndelsDict,
                sample2IndelsSignature2NumberofMutationsDict,
                simNum2Type2SignalArrayDict,
                simNum2Type2CountArrayDict,
                simNum2Sample2Type2SignalArrayDict,
                simNum2Sample2Type2CountArrayDict,
                indelsSignature2PropertiesListDict,
                AGGREGATEDINDELS,
                plusorMinus,
                mycolumns) for row in indels_df[mycolumns].values]

        # Fill for dinucs
        if ((dinucs_df is not None) and (not dinucs_df.empty)):
            columnNamesList = list(dinucs_df.columns.values)
            # We assume that after the column named 'Mutation' there are the signature columns in tab separated way.
            # Last column is the simulation number
            mutationtIndex = columnNamesList.index(MUTATION)
            signatures = columnNamesList[(mutationtIndex + 1):-1]
            mycolumns = [SAMPLE, START, SIMULATION_NUMBER]
            mycolumns.extend(signatures)
            [fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated_using_pyBigWig_using_list_comp(
                row,
                chrLong,
                library_file,
                chrBasedSignalArray,
                library_file_type,
                signal_index,
                my_upperBound,
                maximum_chrom_size,
                sample2NumberofDinucsDict,
                sample2DinucsSignature2NumberofMutationsDict,
                simNum2Type2SignalArrayDict,
                simNum2Type2CountArrayDict,
                simNum2Sample2Type2SignalArrayDict,
                simNum2Sample2Type2CountArrayDict,
                dinucsSignature2PropertiesListDict,
                AGGREGATEDDINUCS,
                plusorMinus,
                mycolumns) for row in dinucs_df[mycolumns].values]
            ###############################################################################
            ################### Fill signal and count array ends ##########################
            ###############################################################################

    if (library_file is not None):
        library_file.close()

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


# Using pyBigWig for bigBed and bigWig files ends
# Using bed files prepared on the fly ends
########################################################################################

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
    sample2SubsSignature2NumberofMutationsDict =  inputList[8]
    sample2IndelsSignature2NumberofMutationsDict = inputList[9]
    sample2DinucsSignature2NumberofMutationsDict = inputList[10]
    subsSignature2PropertiesListDict = inputList[11]
    indelsSignature2PropertiesListDict = inputList[12]
    dinucsSignature2PropertiesListDict = inputList[13]

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
                           sample2Signature2NumberofMutationsDict=sample2SubsSignature2NumberofMutationsDict,
                            simNum2Type2SignalArrayDict=simNum2Type2SignalArrayDict,
                            simNum2Type2CountArrayDict = simNum2Type2CountArrayDict,
                            simNum2Sample2Type2SignalArrayDict = simNum2Sample2Type2SignalArrayDict,
                            simNum2Sample2Type2CountArrayDict = simNum2Sample2Type2CountArrayDict,
                            signature2PropertiesListDict=subsSignature2PropertiesListDict,
                           type = AGGREGATEDSUBSTITUTIONS,
                           axis=1)

    #Fill for indels
    if ((indels_df is not None) and (not indels_df.empty)):
        indels_df.apply(fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated,
                                 nucleosome_array=chrbased_nucleosome_signal_array,
                                 maximum_chrom_size=maximum_chrom_size,
                                 sample2NumberofMutationsDict = sample2NumberofIndelsDict,
                                 sample2Signature2NumberofMutationsDict = sample2IndelsSignature2NumberofMutationsDict,
                                simNum2Type2SignalArrayDict=simNum2Type2SignalArrayDict,
                                simNum2Type2CountArrayDict=simNum2Type2CountArrayDict,
                                simNum2Sample2Type2SignalArrayDict=simNum2Sample2Type2SignalArrayDict,
                                simNum2Sample2Type2CountArrayDict=simNum2Sample2Type2CountArrayDict,
                                signature2PropertiesListDict=indelsSignature2PropertiesListDict,
                                type=AGGREGATEDINDELS,
                                 axis=1)

    #Fill for dinucs
    if ((dinucs_df is not None) and (not dinucs_df.empty)):
        dinucs_df.apply(fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated,
            nucleosome_array=chrbased_nucleosome_signal_array,
            maximum_chrom_size=maximum_chrom_size,
            sample2NumberofMutationsDict=sample2NumberofDinucsDict,
            sample2Signature2NumberofMutationsDict=sample2DinucsSignature2NumberofMutationsDict,
          simNum2Type2SignalArrayDict=simNum2Type2SignalArrayDict,
          simNum2Type2CountArrayDict=simNum2Type2CountArrayDict,
          simNum2Sample2Type2SignalArrayDict=simNum2Sample2Type2SignalArrayDict,
          simNum2Sample2Type2CountArrayDict=simNum2Sample2Type2CountArrayDict,
        signature2PropertiesListDict=dinucsSignature2PropertiesListDict,
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
def nucleosome_occupancy_analysis_all_chroms_parallel(chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,nucleosomeFilename,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict):

    ##########################################################################
    sample2NumberofSubsDict = getSample2NumberofSubsDict(outputDir,jobname)
    sample2NumberofIndelsDict = getSample2NumberofIndelsDict(outputDir,jobname)
    sample2NumberofDinucsDict = getDictionary(outputDir,jobname, Sample2NumberofDinucsDictFilename)

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
            inputList.append(sample2SubsSignature2NumberofMutationsDict)
            inputList.append(sample2IndelsSignature2NumberofMutationsDict)
            inputList.append(sample2DinucsSignature2NumberofMutationsDict)
            inputList.append(subsSignature2PropertiesListDict)
            inputList.append(indelsSignature2PropertiesListDict)
            inputList.append(dinucsSignature2PropertiesListDict)

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
class Accumulator:
    def __init__(self):
        # self.signal = np.zeros((10000,), dtype=float)
        # self.count = np.zeros((10000,), dtype=int)
        self.simNum2Type2AccumulatedSignalArrayDict = {}
        self.simNum2Type2AccumulatedCountArrayDict = {}
        self.simNum2Sample2Type2AccumulatedSignalArrayDict = {}
        self.simNum2Sample2Type2AccumulatedCountArrayDict = {}

    def on_result(self, result):
        # self.signal += result[0]
        # self.count += result[1]
        simNum2Type2SignalArrayDict = result[0]
        simNum2Type2CountArrayDict = result[1]
        simNum2Sample2Type2SignalArrayDict = result[2]
        simNum2Sample2Type2CountArrayDict = result[3]

        accumulateSimulationBasedTypeBasedArrays(self.simNum2Type2AccumulatedSignalArrayDict, simNum2Type2SignalArrayDict)
        accumulateSimulationBasedTypeBasedArrays(self.simNum2Type2AccumulatedCountArrayDict, simNum2Type2CountArrayDict)
        accumulateSimulationBasedSampleBasedTypeBasedArrays(self.simNum2Sample2Type2AccumulatedSignalArrayDict,simNum2Sample2Type2SignalArrayDict)
        accumulateSimulationBasedSampleBasedTypeBasedArrays(self.simNum2Sample2Type2AccumulatedCountArrayDict,simNum2Sample2Type2CountArrayDict)
########################################################################################



########################################################################################
def fillInputList(outputDir, jobname,chrLong,simNum,chromSizesDict,library_file_with_path,library_file_type,library_file_df,library_file_df_grouped,
                  sample2NumberofSubsDict,sample2NumberofIndelsDict,sample2NumberofDinucsDict,
                  sample2SubsSignature2NumberofMutationsDict,sample2IndelsSignature2NumberofMutationsDict,sample2DinucsSignature2NumberofMutationsDict,
                  subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,plusorMinus):

    # print('Debug We are in fillInputList: chrLong:%s simNum:%d' %(chrLong,simNum))

    inputList=[]
    chromSize = chromSizesDict[chrLong]
    chrBasedSignalArray = None

    #################################################################################################################
    # if (os.path.exists(chrBasedSignalNucleosmeFile) and os.path.exists(chrBasedCountNucleosmeFile)):
    if (os.path.exists(library_file_with_path)):
        # TODO: This is specific to our data right now
        # Nucleosomes have chrM
        # SinglePointMutations and Indels have chrMT
        chrLong_for_mutations_data = chrLong
        if (chrLong == 'chrM'):
            chrLong_for_mutations_data = 'chrMT'

        if (((library_file_type == BED) or (library_file_type == NARROWPEAK)) and (library_file_df is not None)):
            chrom_based_library_df = library_file_df_grouped.get_group(chrLong)
            # chrBasedSignalArray and library_file_df  signal column is of type np.float32
            chrBasedSignalArray = np.zeros(chromSize, dtype=np.float32)
            # TODO Can we fill chrBasedSignalArray faster?
            # chrom_based_library_df.apply(updateChrBasedSignalArray, chrBasedSignalArray=chrBasedSignalArray, axis=1)
            [fillNumpyArray(start, end, signal, chrBasedSignalArray) for start, end, signal in
             zip(chrom_based_library_df['start'], chrom_based_library_df['end'], chrom_based_library_df['signal'])]

        chrBased_spms_df = readChrBasedSubsDF(outputDir, jobname, chrLong_for_mutations_data, SUBS, simNum)
        chrBased_indels_df = readChrBasedMutationsDF(outputDir, jobname, chrLong_for_mutations_data, INDELS, simNum)
        chrBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong_for_mutations_data, DINUCS, simNum)

        inputList.append(chrLong)
        inputList.append(library_file_with_path)
        inputList.append(library_file_type)
        inputList.append(chrBasedSignalArray)
        inputList.append(chrBased_spms_df)
        inputList.append(chrBased_indels_df)
        inputList.append(chrBased_dinucs_df)
        inputList.append(chromSize)
        inputList.append(sample2NumberofSubsDict)
        inputList.append(sample2NumberofIndelsDict)
        inputList.append(sample2NumberofDinucsDict)
        inputList.append(sample2SubsSignature2NumberofMutationsDict)
        inputList.append(sample2IndelsSignature2NumberofMutationsDict)
        inputList.append(sample2DinucsSignature2NumberofMutationsDict)
        inputList.append(subsSignature2PropertiesListDict)
        inputList.append(indelsSignature2PropertiesListDict)
        inputList.append(dinucsSignature2PropertiesListDict)
        inputList.append(plusorMinus)

    #################################################################################################################
    return inputList
########################################################################################

########################################################################################
#Using pyBigWig for bigBed and bigWig files starts
#Using bed files prepared on the fly starts
def occupancy_analysis_eachChrSequential_allSimsParallel_using_pyBigWig(occupancy_type,
                                                                        sample_based,
                                                                        plusorMinus,
                                                                        chromSizesDict,
                                                                        chromNamesList,
                                                                        outputDir,
                                                                        jobname,
                                                                        numofSimulations,
                                                                        library_file_with_path,
                                                                        library_file_memo,
                                                                        subsSignature2PropertiesListDict,
                                                                        indelsSignature2PropertiesListDict,
                                                                        dinucsSignature2PropertiesListDict):

    if sample_based:
        ##########################################################################
        sample2NumberofSubsDict = getSample2NumberofSubsDict(outputDir,jobname)
        sample2NumberofIndelsDict = getSample2NumberofIndelsDict(outputDir,jobname)
        sample2NumberofDinucsDict = getDictionary(outputDir,jobname, Sample2NumberofDinucsDictFilename)

        sample2SubsSignature2NumberofMutationsDict = getSample2SubsSignature2NumberofMutationsDict(outputDir,jobname)
        sample2IndelsSignature2NumberofMutationsDict = getSample2IndelsSignature2NumberofMutationsDict(outputDir,jobname)
        sample2DinucsSignature2NumberofMutationsDict = getDictionary(outputDir, jobname,Sample2DinucsSignature2NumberofMutationsDictFilename)
        ##########################################################################
    else:
        ##########################################################################
        sample2NumberofSubsDict = {}
        sample2NumberofIndelsDict = {}
        sample2NumberofDinucsDict = {}

        sample2SubsSignature2NumberofMutationsDict ={}
        sample2IndelsSignature2NumberofMutationsDict = {}
        sample2DinucsSignature2NumberofMutationsDict = {}
        ##########################################################################

    ##########################################################################
    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)
    ##########################################################################

    simNum2Type2AccumulatedSignalArrayDict = {}
    simNum2Type2AccumulatedCountArrayDict = {}
    simNum2Sample2Type2AccumulatedSignalArrayDict = {}
    simNum2Sample2Type2AccumulatedCountArrayDict = {}

    ##############################################################
    #What is the type of the signal_file_with_path?
    #If it is a bed file read signal_file_with_path here
    file_extension = os.path.basename(library_file_with_path).split('.')[-1]

    library_file_df=None
    library_file_df_grouped=None

    if ((file_extension.lower()=='bigwig') or (file_extension.lower()=='bw')):
        library_file_type=BIGWIG
    elif ((file_extension.lower()=='bigbed') or (file_extension.lower()=='bb')):
        library_file_type=BIGBED
    elif (file_extension.lower()=='bed'):
        library_file_type=BED
        library_file_df=readFileInBEDFormat(library_file_with_path)
        library_file_df_grouped=library_file_df.groupby(chrom)
    elif (file_extension.lower()=='wig'):
        library_file_type=WIG
    elif (file_extension.lower()=='narrowpeak'):
        library_file_type = NARROWPEAK
        library_file_df = readFileInNarrowPeakFormat(library_file_with_path)
        library_file_df_grouped = library_file_df.groupby(chrom)
    else:
        library_file_type=LIBRARY_FILE_TYPE_OTHER
    ##############################################################

    # using a list
    # results = []

    #Using pool.apply_async
    # accumulator = Accumulator()

    # poolInputList = []

    # #October 18, 2019 starts
    # poolInputList=[fillInputList(outputDir, jobname,chrLong,simNum,chromSizesDict,library_file_with_path,library_file_type,library_file_df,library_file_df_grouped,
    #                             sample2NumberofSubsDict,sample2NumberofIndelsDict,sample2NumberofDinucsDict,
    #                             sample2SubsSignature2NumberofMutationsDict,sample2IndelsSignature2NumberofMutationsDict,sample2DinucsSignature2NumberofMutationsDict,
    #                             subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,plusorMinus)  for chrLong in chromNamesList for simNum in range(0,numofSimulations+1)]
    # #October 18, 2019 ends


    ###################################################################################
    ##################  For all chromsomes sequential starts ##########################
    ###################################################################################
    for chrLong in chromNamesList:

        # print('Debug we are in chrLong: %s' %(chrLong))
        poolInputList = []

        # chromSize = chromSizesDict[chrLong]
        # chrBasedSignalArray=None
        #
        # #################################################################################################################
        # # if (os.path.exists(chrBasedSignalNucleosmeFile) and os.path.exists(chrBasedCountNucleosmeFile)):
        # if (os.path.exists(library_file_with_path)):
        #     #TODO: This is specific to our data right now
        #     #Nucleosomes have chrM
        #     #SinglePointMutations and Indels have chrMT
        #     chrLong_for_mutations_data = chrLong
        #     if (chrLong=='chrM'):
        #         chrLong_for_mutations_data='chrMT'
        #
        #     #chrBasedSignalArray will be filled when library_file_type is BED or NARROWPEAK otherwise it will be None
        #     if (((library_file_type==BED) or (library_file_type==NARROWPEAK)) and (library_file_df is not None)):
        #          chrom_based_library_df=library_file_df_grouped.get_group(chrLong)
        #          #chrBasedSignalArray and library_file_df  signal column is of type np.float32
        #          chrBasedSignalArray = np.zeros(chromSize,dtype=np.float32)
        #          #TODO Can we fill chrBasedSignalArray faster?
        #          # chrom_based_library_df.apply(updateChrBasedSignalArray, chrBasedSignalArray=chrBasedSignalArray, axis=1)
        #          [fillNumpyArray(start,end,signal,chrBasedSignalArray) for start,end,signal in zip(chrom_based_library_df['start'],chrom_based_library_df['end'],chrom_based_library_df['signal'])]
        #
        #     ###################################################################################
        #     ##################  For all simulations in parallel starts ########################
        #     ###################################################################################
        #     for simNum in range(0,numofSimulations+1):
        #         print('simNum: %d' % (simNum))
        #         inputList = []
        #
        #         # READ CHRBASED SINGLE POINT MUTATIONS
        #         # READ CHRBASED INDELS
        #         # READ CHRBASED_DINUCS
        #         chrBased_spms_df = readChrBasedSubsDF(outputDir, jobname, chrLong_for_mutations_data, SUBS, simNum)
        #         chrBased_indels_df = readChrBasedMutationsDF(outputDir, jobname, chrLong_for_mutations_data, INDELS,simNum)
        #         chrBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong_for_mutations_data, DINUCS,simNum)
        #
        #         # ########## For information starts ##########
        #         # if (chrBased_spms_df is not None):
        #         #     print('For sim%d chromosome %s  -- chrBased_spms_df: %d in Bytes %f in GigaBytes' %(simNum,chrLong_for_mutations_data,sys.getsizeof(chrBased_spms_df),sys.getsizeof(chrBased_spms_df)/GIGABYTE_IN_BYTES))
        #         #
        #         # if (chrBased_indels_df is not None):
        #         #     print('For sim%d chromosome %s  -- chrBased_indels_df: %d in Bytes %f in GigaBytes' %(simNum,chrLong_for_mutations_data,sys.getsizeof(chrBased_indels_df),sys.getsizeof(chrBased_indels_df)/GIGABYTE_IN_BYTES))
        #         #
        #         # if (chrBased_dinucs_df is not None):
        #         #     print('For sim%d chromosome %s  -- chrBased_dinucss_df: %d in Bytes %f in GigaBytes' %(simNum,chrLong_for_mutations_data,sys.getsizeof(chrBased_dinucs_df), sys.getsizeof(chrBased_dinucs_df)/GIGABYTE_IN_BYTES))
        #         # ########## For information ends ##########
        #
        #         inputList.append(chrLong)
        #         inputList.append(library_file_with_path)
        #         inputList.append(library_file_type)
        #         inputList.append(chrBasedSignalArray)
        #         inputList.append(chrBased_spms_df)
        #         inputList.append(chrBased_indels_df)
        #         inputList.append(chrBased_dinucs_df)
        #         inputList.append(chromSize)
        #         inputList.append(sample2NumberofSubsDict)
        #         inputList.append(sample2NumberofIndelsDict)
        #         inputList.append(sample2NumberofDinucsDict)
        #         inputList.append(sample2SubsSignature2NumberofMutationsDict)
        #         inputList.append(sample2IndelsSignature2NumberofMutationsDict)
        #         inputList.append(sample2DinucsSignature2NumberofMutationsDict)
        #         inputList.append(subsSignature2PropertiesListDict)
        #         inputList.append(indelsSignature2PropertiesListDict)
        #         inputList.append(dinucsSignature2PropertiesListDict)
        #         inputList.append(plusorMinus)
        #         poolInputList.append(inputList)
        #
        #         #using a list
        #         # results.append(pool.apply_async(fillSignalArrayAndCountArrays_using_pyBigWig_using_apply_async,(chrLong,library_file_with_path,library_file_type,chrBasedSignalArray,chrBased_spms_df,chrBased_indels_df,chrBased_dinucs_df,chromSize,sample2NumberofSubsDict,sample2NumberofIndelsDict,sample2NumberofDinucsDict,sample2SubsSignature2NumberofMutationsDict,sample2IndelsSignature2NumberofMutationsDict,sample2DinucsSignature2NumberofMutationsDict,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,plusorMinus,)))
        #
        #         #using callback
        #         # pool.apply_async(fillSignalArrayAndCountArrays_using_pyBigWig_using_apply_async, (chrLong, library_file_with_path, library_file_type, chrBasedSignalArray, chrBased_spms_df, chrBased_indels_df, chrBased_dinucs_df, chromSize, sample2NumberofSubsDict, sample2NumberofIndelsDict,
        #         # sample2NumberofDinucsDict, sample2SubsSignature2NumberofMutationsDict, sample2IndelsSignature2NumberofMutationsDict, sample2DinucsSignature2NumberofMutationsDict, subsSignature2PropertiesListDict, indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict, plusorMinus,),callback=accumulator.on_result)
        #     ###################################################################################
        #     ##################  For all simulations in parallel ends ##########################
        #     ###################################################################################

        # poolInputList=[fillInputList(outputDir, jobname,chrLong,simNum,chromSizesDict,library_file_with_path,library_file_type,library_file_df,library_file_df_grouped,
        #                             sample2NumberofSubsDict,sample2NumberofIndelsDict,sample2NumberofDinucsDict,
        #                             sample2SubsSignature2NumberofMutationsDict,sample2IndelsSignature2NumberofMutationsDict,sample2DinucsSignature2NumberofMutationsDict,
        #                             subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,plusorMinus) for simNum in range(0,numofSimulations+1)]

        ########################################################################
        # Provides list of  lists (each sub list comes from a simulation)
        # allSimulations_SignalArrayAndCountArrayList_List = pool.map(fillSignalArrayAndCountArrays_using_pyBigWig,poolInputList)

        # chunksize = len(poolInputList) // numofProcesses
        # if chunksize < 1:
        #     chunksize = 1
        # print('chrLong:%s len(poolInputList):%d numofProcesses:%d chunksize:%d' %(chrLong,len(poolInputList),numofProcesses,chunksize))

        # allSimulations_SignalArrayAndCountArrayList_List = pool.map(fillSignalArrayAndCountArrays_using_pyBigWig,poolInputList, chunksize=chunksize)
        ########################################################################

        ########################################################################
        chunksize = max(1,(numofSimulations//numofProcesses))
        # print('Debug we are in chrLong:%s chunkSize:%d' %(chrLong,chunksize))

        for simulatonBased_SignalArrayAndCountArrayList in pool.imap_unordered(fillSignalArrayAndCountArrays_using_pyBigWig,(fillInputList(outputDir, jobname,chrLong,simNum,chromSizesDict,library_file_with_path,library_file_type,library_file_df,library_file_df_grouped,
                                    sample2NumberofSubsDict,sample2NumberofIndelsDict,sample2NumberofDinucsDict,
                                    sample2SubsSignature2NumberofMutationsDict,sample2IndelsSignature2NumberofMutationsDict,sample2DinucsSignature2NumberofMutationsDict,
                                    subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,plusorMinus) for simNum in range(0,numofSimulations+1)), chunksize=chunksize):

            # print('Debug We are in accumulate for chrLong:%s' %(chrLong))

            simNum2Type2SignalArrayDict = simulatonBased_SignalArrayAndCountArrayList[0]
            simNum2Type2CountArrayDict = simulatonBased_SignalArrayAndCountArrayList[1]
            simNum2Sample2Type2SignalArrayDict = simulatonBased_SignalArrayAndCountArrayList[2]
            simNum2Sample2Type2CountArrayDict = simulatonBased_SignalArrayAndCountArrayList[3]

            accumulateSimulationBasedTypeBasedArrays(simNum2Type2AccumulatedSignalArrayDict, simNum2Type2SignalArrayDict)
            accumulateSimulationBasedTypeBasedArrays(simNum2Type2AccumulatedCountArrayDict, simNum2Type2CountArrayDict)
            accumulateSimulationBasedSampleBasedTypeBasedArrays(simNum2Sample2Type2AccumulatedSignalArrayDict,simNum2Sample2Type2SignalArrayDict)
            accumulateSimulationBasedSampleBasedTypeBasedArrays(simNum2Sample2Type2AccumulatedCountArrayDict,simNum2Sample2Type2CountArrayDict)
        ########################################################################


        # #####################################################################################################
        # ######################### Accumulate right in the left starts  ######################################
        # #####################################################################################################
        # for simulatonBased_SignalArrayAndCountArrayList in allSimulations_SignalArrayAndCountArrayList_List:
        #     simNum2Type2SignalArrayDict = simulatonBased_SignalArrayAndCountArrayList[0]
        #     simNum2Type2CountArrayDict = simulatonBased_SignalArrayAndCountArrayList[1]
        #     simNum2Sample2Type2SignalArrayDict = simulatonBased_SignalArrayAndCountArrayList[2]
        #     simNum2Sample2Type2CountArrayDict = simulatonBased_SignalArrayAndCountArrayList[3]
        #
        #     accumulateSimulationBasedTypeBasedArrays(simNum2Type2AccumulatedSignalArrayDict, simNum2Type2SignalArrayDict)
        #     accumulateSimulationBasedTypeBasedArrays(simNum2Type2AccumulatedCountArrayDict, simNum2Type2CountArrayDict)
        #     accumulateSimulationBasedSampleBasedTypeBasedArrays(simNum2Sample2Type2AccumulatedSignalArrayDict,simNum2Sample2Type2SignalArrayDict)
        #     accumulateSimulationBasedSampleBasedTypeBasedArrays(simNum2Sample2Type2AccumulatedCountArrayDict,simNum2Sample2Type2CountArrayDict)
        # #####################################################################################################
        # ######################### Accumulate right in the left ends  ########################################
        # #####################################################################################################

    ###################################################################################
    ##################  For all chromsomes sequential ends ############################
    ###################################################################################

    # using a list
    # while results:
    #     for result in results[:]:
    #         if result.ready():
    #             # print('{} is ready'.format(r))
    #             simNum2Type2SignalArrayDict = result.get()[0]
    #             simNum2Type2CountArrayDict = result.get()[1]
    #             simNum2Sample2Type2SignalArrayDict = result.get()[2]
    #             simNum2Sample2Type2CountArrayDict = result.get()[3]
    #             accumulateSimulationBasedTypeBasedArrays(simNum2Type2AccumulatedSignalArrayDict, simNum2Type2SignalArrayDict)
    #             accumulateSimulationBasedTypeBasedArrays(simNum2Type2AccumulatedCountArrayDict, simNum2Type2CountArrayDict)
    #             accumulateSimulationBasedSampleBasedTypeBasedArrays(simNum2Sample2Type2AccumulatedSignalArrayDict,simNum2Sample2Type2SignalArrayDict)
    #             accumulateSimulationBasedSampleBasedTypeBasedArrays(simNum2Sample2Type2AccumulatedCountArrayDict,simNum2Sample2Type2CountArrayDict)
    #             results.remove(result)



    # ########################################################################
    # #Provides list of  lists (each sub list comes from a simulation)
    # # allSimulations_SignalArrayAndCountArrayList_List = pool.map(fillSignalArrayAndCountArrays_using_pyBigWig,poolInputList)
    # chunksize=len(poolInputList)//numofProcesses
    # if chunksize<1:
    #     chunksize=1
    # allSimulations_SignalArrayAndCountArrayList_List = pool.imap(fillSignalArrayAndCountArrays_using_pyBigWig,poolInputList,chunksize=chunksize)
    # ########################################################################
    #
    # #####################################################################################################
    # ######################### Accumulate right in the left starts  ######################################
    # #####################################################################################################
    # for simulatonBased_SignalArrayAndCountArrayList in  allSimulations_SignalArrayAndCountArrayList_List:
    #     simNum2Type2SignalArrayDict = simulatonBased_SignalArrayAndCountArrayList[0]
    #     simNum2Type2CountArrayDict = simulatonBased_SignalArrayAndCountArrayList[1]
    #     simNum2Sample2Type2SignalArrayDict = simulatonBased_SignalArrayAndCountArrayList[2]
    #     simNum2Sample2Type2CountArrayDict = simulatonBased_SignalArrayAndCountArrayList[3]
    #
    #     accumulateSimulationBasedTypeBasedArrays(simNum2Type2AccumulatedSignalArrayDict,simNum2Type2SignalArrayDict)
    #     accumulateSimulationBasedTypeBasedArrays(simNum2Type2AccumulatedCountArrayDict, simNum2Type2CountArrayDict)
    #     accumulateSimulationBasedSampleBasedTypeBasedArrays(simNum2Sample2Type2AccumulatedSignalArrayDict,simNum2Sample2Type2SignalArrayDict)
    #     accumulateSimulationBasedSampleBasedTypeBasedArrays(simNum2Sample2Type2AccumulatedCountArrayDict,simNum2Sample2Type2CountArrayDict)
    #
    # #####################################################################################################
    # ######################### Accumulate right in the left ends  ########################################
    # #####################################################################################################

    ###################################################################################
    ##################  For all chromsomes sequential ends ############################
    ###################################################################################

    ################################
    pool.close()
    pool.join()
    ################################

    writeSimulationBasedAverageNucleosomeOccupancy(occupancy_type,
                                                   sample_based,
                                                   plusorMinus,
                                                   simNum2Type2AccumulatedSignalArrayDict,
                                                   simNum2Type2AccumulatedCountArrayDict,
                                                   simNum2Sample2Type2AccumulatedSignalArrayDict,
                                                   simNum2Sample2Type2AccumulatedCountArrayDict,
                                                   outputDir, jobname,library_file_memo)


#Using pyBigWig for bigBed and bigWig files ends
#Using bed files prepared on the fly ends
########################################################################################



########################################################################################
#For each chromosome sequential starts
#type can be NUCLEOSOME or HISTONE_MODIFICATION
def occupancy_analysis_eachChrSequential_allSimsParallel(chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,nucleosomeFilename,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,occupancy_type):

    ##########################################################################
    sample2NumberofSubsDict = getSample2NumberofSubsDict(outputDir,jobname)
    sample2NumberofIndelsDict = getSample2NumberofIndelsDict(outputDir,jobname)
    sample2NumberofDinucsDict = getDictionary(outputDir,jobname, Sample2NumberofDinucsDictFilename)

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
        nucleosomeFilenameWoExtension = os.path.splitext(os.path.basename(nucleosomeFilename))[0]

        ##############################################################
        signalArrayFilename = '%s_signal_%s.npy' % (chrLong, nucleosomeFilenameWoExtension)

        if (occupancy_type==NUCLEOSOMEOCCUPANCY):
            #For Nucleosome Occupancy
            chrBasedSignalFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB,NUCLEOSOME, CHRBASED, signalArrayFilename)
        elif (occupancy_type==EPIGENOMICSOCCUPANCY):
            #For HM
            chrBasedSignalFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, EPIGENOMICS, CHRBASED, signalArrayFilename)

        #################################################################################################################
        # if (os.path.exists(chrBasedSignalNucleosmeFile) and os.path.exists(chrBasedCountNucleosmeFile)):
        if (os.path.exists(chrBasedSignalFile)):
            chrbased_signal_array = np.load(chrBasedSignalFile)
            print('chromosome %s  -- signal_array_npy: %d in bytes %f in GB' % (chrLong,sys.getsizeof(chrbased_signal_array), sys.getsizeof(chrbased_signal_array)/GIGABYTE_IN_BYTES))

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

                inputList.append(chrbased_signal_array)
                inputList.append(chrBased_spms_df)
                inputList.append(chrBased_indels_df)
                inputList.append(chrBased_dinucs_df)
                inputList.append(chromSize)
                inputList.append(sample2NumberofSubsDict)
                inputList.append(sample2NumberofIndelsDict)
                inputList.append(sample2NumberofDinucsDict)
                inputList.append(sample2SubsSignature2NumberofMutationsDict)
                inputList.append(sample2IndelsSignature2NumberofMutationsDict)
                inputList.append(sample2DinucsSignature2NumberofMutationsDict)
                inputList.append(subsSignature2PropertiesListDict)
                inputList.append(indelsSignature2PropertiesListDict)
                inputList.append(dinucsSignature2PropertiesListDict)
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

    writeSimulationBasedAverageNucleosomeOccupancy(occupancy_type,
                                                   simNum2Type2AccumulatedSignalArrayDict,
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