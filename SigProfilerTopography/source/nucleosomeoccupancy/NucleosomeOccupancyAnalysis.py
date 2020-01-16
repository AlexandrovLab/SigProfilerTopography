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

import time
import sys
import multiprocessing
import os
import pandas as pd
import numpy as np
import math

from SigProfilerTopography.source.commons.TopographyCommons import memory_usage
from SigProfilerTopography.source.commons.TopographyCommons import readChrBasedMutationsDF
from SigProfilerTopography.source.commons.TopographyCommons import func_addSignal
from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofSubsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofIndelsDict
from SigProfilerTopography.source.commons.TopographyCommons import getDictionary
from SigProfilerTopography.source.commons.TopographyCommons import accumulateSimulationBasedTypeBasedArrays
from SigProfilerTopography.source.commons.TopographyCommons import accumulateSimulationBasedSampleBasedTypeBasedArrays

from SigProfilerTopography.source.commons.TopographyCommons import getSample2SubsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2IndelsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import writeSimulationBasedAverageNucleosomeOccupancy

from SigProfilerTopography.source.commons.TopographyCommons import SUBS
from SigProfilerTopography.source.commons.TopographyCommons import INDELS
from SigProfilerTopography.source.commons.TopographyCommons import DINUCS
from SigProfilerTopography.source.commons.TopographyCommons import GIGABYTE_IN_BYTES

from SigProfilerTopography.source.commons.TopographyCommons import EPIGENOMICSOCCUPANCY
from SigProfilerTopography.source.commons.TopographyCommons import NUCLEOSOMEOCCUPANCY
from SigProfilerTopography.source.commons.TopographyCommons import ONE_DIRECTORY_UP
from SigProfilerTopography.source.commons.TopographyCommons import LIB
from SigProfilerTopography.source.commons.TopographyCommons import DATA
from SigProfilerTopography.source.commons.TopographyCommons import EPIGENOMICS
from SigProfilerTopography.source.commons.TopographyCommons import NUCLEOSOME
from SigProfilerTopography.source.commons.TopographyCommons import CHRBASED

from SigProfilerTopography.source.commons.TopographyCommons import current_abs_path

from SigProfilerTopography.source.commons.TopographyCommons import BED
from SigProfilerTopography.source.commons.TopographyCommons import NARROWPEAK
from SigProfilerTopography.source.commons.TopographyCommons import BIGBED
from SigProfilerTopography.source.commons.TopographyCommons import BIGWIG
from SigProfilerTopography.source.commons.TopographyCommons import WIG
from SigProfilerTopography.source.commons.TopographyCommons import LIBRARY_FILE_TYPE_OTHER

from SigProfilerTopography.source.commons.TopographyCommons import BED_6PLUS4
from SigProfilerTopography.source.commons.TopographyCommons import BED_9PLUS2

from SigProfilerTopography.source.commons.TopographyCommons import AGGREGATEDSUBSTITUTIONS
from SigProfilerTopography.source.commons.TopographyCommons import AGGREGATEDINDELS
from SigProfilerTopography.source.commons.TopographyCommons import AGGREGATEDDINUCS

from SigProfilerTopography.source.commons.TopographyCommons import SAMPLE
from SigProfilerTopography.source.commons.TopographyCommons import START
from SigProfilerTopography.source.commons.TopographyCommons import SIMULATION_NUMBER

from SigProfilerTopography.source.commons.TopographyCommons import Sample2NumberofDinucsDictFilename
from SigProfilerTopography.source.commons.TopographyCommons import Sample2DinucsSignature2NumberofMutationsDictFilename

from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC
from SigProfilerTopography.source.commons.TopographyCommons import USING_IMAP_UNORDERED

from SigProfilerTopography.source.nucleosomeoccupancy.ChrBasedSignalArrays import readBEDandWriteChromBasedSignalArrays
from SigProfilerTopography.source.nucleosomeoccupancy.ChrBasedSignalArrays import readWig_with_fixedStep_variableStep_writeChrBasedSignalArrays
from SigProfilerTopography.source.nucleosomeoccupancy.ChrBasedSignalArrays import readWig_write_derived_from_bedgraph
from SigProfilerTopography.source.nucleosomeoccupancy.ChrBasedSignalArrays import readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArraysInParallel
from SigProfilerTopography.source.nucleosomeoccupancy.ChrBasedSignalArrays import readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArraysSequentially

from SigProfilerTopography.source.commons.TopographyCommons import decideFileType

##############################################################################################################
#main function
def occupancyAnalysis(genome,computationType,occupancy_type,sample_based,plusorMinus,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,library_file_with_path,library_file_memo,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,verbose):

    print('\n#################################################################################')
    print('--- %s Analysis starts' %(occupancy_type))
    print('--- Computation Type:%s' % (computationType))
    print('--- Occupancy Type:%s' % (occupancy_type))
    print('--- Library file with path: %s\n' %library_file_with_path)

    #Using pyBigWig for BigWig and BigBed files if operating system is not Windows and using_pyBigWig is set to True
    #Using Bed files preparing chr based signal array runtime
    occupancy_analysis(genome,computationType,occupancy_type,sample_based,plusorMinus,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,library_file_with_path,library_file_memo,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,verbose)
    print('--- %s Analysis ends' %(occupancy_type))
    print('#################################################################################\n')
##############################################################################################################



########################################################################################
# November 1, 2019
# Just fill the list
def fillInputList(occupancy_type,outputDir,jobname,chrLong,simNum,chromSizesDict,library_file_with_path,library_file_type,sample2NumberofSubsDict,sample2NumberofIndelsDict,sample2NumberofDinucsDict,
                      sample2SubsSignature2NumberofMutationsDict,sample2IndelsSignature2NumberofMutationsDict,sample2DinucsSignature2NumberofMutationsDict,
                      subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,plusorMinus,sample_based,verbose):

    if verbose: print('FillInputList: Worker pid %s current_mem_usage %.2f (mb) chrLong:%s simNum:%d' % (str(os.getpid()), memory_usage(),chrLong,simNum))
    inputList=[]

    inputList.append(occupancy_type)
    inputList.append(outputDir)
    inputList.append(jobname)
    inputList.append(chrLong)
    inputList.append(simNum)
    inputList.append(chromSizesDict)
    inputList.append(library_file_with_path)
    inputList.append(library_file_type)
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
    inputList.append(sample_based)
    inputList.append(verbose)
    # print('FillInputList: Worker pid %s maximum memory usage %.2f (mb)' % (str(os.getpid()), current_mem_usage()))
    return inputList
########################################################################################


########################################################################################
#Updated JAN 9, 2020
#November 26 2019
#Called from apply_async
def combined_prepare_chrbased_data_fill_signal_count_arrays(occupancy_type,
        outputDir,
        jobname, chrLong, simNum, chromSizesDict,
        library_file_with_path,
        library_file_type,
        sample2NumberofSubsDict, sample2NumberofIndelsDict, sample2NumberofDinucsDict,
        sample2SubsSignature2NumberofMutationsDict, sample2IndelsSignature2NumberofMutationsDict,
        sample2DinucsSignature2NumberofMutationsDict,
        subsSignature2PropertiesListDict, indelsSignature2PropertiesListDict, dinucsSignature2PropertiesListDict,
        plusorMinus, sample_based,verbose):

    if verbose: print('\tWorker pid %s memory_usage in %.2f MB START chrLong:%s simNum:%d' %(str(os.getpid()), memory_usage(), chrLong, simNum))
    # 1st part  Prepare chr based mutations dataframes
    maximum_chrom_size = chromSizesDict[chrLong]

    ##############################################################
    simNum2Type2SignalArrayDict = {}
    simNum2Type2CountArrayDict = {}
    simNum2Sample2Type2SignalArrayDict = {}
    simNum2Sample2Type2CountArrayDict = {}
    ##############################################################

    ##############################################################
    chrBasedSignalArray = None #Will be filled from chrBasedSignal files if they exists
    library_file_opened_by_pyBigWig = None #Will be filled by pyBigWig from bigWig or bigBed
    my_upperBound = None
    signal_index = None
    ##############################################################

    #################################################################################################################
    # If library file does not exists there is no library file to use and fill the signal and count arrays
    # Nucleosomes have chrM
    # SinglePointMutations and Indels have chrMT
    chrLong_for_mutations_data = chrLong
    if (chrLong == 'chrM'):
        chrLong_for_mutations_data = 'chrMT'

    chrBased_subs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong_for_mutations_data, SUBS, simNum)
    chrBased_indels_df = readChrBasedMutationsDF(outputDir, jobname, chrLong_for_mutations_data, INDELS, simNum)
    chrBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong_for_mutations_data, DINUCS, simNum)
    if verbose: print('\tWorker pid %s memory_usage in %.2f MB Check1 Read Signal Array and Dataframes chrLong:%s simNum:%d' % (str(os.getpid()), memory_usage(), chrLong, simNum))
    if verbose: print('\tWorker pid %s -- signal_array_npy: %f in GB -- subs_df: %f in GB -- indels_df: %f in GB -- dinucs_df: %f in GB -- chrLong:%s simNum:%d' % (
            str(os.getpid()),
            sys.getsizeof(chrBasedSignalArray) / GIGABYTE_IN_BYTES,
            sys.getsizeof(chrBased_subs_df) / GIGABYTE_IN_BYTES,
            sys.getsizeof(chrBased_indels_df) / GIGABYTE_IN_BYTES,
            sys.getsizeof(chrBased_dinucs_df) / GIGABYTE_IN_BYTES,
            chrLong, simNum))
    #################################################################################################################

    #################################################################################################################
    libraryFilenameWoExtension = os.path.splitext(os.path.basename(library_file_with_path))[0]
    signalArrayFilename = '%s_signal_%s.npy' % (chrLong, libraryFilenameWoExtension)
    if (occupancy_type== EPIGENOMICSOCCUPANCY):
        chrBasedSignalFile = os.path.join(outputDir,jobname,DATA,occupancy_type,LIB,CHRBASED,signalArrayFilename)
    elif (occupancy_type==NUCLEOSOMEOCCUPANCY):
        chrBasedSignalFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME, CHRBASED,signalArrayFilename)

    #Downloaded or created runtime
    if (os.path.exists(chrBasedSignalFile)):
        chrBasedSignalArray = np.load(chrBasedSignalFile, mmap_mode='r')
        # For testing purposes
        # chrBasedSignalArray = np.load(chrBasedSignalFile,mmap_mode=None)
        # chrBasedSignalArray = np.random.uniform(low=0.0, high=13.3, size=(maximum_chrom_size,))

    #If library_file_with_path is abs path and library_file_type is BIGWIG or BIGBED
    #For nucleosome_biosample==GM12878 or nucleosome_biosample==K562 library_file_with_path is only filename with extension
    if os.path.isabs(library_file_with_path):
        # Comment below to make it run in windows
        if (library_file_type == BIGWIG):
            try:
                import pyBigWig
                library_file_opened_by_pyBigWig = pyBigWig.open(library_file_with_path)
                if chrLong in library_file_opened_by_pyBigWig.chroms():
                    maximum_chrom_size = library_file_opened_by_pyBigWig.chroms()[chrLong]
                # For BigWig Files information in header is correct
                if ('sumData' in library_file_opened_by_pyBigWig.header()) and ('nBasesCovered' in library_file_opened_by_pyBigWig.header()):
                    my_mean = library_file_opened_by_pyBigWig.header()['sumData'] / library_file_opened_by_pyBigWig.header()['nBasesCovered']
                    std_dev = (library_file_opened_by_pyBigWig.header()['sumSquared'] - 2 * my_mean * library_file_opened_by_pyBigWig.header()['sumData'] +
                               library_file_opened_by_pyBigWig.header()['nBasesCovered'] * my_mean * my_mean) ** (0.5) / (
                            library_file_opened_by_pyBigWig.header()['nBasesCovered'] ** (0.5))
                    # Scientific definition of outlier
                    my_upperBound = std_dev * 3
                else:
                    # Undefined
                    my_upperBound = np.iinfo(np.int16).max
            except:
                print('Exception %s' %library_file_with_path)

        elif (library_file_type == BIGBED):
            try:
                import pyBigWig
                library_file_opened_by_pyBigWig = pyBigWig.open(library_file_with_path)
                if BED_6PLUS4 in str(library_file_opened_by_pyBigWig.SQL()):
                    signal_index = 3
                elif BED_9PLUS2 in str(library_file_opened_by_pyBigWig.SQL()):
                    signal_index = 7
                if chrLong in library_file_opened_by_pyBigWig.chroms():
                    # For BigBed Files information in header is not meaningful
                    maximum_chrom_size = library_file_opened_by_pyBigWig.chroms()[chrLong]
                    my_mean = np.mean([float(entry[2].split('\t')[signal_index]) for entry in
                                       library_file_opened_by_pyBigWig.entries(chrLong, 0, maximum_chrom_size)])
                    # Not scientific definition of outlier
                    my_upperBound = my_mean * 10
                else:
                    # Undefined
                    my_upperBound = np.iinfo(np.int16).max
            except:
                print('Exception %s' %library_file_with_path)

    #################################################################################################################

    #################################################################################################################
    if ((chrBasedSignalArray is not None) or ((library_file_opened_by_pyBigWig is not None) and (chrLong in library_file_opened_by_pyBigWig.chroms()))):
        ######################################################## #######################
        ################### Fill signal and count array starts ########################
        ###############################################################################

        if verbose: print('\tWorker pid %s memory_usage in %.2f MB Check2_1 Dinucs Start chrLong:%s simNum:%d' % (str(os.getpid()), memory_usage(), chrLong, simNum))
        # Fill for dinucs
        if ((chrBased_dinucs_df is not None) and (not chrBased_dinucs_df.empty)):
            df_columns = list(chrBased_dinucs_df.columns.values)
            [fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated_using_pyBigWig_using_list_comp(
                row,
                chrLong,
                library_file_opened_by_pyBigWig,
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
                sample_based,
                df_columns) for row in chrBased_dinucs_df[df_columns].values]
        if verbose: print('\tWorker pid %s memory_usage in %.2f MB Check2_2 Dinucs End chrLong:%s simNum:%d' % (str(os.getpid()), memory_usage(), chrLong, simNum))

        if verbose: print('\tWorker pid %s memory_usage in %.2f MB Check3_1 Indels Start chrLong:%s simNum:%d' % (str(os.getpid()), memory_usage(), chrLong, simNum))
        # Fill for indels
        if ((chrBased_indels_df is not None) and (not chrBased_indels_df.empty)):
            df_columns = list(chrBased_indels_df.columns.values)
            [fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated_using_pyBigWig_using_list_comp(
                row,
                chrLong,
                library_file_opened_by_pyBigWig,
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
                sample_based,
                df_columns) for row in chrBased_indels_df[df_columns].values]
        if verbose: print('\tWorker pid %s memory_usage in %.2f MB Check3_2 Indels End chrLong:%s simNum:%d' % (str(os.getpid()), memory_usage(), chrLong, simNum))

        if verbose: print('\tWorker pid %s memory_usage in %.2f MB Check4_1 Subs Start chrLong:%s simNum:%d' % (str(os.getpid()), memory_usage(), chrLong, simNum))
        # Fill for single point mutations
        if ((chrBased_subs_df is not None) and (not chrBased_subs_df.empty)):
            list_of_dfs = None
            df_columns = list(chrBased_subs_df.columns.values)

            # 1 MB 1024*1024= 1048576 B
            size_in_mbs = sys.getsizeof(chrBased_subs_df) / 1048576
            if verbose: print('\tWorker pid %s ##################### subs_df: %f in MB chrLong:%s simNum:%d' % (str(os.getpid()),size_in_mbs, chrLong, simNum))
            max_size_in_mbs = 50
            if (size_in_mbs > max_size_in_mbs):
                numberofSplits = math.ceil(size_in_mbs / max_size_in_mbs)
                if verbose: print('\tWorker pid %s numberofSplits: %d chrLong:%s simNum:%d' % (str(os.getpid()),numberofSplits, chrLong, simNum))
                list_of_dfs = np.array_split(chrBased_subs_df, numberofSplits)

            # This is 3X-4X faster with almost same memory usage
            start_time = time.time()

            if list_of_dfs is not None:
                for part_index, part_df in enumerate(list_of_dfs, 1):
                    [fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated_using_pyBigWig_using_list_comp(
                        row,
                        chrLong,
                        library_file_opened_by_pyBigWig,
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
                        sample_based,
                        df_columns) for row in part_df[df_columns].values]
            else:
                [fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated_using_pyBigWig_using_list_comp(
                    row,
                    chrLong,
                    library_file_opened_by_pyBigWig,
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
                    sample_based,
                    df_columns) for row in chrBased_subs_df[df_columns].values]

        if verbose: print('\tWorker pid %s memory_usage in %.2f MB Check4_2 Subs End chrLong:%s simNum:%d' % (str(os.getpid()), memory_usage(), chrLong, simNum))
        ###############################################################################
        ################### Fill signal and count array ends ##########################
        ###############################################################################

    if (library_file_opened_by_pyBigWig is not None):
        library_file_opened_by_pyBigWig.close()

    if verbose: print('\tWorker pid %s memory_usage in %.2f MB END  chrLong:%s simNum:%d' % (str(os.getpid()), memory_usage(), chrLong, simNum))
    if verbose: print('----->\tWorker pid %s took %f seconds chrLong:%s simNum:%d' % (str(os.getpid()), (time.time() - start_time), chrLong, simNum))
    ###############################################################################
    ################### Return  starts ############################################
    ###############################################################################
    # Initialzie the list, you will return this list
    SignalArrayAndCountArrayList = []
    SignalArrayAndCountArrayList.append(simNum2Type2SignalArrayDict)
    SignalArrayAndCountArrayList.append(simNum2Type2CountArrayDict)
    SignalArrayAndCountArrayList.append(simNum2Sample2Type2SignalArrayDict)
    SignalArrayAndCountArrayList.append(simNum2Sample2Type2CountArrayDict)

    return SignalArrayAndCountArrayList
########################################################################################


########################################################################################
#Using list comprehension
#September 18, 2019
#You need to send mutation_row[START], mutation_row[SAMPLE], mutation_row[SIMULATION_NUMBER], and mutation_row[signature]
def fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated_using_pyBigWig_using_list_comp(
        row,
        chrLong,
        library_file_opened_by_pyBigWig,
        chrBasedSignalArray,
        library_file_type,
        signal_index,
        my_upperBound,
        maximum_chrom_size,
        sample2NumberofMutationsDict,
        sample2Signature2NumberofMutationsDict,
        simNum2Type2SignalArrayDict,
        simNum2Type2CountArrayDict,
        simNum2Sample2Type2SignalArrayDict,
        simNum2Sample2Type2CountArrayDict,
        signature2PropertiesListDict,
        my_type,
        plusOrMinus,
        sample_based,
        df_columns):

    window_array=None
    windowSize=plusOrMinus*2+1

    # df_columns: ['Sample', 'Chrom', 'Start', 'MutationLong', 'PyramidineStrand', 'TranscriptionStrand', 'Mutation',
    #              'SBS1', 'SBS2', 'SBS3', 'SBS4', 'SBS5', 'SBS6', 'SBS7a', 'SBS7b', 'SBS7c', 'SBS7d', 'SBS8', 'SBS9',
    #              'SBS10a', 'SBS10b', 'SBS11', 'SBS12', 'SBS13', 'SBS14', 'SBS15', 'SBS16', 'SBS17a', 'SBS17b', 'SBS18',
    #              'SBS19', 'SBS20', 'SBS21', 'SBS22', 'SBS23', 'SBS24', 'SBS25', 'SBS26', 'SBS27', 'SBS28', 'SBS29',
    #              'SBS30', 'SBS31', 'SBS32', 'SBS33', 'SBS34', 'SBS35', 'SBS36', 'SBS37', 'SBS38', 'SBS39', 'SBS40',
    #              'SBS41', 'SBS42', 'SBS43', 'SBS44', 'SBS45', 'SBS46', 'SBS47', 'SBS48', 'SBS49', 'SBS50', 'SBS51',
    #              'SBS52', 'SBS53', 'SBS54', 'SBS55', 'SBS56', 'SBS57', 'SBS58', 'SBS59', 'SBS60', 'Simulation_Number']

    indexofSample = df_columns.index(SAMPLE)
    indexofStart = df_columns.index(START)
    indexofSimulationNumber = df_columns.index(SIMULATION_NUMBER)

    mutation_row_sample = row[indexofSample]
    mutation_row_start = row[indexofStart]
    mutation_row_simulation_number = row[indexofSimulationNumber]

    # mutation_row_sample=row[0]
    # mutation_row_start=row[1]
    # mutation_row_simulation_number=row[2]

    #Get or fill window_array using Case1, Case2, and Case3
    # Case 1: start is very close to the chromosome start
    if (mutation_row_start<plusOrMinus):
        # print('Case 1: start is very close to the chromosome start --- mutation[Start]:%d' %(mutation_row_start))
        #Faster
        if (chrBasedSignalArray is not None):
            window_array = chrBasedSignalArray[0:(mutation_row_start + plusOrMinus + 1)]
            window_array = np.pad(window_array, (plusOrMinus - mutation_row_start, 0), 'constant', constant_values=(0, 0))

        # #Slower but uses less memory solution
        # elif (chrom_based_signal_df is not None):
        #     start=0
        #     end=mutation_row_start+plusOrMinus+1
        #     overlaps_df=chrom_based_signal_df[((start <= chrom_based_signal_df['end']) & (chrom_based_signal_df['start'] <= end))]
        #     if overlaps_df is not None:
        #         window_array = np.zeros((windowSize,), dtype=np.float32)
        #         #Overlap
        #         #chr1  240106575  240107544  3.78085
        #         [(func_addSignal(window_array, overlap[1], overlap[2], np.float32(overlap[3]),mutation_row_start, plusOrMinus)) for overlap in overlaps_df.values]

        elif (library_file_type==BIGWIG):
            #Important: The bigWig format does not support overlapping intervals.
            window_array=library_file_opened_by_pyBigWig.values(chrLong,0,(mutation_row_start+plusOrMinus+1),numpy=True)
            # How do you handle outliers?
            window_array[np.isnan(window_array)] = 0
            window_array[window_array>my_upperBound]=my_upperBound
            window_array = np.pad(window_array, (plusOrMinus - mutation_row_start, 0), 'constant',constant_values=(0, 0))

        elif (library_file_type==BIGBED):
            #We assume that in the 7th column there is signal data
            list_of_entries=library_file_opened_by_pyBigWig.entries(chrLong,0,(mutation_row_start+plusOrMinus+1))
            if list_of_entries is not None:
                window_array = np.zeros((windowSize,),dtype=np.float32)
                # We did not handle outliers for BigBed files.

                #From DNA methylation get the 7th
                # library_file_bed_format==BED_6PLUS4):
                # (713235, 713435, 'Peak_40281\t15\t.\t3.48949\t5.67543\t3.79089\t158')
                #signal_index=3
                #library_file_bed_format==BED_9PLUS2):
                #[(10810, 10811, 'MCF7_NoStarve_B1__GC_\t3\t+\t10810\t10811\t255,0,0\t3\t100'), (10812, 10813, 'MCF7_NoStarve_B1__GC_\t3\t+\t10812\t10813\t255,0,0\t3\t100'), (10815, 10816, 'MCF7_NoStarve_B1__GC_\t3\t+\t10815\t10816\t0,255,0\t3\t0')]
                #signal_index=7
                [(func_addSignal(window_array, entry[0], entry[1], np.float32(entry[2].split()[signal_index]),mutation_row_start, plusOrMinus) if len(entry) >= 3 else (func_addSignal(window_array, entry[0], entry[1], 1, mutation_row_start, plusOrMinus))) for entry in list_of_entries]

    # Case 2: start is very close to the chromosome end
    elif (mutation_row_start+plusOrMinus+1 > maximum_chrom_size):
        # print('Case2: start is very close to the chromosome end ---  mutation[Start]:%d' %(mutation_row_start))
        if ((chrBasedSignalArray is not None)):
            window_array = chrBasedSignalArray[(mutation_row_start-plusOrMinus):maximum_chrom_size]
            window_array = np.pad(window_array, (0,mutation_row_start+plusOrMinus-maximum_chrom_size+1),'constant',constant_values=(0,0))

        # elif (chrom_based_signal_df is not None):
        #     start=mutation_row_start-plusOrMinus
        #     end=maximum_chrom_size
        #     overlaps_df=chrom_based_signal_df[((start <= chrom_based_signal_df['end']) & (chrom_based_signal_df['start'] <= end))]
        #     if overlaps_df is not None:
        #         window_array = np.zeros((windowSize,), dtype=np.float32)
        #         #Overlap
        #         #chr1  240106575  240107544  3.78085
        #         [(func_addSignal(window_array, overlap[1], overlap[2], np.float32(overlap[3]),mutation_row_start, plusOrMinus)) for overlap in overlaps_df.values]

        elif (library_file_type==BIGWIG):
            #Important: The bigWig format does not support overlapping intervals.
            window_array = library_file_opened_by_pyBigWig.values(chrLong,(mutation_row_start-plusOrMinus),maximum_chrom_size,numpy=True)
            # How do you handle outliers?
            window_array[np.isnan(window_array)] = 0
            window_array[window_array>my_upperBound]=my_upperBound
            window_array = np.pad(window_array, (0,mutation_row_start+plusOrMinus-maximum_chrom_size+1),'constant',constant_values=(0,0))

        elif (library_file_type==BIGBED):
            # print('Case2 Debug Sep 5, 2019 %s mutation_row[START]:%d mutation_row[START]-plusOrMinus:%d maximum_chrom_size:%d' %(chrLong,mutation_row[START],mutation_row[START]-plusOrMinus,maximum_chrom_size))
            if ((mutation_row_start-plusOrMinus)<maximum_chrom_size):
                list_of_entries=library_file_opened_by_pyBigWig.entries(chrLong,(mutation_row_start-plusOrMinus),maximum_chrom_size)
                if list_of_entries is not None:
                    window_array = np.zeros((windowSize,),dtype=np.float32)
                    # We did not handle outliers for BigBed files.
                    [(func_addSignal(window_array, entry[0], entry[1], np.float32(entry[2].split()[signal_index]),mutation_row_start,plusOrMinus) if len(entry) >= 3 else (func_addSignal(window_array, entry[0], entry[1],1, mutation_row_start,plusOrMinus))) for entry in list_of_entries]


    #Case 3: No problem
    else:
        if (chrBasedSignalArray is not None):
            window_array = chrBasedSignalArray[(mutation_row_start-plusOrMinus):(mutation_row_start+plusOrMinus+1)]

        # elif (chrom_based_signal_df is not None):
        #     start=mutation_row_start-plusOrMinus
        #     end=mutation_row_start+plusOrMinus+1
        #     overlaps_df=chrom_based_signal_df[((start <= chrom_based_signal_df['end']) & (chrom_based_signal_df['start'] <= end))]
        #     if overlaps_df is not None:
        #         window_array = np.zeros((windowSize,), dtype=np.float32)
        #         #Overlap
        #         #chr1  240106575  240107544  3.78085
        #         [(func_addSignal(window_array, overlap[1], overlap[2], np.float32(overlap[3]),mutation_row_start, plusOrMinus)) for overlap in overlaps_df.values]

        elif (library_file_type==BIGWIG):
            #Important: You have to go over intervals if there are overlapping intervals.
            window_array = library_file_opened_by_pyBigWig.values(chrLong, (mutation_row_start-plusOrMinus), (mutation_row_start+plusOrMinus+1),numpy=True)
            #How do you handle outliers?
            window_array[np.isnan(window_array)] = 0
            window_array[window_array>my_upperBound]=my_upperBound

        elif (library_file_type==BIGBED):
            # print('Case3 Debug Sep 5, 2019 %s mutation_row[START]:%d mutation_row[START]-plusOrMinus:%d mutation_row[START]+plusOrMinus+1:%d' %(chrLong,mutation_row[START],mutation_row[START]-plusOrMinus,mutation_row[START]+plusOrMinus+1))
            if ((mutation_row_start+plusOrMinus+1)<=maximum_chrom_size):
                list_of_entries=library_file_opened_by_pyBigWig.entries(chrLong, (mutation_row_start-plusOrMinus), (mutation_row_start+plusOrMinus+1))
                if list_of_entries is not None:
                    window_array = np.zeros((windowSize,),dtype=np.float32)
                    # We did not handle outliers for BigBed files.
                    [(func_addSignal(window_array, entry[0], entry[1], np.float32(entry[2].split()[signal_index]),mutation_row_start,plusOrMinus) if len(entry) >= 3 else (func_addSignal(window_array, entry[0], entry[1],1, mutation_row_start,plusOrMinus))) for entry in list_of_entries]

    ##########################################################

    #Get the sample at this mutation_row
    sample = mutation_row_sample
    simulationNumber= mutation_row_simulation_number

    #####################################################
    if simulationNumber not in simNum2Type2SignalArrayDict:
        simNum2Type2SignalArrayDict[simulationNumber] = {}
        simNum2Type2CountArrayDict[simulationNumber] = {}

    type2SignalArrayDict = simNum2Type2SignalArrayDict[simulationNumber]
    type2CountArrayDict =  simNum2Type2CountArrayDict[simulationNumber]
    #####################################################

    #####################################################
    if sample_based:
        if simulationNumber not in simNum2Sample2Type2SignalArrayDict:
            simNum2Sample2Type2SignalArrayDict[simulationNumber] = {}
            simNum2Sample2Type2CountArrayDict[simulationNumber] = {}
        sample2Type2SignalArrayDict = simNum2Sample2Type2SignalArrayDict[simulationNumber]
        sample2Type2CountArrayDict = simNum2Sample2Type2CountArrayDict[simulationNumber]
    #####################################################

    #Fill dictionaries uisng window_array
    if (window_array is not None) and (np.any(window_array)):
        #TODO: Is there a faster way than using for loop?
        ################# Signatures starts #######################
        #mutation_row[signature] mutation probability for that signature
        #signature2PropertiesListDict[signature][0] cutoff probability for that signature
        for signature in signature2PropertiesListDict:
            indexofSignature = df_columns.index(signature)
            mutation_row_signature = row[indexofSignature]

            if (mutation_row_signature >= float(signature2PropertiesListDict[signature][0])):
                if (signature in type2SignalArrayDict):
                    type2SignalArrayDict[signature] += window_array
                    type2CountArrayDict[signature] += (window_array>0)
                else:
                    type2SignalArrayDict[signature] = np.zeros(windowSize)
                    type2CountArrayDict[signature] = np.zeros(windowSize, dtype=int)
                    type2SignalArrayDict[signature] += window_array
                    type2CountArrayDict[signature] += (window_array>0)

                ####################################################
                if sample_based:
                    if (sample in sample2Signature2NumberofMutationsDict) and (signature in sample2Signature2NumberofMutationsDict[sample]):
                        if sample in sample2Type2SignalArrayDict:
                            if signature in sample2Type2SignalArrayDict[sample]:
                                sample2Type2SignalArrayDict[sample][signature] += window_array
                                sample2Type2CountArrayDict[sample][signature] += (window_array>0)
                            else:
                                sample2Type2SignalArrayDict[sample][signature] = np.zeros(windowSize)
                                sample2Type2CountArrayDict[sample][signature] = np.zeros(windowSize, dtype=int)
                                sample2Type2SignalArrayDict[sample][signature] += window_array
                                sample2Type2CountArrayDict[sample][signature] += (window_array > 0)

                        else:
                            sample2Type2SignalArrayDict[sample] = {}
                            sample2Type2CountArrayDict[sample] = {}
                            sample2Type2SignalArrayDict[sample][signature] = np.zeros(windowSize)
                            sample2Type2CountArrayDict[sample][signature] = np.zeros(windowSize, dtype=int)
                            sample2Type2SignalArrayDict[sample][signature] += window_array
                            sample2Type2CountArrayDict[sample][signature] += (window_array > 0)
                ####################################################
        ################# Signatures ends #########################

        ######################################################################
        if my_type in type2SignalArrayDict:
            type2SignalArrayDict[my_type] += window_array
            type2CountArrayDict[my_type] += (window_array > 0)
        else:
            type2SignalArrayDict[my_type] = np.zeros(windowSize)
            type2CountArrayDict[my_type] = np.zeros(windowSize, dtype=int)
            type2SignalArrayDict[my_type] += window_array
            type2CountArrayDict[my_type] += (window_array > 0)

        if sample_based:
            if (sample in sample2NumberofMutationsDict):
                if sample in sample2Type2SignalArrayDict:
                    if my_type in sample2Type2SignalArrayDict[sample]:
                        sample2Type2SignalArrayDict[sample][my_type] += window_array
                        sample2Type2CountArrayDict[sample][my_type] += (window_array > 0)
                    else:
                        sample2Type2SignalArrayDict[sample][my_type] = np.zeros(windowSize)
                        sample2Type2CountArrayDict[sample][my_type] = np.zeros(windowSize, dtype=int)
                        sample2Type2SignalArrayDict[sample][my_type] += window_array
                        sample2Type2CountArrayDict[sample][my_type] += (window_array > 0)
                else:
                    sample2Type2SignalArrayDict[sample] = {}
                    sample2Type2CountArrayDict[sample] = {}
                    sample2Type2SignalArrayDict[sample][my_type] = np.zeros(windowSize)
                    sample2Type2CountArrayDict[sample][my_type] = np.zeros(windowSize, dtype=int)
                    sample2Type2SignalArrayDict[sample][my_type] += window_array
                    sample2Type2CountArrayDict[sample][my_type] += (window_array > 0)
        ######################################################################


########################################################################################



########################################################################################
#Updated JAN 9, 2020
# November 1, 2019
#Called from  imap_unordered
def combined_prepare_chrbased_data_fill_signal_count_arrays_using_inputList(inputList):
    occupancy_type=inputList[0]
    outputDir=inputList[1]
    jobname=inputList[2]
    chrLong=inputList[3]
    simNum=inputList[4]
    chromSizesDict=inputList[5]
    library_file_with_path=inputList[6]
    library_file_type=inputList[7]
    sample2NumberofSubsDict=inputList[8]
    sample2NumberofIndelsDict=inputList[9]
    sample2NumberofDinucsDict=inputList[10]
    sample2SubsSignature2NumberofMutationsDict=inputList[11]
    sample2IndelsSignature2NumberofMutationsDict=inputList[12]
    sample2DinucsSignature2NumberofMutationsDict=inputList[13]
    subsSignature2PropertiesListDict=inputList[14]
    indelsSignature2PropertiesListDict=inputList[15]
    dinucsSignature2PropertiesListDict=inputList[16]
    plusorMinus=inputList[17]
    sample_based=inputList[18]
    verbose=inputList[19]

    return combined_prepare_chrbased_data_fill_signal_count_arrays(occupancy_type,outputDir, jobname, chrLong, simNum, chromSizesDict,
                                        library_file_with_path,
                                        library_file_type,
                                        sample2NumberofSubsDict, sample2NumberofIndelsDict, sample2NumberofDinucsDict,
                                        sample2SubsSignature2NumberofMutationsDict, sample2IndelsSignature2NumberofMutationsDict,sample2DinucsSignature2NumberofMutationsDict,
                                        subsSignature2PropertiesListDict, indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,
                                        plusorMinus, sample_based,verbose)
########################################################################################


########################################################################################
#Using pyBigWig for bigBed and bigWig files starts Optional for unix, linux
#Using chrBasedSignalArrays for big files
#Using dataframes for small bed files
def occupancy_analysis(genome,
                       computation_type,
                        occupancy_type,
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
                        dinucsSignature2PropertiesListDict,
                       verbose):

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
    # If chunksize is 1, maxtasksperchild=x will call the function x times in each process,
    # but if chunksize is y, it will call the function x*y times in each process.
    # Setting maxtaskperchild to 1 would restart each process in your pool after it processed a single task, which is the most aggressive setting you could use to free any leaked resources.
    numofProcesses = multiprocessing.cpu_count()
    # pool = multiprocessing.Pool(numofProcesses, maxtasksperchild=1)
    ##########################################################################

    simNum2Type2AccumulatedSignalArrayDict = {}
    simNum2Type2AccumulatedCountArrayDict = {}
    simNum2Sample2Type2AccumulatedSignalArrayDict = {}
    simNum2Sample2Type2AccumulatedCountArrayDict = {}

    ##############################################################
    #What is the type of the signal_file_with_path?
    #If it is a bed file read signal_file_with_path here
    file_extension = os.path.splitext(os.path.basename(library_file_with_path))[1]

    quantileValue = 0.97
    remove_outliers = False

    if ((file_extension.lower()=='.bigwig') or (file_extension.lower()=='.bw')):
        library_file_type=BIGWIG
        #if chrBasedSignalArrays does not exist we will use pyBigWig if installed and we will not create chrBasedSignalArrays but use BigWig file opened by pyBigWig to fill windowArray
    elif ((file_extension.lower()=='.bigbed') or (file_extension.lower()=='.bb')):
        library_file_type=BIGBED
        #if chrBasedSignalArrays does not exist we will use pyBigWig if installed and we will not create chrBasedSignalArrays but use BigBed file opened by pyBigWig to fill windowArray
    elif (file_extension.lower()=='.bed'):
        library_file_type=BED
        readBEDandWriteChromBasedSignalArrays(outputDir, jobname, genome, library_file_with_path, occupancy_type,quantileValue,remove_outliers)
    elif ((file_extension.lower()=='.narrowpeak') or (file_extension.lower()=='.np')):
        library_file_type=NARROWPEAK
        readBEDandWriteChromBasedSignalArrays(outputDir, jobname, genome, library_file_with_path, occupancy_type,quantileValue,remove_outliers)
    elif (file_extension.lower()=='.wig'):
        library_file_type=WIG
        #For inhouse preparation
        #readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArraysSequentially(genome, quantileValue,library_file_with_path)
        #readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArraysInParallel(genome, quantileValue,library_file_with_path)
        BEDGRAPH=decideFileType(library_file_with_path)
        if BEDGRAPH==True:
            readWig_write_derived_from_bedgraph(outputDir, jobname, genome, library_file_with_path,occupancy_type)
        else:
            readWig_with_fixedStep_variableStep_writeChrBasedSignalArrays(outputDir, jobname, genome, library_file_with_path,occupancy_type,quantileValue,remove_outliers)
    else:
        library_file_type=LIBRARY_FILE_TYPE_OTHER
    ##############################################################


    ###################################################################################
    ##################  USING IMAP UNORDERED starts ###################################
    ###################################################################################
    if (computation_type==USING_IMAP_UNORDERED):

        #original
        sim_nums = range(0, numofSimulations+1)
        sim_num_chr_tuples=((sim_num,chrLong) for sim_num in sim_nums for chrLong in chromNamesList)

        # test
        # sim_nums = range(0,10)
        # chromNamesList=['chr1']
        # sim_num_chr_tuples=((sim_num,chrLong) for sim_num in sim_nums for chrLong in chromNamesList)

        #on tscc-11-9 all works, on tscc-11-17 it does not work
        pool = multiprocessing.Pool(numofProcesses)

        # works in laptop, works in login node in tscc, submits each task to different process, but not hpc on tscc
        # pool = multiprocessing.Pool(numofProcesses, maxtasksperchild=1000)

        # This worked. Sends each job to a new process. Works for laptop, login node,  tscc hpc
        # pool = multiprocessing.Pool(numofProcesses,maxtasksperchild=1)

        #imap uses iterable as input
        # Note that map may cause high memory usage for very long iterables.
        # Consider using imap() or imap_unordered() with explicit chunksize option for better efficiency.
        # This method chops the iterable into a number of chunks which it submits to the process pool as separate tasks.
        # The (approximate) size of these chunks can be specified by setting chunksize to a positive integer.
        # default chunksize=1

        for simulatonBased_SignalArrayAndCountArrayList in pool.imap_unordered(combined_prepare_chrbased_data_fill_signal_count_arrays_using_inputList, (fillInputList(occupancy_type,outputDir,jobname,chrLong,simNum,chromSizesDict,library_file_with_path,library_file_type,
                                                            sample2NumberofSubsDict,sample2NumberofIndelsDict,sample2NumberofDinucsDict,sample2SubsSignature2NumberofMutationsDict,sample2IndelsSignature2NumberofMutationsDict,sample2DinucsSignature2NumberofMutationsDict,
                                                            subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,plusorMinus,sample_based,verbose) for simNum,chrLong in sim_num_chr_tuples),chunksize=1):

                simNum2Type2SignalArrayDict = simulatonBased_SignalArrayAndCountArrayList[0]
                simNum2Type2CountArrayDict = simulatonBased_SignalArrayAndCountArrayList[1]
                simNum2Sample2Type2SignalArrayDict = simulatonBased_SignalArrayAndCountArrayList[2]
                simNum2Sample2Type2CountArrayDict = simulatonBased_SignalArrayAndCountArrayList[3]

                keys = simNum2Type2SignalArrayDict.keys()
                if verbose: print('Accumulate: Worker pid %s current_mem_usage %.2f (mb) simNum:%s' % (str(os.getpid()),memory_usage(),keys))

                accumulateSimulationBasedTypeBasedArrays(simNum2Type2AccumulatedSignalArrayDict, simNum2Type2SignalArrayDict)
                accumulateSimulationBasedTypeBasedArrays(simNum2Type2AccumulatedCountArrayDict, simNum2Type2CountArrayDict)
                accumulateSimulationBasedSampleBasedTypeBasedArrays(simNum2Sample2Type2AccumulatedSignalArrayDict,simNum2Sample2Type2SignalArrayDict)
                accumulateSimulationBasedSampleBasedTypeBasedArrays(simNum2Sample2Type2AccumulatedCountArrayDict,simNum2Sample2Type2CountArrayDict)

        ################################
        pool.close()
        pool.join()
        ################################

    ###################################################################################
    ##################  USING IMAP UNORDERED ends #####################################
    ###################################################################################

    ###################################################################################
    ##################  USING APPLY ASYNC starts ######################################
    ###################################################################################
    elif (computation_type == USING_APPLY_ASYNC):

        sim_nums = range(0, numofSimulations+1)
        sim_num_chr_tuples=((sim_num,chrLong) for sim_num in sim_nums for chrLong in chromNamesList)

        pool = multiprocessing.Pool(numofProcesses)

        #########################################################################################
        def accumulate_apply_async_result(simulatonBased_SignalArrayAndCountArrayList):
            simNum2Type2SignalArrayDict = simulatonBased_SignalArrayAndCountArrayList[0]
            simNum2Type2CountArrayDict = simulatonBased_SignalArrayAndCountArrayList[1]
            simNum2Sample2Type2SignalArrayDict = simulatonBased_SignalArrayAndCountArrayList[2]
            simNum2Sample2Type2CountArrayDict = simulatonBased_SignalArrayAndCountArrayList[3]

            keys = simNum2Type2SignalArrayDict.keys()
            if verbose: print('Worker pid %s memory_usage %.2f MB ACCUMULATE simNum:%s' % (str(os.getpid()), memory_usage(), keys))

            accumulateSimulationBasedTypeBasedArrays(simNum2Type2AccumulatedSignalArrayDict,simNum2Type2SignalArrayDict)
            accumulateSimulationBasedTypeBasedArrays(simNum2Type2AccumulatedCountArrayDict, simNum2Type2CountArrayDict)
            accumulateSimulationBasedSampleBasedTypeBasedArrays(simNum2Sample2Type2AccumulatedSignalArrayDict,simNum2Sample2Type2SignalArrayDict)
            accumulateSimulationBasedSampleBasedTypeBasedArrays(simNum2Sample2Type2AccumulatedCountArrayDict,simNum2Sample2Type2CountArrayDict)
        #########################################################################################

        for simNum, chrLong in sim_num_chr_tuples:

            pool.apply_async(combined_prepare_chrbased_data_fill_signal_count_arrays, (occupancy_type,
                          outputDir, jobname, chrLong, simNum, chromSizesDict,
                          library_file_with_path,
                          library_file_type,
                          sample2NumberofSubsDict, sample2NumberofIndelsDict, sample2NumberofDinucsDict,
                          sample2SubsSignature2NumberofMutationsDict, sample2IndelsSignature2NumberofMutationsDict,
                          sample2DinucsSignature2NumberofMutationsDict,
                          subsSignature2PropertiesListDict, indelsSignature2PropertiesListDict,
                          dinucsSignature2PropertiesListDict, plusorMinus, sample_based,verbose), callback=accumulate_apply_async_result)


        pool.close()
        pool.join()
    ###################################################################################
    ##################  USING APPLY ASYNC ends ########################################
    ###################################################################################

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