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
import traceback

from SigProfilerTopography.source.commons.TopographyCommons import memory_usage
from SigProfilerTopography.source.commons.TopographyCommons import readChrBasedMutationsDF
from SigProfilerTopography.source.commons.TopographyCommons import func_addSignal
from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofSubsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofIndelsDict
from SigProfilerTopography.source.commons.TopographyCommons import getDictionary

from SigProfilerTopography.source.commons.TopographyCommons import getSample2SubsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2IndelsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import writeSimulationBasedAverageNucleosomeOccupancyUsingNumpyArray

from SigProfilerTopography.source.commons.TopographyCommons import TYPE
from SigProfilerTopography.source.commons.TopographyCommons import SUBS
from SigProfilerTopography.source.commons.TopographyCommons import INDELS
from SigProfilerTopography.source.commons.TopographyCommons import DINUCS
from SigProfilerTopography.source.commons.TopographyCommons import MEGABYTE_IN_BYTES

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
from SigProfilerTopography.source.commons.TopographyCommons import BEDGRAPH
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

from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM
from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT

from SigProfilerTopography.source.occupancy.ChrBasedSignalArrays import readBEDandWriteChromBasedSignalArrays
from SigProfilerTopography.source.occupancy.ChrBasedSignalArrays import readWig_with_fixedStep_variableStep_writeChrBasedSignalArrays
from SigProfilerTopography.source.occupancy.ChrBasedSignalArrays import readWig_write_derived_from_bedgraph

from SigProfilerTopography.source.occupancy.ChrBasedSignalArrays import readWig_write_derived_from_bedgraph_using_pool_chunks
from SigProfilerTopography.source.occupancy.ChrBasedSignalArrays import readWig_write_derived_from_bedgraph_using_pool_read_all

from SigProfilerTopography.source.commons.TopographyCommons import decideFileType
from SigProfilerTopography.source.commons.TopographyCommons import get_chrBased_simBased_combined_df_split
from SigProfilerTopography.source.commons.TopographyCommons import get_chrBased_simBased_combined_df
from SigProfilerTopography.source.commons.TopographyCommons import get_chrBased_simBased_dfs

from SigProfilerTopography.source.commons.TopographyCommons import MISSING_SIGNAL

from SigProfilerTopography.source.commons.TopographyCommons import SIGPROFILERTOPOGRAPHY_DEFAULT_FILES
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_ATAC_SEQ_OCCUPANCY_FILE

########################################################################################
# April 27, 2020
# requires chrBased_simBased_combined_df_split which can be real split or whole in fact
# This is common for pool.imap_unordered and pool.apply_async variations
def chrbased_data_fill_signal_count_arrays_for_all_mutations(occupancy_type,
                                                             occupancy_calculation_type,
                                                             outputDir,
                                                             jobname,
                                                             chrLong,
                                                             simNum,
                                                             chrBased_simBased_subs_df,
                                                             chrBased_simBased_dinucs_df,
                                                             chrBased_simBased_indels_df,
                                                             chromSizesDict,
                                                             library_file_with_path,
                                                             library_file_type,
                                                             ordered_sbs_signatures,
                                                             ordered_dbs_signatures,
                                                             ordered_id_signatures,
                                                             ordered_sbs_signatures_cutoffs,
                                                             ordered_dbs_signatures_cutoffs,
                                                             ordered_id_signatures_cutoffs,
                                                             plusorMinus,
                                                             verbose):

    if verbose: print('\tVerbose %s Worker pid %s memory_usage in %.2f MB START chrLong:%s simNum:%d' %(occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum))
    # 1st part  Prepare chr based mutations dataframes
    maximum_chrom_size = chromSizesDict[chrLong]
    start_time = time.time()

    ##############################################################
    chrBasedSignalArray = None #Will be filled from chrBasedSignal files if they exists
    library_file_opened_by_pyBigWig = None #Will be filled by pyBigWig from bigWig or bigBed
    my_upperBound = None
    signal_index = None
    ##############################################################

    if verbose: print('\tVerbose %s Worker pid %s memory_usage in %.2f MB Check1 Read Signal Array and Dataframes chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum))
    # if verbose: print('\tVerbose %s Worker pid %s -- signal_array_npy: %f in MB -- chrBased_simBased_combined_df_split: %f in MB -- chrLong:%s simNum:%d' % (
    #         occupancy_type,
    #         str(os.getpid()),
    #         sys.getsizeof(chrBasedSignalArray) / MEGABYTE_IN_BYTES,
    #         sys.getsizeof(chrBased_simBased_combined_df_split) / MEGABYTE_IN_BYTES,
    #         chrLong, simNum))
    #################################################################################################################

    #################################################################################################################
    libraryFilenameWoExtension = os.path.splitext(os.path.basename(library_file_with_path))[0]
    signalArrayFilename = '%s_signal_%s.npy' % (chrLong, libraryFilenameWoExtension)
    if (occupancy_type==NUCLEOSOMEOCCUPANCY):
        chrBasedSignalFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME, CHRBASED, signalArrayFilename)
    elif (occupancy_type== EPIGENOMICSOCCUPANCY):
        if os.path.basename(library_file_with_path) == DEFAULT_ATAC_SEQ_OCCUPANCY_FILE:
            chrBasedSignalFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, EPIGENOMICS, CHRBASED, signalArrayFilename)
        else:
            chrBasedSignalFile = os.path.join(outputDir,jobname,DATA,occupancy_type,LIB,CHRBASED,signalArrayFilename)
    else:
        #It can be EPIGENOMICSOCCUPANCY or user provided name e.g.: Epigenomics_ATAC_ENCFF317TWD
        chrBasedSignalFile = os.path.join(outputDir,jobname,DATA,occupancy_type,LIB,CHRBASED,signalArrayFilename)

    #Downloaded or created runtime
    if (os.path.exists(chrBasedSignalFile)):
        #Can this cause to deep sleep of processes?
        # chrBasedSignalArray = np.load(chrBasedSignalFile, mmap_mode='r')
        chrBasedSignalArray = np.load(chrBasedSignalFile)

    #If library_file_with_path is abs path and library_file_type is BIGWIG or BIGBED
    #For nucleosome_biosample==GM12878 or nucleosome_biosample==K562 library_file_with_path is only filename with extension, it is not absolute path
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
                    my_upperBound = my_mean + std_dev * 3
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

    ###############################################################################
    ################################ Initialization ###############################
    ###############################################################################
    # Add one more row for the aggregated analysis
    subsSignature_accumulated_signal_np_array = np.zeros((ordered_sbs_signatures.size + 1, plusorMinus * 2 + 1))
    dinucsSignature_accumulated_signal_np_array = np.zeros((ordered_dbs_signatures.size + 1, plusorMinus * 2 + 1))
    indelsSignature_accumulated_signal_np_array = np.zeros((ordered_id_signatures.size + 1, plusorMinus * 2 + 1))

    # Add one more row for the aggregated analysis
    subsSignature_accumulated_count_np_array = np.zeros((ordered_sbs_signatures.size + 1, plusorMinus * 2 + 1))
    dinucsSignature_accumulated_count_np_array = np.zeros((ordered_dbs_signatures.size + 1, plusorMinus * 2 + 1))
    indelsSignature_accumulated_count_np_array = np.zeros((ordered_id_signatures.size + 1, plusorMinus * 2 + 1))
    ###############################################################################
    ################################ Initialization ###############################
    ###############################################################################

    #################################################################################################################
    if ((chrBasedSignalArray is not None) or ((library_file_opened_by_pyBigWig is not None) and (chrLong in library_file_opened_by_pyBigWig.chroms()))):
        ######################################################## #######################
        ################### Fill signal and count array starts ########################
        ###############################################################################
        if verbose: print('\tVerbose %s Worker pid %s memory_usage in %.2f MB Check2_1 Start chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum))

        #For Subs
        if ((chrBased_simBased_subs_df is not None) and (not chrBased_simBased_subs_df.empty)):

            #df_columns is a numpy array
            df_columns = chrBased_simBased_subs_df.columns.values
            subsSignatures_mask_array = np.isin(df_columns,ordered_sbs_signatures)

            #July 25, 2020
            [fillSignalArrayAndCountArray_using_list_comp(
                row,
                chrLong,
                library_file_opened_by_pyBigWig,
                chrBasedSignalArray,
                library_file_type,
                signal_index,
                my_upperBound,
                maximum_chrom_size,
                ordered_sbs_signatures_cutoffs,
                subsSignatures_mask_array,
                subsSignature_accumulated_signal_np_array,
                subsSignature_accumulated_count_np_array,
                plusorMinus,
                df_columns,
                occupancy_calculation_type) for row in chrBased_simBased_subs_df[df_columns].values]

            # print('DEBUG SUBS chrLong:%s simNum:%d ENDED' %(chrLong,simNum),flush=True)

        #For Dinusc
        if ((chrBased_simBased_dinucs_df is not None) and (not chrBased_simBased_dinucs_df.empty)):

            #df_columns is a numpy array
            df_columns = chrBased_simBased_dinucs_df.columns.values
            dinucsSignatures_mask_array = np.isin(df_columns,ordered_dbs_signatures)

            #July 25, 2020
            [fillSignalArrayAndCountArray_using_list_comp(
                row,
                chrLong,
                library_file_opened_by_pyBigWig,
                chrBasedSignalArray,
                library_file_type,
                signal_index,
                my_upperBound,
                maximum_chrom_size,
                ordered_dbs_signatures_cutoffs,
                dinucsSignatures_mask_array,
                dinucsSignature_accumulated_signal_np_array,
                dinucsSignature_accumulated_count_np_array,
                plusorMinus,
                df_columns,
                occupancy_calculation_type) for row in chrBased_simBased_dinucs_df[df_columns].values]

            # print('DEBUG DINUCS chrLong:%s simNum:%d ENDED' %(chrLong,simNum),flush=True)

        #For Indels
        if ((chrBased_simBased_indels_df is not None) and (not chrBased_simBased_indels_df.empty)):

            #df_columns is a numpy array
            df_columns = chrBased_simBased_indels_df.columns.values
            indelsSignatures_mask_array = np.isin(df_columns,ordered_id_signatures)

            #July 25, 2020
            [fillSignalArrayAndCountArray_using_list_comp(
                row,
                chrLong,
                library_file_opened_by_pyBigWig,
                chrBasedSignalArray,
                library_file_type,
                signal_index,
                my_upperBound,
                maximum_chrom_size,
                ordered_id_signatures_cutoffs,
                indelsSignatures_mask_array,
                indelsSignature_accumulated_signal_np_array,
                indelsSignature_accumulated_count_np_array,
                plusorMinus,
                df_columns,
                occupancy_calculation_type) for row in chrBased_simBased_indels_df[df_columns].values]

            # print('DEBUG INDELS chrLong:%s simNum:%d ENDED' %(chrLong,simNum),flush=True)

        if verbose: print('\tVerbose %s Worker pid %s memory_usage in %.2f MB Check2_2 End chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum))
        ###############################################################################
        ################### Fill signal and count array ends ##########################
        ###############################################################################

    if (library_file_opened_by_pyBigWig is not None):
        library_file_opened_by_pyBigWig.close()

    if verbose: print('\tVerbose %s Worker pid %s memory_usage in %.2f MB END  chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum))
    if verbose: print('\tVerbose %s Worker pid %s took %f seconds chrLong:%s simNum:%d\n' % (occupancy_type,str(os.getpid()), (time.time() - start_time), chrLong, simNum))
    ###############################################################################
    ################### Return  starts ############################################
    ###############################################################################

    # Initialzie the list, you will return this list
    SignalArrayAndCountArrayList = []

    #new way
    SignalArrayAndCountArrayList.append(chrLong)
    SignalArrayAndCountArrayList.append(simNum)
    SignalArrayAndCountArrayList.append(subsSignature_accumulated_signal_np_array)
    SignalArrayAndCountArrayList.append(dinucsSignature_accumulated_signal_np_array)
    SignalArrayAndCountArrayList.append(indelsSignature_accumulated_signal_np_array)
    SignalArrayAndCountArrayList.append(subsSignature_accumulated_count_np_array)
    SignalArrayAndCountArrayList.append(dinucsSignature_accumulated_count_np_array)
    SignalArrayAndCountArrayList.append(indelsSignature_accumulated_count_np_array)

    return SignalArrayAndCountArrayList
########################################################################################


########################################################################################
# May 5, 2020
# For apply_async
# Read chromBased simBased combined mutations df in the process
def chrbased_data_fill_signal_count_arrays_for_all_mutations_read_mutations(
        # input_list):
        occupancy_type,
        occupancy_calculation_type,
        outputDir,
        jobname,
        chrLong,
        simNum,
        chromSizesDict,
        library_file_with_path,
        library_file_type,
        ordered_sbs_signatures,
        ordered_dbs_signatures,
        ordered_id_signatures,
        ordered_sbs_signatures_cutoffs,
        ordered_dbs_signatures_cutoffs,
        ordered_id_signatures_cutoffs,
        plusorMinus,
        verbose):

    # occupancy_type=input_list[0]
    # occupancy_calculation_type=input_list[1]
    # outputDir=input_list[2]
    # jobname=input_list[3]
    # chrLong=input_list[4]
    # simNum=input_list[5]
    # chromSizesDict=input_list[6]
    # library_file_with_path=input_list[7]
    # library_file_type=input_list[8]
    # ordered_sbs_signatures=input_list[9]
    # ordered_dbs_signatures=input_list[10]
    # ordered_id_signatures=input_list[11]
    # ordered_sbs_signatures_cutoffs=input_list[12]
    # ordered_dbs_signatures_cutoffs=input_list[13]
    # ordered_id_signatures_cutoffs=input_list[14]
    # plusorMinus=input_list[15]
    # verbose=input_list[16]

    try:
        # chrBased_simBased_mutations_df = get_chrBased_simBased_combined_df(outputDir, jobname, chrLong, simNum)
        chrBased_simBased_subs_df, chrBased_simBased_dinucs_df, chrBased_simBased_indels_df = get_chrBased_simBased_dfs(outputDir, jobname, chrLong, simNum)

        return chrbased_data_fill_signal_count_arrays_for_all_mutations(occupancy_type,
                                                                        occupancy_calculation_type,
                                                                        outputDir,
                                                                        jobname,
                                                                        chrLong,
                                                                        simNum,
                                                                        chrBased_simBased_subs_df,
                                                                        chrBased_simBased_dinucs_df,
                                                                        chrBased_simBased_indels_df,
                                                                        chromSizesDict,
                                                                        library_file_with_path,
                                                                        library_file_type,
                                                                        ordered_sbs_signatures,
                                                                        ordered_dbs_signatures,
                                                                        ordered_id_signatures,
                                                                        ordered_sbs_signatures_cutoffs,
                                                                        ordered_dbs_signatures_cutoffs,
                                                                        ordered_id_signatures_cutoffs,
                                                                        plusorMinus,
                                                                        verbose)
    except Exception as e:
        print("Exception: %s" % (e))
        # print(e)
        # traceback.print_exc()
        # raise e
########################################################################################

########################################################################################
#May 19, 2020
# For apply_async split using poolInputList
# Read chromBased simBased combined mutations df split in the process
def chrbased_data_fill_signal_count_arrays_for_all_mutations_read_mutations_split(occupancy_type,
                                                                                  occupancy_calculation_type,
                                                                                  outputDir,
                                                                                  jobname,
                                                                                  chrLong,
                                                                                  simNum,
                                                                                  splitIndex,
                                                                                  chromSizesDict,
                                                                                  library_file_with_path,
                                                                                  library_file_type,
                                                                                  ordered_sbs_signatures,
                                                                                  ordered_dbs_signatures,
                                                                                  ordered_id_signatures,
                                                                                  ordered_sbs_signatures_cutoffs,
                                                                                  ordered_dbs_signatures_cutoffs,
                                                                                  ordered_id_signatures_cutoffs,
                                                                                  plusorMinus,
                                                                                  verbose):

    chrBased_simBased_combined_df_split = get_chrBased_simBased_combined_df_split(outputDir, jobname, chrLong, simNum,splitIndex)

    return chrbased_data_fill_signal_count_arrays_for_all_mutations_for_df_split(occupancy_type,
                                                                    occupancy_calculation_type,
                                                                    outputDir,
                                                                    jobname,
                                                                    chrLong,
                                                                    simNum,
                                                                    chrBased_simBased_combined_df_split,
                                                                    chromSizesDict,
                                                                    library_file_with_path,
                                                                    library_file_type,
                                                                    ordered_sbs_signatures,
                                                                    ordered_dbs_signatures,
                                                                    ordered_id_signatures,
                                                                    ordered_sbs_signatures_cutoffs,
                                                                    ordered_dbs_signatures_cutoffs,
                                                                    ordered_id_signatures_cutoffs,
                                                                    plusorMinus,
                                                                    verbose)
########################################################################################


########################################################################################
# For df_split
# April 27, 2020
# requires chrBased_simBased_combined_df_split which can be real split or whole in fact
# This is common for pool.imap_unordered and pool.apply_async variations
def chrbased_data_fill_signal_count_arrays_for_all_mutations_for_df_split(occupancy_type,
                                                             occupancy_calculation_type,
                                                             outputDir,
                                                             jobname,
                                                             chrLong,
                                                             simNum,
                                                             chrBased_simBased_combined_df_split,
                                                             chromSizesDict,
                                                             library_file_with_path,
                                                             library_file_type,
                                                             ordered_sbs_signatures,
                                                             ordered_dbs_signatures,
                                                             ordered_id_signatures,
                                                             ordered_sbs_signatures_cutoffs,
                                                             ordered_dbs_signatures_cutoffs,
                                                             ordered_id_signatures_cutoffs,
                                                             plusorMinus,
                                                             verbose):

    if verbose: print('\tVerbose %s Worker pid %s memory_usage in %.2f MB START chrLong:%s simNum:%d' %(occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum))
    # 1st part  Prepare chr based mutations dataframes
    maximum_chrom_size = chromSizesDict[chrLong]
    start_time = time.time()

    ##############################################################
    chrBasedSignalArray = None #Will be filled from chrBasedSignal files if they exists
    library_file_opened_by_pyBigWig = None #Will be filled by pyBigWig from bigWig or bigBed
    my_upperBound = None
    signal_index = None
    ##############################################################

    if (chrBased_simBased_combined_df_split is not None) and verbose:
        print('\tVerbose %s Worker pid %s chrBased_mutations_df(%d,%d) ' %(occupancy_type,str(os.getpid()),chrBased_simBased_combined_df_split.shape[0],chrBased_simBased_combined_df_split.shape[1]))

    if verbose: print('\tVerbose %s Worker pid %s memory_usage in %.2f MB Check1 Read Signal Array and Dataframes chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum))
    if verbose: print('\tVerbose %s Worker pid %s -- signal_array_npy: %f in MB -- chrBased_simBased_combined_df_split: %f in MB -- chrLong:%s simNum:%d' % (
            occupancy_type,
            str(os.getpid()),
            sys.getsizeof(chrBasedSignalArray) / MEGABYTE_IN_BYTES,
            sys.getsizeof(chrBased_simBased_combined_df_split) / MEGABYTE_IN_BYTES,
            chrLong, simNum))
    #################################################################################################################

    #################################################################################################################
    libraryFilenameWoExtension = os.path.splitext(os.path.basename(library_file_with_path))[0]
    signalArrayFilename = '%s_signal_%s.npy' % (chrLong, libraryFilenameWoExtension)
    if (occupancy_type==NUCLEOSOMEOCCUPANCY):
        chrBasedSignalFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME, CHRBASED,signalArrayFilename)
    elif (occupancy_type== EPIGENOMICSOCCUPANCY):
        chrBasedSignalFile = os.path.join(outputDir,jobname,DATA,occupancy_type,LIB,CHRBASED,signalArrayFilename)
    else:
        #It can be EPIGENOMICSOCCUPANCY or user provided name e.g.: Epigenomics_ATAC_ENCFF317TWD
        chrBasedSignalFile = os.path.join(outputDir,jobname,DATA,occupancy_type,LIB,CHRBASED,signalArrayFilename)

    #Downloaded or created runtime
    if (os.path.exists(chrBasedSignalFile)):
        #Can this cause to deep sleep of processes?
        # chrBasedSignalArray = np.load(chrBasedSignalFile, mmap_mode='r')
        chrBasedSignalArray = np.load(chrBasedSignalFile)

    #If library_file_with_path is abs path and library_file_type is BIGWIG or BIGBED
    #For nucleosome_biosample==GM12878 or nucleosome_biosample==K562 library_file_with_path is only filename with extension, it is not absolute path
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
                    my_upperBound = my_mean + std_dev * 3
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
        if verbose: print('\tVerbose %s Worker pid %s memory_usage in %.2f MB Check2_1 Start chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum))
        if ((chrBased_simBased_combined_df_split is not None) and (not chrBased_simBased_combined_df_split.empty)):

            #df_columns is a numpy array
            df_columns = chrBased_simBased_combined_df_split.columns.values

            ###############################################################################
            ################################ Initialization ###############################
            ###############################################################################

            subsSignatures_mask_array = np.isin(df_columns,ordered_sbs_signatures)
            dinucsSignatures_mask_array = np.isin(df_columns,ordered_dbs_signatures)
            indelsSignatures_mask_array = np.isin(df_columns,ordered_id_signatures)

            #Add one more row for the aggregated analysis
            subsSignature_accumulated_signal_np_array=np.zeros((ordered_sbs_signatures.size+1,plusorMinus*2+1))
            dinucsSignature_accumulated_signal_np_array=np.zeros((ordered_dbs_signatures.size+1,plusorMinus*2+1))
            indelsSignature_accumulated_signal_np_array=np.zeros((ordered_id_signatures.size+1,plusorMinus*2+1))

            #Add one more row for the aggregated analysis
            subsSignature_accumulated_count_np_array=np.zeros((ordered_sbs_signatures.size+1,plusorMinus*2+1))
            dinucsSignature_accumulated_count_np_array=np.zeros((ordered_dbs_signatures.size+1,plusorMinus*2+1))
            indelsSignature_accumulated_count_np_array=np.zeros((ordered_id_signatures.size+1,plusorMinus*2+1))
            ###############################################################################
            ################################ Initialization ###############################
            ###############################################################################

            #July 25, 2020
            [fillSignalArrayAndCountArray_using_list_comp_for_df_split(
                row,
                chrLong,
                library_file_opened_by_pyBigWig,
                chrBasedSignalArray,
                library_file_type,
                signal_index,
                my_upperBound,
                maximum_chrom_size,
                ordered_sbs_signatures_cutoffs,
                ordered_dbs_signatures_cutoffs,
                ordered_id_signatures_cutoffs,
                subsSignatures_mask_array,
                dinucsSignatures_mask_array,
                indelsSignatures_mask_array,
                subsSignature_accumulated_signal_np_array,
                dinucsSignature_accumulated_signal_np_array,
                indelsSignature_accumulated_signal_np_array,
                subsSignature_accumulated_count_np_array,
                dinucsSignature_accumulated_count_np_array,
                indelsSignature_accumulated_count_np_array,
                plusorMinus,
                df_columns,
                occupancy_calculation_type) for row in chrBased_simBased_combined_df_split[df_columns].values]


        if verbose: print('\tVerbose %s Worker pid %s memory_usage in %.2f MB Check2_2 End chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum))
        ###############################################################################
        ################### Fill signal and count array ends ##########################
        ###############################################################################

    if (library_file_opened_by_pyBigWig is not None):
        library_file_opened_by_pyBigWig.close()

    if verbose: print('\tVerbose %s Worker pid %s memory_usage in %.2f MB END  chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum))
    if verbose: print('\tVerbose %s Worker pid %s took %f seconds chrLong:%s simNum:%d\n' % (occupancy_type,str(os.getpid()), (time.time() - start_time), chrLong, simNum))
    ###############################################################################
    ################### Return  starts ############################################
    ###############################################################################

    # Initialzie the list, you will return this list
    SignalArrayAndCountArrayList = []

    #new way
    SignalArrayAndCountArrayList.append(simNum)
    SignalArrayAndCountArrayList.append(subsSignature_accumulated_signal_np_array)
    SignalArrayAndCountArrayList.append(dinucsSignature_accumulated_signal_np_array)
    SignalArrayAndCountArrayList.append(indelsSignature_accumulated_signal_np_array)
    SignalArrayAndCountArrayList.append(subsSignature_accumulated_count_np_array)
    SignalArrayAndCountArrayList.append(dinucsSignature_accumulated_count_np_array)
    SignalArrayAndCountArrayList.append(indelsSignature_accumulated_count_np_array)

    return SignalArrayAndCountArrayList
########################################################################################


########################################################################################
#For df_split
#July 25, 2020, Vectorization
def fillSignalArrayAndCountArray_using_list_comp_for_df_split(
        row,
        chrLong,
        library_file_opened_by_pyBigWig,
        chrBasedSignalArray,
        library_file_type,
        signal_index,
        my_upperBound,
        maximum_chrom_size,
        ordered_sbs_signatures_cutoffs,
        ordered_dbs_signatures_cutoffs,
        ordered_id_signatures_cutoffs,
        subsSignatures_mask_array,
        dinucsSignatures_mask_array,
        indelsSignatures_mask_array,
        subsSignature_accumulated_signal_np_array,
        dinucsSignature_accumulated_signal_np_array,
        indelsSignature_accumulated_signal_np_array,
        subsSignature_accumulated_count_np_array,
        dinucsSignature_accumulated_count_np_array,
        indelsSignature_accumulated_count_np_array,
        plusOrMinus,
        df_columns,
        occupancy_calculation_type):

    indexofType = np.where(df_columns == TYPE)[0][0]
    indexofStart = np.where(df_columns == START)[0][0]

    mutation_row_type = row[indexofType]
    mutation_row_start = row[indexofStart]

    ###########################################
    if mutation_row_type == SUBS:
        accumulated_signal_np_array=subsSignature_accumulated_signal_np_array
        accumulated_count_np_array=subsSignature_accumulated_count_np_array
        cutoffs=ordered_sbs_signatures_cutoffs
        signatures_mask_array=subsSignatures_mask_array
    elif mutation_row_type == DINUCS:
        accumulated_signal_np_array = dinucsSignature_accumulated_signal_np_array
        accumulated_count_np_array = dinucsSignature_accumulated_count_np_array
        cutoffs=ordered_dbs_signatures_cutoffs
        signatures_mask_array=dinucsSignatures_mask_array
    elif mutation_row_type == INDELS:
        accumulated_signal_np_array=indelsSignature_accumulated_signal_np_array
        accumulated_count_np_array=indelsSignature_accumulated_count_np_array
        cutoffs=ordered_id_signatures_cutoffs
        signatures_mask_array=indelsSignatures_mask_array
    ###########################################

    window_array=None
    windowSize=plusOrMinus*2+1

    # df_columns 'numpy.ndarray'
    # df_columns: ['Sample', 'Chrom', 'Start', 'MutationLong', 'PyramidineStrand', 'TranscriptionStrand', 'Mutation',
    #              'SBS1', 'SBS2', 'SBS3', 'SBS4', 'SBS5', 'SBS6', 'SBS7a', 'SBS7b', 'SBS7c', 'SBS7d', 'SBS8', 'SBS9',
    #              'SBS10a', 'SBS10b', 'SBS11', 'SBS12', 'SBS13', 'SBS14', 'SBS15', 'SBS16', 'SBS17a', 'SBS17b', 'SBS18',
    #              'SBS19', 'SBS20', 'SBS21', 'SBS22', 'SBS23', 'SBS24', 'SBS25', 'SBS26', 'SBS27', 'SBS28', 'SBS29',
    #              'SBS30', 'SBS31', 'SBS32', 'SBS33', 'SBS34', 'SBS35', 'SBS36', 'SBS37', 'SBS38', 'SBS39', 'SBS40',
    #              'SBS41', 'SBS42', 'SBS43', 'SBS44', 'SBS45', 'SBS46', 'SBS47', 'SBS48', 'SBS49', 'SBS50', 'SBS51',
    #              'SBS52', 'SBS53', 'SBS54', 'SBS55', 'SBS56', 'SBS57', 'SBS58', 'SBS59', 'SBS60', 'Simulation_Number',
    #              'Type', 'Ref', 'Alt', 'Length', 'ID1', 'ID2', 'ID3', 'ID4', 'ID5', 'ID6', 'ID7', 'ID8', 'ID9', 'ID10',
    #              'ID11', 'ID12', 'ID13', 'ID14', 'ID15', 'ID16', 'ID17', 'DBS1', 'DBS2', 'DBS3', 'DBS4', 'DBS5', 'DBS6',
    #              'DBS7', 'DBS8', 'DBS9', 'DBS10', 'DBS11']

    #Get or fill window_array using Case1, Case2, and Case3
    # Case 1: start is very close to the chromosome start
    if (mutation_row_start<plusOrMinus):
        # print('Case 1: start is very close to the chromosome start --- mutation[Start]:%d' %(mutation_row_start))
        #Faster
        if (chrBasedSignalArray is not None):
            window_array = chrBasedSignalArray[0:(mutation_row_start + plusOrMinus + 1)]
            window_array = np.pad(window_array, (plusOrMinus - mutation_row_start, 0), 'constant', constant_values=(0, 0))

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

    ##########################################################
    #Get the sample at this mutation_row
    # sample = mutation_row_sample
    # simulationNumber= mutation_row_simulation_number

    #September 18, 2020 NO SIGNAL caseis added
    #Vectorize July 25, 2020
    #Fill numpy arrays using window_array
    if (window_array is not None) and (np.any(window_array)):
        probabilities = row[signatures_mask_array]
        threshold_mask_array = np.greater_equal(probabilities, cutoffs)

        #Convert True into 1, and False into 0
        mask_array = threshold_mask_array.astype(int)

        #Add 1 for the Aggregated analysis to the mask array
        mask_array = np.append(mask_array, 1)

        #Add one more dimension to window_array and mask_array
        window_array_1x2001=np.array([window_array])
        mask_array_1xnumofsignatures=np.array([mask_array])

        to_be_accumulated_array = mask_array_1xnumofsignatures.T * window_array_1x2001
        accumulated_signal_np_array += to_be_accumulated_array

        #default
        if occupancy_calculation_type==MISSING_SIGNAL:
            accumulated_count_np_array += (to_be_accumulated_array>0)
        else:
            accumulated_count_np_array += 1
    ##########################################################

########################################################################################


########################################################################################
#Dec 17, 2020
def get_window_array(mutation_row_start,
                     plusOrMinus,
                     chrLong,
                     chrBasedSignalArray,
                     library_file_type,
                     library_file_opened_by_pyBigWig,
                     signal_index,
                     my_upperBound,
                     maximum_chrom_size):

    window_array=None
    windowSize=plusOrMinus*2+1

    # df_columns 'numpy.ndarray'
    # df_columns: ['Sample', 'Chrom', 'Start', 'MutationLong', 'PyramidineStrand', 'TranscriptionStrand', 'Mutation',
    #              'SBS1', 'SBS2', 'SBS3', 'SBS4', 'SBS5', 'SBS6', 'SBS7a', 'SBS7b', 'SBS7c', 'SBS7d', 'SBS8', 'SBS9',
    #              'SBS10a', 'SBS10b', 'SBS11', 'SBS12', 'SBS13', 'SBS14', 'SBS15', 'SBS16', 'SBS17a', 'SBS17b', 'SBS18',
    #              'SBS19', 'SBS20', 'SBS21', 'SBS22', 'SBS23', 'SBS24', 'SBS25', 'SBS26', 'SBS27', 'SBS28', 'SBS29',
    #              'SBS30', 'SBS31', 'SBS32', 'SBS33', 'SBS34', 'SBS35', 'SBS36', 'SBS37', 'SBS38', 'SBS39', 'SBS40',
    #              'SBS41', 'SBS42', 'SBS43', 'SBS44', 'SBS45', 'SBS46', 'SBS47', 'SBS48', 'SBS49', 'SBS50', 'SBS51',
    #              'SBS52', 'SBS53', 'SBS54', 'SBS55', 'SBS56', 'SBS57', 'SBS58', 'SBS59', 'SBS60', 'Simulation_Number',
    #              'Type', 'Ref', 'Alt', 'Length', 'ID1', 'ID2', 'ID3', 'ID4', 'ID5', 'ID6', 'ID7', 'ID8', 'ID9', 'ID10',
    #              'ID11', 'ID12', 'ID13', 'ID14', 'ID15', 'ID16', 'ID17', 'DBS1', 'DBS2', 'DBS3', 'DBS4', 'DBS5', 'DBS6',
    #              'DBS7', 'DBS8', 'DBS9', 'DBS10', 'DBS11']

    #Get or fill window_array using Case1, Case2, and Case3
    # Case 1: start is very close to the chromosome start
    if (mutation_row_start<plusOrMinus):
        # print('Case 1: start is very close to the chromosome start --- mutation[Start]:%d' %(mutation_row_start))
        #Faster
        if (chrBasedSignalArray is not None):
            window_array = chrBasedSignalArray[0:(mutation_row_start + plusOrMinus + 1)]
            window_array = np.pad(window_array, (plusOrMinus - mutation_row_start, 0), 'constant', constant_values=(0, 0))

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

    return window_array
########################################################################################

########################################################################################
def accumulate_arrays(row,
                      window_array,
                      signatures_mask_array,
                      ordered_signatures_cutoffs,
                      occupancy_calculation_type,
                      accumulated_signal_np_array,
                      accumulated_count_np_array):

    ##########################################################
    #September 18, 2020 NO SIGNAL case is added
    #Vectorize July 25, 2020
    #Fill numpy arrays using window_array
    if (window_array is not None) and (np.any(window_array)):
        probabilities = row[signatures_mask_array]
        threshold_mask_array = np.greater_equal(probabilities, ordered_signatures_cutoffs)

        #Convert True into 1, and False into 0
        mask_array = threshold_mask_array.astype(int)

        #Add 1 for the Aggregated analysis to the mask array
        mask_array = np.append(mask_array, 1)

        #Add one more dimension to window_array and mask_array
        window_array_1x2001=np.array([window_array])
        mask_array_1xnumofsignatures=np.array([mask_array])

        to_be_accumulated_array = mask_array_1xnumofsignatures.T * window_array_1x2001
        accumulated_signal_np_array += to_be_accumulated_array

        #default
        if occupancy_calculation_type==MISSING_SIGNAL:
            accumulated_count_np_array += (to_be_accumulated_array>0)
        else:
            accumulated_count_np_array += 1
    ##########################################################

########################################################################################

########################################################################################
#July 25, 2020, Vectorization
def fillSignalArrayAndCountArray_using_list_comp(
        row,
        chrLong,
        library_file_opened_by_pyBigWig,
        chrBasedSignalArray,
        library_file_type,
        signal_index,
        my_upperBound,
        maximum_chrom_size,
        ordered_signatures_cutoffs,
        signatures_mask_array,
        accumulated_signal_np_array,
        accumulated_count_np_array,
        plusOrMinus,
        df_columns,
        occupancy_calculation_type):

    indexofStart = np.where(df_columns == START)[0][0]
    mutation_row_start = row[indexofStart]

    window_array=get_window_array(mutation_row_start,
                     plusOrMinus,
                     chrLong,
                     chrBasedSignalArray,
                     library_file_type,
                     library_file_opened_by_pyBigWig,
                     signal_index,
                     my_upperBound,
                     maximum_chrom_size)

    accumulate_arrays(row,
                      window_array,
                      signatures_mask_array,
                      ordered_signatures_cutoffs,
                      occupancy_calculation_type,
                      accumulated_signal_np_array,
                      accumulated_count_np_array)
########################################################################################


#######################################################
#October 13, 2020
def check_download_chrbased_npy_atac_seq_files(outputDir,jobname,occupancy_type,atac_seq_file,chromNamesList):
    current_abs_path = os.path.dirname(os.path.abspath(__file__))
    # print(current_abs_path)

    os.makedirs(os.path.join(outputDir,jobname,DATA,occupancy_type,LIB,CHRBASED),exist_ok=True)
    chrombased_npy_path = os.path.join(outputDir,jobname,DATA,occupancy_type,LIB,CHRBASED)
    # print(chrombased_npy_path)

    if os.path.isabs(chrombased_npy_path):
        # print('%s an absolute path.' %(chrombased_npy_path))
        os.chdir(chrombased_npy_path)

        atac_seq_filename_wo_extension = os.path.splitext(os.path.basename(atac_seq_file))[0]

        for chrLong in chromNamesList:
            filename = '%s_signal_%s.npy' % (chrLong, atac_seq_filename_wo_extension)

            chrbased_npy_array_path = os.path.join(chrombased_npy_path, filename)
            if not os.path.exists(chrbased_npy_array_path):
                print('Does not exists: %s' % (chrbased_npy_array_path))
                try:
                    print('Downloading %s under %s' % (filename, chrbased_npy_array_path))

                    # wget -c Continue getting a partially-downloaded file
                    # wget -nc  If a file is downloaded more than once in the same directory, the local file will be clobbered, or overwritten
                    # cmd="bash -c '" + 'wget -r -l1 -c -nc --no-parent -nd -P ' + chrombased_npy_path + ' ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerTopography/lib/nucleosome/chrbased/' + filename + "'"

                    # -r When included, the wget will recursively traverse subdirectories in order to obtain all content.
                    # -l1 Limit recursion depth to a specific number of levels, by setting the <#> variable to the desired number.
                    # -c option to resume a download
                    # -nc, --no-clobber If a file is downloaded more than once in the same directory, Wget's behavior depends on a few options, including -nc.  In certain cases, the local file will be clobbered, or overwritten, upon repeated download.  In other cases it will be preserved.
                    # -np, --no-parent Do not ever ascend to the parent directory when retrieving recursively.  This is a useful option, since it guarantees that only the files below a certain hierarchy will be downloaded.
                    # -nd, --no-directories When included, directories will not be created. All files captured in the wget will be copied directly in to the active directory
                    cmd = "bash -c '" + 'wget -r -l1 -c -nc --no-parent -nd ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerTopography/lib/epigenomics/chrbased/' + filename + "'"
                    print("cmd: %s" %cmd)
                    os.system(cmd)
                except:
                    # print("The UCSD ftp site is not responding...pulling from sanger ftp now.")
                    print("The UCSD ftp site is not responding...")

    else:
        #It has to be an absolute path
        print('%s is not an absolute path.' %(chrombased_npy_path))

    #go back
    os.chdir(current_abs_path)
#######################################################



########################################################################################
#main function
#Using pyBigWig for bigBed and bigWig files starts Optional for unix, linux
#Using chrBasedSignalArrays for big files
#Using dataframes for small bed files
def occupancyAnalysis(genome,
                       computation_type,
                        occupancy_type,
                        occupancy_calculation_type,
                        sample_based,
                        plusorMinus,
                        chromSizesDict,
                        chromNamesList,
                        outputDir,
                        jobname,
                        numofSimulations,
                        job_tuples,
                        library_file_with_path,
                        library_file_memo,
                        ordered_sbs_signatures,
                        ordered_dbs_signatures,
                        ordered_id_signatures,
                        ordered_sbs_signatures_cutoffs,
                        ordered_dbs_signatures_cutoffs,
                        ordered_id_signatures_cutoffs,
                        remove_outliers,
                        quantileValue,
                        verbose):

    print('\n#################################################################################')
    print('--- %s Analysis starts' %(occupancy_type))
    print('--- Computation Type:%s' % (computation_type))
    print('--- Occupancy Type:%s' % (occupancy_type))
    print('--- Library file with path: %s\n' %library_file_with_path)

    #Using pyBigWig for BigWig and BigBed files if you can import pyBigWig (linux only) otherwise no
    #By the way pyBigWig can be imported in unix, linux like os not available in windows
    #Using HM and CTCF bed files preparing chr based signal array runtime
    #Using ATAC-seq wig files preparing chr based signal array runtime

    ##########################################################################
    # If chunksize is 1, maxtasksperchild=x will call the function x times in each process,
    # but if chunksize is y, it will call the function x*y times in each process.
    # Setting maxtaskperchild to 1 would restart each process in your pool after it processed a single task, which is the most aggressive setting you could use to free any leaked resources.
    # numofProcesses = multiprocessing.cpu_count()
    # pool = multiprocessing.Pool(numofProcesses, maxtasksperchild=1)
    ##########################################################################

    ##########################################################################
    #July 26, 2020
    #For Vectorization
    #These are used in writing tables
    subsSignatures = np.append(ordered_sbs_signatures, AGGREGATEDSUBSTITUTIONS)
    dinucsSignatures = np.append(ordered_dbs_signatures, AGGREGATEDDINUCS)
    indelsSignatures = np.append(ordered_id_signatures, AGGREGATEDINDELS)

    allSims_subsSignature_accumulated_signal_np_array = np.zeros((numofSimulations+1,ordered_sbs_signatures.size+1, plusorMinus * 2 + 1))
    allSims_dinucsSignature_accumulated_signal_np_array = np.zeros((numofSimulations+1,ordered_dbs_signatures.size+1, plusorMinus * 2 + 1))
    allSims_indelsSignature_accumulated_signal_np_array = np.zeros((numofSimulations+1,ordered_id_signatures.size+1, plusorMinus * 2 + 1))

    allSims_subsSignature_accumulated_count_np_array = np.zeros((numofSimulations+1,ordered_sbs_signatures.size+1, plusorMinus * 2 + 1))
    allSims_dinucsSignature_accumulated_count_np_array = np.zeros((numofSimulations+1,ordered_dbs_signatures.size+1, plusorMinus * 2 + 1))
    allSims_indelsSignature_accumulated_count_np_array = np.zeros((numofSimulations+1,ordered_id_signatures.size+1, plusorMinus * 2 + 1))
    ##########################################################################

    ##############################################################
    # This code reads the file abd preaore chrbased signal files
    # If file is in default files,chr based signal files are downloaded from ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerTopography/lib/
    # No need for preparing here
    library_file_type=None

    if (os.path.basename(library_file_with_path) not in SIGPROFILERTOPOGRAPHY_DEFAULT_FILES):
        #What is the type of the signal_file_with_path?
        #If it is a bed file read signal_file_with_path here
        file_extension = os.path.splitext(os.path.basename(library_file_with_path))[1]

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
            isFileTypeBEDGRAPH=decideFileType(library_file_with_path)
            if isFileTypeBEDGRAPH:
                if verbose: start_time = time.time()
                #Read by chunks
                # readWig_write_derived_from_bedgraph_using_pool_chunks(outputDir, jobname, genome, library_file_with_path,occupancy_type,remove_outliers,verbose)
                #Read at once
                readWig_write_derived_from_bedgraph_using_pool_read_all(outputDir, jobname, genome, library_file_with_path, occupancy_type,remove_outliers,verbose,quantileValue)
                if verbose: print('\tVerbose Read wig file and write chrbased arrays took %f seconds' %((time.time() - start_time)))

                #For 6 GB ATAC-seq file using pool took 8 min whereas without pool took 16 min.
                # start_time = time.time()
                # readWig_write_derived_from_bedgraph(outputDir, jobname, genome, library_file_with_path,occupancy_type,verbose)
                # print('Without pool Took %f seconds' %((time.time() - start_time)))
            else:
                readWig_with_fixedStep_variableStep_writeChrBasedSignalArrays(outputDir, jobname, genome, library_file_with_path,occupancy_type,quantileValue,remove_outliers)
        elif (file_extension.lower()=='.bedgraph'):
            library_file_type = BEDGRAPH
            readWig_write_derived_from_bedgraph_using_pool_read_all(outputDir, jobname, genome, library_file_with_path,occupancy_type, remove_outliers, verbose, quantileValue)
        else:
            library_file_type=LIBRARY_FILE_TYPE_OTHER

    # elif os.path.basename(library_file_with_path)==DEFAULT_ATAC_SEQ_OCCUPANCY_FILE:
    #     check_download_chrbased_npy_atac_seq_files(outputDir, jobname, occupancy_type,DEFAULT_ATAC_SEQ_OCCUPANCY_FILE,chromNamesList)
    ##############################################################


    #########################################################################################
    def accumulate_apply_async_result_vectorization(simulatonBased_SignalArrayAndCountArrayList):
        try:
            chrLong=simulatonBased_SignalArrayAndCountArrayList[0]
            simNum=simulatonBased_SignalArrayAndCountArrayList[1]
            subsSignature_accumulated_signal_np_array=simulatonBased_SignalArrayAndCountArrayList[2]
            dinucsSignature_accumulated_signal_np_array=simulatonBased_SignalArrayAndCountArrayList[3]
            indelsSignature_accumulated_signal_np_array=simulatonBased_SignalArrayAndCountArrayList[4]
            subsSignature_accumulated_count_np_array=simulatonBased_SignalArrayAndCountArrayList[5]
            dinucsSignature_accumulated_count_np_array=simulatonBased_SignalArrayAndCountArrayList[6]
            indelsSignature_accumulated_count_np_array=simulatonBased_SignalArrayAndCountArrayList[7]

            #Accumulation
            allSims_subsSignature_accumulated_signal_np_array[simNum] += subsSignature_accumulated_signal_np_array
            allSims_dinucsSignature_accumulated_signal_np_array[simNum] += dinucsSignature_accumulated_signal_np_array
            allSims_indelsSignature_accumulated_signal_np_array[simNum] += indelsSignature_accumulated_signal_np_array

            allSims_subsSignature_accumulated_count_np_array[simNum] += subsSignature_accumulated_count_np_array
            allSims_dinucsSignature_accumulated_count_np_array[simNum] += dinucsSignature_accumulated_count_np_array
            allSims_indelsSignature_accumulated_count_np_array[simNum] += indelsSignature_accumulated_count_np_array
            # print('ACCUMULATION chrLong:%s simNum:%d ENDS' %(chrLong,simNum))
        except Exception as e:
            print("Exception: %s" % (e))
    #########################################################################################

    jobs=[]

    ##################################################################################
    if (computation_type == USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM):
        print(USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM,flush=True)

        sim_nums = range(0, numofSimulations + 1)
        sim_num_chr_tuples = ((sim_num, chrLong) for sim_num in sim_nums for chrLong in chromNamesList)

        numofProcesses = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=numofProcesses)

        for simNum, chrLong in sim_num_chr_tuples:
            jobs.append(pool.apply_async(chrbased_data_fill_signal_count_arrays_for_all_mutations_read_mutations,
                             args=(occupancy_type,
                                   occupancy_calculation_type,
                                   outputDir,
                                   jobname,
                                   chrLong,
                                   simNum,
                                   chromSizesDict,
                                   library_file_with_path,
                                   library_file_type,
                                   ordered_sbs_signatures,
                                   ordered_dbs_signatures,
                                   ordered_id_signatures,
                                   ordered_sbs_signatures_cutoffs,
                                   ordered_dbs_signatures_cutoffs,
                                   ordered_id_signatures_cutoffs,
                                   plusorMinus,
                                   verbose,),
                             callback=accumulate_apply_async_result_vectorization))

        pool.close()
        pool.join()
    ##################################################################################


    ##################################################################################
    elif (computation_type==USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT):
        print(USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT, flush=True)

        numofProcesses = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=numofProcesses)

        for chrLong, simNum, splitIndex in job_tuples:
            jobs.append(pool.apply_async(chrbased_data_fill_signal_count_arrays_for_all_mutations_read_mutations_split,
                                         args=(occupancy_type,
                                               occupancy_calculation_type,
                                               outputDir,
                                               jobname,
                                               chrLong,
                                               simNum,
                                               splitIndex,
                                               chromSizesDict,
                                               library_file_with_path,
                                               library_file_type,
                                               ordered_sbs_signatures,
                                               ordered_dbs_signatures,
                                               ordered_id_signatures,
                                               ordered_sbs_signatures_cutoffs,
                                               ordered_dbs_signatures_cutoffs,
                                               ordered_id_signatures_cutoffs,
                                               plusorMinus,
                                               verbose,),
                                         callback=accumulate_apply_async_result_vectorization))
        pool.close()
        pool.join()
    ##################################################################################


    ##################################################################################
    #July 26, 2020, For Vectorization
    writeSimulationBasedAverageNucleosomeOccupancyUsingNumpyArray(occupancy_type,
                                                   sample_based,
                                                   plusorMinus,
                                                   subsSignatures,
                                                   dinucsSignatures,
                                                   indelsSignatures,
                                                   allSims_subsSignature_accumulated_signal_np_array,
                                                   allSims_dinucsSignature_accumulated_signal_np_array,
                                                   allSims_indelsSignature_accumulated_signal_np_array,
                                                   allSims_subsSignature_accumulated_count_np_array,
                                                   allSims_dinucsSignature_accumulated_count_np_array,
                                                   allSims_indelsSignature_accumulated_count_np_array,
                                                   outputDir,
                                                   jobname,
                                                   numofSimulations,
                                                   library_file_memo)
    ##################################################################################

    print('--- %s Analysis ends' %(occupancy_type))
    print('#################################################################################\n')

########################################################################################
