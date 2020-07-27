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

from SigProfilerTopography.source.commons.TopographyCommons import accumulateChrBasedSimBasedSplitBasedArrays
from SigProfilerTopography.source.commons.TopographyCommons import accumulateChrBasedSimBasedSplitBasedSampleBasedArrays

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

from SigProfilerTopography.source.commons.TopographyCommons import USING_IMAP_UNORDERED

from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM
from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT
from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT_USING_POOL_INPUT_LIST


from SigProfilerTopography.source.commons.TopographyCommons import MAXIMUM_NUMBER_JOBS_IN_THE_POOL_AT_ONCE

from SigProfilerTopography.source.nucleosomeoccupancy.ChrBasedSignalArrays import readBEDandWriteChromBasedSignalArrays
from SigProfilerTopography.source.nucleosomeoccupancy.ChrBasedSignalArrays import readWig_with_fixedStep_variableStep_writeChrBasedSignalArrays
from SigProfilerTopography.source.nucleosomeoccupancy.ChrBasedSignalArrays import readWig_write_derived_from_bedgraph

from SigProfilerTopography.source.nucleosomeoccupancy.ChrBasedSignalArrays import readWig_write_derived_from_bedgraph_using_pool_chunks
from SigProfilerTopography.source.nucleosomeoccupancy.ChrBasedSignalArrays import readWig_write_derived_from_bedgraph_using_pool_read_all

from SigProfilerTopography.source.nucleosomeoccupancy.ChrBasedSignalArrays import readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArraysInParallel
from SigProfilerTopography.source.nucleosomeoccupancy.ChrBasedSignalArrays import readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArraysSequentially

from SigProfilerTopography.source.commons.TopographyCommons import decideFileType
from SigProfilerTopography.source.commons.TopographyCommons import get_chrBased_simBased_combined_df_split
from SigProfilerTopography.source.commons.TopographyCommons import get_chrBased_simBased_combined_chunks_df
from SigProfilerTopography.source.commons.TopographyCommons import get_chrBased_simBased_combined_df

##############################################################################################################
#main function
def occupancyAnalysis(genome,computationType,occupancy_type,sample_based,plusorMinus,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,job_tuples,library_file_with_path,library_file_memo,subsSignature_cutoff_numberofmutations_averageprobability_df,indelsSignature_cutoff_numberofmutations_averageprobability_df,dinucsSignature_cutoff_numberofmutations_averageprobability_df,remove_outliers,quantileValue,verbose):

    print('\n#################################################################################')
    print('--- %s Analysis starts' %(occupancy_type))
    print('--- Computation Type:%s' % (computationType))
    print('--- Occupancy Type:%s' % (occupancy_type))
    print('--- Library file with path: %s\n' %library_file_with_path)

    #Using pyBigWig for BigWig and BigBed files if you can import pyBigWig (linux only) otherwise no
    #By the way pyBigWig can be imported in unix, linux like os not available in windows
    #Using HM and CTCF bed files preparing chr based signal array runtime
    #Using ATAC-seq wig files preparing chr based signal array runtime
    occupancy_analysis(genome,computationType,occupancy_type,sample_based,plusorMinus,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,job_tuples,library_file_with_path,library_file_memo,subsSignature_cutoff_numberofmutations_averageprobability_df,indelsSignature_cutoff_numberofmutations_averageprobability_df,dinucsSignature_cutoff_numberofmutations_averageprobability_df,remove_outliers,quantileValue,verbose)
    print('--- %s Analysis ends' %(occupancy_type))
    print('#################################################################################\n')
##############################################################################################################


########################################################################################
# April 27, 2020
# requires chrBased_simBased_combined_df_split which can be real split or whole in fact
# This is common for pool.imap_unordered and pool.apply_async variations
def chrbased_data_fill_signal_count_arrays_for_all_mutations(occupancy_type,
                                                             outputDir,
                                                             jobname,
                                                             chrLong,
                                                             simNum,
                                                             chrBased_simBased_combined_df_split,
                                                             chromSizesDict,
                                                             library_file_with_path,
                                                             library_file_type,
                                                             sample2NumberofSubsDict,
                                                             sample2SubsSignature2NumberofMutationsDict,
                                                             subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                             indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                                             dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                                             plusorMinus,
                                                             sample_based,
                                                             verbose):

    if verbose: print('\tVerbose %s Worker pid %s memory_usage in %.2f MB START chrLong:%s simNum:%d' %(occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum))
    # 1st part  Prepare chr based mutations dataframes
    maximum_chrom_size = chromSizesDict[chrLong]
    start_time = time.time()

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
            df_columns = list(chrBased_simBased_combined_df_split.columns.values)
            #Initialize
            subsSignatures = subsSignature_cutoff_numberofmutations_averageprobability_df['signature'].values
            dinucsSignatures = dinucsSignature_cutoff_numberofmutations_averageprobability_df['signature'].values
            indelsSignatures = indelsSignature_cutoff_numberofmutations_averageprobability_df['signature'].values

            subsSignatures_cutoffs = subsSignature_cutoff_numberofmutations_averageprobability_df['cutoff'].values
            dinucsSignatures_cutoffs = dinucsSignature_cutoff_numberofmutations_averageprobability_df['cutoff'].values
            indelsSignatures_cutoffs = indelsSignature_cutoff_numberofmutations_averageprobability_df['cutoff'].values

            subsSignatures_mask_array = np.array([True if column_name in subsSignatures else False for column_name in df_columns])
            dinucsSignatures_mask_array = np.array([True if column_name in dinucsSignatures else False for column_name in df_columns])
            indelsSignatures_mask_array = np.array([True if column_name in indelsSignatures else False for column_name in df_columns])

            #Add one more row for the aggregated analysis
            subsSignature_accumulated_signal_np_array=np.zeros((subsSignatures.size+1,plusorMinus*2+1))
            dinucsSignature_accumulated_signal_np_array=np.zeros((dinucsSignatures.size+1,plusorMinus*2+1))
            indelsSignature_accumulated_signal_np_array=np.zeros((indelsSignatures.size+1,plusorMinus*2+1))

            #Add one more row for the aggregated analysis
            subsSignature_accumulated_count_np_array=np.zeros((subsSignatures.size+1,plusorMinus*2+1))
            dinucsSignature_accumulated_count_np_array=np.zeros((dinucsSignatures.size+1,plusorMinus*2+1))
            indelsSignature_accumulated_count_np_array=np.zeros((indelsSignatures.size+1,plusorMinus*2+1))

            #July 25, 2020
            [fillSignalArrayAndCountArray_using_list_comp_vectorization(
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
                subsSignatures_cutoffs,
                dinucsSignatures_cutoffs,
                indelsSignatures_cutoffs,
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
                sample_based,
                df_columns) for row in chrBased_simBased_combined_df_split[df_columns].values]


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
# May 5, 2020
# For apply_async
# Read chromBased simBased combined mutations df in the process
def chrbased_data_fill_signal_count_arrays_for_all_mutations_read_mutations(occupancy_type,
                                                                            outputDir,
                                                                            jobname,
                                                                            chrLong,
                                                                            simNum,
                                                                            chromSizesDict,
                                                                            library_file_with_path,
                                                                            library_file_type,
                                                                            sample2NumberofSubsDict,
                                                                            sample2SubsSignature2NumberofMutationsDict,
                                                                            subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                            indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                            dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                            plusorMinus,
                                                                            sample_based,
                                                                            verbose):

    chrBased_simBased_mutations_df = get_chrBased_simBased_combined_df(outputDir, jobname, chrLong, simNum)

    return chrbased_data_fill_signal_count_arrays_for_all_mutations(occupancy_type,
                                                                    outputDir,
                                                                    jobname,
                                                                    chrLong,
                                                                    simNum,
                                                                    chrBased_simBased_mutations_df,
                                                                    chromSizesDict,
                                                                    library_file_with_path,
                                                                    library_file_type,
                                                                    sample2NumberofSubsDict,
                                                                    sample2SubsSignature2NumberofMutationsDict,
                                                                    subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                    indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                    dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                    plusorMinus,
                                                                    sample_based,
                                                                    verbose)
########################################################################################

########################################################################################
#May 19, 2020
# For apply_async split using poolInputList
# Read chromBased simBased combined mutations df split in the process
def chrbased_data_fill_signal_count_arrays_for_all_mutations_read_mutations_split(occupancy_type, outputDir, jobname, chrLong, simNum, splitIndex,
                                                                                  chromSizesDict, library_file_with_path,
                                                                                  library_file_type, sample2NumberofSubsDict, sample2SubsSignature2NumberofMutationsDict,
                                                                                  subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                                  indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                                  dinucsSignature_cutoff_numberofmutations_averageprobability_df, plusorMinus, sample_based, verbose):

    chrBased_simBased_combined_df_split = get_chrBased_simBased_combined_df_split(outputDir, jobname, chrLong, simNum,splitIndex)

    return chrbased_data_fill_signal_count_arrays_for_all_mutations(occupancy_type, outputDir, jobname, chrLong, simNum,
                                                                    chrBased_simBased_combined_df_split, chromSizesDict,
                                                                    library_file_with_path, library_file_type,
                                                                    sample2NumberofSubsDict,
                                                                    sample2SubsSignature2NumberofMutationsDict,
                                                                    subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                    indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                    dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                    plusorMinus, sample_based, verbose)
########################################################################################


########################################################################################
#July 25, 2020, Vectorization
def fillSignalArrayAndCountArray_using_list_comp_vectorization(
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
        subsSignatures_cutoffs,
        dinucsSignatures_cutoffs,
        indelsSignatures_cutoffs,
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
        sample_based,
        df_columns):

    indexofType = df_columns.index(TYPE)
    mutation_row_type = row[indexofType]

    ###########################################
    if mutation_row_type == SUBS:
        accumulated_signal_np_array=subsSignature_accumulated_signal_np_array
        accumulated_count_np_array=subsSignature_accumulated_count_np_array
        cutoffs=subsSignatures_cutoffs
        signatures_mask_array=subsSignatures_mask_array
    elif mutation_row_type == DINUCS:
        accumulated_signal_np_array = dinucsSignature_accumulated_signal_np_array
        accumulated_count_np_array = dinucsSignature_accumulated_count_np_array
        cutoffs=dinucsSignatures_cutoffs
        signatures_mask_array=dinucsSignatures_mask_array
    elif mutation_row_type == INDELS:
        accumulated_signal_np_array=indelsSignature_accumulated_signal_np_array
        accumulated_count_np_array=indelsSignature_accumulated_count_np_array
        cutoffs=indelsSignatures_cutoffs
        signatures_mask_array=indelsSignatures_mask_array
    ###########################################

    window_array=None
    windowSize=plusOrMinus*2+1

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

    #Get the sample at this mutation_row
    # sample = mutation_row_sample
    # simulationNumber= mutation_row_simulation_number

    #Vectorize July 25, 2020
    #Fill dictionaries using window_array
    if (window_array is not None) and (np.any(window_array)):
        probabilities = row[signatures_mask_array]
        threshold_mask_array = np.greater(probabilities, cutoffs)

        #Convert True into 1, and False into 0
        mask_array = threshold_mask_array.astype(int)

        #Add 1 for the aggregated analysis to the mask array
        mask_array = np.append(mask_array, 1)

        #Add one more dimension to window_array and mask_array
        window_array_1x2001=np.array([window_array])
        mask_array_1xnumofsignatures=np.array([mask_array])

        to_be_accumulated_array = mask_array_1xnumofsignatures.T * window_array_1x2001
        accumulated_signal_np_array += to_be_accumulated_array
        accumulated_count_np_array += (to_be_accumulated_array>0)
########################################################################################


########################################################################################
#TODO To be depreceated
#April 6, 2020
def fillSignalArrayAndCountArray_using_list_comp(
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
        subsSignature_cutoff_numberofmutations_averageprobability_df,
        indelsSignature_cutoff_numberofmutations_averageprobability_df,
        dinucsSignature_cutoff_numberofmutations_averageprobability_df,
        plusOrMinus,
        sample_based,
        df_columns):

    indexofType = df_columns.index(TYPE)
    mutation_row_type = row[indexofType]

    ###########################################
    if mutation_row_type==SUBS:
        my_type=AGGREGATEDSUBSTITUTIONS
        signature_cutoff_numberofmutations_averageprobability_df=subsSignature_cutoff_numberofmutations_averageprobability_df
    elif mutation_row_type == INDELS:
        my_type=AGGREGATEDINDELS
        signature_cutoff_numberofmutations_averageprobability_df=indelsSignature_cutoff_numberofmutations_averageprobability_df
    elif mutation_row_type == DINUCS:
        my_type=AGGREGATEDDINUCS
        signature_cutoff_numberofmutations_averageprobability_df=dinucsSignature_cutoff_numberofmutations_averageprobability_df
    ###########################################

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

    #Fill dictionaries using window_array
    if (window_array is not None) and (np.any(window_array)):
        #TODO: Is there a faster way than using for loop? Yes, there is a way.
        ################# Signatures starts #######################
        #mutation_row[signature] mutation probability for that signature
        #signature2PropertiesListDict[signature][0] cutoff probability for that signature

        for signature in signature_cutoff_numberofmutations_averageprobability_df['signature'].unique():
            indexofSignature = df_columns.index(signature)
            mutation_row_signature = row[indexofSignature]

            cutoff=float(signature_cutoff_numberofmutations_averageprobability_df[signature_cutoff_numberofmutations_averageprobability_df['signature']==signature]['cutoff'].values[0])

            if (mutation_row_signature >= cutoff):
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



################################
#May 19, 2020
def my_initializer(my_occupancy_type,
                   my_outputDir,
                   my_jobname,
                   my_chromSizesDict,
                   my_library_file_with_path,
                   my_library_file_type,
                   my_sample2NumberofSubsDict,
                   my_sample2SubsSignature2NumberofMutationsDict,
                   my_subsSignature_cutoff_numberofmutations_averageprobability_df,
                   my_indelsSignature_cutoff_numberofmutations_averageprobability_df,
                   my_dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                   my_plusorMinus,
                   my_sample_based,
                   my_verbose):

    global occupancy_type
    global outputDir
    global jobname
    global chromSizesDict
    global library_file_with_path
    global library_file_type
    global sample2NumberofSubsDict
    global sample2SubsSignature2NumberofMutationsDict
    global subsSignature_cutoff_numberofmutations_averageprobability_df
    global indelsSignature_cutoff_numberofmutations_averageprobability_df
    global dinucsSignature_cutoff_numberofmutations_averageprobability_df
    global plusorMinus
    global sample_based
    global verbose

    occupancy_type = my_occupancy_type
    outputDir=my_outputDir
    jobname=my_jobname
    chromSizesDict=my_chromSizesDict
    library_file_with_path =my_library_file_with_path
    library_file_type = my_library_file_type
    sample2NumberofSubsDict= my_sample2NumberofSubsDict
    sample2SubsSignature2NumberofMutationsDict=my_sample2SubsSignature2NumberofMutationsDict
    subsSignature_cutoff_numberofmutations_averageprobability_df=my_subsSignature_cutoff_numberofmutations_averageprobability_df
    indelsSignature_cutoff_numberofmutations_averageprobability_df=my_indelsSignature_cutoff_numberofmutations_averageprobability_df
    dinucsSignature_cutoff_numberofmutations_averageprobability_df=my_dinucsSignature_cutoff_numberofmutations_averageprobability_df
    plusorMinus=my_plusorMinus
    sample_based=my_sample_based
    verbose=my_verbose
################################

################################
#May 19, 2020
def chrbased_data_fill_signal_count_arrays_for_all_mutations_imap_unordered(tuple):
    chrLong, simNum, chunk =tuple
    # print('MONITOR %s simNum:%d chunk' %(chrLong,simNum),flush=True)

    return chrbased_data_fill_signal_count_arrays_for_all_mutations(occupancy_type, outputDir, jobname, chrLong, simNum,
                                                                    chunk, chromSizesDict,
                                                                    library_file_with_path, library_file_type,
                                                                    sample2NumberofSubsDict,
                                                                    sample2SubsSignature2NumberofMutationsDict,
                                                                    subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                    indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                    dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                    plusorMinus, sample_based, verbose)

################################


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
                        job_tuples,
                        library_file_with_path,
                        library_file_memo,
                        subsSignature_cutoff_numberofmutations_averageprobability_df,
                        indelsSignature_cutoff_numberofmutations_averageprobability_df,
                        dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                        remove_outliers,
                        quantileValue,
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
    # numofProcesses = multiprocessing.cpu_count()
    # pool = multiprocessing.Pool(numofProcesses, maxtasksperchild=1)
    ##########################################################################


    ##########################################################################
    #July 26, 2020
    #For Vectorization
    subsSignatures = subsSignature_cutoff_numberofmutations_averageprobability_df['signature'].values
    dinucsSignatures = dinucsSignature_cutoff_numberofmutations_averageprobability_df['signature'].values
    indelsSignatures = indelsSignature_cutoff_numberofmutations_averageprobability_df['signature'].values

    subsSignatures = np.append(subsSignatures, AGGREGATEDSUBSTITUTIONS)
    dinucsSignatures = np.append(dinucsSignatures, AGGREGATEDDINUCS)
    indelsSignatures = np.append(indelsSignatures, AGGREGATEDINDELS)

    allSims_subsSignature_accumulated_signal_np_array = np.zeros((numofSimulations+1,subsSignatures.size, plusorMinus * 2 + 1))
    allSims_dinucsSignature_accumulated_signal_np_array = np.zeros((numofSimulations+1,dinucsSignatures.size, plusorMinus * 2 + 1))
    allSims_indelsSignature_accumulated_signal_np_array = np.zeros((numofSimulations+1,indelsSignatures.size, plusorMinus * 2 + 1))

    allSims_subsSignature_accumulated_count_np_array = np.zeros((numofSimulations+1,subsSignatures.size, plusorMinus * 2 + 1))
    allSims_dinucsSignature_accumulated_count_np_array = np.zeros((numofSimulations+1,dinucsSignatures.size, plusorMinus * 2 + 1))
    allSims_indelsSignature_accumulated_count_np_array = np.zeros((numofSimulations+1,indelsSignatures.size, plusorMinus * 2 + 1))
    ##########################################################################

    ##############################################################
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
    ##############################################################


    #########################################################################################
    def accumulate_apply_async_result_vectorization(simulatonBased_SignalArrayAndCountArrayList):
        simNum=simulatonBased_SignalArrayAndCountArrayList[0]
        subsSignature_accumulated_signal_np_array=simulatonBased_SignalArrayAndCountArrayList[1]
        dinucsSignature_accumulated_signal_np_array=simulatonBased_SignalArrayAndCountArrayList[2]
        indelsSignature_accumulated_signal_np_array=simulatonBased_SignalArrayAndCountArrayList[3]
        subsSignature_accumulated_count_np_array=simulatonBased_SignalArrayAndCountArrayList[4]
        dinucsSignature_accumulated_count_np_array=simulatonBased_SignalArrayAndCountArrayList[5]
        indelsSignature_accumulated_count_np_array=simulatonBased_SignalArrayAndCountArrayList[6]

        print('simNum:%s ACCUMULATION' %(simNum))
        #Accumulation
        allSims_subsSignature_accumulated_signal_np_array[simNum] += subsSignature_accumulated_signal_np_array
        allSims_dinucsSignature_accumulated_signal_np_array[simNum] += dinucsSignature_accumulated_signal_np_array
        allSims_indelsSignature_accumulated_signal_np_array[simNum] += indelsSignature_accumulated_signal_np_array

        allSims_subsSignature_accumulated_count_np_array[simNum] += subsSignature_accumulated_count_np_array
        allSims_dinucsSignature_accumulated_count_np_array[simNum] += dinucsSignature_accumulated_count_np_array
        allSims_indelsSignature_accumulated_count_np_array[simNum] += indelsSignature_accumulated_count_np_array
    #########################################################################################


    ##################################################################################
    if (computation_type == USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM):
        print(USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM,flush=True)

        sim_nums = range(0, numofSimulations + 1)
        sim_num_chr_tuples = ((sim_num, chrLong) for sim_num in sim_nums for chrLong in chromNamesList)

        ################################
        numofProcesses = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=numofProcesses)
        ################################

        ################################
        jobs=[]
        ################################

        ################################
        for simNum, chrLong in sim_num_chr_tuples:
            jobs.append(pool.apply_async(chrbased_data_fill_signal_count_arrays_for_all_mutations_read_mutations,
                                         args=(occupancy_type,
                                               outputDir,
                                               jobname,
                                               chrLong,
                                               simNum,
                                               chromSizesDict,
                                               library_file_with_path,
                                               library_file_type,
                                               sample2NumberofSubsDict,
                                               sample2SubsSignature2NumberofMutationsDict,
                                               subsSignature_cutoff_numberofmutations_averageprobability_df,
                                               indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                               dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                               plusorMinus,
                                               sample_based,
                                               verbose,),
                                         callback=accumulate_apply_async_result_vectorization))

            print('MONITOR %s %d len(jobs):%d' %(chrLong,simNum,len(jobs)),flush=True)
        ################################

        ##############################################################################
        # wait for all jobs to finish
        for job in jobs:
            if verbose: print('\tVerbose %s Worker pid %s job.get():%s ' % (occupancy_type, str(os.getpid()), job.get()))
        ##############################################################################

        ################################
        pool.close()
        pool.join()
        ################################
    ##################################################################################


    ##################################################################################
    elif (computation_type==USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT):
        print(USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT, flush=True)

        ################################
        numofProcesses = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=numofProcesses)
        ################################

        ################################
        jobs=[]
        ################################

        ################################
        for chrLong, simNum, splitIndex in job_tuples:
            jobs.append(pool.apply_async(chrbased_data_fill_signal_count_arrays_for_all_mutations_read_mutations_split,
                                         args=(occupancy_type,
                                               outputDir,
                                               jobname,
                                               chrLong,
                                               simNum,
                                               splitIndex,
                                               chromSizesDict,
                                               library_file_with_path,
                                               library_file_type,
                                               sample2NumberofSubsDict,
                                               sample2SubsSignature2NumberofMutationsDict,
                                               subsSignature_cutoff_numberofmutations_averageprobability_df,
                                               indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                               dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                               plusorMinus,
                                               sample_based,
                                               verbose,),
                                         callback=accumulate_apply_async_result_vectorization))

            print('MONITOR %s %d len(jobs):%d' % (chrLong, simNum, len(jobs)), flush=True)
        ################################

        ##############################################################################
        # wait for all jobs to finish
        for job in jobs:
            if verbose: print('\tVerbose %s Worker pid %s job.get():%s ' % (occupancy_type, str(os.getpid()), job.get()))
        ##############################################################################

        ################################
        pool.close()
        pool.join()
        ################################
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


########################################################################################