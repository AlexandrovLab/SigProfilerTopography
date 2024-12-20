# !/usr/bin/env python3

# Author: burcakotlu

# Contact: burcakotlu@eng.ucsd.edu

# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018-2020 Burcak Otlu

#############################################################
# This python code generates normalized mutation density data for different analyses types which are
# Using all mutations (all chromosomes, all samples, signatures) which leads to Aggregated Substitutions
# Using all mutations but in signature based manner leads to Signature Based
# Using all mutations in signature based and sample based manner leads to Signature based and Sample Based
# Using all indels (all chromosomes, all samples) which leads to Aggregated Indels
# Using all indels (diving indels of length >=3 as microhomology-mediated indels and indels of length <3 as repeat-mediated  indels)

## In this python code
## I read the wavelet-transformed signal wig across whole genome
## Sort the signal in descending order
## Divide the data into 10 equal deciles
## Return 10 deciles as data frames
## Read the mutation file
## Group by mutation file by chromosome
## Group by decile file by chromosome
## For each chromosome (Done in the loop sequentially)
## Make a list of input
## chromosome based mutation file mutation file
## chromosome based decile files coming from 10 different deciles (Done in parallel for all of the 10 deciles)
## In the function I will generate interval tree from decile file
## Overlap the mutations with the interval tree
## Return the total number of mutations that overlap between mutation file and decile file
## Get chromosome based mutation file
## Get chromosome based decile file
## Accumulate the total number of mutations that overlap for each decile file
#############################################################

#############################################################
# Constraints, Thresholds
# Please note that for sample based and signature based replication time analysis
# We consider samples with at least 3000 mutations at total in all deciles.
#############################################################

import multiprocessing
import sys
import os
import math
import numpy as np
import pandas as pd
import pyranges as pr

from SigProfilerTopography.source.commons.TopographyCommons import DATA
from SigProfilerTopography.source.commons.TopographyCommons import REPLICATIONTIME
from SigProfilerTopography.source.commons.TopographyCommons import SAMPLES
from SigProfilerTopography.source.commons.TopographyCommons import SIGNATUREBASED

from SigProfilerTopography.source.commons.TopographyCommons import NUMOFBASES

from SigProfilerTopography.source.commons.TopographyCommons import CHROM
from SigProfilerTopography.source.commons.TopographyCommons import START
from SigProfilerTopography.source.commons.TopographyCommons import END
from SigProfilerTopography.source.commons.TopographyCommons import SIGNAL

from SigProfilerTopography.source.commons.TopographyCommons import SAMPLE
from SigProfilerTopography.source.commons.TopographyCommons import LENGTH
from SigProfilerTopography.source.commons.TopographyCommons import SIMULATION_NUMBER

from SigProfilerTopography.source.commons.TopographyCommons import TYPE
from SigProfilerTopography.source.commons.TopographyCommons import SUBS
from SigProfilerTopography.source.commons.TopographyCommons import INDELS
from SigProfilerTopography.source.commons.TopographyCommons import DINUCS

from SigProfilerTopography.source.commons.TopographyCommons import AGGREGATEDSUBSTITUTIONS
from SigProfilerTopography.source.commons.TopographyCommons import AGGREGATEDINDELS
from SigProfilerTopography.source.commons.TopographyCommons import AGGREGATEDDINUCS

from SigProfilerTopography.source.commons.TopographyCommons import MICROHOMOLOGY
from SigProfilerTopography.source.commons.TopographyCommons import REPEAT

from SigProfilerTopography.source.commons.TopographyCommons import MUTATION
from SigProfilerTopography.source.commons.TopographyCommons import REPLICATIONTIMECOLUMN

from SigProfilerTopography.source.commons.TopographyCommons import SBS
from SigProfilerTopography.source.commons.TopographyCommons import DBS
from SigProfilerTopography.source.commons.TopographyCommons import ID

from SigProfilerTopography.source.commons.TopographyCommons import write_chr_based_mutations_df
from SigProfilerTopography.source.commons.TopographyCommons import readWig_with_fixedStep_variableStep
from SigProfilerTopography.source.commons.TopographyCommons import memory_usage

from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM
from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT

from SigProfilerTopography.source.commons.TopographyCommons import Sample2NumberofDinucsDictFilename
from SigProfilerTopography.source.commons.TopographyCommons import Sample2DinucsSignature2NumberofMutationsDictFilename

from SigProfilerTopography.source.commons.TopographyCommons import getDictionary
from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofSubsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofIndelsDict

from SigProfilerTopography.source.commons.TopographyCommons import getSample2SubsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2IndelsSignature2NumberofMutationsDict

from SigProfilerTopography.source.commons.TopographyCommons import getChrShort

from SigProfilerTopography.source.commons.TopographyCommons import get_chrBased_simBased_combined_df_split

from SigProfilerTopography.source.commons.TopographyCommons import  MEGABYTE_IN_BYTES
from SigProfilerTopography.source.commons.TopographyCommons import  isFileTypeBedGraph

from SigProfilerTopography.source.commons.TopographyCommons import  get_chrBased_simBased_combined_df
from SigProfilerTopography.source.commons.TopographyCommons import  get_chrBased_simBased_dfs

# June 6, 2021
mutations_in_early_replicating_regions_df = pd.DataFrame()
mutations_in_late_replicating_regions_df = pd.DataFrame()


##################################################################
# Please note that this dictionary is copied from .../SigProfilerMatrixGenerator/SigProfilerMatrixGenerator/scripts/SigProfilerMatrixGeneratorFunc.py line 102
# it is defined in a fuction called SigProfilerMatrixGeneratorFunc
# If this dictionary is updated in SigProfilerMatrixGeneratorFunc.py, it has to be updated in ReplicationTimeAnalysis.py
# Provides the reference file conversion from binary to base information
tsb_ref = {
    0: ["N", "A"],
    1: ["N", "C"],
    2: ["N", "G"],
    3: ["N", "T"],
    4: ["T", "A"],
    5: ["T", "C"],
    6: ["T", "G"],
    7: ["T", "T"],
    8: ["U", "A"],
    9: ["U", "C"],
    10: ["U", "G"],
    11: ["U", "T"],
    12: ["B", "A"],
    13: ["B", "C"],
    14: ["B", "G"],
    15: ["B", "T"],
    16: ["N", "N"],
    17: ["T", "N"],
    18: ["U", "N"],
    19: ["B", "N"],
}
##################################################################



# Higher the replication time signal earlier the replication is
# Regions with high values indicate domains of early replication where initiaion occurs earlier in S-phase or early in a higher proportion of cells.
def read_repli_seq_time_data(genome, chromNamesList, repliseqDataFilename, matrix_generator_path, verbose, log_file):

    ###################################################################
    ############### Read RepliSeq Time data starts ####################
    ###################################################################

    filetype_BEDGRAPH = isFileTypeBedGraph(repliseqDataFilename)
    file_extension = os.path.splitext(os.path.basename(repliseqDataFilename))[1]

    if file_extension.lower() == '.bedgraph':
        column_names = [CHROM, START, END, SIGNAL]
        replication_time_interval_version_df = pd.read_csv(repliseqDataFilename, sep='\t', header=None, comment='#',
                                                           names=column_names,
                                                           dtype={CHROM: 'string', START: np.int32, END: np.int32, SIGNAL: np.float32}) # legacy category
    elif filetype_BEDGRAPH:
        column_names = [CHROM, START, END, SIGNAL]
        replication_time_interval_version_df = pd.read_csv(repliseqDataFilename, sep='\t', header=None, comment='#',
                                                           names=column_names,
                                                           dtype={CHROM: 'string', START: np.int32, END: np.int32, SIGNAL: np.float32}) # legacy category
    else:
        replication_time_interval_version_df = readWig_with_fixedStep_variableStep(repliseqDataFilename)

    chrNamesInReplicationTimeDataArray = replication_time_interval_version_df[CHROM].unique()
    log_out = open(log_file, 'a')
    print('Before --- Chromosome names in replication time signal data: %s replication_time_interval_version_df.shape(%d,%d)\n' %(chrNamesInReplicationTimeDataArray,replication_time_interval_version_df.shape[0],replication_time_interval_version_df.shape[1]), file=log_out)

    # Remove rows with chromosomes that are not in chromNamesList
    replication_time_interval_version_df = replication_time_interval_version_df[replication_time_interval_version_df[CHROM].isin(chromNamesList)]

    chrNamesInReplicationTimeDataArray = replication_time_interval_version_df[CHROM].unique()
    print('After considering only chromosomes in chromNamesList --- Chromosome names in replication time signal data: %s replication_time_interval_version_df.shape(%d,%d)\n' %(chrNamesInReplicationTimeDataArray,replication_time_interval_version_df.shape[0],replication_time_interval_version_df.shape[1]), file=log_out)
    log_out.close()

    # Augment wavelet_processed_df with numberofAttributableBases
    # each row has numberofAttributableBases
    wavelet_processed_augmented_df = augment(genome,
                                             replication_time_interval_version_df,
                                             matrix_generator_path,
                                             log_file,
                                             verbose)

    # Return 10 deciles: the first decile is the earliest one and the tenth decile is the latest one

    # Sort in descending order
    # Higher the replication time signal earlier the replication is
    wavelet_processed_augmented_df.sort_values(SIGNAL, ascending=False, inplace=True)

    # Split wavelet_processed_augmented_df into 10 deciles
    deciles_df_list = np.array_split(wavelet_processed_augmented_df,10)
    # print('Number of decile:%d' %len(deciles))
    # deciles_df_list is a list and each decile_df is a dataframe <class 'pandas.core.frame.DataFrame'>
    # The first decile_df is the earliest one
    # The last decile_df is the latest one
    # type(deciles_df_list) --> <class 'list'>

    return chrNamesInReplicationTimeDataArray, deciles_df_list
    ###################################################################
    ############### Read RepliSeq Time data ends ######################
    ###################################################################


def addNumofAttributableBasesColumnForApplyAsync(chrLong, chrBased_wavelet_processed_df_group, chrbased_file_path, log_file, verbose):
    if verbose:
        log_out = open(log_file, 'a')
        print('\tVerbose Worker pid %s %s Before Augment with number of attributable bases memory_usage %.2f MB' % (str(os.getpid()), chrLong, memory_usage()), file=log_out)
        log_out.close()

    with open(chrbased_file_path, "rb") as f2:
        chrom_string = f2.read()

    # old way
    # resulting_df = chrBased_wavelet_processed_df_group.apply(getNumberofAttributableBasesUsingMatrixGeneratorGenome,
    #                                                          chrom_string=chrom_string, axis= 1)
    # if (len(chrBased_wavelet_processed_df_group) != len(resulting_df)):
    #     log_out = open(log_file, 'a')
    #     print('There is a situation: len(chrBased_wavelet_processed_df_group) is not equal  to len(resulting_df)', file=log_out)
    #     log_out.close()
    # chrBased_wavelet_processed_df_group[NUMOFBASES] = resulting_df

    # new way
    # chrom_based_byte_array is of class bytes
    # contains only integers from 0 to 19
    # have look at tsb_ref
    # 16: ["N", "N"],
    # 17: ["T", "N"],
    # 18: ["U", "N"],
    # 19: ["B", "N"],
    # convert bytes to integers
    chrom_based_byte_array = np.frombuffer(chrom_string, dtype=np.int8)  # Convert bytes to integers

    num_of_attributable_bases = np.hstack([np.sum(chrom_based_byte_array[start:end+1] < 16) for start, end in zip(chrBased_wavelet_processed_df_group.iloc[:,1].values,
                                                                             chrBased_wavelet_processed_df_group.iloc[:,2].values)])

    chrBased_wavelet_processed_df_group[NUMOFBASES] = num_of_attributable_bases

    if (len(chrBased_wavelet_processed_df_group) != len(num_of_attributable_bases)):
        log_out = open(log_file, 'a')
        print('There is a situation: len(chrBased_wavelet_processed_df_group) is not equal  to len(all_indexes)', file=log_out)
        log_out.close()

    if verbose:
        log_out = open(log_file, 'a')
        print('\tVerbose Worker pid %s %s After Augment with number of attributable bases memory_usage %.2f MB' % (str(os.getpid()),chrLong, memory_usage()), file=log_out)
        log_out.close()

    return (chrLong, chrBased_wavelet_processed_df_group)

# sample based
# numpy array
# Main engine function
def search_for_each_mutation_using_list_comprehension_sample_based_using_numpy_array(mutation_row,
                                                                        mutation_type,
                                                                        chrBasedReplicationTimeDataArrayWithDecileIndex,
                                                                        ordered_signatures_cutoffs,
                                                                        signatures_mask_array,
                                                                        sample_signature_decile_index_accumulated_np_array,
                                                                        sample_2_sample_index_dict,
                                                                        discreet_mode,
                                                                        default_cutoff,
                                                                        df_columns):
    decile_index = None

    # df_columns: numpy array
    indexofStart = np.where(df_columns == START) [0][0]
    start = mutation_row[indexofStart]

    indexofSample = np.where(df_columns == SAMPLE) [0][0]
    sample = mutation_row[indexofSample]

    # For single point mutations start and end are the same, therefore we need to add 1 to mutation_row[END]
    # end is exclusive for subs, indels and dinucs provided by readChrBased methods

    # Simulation Number
    index_of_sim_num = np.where(df_columns == 'Simulation_Number')[0][0]
    sim_num = mutation_row[index_of_sim_num]

    if mutation_type == SUBS:
        end = start + 1
        # column_names = ['Sample', 'Chrom', 'Start', 'MutationLong']
        # num_of_rows = 4
    elif mutation_type == DINUCS:
        end = start + 2
        # column_names = ['Sample', 'Chrom', 'Start', 'MutationLong']
        # num_of_rows = 4
    elif mutation_type == INDELS:
        indexofLength = np.where(df_columns == LENGTH)
        length=mutation_row[indexofLength]
        end = start + int(length)
        # column_names = ['Sample', 'Chrom', 'Start', 'MutationLong', 'Ref', 'Alt', 'Length']
        # num_of_rows = 7

    slicedArray = chrBasedReplicationTimeDataArrayWithDecileIndex[start:end]

    # np.nonzero returns the indices of the elements that are non-zero.
    # np.unique finds the unique elements of an array returns ndarray the sorted unique values.
    # np.nditer efficient multi-dimensional iterator object to iterate over arrays.
    # We aim to get rid of zero if any exists in slicedArray.
    uniqueIndexesArray = np.unique(slicedArray[np.nonzero(slicedArray)])

    decile_index_array=np.zeros((10), dtype=int)

    # decile 10  will be accumulated in index 9
    # decile 1 will be accumulated in index 0
    # Therefore uniqueIndexesArray minus 1
    if (uniqueIndexesArray.size > 0):

        uniqueIndexesArray -= 1
        decile_index_array[uniqueIndexesArray] = 1

        # for mutation annotation
        if (decile_index_array[0]):
            # my_dict = {a: b for a, b in zip(df_columns, mutation_row)}
            # list_of_rep_timing_decile1_dict.append(my_dict)
            # list_of_repli_seq_timing_for_mutations.append(np.append(mutation_row[:num_of_rows], 1))
            decile_index = 1
        if (decile_index_array[1]):
            # my_dict = {a: b for a, b in zip(df_columns, mutation_row)}
            # list_of_rep_timing_decile2_dict.append(my_dict)
            # list_of_repli_seq_timing_for_mutations.append(np.append(mutation_row[:num_of_rows], 2))
            decile_index = 2
        if (decile_index_array[2]):
            # my_dict = {a: b for a, b in zip(df_columns, mutation_row)}
            # list_of_rep_timing_decile3_dict.append(my_dict)
            # list_of_repli_seq_timing_for_mutations.append(np.append(mutation_row[:num_of_rows], 3))
            decile_index = 3
        if (decile_index_array[3]):
            # my_dict = {a: b for a, b in zip(df_columns, mutation_row)}
            # list_of_rep_timing_decile4_dict.append(my_dict)
            # list_of_repli_seq_timing_for_mutations.append(np.append(mutation_row[:num_of_rows], 4))
            decile_index = 4
        if (decile_index_array[4]):
            # my_dict = {a: b for a, b in zip(df_columns, mutation_row)}
            # list_of_rep_timing_decile5_dict.append(my_dict)
            # list_of_repli_seq_timing_for_mutations.append(np.append(mutation_row[:num_of_rows], 5))
            decile_index = 5
        if (decile_index_array[5]):
            decile_index = 6
        if (decile_index_array[6]):
            decile_index = 7
        if (decile_index_array[7]):
            decile_index = 8
        if (decile_index_array[8]):
            decile_index = 9
        if (decile_index_array[9]):
            decile_index = 10

        probabilities = mutation_row[signatures_mask_array]

        # new code starts
        if discreet_mode:
            # Discreet way 1 or 0
            # Convert True into 1, and False into 0
            threshold_mask_array = np.greater_equal(probabilities, ordered_signatures_cutoffs)
            mask_array = threshold_mask_array.astype(int)
            if mutation_type == INDELS:
                if (length >= 3):
                    # Order is important
                    # MICROHOMOLOGY --> 1
                    # REPEAT --> 0
                    # AGGREGATEDINDELS --> 1
                    mask_array = np.append(mask_array, [1, 0, 1])
                else:
                    # Order is important
                    # MICROHOMOLOGY --> 0
                    # REPEAT --> 1
                    # AGGREGATEDINDELS --> 1
                    mask_array = np.append(mask_array, [0, 1, 1])

            else:
                # Add 1 for the aggregated analysis to the mask array
                # For SUBS and DINUCS
                mask_array = np.append(mask_array, 1)

        else:
            # new way
            probabilities[probabilities < default_cutoff ] = 0
            mask_array = np.array(probabilities).astype(float)

            # old way
            # mask_array = np.array(probabilities).astype(float)

            if mutation_type == INDELS:
                if (length >= 3):
                    # Order is important
                    # MICROHOMOLOGY --> 1
                    # REPEAT --> 0
                    # AGGREGATEDINDELS --> 1
                    mask_array = np.append(mask_array, [1.0, 0.0, 1.0])
                else:
                    # Order is important
                    # MICROHOMOLOGY --> 0
                    # REPEAT --> 1
                    # AGGREGATEDINDELS --> 1
                    mask_array = np.append(mask_array, [0.0, 1.0, 1.0])

            else:
                # Add 1 for the aggregated analysis to the mask array
                # For SUBS and DINUCS
                mask_array = np.append(mask_array, 1.0)
        # new code ends

        # Add 1 more dimension to the arrays
        decile_index_array_1x10 = np.expand_dims(decile_index_array, axis=0)
        mask_array_1xnumofsignatures = np.expand_dims(mask_array, axis=0)

        signatures_decile_index_np_array = mask_array_1xnumofsignatures.T * decile_index_array_1x10

        sample_index = sample_2_sample_index_dict[sample]

        # old code
        # signature_decile_index_accumulated_np_array += signatures_decile_index_np_array
        sample_signature_decile_index_accumulated_np_array[sample_index] += signatures_decile_index_np_array

    else:
        decile_index = -1 # Unknown

    return decile_index


# Using numpy array
# Main engine function
def search_for_each_mutation_using_list_comprehension_using_numpy_array(mutation_row,
                                                                        mutation_type,
                                                                        chrBasedReplicationTimeDataArrayWithDecileIndex,
                                                                        ordered_signatures_cutoffs,
                                                                        signatures_mask_array,
                                                                        signature_decile_index_accumulated_np_array,
                                                                        discreet_mode,
                                                                        default_cutoff,
                                                                        list_of_rep_timing_decile1_dict,
                                                                        list_of_rep_timing_decile2_dict,
                                                                        list_of_rep_timing_decile3_dict,
                                                                        list_of_rep_timing_decile4_dict,
                                                                        list_of_rep_timing_decile5_dict,
                                                                        list_of_rep_timing_decile6_dict,
                                                                        list_of_rep_timing_decile7_dict,
                                                                        list_of_rep_timing_decile8_dict,
                                                                        list_of_rep_timing_decile9_dict,
                                                                        list_of_rep_timing_decile10_dict,
                                                                        df_columns):


    # df_columns: numpy array
    indexofStart = np.where(df_columns == START) [0][0]
    start = mutation_row[indexofStart]

    # For single point mutations start and end are the same, therefore we need to add 1 to mutation_row[END]
    # end is exclusive for subs, indels and dinucs provided by readChrBased methods

    # Simulation Number
    index_of_sim_num = np.where(df_columns == 'Simulation_Number')[0][0]
    sim_num = mutation_row[index_of_sim_num]

    if mutation_type == SUBS:
        end = start + 1
    elif mutation_type == DINUCS:
        end = start + 2
    elif mutation_type == INDELS:
        indexofLength = np.where(df_columns == LENGTH)
        length = mutation_row[indexofLength]
        end = start + int(length)

    slicedArray = chrBasedReplicationTimeDataArrayWithDecileIndex[start:end]

    # np.nonzero returns the indices of the elements that are non-zero.
    # np.unique finds the unique elements of an array returns ndarray the sorted unique values.
    # np.nditer efficient multi-dimensional iterator object to iterate over arrays.
    # We aim to get rid of zero if any exists in slicedArray.
    uniqueIndexesArray = np.unique(slicedArray[np.nonzero(slicedArray)])

    decile_index_array = np.zeros((10), dtype=int)

    # decile 10  will be accumulated in index 9
    # decile 1 will be accumulated in index 0
    # Therefore uniqueIndexesArray minus 1
    if (uniqueIndexesArray.size > 0):
        uniqueIndexesArray -= 1
        decile_index_array[uniqueIndexesArray] = 1

        if (sim_num == 0) and (decile_index_array[0]):
            my_dict = {a: b for a, b in zip(df_columns, mutation_row)}
            list_of_rep_timing_decile1_dict.append(my_dict)
        if (sim_num == 0) and (decile_index_array[1]):
            my_dict = {a: b for a, b in zip(df_columns, mutation_row)}
            list_of_rep_timing_decile2_dict.append(my_dict)
        if (sim_num == 0) and (decile_index_array[2]):
            my_dict = {a: b for a, b in zip(df_columns, mutation_row)}
            list_of_rep_timing_decile3_dict.append(my_dict)
        if (sim_num == 0) and (decile_index_array[3]):
            my_dict = {a: b for a, b in zip(df_columns, mutation_row)}
            list_of_rep_timing_decile4_dict.append(my_dict)
        if (sim_num == 0) and (decile_index_array[4]):
            my_dict = {a: b for a, b in zip(df_columns, mutation_row)}
            list_of_rep_timing_decile5_dict.append(my_dict)
        if (sim_num == 0) and (decile_index_array[5]):
            my_dict = {a: b for a, b in zip(df_columns, mutation_row)}
            list_of_rep_timing_decile6_dict.append(my_dict)
        if (sim_num == 0) and (decile_index_array[6]):
            my_dict = {a: b for a, b in zip(df_columns, mutation_row)}
            list_of_rep_timing_decile7_dict.append(my_dict)
        if (sim_num == 0) and (decile_index_array[7]):
            my_dict = {a: b for a, b in zip(df_columns, mutation_row)}
            list_of_rep_timing_decile8_dict.append(my_dict)
        if (sim_num == 0) and (decile_index_array[8]):
            my_dict = {a: b for a, b in zip(df_columns, mutation_row)}
            list_of_rep_timing_decile9_dict.append(my_dict)
        if (sim_num == 0) and (decile_index_array[9]):
            my_dict = {a: b for a, b in zip(df_columns, mutation_row)}
            list_of_rep_timing_decile10_dict.append(my_dict)

        probabilities = mutation_row[signatures_mask_array]

        # new code starts
        if discreet_mode:
            # Discreet way 1 or 0
            # Convert True into 1, and False into 0
            threshold_mask_array = np.greater_equal(probabilities, ordered_signatures_cutoffs)
            mask_array = threshold_mask_array.astype(int)
            if mutation_type == INDELS:
                if (length >= 3):
                    # Order is important
                    # MICROHOMOLOGY --> 1
                    # REPEAT --> 0
                    # AGGREGATEDINDELS --> 1
                    mask_array = np.append(mask_array, [1, 0, 1])
                else:
                    # Order is important
                    # MICROHOMOLOGY --> 0
                    # REPEAT --> 1
                    # AGGREGATEDINDELS --> 1
                    mask_array = np.append(mask_array, [0, 1, 1])

            else:
                # Add 1 for the aggregated analysis to the mask array
                # For SUBS and DINUCS
                mask_array = np.append(mask_array, 1)

        else:
            # new way
            probabilities[probabilities < default_cutoff ] = 0
            mask_array = np.array(probabilities).astype(float)

            if mutation_type == INDELS:
                if (length >= 3):
                    # Order is important
                    # MICROHOMOLOGY --> 1
                    # REPEAT --> 0
                    # AGGREGATEDINDELS --> 1
                    mask_array = np.append(mask_array, [1.0, 0.0, 1.0])
                else:
                    # Order is important
                    # MICROHOMOLOGY --> 0
                    # REPEAT --> 1
                    # AGGREGATEDINDELS --> 1
                    mask_array = np.append(mask_array, [0.0, 1.0, 1.0])

            else:
                # Add 1 for the aggregated analysis to the mask array
                # For SUBS and DINUCS
                mask_array = np.append(mask_array, 1.0)

        # Add 1 more dimension to the arrays
        decile_index_array_1x10 = np.expand_dims(decile_index_array, axis=0)
        mask_array_1xnumofsignatures = np.expand_dims(mask_array, axis=0)

        signatures_decile_index_np_array = mask_array_1xnumofsignatures.T * decile_index_array_1x10
        signature_decile_index_accumulated_np_array += signatures_decile_index_np_array


# For df_split
# Using numpy array
def search_for_each_mutation_using_list_comprehension_using_numpy_array_for_df_split(mutation_row,
                                                                        chrBasedReplicationTimeDataArrayWithDecileIndex,
                                                                        ordered_sbs_signatures_cutoffs,
                                                                        ordered_dbs_signatures_cutoffs,
                                                                        ordered_id_signatures_cutoffs,
                                                                        df_columns_subs_signatures_mask_array,
                                                                        df_columns_dinucs_signatures_mask_array,
                                                                        df_columns_indels_signatures_mask_array,
                                                                        subs_signature_decile_index_accumulated_np_array,
                                                                        dinucs_signature_decile_index_accumulated_np_array,
                                                                        indels_signature_decile_index_accumulated_np_array,
                                                                        discreet_mode,
                                                                        df_columns):

    # df_columns: numpy array
    indexofType = np.where(df_columns == TYPE)[0][0]
    indexofSample = np.where(df_columns == SAMPLE)[0][0]
    indexofStart = np.where(df_columns == START) [0][0]
    indexofSimulationNumber = np.where(df_columns == SIMULATION_NUMBER)[0][0]

    mutation_row_type = mutation_row[indexofType]
    sample = mutation_row[indexofSample]
    start = mutation_row[indexofStart]
    mutation_row_sim_num = mutation_row[indexofSimulationNumber]

    # For single point mutations start and end are the same, therefore we need to add 1 to mutation_row[END]
    # end is exclusive for subs, indels and dinucs provided by readChrBased methods

    if mutation_row_type==SUBS:
        end = start + 1
        cutoffs=ordered_sbs_signatures_cutoffs
        signatures_mask_array=df_columns_subs_signatures_mask_array
        signature_decile_index_accumulated_np_array = subs_signature_decile_index_accumulated_np_array
    elif mutation_row_type == DINUCS:
        end = start + 2
        cutoffs=ordered_dbs_signatures_cutoffs
        signatures_mask_array=df_columns_dinucs_signatures_mask_array
        signature_decile_index_accumulated_np_array = dinucs_signature_decile_index_accumulated_np_array
    elif mutation_row_type == INDELS:
        indexofLength = np.where(df_columns == LENGTH)
        length=mutation_row[indexofLength]
        end = start + int(length)
        cutoffs=ordered_id_signatures_cutoffs
        signatures_mask_array=df_columns_indels_signatures_mask_array
        signature_decile_index_accumulated_np_array = indels_signature_decile_index_accumulated_np_array

    slicedArray = chrBasedReplicationTimeDataArrayWithDecileIndex[start:end]

    # np.nonzero returns the indices of the elements that are non-zero.
    # np.unique finds the unique elements of an array returns ndarray the sorted unique values.
    # np.nditer efficient multi-dimensional iterator object to iterate over arrays.
    # We aim to get rid of zero if any exists in slicedArray.
    uniqueIndexesArray = np.unique(slicedArray[np.nonzero(slicedArray)])

    decile_index_array=np.zeros((10),dtype=int)

    #decile 10  will be accumulated in index 9
    #decile 1 will be accumulated in index 0
    #Therefore uniqueIndexesArray minus 1
    if (uniqueIndexesArray.size>0):
        uniqueIndexesArray -=1
        decile_index_array[uniqueIndexesArray]=1

        probabilities = mutation_row[signatures_mask_array]

        if discreet_mode:
            threshold_mask_array = np.greater_equal(probabilities, cutoffs)

            # Convert True into 1, and False into 0
            mask_array = threshold_mask_array.astype(int)
            if  mutation_row_type == INDELS:
                if (length >= 3):
                    # Order is important
                    # MICROHOMOLOGY --> 1
                    # REPEAT --> 0
                    # AGGREGATEDINDELS --> 1
                    mask_array = np.append(mask_array, [1, 0, 1])
                else:
                    # Order is important
                    # MICROHOMOLOGY --> 0
                    # REPEAT --> 1
                    # AGGREGATEDINDELS --> 1
                    mask_array = np.append(mask_array, [0, 1, 1])

            else:
                # Add 1 for the aggregated analysis to the mask array
                # For SUBS and DINUCS
                mask_array = np.append(mask_array, 1)

        else:
            mask_array = np.array(probabilities).astype(float)
            if mutation_row_type == INDELS:
                if (length >= 3):
                    # Order is important
                    # MICROHOMOLOGY --> 1
                    # REPEAT --> 0
                    # AGGREGATEDINDELS --> 1
                    mask_array = np.append(mask_array, [1.0, 0.0, 1.0])
                else:
                    # Order is important
                    # MICROHOMOLOGY --> 0
                    # REPEAT --> 1
                    # AGGREGATEDINDELS --> 1
                    mask_array = np.append(mask_array, [0.0, 1.0, 1.0])

            else:
                # Add 1 for the aggregated analysis to the mask array
                # For SUBS and DINUCS
                mask_array = np.append(mask_array, 1.0)


        #Add 1 more dimension to the arrays
        mask_array_1xnumofsignatures = np.expand_dims(mask_array, axis=0)
        decile_index_array_1x10 = np.expand_dims(decile_index_array, axis=0)

        signatures_decile_index_np_array = mask_array_1xnumofsignatures.T * decile_index_array_1x10
        signature_decile_index_accumulated_np_array += signatures_decile_index_np_array


# sample based
def searchforAllMutations_sample_based_using_numpy_array(chrLong,
                                                         sim_num,
                                                         chrBased_simBased_subs_df,
                                                         chrBased_simBased_dinucs_df,
                                                         chrBased_simBased_indels_df,
                                                         chrBasedReplicationTimeDataArrayWithDecileIndex,
                                                         ordered_sbs_signatures,
                                                         ordered_dbs_signatures,
                                                         ordered_id_signatures,
                                                         ordered_sbs_signatures_cutoffs,
                                                         ordered_dbs_signatures_cutoffs,
                                                         ordered_id_signatures_cutoffs,
                                                         all_samples_list,
                                                         sample_2_sample_index_dict,
                                                         discreet_mode,
                                                         default_cutoff):

    chrBased_simBased_subs_replication_timing_list = None
    chrBased_simBased_doublets_replication_timing_list = None
    chrBased_simBased_indels_replication_timing_list = None

    number_of_sbs_signatures = ordered_sbs_signatures.size
    number_of_dbs_signatures = ordered_dbs_signatures.size
    number_of_id_signatures = ordered_id_signatures.size

    # Add one more row for the Aggregated analysis, there are 10 deciles
    # Add three more rows for the Microhomology, Repeat Mediated and Aggregated analysis, there are 10 deciles
    samples_subs_signature_decile_index_accumulated_np_array = np.zeros((len(all_samples_list), number_of_sbs_signatures + 1, 10), dtype=float) # legacy int
    samples_dinucs_signature_decile_index_accumulated_np_array = np.zeros((len(all_samples_list),number_of_dbs_signatures + 1, 10), dtype=float) # legacy int
    samples_indels_signature_decile_index_accumulated_np_array = np.zeros((len(all_samples_list),number_of_id_signatures + 3, 10), dtype=float) # # legacy int

    # SUBS
    if ((chrBased_simBased_subs_df is not None) and (not chrBased_simBased_subs_df.empty)):
        df_columns = chrBased_simBased_subs_df.columns.values

        df_columns_subs_signatures_mask_array = np.isin(df_columns, ordered_sbs_signatures)

        chrBased_simBased_subs_replication_timing_list = [search_for_each_mutation_using_list_comprehension_sample_based_using_numpy_array(mutation_row,
                                                                            SUBS,
                                                                            chrBasedReplicationTimeDataArrayWithDecileIndex,
                                                                            ordered_sbs_signatures_cutoffs,
                                                                            df_columns_subs_signatures_mask_array,
                                                                            samples_subs_signature_decile_index_accumulated_np_array,
                                                                            sample_2_sample_index_dict,
                                                                            discreet_mode,
                                                                            default_cutoff,
                                                                            df_columns) for mutation_row in chrBased_simBased_subs_df.values]

    # DINUCS
    if ((chrBased_simBased_dinucs_df is not None) and (not chrBased_simBased_dinucs_df.empty)):
        df_columns = chrBased_simBased_dinucs_df.columns.values

        df_columns_dinucs_signatures_mask_array = np.isin(df_columns, ordered_dbs_signatures)

        chrBased_simBased_doublets_replication_timing_list = [search_for_each_mutation_using_list_comprehension_sample_based_using_numpy_array(mutation_row,
                                                                            DINUCS,
                                                                            chrBasedReplicationTimeDataArrayWithDecileIndex,
                                                                            ordered_dbs_signatures_cutoffs,
                                                                            df_columns_dinucs_signatures_mask_array,
                                                                            samples_dinucs_signature_decile_index_accumulated_np_array,
                                                                            sample_2_sample_index_dict,
                                                                            discreet_mode,
                                                                            default_cutoff,
                                                                            df_columns) for mutation_row in chrBased_simBased_dinucs_df.values]

    # INDELS
    if ((chrBased_simBased_indels_df is not None) and (not chrBased_simBased_indels_df.empty)):
        df_columns = chrBased_simBased_indels_df.columns.values

        df_columns_indels_signatures_mask_array = np.isin(df_columns, ordered_id_signatures)

        chrBased_simBased_indels_replication_timing_list = [search_for_each_mutation_using_list_comprehension_sample_based_using_numpy_array(mutation_row,
                                                                                                                                             INDELS,
                                                                                                                                             chrBasedReplicationTimeDataArrayWithDecileIndex,
                                                                                                                                             ordered_id_signatures_cutoffs,
                                                                                                                                             df_columns_indels_signatures_mask_array,
                                                                                                                                             samples_indels_signature_decile_index_accumulated_np_array,
                                                                                                                                             sample_2_sample_index_dict,
                                                                                                                                             discreet_mode,
                                                                                                                                             default_cutoff,
                                                                                                                                             df_columns) for mutation_row in chrBased_simBased_indels_df.values]


    chrBased_simBased_subs_replication_timing_array = np.array(chrBased_simBased_subs_replication_timing_list)
    chrBased_simBased_doublets_replication_timing_array = np.array(chrBased_simBased_doublets_replication_timing_list)
    chrBased_simBased_indels_replication_timing_array = np.array(chrBased_simBased_indels_replication_timing_list)

    if chrBased_simBased_subs_df is not None:
        chrBased_simBased_subs_df[REPLICATIONTIMECOLUMN] = chrBased_simBased_subs_replication_timing_array

        # Insert ReplicationStrand column just before Mutation column
        columns = chrBased_simBased_subs_df.columns.tolist()
        mutation_index = columns.index(MUTATION)
        replication_time_index = columns.index(REPLICATIONTIMECOLUMN)
        if replication_time_index > mutation_index:
            columns.insert(mutation_index, columns.pop(columns.index(REPLICATIONTIMECOLUMN)))
        columns.pop(columns.index(SIMULATION_NUMBER))
        chrBased_simBased_subs_df = chrBased_simBased_subs_df[columns]

    if chrBased_simBased_dinucs_df is not None:
        chrBased_simBased_dinucs_df[REPLICATIONTIMECOLUMN] = chrBased_simBased_doublets_replication_timing_array

        # Insert ReplicationStrand column just before Mutation column
        columns = chrBased_simBased_dinucs_df.columns.tolist()
        mutation_index = columns.index(MUTATION)
        replication_time_index = columns.index(REPLICATIONTIMECOLUMN)
        if replication_time_index > mutation_index:
            columns.insert(mutation_index, columns.pop(columns.index(REPLICATIONTIMECOLUMN)))
        columns.pop(columns.index(SIMULATION_NUMBER))
        chrBased_simBased_dinucs_df = chrBased_simBased_dinucs_df[columns]

    if chrBased_simBased_indels_df is not None:
        chrBased_simBased_indels_df[REPLICATIONTIMECOLUMN] = chrBased_simBased_indels_replication_timing_array

        # Insert ReplicationStrand column just before Mutation column
        columns = chrBased_simBased_indels_df.columns.tolist()
        mutation_index = columns.index(MUTATION)
        replication_strand_index = columns.index(REPLICATIONTIMECOLUMN)
        if replication_strand_index > mutation_index:
            columns.insert(mutation_index, columns.pop(columns.index(REPLICATIONTIMECOLUMN)))
        columns.pop(columns.index(SIMULATION_NUMBER))
        chrBased_simBased_indels_df = chrBased_simBased_indels_df[columns]


    return  chrLong, \
            sim_num, \
            chrBased_simBased_subs_df, \
            chrBased_simBased_dinucs_df, \
            chrBased_simBased_indels_df, \
            samples_subs_signature_decile_index_accumulated_np_array, \
            samples_dinucs_signature_decile_index_accumulated_np_array, \
            samples_indels_signature_decile_index_accumulated_np_array


# using pyranges
def find_overlaps_core_enhanced(chrLong,
                                sim_num,
                                mutation_type,
                                repli_seq_intervals_df,
                                mutations_intervals_df,
                                ordered_signatures_cutoffs,
                                signatures_mask_array,
                                signature_num_of_mutations_in_bins_arr,
                                discreet_mode,
                                default_cutoff,
                                df_columns):

    # signatures_mask_array: [False False False False False False False  True  True False  True False
    #   True  True False False]
    signature_probabilities = mutations_intervals_df.loc [:, signatures_mask_array]
    # signature_probabilities = mutations_intervals_df.loc[:, mutations_intervals_df.columns.isin(['B', 'C'])]

    if discreet_mode:
        # Discreet way 1 or 0
        # Convert True into 1, and False into 0
        # mask array: rows -> number of mutations, columns -> number of signatures
        threshold_mask_array = np.greater_equal(signature_probabilities, ordered_signatures_cutoffs)
        signatures_mask_array = threshold_mask_array.astype(int)

        # Add 1 column for the aggregated analysis to the mask array
        # Create a new array with the desired shape, filled with ones
        # rows: number of mutations, columns: 1 column only
        aggregated_column = np.ones((signatures_mask_array.shape[0], 1), dtype=int)

        # Create a new array by horizontally stacking the original array and the new column
        signatures_aggregated_mask_array = np.hstack((signatures_mask_array, aggregated_column))

    else:
        signature_probabilities[signature_probabilities < default_cutoff] = 0
        signatures_mask_array = np.array(signature_probabilities).astype(float)

        # Add 1 for the aggregated analysis to the mask array
        # Create a new array with the desired shape, filled with ones
        aggregated_column = np.ones((signatures_mask_array.shape[0], 1.0))

        # Create a new array by horizontally stacking the original array and the new column
        signatures_aggregated_mask_array = np.hstack((signatures_mask_array, aggregated_column))

    # Add a new column with values from 0 to n-1 where n is the number of rows
    repli_seq_intervals_df['repli_seq_index'] = range(len(repli_seq_intervals_df))
    mutations_intervals_df['mutation_index'] = range(len(mutations_intervals_df))

    repli_seq_ranges = pr.PyRanges(repli_seq_intervals_df)

    mutations_intervals_df.rename(columns={'Chrom': 'Chromosome'}, inplace=True)
    mutations_intervals_df['Chromosome'] = 'chr' + mutations_intervals_df['Chromosome']
    mutations_ranges = pr.PyRanges(mutations_intervals_df)

    overlapping_intervals = repli_seq_ranges.join(mutations_ranges)

    # another way to add ids
    # add ids
    # gr.id1 = np.arange(len(gr))
    # gr2.id2 = np.arange(len(gr))

    # bins_2d_arr: rows: number of overlapping intervals
    # bins_2d_arr: columns: 10 for 10 bins
    # bins_2d_arr = np.eye(10)[repli_seq_bins[overlapping_indices[:,0]]]
    bins_2d_arr = np.eye(10)[repli_seq_intervals_df.loc[overlapping_intervals.df['repli_seq_index'].values]['Bin'].values.astype(int)]

    # signatures_aggregated_mask_array = signatures_aggregated_mask_array[overlapping_indices[:, 1]]
    signatures_aggregated_mask_array = signatures_aggregated_mask_array[overlapping_intervals.df['mutation_index'].values]

    # signatures_decile_index_np_array = mask_array.T * bins_2d_arr
    signatures_decile_index_np_array = np.dot(signatures_aggregated_mask_array.T, bins_2d_arr)
    signature_num_of_mutations_in_bins_arr += signatures_decile_index_np_array

    return signature_num_of_mutations_in_bins_arr


# using pyranges
def find_overlaps_enhanced(chrLong,
                           sim_num,
                           chrBased_simBased_subs_df,
                           mutation_type,
                           chr_based_repli_seq_array,
                           ordered_signatures,
                           ordered_signatures_cutoffs,
                           signature_sample_mutation_type_count_arr,
                           discreet_mode,
                           default_cutoff):

    if mutation_type == SUBS:
        chrBased_simBased_subs_df['End'] = chrBased_simBased_subs_df['Start'] + 1
    elif mutation_type == DINUCS:
        chrBased_simBased_subs_df['End'] = chrBased_simBased_subs_df['Start'] + 2
    elif mutation_type == INDELS:
        chrBased_simBased_subs_df['End'] = chrBased_simBased_subs_df['Start'] + chrBased_simBased_subs_df['Length']

    df_columns = chrBased_simBased_subs_df.columns.values
    df_columns_signatures_mask_array = np.isin(df_columns, ordered_signatures)

    # Sample  Chrom   Start   MutationLong    PyramidineStrand        TranscriptionStrand     Mutation        SBS1    SBS2    SBS3    SBS5    SBS8    SBS13   SBS40
    # PD3851a 1       809687  N:TC[C>G]TC     -1      N       C>G     2.6532282126713807e-15  0.0     0.0     0.1815425900499847      0.0     0.0     0.8184574099500127
    # PD3851a 1       121101678       N:TC[C>G]TA     -1      N       C>G     2.6532282126713807e-15  0.0     0.0     0.1815425900499847      0.0     0.0     0.8184574099500127
    # PD3851a 1       819245  T:GG[C>A]TA     -1      T       C>A     3.4644592534890214e-15  0.0     0.0     0.16964101113950414     0.0     0.0     0.8303589888604923
    # PD3851a 1       28195262        N:AG[C>A]TG     1       N       C>A     3.4644592534890214e-15  0.0     0.0     0.16964101113950414     0.0     0.0     0.8303589888604923

    # chr_based_mutations_arr = chrBased_simBased_subs_df.iloc[:, 1:].values
    # chr_based_mutations_arr = chrBased_simBased_subs_df.values

    chr_based_repli_seq_df = pd.DataFrame(data=chr_based_repli_seq_array,
                                         columns= ['Chromosome', 'Start', 'End', 'Signal', 'Bin'])

    return find_overlaps_core_enhanced(chrLong,
                                       sim_num,
                                       mutation_type,
                                       chr_based_repli_seq_df,
                                       chrBased_simBased_subs_df,
                                       ordered_signatures_cutoffs,
                                       df_columns_signatures_mask_array,
                                       signature_sample_mutation_type_count_arr,
                                       discreet_mode,
                                       default_cutoff,
                                       df_columns)




# using pyranges
def distribute_mutations_into_replication_timing_bins_enhanced(chrLong,
                                                      chromSize,
                                                      sim_num,
                                                      chr_based_repli_seq_array,
                                                      chrBased_simBased_subs_df,
                                                      chrBased_simBased_dinucs_df,
                                                      chrBased_simBased_indels_df,
                                                      ordered_sbs_signatures,
                                                      ordered_dbs_signatures,
                                                      ordered_id_signatures,
                                                      ordered_sbs_signatures_cutoffs,
                                                      ordered_dbs_signatures_cutoffs,
                                                      ordered_id_signatures_cutoffs,
                                                      discreet_mode,
                                                      default_cutoff):

    number_of_sbs_signatures = ordered_sbs_signatures.size
    number_of_dbs_signatures = ordered_dbs_signatures.size
    number_of_id_signatures = ordered_id_signatures.size

    # Add one more row for the Aggregated analysis, there are 10 deciles
    sbs_signature_sample_mutation_type_count_arr = np.zeros((number_of_sbs_signatures + 1, 10), dtype=float) # legacy int
    dbs_signature_sample_mutation_type_count_arr = np.zeros((number_of_dbs_signatures + 1, 10), dtype=float) # legacy int
    id_signature_sample_mutation_type_count_arr = np.zeros((number_of_id_signatures + 1, 10), dtype=float) # # legacy int

    # SUBS
    if ((chrBased_simBased_subs_df is not None) and (not chrBased_simBased_subs_df.empty)):
        sbs_signature_sample_mutation_type_count_arr = find_overlaps_enhanced(chrLong,
                                                                              sim_num,
                                                                              chrBased_simBased_subs_df,
                                                                              SUBS,
                                                                              chr_based_repli_seq_array,
                                                                              ordered_sbs_signatures,
                                                                              ordered_sbs_signatures_cutoffs,
                                                                              sbs_signature_sample_mutation_type_count_arr,
                                                                              discreet_mode,
                                                                              default_cutoff)



    # DOUBLETS
    if ((chrBased_simBased_dinucs_df is not None) and (not chrBased_simBased_dinucs_df.empty)):
        dbs_signature_sample_mutation_type_count_arr = find_overlaps_enhanced(chrLong,
                                                                              sim_num,
                                                                              chrBased_simBased_dinucs_df,
                                                                              DINUCS,
                                                                              chr_based_repli_seq_array,
                                                                              ordered_dbs_signatures,
                                                                              ordered_dbs_signatures_cutoffs,
                                                                              dbs_signature_sample_mutation_type_count_arr,
                                                                              discreet_mode,
                                                                              default_cutoff)

    # INDELS
    if ((chrBased_simBased_indels_df is not None) and (not chrBased_simBased_indels_df.empty)):
        id_signature_sample_mutation_type_count_arr = find_overlaps_enhanced(chrLong,
                                                                             sim_num,
                                                                             chrBased_simBased_indels_df,
                                                                              INDELS,
                                                                              chr_based_repli_seq_array,
                                                                              ordered_id_signatures,
                                                                              ordered_id_signatures_cutoffs,
                                                                              id_signature_sample_mutation_type_count_arr,
                                                                              discreet_mode,
                                                                              default_cutoff)

    return sim_num, \
           sbs_signature_sample_mutation_type_count_arr, \
           dbs_signature_sample_mutation_type_count_arr, \
           id_signature_sample_mutation_type_count_arr



def searchforAllMutations_using_numpy_array(sim_num,
                                            chrBased_simBased_subs_df,
                                            chrBased_simBased_dinucs_df,
                                            chrBased_simBased_indels_df,
                                            chrBasedReplicationTimeDataArrayWithDecileIndex,
                                            ordered_sbs_signatures,
                                            ordered_dbs_signatures,
                                            ordered_id_signatures,
                                            ordered_sbs_signatures_cutoffs,
                                            ordered_dbs_signatures_cutoffs,
                                            ordered_id_signatures_cutoffs,
                                            discreet_mode,
                                            default_cutoff):

    number_of_sbs_signatures = ordered_sbs_signatures.size
    number_of_dbs_signatures = ordered_dbs_signatures.size
    number_of_id_signatures = ordered_id_signatures.size

    # Add one more row for the Aggregated analysis, there are 10 deciles
    # Add three more rows for the Microhomology, Repeat Mediated and Aggregated analysis, there are 10 deciles
    subs_signature_decile_index_accumulated_np_array = np.zeros((number_of_sbs_signatures + 1, 10), dtype=float) # legacy int
    dinucs_signature_decile_index_accumulated_np_array = np.zeros((number_of_dbs_signatures + 1, 10), dtype=float) # legacy int
    indels_signature_decile_index_accumulated_np_array = np.zeros((number_of_id_signatures + 3, 10), dtype=float) # # legacy int

    list_of_rep_timing_decile1_dict = []
    list_of_rep_timing_decile2_dict = []
    list_of_rep_timing_decile3_dict = []
    list_of_rep_timing_decile4_dict = []
    list_of_rep_timing_decile5_dict = []
    list_of_rep_timing_decile6_dict = []
    list_of_rep_timing_decile7_dict = []
    list_of_rep_timing_decile8_dict = []
    list_of_rep_timing_decile9_dict = []
    list_of_rep_timing_decile10_dict = []

    # SUBS
    if ((chrBased_simBased_subs_df is not None) and (not chrBased_simBased_subs_df.empty)):
        df_columns = chrBased_simBased_subs_df.columns.values

        df_columns_subs_signatures_mask_array = np.isin(df_columns, ordered_sbs_signatures)

        [search_for_each_mutation_using_list_comprehension_using_numpy_array(mutation_row,
                                                                            SUBS,
                                                                            chrBasedReplicationTimeDataArrayWithDecileIndex,
                                                                            ordered_sbs_signatures_cutoffs,
                                                                            df_columns_subs_signatures_mask_array,
                                                                            subs_signature_decile_index_accumulated_np_array,
                                                                            discreet_mode,
                                                                            default_cutoff,
                                                                            list_of_rep_timing_decile1_dict,
                                                                            list_of_rep_timing_decile2_dict,
                                                                            list_of_rep_timing_decile3_dict,
                                                                            list_of_rep_timing_decile4_dict,
                                                                            list_of_rep_timing_decile5_dict,
                                                                            list_of_rep_timing_decile6_dict,
                                                                            list_of_rep_timing_decile7_dict,
                                                                            list_of_rep_timing_decile8_dict,
                                                                            list_of_rep_timing_decile9_dict,
                                                                            list_of_rep_timing_decile10_dict,
                                                                            df_columns) for mutation_row in chrBased_simBased_subs_df.values]

    # DINUCS
    if ((chrBased_simBased_dinucs_df is not None) and (not chrBased_simBased_dinucs_df.empty)):
        df_columns = chrBased_simBased_dinucs_df.columns.values

        df_columns_dinucs_signatures_mask_array = np.isin(df_columns, ordered_dbs_signatures)

        [search_for_each_mutation_using_list_comprehension_using_numpy_array(mutation_row,
                                                                            DINUCS,
                                                                            chrBasedReplicationTimeDataArrayWithDecileIndex,
                                                                            ordered_dbs_signatures_cutoffs,
                                                                            df_columns_dinucs_signatures_mask_array,
                                                                            dinucs_signature_decile_index_accumulated_np_array,
                                                                            discreet_mode,
                                                                            default_cutoff,
                                                                            list_of_rep_timing_decile1_dict,
                                                                            list_of_rep_timing_decile2_dict,
                                                                            list_of_rep_timing_decile3_dict,
                                                                            list_of_rep_timing_decile4_dict,
                                                                            list_of_rep_timing_decile5_dict,
                                                                            list_of_rep_timing_decile6_dict,
                                                                            list_of_rep_timing_decile7_dict,
                                                                            list_of_rep_timing_decile8_dict,
                                                                            list_of_rep_timing_decile9_dict,
                                                                            list_of_rep_timing_decile10_dict,
                                                                            df_columns) for mutation_row in chrBased_simBased_dinucs_df.values]

    # INDELS
    if ((chrBased_simBased_indels_df is not None) and (not chrBased_simBased_indels_df.empty)):
        df_columns = chrBased_simBased_indels_df.columns.values

        df_columns_indels_signatures_mask_array = np.isin(df_columns, ordered_id_signatures)

        [search_for_each_mutation_using_list_comprehension_using_numpy_array(mutation_row,
                                                                             INDELS,
                                                                            chrBasedReplicationTimeDataArrayWithDecileIndex,
                                                                            ordered_id_signatures_cutoffs,
                                                                            df_columns_indels_signatures_mask_array,
                                                                            indels_signature_decile_index_accumulated_np_array,
                                                                            discreet_mode,
                                                                            default_cutoff,
                                                                             list_of_rep_timing_decile1_dict,
                                                                             list_of_rep_timing_decile2_dict,
                                                                             list_of_rep_timing_decile3_dict,
                                                                             list_of_rep_timing_decile4_dict,
                                                                             list_of_rep_timing_decile5_dict,
                                                                             list_of_rep_timing_decile6_dict,
                                                                             list_of_rep_timing_decile7_dict,
                                                                             list_of_rep_timing_decile8_dict,
                                                                             list_of_rep_timing_decile9_dict,
                                                                             list_of_rep_timing_decile10_dict,
                                                                             df_columns) for mutation_row in chrBased_simBased_indels_df.values]

    return sim_num, \
           subs_signature_decile_index_accumulated_np_array, \
           dinucs_signature_decile_index_accumulated_np_array, \
           indels_signature_decile_index_accumulated_np_array,\
           list_of_rep_timing_decile1_dict,\
           list_of_rep_timing_decile2_dict,\
           list_of_rep_timing_decile3_dict,\
           list_of_rep_timing_decile4_dict,\
           list_of_rep_timing_decile5_dict,\
           list_of_rep_timing_decile6_dict,\
           list_of_rep_timing_decile7_dict,\
           list_of_rep_timing_decile8_dict,\
           list_of_rep_timing_decile9_dict,\
           list_of_rep_timing_decile10_dict


# For df split
# August 1, 2020
# Using numpy arrays
def searchforAllMutations_using_numpy_array_for_df_split(sim_num,
                                            chrBased_simBased_combined_df,
                                            chrBasedReplicationTimeDataArrayWithDecileIndex,
                                            ordered_sbs_signatures,
                                            ordered_dbs_signatures,
                                            ordered_id_signatures,
                                            ordered_sbs_signatures_cutoffs,
                                            ordered_dbs_signatures_cutoffs,
                                            ordered_id_signatures_cutoffs,
                                            discreet_mode):

    df_columns = chrBased_simBased_combined_df.columns.values

    ###############################################################################
    ################################ Initialization ###############################
    ###############################################################################
    df_columns_subs_signatures_mask_array = np.isin(df_columns, ordered_sbs_signatures)
    df_columns_dinucs_signatures_mask_array = np.isin(df_columns, ordered_dbs_signatures)
    df_columns_indels_signatures_mask_array = np.isin(df_columns, ordered_id_signatures)

    # Add one more row for the Aggregated analysis, there are 10 deciles
    # Add three more rows for the Microhomology, Repeat Mediated and Aggregated analysis, there are 10 deciles
    subs_signature_decile_index_accumulated_np_array = np.zeros((ordered_sbs_signatures.size + 1, 10),dtype=int)
    dinucs_signature_decile_index_accumulated_np_array = np.zeros((ordered_dbs_signatures.size + 1, 10),dtype=int)
    indels_signature_decile_index_accumulated_np_array = np.zeros((ordered_id_signatures.size + 3, 10),dtype=int)

    ###############################################################################
    ################################ Initialization ###############################
    ###############################################################################

    if ((chrBased_simBased_combined_df is not None) and (not chrBased_simBased_combined_df.empty)):

        [search_for_each_mutation_using_list_comprehension_using_numpy_array_for_df_split(mutation_row,
                                                                            chrBasedReplicationTimeDataArrayWithDecileIndex,
                                                                            ordered_sbs_signatures_cutoffs,
                                                                            ordered_dbs_signatures_cutoffs,
                                                                            ordered_id_signatures_cutoffs,
                                                                            df_columns_subs_signatures_mask_array,
                                                                            df_columns_dinucs_signatures_mask_array,
                                                                            df_columns_indels_signatures_mask_array,
                                                                            subs_signature_decile_index_accumulated_np_array,
                                                                            dinucs_signature_decile_index_accumulated_np_array,
                                                                            indels_signature_decile_index_accumulated_np_array,
                                                                            discreet_mode,
                                                                            df_columns) for mutation_row in chrBased_simBased_combined_df.values]

    return sim_num, subs_signature_decile_index_accumulated_np_array, dinucs_signature_decile_index_accumulated_np_array, indels_signature_decile_index_accumulated_np_array


# Please notice that replication time data are not overlappig data therefore only setting one decileIndex will be correct.
# e.g.:
# start end
#   10 1009
# 1010 2009
def fill_array(chrBased_replicationtimedata_row, chrBasedDecileIndexArray, decileIndex):
    start = chrBased_replicationtimedata_row[START]
    end = chrBased_replicationtimedata_row[END] + 1
    chrBasedDecileIndexArray[start:end] = decileIndex


# Explanation of this function
# This function is called for each chromosome
# For each chromosome, a numpy array is filled.
# If there is a chromosome locus with a decile index 8 let's say, in the array that locus is filled with 8.
# Decile index can be between 1-10.
# In this function, each available interval with an index is filled in the corresponding numpy array with chrom size
def fillChrBasedReplicationTimeNPArray(chrLong, chromSize, chrBased_grouped_decile_df_list):
    # We can set the starting index as 1 in builtin function enumerate
    # First chrBased_grouped_decile has index of 1
    # Last chrBased_grouped_decile has index of 10

    # int8	Byte (-128 to 127)
    chrBasedReplicationTimeDataArrayWithDecileIndex = np.zeros(chromSize, dtype=np.int8)

    # First decileIndex is 1, last decile index is 10.
    for decileIndex, chrBased_grouped_decile_df in enumerate(chrBased_grouped_decile_df_list,1):
        # Solution to keyError
        for name, chrBased_replicationtimedata_df in chrBased_grouped_decile_df:
            if (chrLong == name) and (chrBased_replicationtimedata_df is not None) and  (not chrBased_replicationtimedata_df.empty):

                # type(chrBased_decile) -> DataFrame
                chrBased_replicationtimedata_df.apply(fill_array,
                                                      chrBasedDecileIndexArray=chrBasedReplicationTimeDataArrayWithDecileIndex,
                                                      decileIndex=decileIndex,
                                                      axis=1)


    return chrBasedReplicationTimeDataArrayWithDecileIndex


# Using numpy array
def combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray_simbased_chrbased_splitbased_using_numpy_array(outputDir,
                                       jobname,
                                       chrLong,
                                       chromSize,
                                       simNum,
                                       splitIndex,
                                       chrBased_grouped_decile_df_list,
                                       ordered_sbs_signatures,
                                       ordered_dbs_signatures,
                                       ordered_id_signatures,
                                       ordered_sbs_signatures_cutoffs,
                                       ordered_dbs_signatures_cutoffs,
                                       ordered_id_signatures_cutoffs,
                                       discreet_mode):

    chrBased_simBased_combined_df_split = get_chrBased_simBased_combined_df_split(outputDir, jobname, chrLong, simNum, splitIndex)

    return combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray_using_numpy_array_for_df_split(chrLong,
                                                                                                chromSize,
                                                                                                simNum,
                                                                                                chrBased_grouped_decile_df_list,
                                                                                                chrBased_simBased_combined_df_split,
                                                                                                ordered_sbs_signatures,
                                                                                                ordered_dbs_signatures,
                                                                                                ordered_id_signatures,
                                                                                                ordered_sbs_signatures_cutoffs,
                                                                                                ordered_dbs_signatures_cutoffs,
                                                                                                ordered_id_signatures_cutoffs,
                                                                                                discreet_mode)


# sample based
# using numpy arrays
def combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray_simbased_chrbased_sample_based_using_numpy_array(outputDir,
                                                                                                           jobname,
                                                                                                           chrLong,
                                                                                                           chromSize,
                                                                                                           sim_num,
                                                                                                           samples_of_interest,
                                                                                                           chrBased_grouped_decile_df_list,
                                                                                                           ordered_sbs_signatures,
                                                                                                           ordered_dbs_signatures,
                                                                                                           ordered_id_signatures,
                                                                                                           ordered_sbs_signatures_cutoffs,
                                                                                                           ordered_dbs_signatures_cutoffs,
                                                                                                           ordered_id_signatures_cutoffs,
                                                                                                           all_samples_list,
                                                                                                           sample_2_sample_index_dict,
                                                                                                           discreet_mode,
                                                                                                           default_cutoff):

    # To reduce memory usage
    chrBased_simBased_subs_df, chrBased_simBased_dinucs_df, chrBased_simBased_indels_df = get_chrBased_simBased_dfs(outputDir, jobname, chrLong, sim_num)

    # filter chrbased_df for samples_of_interest
    if samples_of_interest is not None:
        if chrBased_simBased_subs_df is not None:
            chrBased_simBased_subs_df = chrBased_simBased_subs_df[chrBased_simBased_subs_df['Sample'].isin(samples_of_interest)]

        if chrBased_simBased_dinucs_df is not None:
            chrBased_simBased_dinucs_df = chrBased_simBased_dinucs_df[chrBased_simBased_dinucs_df['Sample'].isin(samples_of_interest)]

        if chrBased_simBased_indels_df is not None:
            chrBased_simBased_indels_df = chrBased_simBased_indels_df[chrBased_simBased_indels_df['Sample'].isin(samples_of_interest)]

    return combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray_sample_based_using_numpy_array(chrLong,
                                                                                                chromSize,
                                                                                                sim_num,
                                                                                                chrBased_grouped_decile_df_list,
                                                                                                chrBased_simBased_subs_df,
                                                                                                chrBased_simBased_dinucs_df,
                                                                                                chrBased_simBased_indels_df,
                                                                                                ordered_sbs_signatures,
                                                                                                ordered_dbs_signatures,
                                                                                                ordered_id_signatures,
                                                                                                ordered_sbs_signatures_cutoffs,
                                                                                                ordered_dbs_signatures_cutoffs,
                                                                                                ordered_id_signatures_cutoffs,
                                                                                                all_samples_list,
                                                                                                sample_2_sample_index_dict,
                                                                                                discreet_mode,
                                                                                                default_cutoff)


# using pyranges
def mediator_enhanced(outputDir,
                      jobname,
                      chrLong,
                      chromSize,
                      sim_num,
                      samples_of_interest,
                      sorted_repli_seq_arr,
                      ordered_sbs_signatures,
                      ordered_dbs_signatures,
                      ordered_id_signatures,
                      ordered_sbs_signatures_cutoffs,
                      ordered_dbs_signatures_cutoffs,
                      ordered_id_signatures_cutoffs,
                      discreet_mode,
                      default_cutoff):

    # To reduce memory usage
    chrBased_simBased_subs_df, \
    chrBased_simBased_dinucs_df, \
    chrBased_simBased_indels_df = get_chrBased_simBased_dfs(outputDir, jobname, chrLong, sim_num)

    # 0:chr, 1:start, 2:end, 3:signal, 4:bin
    chr_based_repli_seq_array = sorted_repli_seq_arr[sorted_repli_seq_arr[:, 0] == chrLong]

    # filter chrbased_df for samples_of_interest
    if samples_of_interest is not None:
        if chrBased_simBased_subs_df is not None:
            chrBased_simBased_subs_df = chrBased_simBased_subs_df[chrBased_simBased_subs_df['Sample'].isin(samples_of_interest)]

        if chrBased_simBased_dinucs_df is not None:
            chrBased_simBased_dinucs_df = chrBased_simBased_dinucs_df[chrBased_simBased_dinucs_df['Sample'].isin(samples_of_interest)]

        if chrBased_simBased_indels_df is not None:
            chrBased_simBased_indels_df = chrBased_simBased_indels_df[chrBased_simBased_indels_df['Sample'].isin(samples_of_interest)]

    return distribute_mutations_into_replication_timing_bins_enhanced(chrLong,
                                                            chromSize,
                                                            sim_num,
                                                            chr_based_repli_seq_array,
                                                            chrBased_simBased_subs_df,
                                                            chrBased_simBased_dinucs_df,
                                                            chrBased_simBased_indels_df,
                                                            ordered_sbs_signatures,
                                                            ordered_dbs_signatures,
                                                            ordered_id_signatures,
                                                            ordered_sbs_signatures_cutoffs,
                                                            ordered_dbs_signatures_cutoffs,
                                                            ordered_id_signatures_cutoffs,
                                                            discreet_mode,
                                                            default_cutoff)


# using numpy arrays
def combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray_simbased_chrbased_using_numpy_array(outputDir,
                                                                                                           jobname,
                                                                                                           chrLong,
                                                                                                           chromSize,
                                                                                                           sim_num,
                                                                                                           samples_of_interest,
                                                                                                           chrBased_grouped_decile_df_list,
                                                                                                           ordered_sbs_signatures,
                                                                                                           ordered_dbs_signatures,
                                                                                                           ordered_id_signatures,
                                                                                                           ordered_sbs_signatures_cutoffs,
                                                                                                           ordered_dbs_signatures_cutoffs,
                                                                                                           ordered_id_signatures_cutoffs,
                                                                                                           discreet_mode,
                                                                                                           default_cutoff):

    # To reduce memory usage
    chrBased_simBased_subs_df, chrBased_simBased_dinucs_df, chrBased_simBased_indels_df = get_chrBased_simBased_dfs(outputDir, jobname, chrLong, sim_num)

    # filter chrbased_df for samples_of_interest
    if samples_of_interest is not None:
        if chrBased_simBased_subs_df is not None:
            chrBased_simBased_subs_df = chrBased_simBased_subs_df[chrBased_simBased_subs_df['Sample'].isin(samples_of_interest)]

        if chrBased_simBased_dinucs_df is not None:
            chrBased_simBased_dinucs_df = chrBased_simBased_dinucs_df[chrBased_simBased_dinucs_df['Sample'].isin(samples_of_interest)]

        if chrBased_simBased_indels_df is not None:
            chrBased_simBased_indels_df = chrBased_simBased_indels_df[chrBased_simBased_indels_df['Sample'].isin(samples_of_interest)]

    return combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray_using_numpy_array(chrLong,
                                                                                                chromSize,
                                                                                                sim_num,
                                                                                                chrBased_grouped_decile_df_list,
                                                                                                chrBased_simBased_subs_df,
                                                                                                chrBased_simBased_dinucs_df,
                                                                                                chrBased_simBased_indels_df,
                                                                                                ordered_sbs_signatures,
                                                                                                ordered_dbs_signatures,
                                                                                                ordered_id_signatures,
                                                                                                ordered_sbs_signatures_cutoffs,
                                                                                                ordered_dbs_signatures_cutoffs,
                                                                                                ordered_id_signatures_cutoffs,
                                                                                                discreet_mode,
                                                                                                default_cutoff)

# sample_based
def combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray_sample_based_using_numpy_array(chrLong,
                                                                                        chromSize,
                                                                                        sim_num,
                                                                                        chrBased_grouped_decile_df_list,
                                                                                        chrBased_simBased_subs_df,
                                                                                        chrBased_simBased_dinucs_df,
                                                                                        chrBased_simBased_indels_df,
                                                                                        ordered_sbs_signatures,
                                                                                        ordered_dbs_signatures,
                                                                                        ordered_id_signatures,
                                                                                        ordered_sbs_signatures_cutoffs,
                                                                                        ordered_dbs_signatures_cutoffs,
                                                                                        ordered_id_signatures_cutoffs,
                                                                                        all_samples_list,
                                                                                        sample_2_sample_index_dict,
                                                                                        discreet_mode,
                                                                                        default_cutoff):

    # Fill replication time numpy array
    chrBasedReplicationTimeDataArrayWithDecileIndex = fillChrBasedReplicationTimeNPArray(chrLong,
                                                                                         chromSize,
                                                                                         chrBased_grouped_decile_df_list)

    return searchforAllMutations_sample_based_using_numpy_array(chrLong,
                                                                sim_num,
                                                                chrBased_simBased_subs_df,
                                                                chrBased_simBased_dinucs_df,
                                                                chrBased_simBased_indels_df,
                                                                chrBasedReplicationTimeDataArrayWithDecileIndex,
                                                                ordered_sbs_signatures,
                                                                ordered_dbs_signatures,
                                                                ordered_id_signatures,
                                                                ordered_sbs_signatures_cutoffs,
                                                                ordered_dbs_signatures_cutoffs,
                                                                ordered_id_signatures_cutoffs,
                                                                all_samples_list,
                                                                sample_2_sample_index_dict,
                                                                discreet_mode,
                                                                default_cutoff)




# using numpy arrays
def combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray_using_numpy_array(chrLong,
                                                                                        chromSize,
                                                                                        sim_num,
                                                                                        chrBased_grouped_decile_df_list,
                                                                                        chrBased_simBased_subs_df,
                                                                                        chrBased_simBased_dinucs_df,
                                                                                        chrBased_simBased_indels_df,
                                                                                        ordered_sbs_signatures,
                                                                                        ordered_dbs_signatures,
                                                                                        ordered_id_signatures,
                                                                                        ordered_sbs_signatures_cutoffs,
                                                                                        ordered_dbs_signatures_cutoffs,
                                                                                        ordered_id_signatures_cutoffs,
                                                                                        discreet_mode,
                                                                                        default_cutoff):


    # Fill replication time numpy array
    chrBasedReplicationTimeDataArrayWithDecileIndex = fillChrBasedReplicationTimeNPArray(chrLong, chromSize, chrBased_grouped_decile_df_list)

    return searchforAllMutations_using_numpy_array(sim_num,
                                                   chrBased_simBased_subs_df,
                                                   chrBased_simBased_dinucs_df,
                                                   chrBased_simBased_indels_df,
                                                   chrBasedReplicationTimeDataArrayWithDecileIndex,
                                                   ordered_sbs_signatures,
                                                   ordered_dbs_signatures,
                                                   ordered_id_signatures,
                                                   ordered_sbs_signatures_cutoffs,
                                                   ordered_dbs_signatures_cutoffs,
                                                   ordered_id_signatures_cutoffs,
                                                   discreet_mode,
                                                   default_cutoff)


# For df_split
# August 1, 2020
# using numpy arrays
def combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray_using_numpy_array_for_df_split(chrLong,
                                                                                        chromSize,
                                                                                        sim_num,
                                                                                        chrBased_grouped_decile_df_list,
                                                                                        chrBased_simBased_combined_df,
                                                                                        ordered_sbs_signatures,
                                                                                        ordered_dbs_signatures,
                                                                                        ordered_id_signatures,
                                                                                        ordered_sbs_signatures_cutoffs,
                                                                                        ordered_dbs_signatures_cutoffs,
                                                                                        ordered_id_signatures_cutoffs,
                                                                                        discreet_mode):

    # Fill replication time numpy array
    chrBasedReplicationTimeDataArrayWithDecileIndex = fillChrBasedReplicationTimeNPArray(chrLong, chromSize, chrBased_grouped_decile_df_list)

    return searchforAllMutations_using_numpy_array_for_df_split(sim_num,
                                                   chrBased_simBased_combined_df,
                                                   chrBasedReplicationTimeDataArrayWithDecileIndex,
                                                   ordered_sbs_signatures,
                                                   ordered_dbs_signatures,
                                                   ordered_id_signatures,
                                                   ordered_sbs_signatures_cutoffs,
                                                   ordered_dbs_signatures_cutoffs,
                                                   ordered_id_signatures_cutoffs,
                                                   discreet_mode)


def getNormalizedMutationDensityList(mutationDensityDict):
    densities = mutationDensityDict.values()
    maxDensity = max(densities)
    if maxDensity>0:
        normalizedMutationDensities = [x/maxDensity for x in densities]
    else:
        normalizedMutationDensities = densities

    return normalizedMutationDensities



# March 22, 2019 starts
def getMutationDensityDictNewVersion(decileIndex2NumberofAttributableBasesDict,decileBasedAllChrAccumulatedCountDict):
    decileBasedMutationDensityDict = {}
    numberofMutations = 0

    numberofMutationsList = []
    numberofAttributableBasesList = []

    for decileIndex in decileIndex2NumberofAttributableBasesDict:
        if (decileIndex in decileBasedAllChrAccumulatedCountDict):
            count = decileBasedAllChrAccumulatedCountDict[decileIndex]
            numberofMutations += count
            numofAttBases = decileIndex2NumberofAttributableBasesDict[decileIndex]
            mutationDensity = float(count) / numofAttBases
            decileBasedMutationDensityDict[decileIndex] = mutationDensity
            numberofMutationsList.append(count)
            numberofAttributableBasesList.append(numofAttBases)

        else:
            decileBasedMutationDensityDict[decileIndex] = 0
            numberofMutationsList.append(0)
            numberofAttributableBasesList.append(0)


    return numberofMutations, decileBasedMutationDensityDict, numberofMutationsList, numberofAttributableBasesList

# Feb5, 2021
def getNumberofAttributableBases(decile_df_list):
    numberofAttributableBasesList = []

    # Modifiled as enumerate(deciles,1) formerly it was enumerate(deciles,0)
    for i, decile_df in enumerate(decile_df_list,1):
        numofAttBases = decile_df[NUMOFBASES].sum()
        numberofAttributableBasesList.append(numofAttBases)

    return numberofAttributableBasesList


# num_of_attributable_bases_list
def getMutationDensityDictUsingNumpyArray_enhanced(num_of_attributable_bases_list, decile_counts_np_array):
    numberofMutations = 0
    decileBasedMutationDensityDict = {}
    numberofMutationsList = []

    # Modifiled as enumerate(deciles,1) formerly it was enumerate(deciles,0)
    for i, numofAttBases in enumerate(num_of_attributable_bases_list, 1):
        count = decile_counts_np_array[i-1]
        numberofMutations += count
        mutationDensity = float(count)/numofAttBases
        decileBasedMutationDensityDict[i] = mutationDensity
        numberofMutationsList.append(count)

        # decileBasedMutationDensityDict[i] = 0
        # numberofMutationsList.append(0)
        # numberofAttributableBasesList.append(0)
    return numberofMutations, decileBasedMutationDensityDict, numberofMutationsList


# Feb5, 2021
def getMutationDensityDictUsingNumpyArray(decile_df_list, decile_counts_np_array):
    numberofMutations = 0
    decileBasedMutationDensityDict = {}
    numberofMutationsList = []

    # Modifiled as enumerate(deciles,1) formerly it was enumerate(deciles,0)
    for i, decile_df in enumerate(decile_df_list,1):
        count = decile_counts_np_array[i-1]
        numberofMutations += count
        numofAttBases = decile_df[NUMOFBASES].sum()
        mutationDensity = float(count)/numofAttBases
        decileBasedMutationDensityDict[i] = mutationDensity
        numberofMutationsList.append(count)

        # decileBasedMutationDensityDict[i] = 0
        # numberofMutationsList.append(0)
        # numberofAttributableBasesList.append(0)
    return numberofMutations, decileBasedMutationDensityDict, numberofMutationsList

def getMutationDensityDict(decile_df_list,decileBasedAllChrAccumulatedCountDict):
    decileBasedMutationDensityDict = {}
    numberofMutations = 0

    numberofMutationsList = []
    numberofAttributableBasesList = []

    #Modifiled as enumerate(deciles,1) formerly it was enumerate(deciles,0)
    for i,decile_df in enumerate(decile_df_list,1):
        if (i in decileBasedAllChrAccumulatedCountDict):
            count = decileBasedAllChrAccumulatedCountDict[i]
            numberofMutations += count
            numofAttBases = decile_df[NUMOFBASES].sum()
            mutationDensity=float(count)/numofAttBases
            decileBasedMutationDensityDict[i] = mutationDensity
            numberofMutationsList.append(count)
            numberofAttributableBasesList.append(numofAttBases)
        else:
            decileBasedMutationDensityDict[i] = 0
            numberofMutationsList.append(0)
            numberofAttributableBasesList.append(0)

        # print('decile: %d numofAttBases: %d' %(i,numofAttBases))

    return numberofMutations, decileBasedMutationDensityDict, numberofMutationsList, numberofAttributableBasesList


# sample based
def calculateCountsForMutationsFillingReplicationTimeNPArrayRuntime_sample_based_using_numpy_array(computationType,
        outputDir,
        jobname,
        numofSimulations,
        samples_of_interest,
        all_samples_list,
        sample_2_sample_index_dict,
        job_tuples,
        chromSizesDict,
        chromNamesList,
        chrBased_grouped_decile_df_list,
        ordered_sbs_signatures,
        ordered_dbs_signatures,
        ordered_id_signatures,
        ordered_sbs_signatures_cutoffs,
        ordered_dbs_signatures_cutoffs,
        ordered_id_signatures_cutoffs,
        discreet_mode,
        default_cutoff,
        parallel_mode):

    number_of_sbs_signatures = ordered_sbs_signatures.size
    number_of_dbs_signatures = ordered_dbs_signatures.size
    number_of_id_signatures = ordered_id_signatures.size

    # Samples are considered
    # +1 for Aggregated
    all_sims_samples_subs_signature_decile_index_accumulated_np_array=np.zeros((numofSimulations+1, len(all_samples_list), number_of_sbs_signatures+1, 10), dtype=float) # int
    all_sims_samples_dinucs_signature_decile_index_accumulated_np_array=np.zeros((numofSimulations+1, len(all_samples_list), number_of_dbs_signatures+1, 10), dtype=float) # int
    # +3 for Microhomology, Repeat, Aggregated
    all_sims_samples_indels_signature_decile_index_accumulated_np_array=np.zeros((numofSimulations+1, len(all_samples_list), number_of_id_signatures+3, 10), dtype=float) # int

    def accumulate_np_arrays(result_tuple):
        chrLong = result_tuple[0]
        sim_num = result_tuple[1]
        chrBased_simBased_subs_df = result_tuple[2]
        chrBased_simBased_dinucs_df = result_tuple[3]
        chrBased_simBased_indels_df = result_tuple[4]

        # save chrom based sim based files
        if chrBased_simBased_subs_df is not None:
            write_chr_based_mutations_df(outputDir, jobname, chrLong, SBS, sim_num, chrBased_simBased_subs_df)

        if chrBased_simBased_dinucs_df is not None:
            write_chr_based_mutations_df(outputDir, jobname, chrLong, DBS, sim_num, chrBased_simBased_dinucs_df)

        if chrBased_simBased_indels_df is not None:
            write_chr_based_mutations_df(outputDir, jobname, chrLong, ID, sim_num, chrBased_simBased_indels_df)

        samples_subs_signature_decile_index_accumulated_np_array = result_tuple[5]
        samples_dinucs_signature_decile_index_accumulated_np_array = result_tuple[6]
        samples_indels_signature_decile_index_accumulated_np_array = result_tuple[7]

        # print('MONITOR ACCUMULATE', flush=True)

        # accumulate
        all_sims_samples_subs_signature_decile_index_accumulated_np_array[sim_num] += samples_subs_signature_decile_index_accumulated_np_array
        all_sims_samples_dinucs_signature_decile_index_accumulated_np_array[sim_num] += samples_dinucs_signature_decile_index_accumulated_np_array
        all_sims_samples_indels_signature_decile_index_accumulated_np_array[sim_num] += samples_indels_signature_decile_index_accumulated_np_array


    sim_nums = range(0, numofSimulations + 1)
    sim_num_chr_tuples = ((sim_num, chrLong) for sim_num in sim_nums for chrLong in chromNamesList)

    if parallel_mode:
        jobs = []

        # Using numpy array
        if (computationType == USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM):
            numofProcesses = multiprocessing.cpu_count()
            pool = multiprocessing.Pool(processes=numofProcesses)

            for simNum, chrLong in sim_num_chr_tuples:
                chromSize = chromSizesDict[chrLong]

                jobs.append(pool.apply_async(combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray_simbased_chrbased_sample_based_using_numpy_array,
                                     args=(outputDir,
                                           jobname,
                                           chrLong,
                                           chromSize,
                                           simNum,
                                           samples_of_interest,
                                           chrBased_grouped_decile_df_list,
                                           ordered_sbs_signatures,
                                           ordered_dbs_signatures,
                                           ordered_id_signatures,
                                           ordered_sbs_signatures_cutoffs,
                                           ordered_dbs_signatures_cutoffs,
                                           ordered_id_signatures_cutoffs,
                                           all_samples_list,
                                           sample_2_sample_index_dict,
                                           discreet_mode,
                                           default_cutoff,),
                                     callback=accumulate_np_arrays))

            pool.close()
            pool.join()

    else:
        # Sequential mode for testing, debugging and profiling purposes
        for simNum, chrLong in sim_num_chr_tuples:
            chromSize = chromSizesDict[chrLong]
            result_tuple = combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray_simbased_chrbased_sample_based_using_numpy_array(outputDir,
                      jobname,
                      chrLong,
                      chromSize,
                      simNum,
                      samples_of_interest,
                      chrBased_grouped_decile_df_list,
                      ordered_sbs_signatures,
                      ordered_dbs_signatures,
                      ordered_id_signatures,
                      ordered_sbs_signatures_cutoffs,
                      ordered_dbs_signatures_cutoffs,
                      ordered_id_signatures_cutoffs,
                      all_samples_list,
                      sample_2_sample_index_dict,
                      discreet_mode,
                      default_cutoff)

            accumulate_np_arrays(result_tuple)

    return all_sims_samples_subs_signature_decile_index_accumulated_np_array, \
           all_sims_samples_dinucs_signature_decile_index_accumulated_np_array, \
           all_sims_samples_indels_signature_decile_index_accumulated_np_array



# Using numpy array
def calculateCountsForMutationsFillingReplicationTimeNPArrayRuntime_using_numpy_array(computationType,
        outputDir,
        jobname,
        numofSimulations,
        samples_of_interest,
        job_tuples,
        chromSizesDict,
        chromNamesList,
        chrBased_grouped_decile_df_list,
        ordered_sbs_signatures,
        ordered_dbs_signatures,
        ordered_id_signatures,
        ordered_sbs_signatures_cutoffs,
        ordered_dbs_signatures_cutoffs,
        ordered_id_signatures_cutoffs,
        discreet_mode,
        default_cutoff,
        parallel_mode):

    number_of_sbs_signatures = ordered_sbs_signatures.size
    number_of_dbs_signatures = ordered_dbs_signatures.size
    number_of_id_signatures = ordered_id_signatures.size

    # +1 for Aggregated
    all_sims_subs_signature_decile_index_accumulated_np_array=np.zeros((numofSimulations+1, number_of_sbs_signatures+1, 10), dtype=float) # int
    all_sims_dinucs_signature_decile_index_accumulated_np_array=np.zeros((numofSimulations+1, number_of_dbs_signatures+1, 10), dtype=float) # int
    # +3 for Microhomology, Repeat, Aggregated
    all_sims_indels_signature_decile_index_accumulated_np_array=np.zeros((numofSimulations+1, number_of_id_signatures+3, 10), dtype=float) # int

    list_of_all_decile1_dicts = []
    list_of_all_decile2_dicts = []
    list_of_all_decile3_dicts = []
    list_of_all_decile4_dicts = []
    list_of_all_decile5_dicts = []
    list_of_all_decile6_dicts = []
    list_of_all_decile7_dicts = []
    list_of_all_decile8_dicts = []
    list_of_all_decile9_dicts = []
    list_of_all_decile10_dicts = []

    def accumulate_np_arrays(result_tuple):
        sim_num = result_tuple[0]
        subs_signature_decile_index_accumulated_np_array = result_tuple[1]
        dinucs_signature_decile_index_accumulated_np_array = result_tuple[2]
        indels_signature_decile_index_accumulated_np_array = result_tuple[3]

        # for mutation annotation
        list_of_rep_timing_decile1_dict = result_tuple[4]
        list_of_rep_timing_decile2_dict = result_tuple[5]
        list_of_rep_timing_decile3_dict = result_tuple[6]
        list_of_rep_timing_decile4_dict = result_tuple[7]
        list_of_rep_timing_decile5_dict = result_tuple[8]
        list_of_rep_timing_decile6_dict = result_tuple[9]
        list_of_rep_timing_decile7_dict  = result_tuple[10]
        list_of_rep_timing_decile8_dict = result_tuple[11]
        list_of_rep_timing_decile9_dict = result_tuple[12]
        list_of_rep_timing_decile10_dict = result_tuple[13]

        # print('MONITOR ACCUMULATE', flush=True)

        # accumulate
        all_sims_subs_signature_decile_index_accumulated_np_array[sim_num] += subs_signature_decile_index_accumulated_np_array
        all_sims_dinucs_signature_decile_index_accumulated_np_array[sim_num] += dinucs_signature_decile_index_accumulated_np_array
        all_sims_indels_signature_decile_index_accumulated_np_array[sim_num] += indels_signature_decile_index_accumulated_np_array

        list_of_all_decile1_dicts.extend(list_of_rep_timing_decile1_dict)
        list_of_all_decile2_dicts.extend(list_of_rep_timing_decile2_dict)
        list_of_all_decile3_dicts.extend(list_of_rep_timing_decile3_dict)
        list_of_all_decile4_dicts.extend(list_of_rep_timing_decile4_dict)
        list_of_all_decile5_dicts.extend(list_of_rep_timing_decile5_dict)
        list_of_all_decile6_dicts.extend(list_of_rep_timing_decile6_dict)
        list_of_all_decile7_dicts.extend(list_of_rep_timing_decile7_dict)
        list_of_all_decile8_dicts.extend(list_of_rep_timing_decile8_dict)
        list_of_all_decile9_dicts.extend(list_of_rep_timing_decile9_dict)
        list_of_all_decile10_dicts.extend(list_of_rep_timing_decile10_dict)

    sim_nums = range(0, numofSimulations + 1)
    sim_num_chr_tuples = ((sim_num, chrLong) for sim_num in sim_nums for chrLong in chromNamesList)

    if parallel_mode:
        jobs = []

        # Using numpy array
        if (computationType == USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM):
            numofProcesses = multiprocessing.cpu_count()
            pool = multiprocessing.Pool(processes=numofProcesses)

            for simNum, chrLong in sim_num_chr_tuples:
                chromSize = chromSizesDict[chrLong]

                jobs.append(pool.apply_async(combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray_simbased_chrbased_using_numpy_array,
                                     args=(outputDir,
                                           jobname,
                                           chrLong,
                                           chromSize,
                                           simNum,
                                           samples_of_interest,
                                           chrBased_grouped_decile_df_list,
                                           ordered_sbs_signatures,
                                           ordered_dbs_signatures,
                                           ordered_id_signatures,
                                           ordered_sbs_signatures_cutoffs,
                                           ordered_dbs_signatures_cutoffs,
                                           ordered_id_signatures_cutoffs,
                                           discreet_mode,
                                           default_cutoff,),
                                     callback=accumulate_np_arrays))

            pool.close()
            pool.join()

        elif (computationType == USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT):

            numofProcesses = multiprocessing.cpu_count()
            pool = multiprocessing.Pool(processes=numofProcesses)

            for chrLong, simNum, splitIndex in job_tuples:
                chromSize = chromSizesDict[chrLong]

                jobs.append(pool.apply_async(combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray_simbased_chrbased_splitbased_using_numpy_array,
                                     args=(outputDir,
                                           jobname,
                                           chrLong,
                                           chromSize,
                                           simNum,
                                           splitIndex,
                                           chrBased_grouped_decile_df_list,
                                           ordered_sbs_signatures,
                                           ordered_dbs_signatures,
                                           ordered_id_signatures,
                                           ordered_sbs_signatures_cutoffs,
                                           ordered_dbs_signatures_cutoffs,
                                           ordered_id_signatures_cutoffs,
                                           discreet_mode,),
                                     callback=accumulate_np_arrays))

            pool.close()
            pool.join()

    else:
        # Sequential mode for testing, debugging and profiling purposes
        for simNum, chrLong in sim_num_chr_tuples:
            chromSize = chromSizesDict[chrLong]
            result_tuple = combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray_simbased_chrbased_using_numpy_array(outputDir,
                      jobname,
                      chrLong,
                      chromSize,
                      simNum,
                      samples_of_interest,
                      chrBased_grouped_decile_df_list,
                      ordered_sbs_signatures,
                      ordered_dbs_signatures,
                      ordered_id_signatures,
                      ordered_sbs_signatures_cutoffs,
                      ordered_dbs_signatures_cutoffs,
                      ordered_id_signatures_cutoffs,
                      discreet_mode,
                      default_cutoff)

            accumulate_np_arrays(result_tuple)

    return all_sims_subs_signature_decile_index_accumulated_np_array, \
           all_sims_dinucs_signature_decile_index_accumulated_np_array, \
           all_sims_indels_signature_decile_index_accumulated_np_array,\
           list_of_all_decile1_dicts,\
           list_of_all_decile2_dicts,\
           list_of_all_decile3_dicts,\
           list_of_all_decile4_dicts,\
           list_of_all_decile5_dicts,\
           list_of_all_decile6_dicts,\
           list_of_all_decile7_dicts,\
           list_of_all_decile8_dicts,\
           list_of_all_decile9_dicts,\
           list_of_all_decile10_dicts



# Augment wavelet_processed_df with numberofAttributableBases
def augment(genome, wavelet_processed_df, matrix_generator_path, log_file, verbose):

    # numofProcesses = multiprocessing.cpu_count()
    # pool = multiprocessing.Pool(numofProcesses)

    # Augment for each chromosome
    chrBased_wavelet_processed_df_groups = wavelet_processed_df.groupby(CHROM)

    frames = []

    # tuple contains (chrLong,chrBased_wavelet_processed_df_group)
    def accumulate_apply_async_result(tuple):
        chrBased_augmented_df = tuple[1]
        frames.append(chrBased_augmented_df)

    for chrLong, chrBased_wavelet_processed_df_group in chrBased_wavelet_processed_df_groups:
        chrShort = getChrShort(chrLong)
        chrbased_filename = chrShort + ".txt"
        chrbased_file_path = os.path.join(matrix_generator_path, 'references', 'chromosomes', 'tsb', genome, chrbased_filename)
        if os.path.exists(chrbased_file_path):
            accumulate_apply_async_result(addNumofAttributableBasesColumnForApplyAsync(chrLong,
                                                         chrBased_wavelet_processed_df_group,
                                                         chrbased_file_path,
                                                         log_file,
                                                         verbose))

            # pool.apply_async(addNumofAttributableBasesColumnForApplyAsync, (chrLong,
            #                                                                 chrBased_wavelet_processed_df_group,
            #                                                                 chrbased_file_path,
            #                                                                 log_file,
            #                                                                 verbose),
            #                  callback=accumulate_apply_async_result)

    # pool.close()
    # pool.join()

    augment_df = pd.concat(frames, ignore_index=True)

    return augment_df




def write_mutations_decile_based_replication_timing(outputDir,
                                                    jobname,
                                                    list_of_all_decile1_dicts,
                                                    list_of_all_decile2_dicts,
                                                    list_of_all_decile3_dicts,
                                                    list_of_all_decile4_dicts,
                                                    list_of_all_decile5_dicts,
                                                    list_of_all_decile6_dicts,
                                                    list_of_all_decile7_dicts,
                                                    list_of_all_decile8_dicts,
                                                    list_of_all_decile9_dicts,
                                                    list_of_all_decile10_dicts):

    os.makedirs(os.path.join(outputDir, jobname, DATA, REPLICATIONTIME), exist_ok=True)

    for i in range(1,11):
        filename = 'Mutations_decile%d_replicating_timing.txt' %i

        if i == 1:
            list_of_all_decile_dicts = list_of_all_decile1_dicts
        elif i == 2:
            list_of_all_decile_dicts = list_of_all_decile2_dicts
        elif i == 3:
            list_of_all_decile_dicts = list_of_all_decile3_dicts
        elif i == 4:
            list_of_all_decile_dicts = list_of_all_decile4_dicts
        elif i == 5:
            list_of_all_decile_dicts = list_of_all_decile5_dicts
        elif i == 6:
            list_of_all_decile_dicts = list_of_all_decile6_dicts
        elif i == 7:
            list_of_all_decile_dicts = list_of_all_decile7_dicts
        elif i == 8:
            list_of_all_decile_dicts = list_of_all_decile8_dicts
        elif i == 9:
            list_of_all_decile_dicts = list_of_all_decile9_dicts
        elif i == 10:
            list_of_all_decile_dicts = list_of_all_decile10_dicts

        mutations_decile1_replicating_timing_df = pd.DataFrame(list_of_all_decile_dicts)
        mutations_decile1_replicating_timing_df.to_csv(os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, filename),
                                                       sep='\t', header=True, index=False)


def get_chrom_based_number_of_attributable_bases_in_each_bin_enhanced(chrom_based_byte_array, repli_seq_intervals):

    # Initialize an empty 2D array of shape (1, 10)
    bin_num_of_attributable_bases_array_2d_arr = np.zeros((1, 10), dtype=int)

    # Get the unique bin numbers in the 4th column
    # 0:chr, 1:start, 2:end, 3:signal, 4:bin
    bin_numbers = np.unique(repli_seq_intervals[:, 4])

    for bin_number in bin_numbers:
        bin_based_repli_seq_intervals = repli_seq_intervals[repli_seq_intervals[:,4] == bin_number]
        all_indexes = np.hstack([np.arange(start, end+1) for start, end in zip(bin_based_repli_seq_intervals[:,1],
                                                                                 bin_based_repli_seq_intervals[:,2])])

        # 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 means As, Cs, Gs, and Ts
        # 16 17 18 19 means Ns
        num_of_attributable_bases = np.sum(chrom_based_byte_array[all_indexes] < 16)
        # Initialize an empty 1D array of shape (1, 10)
        bin_num_of_attributable_bases_array_2d_arr[0, bin_number] = num_of_attributable_bases

    # Get the column indices and values from the 2D array
    # column_indices = bin_num_of_attributable_bases_array[:, 0]  # bin index
    # values = bin_num_of_attributable_bases_array[:, 1]  #  number of attributable bases

    # Accumulate values into the corresponding columns using advanced indexing and broadcasting
    # bin_num_of_attributable_bases_array_2d_arr[0, column_indices] += values

    # Accumulate values into the corresponding columns using advanced indexing and broadcasting
    # np.add.at(bin_num_of_attributable_bases_array_2d_arr, (0, column_indices), values)

    return bin_num_of_attributable_bases_array_2d_arr


def get_number_of_attributable_bases_in_each_bin_enhanced(genome, sorted_repli_seq_arr, matrix_generator_path):
    # 0: chr, 1:start, 2:end, 3:signal 4:bin
    # sorted_repli_seq_arr

    # Get the unique values in the first column
    chromosomes = np.unique(sorted_repli_seq_arr[:,0])

    # Initialize an empty 1D array of shape (1, 10)
    all_chroms_bin_num_of_attributable_bases_array_2d_arr = np.zeros((1, 10), dtype=int)

    for chrLong in chromosomes:
        chrom_based_repli_seq_intervals = sorted_repli_seq_arr[sorted_repli_seq_arr[:,0] == chrLong]

        chrShort = getChrShort(chrLong)
        chrbased_filename = chrShort + ".txt"
        chrbased_file_path = os.path.join(matrix_generator_path, 'references', 'chromosomes', 'tsb', genome, chrbased_filename)

        with open(chrbased_file_path, "rb") as f2:
            # chrom_based_bytes_object is of class bytes
            chrom_based_bytes_object = f2.read()

        # chrom_based_byte_array is of class bytes
        # contains only integers from 0 to 19
        # have look at tsb_ref
        # 16: ["N", "N"],
        # 17: ["T", "N"],
        # 18: ["U", "N"],
        # 19: ["B", "N"],
        # convert bytes to integers
        chrom_based_byte_array = np.frombuffer(chrom_based_bytes_object, dtype=np.int8)  # Convert bytes to integers
        bin_num_of_attributable_bases_array_2d_arr = get_chrom_based_number_of_attributable_bases_in_each_bin_enhanced(
            chrom_based_byte_array,
            chrom_based_repli_seq_intervals)
        all_chroms_bin_num_of_attributable_bases_array_2d_arr += bin_num_of_attributable_bases_array_2d_arr

    return all_chroms_bin_num_of_attributable_bases_array_2d_arr

def writeReplicationTimeDataUsingNumpyArray_enhanced(outputDir,
                                            jobname,
                                            all_chroms_bin_num_of_attributable_bases_array_2d_arr,
                                            subs_signatures,
                                            dinucs_signatures,
                                            indels_signatures,
                                            all_sims_subs_signature_decile_index_accumulated_np_array,
                                            all_sims_dinucs_signature_decile_index_accumulated_np_array,
                                            all_sims_indels_signature_decile_index_accumulated_np_array):

    num_of_attributable_bases_list = all_chroms_bin_num_of_attributable_bases_array_2d_arr.flatten().tolist()

    os.makedirs(os.path.join(outputDir, jobname, DATA, REPLICATIONTIME), exist_ok=True)

    my_list=[(SUBS, subs_signatures,all_sims_subs_signature_decile_index_accumulated_np_array),
             (DINUCS, dinucs_signatures, all_sims_dinucs_signature_decile_index_accumulated_np_array),
             (INDELS, indels_signatures, all_sims_indels_signature_decile_index_accumulated_np_array)]

    # Write Number of Attributable Bases List
    if (all_chroms_bin_num_of_attributable_bases_array_2d_arr is not None):
        numberofAttributabelBasesFilename = 'NumberofAttributableBases.txt'
        numberofAttributabelBasesFilePath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, numberofAttributabelBasesFilename)

        # Number of Attributable Bases
        with open(numberofAttributabelBasesFilePath, 'w') as file:
            for numberofAttributabelBases in num_of_attributable_bases_list:
                file.write(str(numberofAttributabelBases) + ' ')
            file.write('\n')

    for my_tuple in my_list:
        my_type, signatures, all_sims_signature_decile_index_accumulated_np_array = my_tuple
        num_of_sims, num_of_signatures_with_extras, num_of_deciles = all_sims_signature_decile_index_accumulated_np_array.shape

        for sim_index in range(0,num_of_sims):
            for signature_index in range(0,num_of_signatures_with_extras):
                decile_counts_np_array = all_sims_signature_decile_index_accumulated_np_array[sim_index,signature_index]

                # Order is important
                # -1 contains AGGREGATEDSUBSTITUTIONS, AGGREGATEDDINUCS, AGGREGATEDINDELS for my_type==SUBS, my_type==DINUCS, my_type==INDELS respectively
                # -2 contains REPEAT for  my_type==INDELS
                # -3 contains MICROHOMOLOGY for  my_type==INDELS
                if signature_index < signatures.size:
                    signature = signatures[signature_index]
                elif signature_index == (num_of_signatures_with_extras-1) and my_type == SUBS:
                    signature = AGGREGATEDSUBSTITUTIONS
                elif signature_index == (num_of_signatures_with_extras-1) and my_type == DINUCS:
                    signature = AGGREGATEDDINUCS
                elif signature_index == (num_of_signatures_with_extras-1) and my_type == INDELS:
                    signature = AGGREGATEDINDELS

                if (sim_index==0):
                    normalizedMutationDensityFilename = '%s_NormalizedMutationDensity.txt' %(signature)
                    numberofMutationsFilename = '%s_NumberofMutations.txt' %(signature)
                else:
                    normalizedMutationDensityFilename = '%s_sim%d_NormalizedMutationDensity.txt' %(signature,sim_index)
                    numberofMutationsFilename = '%s_sim%d_NumberofMutations.txt' %(signature,sim_index)

                if (signature == AGGREGATEDSUBSTITUTIONS) or (signature == AGGREGATEDINDELS) or (signature == AGGREGATEDDINUCS) or (signature == MICROHOMOLOGY) or (signature == REPEAT):
                    os.makedirs(os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, signature), exist_ok=True)
                    normalizedMutationDensityFilePath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, signature,normalizedMutationDensityFilename)
                    numberofMutationsFilePath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, signature,numberofMutationsFilename)
                else:
                    os.makedirs(os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, SIGNATUREBASED), exist_ok=True)
                    normalizedMutationDensityFilePath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, SIGNATUREBASED,normalizedMutationDensityFilename)
                    numberofMutationsFilePath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, SIGNATUREBASED,numberofMutationsFilename)

                numberofMutations, mutationDensityDict, numberofMutationsList = getMutationDensityDictUsingNumpyArray_enhanced(num_of_attributable_bases_list,
                                                                                                                               decile_counts_np_array)
                normalizedMutationDensityList = getNormalizedMutationDensityList(mutationDensityDict)

                # Normalized Mutation Density
                with open(normalizedMutationDensityFilePath, 'w') as file:
                    for normalizedMutationDensity in normalizedMutationDensityList:
                        file.write(str(normalizedMutationDensity) + ' ')
                    file.write('\n')

                # Number of Mutations
                with open(numberofMutationsFilePath, 'w') as file:
                    for numberofMutations in numberofMutationsList:
                        file.write(str(numberofMutations) + ' ')
                    file.write('\n')



# August 3, 2020
# Using numpy array
# decile_df_list is RepliSeq input file dependent
def writeReplicationTimeDataUsingNumpyArray(outputDir,
                                            jobname,
                                            decile_df_list,
                                            subs_signatures,
                                            dinucs_signatures,
                                            indels_signatures,
                                            all_sims_subs_signature_decile_index_accumulated_np_array,
                                            all_sims_dinucs_signature_decile_index_accumulated_np_array,
                                            all_sims_indels_signature_decile_index_accumulated_np_array):

    os.makedirs(os.path.join(outputDir, jobname, DATA, REPLICATIONTIME), exist_ok=True)

    my_list=[(SUBS, subs_signatures, all_sims_subs_signature_decile_index_accumulated_np_array),
             (DINUCS, dinucs_signatures, all_sims_dinucs_signature_decile_index_accumulated_np_array),
             (INDELS, indels_signatures, all_sims_indels_signature_decile_index_accumulated_np_array)]

    # Write Number of Attributable Bases List
    if (decile_df_list is not None):
        numberofAttributableBasesList = getNumberofAttributableBases(decile_df_list)
        numberofAttributabelBasesFilename = 'NumberofAttributableBases.txt'
        numberofAttributabelBasesFilePath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, numberofAttributabelBasesFilename)

        # Number of Attributable Bases
        with open(numberofAttributabelBasesFilePath, 'w') as file:
            for numberofAttributabelBases in numberofAttributableBasesList:
                file.write(str(numberofAttributabelBases) + ' ')
            file.write('\n')

    for my_tuple in my_list:
        my_type, signatures, all_sims_signature_decile_index_accumulated_np_array = my_tuple
        num_of_sims, num_of_signatures_with_extras, num_of_deciles = all_sims_signature_decile_index_accumulated_np_array.shape

        for sim_index in range(0,num_of_sims):
            for signature_index in range(0,num_of_signatures_with_extras):
                decile_counts_np_array = all_sims_signature_decile_index_accumulated_np_array[sim_index,signature_index]

                # Order is important
                # -1 contains AGGREGATEDSUBSTITUTIONS, AGGREGATEDDINUCS, AGGREGATEDINDELS for my_type==SUBS, my_type==DINUCS, my_type==INDELS respectively
                # -2 contains REPEAT for  my_type==INDELS
                # -3 contains MICROHOMOLOGY for  my_type==INDELS
                if signature_index < signatures.size:
                    signature = signatures[signature_index]
                elif signature_index == (num_of_signatures_with_extras-1) and my_type == SUBS:
                    signature = AGGREGATEDSUBSTITUTIONS
                elif signature_index == (num_of_signatures_with_extras-1) and my_type == DINUCS:
                    signature = AGGREGATEDDINUCS
                elif signature_index == (num_of_signatures_with_extras-1) and my_type == INDELS:
                    signature = AGGREGATEDINDELS
                elif signature_index == (num_of_signatures_with_extras-2) and my_type == INDELS:
                    signature = REPEAT
                elif signature_index == (num_of_signatures_with_extras-3) and my_type == INDELS:
                    signature = MICROHOMOLOGY

                if (sim_index == 0):
                    normalizedMutationDensityFilename = '%s_NormalizedMutationDensity.txt' %(signature)
                    numberofMutationsFilename = '%s_NumberofMutations.txt' %(signature)
                else:
                    normalizedMutationDensityFilename = '%s_sim%d_NormalizedMutationDensity.txt' %(signature,sim_index)
                    numberofMutationsFilename = '%s_sim%d_NumberofMutations.txt' %(signature,sim_index)

                if (signature == AGGREGATEDSUBSTITUTIONS) or (signature == AGGREGATEDINDELS) or (signature == AGGREGATEDDINUCS) or (signature == MICROHOMOLOGY) or (signature == REPEAT):
                    os.makedirs(os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, signature), exist_ok=True)
                    normalizedMutationDensityFilePath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, signature,normalizedMutationDensityFilename)
                    numberofMutationsFilePath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, signature,numberofMutationsFilename)
                else:
                    os.makedirs(os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, SIGNATUREBASED), exist_ok=True)
                    normalizedMutationDensityFilePath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, SIGNATUREBASED,normalizedMutationDensityFilename)
                    numberofMutationsFilePath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, SIGNATUREBASED,numberofMutationsFilename)

                if (decile_df_list is not None):
                    numberofMutations, mutationDensityDict, numberofMutationsList = getMutationDensityDictUsingNumpyArray(decile_df_list, decile_counts_np_array)
                    normalizedMutationDensityList = getNormalizedMutationDensityList(mutationDensityDict)

                    # Normalized Mutation Density
                    with open(normalizedMutationDensityFilePath, 'w') as file:
                        for normalizedMutationDensity in normalizedMutationDensityList:
                            file.write(str(normalizedMutationDensity) + ' ')
                        file.write('\n')

                    # Number of Mutations
                    with open(numberofMutationsFilePath, 'w') as file:
                        for numberofMutations in numberofMutationsList:
                            file.write(str(numberofMutations) + ' ')
                        file.write('\n')


# using pyranges
def find_num_of_mutations_in_bins_enhanced(computationType,
            outputDir,
            jobname,
            numofSimulations,
            samples_of_interest,
            job_tuples,
            chromSizesDict,
            chromNamesList,
            sorted_repli_seq_arr,
            ordered_sbs_signatures,
            ordered_dbs_signatures,
            ordered_id_signatures,
            ordered_sbs_signatures_cutoffs,
            ordered_dbs_signatures_cutoffs,
            ordered_id_signatures_cutoffs,
            discreet_mode,
            default_cutoff,
            parallel_mode):

    number_of_sbs_signatures = ordered_sbs_signatures.size
    number_of_dbs_signatures = ordered_dbs_signatures.size
    number_of_id_signatures = ordered_id_signatures.size

    # +1 for Aggregated
    all_sims_subs_signature_decile_index_accumulated_np_array = np.zeros(
        (numofSimulations + 1, number_of_sbs_signatures + 1, 10), dtype=float)  # int
    all_sims_dinucs_signature_decile_index_accumulated_np_array = np.zeros(
        (numofSimulations + 1, number_of_dbs_signatures + 1, 10), dtype=float)  # int
    all_sims_indels_signature_decile_index_accumulated_np_array = np.zeros(
        (numofSimulations + 1, number_of_id_signatures + 1, 10), dtype=float)  # int

    def accumulate_np_arrays(result_tuple):
        sim_num = result_tuple[0]
        subs_signature_decile_index_accumulated_np_array = result_tuple[1]
        dinucs_signature_decile_index_accumulated_np_array = result_tuple[2]
        indels_signature_decile_index_accumulated_np_array = result_tuple[3]

        # print('MONITOR ACCUMULATE', flush=True)

        # accumulate
        all_sims_subs_signature_decile_index_accumulated_np_array[sim_num] += subs_signature_decile_index_accumulated_np_array
        all_sims_dinucs_signature_decile_index_accumulated_np_array[sim_num] += dinucs_signature_decile_index_accumulated_np_array
        all_sims_indels_signature_decile_index_accumulated_np_array[sim_num] += indels_signature_decile_index_accumulated_np_array

    sim_nums = range(0, numofSimulations + 1)
    sim_num_chr_tuples = ((sim_num, chrLong) for sim_num in sim_nums for chrLong in chromNamesList)

    if parallel_mode:
        jobs = []

        # Using numpy array
        if (computationType == USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM):
            numofProcesses = multiprocessing.cpu_count()
            pool = multiprocessing.Pool(processes=numofProcesses)

            for simNum, chrLong in sim_num_chr_tuples:
                chromSize = chromSizesDict[chrLong]

                jobs.append(pool.apply_async(
                    mediator_enhanced,
                    args=(outputDir,
                          jobname,
                          chrLong,
                          chromSize,
                          simNum,
                          samples_of_interest,
                          sorted_repli_seq_arr,
                          ordered_sbs_signatures,
                          ordered_dbs_signatures,
                          ordered_id_signatures,
                          ordered_sbs_signatures_cutoffs,
                          ordered_dbs_signatures_cutoffs,
                          ordered_id_signatures_cutoffs,
                          discreet_mode,
                          default_cutoff,),
                    callback=accumulate_np_arrays))

            pool.close()
            pool.join()

    else:

        # Sequential mode for testing, debugging and profiling purposes
        for simNum, chrLong in sim_num_chr_tuples:
            chromSize = chromSizesDict[chrLong]

            result_tuple = mediator_enhanced(
                outputDir,
                jobname,
                chrLong,
                chromSize,
                simNum,
                samples_of_interest,
                sorted_repli_seq_arr,
                ordered_sbs_signatures,
                ordered_dbs_signatures,
                ordered_id_signatures,
                ordered_sbs_signatures_cutoffs,
                ordered_dbs_signatures_cutoffs,
                ordered_id_signatures_cutoffs,
                discreet_mode,
                default_cutoff)

            accumulate_np_arrays(result_tuple)

    return all_sims_subs_signature_decile_index_accumulated_np_array, \
           all_sims_dinucs_signature_decile_index_accumulated_np_array, \
           all_sims_indels_signature_decile_index_accumulated_np_array,



# main function using pyranges
# not preferred since it is not better/faster than the old method
def replicationTimeAnalysis_enhanced(computationType,
                            sample_based,
                            genome,
                            chromSizesDict,
                            chromNamesList,
                            outputDir,
                            jobname,
                            numofSimulations,
                            samples_of_interest,
                            all_samples_list,
                            job_tuples,
                            repliseqDataFilename,
                            ordered_sbs_signatures_with_cutoffs,
                            ordered_dbs_signatures_with_cutoffs,
                            ordered_id_signatures_with_cutoffs,
                            ordered_sbs_signatures_cutoffs,
                            ordered_dbs_signatures_cutoffs,
                            ordered_id_signatures_cutoffs,
                            discreet_mode,
                            default_cutoff,
                            parallel_mode,
                            log_file,
                            verbose,
                            matrix_generator_path):

    sbs_signatures = ordered_sbs_signatures_with_cutoffs
    dbs_signatures = ordered_dbs_signatures_with_cutoffs
    id_signatures = ordered_id_signatures_with_cutoffs

    log_out = open(log_file, 'a')
    print('\n#################################################################################', file=log_out)
    print('--- Replication Timing Analysis starts', file=log_out)
    print('--- Replication Timing Analyis Computation Type:%s' % (computationType), file=log_out)
    log_out.close()

    #########################################################################
    # Analysis Type can be
    # AggregatedSubstitutions: All in one
    # AggregatedIndels : All in one
    # AggregatedDinucs : All in one
    # IndelsBased : Microhomology, Repeat
    # SignatureBased: Subs Signatures Sig1, Sig2, ... and
    # Indels Signatures  ID1, ID2, ..., ... and
    # Dinucs Signatures  DBS1, DBS2, ...
    #########################################################################

    # bed files with 4 columns only
    # repli_seq_df = pd.read_csv(repliseqDataFilename, sep='\t', names=['chr', 'start', 'end', 'signal'])

    # Chrom      Start        End     Signal
    repli_seq_df = readWig_with_fixedStep_variableStep(repliseqDataFilename)

    # add bin column
    # valid values are 0, 1, ..., 9
    repli_seq_df['bin'] = -1 # -1 dummy value

    # create numpy array from repli-seq file
    repli_seq_arr = repli_seq_df.iloc[:,:].values

    # sort the array w.r.t. signal in descending order
    # 0: chr, 1:start, 2:end, 3:signal 4:bin
    sorted_repli_seq_arr = repli_seq_arr[np.argsort(repli_seq_arr[:, 3])[::-1]]

    # divide into 10 bins
    num_of_parts = 10
    indices = np.linspace(0, len(sorted_repli_seq_arr), num_of_parts+1, dtype=int)

    # update bin indices
    # 0: chr, 1:start, 2:end, 3:signal, 4:bin
    for i in range(num_of_parts):
        start_index = indices[i]
        end_index = indices[i+1]
        sorted_repli_seq_arr[start_index:end_index, 4] = i

    all_chroms_bin_num_of_attributable_bases_array_2d_arr = get_number_of_attributable_bases_in_each_bin_enhanced(genome,
                                                                                                                  sorted_repli_seq_arr,
                                                                                                                  matrix_generator_path)

    sim_sbs_signature_sample_mutation_type_count_arr, \
    sim_dbs_signature_sample_mutation_type_count_arr, \
    sim_id_signature_sample_mutation_type_count_array = find_num_of_mutations_in_bins_enhanced(computationType,
            outputDir,
            jobname,
            numofSimulations,
            samples_of_interest,
            job_tuples,
            chromSizesDict,
            chromNamesList,
            sorted_repli_seq_arr,
            ordered_sbs_signatures_with_cutoffs,
            ordered_dbs_signatures_with_cutoffs,
            ordered_id_signatures_with_cutoffs,
            ordered_sbs_signatures_cutoffs,
            ordered_dbs_signatures_cutoffs,
            ordered_id_signatures_cutoffs,
            discreet_mode,
            default_cutoff,
            parallel_mode)

    writeReplicationTimeDataUsingNumpyArray_enhanced(outputDir,
                                                jobname,
                                                all_chroms_bin_num_of_attributable_bases_array_2d_arr,
                                                sbs_signatures,
                                                dbs_signatures,
                                                id_signatures,
                                                sim_sbs_signature_sample_mutation_type_count_arr,
                                                sim_dbs_signature_sample_mutation_type_count_arr,
                                                sim_id_signature_sample_mutation_type_count_array)


    log_out = open(log_file, 'a')
    print('--- Replication Timing Analysis ends', file=log_out)
    print('#################################################################################\n', file=log_out)
    log_out.close()


# main function
def replication_time_analysis(computationType,
                            sample_based,
                            genome,
                            chromSizesDict,
                            chromNamesList,
                            outputDir,
                            jobname,
                            numofSimulations,
                            samples_of_interest,
                            all_samples_list,
                            job_tuples,
                            repli_seq_data_file_name,
                            ordered_sbs_signatures_with_cutoffs,
                            ordered_dbs_signatures_with_cutoffs,
                            ordered_id_signatures_with_cutoffs,
                            ordered_sbs_signatures_cutoffs,
                            ordered_dbs_signatures_cutoffs,
                            ordered_id_signatures_cutoffs,
                            discreet_mode,
                            default_cutoff,
                            parallel_mode,
                            log_file,
                            verbose,
                            matrix_generator_path):

    sbs_signatures = ordered_sbs_signatures_with_cutoffs
    dbs_signatures = ordered_dbs_signatures_with_cutoffs
    id_signatures = ordered_id_signatures_with_cutoffs

    log_out = open(log_file, 'a')
    print('\n#################################################################################', file=log_out)
    print('--- Replication Timing Analysis starts', file=log_out)
    print('--- Replication Timing Analyis Computation Type:%s' % (computationType), file=log_out)
    log_out.close()

    #########################################################################
    # Analysis Type can be
    # AggregatedSubstitutions: All in one
    # AggregatedIndels : All in one
    # AggregatedDinucs : All in one
    # IndelsBased : Microhomology, Repeat
    # SignatureBased: Subs Signatures Sig1, Sig2, ... and Indels Signatures  ID1, ID2, ..., ... and Dinucs Signatures  DBS1, DBS2, ...

    # We know the indels type we are interested in.
    # Microhomology indels --- len(indels) >= 3
    # Repeat indels --- len(indels) < 3
    #########################################################################

    # Fill replication np arrays during runtime

    ###################################################################
    ############### Read RepliSeq Time data starts ####################
    ###################################################################
    # Fist decile_df in decile_df_list contains the intervals that are replicated the earliest.
    # Last decile_df in decile_df_list contains the intervals that are replicated the latest.
    # Please note that each decile_df contains intervals from all chroms (mixed chroms)
    # What is the type of deciles? Deciles is a list of dataframes.
    # if verbose: print('\tVerbose Worker pid %s READ Repliseq DATA STARTS  %s MB' % (str(os.getpid()), memory_usage()))
    # if verbose: print('\tVerbose READ Repliseq DATA STARTS')

    # formerly downloaded matrix generator reference genome is being used.
    # decile_df_list is a list, each decile_df is a dataframe
    # first decile_df with index=1 contains the intervals that are replicated the earliest
    # tenth (last) decile_df with index=10 contains the intervals that are replicated the latest
    chrNamesInReplicationTimeDataArray, decile_df_list = read_repli_seq_time_data(genome,
                                                                              chromNamesList,
                                                                              repli_seq_data_file_name,
                                                                              matrix_generator_path,
                                                                              verbose,
                                                                              log_file)

    # if verbose: print('\tVerbose Worker pid %s READ Repliseq DATA ENDS  %s MB' % (str(os.getpid()), memory_usage()))
    # if verbose: print('\tVerbose READ Repliseq DATA ENDS')

    # Get chrBased grouped deciles
    chrBased_grouped_decile_df_list = []
    # Get this chrLong of each decile
    # The first decile is the earliest one with index 1
    # The last decile is the latest one with index 10
    # We provide these indexes later in the fillChrBasedDecileIndexArray function by enumerate
    for decile_df in decile_df_list:
        chrBased_grouped_decile_df = decile_df.groupby(CHROM)
        chrBased_grouped_decile_df_list.append(chrBased_grouped_decile_df)
    ###################################################################
    ############### Read RepliSeq Time data ends ######################
    ###################################################################

    #######################################################################################################
    ################################### Replication Time Data Analysis starts #############################
    #######################################################################################################

    # old method keep it for further guidance for sample based
    # writeReplicationTimeData(outputDir,jobname,sample_based,decile_df_list,None,simNum2Type2DecileBasedAllChrAccumulatedCountDict,simNum2Sample2Type2DecileBasedAllChrAccumulatedCountDict)

    if sample_based:
        sample_2_sample_index_dict = {sample: ind for ind, sample in enumerate(all_samples_list, 0)}
        sample_index_2_sample_dict = {ind: sample for ind, sample in enumerate(all_samples_list, 0)}

        # log_out = open(log_file, 'a')
        # print('\n#################################################################################', file=log_out)
        # print('--- sample_2_sample_index_dict:', sample_2_sample_index_dict, file=log_out)
        # print('--- sample_index_2_sample_dict:', sample_index_2_sample_dict, file=log_out)
        # log_out.close()

        # Ordered signatures will only have signatures since later on, they are used in filtering mutation row columns
        all_sims_samples_subs_signature_decile_index_accumulated_np_array, \
        all_sims_samples_dinucs_signature_decile_index_accumulated_np_array, \
        all_sims_samples_indels_signature_decile_index_accumulated_np_array = calculateCountsForMutationsFillingReplicationTimeNPArrayRuntime_sample_based_using_numpy_array(
            computationType,
            outputDir,
            jobname,
            numofSimulations,
            samples_of_interest,
            all_samples_list,
            sample_2_sample_index_dict,
            job_tuples,
            chromSizesDict,
            chromNamesList,
            chrBased_grouped_decile_df_list,
            ordered_sbs_signatures_with_cutoffs,
            ordered_dbs_signatures_with_cutoffs,
            ordered_id_signatures_with_cutoffs,
            ordered_sbs_signatures_cutoffs,
            ordered_dbs_signatures_cutoffs,
            ordered_id_signatures_cutoffs,
            discreet_mode,
            default_cutoff,
            parallel_mode)

        # log_out = open(log_file, 'a')

        # Sum across all samples
        all_sims_subs_signature_decile_index_accumulated_np_array = np.sum(all_sims_samples_subs_signature_decile_index_accumulated_np_array, axis=1)
        all_sims_dinucs_signature_decile_index_accumulated_np_array = np.sum(all_sims_samples_dinucs_signature_decile_index_accumulated_np_array, axis=1)
        all_sims_indels_signature_decile_index_accumulated_np_array = np.sum(all_sims_samples_indels_signature_decile_index_accumulated_np_array, axis=1)

    else:
        # Ordered signatures will only have signatures since later on, they are used in filtering mutation row columns
        all_sims_subs_signature_decile_index_accumulated_np_array, \
        all_sims_dinucs_signature_decile_index_accumulated_np_array, \
        all_sims_indels_signature_decile_index_accumulated_np_array, \
        list_of_all_decile1_dicts, \
        list_of_all_decile2_dicts, \
        list_of_all_decile3_dicts, \
        list_of_all_decile4_dicts, \
        list_of_all_decile5_dicts, \
        list_of_all_decile6_dicts, \
        list_of_all_decile7_dicts, \
        list_of_all_decile8_dicts, \
        list_of_all_decile9_dicts, \
        list_of_all_decile10_dicts = calculateCountsForMutationsFillingReplicationTimeNPArrayRuntime_using_numpy_array(
            computationType,
            outputDir,
            jobname,
            numofSimulations,
            samples_of_interest,
            job_tuples,
            chromSizesDict,
            chromNamesList,
            chrBased_grouped_decile_df_list,
            ordered_sbs_signatures_with_cutoffs,
            ordered_dbs_signatures_with_cutoffs,
            ordered_id_signatures_with_cutoffs,
            ordered_sbs_signatures_cutoffs,
            ordered_dbs_signatures_cutoffs,
            ordered_id_signatures_cutoffs,
            discreet_mode,
            default_cutoff,
            parallel_mode)

    writeReplicationTimeDataUsingNumpyArray(outputDir,
                                                jobname,
                                                decile_df_list,
                                                sbs_signatures,
                                                dbs_signatures,
                                                id_signatures,
                                                all_sims_subs_signature_decile_index_accumulated_np_array,
                                                all_sims_dinucs_signature_decile_index_accumulated_np_array,
                                                all_sims_indels_signature_decile_index_accumulated_np_array)

    # # for mutation annotation
    # write_mutations_decile_based_replication_timing(outputDir,
    #                                                 jobname,
    #                                                 list_of_all_decile1_dicts,
    #                                                 list_of_all_decile2_dicts,
    #                                                 list_of_all_decile3_dicts,
    #                                                 list_of_all_decile4_dicts,
    #                                                 list_of_all_decile5_dicts,
    #                                                 list_of_all_decile6_dicts,
    #                                                 list_of_all_decile7_dicts,
    #                                                 list_of_all_decile8_dicts,
    #                                                 list_of_all_decile9_dicts,
    #                                                 list_of_all_decile10_dicts)


    #######################################################################################################
    ################################### Replication Time Data Analysis ends ###############################
    #######################################################################################################

    log_out = open(log_file, 'a')
    print('--- Replication Timing Analysis ends', file=log_out)
    print('#################################################################################\n', file=log_out)
    log_out.close()