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
import pyranges as pr
import pandas as pd
import numpy as np
import math
import traceback
import pyBigWig

from typing import List, Tuple

from SigProfilerTopography.source.commons.TopographyCommons import memory_usage
from SigProfilerTopography.source.commons.TopographyCommons import readChrBasedMutationsDF
from SigProfilerTopography.source.commons.TopographyCommons import func_addSignal
from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofSubsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofIndelsDict
from SigProfilerTopography.source.commons.TopographyCommons import getDictionary

from SigProfilerTopography.source.commons.TopographyCommons import getSample2SubsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2IndelsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import writeSimulationBasedAverageOccupancyUsingNumpyArray

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
from SigProfilerTopography.source.occupancy.ChrBasedSignalArrays import readFileInBEDFormat

from SigProfilerTopography.source.commons.TopographyCommons import CHROM
from SigProfilerTopography.source.commons.TopographyCommons import START
from SigProfilerTopography.source.commons.TopographyCommons import END
from SigProfilerTopography.source.commons.TopographyCommons import SIGNAL

from SigProfilerTopography.source.occupancy.ChrBasedSignalArrays import readWig_write_derived_from_bedgraph_using_pool_chunks
from SigProfilerTopography.source.occupancy.ChrBasedSignalArrays import readWig_write_derived_from_bedgraph_using_pool_read_all

from SigProfilerTopography.source.commons.TopographyCommons import isFileTypeBedGraph
from SigProfilerTopography.source.commons.TopographyCommons import get_chrBased_simBased_combined_df_split
from SigProfilerTopography.source.commons.TopographyCommons import get_chrBased_simBased_combined_df
from SigProfilerTopography.source.commons.TopographyCommons import get_chrBased_simBased_dfs

from SigProfilerTopography.source.commons.TopographyCommons import MISSING_SIGNAL

from SigProfilerTopography.source.commons.TopographyCommons import SIGPROFILERTOPOGRAPHY_DEFAULT_FILES
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_ATAC_SEQ_OCCUPANCY_FILE
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_ATAC_SEQ_GRCh38_OCCUPANCY_FILE
from SigProfilerTopography.source.commons.TopographyCommons import ENCFF575PMI_mm10_embryonic_facial_prominence_ATAC_seq

pd.set_option('display.max_columns', None)
pd.set_option('display.width', 5000)

def chrbased_data_fill_signal_count_arrays_for_all_mutations_memory_efficient(occupancy_type,
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
                                                             plus_minus,
                                                             discreet_mode,
                                                             default_cutoff,
                                                             log_file,
                                                             verbose):


    if verbose:
        log_out = open(log_file, 'a')
        print('\tVerbose %s Worker pid %s memory_usage in %.2f MB START chrLong:%s simNum:%d' %(occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum), file=log_out)
        log_out.close()

    # 1st part  Prepare chr based mutations dataframes
    maximum_chrom_size = chromSizesDict[chrLong]
    start_time = time.time()
    array_length = 2*plus_minus+1

    # increase batch size for big data
    # decrease batch size for small data
    batch_size = 10000  # Adjust batch size as needed

    if verbose:
        log_out = open(log_file, 'a')
        print('\tVerbose %s Worker pid %s memory_usage in %.2f MB Check1 Read Signal Array and Dataframes chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum), file=log_out)
        log_out.close()

    number_of_sbs_signatures = ordered_sbs_signatures.size
    number_of_dbs_signatures = ordered_dbs_signatures.size
    number_of_id_signatures = ordered_id_signatures.size

    ###############################################################################
    ################################ Initialization ###############################
    ###############################################################################
    # Add one more row for the aggregated analysis
    subsSignature_accumulated_signal_np_array = np.zeros((number_of_sbs_signatures + 1, plus_minus * 2 + 1)) # dtype=float
    dinucsSignature_accumulated_signal_np_array = np.zeros((number_of_dbs_signatures + 1, plus_minus * 2 + 1)) # dtype=float
    indelsSignature_accumulated_signal_np_array = np.zeros((number_of_id_signatures + 1, plus_minus * 2 + 1)) # dtype=float

    # Add one more row for the aggregated analysis
    subsSignature_accumulated_count_np_array = np.zeros((number_of_sbs_signatures + 1, plus_minus * 2 + 1))
    dinucsSignature_accumulated_count_np_array = np.zeros((number_of_dbs_signatures + 1, plus_minus * 2 + 1))
    indelsSignature_accumulated_count_np_array = np.zeros((number_of_id_signatures + 1, plus_minus * 2 + 1))
    ###############################################################################
    ################################ Initialization ###############################
    ###############################################################################

    ######################################################## ######################
    ################### Fill signal and count array starts ########################
    ###############################################################################

    if verbose:
        log_out = open(log_file, 'a')
        print('\tVerbose %s Worker pid %s memory_usage in %.2f MB Check2_1 Start chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum), file=log_out)
        log_out.close()

    # For subs
    if ((chrBased_simBased_subs_df is not None) and (not chrBased_simBased_subs_df.empty)):

        sorted_chrBased_simBased_subs_df = chrBased_simBased_subs_df.sort_values(by='Start', ascending=True)
        is_sorted_ascending = sorted_chrBased_simBased_subs_df['Start'].is_monotonic_increasing
        print('DEBUG0 is_sorted_ascending for subs chrLong:', chrLong, 'simNum:', simNum, 'is_sorted_ascending:', is_sorted_ascending)

        # df_columns is a numpy array
        df_columns = chrBased_simBased_subs_df.columns.values

        # subsSignatures_mask_array is a boolean array
        # ordered_sbs_signatures is the names of the signatures we are interested in which do not necesaarily contain all of the signatures
        subsSignatures_mask_array = np.isin(df_columns, ordered_sbs_signatures)

        # df_columns: ['Sample' 'Chrom' 'Start' 'MutationLong' 'PyramidineStrand'
        #              'TranscriptionStrand' 'ReplicationTime' 'ReplicationStrand' 'Mutation'
        #              'SBS1' 'SBS2' 'SBS3' 'SBS4' 'SBS5' 'SBS8' 'SBS10a' 'SBS10b' 'SBS12'
        #              'SBS13' 'SBS14' 'SBS15' 'SBS16' 'SBS17a' 'SBS17b' 'SBS18' 'SBS20' 'SBS22'
        #              'SBS23' 'SBS25' 'SBS28' 'SBS30' 'SBS33' 'SBS34' 'SBS39' 'SBS40' 'SBS44'
        #              'SBS288P' 'Simulation_Number']

        mutations_df = sorted_chrBased_simBased_subs_df[[CHROM, START]]
        mutations_df[CHROM] = "chr" + mutations_df[CHROM]
        mutations_df[END] = mutations_df[START] + plus_minus
        mutations_df[START] = mutations_df[START] - plus_minus

        # probabilities rows: num_of_mutations, columns: num_of_signatures_of_interest, cells: probabilities
        probabilities = chrBased_simBased_subs_df[df_columns[subsSignatures_mask_array]]

        if discreet_mode:
            # Discreet way 1 or 0
            # Convert True into 1, and False into 0
            # threshold_mask_array rows: num_of_mutations, columns: num_of_signatures_of_interest, cells: True or False
            threshold_mask_array = np.greater_equal(probabilities, ordered_sbs_signatures_cutoffs)
            # mask_array rows: num_of_mutations, columns: num_of_signatures_of_interest, cells: 1 or 0
            # if a mutation is assigned to a signature, then corresponding cell is set to 1, otherwise 0
            # if a mutation is not assigned to any signature, then all cells in that row are set to 0
            # Case1: only one cell in a row can be set to 1, the rest of the row is all zeros.
            # Case2: all cells are zero
            mask_array = threshold_mask_array.astype(int)  # num_of_mutations * num_of_signatures_of_interest
        else:
            probabilities[probabilities < default_cutoff] = 0
            # Case1: only one cell in a row is greater than or equal to default cutoff, the rest of the row is all zeros.
            # Case2: only two cells in a row are greater than or equal to default cutoff (e.g.: 0.5 0.5), the rest of the row is all zeros.
            # Case3: all cells are zero
            mask_array = np.array(probabilities).astype(float)

        if mask_array.size == 0:
            if probabilities.shape[0] > 0:
                mask_array = np.zeros((probabilities.shape[0],1), dtype=float)
        else:
            mask_array = mask_array.to_numpy()

        # Get the indices of the maximum values for each row
        # axis=1 refers to the rows (finding the maximum value in each row).
        # case1: all cells in a row are zeros, max_index for that row will be 0
        # case2: The first cell is 1, the rest is zero, max_index for that row will be 0
        # to distinguish the case1 from case2, we need to keep mask_array and set max_index to -1 if all cells in a row are zeros.
        max_indices = np.argmax(mask_array, axis=1)

        # if there is no signature assigned then assign -1
        # All zeros in row, means that there is no signature assigned
        # Assigning 0 will be misleading meaning it is assigned to the first signature in ordered_signatures
        max_indices[np.all(mask_array == 0, axis=1)] = -1

        # add Signature column to mutations_df
        mutations_df['Signature'] = max_indices

        # rename column name from Chrom to
        mutations_df.rename(columns={'Chrom': 'Chromosome'}, inplace=True)

        num_of_mutations = len(mutations_df)
        current_start_index = 0  # This tracks where the current batch begins

        # Open the BigWig file (using context manager for automatic closing)
        try:
            with pyBigWig.open(library_file_with_path) as bw:

                # Process mutations in batches
                while current_start_index < num_of_mutations:
                    current_end_index = min(current_start_index + batch_size, num_of_mutations)
                    print('DEBUG1')

                    batch_mutations_df = mutations_df.iloc[current_start_index:current_end_index]
                    batch_mutations_pr = pr.PyRanges(batch_mutations_df)

                    # numpy arrays
                    current_batch_starts = batch_mutations_df[START].values
                    current_batch_ends = batch_mutations_df[END].values

                    max_chrom_size = bw.chroms().get(chrLong)
                    mutation_start_pos_min = np.clip(np.min(current_batch_starts), a_min=0, a_max=max_chrom_size)
                    mutation_end_pos_max = np.clip(np.max(current_batch_ends), a_min=0, a_max=max_chrom_size)

                    try:
                        # Fetch all overlapping intervals in one go (s, e, v)
                        # s is 0-based start, e is 1-based exclusive end
                        intervals: List[Tuple[int, int, float]] = bw.intervals(
                            chrLong,
                            mutation_start_pos_min,
                            mutation_end_pos_max
                        )

                    except Exception as e:
                        # Handle cases where the chromosome might not exist in the BigWig file
                        print(f"Warning: SBS Could not fetch data for {chrLong} simNum:{simNum}:{mutation_start_pos_min} and {mutation_end_pos_max}. Error: {e}")

                    # create library_intervals_pr
                    if intervals is not None:
                        data_array = np.array(intervals)

                        # 0-indexed column access: start, end, signal
                        starts = data_array[:, 0].astype(int)
                        ends = data_array[:, 1].astype(int)

                        if library_file_type == BIGWIG:
                            # signals = np.array([value for start, end, value in values_list])
                            signals = data_array[:, 2].astype(float)

                        chromosomes = np.full(len(starts), chrLong, dtype='<U10')

                        # avoiding pandas dataframe to reduce memory usage
                        # doesn't work because it gives no empty() attribute error.
                        # Create pandas dataframe.
                        chrBased_library_df = pd.DataFrame({'Chromosome': chromosomes,
                                                            'Start': starts,
                                                            'End': ends,
                                                            'Signal': signals})

                        chrBased_library_df['Signal'] = chrBased_library_df['Signal'].astype(float)
                        library_intervals_pr = pr.PyRanges(chrBased_library_df)

                        # combine information from two PyRanges objects based on overlapping intervals, by default inner
                        try:
                            joined = batch_mutations_pr.join(library_intervals_pr)

                        except Exception as e:
                            print(f"Warning: SBS Could not join data for {chrLong} simNum:{simNum}. Error: {e}", flush=True)

                        if len(joined) > 0:
                            print("DEBUG2", chrLong, "simNum:", simNum, "mutation_start_pos_min:", mutation_start_pos_min, "mutation_end_pos_max:", mutation_end_pos_max,  "len(joined):", len(joined))
                            overlap_starts = joined.df[['Start', 'Start_b']].max(axis=1).values - joined.df['Start'].values
                            overlap_ends = joined.df[['End', 'End_b']].min(axis=1).values - joined.df['Start'].values
                            overlap_signals = joined.df['Signal'].values
                            overlap_signatures = joined.df['Signature'].values

                            # 1st way
                            # Part1, for mutations that are not assigned to any signature
                            # update aggregated signal and count arrays
                            mask = overlap_signatures == -1
                            current_starts = overlap_starts[mask]
                            current_ends = overlap_ends[mask]
                            current_signals = overlap_signals[mask]

                            # vectorized np.add.at starts
                            # 1. Calculate total number of base-pairs to update
                            total_bp = np.sum(current_ends - current_starts)
                            print("DEBUG PART1 for vectorized_np_add_at", chrLong, "simNum:", simNum, "total_bp:", total_bp, "total_bp.shape:", total_bp.shape)

                            # 2. Pre-allocate arrays for vectorized update
                            # The 'indices' array holds the positional index (column) for the np.add.at call
                            global_col_indices = np.empty(total_bp, dtype=np.int32)
                            print("DEBUG PART1 for vectorized_np_add_at", chrLong, "simNum:", simNum, "global_col_indices.shape:", global_col_indices.shape)

                            # The 'signatures' array holds the signature index (row) for the np.add.at call
                            global_row_indices = np.empty(total_bp, dtype=np.int8)  # Assuming signature is small integer
                            print("DEBUG PART1 for vectorized_np_add_at", chrLong, "simNum:", simNum, "global_row_indices.shape:", global_row_indices.shape)

                            # The 'signals' array holds the signal value to add
                            global_signals = np.empty(total_bp, dtype=np.float32)
                            global_counts = np.ones(total_bp, dtype=np.int8)
                            print("DEBUG PART1 for vectorized_np_add_at", chrLong, "simNum:", simNum, "global_signals.shape:", global_signals.shape)

                            current_idx = 0
                            for start, end, signal in zip(current_starts,
                                                                     current_ends,
                                                                     current_signals):
                                size = end - start

                                # 3. Fill the arrays
                                # Positional indices (columns): 0, 1, 2, ..., size-1, 0, 1, 2, ...
                                global_col_indices[current_idx: current_idx + size] = np.arange(start, end)

                                # Signal values: signal, signal, ...
                                global_signals[current_idx: current_idx + size] = signal

                                current_idx += size

                            # 5. Perform the single, vectorized np.add.at call for aggregated (-1) update
                            # The aggregated row index is always -1 (or the last row, assuming -1 maps to the last index)
                            # NOTE: Ensure your array indexing correctly maps -1 to the aggregated row index (e.g., array.shape[0]-1)
                            aggregated_row_index = subsSignature_accumulated_signal_np_array.shape[0] - 1
                            global_aggregated_row_indices = np.full(total_bp, aggregated_row_index, dtype=np.int8)

                            np.add.at(
                                subsSignature_accumulated_signal_np_array,
                                (global_aggregated_row_indices, global_col_indices),
                                global_signals
                            )

                            # Repeat the same process for subsSignature_accumulated_count_np_array
                            np.add.at(
                                subsSignature_accumulated_count_np_array,
                                (global_aggregated_row_indices, global_col_indices),
                                global_counts
                            )

                            # end == plus_minus * 2 condition handling
                            # 1. Create a mask for the conditional update
                            end_point_mask_part2 = current_ends == (plus_minus * 2)
                            if np.any(end_point_mask_part2):
                                # Select the relevant data points
                                update_ends = current_ends[end_point_mask_part2]  # The column index
                                update_signals = current_signals[end_point_mask_part2]  # The signal value

                                # The row index for aggregation
                                aggregated_row_index = subsSignature_accumulated_signal_np_array.shape[0] - 1
                                # The count value is always 1
                                update_counts = np.ones_like(update_ends, dtype=np.int8) # or 1

                                # Update aggregated signal
                                np.add.at(
                                    subsSignature_accumulated_signal_np_array,
                                    (aggregated_row_index, update_ends),
                                    update_signals
                                )

                                # Update aggregated count
                                np.add.at(
                                    subsSignature_accumulated_count_np_array,
                                    (aggregated_row_index, update_ends),
                                    update_counts
                                )
                            # vectorized np.add.at ends

                            # # former np.add.at starts
                            # for start, end, signal in zip(current_starts, current_ends, current_signals):
                            #     # Indices to update for the current interval
                            #     indices = np.arange(start, end)
                            #
                            #     # Update aggregated index (-1)
                            #     np.add.at(subsSignature_accumulated_signal_np_array, (-1, indices), signal)
                            #     np.add.at(subsSignature_accumulated_count_np_array, (-1, indices), 1)
                            #
                            #     # subsSignature_accumulated_signal_np_array[-1, start:end] += signal  # for aggregated
                            #     # we may need to handle probability case, instead if 1, we may consider probability
                            #     # or we may leave as it is, since probabilities are considered in mask array
                            #     # subsSignature_accumulated_count_np_array[-1, start:end] += 1  # for aggregated
                            #
                            #     # starts: [  0  25  56 ... 199 282 422]  --> bigWig
                            #     # ends: [  25   56  118 ...  282  422 2000] --> bigWig
                            #     # signals: [0.08958    0.10479    0.12298    ... 0.36173999 0.19921    0.36173999] --> bigWig
                            #     if end == plus_minus * 2:
                            #         np.add.at(subsSignature_accumulated_signal_np_array, (-1, end), signal)
                            #         np.add.at(subsSignature_accumulated_count_np_array, (-1, end), 1)
                            #
                            #         # subsSignature_accumulated_signal_np_array[-1, end] += signal
                            #         # subsSignature_accumulated_count_np_array[-1, end] += 1
                            # # former np.add.at ends

                            # Part2, for mutations that are assigned to a signature
                            # update signature and aggregated signal and count arrays
                            mask = overlap_signatures != -1
                            current_starts = overlap_starts[mask]
                            current_ends = overlap_ends[mask]
                            current_signals = overlap_signals[mask]
                            current_signatures = overlap_signatures[mask]

                            # for vectorized np.add.at starts
                            # 1. Calculate total number of base-pairs to update
                            total_bp = np.sum(current_ends - current_starts)
                            print("DEBUG PART2 for vectorized_np_add_at", chrLong, "simNum:", simNum, "total_bp:", total_bp, "total_bp.shape:", total_bp.shape)

                            # 2. Pre-allocate arrays for vectorized update
                            # The 'indices' array holds the positional index (column) for the np.add.at call
                            global_col_indices = np.empty(total_bp, dtype=np.int32)
                            print("DEBUG PART2 for vectorized_np_add_at", chrLong, "simNum:", simNum, "global_col_indices.shape:", global_col_indices.shape)

                            # The 'signatures' array holds the signature index (row) for the np.add.at call
                            global_row_indices = np.empty(total_bp, dtype=np.int8)  # Assuming signature is small integer
                            print("DEBUG PART2 for vectorized_np_add_at", chrLong, "simNum:", simNum, "global_row_indices.shape:", global_row_indices.shape)

                            # The 'signals' array holds the signal value to add
                            global_signals = np.empty(total_bp, dtype=np.float32)
                            global_counts = np.ones(total_bp, dtype=np.int8)
                            print("DEBUG PART2 for vectorized_np_add_at", chrLong, "simNum:", simNum, "global_signals.shape:", global_signals.shape)

                            current_idx = 0
                            for start, end, signal, signature in zip(current_starts,
                                                                     current_ends,
                                                                     current_signals,
                                                                     current_signatures):
                                size = end - start

                                # 3. Fill the arrays
                                # Positional indices (columns): 0, 1, 2, ..., size-1, 0, 1, 2, ...
                                global_col_indices[current_idx: current_idx + size] = np.arange(start, end)

                                # Signature indices (rows): signature, signature, ...
                                global_row_indices[current_idx: current_idx + size] = signature

                                # Signal values: signal, signal, ...
                                global_signals[current_idx: current_idx + size] = signal

                                current_idx += size

                            # 4. Perform the single, vectorized np.add.at call for signature-specific update
                            np.add.at(
                                subsSignature_accumulated_signal_np_array,
                                (global_row_indices, global_col_indices),
                                global_signals
                            )

                            # 5. Perform the single, vectorized np.add.at call for aggregated (-1) update
                            # The aggregated row index is always -1 (or the last row, assuming -1 maps to the last index)
                            # NOTE: Ensure your array indexing correctly maps -1 to the aggregated row index (e.g., array.shape[0]-1)
                            aggregated_row_index = subsSignature_accumulated_signal_np_array.shape[0] - 1
                            global_aggregated_row_indices = np.full(total_bp, aggregated_row_index, dtype=np.int8)

                            np.add.at(
                                subsSignature_accumulated_signal_np_array,
                                (global_aggregated_row_indices, global_col_indices),
                                global_signals
                            )

                            # Repeat the same process for subsSignature_accumulated_count_np_array
                            np.add.at(
                                subsSignature_accumulated_count_np_array,
                                (global_row_indices, global_col_indices),
                                global_counts
                            )

                            np.add.at(
                                subsSignature_accumulated_count_np_array,
                                (global_aggregated_row_indices, global_col_indices),
                                global_counts
                            )

                            # end == plus_minus * 2 condition handling
                            # 1. Create a mask for the conditional update
                            end_point_mask_part2 = current_ends == (plus_minus * 2)

                            if np.any(end_point_mask_part2):
                                # Select the relevant data points
                                update_ends = current_ends[end_point_mask_part2]  # The column index
                                update_signals = current_signals[end_point_mask_part2]  # The signal value
                                update_signatures = current_signatures[end_point_mask_part2]  # The signature row index

                                # The row index for aggregation
                                aggregated_row_index = subsSignature_accumulated_signal_np_array.shape[0] - 1
                                # The count value is always 1
                                update_counts = np.ones_like(update_ends, dtype=np.int8) # 1

                                # Update specific signature signal
                                np.add.at(
                                    subsSignature_accumulated_signal_np_array,
                                    (update_signatures, update_ends),
                                    update_signals
                                )

                                # Update aggregated signal
                                np.add.at(
                                    subsSignature_accumulated_signal_np_array,
                                    (aggregated_row_index, update_ends),
                                    update_signals
                                )

                                # Update specific signature count
                                np.add.at(
                                    subsSignature_accumulated_count_np_array,
                                    (update_signatures, update_ends),
                                    update_counts
                                )

                                # Update aggregated count
                                np.add.at(
                                    subsSignature_accumulated_count_np_array,
                                    (aggregated_row_index, update_ends),
                                    update_counts
                                )
                            # for vectorized np.add.at ends

                            # # former np.add.at starts
                            # # fill arrays
                            # for start, end, signal, signature in zip(current_starts, current_ends, current_signals, current_signatures):
                            #
                            #     # Indices to update for the current interval
                            #     indices = np.arange(start, end)
                            #
                            #     # Update specific signature index
                            #     # subsSignature_accumulated_signal_np_array[signature, start:end] += signal
                            #     # np.add.at(array_to_update, (row_indices, col_indices), value_to_add)
                            #     np.add.at(subsSignature_accumulated_signal_np_array, (signature, indices), signal)
                            #     np.add.at(subsSignature_accumulated_count_np_array, (signature, indices), 1)
                            #
                            #     # Update aggregated index (-1)
                            #     np.add.at(subsSignature_accumulated_signal_np_array, (-1, indices), signal)
                            #     np.add.at(subsSignature_accumulated_count_np_array, (-1, indices), 1)
                            #
                            #     # subsSignature_accumulated_signal_np_array[signature, start:end] += signal
                            #     # subsSignature_accumulated_signal_np_array[-1, start:end] += signal  # for aggregated
                            #
                            #     # subsSignature_accumulated_count_np_array[signature, start:end] += 1
                            #     # subsSignature_accumulated_count_np_array[-1, start:end] += 1  # for aggregated
                            #
                            #     if end == plus_minus * 2:
                            #         # Update specific signature index
                            #         np.add.at(subsSignature_accumulated_signal_np_array, (signature, end), signal)
                            #         np.add.at(subsSignature_accumulated_count_np_array, (signature, end), 1)
                            #
                            #         # Update aggregated index (-1)
                            #         np.add.at(subsSignature_accumulated_signal_np_array, (-1, end), signal)
                            #         np.add.at(subsSignature_accumulated_count_np_array, (-1, end), 1)
                            #
                            #         # subsSignature_accumulated_signal_np_array[signature, end] += signal
                            #         # subsSignature_accumulated_signal_np_array[-1, end] += signal  # for aggregated
                            #
                            #         # subsSignature_accumulated_count_np_array[signature, end] += 1
                            #         # subsSignature_accumulated_count_np_array[-1, end] += 1  # for aggregated
                            # # former np.add.at ends

                    # Move the index to the next batch
                    current_start_index = current_end_index

        except Exception as e:
            print(f"A critical error occurred for SBS: {e}")

        print('DEBUG3 Finished SBS chrLong:', chrLong, 'simNum:', simNum, 'time taken in seconds:', time.time() - start_time, flush=True)


    # # For Dinucs
    # if ((chrBased_simBased_dinucs_df is not None) and (not chrBased_simBased_dinucs_df.empty)):
    #
    #     # df_columns is a numpy array
    #     df_columns = chrBased_simBased_dinucs_df.columns.values
    #
    #     dinucsSignatures_mask_array = np.isin(df_columns, ordered_dbs_signatures)
    #
    #     mutations_df = chrBased_simBased_dinucs_df[[CHROM, START]]
    #     mutations_df[CHROM] = "chr" + mutations_df[CHROM]
    #     mutations_df[END] = mutations_df[START] + plus_minus
    #     mutations_df[START] = mutations_df[START] - plus_minus
    #
    #     probabilities = chrBased_simBased_dinucs_df[df_columns[dinucsSignatures_mask_array]]
    #
    #     if discreet_mode:
    #         # Discreet way 1 or 0
    #         # Convert True into 1, and False into 0
    #         threshold_mask_array = np.greater_equal(probabilities, ordered_dbs_signatures_cutoffs)
    #         mask_array = threshold_mask_array.astype(int)  # num_of_mutations * num_of_signatures_of_interest
    #     else:
    #         probabilities[probabilities < default_cutoff] = 0
    #         mask_array = np.array(probabilities).astype(float)
    #
    #     if mask_array.size == 0:
    #         if probabilities.shape[0] > 0:
    #             mask_array = np.zeros((probabilities.shape[0],1), dtype=float)
    #     else:
    #         mask_array = mask_array.to_numpy()
    #
    #     # Get the indices of the maximum values for each row
    #     max_indices = np.argmax(mask_array, axis=1)
    #     max_indices[np.all(mask_array == 0, axis=1)] = -1
    #
    #     # add Signature column to mutations_df
    #     mutations_df['Signature'] = max_indices
    #
    #     # Open the BigWig file (using context manager for automatic closing)
    #     try:
    #         with pyBigWig.open(library_file_with_path) as bw:
    #
    #             for _, row in mutations_df.iterrows():
    #                 chrom = row[CHROM]
    #                 window_start_pos = row[START]
    #                 window_end_pos = row[END]  # Exclusive end for pyBigWig
    #                 signature = row['Signature']
    #
    #                 # Genomic coordinate corresponding to index 0 of the local array (B)
    #                 array_base_pos = window_start_pos
    #
    #                 try:
    #                     # Fetch all overlapping intervals in one go (s, e, v)
    #                     # s is 0-based start, e is 1-based exclusive end
    #                     intervals: List[Tuple[int, int, float]] = bw.intervals(
    #                         chrom,
    #                         window_start_pos,
    #                         window_end_pos
    #                     )
    #                 except Exception as e:
    #                     # Handle cases where the chromosome might not exist in the BigWig file
    #                     print(f"Warning: DBS Could not fetch data for {chrom}:{window_start_pos} and {window_end_pos}. Error: {e}")
    #                     continue
    #
    #                 if not intervals:
    #                     # print(f"No intervals found for {chrom}:{mutation_pos}")
    #                     continue
    #
    #                 # Iterate through the fetched intervals and apply the vectorized update
    #                 for interval_start, interval_end, value in intervals:
    #                     # interval_start (s): 0-based start of BigWig interval
    #                     # interval_end (e): 1-based exclusive end of BigWig interval
    #                     # array_base_pos (B): 1-based position of array index 0
    #
    #                     # 1. Calculate the start index within the total_arrays
    #                     # Array index = interval_position - array_base_pos
    #                     # The first affected index is max(0, interval_start - array_base_pos)
    #                     # We ensure it's not negative and not outside the array bounds [0, 2001)
    #                     array_start_index = max(window_start_pos, interval_start) - array_base_pos
    #                     # starts = joined.df[['Start', 'Start_b']].max(axis=1).values - joined.df['Start'].values
    #
    #                     # 2. Calculate the end index within the total_arrays (exclusive)
    #                     # The last affected index is interval_end - array_base_pos
    #                     # We use min(array_length, ...) to ensure it stays within bounds
    #                     array_end_index = min(window_end_pos, interval_end) - array_base_pos
    #                     # ends = joined.df[['End', 'End_b']].min(axis=1).values - joined.df['Start'].values
    #
    #                     # --- Vectorized Update ---
    #                     # Update the slices directly. This is the fast part.
    #                     if signature == -1 and array_start_index < array_end_index:
    #                         # Add the signal value to the signals array slice
    #                         dinucsSignature_accumulated_signal_np_array[-1, array_start_index:array_end_index] += value
    #
    #                         # Add 1 to the counts array slice
    #                         dinucsSignature_accumulated_count_np_array[-1, array_start_index:array_end_index] += 1
    #
    #                         if array_end_index == plus_minus * 2:
    #                             dinucsSignature_accumulated_signal_np_array[-1, array_end_index] += value
    #                             dinucsSignature_accumulated_count_np_array[-1, array_end_index] += 1
    #
    #                     elif array_start_index < array_end_index:
    #                         # Update the specific signature row
    #                         dinucsSignature_accumulated_signal_np_array[signature, array_start_index:array_end_index] += value
    #                         dinucsSignature_accumulated_count_np_array[signature, array_start_index:array_end_index] += 1
    #
    #                         # Also update the aggregated row
    #                         dinucsSignature_accumulated_signal_np_array[-1, array_start_index:array_end_index] += value
    #                         dinucsSignature_accumulated_count_np_array[-1, array_start_index:array_end_index] += 1
    #
    #                         if array_end_index == plus_minus * 2:
    #                             dinucsSignature_accumulated_signal_np_array[signature, array_end_index] += value
    #                             dinucsSignature_accumulated_signal_np_array[-1, array_end_index] += value  # for aggregated
    #
    #                             dinucsSignature_accumulated_count_np_array[signature, array_end_index] += 1
    #                             dinucsSignature_accumulated_count_np_array[-1, array_end_index] += 1  # for aggregated
    #
    #
    #     except Exception as e:
    #         print(f"A critical error occurred for DBS: {e}")

    # # For Indels
    # if ((chrBased_simBased_indels_df is not None) and (not chrBased_simBased_indels_df.empty)):
    #
    #     # df_columns is a numpy array
    #     df_columns = chrBased_simBased_indels_df.columns.values
    #
    #     indelsSignatures_mask_array = np.isin(df_columns, ordered_id_signatures)
    #
    #     mutations_df = chrBased_simBased_indels_df[[CHROM, START]]
    #     mutations_df[CHROM] = "chr" + mutations_df[CHROM]
    #     mutations_df[END] = mutations_df[START] + plus_minus
    #     mutations_df[START] = mutations_df[START] - plus_minus
    #
    #     probabilities = chrBased_simBased_indels_df[df_columns[indelsSignatures_mask_array]]
    #
    #     if discreet_mode:
    #         # Discreet way 1 or 0
    #         # Convert True into 1, and False into 0
    #         threshold_mask_array = np.greater_equal(probabilities, ordered_id_signatures_cutoffs)
    #         mask_array = threshold_mask_array.astype(int)  # num_of_mutations * num_of_signatures_of_interest
    #     else:
    #         probabilities[probabilities < default_cutoff] = 0
    #         mask_array = np.array(probabilities).astype(float)
    #
    #     if mask_array.size == 0:
    #         if probabilities.shape[0] > 0:
    #             mask_array = np.zeros((probabilities.shape[0],1), dtype=float)
    #     else:
    #         mask_array = mask_array.to_numpy()
    #
    #     # Get the indices of the maximum values for each row
    #     max_indices = np.argmax(mask_array, axis=1)
    #     max_indices[np.all(mask_array == 0, axis=1)] = -1
    #
    #     # add Signature column to mutations_df
    #     mutations_df['Signature'] = max_indices
    #
    #     # Open the BigWig file (using context manager for automatic closing)
    #     try:
    #         with pyBigWig.open(library_file_with_path) as bw:
    #
    #             for _, row in mutations_df.iterrows():
    #                 chrom = row[CHROM]
    #                 window_start_pos = row[START]
    #                 window_end_pos = row[END]  # Exclusive end for pyBigWig
    #                 signature = row['Signature']
    #
    #                 # Genomic coordinate corresponding to index 0 of the local array (B)
    #                 array_base_pos = window_start_pos
    #
    #                 try:
    #                     # Fetch all overlapping intervals in one go (s, e, v)
    #                     # s is 0-based start, e is 1-based exclusive end
    #                     intervals: List[Tuple[int, int, float]] = bw.intervals(
    #                         chrom,
    #                         window_start_pos,
    #                         window_end_pos
    #                     )
    #                 except Exception as e:
    #                     # Handle cases where the chromosome might not exist in the BigWig file
    #                     print(f"Warning: ID Could not fetch data for {chrom}:{window_start_pos} and {window_end_pos}. Error: {e}")
    #                     continue
    #
    #                 if not intervals:
    #                     # print(f"No intervals found for {chrom}:{mutation_pos}")
    #                     continue
    #
    #                 # Iterate through the fetched intervals and apply the vectorized update
    #                 for interval_start, interval_end, value in intervals:
    #                     # interval_start (s): 0-based start of BigWig interval
    #                     # interval_end (e): 1-based exclusive end of BigWig interval
    #                     # array_base_pos (B): 1-based position of array index 0
    #
    #                     # 1. Calculate the start index within the total_arrays
    #                     # Array index = interval_position - array_base_pos
    #                     # The first affected index is max(0, interval_start - array_base_pos)
    #                     # We ensure it's not negative and not outside the array bounds [0, 2001)
    #                     array_start_index = max(window_start_pos, interval_start) - array_base_pos
    #                     # starts = joined.df[['Start', 'Start_b']].max(axis=1).values - joined.df['Start'].values
    #
    #                     # 2. Calculate the end index within the total_arrays (exclusive)
    #                     # The last affected index is interval_end - array_base_pos
    #                     # We use min(array_length, ...) to ensure it stays within bounds
    #                     array_end_index = min(window_end_pos, interval_end) - array_base_pos
    #                     # ends = joined.df[['End', 'End_b']].min(axis=1).values - joined.df['Start'].values
    #
    #                     # --- Vectorized Update ---
    #                     # Update the slices directly. This is the fast part.
    #                     if signature == -1 and array_start_index < array_end_index:
    #                         # Add the signal value to the signals array slice
    #                         indelsSignature_accumulated_signal_np_array[-1, array_start_index:array_end_index] += value
    #
    #                         # Add 1 to the counts array slice
    #                         indelsSignature_accumulated_count_np_array[-1, array_start_index:array_end_index] += 1
    #
    #                         if array_end_index == plus_minus * 2:
    #                             indelsSignature_accumulated_signal_np_array[-1, array_end_index] += value
    #                             indelsSignature_accumulated_count_np_array[-1, array_end_index] += 1
    #
    #                     elif array_start_index < array_end_index:
    #                         # Update the specific signature row
    #                         indelsSignature_accumulated_signal_np_array[signature, array_start_index:array_end_index] += value
    #                         indelsSignature_accumulated_count_np_array[signature, array_start_index:array_end_index] += 1
    #
    #                         # Also update the aggregated row
    #                         indelsSignature_accumulated_signal_np_array[-1, array_start_index:array_end_index] += value
    #                         indelsSignature_accumulated_count_np_array[-1, array_start_index:array_end_index] += 1
    #
    #                         if array_end_index == plus_minus * 2:
    #                             indelsSignature_accumulated_signal_np_array[signature, array_end_index] += value
    #                             indelsSignature_accumulated_count_np_array[-1, array_end_index] += value  # for aggregated
    #
    #                             indelsSignature_accumulated_signal_np_array[signature, array_end_index] += 1
    #                             indelsSignature_accumulated_count_np_array[-1, array_end_index] += 1  # for aggregated
    #
    #
    #     except Exception as e:
    #         print(f"A critical error occurred for ID: {e}")

    if verbose:
        log_out = open(log_file, 'a')
        print('\tVerbose %s Worker pid %s memory_usage in %.2f MB Check2_2 End chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum), file=log_out)
        log_out.close()
    ###############################################################################
    ################### Fill signal and count array ends ##########################
    ###############################################################################

    if verbose:
        log_out = open(log_file, 'a')
        print('\tVerbose %s Worker pid %s memory_usage in %.2f MB END  chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum), file=log_out)
        print('\tVerbose %s Worker pid %s took %f seconds chrLong:%s simNum:%d\n' % (occupancy_type,str(os.getpid()), (time.time() - start_time), chrLong, simNum), file=log_out)
        log_out.close()


    # Initialzie the list, you will return this list
    SignalArrayAndCountArrayList = []

    SignalArrayAndCountArrayList.append(chrLong)
    SignalArrayAndCountArrayList.append(simNum)
    SignalArrayAndCountArrayList.append(subsSignature_accumulated_signal_np_array)
    SignalArrayAndCountArrayList.append(dinucsSignature_accumulated_signal_np_array)
    SignalArrayAndCountArrayList.append(indelsSignature_accumulated_signal_np_array)
    SignalArrayAndCountArrayList.append(subsSignature_accumulated_count_np_array)
    SignalArrayAndCountArrayList.append(dinucsSignature_accumulated_count_np_array)
    SignalArrayAndCountArrayList.append(indelsSignature_accumulated_count_np_array)

    return SignalArrayAndCountArrayList


def chrbased_data_fill_signal_count_arrays_for_all_mutations_using_pyranges(occupancy_type,
                                                             chrLong,
                                                             simNum,
                                                             chrBased_simBased_subs_df,
                                                             chrBased_simBased_dinucs_df,
                                                             chrBased_simBased_indels_df,
                                                             chromSizesDict,
                                                             chrBased_library_df,
                                                             ordered_sbs_signatures,
                                                             ordered_dbs_signatures,
                                                             ordered_id_signatures,
                                                             ordered_sbs_signatures_cutoffs,
                                                             ordered_dbs_signatures_cutoffs,
                                                             ordered_id_signatures_cutoffs,
                                                             plus_minus,
                                                             discreet_mode,
                                                             default_cutoff,
                                                             log_file,
                                                             verbose):


    if verbose:
        log_out = open(log_file, 'a')
        print('\tVerbose %s Worker pid %s memory_usage in %.2f MB START chrLong:%s simNum:%d' %(occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum), file=log_out)
        log_out.close()

    # 1st part  Prepare chr based mutations dataframes
    maximum_chrom_size = chromSizesDict[chrLong]
    start_time = time.time()

    # pyranges requires dataframe to have all the columns Chromosome, Start and End.
    # chrBased_library_df.columns.values: ['Chrom' 'Start' 'End' 'Signal']
    if chrBased_library_df is not None:
        chrBased_library_df.rename(columns={'Chrom':'Chromosome'}, inplace=True)

    # # CHROM, START, END, SIGNAL
    # library_intervals_pr = pr.PyRanges(chromosomes = chrBased_library_df[CHROM],
    #                                     starts = chrBased_library_df[START],
    #                                     ends=chrBased_library_df[END],
    #                                     signals = chrBased_library_df[SIGNAL])

    library_intervals_pr = pr.PyRanges(chrBased_library_df)


    if verbose:
        log_out = open(log_file, 'a')
        print('\tVerbose %s Worker pid %s memory_usage in %.2f MB Check1 Read Signal Array and Dataframes chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum), file=log_out)
        log_out.close()

    number_of_sbs_signatures = ordered_sbs_signatures.size
    number_of_dbs_signatures = ordered_dbs_signatures.size
    number_of_id_signatures = ordered_id_signatures.size

    ###############################################################################
    ################################ Initialization ###############################
    ###############################################################################
    # Add one more row for the aggregated analysis
    subsSignature_accumulated_signal_np_array = np.zeros((number_of_sbs_signatures + 1, plus_minus * 2 + 1)) # dtype=float
    dinucsSignature_accumulated_signal_np_array = np.zeros((number_of_dbs_signatures + 1, plus_minus * 2 + 1)) # dtype=float
    indelsSignature_accumulated_signal_np_array = np.zeros((number_of_id_signatures + 1, plus_minus * 2 + 1)) # dtype=float

    # Add one more row for the aggregated analysis
    subsSignature_accumulated_count_np_array = np.zeros((number_of_sbs_signatures + 1, plus_minus * 2 + 1))
    dinucsSignature_accumulated_count_np_array = np.zeros((number_of_dbs_signatures + 1, plus_minus * 2 + 1))
    indelsSignature_accumulated_count_np_array = np.zeros((number_of_id_signatures + 1, plus_minus * 2 + 1))
    ###############################################################################
    ################################ Initialization ###############################
    ###############################################################################

    if (library_intervals_pr is not None):
        ######################################################## ######################
        ################### Fill signal and count array starts ########################
        ###############################################################################

        if verbose:
            log_out = open(log_file, 'a')
            print('\tVerbose %s Worker pid %s memory_usage in %.2f MB Check2_1 Start chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum), file=log_out)
            log_out.close()

        # For subs
        if ((chrBased_simBased_subs_df is not None) and (not chrBased_simBased_subs_df.empty)):

            # df_columns is a numpy array
            df_columns = chrBased_simBased_subs_df.columns.values

            # subsSignatures_mask_array is a boolean array
            # ordered_sbs_signatures is the names of the signatures we are interested in which do not necesaarily contain all of the signatures
            subsSignatures_mask_array = np.isin(df_columns, ordered_sbs_signatures)

            # df_columns: ['Sample' 'Chrom' 'Start' 'MutationLong' 'PyramidineStrand'
            #              'TranscriptionStrand' 'ReplicationTime' 'ReplicationStrand' 'Mutation'
            #              'SBS1' 'SBS2' 'SBS3' 'SBS4' 'SBS5' 'SBS8' 'SBS10a' 'SBS10b' 'SBS12'
            #              'SBS13' 'SBS14' 'SBS15' 'SBS16' 'SBS17a' 'SBS17b' 'SBS18' 'SBS20' 'SBS22'
            #              'SBS23' 'SBS25' 'SBS28' 'SBS30' 'SBS33' 'SBS34' 'SBS39' 'SBS40' 'SBS44'
            #              'SBS288P' 'Simulation_Number']

            mutations_df = chrBased_simBased_subs_df[[CHROM, START]]
            mutations_df[CHROM] = "chr" + mutations_df[CHROM]
            mutations_df[END] = mutations_df[START] + plus_minus
            mutations_df[START] = mutations_df[START] - plus_minus

            # probabilities rows: num_of_mutations, columns: num_of_signatures_of_interest, cells: probabilities
            probabilities = chrBased_simBased_subs_df[df_columns[subsSignatures_mask_array]]

            if discreet_mode:
                # Discreet way 1 or 0
                # Convert True into 1, and False into 0
                # threshold_mask_array rows: num_of_mutations, columns: num_of_signatures_of_interest, cells: True or False
                threshold_mask_array = np.greater_equal(probabilities, ordered_sbs_signatures_cutoffs)
                # mask_array rows: num_of_mutations, columns: num_of_signatures_of_interest, cells: 1 or 0
                # if a mutation is assigned to a signature, then corresponding cell is set to 1, otherwise 0
                # if a mutation is not assigned to any signature, then all cells in that row are set to 0
                # Case1: only one cell in a row can be set to 1, the rest of the row is all zeros.
                # Case2: all cells are zero
                mask_array = threshold_mask_array.astype(int)  # num_of_mutations * num_of_signatures_of_interest
            else:
                probabilities[probabilities < default_cutoff] = 0
                # Case1: only one cell in a row is greater than or equal to default cutoff, the rest of the row is all zeros.
                # Case2: only two cells in a row are greater than or equal to default cutoff (e.g.: 0.5 0.5), the rest of the row is all zeros.
                # Case3: all cells are zero
                mask_array = np.array(probabilities).astype(float)

            if mask_array.size == 0:
                if probabilities.shape[0] > 0:
                    mask_array = np.zeros((probabilities.shape[0],1), dtype=float)
            else:
                mask_array = mask_array.to_numpy()

            # Get the indices of the maximum values for each row
            # axis=1 refers to the rows (finding the maximum value in each row).
            # case1: all cells in a row are zeros, max_index for that row will be 0
            # case2: The first cell is 1, the rest is zero, max_index for that row will be 0
            # to distinguish the case1 from case2, we need to keep mask_array and set max_index to -1 if all cells in a row are zeros.
            max_indices = np.argmax(mask_array, axis=1)

            # if there is no signature assigned then assign -1
            # All zeros in row, means that there is no signature assigned
            # Assigning 0 will be misleading meaning it is assigned to the first signature in ordered_signatures
            max_indices[np.all(mask_array == 0, axis=1)] = -1

            # add Signature column to mutations_df
            mutations_df['Signature'] = max_indices

            # rename column name from Chrom to
            mutations_df.rename(columns={'Chrom':'Chromosome'}, inplace=True)

            mutations_pr = pr.PyRanges(mutations_df)

            # combine information from two PyRanges objects based on overlapping intervals, by default inner
            joined = mutations_pr.join(library_intervals_pr)

            # Start and End holds the mutation positions as mutation_start - 1000 and mutation_start + 1000
            # Start_b and End_b holds the library positions
            # therefore, we substract Start
            # we convert the overlap coordinates into [0, 2000] array coordinates.

            if len(joined) > 0:
                starts = joined.df[['Start', 'Start_b']].max(axis=1).values - joined.df['Start'].values
                ends = joined.df[['End', 'End_b']].min(axis=1).values -  joined.df['Start'].values
                signals = joined.df['Signal'].values
                signatures = joined.df['Signature'].values

                # 1st way
                # Part1, for mutations that are not assigned to any signature
                # update aggregated signal and count arrays
                mask = signatures == -1
                current_starts = starts[mask]
                current_ends = ends[mask]
                current_signals = signals[mask]

                for start, end, signal in zip(current_starts, current_ends, current_signals):
                    subsSignature_accumulated_signal_np_array[-1, start:end ] += signal # for aggregated

                    # we may need to handle probability case, instead if 1, we may consider probability
                    # or we may leave as it is, since probabilities are considered in mask array
                    subsSignature_accumulated_count_np_array[-1, start:end ] += 1  # for aggregated

                    # starts: [  0  25  56 ... 199 282 422]  --> bigWig
                    # ends: [  25   56  118 ...  282  422 2000] --> bigWig
                    # signals: [0.08958    0.10479    0.12298    ... 0.36173999 0.19921    0.36173999] --> bigWig
                    if end == plus_minus * 2:
                        subsSignature_accumulated_signal_np_array[-1, end] += signal
                        subsSignature_accumulated_count_np_array[-1, end] += 1

                # Part2, for mutations that are assigned to a signature
                # update signature and aggregated signal and count arrays
                mask = signatures != -1
                current_starts = starts[mask]
                current_ends = ends[mask]
                current_signals = signals[mask]
                current_signatures = signatures[mask]

                # fill arrays
                # np.add.at can be tried and compared
                for start, end, signal, signature in zip(current_starts, current_ends, current_signals, current_signatures):
                    subsSignature_accumulated_signal_np_array[signature, start:end ] += signal
                    subsSignature_accumulated_signal_np_array[-1, start:end ] += signal # for aggregated

                    subsSignature_accumulated_count_np_array[signature, start:end ] += 1
                    subsSignature_accumulated_count_np_array[-1, start:end ] += 1  # for aggregated

                    if end == plus_minus * 2:
                        subsSignature_accumulated_signal_np_array[signature, end] += signal
                        subsSignature_accumulated_signal_np_array[-1, end] += signal  # for aggregated

                        subsSignature_accumulated_count_np_array[signature, end] += 1
                        subsSignature_accumulated_count_np_array[-1, end] += 1  # for aggregated

                # # 2nd way
                # # Ensure starts, ends, signals, and signatures are NumPy arrays)
                #
                # # Pre-calculate the maximum end position to avoid repeated calls to max()
                # max_end = plus_minus * 2
                #
                # # Optimized accumulation for each signature:
                # for signature in np.unique(signatures):
                #     mask = signatures == signature
                #     current_starts = starts[mask]
                #     current_ends = ends[mask]
                #     current_signals = signals[mask]
                #
                #     # Create a temporary array to store the updates for this signature
                #     updates_signal = np.zeros((max_end + 1,), dtype=subsSignature_accumulated_signal_np_array.dtype)
                #     updates_count = np.zeros((max_end + 1,), dtype=subsSignature_accumulated_count_np_array.dtype)
                #
                #     for start, end, signal in zip(current_starts, current_ends, current_signals):
                #         updates_signal[start:end + 1] += signal
                #         updates_count[start:end + 1] += 1
                #
                #     subsSignature_accumulated_signal_np_array[signature, :max_end + 1] += updates_signal
                #     subsSignature_accumulated_count_np_array[signature, :max_end + 1] += updates_count
                #
                #     # here might be problematic since we consider only mutations that are assigned to a signature
                #     # what about mutations that are not assigned to any signature?
                #     subsSignature_accumulated_signal_np_array[-1, :max_end + 1] += updates_signal  # Aggregate
                #     subsSignature_accumulated_count_np_array[-1, :max_end + 1] += updates_count  # Aggregate

        # For Dinucs
        if ((chrBased_simBased_dinucs_df is not None) and (not chrBased_simBased_dinucs_df.empty)):

            # df_columns is a numpy array
            df_columns = chrBased_simBased_dinucs_df.columns.values

            dinucsSignatures_mask_array = np.isin(df_columns, ordered_dbs_signatures)

            mutations_df = chrBased_simBased_dinucs_df[[CHROM, START]]
            mutations_df[CHROM] = "chr" + mutations_df[CHROM]
            mutations_df[END] = mutations_df[START] + plus_minus
            mutations_df[START] = mutations_df[START] - plus_minus

            probabilities = chrBased_simBased_dinucs_df[df_columns[dinucsSignatures_mask_array]]

            if discreet_mode:
                # Discreet way 1 or 0
                # Convert True into 1, and False into 0
                threshold_mask_array = np.greater_equal(probabilities, ordered_dbs_signatures_cutoffs)
                mask_array = threshold_mask_array.astype(int)  # num_of_mutations * num_of_signatures_of_interest
            else:
                probabilities[probabilities < default_cutoff] = 0
                mask_array = np.array(probabilities).astype(float)

            if mask_array.size == 0:
                if probabilities.shape[0] > 0:
                    mask_array = np.zeros((probabilities.shape[0],1), dtype=float)
            else:
                mask_array = mask_array.to_numpy()

            # Get the indices of the maximum values for each row
            max_indices = np.argmax(mask_array, axis=1)
            max_indices[np.all(mask_array == 0, axis=1)] = -1

            # add Signature column to mutations_df
            mutations_df['Signature'] = max_indices

            # rename column name from Chrom to
            mutations_df.rename(columns={'Chrom':'Chromosome'}, inplace=True)

            mutations_pr = pr.PyRanges(mutations_df)

            joined = mutations_pr.join(library_intervals_pr)

            if len(joined) > 0:
                starts = joined.df[['Start', 'Start_b']].max(axis=1).values - joined.df['Start'].values
                ends = joined.df[['End', 'End_b']].min(axis=1).values -  joined.df['Start'].values
                signals = joined.df['Signal'].values
                signatures = joined.df['Signature'].values

                # 1st way
                # Part1, for mutations that are not assigned to any signature
                # update aggregated signal and count arrays
                mask = signatures == -1
                current_starts = starts[mask]
                current_ends = ends[mask]
                current_signals = signals[mask]

                for start, end, signal in zip(current_starts, current_ends, current_signals):
                    dinucsSignature_accumulated_signal_np_array[-1, start:end ] += signal # for aggregated
                    dinucsSignature_accumulated_count_np_array[-1, start:end ] += 1  # for aggregated

                    if end == plus_minus * 2:
                        dinucsSignature_accumulated_signal_np_array[-1, end] += signal
                        dinucsSignature_accumulated_count_np_array[-1, end] += 1

                # Part2, for mutations that are assigned to a signature
                # update signature and aggregated signal and count arrays
                mask = signatures != -1
                current_starts = starts[mask]
                current_ends = ends[mask]
                current_signals = signals[mask]
                current_signatures = signatures[mask]

                for start, end, signal, signature in zip(current_starts, current_ends, current_signals, current_signatures):
                    dinucsSignature_accumulated_signal_np_array[signature, start:end ] += signal
                    dinucsSignature_accumulated_signal_np_array[-1, start:end ] += signal # for aggregated

                    dinucsSignature_accumulated_count_np_array[signature, start:end ] += 1
                    dinucsSignature_accumulated_count_np_array[-1, start:end ] += 1  # for aggregated

                    if end == plus_minus * 2:
                        dinucsSignature_accumulated_signal_np_array[signature, end] += signal
                        dinucsSignature_accumulated_signal_np_array[-1, end] += signal  # for aggregated

                        dinucsSignature_accumulated_count_np_array[signature, end] += 1
                        dinucsSignature_accumulated_count_np_array[-1, end] += 1  # for aggregated

                # # 2nd way
                # # Ensure starts, ends, signals, and signatures are NumPy arrays)
                #
                # # Pre-calculate the maximum end position to avoid repeated calls to max()
                # max_end = plus_minus * 2
                #
                # # Optimized accumulation for each signature:
                # for signature in np.unique(signatures):
                #     mask = signatures == signature
                #     current_starts = starts[mask]
                #     current_ends = ends[mask]
                #     current_signals = signals[mask]
                #
                #     # Create a temporary array to store the updates for this signature
                #     updates_signal = np.zeros((max_end + 1,), dtype=dinucsSignature_accumulated_signal_np_array.dtype)
                #     updates_count = np.zeros((max_end + 1,), dtype=dinucsSignature_accumulated_count_np_array.dtype)
                #
                #     for start, end, signal in zip(current_starts, current_ends, current_signals):
                #         updates_signal[start:end + 1] += signal
                #         updates_count[start:end + 1] += 1
                #
                #     dinucsSignature_accumulated_signal_np_array[signature, :max_end + 1] += updates_signal
                #     dinucsSignature_accumulated_count_np_array[signature, :max_end + 1] += updates_count
                #
                #     dinucsSignature_accumulated_signal_np_array[-1, :max_end + 1] += updates_signal  # Aggregate
                #     dinucsSignature_accumulated_count_np_array[-1, :max_end + 1] += updates_count  # Aggregate

        # For Indels
        if ((chrBased_simBased_indels_df is not None) and (not chrBased_simBased_indels_df.empty)):

            # df_columns is a numpy array
            df_columns = chrBased_simBased_indels_df.columns.values

            indelsSignatures_mask_array = np.isin(df_columns, ordered_id_signatures)

            mutations_df = chrBased_simBased_indels_df[[CHROM, START]]
            mutations_df[CHROM] = "chr" + mutations_df[CHROM]
            mutations_df[END] = mutations_df[START] + plus_minus
            mutations_df[START] = mutations_df[START] - plus_minus

            probabilities = chrBased_simBased_indels_df[df_columns[indelsSignatures_mask_array]]

            if discreet_mode:
                # Discreet way 1 or 0
                # Convert True into 1, and False into 0
                threshold_mask_array = np.greater_equal(probabilities, ordered_id_signatures_cutoffs)
                mask_array = threshold_mask_array.astype(int)  # num_of_mutations * num_of_signatures_of_interest
            else:
                probabilities[probabilities < default_cutoff] = 0
                mask_array = np.array(probabilities).astype(float)

            if mask_array.size == 0:
                if probabilities.shape[0] > 0:
                    mask_array = np.zeros((probabilities.shape[0],1), dtype=float)
            else:
                mask_array = mask_array.to_numpy()

            # Get the indices of the maximum values for each row
            max_indices = np.argmax(mask_array, axis=1)
            max_indices[np.all(mask_array == 0, axis=1)] = -1

            # add Signature column to mutations_df
            mutations_df['Signature'] = max_indices

            # rename column name from Chrom to
            mutations_df.rename(columns={'Chrom':'Chromosome'}, inplace=True)

            mutations_pr = pr.PyRanges(mutations_df)

            joined = mutations_pr.join(library_intervals_pr)

            if len(joined) > 0:
                starts = joined.df[['Start', 'Start_b']].max(axis=1).values - joined.df['Start'].values
                ends = joined.df[['End', 'End_b']].min(axis=1).values -  joined.df['Start'].values
                signals = joined.df['Signal'].values
                signatures = joined.df['Signature'].values

                # 1st way
                # Part1, for mutations that are not assigned to any signature
                # update aggregated signal and count arrays
                mask = signatures == -1
                current_starts = starts[mask]
                current_ends = ends[mask]
                current_signals = signals[mask]

                for start, end, signal in zip(current_starts, current_ends, current_signals):
                    indelsSignature_accumulated_signal_np_array[-1, start:end ] += signal # for aggregated
                    indelsSignature_accumulated_count_np_array[-1, start:end ] += 1  # for aggregated

                    if end == plus_minus * 2:
                        indelsSignature_accumulated_signal_np_array[-1, end] += signal
                        indelsSignature_accumulated_count_np_array[-1, end] += 1

                # Part2, for mutations that are assigned to a signature
                # update signature and aggregated signal and count arrays
                mask = signatures != -1
                current_starts = starts[mask]
                current_ends = ends[mask]
                current_signals = signals[mask]
                current_signatures = signatures[mask]

                for start, end, signal, signature in zip(current_starts, current_ends, current_signals, current_signatures):
                    indelsSignature_accumulated_signal_np_array[signature, start:end ] += signal
                    indelsSignature_accumulated_signal_np_array[-1, start:end ] += signal # for aggregated

                    indelsSignature_accumulated_count_np_array[signature, start:end ] += 1
                    indelsSignature_accumulated_count_np_array[-1, start:end ] += 1  # for aggregated

                    if end == plus_minus * 2:
                        indelsSignature_accumulated_signal_np_array[signature, end] += signal
                        indelsSignature_accumulated_signal_np_array[-1, end] += signal  # for aggregated

                        indelsSignature_accumulated_count_np_array[signature, end] += 1
                        indelsSignature_accumulated_count_np_array[-1, end] += 1  # for aggregated

                # # 2nd way
                # # Ensure starts, ends, signals, and signatures are NumPy arrays)
                #
                # # Pre-calculate the maximum end position to avoid repeated calls to max()
                # max_end = plus_minus * 2
                #
                # # Optimized accumulation for each signature:
                # for signature in np.unique(signatures):
                #     mask = signatures == signature
                #     current_starts = starts[mask]
                #     current_ends = ends[mask]
                #     current_signals = signals[mask]
                #
                #     # Create a temporary array to store the updates for this signature
                #     updates_signal = np.zeros((max_end + 1,), dtype=indelsSignature_accumulated_signal_np_array.dtype)
                #     updates_count = np.zeros((max_end + 1,), dtype=indelsSignature_accumulated_count_np_array.dtype)
                #
                #     for start, end, signal in zip(current_starts, current_ends, current_signals):
                #         updates_signal[start:end + 1] += signal
                #         updates_count[start:end + 1] += 1
                #
                #     indelsSignature_accumulated_signal_np_array[signature, :max_end + 1] += updates_signal
                #     indelsSignature_accumulated_count_np_array[signature, :max_end + 1] += updates_count
                #
                #     indelsSignature_accumulated_signal_np_array[-1, :max_end + 1] += updates_signal  # Aggregate
                #     indelsSignature_accumulated_count_np_array[-1, :max_end + 1] += updates_count  # Aggregate

        if verbose:
            log_out = open(log_file, 'a')
            print('\tVerbose %s Worker pid %s memory_usage in %.2f MB Check2_2 End chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum), file=log_out)
            log_out.close()
        ###############################################################################
        ################### Fill signal and count array ends ##########################
        ###############################################################################

    if verbose:
        log_out = open(log_file, 'a')
        print('\tVerbose %s Worker pid %s memory_usage in %.2f MB END  chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum), file=log_out)
        print('\tVerbose %s Worker pid %s took %f seconds chrLong:%s simNum:%d\n' % (occupancy_type,str(os.getpid()), (time.time() - start_time), chrLong, simNum), file=log_out)
        log_out.close()


    # Initialzie the list, you will return this list
    SignalArrayAndCountArrayList = []

    SignalArrayAndCountArrayList.append(chrLong)
    SignalArrayAndCountArrayList.append(simNum)
    SignalArrayAndCountArrayList.append(subsSignature_accumulated_signal_np_array)
    SignalArrayAndCountArrayList.append(dinucsSignature_accumulated_signal_np_array)
    SignalArrayAndCountArrayList.append(indelsSignature_accumulated_signal_np_array)
    SignalArrayAndCountArrayList.append(subsSignature_accumulated_count_np_array)
    SignalArrayAndCountArrayList.append(dinucsSignature_accumulated_count_np_array)
    SignalArrayAndCountArrayList.append(indelsSignature_accumulated_count_np_array)

    return SignalArrayAndCountArrayList



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
                                                             remove_outliers,
                                                             plusorMinus,
                                                             discreet_mode,
                                                             default_cutoff,
                                                             log_file,
                                                             verbose):


    if verbose:
        log_out = open(log_file, 'a')
        print('\tVerbose %s Worker pid %s memory_usage in %.2f MB START chrLong:%s simNum:%d' %(occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum), file=log_out)
        log_out.close()

    # 1st part  Prepare chr based mutations dataframes
    maximum_chrom_size = chromSizesDict[chrLong]
    start_time = time.time()

    chrBasedSignalArray = None # Will be filled from chrBasedSignal files if they exists
    library_file_opened_by_pyBigWig = None # Will be filled by pyBigWig from bigWig or bigBed
    my_upperBound = np.iinfo(np.int32).max
    signal_index = None

    if verbose:
        log_out = open(log_file, 'a')
        print('\tVerbose %s Worker pid %s memory_usage in %.2f MB Check1 Read Signal Array and Dataframes chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum), file=log_out)
        log_out.close()

    libraryFilenameWoExtension = os.path.splitext(os.path.basename(library_file_with_path))[0]
    signalArrayFilename = '%s_signal_%s.npy' % (chrLong, libraryFilenameWoExtension)
    if (occupancy_type == NUCLEOSOMEOCCUPANCY):
        chrBasedSignalFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME, CHRBASED, signalArrayFilename)
    elif (occupancy_type == EPIGENOMICSOCCUPANCY):
        if ((os.path.basename(library_file_with_path) == DEFAULT_ATAC_SEQ_OCCUPANCY_FILE) or \
                (os.path.basename(library_file_with_path) == DEFAULT_ATAC_SEQ_GRCh38_OCCUPANCY_FILE) or \
                (os.path.basename(library_file_with_path) == ENCFF575PMI_mm10_embryonic_facial_prominence_ATAC_seq)):
            chrBasedSignalFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, EPIGENOMICS, CHRBASED, signalArrayFilename)
        else:
            chrBasedSignalFile = os.path.join(outputDir,jobname, DATA, occupancy_type, LIB, CHRBASED, signalArrayFilename)

    else:
        # It can be EPIGENOMICSOCCUPANCY or user provided name e.g.: Epigenomics_ATAC_ENCFF317TWD
        chrBasedSignalFile = os.path.join(outputDir, jobname, DATA, occupancy_type, LIB, CHRBASED, signalArrayFilename)

    # Downloaded or created runtime
    if (os.path.exists(chrBasedSignalFile)):
        # Can this cause to deep sleep of processes?
        # chrBasedSignalArray = np.load(chrBasedSignalFile, mmap_mode='r')
        chrBasedSignalArray = np.load(chrBasedSignalFile)

    # If library_file_with_path is abs path and library_file_type is BIGWIG or BIGBED
    # For nucleosome_biosample==GM12878 or nucleosome_biosample==K562 library_file_with_path is only filename with extension, it is not absolute path
    if os.path.isabs(library_file_with_path):

        # Comment below to make it run in windows
        if (library_file_type == BIGWIG):
            try:
                library_file_opened_by_pyBigWig = pyBigWig.open(library_file_with_path)
                if (library_file_opened_by_pyBigWig is not None) and (chrLong in library_file_opened_by_pyBigWig.chroms()):
                    maximum_chrom_size = library_file_opened_by_pyBigWig.chroms()[chrLong]

                if (library_file_opened_by_pyBigWig is not None) and remove_outliers:

                    # For BigWig Files information in header is correct
                    # Using mean and standard deviation to define outlier
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

                    # # Using IQR better when normal distribution is not satisfied
                    # chrom_length = library_file_opened_by_pyBigWig.chroms(chrLong)
                    # values = library_file_opened_by_pyBigWig.values(chrLong, 0, chrom_length)
                    #
                    # # Step 2: Convert the list to a NumPy array for efficient calculation.
                    # # The numpy.nanpercentile() function will automatically ignore NaN values.
                    # data = np.array(values)
                    #
                    # # Step 3: Calculate Q1, Q3, and the Interquartile Range (IQR).
                    # q1 = np.nanpercentile(data, 25)
                    # q3 = np.nanpercentile(data, 75)
                    # iqr = q3 - q1
                    #
                    # # Step 4: Compute the upper outlier boundary.
                    # upper_outlier_bound = q3 + 1.5 * iqr
                    #
                    # # Step 5: Count how many items in the dataset are greater than the boundary.
                    # # Filter out nan values first to avoid incorrect comparisons.
                    # clean_data = data[~np.isnan(data)]
                    # num_outliers = np.sum(clean_data > upper_outlier_bound)
                    #
                    # print(f"Statistics for {chrLong} from {library_file_with_path}:")
                    # print(f"  Q1: {q1}")
                    # print(f"  Q3: {q3}")
                    # print(f"  IQR: {iqr}")
                    # print(f"  Upper Outlier Boundary (Q3 + 1.5 * IQR): {upper_outlier_bound}")
                    # print(f"  Number of items greater than the boundary: {num_outliers}")
                    # my_upperBound = upper_outlier_bound


            except:
                log_out = open(log_file, 'a')
                print('Exception %s' %library_file_with_path, file=log_out)
                log_out.close()

        elif (library_file_type == BIGBED):

            try:
                library_file_opened_by_pyBigWig = pyBigWig.open(library_file_with_path)
                if (library_file_opened_by_pyBigWig is not None) and (BED_6PLUS4 in str(library_file_opened_by_pyBigWig.SQL())):
                    signal_index = 3
                elif (library_file_opened_by_pyBigWig is not None) and (BED_9PLUS2 in str(library_file_opened_by_pyBigWig.SQL())):
                    signal_index = 7

                if (library_file_opened_by_pyBigWig is not None) and remove_outliers:
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
                log_out = open(log_file, 'a')
                print('Exception %s' %library_file_with_path, file=log_out)
                log_out.close()

    number_of_sbs_signatures = ordered_sbs_signatures.size
    number_of_dbs_signatures = ordered_dbs_signatures.size
    number_of_id_signatures = ordered_id_signatures.size

    ###############################################################################
    ################################ Initialization ###############################
    ###############################################################################
    # Add one more row for the aggregated analysis
    subsSignature_accumulated_signal_np_array = np.zeros((number_of_sbs_signatures + 1, plusorMinus * 2 + 1)) # dtype=float
    dinucsSignature_accumulated_signal_np_array = np.zeros((number_of_dbs_signatures + 1, plusorMinus * 2 + 1)) # dtype=float
    indelsSignature_accumulated_signal_np_array = np.zeros((number_of_id_signatures + 1, plusorMinus * 2 + 1)) # dtype=float

    # Add one more row for the aggregated analysis
    subsSignature_accumulated_count_np_array = np.zeros((number_of_sbs_signatures + 1, plusorMinus * 2 + 1))
    dinucsSignature_accumulated_count_np_array = np.zeros((number_of_dbs_signatures + 1, plusorMinus * 2 + 1))
    indelsSignature_accumulated_count_np_array = np.zeros((number_of_id_signatures + 1, plusorMinus * 2 + 1))
    ###############################################################################
    ################################ Initialization ###############################
    ###############################################################################

    if ((chrBasedSignalArray is not None) or ((library_file_opened_by_pyBigWig is not None) and (chrLong in library_file_opened_by_pyBigWig.chroms()))):
        ######################################################## #######################
        ################### Fill signal and count array starts ########################
        ###############################################################################
        if verbose:
            log_out = open(log_file, 'a')
            print('\tVerbose %s Worker pid %s memory_usage in %.2f MB Check2_1 Start chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum), file=log_out)
            log_out.close()

        # For subs
        if ((chrBased_simBased_subs_df is not None) and (not chrBased_simBased_subs_df.empty)):

            # df_columns is a numpy array
            df_columns = chrBased_simBased_subs_df.columns.values

            subsSignatures_mask_array = np.isin(df_columns, ordered_sbs_signatures)

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
                discreet_mode,
                default_cutoff,
                df_columns,
                occupancy_calculation_type) for row in chrBased_simBased_subs_df[df_columns].values]

        # For Dinucs
        if ((chrBased_simBased_dinucs_df is not None) and (not chrBased_simBased_dinucs_df.empty)):

            # df_columns is a numpy array
            df_columns = chrBased_simBased_dinucs_df.columns.values

            dinucsSignatures_mask_array = np.isin(df_columns, ordered_dbs_signatures)

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
                discreet_mode,
                default_cutoff,
                df_columns,
                occupancy_calculation_type) for row in chrBased_simBased_dinucs_df[df_columns].values]

        # For Indels
        if ((chrBased_simBased_indels_df is not None) and (not chrBased_simBased_indels_df.empty)):

            # df_columns is a numpy array
            df_columns = chrBased_simBased_indels_df.columns.values

            indelsSignatures_mask_array = np.isin(df_columns, ordered_id_signatures)

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
                discreet_mode,
                default_cutoff,
                df_columns,
                occupancy_calculation_type) for row in chrBased_simBased_indels_df[df_columns].values]

        if verbose:
            log_out = open(log_file, 'a')
            print('\tVerbose %s Worker pid %s memory_usage in %.2f MB Check2_2 End chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum), file=log_out)
            log_out.close()
        ###############################################################################
        ################### Fill signal and count array ends ##########################
        ###############################################################################

    if (library_file_opened_by_pyBigWig is not None):
        library_file_opened_by_pyBigWig.close()

    if verbose:
        log_out = open(log_file, 'a')
        print('\tVerbose %s Worker pid %s memory_usage in %.2f MB END  chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum), file=log_out)
        print('\tVerbose %s Worker pid %s took %f seconds chrLong:%s simNum:%d\n' % (occupancy_type,str(os.getpid()), (time.time() - start_time), chrLong, simNum), file=log_out)
        log_out.close()

    # Initialzie the list, you will return this list
    SignalArrayAndCountArrayList = []

    SignalArrayAndCountArrayList.append(chrLong)
    SignalArrayAndCountArrayList.append(simNum)
    SignalArrayAndCountArrayList.append(subsSignature_accumulated_signal_np_array)
    SignalArrayAndCountArrayList.append(dinucsSignature_accumulated_signal_np_array)
    SignalArrayAndCountArrayList.append(indelsSignature_accumulated_signal_np_array)
    SignalArrayAndCountArrayList.append(subsSignature_accumulated_count_np_array)
    SignalArrayAndCountArrayList.append(dinucsSignature_accumulated_count_np_array)
    SignalArrayAndCountArrayList.append(indelsSignature_accumulated_count_np_array)

    return SignalArrayAndCountArrayList

def read_chrom_based_bigwig_file(library_file_with_path, library_file_type, chrLong):

    chrBased_library_df = None

    if os.path.exists(library_file_with_path):
        bw = pyBigWig.open(library_file_with_path)
        try:
            if (bw is not None) and (chrLong in bw.chroms()):
                chrom_length = bw.chroms(chrLong)

                if library_file_type == BIGBED:
                    values_list = bw.entries(chrLong, 0, chrom_length) # for bigBed files

                elif library_file_type == BIGWIG:
                    # Get all intervals at once. Intervals are returned as a list of tuples (start, end, value).
                    values_list = bw.intervals(chrLong, 0, chrom_length) # for bigWig files

                if values_list is not None:
                    data_array = np.array(values_list)

                    # Convert the list of tuples to numpy arrays.
                    # starts = np.array([start for start, end, value in values_list])
                    # ends = np.array([end for start, end, value in values_list])

                    # 0-indexed column access: start, end, signal
                    starts = data_array[:, 0].astype(int)
                    ends = data_array[:, 1].astype(int)

                    if library_file_type == BIGWIG:
                        # signals = np.array([value for start, end, value in values_list])
                        signals = data_array[:, 2].astype(float)
                    elif library_file_type == BIGBED:
                        signals = np.array([value.split("\t")[3] for start, end, value in values_list])
                        # starts: [0] ends: [1025] signals: [59.6297] --> bigBed
                        ends = ends + 1

                    chromosomes = np.full(len(starts), chrLong, dtype='<U10')

                    # Create pandas dataframe.
                    chrBased_library_df = pd.DataFrame({'Chromosome': chromosomes,
                                                        'Start': starts,
                                                        'End': ends,
                                                        'Signal': signals})

        except KeyError:
            print(f"Chromosome '{chrLong}' not found in the BigWig file.")
        except RuntimeError as e:
            print(f"Error processing {chrLong}: {e}")
        finally:
            bw.close()

    return chrBased_library_df

def chrbased_data_fill_signal_count_arrays_for_all_mutations_read_mutations_memory_efficient(occupancy_type,
        occupancy_calculation_type,
        outputDir,
        jobname,
        chrLong,
        simNum,
        samples_of_interest,
        chromSizesDict,
        library_grouped_df,
        chrBased_library_df,
        library_file_with_path,
        library_file_type,
        ordered_sbs_signatures,
        ordered_dbs_signatures,
        ordered_id_signatures,
        ordered_sbs_signatures_cutoffs,
        ordered_dbs_signatures_cutoffs,
        ordered_id_signatures_cutoffs,
        plusorMinus,
        discreet_mode,
        default_cutoff,
        log_file,
        verbose):

    try:
        chrBased_simBased_subs_df, chrBased_simBased_dinucs_df, chrBased_simBased_indels_df = get_chrBased_simBased_dfs(outputDir, jobname, chrLong, simNum)

        # filter chrbased_df for samples_of_interest
        if samples_of_interest is not None:
            if chrBased_simBased_subs_df is not None:
                chrBased_simBased_subs_df = chrBased_simBased_subs_df[chrBased_simBased_subs_df['Sample'].isin(samples_of_interest)]

            if chrBased_simBased_dinucs_df is not None:
                chrBased_simBased_dinucs_df = chrBased_simBased_dinucs_df[chrBased_simBased_dinucs_df['Sample'].isin(samples_of_interest)]

            if chrBased_simBased_indels_df is not None:
                chrBased_simBased_indels_df = chrBased_simBased_indels_df[chrBased_simBased_indels_df['Sample'].isin(samples_of_interest)]

        return chrbased_data_fill_signal_count_arrays_for_all_mutations_memory_efficient(occupancy_type,
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
                                                                        discreet_mode,
                                                                        default_cutoff,
                                                                        log_file,
                                                                        verbose)
    except Exception as e:
        log_out = open(log_file, 'a')
        print("There is exception: %s --- chromosome: %s simulation_number: %s" % (e, chrLong, simNum), file=log_out)
        log_out.close()


def chrbased_data_fill_signal_count_arrays_for_all_mutations_read_mutations_using_pyranges(occupancy_type,
        occupancy_calculation_type,
        outputDir,
        jobname,
        chrLong,
        simNum,
        samples_of_interest,
        chromSizesDict,
        library_grouped_df,
        chrBased_library_df,
        library_file_with_path,
        library_file_type,
        ordered_sbs_signatures,
        ordered_dbs_signatures,
        ordered_id_signatures,
        ordered_sbs_signatures_cutoffs,
        ordered_dbs_signatures_cutoffs,
        ordered_id_signatures_cutoffs,
        plusorMinus,
        discreet_mode,
        default_cutoff,
        log_file,
        verbose):

    try:
        chrBased_simBased_subs_df, chrBased_simBased_dinucs_df, chrBased_simBased_indels_df = get_chrBased_simBased_dfs(outputDir, jobname, chrLong, simNum)

        if library_grouped_df is not None:
            chrBased_library_df = library_grouped_df.get_group(chrLong)

        if (library_file_type == BIGWIG) or (library_file_type == BIGBED):
            # read the bigWig file
            chrBased_library_df = read_chrom_based_bigwig_file(library_file_with_path, library_file_type, chrLong)
            if chrBased_library_df is not None:
                chrBased_library_df['Signal'] = chrBased_library_df['Signal'].astype(float)

        # filter chrbased_df for samples_of_interest
        if samples_of_interest is not None:
            if chrBased_simBased_subs_df is not None:
                chrBased_simBased_subs_df = chrBased_simBased_subs_df[chrBased_simBased_subs_df['Sample'].isin(samples_of_interest)]

            if chrBased_simBased_dinucs_df is not None:
                chrBased_simBased_dinucs_df = chrBased_simBased_dinucs_df[chrBased_simBased_dinucs_df['Sample'].isin(samples_of_interest)]

            if chrBased_simBased_indels_df is not None:
                chrBased_simBased_indels_df = chrBased_simBased_indels_df[chrBased_simBased_indels_df['Sample'].isin(samples_of_interest)]

        return chrbased_data_fill_signal_count_arrays_for_all_mutations_using_pyranges(occupancy_type,
                                                                        chrLong,
                                                                        simNum,
                                                                        chrBased_simBased_subs_df,
                                                                        chrBased_simBased_dinucs_df,
                                                                        chrBased_simBased_indels_df,
                                                                        chromSizesDict,
                                                                        chrBased_library_df,
                                                                        ordered_sbs_signatures,
                                                                        ordered_dbs_signatures,
                                                                        ordered_id_signatures,
                                                                        ordered_sbs_signatures_cutoffs,
                                                                        ordered_dbs_signatures_cutoffs,
                                                                        ordered_id_signatures_cutoffs,
                                                                        plusorMinus,
                                                                        discreet_mode,
                                                                        default_cutoff,
                                                                        log_file,
                                                                        verbose)
    except Exception as e:
        log_out = open(log_file, 'a')
        print("There is exception: %s --- chromosome: %s simulation_number: %s" % (e, chrLong, simNum), file=log_out)
        log_out.close()



# For apply_async
# Read chromBased simBased combined mutations df in the process
def chrbased_data_fill_signal_count_arrays_for_all_mutations_read_mutations(
        occupancy_type,
        occupancy_calculation_type,
        outputDir,
        jobname,
        chrLong,
        simNum,
        samples_of_interest,
        chromSizesDict,
        library_file_with_path,
        library_file_type,
        ordered_sbs_signatures,
        ordered_dbs_signatures,
        ordered_id_signatures,
        ordered_sbs_signatures_cutoffs,
        ordered_dbs_signatures_cutoffs,
        ordered_id_signatures_cutoffs,
        remove_outliers,
        plusorMinus,
        discreet_mode,
        default_cutoff,
        log_file,
        verbose):

    try:
        chrBased_simBased_subs_df, chrBased_simBased_dinucs_df, chrBased_simBased_indels_df = get_chrBased_simBased_dfs(outputDir, jobname, chrLong, simNum)

        # filter chrbased_df for samples_of_interest
        if samples_of_interest is not None:
            if chrBased_simBased_subs_df is not None:
                chrBased_simBased_subs_df = chrBased_simBased_subs_df[chrBased_simBased_subs_df['Sample'].isin(samples_of_interest)]

            if chrBased_simBased_dinucs_df is not None:
                chrBased_simBased_dinucs_df = chrBased_simBased_dinucs_df[chrBased_simBased_dinucs_df['Sample'].isin(samples_of_interest)]

            if chrBased_simBased_indels_df is not None:
                chrBased_simBased_indels_df = chrBased_simBased_indels_df[chrBased_simBased_indels_df['Sample'].isin(samples_of_interest)]


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
                                                                        remove_outliers,
                                                                        plusorMinus,
                                                                        discreet_mode,
                                                                        default_cutoff,
                                                                        log_file,
                                                                        verbose)
    except Exception as e:
        log_out = open(log_file, 'a')
        print("Exception: %s" % (e), file=log_out)
        log_out.close()



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
                                                                                  discreet_mode,
                                                                                  default_cutoff,
                                                                                  log_file,
                                                                                  verbose):

    chrBased_simBased_combined_df_split = get_chrBased_simBased_combined_df_split(outputDir, jobname, chrLong, simNum, splitIndex)

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
                                                                    discreet_mode,
                                                                    default_cutoff,
                                                                    log_file,
                                                                    verbose)


# For df_split
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
                                                             discreet_mode,
                                                             default_cutoff,
                                                             log_file,
                                                             verbose):


    if verbose:
        log_out = open(log_file, 'a')
        print('\tVerbose %s Worker pid %s memory_usage in %.2f MB START chrLong:%s simNum:%d' %(occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum), file=log_out)
        log_out.close()

    # 1st part  Prepare chr based mutations dataframes
    maximum_chrom_size = chromSizesDict[chrLong]
    start_time = time.time()

    chrBasedSignalArray = None #Will be filled from chrBasedSignal files if they exists
    library_file_opened_by_pyBigWig = None #Will be filled by pyBigWig from bigWig or bigBed
    my_upperBound = None
    signal_index = None

    if (chrBased_simBased_combined_df_split is not None) and verbose:
        log_out = open(log_file, 'a')
        print('\tVerbose %s Worker pid %s chrBased_mutations_df(%d,%d) ' %(occupancy_type,str(os.getpid()),chrBased_simBased_combined_df_split.shape[0],chrBased_simBased_combined_df_split.shape[1]), file=log_out)
        print('\tVerbose %s Worker pid %s memory_usage in %.2f MB Check1 Read Signal Array and Dataframes chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum), file=log_out)
        print('\tVerbose %s Worker pid %s -- signal_array_npy: %f in MB -- chrBased_simBased_combined_df_split: %f in MB -- chrLong:%s simNum:%d' % (
            occupancy_type,
            str(os.getpid()),
            sys.getsizeof(chrBasedSignalArray) / MEGABYTE_IN_BYTES,
            sys.getsizeof(chrBased_simBased_combined_df_split) / MEGABYTE_IN_BYTES,
            chrLong, simNum), file=log_out)
        log_out.close()

    libraryFilenameWoExtension = os.path.splitext(os.path.basename(library_file_with_path))[0]
    signalArrayFilename = '%s_signal_%s.npy' % (chrLong, libraryFilenameWoExtension)
    if (occupancy_type==NUCLEOSOMEOCCUPANCY):
        chrBasedSignalFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME, CHRBASED,signalArrayFilename)
    elif (occupancy_type== EPIGENOMICSOCCUPANCY):
        chrBasedSignalFile = os.path.join(outputDir,jobname,DATA,occupancy_type,LIB,CHRBASED,signalArrayFilename)
    else:
        #It can be EPIGENOMICSOCCUPANCY or user provided name e.g.: Epigenomics_ATAC_ENCFF317TWD
        chrBasedSignalFile = os.path.join(outputDir,jobname,DATA,occupancy_type,LIB,CHRBASED,signalArrayFilename)

    # Downloaded or created runtime
    if (os.path.exists(chrBasedSignalFile)):
        #Can this cause to deep sleep of processes?
        # chrBasedSignalArray = np.load(chrBasedSignalFile, mmap_mode='r')
        chrBasedSignalArray = np.load(chrBasedSignalFile)

    # If library_file_with_path is abs path and library_file_type is BIGWIG or BIGBED
    # For nucleosome_biosample==GM12878 or nucleosome_biosample==K562 library_file_with_path is only filename with extension, it is not absolute path
    if os.path.isabs(library_file_with_path):

        # Comment below to make it run in windows
        if (library_file_type == BIGWIG):
            try:
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
                log_out = open(log_file, 'a')
                print('Exception %s' %library_file_with_path, file=log_out)
                log_out.close()

        elif (library_file_type == BIGBED):
            try:
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
                log_out = open(log_file, 'a')
                print('Exception %s' %library_file_with_path, file=log_out)
                log_out.close()


    if ((chrBasedSignalArray is not None) or ((library_file_opened_by_pyBigWig is not None) and (chrLong in library_file_opened_by_pyBigWig.chroms()))):
        ######################################################## #######################
        ################### Fill signal and count array starts ########################
        ###############################################################################
        if verbose:
            log_out = open(log_file, 'a')
            print('\tVerbose %s Worker pid %s memory_usage in %.2f MB Check2_1 Start chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum), file=log_out)
            log_out.close()

        if ((chrBased_simBased_combined_df_split is not None) and (not chrBased_simBased_combined_df_split.empty)):

            # df_columns is a numpy array
            df_columns = chrBased_simBased_combined_df_split.columns.values

            ###############################################################################
            ################################ Initialization ###############################
            ###############################################################################
            subsSignatures_mask_array = np.isin(df_columns,ordered_sbs_signatures)
            dinucsSignatures_mask_array = np.isin(df_columns,ordered_dbs_signatures)
            indelsSignatures_mask_array = np.isin(df_columns,ordered_id_signatures)

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
                discreet_mode,
                default_cutoff,
                df_columns,
                occupancy_calculation_type) for row in chrBased_simBased_combined_df_split[df_columns].values]

        if verbose:
            log_out = open(log_file, 'a')
            print('\tVerbose %s Worker pid %s memory_usage in %.2f MB Check2_2 End chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum), file=log_out)
            log_out.close()
        ###############################################################################
        ################### Fill signal and count array ends ##########################
        ###############################################################################

    if (library_file_opened_by_pyBigWig is not None):
        library_file_opened_by_pyBigWig.close()

    if verbose:
        log_out = open(log_file, 'a')
        print('\tVerbose %s Worker pid %s memory_usage in %.2f MB END  chrLong:%s simNum:%d' % (occupancy_type,str(os.getpid()), memory_usage(), chrLong, simNum), file=log_out)
        print('\tVerbose %s Worker pid %s took %f seconds chrLong:%s simNum:%d\n' % (occupancy_type,str(os.getpid()), (time.time() - start_time), chrLong, simNum), file=log_out)
        log_out.close()

    ###############################################################################
    ################### Return  starts ############################################
    ###############################################################################

    # Initialzie the list, you will return this list
    SignalArrayAndCountArrayList = []

    # new way
    SignalArrayAndCountArrayList.append(chrLong)
    SignalArrayAndCountArrayList.append(simNum)
    SignalArrayAndCountArrayList.append(subsSignature_accumulated_signal_np_array)
    SignalArrayAndCountArrayList.append(dinucsSignature_accumulated_signal_np_array)
    SignalArrayAndCountArrayList.append(indelsSignature_accumulated_signal_np_array)
    SignalArrayAndCountArrayList.append(subsSignature_accumulated_count_np_array)
    SignalArrayAndCountArrayList.append(dinucsSignature_accumulated_count_np_array)
    SignalArrayAndCountArrayList.append(indelsSignature_accumulated_count_np_array)

    return SignalArrayAndCountArrayList


# For df_split
# Vectorization
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
        discreet_mode,
        default_cutoff,
        df_columns,
        occupancy_calculation_type):

    indexofType = np.where(df_columns == TYPE)[0][0]
    indexofStart = np.where(df_columns == START)[0][0]

    mutation_row_type = row[indexofType]
    mutation_row_start = row[indexofStart]

    if mutation_row_type == SUBS:
        accumulated_signal_np_array = subsSignature_accumulated_signal_np_array
        accumulated_count_np_array = subsSignature_accumulated_count_np_array
        cutoffs = ordered_sbs_signatures_cutoffs
        signatures_mask_array = subsSignatures_mask_array
    elif mutation_row_type == DINUCS:
        accumulated_signal_np_array = dinucsSignature_accumulated_signal_np_array
        accumulated_count_np_array = dinucsSignature_accumulated_count_np_array
        cutoffs = ordered_dbs_signatures_cutoffs
        signatures_mask_array = dinucsSignatures_mask_array
    elif mutation_row_type == INDELS:
        accumulated_signal_np_array = indelsSignature_accumulated_signal_np_array
        accumulated_count_np_array = indelsSignature_accumulated_count_np_array
        cutoffs = ordered_id_signatures_cutoffs
        signatures_mask_array = indelsSignatures_mask_array

    window_array = None
    windowSize = plusOrMinus * 2 + 1

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

    # Get or fill window_array using Case1, Case2, and Case3
    # Case 1: start is very close to the chromosome start
    if (mutation_row_start < plusOrMinus):
        # print('Case 1: start is very close to the chromosome start --- mutation[Start]:%d' %(mutation_row_start))
        #Faster
        if (chrBasedSignalArray is not None):
            window_array = chrBasedSignalArray[0:(mutation_row_start + plusOrMinus + 1)]
            window_array = np.pad(window_array, (plusOrMinus - mutation_row_start, 0), 'constant', constant_values=(0, 0))

        elif (library_file_type == BIGWIG):
            #Important: The bigWig format does not support overlapping intervals.
            window_array=library_file_opened_by_pyBigWig.values(chrLong,0,(mutation_row_start+plusOrMinus+1),numpy=True)
            # How do you handle outliers?
            window_array[np.isnan(window_array)] = 0
            window_array[window_array>my_upperBound]=my_upperBound
            window_array = np.pad(window_array, (plusOrMinus - mutation_row_start, 0), 'constant',constant_values=(0, 0))

        elif (library_file_type == BIGBED):
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

        elif (library_file_type == BIGWIG):
            # Important: The bigWig format does not support overlapping intervals.
            window_array = library_file_opened_by_pyBigWig.values(chrLong,(mutation_row_start-plusOrMinus),maximum_chrom_size,numpy=True)
            # How do you handle outliers?
            window_array[np.isnan(window_array)] = 0
            window_array[window_array>my_upperBound]=my_upperBound
            window_array = np.pad(window_array, (0,mutation_row_start+plusOrMinus-maximum_chrom_size+1),'constant',constant_values=(0,0))

        elif (library_file_type == BIGBED):
            # print('Case2 %s mutation_row[START]:%d mutation_row[START]-plusOrMinus:%d maximum_chrom_size:%d' %(chrLong,mutation_row[START],mutation_row[START]-plusOrMinus,maximum_chrom_size))
            if ((mutation_row_start-plusOrMinus)<maximum_chrom_size):
                list_of_entries=library_file_opened_by_pyBigWig.entries(chrLong,(mutation_row_start-plusOrMinus),maximum_chrom_size)
                if list_of_entries is not None:
                    window_array = np.zeros((windowSize,),dtype=np.float32)
                    # We did not handle outliers for BigBed files.
                    [(func_addSignal(window_array, entry[0], entry[1], np.float32(entry[2].split()[signal_index]),mutation_row_start,plusOrMinus) if len(entry) >= 3 else (func_addSignal(window_array, entry[0], entry[1],1, mutation_row_start,plusOrMinus))) for entry in list_of_entries]

    # Case 3: No problem
    else:
        if (chrBasedSignalArray is not None):
            window_array = chrBasedSignalArray[(mutation_row_start-plusOrMinus):(mutation_row_start+plusOrMinus+1)]

        elif (library_file_type == BIGWIG):
            # Important: You have to go over intervals if there are overlapping intervals.
            window_array = library_file_opened_by_pyBigWig.values(chrLong, (mutation_row_start-plusOrMinus), (mutation_row_start+plusOrMinus+1),numpy=True)
            # How do you handle outliers?
            window_array[np.isnan(window_array)] = 0
            window_array[window_array>my_upperBound]=my_upperBound

        elif (library_file_type == BIGBED):
            # print('Case3 %s mutation_row[START]:%d mutation_row[START]-plusOrMinus:%d mutation_row[START]+plusOrMinus+1:%d' %(chrLong,mutation_row[START],mutation_row[START]-plusOrMinus,mutation_row[START]+plusOrMinus+1))
            if ((mutation_row_start+plusOrMinus+1)<=maximum_chrom_size):
                list_of_entries=library_file_opened_by_pyBigWig.entries(chrLong, (mutation_row_start-plusOrMinus), (mutation_row_start+plusOrMinus+1))
                if list_of_entries is not None:
                    window_array = np.zeros((windowSize,),dtype=np.float32)
                    # We did not handle outliers for BigBed files.
                    [(func_addSignal(window_array, entry[0], entry[1], np.float32(entry[2].split()[signal_index]),mutation_row_start,plusOrMinus) if len(entry) >= 3 else (func_addSignal(window_array, entry[0], entry[1],1, mutation_row_start,plusOrMinus))) for entry in list_of_entries]


    # Get the sample at this mutation_row
    # sample = mutation_row_sample
    # simulationNumber= mutation_row_simulation_number

    # NO SIGNAL case is added
    # Vectorize
    # Fill numpy arrays using window_array
    if (window_array is not None) and (np.any(window_array)):
        probabilities = row[signatures_mask_array]

        if discreet_mode:
            threshold_mask_array = np.greater_equal(probabilities, cutoffs)
            #Convert True into 1, and False into 0
            mask_array = threshold_mask_array.astype(int)
            #Add 1 for the Aggregated analysis to the mask array
            mask_array = np.append(mask_array, 1)
        else:
            probabilities[probabilities < default_cutoff] = 0
            mask_array = np.array(probabilities).astype(float)
            # Add 1 for the Aggregated analysis to the mask array
            mask_array = np.append(mask_array, 1.0)

        # Add one more dimension to window_array and mask_array
        window_array_1x2001=np.array([window_array])
        mask_array_1xnumofsignatures=np.array([mask_array])

        to_be_accumulated_signal_array = mask_array_1xnumofsignatures.T * window_array_1x2001
        accumulated_signal_np_array += to_be_accumulated_signal_array

        # count is prob * 1 case otherwise y-axis range becomes quite low
        to_be_accumulated_count_array = mask_array_1xnumofsignatures.T * (window_array_1x2001 > 0).astype(int) # legacy code

        # default
        if occupancy_calculation_type == MISSING_SIGNAL:
            # accumulated_count_np_array += (to_be_accumulated_signal_array>0)
            accumulated_count_np_array += to_be_accumulated_count_array
        else:
            accumulated_count_np_array += 1



def get_window_array(mutation_row_start,
                     plusOrMinus,
                     chrLong,
                     chrBasedSignalArray,
                     library_file_type,
                     library_file_opened_by_pyBigWig,
                     signal_index,
                     my_upperBound,
                     maximum_chrom_size):

    window_array = None
    windowSize = plusOrMinus*2+1

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

    # Get or fill window_array using Case1, Case2, and Case3
    # Case 1: start is very close to the chromosome start
    if (mutation_row_start < plusOrMinus):
        # print('Case 1: start is very close to the chromosome start --- mutation[Start]:%d' %(mutation_row_start))
        # Faster
        if (chrBasedSignalArray is not None):
            window_array = chrBasedSignalArray[0:(mutation_row_start + plusOrMinus + 1)]
            window_array = np.pad(window_array, (plusOrMinus - mutation_row_start, 0), 'constant', constant_values=(0, 0))

        elif (library_file_type == BIGWIG):
            # Important: The bigWig format does not support overlapping intervals.
            # values() method return all the values between start and end, of size end-start maybe end-start+1
            # intervals() method return all the intervals between start and end, of size variable

            try:
                if ((mutation_row_start+plusOrMinus+1) <= maximum_chrom_size):
                    window_array = library_file_opened_by_pyBigWig.values(chrLong,0,(mutation_row_start+plusOrMinus+1), numpy=True)
            except Exception as e:
                print('Exception in get_window_array %s %s mutation_row_start: %d maximum_chrom_size: %d' %(str(e), chrLong, mutation_row_start, maximum_chrom_size))

            if (window_array is None) or np.sum(np.isnan(window_array))>0:
                return None

            # How do you handle outliers?
            window_array[np.isnan(window_array)] = 0
            window_array[window_array>my_upperBound] = my_upperBound
            window_array = np.pad(window_array, (plusOrMinus - mutation_row_start, 0), 'constant',constant_values=(0, 0))

        elif (library_file_type == BIGBED):
            # We assume that in the 7th column there is signal data
            list_of_entries = library_file_opened_by_pyBigWig.entries(chrLong,0,(mutation_row_start+plusOrMinus+1))

            if list_of_entries is not None:
                window_array = np.zeros((windowSize,),dtype=np.float32)
                # We did not handle outliers for BigBed files.

                # From DNA methylation get the 7th

                # library_file_bed_format==BED_6PLUS4):
                # (713235, 713435, 'Peak_40281\t15\t.\t3.48949\t5.67543\t3.79089\t158')
                # signal_index=3

                # library_file_bed_format==BED_9PLUS2):
                # [(10810, 10811, 'MCF7_NoStarve_B1__GC_\t3\t+\t10810\t10811\t255,0,0\t3\t100'), (10812, 10813, 'MCF7_NoStarve_B1__GC_\t3\t+\t10812\t10813\t255,0,0\t3\t100'), (10815, 10816, 'MCF7_NoStarve_B1__GC_\t3\t+\t10815\t10816\t0,255,0\t3\t0')]
                # signal_index=7
                [(func_addSignal(window_array, entry[0], entry[1], np.float32(entry[2].split()[signal_index]),mutation_row_start, plusOrMinus) if len(entry) >= 3 else (func_addSignal(window_array, entry[0], entry[1], 1, mutation_row_start, plusOrMinus))) for entry in list_of_entries]

    # Case 2: start is very close to the chromosome end
    elif (mutation_row_start+plusOrMinus+1 > maximum_chrom_size):
        # print('Case2: start is very close to the chromosome end ---  mutation[Start]:%d' %(mutation_row_start))
        if ((chrBasedSignalArray is not None)):
            window_array = chrBasedSignalArray[(mutation_row_start-plusOrMinus):maximum_chrom_size]
            window_array = np.pad(window_array, (0,mutation_row_start+plusOrMinus-maximum_chrom_size+1),'constant',constant_values=(0,0))

        elif (library_file_type == BIGWIG):
            #Important: The bigWig format does not support overlapping intervals.
            try:
                if ((mutation_row_start-plusOrMinus) <= maximum_chrom_size):
                    window_array = library_file_opened_by_pyBigWig.values(chrLong,(mutation_row_start-plusOrMinus),maximum_chrom_size, numpy=True)
            except Exception as e:
                print('Exception in get_window_array %s %s mutation_row_start: %d maximum_chrom_size: %d' %(str(e), chrLong, mutation_row_start, maximum_chrom_size))

            if (window_array is None) or np.sum(np.isnan(window_array))>0:
                return None

            # How do you handle outliers?
            window_array[np.isnan(window_array)] = 0
            window_array[window_array>my_upperBound] = my_upperBound
            window_array = np.pad(window_array, (0,mutation_row_start+plusOrMinus-maximum_chrom_size+1),'constant',constant_values=(0,0))

        elif (library_file_type == BIGBED):
            # print('Case2 %s mutation_row[START]:%d mutation_row[START]-plusOrMinus:%d maximum_chrom_size:%d' %(chrLong,mutation_row[START],mutation_row[START]-plusOrMinus,maximum_chrom_size))
            if ((mutation_row_start-plusOrMinus)<maximum_chrom_size):
                list_of_entries = library_file_opened_by_pyBigWig.entries(chrLong,(mutation_row_start-plusOrMinus),maximum_chrom_size)

                if list_of_entries is not None:
                    window_array = np.zeros((windowSize,),dtype=np.float32)
                    # We did not handle outliers for BigBed files.
                    [(func_addSignal(window_array, entry[0], entry[1], np.float32(entry[2].split()[signal_index]),mutation_row_start,plusOrMinus) if len(entry) >= 3 else (func_addSignal(window_array, entry[0], entry[1],1, mutation_row_start,plusOrMinus))) for entry in list_of_entries]

    # Case 3: No problem
    else:
        if (chrBasedSignalArray is not None):
            window_array = chrBasedSignalArray[(mutation_row_start-plusOrMinus):(mutation_row_start+plusOrMinus+1)]

        elif (library_file_type == BIGWIG):
            # Important: You have to go over intervals if there are overlapping intervals.
            try:
                if ((mutation_row_start-plusOrMinus) <= maximum_chrom_size) and ((mutation_row_start+plusOrMinus+1) <= maximum_chrom_size):
                    window_array = library_file_opened_by_pyBigWig.values(chrLong, (mutation_row_start-plusOrMinus), (mutation_row_start+plusOrMinus+1), numpy=True)
            except Exception as e:
                print('Exception in get_window_array %s %s mutation_row_start: %d maximum_chrom_size: %d' %(str(e), chrLong, mutation_row_start, maximum_chrom_size))

            if (window_array is None) or np.sum(np.isnan(window_array))>0:
                return None

            # How do you handle outliers?
            window_array[np.isnan(window_array)] = 0
            window_array[window_array>my_upperBound] = my_upperBound

        elif (library_file_type == BIGBED):
            # print('Case3 %s mutation_row[START]:%d mutation_row[START]-plusOrMinus:%d mutation_row[START]+plusOrMinus+1:%d' %(chrLong,mutation_row[START],mutation_row[START]-plusOrMinus,mutation_row[START]+plusOrMinus+1))
            if ((mutation_row_start+plusOrMinus+1) <= maximum_chrom_size):
                list_of_entries = library_file_opened_by_pyBigWig.entries(chrLong, (mutation_row_start-plusOrMinus), (mutation_row_start+plusOrMinus+1))

                if list_of_entries is not None:
                    window_array = np.zeros((windowSize,),dtype=np.float32)
                    # We did not handle outliers for BigBed files.
                    [(func_addSignal(window_array, entry[0], entry[1], np.float32(entry[2].split()[signal_index]),mutation_row_start,plusOrMinus) if len(entry) >= 3 else (func_addSignal(window_array, entry[0], entry[1],1, mutation_row_start,plusOrMinus))) for entry in list_of_entries]

    return window_array

# Main engine function for accumulate signals and counts
def accumulate_arrays(row,
                      window_array,
                      signatures_mask_array,
                      ordered_signatures_cutoffs,
                      occupancy_calculation_type,
                      accumulated_signal_np_array,
                      accumulated_count_np_array,
                      discreet_mode,
                      default_cutoff):

    # September 18, 2020 NO SIGNAL case is added
    # Vectorize July 25, 2020
    # Fill numpy arrays using window_array
    if (window_array is not None) and (np.any(window_array)):
        probabilities = row[signatures_mask_array]

        if discreet_mode:
            # Discreet way 1 or 0
            # Convert True into 1, and False into 0
            threshold_mask_array = np.greater_equal(probabilities, ordered_signatures_cutoffs)
            mask_array = threshold_mask_array.astype(int)
            # Add 1 for the Aggregated analysis to the mask array
            mask_array = np.append(mask_array, 1) # For discreet 1
        else:
            probabilities[probabilities < default_cutoff ] = 0
            mask_array = np.array(probabilities).astype(float)
            # Add 1 for the Aggregated analysis to the mask array
            mask_array = np.append(mask_array, 1.0)

        # Add one more dimension to mask_array: (num_of_signatures,) --> (1, num_of_signatures)
        mask_array_1xnumofsignatures = np.expand_dims(mask_array, axis=0)

        # Add one more dimension to window_array: (2001,) --> (1, 2001)
        window_array_1x2001 = np.expand_dims(window_array, axis=0)

        # to_be_accumulated_signal_array: (num_of_signatures, 1) * (1,2001) = (num_of_signatures, 2001)
        to_be_accumulated_signal_array = mask_array_1xnumofsignatures.T * window_array_1x2001
        accumulated_signal_np_array += to_be_accumulated_signal_array

        # question: if signal > 0, is count 1 or prob * 1 ?
        # count is taken 1, otherwise with prob * 1, y-axis range becomes quite low
        # if signal is greater than 0, count is 1, in fact count can be more than 1 if there are multiple overlaps
        to_be_accumulated_count_array = mask_array_1xnumofsignatures.T * (window_array_1x2001 > 0).astype(int)

        # Default
        if occupancy_calculation_type == MISSING_SIGNAL:
            accumulated_count_np_array += to_be_accumulated_count_array
            # For debugging uncomment the following line
            # print('to_be_accumulated_signal_array:', to_be_accumulated_signal_array,
            #       'to_be_accumulated_signal_array.shape:', to_be_accumulated_signal_array.shape,
            #       'np.sum(to_be_accumulated_signal_array):', np.sum(to_be_accumulated_signal_array),
            #       'np.count_nonzero(to_be_accumulated_signal_array):', np.count_nonzero(to_be_accumulated_signal_array),
            #       'to_be_accumulated_count_array:', to_be_accumulated_count_array,
            #       'to_be_accumulated_count_array.shape:', to_be_accumulated_count_array.shape,
            #       'np.sum(to_be_accumulated_count_array):', np.sum(to_be_accumulated_count_array),
            #       'np.count_nonzero(to_be_accumulated_count_array):', np.count_nonzero(to_be_accumulated_count_array),
            #       'equal:', np.array_equal(to_be_accumulated_signal_array, to_be_accumulated_count_array))
        else:
            # this part is not tested
            accumulated_count_np_array += 1




# Vectorization
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
        discreet_mode,
        default_cutoff,
        df_columns,
        occupancy_calculation_type):

    indexofStart = np.where(df_columns == START)[0][0]
    mutation_row_start = row[indexofStart]

    # added for debug starts
    simNum = np.where(df_columns == SIMULATION_NUMBER)[0][0]
    sim_num = row[simNum]
    # added for debug ends

    window_array = get_window_array(mutation_row_start,
                     plusOrMinus,
                     chrLong,
                     chrBasedSignalArray,
                     library_file_type,
                     library_file_opened_by_pyBigWig,
                     signal_index,
                     my_upperBound,
                     maximum_chrom_size)

    if window_array is not None:
        accumulate_arrays(row,
                          window_array,
                          signatures_mask_array,
                          ordered_signatures_cutoffs,
                          occupancy_calculation_type,
                          accumulated_signal_np_array,
                          accumulated_count_np_array,
                          discreet_mode,
                          default_cutoff)




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

    # go back
    os.chdir(current_abs_path)


def occupancy_analysis_memory_efficient(genome,
                      computation_type,
                      occupancy_type,
                      occupancy_calculation_type,
                      sample_based,
                      plus_minus_epigenomics,
                      chromSizesDict,
                      chromNamesList,
                      outputDir,
                      jobname,
                      numofSimulations,
                      samples_of_interest,
                      library_file_with_path,
                      library_file_memo,
                      ordered_sbs_signatures,
                      ordered_dbs_signatures,
                      ordered_id_signatures,
                      ordered_sbs_signatures_cutoffs,
                      ordered_dbs_signatures_cutoffs,
                      ordered_id_signatures_cutoffs,
                      remove_outliers,
                      quantile_value,
                      discreet_mode,
                      default_cutoff,
                      parallel_mode,
                      log_file,
                      verbose):

    log_out = open(log_file, 'a')
    if (os.path.basename(library_file_with_path) not in SIGPROFILERTOPOGRAPHY_DEFAULT_FILES) and (not os.path.exists(library_file_with_path)):
        print('There is no such file under %s' %(library_file_with_path), file=log_out)

    print('\n#################################################################################', file=log_out)
    print('--- %s Analysis starts' %(occupancy_type), file=log_out)
    print('--- Computation Type:%s' % (computation_type), file=log_out)
    print('--- Occupancy Type:%s' % (occupancy_type), file=log_out)
    print('--- Library file with path: %s\n' %library_file_with_path, file=log_out)
    log_out.close()

    # Using pyBigWig for BigWig and BigBed files if you can import pyBigWig (linux only) otherwise no
    # By the way pyBigWig can be imported in unix, linux like os not available in windows
    # Using HM and CTCF bed files preparing chr based signal array runtime
    # Using ATAC-seq wig files preparing chr based signal array runtime

    ##########################################################################
    # If chunksize is 1, maxtasksperchild=x will call the function x times in each process,
    # but if chunksize is y, it will call the function x*y times in each process.
    # Setting maxtaskperchild to 1 would restart each process in your pool after it processed a single task, which is the most aggressive setting you could use to free any leaked resources.
    # numofProcesses = multiprocessing.cpu_count()
    # pool = multiprocessing.Pool(numofProcesses, maxtasksperchild=1)
    ##########################################################################

    number_of_sbs_signatures = ordered_sbs_signatures.size
    number_of_dbs_signatures = ordered_dbs_signatures.size
    number_of_id_signatures = ordered_id_signatures.size

    # at index=0 we have the first signature of the ordered_XX_signatures
    # at index=1 we have the second signature of the ordered_XX_signatures
    # at index=-1 we have the aggregated mutations
    subsSignatures = np.append(ordered_sbs_signatures, AGGREGATEDSUBSTITUTIONS)
    dinucsSignatures = np.append(ordered_dbs_signatures, AGGREGATEDDINUCS)
    indelsSignatures = np.append(ordered_id_signatures, AGGREGATEDINDELS)

    # For Vectorization
    # These are used in writing tables
    allSims_subsSignature_accumulated_signal_np_array = np.zeros((numofSimulations+1, number_of_sbs_signatures + 1, plus_minus_epigenomics * 2 + 1))
    allSims_dinucsSignature_accumulated_signal_np_array = np.zeros((numofSimulations+1, number_of_dbs_signatures + 1, plus_minus_epigenomics * 2 + 1))
    allSims_indelsSignature_accumulated_signal_np_array = np.zeros((numofSimulations+1, number_of_id_signatures + 1, plus_minus_epigenomics * 2 + 1))

    allSims_subsSignature_accumulated_count_np_array = np.zeros((numofSimulations+1, number_of_sbs_signatures + 1, plus_minus_epigenomics * 2 + 1))
    allSims_dinucsSignature_accumulated_count_np_array = np.zeros((numofSimulations+1, number_of_dbs_signatures + 1, plus_minus_epigenomics * 2 + 1))
    allSims_indelsSignature_accumulated_count_np_array = np.zeros((numofSimulations+1, number_of_id_signatures+1, plus_minus_epigenomics * 2 + 1))

    # This code reads the file and prepare chrbased signal files
    # If file is in default files,chr based signal files are downloaded from ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerTopography/lib/
    # No need for preparing here
    library_file_type = None
    chrBased_library_df = None
    library_grouped_df = None

    if (os.path.basename(library_file_with_path) not in SIGPROFILERTOPOGRAPHY_DEFAULT_FILES):
        # What is the type of the signal_file_with_path?
        # If it is a bed file read signal_file_with_path here
        file_extension = os.path.splitext(os.path.basename(library_file_with_path))[1]

        if ((file_extension.lower() == '.bigwig') or (file_extension.lower() == '.bw')):
            library_file_type = BIGWIG
        elif ((file_extension.lower() == '.bigbed') or (file_extension.lower() == '.bb')):
            library_file_type = BIGBED
        elif (file_extension.lower() == '.bed'):
            library_file_type = BED

            if os.path.exists(library_file_with_path):
                bedfilename, bedfile_extention = os.path.splitext(os.path.basename(library_file_with_path))

                if (bedfile_extention.lower() == '.bed' or bedfile_extention.lower() == '.np' or bedfile_extention.lower() == '.narrowpeak'):
                    discard_signal = False
                    library_file_df, max_signal, min_signal = readFileInBEDFormat(library_file_with_path, discard_signal, log_file)

                    # Eliminate Outliers
                    if ((remove_outliers == True) and (quantile_value < 1.0)):
                        # remove the outliers
                        q = library_file_df[SIGNAL].quantile(quantile_value)
                        log_out = open(log_file, 'a')
                        print('Signal greater than %f is removed for quantile: %f' % (q, quantile_value), file=log_out)
                        print('before outlier removal number of rows: %d' % (library_file_df.shape[0]), file=log_out)
                        library_file_df = library_file_df[library_file_df[SIGNAL] < q]
                        print('after outlier removal number of rows: %d' % (library_file_df.shape[0]), file=log_out)
                        log_out.close()

                    library_grouped_df = library_file_df.groupby(CHROM)

        elif ((file_extension.lower() == '.narrowpeak') or (file_extension.lower() == '.np')):
            library_file_type = NARROWPEAK
            readBEDandWriteChromBasedSignalArrays(outputDir, jobname, genome, library_file_with_path, occupancy_type, quantile_value, remove_outliers, log_file)
        elif (file_extension.lower() == '.wig'):
            library_file_type = WIG
            # For inhouse preparation
            # readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArraysSequentially(genome, quantileValue,library_file_with_path)
            # readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArraysInParallel(genome, quantileValue,library_file_with_path)
            filetype_BEDGRAPH = isFileTypeBedGraph(library_file_with_path)
            if filetype_BEDGRAPH:
                if verbose: start_time = time.time()
                # Read by chunks
                # readWig_write_derived_from_bedgraph_using_pool_chunks(outputDir, jobname, genome, library_file_with_path,occupancy_type,remove_outliers,verbose)
                # Read at once
                readWig_write_derived_from_bedgraph_using_pool_read_all(outputDir, jobname, genome, library_file_with_path, occupancy_type,remove_outliers,verbose,quantile_value, log_file)
                if verbose:
                    log_out = open(log_file, 'a')
                    print('\tVerbose Read wig file and write chrbased arrays took %f seconds' %((time.time() - start_time)), file=log_out)
                    log_out.close()

                # For 6 GB ATAC-seq file using pool took 8 min whereas without pool it took 16 min.
                # start_time = time.time()
                # readWig_write_derived_from_bedgraph(outputDir, jobname, genome, library_file_with_path,occupancy_type,verbose)
                # print('Without pool Took %f seconds' %((time.time() - start_time)))
            else:
                readWig_with_fixedStep_variableStep_writeChrBasedSignalArrays(outputDir, jobname, genome, library_file_with_path, occupancy_type, quantile_value, remove_outliers, log_file)
        elif (file_extension.lower() == '.bedgraph'):
            library_file_type = BEDGRAPH
            readWig_write_derived_from_bedgraph_using_pool_read_all(outputDir, jobname, genome, library_file_with_path,occupancy_type, remove_outliers, verbose, quantile_value, log_file)
        else:
            library_file_type = LIBRARY_FILE_TYPE_OTHER


    def accumulate_apply_async_result_vectorization(simulatonBased_SignalArrayAndCountArrayList):
        try:
            if simulatonBased_SignalArrayAndCountArrayList is not None:
                chrLong = simulatonBased_SignalArrayAndCountArrayList[0]
                simNum = simulatonBased_SignalArrayAndCountArrayList[1]
                subsSignature_accumulated_signal_np_array = simulatonBased_SignalArrayAndCountArrayList[2]
                dinucsSignature_accumulated_signal_np_array = simulatonBased_SignalArrayAndCountArrayList[3]
                indelsSignature_accumulated_signal_np_array = simulatonBased_SignalArrayAndCountArrayList[4]
                subsSignature_accumulated_count_np_array = simulatonBased_SignalArrayAndCountArrayList[5]
                dinucsSignature_accumulated_count_np_array = simulatonBased_SignalArrayAndCountArrayList[6]
                indelsSignature_accumulated_count_np_array = simulatonBased_SignalArrayAndCountArrayList[7]

                # Accumulation
                allSims_subsSignature_accumulated_signal_np_array[simNum] += subsSignature_accumulated_signal_np_array
                allSims_dinucsSignature_accumulated_signal_np_array[simNum] += dinucsSignature_accumulated_signal_np_array
                allSims_indelsSignature_accumulated_signal_np_array[simNum] += indelsSignature_accumulated_signal_np_array

                allSims_subsSignature_accumulated_count_np_array[simNum] += subsSignature_accumulated_count_np_array
                allSims_dinucsSignature_accumulated_count_np_array[simNum] += dinucsSignature_accumulated_count_np_array
                allSims_indelsSignature_accumulated_count_np_array[simNum] += indelsSignature_accumulated_count_np_array
                # print('ACCUMULATION chrLong:%s simNum:%d ENDS' %(chrLong,simNum))

        except Exception as e:
            print("Exception in accumulate_apply_async_result_vectorization function: %s" % (e))

    sim_nums = range(0, numofSimulations + 1)
    sim_num_chr_tuples = ((sim_num, chrLong) for sim_num in sim_nums for chrLong in chromNamesList)

    if parallel_mode:
        jobs = []

        if (computation_type == USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM):
            numofProcesses = multiprocessing.cpu_count()
            pool = multiprocessing.Pool(processes=numofProcesses)

            for simNum, chrLong in sim_num_chr_tuples:
                jobs.append(pool.apply_async(chrbased_data_fill_signal_count_arrays_for_all_mutations_read_mutations_memory_efficient,
                                             args=(occupancy_type,
                                                   occupancy_calculation_type,
                                                   outputDir,
                                                   jobname,
                                                   chrLong,
                                                   simNum,
                                                   samples_of_interest,
                                                   chromSizesDict,
                                                   library_grouped_df,
                                                   chrBased_library_df,
                                                   library_file_with_path,
                                                   library_file_type,
                                                   ordered_sbs_signatures,
                                                   ordered_dbs_signatures,
                                                   ordered_id_signatures,
                                                   ordered_sbs_signatures_cutoffs,
                                                   ordered_dbs_signatures_cutoffs,
                                                   ordered_id_signatures_cutoffs,
                                                   plus_minus_epigenomics,
                                                   discreet_mode,
                                                   default_cutoff,
                                                   log_file,
                                                   verbose,),
                                             callback=accumulate_apply_async_result_vectorization))

            pool.close()
            pool.join()

    else:
        # Sequential mode for profiling, debugging and testing purposes
        for simNum, chrLong in sim_num_chr_tuples:
            simulatonBased_SignalArrayAndCountArrayList = chrbased_data_fill_signal_count_arrays_for_all_mutations_read_mutations_memory_efficient(occupancy_type,
                                                   occupancy_calculation_type,
                                                   outputDir,
                                                   jobname,
                                                   chrLong,
                                                   simNum,
                                                   samples_of_interest,
                                                   chromSizesDict,
                                                   library_grouped_df,
                                                   chrBased_library_df,
                                                   library_file_with_path,
                                                   library_file_type,
                                                   ordered_sbs_signatures,
                                                   ordered_dbs_signatures,
                                                   ordered_id_signatures,
                                                   ordered_sbs_signatures_cutoffs,
                                                   ordered_dbs_signatures_cutoffs,
                                                   ordered_id_signatures_cutoffs,
                                                   plus_minus_epigenomics,
                                                   discreet_mode,
                                                   default_cutoff,
                                                   log_file,
                                                   verbose)

            accumulate_apply_async_result_vectorization(simulatonBased_SignalArrayAndCountArrayList)


    # Same for parallel or sequential run
    writeSimulationBasedAverageOccupancyUsingNumpyArray(occupancy_type,
                                                   sample_based,
                                                   plus_minus_epigenomics,
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

    log_out = open(log_file, 'a')
    print('--- %s Analysis ends' %(occupancy_type), file=log_out)
    print('#################################################################################\n', file=log_out)
    log_out.close()



def occupancy_analysis_using_pyranges(genome,
                      computation_type,
                      occupancy_type,
                      occupancy_calculation_type,
                      sample_based,
                      plus_minus_epigenomics,
                      chromSizesDict,
                      chromNamesList,
                      outputDir,
                      jobname,
                      numofSimulations,
                      samples_of_interest,
                      library_file_with_path,
                      library_file_memo,
                      ordered_sbs_signatures,
                      ordered_dbs_signatures,
                      ordered_id_signatures,
                      ordered_sbs_signatures_cutoffs,
                      ordered_dbs_signatures_cutoffs,
                      ordered_id_signatures_cutoffs,
                      remove_outliers,
                      quantile_value,
                      discreet_mode,
                      default_cutoff,
                      parallel_mode,
                      log_file,
                      verbose):

    log_out = open(log_file, 'a')
    if (os.path.basename(library_file_with_path) not in SIGPROFILERTOPOGRAPHY_DEFAULT_FILES) and (not os.path.exists(library_file_with_path)):
        print('There is no such file under %s' %(library_file_with_path), file=log_out)

    print('\n#################################################################################', file=log_out)
    print('--- %s Analysis starts' %(occupancy_type), file=log_out)
    print('--- Computation Type:%s' % (computation_type), file=log_out)
    print('--- Occupancy Type:%s' % (occupancy_type), file=log_out)
    print('--- Library file with path: %s\n' %library_file_with_path, file=log_out)
    log_out.close()

    # Using pyBigWig for BigWig and BigBed files if you can import pyBigWig (linux only) otherwise no
    # By the way pyBigWig can be imported in unix, linux like os not available in windows
    # Using HM and CTCF bed files preparing chr based signal array runtime
    # Using ATAC-seq wig files preparing chr based signal array runtime

    ##########################################################################
    # If chunksize is 1, maxtasksperchild=x will call the function x times in each process,
    # but if chunksize is y, it will call the function x*y times in each process.
    # Setting maxtaskperchild to 1 would restart each process in your pool after it processed a single task, which is the most aggressive setting you could use to free any leaked resources.
    # numofProcesses = multiprocessing.cpu_count()
    # pool = multiprocessing.Pool(numofProcesses, maxtasksperchild=1)
    ##########################################################################

    number_of_sbs_signatures = ordered_sbs_signatures.size
    number_of_dbs_signatures = ordered_dbs_signatures.size
    number_of_id_signatures = ordered_id_signatures.size

    # at index=0 we have the first signature of the ordered_XX_signatures
    # at index=1 we have the second signature of the ordered_XX_signatures
    # at index=-1 we have the aggregated mutations
    subsSignatures = np.append(ordered_sbs_signatures, AGGREGATEDSUBSTITUTIONS)
    dinucsSignatures = np.append(ordered_dbs_signatures, AGGREGATEDDINUCS)
    indelsSignatures = np.append(ordered_id_signatures, AGGREGATEDINDELS)

    # For Vectorization
    # These are used in writing tables
    allSims_subsSignature_accumulated_signal_np_array = np.zeros((numofSimulations+1, number_of_sbs_signatures + 1, plus_minus_epigenomics * 2 + 1))
    allSims_dinucsSignature_accumulated_signal_np_array = np.zeros((numofSimulations+1, number_of_dbs_signatures + 1, plus_minus_epigenomics * 2 + 1))
    allSims_indelsSignature_accumulated_signal_np_array = np.zeros((numofSimulations+1, number_of_id_signatures + 1, plus_minus_epigenomics * 2 + 1))

    allSims_subsSignature_accumulated_count_np_array = np.zeros((numofSimulations+1, number_of_sbs_signatures + 1, plus_minus_epigenomics * 2 + 1))
    allSims_dinucsSignature_accumulated_count_np_array = np.zeros((numofSimulations+1, number_of_dbs_signatures + 1, plus_minus_epigenomics * 2 + 1))
    allSims_indelsSignature_accumulated_count_np_array = np.zeros((numofSimulations+1, number_of_id_signatures+1, plus_minus_epigenomics * 2 + 1))

    # This code reads the file and prepare chrbased signal files
    # If file is in default files,chr based signal files are downloaded from ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerTopography/lib/
    # No need for preparing here
    library_file_type = None
    chrBased_library_df = None
    library_grouped_df = None

    if (os.path.basename(library_file_with_path) not in SIGPROFILERTOPOGRAPHY_DEFAULT_FILES):
        # What is the type of the signal_file_with_path?
        # If it is a bed file read signal_file_with_path here
        file_extension = os.path.splitext(os.path.basename(library_file_with_path))[1]

        if ((file_extension.lower() == '.bigwig') or (file_extension.lower() == '.bw')):
            library_file_type = BIGWIG
            # if chrBasedSignalArrays does not exist we will use pyBigWig if installed and we will not create chrBasedSignalArrays but use BigWig file opened by pyBigWig to fill windowArray
        elif ((file_extension.lower() == '.bigbed') or (file_extension.lower() == '.bb')):
            library_file_type = BIGBED
            # if chrBasedSignalArrays does not exist we will use pyBigWig if installed and we will not create chrBasedSignalArrays but use BigBed file opened by pyBigWig to fill windowArray
        elif (file_extension.lower() == '.bed'):
            library_file_type = BED

            if os.path.exists(library_file_with_path):
                bedfilename, bedfile_extention = os.path.splitext(os.path.basename(library_file_with_path))

                if (bedfile_extention.lower() == '.bed' or bedfile_extention.lower() == '.np' or bedfile_extention.lower() == '.narrowpeak'):
                    discard_signal = False
                    library_file_df, max_signal, min_signal = readFileInBEDFormat(library_file_with_path, discard_signal, log_file)

                    # Eliminate Outliers
                    if ((remove_outliers == True) and (quantile_value < 1.0)):
                        # remove the outliers
                        q = library_file_df[SIGNAL].quantile(quantile_value)
                        log_out = open(log_file, 'a')
                        print('Signal greater than %f is removed for quantile: %f' % (q, quantile_value), file=log_out)
                        print('before outlier removal number of rows: %d' % (library_file_df.shape[0]), file=log_out)
                        library_file_df = library_file_df[library_file_df[SIGNAL] < q]
                        print('after outlier removal number of rows: %d' % (library_file_df.shape[0]), file=log_out)
                        log_out.close()

                    library_grouped_df = library_file_df.groupby(CHROM)

        elif ((file_extension.lower() == '.narrowpeak') or (file_extension.lower() == '.np')):
            library_file_type = NARROWPEAK
            readBEDandWriteChromBasedSignalArrays(outputDir, jobname, genome, library_file_with_path, occupancy_type, quantile_value, remove_outliers, log_file)
        elif (file_extension.lower() == '.wig'):
            library_file_type = WIG
            # For inhouse preparation
            # readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArraysSequentially(genome, quantileValue,library_file_with_path)
            # readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArraysInParallel(genome, quantileValue,library_file_with_path)
            filetype_BEDGRAPH = isFileTypeBedGraph(library_file_with_path)
            if filetype_BEDGRAPH:
                if verbose: start_time = time.time()
                # Read by chunks
                # readWig_write_derived_from_bedgraph_using_pool_chunks(outputDir, jobname, genome, library_file_with_path,occupancy_type,remove_outliers,verbose)
                # Read at once
                readWig_write_derived_from_bedgraph_using_pool_read_all(outputDir, jobname, genome, library_file_with_path, occupancy_type,remove_outliers,verbose,quantile_value, log_file)
                if verbose:
                    log_out = open(log_file, 'a')
                    print('\tVerbose Read wig file and write chrbased arrays took %f seconds' %((time.time() - start_time)), file=log_out)
                    log_out.close()

                # For 6 GB ATAC-seq file using pool took 8 min whereas without pool it took 16 min.
                # start_time = time.time()
                # readWig_write_derived_from_bedgraph(outputDir, jobname, genome, library_file_with_path,occupancy_type,verbose)
                # print('Without pool Took %f seconds' %((time.time() - start_time)))
            else:
                readWig_with_fixedStep_variableStep_writeChrBasedSignalArrays(outputDir, jobname, genome, library_file_with_path, occupancy_type, quantile_value, remove_outliers, log_file)
        elif (file_extension.lower() == '.bedgraph'):
            library_file_type = BEDGRAPH
            readWig_write_derived_from_bedgraph_using_pool_read_all(outputDir, jobname, genome, library_file_with_path,occupancy_type, remove_outliers, verbose, quantile_value, log_file)
        else:
            library_file_type = LIBRARY_FILE_TYPE_OTHER


    def accumulate_apply_async_result_vectorization(simulatonBased_SignalArrayAndCountArrayList):
        try:
            if simulatonBased_SignalArrayAndCountArrayList is not None:
                chrLong = simulatonBased_SignalArrayAndCountArrayList[0]
                simNum = simulatonBased_SignalArrayAndCountArrayList[1]
                subsSignature_accumulated_signal_np_array = simulatonBased_SignalArrayAndCountArrayList[2]
                dinucsSignature_accumulated_signal_np_array = simulatonBased_SignalArrayAndCountArrayList[3]
                indelsSignature_accumulated_signal_np_array = simulatonBased_SignalArrayAndCountArrayList[4]
                subsSignature_accumulated_count_np_array = simulatonBased_SignalArrayAndCountArrayList[5]
                dinucsSignature_accumulated_count_np_array = simulatonBased_SignalArrayAndCountArrayList[6]
                indelsSignature_accumulated_count_np_array = simulatonBased_SignalArrayAndCountArrayList[7]

                # Accumulation
                allSims_subsSignature_accumulated_signal_np_array[simNum] += subsSignature_accumulated_signal_np_array
                allSims_dinucsSignature_accumulated_signal_np_array[simNum] += dinucsSignature_accumulated_signal_np_array
                allSims_indelsSignature_accumulated_signal_np_array[simNum] += indelsSignature_accumulated_signal_np_array

                allSims_subsSignature_accumulated_count_np_array[simNum] += subsSignature_accumulated_count_np_array
                allSims_dinucsSignature_accumulated_count_np_array[simNum] += dinucsSignature_accumulated_count_np_array
                allSims_indelsSignature_accumulated_count_np_array[simNum] += indelsSignature_accumulated_count_np_array
                # print('ACCUMULATION chrLong:%s simNum:%d ENDS' %(chrLong,simNum))

        except Exception as e:
            print("Exception in accumulate_apply_async_result_vectorization function: %s" % (e))

    sim_nums = range(0, numofSimulations + 1)
    sim_num_chr_tuples = ((sim_num, chrLong) for sim_num in sim_nums for chrLong in chromNamesList)

    if parallel_mode:
        jobs = []

        if (computation_type == USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM):
            numofProcesses = multiprocessing.cpu_count()
            pool = multiprocessing.Pool(processes=numofProcesses)

            for simNum, chrLong in sim_num_chr_tuples:
                jobs.append(pool.apply_async(chrbased_data_fill_signal_count_arrays_for_all_mutations_read_mutations_using_pyranges,
                                             args=(occupancy_type,
                                                   occupancy_calculation_type,
                                                   outputDir,
                                                   jobname,
                                                   chrLong,
                                                   simNum,
                                                   samples_of_interest,
                                                   chromSizesDict,
                                                   library_grouped_df,
                                                   chrBased_library_df,
                                                   library_file_with_path,
                                                   library_file_type,
                                                   ordered_sbs_signatures,
                                                   ordered_dbs_signatures,
                                                   ordered_id_signatures,
                                                   ordered_sbs_signatures_cutoffs,
                                                   ordered_dbs_signatures_cutoffs,
                                                   ordered_id_signatures_cutoffs,
                                                   plus_minus_epigenomics,
                                                   discreet_mode,
                                                   default_cutoff,
                                                   log_file,
                                                   verbose,),
                                             callback=accumulate_apply_async_result_vectorization))

            pool.close()
            pool.join()

    else:
        # Sequential mode for profiling, debugging and testing purposes
        for simNum, chrLong in sim_num_chr_tuples:
            simulatonBased_SignalArrayAndCountArrayList = chrbased_data_fill_signal_count_arrays_for_all_mutations_read_mutations_using_pyranges(occupancy_type,
                                                   occupancy_calculation_type,
                                                   outputDir,
                                                   jobname,
                                                   chrLong,
                                                   simNum,
                                                   samples_of_interest,
                                                   chromSizesDict,
                                                   library_grouped_df,
                                                   chrBased_library_df,
                                                   library_file_with_path,
                                                   library_file_type,
                                                   ordered_sbs_signatures,
                                                   ordered_dbs_signatures,
                                                   ordered_id_signatures,
                                                   ordered_sbs_signatures_cutoffs,
                                                   ordered_dbs_signatures_cutoffs,
                                                   ordered_id_signatures_cutoffs,
                                                   plus_minus_epigenomics,
                                                   discreet_mode,
                                                   default_cutoff,
                                                   log_file,
                                                   verbose)

            accumulate_apply_async_result_vectorization(simulatonBased_SignalArrayAndCountArrayList)


    # Same for parallel or sequential run
    writeSimulationBasedAverageOccupancyUsingNumpyArray(occupancy_type,
                                                   sample_based,
                                                   plus_minus_epigenomics,
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

    log_out = open(log_file, 'a')
    print('--- %s Analysis ends' %(occupancy_type), file=log_out)
    print('#################################################################################\n', file=log_out)
    log_out.close()


# main function
# Using pyBigWig for bigBed and bigWig files starts Optional for unix, linux
# Using chrBasedSignalArrays for big files
# Using dataframes for small bed files
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
                        samples_of_interest,
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
                        discreet_mode,
                        default_cutoff,
                        parallel_mode,
                        log_file,
                        verbose):

    log_out = open(log_file, 'a')
    if (os.path.basename(library_file_with_path) not in SIGPROFILERTOPOGRAPHY_DEFAULT_FILES) and (not os.path.exists(library_file_with_path)):
        print('There is no such file under %s' %(library_file_with_path), file=log_out)

    print('\n#################################################################################', file=log_out)
    print('--- %s Analysis starts' %(occupancy_type), file=log_out)
    print('--- Computation Type:%s' % (computation_type), file=log_out)
    print('--- Occupancy Type:%s' % (occupancy_type), file=log_out)
    print('--- Library file with path: %s\n' %library_file_with_path, file=log_out)
    log_out.close()

    # Using pyBigWig for BigWig and BigBed files if you can import pyBigWig (linux only) otherwise no
    # By the way pyBigWig can be imported in unix, linux like os not available in windows
    # Using HM and CTCF bed files preparing chr based signal array runtime
    # Using ATAC-seq wig files preparing chr based signal array runtime

    ##########################################################################
    # If chunksize is 1, maxtasksperchild=x will call the function x times in each process,
    # but if chunksize is y, it will call the function x*y times in each process.
    # Setting maxtaskperchild to 1 would restart each process in your pool after it processed a single task, which is the most aggressive setting you could use to free any leaked resources.
    # numofProcesses = multiprocessing.cpu_count()
    # pool = multiprocessing.Pool(numofProcesses, maxtasksperchild=1)
    ##########################################################################

    number_of_sbs_signatures = ordered_sbs_signatures.size
    number_of_dbs_signatures = ordered_dbs_signatures.size
    number_of_id_signatures = ordered_id_signatures.size

    subsSignatures = np.append(ordered_sbs_signatures, AGGREGATEDSUBSTITUTIONS)
    dinucsSignatures = np.append(ordered_dbs_signatures, AGGREGATEDDINUCS)
    indelsSignatures = np.append(ordered_id_signatures, AGGREGATEDINDELS)

    # For Vectorization
    # These are used in writing tables
    allSims_subsSignature_accumulated_signal_np_array = np.zeros((numofSimulations+1, number_of_sbs_signatures + 1, plusorMinus * 2 + 1))
    allSims_dinucsSignature_accumulated_signal_np_array = np.zeros((numofSimulations+1, number_of_dbs_signatures + 1, plusorMinus * 2 + 1))
    allSims_indelsSignature_accumulated_signal_np_array = np.zeros((numofSimulations+1, number_of_id_signatures + 1, plusorMinus * 2 + 1))

    allSims_subsSignature_accumulated_count_np_array = np.zeros((numofSimulations+1, number_of_sbs_signatures + 1, plusorMinus * 2 + 1))
    allSims_dinucsSignature_accumulated_count_np_array = np.zeros((numofSimulations+1, number_of_dbs_signatures + 1, plusorMinus * 2 + 1))
    allSims_indelsSignature_accumulated_count_np_array = np.zeros((numofSimulations+1, number_of_id_signatures+1, plusorMinus * 2 + 1))

    # This code reads the file and prepare chrbased signal files
    # If file is in default files,chr based signal files are downloaded from ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerTopography/lib/
    # No need for preparing here
    library_file_type = None

    if (os.path.basename(library_file_with_path) not in SIGPROFILERTOPOGRAPHY_DEFAULT_FILES):
        # What is the type of the signal_file_with_path?
        # If it is a bed file read signal_file_with_path here
        file_extension = os.path.splitext(os.path.basename(library_file_with_path))[1]

        if ((file_extension.lower() == '.bigwig') or (file_extension.lower() == '.bw')):
            library_file_type = BIGWIG
            # if chrBasedSignalArrays does not exist we will use pyBigWig if installed and we will not create chrBasedSignalArrays but use BigWig file opened by pyBigWig to fill windowArray
        elif ((file_extension.lower() == '.bigbed') or (file_extension.lower() == '.bb')):
            library_file_type = BIGBED
            # if chrBasedSignalArrays does not exist we will use pyBigWig if installed and we will not create chrBasedSignalArrays but use BigBed file opened by pyBigWig to fill windowArray
        elif (file_extension.lower() == '.bed'):
            library_file_type = BED
            readBEDandWriteChromBasedSignalArrays(outputDir, jobname, genome, library_file_with_path, occupancy_type, quantileValue, remove_outliers, log_file)
        elif ((file_extension.lower() == '.narrowpeak') or (file_extension.lower() == '.np')):
            library_file_type = NARROWPEAK
            readBEDandWriteChromBasedSignalArrays(outputDir, jobname, genome, library_file_with_path, occupancy_type, quantileValue, remove_outliers, log_file)
        elif (file_extension.lower() == '.wig'):
            library_file_type = WIG
            # For inhouse preparation
            # readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArraysSequentially(genome, quantileValue,library_file_with_path)
            # readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArraysInParallel(genome, quantileValue,library_file_with_path)
            filetype_BEDGRAPH = isFileTypeBedGraph(library_file_with_path)
            if filetype_BEDGRAPH:
                if verbose: start_time = time.time()
                # Read by chunks
                # readWig_write_derived_from_bedgraph_using_pool_chunks(outputDir, jobname, genome, library_file_with_path,occupancy_type,remove_outliers,verbose)
                # Read at once
                readWig_write_derived_from_bedgraph_using_pool_read_all(outputDir, jobname, genome, library_file_with_path, occupancy_type,remove_outliers,verbose,quantileValue, log_file)
                if verbose:
                    log_out = open(log_file, 'a')
                    print('\tVerbose Read wig file and write chrbased arrays took %f seconds' %((time.time() - start_time)), file=log_out)
                    log_out.close()

                # For 6 GB ATAC-seq file using pool took 8 min whereas without pool it took 16 min.
                # start_time = time.time()
                # readWig_write_derived_from_bedgraph(outputDir, jobname, genome, library_file_with_path,occupancy_type,verbose)
                # print('Without pool Took %f seconds' %((time.time() - start_time)))
            else:
                readWig_with_fixedStep_variableStep_writeChrBasedSignalArrays(outputDir, jobname, genome, library_file_with_path, occupancy_type, quantileValue, remove_outliers, log_file)
        elif (file_extension.lower() == '.bedgraph'):
            library_file_type = BEDGRAPH
            readWig_write_derived_from_bedgraph_using_pool_read_all(outputDir, jobname, genome, library_file_with_path,occupancy_type, remove_outliers, verbose, quantileValue, log_file)
        else:
            library_file_type = LIBRARY_FILE_TYPE_OTHER


    def accumulate_apply_async_result_vectorization(simulatonBased_SignalArrayAndCountArrayList):
        try:
            chrLong = simulatonBased_SignalArrayAndCountArrayList[0]
            simNum = simulatonBased_SignalArrayAndCountArrayList[1]
            subsSignature_accumulated_signal_np_array = simulatonBased_SignalArrayAndCountArrayList[2]
            dinucsSignature_accumulated_signal_np_array = simulatonBased_SignalArrayAndCountArrayList[3]
            indelsSignature_accumulated_signal_np_array = simulatonBased_SignalArrayAndCountArrayList[4]
            subsSignature_accumulated_count_np_array = simulatonBased_SignalArrayAndCountArrayList[5]
            dinucsSignature_accumulated_count_np_array = simulatonBased_SignalArrayAndCountArrayList[6]
            indelsSignature_accumulated_count_np_array = simulatonBased_SignalArrayAndCountArrayList[7]

            # Accumulation
            allSims_subsSignature_accumulated_signal_np_array[simNum] += subsSignature_accumulated_signal_np_array
            allSims_dinucsSignature_accumulated_signal_np_array[simNum] += dinucsSignature_accumulated_signal_np_array
            allSims_indelsSignature_accumulated_signal_np_array[simNum] += indelsSignature_accumulated_signal_np_array

            allSims_subsSignature_accumulated_count_np_array[simNum] += subsSignature_accumulated_count_np_array
            allSims_dinucsSignature_accumulated_count_np_array[simNum] += dinucsSignature_accumulated_count_np_array
            allSims_indelsSignature_accumulated_count_np_array[simNum] += indelsSignature_accumulated_count_np_array
            # print('ACCUMULATION chrLong:%s simNum:%d ENDS' %(chrLong,simNum))

        except Exception as e:
            print("Exception in accumulate_apply_async_result_vectorization function: %s" % (e))

    sim_nums = range(0, numofSimulations + 1)
    sim_num_chr_tuples = ((sim_num, chrLong) for sim_num in sim_nums for chrLong in chromNamesList)

    if parallel_mode:
        jobs = []

        if (computation_type == USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM):
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
                                                   samples_of_interest,
                                                   chromSizesDict,
                                                   library_file_with_path,
                                                   library_file_type,
                                                   ordered_sbs_signatures,
                                                   ordered_dbs_signatures,
                                                   ordered_id_signatures,
                                                   ordered_sbs_signatures_cutoffs,
                                                   ordered_dbs_signatures_cutoffs,
                                                   ordered_id_signatures_cutoffs,
                                                   remove_outliers,
                                                   plusorMinus,
                                                   discreet_mode,
                                                   default_cutoff,
                                                   log_file,
                                                   verbose,),
                                             callback=accumulate_apply_async_result_vectorization))

            pool.close()
            pool.join()

        elif (computation_type == USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT):
            numofProcesses = multiprocessing.cpu_count()
            pool = multiprocessing.Pool(processes=numofProcesses)

            for chrLong, simNum, splitIndex in job_tuples:
                jobs.append(
                    pool.apply_async(chrbased_data_fill_signal_count_arrays_for_all_mutations_read_mutations_split,
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
                                           discreet_mode,
                                           default_cutoff,
                                           log_file,
                                           verbose,),
                                     callback=accumulate_apply_async_result_vectorization))
            pool.close()
            pool.join()

    else:
        # Sequential mode for profiling, debugging and testing purposes
        for simNum, chrLong in sim_num_chr_tuples:
            simulatonBased_SignalArrayAndCountArrayList = chrbased_data_fill_signal_count_arrays_for_all_mutations_read_mutations(occupancy_type,
                                   occupancy_calculation_type,
                                   outputDir,
                                   jobname,
                                   chrLong,
                                   simNum,
                                   samples_of_interest,
                                   chromSizesDict,
                                   library_file_with_path,
                                   library_file_type,
                                   ordered_sbs_signatures,
                                   ordered_dbs_signatures,
                                   ordered_id_signatures,
                                   ordered_sbs_signatures_cutoffs,
                                   ordered_dbs_signatures_cutoffs,
                                   ordered_id_signatures_cutoffs,
                                   remove_outliers,
                                   plusorMinus,
                                   discreet_mode,
                                   default_cutoff,
                                   log_file,
                                   verbose)

            accumulate_apply_async_result_vectorization(simulatonBased_SignalArrayAndCountArrayList)

    # Same for parallel or sequential run
    writeSimulationBasedAverageOccupancyUsingNumpyArray(occupancy_type,
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

    log_out = open(log_file, 'a')
    print('--- %s Analysis ends' %(occupancy_type), file=log_out)
    print('#################################################################################\n', file=log_out)
    log_out.close()