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


# Version2
# This version use np.arrays
# Right now replication strand bias analysis works for single point mutations and signatures.
# This python code analyses the Replication Strand Bias

import multiprocessing
import numpy as np
import pandas as pd
import os

from functools import reduce
from scipy.signal import find_peaks

from SigProfilerTopography.source.commons.TopographyCommons import CHROM
from SigProfilerTopography.source.commons.TopographyCommons import START
from SigProfilerTopography.source.commons.TopographyCommons import END
from SigProfilerTopography.source.commons.TopographyCommons import SIGNAL

from SigProfilerTopography.source.commons.TopographyCommons import PYRAMIDINESTRAND
from SigProfilerTopography.source.commons.TopographyCommons import SAMPLE

from SigProfilerTopography.source.commons.TopographyCommons import TYPE

from SigProfilerTopography.source.commons.TopographyCommons import SUBS
from SigProfilerTopography.source.commons.TopographyCommons import INDELS
from SigProfilerTopography.source.commons.TopographyCommons import DINUCS

from SigProfilerTopography.source.commons.TopographyCommons import MUTATION
from SigProfilerTopography.source.commons.TopographyCommons import MUTATIONLONG
from SigProfilerTopography.source.commons.TopographyCommons import LENGTH

from SigProfilerTopography.source.commons.TopographyCommons import LEADING
from SigProfilerTopography.source.commons.TopographyCommons import LAGGING
from SigProfilerTopography.source.commons.TopographyCommons import REPLICATIONSTRANDBIAS

from SigProfilerTopography.source.commons.TopographyCommons import DATA
from SigProfilerTopography.source.commons.TopographyCommons import LIB
from SigProfilerTopography.source.commons.TopographyCommons import CHRBASED

from SigProfilerTopography.source.commons.TopographyCommons import readWig_with_fixedStep_variableStep

from SigProfilerTopography.source.occupancy.ChrBasedSignalArrays import readFileInBEDFormat
from SigProfilerTopography.source.commons.TopographyCommons import memory_usage

from SigProfilerTopography.source.commons.TopographyCommons import write_type_strand_bias_np_array_as_dataframe
from SigProfilerTopography.source.commons.TopographyCommons import write_signature_mutation_type_strand_bias_np_array_as_dataframe
from SigProfilerTopography.source.commons.TopographyCommons import write_sbs_signature_sbs96_mutation_type_replication_strand_bias
from SigProfilerTopography.source.commons.TopographyCommons import write_sample_based_strand1_strand2_as_dataframe

from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM
from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT

from SigProfilerTopography.source.commons.TopographyCommons import get_chrBased_simBased_combined_df_split
from SigProfilerTopography.source.commons.TopographyCommons import get_chrBased_simBased_combined_df
from SigProfilerTopography.source.commons.TopographyCommons import get_chrBased_simBased_dfs

from SigProfilerTopography.source.commons.TopographyCommons import isFileTypeBedGraph

from SigProfilerTopography.source.commons.TopographyCommons import C2A
from SigProfilerTopography.source.commons.TopographyCommons import C2G
from SigProfilerTopography.source.commons.TopographyCommons import C2T
from SigProfilerTopography.source.commons.TopographyCommons import T2A
from SigProfilerTopography.source.commons.TopographyCommons import T2C
from SigProfilerTopography.source.commons.TopographyCommons import T2G

# For Supp Fig2B
CHR10_THRESHOLD_START = 16400000
CHR10_THRESHOLD_END = 26400000

# For Supp Fig2A
CHR20_START = 36260000
CHR20_END = 36830000

# FOR FINDING TRANSITION ZONES (LEADING or LAGGING)
# THRESHOLD_CONSECUTIVE_LONG_STRETCH_LENGTH= 250000 #used in Supp Fig2B
# THRESHOLD_CONSECUTIVE_LONG_STRETCH_LENGTH= 150000
THRESHOLD_CONSECUTIVE_LONG_STRETCH_LENGTH= 10000

# THRESHOLD_DISCARD_LATEST_TRANSITION_ZONE = 100000 #used in Supp Fig2B
THRESHOLD_DISCARD_LATEST_TRANSITION_ZONE = 25000
# THRESHOLD_LATEST_TRANSITION_ZONE = 0

def check_for_consecutive_regions_with_same_signed_slope(chrLong,
                                                         start, # can be repli_seq_start or peak_midpoint or valley_midpoint
                                                         end, # can be repli_seq_end or peak_midpoint or valley_midpoint
                                                         chrBased_repli_seq_signal_df,
                                                         info='No info'):
    transition_zone_list = []

    # Get the repli_seq signals within the given start and end
    sub_repli_seq_df = chrBased_repli_seq_signal_df[( (chrBased_repli_seq_signal_df[START] + chrBased_repli_seq_signal_df[END])//2 >= start) &
                                             ((chrBased_repli_seq_signal_df[START] + chrBased_repli_seq_signal_df[END])//2 <= end)]

    sub_repli_seq_df['signal_diff'] = sub_repli_seq_df['Signal'].diff() # signal[i+1] - signal[i]
    signal_diff_arr = sub_repli_seq_df['signal_diff'].values

    # Are these all nans?
    all_nans = np.all(np.isnan(signal_diff_arr))

    # Remove nans
    signal_diff_arr = signal_diff_arr[~np.isnan(signal_diff_arr)]

    if all_nans:
        direction_sign = 0
        if info == 'up_to_peak':
            direction_sign = +1
        elif info == 'down_to_valley':
            direction_sign = -1
        # For the last one if there is no other sign it will be 0
    elif np.all(signal_diff_arr>0):
        direction_sign = +1
    elif np.all(signal_diff_arr<0):
        direction_sign = -1

    transition_zone_list.append((chrLong, start, end, direction_sign, end-start))

    return transition_zone_list


# # deprecated
# def check_for_same_signed_slope_between_consecutive_peak_and_valley(chrLong,
#                                                                     peak_or_valley_start,
#                                                                     peak_or_valley_end,
#                                                                     chrBased_repli_seq_signal_df):
#
#     transition_zone_list = []
#
#     # Get the repli-seq signals within peak_or_valley_start and peak_or_valley_end
#     subset_df = chrBased_repli_seq_signal_df[(chrBased_repli_seq_signal_df[START] >= peak_or_valley_start) &
#                                              (chrBased_repli_seq_signal_df[END] <= peak_or_valley_end)]
#
#     consecutive_length = 0
#     former_row = None
#     former_slope_direction = None
#     start = peak_or_valley_start
#
#     for index, row in subset_df.iterrows():
#         region_length = row.get(END) - row.get(START)
#
#         if former_row is None:
#             # We read the row for the first time
#             former_row = row
#             consecutive_length += region_length
#         else:
#             slope = (row.get(SIGNAL) - former_row.get(SIGNAL)) / region_length
#             former_row = row
#
#             if (former_slope_direction is None):
#                 former_slope_direction = np.sign(slope)
#                 consecutive_length += region_length
#             elif (former_slope_direction == np.sign(slope)):
#                 consecutive_length += region_length
#             else:
#                 # They have different signs, there is a sign transition
#                 if (consecutive_length >= THRESHOLD_CONSECUTIVE_LONG_STRETCH_LENGTH):
#                     transition_zone_list.append((chrLong, start, (row.get(START) + row.get(END))//2, former_slope_direction, consecutive_length))
#
#                 # initialize and start again
#                 consecutive_length = region_length
#                 start = (row.get(START) + row.get(END)) // 2
#                 former_row = row
#                 former_slope_direction = np.sign(slope)
#                 continue
#
#             former_slope_direction = np.sign(slope)
#
#     # This is for the last probable transition zone.
#     if (consecutive_length >= THRESHOLD_CONSECUTIVE_LONG_STRETCH_LENGTH):
#         transition_zone_list.append((chrLong, start, (row.get(START) + row.get(END)) // 2, former_slope_direction, consecutive_length))
#
#     return transition_zone_list



# chr10_subset_wavelet_processed_df
#           chr     start       end   signal
# 265577  chr10  16400500  16401499  24.9438
# 265578  chr10  16401500  16402499  24.9585

# valleys_peaks_df
#         chr     start       end    type
# 415     chr10  16454500  16455500    Peak
# 415  chr10  16528500  16529500  Valley

def find_long_stretches_of_consistent_transition_zones(chrLong,
                                                       start,
                                                       end,
                                                       chrBased_repli_seq_signal_df,
                                                       valleys_peaks_df):

    transition_zones_list = []

    for index, row in valleys_peaks_df.iterrows():
        if (row['type'] == 'Peak'):
            peak_start = row[START]
            peak_end = row[END]
            peak_midpoint = (peak_start + peak_end) // 2

            if (peak_midpoint > start):
                found_zone = check_for_consecutive_regions_with_same_signed_slope(chrLong,
                                                                     start,
                                                                     peak_midpoint,
                                                                     chrBased_repli_seq_signal_df,
                                                                     info='up_to_peak')

                transition_zones_list.extend(found_zone)

            start = peak_midpoint

        elif (row['type'] == 'Valley'):
            valley_start = row[START]
            valley_end = row[END]
            valley_midpoint = (valley_start + valley_end) // 2

            # New Way
            if valley_midpoint > start:
                found_zone = check_for_consecutive_regions_with_same_signed_slope(chrLong,
                                                                     start,
                                                                     valley_midpoint,
                                                                     chrBased_repli_seq_signal_df,
                                                                     info='down_to_valley')

                transition_zones_list.extend(found_zone)
            start = valley_midpoint

            # # Old Way
            # # This is something special to valley
            # # To remain conservative in downstream assignments,
            # # we removed the last 25kb of the latest zones of the replicating domains
            # new_valley_start1 = valley_midpoint - THRESHOLD_DISCARD_LATEST_TRANSITION_ZONE
            # new_valley_start2 = valley_midpoint + THRESHOLD_DISCARD_LATEST_TRANSITION_ZONE
            #
            # if (new_valley_start1 > start):
            #     found_zone = check_for_consecutive_regions_with_same_signed_slope(chrLong,
            #                                                          start,
            #                                                          new_valley_start1,
            #                                                          chrBased_repli_seq_signal_df,
            #                                                          file2,
            #                                                          info='down_to_valley')
            #
            #     transition_zones_list.extend(found_zone)
            #
            # # bypass the genome region between newValleyStart1 and newValleyStart2
            # start = new_valley_start2

    # For the last interval
    if (end > start):
        found_zone = check_for_consecutive_regions_with_same_signed_slope(chrLong,
                                                             start,
                                                             end,
                                                             chrBased_repli_seq_signal_df,
                                                             info='from start to end')

        transition_zones_list.extend(found_zone)

    return transition_zones_list


# We assume that there are no overlapping intervals with positive and negative slopes.
# To test it have one array for positive slope fill with 1s
#                one array for negative slope fill with -2a
#                add them if you habe any -1 that means that you contradict this assumption.
def fill_replication_strand_array(replication_strand_row,chrBased_replication_array):
    # e.g.: replicationStrand_row
    # chr chrX
    # start   154861998
    # end 155096999
    # slope_direction  1 (1 means leading strand -1 means lagging strand on positive strand)
    # length  235000

    # labels = ['chr', 'start', 'end', 'slope_direction', 'length']
    chrBased_replication_array[replication_strand_row['start']:replication_strand_row['end']] = replication_strand_row['slope_direction']

# Using numpy arrays
#   if mutationPyramidineStrand and slope have the same sign increase LEADING STRAND count
#   else mutationPyramidineStrand and slope have the opposite sign increase LAGGING STRAND count
# sample_based for further usage
def searchAllMutationOnReplicationStrandArray_using_list_comprehension_using_numpy_array(
        mutation_row,
        my_type,
        chrBasedReplicationArray,
        SBS6_mutation_types_np_array,
        SBS96_mutation_types_np_array,
        ordered_signatures_cutoffs,
        df_columns_signatures_mask_array,
        SBS6_mutation_types_default_zeros_array,
        SBS96_mutation_types_default_zeros_array,
        subs_signatures_default_zeros_array,
        dinucs_signatures_default_zeros_array,
        indels_signatures_default_zeros_array,
        subs_signatures_SBS6_mutation_types_default_zeros_array,
        subs_signatures_SBS96_mutation_types_default_zeros_array,
        all_types_leading_np_array,
        all_types_lagging_np_array,
        subs_signature_SBS6_mutation_type_leading_np_array,
        subs_signature_SBS6_mutation_type_lagging_np_array,
        subs_signature_SBS96_mutation_type_leading_np_array,
        subs_signature_SBS96_mutation_type_lagging_np_array,
        all_samples_all_types_leading_np_array,
        all_samples_all_types_lagging_np_array,
        all_samples_subs_signature_mutation_type_leading_np_array,
        all_samples_subs_signature_mutation_type_lagging_np_array,
        sample_based,
        all_samples_np_array,
        discreet_mode,
        default_cutoff,
        df_columns):

    if sample_based:
        indexofSample = np.where(df_columns == SAMPLE)[0][0]
        mutation_sample = mutation_row[indexofSample]
        sample_index = np.where(all_samples_np_array == mutation_sample)[0][0]

    indexofStart = np.where(df_columns == START)[0][0]
    start = mutation_row[indexofStart]

    indexofPyrimidineStrand = np.where(df_columns == PYRAMIDINESTRAND)[0][0]
    pyramidineStrand = mutation_row[indexofPyrimidineStrand]

    subs_signature_SBS6_mutation_type_mask_array = subs_signatures_SBS6_mutation_types_default_zeros_array
    subs_signature_SBS96_mutation_type_mask_array = subs_signatures_SBS96_mutation_types_default_zeros_array

    probabilities = mutation_row[df_columns_signatures_mask_array]

    if (my_type == SUBS):
        end = start+1
        # e.g.: C>A

        indexofMutation = np.where(df_columns == MUTATION)[0][0]
        SBS6_mutation_type = mutation_row[indexofMutation]

        index_of_mutation_long = np.where(df_columns == MUTATIONLONG)[0][0]

        # e.g.: T:AA[C>A]AA
        mutation_type_long = mutation_row[index_of_mutation_long]
        SBS96_mutation_type = mutation_type_long[3:10]

        # six_mutation_types_mask_array.shape (6,)
        SBS6_mutation_types_mask_array = np.where(SBS6_mutation_types_np_array == SBS6_mutation_type, 1, 0)
        SBS96_mutation_types_mask_array = np.where(SBS96_mutation_types_np_array == SBS96_mutation_type, 1, 0)

        if discreet_mode:
            # Convert True into 1, and False into 0
            # subs_signatures_mask_array.shape (num_of_subs_signatures,)
            threshold_mask_array = np.greater_equal(probabilities, ordered_signatures_cutoffs)
            subs_signatures_mask_array = threshold_mask_array.astype(int)
        else:
            # new way
            probabilities[probabilities < default_cutoff ] = 0
            subs_signatures_mask_array = np.array(probabilities).astype(float)
            # old way:
            # subs_signatures_mask_array = np.array(probabilities).astype(float)

        # Concetanate
        all_types_mask_array = np.concatenate((SBS6_mutation_types_mask_array,
                                              SBS96_mutation_types_mask_array,
                                              subs_signatures_mask_array, # SUBS
                                              dinucs_signatures_default_zeros_array,
                                              indels_signatures_default_zeros_array), axis=None)

        # multiply subs_signatures_mask_array times six_mutation_types_mask_array
        # Add one more dimension to subs_signatures_mask_array and six_mutation_types_mask_array

        # subs_signatures_mask_array_2d.shape (1,num_of_subs_signatures)
        subs_signatures_mask_array_2d = np.expand_dims(subs_signatures_mask_array, axis=0)

        # six_mutation_types_mask_array_2d.shape (1,6)
        SBS6_mutation_types_mask_array_2d = np.expand_dims(SBS6_mutation_types_mask_array, axis=0)

        # SBS96_mutation_types_mask_array.shape (96,) --> (1,96)
        SBS96_mutation_types_mask_array_2d = np.expand_dims(SBS96_mutation_types_mask_array, axis=0)

        # to_be_accumulated_array.shape (num_of_subs_signatures,6)
        subs_signature_SBS6_mutation_type_mask_array = subs_signatures_mask_array_2d.T * SBS6_mutation_types_mask_array_2d

        # to_be_accumulated_array.shape (num_of_subs_signatures,96)
        subs_signature_SBS96_mutation_type_mask_array = subs_signatures_mask_array_2d.T * SBS96_mutation_types_mask_array_2d

    elif (my_type == DINUCS):
        end = start+2

        if discreet_mode:
            # Convert True into 1, and False into 0
            threshold_mask_array = np.greater_equal(probabilities, ordered_signatures_cutoffs)
            dinucs_signatures_mask_array = threshold_mask_array.astype(int)
        else:
            # new way
            probabilities[probabilities < default_cutoff ] = 0
            dinucs_signatures_mask_array = np.array(probabilities).astype(float)
            # old way
            # dinucs_signatures_mask_array = np.array(probabilities).astype(float)

        # Concetanate
        all_types_mask_array= np.concatenate((SBS6_mutation_types_default_zeros_array,
                                              SBS96_mutation_types_default_zeros_array,
                                              subs_signatures_default_zeros_array,
                                              dinucs_signatures_mask_array, # DINUCS
                                              indels_signatures_default_zeros_array), axis=None)

    elif (my_type == INDELS):
        indexofLength = np.where(df_columns == LENGTH)[0][0]
        end = start + int(mutation_row[indexofLength])

        if discreet_mode:
            # Convert True into 1, and False into 0
            threshold_mask_array = np.greater_equal(probabilities, ordered_signatures_cutoffs)
            indels_signatures_mask_array = threshold_mask_array.astype(int)
        else:
            probabilities[probabilities < default_cutoff ] = 0
            indels_signatures_mask_array = np.array(probabilities).astype(float)
            # old way:
            # indels_signatures_mask_array = np.array(probabilities).astype(float)

        # Concetanate
        all_types_mask_array= np.concatenate((SBS6_mutation_types_default_zeros_array,
                                              SBS96_mutation_types_default_zeros_array,
                                              subs_signatures_default_zeros_array,
                                              dinucs_signatures_default_zeros_array,
                                              indels_signatures_mask_array), axis=None)


    # if there is overlap with chrBasedReplicationArray
    slicedArray = chrBasedReplicationArray[int(start):int(end)]

    if (np.any(slicedArray)):
        #It must be full with at most -1 and +1
        uniqueValueArray = np.unique(slicedArray[np.nonzero(slicedArray)])

        # I expect the value of 1 (LEADING on the positive strand) or -1 (LAGGING on the positive strand) so size must be one.
        if (uniqueValueArray.size == 1):
            for uniqueValue in np.nditer(uniqueValueArray):
                # type(decileIndex) is numpy.ndarray
                slope = int(uniqueValue)
                #They have the same sign, multiplication (1,1) (-1,-1) must be 1
                if (slope*pyramidineStrand > 0):
                    all_types_leading_np_array += all_types_mask_array
                    subs_signature_SBS6_mutation_type_leading_np_array += subs_signature_SBS6_mutation_type_mask_array
                    subs_signature_SBS96_mutation_type_leading_np_array += subs_signature_SBS96_mutation_type_mask_array
                # They have the opposite sign, multiplication(1,-1) (-1,-)  must be -1
                elif (slope*pyramidineStrand < 0):
                    all_types_lagging_np_array += all_types_mask_array
                    subs_signature_SBS6_mutation_type_lagging_np_array += subs_signature_SBS6_mutation_type_mask_array
                    subs_signature_SBS96_mutation_type_lagging_np_array += subs_signature_SBS96_mutation_type_mask_array

        elif ((uniqueValueArray.size==2) and (pyramidineStrand!=0)):
            # Increment both LEADING and LAGGING
            all_types_leading_np_array += all_types_mask_array
            all_types_lagging_np_array += all_types_mask_array
            subs_signature_SBS6_mutation_type_leading_np_array += subs_signature_SBS6_mutation_type_mask_array
            subs_signature_SBS6_mutation_type_lagging_np_array += subs_signature_SBS6_mutation_type_mask_array
            subs_signature_SBS96_mutation_type_leading_np_array += subs_signature_SBS96_mutation_type_mask_array
            subs_signature_SBS96_mutation_type_lagging_np_array += subs_signature_SBS96_mutation_type_mask_array
        elif (uniqueValueArray.size>2):
            print('There is a situation!!!')
        else:
            print('There is a situation!!!')

    if sample_based:
        if (np.any(slicedArray)):
            # It must be full with at most -1 and +1
            uniqueValueArray = np.unique(slicedArray[np.nonzero(slicedArray)])

            # I expect the value of 1 (LEADING on the positive strand) or -1 (LAGGING on the positive strand) so size must be one.
            if (uniqueValueArray.size == 1):
                for uniqueValue in np.nditer(uniqueValueArray):
                    # type(decileIndex) is numpy.ndarray
                    slope = int(uniqueValue)
                    # They have the same sign, multiplication (1,1) (-1,-1) must be 1
                    if (slope * pyramidineStrand > 0):
                        all_samples_all_types_leading_np_array[sample_index] += all_types_mask_array
                        all_samples_subs_signature_mutation_type_leading_np_array[sample_index] += subs_signature_SBS6_mutation_type_mask_array
                    # They have the opposite sign, multiplication(1,-1) (-1,-)  must be -1
                    elif (slope * pyramidineStrand < 0):
                        all_samples_all_types_lagging_np_array[sample_index] += all_types_mask_array
                        all_samples_subs_signature_mutation_type_lagging_np_array[sample_index] += subs_signature_SBS6_mutation_type_mask_array
            elif ((uniqueValueArray.size == 2) and (pyramidineStrand != 0)):
                # Increment both LEADING and LAGGING
                all_samples_all_types_leading_np_array[sample_index] += all_types_mask_array
                all_samples_all_types_lagging_np_array[sample_index] += all_types_mask_array
                all_samples_subs_signature_mutation_type_leading_np_array[sample_index] += subs_signature_SBS6_mutation_type_mask_array
                all_samples_subs_signature_mutation_type_lagging_np_array[sample_index] += subs_signature_SBS6_mutation_type_mask_array
            elif (uniqueValueArray.size > 2):
                print('There is a situation!!!')
            else:
                print('There is a situation!!!')



# For df split
# July 28, 2020
# Using numpy arrays
#   if mutationPyramidineStrand and slope have the same sign increase LEADING STRAND count
#   else mutationPyramidineStrand and slope have the opposite sign increase LAGGING STRAND count
# sample_based for further usage
def searchAllMutationOnReplicationStrandArray_using_list_comprehension_using_numpy_array_for_df_split(
        mutation_row,
        chrBasedReplicationArray,
        six_mutation_types_np_array,
        ordered_sbs_signatures_cutoffs,
        ordered_dbs_signatures_cutoffs,
        ordered_id_signatures_cutoffs,
        df_columns_subs_signatures_mask_array,
        df_columns_dinucs_signatures_mask_array,
        df_columns_indels_signatures_mask_array,
        six_mutation_types_default_zeros_array,
        subs_signatures_default_zeros_array,
        dinucs_signatures_default_zeros_array,
        indels_signatures_default_zeros_array,
        subs_signatures_mutation_types_default_zeros_array,
        all_types_leading_np_array,
        all_types_lagging_np_array,
        subs_signature_mutation_type_leading_np_array,
        subs_signature_mutation_type_lagging_np_array,
        sample_based,
        discreet_mode,
        df_columns):

    indexofStart = np.where(df_columns == START)[0][0]
    indexofPyrimidineStrand = np.where(df_columns == PYRAMIDINESTRAND)[0][0]
    indexofSample = np.where(df_columns == SAMPLE)[0][0]
    indexofType = np.where(df_columns == TYPE)[0][0]

    start = mutation_row[indexofStart]
    pyramidineStrand = mutation_row[indexofPyrimidineStrand]
    sample = mutation_row[indexofSample]
    my_type=mutation_row[indexofType]

    mutationType = None
    subs_signature_mutation_type_mask_array=subs_signatures_mutation_types_default_zeros_array

    #############################################################################################################
    if(my_type==SUBS):
        end = start+1
        #e.g.: C>A
        indexofMutation = np.where(df_columns == MUTATION)[0][0]
        mutationType = mutation_row[indexofMutation]

        #six_mutation_types_mask_array.shape (6,)
        six_mutation_types_mask_array= np.where(six_mutation_types_np_array == mutationType, 1, 0)

        probabilities = mutation_row[df_columns_subs_signatures_mask_array]

        if discreet_mode:
            threshold_mask_array = np.greater_equal(probabilities, ordered_sbs_signatures_cutoffs)
            # Convert True into 1, and False into 0
            # subs_signatures_mask_array.shape (num_of_subs_signatures,)
            subs_signatures_mask_array = threshold_mask_array.astype(int)
        else:
            subs_signatures_mask_array = np.array(probabilities).astype(float)

        #Concetanate
        all_types_mask_array= np.concatenate((six_mutation_types_mask_array,
                                              subs_signatures_mask_array,
                                              dinucs_signatures_default_zeros_array,
                                              indels_signatures_default_zeros_array), axis=None)

        # Add one more dimension to subs_signatures_mask_array and six_mutation_types_mask_array
        # subs_signatures_mask_array_2d.shape (1,num_of_subs_signatures)
        subs_signatures_mask_array_2d = np.expand_dims(subs_signatures_mask_array, axis=0)

        # six_mutation_types_mask_array_2d.shape (1,6)
        six_mutation_types_mask_array_2d = np.expand_dims(six_mutation_types_mask_array, axis=0)

        # multiply subs_signatures_mask_array times six_mutation_types_mask_array
        subs_signature_mutation_type_mask_array = subs_signatures_mask_array_2d.T * six_mutation_types_mask_array_2d

    elif (my_type==DINUCS):
        end = start+2

        probabilities = mutation_row[df_columns_dinucs_signatures_mask_array]

        if discreet_mode:
            threshold_mask_array = np.greater_equal(probabilities, ordered_dbs_signatures_cutoffs)
            # Convert True into 1, and False into 0
            dinucs_signatures_mask_array = threshold_mask_array.astype(int)
        else:
            dinucs_signatures_mask_array = np.array(probabilities).astype(float)

        #Concetanate
        all_types_mask_array= np.concatenate((six_mutation_types_default_zeros_array,
                                              subs_signatures_default_zeros_array,
                                              dinucs_signatures_mask_array,
                                              indels_signatures_default_zeros_array), axis=None)

    elif (my_type==INDELS):
        indexofLength = np.where(df_columns == LENGTH)[0][0]
        end = start+int(mutation_row[indexofLength])

        probabilities = mutation_row[df_columns_indels_signatures_mask_array]

        if discreet_mode:
            threshold_mask_array = np.greater_equal(probabilities, ordered_id_signatures_cutoffs)
            # Convert True into 1, and False into 0
            indels_signatures_mask_array = threshold_mask_array.astype(int)
        else:
            indels_signatures_mask_array = np.array(probabilities).astype(float)

        #Concetanate
        all_types_mask_array= np.concatenate((six_mutation_types_default_zeros_array,
                                              subs_signatures_default_zeros_array,
                                              dinucs_signatures_default_zeros_array,
                                              indels_signatures_mask_array), axis=None)

    #############################################################################################################

    #############################################################################################################
    #if there is overlap with chrBasedReplicationArray
    slicedArray = chrBasedReplicationArray[int(start):int(end)]

    if (np.any(slicedArray)):
        #It must be full with at most -1 and +1
        uniqueValueArray = np.unique(slicedArray[np.nonzero(slicedArray)])

        # I expect the value of 1 (LEADING on the positive strand) or -1 (LAGGING on the positive strand) so size must be one.
        if (uniqueValueArray.size == 1):
            for uniqueValue in np.nditer(uniqueValueArray):
                # type(decileIndex) is numpy.ndarray
                slope = int(uniqueValue)
                #They have the same sign, multiplication (1,1) (-1,-1) must be 1
                if (slope*pyramidineStrand > 0):
                    all_types_leading_np_array += all_types_mask_array
                    subs_signature_mutation_type_leading_np_array += subs_signature_mutation_type_mask_array
                # They have the opposite sign, multiplication(1,-1) (-1,-)  must be -1
                elif (slope*pyramidineStrand < 0):
                    all_types_lagging_np_array += all_types_mask_array
                    subs_signature_mutation_type_lagging_np_array += subs_signature_mutation_type_mask_array
        elif ((uniqueValueArray.size==2) and (pyramidineStrand!=0)):
            #Increment both LEADING and LAGGING
            all_types_leading_np_array += all_types_mask_array
            all_types_lagging_np_array += all_types_mask_array
            subs_signature_mutation_type_leading_np_array += subs_signature_mutation_type_mask_array
            subs_signature_mutation_type_lagging_np_array += subs_signature_mutation_type_mask_array
        elif (uniqueValueArray.size>2):
            print('There is a situation!!!')
        else:
            print('There is a situation!!!')



# This code checks whether valleys and peaks are one after another, not two consecutive elements are both valley and peak.
def check_for_validness(chrBased_valleys_peaks_df):
    former_row_type = None

    for index, row in chrBased_valleys_peaks_df.iterrows():
        if former_row_type is None:
            former_row_type = row['type']
        elif (row['type'] == former_row_type):
            print('DEBUG There is a problem in repli-seq valleys and peaks.')
            return False
        else:
            former_row_type = row['type']

    return True


def get_chr_based_replication_strand_array_for_callback(chrLong, chromSize, repliseq_signal_df, valleys_df, peaks_df):
    # global THRESHOLD_CONSECUTIVE_LONG_STRETCH_LENGTH
    # most_frequent_length = (repliseq_signal_df[END] - repliseq_signal_df[START]).value_counts().idxmax()
    # THRESHOLD_CONSECUTIVE_LONG_STRETCH_LENGTH = most_frequent_length * 10
    # print('DEBUG most_frequent_length:', most_frequent_length, 'THRESHOLD_CONSECUTIVE_LONG_STRETCH_LENGTH:', THRESHOLD_CONSECUTIVE_LONG_STRETCH_LENGTH)

    chrBased_replication_array = get_chr_based_replication_strand_array(chrLong,
                                                                        chromSize,
                                                                        repliseq_signal_df,
                                                                        valleys_df,
                                                                        peaks_df)
    return (chrLong, chrBased_replication_array)


def get_chr_based_replication_strand_array(chrLong, chromSize, repliseq_signal_df, valleys_df, peaks_df):

    # Read chrBased_repli_seq_signal_df
    chrBased_repli_seq_signal_df = repliseq_signal_df[repliseq_signal_df[CHROM] == chrLong]

    chrBased_valleys_df = valleys_df[valleys_df[CHROM] == chrLong].copy()
    chrBased_valleys_df['type'] = 'Valley'
    chrBased_valleys_df.astype(dtype={START: int, END: int})

    chrBased_peaks_df = peaks_df[peaks_df[CHROM] == chrLong].copy()
    chrBased_peaks_df['type'] = 'Peak'
    chrBased_peaks_df.astype(dtype={START: int, END: int})

    # Concat valleys and peaks vertically
    chrBased_valleys_peaks_df = pd.concat([chrBased_valleys_df, chrBased_peaks_df], axis=0)

    # Sort valleys and peaks in ascending order
    chrBased_valleys_peaks_df.sort_values(START, inplace=True)

    # start the index from zero
    chrBased_valleys_peaks_df.reset_index(drop=True, inplace=True)

    if ((chrBased_repli_seq_signal_df is not None) and
            (not chrBased_repli_seq_signal_df.empty) and
            (check_for_validness(chrBased_valleys_peaks_df))):
        chrBased_replication_array = fill_chr_based_replication_strand_array(chrLong,
                                                                             chromSize,
                                                                             chrBased_repli_seq_signal_df,
                                                                             chrBased_valleys_peaks_df)

        return chrBased_replication_array
    else:
        return None



def fill_chr_based_replication_strand_array(chrLong,
                                            chromSize,
                                            chrBased_repli_seq_signal_df,
                                            chrBased_valleys_peaks_df):

    # +1 means leading strand, -1 means lagging strand
    # we will fill this array using smoothedSignal, peaks and valleys for each chromosome
    chrBased_replication_array = np.zeros(chromSize, dtype=np.int8)

    start_column_index = chrBased_repli_seq_signal_df.columns.get_loc(START)
    end_column_index = chrBased_repli_seq_signal_df.columns.get_loc(END)

    start = chrBased_repli_seq_signal_df.iloc[0, start_column_index]  # get the first row start
    end = chrBased_repli_seq_signal_df.iloc[-1, end_column_index]  # get the last row end

    # To remain conservative in downstream assignments,
    # we removed the last 25kb of the latest zones of the replicating domains
    start += 25000
    end -= 25000

    # Step1 Find the transition zones
    chrBased_transition_zones_list = find_long_stretches_of_consistent_transition_zones(chrLong,
                                                                                        start,
                                                                                        end,
                                                                                        chrBased_repli_seq_signal_df,
                                                                                        chrBased_valleys_peaks_df)

    labels = ['chr', 'start', 'end', 'slope_direction', 'length']
    chrBased_transition_zones_df = pd.DataFrame.from_records(chrBased_transition_zones_list, columns=labels)

    # Step2 Fill the replication array using transition zones
    chrBased_transition_zones_df.apply(fill_replication_strand_array,
                                       chrBased_replication_array=chrBased_replication_array,
                                       axis=1)

    return chrBased_replication_array

def find_repli_seq_peaks_valleys_using_scipy(repli_seq_df):

    grouped = repli_seq_df.groupby(CHROM)

    valleys_df_list = []
    peaks_df_list = []

    for chrom, chrBased_replicationtimedata_df in grouped:

        # sort w.r.t. start
        chrBased_replicationtimedata_df = chrBased_replicationtimedata_df.sort_values(by=[START], ascending=True)

        # reset_index
        chrBased_replicationtimedata_df = chrBased_replicationtimedata_df.reset_index()

        # Find local minima and maxima using peak detection
        maxima_indices, _ = find_peaks(chrBased_replicationtimedata_df[SIGNAL], prominence=0) # 0
        minima_indices, _ = find_peaks(-(chrBased_replicationtimedata_df[SIGNAL]), prominence=0) # 0

        # Extract local minima and maxima from indices
        local_maxima_df = chrBased_replicationtimedata_df.iloc[maxima_indices]
        local_minima_df = chrBased_replicationtimedata_df.iloc[minima_indices]

        local_maxima_df = local_maxima_df[[CHROM, START, END]]
        local_minima_df = local_minima_df[[CHROM, START, END]]

        peaks_df_list.append(local_maxima_df)
        valleys_df_list.append(local_minima_df)

    # Vertically concatenate dataframes using reduce
    all_peaks_df = reduce(lambda x, y: pd.concat([x, y], axis=0), peaks_df_list)
    all_valleys_df = reduce(lambda x, y: pd.concat([x, y], axis=0), valleys_df_list)

    return all_peaks_df, all_valleys_df


def read_repliseq_dataframes(smoothedWaveletRepliseqDataFilename,
                             valleysBEDFilename,
                             peaksBEDFilename,
                             chromNamesList,
                             log_file):

    # Read the Smoothed Wavelet Replication Time Signal
    file_extension = os.path.splitext(os.path.basename(smoothedWaveletRepliseqDataFilename))[1]
    if (file_extension.lower() == '.wig'):
        filetype_BEDGRAPH = isFileTypeBedGraph(smoothedWaveletRepliseqDataFilename)
        if filetype_BEDGRAPH:
            repliseq_wavelet_signal_df = pd.read_csv(smoothedWaveletRepliseqDataFilename, sep='\t', comment='#',
                                                     header=None, names=[CHROM, START, END, SIGNAL])
        else:
            repliseq_wavelet_signal_df = readWig_with_fixedStep_variableStep(smoothedWaveletRepliseqDataFilename)
    elif (file_extension.lower() == '.bedgraph'):
        repliseq_wavelet_signal_df = pd.read_csv(smoothedWaveletRepliseqDataFilename, sep='\t', comment='#',
                                                 header=None, names=[CHROM, START, END, SIGNAL])
    elif (file_extension.lower() == '.bed'):
        repliseq_wavelet_signal_df = pd.read_csv(smoothedWaveletRepliseqDataFilename, sep='\t', comment='#',
                                                 header=None, names=[CHROM, START, END, SIGNAL])

    # Remove rows with chromosomes that are not in chromNamesList
    repliseq_wavelet_signal_df = repliseq_wavelet_signal_df[repliseq_wavelet_signal_df[CHROM].isin(chromNamesList)]

    if valleysBEDFilename is not None:
        # Read the Valleys
        discard_signal = True
        valleys_df = readFileInBEDFormat(valleysBEDFilename, discard_signal, log_file)
        valleys_df[END] = valleys_df[END] - 1

    if peaksBEDFilename is not None:
        # Read the Peaks
        discard_signal = True
        peaks_df = readFileInBEDFormat(peaksBEDFilename, discard_signal, log_file)
        peaks_df[END] = peaks_df[END] - 1

    if (valleysBEDFilename is None) and (peaksBEDFilename is None):
        peaks_df, valleys_df = find_repli_seq_peaks_valleys_using_scipy(repliseq_wavelet_signal_df)

    # Remove rows with chromosomes that are not in chromNamesList
    valleys_df = valleys_df[valleys_df[CHROM].isin(chromNamesList)]
    peaks_df = peaks_df[peaks_df[CHROM].isin(chromNamesList)]

    log_out = open(log_file, 'a')
    repliseq_wavelet_signal_df_unique_chr_names = repliseq_wavelet_signal_df[CHROM].unique()
    valleys_df_unique_chr_names = valleys_df[CHROM].unique()
    peaks_df_unique_chr_names = peaks_df[CHROM].unique()

    print('After considering only chromosomes in chromNamesList --- ' +
          'Chromosome names in repliseq_wavelet_signal_df: %s repliseq_wavelet_signal_df.shape(%d,%d)\n'
          %(repliseq_wavelet_signal_df_unique_chr_names,
            repliseq_wavelet_signal_df.shape[0],
            repliseq_wavelet_signal_df.shape[1]),
          file=log_out)

    print('After considering only chromosomes in chromNamesList --- ' +
          'Chromosome names in valleys_df: %s valleys_df.shape(%d,%d)\n'
          %(valleys_df_unique_chr_names, valleys_df.shape[0], valleys_df.shape[1]), file=log_out)

    print('After considering only chromosomes in chromNamesList --- ' +
          'Chromosome names in peaks_df: %s peaks_df.shape(%d,%d)\n'
          %(peaks_df_unique_chr_names, peaks_df.shape[0], peaks_df.shape[1]), file=log_out)

    log_out.close()

    return repliseq_wavelet_signal_df, valleys_df, peaks_df


# For df_split
# Using Numpy Arrays
# Search for all mutatitions
def searchAllMutationsOnReplicationStrandArray_for_df_split(chrBased_simBased_combined_df_split,
                                            chrBased_replication_array,
                                            sim_num,
                                            six_mutation_types_np_array,
                                            ordered_sbs_signatures,
                                            ordered_dbs_signatures,
                                            ordered_id_signatures,
                                            ordered_sbs_signatures_cutoffs,
                                            ordered_dbs_signatures_cutoffs,
                                            ordered_id_signatures_cutoffs,
                                            all_types_leading_np_array,
                                            all_types_lagging_np_array,
                                            subs_signature_mutation_type_leading_np_array,
                                            subs_signature_mutation_type_lagging_np_array,
                                            sample_based,
                                            discreet_mode,
                                            log_file,
                                            verbose):


    if ((chrBased_simBased_combined_df_split is not None) and (not chrBased_simBased_combined_df_split.empty)):
        if verbose:
            log_out = open(log_file, 'a')
            print('\tVerbose Worker pid %s SBS searchMutationd_comOnReplicationStrandArray_simulations_integrated starts %s MB' % (str(os.getpid()), memory_usage()), file=log_out)
            log_out.close()

        # df_columns numpy array
        df_columns = chrBased_simBased_combined_df_split.columns.values

        ###############################################################################
        ################################ Initialization ###############################
        ###############################################################################
        df_columns_subs_signatures_mask_array = np.isin(df_columns,ordered_sbs_signatures)
        df_columns_dinucs_signatures_mask_array = np.isin(df_columns,ordered_dbs_signatures)
        df_columns_indels_signatures_mask_array = np.isin(df_columns,ordered_id_signatures)

        six_mutation_types_default_zeros_array= np.zeros(six_mutation_types_np_array.size) # dtype=int
        subs_signatures_default_zeros_array = np.zeros(ordered_sbs_signatures.size) # dtype=int
        dinucs_signatures_default_zeros_array = np.zeros(ordered_dbs_signatures.size) # dtype=int
        indels_signatures_default_zeros_array = np.zeros(ordered_id_signatures.size) # dtype=int
        subs_signatures_mutation_types_default_zeros_array = np.zeros((ordered_sbs_signatures.size, six_mutation_types_np_array.size)) # dtype=int
        ###############################################################################
        ################################ Initialization ###############################
        ###############################################################################

        # In list comprehesion, mutation_row becomes <class 'numpy.ndarray'>
        # In apply, mutation_row becomes <class 'pandas.core.series.Series'>
        # Therefore, list comprehesion is adopted.
        [searchAllMutationOnReplicationStrandArray_using_list_comprehension_using_numpy_array_for_df_split(mutation_row,
                                                                            chrBased_replication_array,
                                                                            six_mutation_types_np_array,
                                                                            ordered_sbs_signatures_cutoffs,
                                                                            ordered_dbs_signatures_cutoffs,
                                                                            ordered_id_signatures_cutoffs,
                                                                            df_columns_subs_signatures_mask_array,
                                                                            df_columns_dinucs_signatures_mask_array,
                                                                            df_columns_indels_signatures_mask_array,
                                                                            six_mutation_types_default_zeros_array,
                                                                            subs_signatures_default_zeros_array,
                                                                            dinucs_signatures_default_zeros_array,
                                                                            indels_signatures_default_zeros_array,
                                                                            subs_signatures_mutation_types_default_zeros_array,
                                                                            all_types_leading_np_array,
                                                                            all_types_lagging_np_array,
                                                                            subs_signature_mutation_type_leading_np_array,
                                                                            subs_signature_mutation_type_lagging_np_array,
                                                                            sample_based,
                                                                            discreet_mode,
                                                                            df_columns) for mutation_row in chrBased_simBased_combined_df_split.values]


        if verbose:
            log_out = open(log_file, 'a')
            print('\tVerbose Worker pid %s SBS searchMutationOnReplicationStrandArray_simulations_integrated ends %s MB' % (str(os.getpid()), memory_usage()), file=log_out)
            log_out.close()

    return(sim_num,
           all_types_leading_np_array,
           all_types_lagging_np_array,
           subs_signature_mutation_type_leading_np_array,
           subs_signature_mutation_type_lagging_np_array)


# Using Numpy Arrays
# Search for all mutations
def searchAllMutationsOnReplicationStrandArray(chrBased_simBased_subs_df,
                                            chrBased_simBased_dinucs_df,
                                            chrBased_simBased_indels_df,
                                            chrBased_replication_array,
                                            sim_num,
                                            SBS6_mutation_types_np_array,
                                            SBS96_mutation_types_np_array,
                                            ordered_sbs_signatures,
                                            ordered_dbs_signatures,
                                            ordered_id_signatures,
                                            ordered_sbs_signatures_cutoffs,
                                            ordered_dbs_signatures_cutoffs,
                                            ordered_id_signatures_cutoffs,
                                            all_types_leading_np_array,
                                            all_types_lagging_np_array,
                                            subs_signature_SBS6_mutation_type_leading_np_array,
                                            subs_signature_SBS6_mutation_type_lagging_np_array,
                                            subs_signature_SBS96_mutation_type_leading_np_array,
                                            subs_signature_SBS96_mutation_type_lagging_np_array,
                                            all_samples_all_types_leading_np_array,
                                            all_samples_all_types_lagging_np_array,
                                            all_samples_subs_signature_mutation_type_leading_np_array,
                                            all_samples_subs_signature_mutation_type_lagging_np_array,
                                            sample_based,
                                            all_samples_np_array,
                                            discreet_mode,
                                            default_cutoff,
                                            log_file,
                                            verbose):

    number_of_sbs_signatures = ordered_sbs_signatures.size
    number_of_dbs_signatures = ordered_dbs_signatures.size
    number_of_id_signatures = ordered_id_signatures.size

    ###############################################################################
    ################################ Initialization ###############################
    ###############################################################################
    SBS6_mutation_types_default_zeros_array = np.zeros(SBS6_mutation_types_np_array.size) # dtype=int
    SBS96_mutation_types_default_zeros_array = np.zeros(SBS96_mutation_types_np_array.size) # dtype=int
    subs_signatures_default_zeros_array = np.zeros(number_of_sbs_signatures) # dtype=int
    dinucs_signatures_default_zeros_array = np.zeros(number_of_dbs_signatures) # dtype=int
    indels_signatures_default_zeros_array = np.zeros(number_of_id_signatures) # dtype=int
    subs_signatures_SBS6_mutation_types_default_zeros_array = np.zeros((number_of_sbs_signatures, SBS6_mutation_types_np_array.size)) # dtype=int
    subs_signatures_SBS96_mutation_types_default_zeros_array = np.zeros((number_of_sbs_signatures, SBS96_mutation_types_np_array.size)) # dtype=int
    ###############################################################################
    ################################ Initialization ###############################
    ###############################################################################

    # SUBS
    if ((chrBased_simBased_subs_df is not None) and (not chrBased_simBased_subs_df.empty)):
        if verbose:
            log_out = open(log_file, 'a')
            print('\tVerbose Worker pid %s SBS searchMutationd_comOnReplicationStrandArray_simulations_integrated starts %s MB' % (str(os.getpid()), memory_usage()), file=log_out)
            log_out.close()

        # df_columns numpy array
        df_columns = chrBased_simBased_subs_df.columns.values
        df_columns_subs_signatures_mask_array = np.isin(df_columns, ordered_sbs_signatures)

        # In list comprehesion, mutation_row becomes <class 'numpy.ndarray'>
        # In apply, mutation_row becomes <class 'pandas.core.series.Series'>
        # Therefore, list comprehesion is adopted.
        [searchAllMutationOnReplicationStrandArray_using_list_comprehension_using_numpy_array(mutation_row,
                                                                            SUBS,
                                                                            chrBased_replication_array,
                                                                            SBS6_mutation_types_np_array,
                                                                            SBS96_mutation_types_np_array,
                                                                            ordered_sbs_signatures_cutoffs,
                                                                            df_columns_subs_signatures_mask_array,
                                                                            SBS6_mutation_types_default_zeros_array,
                                                                            SBS96_mutation_types_default_zeros_array,
                                                                            subs_signatures_default_zeros_array,
                                                                            dinucs_signatures_default_zeros_array,
                                                                            indels_signatures_default_zeros_array,
                                                                            subs_signatures_SBS6_mutation_types_default_zeros_array,
                                                                            subs_signatures_SBS96_mutation_types_default_zeros_array,
                                                                            all_types_leading_np_array,
                                                                            all_types_lagging_np_array,
                                                                            subs_signature_SBS6_mutation_type_leading_np_array,
                                                                            subs_signature_SBS6_mutation_type_lagging_np_array,
                                                                            subs_signature_SBS96_mutation_type_leading_np_array,
                                                                            subs_signature_SBS96_mutation_type_lagging_np_array,
                                                                            all_samples_all_types_leading_np_array,
                                                                            all_samples_all_types_lagging_np_array,
                                                                            all_samples_subs_signature_mutation_type_leading_np_array,
                                                                            all_samples_subs_signature_mutation_type_lagging_np_array,
                                                                            sample_based,
                                                                            all_samples_np_array,
                                                                            discreet_mode,
                                                                            default_cutoff,
                                                                            df_columns) for mutation_row in chrBased_simBased_subs_df.values]

    # DINUCS
    if ((chrBased_simBased_dinucs_df is not None) and (not chrBased_simBased_dinucs_df.empty)):
        if verbose:
            log_out = open(log_file, 'a')
            print('\tVerbose Worker pid %s SBS searchMutationd_comOnReplicationStrandArray_simulations_integrated starts %s MB' % (str(os.getpid()), memory_usage()), file=log_out)
            log_out.close()

        # df_columns numpy array
        df_columns = chrBased_simBased_dinucs_df.columns.values
        df_columns_dinucs_signatures_mask_array = np.isin(df_columns, ordered_dbs_signatures)

        # In list comprehesion, mutation_row becomes <class 'numpy.ndarray'>
        # In apply, mutation_row becomes <class 'pandas.core.series.Series'>
        # Therefore, list comprehesion is adopted.
        [searchAllMutationOnReplicationStrandArray_using_list_comprehension_using_numpy_array(mutation_row,
                                                                                              DINUCS,
                                                                                              chrBased_replication_array,
                                                                                              SBS6_mutation_types_np_array,
                                                                                              SBS96_mutation_types_np_array,
                                                                                              ordered_dbs_signatures_cutoffs,
                                                                                              df_columns_dinucs_signatures_mask_array,
                                                                                              SBS6_mutation_types_default_zeros_array,
                                                                                              SBS96_mutation_types_default_zeros_array,
                                                                                              subs_signatures_default_zeros_array,
                                                                                              dinucs_signatures_default_zeros_array,
                                                                                              indels_signatures_default_zeros_array,
                                                                                              subs_signatures_SBS6_mutation_types_default_zeros_array,
                                                                                              subs_signatures_SBS96_mutation_types_default_zeros_array,
                                                                                              all_types_leading_np_array,
                                                                                              all_types_lagging_np_array,
                                                                                              subs_signature_SBS6_mutation_type_leading_np_array,
                                                                                              subs_signature_SBS6_mutation_type_lagging_np_array,
                                                                                              subs_signature_SBS96_mutation_type_leading_np_array,
                                                                                              subs_signature_SBS96_mutation_type_lagging_np_array,
                                                                                              all_samples_all_types_leading_np_array,
                                                                                              all_samples_all_types_lagging_np_array,
                                                                                              all_samples_subs_signature_mutation_type_leading_np_array,
                                                                                              all_samples_subs_signature_mutation_type_lagging_np_array,
                                                                                              sample_based,
                                                                                              all_samples_np_array,
                                                                                              discreet_mode,
                                                                                              default_cutoff,
                                                                                              df_columns) for mutation_row in chrBased_simBased_dinucs_df.values]

    # INDELS
    if ((chrBased_simBased_indels_df is not None) and (not chrBased_simBased_indels_df.empty)):
        if verbose:
            log_out = open(log_file, 'a')
            print('\tVerbose Worker pid %s SBS searchMutationd_comOnReplicationStrandArray_simulations_integrated starts %s MB' % (str(os.getpid()), memory_usage()), file=log_out)
            log_out.close()

        # df_columns numpy array
        df_columns = chrBased_simBased_indels_df.columns.values
        df_columns_indels_signatures_mask_array = np.isin(df_columns, ordered_id_signatures)

        # In list comprehesion, mutation_row becomes <class 'numpy.ndarray'>
        # In apply, mutation_row becomes <class 'pandas.core.series.Series'>
        # Therefore, list comprehesion is adopted.
        [searchAllMutationOnReplicationStrandArray_using_list_comprehension_using_numpy_array(mutation_row,
                                                                                              INDELS,
                                                                                              chrBased_replication_array,
                                                                                              SBS6_mutation_types_np_array,
                                                                                              SBS96_mutation_types_np_array,
                                                                                              ordered_id_signatures_cutoffs,
                                                                                              df_columns_indels_signatures_mask_array,
                                                                                              SBS6_mutation_types_default_zeros_array,
                                                                                              SBS96_mutation_types_default_zeros_array,
                                                                                              subs_signatures_default_zeros_array,
                                                                                              dinucs_signatures_default_zeros_array,
                                                                                              indels_signatures_default_zeros_array,
                                                                                              subs_signatures_SBS6_mutation_types_default_zeros_array,
                                                                                              subs_signatures_SBS96_mutation_types_default_zeros_array,
                                                                                              all_types_leading_np_array,
                                                                                              all_types_lagging_np_array,
                                                                                              subs_signature_SBS6_mutation_type_leading_np_array,
                                                                                              subs_signature_SBS6_mutation_type_lagging_np_array,
                                                                                              subs_signature_SBS96_mutation_type_leading_np_array,
                                                                                              subs_signature_SBS96_mutation_type_lagging_np_array,
                                                                                              all_samples_all_types_leading_np_array,
                                                                                              all_samples_all_types_lagging_np_array,
                                                                                              all_samples_subs_signature_mutation_type_leading_np_array,
                                                                                              all_samples_subs_signature_mutation_type_lagging_np_array,
                                                                                              sample_based,
                                                                                              all_samples_np_array,
                                                                                              discreet_mode,
                                                                                              default_cutoff,
                                                                                              df_columns) for mutation_row in chrBased_simBased_indels_df.values]

        if verbose:
            log_out = open(log_file, 'a')
            print('\tVerbose Worker pid %s SBS searchMutationOnReplicationStrandArray_simulations_integrated ends %s MB' % (str(os.getpid()), memory_usage()), file=log_out)
            log_out.close()

    return(sim_num,
           all_types_leading_np_array,
           all_types_lagging_np_array,
           subs_signature_SBS6_mutation_type_leading_np_array,
           subs_signature_SBS6_mutation_type_lagging_np_array,
           subs_signature_SBS96_mutation_type_leading_np_array,
           subs_signature_SBS96_mutation_type_lagging_np_array,
           all_samples_all_types_leading_np_array,
           all_samples_all_types_lagging_np_array,
           all_samples_subs_signature_mutation_type_leading_np_array,
           all_samples_subs_signature_mutation_type_lagging_np_array)


def searchAllMutationsOnReplicationStrandArray_simbased_chrombased_splitbased(outputDir,
                                                                              jobname,
                                                                              chrLong,
                                                                              simNum,
                                                                              splitIndex,
                                                                              six_mutation_types_np_array,
                                                                              ordered_sbs_signatures,
                                                                              ordered_dbs_signatures,
                                                                              ordered_id_signatures,
                                                                              ordered_sbs_signatures_cutoffs,
                                                                              ordered_dbs_signatures_cutoffs,
                                                                              ordered_id_signatures_cutoffs,
                                                                              all_types_np_array_size,
                                                                              sample_based,
                                                                              discreet_mode,
                                                                              log_file,
                                                                              verbose):

    chr_based_replication_time_file_name = '%s_replication_time.npy' % (chrLong)
    chr_based_replication_time_file_path = os.path.join(outputDir, jobname, DATA, REPLICATIONSTRANDBIAS, LIB, CHRBASED,chr_based_replication_time_file_name)

    # Initialization
    all_types_leading_np_array = np.zeros((all_types_np_array_size)) # dtype=int
    all_types_lagging_np_array = np.zeros((all_types_np_array_size)) # dtype=int
    subs_signature_mutation_type_leading_np_array = np.zeros((ordered_sbs_signatures.size, six_mutation_types_np_array.size)) # dtype=int
    subs_signature_mutation_type_lagging_np_array = np.zeros((ordered_sbs_signatures.size, six_mutation_types_np_array.size)) # dtype=int

    if (os.path.exists(chr_based_replication_time_file_path)):
        chrBased_replication_array = np.load(chr_based_replication_time_file_path)
    else:
        chrBased_replication_array=None

    if chrBased_replication_array is not None:
        chrBased_simBased_combined_df_split = get_chrBased_simBased_combined_df_split(outputDir, jobname, chrLong,simNum, splitIndex)

        return searchAllMutationsOnReplicationStrandArray_for_df_split(chrBased_simBased_combined_df_split,
                                                          chrBased_replication_array,
                                                          simNum,
                                                          six_mutation_types_np_array,
                                                          ordered_sbs_signatures,
                                                          ordered_dbs_signatures,
                                                          ordered_id_signatures,
                                                          ordered_sbs_signatures_cutoffs,
                                                          ordered_dbs_signatures_cutoffs,
                                                          ordered_id_signatures_cutoffs,
                                                          all_types_leading_np_array,
                                                          all_types_lagging_np_array,
                                                          subs_signature_mutation_type_leading_np_array,
                                                          subs_signature_mutation_type_lagging_np_array,
                                                          sample_based,
                                                          discreet_mode,
                                                          log_file,
                                                          verbose)
    else:
        return (simNum,
                all_types_leading_np_array,
                all_types_lagging_np_array,
                subs_signature_mutation_type_leading_np_array,
                subs_signature_mutation_type_lagging_np_array)


# Read chromBased and simBased combined (SBS, DBS and ID) dataframe in the process
def searchAllMutationsOnReplicationStrandArray_simbased_chrombased(outputDir,
                                                                   jobname,
                                                                   chrLong,
                                                                   simNum,
                                                                   samples_of_interest,
                                                                   SBS6_mutation_types_np_array,
                                                                   SBS96_mutation_types_np_array,
                                                                   ordered_sbs_signatures,
                                                                   ordered_dbs_signatures,
                                                                   ordered_id_signatures,
                                                                   ordered_sbs_signatures_cutoffs,
                                                                   ordered_dbs_signatures_cutoffs,
                                                                   ordered_id_signatures_cutoffs,
                                                                   all_types_np_array_size,
                                                                   sample_based,
                                                                   all_samples_np_array,
                                                                   discreet_mode,
                                                                   default_cutoff,
                                                                   log_file,
                                                                   verbose):

    chr_based_replication_time_file_name = '%s_replication_time.npy' % (chrLong)
    chr_based_replication_time_file_path = os.path.join(outputDir, jobname, DATA, REPLICATIONSTRANDBIAS, LIB, CHRBASED,chr_based_replication_time_file_name)

    number_of_sbs_signatures = ordered_sbs_signatures.size

    # Initialization
    all_types_leading_np_array = np.zeros((all_types_np_array_size)) # dtype=int
    all_types_lagging_np_array = np.zeros((all_types_np_array_size)) # dtype=int
    subs_signature_SBS6_mutation_type_leading_np_array = np.zeros((number_of_sbs_signatures, SBS6_mutation_types_np_array.size)) # dtype=int
    subs_signature_SBS6_mutation_type_lagging_np_array = np.zeros((number_of_sbs_signatures, SBS6_mutation_types_np_array.size)) # dtype=int
    subs_signature_SBS96_mutation_type_leading_np_array = np.zeros((number_of_sbs_signatures, SBS96_mutation_types_np_array.size)) # dtype=int
    subs_signature_SBS96_mutation_type_lagging_np_array = np.zeros((number_of_sbs_signatures, SBS96_mutation_types_np_array.size)) # dtype=int

    all_samples_np_array_size = all_samples_np_array.size
    all_samples_all_types_leading_np_array = np.zeros((all_samples_np_array_size, all_types_np_array_size)) # dtype=int
    all_samples_all_types_lagging_np_array = np.zeros((all_samples_np_array_size, all_types_np_array_size)) # dtype=int
    all_samples_subs_signature_mutation_type_leading_np_array = np.zeros((all_samples_np_array_size, number_of_sbs_signatures, SBS6_mutation_types_np_array.size)) # dtype=int
    all_samples_subs_signature_mutation_type_lagging_np_array = np.zeros((all_samples_np_array_size, number_of_sbs_signatures, SBS6_mutation_types_np_array.size)) # dtype=int

    if (os.path.exists(chr_based_replication_time_file_path)):
        chrBased_replication_array = np.load(chr_based_replication_time_file_path)
    else:
        chrBased_replication_array=None

    if chrBased_replication_array is not None:
        chrBased_simBased_subs_df, chrBased_simBased_dinucs_df, chrBased_simBased_indels_df = get_chrBased_simBased_dfs(outputDir, jobname, chrLong, simNum)

        # filter chrbased_df for samples_of_interest
        if samples_of_interest is not None:

            if chrBased_simBased_subs_df is not None:
                chrBased_simBased_subs_df = chrBased_simBased_subs_df[chrBased_simBased_subs_df['Sample'].isin(samples_of_interest)]

            if chrBased_simBased_dinucs_df is not None:
                chrBased_simBased_dinucs_df = chrBased_simBased_dinucs_df[chrBased_simBased_dinucs_df['Sample'].isin(samples_of_interest)]

            if chrBased_simBased_indels_df is not None:
                chrBased_simBased_indels_df = chrBased_simBased_indels_df[chrBased_simBased_indels_df['Sample'].isin(samples_of_interest)]

        return searchAllMutationsOnReplicationStrandArray(chrBased_simBased_subs_df,
                                                          chrBased_simBased_dinucs_df,
                                                          chrBased_simBased_indels_df,
                                                          chrBased_replication_array,
                                                          simNum,
                                                          SBS6_mutation_types_np_array,
                                                          SBS96_mutation_types_np_array,
                                                          ordered_sbs_signatures,
                                                          ordered_dbs_signatures,
                                                          ordered_id_signatures,
                                                          ordered_sbs_signatures_cutoffs,
                                                          ordered_dbs_signatures_cutoffs,
                                                          ordered_id_signatures_cutoffs,
                                                          all_types_leading_np_array,
                                                          all_types_lagging_np_array,
                                                          subs_signature_SBS6_mutation_type_leading_np_array,
                                                          subs_signature_SBS6_mutation_type_lagging_np_array,
                                                          subs_signature_SBS96_mutation_type_leading_np_array,
                                                          subs_signature_SBS96_mutation_type_lagging_np_array,
                                                          all_samples_all_types_leading_np_array,
                                                          all_samples_all_types_lagging_np_array,
                                                          all_samples_subs_signature_mutation_type_leading_np_array,
                                                          all_samples_subs_signature_mutation_type_lagging_np_array,
                                                          sample_based,
                                                          all_samples_np_array,
                                                          discreet_mode,
                                                          default_cutoff,
                                                          log_file,
                                                          verbose)
    else:
        return (simNum,
                all_types_leading_np_array,
                all_types_lagging_np_array,
                subs_signature_SBS6_mutation_type_leading_np_array,
                subs_signature_SBS6_mutation_type_lagging_np_array,
                subs_signature_SBS96_mutation_type_leading_np_array,
                subs_signature_SBS96_mutation_type_lagging_np_array,
                all_samples_all_types_leading_np_array,
                all_samples_all_types_lagging_np_array,
                all_samples_subs_signature_mutation_type_leading_np_array,
                all_samples_subs_signature_mutation_type_lagging_np_array)


def read_create_write_replication_time_array_in_parallel(outputDir,
                                                         jobname,
                                                         chromNamesList,
                                                         chromSizesDict,
                                                         smoothedWaveletRepliseqDataFilename,
                                                         valleysBEDFilename,
                                                         peaksBEDFilename,
                                                         verbose,
                                                         log_file):

    repliseq_signal_df, valleys_df, peaks_df = read_repliseq_dataframes(smoothedWaveletRepliseqDataFilename,
                                                                        valleysBEDFilename,
                                                                        peaksBEDFilename,
                                                                        chromNamesList,
                                                                        log_file)

    def write_chrom_based_replication_array(result_tuple):
        chrLong = result_tuple[0]
        chrBased_replication_array = result_tuple[1]
        if (chrBased_replication_array is not None):
            os.makedirs(os.path.join(outputDir, jobname, DATA, REPLICATIONSTRANDBIAS, LIB, CHRBASED), exist_ok=True)
            # File name without extension
            chr_based_replication_time_file_name = '%s_replication_time' %(chrLong)
            chr_based_replication_time_file_path = os.path.join(outputDir,jobname,DATA,REPLICATIONSTRANDBIAS,LIB,CHRBASED,chr_based_replication_time_file_name)
            np.save(chr_based_replication_time_file_path, chrBased_replication_array)


    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=numofProcesses)

    jobs = []

    for chrLong in chromNamesList:
        chromSize = chromSizesDict[chrLong]
        jobs.append(pool.apply_async(get_chr_based_replication_strand_array_for_callback,
                                 args=(chrLong,
                                       chromSize,
                                       repliseq_signal_df,
                                       valleys_df,
                                       peaks_df,),
                                 callback=write_chrom_based_replication_array))

    log_out = open(log_file, 'a')

    # wait for all jobs to finish
    for job in jobs:
        if verbose: print('\tVerbose Write Chrom Based Replication Time Array for Replicatio Strand Bias Analysis Worker pid %s job.get():%s ' % (str(os.getpid()), job.get()), file=log_out)

    log_out.close()

    pool.close()
    pool.join()


########################################################################
# pool.imap_unordered : Fills poolInputList (jobs), sends poolInputList and accumulates results one by one.
# Slowest one but completes without memory problem.
# Can be updated and tested to read chrom based sim based mutations data and chrom based replication array in the worker process.
#
# pool.apply_async:  USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT
#
# pool.apply_async:  USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT_USING_POOL_INPUT_LIST
# You fill your pool input list therefore you decide how many jobs  to send at once.
# Faster than imap_unordered. Low memory usage.
# Can be updated to read chrom based sim based mutations data and chrom based replication array in the worker process.
# If USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM does not end because of memory error, this can be used.
# All 28/28 processes are running. When the jobs in the pool input list are finishing some processes waits others to finish.
#
# pool.apply_async:  USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM
# For each possible (chrLong,simNum) couple read the data and array on the worker process
# Fastest, consumes more memory than others. 22/28 processes are running. For Combined_PACWG_nonPCAWG Skin_Melanoma after 1 hour all 28/28 running.
def replicationStrandBiasAnalysis(outputDir,
                                  jobname,
                                  numofSimulations,
                                  samples_of_interest,
                                  job_tuples,
                                  sample_based,
                                  all_samples_np_array,
                                  chromSizesDict,
                                  chromNamesList,
                                  computation_type,
                                  smoothedWaveletRepliseqDataFilename,
                                  valleysBEDFilename,
                                  peaksBEDFilename,
                                  ordered_sbs_signatures,
                                  ordered_dbs_signatures,
                                  ordered_id_signatures,
                                  ordered_sbs_signatures_cutoffs,
                                  ordered_dbs_signatures_cutoffs,
                                  ordered_id_signatures_cutoffs,
                                  discreet_mode,
                                  default_cutoff,
                                  parallel_mode,
                                  log_file,
                                  verbose):

    log_out = open(log_file, 'a')
    print('\n#################################################################################', file=log_out)
    print('--- Replication Strand Asymmetry starts', file=log_out)
    log_out.close()

    read_create_write_replication_time_array_in_parallel(outputDir,
                                                         jobname,
                                                         chromNamesList,
                                                         chromSizesDict,
                                                         smoothedWaveletRepliseqDataFilename,
                                                         valleysBEDFilename,
                                                         peaksBEDFilename,
                                                         verbose,
                                                         log_file)

    SBS6_mutation_types_np_array = np.array([C2A, C2G, C2T, T2A, T2C, T2G])

    nucleotides = ['A', 'C', 'G', 'T']
    mutations = ['[C>A]', '[C>G]', '[C>T]', '[T>A]', '[T>C]', '[T>G]']
    SBS96_mutation_types_np_array = np.array([nucleotide_left + middle + nucleotide_right
                                                for nucleotide_left in nucleotides
                                                    for nucleotide_right in nucleotides
                                                        for middle in mutations])

    all_types_np_array = np.concatenate((SBS6_mutation_types_np_array,
                                         SBS96_mutation_types_np_array,
                                         ordered_sbs_signatures,
                                         ordered_dbs_signatures,
                                         ordered_id_signatures), axis=None)

    number_of_sbs_signatures = ordered_sbs_signatures.size
    sbs_signatures = ordered_sbs_signatures

    all_types_np_array_size = all_types_np_array.size

    # Initialization
    all_sims_all_types_leading_np_array = np.zeros((numofSimulations+1, all_types_np_array_size)) # dtype=int
    all_sims_all_types_lagging_np_array = np.zeros((numofSimulations+1, all_types_np_array_size)) # dtype=int
    all_sims_subs_signature_SBS6_mutation_type_leading_np_array = np.zeros((numofSimulations+1, number_of_sbs_signatures, SBS6_mutation_types_np_array.size)) # dtype=int
    all_sims_subs_signature_SBS6_mutation_type_lagging_np_array= np.zeros((numofSimulations+1, number_of_sbs_signatures, SBS6_mutation_types_np_array.size)) # dtype=int
    all_sims_subs_signature_SBS96_mutation_type_leading_np_array = np.zeros((numofSimulations+1, number_of_sbs_signatures, SBS96_mutation_types_np_array.size)) # dtype=int
    all_sims_subs_signature_SBS96_mutation_type_lagging_np_array= np.zeros((numofSimulations+1, number_of_sbs_signatures, SBS96_mutation_types_np_array.size)) # dtype=int

    if sample_based:
        # Initialization for accumulated arrays
        all_samples_np_array_size = all_samples_np_array.size
        all_sims_all_samples_all_types_leading_np_array = np.zeros((numofSimulations + 1, all_samples_np_array_size, all_types_np_array_size)) # dtype=int
        all_sims_all_samples_all_types_lagging_np_array = np.zeros((numofSimulations + 1, all_samples_np_array_size, all_types_np_array_size)) # dtype=int
        all_sims_all_samples_subs_signature_mutation_type_leading_np_array = np.zeros((numofSimulations + 1, all_samples_np_array_size, number_of_sbs_signatures, SBS6_mutation_types_np_array.size)) # dtype=int
        all_sims_all_samples_subs_signature_mutation_type_lagging_np_array = np.zeros((numofSimulations + 1, all_samples_np_array_size, number_of_sbs_signatures, SBS6_mutation_types_np_array.size)) # dtype=int

    # Accumulate Numpy Arrays
    def accumulate_np_arrays(result_tuple):
        sim_num = result_tuple[0]
        all_types_leading_np_array = result_tuple[1]
        all_types_lagging_np_array = result_tuple[2]
        subs_signature_SBS6_mutation_type_leading_np_array = result_tuple[3]
        subs_signature_SBS6_mutation_type_lagging_np_array = result_tuple[4]
        subs_signature_SBS96_mutation_type_leading_np_array = result_tuple[5]
        subs_signature_SBS96_mutation_type_lagging_np_array = result_tuple[6]

        all_sims_all_types_leading_np_array[sim_num] += all_types_leading_np_array
        all_sims_all_types_lagging_np_array[sim_num] += all_types_lagging_np_array
        all_sims_subs_signature_SBS6_mutation_type_leading_np_array[sim_num] += subs_signature_SBS6_mutation_type_leading_np_array
        all_sims_subs_signature_SBS6_mutation_type_lagging_np_array[sim_num] += subs_signature_SBS6_mutation_type_lagging_np_array
        all_sims_subs_signature_SBS96_mutation_type_leading_np_array[sim_num] += subs_signature_SBS96_mutation_type_leading_np_array
        all_sims_subs_signature_SBS96_mutation_type_lagging_np_array[sim_num] += subs_signature_SBS96_mutation_type_lagging_np_array
        # print('MONITOR ACCUMULATE', flush=True)

        if sample_based:
            all_samples_all_types_leading_np_array = result_tuple[5]
            all_samples_all_types_lagging_np_array = result_tuple[6]
            all_samples_subs_signature_mutation_type_leading_np_array = result_tuple[7]
            all_samples_subs_signature_mutation_type_lagging_np_array = result_tuple[8]

            all_sims_all_samples_all_types_leading_np_array[sim_num] += all_samples_all_types_leading_np_array
            all_sims_all_samples_all_types_lagging_np_array[sim_num] += all_samples_all_types_lagging_np_array
            all_sims_all_samples_subs_signature_mutation_type_leading_np_array[sim_num] += all_samples_subs_signature_mutation_type_leading_np_array
            all_sims_all_samples_subs_signature_mutation_type_lagging_np_array[sim_num] += all_samples_subs_signature_mutation_type_lagging_np_array


    sim_nums = range(0, numofSimulations + 1)
    sim_num_chr_tuples = ((sim_num, chrLong) for sim_num in sim_nums for chrLong in chromNamesList)

    if parallel_mode:
        jobs = []

        # April 30, 2020
        # read chrom based sim based mutations data and chrom based replication time data in each worker process
        if (computation_type == USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM):

            numofProcesses = multiprocessing.cpu_count()
            pool = multiprocessing.Pool(processes=numofProcesses)

            for simNum, chrLong in sim_num_chr_tuples:
                jobs.append(pool.apply_async(searchAllMutationsOnReplicationStrandArray_simbased_chrombased,
                                        args=(outputDir,
                                              jobname,
                                              chrLong,
                                              simNum,
                                              samples_of_interest,
                                              SBS6_mutation_types_np_array,
                                              SBS96_mutation_types_np_array,
                                              ordered_sbs_signatures,
                                              ordered_dbs_signatures,
                                              ordered_id_signatures,
                                              ordered_sbs_signatures_cutoffs,
                                              ordered_dbs_signatures_cutoffs,
                                              ordered_id_signatures_cutoffs,
                                              all_types_np_array_size,
                                              sample_based,
                                              all_samples_np_array,
                                              discreet_mode,
                                              default_cutoff,
                                              log_file,
                                              verbose,),
                                        callback=accumulate_np_arrays))
            pool.close()
            pool.join()

        elif (computation_type == USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT):
            numofProcesses = multiprocessing.cpu_count()
            pool = multiprocessing.Pool(processes=numofProcesses)

            for chrLong, simNum, splitIndex in job_tuples:
                jobs.append(pool.apply_async(searchAllMutationsOnReplicationStrandArray_simbased_chrombased_splitbased,
                                        args=(outputDir,
                                              jobname,
                                              chrLong,
                                              simNum,
                                              splitIndex,
                                              SBS6_mutation_types_np_array,
                                              ordered_sbs_signatures,
                                              ordered_dbs_signatures,
                                              ordered_id_signatures,
                                              ordered_sbs_signatures_cutoffs,
                                              ordered_dbs_signatures_cutoffs,
                                              ordered_id_signatures_cutoffs,
                                              all_types_np_array_size,
                                              sample_based,
                                              discreet_mode,
                                              log_file,
                                              verbose,),
                                        callback=accumulate_np_arrays))

            pool.close()
            pool.join()

    else:
        # Sequential_mode for testing, debugging and profiling purposes
        for simNum, chrLong in sim_num_chr_tuples:
            result_tuple = searchAllMutationsOnReplicationStrandArray_simbased_chrombased(outputDir,
                                               jobname,
                                               chrLong,
                                               simNum,
                                               samples_of_interest,
                                               SBS6_mutation_types_np_array,
                                               SBS96_mutation_types_np_array,
                                               ordered_sbs_signatures,
                                               ordered_dbs_signatures,
                                               ordered_id_signatures,
                                               ordered_sbs_signatures_cutoffs,
                                               ordered_dbs_signatures_cutoffs,
                                               ordered_id_signatures_cutoffs,
                                               all_types_np_array_size,
                                               sample_based,
                                               all_samples_np_array,
                                               discreet_mode,
                                               default_cutoff,
                                               log_file,
                                               verbose)
            accumulate_np_arrays(result_tuple)

    ############################################################################################################
    #####################################       Output starts      #############################################
    ############################################################################################################

    strand_bias = REPLICATIONSTRANDBIAS
    replication_strands = [LAGGING,LEADING]

    SBS6_np_arrays_list =[  all_sims_subs_signature_SBS6_mutation_type_lagging_np_array,
                            all_sims_subs_signature_SBS6_mutation_type_leading_np_array]

    # Write files for real data
    # For each signature all_sims_subs_signature_SBS96_mutation_type_lagging_np_array
    # For each signature all_sims_subs_signature_SBS96_mutation_type_leading_np_array
    write_sbs_signature_sbs96_mutation_type_replication_strand_bias(all_sims_subs_signature_SBS96_mutation_type_lagging_np_array[0],
                                                                    all_sims_subs_signature_SBS96_mutation_type_leading_np_array[0],
                                                                    SBS96_mutation_types_np_array,
                                                                    sbs_signatures,
                                                                    strand_bias,
                                                                    outputDir,
                                                                    jobname)

    write_signature_mutation_type_strand_bias_np_array_as_dataframe(SBS6_np_arrays_list,
                                                                    SBS6_mutation_types_np_array,
                                                                    sbs_signatures,
                                                                    strand_bias,
                                                                    replication_strands,
                                                                    outputDir,
                                                                    jobname)

    all_sims_all_types_strand_np_arrays_list =[all_sims_all_types_lagging_np_array,
                                               all_sims_all_types_leading_np_array]

    write_type_strand_bias_np_array_as_dataframe(all_sims_all_types_strand_np_arrays_list,
                                                all_types_np_array,
                                                strand_bias,
                                                replication_strands,
                                                outputDir,
                                                jobname)


    if sample_based:
        write_sample_based_strand1_strand2_as_dataframe(outputDir,
                                                        jobname,
                                                        numofSimulations,
                                                        strand_bias,
                                                        all_samples_np_array,
                                                        all_types_np_array,
                                                        all_sims_all_samples_all_types_lagging_np_array,
                                                        all_sims_all_samples_all_types_leading_np_array)

    ############################################################################################################
    #####################################       Output ends      ###############################################
    ############################################################################################################

    log_out = open(log_file, 'a')
    print('--- Replication Strand Asymmetry Analysis ends', file=log_out)
    print('#################################################################################\n', file=log_out)
    log_out.close()




