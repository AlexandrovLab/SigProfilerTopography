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
# This version use np.arrays
# Right now transcription strand bias analysis works for single point mutations and signatures.
# Constraints, Thresholds
# Please note that for sample based transcription strand bias analysis
# We consider samples with at least 1000 mutations both on transcribed and non-transcribed strands.
#############################################################

#############################################################
# What is transcription strand bias?
# It is the ratio of = (number of mutations on transcribed strand) / (number of mutations on un-transcribed strand)
#############################################################

import multiprocessing
import numpy as np
import os

from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIPTIONSTRAND
from SigProfilerTopography.source.commons.TopographyCommons import SAMPLE
from SigProfilerTopography.source.commons.TopographyCommons import MUTATION
from SigProfilerTopography.source.commons.TopographyCommons import MUTATIONLONG

from SigProfilerTopography.source.commons.TopographyCommons import TYPE

from SigProfilerTopography.source.commons.TopographyCommons import SUBS
from SigProfilerTopography.source.commons.TopographyCommons import INDELS
from SigProfilerTopography.source.commons.TopographyCommons import DINUCS

from SigProfilerTopography.source.commons.TopographyCommons import BIDIRECTIONAL
from SigProfilerTopography.source.commons.TopographyCommons import NONTRANSCRIBED
from SigProfilerTopography.source.commons.TopographyCommons import QUESTIONABLE
from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIBED
from SigProfilerTopography.source.commons.TopographyCommons import UNTRANSCRIBED

from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIPTIONSTRANDBIAS

from SigProfilerTopography.source.commons.TopographyCommons import write_type_strand_bias_np_array_as_dataframe
from SigProfilerTopography.source.commons.TopographyCommons import write_sample_type_strand_bias_np_array_as_dataframe

from SigProfilerTopography.source.commons.TopographyCommons import write_signature_mutation_type_strand_bias_np_array_as_dataframe
from SigProfilerTopography.source.commons.TopographyCommons import write_sample_signature_mutation_type_strand_np_array_as_dataframe

from SigProfilerTopography.source.commons.TopographyCommons import get_chrBased_simBased_combined_df_split
from SigProfilerTopography.source.commons.TopographyCommons import get_chrBased_simBased_combined_df
from SigProfilerTopography.source.commons.TopographyCommons import get_chrBased_simBased_dfs

from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM
from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT

from SigProfilerTopography.source.commons.TopographyCommons import memory_usage

from SigProfilerTopography.source.commons.TopographyCommons import C2A
from SigProfilerTopography.source.commons.TopographyCommons import C2G
from SigProfilerTopography.source.commons.TopographyCommons import C2T
from SigProfilerTopography.source.commons.TopographyCommons import T2A
from SigProfilerTopography.source.commons.TopographyCommons import T2C
from SigProfilerTopography.source.commons.TopographyCommons import T2G

from SigProfilerTopography.source.commons.TopographyCommons import SBS6_mutation_types_np_array
from SigProfilerTopography.source.commons.TopographyCommons import SBS96_mutation_types_np_array
from SigProfilerTopography.source.commons.TopographyCommons import DBS78_mutation_types_np_array
from SigProfilerTopography.source.commons.TopographyCommons import ID83_mutation_types_np_array

# For df_split
# April 24, 2020
# Updated July 29, 2020
def search_all_mutations_on_transcription_strand_array_using_list_comprehension_using_numpy_array_for_df_split(mutation_row,
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
                                                                            all_types_transcribed_np_array,
                                                                            all_types_untranscribed_np_array,
                                                                            all_types_nontranscribed_np_array,
                                                                            subs_signature_mutation_type_transcribed_np_array,
                                                                            subs_signature_mutation_type_untranscribed_np_array,
                                                                            subs_signature_mutation_type_nontranscribed_np_array,
                                                                            discreet_mode,
                                                                            df_columns):

    indexofTranscriptionStrand = np.where(df_columns == TRANSCRIPTIONSTRAND)[0][0]
    indexofType = np.where(df_columns == TYPE)[0][0]

    mutationTranscriptionStrand = mutation_row[indexofTranscriptionStrand]
    my_type = mutation_row[indexofType]

    subs_signatures_mutation_types_mask_array=subs_signatures_mutation_types_default_zeros_array

    if (my_type==SUBS):
        #e.g.: C>A
        indexofMutation = np.where(df_columns == MUTATION)[0][0]
        mutationType = mutation_row[indexofMutation]

        six_mutation_types_mask_array= np.where(six_mutation_types_np_array == mutationType, 1, 0)

        probabilities = mutation_row[df_columns_subs_signatures_mask_array]

        if discreet_mode:
            threshold_mask_array = np.greater_equal(probabilities, ordered_sbs_signatures_cutoffs)
            # Convert True into 1, and False into 0
            subs_signatures_mask_array = threshold_mask_array.astype(int)
        else:
            subs_signatures_mask_array = np.array(probabilities).astype(float)

        # Concetanate
        all_types_mask_array= np.concatenate((six_mutation_types_mask_array,
                                              subs_signatures_mask_array,
                                              dinucs_signatures_default_zeros_array,
                                              indels_signatures_default_zeros_array), axis=None)

        # Add one more dimension to subs_signatures_mask_array and six_mutation_types_mask_array
        subs_signatures_mask_array_2d = np.expand_dims(subs_signatures_mask_array, axis=0)
        six_mutation_types_mask_array_2d = np.expand_dims(six_mutation_types_mask_array, axis=0)

        # multiply subs_signatures_mask_array times six_mutation_types_mask_array
        subs_signatures_mutation_types_mask_array = subs_signatures_mask_array_2d.T * six_mutation_types_mask_array_2d
    elif (my_type == DINUCS):
        probabilities = mutation_row[df_columns_dinucs_signatures_mask_array]

        if discreet_mode:
            threshold_mask_array = np.greater_equal(probabilities, ordered_dbs_signatures_cutoffs)
            # Convert True into 1, and False into 0
            dinucs_signatures_mask_array = threshold_mask_array.astype(int)
        else:
            dinucs_signatures_mask_array = np.array(probabilities).astype(float)

        # Concetanate
        all_types_mask_array= np.concatenate((six_mutation_types_default_zeros_array,
                                              subs_signatures_default_zeros_array,
                                              dinucs_signatures_mask_array,
                                              indels_signatures_default_zeros_array), axis=None)

    elif (my_type == INDELS):
        probabilities = mutation_row[df_columns_indels_signatures_mask_array]

        if discreet_mode:
            threshold_mask_array = np.greater_equal(probabilities, ordered_id_signatures_cutoffs)
            # Convert True into 1, and False into 0
            indels_signatures_mask_array = threshold_mask_array.astype(int)
        else:
            indels_signatures_mask_array = np.array(probabilities).astype(float)

        # Concetanate
        all_types_mask_array= np.concatenate((six_mutation_types_default_zeros_array,
                                              subs_signatures_default_zeros_array,
                                              dinucs_signatures_default_zeros_array,
                                              indels_signatures_mask_array), axis=None)

    # Values on TranscriptionStrand column
    # N --> Non-transcribed
    # T --> Transcribed
    # U --> Untranscribed
    # Q --> Question Not known

    if (mutationTranscriptionStrand == 'U'):
        all_types_untranscribed_np_array += all_types_mask_array
        subs_signature_mutation_type_untranscribed_np_array += subs_signatures_mutation_types_mask_array

    elif (mutationTranscriptionStrand == 'T'):
        all_types_transcribed_np_array += all_types_mask_array
        subs_signature_mutation_type_transcribed_np_array += subs_signatures_mutation_types_mask_array

    elif (mutationTranscriptionStrand == 'B'):
        all_types_untranscribed_np_array += all_types_mask_array
        all_types_transcribed_np_array += all_types_mask_array
        subs_signature_mutation_type_untranscribed_np_array += subs_signatures_mutation_types_mask_array
        subs_signature_mutation_type_transcribed_np_array += subs_signatures_mutation_types_mask_array

    elif (mutationTranscriptionStrand == 'N'):
        all_types_nontranscribed_np_array += all_types_mask_array
        subs_signature_mutation_type_nontranscribed_np_array += subs_signatures_mutation_types_mask_array



def search_all_mutations_on_transcription_strand_array_using_list_comprehension_using_numpy_array(mutation_row,
                                                                            my_type,
                                                                            sample_based,
                                                                            all_samples_np_array,
                                                                            ordered_signatures_cutoffs,
                                                                            df_columns_signatures_mask_array,
                                                                            sbs_signatures_np_array,
                                                                            dbs_signatures_np_array,
                                                                            id_signatures_np_array,
                                                                            SBS96_mutation_types_np_array,
                                                                            DBS78_mutation_types_np_array,
                                                                            ID83_mutation_types_np_array,
                                                                            all_mutation_types_np_array,
                                                                            sample_mutation_type_strand_np_array,
                                                                            sample_sbs_signature_mutation_type_strand_np_array,
                                                                            sample_dbs_signature_mutation_type_strand_np_array,
                                                                            sample_id_signature_mutation_type_strand_np_array,
                                                                            discreet_mode,
                                                                            default_cutoff,
                                                                            df_columns):

    sbs_signature_index = None
    dbs_signature_index = None
    id_signature_index = None

    if sample_based:
        indexofSample = np.where(df_columns == SAMPLE)[0][0]
        mutation_sample = mutation_row[indexofSample]
        sample_index = np.where(all_samples_np_array == mutation_sample)[0][0]

    indexofTranscriptionStrand = np.where(df_columns == TRANSCRIPTIONSTRAND)[0][0]
    mutationTranscriptionStrand = mutation_row[indexofTranscriptionStrand]

    probabilities = mutation_row[df_columns_signatures_mask_array]

    ##########################################
    if (my_type == SUBS):
        # e.g.: C>A
        index_of_mutation_long = np.where(df_columns == MUTATIONLONG)[0][0]

        # e.g.: T:AA[C>A]AA
        mutation_type_long = mutation_row[index_of_mutation_long]
        SBS96_mutation_type = mutation_type_long[3:10]
        SBS96_mutation_type_index = np.where(SBS96_mutation_types_np_array == SBS96_mutation_type)[0][0]
        mutation_type_index = np.where(all_mutation_types_np_array == SBS96_mutation_type)[0][0]

        if discreet_mode:
            # Convert True into 1, and False into 0
            threshold_mask_array = np.greater_equal(probabilities, ordered_signatures_cutoffs)
            subs_signatures_mask_array = threshold_mask_array.astype(int)
        else:
            # new way
            probabilities[probabilities < default_cutoff] = 0
            subs_signatures_mask_array = np.array(probabilities).astype(float)
            # old way
            # subs_signatures_mask_array = np.array(probabilities).astype(float)

        if np.sum(subs_signatures_mask_array) > 0:
            sbs_signature = sbs_signatures_np_array[subs_signatures_mask_array == 1][0] # e.g. ['SBS2']
            sbs_signature_index = np.where(sbs_signatures_np_array == sbs_signature)[0][0]

    elif (my_type == DINUCS):
        index_of_mutation_long = np.where(df_columns == MUTATIONLONG)[0][0]

        # T:T[TC>AA]A
        mutation_type_long = mutation_row[index_of_mutation_long]
        DBS78_mutation_type = mutation_type_long[4:9]
        DBS78_mutation_type_index = np.where(DBS78_mutation_types_np_array == DBS78_mutation_type)[0][0]
        mutation_type_index = np.where(all_mutation_types_np_array == DBS78_mutation_type)[0][0]

        if discreet_mode:
            threshold_mask_array = np.greater_equal(probabilities, ordered_signatures_cutoffs)
            # Convert True into 1, and False into 0
            dinucs_signatures_mask_array = threshold_mask_array.astype(int)
        else:
            # new way
            probabilities[probabilities < default_cutoff] = 0
            dinucs_signatures_mask_array = np.array(probabilities).astype(float)
            # old way
            # dinucs_signatures_mask_array = np.array(probabilities).astype(float)

        if np.sum(dinucs_signatures_mask_array) > 0:
            dbs_signature = dbs_signatures_np_array[dinucs_signatures_mask_array == 1][0] # e.g. ['DBS2']
            dbs_signature_index = np.where(dbs_signatures_np_array == dbs_signature)[0][0]

    elif (my_type == INDELS):
        index_of_mutation_long = np.where(df_columns == MUTATIONLONG)[0][0]

        # T:1:Del:T:5
        mutation_type_long = mutation_row[index_of_mutation_long]
        ID83_mutation_type = mutation_type_long[2:]
        ID83_mutation_type_index = np.where(ID83_mutation_types_np_array == ID83_mutation_type)[0][0]
        mutation_type_index = np.where(all_mutation_types_np_array == ID83_mutation_type)[0][0]

        if discreet_mode:
            threshold_mask_array = np.greater_equal(probabilities, ordered_signatures_cutoffs)
            # Convert True into 1, and False into 0
            indels_signatures_mask_array = threshold_mask_array.astype(int)
        else:
            # new way
            probabilities[probabilities < default_cutoff] = 0
            indels_signatures_mask_array = np.array(probabilities).astype(float)
            # old way
            # indels_signatures_mask_array = np.array(probabilities).astype(float)

        if np.sum(indels_signatures_mask_array) > 0:
            id_signature = id_signatures_np_array[indels_signatures_mask_array == 1][0] # e.g. ['ID2']
            id_signature_index = np.where(id_signatures_np_array == id_signature)[0][0]

    # Values on TranscriptionStrand column
    # N --> Non-transcribed (Intergenic)
    # T --> Transcribed (Genic)
    # U --> Untranscribed (Genic)
    # B --> Both: Transcribed and Untranscribed
    # Q --> Question not known

    # Strand Index
    # Bidirectional no index is needed. Accumulated in Transcribed and Untranscribed
    # TRANSCRIBED --> 0
    # UNTRANSCRIBED --> 1
    # NONTRANSCRIBED --> 2
    # QUESTIONABLE --> 3

    if sample_based:
        if (mutationTranscriptionStrand == 'B'):
            sample_mutation_type_strand_np_array[sample_index][mutation_type_index][0] += 1
            sample_mutation_type_strand_np_array[sample_index][mutation_type_index][1] += 1

            if sbs_signature_index is not None:  # if not None
                sample_sbs_signature_mutation_type_strand_np_array[sample_index][sbs_signature_index][SBS96_mutation_type_index][0] += 1
                sample_sbs_signature_mutation_type_strand_np_array[sample_index][sbs_signature_index][SBS96_mutation_type_index][1] += 1

            if dbs_signature_index is not None:  # if not None
                sample_dbs_signature_mutation_type_strand_np_array[sample_index][dbs_signature_index][DBS78_mutation_type_index][0] += 1
                sample_dbs_signature_mutation_type_strand_np_array[sample_index][dbs_signature_index][DBS78_mutation_type_index][1] += 1

            if id_signature_index is not None:  # if not None
                sample_id_signature_mutation_type_strand_np_array[sample_index][id_signature_index][ID83_mutation_type_index][0] += 1
                sample_id_signature_mutation_type_strand_np_array[sample_index][id_signature_index][ID83_mutation_type_index][1] += 1

        elif (mutationTranscriptionStrand == 'N'):
            sample_mutation_type_strand_np_array[sample_index][mutation_type_index][2] += 1

            if sbs_signature_index is not None:  # if not None
                sample_sbs_signature_mutation_type_strand_np_array[sample_index][sbs_signature_index][SBS96_mutation_type_index][2] += 1

            if dbs_signature_index is not None:  # if not None
                sample_dbs_signature_mutation_type_strand_np_array[sample_index][dbs_signature_index][DBS78_mutation_type_index][2] += 1

            if id_signature_index is not None:  # if not None
                sample_id_signature_mutation_type_strand_np_array[sample_index][id_signature_index][ID83_mutation_type_index][2] += 1

        elif (mutationTranscriptionStrand == 'Q'):
            sample_mutation_type_strand_np_array[sample_index][mutation_type_index][3] += 1

            if sbs_signature_index is not None:  # if not None
                sample_sbs_signature_mutation_type_strand_np_array[sample_index][sbs_signature_index][SBS96_mutation_type_index][3] += 1

            if dbs_signature_index is not None:  # if not None
                sample_dbs_signature_mutation_type_strand_np_array[sample_index][dbs_signature_index][DBS78_mutation_type_index][3] += 1

            if id_signature_index is not None:  # if not None
                sample_id_signature_mutation_type_strand_np_array[sample_index][id_signature_index][ID83_mutation_type_index][3] += 1

        elif (mutationTranscriptionStrand == 'T'):
            sample_mutation_type_strand_np_array[sample_index][mutation_type_index][0] += 1

            if sbs_signature_index is not None:  # if not None
                sample_sbs_signature_mutation_type_strand_np_array[sample_index][sbs_signature_index][SBS96_mutation_type_index][0] += 1

            if dbs_signature_index is not None:  # if not None
                sample_dbs_signature_mutation_type_strand_np_array[sample_index][dbs_signature_index][DBS78_mutation_type_index][0] += 1

            if id_signature_index is not None:  # if not None
                sample_id_signature_mutation_type_strand_np_array[sample_index][id_signature_index][ID83_mutation_type_index][0] += 1

        if (mutationTranscriptionStrand == 'U'):
            sample_mutation_type_strand_np_array[sample_index][mutation_type_index][1] += 1

            if sbs_signature_index is not None:  # if not None
                sample_sbs_signature_mutation_type_strand_np_array[sample_index][sbs_signature_index][SBS96_mutation_type_index][1] += 1

            if dbs_signature_index is not None:  # if not None
                sample_dbs_signature_mutation_type_strand_np_array[sample_index][dbs_signature_index][DBS78_mutation_type_index][1] += 1

            if id_signature_index is not None:  # if not None
                sample_id_signature_mutation_type_strand_np_array[sample_index][id_signature_index][ID83_mutation_type_index][1] += 1


# For df_split
def searchAllMutations_for_df_split(chrBased_simBased_combined_df,
            sim_num,
            six_mutation_types_np_array,
            ordered_sbs_signatures,
            ordered_dbs_signatures,
            ordered_id_signatures,
            ordered_sbs_signatures_cutoffs,
            ordered_dbs_signatures_cutoffs,
            ordered_id_signatures_cutoffs,
            all_types_transcribed_np_array,
            all_types_untranscribed_np_array,
            all_types_nontranscribed_np_array,
            subs_signature_mutation_type_transcribed_np_array,
            subs_signature_mutation_type_untranscribed_np_array,
            subs_signature_mutation_type_nontranscribed_np_array,
            discreet_mode,
            log_file,
            verbose):

    if ((chrBased_simBased_combined_df is not None) and (not chrBased_simBased_combined_df.empty)):
        if verbose:
            log_out = open(log_file, 'a')
            print('Worker pid %s searchMutationUsingTranscriptionStrandColumn_simulations_integrated starts %s MB' % (str(os.getpid()), memory_usage()), file=log_out)
            log_out.close()

        # df_columns numpy array
        df_columns = chrBased_simBased_combined_df.columns.values

        ###############################################################################
        ################################ Initialization ###############################
        ###############################################################################
        df_columns_subs_signatures_mask_array = np.isin(df_columns,ordered_sbs_signatures)
        df_columns_dinucs_signatures_mask_array=np.isin(df_columns,ordered_dbs_signatures)
        df_columns_indels_signatures_mask_array= np.isin(df_columns,ordered_id_signatures)

        six_mutation_types_default_zeros_array= np.zeros(six_mutation_types_np_array.size) # dtype=int
        subs_signatures_default_zeros_array = np.zeros(ordered_sbs_signatures.size) # dtype=int
        dinucs_signatures_default_zeros_array = np.zeros(ordered_dbs_signatures.size) # dtype=int
        indels_signatures_default_zeros_array = np.zeros(ordered_id_signatures.size) # dtype=int
        subs_signatures_mutation_types_default_zeros_array= np.zeros((ordered_sbs_signatures.size,six_mutation_types_np_array.size)) # dtype=int
        ###############################################################################
        ################################ Initialization ###############################
        ###############################################################################


        # Search for each row using Numpy Array
        [search_all_mutations_on_transcription_strand_array_using_list_comprehension_using_numpy_array_for_df_split(mutation_row,
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
                                                                            all_types_transcribed_np_array,
                                                                            all_types_untranscribed_np_array,
                                                                            all_types_nontranscribed_np_array,
                                                                            subs_signature_mutation_type_transcribed_np_array,
                                                                            subs_signature_mutation_type_untranscribed_np_array,
                                                                            subs_signature_mutation_type_nontranscribed_np_array,
                                                                            discreet_mode,
                                                                            df_columns) for mutation_row in chrBased_simBased_combined_df.values]


        if verbose:
            log_out = open(log_file, 'a')
            print('Worker pid %s searchMutationUsingTranscriptionStrandColumn_simulations_integrated ends %s MB' % (str(os.getpid()), memory_usage()), file=log_out)
            log_out.close()

    return (sim_num,
            all_types_transcribed_np_array,
            all_types_untranscribed_np_array,
            all_types_nontranscribed_np_array,
            subs_signature_mutation_type_transcribed_np_array,
            subs_signature_mutation_type_untranscribed_np_array,
            subs_signature_mutation_type_nontranscribed_np_array)



def searchAllMutations(chrBased_simBased_subs_df,
            chrBased_simBased_dinucs_df,
            chrBased_simBased_indels_df,
            sim_num,
            SBS96_mutation_types_np_array,
            DBS78_mutation_types_np_array,
            ID83_mutation_types_np_array,
            all_mutation_types_np_array,
            sample_based,
            all_samples_np_array,
            ordered_sbs_signatures_np_array,
            ordered_dbs_signatures_np_array,
            ordered_id_signatures_np_array,
            ordered_sbs_signatures_cutoffs,
            ordered_dbs_signatures_cutoffs,
            ordered_id_signatures_cutoffs,
            sample_mutation_type_strand_np_array,
            sample_sbs_signature_mutation_type_strand_np_array,
            sample_dbs_signature_mutation_type_strand_np_array,
            sample_id_signature_mutation_type_strand_np_array,
            discreet_mode,
            default_cutoff,
            log_file,
            verbose):

    number_of_sbs_signatures = ordered_sbs_signatures_np_array.size

    # SUBS
    if ((chrBased_simBased_subs_df is not None) and (not chrBased_simBased_subs_df.empty)):
        if verbose:
            log_out = open(log_file, 'a')
            print('Worker pid %s searchMutationUsingTranscriptionStrandColumn_simulations_integrated starts %s MB' % (str(os.getpid()), memory_usage()), file=log_out)
            log_out.close()

        # df_columns numpy array
        df_columns = chrBased_simBased_subs_df.columns.values

        df_columns_subs_signatures_mask_array = np.isin(df_columns, ordered_sbs_signatures_np_array)

        # Search for each row using Numpy Array
        [search_all_mutations_on_transcription_strand_array_using_list_comprehension_using_numpy_array(mutation_row,
                                                                            SUBS,
                                                                            sample_based,
                                                                            all_samples_np_array,
                                                                            ordered_sbs_signatures_cutoffs,
                                                                            df_columns_subs_signatures_mask_array,
                                                                            ordered_sbs_signatures_np_array,
                                                                            ordered_dbs_signatures_np_array,
                                                                            ordered_id_signatures_np_array,
                                                                            SBS96_mutation_types_np_array,
                                                                            DBS78_mutation_types_np_array,
                                                                            ID83_mutation_types_np_array,
                                                                            all_mutation_types_np_array,
                                                                            sample_mutation_type_strand_np_array,
                                                                            sample_sbs_signature_mutation_type_strand_np_array,
                                                                            sample_dbs_signature_mutation_type_strand_np_array,
                                                                            sample_id_signature_mutation_type_strand_np_array,
                                                                            discreet_mode,
                                                                            default_cutoff,
                                                                            df_columns) for mutation_row in chrBased_simBased_subs_df.values]


    # DINUCS
    if ((chrBased_simBased_dinucs_df is not None) and (not chrBased_simBased_dinucs_df.empty)):
        if verbose:
            log_out = open(log_file, 'a')
            print('Worker pid %s searchMutationUsingTranscriptionStrandColumn_simulations_integrated starts %s MB' % (str(os.getpid()), memory_usage()), file=log_out)
            log_out.close()

        # df_columns numpy array
        df_columns = chrBased_simBased_dinucs_df.columns.values

        df_columns_dinucs_signatures_mask_array = np.isin(df_columns, ordered_dbs_signatures_np_array)

        # Search for each row using Numpy Array
        [search_all_mutations_on_transcription_strand_array_using_list_comprehension_using_numpy_array(mutation_row,
                                                                                              DINUCS,
                                                                                              sample_based,
                                                                                              all_samples_np_array,
                                                                                              ordered_dbs_signatures_cutoffs,
                                                                                              df_columns_dinucs_signatures_mask_array,
                                                                                            ordered_sbs_signatures_np_array,
                                                                                            ordered_dbs_signatures_np_array,
                                                                                            ordered_id_signatures_np_array,
                                                                                            SBS96_mutation_types_np_array,
                                                                                            DBS78_mutation_types_np_array,
                                                                                            ID83_mutation_types_np_array,
                                                                                            all_mutation_types_np_array,
                                                                                            sample_mutation_type_strand_np_array,
                                                                                            sample_sbs_signature_mutation_type_strand_np_array,
                                                                                            sample_dbs_signature_mutation_type_strand_np_array,
                                                                                            sample_id_signature_mutation_type_strand_np_array,
                                                                                              discreet_mode,
                                                                                              default_cutoff,
                                                                                              df_columns) for mutation_row in chrBased_simBased_dinucs_df.values]

    # INDELS
    if ((chrBased_simBased_indels_df is not None) and (not chrBased_simBased_indels_df.empty)):
        if verbose:
            log_out = open(log_file, 'a')
            print('Worker pid %s searchMutationUsingTranscriptionStrandColumn_simulations_integrated starts %s MB' % (str(os.getpid()), memory_usage()), file=log_out)
            log_out.close()

        # df_columns numpy array
        df_columns = chrBased_simBased_indels_df.columns.values

        df_columns_indels_signatures_mask_array = np.isin(df_columns, ordered_id_signatures_np_array)

        # Search for each row using Numpy Array
        [search_all_mutations_on_transcription_strand_array_using_list_comprehension_using_numpy_array(mutation_row,
                                                                                              INDELS,
                                                                                              sample_based,
                                                                                              all_samples_np_array,
                                                                                              ordered_id_signatures_cutoffs,
                                                                                              df_columns_indels_signatures_mask_array,
                                                                                           ordered_sbs_signatures_np_array,
                                                                                           ordered_dbs_signatures_np_array,
                                                                                           ordered_id_signatures_np_array,
                                                                                           SBS96_mutation_types_np_array,
                                                                                           DBS78_mutation_types_np_array,
                                                                                           ID83_mutation_types_np_array,
                                                                                           all_mutation_types_np_array,
                                                                                            sample_mutation_type_strand_np_array,
                                                                                            sample_sbs_signature_mutation_type_strand_np_array,
                                                                                            sample_dbs_signature_mutation_type_strand_np_array,
                                                                                            sample_id_signature_mutation_type_strand_np_array,
                                                                                              discreet_mode,
                                                                                              default_cutoff,
                                                                                              df_columns) for mutation_row in chrBased_simBased_indels_df.values]



    return (sim_num, #0
            sample_mutation_type_strand_np_array,  # 1
            sample_sbs_signature_mutation_type_strand_np_array,  # 2
            sample_dbs_signature_mutation_type_strand_np_array,  # 3
            sample_id_signature_mutation_type_strand_np_array)  # 4


def searchAllMutations_simbased_chrombased_splitbased(outputDir,
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
                                          discreet_mode,
                                          log_file,
                                          verbose):

    # Initialization
    all_types_transcribed_np_array = np.zeros((all_types_np_array_size)) # dtype=int
    all_types_untranscribed_np_array = np.zeros((all_types_np_array_size)) # dtype=int
    all_types_nontranscribed_np_array = np.zeros((all_types_np_array_size)) # dtype=int
    subs_signature_mutation_type_transcribed_np_array = np.zeros((ordered_sbs_signatures.size, six_mutation_types_np_array.size)) # dtype=int
    subs_signature_mutation_type_untranscribed_np_array = np.zeros((ordered_sbs_signatures.size, six_mutation_types_np_array.size)) # dtype=int
    subs_signature_mutation_type_nontranscribed_np_array = np.zeros((ordered_sbs_signatures.size, six_mutation_types_np_array.size)) # dtype=int

    chrBased_simBased_combined_df_split = get_chrBased_simBased_combined_df_split(outputDir, jobname, chrLong, simNum, splitIndex)

    return searchAllMutations_for_df_split(chrBased_simBased_combined_df_split,
                              simNum,
                              six_mutation_types_np_array,
                              ordered_sbs_signatures,
                              ordered_dbs_signatures,
                              ordered_id_signatures,
                              ordered_sbs_signatures_cutoffs,
                              ordered_dbs_signatures_cutoffs,
                              ordered_id_signatures_cutoffs,
                              all_types_transcribed_np_array,
                              all_types_untranscribed_np_array,
                              all_types_nontranscribed_np_array,
                              subs_signature_mutation_type_transcribed_np_array,
                              subs_signature_mutation_type_untranscribed_np_array,
                              subs_signature_mutation_type_nontranscribed_np_array,
                              discreet_mode,
                              log_file,
                              verbose)



def searchAllMutations_simbased_chrombased(outputDir,
                                        jobname,
                                        chrLong,
                                        simNum,
                                        samples_of_interest,
                                        sample_based,
                                        all_samples_np_array,
                                        SBS96_mutation_types_np_array,
                                        DBS78_mutation_types_np_array,
                                        ID83_mutation_types_np_array,
                                        all_mutation_types_np_array,
                                        all_mutation_types_np_array_size,
                                        ordered_sbs_signatures_np_array,
                                        ordered_dbs_signatures_np_array,
                                        ordered_id_signatures_np_array,
                                        ordered_sbs_signatures_cutoffs_np_array,
                                        ordered_dbs_signatures_cutoffs_np_array,
                                        ordered_id_signatures_cutoffs_np_array,
                                        discreet_mode,
                                        default_cutoff,
                                        log_file,
                                        verbose):

    number_of_sbs_signatures = ordered_sbs_signatures_np_array.size
    number_of_dbs_signatures = ordered_dbs_signatures_np_array.size
    number_of_id_signatures = ordered_id_signatures_np_array.size

    all_samples_np_array_size = all_samples_np_array.size

    # Initialization for each chr and sim
    # BIDIRECTIONAL --> No index is needed. Will be accumulated in Transcribed and Untranscribed
    # TRANSCRIBED --> 0
    # UNTRANSCRIBED --> 1
    # NONTRANSCRIBED --> 2
    # QUESTIONABLE --> 3
    transcription_strands = [TRANSCRIBED,
                             UNTRANSCRIBED,
                             NONTRANSCRIBED,
                             QUESTIONABLE]
    sample_mutation_type_strand_np_array = np.zeros((all_samples_np_array_size, all_mutation_types_np_array_size , len(transcription_strands)))  # dtype=int
    sample_sbs_signature_mutation_type_strand_np_array = np.zeros((all_samples_np_array_size, number_of_sbs_signatures, SBS96_mutation_types_np_array.size, len(transcription_strands))) # dtype=int
    sample_dbs_signature_mutation_type_strand_np_array = np.zeros((all_samples_np_array_size, number_of_dbs_signatures, DBS78_mutation_types_np_array.size, len(transcription_strands))) # dtype=int
    sample_id_signature_mutation_type_strand_np_array = np.zeros((all_samples_np_array_size, number_of_id_signatures, ID83_mutation_types_np_array.size, len(transcription_strands))) # dtype=int

    chrBased_simBased_subs_df, chrBased_simBased_dinucs_df, chrBased_simBased_indels_df = get_chrBased_simBased_dfs(outputDir, jobname, chrLong, simNum)

    # filter chrbased_df for samples_of_interest
    if samples_of_interest is not None:
        if chrBased_simBased_subs_df is not None:
            chrBased_simBased_subs_df = chrBased_simBased_subs_df[chrBased_simBased_subs_df['Sample'].isin(samples_of_interest)]

        if chrBased_simBased_dinucs_df is not None:
            chrBased_simBased_dinucs_df = chrBased_simBased_dinucs_df[chrBased_simBased_dinucs_df['Sample'].isin(samples_of_interest)]

        if chrBased_simBased_indels_df is not None:
            chrBased_simBased_indels_df = chrBased_simBased_indels_df[chrBased_simBased_indels_df['Sample'].isin(samples_of_interest)]

    return searchAllMutations(chrBased_simBased_subs_df,
                               chrBased_simBased_dinucs_df,
                               chrBased_simBased_indels_df,
                               simNum,
                              SBS96_mutation_types_np_array,
                              DBS78_mutation_types_np_array,
                              ID83_mutation_types_np_array,
                              all_mutation_types_np_array,
                              sample_based,
                               all_samples_np_array,
                               ordered_sbs_signatures_np_array,
                               ordered_dbs_signatures_np_array,
                               ordered_id_signatures_np_array,
                               ordered_sbs_signatures_cutoffs_np_array,
                               ordered_dbs_signatures_cutoffs_np_array,
                               ordered_id_signatures_cutoffs_np_array,
                              sample_mutation_type_strand_np_array,
                              sample_sbs_signature_mutation_type_strand_np_array,
                              sample_dbs_signature_mutation_type_strand_np_array,
                              sample_id_signature_mutation_type_strand_np_array,
                               discreet_mode,
                               default_cutoff,
                               log_file,
                               verbose)


# main function
def transcription_strand_bias_analysis(outputDir,
                                    jobname,
                                    numofSimulations,
                                    samples_of_interest,
                                    job_tuples,
                                    sample_based,
                                    all_samples_np_array,
                                    computationType,
                                    chromNamesList,
                                    ordered_sbs_signatures_np_array,
                                    ordered_dbs_signatures_np_array,
                                    ordered_id_signatures_np_array,
                                    ordered_sbs_signatures_cutoffs_np_array,
                                    ordered_dbs_signatures_cutoffs_np_array,
                                    ordered_id_signatures_cutoffs_np_array,
                                    discreet_mode,
                                    default_cutoff,
                                    parallel_mode,
                                    log_file,
                                    verbose):

    log_out = open(log_file, 'a')
    print('\n#################################################################################', file=log_out)
    print('--- Transcription Strand Asymmetry Analysis starts', file=log_out)
    log_out.close()

    all_mutation_types_np_array = np.concatenate((SBS96_mutation_types_np_array,
                                                  DBS78_mutation_types_np_array,
                                                  ID83_mutation_types_np_array), axis=None)

    all_mutation_types_np_array_size = all_mutation_types_np_array.size

    transcription_strands = [TRANSCRIBED,
                             UNTRANSCRIBED,
                             NONTRANSCRIBED,
                             QUESTIONABLE
                            # BIDIRECTIONAL, will be stored as TRANSCRIBED and UNTRANSCRIBED
                            ]

    six_mutation_types_np_array = np.array([C2A, C2G, C2T, T2A, T2C, T2G])

    all_types_np_array = np.concatenate((SBS6_mutation_types_np_array,
                                         SBS96_mutation_types_np_array,
                                         DBS78_mutation_types_np_array,
                                         ID83_mutation_types_np_array,
                                         ordered_sbs_signatures_np_array,
                                         ordered_dbs_signatures_np_array,
                                         ordered_id_signatures_np_array), axis=None)

    all_types_np_array_size = all_types_np_array.size

    number_of_sbs_signatures = ordered_sbs_signatures_np_array.size
    number_of_dbs_signatures = ordered_dbs_signatures_np_array.size
    number_of_id_signatures = ordered_id_signatures_np_array.size

    if sample_based:
        all_samples_np_array_size = all_samples_np_array.size

        # Initialization for accumulated arrays
        all_sims_sample_mutation_type_strand_np_array = np.zeros((numofSimulations + 1,
                                                                  all_samples_np_array_size,
                                                                  all_mutation_types_np_array_size,
                                                                  len(transcription_strands)))  # dtype=int

        all_sims_sample_sbs_signature_mutation_type_strand_np_array = np.zeros((numofSimulations + 1,
                                                                                all_samples_np_array_size,
                                                                                number_of_sbs_signatures,
                                                                                SBS96_mutation_types_np_array.size,
                                                                                len(transcription_strands)))  # dtype=int

        all_sims_sample_dbs_signature_mutation_type_strand_np_array = np.zeros((numofSimulations + 1,
                                                                                all_samples_np_array_size,
                                                                                number_of_dbs_signatures,
                                                                                DBS78_mutation_types_np_array.size,
                                                                                len(transcription_strands)))  # dtype=int

        all_sims_sample_id_signature_mutation_type_strand_np_array = np.zeros((numofSimulations + 1,
                                                                               all_samples_np_array_size,
                                                                               number_of_id_signatures,
                                                                               ID83_mutation_types_np_array.size,
                                                                               len(transcription_strands)))  # dtype=int

    # Accumulate Numpy Arrays
    # July 27, 2020
    def accumulate_np_arrays(result_tuple):
        sim_num = result_tuple[0]

        # print('MONITOR ACCUMULATE', flush=True)
        if sample_based:
            sample_mutation_type_strand_np_array = result_tuple[1]
            sample_sbs_signature_mutation_type_strand_np_array = result_tuple[2]
            sample_dbs_signature_mutation_type_strand_np_array = result_tuple[3]
            sample_id_signature_mutation_type_strand_np_array = result_tuple[4]

            all_sims_sample_mutation_type_strand_np_array[sim_num] += sample_mutation_type_strand_np_array
            all_sims_sample_sbs_signature_mutation_type_strand_np_array[sim_num] += sample_sbs_signature_mutation_type_strand_np_array
            all_sims_sample_dbs_signature_mutation_type_strand_np_array[sim_num] += sample_dbs_signature_mutation_type_strand_np_array
            all_sims_sample_id_signature_mutation_type_strand_np_array[sim_num] += sample_id_signature_mutation_type_strand_np_array

    sim_nums = range(0, numofSimulations + 1)
    sim_num_chr_tuples = ((sim_num, chrLong) for sim_num in sim_nums for chrLong in chromNamesList)

    if parallel_mode:
        jobs = []

        # Read the chrom based sim based mutations data in the worker process
        if (computationType==USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM):

            numofProcesses = multiprocessing.cpu_count()
            pool = multiprocessing.Pool(processes=numofProcesses)

            for simNum, chrLong in sim_num_chr_tuples:
                jobs.append(pool.apply_async(searchAllMutations_simbased_chrombased,
                                        args=(outputDir,
                                              jobname,
                                              chrLong,
                                              simNum,
                                              samples_of_interest,
                                              sample_based,
                                              all_samples_np_array,
                                              SBS96_mutation_types_np_array,
                                              DBS78_mutation_types_np_array,
                                              ID83_mutation_types_np_array,
                                              all_mutation_types_np_array,
                                              all_mutation_types_np_array_size,
                                              ordered_sbs_signatures_np_array,
                                              ordered_dbs_signatures_np_array,
                                              ordered_id_signatures_np_array,
                                              ordered_sbs_signatures_cutoffs_np_array,
                                              ordered_dbs_signatures_cutoffs_np_array,
                                              ordered_id_signatures_cutoffs_np_array,
                                              discreet_mode,
                                              default_cutoff,
                                              log_file,
                                              verbose,),
                                        callback=accumulate_np_arrays))

            pool.close()
            pool.join()

        elif (computationType == USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT):
            numofProcesses = multiprocessing.cpu_count()
            pool = multiprocessing.Pool(processes=numofProcesses)

            for chrLong, simNum, splitIndex in job_tuples:
                jobs.append(pool.apply_async(searchAllMutations_simbased_chrombased_splitbased,
                                        args=(outputDir,
                                              jobname,
                                              chrLong,
                                              simNum,
                                              splitIndex,
                                              six_mutation_types_np_array,
                                              ordered_sbs_signatures_np_array,
                                              ordered_dbs_signatures_np_array,
                                              ordered_id_signatures_np_array,
                                              ordered_sbs_signatures_cutoffs_np_array,
                                              ordered_dbs_signatures_cutoffs_np_array,
                                              ordered_id_signatures_cutoffs_np_array,
                                              all_types_np_array_size,
                                              discreet_mode,
                                              log_file,
                                              verbose,),
                                        callback=accumulate_np_arrays))
            pool.close()
            pool.join()

    else:
        for simNum, chrLong in sim_num_chr_tuples:
            result_tuple = searchAllMutations_simbased_chrombased(outputDir,
                                               jobname,
                                               chrLong,
                                               simNum,
                                               samples_of_interest,
                                               sample_based,
                                               all_samples_np_array,
                                               SBS96_mutation_types_np_array,
                                               DBS78_mutation_types_np_array,
                                               ID83_mutation_types_np_array,
                                               all_mutation_types_np_array,
                                               all_mutation_types_np_array_size,
                                               ordered_sbs_signatures_np_array,
                                               ordered_dbs_signatures_np_array,
                                               ordered_id_signatures_np_array,
                                               ordered_sbs_signatures_cutoffs_np_array,
                                               ordered_dbs_signatures_cutoffs_np_array,
                                               ordered_id_signatures_cutoffs_np_array,
                                               discreet_mode,
                                               default_cutoff,
                                               log_file,
                                               verbose)

            accumulate_np_arrays(result_tuple)

    #################################################################################################################
    ##########################################      Output starts      ##############################################
    #################################################################################################################
    strand_bias = TRANSCRIPTIONSTRANDBIAS

    transcription_strands = [TRANSCRIBED,
                             UNTRANSCRIBED,
                             NONTRANSCRIBED]

    zipped_arrays = zip([all_sims_sample_sbs_signature_mutation_type_strand_np_array,
                         all_sims_sample_dbs_signature_mutation_type_strand_np_array,
                         all_sims_sample_id_signature_mutation_type_strand_np_array],
                        [SBS96_mutation_types_np_array,
                         DBS78_mutation_types_np_array,
                         ID83_mutation_types_np_array],
                        [ordered_sbs_signatures_np_array,
                         ordered_dbs_signatures_np_array,
                         ordered_id_signatures_np_array])

    write_signature_mutation_type_strand_bias_np_array_as_dataframe(zipped_arrays,
                                                                    strand_bias,
                                                                    transcription_strands,
                                                                    outputDir,
                                                                    jobname)

    write_type_strand_bias_np_array_as_dataframe(all_sims_sample_mutation_type_strand_np_array,
                                                       all_sims_sample_sbs_signature_mutation_type_strand_np_array,
                                                       all_sims_sample_dbs_signature_mutation_type_strand_np_array,
                                                       all_sims_sample_id_signature_mutation_type_strand_np_array,
                                                       SBS6_mutation_types_np_array,
                                                       SBS96_mutation_types_np_array,
                                                       ordered_sbs_signatures_np_array,
                                                       ordered_dbs_signatures_np_array,
                                                       ordered_id_signatures_np_array,
                                                       all_types_np_array,
                                                       strand_bias,
                                                       transcription_strands,
                                                       outputDir,
                                                       jobname)


    if sample_based:
        sample_based_zipped_arrays = zip([all_sims_sample_sbs_signature_mutation_type_strand_np_array,
                                          all_sims_sample_dbs_signature_mutation_type_strand_np_array,
                                          all_sims_sample_id_signature_mutation_type_strand_np_array],
                                         [SBS96_mutation_types_np_array,
                                          DBS78_mutation_types_np_array,
                                          ID83_mutation_types_np_array],
                                         [ordered_sbs_signatures_np_array,
                                          ordered_dbs_signatures_np_array,
                                          ordered_id_signatures_np_array])

        write_sample_signature_mutation_type_strand_np_array_as_dataframe(outputDir,
                                                                          jobname,
                                                                          strand_bias,
                                                                          transcription_strands,
                                                                          all_samples_np_array,
                                                                          sample_based_zipped_arrays)

        write_sample_type_strand_bias_np_array_as_dataframe(outputDir,
                        jobname,
                        strand_bias,
                        transcription_strands,
                        all_samples_np_array,
                        all_types_np_array,
                        SBS6_mutation_types_np_array,
                        SBS96_mutation_types_np_array,
                        DBS78_mutation_types_np_array,
                        ID83_mutation_types_np_array,
                        ordered_sbs_signatures_np_array,
                        ordered_dbs_signatures_np_array,
                        ordered_id_signatures_np_array,
                        all_sims_sample_mutation_type_strand_np_array,
                        all_sims_sample_sbs_signature_mutation_type_strand_np_array,
                        all_sims_sample_dbs_signature_mutation_type_strand_np_array,
                        all_sims_sample_id_signature_mutation_type_strand_np_array)


    #################################################################################################################
    ##########################################      Output ends      ################################################
    #################################################################################################################

    log_out = open(log_file, 'a')
    print('--- Transcription Strand Asymmetry Analysis ends', file=log_out)
    print('#################################################################################\n', file=log_out)
    log_out.close()
