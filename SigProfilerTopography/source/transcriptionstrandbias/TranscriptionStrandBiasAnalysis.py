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

from SigProfilerTopography.source.commons.TopographyCommons import TYPE

from SigProfilerTopography.source.commons.TopographyCommons import SUBS
from SigProfilerTopography.source.commons.TopographyCommons import INDELS
from SigProfilerTopography.source.commons.TopographyCommons import DINUCS

from SigProfilerTopography.source.commons.TopographyCommons import UNTRANSCRIBED_STRAND
from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIBED_STRAND
from SigProfilerTopography.source.commons.TopographyCommons import NONTRANSCRIBED_STRAND
from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIPTIONSTRANDBIAS

from SigProfilerTopography.source.commons.TopographyCommons import write_signature_mutation_type_strand_bias_np_array_as_dataframe
from SigProfilerTopography.source.commons.TopographyCommons import write_type_strand_bias_np_array_as_dataframe
from SigProfilerTopography.source.commons.TopographyCommons import write_sample_based_strand1_strand2_as_dataframe

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

########################################################################
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
                                                                            df_columns):

    indexofTranscriptionStrand = np.where(df_columns == TRANSCRIPTIONSTRAND)[0][0]
    indexofType = np.where(df_columns == TYPE)[0][0]

    mutationTranscriptionStrand = mutation_row[indexofTranscriptionStrand]
    my_type = mutation_row[indexofType]

    subs_signatures_mutation_types_mask_array=subs_signatures_mutation_types_default_zeros_array

    ##########################################
    if (my_type==SUBS):
        #e.g.: C>A
        indexofMutation = np.where(df_columns == MUTATION)[0][0]
        mutationType = mutation_row[indexofMutation]

        six_mutation_types_mask_array= np.where(six_mutation_types_np_array == mutationType, 1, 0)

        probabilities = mutation_row[df_columns_subs_signatures_mask_array]
        threshold_mask_array = np.greater_equal(probabilities, ordered_sbs_signatures_cutoffs)
        # Convert True into 1, and False into 0
        subs_signatures_mask_array = threshold_mask_array.astype(int)

        #Concetanate
        all_types_mask_array= np.concatenate((six_mutation_types_mask_array, subs_signatures_mask_array, dinucs_signatures_default_zeros_array, indels_signatures_default_zeros_array), axis=None)

        #multiply subs_signatures_mask_array times six_mutation_types_mask_array
        # Add one more dimension to subs_signatures_mask_array and six_mutation_types_mask_array
        subs_signatures_mask_array_2d = np.array([subs_signatures_mask_array])
        six_mutation_types_mask_array_2d = np.array([six_mutation_types_mask_array])
        subs_signatures_mutation_types_mask_array = subs_signatures_mask_array_2d.T * six_mutation_types_mask_array_2d
    elif (my_type == DINUCS):
        probabilities = mutation_row[df_columns_dinucs_signatures_mask_array]
        threshold_mask_array = np.greater_equal(probabilities, ordered_dbs_signatures_cutoffs)
        # Convert True into 1, and False into 0
        dinucs_signatures_mask_array = threshold_mask_array.astype(int)

        #Concetanate
        all_types_mask_array= np.concatenate((six_mutation_types_default_zeros_array, subs_signatures_default_zeros_array, dinucs_signatures_mask_array, indels_signatures_default_zeros_array), axis=None)

    elif (my_type == INDELS):
        probabilities = mutation_row[df_columns_indels_signatures_mask_array]
        threshold_mask_array = np.greater_equal(probabilities, ordered_id_signatures_cutoffs)
        # Convert True into 1, and False into 0
        indels_signatures_mask_array = threshold_mask_array.astype(int)

        #Concetanate
        all_types_mask_array= np.concatenate((six_mutation_types_default_zeros_array, subs_signatures_default_zeros_array, dinucs_signatures_default_zeros_array, indels_signatures_mask_array), axis=None)
    ##########################################

    #Values on TranscriptionStrand column
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

########################################################################


########################################################################
#April 24, 2020
#Updated July 29, 2020
def search_all_mutations_on_transcription_strand_array_using_list_comprehension_using_numpy_array(mutation_row,
                                                                            my_type,
                                                                            sample_based,
                                                                            all_samples_np_array,
                                                                            six_mutation_types_np_array,
                                                                            ordered_signatures_cutoffs,
                                                                            df_columns_signatures_mask_array,
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
                                                                            all_samples_all_types_transcribed_np_array,
                                                                            all_samples_all_types_untranscribed_np_array,
                                                                            all_samples_all_types_nontranscribed_np_array,
                                                                            all_samples_subs_signature_mutation_type_transcribed_np_array,
                                                                            all_samples_subs_signature_mutation_type_untranscribed_np_array,
                                                                            all_samples_subs_signature_mutation_type_nontranscribed_np_array,
                                                                            df_columns):

    if sample_based:
        indexofSample = np.where(df_columns == SAMPLE)[0][0]
        mutation_sample = mutation_row[indexofSample]
        sample_index=np.where(all_samples_np_array == mutation_sample)[0][0]

    indexofTranscriptionStrand = np.where(df_columns == TRANSCRIPTIONSTRAND)[0][0]
    mutationTranscriptionStrand = mutation_row[indexofTranscriptionStrand]

    subs_signatures_mutation_types_mask_array=subs_signatures_mutation_types_default_zeros_array

    probabilities = mutation_row[df_columns_signatures_mask_array]
    threshold_mask_array = np.greater_equal(probabilities, ordered_signatures_cutoffs)


    ##########################################
    if (my_type==SUBS):
        #e.g.: C>A
        indexofMutation = np.where(df_columns == MUTATION)[0][0]
        mutationType = mutation_row[indexofMutation]

        six_mutation_types_mask_array= np.where(six_mutation_types_np_array == mutationType, 1, 0)

        # Convert True into 1, and False into 0
        subs_signatures_mask_array = threshold_mask_array.astype(int)

        #Concetanate
        all_types_mask_array= np.concatenate((six_mutation_types_mask_array, subs_signatures_mask_array, dinucs_signatures_default_zeros_array, indels_signatures_default_zeros_array), axis=None)

        #multiply subs_signatures_mask_array times six_mutation_types_mask_array
        # Add one more dimension to subs_signatures_mask_array and six_mutation_types_mask_array
        subs_signatures_mask_array_2d = np.array([subs_signatures_mask_array])
        six_mutation_types_mask_array_2d = np.array([six_mutation_types_mask_array])
        subs_signatures_mutation_types_mask_array = subs_signatures_mask_array_2d.T * six_mutation_types_mask_array_2d
    elif (my_type == DINUCS):
        # Convert True into 1, and False into 0
        dinucs_signatures_mask_array = threshold_mask_array.astype(int)

        #Concetanate
        all_types_mask_array= np.concatenate((six_mutation_types_default_zeros_array, subs_signatures_default_zeros_array, dinucs_signatures_mask_array, indels_signatures_default_zeros_array), axis=None)

    elif (my_type == INDELS):
        # Convert True into 1, and False into 0
        indels_signatures_mask_array = threshold_mask_array.astype(int)

        #Concetanate
        all_types_mask_array= np.concatenate((six_mutation_types_default_zeros_array, subs_signatures_default_zeros_array, dinucs_signatures_default_zeros_array, indels_signatures_mask_array), axis=None)
    ##########################################

    #Values on TranscriptionStrand column
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

    if sample_based:
        if (mutationTranscriptionStrand == 'U'):
            # Jan 21, 2021
            all_samples_all_types_untranscribed_np_array[sample_index] += all_types_mask_array
            all_samples_subs_signature_mutation_type_untranscribed_np_array[sample_index] += subs_signatures_mutation_types_mask_array

        elif (mutationTranscriptionStrand == 'T'):
            # Jan 21, 2021
            all_samples_all_types_transcribed_np_array[sample_index] += all_types_mask_array
            all_samples_subs_signature_mutation_type_transcribed_np_array[sample_index] += subs_signatures_mutation_types_mask_array

        elif (mutationTranscriptionStrand == 'B'):
            # Jan 21, 2021
            all_samples_all_types_untranscribed_np_array[sample_index] += all_types_mask_array
            all_samples_subs_signature_mutation_type_untranscribed_np_array[sample_index] += subs_signatures_mutation_types_mask_array
            all_samples_all_types_transcribed_np_array[sample_index] += all_types_mask_array
            all_samples_subs_signature_mutation_type_transcribed_np_array[sample_index] += subs_signatures_mutation_types_mask_array

        elif (mutationTranscriptionStrand == 'N'):
            # Jan 21, 2021
            all_samples_all_types_nontranscribed_np_array[sample_index] += all_types_mask_array
            all_samples_subs_signature_mutation_type_nontranscribed_np_array[sample_index] += subs_signatures_mutation_types_mask_array


########################################################################


########################################################################
# For df_split
# April 30, 2020
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
            verbose):

    ################################################################################
    if ((chrBased_simBased_combined_df is not None) and (not chrBased_simBased_combined_df.empty)):
        if verbose: print('Worker pid %s searchMutationUsingTranscriptionStrandColumn_simulations_integrated starts %s MB' % (str(os.getpid()), memory_usage()))

        #df_columns numpy array
        df_columns = chrBased_simBased_combined_df.columns.values

        ###############################################################################
        ################################ Initialization ###############################
        ###############################################################################
        df_columns_subs_signatures_mask_array = np.isin(df_columns,ordered_sbs_signatures)
        df_columns_dinucs_signatures_mask_array=np.isin(df_columns,ordered_dbs_signatures)
        df_columns_indels_signatures_mask_array= np.isin(df_columns,ordered_id_signatures)

        six_mutation_types_default_zeros_array= np.zeros(six_mutation_types_np_array.size,dtype=int)
        subs_signatures_default_zeros_array = np.zeros(ordered_sbs_signatures.size,dtype=int)
        dinucs_signatures_default_zeros_array = np.zeros(ordered_dbs_signatures.size,dtype=int)
        indels_signatures_default_zeros_array = np.zeros(ordered_id_signatures.size,dtype=int)
        subs_signatures_mutation_types_default_zeros_array= np.zeros((ordered_sbs_signatures.size,six_mutation_types_np_array.size),dtype=int)
        ###############################################################################
        ################################ Initialization ###############################
        ###############################################################################

        #####################################################################################
        #July 28, 2020
        #Search for each row using Numpy Array
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
                                                                            df_columns) for mutation_row in chrBased_simBased_combined_df.values]
        #####################################################################################


        if verbose: print('Worker pid %s searchMutationUsingTranscriptionStrandColumn_simulations_integrated ends %s MB' % (str(os.getpid()), memory_usage()))
    ################################################################################

    return (sim_num,
            all_types_transcribed_np_array,
            all_types_untranscribed_np_array,
            all_types_nontranscribed_np_array,
            subs_signature_mutation_type_transcribed_np_array,
            subs_signature_mutation_type_untranscribed_np_array,
            subs_signature_mutation_type_nontranscribed_np_array)
########################################################################



########################################################################
# April 30, 2020
def searchAllMutations(chrBased_simBased_subs_df,
            chrBased_simBased_dinucs_df,
            chrBased_simBased_indels_df,
            sim_num,
            sample_based,
            all_samples_np_array,
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
            all_samples_all_types_transcribed_np_array,
            all_samples_all_types_untranscribed_np_array,
            all_samples_all_types_nontranscribed_np_array,
            all_samples_subs_signature_mutation_type_transcribed_np_array,
            all_samples_subs_signature_mutation_type_untranscribed_np_array,
            all_samples_subs_signature_mutation_type_nontranscribed_np_array,
            verbose):

    ###############################################################################
    ################################ Initialization ###############################
    ###############################################################################
    six_mutation_types_default_zeros_array = np.zeros(six_mutation_types_np_array.size, dtype=int)
    subs_signatures_default_zeros_array = np.zeros(ordered_sbs_signatures.size, dtype=int)
    dinucs_signatures_default_zeros_array = np.zeros(ordered_dbs_signatures.size, dtype=int)
    indels_signatures_default_zeros_array = np.zeros(ordered_id_signatures.size, dtype=int)
    subs_signatures_mutation_types_default_zeros_array = np.zeros((ordered_sbs_signatures.size, six_mutation_types_np_array.size), dtype=int)
    ###############################################################################
    ################################ Initialization ###############################
    ###############################################################################

    ################################################################################
    #SUBS
    if ((chrBased_simBased_subs_df is not None) and (not chrBased_simBased_subs_df.empty)):
        if verbose: print('Worker pid %s searchMutationUsingTranscriptionStrandColumn_simulations_integrated starts %s MB' % (str(os.getpid()), memory_usage()))

        #df_columns numpy array
        df_columns = chrBased_simBased_subs_df.columns.values
        df_columns_subs_signatures_mask_array = np.isin(df_columns, ordered_sbs_signatures)

        #Search for each row using Numpy Array
        [search_all_mutations_on_transcription_strand_array_using_list_comprehension_using_numpy_array(mutation_row,
                                                                            SUBS,
                                                                            sample_based,
                                                                            all_samples_np_array,
                                                                            six_mutation_types_np_array,
                                                                            ordered_sbs_signatures_cutoffs,
                                                                            df_columns_subs_signatures_mask_array,
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
                                                                            all_samples_all_types_transcribed_np_array,
                                                                            all_samples_all_types_untranscribed_np_array,
                                                                            all_samples_all_types_nontranscribed_np_array,
                                                                            all_samples_subs_signature_mutation_type_transcribed_np_array,
                                                                            all_samples_subs_signature_mutation_type_untranscribed_np_array,
                                                                            all_samples_subs_signature_mutation_type_nontranscribed_np_array,
                                                                            df_columns) for mutation_row in chrBased_simBased_subs_df.values]


    #DINUCS
    if ((chrBased_simBased_dinucs_df is not None) and (not chrBased_simBased_dinucs_df.empty)):
        if verbose: print('Worker pid %s searchMutationUsingTranscriptionStrandColumn_simulations_integrated starts %s MB' % (str(os.getpid()), memory_usage()))

        # df_columns numpy array
        df_columns = chrBased_simBased_dinucs_df.columns.values
        df_columns_dinucs_signatures_mask_array = np.isin(df_columns, ordered_dbs_signatures)

        # Search for each row using Numpy Array
        [search_all_mutations_on_transcription_strand_array_using_list_comprehension_using_numpy_array(mutation_row,
                                                                                              DINUCS,
                                                                                              sample_based,
                                                                                              all_samples_np_array,
                                                                                              six_mutation_types_np_array,
                                                                                              ordered_dbs_signatures_cutoffs,
                                                                                              df_columns_dinucs_signatures_mask_array,
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
                                                                                              all_samples_all_types_transcribed_np_array,
                                                                                              all_samples_all_types_untranscribed_np_array,
                                                                                              all_samples_all_types_nontranscribed_np_array,
                                                                                              all_samples_subs_signature_mutation_type_transcribed_np_array,
                                                                                              all_samples_subs_signature_mutation_type_untranscribed_np_array,
                                                                                              all_samples_subs_signature_mutation_type_nontranscribed_np_array,
                                                                                              df_columns) for mutation_row in chrBased_simBased_dinucs_df.values]

    # INDELS
    if ((chrBased_simBased_indels_df is not None) and (not chrBased_simBased_indels_df.empty)):
        if verbose: print('Worker pid %s searchMutationUsingTranscriptionStrandColumn_simulations_integrated starts %s MB' % (str(os.getpid()), memory_usage()))

        # df_columns numpy array
        df_columns = chrBased_simBased_indels_df.columns.values
        df_columns_indels_signatures_mask_array = np.isin(df_columns, ordered_id_signatures)

        # Search for each row using Numpy Array
        [search_all_mutations_on_transcription_strand_array_using_list_comprehension_using_numpy_array(mutation_row,
                                                                                              INDELS,
                                                                                              sample_based,
                                                                                              all_samples_np_array,
                                                                                              six_mutation_types_np_array,
                                                                                              ordered_id_signatures_cutoffs,
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
                                                                                              all_samples_all_types_transcribed_np_array,
                                                                                              all_samples_all_types_untranscribed_np_array,
                                                                                              all_samples_all_types_nontranscribed_np_array,
                                                                                              all_samples_subs_signature_mutation_type_transcribed_np_array,
                                                                                              all_samples_subs_signature_mutation_type_untranscribed_np_array,
                                                                                              all_samples_subs_signature_mutation_type_nontranscribed_np_array,
                                                                                              df_columns) for mutation_row in chrBased_simBased_indels_df.values]
    ################################################################################

    return (sim_num,
            all_types_transcribed_np_array,
            all_types_untranscribed_np_array,
            all_types_nontranscribed_np_array,
            subs_signature_mutation_type_transcribed_np_array,
            subs_signature_mutation_type_untranscribed_np_array,
            subs_signature_mutation_type_nontranscribed_np_array,
            all_samples_all_types_transcribed_np_array,
            all_samples_all_types_untranscribed_np_array,
            all_samples_all_types_nontranscribed_np_array,
            all_samples_subs_signature_mutation_type_transcribed_np_array,
            all_samples_subs_signature_mutation_type_untranscribed_np_array,
            all_samples_subs_signature_mutation_type_nontranscribed_np_array)
########################################################################


########################################################################
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
                                           verbose):

    #Initialization
    all_types_transcribed_np_array = np.zeros((all_types_np_array_size), dtype=int)
    all_types_untranscribed_np_array = np.zeros((all_types_np_array_size), dtype=int)
    all_types_nontranscribed_np_array = np.zeros((all_types_np_array_size), dtype=int)
    subs_signature_mutation_type_transcribed_np_array = np.zeros((ordered_sbs_signatures.size, six_mutation_types_np_array.size),dtype=int)
    subs_signature_mutation_type_untranscribed_np_array = np.zeros((ordered_sbs_signatures.size, six_mutation_types_np_array.size),dtype=int)
    subs_signature_mutation_type_nontranscribed_np_array = np.zeros((ordered_sbs_signatures.size, six_mutation_types_np_array.size),dtype=int)

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
                              verbose)
########################################################################

########################################################################
def searchAllMutations_simbased_chrombased(outputDir,
                                        jobname,
                                        chrLong,
                                        simNum,
                                        sample_based,
                                        all_samples_np_array,
                                        six_mutation_types_np_array,
                                        ordered_sbs_signatures,
                                        ordered_dbs_signatures,
                                        ordered_id_signatures,
                                        ordered_sbs_signatures_cutoffs,
                                        ordered_dbs_signatures_cutoffs,
                                        ordered_id_signatures_cutoffs,
                                        all_types_np_array_size,
                                        verbose):

    #Initialization for each chr and sim
    all_types_transcribed_np_array = np.zeros((all_types_np_array_size), dtype=int)
    all_types_untranscribed_np_array = np.zeros((all_types_np_array_size), dtype=int)
    all_types_nontranscribed_np_array = np.zeros((all_types_np_array_size), dtype=int)
    subs_signature_mutation_type_transcribed_np_array = np.zeros((ordered_sbs_signatures.size, six_mutation_types_np_array.size),dtype=int)
    subs_signature_mutation_type_untranscribed_np_array = np.zeros((ordered_sbs_signatures.size, six_mutation_types_np_array.size),dtype=int)
    subs_signature_mutation_type_nontranscribed_np_array = np.zeros((ordered_sbs_signatures.size, six_mutation_types_np_array.size),dtype=int)

    chrBased_simBased_subs_df, chrBased_simBased_dinucs_df, chrBased_simBased_indels_df = get_chrBased_simBased_dfs(outputDir, jobname, chrLong, simNum)

    if sample_based:
        #Jan 21, 2021
        all_samples_np_array_size=all_samples_np_array.size
        all_samples_all_types_transcribed_np_array = np.zeros((all_samples_np_array_size,all_types_np_array_size), dtype=int)
        all_samples_all_types_untranscribed_np_array = np.zeros((all_samples_np_array_size,all_types_np_array_size), dtype=int)
        all_samples_all_types_nontranscribed_np_array = np.zeros((all_samples_np_array_size,all_types_np_array_size), dtype=int)
        all_samples_subs_signature_mutation_type_transcribed_np_array = np.zeros((all_samples_np_array_size,ordered_sbs_signatures.size, six_mutation_types_np_array.size),dtype=int)
        all_samples_subs_signature_mutation_type_untranscribed_np_array = np.zeros((all_samples_np_array_size,ordered_sbs_signatures.size, six_mutation_types_np_array.size),dtype=int)
        all_samples_subs_signature_mutation_type_nontranscribed_np_array = np.zeros((all_samples_np_array_size,ordered_sbs_signatures.size, six_mutation_types_np_array.size),dtype=int)
    else:
        #Jan 21, 2021
        all_samples_all_types_transcribed_np_array = None
        all_samples_all_types_untranscribed_np_array = None
        all_samples_all_types_nontranscribed_np_array = None
        all_samples_subs_signature_mutation_type_transcribed_np_array = None
        all_samples_subs_signature_mutation_type_untranscribed_np_array = None
        all_samples_subs_signature_mutation_type_nontranscribed_np_array = None

    return  searchAllMutations(chrBased_simBased_subs_df,
                               chrBased_simBased_dinucs_df,
                               chrBased_simBased_indels_df,
                               simNum,
                               sample_based,
                               all_samples_np_array,
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
                               all_samples_all_types_transcribed_np_array,
                               all_samples_all_types_untranscribed_np_array,
                               all_samples_all_types_nontranscribed_np_array,
                               all_samples_subs_signature_mutation_type_transcribed_np_array,
                               all_samples_subs_signature_mutation_type_untranscribed_np_array,
                               all_samples_subs_signature_mutation_type_nontranscribed_np_array,
                               verbose)
########################################################################


########################################################################
#main function
def transcriptionStrandBiasAnalysis(outputDir,
                                    jobname,
                                    numofSimulations,
                                    job_tuples,
                                    sample_based,
                                    all_samples_np_array,
                                    computationType,
                                    chromNamesList,
                                    ordered_sbs_signatures,
                                    ordered_dbs_signatures,
                                    ordered_id_signatures,
                                    ordered_sbs_signatures_cutoffs,
                                    ordered_dbs_signatures_cutoffs,
                                    ordered_id_signatures_cutoffs,
                                    verbose):

    print('\n#################################################################################')
    print('--- TranscriptionStrandBias Analysis starts')

    #########################################################################################
    six_mutation_types_np_array = np.array([C2A, C2G, C2T, T2A, T2C, T2G])
    all_types_np_array = np.concatenate((six_mutation_types_np_array,ordered_sbs_signatures,ordered_dbs_signatures,ordered_id_signatures), axis=None)
    all_types_np_array_size = six_mutation_types_np_array.size + ordered_sbs_signatures.size + ordered_dbs_signatures.size + ordered_id_signatures.size
    #########################################################################################

    #########################################################################################
    # Initialization for accumulated arrays
    all_sims_all_types_transcribed_np_array = np.zeros((numofSimulations+1,all_types_np_array_size),dtype=int)
    all_sims_all_types_untranscribed_np_array = np.zeros((numofSimulations+1,all_types_np_array_size),dtype=int)
    all_sims_all_types_nontranscribed_np_array = np.zeros((numofSimulations+1,all_types_np_array_size),dtype=int)
    all_sims_subs_signature_mutation_type_transcribed_np_array = np.zeros((numofSimulations+1,ordered_sbs_signatures.size, six_mutation_types_np_array.size),dtype=int)
    all_sims_subs_signature_mutation_type_untranscribed_np_array= np.zeros((numofSimulations+1,ordered_sbs_signatures.size, six_mutation_types_np_array.size),dtype=int)
    all_sims_subs_signature_mutation_type_nontranscribed_np_array= np.zeros((numofSimulations+1,ordered_sbs_signatures.size, six_mutation_types_np_array.size),dtype=int)
    #########################################################################################

    #########################################################################################
    if sample_based:
        # Initialization for accumulated arrays
        all_samples_np_array_size=all_samples_np_array.size
        all_sims_all_samples_all_types_transcribed_np_array = np.zeros((numofSimulations+1,all_samples_np_array_size, all_types_np_array_size),dtype=int)
        all_sims_all_samples_all_types_untranscribed_np_array = np.zeros((numofSimulations+1,all_samples_np_array_size,all_types_np_array_size),dtype=int)
        all_sims_all_samples_all_types_nontranscribed_np_array = np.zeros((numofSimulations+1,all_samples_np_array_size, all_types_np_array_size),dtype=int)
        all_sims_all_samples_subs_signature_mutation_type_transcribed_np_array = np.zeros((numofSimulations+1,all_samples_np_array_size,ordered_sbs_signatures.size, six_mutation_types_np_array.size),dtype=int)
        all_sims_all_samples_subs_signature_mutation_type_untranscribed_np_array= np.zeros((numofSimulations+1,all_samples_np_array_size,ordered_sbs_signatures.size, six_mutation_types_np_array.size),dtype=int)
        all_sims_all_samples_subs_signature_mutation_type_nontranscribed_np_array= np.zeros((numofSimulations+1,all_samples_np_array_size,ordered_sbs_signatures.size, six_mutation_types_np_array.size),dtype=int)
    #########################################################################################

    #########################################################################################
    # Accumulate Numpy Arrays
    # July 27, 2020
    def accumulate_np_arrays(result_tuple):
        sim_num = result_tuple[0]
        all_types_transcribed_np_array = result_tuple[1]
        all_types_untranscribed_np_array = result_tuple[2]
        all_types_nontranscribed_np_array = result_tuple[3]
        subs_signature_mutation_type_transcribed_np_array = result_tuple[4]
        subs_signature_mutation_type_untranscribed_np_array = result_tuple[5]
        subs_signature_mutation_type_nontranscribed_np_array = result_tuple[6]

        all_sims_all_types_transcribed_np_array[sim_num] += all_types_transcribed_np_array
        all_sims_all_types_untranscribed_np_array[sim_num] += all_types_untranscribed_np_array
        all_sims_all_types_nontranscribed_np_array[sim_num] += all_types_nontranscribed_np_array
        all_sims_subs_signature_mutation_type_transcribed_np_array[sim_num] += subs_signature_mutation_type_transcribed_np_array
        all_sims_subs_signature_mutation_type_untranscribed_np_array[sim_num] += subs_signature_mutation_type_untranscribed_np_array
        all_sims_subs_signature_mutation_type_nontranscribed_np_array[sim_num] += subs_signature_mutation_type_nontranscribed_np_array

        # print('MONITOR ACCUMULATE', flush=True)
        if sample_based:
            #Jan 21, 2021
            all_samples_all_types_transcribed_np_array = result_tuple[7]
            all_samples_all_types_untranscribed_np_array = result_tuple[8]
            all_samples_all_types_nontranscribed_np_array = result_tuple[9]
            all_samples_subs_signature_mutation_type_transcribed_np_array = result_tuple[10]
            all_samples_subs_signature_mutation_type_untranscribed_np_array = result_tuple[11]
            all_samples_subs_signature_mutation_type_nontranscribed_np_array = result_tuple[12]

            all_sims_all_samples_all_types_transcribed_np_array[sim_num] += all_samples_all_types_transcribed_np_array
            all_sims_all_samples_all_types_untranscribed_np_array[sim_num] += all_samples_all_types_untranscribed_np_array
            all_sims_all_samples_all_types_nontranscribed_np_array[sim_num] += all_samples_all_types_nontranscribed_np_array
            all_sims_all_samples_subs_signature_mutation_type_transcribed_np_array[sim_num] += all_samples_subs_signature_mutation_type_transcribed_np_array
            all_sims_all_samples_subs_signature_mutation_type_untranscribed_np_array[sim_num] += all_samples_subs_signature_mutation_type_untranscribed_np_array
            all_sims_all_samples_subs_signature_mutation_type_nontranscribed_np_array[sim_num] += all_samples_subs_signature_mutation_type_nontranscribed_np_array
    #########################################################################################

    jobs = []

    ###############################################################################
    # April 30, 2020
    # Read the chrom based sim based mutations data in the worker process
    if (computationType==USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM):
        sim_nums = range(0, numofSimulations + 1)
        sim_num_chr_tuples = ((sim_num, chrLong) for sim_num in sim_nums for chrLong in chromNamesList)

        numofProcesses = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=numofProcesses)

        for simNum, chrLong in sim_num_chr_tuples:
            jobs.append(pool.apply_async(searchAllMutations_simbased_chrombased,
                                    args=(outputDir,
                                          jobname,
                                          chrLong,
                                          simNum,
                                          sample_based,
                                          all_samples_np_array,
                                          six_mutation_types_np_array,
                                          ordered_sbs_signatures,
                                          ordered_dbs_signatures,
                                          ordered_id_signatures,
                                          ordered_sbs_signatures_cutoffs,
                                          ordered_dbs_signatures_cutoffs,
                                          ordered_id_signatures_cutoffs,
                                          all_types_np_array_size,
                                          verbose,),
                                    callback=accumulate_np_arrays))

        pool.close()
        pool.join()
    ###############################################################################

    ###############################################################################
    elif (computationType==USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT):
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
                                          ordered_sbs_signatures,
                                          ordered_dbs_signatures,
                                          ordered_id_signatures,
                                          ordered_sbs_signatures_cutoffs,
                                          ordered_dbs_signatures_cutoffs,
                                          ordered_id_signatures_cutoffs,
                                          all_types_np_array_size,
                                          verbose,),
                                    callback=accumulate_np_arrays))
        pool.close()
        pool.join()
    ###############################################################################

    #################################################################################################################
    ##########################################      Output starts      ##############################################
    #################################################################################################################
    strand_bias = TRANSCRIPTIONSTRANDBIAS

    transcription_strands = [TRANSCRIBED_STRAND,
                             UNTRANSCRIBED_STRAND,
                             NONTRANSCRIBED_STRAND]

    all_sims_subs_signature_mutation_type_strand_np_arrays_list =[all_sims_subs_signature_mutation_type_transcribed_np_array,
                     all_sims_subs_signature_mutation_type_untranscribed_np_array,
                     all_sims_subs_signature_mutation_type_nontranscribed_np_array]


    write_signature_mutation_type_strand_bias_np_array_as_dataframe(all_sims_subs_signature_mutation_type_strand_np_arrays_list,
                                                                    six_mutation_types_np_array,
                                                                    ordered_sbs_signatures,
                                                                    strand_bias,
                                                                    transcription_strands,
                                                                    outputDir,
                                                                    jobname)

    all_sims_all_types_strand_np_arrays_list =[all_sims_all_types_transcribed_np_array,
                                            all_sims_all_types_untranscribed_np_array,
                                            all_sims_all_types_nontranscribed_np_array]


    write_type_strand_bias_np_array_as_dataframe(all_sims_all_types_strand_np_arrays_list,
                                                all_types_np_array,
                                                strand_bias,
                                                transcription_strands,
                                                outputDir,
                                                jobname)

    if sample_based:
        write_sample_based_strand1_strand2_as_dataframe(outputDir,
                                                        jobname,
                                                        numofSimulations,
                                                        strand_bias,
                                                        all_samples_np_array,
                                                        all_types_np_array,
                                                        all_sims_all_samples_all_types_transcribed_np_array,
                                                        all_sims_all_samples_all_types_untranscribed_np_array)

    #################################################################################################################
    ##########################################      Output ends      ################################################
    #################################################################################################################

    print('--- TranscriptionStrandBias Analysis ends')
    print('#################################################################################\n')

########################################################################
