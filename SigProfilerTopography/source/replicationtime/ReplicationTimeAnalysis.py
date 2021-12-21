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
from SigProfilerTopography.source.commons.TopographyCommons import  decideFileType

from SigProfilerTopography.source.commons.TopographyCommons import  get_chrBased_simBased_combined_df
from SigProfilerTopography.source.commons.TopographyCommons import  get_chrBased_simBased_dfs

# June 6, 2021
mutations_in_early_replicating_regions_df = pd.DataFrame()
mutations_in_late_replicating_regions_df = pd.DataFrame()

##################################################################
# Please note that this dictionary is copied from .../SigProfilerMatrixGenerator/SigProfilerMatrixGenerator/scripts/SigProfilerMatrixGeneratorFunc.py line 101
# If this dictionary is updated in SigProfilerMatrixGeneratorFunc.py, it has to be updated in ReplicationTimeAnalysis.py
# Provides the reference file conversion from binary to base information
tsb_ref = {0: ['N', 'A'], 1: ['N', 'C'], 2: ['N', 'G'], 3: ['N', 'T'],
           4: ['T', 'A'], 5: ['T', 'C'], 6: ['T', 'G'], 7: ['T', 'T'],
           8: ['U', 'A'], 9: ['U', 'C'], 10: ['U', 'G'], 11: ['U', 'T'],
           12: ['B', 'A'], 13: ['B', 'C'], 14: ['B', 'G'], 15: ['B', 'T'],
           16: ['N', 'N'], 17: ['T', 'N'], 18: ['U', 'N'], 19: ['B', 'N']}
##################################################################


##################################################################
#Higher the replication time signal earlier the replication is
#Regions with high values indicate domains of early replication where initiaion occurs earlier in S-phase or early in a higher proportion of cells.
def readRepliSeqTimeData(genome,chromNamesList,repliseqDataFilename,matrix_generator_path,verbose):

    ###################################################################
    ############### Read RepliSeq Time data starts ####################
    ###################################################################

    isFileTypeBEDGRAPH = decideFileType(repliseqDataFilename)

    if isFileTypeBEDGRAPH:
        column_names = [CHROM, START, END, SIGNAL]
        replication_time_interval_version_df = pd.read_csv(repliseqDataFilename, sep='\t', header=None, comment='#', names=column_names,dtype={CHROM: 'category', START: np.int32, END: np.int32, SIGNAL: np.float32})
    else:
        #JAN 7, 2020
        replication_time_interval_version_df = readWig_with_fixedStep_variableStep(repliseqDataFilename)

    chrNamesInReplicationTimeDataArray = replication_time_interval_version_df[CHROM].unique()
    print('Before --- Chromosome names in replication time signal data: %s replication_time_interval_version_df.shape(%d,%d)\n' %(chrNamesInReplicationTimeDataArray,replication_time_interval_version_df.shape[0],replication_time_interval_version_df.shape[1]))

    #Remove rows with chromosomes that are not in chromNamesList
    replication_time_interval_version_df=replication_time_interval_version_df[replication_time_interval_version_df[CHROM].isin(chromNamesList)]

    chrNamesInReplicationTimeDataArray = replication_time_interval_version_df[CHROM].unique()
    print('After considering only chromosomes in chromNamesList --- Chromosome names in replication time signal data: %s replication_time_interval_version_df.shape(%d,%d)\n' %(chrNamesInReplicationTimeDataArray,replication_time_interval_version_df.shape[0],replication_time_interval_version_df.shape[1]))

    #Augment wavelet_processed_df with numberofAttributableBases
    wavelet_processed_augmented_df = augment(genome,replication_time_interval_version_df,matrix_generator_path,verbose)

    #Return 10 deciles: the first decile is the earliest one and the tenth decile is the latest one

    #Sort in descending order
    #Higher the replication time signal earlier the replication is
    wavelet_processed_augmented_df.sort_values(SIGNAL, ascending=False, inplace=True)


    # print('############ after sort wavelet_processed_augmented_df ###################')
    # print(wavelet_processed_augmented_df.head())
    # print('############ after sort wavelet_processed_augmented_df ###################')

    #Split wavelet_processed_augmented_df into 10 deciles
    deciles_df_list = np.array_split(wavelet_processed_augmented_df,10)
    # print('Number of decile:%d' %len(deciles))
    #deciles is a list and each decile is a dataframe <class 'pandas.core.frame.DataFrame'>
    #The first decile is the earliest one
    #The last decile is the latest one
    # print('type(deciles_df_list):%s' %type(deciles_df_list))
    # type(deciles_df_list) --> <class 'list'>

    # ############################################################
    # #For information
    # totalNumberofIntervals = 0
    # for decileIndex, decile_df in enumerate(deciles_df_list,1):
    #     # print('decileIndex: %d' %decileIndex)
    #     # print('type(decile_df)')
    #     # print(type(decile_df))
    #     # type(decile_df) --> <class 'pandas.core.frame.DataFrame'>
    #     totalNumberofIntervals += len(decile_df)
    # print('totalNumberofIntervals in all deciles : %d' %totalNumberofIntervals)
    # ############################################################

    return chrNamesInReplicationTimeDataArray, deciles_df_list
    ###################################################################
    ############### Read RepliSeq Time data ends ######################
    ###################################################################

##################################################################

##################################################################
#Jan 13, 2020
def getNumberofAttributableBasesUsingMatrixGeneratorGenome(wavelet_row,chrom_string):
    start =wavelet_row[1]
    end = wavelet_row[2]

    #old way
    #In my code ends are inclusive
    #twobitreader uses ends exlusive
    # seq_old_way = chrBasedGenome.get_slice(start, end+1)

    seq = ''
    for i in range(start,end+1,1):
        seq += tsb_ref[chrom_string[i - 1]][1]

    numofAttributableBases = seq.count('A') + seq.count('T') + seq.count('G') + seq.count('C') + seq.count('a') + seq.count('t') + seq.count('g') + seq.count('c')
    # print('######### debug starts ##############')
    # print(wavelet_row)
    # print('len(seq):%d' %len(seq))
    # print('numofAttributableBases:%d' %numofAttributableBases)
    # print('##########  debug ends #############')
    return numofAttributableBases
##################################################################


##################################################################
def addNumofAttributableBasesColumnForApplyAsync(chrLong,chrBased_wavelet_processed_df_group,chrbased_file_path,verbose):
    if verbose: print('\tVerbose Worker pid %s %s Before Augment with number of attributable bases memory_usage %.2f MB' % (str(os.getpid()), chrLong, memory_usage()))

    # 1st way Slower. Not tested in terms fo final outcome/ results.
    # chrom_string = np.memmap(chrbased_file_path, dtype=np.byte, mode='r')

    # 2nd way Faster than 1st way
    with open(chrbased_file_path, "rb") as f2:
        chrom_string = f2.read()

    resulting_df = chrBased_wavelet_processed_df_group.apply(getNumberofAttributableBasesUsingMatrixGeneratorGenome,chrom_string=chrom_string, axis= 1)

    if (len(chrBased_wavelet_processed_df_group)!=len(resulting_df)):
        print('There is a situation: len(chrBased_wavelet_processed_df_group) is not equal  to len(resulting_df)')

    chrBased_wavelet_processed_df_group[NUMOFBASES] = resulting_df

    if verbose: print('\tVerbose Worker pid %s %s After Augment with number of attributable bases memory_usage %.2f MB' % (str(os.getpid()),chrLong, memory_usage()))

    return (chrLong,chrBased_wavelet_processed_df_group)
##################################################################


##################################################################
# August 1, 2020
# Using numpy array
# Main engine function
def search_for_each_mutation_using_list_comprehension_using_numpy_array(mutation_row,
                                                                        mutation_type,
                                                                        chrBasedReplicationTimeDataArrayWithDecileIndex,
                                                                        ordered_signatures_cutoffs,
                                                                        signatures_mask_array,
                                                                        signature_decile_index_accumulated_np_array,
                                                                        is_discreet,
                                                                        df_columns):
    # df_columns: numpy array
    indexofStart = np.where(df_columns == START) [0][0]
    start = mutation_row[indexofStart]

    # For single point mutations start and end are the same, therefore we need to add 1 to mutation_row[END]
    # end is exclusive for subs, indels and dinucs provided by readChrBased methods

    #Simulation Number
    index_of_sim_num = np.where(df_columns == 'Simulation_Number')[0][0]
    sim_num = mutation_row[index_of_sim_num]

    if mutation_type==SUBS:
        end = start + 1
    elif mutation_type==DINUCS:
        end = start + 2
    elif mutation_type==INDELS:
        indexofLength = np.where(df_columns == LENGTH)
        length=mutation_row[indexofLength]
        end = start + int(length)

    slicedArray = chrBasedReplicationTimeDataArrayWithDecileIndex[start:end]

    # np.nonzero returns the indices of the elements that are non-zero.
    # np.unique finds the unique elements of an array returns ndarray the sorted unique values.
    # np.nditer efficient multi-dimensional iterator object to iterate over arrays.
    # We aim to get rid of zero if any exists in slicedArray.
    uniqueIndexesArray = np.unique(slicedArray[np.nonzero(slicedArray)])

    decile_index_array=np.zeros((10), dtype=int)

    #decile 10  will be accumulated in index 9
    #decile 1 will be accumulated in index 0
    #Therefore uniqueIndexesArray minus 1
    if (uniqueIndexesArray.size>0):
        uniqueIndexesArray -= 1
        decile_index_array[uniqueIndexesArray] = 1

        probabilities = mutation_row[signatures_mask_array]

        # # old code starts
        # threshold_mask_array = np.greater_equal(probabilities, ordered_signatures_cutoffs)
        #
        # # Convert True into 1, and False into 0
        # mask_array = threshold_mask_array.astype(int)
        # if  mutation_type == INDELS:
        #     if (length >= 3):
        #         # Order is important
        #         # MICROHOMOLOGY --> 1
        #         # REPEAT --> 0
        #         # AGGREGATEDINDELS --> 1
        #         mask_array = np.append(mask_array, [1, 0, 1])
        #     else:
        #         # Order is important
        #         # MICROHOMOLOGY --> 0
        #         # REPEAT --> 1
        #         # AGGREGATEDINDELS --> 1
        #         mask_array = np.append(mask_array, [0, 1, 1])
        #
        # else:
        #     # Add 1 for the aggregated analysis to the mask array
        #     # For SUBS and DINUCS
        #     mask_array = np.append(mask_array, 1)
        # # old code ends

        # new code starts
        if is_discreet:
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
        # new code ends

        # Add 1 more dimension to the arrays
        decile_index_array_1x10 = np.expand_dims(decile_index_array, axis=0)
        mask_array_1xnumofsignatures = np.expand_dims(mask_array, axis=0)

        signatures_decile_index_np_array = mask_array_1xnumofsignatures.T * decile_index_array_1x10
        signature_decile_index_accumulated_np_array += signatures_decile_index_np_array
##################################################################



##################################################################
# For df_split
# August 1, 2020
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
                                                                        sample_based,
                                                                        is_disceet,
                                                                        verbose,
                                                                        df_columns):

    # df_columns: numpy array
    indexofType = np.where(df_columns == TYPE)[0][0]
    indexofSample = np.where(df_columns == SAMPLE)[0][0]
    indexofStart = np.where(df_columns == START) [0][0]
    indexofSimulationNumber = np.where(df_columns == SIMULATION_NUMBER)[0][0]

    ###########################################
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
    ###########################################

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

        if is_disceet:
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

##################################################################



##################################################################
# DEC 18, 2020
def searchforAllMutations_using_numpy_array(sim_num,
                                            chrBased_simBased_subs_df,
                                            chrBased_simBased_dinucs_df,
                                            chrBased_simBased_indels_df,
                                            chrBasedReplicationTimeDataArrayWithDecileIndex,
                                            ordered_all_sbs_signatures_array,
                                            ordered_all_dbs_signatures_array,
                                            ordered_all_id_signatures_array,
                                            ordered_sbs_signatures,
                                            ordered_dbs_signatures,
                                            ordered_id_signatures,
                                            ordered_sbs_signatures_cutoffs,
                                            ordered_dbs_signatures_cutoffs,
                                            ordered_id_signatures_cutoffs,
                                            is_discreet,
                                            verbose):

    if is_discreet:
        number_of_sbs_signatures = ordered_sbs_signatures.size
        number_of_dbs_signatures = ordered_dbs_signatures.size
        number_of_id_signatures = ordered_id_signatures.size
    else:
        number_of_sbs_signatures = ordered_all_sbs_signatures_array.size
        number_of_dbs_signatures = ordered_all_dbs_signatures_array.size
        number_of_id_signatures = ordered_all_id_signatures_array.size


    # Add one more row for the Aggregated analysis, there are 10 deciles
    # Add three more rows for the Microhomology, Repeat Mediated and Aggregated analysis, there are 10 deciles
    subs_signature_decile_index_accumulated_np_array = np.zeros((number_of_sbs_signatures + 1, 10), dtype=float) # legacy int
    dinucs_signature_decile_index_accumulated_np_array = np.zeros((number_of_dbs_signatures + 1, 10), dtype=float) # legacy int
    indels_signature_decile_index_accumulated_np_array = np.zeros((number_of_id_signatures + 3, 10), dtype=float) # # legacy int

    #SUBS
    if ((chrBased_simBased_subs_df is not None) and (not chrBased_simBased_subs_df.empty)):
        df_columns = chrBased_simBased_subs_df.columns.values

        if is_discreet:
            df_columns_subs_signatures_mask_array = np.isin(df_columns, ordered_sbs_signatures)
        else:
            df_columns_subs_signatures_mask_array = np.isin(df_columns, ordered_all_sbs_signatures_array)

        [search_for_each_mutation_using_list_comprehension_using_numpy_array(mutation_row,
                                                                             SUBS,
                                                                            chrBasedReplicationTimeDataArrayWithDecileIndex,
                                                                            ordered_sbs_signatures_cutoffs,
                                                                            df_columns_subs_signatures_mask_array,
                                                                            subs_signature_decile_index_accumulated_np_array,
                                                                            is_discreet,
                                                                            df_columns) for mutation_row in chrBased_simBased_subs_df.values]

    #DINUCS
    if ((chrBased_simBased_dinucs_df is not None) and (not chrBased_simBased_dinucs_df.empty)):
        df_columns = chrBased_simBased_dinucs_df.columns.values

        if is_discreet:
            df_columns_dinucs_signatures_mask_array = np.isin(df_columns, ordered_dbs_signatures)
        else:
            df_columns_dinucs_signatures_mask_array = np.isin(df_columns, ordered_all_dbs_signatures_array)

        [search_for_each_mutation_using_list_comprehension_using_numpy_array(mutation_row,
                                                                             DINUCS,
                                                                            chrBasedReplicationTimeDataArrayWithDecileIndex,
                                                                            ordered_dbs_signatures_cutoffs,
                                                                            df_columns_dinucs_signatures_mask_array,
                                                                            dinucs_signature_decile_index_accumulated_np_array,
                                                                            is_discreet,
                                                                            df_columns) for mutation_row in chrBased_simBased_dinucs_df.values]

    #INDELS
    if ((chrBased_simBased_indels_df is not None) and (not chrBased_simBased_indels_df.empty)):
        df_columns = chrBased_simBased_indels_df.columns.values

        if is_discreet:
            df_columns_indels_signatures_mask_array = np.isin(df_columns, ordered_id_signatures)
        else:
            df_columns_indels_signatures_mask_array = np.isin(df_columns, ordered_all_id_signatures_array)


        [search_for_each_mutation_using_list_comprehension_using_numpy_array(mutation_row,
                                                                             INDELS,
                                                                            chrBasedReplicationTimeDataArrayWithDecileIndex,
                                                                            ordered_id_signatures_cutoffs,
                                                                            df_columns_indels_signatures_mask_array,
                                                                            indels_signature_decile_index_accumulated_np_array,
                                                                            is_discreet,
                                                                            df_columns) for mutation_row in chrBased_simBased_indels_df.values]

    return sim_num, \
           subs_signature_decile_index_accumulated_np_array, \
           dinucs_signature_decile_index_accumulated_np_array, \
           indels_signature_decile_index_accumulated_np_array
##################################################################


##################################################################
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
                                            sample_based,
                                            is_disceet,
                                            verbose):

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
                                                                            sample_based,
                                                                            is_disceet,
                                                                            verbose,
                                                                            df_columns) for mutation_row in chrBased_simBased_combined_df.values]

    return sim_num, subs_signature_decile_index_accumulated_np_array, dinucs_signature_decile_index_accumulated_np_array, indels_signature_decile_index_accumulated_np_array
##################################################################


##################################################################
#Please notice that replication time data are not overlappig data therefore only setting one decileIndex will be correct.
# e.g.:
# start end
#   10 1009
# 1010 2009
def fillArray(chrBased_replicationtimedata_row,chrBasedDecileIndexArray,decileIndex):
    start= chrBased_replicationtimedata_row[START]
    end = chrBased_replicationtimedata_row[END] + 1
    chrBasedDecileIndexArray[start:end] = decileIndex
##################################################################

##################################################################
#Explanation of this function
#This function is called for each chromosome
#For each chromosome, a numpy array is filled.
#If there is a chromosome locus with a decile index 8 let's say, in the array that locus is filled with 8.
#Decile index can be between 1-10.
#In this function, each available interval with an index is filled in the corresponding numpy array with chrom size
def fillChrBasedReplicationTimeNPArray(chrLong,chromSize,chrBased_grouped_decile_df_list):
    #We can set the starting index as 1 in builtin function enumerate
    #First chrBased_grouped_decile has index of 1
    #Last chrBased_grouped_decile has index of 10

    # int8	Byte (-128 to 127)
    chrBasedReplicationTimeDataArrayWithDecileIndex = np.zeros(chromSize, dtype=np.int8)

    #First decileIndex is 1, last decile index is 10.
    for decileIndex, chrBased_grouped_decile_df in enumerate(chrBased_grouped_decile_df_list,1):
        # Solution to keyError
        for name, chrBased_replicationtimedata_df in chrBased_grouped_decile_df:
            if (chrLong==name) and (chrBased_replicationtimedata_df is not None) and  (not chrBased_replicationtimedata_df.empty):
                #what is chrBased_decile's type? DataFrame
                chrBased_replicationtimedata_df.apply(fillArray,chrBasedDecileIndexArray=chrBasedReplicationTimeDataArrayWithDecileIndex,decileIndex=decileIndex, axis=1)

        # # Please note that although key exists if dataframe to be returned is empty dataframe it gives KeyError
        # if chrLong in chrBased_grouped_decile_df.groups.keys():
        #     if chrLong in chrBased_grouped_decile_df.groups:
        #         chrBased_replicationtimedata_df = chrBased_grouped_decile_df.get_group(chrLong)
        #         print('DEBUG %s decileIndex:%d chrBased_replicationtimedata_df.shape(%d,%d)' %(chrLong,decileIndex,chrBased_replicationtimedata_df.shape[0],chrBased_replicationtimedata_df.shape[1]))
        #     else:
        #         chrBased_replicationtimedata_df=None
        #     #what is chrBased_decile's type? DataFrame
        #     if ((chrBased_replicationtimedata_df is not None) and (not chrBased_replicationtimedata_df.empty)):
        #         chrBased_replicationtimedata_df.apply(fillArray,chrBasedDecileIndexArray=chrBasedReplicationTimeDataArrayWithDecileIndex,decileIndex=decileIndex, axis=1)

    return chrBasedReplicationTimeDataArrayWithDecileIndex
##################################################################


##################################################################
# August 2, 2020
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
                                       sample_based,
                                       is_discreet,
                                       verbose):

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
                                                                                                sample_based,
                                                                                                is_discreet,
                                                                                                verbose)

##################################################################

##################################################################
# August 1, 2020
# using numpy arrays
def combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray_simbased_chrbased_using_numpy_array(outputDir,
                                                                                                           jobname,
                                                                                                           chrLong,
                                                                                                           chromSize,
                                                                                                           sim_num,
                                                                                                           chrBased_grouped_decile_df_list,
                                                                                                           ordered_all_sbs_signatures_array,
                                                                                                           ordered_all_dbs_signatures_array,
                                                                                                           ordered_all_id_signatures_array,
                                                                                                           ordered_sbs_signatures,
                                                                                                           ordered_dbs_signatures,
                                                                                                           ordered_id_signatures,
                                                                                                           ordered_sbs_signatures_cutoffs,
                                                                                                           ordered_dbs_signatures_cutoffs,
                                                                                                           ordered_id_signatures_cutoffs,
                                                                                                           is_discreet,
                                                                                                           verbose):

    #DEC18, 2020 To reduce memory usage
    chrBased_simBased_subs_df, chrBased_simBased_dinucs_df, chrBased_simBased_indels_df = get_chrBased_simBased_dfs(outputDir, jobname, chrLong, sim_num)

    return combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray_using_numpy_array(chrLong,
                                                                                                chromSize,
                                                                                                sim_num,
                                                                                                chrBased_grouped_decile_df_list,
                                                                                                chrBased_simBased_subs_df,
                                                                                                chrBased_simBased_dinucs_df,
                                                                                                chrBased_simBased_indels_df,
                                                                                                ordered_all_sbs_signatures_array,
                                                                                                ordered_all_dbs_signatures_array,
                                                                                                ordered_all_id_signatures_array,
                                                                                                ordered_sbs_signatures,
                                                                                                ordered_dbs_signatures,
                                                                                                ordered_id_signatures,
                                                                                                ordered_sbs_signatures_cutoffs,
                                                                                                ordered_dbs_signatures_cutoffs,
                                                                                                ordered_id_signatures_cutoffs,
                                                                                                is_discreet,
                                                                                                verbose)

##################################################################


##################################################################
#Dec 18, 2020
# using numpy arrays
def combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray_using_numpy_array(chrLong,
                                                                                        chromSize,
                                                                                        sim_num,
                                                                                        chrBased_grouped_decile_df_list,
                                                                                        chrBased_simBased_subs_df,
                                                                                        chrBased_simBased_dinucs_df,
                                                                                        chrBased_simBased_indels_df,
                                                                                        ordered_all_sbs_signatures_array,
                                                                                        ordered_all_dbs_signatures_array,
                                                                                        ordered_all_id_signatures_array,
                                                                                        ordered_sbs_signatures,
                                                                                        ordered_dbs_signatures,
                                                                                        ordered_id_signatures,
                                                                                        ordered_sbs_signatures_cutoffs,
                                                                                        ordered_dbs_signatures_cutoffs,
                                                                                        ordered_id_signatures_cutoffs,
                                                                                        is_discreet,
                                                                                        verbose):

    #Fill replication time numpy array
    chrBasedReplicationTimeDataArrayWithDecileIndex = fillChrBasedReplicationTimeNPArray(chrLong, chromSize, chrBased_grouped_decile_df_list)

    return searchforAllMutations_using_numpy_array(sim_num,
                                                   chrBased_simBased_subs_df,
                                                   chrBased_simBased_dinucs_df,
                                                   chrBased_simBased_indels_df,
                                                   chrBasedReplicationTimeDataArrayWithDecileIndex,
                                                   ordered_all_sbs_signatures_array,
                                                   ordered_all_dbs_signatures_array,
                                                   ordered_all_id_signatures_array,
                                                   ordered_sbs_signatures,
                                                   ordered_dbs_signatures,
                                                   ordered_id_signatures,
                                                   ordered_sbs_signatures_cutoffs,
                                                   ordered_dbs_signatures_cutoffs,
                                                   ordered_id_signatures_cutoffs,
                                                   is_discreet,
                                                   verbose)
##################################################################


##################################################################
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
                                                                                        sample_based,
                                                                                        is_discreet,
                                                                                        verbose):

    #Fill replication time numpy array
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
                                                   sample_based,
                                                   is_discreet,
                                                   verbose)
##################################################################


##################################################################
def getNormalizedMutationDensityList(mutationDensityDict):
    densities = mutationDensityDict.values()
    maxDensity = max(densities)
    if maxDensity>0:
        normalizedMutationDensities = [x/maxDensity for x in densities]
    else:
        normalizedMutationDensities = densities

    return normalizedMutationDensities
##################################################################



##################################################################
#March 22, 2019 starts
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
#March 22, 2019 starts
##################################################################

##################################################################
# Feb5, 2021
def getNumberofAttributableBases(decile_df_list):
    numberofAttributableBasesList = []

    #Modifiled as enumerate(deciles,1) formerly it was enumerate(deciles,0)
    for i,decile_df in enumerate(decile_df_list,1):
        numofAttBases = decile_df[NUMOFBASES].sum()
        numberofAttributableBasesList.append(numofAttBases)

    return numberofAttributableBasesList
##################################################################

##################################################################
# Feb5, 2021
def getMutationDensityDictUsingNumpyArray(decile_df_list,decile_counts_np_array):
    numberofMutations = 0
    decileBasedMutationDensityDict = {}
    numberofMutationsList = []

    #Modifiled as enumerate(deciles,1) formerly it was enumerate(deciles,0)
    for i,decile_df in enumerate(decile_df_list,1):
        count = decile_counts_np_array[i-1]
        numberofMutations += count
        numofAttBases = decile_df[NUMOFBASES].sum()
        mutationDensity=float(count)/numofAttBases
        decileBasedMutationDensityDict[i] = mutationDensity
        numberofMutationsList.append(count)

        # decileBasedMutationDensityDict[i] = 0
        # numberofMutationsList.append(0)
        # numberofAttributableBasesList.append(0)
    return numberofMutations, decileBasedMutationDensityDict, numberofMutationsList
##################################################################

##################################################################
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
##################################################################


##################################################################
# August 2, 2020
# Using numpy array
def calculateCountsForMutationsFillingReplicationTimeNPArrayRuntime_using_numpy_array(computationType,
        outputDir,
        jobname,
        numofSimulations,
        job_tuples,
        chromSizesDict,
        chromNamesList,
        chrBased_grouped_decile_df_list,
        ordered_all_sbs_signatures_array,
        ordered_all_dbs_signatures_array,
        ordered_all_id_signatures_array,
        ordered_sbs_signatures,
        ordered_dbs_signatures,
        ordered_id_signatures,
        ordered_sbs_signatures_cutoffs,
        ordered_dbs_signatures_cutoffs,
        ordered_id_signatures_cutoffs,
        sample_based,
        is_discreet,
        verbose):

    if is_discreet:
        number_of_sbs_signatures = ordered_sbs_signatures.size
        number_of_dbs_signatures = ordered_dbs_signatures.size
        number_of_id_signatures = ordered_id_signatures.size
    else:
        number_of_sbs_signatures = ordered_all_sbs_signatures_array.size
        number_of_dbs_signatures = ordered_all_dbs_signatures_array.size
        number_of_id_signatures = ordered_all_id_signatures_array.size

    ################################
    # +1 for Aggregated
    all_sims_subs_signature_decile_index_accumulated_np_array=np.zeros((numofSimulations+1, number_of_sbs_signatures+1, 10), dtype=float) # int
    all_sims_dinucs_signature_decile_index_accumulated_np_array=np.zeros((numofSimulations+1, number_of_dbs_signatures+1, 10), dtype=float) # int
    # +3 for Microhomology, Repeat, Aggregated
    all_sims_indels_signature_decile_index_accumulated_np_array=np.zeros((numofSimulations+1, number_of_id_signatures+3, 10), dtype=float) # int
    ################################

    #########################################################################################
    def accumulate_np_arrays(result_tuple):
        sim_num = result_tuple[0]
        subs_signature_decile_index_accumulated_np_array = result_tuple[1]
        dinucs_signature_decile_index_accumulated_np_array = result_tuple[2]
        indels_signature_decile_index_accumulated_np_array = result_tuple[3]

        # print('MONITOR ACCUMULATE', flush=True)

        all_sims_subs_signature_decile_index_accumulated_np_array[sim_num] += subs_signature_decile_index_accumulated_np_array
        all_sims_dinucs_signature_decile_index_accumulated_np_array[sim_num] += dinucs_signature_decile_index_accumulated_np_array
        all_sims_indels_signature_decile_index_accumulated_np_array[sim_num] += indels_signature_decile_index_accumulated_np_array
    #########################################################################################

    jobs = []

    #######################################################################################################################
    # August 1, 2020
    # Using numpy array
    if (computationType == USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM):
        sim_nums = range(0, numofSimulations + 1)
        sim_num_chr_tuples = ((sim_num, chrLong) for sim_num in sim_nums for chrLong in chromNamesList)

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
                                       chrBased_grouped_decile_df_list,
                                       ordered_all_sbs_signatures_array,
                                       ordered_all_dbs_signatures_array,
                                       ordered_all_id_signatures_array,
                                       ordered_sbs_signatures,
                                       ordered_dbs_signatures,
                                       ordered_id_signatures,
                                       ordered_sbs_signatures_cutoffs,
                                       ordered_dbs_signatures_cutoffs,
                                       ordered_id_signatures_cutoffs,
                                       is_discreet,
                                       verbose,),
                                 callback=accumulate_np_arrays))

        pool.close()
        pool.join()
    #######################################################################################################################

    ######################## starts July 22 2020 ###############################
    elif (computationType==USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT):

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
                                       sample_based,
                                       is_discreet,
                                       verbose,),
                                 callback=accumulate_np_arrays))

        pool.close()
        pool.join()
    ####################### ends July 22 2020   ###############################

    return all_sims_subs_signature_decile_index_accumulated_np_array, \
           all_sims_dinucs_signature_decile_index_accumulated_np_array, \
           all_sims_indels_signature_decile_index_accumulated_np_array
##################################################################



##################################################################
# Augment wavelet_processed_df with numberofAttributableBases
def augment(genome,wavelet_processed_df,matrix_generator_path,verbose):

    ################################
    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)
    ################################

    #Augment for each chromosome
    chrBased_wavelet_processed_df_groups = wavelet_processed_df.groupby(CHROM)

    frames=[]

    ####################################################################################
    # tuple contains (chrLong,chrBased_wavelet_processed_df_group)
    def accumulate_apply_async_result(tuple):
        chrBased_augmented_df = tuple[1]
        frames.append(chrBased_augmented_df)
    ####################################################################################

    for chrLong, chrBased_wavelet_processed_df_group in chrBased_wavelet_processed_df_groups:
        chrShort=getChrShort(chrLong)
        chrbased_filename = chrShort + ".txt"
        chrbased_file_path = os.path.join(matrix_generator_path, 'references', 'chromosomes', 'tsb', genome, chrbased_filename)
        if os.path.exists(chrbased_file_path):
            pool.apply_async(addNumofAttributableBasesColumnForApplyAsync, (chrLong,chrBased_wavelet_processed_df_group,chrbased_file_path,verbose), callback=accumulate_apply_async_result)
    #JAN 21, 2020 ends

    ################################
    pool.close()
    pool.join()
    ################################

    augment_df = pd.concat(frames, ignore_index=True)

    return augment_df
##################################################################


##################################################################
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

    ##############################################################
    my_list=[(SUBS,subs_signatures,all_sims_subs_signature_decile_index_accumulated_np_array),
             (DINUCS,dinucs_signatures, all_sims_dinucs_signature_decile_index_accumulated_np_array),
             (INDELS,indels_signatures, all_sims_indels_signature_decile_index_accumulated_np_array)]

    #Write Number of Attributable Bases List
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

                #Order is important
                # -1 contains AGGREGATEDSUBSTITUTIONS, AGGREGATEDDINUCS, AGGREGATEDINDELS for  my_type==SUBS, my_type==DINUCS, my_type==INDELS respectively
                # -2 contains REPEAT for  my_type==INDELS
                # -3 contains MICROHOMOLOGY for  my_type==INDELS
                if signature_index<signatures.size:
                    signature = signatures[signature_index]
                elif signature_index==(num_of_signatures_with_extras-1) and my_type==SUBS:
                    signature=AGGREGATEDSUBSTITUTIONS
                elif signature_index==(num_of_signatures_with_extras-1) and my_type==DINUCS:
                    signature=AGGREGATEDDINUCS
                elif signature_index==(num_of_signatures_with_extras-1) and my_type==INDELS:
                    signature=AGGREGATEDINDELS
                elif signature_index==(num_of_signatures_with_extras-2) and my_type==INDELS:
                    signature=REPEAT
                elif signature_index==(num_of_signatures_with_extras-3) and my_type==INDELS:
                    signature=MICROHOMOLOGY

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

                if (decile_df_list is not None):
                    numberofMutations, mutationDensityDict, numberofMutationsList = getMutationDensityDictUsingNumpyArray(decile_df_list, decile_counts_np_array)
                    normalizedMutationDensityList = getNormalizedMutationDensityList(mutationDensityDict)

                    #Normalized Mutation Density
                    with open(normalizedMutationDensityFilePath, 'w') as file:
                        for normalizedMutationDensity in normalizedMutationDensityList:
                            file.write(str(normalizedMutationDensity) + ' ')
                        file.write('\n')

                    #Number of Mutations
                    with open(numberofMutationsFilePath, 'w') as file:
                        for numberofMutations in numberofMutationsList:
                            file.write(str(numberofMutations) + ' ')
                        file.write('\n')
    ##############################################################

##################################################################

##################################################################
# old method, keep it for further guidance for sample_based
def writeReplicationTimeData(outputDir,jobname,sample_based,decile_df_list,decileIndex2NumberofAttributableBasesDict,simNum2Type2DecileBasedAllChrAccumulatedCountDict,simNum2Sample2Type2DecileBasedAllChrAccumulatedCountDict):
    os.makedirs(os.path.join(outputDir, jobname, DATA, REPLICATIONTIME), exist_ok=True)

    ##############################################################
    for simNum in simNum2Type2DecileBasedAllChrAccumulatedCountDict:
        type2DecileBasedAllChrAccumulatedCountDict = simNum2Type2DecileBasedAllChrAccumulatedCountDict[simNum]
        for type in type2DecileBasedAllChrAccumulatedCountDict:
            decileBasedAllChrAccumulatedCountDict = type2DecileBasedAllChrAccumulatedCountDict[type]

            if (simNum==0):
                normalizedMutationDensityFilename = '%s_NormalizedMutationDensity.txt' %(type)
                numberofMutationsFilename = '%s_NumberofMutations.txt' %(type)
                numberofAttributabelBasesFilename = '%s_NumberofAttributableBases.txt' %(type)
            else:
                normalizedMutationDensityFilename = '%s_sim%d_NormalizedMutationDensity.txt' %(type,simNum)
                numberofMutationsFilename = '%s_sim%d_NumberofMutations.txt' %(type,simNum)
                numberofAttributabelBasesFilename = '%s_sim%d_NumberofAttributableBases.txt' %(type,simNum)

            if (type==AGGREGATEDSUBSTITUTIONS) or (type==AGGREGATEDINDELS) or (type==AGGREGATEDDINUCS) or (type==MICROHOMOLOGY) or (type==REPEAT):
                os.makedirs(os.path.join(outputDir, jobname, DATA, REPLICATIONTIME,type), exist_ok=True)
                normalizedMutationDensityFilePath = os.path.join(outputDir,jobname, DATA, REPLICATIONTIME,type,normalizedMutationDensityFilename)
                numberofMutationsFilePath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME,type,numberofMutationsFilename)
                numberofAttributabelBasesFilePath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME,type,numberofAttributabelBasesFilename)
            else:
                os.makedirs(os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, SIGNATUREBASED), exist_ok=True)
                normalizedMutationDensityFilePath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, SIGNATUREBASED,normalizedMutationDensityFilename)
                numberofMutationsFilePath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, SIGNATUREBASED,numberofMutationsFilename)
                numberofAttributabelBasesFilePath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, SIGNATUREBASED,numberofAttributabelBasesFilename)

            # If decileBasedAllChrAccumulatedCountDict is not empty
            if (decileBasedAllChrAccumulatedCountDict):

                if (decile_df_list is not None):
                    numberofMutations, mutationDensityDict,numberofMutationsList, numberofAttributableBasesList = getMutationDensityDict(decile_df_list,decileBasedAllChrAccumulatedCountDict)
                elif (decileIndex2NumberofAttributableBasesDict is not None):
                    numberofMutations, mutationDensityDict,numberofMutationsList, numberofAttributableBasesList = getMutationDensityDictNewVersion(decileIndex2NumberofAttributableBasesDict, decileBasedAllChrAccumulatedCountDict)

                normalizedMutationDensityList = getNormalizedMutationDensityList(mutationDensityDict)

                #Normalized Mutation Density
                with open(normalizedMutationDensityFilePath, 'w') as file:
                    for normalizedMutationDensity in normalizedMutationDensityList:
                        file.write(str(normalizedMutationDensity) + ' ')
                    file.write('\n')

                #Number of Mutations
                with open(numberofMutationsFilePath, 'w') as file:
                    for numberofMutations in numberofMutationsList:
                        file.write(str(numberofMutations) + ' ')
                    file.write('\n')

                #Number of Attributable Bases
                with open(numberofAttributabelBasesFilePath, 'w') as file:
                    for numberofAttributabelBases in numberofAttributableBasesList:
                        file.write(str(numberofAttributabelBases) + ' ')
                    file.write('\n')
    ##############################################################


    ######################### Sample Based starts #####################
    if sample_based:
        for simNum in simNum2Sample2Type2DecileBasedAllChrAccumulatedCountDict:
            sample2Type2DecileBasedAllChrAccumulatedCountDict= simNum2Sample2Type2DecileBasedAllChrAccumulatedCountDict[simNum]
            for sample in sample2Type2DecileBasedAllChrAccumulatedCountDict:
                os.makedirs(os.path.join(outputDir, jobname, DATA, SAMPLES, sample, REPLICATIONTIME), exist_ok=True)

                for type in sample2Type2DecileBasedAllChrAccumulatedCountDict[sample]:
                    decileBasedAllChrAccumulatedCountDict = sample2Type2DecileBasedAllChrAccumulatedCountDict[sample][type]

                    if (simNum==0):
                        normalizedMutationDensityFilename = '%s_%s_NormalizedMutationDensity.txt' %(sample,type)
                    else:
                        normalizedMutationDensityFilename = '%s_%s_sim%d_NormalizedMutationDensity.txt' %(sample,type,simNum)

                    if (type == AGGREGATEDSUBSTITUTIONS) or (type == AGGREGATEDINDELS) or (type == AGGREGATEDDINUCS) or (type == MICROHOMOLOGY) or (type == REPEAT):
                        os.makedirs(os.path.join(outputDir, jobname, DATA,SAMPLES,sample, REPLICATIONTIME, type), exist_ok=True)
                        normalizedMutationDensityFilePath = os.path.join(outputDir, jobname, DATA, SAMPLES, sample,REPLICATIONTIME,type, normalizedMutationDensityFilename)
                    else:
                        os.makedirs(os.path.join(outputDir, jobname, DATA,SAMPLES,sample, REPLICATIONTIME, SIGNATUREBASED), exist_ok=True)
                        normalizedMutationDensityFilePath = os.path.join(outputDir, jobname, DATA, SAMPLES, sample,REPLICATIONTIME,SIGNATUREBASED, normalizedMutationDensityFilename)

                    #If decileBasedAllChrAccumulatedCountDict is not empty
                    if (decileBasedAllChrAccumulatedCountDict):
                        if (decile_df_list is not None):
                            numberofMutations, mutationDensityDict, numberofMutationsList, numberofAttributableBasesList = getMutationDensityDict(decile_df_list, decileBasedAllChrAccumulatedCountDict)
                        elif (decileIndex2NumberofAttributableBasesDict is not None):
                            numberofMutations, mutationDensityDict, numberofMutationsList, numberofAttributableBasesList = getMutationDensityDictNewVersion(decileIndex2NumberofAttributableBasesDict, decileBasedAllChrAccumulatedCountDict)
                        normalizedMutationDensityList = getNormalizedMutationDensityList(mutationDensityDict)

                        with open(normalizedMutationDensityFilePath, 'w') as file:
                            for normalizedMutationDensity in normalizedMutationDensityList:
                                file.write(str(normalizedMutationDensity) + ' ')
                            file.write('\n')
    ######################### Sample Based ends #######################

##################################################################


##################################################################
#main function
def replicationTimeAnalysis(computationType,
                            sample_based,
                            genome,
                            chromSizesDict,
                            chromNamesList,
                            outputDir,
                            jobname,
                            numofSimulations,
                            job_tuples,
                            repliseqDataFilename,
                            ordered_all_sbs_signatures_array,
                            ordered_all_dbs_signatures_array,
                            ordered_all_id_signatures_array,
                            ordered_sbs_signatures_with_cutoffs,
                            ordered_dbs_signatures_with_cutoffs,
                            ordered_id_signatures_with_cutoffs,
                            ordered_sbs_signatures_cutoffs,
                            ordered_dbs_signatures_cutoffs,
                            ordered_id_signatures_cutoffs,
                            is_discreet,
                            verbose,
                            matrix_generator_path):

    if is_discreet:
        sbs_signatures = ordered_sbs_signatures_with_cutoffs
        dbs_signatures = ordered_dbs_signatures_with_cutoffs
        id_signatures = ordered_id_signatures_with_cutoffs
    else:
        sbs_signatures = ordered_all_sbs_signatures_array
        dbs_signatures = ordered_all_dbs_signatures_array
        id_signatures = ordered_all_id_signatures_array

    print('\n#################################################################################')
    print('--- ReplicationTimeAnalysis starts')
    print('--- Replication Time Analyis Computation Type:%s' % (computationType))

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
    #What is the type of deciles? Deciles is a list of dataframes.
    # if verbose: print('\tVerbose Worker pid %s READ Repliseq DATA STARTS  %s MB' % (str(os.getpid()), memory_usage()))
    if verbose: print('\tVerbose READ Repliseq DATA STARTS')
    #Whole genome is needed here
    #Formerly I was reading 2bit files using twobitreader
    #Now, formerly downloaded matrix generator reference genome is being used.
    chrNamesInReplicationTimeDataArray, decile_df_list = readRepliSeqTimeData(genome,chromNamesList,repliseqDataFilename,matrix_generator_path,verbose)
    # if verbose: print('\tVerbose Worker pid %s READ Repliseq DATA ENDS  %s MB' % (str(os.getpid()), memory_usage()))
    if verbose: print('\tVerbose READ Repliseq DATA ENDS')

    #Get chrBased grouped deciles
    chrBased_grouped_decile_df_list = []
    #Get this chrLong of each decile
    #The first decile is the earliest one with index 1
    #The last decile is the latest one with index 10
    #We provide these indexes later in the fillChrBasedDecileIndexArray function by enumerate
    for decile_df in decile_df_list:
        chrBased_grouped_decile_df = decile_df.groupby(CHROM)
        chrBased_grouped_decile_df_list.append(chrBased_grouped_decile_df)
    ###################################################################
    ############### Read RepliSeq Time data ends ######################
    ###################################################################


    #######################################################################################################
    ################################### Replication Time Data Analysis starts #############################
    #######################################################################################################

    # #old method keep it for further guidance for sample based
    # writeReplicationTimeData(outputDir,jobname,sample_based,decile_df_list,None,simNum2Type2DecileBasedAllChrAccumulatedCountDict,simNum2Sample2Type2DecileBasedAllChrAccumulatedCountDict)

    # Ordered signatures will only have signatures since later on, they are used in filtering mutation row columns
    all_sims_subs_signature_decile_index_accumulated_np_array, \
    all_sims_dinucs_signature_decile_index_accumulated_np_array, \
    all_sims_indels_signature_decile_index_accumulated_np_array = calculateCountsForMutationsFillingReplicationTimeNPArrayRuntime_using_numpy_array(
        computationType,
        outputDir,
        jobname,
        numofSimulations,
        job_tuples,
        chromSizesDict,
        chromNamesList,
        chrBased_grouped_decile_df_list,
        ordered_all_sbs_signatures_array,
        ordered_all_dbs_signatures_array,
        ordered_all_id_signatures_array,
        ordered_sbs_signatures_with_cutoffs,
        ordered_dbs_signatures_with_cutoffs,
        ordered_id_signatures_with_cutoffs,
        ordered_sbs_signatures_cutoffs,
        ordered_dbs_signatures_cutoffs,
        ordered_id_signatures_cutoffs,
        sample_based,
        is_discreet,
        verbose)

    writeReplicationTimeDataUsingNumpyArray(outputDir,
                                            jobname,
                                            decile_df_list,
                                            sbs_signatures,
                                            dbs_signatures,
                                            id_signatures,
                                            all_sims_subs_signature_decile_index_accumulated_np_array,
                                            all_sims_dinucs_signature_decile_index_accumulated_np_array,
                                            all_sims_indels_signature_decile_index_accumulated_np_array)

    #######################################################################################################
    ################################### Replication Time Data Analysis ends ###############################
    #######################################################################################################

    print('--- ReplicationTimeAnalysis ends')
    print('#################################################################################\n')

##################################################################
