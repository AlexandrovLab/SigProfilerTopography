# !/usr/bin/env python3

# Author: burcakotlu

# Contact: burcakotlu@eng.ucsd.edu

from SigProfilerTopography.source.commons.TopographyCommons import AGGREGATEDSUBSTITUTIONS
from SigProfilerTopography.source.commons.TopographyCommons import AGGREGATEDINDELS
from SigProfilerTopography.source.commons.TopographyCommons import AGGREGATEDDINUCS
from SigProfilerTopography.source.commons.TopographyCommons import SIGNATUREBASED

from SigProfilerTopography.source.commons.TopographyCommons import START
from SigProfilerTopography.source.commons.TopographyCommons import TYPE
from SigProfilerTopography.source.commons.TopographyCommons import LENGTH

from SigProfilerTopography.source.commons.TopographyCommons import DATA

from SigProfilerTopography.source.commons.TopographyCommons import SUBS
from SigProfilerTopography.source.commons.TopographyCommons import INDELS
from SigProfilerTopography.source.commons.TopographyCommons import DINUCS

from SigProfilerTopography.source.commons.TopographyCommons import LNCRNA

from SigProfilerTopography.source.commons.TopographyCommons import get_chrBased_simBased_dfs

import pandas as pd
import numpy as np
import os
import multiprocessing

from intervaltree import Interval, IntervalTree

def read_lncRNA_GRCh37_annotation_file(file_path):
    column_names = ['chrom', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    # chr1	HAVANA	gene	29554	31109	.	+	.
    # ID=ENSG00000243485.5;
    # gene_id=ENSG00000243485.5_11;
    # gene_type=lncRNA;
    # gene_name=MIR1302-2HG;
    # level=2;
    # hgnc_id=HGNC:52482;
    # tag=ncRNA_host;
    # havana_gene=OTTHUMG00000000959.2_11;
    # remap_status=full_contig;
    # remap_num_mappings=1;
    # remap_target_status=overlap
    df = pd.read_csv(file_path, header=None, comment='#', sep="\t", names=column_names)
    # print('df.shape:', df.shape) # (294561, 9)

    # Add columns
    df['ID'] = df['attributes'].str.split(';', expand=True)[0].str.split('=', expand=True)[1]
    df['gene_id'] = df['attributes'].str.split(';', expand=True)[1].str.split('=', expand=True)[1]
    df['gene_type'] = df['attributes'].str.split(';', expand=True)[2].str.split('=', expand=True)[1]
    df['gene_name'] = df['attributes'].str.split(';', expand=True)[3].str.split('=', expand=True)[1]

    # print('df.shape:', df.shape) # (294561, 13)
    # print(df.head())
    # print('df.columns.values:', df.columns.values)
    #
    # print('############')
    # print('df[\'chrom\'].unique()', df['chrom'].unique(), df['chrom'].unique().shape)
    # #  ['chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chr10'
    # #  'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19'
    # #  'chr20' 'chr21' 'chr22' 'chrX' 'chrY'] (24,)
    # print('df[\'type\'].unique()', df['type'].unique(), df['type'].unique().shape) # ['gene' 'transcript' 'exon']
    # print('df[\'gene_type\'].unique()', df['gene_type'].unique(), df['gene_type'].unique().shape)
    # # ['lncRNA' 'ENSG00000243485.5_11' 'ENSG00000237613.2_6' ...
    # #  'ENSG00000228786.5_10' 'ENSG00000240450.1_8' 'ENSG00000231141.1_4'] (19966,)
    # print('############')
    #
    # print('############')
    # print(df['ID'])
    # print('############')
    # print(df['type'])
    # print('############')
    # print(df['gene_id'])
    # print('############')
    # print(df['gene_type'])
    # print('############')
    # print(df['gene_name'])
    # print('############')

    # filter only genes
    df = df[df['type']=='gene']
    # print('df.shape:', df.shape)  # (19974, 13)

    # filter some columns
    df = df[['chrom', 'type', 'start', 'end', 'strand', 'gene_type', 'gene_name']]
    # print('df.shape:', df.shape)  # (19974, 6)
    # print(df)
    # print('############')

    return df

def create_intervaltree(df):

    chrom_2_tree_dict = {}

    region_names_array = df['gene_name'].values
    sorted_region_names_array = np.sort(region_names_array)

    # print('region_names_array:', region_names_array)
    # print('region_names_array.shape:', region_names_array.shape, 'region_names_array.size:', region_names_array.size)

    # print('sorted_region_names_array:', sorted_region_names_array)
    # print('sorted_region_names_array.shape:', sorted_region_names_array.shape, 'sorted_region_names_array.size:', sorted_region_names_array.size)


    # test all chroms intervals in one tree
    # tuples = [Interval(start, end, (chrom, strand, gene_name)) for start, end, chrom, strand, gene_type, gene_name in zip(df['start'], df['end'], df['chrom'], df['strand'], df['gene_type'], df['gene_name'])]
    # lncRNA_tree = IntervalTree.from_tuples(tuples)

    # print('############')
    # print('lncRNA_tree:', lncRNA_tree)
    # print('############')
    # print(lncRNA_tree[150000:200000])
    # print('############')
    # print(len(lncRNA_tree[150000:200000]))
    # print('############')
    # print(lncRNA_tree.overlap(150000,200000))
    # print('############')
    # print(len(lncRNA_tree.overlap(150000, 200000)))
    # print('############')
    # print('len(lncRNA_tree):', len(lncRNA_tree))
    # print('############')

    # grouped_df = df.groupby('chrom')
    num_of_lncRNAs = 0

    for chrom, chrom_df  in df.groupby('chrom'):
        num_of_lncRNAs += chrom_df.shape[0]

        chrom_tuples = [Interval(start, end, (chrom, strand, end-start, gene_name))
                        for start, end, chrom, strand, gene_name
                        in zip(chrom_df['start'], chrom_df['end'], chrom_df['chrom'],
                               chrom_df['strand'], chrom_df['gene_name'])]

        chrom_tree = IntervalTree.from_tuples(chrom_tuples)
        chrom_2_tree_dict[chrom] = chrom_tree

    return num_of_lncRNAs, sorted_region_names_array, chrom_2_tree_dict



# Main engine function for for handling overlapping intervals
def accumulate_arrays(row,
                      overlapping_regions_set,
                      ordered_annotated_regions_array,
                      signatures_mask_array,
                      ordered_signatures_cutoffs,
                      accumulated_num_of_hits_np_array,
                      discreet_mode,
                      default_cutoff):

    # Fill numpy arrays using window_array
    if (overlapping_regions_set is not None) and (len(overlapping_regions_set) > 0):
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

        # Interval(start, end, (chrom, strand, gene_name))
        # for each overlap with the mutation, update that annotated region with 1 (from 0 to 1)
        regions_mask_array = np.isin(ordered_annotated_regions_array, np.array([overlap.data[2] for overlap in overlapping_regions_set]))
        annotated_regions_array = regions_mask_array.astype(int)

        # Add one more dimension to annotated_regions_array: (num_of_annotated_regions,) --> (1, num_of_annotated_regions)
        annotated_regions_array_1xnum_of_annotated_regions = np.expand_dims(annotated_regions_array, axis=0)

        to_be_accumulated_num_of_hits_array = mask_array_1xnumofsignatures.T * annotated_regions_array_1xnum_of_annotated_regions
        accumulated_num_of_hits_np_array += to_be_accumulated_num_of_hits_array




# Vectorization
def fillNumofHitsArray_using_list_comp(
        row,
        mutation_type,
        chrBased_tree,
        ordered_annotated_regions_array,
        ordered_signatures_cutoffs,
        signatures_mask_array,
        accumulated_num_of_hits_np_array,
        discreet_mode,
        default_cutoff,
        df_columns):

    # df_columns 'numpy.ndarray'
    # df_columns: ['Sample' 'Chrom' 'Start' 'MutationLong' 'PyramidineStrand'
    #              'TranscriptionStrand' 'Mutation' 'SBS1' 'SBS2' 'SBS3' 'SBS4' 'SBS5'
    #              'SBS8' 'SBS10a' 'SBS10b' 'SBS12' 'SBS13' 'SBS14' 'SBS15' 'SBS16' 'SBS17a'
    #              'SBS17b' 'SBS18' 'SBS20' 'SBS22' 'SBS23' 'SBS25' 'SBS28' 'SBS30' 'SBS33'
    #              'SBS34' 'SBS39' 'SBS40' 'SBS44' 'SBS288P' 'Simulation_Number']


    indexofStart = np.where(df_columns == START)[0][0]
    mutation_row_start = row[indexofStart]

    if mutation_type == SUBS:
        length = 1
    elif mutation_type == DINUCS:
        length = 2
    elif mutation_type == INDELS:
        indexofLength = np.where(df_columns == LENGTH)[0][0]
        length = row[indexofLength]

    overlapping_regions_set = chrBased_tree[mutation_row_start:mutation_row_start+length]

    accumulate_arrays(row,
                      overlapping_regions_set,
                      ordered_annotated_regions_array,
                      signatures_mask_array,
                      ordered_signatures_cutoffs,
                      accumulated_num_of_hits_np_array,
                      discreet_mode,
                      default_cutoff)



# requires chrBased_simBased_combined_df_split which can be real split or whole in fact
# This is common for pool.imap_unordered and pool.apply_async variations
def chrbased_data_fill_num_of_hits_arrays_for_all_mutations(chrBased_tree,
                                                            num_of_annotated_regions,
                                                            ordered_annotated_regions_array,
                                                             chrLong,
                                                             simNum,
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

    ###############################################################################
    ################################ Initialization ###############################
    ###############################################################################
    # Add one more row for the aggregated analysis
    subsSignature_accumulated_num_of_hits_np_array = np.zeros((number_of_sbs_signatures + 1, num_of_annotated_regions)) # dtype=float
    dinucsSignature_accumulated_num_of_hits_np_array = np.zeros((number_of_dbs_signatures + 1, num_of_annotated_regions)) # dtype=float
    indelsSignature_accumulated_num_of_hits_np_array = np.zeros((number_of_id_signatures + 1, num_of_annotated_regions)) # dtype=float
    ###############################################################################
    ################################ Initialization ###############################
    ###############################################################################

    if ((chrBased_tree is not None) and (len(chrBased_tree) > 0)):
        ######################################################## #######################
        ################### Fill signal and count array starts ########################
        ###############################################################################


        # For subs
        if ((chrBased_simBased_subs_df is not None) and (not chrBased_simBased_subs_df.empty)):


            # df_columns is a numpy array
            df_columns = chrBased_simBased_subs_df.columns.values
            subsSignatures_mask_array = np.isin(df_columns, ordered_sbs_signatures)


            [fillNumofHitsArray_using_list_comp(
                row,
                SUBS,
                chrBased_tree,
                ordered_annotated_regions_array,
                ordered_sbs_signatures_cutoffs,
                subsSignatures_mask_array,
                subsSignature_accumulated_num_of_hits_np_array,
                discreet_mode,
                default_cutoff,
                df_columns) for row in chrBased_simBased_subs_df[df_columns].values]


        # For Dinusc
        if ((chrBased_simBased_dinucs_df is not None) and (not chrBased_simBased_dinucs_df.empty)):

            # df_columns is a numpy array
            df_columns = chrBased_simBased_dinucs_df.columns.values
            dinucsSignatures_mask_array = np.isin(df_columns, ordered_dbs_signatures)

            [fillNumofHitsArray_using_list_comp(
                row,
                DINUCS,
                chrBased_tree,
                ordered_annotated_regions_array,
                ordered_dbs_signatures_cutoffs,
                dinucsSignatures_mask_array,
                dinucsSignature_accumulated_num_of_hits_np_array,
                discreet_mode,
                default_cutoff,
                df_columns) for row in chrBased_simBased_dinucs_df[df_columns].values]

        # For Indels
        if ((chrBased_simBased_indels_df is not None) and (not chrBased_simBased_indels_df.empty)):

            # df_columns is a numpy array
            df_columns = chrBased_simBased_indels_df.columns.values
            indelsSignatures_mask_array = np.isin(df_columns, ordered_id_signatures)

            [fillNumofHitsArray_using_list_comp(
                row,
                INDELS,
                chrBased_tree,
                ordered_annotated_regions_array,
                ordered_id_signatures_cutoffs,
                indelsSignatures_mask_array,
                indelsSignature_accumulated_num_of_hits_np_array,
                discreet_mode,
                default_cutoff,
                df_columns) for row in chrBased_simBased_indels_df[df_columns].values]

        ###############################################################################
        ################### Fill signal and count array ends ##########################
        ###############################################################################

    # Initialzie the list, you will return this list
    NumofHitsArrayList = []

    # print('DEBUG8',
    #       chrLong,
    #       'simNum:', simNum,
    #       'subsSignature_accumulated_num_of_hits_np_array.shape', subsSignature_accumulated_num_of_hits_np_array.shape,
    #       'dinucsSignature_accumulated_num_of_hits_np_array.shape', dinucsSignature_accumulated_num_of_hits_np_array.shape,
    #       'indelsSignature_accumulated_num_of_hits_np_array.shape', indelsSignature_accumulated_num_of_hits_np_array.shape)


    NumofHitsArrayList.append(chrLong) # 0
    NumofHitsArrayList.append(simNum) # 1
    NumofHitsArrayList.append(subsSignature_accumulated_num_of_hits_np_array) # 2
    NumofHitsArrayList.append(dinucsSignature_accumulated_num_of_hits_np_array) # 3
    NumofHitsArrayList.append(indelsSignature_accumulated_num_of_hits_np_array) # 4

    return NumofHitsArrayList


# For apply_async
# Read chromBased simBased combined mutations df in the process
def chrbased_fill_num_of_hits_arrays_for_all_mutations_read_mutations(
        chrBased_tree,
        num_of_annotated_regions,
        ordered_annotated_regions_array,
        outputDir,
        jobname,
        chrLong,
        simNum,
        samples_of_interest,
        ordered_sbs_signatures,
        ordered_dbs_signatures,
        ordered_id_signatures,
        ordered_sbs_signatures_cutoffs,
        ordered_dbs_signatures_cutoffs,
        ordered_id_signatures_cutoffs,
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

        # print('DEBUG3',
        #       chrLong,
        #       'simNum:', simNum,
        #       'len(chrBased_tree):', len(chrBased_tree),
        #       'num_of_annotated_regions:', num_of_annotated_regions,
        #       'ordered_annotated_regions_array.shape:', ordered_annotated_regions_array.shape,
        #       'chrBased_simBased_subs_df.shape:', chrBased_simBased_subs_df.shape,
        #       'chrBased_simBased_dinucs_df.shape:', chrBased_simBased_dinucs_df.shape,
        #       'chrBased_simBased_indels_df.shape:', chrBased_simBased_indels_df.shape)

        return chrbased_data_fill_num_of_hits_arrays_for_all_mutations(chrBased_tree,
                                                                       num_of_annotated_regions,
                                                                        ordered_annotated_regions_array,
                                                                        chrLong,
                                                                        simNum,
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

    except Exception as e:
        log_out = open(log_file, 'a')
        print("Exception in chrbased_fill_num_of_hits_arrays_for_all_mutations_read_mutations: %s" % (e), file=log_out)
        log_out.close()


# simulationNumber == 0 means original data
# simulationNumber > 0 means simulation data
def writeNumofHitsFiles(region_type,
                        allMutationsAccumulatedNumofHitsArray,
                        outputDir,
                        jobname,
                        aggregated_mutations_type,
                        simulationNumber):

    os.makedirs(os.path.join(outputDir, jobname, DATA, region_type, aggregated_mutations_type),exist_ok=True)

    if (simulationNumber == 0):
        accumulated_num_of_hits_filename = '%s_AccumulatedNumofHits.txt' %(jobname)
    else:
        accumulated_num_of_hits_filename = '%s_sim%d_AccumulatedNumofHits.txt' %(jobname, simulationNumber)

    accumulatedSignalFilePath = os.path.join(outputDir, jobname,DATA, region_type, aggregated_mutations_type, accumulated_num_of_hits_filename)
    allMutationsAccumulatedNumofHitsArray.tofile(file=accumulatedSignalFilePath, sep="\t", format="%s")

# TODO
def writeSignatureBasedAccumulatedNumofHitsFilesUsingNumpyArray(region_type,
                                                           subsSignatures,
                                                           dinucsSignatures,
                                                           indelsSignatures,
                                                           subsSignature_accumulated_num_of_hits_np_array,
                                                           dinucsSignature_accumulated_num_of_hits_np_array,
                                                           indelsSignature_accumulated_num_of_hits_np_array,
                                                           outputDir,
                                                           jobname,
                                                           simNum):

    os.makedirs(os.path.join(outputDir, jobname, DATA, region_type, SIGNATUREBASED), exist_ok=True)

    all_signatures = [subsSignatures,
                      dinucsSignatures,
                      indelsSignatures]

    all_num_of_hits_arrays = [subsSignature_accumulated_num_of_hits_np_array,
                              dinucsSignature_accumulated_num_of_hits_np_array,
                              indelsSignature_accumulated_num_of_hits_np_array]

    for signatures_index, signatures in enumerate(all_signatures,0):
        number_of_signatures = signatures.size
        num_of_hits_arrays = all_num_of_hits_arrays[signatures_index]

        # Why -1, because last one is aggregated mutations and it is not written here.
        for signature_index in range(0, number_of_signatures-1):
            signature = signatures[signature_index]
            num_of_hits_array = num_of_hits_arrays[signature_index]

            # To provide filename with no space in signature name
            # signatureWithNoSpace = signature.replace(' ','')

            if (simNum == 0):
                accumulated_num_of_hits_filename = '%s_AccumulatedSignalArray.txt' %(signature)
            else:
                accumulated_num_of_hits_filename = '%s_sim%d_AccumulatedSignalArray.txt' %(signature, simNum)

            accumulated_num_of_hits_file_path = os.path.join(outputDir, jobname, DATA, region_type, SIGNATUREBASED, accumulated_num_of_hits_filename)
            num_of_hits_array.tofile(file=accumulated_num_of_hits_file_path, sep="\t",format="%s")



def writeSimulationBasedOverlappingLNCRNAUsingNumpyArray(region_type,
                                                   subsSignatures,
                                                   dinucsSignatures,
                                                   indelsSignatures,
                                                   allSims_subsSignature_accumulated_num_of_hits_np_array,
                                                   allSims_dinucsSignature_accumulated_num_of_hits_np_array,
                                                   allSims_indelsSignature_accumulated_num_of_hits_np_array,
                                                   outputDir,
                                                   jobname,
                                                   numofSimulations):

    os.makedirs(os.path.join(outputDir, jobname, DATA, region_type),exist_ok=True)

    for simNum in range(0, numofSimulations + 1):

        subsSignature_accumulated_num_of_hits_np_array = allSims_subsSignature_accumulated_num_of_hits_np_array[simNum]
        dinucsSignature_accumulated_num_of_hits_np_array =  allSims_dinucsSignature_accumulated_num_of_hits_np_array[simNum]
        indelsSignature_accumulated_num_of_hits_np_array = allSims_indelsSignature_accumulated_num_of_hits_np_array[simNum]

        # Last row contains AGGREGATEDSUBSTITUTIONS
        writeNumofHitsFiles(region_type,
                            subsSignature_accumulated_num_of_hits_np_array[-1],
                            outputDir,
                            jobname,
                            AGGREGATEDSUBSTITUTIONS,
                            simNum)

        # Last row contains AGGREGATEDDINUCS
        writeNumofHitsFiles(region_type,
                            dinucsSignature_accumulated_num_of_hits_np_array[-1],
                            outputDir,
                            jobname,
                            AGGREGATEDDINUCS,
                            simNum)

        # Last row contains AGGREGATEDINDELS
        writeNumofHitsFiles(region_type,
                            indelsSignature_accumulated_num_of_hits_np_array[-1],
                            outputDir,
                            jobname,
                            AGGREGATEDINDELS,
                            simNum)

        # Signatures
        writeSignatureBasedAccumulatedNumofHitsFilesUsingNumpyArray(region_type,
                                                           subsSignatures,
                                                           dinucsSignatures,
                                                           indelsSignatures,
                                                           subsSignature_accumulated_num_of_hits_np_array,
                                                           dinucsSignature_accumulated_num_of_hits_np_array,
                                                           indelsSignature_accumulated_num_of_hits_np_array,
                                                           outputDir,
                                                           jobname,
                                                           simNum)



# main function
def annotated_region_analysis(genome,
                            outputDir,
                            jobname,
                            numofSimulations,
                            region_type,
                            chromNamesList,
                            samples_of_interest,
                            ordered_sbs_signatures_array,
                            ordered_dbs_signatures_array,
                            ordered_id_signatures_array,
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
    print('--- Annotated Region Analysis starts', file=log_out)
    log_out.close()

    number_of_sbs_signatures = ordered_sbs_signatures_array.size
    number_of_dbs_signatures = ordered_dbs_signatures_array.size
    number_of_id_signatures = ordered_id_signatures_array.size

    subsSignatures = np.append(ordered_sbs_signatures_array, AGGREGATEDSUBSTITUTIONS)
    dinucsSignatures = np.append(ordered_dbs_signatures_array, AGGREGATEDDINUCS)
    indelsSignatures = np.append(ordered_id_signatures_array, AGGREGATEDINDELS)

    num_of_annotated_regions = 0
    ordered_annotated_regions_array = None
    chrom_2_tree_dict = None

    if region_type == LNCRNA:
        # TODO genome based
        lncRNA_file_path = os.path.join('/restricted/alexandrov-ddn/users/burcak/data/GENCODE/GRCh37/lncRNA',
                                        'gencode.v43lift37.long_noncoding_RNAs.gff3')
        lncRNA_df = read_lncRNA_GRCh37_annotation_file(lncRNA_file_path)
        num_of_lncRNAs, ordered_lncRNA_regions, chrom_2_tree_dict = create_intervaltree(lncRNA_df)
        num_of_annotated_regions = num_of_lncRNAs
        ordered_annotated_regions_array = ordered_lncRNA_regions

    # For Vectorization
    # These are used in writing tables
    allSims_subsSignature_accumulated_num_of_hits_np_array = np.zeros((numofSimulations+1, number_of_sbs_signatures + 1, num_of_annotated_regions))
    allSims_dinucsSignature_accumulated_num_of_hits_np_array = np.zeros((numofSimulations+1, number_of_dbs_signatures + 1, num_of_annotated_regions))
    allSims_indelsSignature_accumulated_num_of_hits_np_array = np.zeros((numofSimulations+1, number_of_id_signatures + 1, num_of_annotated_regions))

    # accumulate number of hits in the annotated regions comuing from all mutations in all chromosomes
    def accumulate_apply_async_result_vectorization(simulatonBased_num_of_hits_arrays_list):
        try:
            chrLong = simulatonBased_num_of_hits_arrays_list[0]
            simNum = simulatonBased_num_of_hits_arrays_list[1]
            subsSignature_accumulated_num_of_hits_np_array = simulatonBased_num_of_hits_arrays_list[2]
            dinucsSignature_accumulated_num_of_hits_np_array = simulatonBased_num_of_hits_arrays_list[3]
            indelsSignature_accumulated_num_of_hits_np_array = simulatonBased_num_of_hits_arrays_list[4]

            # Accumulation
            allSims_subsSignature_accumulated_num_of_hits_np_array[simNum] += subsSignature_accumulated_num_of_hits_np_array
            allSims_dinucsSignature_accumulated_num_of_hits_np_array[simNum] += dinucsSignature_accumulated_num_of_hits_np_array
            allSims_indelsSignature_accumulated_num_of_hits_np_array[simNum] += indelsSignature_accumulated_num_of_hits_np_array
            # print('ACCUMULATION chrLong:%s simNum:%d ENDS' %(chrLong,simNum))

        except Exception as e:
            print("Exception in accumulate_apply_async_result_vectorization function: %s" % (e))

    sim_nums = range(0, numofSimulations + 1)
    sim_num_chr_tuples = ((sim_num, chrLong) for sim_num in sim_nums for chrLong in chromNamesList)

    if parallel_mode:
        jobs = []
        numofProcesses = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=numofProcesses)

        for simNum, chrLong in sim_num_chr_tuples:
            if chrLong in chrom_2_tree_dict:
                chrBased_tree = chrom_2_tree_dict[chrLong]
                jobs.append(pool.apply_async(chrbased_fill_num_of_hits_arrays_for_all_mutations_read_mutations,
                                             args=(chrBased_tree,
                                                   num_of_annotated_regions,
                                                   ordered_annotated_regions_array,
                                                   outputDir,
                                                   jobname,
                                                   chrLong,
                                                   simNum,
                                                   samples_of_interest,
                                                   ordered_sbs_signatures_array,
                                                   ordered_dbs_signatures_array,
                                                   ordered_id_signatures_array,
                                                   ordered_sbs_signatures_cutoffs,
                                                   ordered_dbs_signatures_cutoffs,
                                                   ordered_id_signatures_cutoffs,
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
            if chrLong in chrom_2_tree_dict:
                chrBased_tree = chrom_2_tree_dict[chrLong]
                simulatonBased_SignalArrayAndCountArrayList = chrbased_fill_num_of_hits_arrays_for_all_mutations_read_mutations(
                                chrBased_tree,
                                num_of_annotated_regions,
                                ordered_annotated_regions_array,
                                outputDir,
                                jobname,
                                chrLong,
                                simNum,
                                samples_of_interest,
                                ordered_sbs_signatures_array,
                                ordered_dbs_signatures_array,
                                ordered_id_signatures_array,
                                ordered_sbs_signatures_cutoffs,
                                ordered_dbs_signatures_cutoffs,
                                ordered_id_signatures_cutoffs,
                                discreet_mode,
                                default_cutoff,
                                log_file,
                                verbose)

                # if simulatonBased_SignalArrayAndCountArrayList is None:
                #     print('DEBUG2', chrLong, 'simulatonBased_SignalArrayAndCountArrayList is None')

                accumulate_apply_async_result_vectorization(simulatonBased_SignalArrayAndCountArrayList)


    # Same for parallel_mode or sequential run
    writeSimulationBasedOverlappingLNCRNAUsingNumpyArray(region_type,
                                                   subsSignatures,
                                                   dinucsSignatures,
                                                   indelsSignatures,
                                                   allSims_subsSignature_accumulated_num_of_hits_np_array,
                                                   allSims_dinucsSignature_accumulated_num_of_hits_np_array,
                                                   allSims_indelsSignature_accumulated_num_of_hits_np_array,
                                                   outputDir,
                                                   jobname,
                                                   numofSimulations)




    log_out = open(log_file, 'a')
    print('--- Annotated Region Analysis ends', file=log_out)
    print('#################################################################################\n', file=log_out)
    log_out.close()
