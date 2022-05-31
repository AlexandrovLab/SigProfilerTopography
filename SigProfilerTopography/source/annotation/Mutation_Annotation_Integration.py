# !/usr/bin/env python3

# Author: burcakotlu

# Contact: burcakotlu@eng.ucsd.edu


"""

This python file integrates somatic mutations annotations with SigProfilerTopography analyses.
This python file aims to answer the following questions:
1) What types of mutations are found from early to late replicating regions?
2) Is there a consistent behaviour among the samples?
3) Is there a pattern for mutations found in early replicating regions? (across all mutations and signature specific)
4) Is there a pattern for mutations found in late replicating regions? (across all mutations and signature specific)
5) Is there a mutational signature specific pattern?
6) Is there a cancer type specific pattern? e.g.: Liver cancer, Lung cancer

"""


from pathlib import Path

import pandas as pd
import os
import numpy as np
from matplotlib import pyplot as plt

CHR_START = 'Chr_Start'

topography_vep_figures = 'topography_vep_figures'
topography_vep_text_files = 'topography_vep_text_files'

def get_samples(inputDir):
    samples = []

    if os.path.exists(inputDir):
        files_list = os.listdir(inputDir)
        samples = [file[:-4] for file in files_list if (file.endswith('.vcf'))]

    return samples


# Somatic mutation annotations integration with SigProfilerTopography replication timing analysis
def mutation_annotation_replication_timing_integration(inputDir,
                                                       outputDir,
                                                       jobname,
                                                       cancer_type = None):

    os.makedirs(os.path.join(outputDir, jobname, topography_vep_text_files), exist_ok=True)
    os.makedirs(os.path.join(outputDir, jobname, topography_vep_figures), exist_ok=True)

    # my colors for 21 colors
    my_colors = ['red', 'royalblue', 'pink', 'black', 'limegreen', 'darkblue', 'gold', 'coral', 'gray', 'springgreen',
                 'mediumvioletred', 'deepskyblue', 'violet', 'black', 'yellow', 'lightgreen', 'orange', 'indianred',
                 'gray', 'green', 'cyan']

    if cancer_type is not None and ((cancer_type == 'Liver-HCC') or (cancer_type == 'Lung-AdenoCA')) :
        vep_files_path = os.path.join('/restricted/alexandrov-group/burcak/data/PCAWG/' + cancer_type + '/vep_files')
    else:
        vep_files_path = os.path.join(Path(inputDir).parents[0], 'vep_files')

    vep_columns = ['Uploaded_variation', 'Location', 'Allele', 'Gene', 'Feature', 'Feature_type', 'Consequence',
                   'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation',
                   'Extra']

    # Collect all possible consequences across all samples and all deciles
    all_consequences = set()

    # Fill for 10 replication timing deciles
    decile_dict = {}

    range_max = 11 # from decile index: 1 to 10

    for decile_index in range(1,range_max):

        # initialize for each decile
        all_samples_dict = {}

        # main topography replication timing data files showing mutations in each replication timing decile
        filename = 'Mutations_decile%s_replicating_timing.txt' %(decile_index)
        filepath = os.path.join(outputDir, jobname, filename)
        mutations_decile_replicating_timing_df = pd.read_csv(filepath, sep='\t', header=0)

        grouped = mutations_decile_replicating_timing_df.groupby('Sample')

        for sample, sample_decile_df in grouped:
            # add a new column to merge with sample vep file
            sample_decile_df[CHR_START] = sample_decile_df['Chrom'].astype(str) + ':' + sample_decile_df['Start'].astype(str)

            # read sample vep file that comes from  vep call
            # sample <- SP97681 in mutations_decile_replicating_timing_df
            # under vep_files_path
            # LICA-FR_SP97681.snv_mnv_filtered.txt
            # LICA-FR_SP97681.snv_mnv_filtered.txt_summary.html
            # LICA-FR_SP97681.indel_filtered.txt
            # LICA-FR_SP97681.indel_filtered.txt_summary.html
            if (cancer_type is not None) and ((cancer_type == 'Liver-HCC') or (cancer_type == 'Lung-AdenoCA')) and os.path.exists(vep_files_path):
                files_list = os.listdir(vep_files_path)
                filenames = [file[:-4] for file in files_list if (sample in file) and (file.endswith('.snv_mnv_filtered.txt'))]

                if len(filenames)>0:
                    filename = filenames[0]
                    sample_vep_filepath = os.path.join(vep_files_path, filename + '.txt')  # TODO change here for Liver and Lung cancer samples
                else:
                    continue

            else:
                # read sample vep file that comes from  vep call
                sample_vep_filepath = os.path.join(vep_files_path,
                                                   sample + '.txt')  # TODO change here for Liver and Lung cancer samples

            # read sample vep file that comes from  vep call
            sample_vep_df = pd.read_csv(sample_vep_filepath, sep='\t', comment='#')
            sample_vep_df.columns = vep_columns

            # merge sample_decile_df and sample_vep_df
            merged_df = pd.merge(sample_decile_df, sample_vep_df, how='left', left_on=[CHR_START], right_on=['Location'])
            merged_grouped = merged_df.groupby('Consequence')

            # decile based and sample based: consequence -> number of consequences
            sample_consequences_dict = {consequence: len(consequence_df) for consequence, consequence_df in merged_grouped}

            # Breakdown consequences
            sample_new_consequences_dict = {}

            for k, v in sample_consequences_dict.items():
                k_list = k.strip().split(',')
                for knew in k_list:
                    if knew not in sample_new_consequences_dict:
                        sample_new_consequences_dict[knew] = v
                    else:
                        sample_new_consequences_dict[knew] += v

            assert sample not in all_samples_dict

            if sample not in all_samples_dict:
                all_samples_dict[sample] = sample_new_consequences_dict

        # after all samples are visited for this decile index
        decile_dict[decile_index] = all_samples_dict

        # after all samples are visited for this decile index
        for sample, sample_dict in all_samples_dict.items():
            sample_consequences = sample_dict.keys()
            all_consequences = all_consequences.union(set(sample_consequences))

    # after all deciles and samples are visited
    # order consequence and create list for each sample
    ordered_all_consequences = sorted(all_consequences)
    print('type(ordered_all_consequences):',  type(ordered_all_consequences),
          'len(ordered_all_consequences):', len(ordered_all_consequences),
          'ordered_all_consequences:', ordered_all_consequences)

    # For all deciles df
    data_for_df_deciles = []
    df_decile_columns = ['Decile'] + ordered_all_consequences

    # Plot png figure for each decile
    for decile_index in range(1,range_max):
        decile_list = [decile_index]
        all_samples_dict = decile_dict[decile_index]

        # For each decile df
        data_for_df = []
        df_columns = ['Sample'] + ordered_all_consequences
        samples = all_samples_dict.keys()

        for sample in samples:
            sample_list= [sample]
            for consequence in ordered_all_consequences:
                if consequence in all_samples_dict[sample]:
                    count = all_samples_dict[sample][consequence]
                else:
                    count = 0
                sample_list.append(count)
            data_for_df.append(sample_list)

        # Number of Mutations for each VEP consequence
        df = pd.DataFrame(data_for_df, columns=df_columns)

        filename = 'Samples_Decile_%s_VEP_Consequences.txt' % (decile_index)
        filepath = os.path.join(outputDir, jobname, topography_vep_text_files, filename)
        df.to_csv(filepath, sep='\t', index=False)

        # Sum number of mutations for each consequence across all samples
        decile_list.extend(df.iloc[:,1:].sum(axis=0).values.tolist())
        data_for_df_deciles.append(decile_list)

        # Percentages for each VEP consequence
        df_perc = df
        df_perc.iloc[:,1:] = df.iloc[:, 1:].div(df.sum(axis=1), axis=0)
        df_perc.iloc[:,1:] = df_perc.iloc[:,1:] * 100
        filename = 'Samples_Decile_%s_VEP_Consequences_Percentages.txt' % (decile_index)
        filepath = os.path.join(outputDir, jobname, topography_vep_text_files, filename)
        df_perc.to_csv(filepath, sep='\t', index=False)

        # Plot data in stack manner of bar type
        fig = df_perc.plot(x='Sample', kind='bar', stacked=True, figsize=(40, 30), color=my_colors, fontsize=20, title='Samples Mutations Consequences Percentages in Decile%s' %(decile_index)).get_figure()

        plt.legend(bbox_to_anchor=(1.01, 1.0), fontsize=20)

        filename = 'Samples_Consequence_Percentages_Decile%s.png' %(decile_index)
        figureFile = os.path.join(outputDir, jobname, topography_vep_figures, filename)

        fig.savefig(figureFile)
        plt.close()

    # Aggregated Mutations
    # Dataframe Deciles Sum of Consequences Across All Samples Text File
    df_deciles = pd.DataFrame(data_for_df_deciles, columns=df_decile_columns)

    filename = 'Deciles_Across_All_Samples_VEP_Consequences.txt'
    filepath = os.path.join(outputDir, jobname, topography_vep_text_files, filename)
    df_deciles.to_csv(filepath, sep='\t', index=False)

    # Plot Deciles on x axis, y axis Sum of Consequences Across All Samples Figure
    # x axis deciles, y vep consequences as stacked bar plots
    fig = df_deciles.plot(x='Decile', kind='bar', stacked=True, figsize=(40, 30), color=my_colors, fontsize=20, title='Decile Mutations VEP Consequences').get_figure()
    plt.legend(bbox_to_anchor=(1.01, 1.0), fontsize=20)

    filename = 'Deciles_Across_All_Samples_VEP_Consequences_Numbers.png'
    figureFile = os.path.join(outputDir, jobname, topography_vep_figures, filename)

    fig.savefig(figureFile)
    plt.close()

    # Aggregated Mutations
    # Dataframe Deciles Percentage of Consequences Across All Samples Text file
    df_deciles_percentages = df_deciles
    df_deciles_percentages.iloc[:, 1:] = df_deciles.iloc[:, 1:].div(df_deciles.sum(axis=1), axis=0)
    df_deciles_percentages.iloc[:, 1:] = df_deciles_percentages.iloc[:, 1:] * 100

    filename = 'Deciles_Across_All_Samples_VEP_Consequences_Percentages.txt'
    filepath = os.path.join(outputDir, jobname, topography_vep_text_files, filename)
    df_deciles_percentages.to_csv(filepath, sep='\t', index=False)

    # Dataframe Deciles Percentage of Consequences Across All Samples Figure
    # # x axis deciles, y percentages of vep consequences as stacked bar plots

    # from test
    # colormap = "tab20"
    # color = my_colors,
    fig = df_deciles_percentages.plot(x='Decile', kind='bar', stacked=True, figsize=(40, 35), color = my_colors,
                                      fontsize=20,
                                      title='Decile Mutations VEP Consequences Percentages').get_figure()

    plt.xlabel('Early to Late Replication Timing', fontsize=30)
    plt.ylabel('Percentage of Mutations Annotations', fontsize=30)

    # put lower left of box to (x=0, y=1)
    # plt.gca().legend(labels, bbox_to_anchor=([0, 1]), loc='lower left', ncol=4, fontsize=30, frameon=False)
    plt.gca().legend(bbox_to_anchor=([0, 1]), loc='lower left', ncol=4, fontsize=30, frameon=False)

    filename = 'Deciles_Across_All_Samples_VEP_Consequences_Percentages.png'
    figureFile = os.path.join(outputDir, jobname, topography_vep_figures, filename)

    fig.savefig(figureFile)
    plt.tight_layout()
    plt.close()

    return ordered_all_consequences


def mutation_annotation_replication_timing_integration_signature_specific(inputDir,
                                                                          outputDir,
                                                                          jobname,
                                                                          ordered_all_consequences,
                                                                          signature_cutoff_df,
                                                                          cancer_type = None):

    os.makedirs(os.path.join(outputDir, jobname, topography_vep_text_files), exist_ok=True)
    os.makedirs(os.path.join(outputDir, jobname, topography_vep_figures), exist_ok=True)

    # my colors for 21 colors
    my_colors = ['red', 'royalblue', 'pink', 'black', 'limegreen', 'darkblue', 'gold', 'coral', 'gray', 'springgreen',
                 'mediumvioletred', 'deepskyblue', 'violet', 'black', 'yellow', 'lightgreen', 'orange', 'indianred',
                 'gray', 'green', 'cyan']


    if cancer_type is not None and ((cancer_type == 'Liver-HCC') or (cancer_type == 'Lung-AdenoCA')):
        vep_files_path = os.path.join('/restricted/alexandrov-group/burcak/data/PCAWG/' + cancer_type + '/vep_files')
    else:
        vep_files_path = os.path.join(Path(inputDir).parents[0], 'vep_files')

    vep_columns = ['Uploaded_variation', 'Location', 'Allele', 'Gene', 'Feature', 'Feature_type', 'Consequence',
                   'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation',
                   'Extra']

    print('subsSignature_cutoff_numberofmutations_averageprobability_df:', signature_cutoff_df)

    signatures = signature_cutoff_df['signature'].unique()

    for signature in signatures:
        cutoff = signature_cutoff_df[signature_cutoff_df['signature']==signature]['cutoff'].values[0]

        # Fill for 10 replication timing deciles
        decile_dict = {}

        range_max = 11 # from decile index: 1 to 10

        for decile_index in range(1,range_max):

            # initialize for each decile
            all_samples_dict = {}

            # main topography replication timing data files showing mutations in each replication timing decile
            filename = 'Mutations_decile%s_replicating_timing.txt' %(decile_index)
            filepath = os.path.join(outputDir, jobname, filename)
            mutations_decile_replicating_timing_df = pd.read_csv(filepath, sep='\t', header=0)

            # filter mutations_decile_replicating_timing_df for cutoff
            mutations_decile_replicating_timing_df= mutations_decile_replicating_timing_df[mutations_decile_replicating_timing_df[signature] >= cutoff]

            grouped = mutations_decile_replicating_timing_df.groupby('Sample')

            for sample, sample_decile_df in grouped:
                # add a new column to merge with sample vep file
                sample_decile_df[CHR_START] = sample_decile_df['Chrom'].astype(str) + ':' + sample_decile_df['Start'].astype(str)

                # read sample vep file that comes from  vep call
                # sample <- SP97681 in mutations_decile_replicating_timing_df
                # under vep_files_path
                # LICA-FR_SP97681.snv_mnv_filtered.txt
                # LICA-FR_SP97681.snv_mnv_filtered.txt_summary.html
                # LICA-FR_SP97681.indel_filtered.txt
                # LICA-FR_SP97681.indel_filtered.txt_summary.html
                if (cancer_type is not None) and ((cancer_type == 'Liver-HCC') or (cancer_type == 'Lung-AdenoCA')) and os.path.exists(vep_files_path):
                    files_list = os.listdir(vep_files_path)
                    filenames = [file[:-4] for file in files_list if (sample in file) and (file.endswith('.snv_mnv_filtered.txt'))]

                    if len(filenames) > 0:
                        filename = filenames[0]
                        sample_vep_filepath = os.path.join(vep_files_path, filename + '.txt')  # TODO change here for Liver and Lung cancer samples
                    else:
                        continue

                else:
                    # read sample vep file that comes from  vep call
                    sample_vep_filepath = os.path.join(vep_files_path, sample + '.txt')  # TODO change here for Liver and Lung cancer samples

                sample_vep_df = pd.read_csv(sample_vep_filepath, sep='\t', comment='#')
                sample_vep_df.columns = vep_columns

                # merge sample_decile_df and sample_vep_df
                merged_df = pd.merge(sample_decile_df, sample_vep_df, how='left', left_on=[CHR_START], right_on=['Location'])
                merged_grouped = merged_df.groupby('Consequence')

                # decile based and sample based: consequence -> number of consequences
                sample_consequences_dict = {consequence: len(consequence_df) for consequence, consequence_df in merged_grouped}

                # Breakdown consequences
                sample_new_consequences_dict = {}

                for k, v in sample_consequences_dict.items():
                    k_list = k.strip().split(',')
                    for knew in k_list:
                        if knew not in sample_new_consequences_dict:
                            sample_new_consequences_dict[knew] = v
                        else:
                            sample_new_consequences_dict[knew] += v

                assert sample not in all_samples_dict

                if sample not in all_samples_dict:
                    all_samples_dict[sample] = sample_new_consequences_dict

            # after all samples are visited for this decile index
            decile_dict[decile_index] = all_samples_dict

        # For all deciles df
        data_for_df_deciles = []
        df_decile_columns = ['Decile'] + ordered_all_consequences

        # Plot png figure for each decile
        for decile_index in range(1,range_max):
            decile_list = [decile_index]
            all_samples_dict = decile_dict[decile_index]

            # For each decile df
            data_for_df = []
            df_columns = ['Sample'] + ordered_all_consequences
            samples = all_samples_dict.keys()

            for sample in samples:
                sample_list= [sample]
                for consequence in ordered_all_consequences:
                    if consequence in all_samples_dict[sample]:
                        count = all_samples_dict[sample][consequence]
                    else:
                        count = 0
                    sample_list.append(count)
                data_for_df.append(sample_list)

            # Number of Mutations for each VEP consequence
            df = pd.DataFrame(data_for_df, columns=df_columns)

            filename = '%s_Samples_Decile_%s_VEP_Consequences.txt' % (signature, decile_index)
            filepath = os.path.join(outputDir, jobname, topography_vep_text_files, filename)
            df.to_csv(filepath, sep='\t', index=False)

            # Sum number of mutations for each consequence across all samples
            decile_list.extend(df.iloc[:,1:].sum(axis=0).values.tolist())
            data_for_df_deciles.append(decile_list)

            # Percentages for each VEP consequence
            df_perc = df
            df_perc.iloc[:,1:] = df.iloc[:, 1:].div(df.sum(axis=1), axis=0)
            df_perc.iloc[:,1:] = df_perc.iloc[:,1:] * 100
            filename = '%s_Samples_Decile_%s_VEP_Consequences_Percentages.txt' % (signature, decile_index)
            filepath = os.path.join(outputDir, jobname, topography_vep_text_files, filename)
            df_perc.to_csv(filepath, sep='\t', index=False)

            # Plot data in stack manner of bar type
            fig = df_perc.plot(x='Sample', kind='bar', stacked=True, figsize=(40, 30), color=my_colors, fontsize=20, title='Samples Mutations Consequences Percentages in Decile%s' %(decile_index)).get_figure()

            plt.legend(bbox_to_anchor=(1.01, 1.0), fontsize=20)

            filename = '%s_Samples_Consequence_Percentages_Decile%s.png' %(signature, decile_index)
            figureFile = os.path.join(outputDir, jobname, topography_vep_figures, filename)

            fig.savefig(figureFile)
            plt.close()

        # Aggregated Mutations
        # Dataframe Deciles Sum of Consequences Across All Samples Text File
        df_deciles = pd.DataFrame(data_for_df_deciles, columns=df_decile_columns)

        filename = '%s_Deciles_Across_All_Samples_VEP_Consequences.txt' %(signature)
        filepath = os.path.join(outputDir, jobname, topography_vep_text_files, filename)
        df_deciles.to_csv(filepath, sep='\t', index=False)

        # Plot Deciles on x axis, y axis Sum of Consequences Across All Samples Figure
        # x axis deciles, y vep consequences as stacked bar plots
        fig = df_deciles.plot(x='Decile', kind='bar', stacked=True, figsize=(40, 30), color=my_colors, fontsize=20, title='Decile Mutations VEP Consequences').get_figure()
        plt.legend(bbox_to_anchor=(1.01, 1.0), fontsize=20)

        filename = '%s_Deciles_Across_All_Samples_VEP_Consequences_Numbers.png' %(signature)
        figureFile = os.path.join(outputDir, jobname, topography_vep_figures, filename)

        fig.savefig(figureFile)
        plt.close()

        # Aggregated Mutations
        # Dataframe Deciles Percentage of Consequences Across All Samples Text file
        df_deciles_percentages = df_deciles
        df_deciles_percentages.iloc[:, 1:] = df_deciles.iloc[:, 1:].div(df_deciles.sum(axis=1), axis=0)
        df_deciles_percentages.iloc[:, 1:] = df_deciles_percentages.iloc[:, 1:] * 100

        filename = '%s_Deciles_Across_All_Samples_VEP_Consequences_Percentages.txt' %(signature)
        filepath = os.path.join(outputDir, jobname, topography_vep_text_files, filename)
        df_deciles_percentages.to_csv(filepath, sep='\t', index=False)

        # Dataframe Deciles Percentage of Consequences Across All Samples Figure
        # # x axis deciles, y percentages of vep consequences as stacked bar plots

        # from test
        # colormap = "tab20"
        # color=my_colors
        fig = df_deciles_percentages.plot(x='Decile', kind='bar', stacked=True, figsize=(40, 35), color = my_colors,
                                          fontsize=20,
                                          title='Decile Mutations VEP Consequences Percentages').get_figure()

        plt.xlabel('Early to Late Replication Timing', fontsize=30)
        plt.ylabel('Percentage of Mutations Annotations', fontsize=30)

        # put lower left of box to (x=0, y=1)
        # plt.gca().legend(labels, bbox_to_anchor=([0, 1]), loc='lower left', ncol=4, fontsize=30, frameon=False)
        plt.gca().legend(bbox_to_anchor=([0, 1]), loc='lower left', ncol=4, fontsize=30, frameon=False)

        filename = '%s_Deciles_Across_All_Samples_VEP_Consequences_Percentages.png' %(signature)
        figureFile = os.path.join(outputDir, jobname, topography_vep_figures, filename)

        fig.savefig(figureFile)
        plt.tight_layout()
        plt.close()
