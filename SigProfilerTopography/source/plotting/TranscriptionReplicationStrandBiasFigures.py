# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018-2020 Burcak Otlu

import os
import numpy as np
import statsmodels.stats.multitest

# import matplotlib
# BACKEND = 'Agg'
# if matplotlib.get_backend().lower() != BACKEND.lower():
#     # If backend is not set properly a call to describe will hang
#     matplotlib.use(BACKEND)

from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import gridspec


import pandas as pd

from SigProfilerTopography.source.commons.TopographyCommons import natural_key
from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIBED_STRAND
from SigProfilerTopography.source.commons.TopographyCommons import UNTRANSCRIBED_STRAND
from SigProfilerTopography.source.commons.TopographyCommons import LAGGING
from SigProfilerTopography.source.commons.TopographyCommons import LEADING

from SigProfilerTopography.source.commons.TopographyCommons import six_mutation_types
from SigProfilerTopography.source.commons.TopographyCommons import STRANDBIAS
from SigProfilerTopography.source.commons.TopographyCommons import DATA
from SigProfilerTopography.source.commons.TopographyCommons import FIGURE
from SigProfilerTopography.source.commons.TopographyCommons import SCATTER_PLOTS
from SigProfilerTopography.source.commons.TopographyCommons import BAR_PLOTS
from SigProfilerTopography.source.commons.TopographyCommons import CIRCLE_PLOTS
from SigProfilerTopography.source.commons.TopographyCommons import CIRCLE_BAR_PLOTS
from SigProfilerTopography.source.commons.TopographyCommons import SAMPLES
from SigProfilerTopography.source.commons.TopographyCommons import TABLES

from SigProfilerTopography.source.commons.TopographyCommons import SUBS_STRAND_BIAS_NUMBER_OF_MUTATIONS_THRESHOLD

from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIPTIONSTRANDBIAS
from SigProfilerTopography.source.commons.TopographyCommons import REPLICATIONSTRANDBIAS

from SigProfilerTopography.source.commons.TopographyCommons import LAGGING_VERSUS_LEADING
from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIBED_VERSUS_UNTRANSCRIBED
from SigProfilerTopography.source.commons.TopographyCommons import GENIC_VERSUS_INTERGENIC

from SigProfilerTopography.source.commons.TopographyCommons import LAGGING_VERSUS_LEADING_P_VALUE
from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIBED_VERSUS_UNTRANSCRIBED_P_VALUE
from SigProfilerTopography.source.commons.TopographyCommons import GENIC_VERSUS_INTERGENIC_P_VALUE

from SigProfilerTopography.source.commons.TopographyCommons import LAGGING_VERSUS_LEADING_Q_VALUE
from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIBED_VERSUS_UNTRANSCRIBED_Q_VALUE
from SigProfilerTopography.source.commons.TopographyCommons import GENIC_VERSUS_INTERGENIC_Q_VALUE

from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIBED_REAL_COUNT
from SigProfilerTopography.source.commons.TopographyCommons import UNTRANSCRIBED_REAL_COUNT

from SigProfilerTopography.source.commons.TopographyCommons import GENIC_REAL_COUNT
from SigProfilerTopography.source.commons.TopographyCommons import INTERGENIC_REAL_COUNT

from SigProfilerTopography.source.commons.TopographyCommons import LAGGING_REAL_COUNT
from SigProfilerTopography.source.commons.TopographyCommons import LEADING_REAL_COUNT

from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIBED_SIMULATIONS_MEAN_COUNT
from SigProfilerTopography.source.commons.TopographyCommons import UNTRANSCRIBED_SIMULATIONS_MEAN_COUNT

from SigProfilerTopography.source.commons.TopographyCommons import GENIC_SIMULATIONS_MEAN_COUNT
from SigProfilerTopography.source.commons.TopographyCommons import INTERGENIC_SIMULATIONS_MEAN_COUNT

from SigProfilerTopography.source.commons.TopographyCommons import LAGGING_SIMULATIONS_MEAN_COUNT
from SigProfilerTopography.source.commons.TopographyCommons import LEADING_SIMULATIONS_MEAN_COUNT

from SigProfilerTopography.source.commons.TopographyCommons import GENIC
from SigProfilerTopography.source.commons.TopographyCommons import INTERGENIC

from SigProfilerTopography.source.commons.TopographyCommons import percentage_numbers
from SigProfilerTopography.source.commons.TopographyCommons import percentage_strings

from SigProfilerTopography.source.commons.TopographyCommons import AT_LEAST_10_PERCENT_DIFF
from SigProfilerTopography.source.commons.TopographyCommons import AT_LEAST_20_PERCENT_DIFF
from SigProfilerTopography.source.commons.TopographyCommons import AT_LEAST_30_PERCENT_DIFF
from SigProfilerTopography.source.commons.TopographyCommons import AT_LEAST_50_PERCENT_DIFF
from SigProfilerTopography.source.commons.TopographyCommons import AT_LEAST_75_PERCENT_DIFF
from SigProfilerTopography.source.commons.TopographyCommons import AT_LEAST_100_PERCENT_DIFF

from SigProfilerTopography.source.commons.TopographyCommons import ID
from SigProfilerTopography.source.commons.TopographyCommons import DBS
from SigProfilerTopography.source.commons.TopographyCommons import SBS_CONTEXTS

from SigProfilerTopography.source.commons.TopographyCommons import PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_TOOL
from SigProfilerTopography.source.commons.TopographyCommons import PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_MANUSCRIPT

from SigProfilerTopography.source.commons.TopographyCommons import EXCEL_FILES
from SigProfilerTopography.source.commons.TopographyCommons import write_excel_file

from SigProfilerTopography.source.commons.TopographyCommons import NUMBER_OF_REQUIRED_MUTATIONS_FOR_STRAND_BIAS_BAR_PLOT

SIGNATURE='signature'
CANCER_TYPE='cancer_type'
MUTATION_TYPE='mutation_type'
TYPE = 'type'
SIGNIFICANT_STRAND='significant_strand'

SIGNIFICANCE_LEVEL=0.05

from SigProfilerTopography.source.commons.TopographyCommons import Table_SubsSignature_Cutoff_NumberofMutations_AverageProbability_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Table_IndelsSignature_Cutoff_NumberofMutations_AverageProbability_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Table_DinucsSignature_Cutoff_NumberofMutations_AverageProbability_Filename

from SigProfilerTopography.source.commons.TopographyCommons import Table_SubsSignature_NumberofMutations_AverageProbability_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Table_DinucsSignature_NumberofMutations_AverageProbability_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Table_IndelsSignature_NumberofMutations_AverageProbability_Filename

from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofSubsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofIndelsDict
from SigProfilerTopography.source.commons.TopographyCommons import Sample2NumberofDinucsDictFilename

from SigProfilerTopography.source.commons.TopographyCommons import getSample2SubsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2IndelsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import Sample2DinucsSignature2NumberofMutationsDictFilename

transcriptionStrands = [TRANSCRIBED_STRAND, UNTRANSCRIBED_STRAND]
genicVersusIntergenicStrands=[GENIC, INTERGENIC]
replicationStrands = [LAGGING, LEADING]

########################################################################
#New way
#For Mutation Types
def plot_mutation_types_transcription_log10_ratio_replication_log_10_ratio_using_dataframes(sample,numberofMutations,
                                                                                        type_transcribed_versus_untranscribed_df,
                                                                                        type_lagging_versus_leading_df,
                                                                                        outputDir, jobname):


    fig = plt.figure(figsize=(8,8), facecolor=None)
    plt.style.use('ggplot')

    # build a rectangle in axes coords
    left, width = .0, 1.
    bottom, height = .0, 1.
    right = left + width
    top = bottom + height

    # This code makes the background white.
    # Always put these statements after plt.figure
    ax = plt.gca()
    ax.set_facecolor('white')
    for edge_i in ['bottom','top','left','right']:
        ax.spines[edge_i].set_edgecolor("black")
        ax.spines[edge_i].set_linewidth(1)
        ax.spines[edge_i].set_bounds(-0.3, 0.3)

    plt.text(0.05, 1.02, 'Leading', ha='center', va='center', transform=ax.transAxes)
    plt.text(0.95, 1.02, 'Lagging', ha='center', va='center', transform=ax.transAxes)
    plt.text((right+0.02),(bottom+top-0.08), 'Transcribed',horizontalalignment='center',verticalalignment='center',rotation='vertical',transform=ax.transAxes)
    plt.text((right+0.02),(bottom+0.1), 'Untranscribed', horizontalalignment='center',verticalalignment='center', rotation='vertical', transform=ax.transAxes)

    if (sample is not None):
        plt.title(sample, fontsize=15, fontweight='bold')

    plt.xlabel('Lagging/leading replication strand\nratio(log10)', fontstyle='normal', fontsize=12, fontweight='bold')
    plt.ylabel('Transcribed/untranscribed strand\nratio(log10)',fontstyle='normal', fontsize=12, fontweight='bold')

    # Put some extra place by xlim if necessary
    plt.xlim(-0.3, 0.3)
    plt.ylim(-0.3, 0.3)

    plt.tick_params(axis='y', which='major', labelsize=10, width=1, length=10)
    plt.tick_params(axis='y', which='minor', labelsize=10, width=1, length=10)
    plt.tick_params(axis='x', which='major', labelsize=10, width=1, length=10)
    plt.tick_params(axis='x', which='minor', labelsize=10, width=1, length=10)

    yticks = [-0.2, -0.1, 0.0, 0.1, 0.2]
    yticklabels = ['-0.2', '-0.1', '0.0', '0.1', '0.2']
    plt.yticks(yticks, yticklabels)

    xticks = [-0.2, -0.1, 0.0, 0.1, 0.2]
    xticklabels = ['-0.2', '-0.1', '0.0', '0.1', '0.2']
    plt.xticks(xticks, xticklabels)

    # type_transcribed_versus_untranscribed_df=type_transcribed_versus_untranscribed_df[['cancer_type', 'type',
    #     'Transcribed_real_count', 'UnTranscribed_real_count', 'Transcribed_mean_sims_count', 'UnTranscribed_mean_sims_count', 'transcribed_versus_untranscribed_p_value','transcribed_versus_untranscribed_q_value',
    #     'Transcribed_real_count.1', 'Transcribed_mean_sims_count.1', 'Transcribed_min_sims_count', 'Transcribed_max_sims_count', 'Transcribed_sims_count_list',
    #     'UnTranscribed_real_count.1', 'UnTranscribed_mean_sims_count.1', 'UnTranscribed_min_sims_count', 'UnTranscribed_max_sims_count', 'UnTranscribed_sims_count_list' ]]
    #
    # type_lagging_versus_leading_df=type_lagging_versus_leading_df[['cancer_type', 'type',
    #     'Lagging_real_count', 'Leading_real_count', 'Lagging_mean_sims_count', 'Leading_mean_sims_count', 'lagging_versus_leading_p_value', 'lagging_versus_leading_q_value',
    #     'Lagging_real_count.1', 'Lagging_mean_sims_count.1', 'Lagging_min_sims_count', 'Lagging_max_sims_count', 'Lagging_sims_count_list',
    #     'Leading_real_count.1', 'Leading_mean_sims_count.1', 'Leading_min_sims_count', 'Leading_max_sims_count', 'Leading_sims_count_list' ]]

    ########################################################################
    transcriptionRatiosDict = {}
    replicationRatiosDict = {}
    for mutationType in six_mutation_types:

        ##################################################################
        transcribed_real_count=0
        untranscribed_real_count=0

        if (type_transcribed_versus_untranscribed_df[type_transcribed_versus_untranscribed_df['type']==mutationType]['Transcribed_real_count'].values.size>0):
            transcribed_real_count= type_transcribed_versus_untranscribed_df[type_transcribed_versus_untranscribed_df['type'] == mutationType]['Transcribed_real_count'].values[0]

        if (type_transcribed_versus_untranscribed_df[type_transcribed_versus_untranscribed_df['type']==mutationType]['UnTranscribed_real_count'].values.size>0):
            untranscribed_real_count= type_transcribed_versus_untranscribed_df[type_transcribed_versus_untranscribed_df['type'] == mutationType]['UnTranscribed_real_count'].values[0]

        if (transcribed_real_count>0 and untranscribed_real_count>0):
            transcriptionRatiosDict[mutationType] = np.log10(transcribed_real_count/untranscribed_real_count)
        ##################################################################

        ##################################################################
        lagging_real_count = 0
        leading_real_count = 0

        if (type_lagging_versus_leading_df[type_lagging_versus_leading_df['type'] == mutationType]['Lagging_real_count'].values.size > 0):
            lagging_real_count = type_lagging_versus_leading_df[type_lagging_versus_leading_df['type'] == mutationType]['Lagging_real_count'].values[0]

        if (type_lagging_versus_leading_df[type_lagging_versus_leading_df['type'] == mutationType]['Leading_real_count'].values.size > 0):
            leading_real_count = type_lagging_versus_leading_df[type_lagging_versus_leading_df['type'] == mutationType]['Leading_real_count'].values[0]

        if (lagging_real_count>0 and leading_real_count>0):
            replicationRatiosDict[mutationType] = np.log10(lagging_real_count/leading_real_count)
        ##################################################################

        ##################################################################
        if (mutationType in replicationRatiosDict) and (mutationType in transcriptionRatiosDict):
            plt.scatter(replicationRatiosDict[mutationType], transcriptionRatiosDict[mutationType], label=mutationType)
        ##################################################################

    ########################################################################

    legend = plt.legend(loc='upper left', frameon=True, fancybox =False,labels=six_mutation_types, bbox_to_anchor=(-0.0095, 1.0095))
    legend.get_frame().set_linewidth(1)

    frame = legend.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('black')

    plt.axvline(x=0.0, color='gray', linestyle='--')
    plt.axhline(y=0.0, color='gray', linestyle='--')

    if sample is None:
        figureName = 'all_mutation_types_%s_scatter_plot.png' %(STRANDBIAS)
        figureFile = os.path.join(outputDir, jobname, FIGURE, STRANDBIAS,SCATTER_PLOTS,figureName)

    else:
        figureName = 'all_mutation_types_%s_%d_%s_scatter_plot.png' %(sample,numberofMutations,STRANDBIAS)
        os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, STRANDBIAS,SCATTER_PLOTS), exist_ok=True)
        figureFile = os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, STRANDBIAS, SCATTER_PLOTS, figureName)

    fig.savefig(figureFile)
    plt.cla()
    plt.close(fig)

########################################################################

########################################################################
#Old way
#For Mutation Types
def plot_ncomms11383_Supp_FigG_AllMutationTypes_TranscriptionLog10Ratio_ReplicationLog10Ratio(sample,numberofMutations,type2TranscriptionStrand2CountDict,type2ReplicationStrand2CountDict,outputDir,jobname):

    fig = plt.figure(figsize=(8,8), facecolor=None)
    plt.style.use('ggplot')

    # build a rectangle in axes coords
    left, width = .0, 1.
    bottom, height = .0, 1.
    right = left + width
    top = bottom + height

    # This code makes the background white.
    # Always put these statements after plt.figure
    ax = plt.gca()
    ax.set_facecolor('white')
    for edge_i in ['bottom','top','left','right']:
        ax.spines[edge_i].set_edgecolor("black")
        ax.spines[edge_i].set_linewidth(1)
        ax.spines[edge_i].set_bounds(-0.3, 0.3)

    plt.text(0.05, 1.02, 'Leading', ha='center', va='center', transform=ax.transAxes)
    plt.text(0.95, 1.02, 'Lagging', ha='center', va='center', transform=ax.transAxes)
    plt.text((right+0.02),(bottom+top-0.08), 'Transcribed',horizontalalignment='center',verticalalignment='center',rotation='vertical',transform=ax.transAxes)
    plt.text((right+0.02),(bottom+0.1), 'Untranscribed', horizontalalignment='center',verticalalignment='center', rotation='vertical', transform=ax.transAxes)

    if (sample is not None):
        plt.title(sample, fontsize=15, fontweight='bold')

    plt.xlabel('Lagging/leading replication strand\nratio(log10)', fontstyle='normal', fontsize=12, fontweight='bold')
    plt.ylabel('Transcribed/untranscribed strand\nratio(log10)',fontstyle='normal', fontsize=12, fontweight='bold')

    # Put some extra place by xlim if necessary
    plt.xlim(-0.3, 0.3)
    plt.ylim(-0.3, 0.3)

    plt.tick_params(axis='y', which='major', labelsize=10, width=1, length=10)
    plt.tick_params(axis='y', which='minor', labelsize=10, width=1, length=10)
    plt.tick_params(axis='x', which='major', labelsize=10, width=1, length=10)
    plt.tick_params(axis='x', which='minor', labelsize=10, width=1, length=10)

    # plt.tick_params(
    #     axis='y',  # changes apply to the x-axis
    #     which='both',  # both major and minor ticks are affected
    #     left='off'  # ticks along the bottom edge are off
    # )

    yticks = [-0.2, -0.1, 0.0, 0.1, 0.2]
    yticklabels = ['-0.2', '-0.1', '0.0', '0.1', '0.2']
    plt.yticks(yticks, yticklabels)

    xticks = [-0.2, -0.1, 0.0, 0.1, 0.2]
    xticklabels = ['-0.2', '-0.1', '0.0', '0.1', '0.2']
    plt.xticks(xticks, xticklabels)

    ########################################################################
    transcriptionRatiosDict = {}
    replicationRatiosDict = {}
    for mutationType in six_mutation_types:
        if (mutationType in type2TranscriptionStrand2CountDict) and (mutationType in type2ReplicationStrand2CountDict):

            if ((TRANSCRIBED_STRAND in type2TranscriptionStrand2CountDict[mutationType]) and (UNTRANSCRIBED_STRAND in type2TranscriptionStrand2CountDict[mutationType])):
                transcriptionRatiosDict[mutationType]= np.log10(type2TranscriptionStrand2CountDict[mutationType][TRANSCRIBED_STRAND]/type2TranscriptionStrand2CountDict[mutationType][UNTRANSCRIBED_STRAND])

            if ((LAGGING in type2ReplicationStrand2CountDict[mutationType]) and (LEADING in type2ReplicationStrand2CountDict[mutationType])):
                replicationRatiosDict[mutationType] = np.log10(type2ReplicationStrand2CountDict[mutationType][LAGGING]/type2ReplicationStrand2CountDict[mutationType][LEADING])

            if (mutationType in replicationRatiosDict) and (mutationType in transcriptionRatiosDict):
                plt.scatter(replicationRatiosDict[mutationType],transcriptionRatiosDict[mutationType], label=mutationType)
    ########################################################################

    legend = plt.legend(loc='upper left', frameon=True, fancybox =False,labels=six_mutation_types, bbox_to_anchor=(-0.0095, 1.0095))
    legend.get_frame().set_linewidth(1)

    frame = legend.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('black')

    plt.axvline(x=0.0, color='gray', linestyle='--')
    plt.axhline(y=0.0, color='gray', linestyle='--')

    if sample is None:
        figureName = 'all_mutation_types_%s_scatter_plot.png' %(STRANDBIAS)
        figureFile = os.path.join(outputDir, jobname, FIGURE, STRANDBIAS,SCATTER_PLOTS,figureName)

    else:
        figureName = 'all_mutation_types_%s_%d_%s_scatter_plot.png' %(sample,numberofMutations,STRANDBIAS)
        os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, STRANDBIAS,SCATTER_PLOTS), exist_ok=True)
        figureFile = os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, STRANDBIAS, SCATTER_PLOTS, figureName)

    fig.savefig(figureFile)
    plt.cla()
    plt.close(fig)
########################################################################

########################################################################
#July 7, 2020
def plot_types_transcription_log10_ratio_replication_log10_ratio_using_dataframes(signatureType,
        sample,
        numberofMutations,
        type_transcribed_versus_untranscribed_df,
        type_lagging_versus_leading_df,
        signature_cutoff_numberofmutations_averageprobability_df,
        outputDir,
        jobname):

    fig = plt.figure(figsize=(8,8), facecolor=None)
    plt.style.use('ggplot')

    # build a rectangle in axes coords
    left, width = .0, 1.
    bottom, height = .0, 1.
    right = left + width
    top = bottom + height

    # This code makes the background white.
    # Always put these statements after plt.figure
    ax = plt.gca()
    ax.set_facecolor('white')
    for edge_i in ['bottom','top','left','right']:
        ax.spines[edge_i].set_edgecolor("black")
        ax.spines[edge_i].set_linewidth(1)
        ax.spines[edge_i].set_bounds(-0.3, 0.3)

    plt.text(0.05, 1.02, 'Leading', ha='center', va='center', transform=ax.transAxes)
    plt.text(0.95, 1.02, 'Lagging', ha='center', va='center', transform=ax.transAxes)
    plt.text((right+0.02),(bottom+top-0.08), 'Transcribed',horizontalalignment='center',verticalalignment='center',rotation='vertical',transform=ax.transAxes)
    plt.text((right+0.02),(bottom+0.1), 'Untranscribed', horizontalalignment='center',verticalalignment='center', rotation='vertical', transform=ax.transAxes)

    if (sample is not None):
        plt.title(sample, fontsize=15, fontweight='bold')

    plt.xlabel('Lagging/leading replication strand\nratio(log10)', fontstyle='normal', fontsize=12, fontweight='bold')
    plt.ylabel('Transcribed/untranscribed strand\nratio(log10)', fontstyle='normal', fontsize=12, fontweight='bold')

    # Put some extra place by xlim if necessary
    plt.xlim(-0.3, 0.3)
    plt.ylim(-0.3, 0.3)

    plt.tick_params(axis='y', which='major', labelsize=10, width=1, length=10)
    plt.tick_params(axis='y', which='minor', labelsize=10, width=1, length=10)
    plt.tick_params(axis='x', which='major', labelsize=10, width=1, length=10)
    plt.tick_params(axis='x', which='minor', labelsize=10, width=1, length=10)

    yticks = [-0.2, -0.1, 0.0, 0.1, 0.2]
    yticklabels = ['-0.2', '-0.1', '0.0', '0.1', '0.2']
    plt.yticks(yticks, yticklabels)

    xticks = [-0.2, -0.1, 0.0, 0.1, 0.2]
    xticklabels = ['-0.2', '-0.1', '0.0', '0.1', '0.2']
    plt.xticks(xticks, xticklabels)

    transcriptionRatiosDict = {}
    replicationRatiosDict = {}

    for signature in signature_cutoff_numberofmutations_averageprobability_df['signature'].unique():
        #################################################################################################
        #First check whether we have this signature or not
        # type_transcribed_versus_untranscribed_df=type_transcribed_versus_untranscribed_df[['cancer_type', 'type',
        #     'Transcribed_real_count', 'UnTranscribed_real_count', 'Transcribed_mean_sims_count', 'UnTranscribed_mean_sims_count', 'transcribed_versus_untranscribed_p_value','transcribed_versus_untranscribed_q_value',
        #     'Transcribed_real_count.1', 'Transcribed_mean_sims_count.1', 'Transcribed_min_sims_count', 'Transcribed_max_sims_count', 'Transcribed_sims_count_list',
        #     'UnTranscribed_real_count.1', 'UnTranscribed_mean_sims_count.1', 'UnTranscribed_min_sims_count', 'UnTranscribed_max_sims_count', 'UnTranscribed_sims_count_list' ]]

        transcribed_real_count=0
        untranscribed_real_count=0

        if (type_transcribed_versus_untranscribed_df[type_transcribed_versus_untranscribed_df['type']==signature]['Transcribed_real_count'].values.size>0):
            transcribed_real_count=type_transcribed_versus_untranscribed_df[type_transcribed_versus_untranscribed_df['type'] == signature]['Transcribed_real_count'].values[0]
        if (type_transcribed_versus_untranscribed_df[type_transcribed_versus_untranscribed_df['type']==signature]['UnTranscribed_real_count'].values.size>0):
            untranscribed_real_count=type_transcribed_versus_untranscribed_df[type_transcribed_versus_untranscribed_df['type'] == signature]['UnTranscribed_real_count'].values[0]

        if (transcribed_real_count+untranscribed_real_count>=SUBS_STRAND_BIAS_NUMBER_OF_MUTATIONS_THRESHOLD):
            transcriptionRatiosDict[signature] = np.log10(transcribed_real_count/untranscribed_real_count)
        #################################################################################################

        #################################################################################################
        # First check whether we have this signature or not
        # type_lagging_versus_leading_df=type_lagging_versus_leading_df[['cancer_type', 'type',
        #     'Lagging_real_count', 'Leading_real_count', 'Lagging_mean_sims_count', 'Leading_mean_sims_count', 'lagging_versus_leading_p_value', 'lagging_versus_leading_q_value',
        #     'Lagging_real_count.1', 'Lagging_mean_sims_count.1', 'Lagging_min_sims_count', 'Lagging_max_sims_count', 'Lagging_sims_count_list',
        #     'Leading_real_count.1', 'Leading_mean_sims_count.1', 'Leading_min_sims_count', 'Leading_max_sims_count', 'Leading_sims_count_list' ]]

        lagging_real_count=0
        leading_real_count = 0

        if (type_lagging_versus_leading_df[type_lagging_versus_leading_df['type']==signature]['Lagging_real_count'].values.size>0):
            lagging_real_count=type_lagging_versus_leading_df[type_lagging_versus_leading_df['type']==signature]['Lagging_real_count'].values[0]

        if (type_lagging_versus_leading_df[type_lagging_versus_leading_df['type']==signature]['Leading_real_count'].values.size>0):
            leading_real_count=type_lagging_versus_leading_df[type_lagging_versus_leading_df['type']==signature]['Leading_real_count'].values[0]

        if (lagging_real_count+leading_real_count>=SUBS_STRAND_BIAS_NUMBER_OF_MUTATIONS_THRESHOLD):
            replicationRatiosDict[signature] = np.log10(lagging_real_count/leading_real_count)
        #################################################################################################

    if (transcriptionRatiosDict and replicationRatiosDict):
        signaturesShownInLegend = []

        for signature in signature_cutoff_numberofmutations_averageprobability_df['signature'].unique():
            if ((signature in replicationRatiosDict.keys()) and (signature in transcriptionRatiosDict.keys())):
                signaturesShownInLegend.append(signature)
                plt.scatter(replicationRatiosDict[signature], transcriptionRatiosDict[signature], label=signature)

        legend = plt.legend(loc='upper left', frameon=True, fancybox=False, labels=signaturesShownInLegend,
                            bbox_to_anchor=(-0.0095, 1.0095))
        legend.get_frame().set_linewidth(1)

        frame = legend.get_frame()
        frame.set_facecolor('white')
        frame.set_edgecolor('black')

        plt.axvline(x=0.0, color='gray', linestyle='--')
        plt.axhline(y=0.0, color='gray', linestyle='--')

        if sample is None:
            figureName = 'all_%s_signatures_%s_scatter_plot.png' % (signatureType, STRANDBIAS)
            figureFile = os.path.join(outputDir, jobname, FIGURE, STRANDBIAS,SCATTER_PLOTS,figureName)
        else:
            figureName = 'all_%s_signatures_%s_%d_%s_scatter_plot.png' % (signatureType, sample, numberofMutations, STRANDBIAS)
            os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, STRANDBIAS,SCATTER_PLOTS), exist_ok=True)
            figureFile = os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, STRANDBIAS, SCATTER_PLOTS, figureName)

        fig.savefig(figureFile)
        plt.cla()
        plt.close(fig)
########################################################################

########################################################################
#May 9, 2018 starts
#For Signatures
def plot_ncomms11383_Supp_FigH_AllSignatures_TranscriptionLog10Ratio_ReplicationLog10Ratio(
        signatureType,
        sample,
        numberofMutations,
        signature2TranscriptionStrand2CountDict,
        signature2ReplicationStrand2CountDict,
        signature_cutoff_numberofmutations_averageprobability_df,
        outputDir,
        jobname):

    fig = plt.figure(figsize=(8,8), facecolor=None)
    plt.style.use('ggplot')

    # build a rectangle in axes coords
    left, width = .0, 1.
    bottom, height = .0, 1.
    right = left + width
    top = bottom + height

    # This code makes the background white.
    # Always put these statements after plt.figure
    ax = plt.gca()
    ax.set_facecolor('white')
    for edge_i in ['bottom','top','left','right']:
        ax.spines[edge_i].set_edgecolor("black")
        ax.spines[edge_i].set_linewidth(1)
        ax.spines[edge_i].set_bounds(-0.3, 0.3)

    plt.text(0.05, 1.02, 'Leading', ha='center', va='center', transform=ax.transAxes)
    plt.text(0.95, 1.02, 'Lagging', ha='center', va='center', transform=ax.transAxes)
    plt.text((right+0.02),(bottom+top-0.08), 'Transcribed',horizontalalignment='center',verticalalignment='center',rotation='vertical',transform=ax.transAxes)
    plt.text((right+0.02),(bottom+0.1), 'Untranscribed', horizontalalignment='center',verticalalignment='center', rotation='vertical', transform=ax.transAxes)

    if (sample is not None):
        plt.title(sample, fontsize=15, fontweight='bold')

    plt.xlabel('Lagging/leading replication strand\nratio(log10)', fontstyle='normal', fontsize=12, fontweight='bold')
    plt.ylabel('Transcribed/untranscribed strand\nratio(log10)', fontstyle='normal', fontsize=12, fontweight='bold')

    # Put some extra place by xlim if necessary
    plt.xlim(-0.3, 0.3)
    plt.ylim(-0.3, 0.3)

    plt.tick_params(axis='y', which='major', labelsize=10, width=1, length=10)
    plt.tick_params(axis='y', which='minor', labelsize=10, width=1, length=10)
    plt.tick_params(axis='x', which='major', labelsize=10, width=1, length=10)
    plt.tick_params(axis='x', which='minor', labelsize=10, width=1, length=10)

    yticks = [-0.2, -0.1, 0.0, 0.1, 0.2]
    yticklabels = ['-0.2', '-0.1', '0.0', '0.1', '0.2']
    plt.yticks(yticks, yticklabels)

    xticks = [-0.2, -0.1, 0.0, 0.1, 0.2]
    xticklabels = ['-0.2', '-0.1', '0.0', '0.1', '0.2']
    plt.xticks(xticks, xticklabels)

    transcriptionRatiosDict = {}
    replicationRatiosDict = {}

    for signature in signature_cutoff_numberofmutations_averageprobability_df['signature'].unique():
        #################################################################################################
        #First check whether we have this signature or not
        if ((signature in signature2TranscriptionStrand2CountDict) and (TRANSCRIBED_STRAND in (signature2TranscriptionStrand2CountDict[signature])) and
             (UNTRANSCRIBED_STRAND in (signature2TranscriptionStrand2CountDict[signature])) ):

            if ((signature2TranscriptionStrand2CountDict[signature][TRANSCRIBED_STRAND]+signature2TranscriptionStrand2CountDict[signature][UNTRANSCRIBED_STRAND]) >= SUBS_STRAND_BIAS_NUMBER_OF_MUTATIONS_THRESHOLD):
                transcriptionRatiosDict[signature]= np.log10(signature2TranscriptionStrand2CountDict[signature][TRANSCRIBED_STRAND]/signature2TranscriptionStrand2CountDict[signature][UNTRANSCRIBED_STRAND])
        #################################################################################################

        #################################################################################################
        # First check whether we have this signature or not
        if ((signature in signature2ReplicationStrand2CountDict) and (LAGGING in (signature2ReplicationStrand2CountDict[signature])) and
                (LEADING in (signature2ReplicationStrand2CountDict[signature]))):

            if ((signature2ReplicationStrand2CountDict[signature][LAGGING]+signature2ReplicationStrand2CountDict[signature][LEADING])>= SUBS_STRAND_BIAS_NUMBER_OF_MUTATIONS_THRESHOLD):
                replicationRatiosDict[signature] = np.log10(signature2ReplicationStrand2CountDict[signature][LAGGING]/signature2ReplicationStrand2CountDict[signature][LEADING])
        #################################################################################################

    if (transcriptionRatiosDict and replicationRatiosDict):
        signaturesShownInLegend = []

        for signature in signature_cutoff_numberofmutations_averageprobability_df['signature'].unique():
            if ((signature in replicationRatiosDict.keys()) and (signature in transcriptionRatiosDict.keys())):
                signaturesShownInLegend.append(signature)
                plt.scatter(replicationRatiosDict[signature], transcriptionRatiosDict[signature], label=signature)

        legend = plt.legend(loc='upper left', frameon=True, fancybox=False, labels=signaturesShownInLegend,
                            bbox_to_anchor=(-0.0095, 1.0095))
        legend.get_frame().set_linewidth(1)

        frame = legend.get_frame()
        frame.set_facecolor('white')
        frame.set_edgecolor('black')

        plt.axvline(x=0.0, color='gray', linestyle='--')
        plt.axhline(y=0.0, color='gray', linestyle='--')

        if sample is None:
            figureName = 'all_%s_signatures_%s_scatter_plot.png' % (signatureType, STRANDBIAS)
            figureFile = os.path.join(outputDir, jobname, FIGURE, STRANDBIAS,SCATTER_PLOTS,figureName)
        else:
            figureName = 'all_%s_signatures_%s_%d_%s_scatter_plot.png' % (
            signatureType, sample, numberofMutations, STRANDBIAS)
            os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, STRANDBIAS,SCATTER_PLOTS), exist_ok=True)
            figureFile = os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, STRANDBIAS, SCATTER_PLOTS, figureName)

        fig.savefig(figureFile)
        plt.cla()
        plt.close(fig)
########################################################################

########################################################################
#MutationTypeBased SampleBased Figures
def plot_ncomms11383_Supp_FigE_MutationTypeBased_AllSamples_TranscriptionLog10Ratio_ReplicationLog10Ratio(
        type2Sample2TranscriptionStrand2CountDict,
        type2Sample2ReplicationStrand2CountDict,
        outputDir,
        jobname,
        isFigureAugmentation):

    mutationType2ColorDict = {'C>A': 'blue', 'C>G':'black', 'C>T':'red', 'T>A':'gray', 'T>C':'green', 'T>G':'pink'}

    transcriptionRatiosDict = {}
    replicationRatiosDict = {}
    for mutationType in six_mutation_types:
        #initialization
        if mutationType not in transcriptionRatiosDict:
            transcriptionRatiosDict[mutationType] = {}
        if mutationType not in replicationRatiosDict:
            replicationRatiosDict[mutationType] = {}

        #Fill the dictionaries
        if mutationType in type2Sample2TranscriptionStrand2CountDict:
            for sample in type2Sample2TranscriptionStrand2CountDict[mutationType].keys():
                if ((TRANSCRIBED_STRAND in type2Sample2TranscriptionStrand2CountDict[mutationType][sample].keys()) and (UNTRANSCRIBED_STRAND in type2Sample2TranscriptionStrand2CountDict[mutationType][sample].keys())):
                    transcriptionRatiosDict[mutationType][sample]= np.log10(type2Sample2TranscriptionStrand2CountDict[mutationType][sample][TRANSCRIBED_STRAND]/type2Sample2TranscriptionStrand2CountDict[mutationType][sample][UNTRANSCRIBED_STRAND])

        if mutationType in type2Sample2ReplicationStrand2CountDict:
            for sample in type2Sample2ReplicationStrand2CountDict[mutationType].keys():
                if ((LAGGING in type2Sample2ReplicationStrand2CountDict[mutationType][sample].keys()) and (LEADING in type2Sample2ReplicationStrand2CountDict[mutationType][sample].keys())):
                    replicationRatiosDict[mutationType][sample] = np.log10(type2Sample2ReplicationStrand2CountDict[mutationType][sample][LAGGING]/type2Sample2ReplicationStrand2CountDict[mutationType][sample][LEADING])

    for mutationType in six_mutation_types:
        fig = plt.figure(figsize=(8, 8), facecolor=None)
        plt.style.use('ggplot')

        # build a rectangle in axes coords
        left, width = .0, 1.
        bottom, height = .0, 1.
        right = left + width
        top = bottom + height

        # This code makes the background white.
        # Always put these statements after plt.figure
        ax = plt.gca()
        ax.set_facecolor('white')
        for edge_i in ['bottom', 'top', 'left', 'right']:
            ax.spines[edge_i].set_edgecolor("black")
            ax.spines[edge_i].set_linewidth(1)
            ax.spines[edge_i].set_bounds(-0.65, 0.65)

        plt.title(mutationType, fontsize=15, fontweight='bold')

        plt.text(0.05, 1.02, 'Leading', ha='center', va='center', transform=ax.transAxes)
        plt.text(0.95, 1.02, 'Lagging', ha='center', va='center', transform=ax.transAxes)
        plt.text((right + 0.02), (bottom + top - 0.08), 'Transcribed', horizontalalignment='center',verticalalignment='center', rotation='vertical', transform=ax.transAxes)
        plt.text((right + 0.02), (bottom + 0.1), 'Untranscribed', horizontalalignment='center',verticalalignment='center', rotation='vertical', transform=ax.transAxes)

        plt.xlabel('Lagging/leading replication strand\nratio(log10)', fontstyle='normal', fontsize=12,fontweight='bold')
        plt.ylabel('Transcribed/untranscribed strand\nratio(log10)', fontstyle='normal', fontsize=12,fontweight='bold')

        # Put some extra place by xlim if necessary
        plt.xlim(-0.65, 0.65)
        plt.ylim(-0.65, 0.65)

        plt.tick_params(axis='y', which='major', labelsize=10, width=1, length=10)
        plt.tick_params(axis='y', which='minor', labelsize=10, width=1, length=10)
        plt.tick_params(axis='x', which='major', labelsize=10, width=1, length=10)
        plt.tick_params(axis='x', which='minor', labelsize=10, width=1, length=10)

        yticks = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6]
        yticklabels = ['-0.6', '-0.4', '-0.2', '0.0', '0.2', '0.4', '0.6']
        plt.yticks(yticks, yticklabels)

        xticks = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6]
        xticklabels = ['-0.6', '-0.4', '-0.2', '0.0', '0.2', '0.4', '0.6']
        plt.xticks(xticks, xticklabels)

        if (mutationType in type2Sample2TranscriptionStrand2CountDict):
            for sample in type2Sample2TranscriptionStrand2CountDict[mutationType].keys():
                if ((sample in replicationRatiosDict[mutationType].keys()) and (sample in transcriptionRatiosDict[mutationType].keys())):
                    plt.scatter(replicationRatiosDict[mutationType][sample],transcriptionRatiosDict[mutationType][sample], facecolor='none', color=mutationType2ColorDict[mutationType])

        plt.axvline(x=0.0, color='gray', linestyle='--')
        plt.axhline(y=0.0, color='gray', linestyle='--')

        if (isFigureAugmentation):
            plt.title(jobname + ' ' + mutationType)

        newMutationType = mutationType.replace('>', '2')

        figureName = newMutationType + '_MutationType_' + STRANDBIAS + '.png'
        figureFile = os.path.join(outputDir,jobname,FIGURE,STRANDBIAS,SCATTER_PLOTS,figureName)
        fig.savefig(figureFile)
        plt.cla()
        plt.close(fig)
########################################################################


########################################################################
#SignatureBased SampleBased Figures
#Sig26 is very different
def plot_ncomms11383_Supp_FigF_SignatureBased_AllSamples_TranscriptionLog10Ratio_ReplicationLog10Ratio(type2Sample2TranscriptionStrand2CountDict,type2Sample2ReplicationStrand2CountDict,signatures,outputDir,jobname,isFigureAugmentation):
    transcriptionRatiosDict = {}
    replicationRatiosDict = {}

    for signature in signatures:
        # initialization
        if signature not in transcriptionRatiosDict:
            transcriptionRatiosDict[signature] = {}
        if signature not in replicationRatiosDict:
            replicationRatiosDict[signature] = {}
        # Fill the dictionaries
        if signature in type2Sample2TranscriptionStrand2CountDict:
            for sample in type2Sample2TranscriptionStrand2CountDict[signature].keys():
                if (UNTRANSCRIBED_STRAND in type2Sample2TranscriptionStrand2CountDict[signature][sample]) and (TRANSCRIBED_STRAND in type2Sample2TranscriptionStrand2CountDict[signature][sample]):
                    transcriptionRatiosDict[signature][sample] = np.log10(type2Sample2TranscriptionStrand2CountDict[signature][sample][TRANSCRIBED_STRAND] /type2Sample2TranscriptionStrand2CountDict[signature][sample][UNTRANSCRIBED_STRAND])
                    # print(signature, sample)
                    # print(signature2Sample2TranscriptionStrand2CountDict[signature][sample][TRANSCRIBED_STRAND])
                    # print(signature2Sample2TranscriptionStrand2CountDict[signature][sample][UNTRANSCRIBED_STRAND])
                    # print(signature,sample,transcriptionRatiosDict[signature][sample])

        if signature in type2Sample2ReplicationStrand2CountDict:
            for sample in type2Sample2ReplicationStrand2CountDict[signature].keys():
                if (LAGGING in type2Sample2ReplicationStrand2CountDict[signature][sample]) and (LEADING in type2Sample2ReplicationStrand2CountDict[signature][sample]):
                    replicationRatiosDict[signature][sample] = np.log10(type2Sample2ReplicationStrand2CountDict[signature][sample][LAGGING] /type2Sample2ReplicationStrand2CountDict[signature][sample][LEADING])

    for signature in signatures:
        if (len(replicationRatiosDict[signature].keys())>0 and len(transcriptionRatiosDict[signature].keys())>0):
            fig = plt.figure(figsize=(8, 8), facecolor=None)
            plt.style.use('ggplot')

            # build a rectangle in axes coords
            left, width = .0, 1.
            bottom, height = .0, 1.
            right = left + width
            top = bottom + height

            # This code makes the background white.
            # Always put these statements after plt.figure
            ax = plt.gca()
            ax.set_facecolor('white')
            for edge_i in ['bottom', 'top', 'left', 'right']:
                ax.spines[edge_i].set_edgecolor("black")
                ax.spines[edge_i].set_linewidth(1)
                ax.spines[edge_i].set_bounds(-0.65, 0.65)

            plt.title(signature, fontsize=15, fontweight='bold')

            plt.text(0.05, 1.02, 'Leading', ha='center', va='center', transform=ax.transAxes)
            plt.text(0.95, 1.02, 'Lagging', ha='center', va='center', transform=ax.transAxes)
            plt.text((right + 0.02), (bottom + top - 0.08), 'Transcribed', horizontalalignment='center',verticalalignment='center', rotation='vertical', transform=ax.transAxes)
            plt.text((right + 0.02), (bottom + 0.1), 'Untranscribed', horizontalalignment='center',verticalalignment='center', rotation='vertical', transform=ax.transAxes)

            plt.xlabel('Lagging/leading replication strand\nratio(log10)', fontstyle='normal', fontsize=12,fontweight='bold')
            plt.ylabel('Transcribed/untranscribed strand\nratio(log10)', fontstyle='normal', fontsize=12,fontweight='bold')

            # Put some extra place by xlim if necessary
            plt.xlim(-0.65, 0.65)
            plt.ylim(-0.65, 0.65)

            plt.tick_params(axis='y', which='major', labelsize=10, width=1, length=10)
            plt.tick_params(axis='y', which='minor', labelsize=10, width=1, length=10)
            plt.tick_params(axis='x', which='major', labelsize=10, width=1, length=10)
            plt.tick_params(axis='x', which='minor', labelsize=10, width=1, length=10)

            yticks = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6]
            yticklabels = ['-0.6', '-0.4', '-0.2', '0.0', '0.2', '0.4', '0.6']
            plt.yticks(yticks, yticklabels)

            xticks = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6]
            xticklabels = ['-0.6', '-0.4', '-0.2', '0.0', '0.2', '0.4', '0.6']
            plt.xticks(xticks, xticklabels)

            for sample in type2Sample2TranscriptionStrand2CountDict[signature].keys():
                if (sample in replicationRatiosDict[signature]) and (sample in transcriptionRatiosDict[signature]):
                    plt.scatter(replicationRatiosDict[signature][sample], transcriptionRatiosDict[signature][sample],facecolor='none',color='green')

            plt.axvline(x=0.0, color='gray', linestyle='--')
            plt.axhline(y=0.0, color='gray', linestyle='--')

            if (isFigureAugmentation):
                plt.title(jobname + ' ' + signature)

            figureName = signature.replace(' ','') + '_Signature_' + STRANDBIAS + '.png'
            figureFile = os.path.join(outputDir,jobname,FIGURE,STRANDBIAS,SCATTER_PLOTS,figureName)
            fig.savefig(figureFile)
            plt.cla()
            plt.close(fig)
########################################################################


##################################################################
#Only this method supports simulations
#key can be a sample or a signature
def plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,key,isKeySample,numberofMutations,N,x_axis_labels,strand1_values,strand2_values,strand1_simulations_median_values ,strand2_simulations_median_values , fdr_bh_adjusted_pvalues, strand1Name, strand2Name, mutationsOrSignatures, color1, color2, figureName, width,plot_mode):
    #Here we can take into difference between strand1_values and strand2_values while deciding on significance

    from matplotlib import rcParams
    rcParams.update({'figure.autolayout': True})

    # the x locations for the groups
    ind = np.arange(N)

    fig, ax = plt.subplots(figsize=(16,10),dpi=300)

    legend=None
    rects1=None
    rects2=None
    rects3=None
    rects4=None

    rects1 = ax.bar(ind, strand1_values, width=width, edgecolor='black', color=color1)
    rects2 = ax.bar(ind + width, strand2_values, width=width, edgecolor='black', color=color2)

    if ((strand1_simulations_median_values is not None) and strand1_simulations_median_values):
        rects3 = ax.bar(ind+ 2*width, strand1_simulations_median_values, width=width, edgecolor='black', color=color1, hatch = '///')
    if ((strand2_simulations_median_values is not None) and strand2_simulations_median_values):
        rects4 = ax.bar(ind +3*width, strand2_simulations_median_values, width=width, edgecolor='black', color=color2, hatch = '///')


    # add some text for labels, title and axes ticks

    ###########################################################################
    if plot_mode==PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_TOOL:
        ax.tick_params(axis='x', labelsize=35)
        ax.tick_params(axis='y', labelsize=35)

        locs, labels = plt.yticks()
        ax.set_ylim(0, locs[-1] + 5000)

        ##############################
        # To make the bar width not too wide
        if len(ind) < 6:
            maxn = 6
            ax.set_xlim(-0.5, maxn - 0.5)
        ##############################

        #Set title
        if key is not None:
            ax.set_title('%s %s vs. %s %s' %(key,strand1Name,strand2Name,mutationsOrSignatures), fontsize=20,fontweight='bold')
        else:
            ax.set_title('%s vs. %s %s' %(strand1Name,strand2Name,mutationsOrSignatures), fontsize=20,fontweight='bold')

        #Set x tick labels
        if len(x_axis_labels) > 6:
            ax.set_xticklabels(x_axis_labels, fontsize=35, rotation=90)
        else:
            ax.set_xticklabels(x_axis_labels, fontsize=35)

        # Set the ylabel
        plt.ylabel('Number of single base substitutions', fontsize=35, fontweight='normal')

        # set the x axis tick locations
        if (numberofSimulations > 0):
            ax.set_xticks(ind + (3 * width) / 2)
            realStrand1Name = 'Real %s' % (strand1Name)
            realStrand2Name = 'Real %s' % (strand2Name)
            simulationsStrand1Name = 'Simulated %s' % (strand1Name)
            simulationsStrand2Name = 'Simulated %s' % (strand2Name)
            if ((rects1 is not None) and (rects2 is not None) and (rects3 is not None) and (rects4 is not None)):
                if ((len(rects1) > 0) and (len(rects2) > 0) and (len(rects3) > 0) and (len(rects4) > 0)):
                    legend = ax.legend((rects1[0], rects2[0], rects3[0], rects4[0]),(realStrand1Name, realStrand2Name, simulationsStrand1Name, simulationsStrand2Name),prop={'size': 25}, ncol=1, loc='best')
        else:
            # Old way with no simulations
            ax.set_xticks(ind + width / 2)
            if ((rects1 is not None) and (rects2 is not None)):
                if ((len(rects1) > 0) and (len(rects2) > 0)):
                    legend = ax.legend((rects1[0], rects2[0]), (strand1Name, strand2Name), prop={'size': 25}, ncol=1, loc='upper right')

    elif plot_mode==PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_MANUSCRIPT:
        # set axis ticks
        # ax.tick_params(axis='both', which='both', length=0)
        ax.tick_params(axis='x', which='both', length=0)
        ax.tick_params(axis='y', which='both', length=0)
        # set axis labels
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax.get_yticklabels(), visible=False)
        if (numberofSimulations > 0):
            realStrand1Name = 'Real %s' % (strand1Name)
            realStrand2Name = 'Real %s' % (strand2Name)
            simulationsStrand1Name = 'Simulated %s' % (strand1Name)
            simulationsStrand2Name = 'Simulated %s' % (strand2Name)

            if ((rects1 is not None) and (rects2 is not None) and (rects3 is not None) and (rects4 is not None)):
                if ((len(rects1) > 0) and (len(rects2) > 0) and (len(rects3) > 0) and (len(rects4) > 0)):
                    #For SigProfilerTopography Overview Figure
                    # legend = ax.legend((rects1[0], rects2[0], rects3[0], rects4[0]),(realStrand1Name, realStrand2Name, simulationsStrand1Name, simulationsStrand2Name),prop={'size': 35}, ncol=2, loc='best')
                    #For Replication Time and Strand Bias Overview Figure
                    legend = ax.legend((rects1[0], rects2[0], rects3[0], rects4[0]),(realStrand1Name, realStrand2Name, simulationsStrand1Name, simulationsStrand2Name),prop={'size': 30}, ncol=1, loc='best')

        else:
            if ((rects1 is not None) and (rects2 is not None)):
                if ((len(rects1) > 0) and (len(rects2) > 0)):
                    legend = ax.legend((rects1[0], rects2[0]), (strand1Name, strand2Name), prop={'size': 35},loc='upper right')
    ###########################################################################

    #To make the barplot background white
    ax.set_facecolor('white')
    #To makes spines black like a rectangle with black stroke
    ax.spines["bottom"].set_color('black')
    ax.spines["left"].set_color('black')
    ax.spines["top"].set_color('black')
    ax.spines["right"].set_color('black')

    if (legend is not None):
        frame = legend.get_frame()
        frame.set_facecolor('white')
        frame.set_edgecolor('black')

    #########################################################################################################
    #Add star above the bars for significant differences between the number of mutations on each strand starts
    # For each bar: Place a label
    if fdr_bh_adjusted_pvalues is not None:
        for fdr_bh_adjusted_pvalue, rect1, rect2 in zip(fdr_bh_adjusted_pvalues,rects1,rects2):
            # Get X and Y placement of label from rect.
            y_value = max(rect1.get_height(),rect2.get_height())
            x_value = rect1.get_x() + rect1.get_width()

            # Number of points between bar and label. Change to your liking.
            space = 3
            # Vertical alignment for positive values
            va = 'bottom'

            # If value of bar is negative: Place label below bar
            if y_value < 0:
                # Invert space to place label below
                space *= -1
                # Vertically align label at top
                va = 'top'

            # Use Y value as label and format number with one decimal place
            label = "{:.1f}".format(y_value)

            # Create annotation
            if ((fdr_bh_adjusted_pvalue is not None) and (fdr_bh_adjusted_pvalue)<=0.0001):
                plt.annotate(
                    '***',  # Use `label` as label
                    (x_value, y_value),  # Place label at end of the bar
                    xytext=(0, space),  # Vertically shift label by `space`
                    textcoords="offset points",  # Interpret `xytext` as offset in points
                    ha='center',  # Horizontally center label
                    va=va,
                    fontsize=20)  # Vertically align label differently for

            elif ((fdr_bh_adjusted_pvalue is not None) and (fdr_bh_adjusted_pvalue)<=0.001):
                plt.annotate(
                    '**',  # Use `label` as label
                    (x_value, y_value),  # Place label at end of the bar
                    xytext=(0, space),  # Vertically shift label by `space`
                    textcoords="offset points",  # Interpret `xytext` as offset in points
                    ha='center',  # Horizontally center label
                    va=va,
                    fontsize=20)  # Vertically align label differently for

            elif ((fdr_bh_adjusted_pvalue is not None) and (fdr_bh_adjusted_pvalue)<=0.05):
                plt.annotate(
                    '*',  # Use `label` as label
                    (x_value, y_value),  # Place label at end of the bar
                    xytext=(0, space),  # Vertically shift label by `space`
                    textcoords="offset points",  # Interpret `xytext` as offset in points
                    ha='center',  # Horizontally center label
                    va=va,
                    fontsize=20)  # Vertically align label differently for

            # positive and negative values.
    #Add star above the bars for significant differences between the number of mutations on each strand ends
    #########################################################################################################

    if (key is None):
        figureName = '%s_bar_plot.png' %(figureName)
        figureFile = os.path.join(outputDir, jobname, FIGURE, STRANDBIAS, BAR_PLOTS, figureName)
    elif (not isKeySample):
        figureName = '%s_%s_bar_plot.png' %(key,figureName)
        figureFile = os.path.join(outputDir, jobname, FIGURE, STRANDBIAS, BAR_PLOTS, figureName)
    else:
        figureName = '%s_%s_%d_bar_plot.png' %(figureName,key,numberofMutations)
        os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES, key, STRANDBIAS, BAR_PLOTS), exist_ok=True)
        figureFile = os.path.join(outputDir, jobname, FIGURE, SAMPLES, key, STRANDBIAS, BAR_PLOTS, figureName)

    fig.savefig(figureFile)
    plt.cla()
    plt.close(fig)
##################################################################


# June 2, 2021
def plot_circle_plot_in_given_axis(ax,
                               percentage_strings,
                               sbs_signature,
                               six_mutation_types,
                               xticklabels_list,
                               signature2mutation_type2strand2percentagedict):

    strand_bias_list=[LAGGING_VERSUS_LEADING, TRANSCRIBED_VERSUS_UNTRANSCRIBED, GENIC_VERSUS_INTERGENIC]

    # make aspect ratio square
    ax.set_aspect(1.0)

    # set title
    title = '%s Strand Bias' %(sbs_signature)
    ax.text(len(percentage_strings) * 3, len(strand_bias_list) + 2.5, title, horizontalalignment='center',fontsize=60, fontweight='bold', fontname='Arial')

    # Colors are from SigProfilerPlotting tool to be consistent
    colors = [[3 / 256, 189 / 256, 239 / 256],
              [1 / 256, 1 / 256, 1 / 256],
              [228 / 256, 41 / 256, 38 / 256],
              [203 / 256, 202 / 256, 202 / 256],
              [162 / 256, 207 / 256, 99 / 256],
              [236 / 256, 199 / 256, 197 / 256]]

    # Put rectangles
    x = 0

    for i in range(0, len(six_mutation_types), 1):
        ax.text((x + (len(percentage_strings) / 2) - 0.75), len(strand_bias_list) + 1.5, six_mutation_types[i],fontsize=55, fontweight='bold', fontname='Arial')
        ax.add_patch(plt.Rectangle((x + .0415, len(strand_bias_list) + 0.75), len(percentage_strings) - (2 * .0415), .5,facecolor=colors[i], clip_on=False))
        ax.add_patch(plt.Rectangle((x, 0), len(percentage_strings), len(strand_bias_list), facecolor=colors[i], zorder=0,alpha=0.25, edgecolor='grey'))
        x += len(percentage_strings)

    # CODE GOES HERE TO CENTER X-AXIS LABELS...
    ax.set_xlim([0, len(six_mutation_types) * len(percentage_strings)])
    ax.set_xticklabels([])
    ax.tick_params(axis='x', which='minor', length=0, labelsize=35)

    # major ticks
    ax.set_xticks(np.arange(0, len(six_mutation_types) * len(percentage_strings), 1))
    # minor ticks
    ax.set_xticks(np.arange(0, len(six_mutation_types) * len(percentage_strings), 1) + 0.5, minor=True)

    ax.set_xticklabels(xticklabels_list, minor=True)
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')

    ax.tick_params(
        axis='x',  # changes apply to the x-axis
        which='major',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False)  # labels along the bottom edge are off

    # CODE GOES HERE TO CENTER Y-AXIS LABELS...
    ax.set_ylim([0, len(strand_bias_list)])
    ax.set_yticklabels([])
    ax.tick_params(axis='y', which='minor', length=0, labelsize=40)

    # major ticks
    ax.set_yticks(np.arange(0, len(strand_bias_list), 1))
    # minor ticks
    ax.set_yticks(np.arange(0, len(strand_bias_list), 1) + 0.5, minor=True)
    ax.set_yticklabels(['', sbs_signature,''], minor=True)  # fontsize

    ax.tick_params(
        axis='y',  # changes apply to the x-axis
        which='major',  # both major and minor ticks are affected
        left=False)  # labels along the bottom edge are off

    # Gridlines based on major ticks
    ax.grid(which='major', color='black', zorder=3)

    # Put the legend
    legend_elements = [
        Line2D([0], [0], marker='o', color='white', label=GENIC, markerfacecolor='cyan', markersize=40),
        Line2D([0], [0], marker='o', color='white', label=INTERGENIC, markerfacecolor='gray', markersize=40),
        Line2D([0], [0], marker='o', color='white', label=TRANSCRIBED_STRAND, markerfacecolor='royalblue',markersize=40),
        Line2D([0], [0], marker='o', color='white', label=UNTRANSCRIBED_STRAND, markerfacecolor='yellowgreen',markersize=40),
        Line2D([0], [0], marker='o', color='white', label=LAGGING, markerfacecolor='indianred', markersize=40),
        Line2D([0], [0], marker='o', color='white', label=LEADING, markerfacecolor='goldenrod', markersize=40)]

    legend = ax.legend(handles=legend_elements, ncol=len(legend_elements), bbox_to_anchor=(0.5, 0), loc='upper center',fontsize=40)
    # legend.get_frame().set_linewidth(1)
    frame = legend.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('black')

    for percentage_diff_index, percentage_string in enumerate(percentage_strings):
        for mutation_type_index, mutation_type in enumerate(six_mutation_types):
            # for row_sbs_signature_index, row_sbs_signature in enumerate(rows_sbs_signatures):
            # strand_bias_list = [TRANSCRIBED_VERSUS_UNTRANSCRIBED, GENIC_VERSUS_INTERGENIC, LAGGING_VERSUS_LEADING]
            for strand_bias_index, strand_bias in enumerate(strand_bias_list):
                if (strand_bias == LAGGING_VERSUS_LEADING):
                    if sbs_signature in signature2mutation_type2strand2percentagedict:
                        if mutation_type in signature2mutation_type2strand2percentagedict[sbs_signature]:
                            lagging_percentage = None
                            leading_percentage = None

                            if (LAGGING in signature2mutation_type2strand2percentagedict[sbs_signature][mutation_type]) and (
                                    signature2mutation_type2strand2percentagedict[sbs_signature][mutation_type][LAGGING][percentage_string] == 1):
                                lagging_percentage = 100
                            if (LEADING in signature2mutation_type2strand2percentagedict[sbs_signature][mutation_type]) and (
                                    signature2mutation_type2strand2percentagedict[sbs_signature][mutation_type][LEADING][percentage_string] == 1):
                                leading_percentage = 100

                            if (lagging_percentage is not None) and (leading_percentage is None):
                                radius = 0.49
                                if (radius > 0):
                                    # print('Plot circle at x=%d y=%d for %s %s %s' % (mutation_type_index * len(percentage_strings) + percentage_diff_index, row_sbs_signature_index, row_sbs_signature,mutation_type, percentage_string))
                                    ax.add_patch(plt.Circle((mutation_type_index * len(percentage_strings) + percentage_diff_index + 0.5,strand_bias_index + 0.5), radius, color='indianred', fill=True))
                            elif (leading_percentage is not None) and (lagging_percentage is None):
                                radius = 0.49
                                if (radius > 0):
                                    # print('Plot circle at x=%d y=%d for %s %s %s' % (mutation_type_index * len(percentage_strings) + percentage_diff_index, row_sbs_signature_index, row_sbs_signature,mutation_type, percentage_string))
                                    ax.add_patch(plt.Circle((mutation_type_index * len(
                                        percentage_strings) + percentage_diff_index + 0.5,
                                                                 strand_bias_index + 0.5), radius,
                                                                color='goldenrod', fill=True))

                            elif (lagging_percentage is not None) and (leading_percentage is not None):
                                radius_lagging = 0.49
                                radius_leading = 0.49
                                if (radius_lagging > radius_leading):
                                    # First lagging
                                    ax.add_patch(plt.Circle((mutation_type_index * len(
                                        percentage_strings) + percentage_diff_index + 0.5,
                                                                 strand_bias_index + 0.5), radius_lagging,
                                                                color='indianred', fill=True))
                                    # Second leading
                                    ax.add_patch(plt.Circle((mutation_type_index * len(
                                        percentage_strings) + percentage_diff_index + 0.5,
                                                                 strand_bias_index + 0.5), radius_leading,
                                                                color='goldenrod', fill=True))
                                else:
                                    # First leading
                                    ax.add_patch(plt.Circle((mutation_type_index * len(
                                        percentage_strings) + percentage_diff_index + 0.5,
                                                                 strand_bias_index + 0.5), radius_leading,
                                                                color='goldenrod', fill=True))
                                    # Second lagging
                                    ax.add_patch(plt.Circle((mutation_type_index * len(
                                        percentage_strings) + percentage_diff_index + 0.5,
                                                                 strand_bias_index + 0.5), radius_lagging,
                                                                color='indianred', fill=True))

                elif (strand_bias == GENIC_VERSUS_INTERGENIC):
                    if sbs_signature in signature2mutation_type2strand2percentagedict:
                        if mutation_type in signature2mutation_type2strand2percentagedict[sbs_signature]:
                            genic_percentage = None
                            intergenic_percentage = None

                            if (GENIC in signature2mutation_type2strand2percentagedict[sbs_signature][
                                mutation_type]) and (
                                    signature2mutation_type2strand2percentagedict[sbs_signature][mutation_type][
                                        GENIC][percentage_string] == 1):
                                genic_percentage = 100
                            if (INTERGENIC in signature2mutation_type2strand2percentagedict[sbs_signature][
                                mutation_type]) and (
                                    signature2mutation_type2strand2percentagedict[sbs_signature][mutation_type][
                                        INTERGENIC][percentage_string] == 1):
                                intergenic_percentage = 100

                            if (genic_percentage is not None) and (intergenic_percentage is None):
                                radius = 0.49
                                if (radius > 0):
                                    # print('Plot circle at x=%d y=%d for %s %s %s' % (mutation_type_index * len(percentage_strings) + percentage_diff_index, row_sbs_signature_index,row_sbs_signature, mutation_type, percentage_string))
                                    ax.add_patch(plt.Circle((mutation_type_index * len(
                                        percentage_strings) + percentage_diff_index + 0.5,
                                                                 strand_bias_index + 0.5), radius, color='cyan',
                                                                fill=True))

                            elif (intergenic_percentage is not None) and (genic_percentage is None):
                                radius = 0.49
                                if (radius > 0):
                                    # print('Plot circle at x=%d y=%d for %s %s %s' % (mutation_type_index * len(percentage_strings) + percentage_diff_index, row_sbs_signature_index,row_sbs_signature, mutation_type, percentage_string))
                                    ax.add_patch(plt.Circle((mutation_type_index * len(
                                        percentage_strings) + percentage_diff_index + 0.5,
                                                                 strand_bias_index + 0.5), radius, color='gray',
                                                                fill=True))

                            elif (genic_percentage is not None) and (intergenic_percentage is not None):
                                radius_genic = 0.49
                                radius_intergenic = 0.49
                                if (radius_genic > radius_intergenic):
                                    # First genic
                                    ax.add_patch(plt.Circle((mutation_type_index * len(
                                        percentage_strings) + percentage_diff_index + 0.5,
                                                                 strand_bias_index + 0.5), radius_genic,
                                                                color='cyan', fill=True))
                                    # Second intergenic
                                    ax.add_patch(plt.Circle((mutation_type_index * len(
                                        percentage_strings) + percentage_diff_index + 0.5,
                                                                 strand_bias_index + 0.5), radius_intergenic,
                                                                color='gray', fill=True))

                                else:
                                    # First intergenic
                                    ax.add_patch(plt.Circle((mutation_type_index * len(
                                        percentage_strings) + percentage_diff_index + 0.5,
                                                                 strand_bias_index + 0.5), radius_intergenic,
                                                                color='gray', fill=True))
                                    # Second genic
                                    ax.add_patch(plt.Circle((mutation_type_index * len(
                                        percentage_strings) + percentage_diff_index + 0.5,
                                                                 strand_bias_index + 0.5), radius_genic,
                                                                color='cyan', fill=True))


                elif (strand_bias == TRANSCRIBED_VERSUS_UNTRANSCRIBED):
                    if sbs_signature in signature2mutation_type2strand2percentagedict:
                        if mutation_type in signature2mutation_type2strand2percentagedict[sbs_signature]:
                            transcribed_percentage = None
                            untranscribed_percentage = None

                            if (TRANSCRIBED_STRAND in signature2mutation_type2strand2percentagedict[sbs_signature][
                                mutation_type]) and (
                                    signature2mutation_type2strand2percentagedict[sbs_signature][mutation_type][
                                        TRANSCRIBED_STRAND][percentage_string] == 1):
                                transcribed_percentage = 100
                            if (UNTRANSCRIBED_STRAND in
                                signature2mutation_type2strand2percentagedict[sbs_signature][mutation_type]) and (
                                    signature2mutation_type2strand2percentagedict[sbs_signature][mutation_type][
                                        UNTRANSCRIBED_STRAND][percentage_string] == 1):
                                untranscribed_percentage = 100

                            if (transcribed_percentage is not None) and (untranscribed_percentage is None):
                                radius = 0.49
                                if (radius > 0):
                                    # print('Plot circle at x=%d y=%d for %s %s %s' % (mutation_type_index * len(percentage_strings) + percentage_diff_index, row_sbs_signature_index,row_sbs_signature, mutation_type, percentage_string))
                                    ax.add_patch(plt.Circle((mutation_type_index * len(
                                        percentage_strings) + percentage_diff_index + 0.5,
                                                                 strand_bias_index + 0.5), radius,
                                                                color='royalblue', fill=True))

                            elif (untranscribed_percentage is not None) and (transcribed_percentage is None):
                                radius = 0.49
                                if (radius > 0):
                                    # print('Plot circle at x=%d y=%d for %s %s %s' % (mutation_type_index * len(percentage_strings) + percentage_diff_index, row_sbs_signature_index,row_sbs_signature, mutation_type, percentage_string))
                                    ax.add_patch(plt.Circle((mutation_type_index * len(
                                        percentage_strings) + percentage_diff_index + 0.5,
                                                                 strand_bias_index + 0.5), radius,
                                                                color='yellowgreen', fill=True))

                            elif (transcribed_percentage is not None) and (untranscribed_percentage is not None):
                                radius_transcribed = 0.49
                                radius_untranscribed = 0.49
                                if (radius_transcribed > radius_untranscribed):
                                    # First transcribed
                                    ax.add_patch(plt.Circle((mutation_type_index * len(
                                        percentage_strings) + percentage_diff_index + 0.5,
                                                                 strand_bias_index + 0.5), radius_transcribed,
                                                                color='royalblue', fill=True))
                                    # Second untranscribed
                                    ax.add_patch(plt.Circle((mutation_type_index * len(
                                        percentage_strings) + percentage_diff_index + 0.5,
                                                                 strand_bias_index + 0.5), radius_untranscribed,
                                                                color='yellowgreen', fill=True))

                                else:
                                    # First untranscribed
                                    ax.add_patch(plt.Circle((mutation_type_index * len(
                                        percentage_strings) + percentage_diff_index + 0.5,
                                                                 strand_bias_index + 0.5), radius_untranscribed,
                                                                color='yellowgreen', fill=True))
                                    # Second transcribed
                                    ax.add_patch(plt.Circle((mutation_type_index * len(
                                        percentage_strings) + percentage_diff_index + 0.5,
                                                                 strand_bias_index + 0.5), radius_transcribed,
                                                                color='royalblue', fill=True))


# June 2, 2021
def plot_strand_bias_figure_with_bar_plots(strand_bias,
                                     strandbias_figures_outputDir,
                                     numberofSimulations,
                                     signature,
                                     N,
                                     x_axis_tick_labels,
                                     y_axis_label,
                                     strand1_values,
                                     strand2_values,
                                     strand1_simulations_median_values,
                                     strand2_simulations_median_values,
                                     fdr_bh_adjusted_pvalues,
                                     strand1Name,
                                     strand2Name,
                                     color1,
                                     color2,
                                     width,
                                     axis_given=None):

    # Here we can take into difference between strand1_values and strand2_values while deciding on significance
    # the x locations for the groups
    ind = np.arange(N)
    if axis_given==None:
        fig, ax = plt.subplots(figsize=(16,10),dpi=100)
    else:
        ax=axis_given

    legend=None
    rects3=None
    rects4=None

    rects1 = ax.bar(ind, strand1_values, width=width, edgecolor='black', color=color1)
    rects2 = ax.bar(ind + width, strand2_values, width=width, edgecolor='black', color=color2)

    if ((strand1_simulations_median_values is not None) and strand1_simulations_median_values):
        rects3 = ax.bar(ind+ 2*width, strand1_simulations_median_values, width=width, edgecolor='black', color=color1, hatch = '///')
    if ((strand2_simulations_median_values is not None) and strand2_simulations_median_values):
        rects4 = ax.bar(ind + 3*width, strand2_simulations_median_values, width=width, edgecolor='black', color=color2, hatch = '///')

    # add some text for labels, title and axes ticks
    ax.tick_params(axis='x', labelsize=35)
    ax.tick_params(axis='y', labelsize=35)

    # Old way for y axis tick labels
    # locs = ax.get_yticks() # numpy.ndarray
    # ax.set_ylim(0, locs[-1] + 5000)

    # New way for y axis ticks and labels
    ymax = max(max(strand1_values), max(strand2_values), max(strand1_simulations_median_values), max(strand2_simulations_median_values))
    y = ymax / 1.025
    ytick_offest = float(y / 3)
    ylabs = [0, ytick_offest, ytick_offest * 2, ytick_offest * 3, ytick_offest * 4]
    ylabels = [0, ytick_offest, ytick_offest * 2, ytick_offest * 3, ytick_offest * 4]

    ylabels = ['{:,}'.format(int(x)) for x in ylabels]
    if len(ylabels[-1]) > 3:
        ylabels_temp = []
        if len(ylabels[-1]) > 7:
            for label in ylabels:
                if len(label) > 7:
                    ylabels_temp.append(label[0:-8] + "m")
                elif len(label) > 3:
                    ylabels_temp.append(label[0:-4] + "k")
                else:
                    ylabels_temp.append(label)
        else:
            for label in ylabels:
                if len(label) > 3:
                    ylabels_temp.append(label[0:-4] + "k")
                else:
                    ylabels_temp.append(label)
        ylabels = ylabels_temp

    ax.set_ylim([0, y])
    ax.set_yticks(ylabs)
    ax.set_yticklabels(ylabels, fontsize=35, fontweight='bold', fontname='Arial')

    # To make the bar width not too wide
    if len(ind) < 6:
        maxn = 6
        ax.set_xlim(-0.5, maxn - 0.5)

    # Set title
    ax.set_title('%s vs. %s' %(strand1Name,strand2Name), fontsize=40, fontweight='bold')

    # Set x tick labels
    if len(x_axis_tick_labels) > 6:
        ax.set_xticklabels(x_axis_tick_labels, fontsize=35, rotation=90)
    else:
        ax.set_xticklabels(x_axis_tick_labels, fontsize=35)

    # Set the ylabel
    if y_axis_label:
        ax.set_ylabel(y_axis_label, fontsize=35, fontweight='normal', labelpad=15)

    # Set the x axis tick locations
    if (numberofSimulations > 0):
        ax.set_xticks(ind + (3 * width) / 2)
        realStrand1Name = 'Real %s' % (strand1Name)
        realStrand2Name = 'Real %s' % (strand2Name)
        simulationsStrand1Name = 'Simulated %s' % (strand1Name)
        simulationsStrand2Name = 'Simulated %s' % (strand2Name)
        if ((rects1 is not None) and (rects2 is not None) and (rects3 is not None) and (rects4 is not None)):
            if ((len(rects1) > 0) and (len(rects2) > 0) and (len(rects3) > 0) and (len(rects4) > 0)):
                legend = ax.legend((rects1[0], rects2[0], rects3[0], rects4[0]),
                                   (realStrand1Name, realStrand2Name, simulationsStrand1Name, simulationsStrand2Name),prop={'size': 25}, ncol=1, loc='best')
    else:
        # Old way with no simulations
        ax.set_xticks(ind + width / 2)
        if ((rects1 is not None) and (rects2 is not None)):
            if ((len(rects1) > 0) and (len(rects2) > 0)):
                legend = ax.legend((rects1[0], rects2[0]), (strand1Name, strand2Name), prop={'size': 25}, ncol=1, loc='upper right')

    #To make the barplot background white
    ax.set_facecolor('white')
    #To makes spines black like a rectangle with black stroke
    ax.spines["bottom"].set_color('black')
    ax.spines["left"].set_color('black')
    ax.spines["top"].set_color('black')
    ax.spines["right"].set_color('black')

    if (legend is not None):
        frame = legend.get_frame()
        frame.set_facecolor('white')
        frame.set_edgecolor('black')

    #Add star above the bars for significant differences between the number of mutations on each strand starts
    # For each bar: Place a label
    if fdr_bh_adjusted_pvalues is not None:
        for fdr_bh_adjusted_pvalue, rect1, rect2 in zip(fdr_bh_adjusted_pvalues,rects1,rects2):
            # Get X and Y placement of label from rect.
            y_value = max(rect1.get_height(),rect2.get_height())
            x_value = rect1.get_x() + rect1.get_width()

            # Number of points between bar and label. Change to your liking.
            space = 3
            # Vertical alignment for positive values
            va = 'bottom'

            # If value of bar is negative: Place label below bar
            if y_value < 0:
                # Invert space to place label below
                space *= -1
                # Vertically align label at top
                va = 'top'

            # Use Y value as label and format number with one decimal place
            label = "{:.1f}".format(y_value)

            # Create annotation
            if ((fdr_bh_adjusted_pvalue is not None) and fdr_bh_adjusted_pvalue<=0.0001):
                ax.annotate(
                    '***',  # Use `label` as label
                    (x_value, y_value),  # Place label at end of the bar
                    xytext=(0, space),  # Vertically shift label by `space`
                    textcoords="offset points",  # Interpret `xytext` as offset in points
                    ha='center',  # Horizontally center label
                    va=va,
                    fontsize=25)  # Vertically align label differently for

            elif ((fdr_bh_adjusted_pvalue is not None) and fdr_bh_adjusted_pvalue<=0.001):
                ax.annotate(
                    '**',  # Use `label` as label
                    (x_value, y_value),  # Place label at end of the bar
                    xytext=(0, space),  # Vertically shift label by `space`
                    textcoords="offset points",  # Interpret `xytext` as offset in points
                    ha='center',  # Horizontally center label
                    va=va,
                    fontsize=25)  # Vertically align label differently for

            elif ((fdr_bh_adjusted_pvalue is not None) and fdr_bh_adjusted_pvalue<=SIGNIFICANCE_LEVEL):
                ax.annotate(
                    '*',  # Use `label` as label
                    (x_value, y_value),  # Place label at end of the bar
                    xytext=(0, space),  # Vertically shift label by `space`
                    textcoords="offset points",  # Interpret `xytext` as offset in points
                    ha='center',  # Horizontally center label
                    va=va,
                    fontsize=25) # Vertically align label differently for

    if axis_given==None:
        filename = '%s_%s_with_bars.png' %(signature,strand_bias)
        figFile = os.path.join(strandbias_figures_outputDir, filename)
        fig.savefig(figFile, dpi=100, bbox_inches="tight")

        plt.cla()
        plt.close(fig)


# June 2, 2021
def plot_bar_plot_in_given_axis(axis,
                                sbs_signature,
                                strand_bias,
                                strands_list,
                                signature_strand1_versus_strand2_df,
                                y_axis_label = None):
    box = axis.get_position()
    axis.set_position([box.x0, box.y0 + 0.125, box.width * 1, box.height * 1], which='both')

    mutation_types = six_mutation_types
    numberofSimulations=100
    width = 0.20

    if strand_bias==LAGGING_VERSUS_LEADING:
        strands=strands_list
        strand1="Lagging_real_count"
        strand2="Leading_real_count"
        strand1_sims="Lagging_mean_sims_count"
        strand2_sims="Leading_mean_sims_count"
        q_value_column_name = "lagging_versus_leading_q_value"
        color1='indianred'
        color2='goldenrod'
    elif strand_bias==TRANSCRIBED_VERSUS_UNTRANSCRIBED:
        strands=strands_list
        strand1="Transcribed_real_count"
        strand2="UnTranscribed_real_count"
        strand1_sims="Transcribed_mean_sims_count"
        strand2_sims="UnTranscribed_mean_sims_count"
        q_value_column_name="transcribed_versus_untranscribed_q_value"
        color1='royalblue'
        color2='yellowgreen'
    elif strand_bias==GENIC_VERSUS_INTERGENIC:
        strands=strands_list
        strand1="genic_real_count"
        strand2="intergenic_real_count"
        strand1_sims="genic_mean_sims_count"
        strand2_sims="intergenic_mean_sims_count"
        q_value_column_name="genic_versus_intergenic_q_value"
        color1='cyan'
        color2='gray'

    groupby_df = signature_strand1_versus_strand2_df.groupby(['signature'])
    group_df = groupby_df.get_group(sbs_signature)

    mutationtype_strand1_real_list=[]
    mutationtype_strand2_real_list=[]
    mutationtype_strand1_sims_mean_list=[]
    mutationtype_strand2_sims_mean_list=[]
    mutationtype_FDR_BH_adjusted_pvalues_list=[]

    for mutation_type in six_mutation_types:
        strand1_real_count=group_df[group_df['mutation_type']==mutation_type][strand1].values[0]
        strand2_real_count=group_df[group_df['mutation_type']==mutation_type][strand2].values[0]
        strand1_sims_count=group_df[group_df['mutation_type']==mutation_type][strand1_sims].values[0]
        strand2_sims_count=group_df[group_df['mutation_type']==mutation_type][strand2_sims].values[0]
        q_value=group_df[group_df['mutation_type']==mutation_type][q_value_column_name].values[0]
        mutationtype_strand1_real_list.append(strand1_real_count)
        mutationtype_strand2_real_list.append(strand2_real_count)
        mutationtype_strand1_sims_mean_list.append(strand1_sims_count)
        mutationtype_strand2_sims_mean_list.append(strand2_sims_count)
        mutationtype_FDR_BH_adjusted_pvalues_list.append(q_value)

    plot_strand_bias_figure_with_bar_plots(strand_bias,
                             None,
                             numberofSimulations,
                             sbs_signature,
                             len(mutation_types),
                             mutation_types,
                             y_axis_label,
                             mutationtype_strand1_real_list,
                             mutationtype_strand2_real_list,
                             mutationtype_strand1_sims_mean_list,
                             mutationtype_strand2_sims_mean_list,
                             mutationtype_FDR_BH_adjusted_pvalues_list,
                             strands[0],
                             strands[1],
                             color1,
                             color2,
                             width,
                             axis_given = axis)


# June 2, 2021
def plot_strand_bias_figure_with_stacked_bar_plots(strand_bias,
                                     strandbias_figures_outputDir,
                                     numberofSimulations,
                                     signature,
                                     N,
                                     x_axis_tick_labels,
                                     y_axis_label,
                                     strand1_values,
                                     strand2_values,
                                     strand1_simulations_median_values,
                                     strand2_simulations_median_values,
                                     fdr_bh_adjusted_pvalues,
                                     strand1Name,
                                     strand2Name,
                                     color1,
                                     color2,
                                     width,
                                     axis_given=None):

    # Replace np.nans with 0
    strand1_values = [0 if np.isnan(x) else x for x in strand1_values]
    strand2_values = [0 if np.isnan(x) else x for x in strand2_values]
    strand1_simulations_median_values = [0 if np.isnan(x) else x for x in strand1_simulations_median_values]
    strand2_simulations_median_values = [0 if np.isnan(x) else x for x in strand2_simulations_median_values]

    # Fill odds_ratio_list
    odds_real_list = []
    odds_sims_list = []
    for a, b in zip(strand1_values, strand2_values):
        odds_real = np.nan
        if b>0:
            odds_real = a/b
        odds_real_list.append(odds_real)

    for x, y in zip(strand1_simulations_median_values, strand2_simulations_median_values):
        odds_sims = np.nan
        if y>0:
            odds_sims = x/y
        odds_sims_list.append(odds_sims)

    odds_ratio_list = [odds_real/odds_sims if odds_sims>0 else np.nan for (odds_real, odds_sims) in zip(odds_real_list,odds_sims_list)]

    # Here we can take into difference between strand1_values and strand2_values while deciding on significance
    # the x locations for the groups
    ind = np.arange(N)
    if axis_given==None:
        fig, ax = plt.subplots(figsize=(16,10),dpi=100)
    else:
        ax=axis_given

    legend=None
    rects3=None
    rects4=None

    rects1 = ax.bar(ind, strand1_values, width=width, edgecolor='black', color=color1)
    rects2 = ax.bar(ind, strand2_values, width=width, edgecolor='black', color=color2, bottom=strand1_values)

    if ((strand1_simulations_median_values is not None) and strand1_simulations_median_values):
        rects3 = ax.bar(ind + width, strand1_simulations_median_values, width=width, edgecolor='black', color=color1, hatch = '///')
    if ((strand2_simulations_median_values is not None) and strand2_simulations_median_values):
        rects4 = ax.bar(ind + width, strand2_simulations_median_values, width=width, edgecolor='black', color=color2, hatch = '///', bottom=strand1_simulations_median_values)

    # Add some text for labels, title and axes ticks
    ax.tick_params(axis='x', labelsize=35)
    ax.tick_params(axis='y', labelsize=35)

    ax.set_ylim(0, 1.1)
    ax.set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=35)

    # To make the bar width not too wide
    if len(ind) < 6:
        maxn = 6
        ax.set_xlim(-0.5, maxn - 0.5)

    # Set title
    stacked_bar_title = 'Real vs. Simulated\nOdds Ratio of %s vs. %s' %(strand1Name, strand2Name)
    ax.set_title(stacked_bar_title, fontsize=40,  fontweight='bold')

    # Set x tick labels
    if len(x_axis_tick_labels) > 6:
        ax.set_xticklabels(x_axis_tick_labels, fontsize=35, rotation=90)
    else:
        ax.set_xticklabels(x_axis_tick_labels, fontsize=35)

    # Set the ylabel
    if y_axis_label:
        ax.set_ylabel(y_axis_label, fontsize=35, fontweight='normal', labelpad=15)

    # Set the x axis tick locations
    if (numberofSimulations > 0):
        ax.set_xticks(ind + (width/2))
        realStrand1Name = 'Real %s' % (strand1Name)
        realStrand2Name = 'Real %s' % (strand2Name)
        simulationsStrand1Name = 'Simulated %s' % (strand1Name)
        simulationsStrand2Name = 'Simulated %s' % (strand2Name)
        # # Let's not have a legend
        # if ((rects1 is not None) and (rects2 is not None) and (rects3 is not None) and (rects4 is not None)):
        #     if ((len(rects1) > 0) and (len(rects2) > 0) and (len(rects3) > 0) and (len(rects4) > 0)):
        #         legend = ax.legend((rects1[0], rects2[0], rects3[0], rects4[0]),
        #                            (realStrand1Name, realStrand2Name, simulationsStrand1Name, simulationsStrand2Name),prop={'size': 25}, ncol=1, loc='best')
    else:
        # Old way with no simulations
        ax.set_xticks(ind + width / 2)
        # if ((rects1 is not None) and (rects2 is not None)):
        #     if ((len(rects1) > 0) and (len(rects2) > 0)):
        #         legend = ax.legend((rects1[0], rects2[0]), (strand1Name, strand2Name), prop={'size': 25}, ncol=1, loc='upper right')

    #To make the barplot background white
    ax.set_facecolor('white')
    #To makes spines black like a rectangle with black stroke
    ax.spines["bottom"].set_color('black')
    ax.spines["left"].set_color('black')
    ax.spines["top"].set_color('black')
    ax.spines["right"].set_color('black')

    if (legend is not None):
        frame = legend.get_frame()
        frame.set_facecolor('white')
        frame.set_edgecolor('black')

    #Add star above the bars for significant differences between the number of mutations on each strand starts
    # For each bar: Place a label
    if odds_ratio_list is not None:
        for odds_ratio, fdr_bh_adjusted_pvalue, rect1, rect2 in zip(odds_ratio_list, fdr_bh_adjusted_pvalues, rects1, rects2):
            # Get X and Y placement of label from rect.
            # y_value = max(rect1.get_height(),rect2.get_height())
            y_value = rect1.get_height() + rect2.get_height()
            x_value = rect1.get_x() + rect1.get_width()

            # Number of points between bar and label. Change to your liking.
            space = 3
            # Vertical alignment for positive values
            va = 'bottom'

            # If value of bar is negative: Place label below bar
            if y_value < 0:
                # Invert space to place label below
                space *= -1
                # Vertically align label at top
                va = 'top'

            # Use Y value as label and format number with one decimal place
            label = "{:.1f}".format(y_value)

            # Create annotation
            if not np.isnan(odds_ratio):
                if ((fdr_bh_adjusted_pvalue is not None) and fdr_bh_adjusted_pvalue<=0.0001):
                    ax.annotate(
                        '%.2f ***' %(odds_ratio),  # Use `label` as label
                        (x_value, y_value),  # Place label at end of the bar
                        xytext=(0, space),  # Vertically shift label by `space`
                        textcoords="offset points",  # Interpret `xytext` as offset in points
                        ha='center',  # Horizontally center label
                        va=va,
                        fontsize=25)  # Vertically align label differently for

                elif ((fdr_bh_adjusted_pvalue is not None) and fdr_bh_adjusted_pvalue<=0.001):
                    ax.annotate(
                        '%.2f **' %(odds_ratio),  # Use `label` as label
                        (x_value, y_value),  # Place label at end of the bar
                        xytext=(0, space),  # Vertically shift label by `space`
                        textcoords="offset points",  # Interpret `xytext` as offset in points
                        ha='center',  # Horizontally center label
                        va=va,
                        fontsize=25)  # Vertically align label differently for

                elif ((fdr_bh_adjusted_pvalue is not None) and fdr_bh_adjusted_pvalue<=SIGNIFICANCE_LEVEL):
                    ax.annotate(
                        '%.2f *' %(odds_ratio),  # Use `label` as label
                        (x_value, y_value),  # Place label at end of the bar
                        xytext=(0, space),  # Vertically shift label by `space`
                        textcoords="offset points",  # Interpret `xytext` as offset in points
                        ha='center',  # Horizontally center label
                        va=va,
                        fontsize=25) # Vertically align label differently for
                else:
                    ax.annotate(
                        '%.2f' %(odds_ratio),  # Use `label` as label
                        (x_value, y_value),  # Place label at end of the bar
                        xytext=(0, space),  # Vertically shift label by `space`
                        textcoords="offset points",  # Interpret `xytext` as offset in points
                        ha='center',  # Horizontally center label
                        va=va,
                        fontsize=25) # Vertically align label differently for

    if axis_given==None:
        filename = '%s_%s_with_bars.png' %(signature,strand_bias)
        figFile = os.path.join(strandbias_figures_outputDir, filename)
        fig.savefig(figFile, dpi=100, bbox_inches="tight")

        plt.cla()
        plt.close(fig)


# June 2, 2021
def plot_stacked_bar_plot_in_given_axis(axis,
                                        sbs_signature,
                                        strand_bias,
                                        strands_list,
                                        signature_strand1_versus_strand2_df,
                                        y_axis_label = None):
    box = axis.get_position()
    axis.set_position([box.x0, box.y0+0.125, box.width * 1, box.height * 1], which='both')

    mutation_types = six_mutation_types
    numberofSimulations=100
    width = 0.20

    if strand_bias==LAGGING_VERSUS_LEADING:
        strands=strands_list
        strand1="Lagging_real_count"
        strand2="Leading_real_count"
        strand1_sims="Lagging_mean_sims_count"
        strand2_sims="Leading_mean_sims_count"
        q_value_column_name = "lagging_versus_leading_q_value"
        color1='indianred'
        color2='goldenrod'
    elif strand_bias==TRANSCRIBED_VERSUS_UNTRANSCRIBED:
        strands=strands_list
        strand1="Transcribed_real_count"
        strand2="UnTranscribed_real_count"
        strand1_sims="Transcribed_mean_sims_count"
        strand2_sims="UnTranscribed_mean_sims_count"
        q_value_column_name="transcribed_versus_untranscribed_q_value"
        color1='royalblue'
        color2='yellowgreen'
    elif strand_bias==GENIC_VERSUS_INTERGENIC:
        strands=strands_list
        strand1="genic_real_count"
        strand2="intergenic_real_count"
        strand1_sims="genic_mean_sims_count"
        strand2_sims="intergenic_mean_sims_count"
        q_value_column_name="genic_versus_intergenic_q_value"
        color1='cyan'
        color2='gray'

    groupby_df = signature_strand1_versus_strand2_df.groupby(['signature'])
    group_df = groupby_df.get_group(sbs_signature)

    mutationtype_strand1_real_list=[]
    mutationtype_strand2_real_list=[]
    mutationtype_strand1_sims_mean_list=[]
    mutationtype_strand2_sims_mean_list=[]
    mutationtype_FDR_BH_adjusted_pvalues_list=[]

    for mutation_type in six_mutation_types:
        strand1_real_count=group_df[group_df['mutation_type']==mutation_type][strand1].values[0]
        strand2_real_count=group_df[group_df['mutation_type']==mutation_type][strand2].values[0]
        strand1_sims_count=group_df[group_df['mutation_type']==mutation_type][strand1_sims].values[0]
        strand2_sims_count=group_df[group_df['mutation_type']==mutation_type][strand2_sims].values[0]
        q_value=group_df[group_df['mutation_type']==mutation_type][q_value_column_name].values[0]
        mutationtype_FDR_BH_adjusted_pvalues_list.append(q_value)

        if (strand1_real_count >= NUMBER_OF_REQUIRED_MUTATIONS_FOR_STRAND_BIAS_BAR_PLOT) or (strand2_real_count >= NUMBER_OF_REQUIRED_MUTATIONS_FOR_STRAND_BIAS_BAR_PLOT):
            mutationtype_strand1_real_list.append(strand1_real_count/(strand1_real_count+strand2_real_count))
            mutationtype_strand2_real_list.append(strand2_real_count/(strand1_real_count+strand2_real_count))
        else:
            mutationtype_strand1_real_list.append(np.nan)
            mutationtype_strand2_real_list.append(np.nan)

        if (strand1_sims_count >= NUMBER_OF_REQUIRED_MUTATIONS_FOR_STRAND_BIAS_BAR_PLOT) or (strand2_sims_count >= NUMBER_OF_REQUIRED_MUTATIONS_FOR_STRAND_BIAS_BAR_PLOT):
            mutationtype_strand1_sims_mean_list.append(strand1_sims_count/(strand1_sims_count+strand2_sims_count))
            mutationtype_strand2_sims_mean_list.append(strand2_sims_count/(strand1_sims_count+strand2_sims_count))
        else:
            mutationtype_strand1_sims_mean_list.append(np.nan)
            mutationtype_strand2_sims_mean_list.append(np.nan)

    plot_strand_bias_figure_with_stacked_bar_plots(strand_bias,
                             None,
                             numberofSimulations,
                             sbs_signature,
                             len(mutation_types),
                             mutation_types,
                             y_axis_label,
                             mutationtype_strand1_real_list,
                             mutationtype_strand2_real_list,
                             mutationtype_strand1_sims_mean_list,
                             mutationtype_strand2_sims_mean_list,
                             mutationtype_FDR_BH_adjusted_pvalues_list,
                             strands[0],
                             strands[1],
                             color1,
                             color2,
                             width,
                             axis_given=axis)



def plot_circle_bar_plots_together(outputDir,
                                jobname,
                                sbs_signature,
                                six_mutation_types,
                                signature2mutation_type2strand2percentagedict,
                               signature_genic_versus_intergenic_df,
                               signature_transcribed_versus_untranscribed_df,
                               signature_lagging_versus_leading_df,
                               genic_vs_intergenic_strands,
                               transcription_strands,
                               replication_strands):

    x_ticklabels_list = percentage_strings * 6
    # fig = plt.figure(figsize=(5 + 1.5 * len(x_ticklabels_list), 30 + 1.5))
    fig = plt.figure(figsize=(5 + 1.5 * len(x_ticklabels_list), 30 + 1.5))
    plt.rc('axes', edgecolor='lightgray')

    width = 6
    height = 6
    width_ratios = [1] * width
    height_ratios = [1] * height
    gs = gridspec.GridSpec(height, width, height_ratios = height_ratios, width_ratios = width_ratios)
    fig.subplots_adjust(hspace=0, wspace=3)

    cirle_plot_axis = plt.subplot(gs[0:2, :])

    genic_vs_intergenic_bar_plot_axis = plt.subplot(gs[2:4, 0:2])
    transcribed_vs_untranscribed_bar_plot_axis = plt.subplot(gs[2:4, 2:4])
    lagging_vs_leading_bar_plot_axis = plt.subplot(gs[2:4, 4:6])

    genic_vs_intergenic_stacked_bar_plot_axis = plt.subplot(gs[4:, 0:2])
    transcribed_vs_untranscribed_stacked_bar_plot_axis = plt.subplot(gs[4:, 2:4])
    lagging_vs_leading_stacked_bar_plot_axis = plt.subplot(gs[4:, 4:6])

    # Circle plot with legends
    plot_circle_plot_in_given_axis(cirle_plot_axis,
                               percentage_strings,
                               sbs_signature,
                               six_mutation_types,
                               x_ticklabels_list,
                               signature2mutation_type2strand2percentagedict)

    # 3 Bar plots side by side
    plot_bar_plot_in_given_axis(genic_vs_intergenic_bar_plot_axis,
                                sbs_signature,
                                GENIC_VERSUS_INTERGENIC,
                                genic_vs_intergenic_strands,
                                signature_genic_versus_intergenic_df,
                                y_axis_label = 'Number of Single Base Substitutions')

    plot_bar_plot_in_given_axis(transcribed_vs_untranscribed_bar_plot_axis,
                                sbs_signature,
                                TRANSCRIBED_VERSUS_UNTRANSCRIBED,
                                transcription_strands,
                                signature_transcribed_versus_untranscribed_df)

    plot_bar_plot_in_given_axis(lagging_vs_leading_bar_plot_axis,
                                sbs_signature,
                                LAGGING_VERSUS_LEADING,
                                replication_strands,
                                signature_lagging_versus_leading_df)

    # 3 Stacked Bar plots side by side
    plot_stacked_bar_plot_in_given_axis(genic_vs_intergenic_stacked_bar_plot_axis,
                                        sbs_signature,
                                        GENIC_VERSUS_INTERGENIC,
                                        genic_vs_intergenic_strands,
                                        signature_genic_versus_intergenic_df,
                                        y_axis_label = 'Ratio of mutations on each strand')

    plot_stacked_bar_plot_in_given_axis(transcribed_vs_untranscribed_stacked_bar_plot_axis,
                                        sbs_signature,
                                        TRANSCRIBED_VERSUS_UNTRANSCRIBED,
                                        transcription_strands,
                                        signature_transcribed_versus_untranscribed_df)

    plot_stacked_bar_plot_in_given_axis(lagging_vs_leading_stacked_bar_plot_axis,
                                        sbs_signature,
                                        LAGGING_VERSUS_LEADING,
                                        replication_strands,
                                        signature_lagging_versus_leading_df)

    # filename = '%s_circle_bar_plot_together_%s.png' % (sbs_signature, str(significance_level).replace('.', '_'))
    filename = '%s_circle_bar_plots.png' % (sbs_signature)
    figurepath = os.path.join(outputDir, jobname, FIGURE, STRANDBIAS, CIRCLE_BAR_PLOTS, filename)
    fig.savefig(figurepath, dpi=100, bbox_inches="tight")

    plt.cla()
    plt.close(fig)



##################################################################
# July 5, 2020 starts
#Key can be signature or sample
def plotBarPlotsUsingDataframes(outputDir,
                 jobname,
                 numberofSimulations,
                 signature_cutoff_numberofmutations_averageprobability_df,
                 isKeySample,
                 existingMutationTypesList,
                 signature_strand1_versus_strand2_df,
                 width,
                 strand1_versus_strand2,
                 strands,
                 color1,
                 color2,
                 title,
                 figureName,
                 plot_mode):

    # signature_strand1_versus_strand2_df column names here
    # ['cancer_type', 'signature', 'mutation_type',
    #  'Transcribed_real_count', 'UnTranscribed_real_count', 'Transcribed_mean_sims_count', 'UnTranscribed_mean_sims_count',
    #  'transcribed_versus_untranscribed_p_value', 'transcribed_versus_untranscribed_q_value',
    #  'Transcribed_real_count.1', 'Transcribed_mean_sims_count.1', 'Transcribed_min_sims_count', 'Transcribed_max_sims_count', 'Transcribed_sims_count_list',
    #  'UnTranscribed_real_count.1', 'UnTranscribed_mean_sims_count.1', 'UnTranscribed_min_sims_count', 'UnTranscribed_max_sims_count', 'UnTranscribed_sims_count_list']

    signatures = signature_strand1_versus_strand2_df['signature'].unique()

    x_axis_labels = existingMutationTypesList
    N = len(x_axis_labels)

    for signature in signatures:

        numberofMutations = int(signature_cutoff_numberofmutations_averageprobability_df[signature_cutoff_numberofmutations_averageprobability_df['signature'] == signature]['number_of_mutations'].values[0])

        mutationtype_strand1_real_list=[]
        mutationtype_strand2_real_list=[]
        mutationtype_strand1_sims_mean_list=[]
        mutationtype_strand2_sims_mean_list=[]
        mutationtype_FDR_BH_adjusted_pvalues_list=[]

        for mutation_type in existingMutationTypesList:
            if (strand1_versus_strand2==TRANSCRIBED_VERSUS_UNTRANSCRIBED):
                strand1_real_count_column_name=TRANSCRIBED_REAL_COUNT
                strand1_sims_mean_count_Column_name=TRANSCRIBED_SIMULATIONS_MEAN_COUNT
                strand2_real_count_column_name=UNTRANSCRIBED_REAL_COUNT
                strand2_sims_mean_count_Column_name=UNTRANSCRIBED_SIMULATIONS_MEAN_COUNT
                q_value_column_name = TRANSCRIBED_VERSUS_UNTRANSCRIBED_Q_VALUE
            elif (strand1_versus_strand2 == GENIC_VERSUS_INTERGENIC):
                strand1_real_count_column_name=GENIC_REAL_COUNT
                strand1_sims_mean_count_Column_name=GENIC_SIMULATIONS_MEAN_COUNT
                strand2_real_count_column_name=INTERGENIC_REAL_COUNT
                strand2_sims_mean_count_Column_name=INTERGENIC_SIMULATIONS_MEAN_COUNT
                q_value_column_name = GENIC_VERSUS_INTERGENIC_Q_VALUE
            elif (strand1_versus_strand2 == LAGGING_VERSUS_LEADING):
                strand1_real_count_column_name=LAGGING_REAL_COUNT
                strand1_sims_mean_count_Column_name=LAGGING_SIMULATIONS_MEAN_COUNT
                strand2_real_count_column_name=LEADING_REAL_COUNT
                strand2_sims_mean_count_Column_name=LEADING_SIMULATIONS_MEAN_COUNT
                q_value_column_name = LAGGING_VERSUS_LEADING_Q_VALUE

            strand1_real_count=0
            strand1_sims_mean_count=0
            strand2_real_count=0
            strand2_sims_mean_count=0
            q_value=None

            if (signature_strand1_versus_strand2_df[(signature_strand1_versus_strand2_df['signature']==signature) & (signature_strand1_versus_strand2_df['mutation_type']==mutation_type)][strand1_real_count_column_name].values.size>0):
                strand1_real_count=signature_strand1_versus_strand2_df[(signature_strand1_versus_strand2_df['signature']==signature) & (signature_strand1_versus_strand2_df['mutation_type']==mutation_type)][strand1_real_count_column_name].values[0]

            if (signature_strand1_versus_strand2_df[(signature_strand1_versus_strand2_df['signature'] == signature) & (signature_strand1_versus_strand2_df['mutation_type'] == mutation_type)][strand1_sims_mean_count_Column_name].values.size>0):
                strand1_sims_mean_count = signature_strand1_versus_strand2_df[(signature_strand1_versus_strand2_df['signature'] == signature) & (signature_strand1_versus_strand2_df['mutation_type'] == mutation_type)][strand1_sims_mean_count_Column_name].values[0]

            if (signature_strand1_versus_strand2_df[(signature_strand1_versus_strand2_df['signature']==signature) & (signature_strand1_versus_strand2_df['mutation_type']==mutation_type)][strand2_real_count_column_name].values.size>0):
                strand2_real_count=signature_strand1_versus_strand2_df[(signature_strand1_versus_strand2_df['signature']==signature) & (signature_strand1_versus_strand2_df['mutation_type']==mutation_type)][strand2_real_count_column_name].values[0]

            if  (signature_strand1_versus_strand2_df[(signature_strand1_versus_strand2_df['signature'] == signature) & (signature_strand1_versus_strand2_df['mutation_type'] == mutation_type)][strand2_sims_mean_count_Column_name].values.size>0):
                strand2_sims_mean_count = signature_strand1_versus_strand2_df[(signature_strand1_versus_strand2_df['signature'] == signature) & (signature_strand1_versus_strand2_df['mutation_type'] == mutation_type)][strand2_sims_mean_count_Column_name].values[0]

            if (signature_strand1_versus_strand2_df[(signature_strand1_versus_strand2_df['signature'] == signature) & (signature_strand1_versus_strand2_df['mutation_type'] == mutation_type)][q_value_column_name].values.size>0):
                q_value = signature_strand1_versus_strand2_df[(signature_strand1_versus_strand2_df['signature'] == signature) & (signature_strand1_versus_strand2_df['mutation_type'] == mutation_type)][q_value_column_name].values[0]

            mutationtype_strand1_real_list.append(strand1_real_count)
            mutationtype_strand1_sims_mean_list.append(strand1_sims_mean_count)
            mutationtype_strand2_real_list.append(strand2_real_count)
            mutationtype_strand2_sims_mean_list.append(strand2_sims_mean_count)
            mutationtype_FDR_BH_adjusted_pvalues_list.append(q_value)

        plotStrandBiasFigureWithBarPlots(outputDir,
                                         jobname,
                                         numberofSimulations,
                                         signature,
                                         isKeySample,
                                         numberofMutations,
                                         N,
                                         x_axis_labels,
                                         mutationtype_strand1_real_list,
                                         mutationtype_strand2_real_list,
                                         mutationtype_strand1_sims_mean_list,
                                         mutationtype_strand2_sims_mean_list,
                                         mutationtype_FDR_BH_adjusted_pvalues_list,
                                         strands[0],
                                         strands[1],
                                         title,
                                         color1,
                                         color2,
                                         figureName,
                                         width,
                                         plot_mode)

# July 5, 2020 ends
##################################################################



###################################################################
# April 20, 2020
# July 4, 2020 starts
# Using dataframes
def transcriptionReplicationStrandBiasFiguresUsingDataframes(outputDir, jobname, numberofSimulations, mutation_types_contexts, strand_bias_list, is_discreet, plot_mode):

    # Initialize these dataframes as empty dataframe
    # We will read these dataframes if there is the corresponding data
    subsSignature_cutoff_numberofmutations_averageprobability_df = pd.DataFrame()
    dinucsSignature_cutoff_numberofmutations_averageprobability_df = pd.DataFrame()
    indelsSignature_cutoff_numberofmutations_averageprobability_df = pd.DataFrame()

    sbs_df = pd.DataFrame()
    dbs_df = pd.DataFrame()
    id_df = pd.DataFrame()

    subsSignatures = np.array([])
    dinucsSignatures = np.array([])
    indelsSignatures = np.array([])

    #######################################################################
    os.makedirs(os.path.join(outputDir, jobname, FIGURE, STRANDBIAS,SCATTER_PLOTS), exist_ok=True)
    os.makedirs(os.path.join(outputDir, jobname, FIGURE, STRANDBIAS,BAR_PLOTS), exist_ok=True)
    os.makedirs(os.path.join(outputDir, jobname, FIGURE, STRANDBIAS,CIRCLE_PLOTS), exist_ok=True)
    os.makedirs(os.path.join(outputDir, jobname, FIGURE, STRANDBIAS,CIRCLE_BAR_PLOTS), exist_ok=True)
    os.makedirs(os.path.join(outputDir, jobname, FIGURE, STRANDBIAS,TABLES), exist_ok=True)
    os.makedirs(os.path.join(outputDir, jobname, FIGURE, STRANDBIAS,EXCEL_FILES), exist_ok=True)
    strandbias_figures_outputDir= os.path.join(outputDir, jobname, FIGURE, STRANDBIAS)
    strandbias_figures_tables_outputDir= os.path.join(outputDir, jobname, FIGURE, STRANDBIAS, TABLES)
    strandbias_figures_excel_files_outputDir= os.path.join(outputDir, jobname, FIGURE, STRANDBIAS, EXCEL_FILES)
    #######################################################################

    ##########################################################################################
    #########################  Read dictionaries related with ################################
    #########################  signatures and samples starts  ################################
    ##########################################################################################
    for mutation_type_context in mutation_types_contexts:
        if (mutation_type_context in SBS_CONTEXTS):
            subsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(os.path.join(outputDir, jobname, DATA,Table_SubsSignature_Cutoff_NumberofMutations_AverageProbability_Filename),sep='\t', header=0,dtype={'cutoff': np.float32,'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})
            subsSignatures = subsSignature_cutoff_numberofmutations_averageprobability_df['signature'].unique()
    if (DBS in mutation_types_contexts):
        dinucsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(os.path.join(outputDir, jobname, DATA,Table_DinucsSignature_Cutoff_NumberofMutations_AverageProbability_Filename),sep='\t', header=0,dtype={'cutoff': np.float32,'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})
        dinucsSignatures = dinucsSignature_cutoff_numberofmutations_averageprobability_df['signature'].unique()
    if (ID in mutation_types_contexts):
        indelsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(os.path.join(outputDir, jobname, DATA,Table_IndelsSignature_Cutoff_NumberofMutations_AverageProbability_Filename),sep='\t', header=0,dtype={'cutoff': np.float32,'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})
        indelsSignatures = indelsSignature_cutoff_numberofmutations_averageprobability_df['signature'].unique()
    ##########################################################################################
    #########################  Read dictionaries related with ################################
    #########################  signatures and samples ends  ##################################
    ##########################################################################################

    if is_discreet:
        sbs_df = subsSignature_cutoff_numberofmutations_averageprobability_df
        dbs_df = dinucsSignature_cutoff_numberofmutations_averageprobability_df
        id_df = indelsSignature_cutoff_numberofmutations_averageprobability_df
    else:
        if os.path.exists(os.path.join(outputDir, jobname, DATA, Table_SubsSignature_NumberofMutations_AverageProbability_Filename)):
            sbs_df = pd.read_csv(os.path.join(outputDir, jobname, DATA, Table_SubsSignature_NumberofMutations_AverageProbability_Filename), sep='\t', header=0, dtype={'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})
            subsSignatures = sbs_df['signature'].unique()
        if os.path.exists(os.path.join(outputDir, jobname, DATA, Table_DinucsSignature_NumberofMutations_AverageProbability_Filename)):
            dbs_df = pd.read_csv(os.path.join(outputDir, jobname, DATA, Table_DinucsSignature_NumberofMutations_AverageProbability_Filename), sep='\t', header=0, dtype={'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})
            dinucsSignatures = dbs_df['signature'].unique()
        if os.path.exists(os.path.join(outputDir, jobname, DATA, Table_IndelsSignature_NumberofMutations_AverageProbability_Filename)):
            id_df = pd.read_csv(os.path.join(outputDir, jobname, DATA, Table_IndelsSignature_NumberofMutations_AverageProbability_Filename), sep='\t', header=0, dtype={'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})
            indelsSignatures = id_df['signature'].unique()


    #######################################################################
    #Step1 Read p_value
    if LAGGING_VERSUS_LEADING in strand_bias_list:
        #Replication Strand Bias
        signature_mutation_type_lagging_versus_leading_table_file_name = 'Signature_Mutation_Type_%s_Strand_Table.txt' % (LAGGING_VERSUS_LEADING)
        signature_mutation_type_lagging_versus_leading_table_filepath = os.path.join(outputDir, jobname, DATA, REPLICATIONSTRANDBIAS,signature_mutation_type_lagging_versus_leading_table_file_name)
        signature_lagging_versus_leading_df = pd.read_csv(signature_mutation_type_lagging_versus_leading_table_filepath, header=0, sep='\t')

        type_lagging_versus_leading_table_file_name = 'Type_%s_Strand_Table.txt' % (LAGGING_VERSUS_LEADING)
        type_lagging_versus_leading_table_filepath = os.path.join(outputDir, jobname, DATA, REPLICATIONSTRANDBIAS,type_lagging_versus_leading_table_file_name)
        type_lagging_versus_leading_df = pd.read_csv(type_lagging_versus_leading_table_filepath, header=0, sep='\t')

    if TRANSCRIBED_VERSUS_UNTRANSCRIBED in strand_bias_list:
        #Transcription Strand Bias
        signature_mutation_type_transcribed_versus_untranscribed_table_file_name = 'Signature_Mutation_Type_%s_Strand_Table.txt' % (TRANSCRIBED_VERSUS_UNTRANSCRIBED)
        signature_mutation_type_transcribed_versus_untranscribed_table_filepath = os.path.join(outputDir, jobname, DATA, TRANSCRIPTIONSTRANDBIAS, signature_mutation_type_transcribed_versus_untranscribed_table_file_name)
        signature_transcribed_versus_untranscribed_df = pd.read_csv(signature_mutation_type_transcribed_versus_untranscribed_table_filepath, header=0, sep='\t')

        type_transcribed_versus_untranscribed_table_file_name = 'Type_%s_Strand_Table.txt' % (TRANSCRIBED_VERSUS_UNTRANSCRIBED)
        type_transcribed_versus_untranscribed_table_filepath = os.path.join(outputDir, jobname, DATA, TRANSCRIPTIONSTRANDBIAS, type_transcribed_versus_untranscribed_table_file_name)
        type_transcribed_versus_untranscribed_df = pd.read_csv(type_transcribed_versus_untranscribed_table_filepath, header=0, sep='\t')

    if GENIC_VERSUS_INTERGENIC in strand_bias_list:
        #Transcription Strand Bias
        signature_mutation_type_genic_versus_intergenic_table_file_name = 'Signature_Mutation_Type_%s_Strand_Table.txt' % (GENIC_VERSUS_INTERGENIC)
        signature_mutation_type_genic_versus_intergenic_table_filepath = os.path.join(outputDir, jobname, DATA, TRANSCRIPTIONSTRANDBIAS, signature_mutation_type_genic_versus_intergenic_table_file_name)
        signature_genic_versus_intergenic_df = pd.read_csv(signature_mutation_type_genic_versus_intergenic_table_filepath, header=0, sep='\t')

        type_genic_versus_intergenic_table_file_name = 'Type_%s_Strand_Table.txt' % (GENIC_VERSUS_INTERGENIC)
        type_genic_versus_intergenic_table_filepath = os.path.join(outputDir, jobname, DATA, TRANSCRIPTIONSTRANDBIAS, type_genic_versus_intergenic_table_file_name)
        type_genic_versus_intergenic_df = pd.read_csv(type_genic_versus_intergenic_table_filepath, header=0, sep='\t')
    #######################################################################

    #######################################################################
    #Step2 Compute q_value
    p_values_list=[]
    element_names=[]

    #Fill p_values_list
    if LAGGING_VERSUS_LEADING in strand_bias_list:
        for index, row in signature_lagging_versus_leading_df.iterrows():
            element_name = (row[CANCER_TYPE], row[SIGNATURE], row[MUTATION_TYPE], LAGGING_VERSUS_LEADING)
            element_names.append(element_name)
            p_values_list.append(row[LAGGING_VERSUS_LEADING_P_VALUE])

        for index, row in type_lagging_versus_leading_df.iterrows():
            element_name=(row[CANCER_TYPE], None, row[TYPE], LAGGING_VERSUS_LEADING)
            element_names.append(element_name)
            p_values_list.append(row[LAGGING_VERSUS_LEADING_P_VALUE])

    if TRANSCRIBED_VERSUS_UNTRANSCRIBED in strand_bias_list:
        for index, row in signature_transcribed_versus_untranscribed_df.iterrows():
            element_name=(row[CANCER_TYPE], row[SIGNATURE], row[MUTATION_TYPE], TRANSCRIBED_VERSUS_UNTRANSCRIBED)
            element_names.append(element_name)
            p_values_list.append(row[TRANSCRIBED_VERSUS_UNTRANSCRIBED_P_VALUE])

        for index, row in type_transcribed_versus_untranscribed_df.iterrows():
            element_name=(row[CANCER_TYPE], None, row[TYPE], TRANSCRIBED_VERSUS_UNTRANSCRIBED)
            element_names.append(element_name)
            p_values_list.append(row[TRANSCRIBED_VERSUS_UNTRANSCRIBED_P_VALUE])

    if GENIC_VERSUS_INTERGENIC in strand_bias_list:
        for index, row in signature_genic_versus_intergenic_df.iterrows():
            element_name = (row[CANCER_TYPE], row[SIGNATURE], row[MUTATION_TYPE], GENIC_VERSUS_INTERGENIC)
            element_names.append(element_name)
            p_values_list.append(row[GENIC_VERSUS_INTERGENIC_P_VALUE])

        for index, row in type_genic_versus_intergenic_df.iterrows():
            element_name=(row[CANCER_TYPE], None, row[TYPE], GENIC_VERSUS_INTERGENIC)
            element_names.append(element_name)
            p_values_list.append(row[GENIC_VERSUS_INTERGENIC_P_VALUE])

    # print('len(p_values_list): %d' %(len(p_values_list)))
    #######################################################################

    #######################################################################
    if ((p_values_list is not None) and p_values_list):
        rejected, all_FDR_BH_adjusted_p_values, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(p_values_list, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

        #Add None q_values
        if LAGGING_VERSUS_LEADING in strand_bias_list:
            signature_lagging_versus_leading_df[LAGGING_VERSUS_LEADING_Q_VALUE] = np.nan
            type_lagging_versus_leading_df[LAGGING_VERSUS_LEADING_Q_VALUE] = np.nan

        if TRANSCRIBED_VERSUS_UNTRANSCRIBED in strand_bias_list:
            signature_transcribed_versus_untranscribed_df[TRANSCRIBED_VERSUS_UNTRANSCRIBED_Q_VALUE] = np.nan
            type_transcribed_versus_untranscribed_df[TRANSCRIBED_VERSUS_UNTRANSCRIBED_Q_VALUE]= np.nan

        if GENIC_VERSUS_INTERGENIC in strand_bias_list:
            signature_genic_versus_intergenic_df[GENIC_VERSUS_INTERGENIC_Q_VALUE]= np.nan
            type_genic_versus_intergenic_df[GENIC_VERSUS_INTERGENIC_Q_VALUE]= np.nan

        #Update q_value
        for element_index, element_name in enumerate(element_names,0):
            (cancer_type, signature, mutation_type, versus_type)=element_name
            q_value=all_FDR_BH_adjusted_p_values[element_index]

            # print('%s %s %s %s -----  q_value: %f' %(cancer_type,signature,mutation_type,versus_type,q_value))

            if (signature is not None) and (versus_type==TRANSCRIBED_VERSUS_UNTRANSCRIBED):
                signature_transcribed_versus_untranscribed_df.loc[(signature_transcribed_versus_untranscribed_df[CANCER_TYPE]==cancer_type) &
                                                   (signature_transcribed_versus_untranscribed_df[SIGNATURE]==signature) &
                                                   (signature_transcribed_versus_untranscribed_df[MUTATION_TYPE]==mutation_type),TRANSCRIBED_VERSUS_UNTRANSCRIBED_Q_VALUE]=q_value

            elif (signature is not None) and (versus_type==GENIC_VERSUS_INTERGENIC):
                signature_genic_versus_intergenic_df.loc[(signature_genic_versus_intergenic_df[CANCER_TYPE]==cancer_type) &
                                                   (signature_genic_versus_intergenic_df[SIGNATURE]==signature) &
                                                   (signature_genic_versus_intergenic_df[MUTATION_TYPE]==mutation_type),GENIC_VERSUS_INTERGENIC_Q_VALUE]=q_value

            elif (signature is not None) and (versus_type==LAGGING_VERSUS_LEADING):
                signature_lagging_versus_leading_df.loc[(signature_lagging_versus_leading_df[CANCER_TYPE]==cancer_type) &
                                                   (signature_lagging_versus_leading_df[SIGNATURE]==signature) &
                                                   (signature_lagging_versus_leading_df[MUTATION_TYPE]==mutation_type),LAGGING_VERSUS_LEADING_Q_VALUE]=q_value


            elif (signature is None) and (versus_type == TRANSCRIBED_VERSUS_UNTRANSCRIBED):
                type_transcribed_versus_untranscribed_df.loc[(type_transcribed_versus_untranscribed_df[CANCER_TYPE] == cancer_type) & (type_transcribed_versus_untranscribed_df[TYPE] == mutation_type),TRANSCRIBED_VERSUS_UNTRANSCRIBED_Q_VALUE] = q_value
            elif (signature is None) and (versus_type == GENIC_VERSUS_INTERGENIC):
                type_genic_versus_intergenic_df.loc[(type_genic_versus_intergenic_df[CANCER_TYPE] == cancer_type) & (type_genic_versus_intergenic_df[TYPE] == mutation_type),GENIC_VERSUS_INTERGENIC_Q_VALUE] = q_value
            elif (signature is None) and (versus_type == LAGGING_VERSUS_LEADING):
                type_lagging_versus_leading_df.loc[(type_lagging_versus_leading_df[CANCER_TYPE] == cancer_type) & (type_lagging_versus_leading_df[TYPE] == mutation_type),LAGGING_VERSUS_LEADING_Q_VALUE] = q_value


        # Reorder columns
        # Write dataframes
        if LAGGING_VERSUS_LEADING in strand_bias_list:
            signature_lagging_versus_leading_df = signature_lagging_versus_leading_df[
                ['cancer_type', 'signature', 'mutation_type',
                 'Lagging_real_count', 'Leading_real_count', 'Lagging_mean_sims_count', 'Leading_mean_sims_count',
                 'lagging_versus_leading_p_value', 'lagging_versus_leading_q_value',
                 'Lagging_real_count.1', 'Lagging_mean_sims_count.1', 'Lagging_min_sims_count',
                 'Lagging_max_sims_count', 'Lagging_sims_count_list',
                 'Leading_real_count.1', 'Leading_mean_sims_count.1', 'Leading_min_sims_count',
                 'Leading_max_sims_count', 'Leading_sims_count_list']]

            type_lagging_versus_leading_df=type_lagging_versus_leading_df[['cancer_type', 'type',
                'Lagging_real_count', 'Leading_real_count', 'Lagging_mean_sims_count', 'Leading_mean_sims_count', 'lagging_versus_leading_p_value', 'lagging_versus_leading_q_value',
                'Lagging_real_count.1', 'Lagging_mean_sims_count.1', 'Lagging_min_sims_count', 'Lagging_max_sims_count', 'Lagging_sims_count_list',
                'Leading_real_count.1', 'Leading_mean_sims_count.1', 'Leading_min_sims_count', 'Leading_max_sims_count', 'Leading_sims_count_list' ]]

            signature_filename = 'Signature_Mutation_Type_%s_Q_Value_Table.txt' % (LAGGING_VERSUS_LEADING)
            signature_filepath = os.path.join(strandbias_figures_tables_outputDir, signature_filename)
            signature_lagging_versus_leading_df.to_csv(signature_filepath, sep='\t', header=True, index=False)

            type_filename = 'Type_%s_Q_Value_Table.txt' % (LAGGING_VERSUS_LEADING)
            type_filepath = os.path.join(strandbias_figures_tables_outputDir, type_filename)
            type_lagging_versus_leading_df.to_csv(type_filepath, sep='\t', header=True, index=False)

        if TRANSCRIBED_VERSUS_UNTRANSCRIBED in strand_bias_list:
            signature_transcribed_versus_untranscribed_df=signature_transcribed_versus_untranscribed_df[['cancer_type', 'signature', 'mutation_type',
                'Transcribed_real_count', 'UnTranscribed_real_count', 'NonTranscribed_real_count',
                'Transcribed_mean_sims_count', 'UnTranscribed_mean_sims_count', 'NonTranscribed_mean_sims_count',
                'transcribed_versus_untranscribed_p_value', 'transcribed_versus_untranscribed_q_value',
                'Transcribed_real_count.1', 'Transcribed_mean_sims_count.1', 'Transcribed_min_sims_count', 'Transcribed_max_sims_count', 'Transcribed_sims_count_list',
                'UnTranscribed_real_count.1', 'UnTranscribed_mean_sims_count.1', 'UnTranscribed_min_sims_count', 'UnTranscribed_max_sims_count', 'UnTranscribed_sims_count_list',
                'NonTranscribed_real_count.1', 'NonTranscribed_mean_sims_count.1', 'NonTranscribed_min_sims_count', 'NonTranscribed_max_sims_count', 'NonTranscribed_sims_count_list']]

            type_transcribed_versus_untranscribed_df=type_transcribed_versus_untranscribed_df[['cancer_type', 'type',
                'Transcribed_real_count', 'UnTranscribed_real_count', 'NonTranscribed_real_count',
                'Transcribed_mean_sims_count', 'UnTranscribed_mean_sims_count', 'NonTranscribed_mean_sims_count',
                'transcribed_versus_untranscribed_p_value', 'transcribed_versus_untranscribed_q_value',
                'Transcribed_real_count.1', 'Transcribed_mean_sims_count.1', 'Transcribed_min_sims_count', 'Transcribed_max_sims_count', 'Transcribed_sims_count_list',
                'UnTranscribed_real_count.1', 'UnTranscribed_mean_sims_count.1', 'UnTranscribed_min_sims_count', 'UnTranscribed_max_sims_count', 'UnTranscribed_sims_count_list',
                'NonTranscribed_real_count.1', 'NonTranscribed_mean_sims_count.1', 'NonTranscribed_min_sims_count', 'NonTranscribed_max_sims_count', 'NonTranscribed_sims_count_list']]

            signature_filename = 'Signature_Mutation_Type_%s_Q_Value_Table.txt' % (TRANSCRIBED_VERSUS_UNTRANSCRIBED)
            signature_filepath = os.path.join(strandbias_figures_tables_outputDir, signature_filename)
            signature_transcribed_versus_untranscribed_df.to_csv(signature_filepath, sep='\t', header=True, index=False)

            type_filename = 'Type_%s_Q_Value_Table.txt' % (TRANSCRIBED_VERSUS_UNTRANSCRIBED)
            type_filepath = os.path.join(strandbias_figures_tables_outputDir, type_filename)
            type_transcribed_versus_untranscribed_df.to_csv(type_filepath, sep='\t', header=True, index=False)

        if GENIC_VERSUS_INTERGENIC in strand_bias_list:
            signature_genic_versus_intergenic_df=signature_genic_versus_intergenic_df[['cancer_type', 'signature', 'mutation_type',
                'genic_real_count', 'intergenic_real_count', 'genic_mean_sims_count', 'intergenic_mean_sims_count', 'genic_versus_intergenic_p_value', 'genic_versus_intergenic_q_value',
                'Transcribed_real_count', 'Transcribed_mean_sims_count', 'Transcribed_min_sims_count', 'Transcribed_max_sims_count', 'Transcribed_sims_count_list',
                'UnTranscribed_real_count', 'UnTranscribed_mean_sims_count', 'UnTranscribed_min_sims_count', 'UnTranscribed_max_sims_count', 'UnTranscribed_sims_count_list',
                'NonTranscribed_real_count', 'NonTranscribed_mean_sims_count', 'NonTranscribed_min_sims_count', 'NonTranscribed_max_sims_count', 'NonTranscribed_sims_count_list' ]]

            type_genic_versus_intergenic_df=type_genic_versus_intergenic_df[['cancer_type', 'type',
                'genic_real_count', 'intergenic_real_count', 'genic_mean_sims_count', 'intergenic_mean_sims_count', 'genic_versus_intergenic_p_value', 'genic_versus_intergenic_q_value',
                'Transcribed_real_count', 'Transcribed_mean_sims_count', 'Transcribed_min_sims_count', 'Transcribed_max_sims_count', 'Transcribed_sims_count_list',
                'UnTranscribed_real_count', 'UnTranscribed_mean_sims_count', 'UnTranscribed_min_sims_count', 'UnTranscribed_max_sims_count', 'UnTranscribed_sims_count_list',
                'NonTranscribed_real_count', 'NonTranscribed_mean_sims_count', 'NonTranscribed_min_sims_count', 'NonTranscribed_max_sims_count', 'NonTranscribed_sims_count_list' ]]

            signature_filename = 'Signature_Mutation_Type_%s_Q_Value_Table.txt' % (GENIC_VERSUS_INTERGENIC)
            signature_filepath = os.path.join(strandbias_figures_tables_outputDir, signature_filename)
            signature_genic_versus_intergenic_df.to_csv(signature_filepath, sep='\t', header=True, index=False)

            type_filename = 'Type_%s_Q_Value_Table.txt' % (GENIC_VERSUS_INTERGENIC)
            type_filepath = os.path.join(strandbias_figures_tables_outputDir, type_filename)
            type_genic_versus_intergenic_df.to_csv(type_filepath, sep='\t', header=True, index=False)
    #######################################################################

    #######################################################################
    # Step3 Filter q-values, Decide significant strand and set 10,20,30,50,75, 100 percent
    # Add Significant Strand
    # Set significant strands
    # Set percentages
    # Write Filtered Q Values dataframes with percentages
    ##################################################################################################################################
    if LAGGING_VERSUS_LEADING in strand_bias_list:
        signature_lagging_versus_leading_filtered_q_value_df = signature_lagging_versus_leading_df[signature_lagging_versus_leading_df[LAGGING_VERSUS_LEADING_Q_VALUE] <= SIGNIFICANCE_LEVEL].copy()
        type_lagging_versus_leading_filtered_q_value_df= type_lagging_versus_leading_df[type_lagging_versus_leading_df[LAGGING_VERSUS_LEADING_Q_VALUE]<= SIGNIFICANCE_LEVEL].copy()

        signature_lagging_versus_leading_filtered_q_value_df[SIGNIFICANT_STRAND] = None
        type_lagging_versus_leading_filtered_q_value_df[SIGNIFICANT_STRAND] = None

        for percentage_string in percentage_strings:
            signature_lagging_versus_leading_filtered_q_value_df[percentage_string] = None
            type_lagging_versus_leading_filtered_q_value_df[percentage_string] = None

        signature_lagging_versus_leading_filtered_q_value_df.loc[(signature_lagging_versus_leading_filtered_q_value_df[LAGGING_REAL_COUNT] >signature_lagging_versus_leading_filtered_q_value_df[LEADING_REAL_COUNT]), SIGNIFICANT_STRAND] = LAGGING
        signature_lagging_versus_leading_filtered_q_value_df.loc[(signature_lagging_versus_leading_filtered_q_value_df[LEADING_REAL_COUNT] >signature_lagging_versus_leading_filtered_q_value_df[LAGGING_REAL_COUNT]), SIGNIFICANT_STRAND] = LEADING
        type_lagging_versus_leading_filtered_q_value_df.loc[(type_lagging_versus_leading_filtered_q_value_df[LAGGING_REAL_COUNT]>type_lagging_versus_leading_filtered_q_value_df[LEADING_REAL_COUNT]), SIGNIFICANT_STRAND]=LAGGING
        type_lagging_versus_leading_filtered_q_value_df.loc[(type_lagging_versus_leading_filtered_q_value_df[LEADING_REAL_COUNT]>type_lagging_versus_leading_filtered_q_value_df[LAGGING_REAL_COUNT]),SIGNIFICANT_STRAND]=LEADING

        for percentage_index, percentage_number in enumerate(percentage_numbers, 0):
            percentage_string = percentage_strings[percentage_index]
            # Set percentages for signature mutation_type
            signature_lagging_versus_leading_filtered_q_value_df.loc[((signature_lagging_versus_leading_filtered_q_value_df[LAGGING_REAL_COUNT] -signature_lagging_versus_leading_filtered_q_value_df[LEADING_REAL_COUNT]) >= (signature_lagging_versus_leading_filtered_q_value_df[LEADING_REAL_COUNT] * percentage_number / 100)), percentage_string] = 1
            signature_lagging_versus_leading_filtered_q_value_df.loc[((signature_lagging_versus_leading_filtered_q_value_df[LEADING_REAL_COUNT] -signature_lagging_versus_leading_filtered_q_value_df[LAGGING_REAL_COUNT]) >= (signature_lagging_versus_leading_filtered_q_value_df[LAGGING_REAL_COUNT] * percentage_number / 100)), percentage_string] = 1
            # Set percentages for type
            type_lagging_versus_leading_filtered_q_value_df.loc[((type_lagging_versus_leading_filtered_q_value_df[LAGGING_REAL_COUNT] -type_lagging_versus_leading_filtered_q_value_df[LEADING_REAL_COUNT]) >= (type_lagging_versus_leading_filtered_q_value_df[LEADING_REAL_COUNT] * percentage_number / 100)), percentage_string] = 1
            type_lagging_versus_leading_filtered_q_value_df.loc[((type_lagging_versus_leading_filtered_q_value_df[LEADING_REAL_COUNT] -type_lagging_versus_leading_filtered_q_value_df[LAGGING_REAL_COUNT]) >= (type_lagging_versus_leading_filtered_q_value_df[LAGGING_REAL_COUNT] * percentage_number / 100)), percentage_string] = 1

        signature_filename = 'Signature_Mutation_Type_%s_Filtered_Q_Value_Percentages_Table.txt' % (LAGGING_VERSUS_LEADING)
        signature_filepath = os.path.join(strandbias_figures_tables_outputDir, signature_filename)
        signature_lagging_versus_leading_filtered_q_value_df.to_csv(signature_filepath, sep='\t', header=True,index=False)

        type_filename = 'Type_%s_Filtered_Q_Value_Percentages_Table.txt' % (LAGGING_VERSUS_LEADING)
        type_filepath = os.path.join(strandbias_figures_tables_outputDir, type_filename)
        type_lagging_versus_leading_filtered_q_value_df.to_csv(type_filepath, sep='\t', header=True, index=False)
    ##################################################################################################################################

    ##################################################################################################################################
    if TRANSCRIBED_VERSUS_UNTRANSCRIBED in strand_bias_list:
        signature_transcribed_versus_untranscribed_filtered_q_value_df = signature_transcribed_versus_untranscribed_df[signature_transcribed_versus_untranscribed_df[TRANSCRIBED_VERSUS_UNTRANSCRIBED_Q_VALUE] <= SIGNIFICANCE_LEVEL].copy()
        type_transcribed_versus_untranscribed_filtered_q_value_df= type_transcribed_versus_untranscribed_df[type_transcribed_versus_untranscribed_df[TRANSCRIBED_VERSUS_UNTRANSCRIBED_Q_VALUE]<= SIGNIFICANCE_LEVEL].copy()

        signature_transcribed_versus_untranscribed_filtered_q_value_df[SIGNIFICANT_STRAND] = None
        type_transcribed_versus_untranscribed_filtered_q_value_df[SIGNIFICANT_STRAND]=None

        for percentage_string in percentage_strings:
            signature_transcribed_versus_untranscribed_filtered_q_value_df[percentage_string]=None
            type_transcribed_versus_untranscribed_filtered_q_value_df[percentage_string] = None

        signature_transcribed_versus_untranscribed_filtered_q_value_df.loc[(signature_transcribed_versus_untranscribed_filtered_q_value_df[TRANSCRIBED_REAL_COUNT] >signature_transcribed_versus_untranscribed_filtered_q_value_df[UNTRANSCRIBED_REAL_COUNT]), SIGNIFICANT_STRAND] = TRANSCRIBED_STRAND
        signature_transcribed_versus_untranscribed_filtered_q_value_df.loc[(signature_transcribed_versus_untranscribed_filtered_q_value_df[UNTRANSCRIBED_REAL_COUNT] >signature_transcribed_versus_untranscribed_filtered_q_value_df[TRANSCRIBED_REAL_COUNT]), SIGNIFICANT_STRAND] = UNTRANSCRIBED_STRAND
        type_transcribed_versus_untranscribed_filtered_q_value_df.loc[(type_transcribed_versus_untranscribed_filtered_q_value_df[TRANSCRIBED_REAL_COUNT] > type_transcribed_versus_untranscribed_filtered_q_value_df[UNTRANSCRIBED_REAL_COUNT]), SIGNIFICANT_STRAND] = TRANSCRIBED_STRAND
        type_transcribed_versus_untranscribed_filtered_q_value_df.loc[(type_transcribed_versus_untranscribed_filtered_q_value_df[UNTRANSCRIBED_REAL_COUNT] > type_transcribed_versus_untranscribed_filtered_q_value_df[TRANSCRIBED_REAL_COUNT]), SIGNIFICANT_STRAND] = UNTRANSCRIBED_STRAND

        for percentage_index, percentage_number in enumerate(percentage_numbers,0):
            percentage_string=percentage_strings[percentage_index]
            # Set percentages for signature mutation_type
            signature_transcribed_versus_untranscribed_filtered_q_value_df.loc[((signature_transcribed_versus_untranscribed_filtered_q_value_df[TRANSCRIBED_REAL_COUNT]-signature_transcribed_versus_untranscribed_filtered_q_value_df[UNTRANSCRIBED_REAL_COUNT]) >= (signature_transcribed_versus_untranscribed_filtered_q_value_df[UNTRANSCRIBED_REAL_COUNT]*percentage_number/100)), percentage_string] = 1
            signature_transcribed_versus_untranscribed_filtered_q_value_df.loc[((signature_transcribed_versus_untranscribed_filtered_q_value_df[UNTRANSCRIBED_REAL_COUNT]-signature_transcribed_versus_untranscribed_filtered_q_value_df[TRANSCRIBED_REAL_COUNT]) >= (signature_transcribed_versus_untranscribed_filtered_q_value_df[TRANSCRIBED_REAL_COUNT]*percentage_number/100)), percentage_string] = 1
            # Set percentages for type
            type_transcribed_versus_untranscribed_filtered_q_value_df.loc[((type_transcribed_versus_untranscribed_filtered_q_value_df[TRANSCRIBED_REAL_COUNT]-type_transcribed_versus_untranscribed_filtered_q_value_df[UNTRANSCRIBED_REAL_COUNT]) >= (type_transcribed_versus_untranscribed_filtered_q_value_df[UNTRANSCRIBED_REAL_COUNT]*percentage_number/100)), percentage_string] = 1
            type_transcribed_versus_untranscribed_filtered_q_value_df.loc[((type_transcribed_versus_untranscribed_filtered_q_value_df[UNTRANSCRIBED_REAL_COUNT]-type_transcribed_versus_untranscribed_filtered_q_value_df[TRANSCRIBED_REAL_COUNT]) >= (type_transcribed_versus_untranscribed_filtered_q_value_df[TRANSCRIBED_REAL_COUNT]*percentage_number/100)), percentage_string] = 1

        signature_filename = 'Signature_Mutation_Type_%s_Filtered_Q_Value_Percentages_Table.txt' % (TRANSCRIBED_VERSUS_UNTRANSCRIBED)
        signature_filepath = os.path.join(strandbias_figures_tables_outputDir, signature_filename)
        signature_transcribed_versus_untranscribed_filtered_q_value_df.to_csv(signature_filepath, sep='\t', header=True, index=False)

        type_filename = 'Type_%s_Filtered_Q_Value_Percentages_Table.txt' % (TRANSCRIBED_VERSUS_UNTRANSCRIBED)
        type_filepath = os.path.join(strandbias_figures_tables_outputDir, type_filename)
        type_transcribed_versus_untranscribed_filtered_q_value_df.to_csv(type_filepath, sep='\t', header=True,index=False)
    ##################################################################################################################################


    ##################################################################################################################################
    if GENIC_VERSUS_INTERGENIC in strand_bias_list:
        signature_genic_versus_intergenic_filtered_q_value_df = signature_genic_versus_intergenic_df[signature_genic_versus_intergenic_df[GENIC_VERSUS_INTERGENIC_Q_VALUE] <= SIGNIFICANCE_LEVEL].copy()
        type_genic_versus_intergenic_filtered_q_value_df= type_genic_versus_intergenic_df[type_genic_versus_intergenic_df[GENIC_VERSUS_INTERGENIC_Q_VALUE]<= SIGNIFICANCE_LEVEL].copy()

        signature_genic_versus_intergenic_filtered_q_value_df[SIGNIFICANT_STRAND] = None
        type_genic_versus_intergenic_filtered_q_value_df[SIGNIFICANT_STRAND] = None

        for percentage_string in percentage_strings:
            signature_genic_versus_intergenic_filtered_q_value_df[percentage_string] = None
            type_genic_versus_intergenic_filtered_q_value_df[percentage_string] = None

        signature_genic_versus_intergenic_filtered_q_value_df.loc[(signature_genic_versus_intergenic_filtered_q_value_df[GENIC_REAL_COUNT] > signature_genic_versus_intergenic_filtered_q_value_df[INTERGENIC_REAL_COUNT]), SIGNIFICANT_STRAND] = GENIC
        signature_genic_versus_intergenic_filtered_q_value_df.loc[(signature_genic_versus_intergenic_filtered_q_value_df[INTERGENIC_REAL_COUNT] > signature_genic_versus_intergenic_filtered_q_value_df[GENIC_REAL_COUNT]),SIGNIFICANT_STRAND] = INTERGENIC
        type_genic_versus_intergenic_filtered_q_value_df.loc[(type_genic_versus_intergenic_filtered_q_value_df[GENIC_REAL_COUNT] > type_genic_versus_intergenic_filtered_q_value_df[INTERGENIC_REAL_COUNT]), SIGNIFICANT_STRAND] = GENIC
        type_genic_versus_intergenic_filtered_q_value_df.loc[(type_genic_versus_intergenic_filtered_q_value_df[INTERGENIC_REAL_COUNT] > type_genic_versus_intergenic_filtered_q_value_df[GENIC_REAL_COUNT]), SIGNIFICANT_STRAND] = INTERGENIC

        #Set percentages
        for percentage_index, percentage_number in enumerate(percentage_numbers,0):
            percentage_string=percentage_strings[percentage_index]
            # Set percentages for signature mutation_type
            signature_genic_versus_intergenic_filtered_q_value_df.loc[((signature_genic_versus_intergenic_filtered_q_value_df[GENIC_REAL_COUNT]-signature_genic_versus_intergenic_filtered_q_value_df[INTERGENIC_REAL_COUNT]) >= (signature_genic_versus_intergenic_filtered_q_value_df[INTERGENIC_REAL_COUNT]*percentage_number/100)), percentage_string] = 1
            signature_genic_versus_intergenic_filtered_q_value_df.loc[((signature_genic_versus_intergenic_filtered_q_value_df[INTERGENIC_REAL_COUNT]-signature_genic_versus_intergenic_filtered_q_value_df[GENIC_REAL_COUNT]) >= (signature_genic_versus_intergenic_filtered_q_value_df[GENIC_REAL_COUNT]*percentage_number/100)), percentage_string] = 1
            # Set percentages for type
            type_genic_versus_intergenic_filtered_q_value_df.loc[((type_genic_versus_intergenic_filtered_q_value_df[GENIC_REAL_COUNT]-type_genic_versus_intergenic_filtered_q_value_df[INTERGENIC_REAL_COUNT]) >= (type_genic_versus_intergenic_filtered_q_value_df[INTERGENIC_REAL_COUNT]*percentage_number/100)), percentage_string] = 1
            type_genic_versus_intergenic_filtered_q_value_df.loc[((type_genic_versus_intergenic_filtered_q_value_df[INTERGENIC_REAL_COUNT]-type_genic_versus_intergenic_filtered_q_value_df[GENIC_REAL_COUNT]) >= (type_genic_versus_intergenic_filtered_q_value_df[GENIC_REAL_COUNT]*percentage_number/100)), percentage_string] = 1

        signature_filename = 'Signature_Mutation_Type_%s_Filtered_Q_Value_Percentages_Table.txt' % (GENIC_VERSUS_INTERGENIC)
        signature_filepath = os.path.join(strandbias_figures_tables_outputDir, signature_filename)
        signature_genic_versus_intergenic_filtered_q_value_df.to_csv(signature_filepath, sep='\t', header=True,index=False)

        type_filename = 'Type_%s_Filtered_Q_Value_Percentages_Table.txt' % (GENIC_VERSUS_INTERGENIC)
        type_filepath = os.path.join(strandbias_figures_tables_outputDir, type_filename)
        type_genic_versus_intergenic_filtered_q_value_df.to_csv(type_filepath, sep='\t', header=True, index=False)
    ##################################################################################################################################

    #######################################################################
    # Write Excel Files
    sheet_list = ['corrected_p_value', 'percentages']
    for strand1_versus_strand2 in strand_bias_list:
        if strand1_versus_strand2==LAGGING_VERSUS_LEADING:
            signatures_df_list=[signature_lagging_versus_leading_df,signature_lagging_versus_leading_filtered_q_value_df]
            types_df_list = [type_lagging_versus_leading_df, type_lagging_versus_leading_filtered_q_value_df]
        elif strand1_versus_strand2==TRANSCRIBED_VERSUS_UNTRANSCRIBED:
            signatures_df_list = [signature_transcribed_versus_untranscribed_df,signature_transcribed_versus_untranscribed_filtered_q_value_df]
            types_df_list = [type_transcribed_versus_untranscribed_df, type_transcribed_versus_untranscribed_filtered_q_value_df]
        elif strand1_versus_strand2==GENIC_VERSUS_INTERGENIC:
            signatures_df_list = [signature_genic_versus_intergenic_df,signature_genic_versus_intergenic_filtered_q_value_df]
            types_df_list = [type_genic_versus_intergenic_df, type_genic_versus_intergenic_filtered_q_value_df]

        signatures_filename="Signatures_Mutation_Types_%s.xlsx" %(strand1_versus_strand2)
        file_name_with_path=os.path.join(strandbias_figures_excel_files_outputDir, signatures_filename)
        write_excel_file(signatures_df_list, sheet_list, file_name_with_path)

        types_filename="Types_%s.xlsx" %(strand1_versus_strand2)
        file_name_with_path=os.path.join(strandbias_figures_excel_files_outputDir, types_filename)
        write_excel_file(types_df_list, sheet_list, file_name_with_path)
    #######################################################################


    #######################################################################
    #Circle plots starts

    #######################################################################
    #Step4 Fill this dictionary
    signature2mutation_type2strand2percentagedict={}

    df_list=[]
    if LAGGING_VERSUS_LEADING in strand_bias_list:
        df_list.append(signature_lagging_versus_leading_filtered_q_value_df)
    if TRANSCRIBED_VERSUS_UNTRANSCRIBED in strand_bias_list:
        df_list.append(signature_transcribed_versus_untranscribed_filtered_q_value_df)
    if GENIC_VERSUS_INTERGENIC in strand_bias_list:
        df_list.append(signature_genic_versus_intergenic_filtered_q_value_df)

    for df in df_list:
        for index, row in df.iterrows():
            cancer_type = row[CANCER_TYPE]
            signature = row[SIGNATURE]
            mutation_type = row[MUTATION_TYPE]
            significant_strand=row[SIGNIFICANT_STRAND]
            percent_10 = row[AT_LEAST_10_PERCENT_DIFF]
            percent_20 = row[AT_LEAST_20_PERCENT_DIFF]
            percent_30 = row[AT_LEAST_30_PERCENT_DIFF]
            percent_50 = row[AT_LEAST_50_PERCENT_DIFF]
            percent_75 = row[AT_LEAST_75_PERCENT_DIFF]
            percent_100 = row[AT_LEAST_100_PERCENT_DIFF]

            if signature in signature2mutation_type2strand2percentagedict:
                if mutation_type in signature2mutation_type2strand2percentagedict[signature]:
                    if significant_strand in signature2mutation_type2strand2percentagedict[signature][mutation_type]:
                        if (percent_10 == 1):
                            signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_10_PERCENT_DIFF]=1
                        if (percent_20 == 1):
                            signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_20_PERCENT_DIFF]=1
                        if (percent_30 == 1):
                            signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_30_PERCENT_DIFF]=1
                        if (percent_50 == 1):
                            signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_50_PERCENT_DIFF]=1
                        if (percent_75 == 1):
                            signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_75_PERCENT_DIFF]=1
                        if (percent_100 == 1):
                            signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_100_PERCENT_DIFF]=1

                    else:
                        signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand]={}
                        signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_10_PERCENT_DIFF] = 0
                        signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_20_PERCENT_DIFF] = 0
                        signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_30_PERCENT_DIFF] = 0
                        signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_50_PERCENT_DIFF] = 0
                        signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_75_PERCENT_DIFF] = 0
                        signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_100_PERCENT_DIFF] = 0

                        if (percent_10 == 1):
                            signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_10_PERCENT_DIFF] = 1
                        if (percent_20 == 1):
                            signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_20_PERCENT_DIFF] = 1
                        if (percent_30 == 1):
                            signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_30_PERCENT_DIFF] = 1
                        if (percent_50 == 1):
                            signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_50_PERCENT_DIFF] = 1
                        if (percent_75 == 1):
                            signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_75_PERCENT_DIFF] = 1
                        if (percent_100 == 1):
                            signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_100_PERCENT_DIFF] = 1
                else:
                    signature2mutation_type2strand2percentagedict[signature][mutation_type] = {}
                    signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand] = {}

                    signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_10_PERCENT_DIFF] = 0
                    signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_20_PERCENT_DIFF] = 0
                    signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_30_PERCENT_DIFF] = 0
                    signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_50_PERCENT_DIFF] = 0
                    signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_75_PERCENT_DIFF] = 0
                    signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_100_PERCENT_DIFF] = 0

                    if (percent_10 == 1):
                        signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_10_PERCENT_DIFF] = 1
                    if (percent_20 == 1):
                        signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_20_PERCENT_DIFF] = 1
                    if (percent_30 == 1):
                        signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_30_PERCENT_DIFF] = 1
                    if (percent_50 == 1):
                        signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_50_PERCENT_DIFF] = 1
                    if (percent_75 == 1):
                        signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_75_PERCENT_DIFF] = 1
                    if (percent_100 == 1):
                        signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_100_PERCENT_DIFF] = 1

            else:
                signature2mutation_type2strand2percentagedict[signature] = {}
                signature2mutation_type2strand2percentagedict[signature][mutation_type] = {}
                signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand] = {}

                signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_10_PERCENT_DIFF] = 0
                signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_20_PERCENT_DIFF] = 0
                signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_30_PERCENT_DIFF] = 0
                signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_50_PERCENT_DIFF] = 0
                signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_75_PERCENT_DIFF] = 0
                signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_100_PERCENT_DIFF] = 0

                if (percent_10 == 1):
                    signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_10_PERCENT_DIFF] = 1
                if (percent_20 == 1):
                    signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_20_PERCENT_DIFF] = 1
                if (percent_30 == 1):
                    signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_30_PERCENT_DIFF] = 1
                if (percent_50 == 1):
                    signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_50_PERCENT_DIFF] = 1
                if (percent_75 == 1):
                    signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_75_PERCENT_DIFF] = 1
                if (percent_100 == 1):
                    signature2mutation_type2strand2percentagedict[signature][mutation_type][significant_strand][AT_LEAST_100_PERCENT_DIFF] = 1
    #######################################################################

    #######################################################################
    #Step4 Fill this dictionary
    type2strand2percentagedict={}

    df_list=[]
    if LAGGING_VERSUS_LEADING in strand_bias_list:
        df_list.append(type_lagging_versus_leading_filtered_q_value_df)
    if TRANSCRIBED_VERSUS_UNTRANSCRIBED in strand_bias_list:
        df_list.append(type_transcribed_versus_untranscribed_filtered_q_value_df)
    if GENIC_VERSUS_INTERGENIC in strand_bias_list:
        df_list.append(type_genic_versus_intergenic_filtered_q_value_df)

    for df in df_list:
        for index, row in df.iterrows():
            cancer_type = row[CANCER_TYPE]
            my_type = row[TYPE]
            significant_strand=row[SIGNIFICANT_STRAND]
            percent_10 = row[AT_LEAST_10_PERCENT_DIFF]
            percent_20 = row[AT_LEAST_20_PERCENT_DIFF]
            percent_30 = row[AT_LEAST_30_PERCENT_DIFF]
            percent_50 = row[AT_LEAST_50_PERCENT_DIFF]
            percent_75 = row[AT_LEAST_75_PERCENT_DIFF]
            percent_100 = row[AT_LEAST_100_PERCENT_DIFF]

            if my_type in type2strand2percentagedict:
                if significant_strand in type2strand2percentagedict[my_type]:
                    if (percent_10 == 1):
                        type2strand2percentagedict[my_type][significant_strand][AT_LEAST_10_PERCENT_DIFF]=1
                    if (percent_20 == 1):
                        type2strand2percentagedict[my_type][significant_strand][AT_LEAST_20_PERCENT_DIFF]=1
                    if (percent_30 == 1):
                        type2strand2percentagedict[my_type][significant_strand][AT_LEAST_30_PERCENT_DIFF]=1
                    if (percent_50 == 1):
                        type2strand2percentagedict[my_type][significant_strand][AT_LEAST_50_PERCENT_DIFF]=1
                    if (percent_75 == 1):
                        type2strand2percentagedict[my_type][significant_strand][AT_LEAST_75_PERCENT_DIFF]=1
                    if (percent_100 == 1):
                        type2strand2percentagedict[my_type][significant_strand][AT_LEAST_100_PERCENT_DIFF]=1

                else:
                    type2strand2percentagedict[my_type][significant_strand]={}
                    type2strand2percentagedict[my_type][significant_strand][AT_LEAST_10_PERCENT_DIFF] = 0
                    type2strand2percentagedict[my_type][significant_strand][AT_LEAST_20_PERCENT_DIFF] = 0
                    type2strand2percentagedict[my_type][significant_strand][AT_LEAST_30_PERCENT_DIFF] = 0
                    type2strand2percentagedict[my_type][significant_strand][AT_LEAST_50_PERCENT_DIFF] = 0
                    type2strand2percentagedict[my_type][significant_strand][AT_LEAST_75_PERCENT_DIFF] = 0
                    type2strand2percentagedict[my_type][significant_strand][AT_LEAST_100_PERCENT_DIFF] = 0

                    if (percent_10 == 1):
                        type2strand2percentagedict[my_type][significant_strand][AT_LEAST_10_PERCENT_DIFF]=1
                    if (percent_20 == 1):
                        type2strand2percentagedict[my_type][significant_strand][AT_LEAST_20_PERCENT_DIFF]=1
                    if (percent_30 == 1):
                        type2strand2percentagedict[my_type][significant_strand][AT_LEAST_30_PERCENT_DIFF]=1
                    if (percent_50 == 1):
                        type2strand2percentagedict[my_type][significant_strand][AT_LEAST_50_PERCENT_DIFF]=1
                    if (percent_75 == 1):
                        type2strand2percentagedict[my_type][significant_strand][AT_LEAST_75_PERCENT_DIFF]=1
                    if (percent_100 == 1):
                        type2strand2percentagedict[my_type][significant_strand][AT_LEAST_100_PERCENT_DIFF]=1
            else:
                type2strand2percentagedict[my_type] = {}
                type2strand2percentagedict[my_type][significant_strand] = {}
                type2strand2percentagedict[my_type][significant_strand][AT_LEAST_10_PERCENT_DIFF] = 0
                type2strand2percentagedict[my_type][significant_strand][AT_LEAST_20_PERCENT_DIFF] = 0
                type2strand2percentagedict[my_type][significant_strand][AT_LEAST_30_PERCENT_DIFF] = 0
                type2strand2percentagedict[my_type][significant_strand][AT_LEAST_50_PERCENT_DIFF] = 0
                type2strand2percentagedict[my_type][significant_strand][AT_LEAST_75_PERCENT_DIFF] = 0
                type2strand2percentagedict[my_type][significant_strand][AT_LEAST_100_PERCENT_DIFF] = 0

                if (percent_10 == 1):
                    type2strand2percentagedict[my_type][significant_strand][AT_LEAST_10_PERCENT_DIFF]=1
                if (percent_20 == 1):
                    type2strand2percentagedict[my_type][significant_strand][AT_LEAST_20_PERCENT_DIFF]=1
                if (percent_30 == 1):
                    type2strand2percentagedict[my_type][significant_strand][AT_LEAST_30_PERCENT_DIFF]=1
                if (percent_50 == 1):
                    type2strand2percentagedict[my_type][significant_strand][AT_LEAST_50_PERCENT_DIFF]=1
                if (percent_75 == 1):
                    type2strand2percentagedict[my_type][significant_strand][AT_LEAST_75_PERCENT_DIFF]=1
                if (percent_100 == 1):
                    type2strand2percentagedict[my_type][significant_strand][AT_LEAST_100_PERCENT_DIFF]=1
    #######################################################################


    #######################################################################
    #Step5 Plot figures
    plot_legend(strandbias_figures_outputDir)

    for strand_bias in strand_bias_list:
        if np.any(subsSignatures):
            plot_six_mutations_sbs_signatures_circle_figures(subsSignatures,
                                        strand_bias,
                                        strandbias_figures_outputDir,
                                        SIGNIFICANCE_LEVEL,
                                        signature2mutation_type2strand2percentagedict,
                                        percentage_strings)
        if np.any(dinucsSignatures):
            plot_dbs_and_id_signatures_circle_figures(DBS,
                                           dinucsSignatures,
                                           strand_bias,
                                           strandbias_figures_outputDir,
                                           SIGNIFICANCE_LEVEL,
                                           type2strand2percentagedict,
                                           percentage_strings)
        if np.any(indelsSignatures):
            plot_dbs_and_id_signatures_circle_figures(ID,
                                           indelsSignatures,
                                           strand_bias,
                                           strandbias_figures_outputDir,
                                           SIGNIFICANCE_LEVEL,
                                           type2strand2percentagedict,
                                           percentage_strings)

    #Circle plots ends
    #######################################################################

    ########################################################################
    ##########################  Part 2 starts ##############################
    ############## Mutation Types Scatter Plots starts #####################
    ############## Signatures Scatter Plots starts #########################
    ########################################################################
    if (TRANSCRIBED_VERSUS_UNTRANSCRIBED in strand_bias_list) and (LAGGING_VERSUS_LEADING in strand_bias_list):
        if ((not type_transcribed_versus_untranscribed_df.empty) and (not type_lagging_versus_leading_df.empty)):
            plot_mutation_types_transcription_log10_ratio_replication_log_10_ratio_using_dataframes(None,None,
                                                                                                    type_transcribed_versus_untranscribed_df,
                                                                                                    type_lagging_versus_leading_df,
                                                                                                    outputDir, jobname)

        if ((not type_transcribed_versus_untranscribed_df.empty) and (not type_lagging_versus_leading_df.empty) and (not sbs_df.empty)):
            plot_types_transcription_log10_ratio_replication_log10_ratio_using_dataframes('subs', None, None,
                                                                                           type_transcribed_versus_untranscribed_df,
                                                                                           type_lagging_versus_leading_df,
                                                                                           sbs_df,
                                                                                           outputDir, jobname)

        if ((not type_transcribed_versus_untranscribed_df.empty) and (not type_lagging_versus_leading_df.empty) and (not dbs_df.empty)):
            plot_types_transcription_log10_ratio_replication_log10_ratio_using_dataframes('dinucs', None, None,
                                                                                           type_transcribed_versus_untranscribed_df,
                                                                                           type_lagging_versus_leading_df,
                                                                                           dbs_df,
                                                                                           outputDir, jobname)

        if ((not type_transcribed_versus_untranscribed_df.empty) and (not type_lagging_versus_leading_df.empty) and (not id_df.empty)):
            plot_types_transcription_log10_ratio_replication_log10_ratio_using_dataframes('indels', None, None,
                                                                                           type_transcribed_versus_untranscribed_df,
                                                                                           type_lagging_versus_leading_df,
                                                                                           id_df,
                                                                                           outputDir, jobname)
    ########################################################################
    ############## Mutation Types Scatter Plots ends #######################
    ############## Signatures Scatter Plots ends ###########################
    ##########################  Part 2 ends ################################
    ########################################################################

    ########################################################################
    ##########################  Part 4 starts ##############################
    ######## Bar plot starts includes sample based bar plots ###############
    ########################################################################
    isKeySample = False
    width = 0.20

    #######################################################
    ################# Plot types starts ###################
    #######################################################
    types_list= [('All Mutations', 'mutationtypes', six_mutation_types),
                 ('All Signatures', 'subs_signatures', subsSignatures),
                 ('All Signatures', 'indels_signatures', indelsSignatures),
                 ('All Signatures', 'dinucs_signatures', dinucsSignatures)]

    for mutationsOrSignatures, sub_figure_name, x_axis_labels in types_list:
        x_axis_labels = sorted(x_axis_labels, key=natural_key)
        N = len(x_axis_labels)

        for strand_bias in strand_bias_list:
            if (strand_bias == TRANSCRIBED_VERSUS_UNTRANSCRIBED):
                type_strand1_versus_strand2_df=type_transcribed_versus_untranscribed_df
                strand1=transcriptionStrands[0]
                strand2=transcriptionStrands[1]
                strand1_real_count_column_name='Transcribed_real_count'
                strand2_real_count_column_name='UnTranscribed_real_count'
                strand1_sims_mean_count_column_name='Transcribed_mean_sims_count'
                strand2_sims_mean_count_column_name='UnTranscribed_mean_sims_count'
                q_value_column_name='transcribed_versus_untranscribed_q_value'
                color1='royalblue'
                color2='yellowgreen'
                figureName='%s_transcription_strand_bias' %(sub_figure_name)

            elif (strand_bias == GENIC_VERSUS_INTERGENIC):
                type_strand1_versus_strand2_df=type_genic_versus_intergenic_df
                strand1=genicVersusIntergenicStrands[0]
                strand2=genicVersusIntergenicStrands[1]
                strand1_real_count_column_name='genic_real_count'
                strand2_real_count_column_name='intergenic_real_count'
                strand1_sims_mean_count_column_name='genic_mean_sims_count'
                strand2_sims_mean_count_column_name='intergenic_mean_sims_count'
                q_value_column_name='genic_versus_intergenic_q_value'
                color1='cyan'
                color2='gray'
                figureName='%s_genic_versus_intergenic_strand_bias' %(sub_figure_name)

            elif (strand_bias == LAGGING_VERSUS_LEADING):
                type_strand1_versus_strand2_df = type_lagging_versus_leading_df
                strand1=replicationStrands[0]
                strand2=replicationStrands[1]
                strand1_real_count_column_name='Lagging_real_count'
                strand2_real_count_column_name='Leading_real_count'
                strand1_sims_mean_count_column_name='Lagging_mean_sims_count'
                strand2_sims_mean_count_column_name='Leading_mean_sims_count'
                q_value_column_name='lagging_versus_leading_q_value'
                color1='indianred'
                color2='goldenrod'
                figureName='%s_replication_strand_bias' %(sub_figure_name)

            types_strand1_real_count_list=[]
            types_strand2_real_count_list=[]
            types_strand1_sims_mean_count_list=[]
            types_strand2_sims_mean_count_list=[]
            types_strand1_versus_strand2_FDR_BH_adjusted_pvalues=[]

            for my_type in x_axis_labels:
                strand1_real_count=0
                strand2_real_count=0
                strand1_sims_mean_count=0
                strand2_sims_mean_count=0
                q_value=None

                if type_strand1_versus_strand2_df[(type_strand1_versus_strand2_df['type']==my_type)][strand1_real_count_column_name].values.size>0:
                    strand1_real_count= type_strand1_versus_strand2_df[(type_strand1_versus_strand2_df['type'] == my_type)][strand1_real_count_column_name].values[0]
                if type_strand1_versus_strand2_df[(type_strand1_versus_strand2_df['type']==my_type)][strand2_real_count_column_name].values.size>0:
                    strand2_real_count= type_strand1_versus_strand2_df[(type_strand1_versus_strand2_df['type'] == my_type)][strand2_real_count_column_name].values[0]
                if type_strand1_versus_strand2_df[(type_strand1_versus_strand2_df['type']==my_type)][strand1_sims_mean_count_column_name].values.size>0:
                    strand1_sims_mean_count= type_strand1_versus_strand2_df[(type_strand1_versus_strand2_df['type'] == my_type)][strand1_sims_mean_count_column_name].values[0]
                if type_strand1_versus_strand2_df[(type_strand1_versus_strand2_df['type']==my_type)][strand2_sims_mean_count_column_name].values.size>0:
                    strand2_sims_mean_count= type_strand1_versus_strand2_df[(type_strand1_versus_strand2_df['type'] == my_type)][strand2_sims_mean_count_column_name].values[0]
                if type_strand1_versus_strand2_df[(type_strand1_versus_strand2_df['type']==my_type)][q_value_column_name].values.size>0:
                    q_value= type_strand1_versus_strand2_df[(type_strand1_versus_strand2_df['type'] == my_type)][q_value_column_name].values[0]

                types_strand1_real_count_list.append(strand1_real_count)
                types_strand2_real_count_list.append(strand2_real_count)
                types_strand1_sims_mean_count_list.append(strand1_sims_mean_count)
                types_strand2_sims_mean_count_list.append(strand2_sims_mean_count)
                types_strand1_versus_strand2_FDR_BH_adjusted_pvalues.append(q_value)

            if ((len(x_axis_labels)>0) and types_strand1_real_count_list and types_strand2_real_count_list and types_strand1_sims_mean_count_list and types_strand2_sims_mean_count_list and (len(types_strand1_versus_strand2_FDR_BH_adjusted_pvalues)>0)):

                if (types_strand1_real_count_list and types_strand2_real_count_list):
                    plotStrandBiasFigureWithBarPlots(outputDir,
                                                     jobname,
                                                     numberofSimulations,
                                                     None,
                                                     isKeySample,
                                                     None,
                                                     N,
                                                     x_axis_labels,
                                                     types_strand1_real_count_list,
                                                     types_strand2_real_count_list,
                                                     types_strand1_sims_mean_count_list,
                                                     types_strand2_sims_mean_count_list,
                                                     types_strand1_versus_strand2_FDR_BH_adjusted_pvalues,
                                                     strand1,strand2,
                                                     mutationsOrSignatures,
                                                     color1, color2,
                                                     figureName,
                                                     width,
                                                     plot_mode)

    #######################################################
    ################# Plot types ends #####################
    #######################################################

    #################################################################
    ########### Plot sub signatures mutation types starts ###########
    #################################################################
    if not sbs_df.empty:
        if TRANSCRIBED_VERSUS_UNTRANSCRIBED in strand_bias_list:
            plotBarPlotsUsingDataframes(outputDir,
                     jobname,
                     numberofSimulations,
                     sbs_df,
                     isKeySample,
                     six_mutation_types,
                     signature_transcribed_versus_untranscribed_df,
                     width,
                     TRANSCRIBED_VERSUS_UNTRANSCRIBED,
                     transcriptionStrands,
                     'royalblue',
                     'yellowgreen',
                     'All Mutations',
                     'mutationtypes_transcription_strand_bias',
                      plot_mode)

        if GENIC_VERSUS_INTERGENIC in strand_bias_list:
            plotBarPlotsUsingDataframes(outputDir,
                     jobname,
                     numberofSimulations,
                     sbs_df,
                     isKeySample,
                     six_mutation_types,
                    signature_genic_versus_intergenic_df,
                     width,
                     GENIC_VERSUS_INTERGENIC,
                    genicVersusIntergenicStrands,
                     'cyan',
                     'gray',
                     'All Mutations',
                     'mutationtypes_genic_versus_intergenic_strand_bias',
                      plot_mode)

        if LAGGING_VERSUS_LEADING in strand_bias_list:
            plotBarPlotsUsingDataframes(outputDir,
                     jobname,
                     numberofSimulations,
                     sbs_df,
                     isKeySample,
                     six_mutation_types,
                    signature_lagging_versus_leading_df,
                     width,
                     LAGGING_VERSUS_LEADING,
                    replicationStrands,
                    'indianred',
                    'goldenrod',
                    'All Mutations',
                    'mutationtypes_replication_strand_bias',
                     plot_mode)
    #################################################################
    ########### Plot sub signatures mutation types ends #############
    #################################################################

    ########################################################################
    ######## Bar plot starts includes sample based bar plots ###############
    ##########################  Part 4 ends ################################
    ########################################################################

    # Circle Bar Plots
    # Plot circle plots and bar plots all together
    # At top ax, circle plots with 3 rows: for genic vs. intergenic, transcribed vs. untranscribed, lagging vs. leading
    # At middle ax, 3 bar plots: for genic vs. intergenic, transcribed vs. untranscribed, lagging vs. leading
    # At below ax, 3 normalized bar plots: for genic vs. intergenic, transcribed vs. untranscribed, lagging vs. leading
    if (TRANSCRIBED_VERSUS_UNTRANSCRIBED in strand_bias_list) and (LAGGING_VERSUS_LEADING in strand_bias_list):
        sbs_signatures = sbs_df['signature'].unique()
        for sbs_signature in sbs_signatures:
            plot_circle_bar_plots_together(outputDir,
                                           jobname,
                                           sbs_signature,
                                            six_mutation_types,
                                           signature2mutation_type2strand2percentagedict,
                                            signature_genic_versus_intergenic_df,
                                            signature_transcribed_versus_untranscribed_df,
                                            signature_lagging_versus_leading_df,
                                            genicVersusIntergenicStrands,
                                            transcriptionStrands,
                                            replicationStrands)

###################################################################


############################################################################################################################
def plot_dbs_and_id_signatures_circle_figures(signature_type,
                                       signatures,
                                       strand_bias,
                                       strandbias_figures_outputDir,
                                       SIGNIFICANCE_LEVEL,
                                       type2strand2percentagedict,
                                       percentage_strings):

    rows_signatures=[]

    #####################################################################
    if strand_bias==LAGGING_VERSUS_LEADING:
        strands=replicationStrands
    elif strand_bias==TRANSCRIBED_VERSUS_UNTRANSCRIBED:
        strands=transcriptionStrands
    elif strand_bias==GENIC_VERSUS_INTERGENIC:
        strands=genicVersusIntergenicStrands
    #####################################################################

    #####################################################################
    #Fill rows_DBS_signatures
    #Fill rows_ID_signatures
    for signature in signatures:
        if signature  in type2strand2percentagedict:
            for strand in strands:
                if strand in type2strand2percentagedict[signature]:
                    for percentage_string in percentage_strings:
                        if percentage_string in type2strand2percentagedict[signature][strand]:
                            print('signature:%s strand:%s percentage_string:%s' %(signature,strand,percentage_string))
                            if signature not in rows_signatures:
                                rows_signatures.append(signature)
    #####################################################################

    #####################################################################
    rows_signatures=sorted(rows_signatures,key=natural_key,reverse=True)
    #####################################################################


    if (len(rows_signatures)>0):
        #####################################################################
        #New plot (width,height)
        fig, ax = plt.subplots(figsize=(5+1.5*len(percentage_strings), 10+1.5*len(rows_signatures)))

        #make aspect ratio square
        ax.set_aspect(1.0)
        #####################################################################

        ######################################################################################################################################
        for percentage_diff_index, percentage_string in enumerate(percentage_strings):
            for row_signature_index, row_signature in enumerate(rows_signatures):
                if (strand_bias==LAGGING_VERSUS_LEADING):
                    if row_signature in type2strand2percentagedict:
                        lagging_percentage=None
                        leading_percentage=None

                        if LAGGING in type2strand2percentagedict[row_signature] and type2strand2percentagedict[row_signature][LAGGING][percentage_string]==1:
                            lagging_percentage = 100
                        if LEADING in type2strand2percentagedict[row_signature] and type2strand2percentagedict[row_signature][LEADING][percentage_string]==1:
                            leading_percentage = 100

                        if (lagging_percentage is not None) and (leading_percentage is None):
                            radius = 0.49
                            if (radius > 0):
                                # print('Plot circle at x=%d y=%d for %s %s' % (percentage_diff_index, row_signature_index, row_signature,percentage_string))
                                circle = plt.Circle((percentage_diff_index + 0.5, row_signature_index + 0.5), radius, color='indianred', fill=True)
                                ax.add_artist(circle)

                        elif (leading_percentage is not None) and (lagging_percentage is None):
                            radius = 0.49
                            if (radius > 0):
                                # print('Plot circle at x=%d y=%d for %s %s' % (percentage_diff_index, row_signature_index, row_signature,percentage_string))
                                circle = plt.Circle((percentage_diff_index + 0.5, row_signature_index + 0.5), radius, color='goldenrod', fill=True)
                                ax.add_artist(circle)

                        elif (lagging_percentage is not None) and (leading_percentage is not None):
                            radius_lagging = 0.49
                            radius_leading = 0.49
                            if (radius_lagging>radius_leading):
                                #First lagging
                                circle = plt.Circle((percentage_diff_index + 0.5, row_signature_index + 0.5), radius_lagging, color='goldenrod', fill=True)
                                ax.add_artist(circle)
                                #Second leading
                                circle = plt.Circle((percentage_diff_index + 0.5, row_signature_index + 0.5), radius_leading, color='goldenrod', fill=True)
                                ax.add_artist(circle)
                            else:
                                #First leading
                                circle = plt.Circle((percentage_diff_index + 0.5, row_signature_index + 0.5), radius_leading, color='goldenrod', fill=True)
                                ax.add_artist(circle)
                                #Second lagging
                                circle = plt.Circle((percentage_diff_index + 0.5, row_signature_index + 0.5), radius_lagging, color='goldenrod', fill=True)
                                ax.add_artist(circle)

                elif (strand_bias==TRANSCRIBED_VERSUS_UNTRANSCRIBED):
                    if row_signature in type2strand2percentagedict:
                        transcribed_percentage=None
                        untranscribed_percentage=None

                        if TRANSCRIBED_STRAND in type2strand2percentagedict[row_signature] and type2strand2percentagedict[row_signature][TRANSCRIBED_STRAND][percentage_string]==1:
                            transcribed_percentage = 100
                        if UNTRANSCRIBED_STRAND in type2strand2percentagedict[row_signature] and type2strand2percentagedict[row_signature][UNTRANSCRIBED_STRAND][percentage_string]==1:
                            untranscribed_percentage = 100

                        if (transcribed_percentage is not None) and (untranscribed_percentage is None):
                            radius = 0.49
                            if (radius > 0):
                                # print('Plot circle at x=%d y=%d for %s %s' % (percentage_diff_index, row_signature_index, row_signature,percentage_string))
                                circle = plt.Circle((percentage_diff_index + 0.5, row_signature_index + 0.5), radius, color='royalblue', fill=True)
                                ax.add_artist(circle)
                        elif (untranscribed_percentage is not None) and (transcribed_percentage is None):
                            radius = 0.49
                            if (radius > 0):
                                # print('Plot circle at x=%d y=%d for %s %s' % (percentage_diff_index, row_signature_index, row_signature,percentage_string))
                                circle = plt.Circle((percentage_diff_index + 0.5, row_signature_index + 0.5), radius, color='yellowgreen', fill=True)
                                ax.add_artist(circle)
                        elif (transcribed_percentage is not None) and (untranscribed_percentage is not None):
                            radius_transcribed = 0.49
                            radius_untranscribed = 0.49
                            if (radius_transcribed>radius_untranscribed):
                                #First transcribed
                                circle = plt.Circle((percentage_diff_index + 0.5, row_signature_index + 0.5), radius_transcribed, color='royalblue', fill=True)
                                ax.add_artist(circle)
                                #Second untranscribed
                                circle = plt.Circle((percentage_diff_index + 0.5, row_signature_index + 0.5), radius_untranscribed, color='yellowgreen', fill=True)
                                ax.add_artist(circle)
                            else:
                                #First untranscribed
                                circle = plt.Circle((percentage_diff_index + 0.5, row_signature_index + 0.5), radius_untranscribed, color='yellowgreen', fill=True)
                                ax.add_artist(circle)
                                #Second transcribed
                                circle = plt.Circle((percentage_diff_index + 0.5, row_signature_index + 0.5), radius_transcribed, color='royalblue', fill=True)
                                ax.add_artist(circle)


                elif (strand_bias==GENIC_VERSUS_INTERGENIC):
                    if row_signature in type2strand2percentagedict:
                        genic_percentage=None
                        intergenic_percentage=None
                        if GENIC in type2strand2percentagedict[row_signature] and type2strand2percentagedict[row_signature][GENIC][percentage_string]==1:
                            genic_percentage = 100
                        if INTERGENIC in type2strand2percentagedict[row_signature] and type2strand2percentagedict[row_signature][INTERGENIC][percentage_string]==1:
                            intergenic_percentage = 100

                        if (genic_percentage is not None) and (intergenic_percentage is None):
                            radius = 0.49
                            if (radius > 0):
                                # print('Plot circle at x=%d y=%d for %s %s' % (percentage_diff_index, row_signature_index, row_signature,percentage_string))
                                circle = plt.Circle((percentage_diff_index + 0.5, row_signature_index + 0.5), radius, color='cyan', fill=True)
                                ax.add_artist(circle)
                        elif (intergenic_percentage is not None) and (genic_percentage is None):
                            radius = 0.49
                            if (radius > 0):
                                # print('Plot circle at x=%d y=%d for %s %s' % (percentage_diff_index, row_signature_index, row_signature,percentage_string))
                                circle = plt.Circle((percentage_diff_index + 0.5, row_signature_index + 0.5), radius, color='gray', fill=True)
                                ax.add_artist(circle)
                        elif (genic_percentage is not None) and (intergenic_percentage is not None):
                            radius_genic = 0.49
                            radius_intergenic = 0.49
                            if (radius_genic>radius_intergenic):
                                #First genic
                                circle = plt.Circle((percentage_diff_index + 0.5, row_signature_index + 0.5), radius_genic, color='cyan', fill=True)
                                ax.add_artist(circle)
                                #Second intergenic
                                circle = plt.Circle((percentage_diff_index + 0.5, row_signature_index + 0.5), radius_intergenic, color='gray', fill=True)
                                ax.add_artist(circle)
                            else:
                                #First untranscribed
                                circle = plt.Circle((percentage_diff_index + 0.5, row_signature_index + 0.5), radius_intergenic, color='gray', fill=True)
                                ax.add_artist(circle)
                                #Second transcribed
                                circle = plt.Circle((percentage_diff_index + 0.5, row_signature_index + 0.5), radius_genic, color='cyan', fill=True)
                                ax.add_artist(circle)
        ######################################################################################################################################


        ##################################################################################
        # CODE GOES HERE TO CENTER X-AXIS LABELS...
        ax.set_xlim([0,len(percentage_strings)])
        ax.set_xticklabels([])

        ax.tick_params(axis='x', which='minor', length=0, labelsize=20)

        #major ticks
        ax.set_xticks(np.arange(0, len(percentage_strings), 1))
        #minor ticks
        ax.set_xticks(np.arange(0, len(percentage_strings), 1)+0.5,minor=True)
        ax.set_xticklabels(percentage_strings,minor=True)

        #Jul 7, 2020
        if strand_bias==LAGGING_VERSUS_LEADING:
            fig.suptitle('Lagging versus Leading Strand Bias', fontsize=30)
        elif strand_bias==TRANSCRIBED_VERSUS_UNTRANSCRIBED:
            fig.suptitle('Transcribed versus Untranscribed Strand Bias', fontsize=30)
        elif strand_bias==GENIC_VERSUS_INTERGENIC:
            fig.suptitle('Genic versus Intergenic Strand Bias', fontsize=30)

        ax.xaxis.set_ticks_position('top')

        plt.tick_params(
            axis='x',  # changes apply to the x-axis
            which='major',  # both major and minor ticks are affected
            bottom=False,  # ticks along the bottom edge are off
            top=False)  # labels along the bottom edge are off
        ##################################################################################


        ##################################################################################
        # CODE GOES HERE TO CENTER Y-AXIS LABELS...
        ax.set_ylim([0,len(rows_signatures)])
        ax.set_yticklabels([])

        ax.tick_params(axis='y', which='minor', length=0, labelsize=30)

        #major ticks
        ax.set_yticks(np.arange(0, len(rows_signatures), 1))
        #minor ticks
        ax.set_yticks(np.arange(0, len(rows_signatures), 1)+0.5,minor=True)
        ax.set_yticklabels(rows_signatures, minor=True)  # fontsize

        plt.tick_params(
            axis='y',  # changes apply to the x-axis
            which='major',  # both major and minor ticks are affected
            left=False)  # labels along the bottom edge are off
        ##################################################################################

        ##################################################################################
        # Gridlines based on major ticks
        ax.grid(which='major', color='black')
        ##################################################################################

        ##################################################################################
        # create the directory if it does not exists
        filename = '%s_Signatures_%s_with_circles_%s.png' % (signature_type,strand_bias,str(SIGNIFICANCE_LEVEL).replace('.','_'))
        figFile = os.path.join(strandbias_figures_outputDir, CIRCLE_PLOTS, filename)
        fig.savefig(figFile)
        fig.tight_layout()

        plt.cla()
        plt.close(fig)
        ##################################################################################

############################################################################################################################



############################################################################################################################
#Plot Legend only
def plot_legend(strandbias_figures_outputDir):

    strand_biases=[TRANSCRIBED_VERSUS_UNTRANSCRIBED, GENIC_VERSUS_INTERGENIC, LAGGING_VERSUS_LEADING]

    for strandbias in strand_biases:
        ##################################################################################
        fig = plt.figure(figsize=(4,1), dpi=300)
        ax = plt.gca()
        plt.axis('off')
        ##################################################################################

        ##################################################################################
        if strandbias==TRANSCRIBED_VERSUS_UNTRANSCRIBED:
            legend_elements = [
                Line2D([0], [0], marker='o', color='white', label=TRANSCRIBED_STRAND, markerfacecolor='royalblue' ,markersize=20),
                Line2D([0], [0], marker='o', color='white', label=UNTRANSCRIBED_STRAND, markerfacecolor='yellowgreen',markersize=20)]
        elif strandbias == GENIC_VERSUS_INTERGENIC:
            legend_elements = [
                Line2D([0], [0], marker='o', color='white', label=GENIC, markerfacecolor='cyan',markersize=20),
                Line2D([0], [0], marker='o', color='white', label=INTERGENIC, markerfacecolor='gray',markersize=20)]
        elif (strandbias==LAGGING_VERSUS_LEADING):
            legend_elements = [
                Line2D([0], [0], marker='o', color='white', label=LAGGING, markerfacecolor='indianred', markersize=20),
                Line2D([0], [0], marker='o', color='white', label=LEADING, markerfacecolor='goldenrod', markersize=20)]

        ax.legend(handles=legend_elements, bbox_to_anchor=(0, 0.5), loc='center left' ,fontsize = 20)
        ##################################################################################

        ##################################################################################
        # create the directory if it does not exists
        filename = 'Legend_%s.png' % (strandbias)
        figFile = os.path.join(strandbias_figures_outputDir, CIRCLE_PLOTS, filename)
        fig.savefig(figFile)
        fig.tight_layout()

        plt.cla()
        plt.close(fig)
        ##################################################################################

############################################################################################################################


############################################################################################################################
#Sep 19, 2020
def plot_six_mutations_sbs_signatures_circle_figures(sbs_signatures,
            strand_bias,
            strandbias_figures_outputDir,
            significance_level,
            signature2mutation_type2strand2percentagedict,
            percentage_strings):

    mutation_types=six_mutation_types

    #####################################################################
    if strand_bias==LAGGING_VERSUS_LEADING:
        strands=replicationStrands
    elif strand_bias==TRANSCRIBED_VERSUS_UNTRANSCRIBED:
        strands=transcriptionStrands
    elif strand_bias==GENIC_VERSUS_INTERGENIC:
        strands=genicVersusIntergenicStrands
    #####################################################################

    #####################################################################
    rows_sbs_signatures=[]

    #Fill rows_sbs_signatures
    for signature in sbs_signatures:
        if signature in signature2mutation_type2strand2percentagedict:
            for mutation_type in signature2mutation_type2strand2percentagedict[signature]:
                for strand in strands:
                    if strand in signature2mutation_type2strand2percentagedict[signature][mutation_type]:
                        for percentage_string in percentage_strings:
                            if (percentage_string in signature2mutation_type2strand2percentagedict[signature][mutation_type][strand]) and (signature2mutation_type2strand2percentagedict[signature][mutation_type][strand][percentage_string]==1):
                                if signature not in rows_sbs_signatures:
                                    rows_sbs_signatures.append(signature)
    #####################################################################

    #####################################################################
    rows_sbs_signatures=sorted(rows_sbs_signatures,key=natural_key,reverse=True)
    #####################################################################

    #####################################################################
    xticklabels_list = percentage_strings * len(mutation_types)
    #####################################################################

    if (len(rows_sbs_signatures)>0):

        #####################################################################
        plot1, panel1 = plt.subplots(figsize=(5+1.5*len(xticklabels_list), 10+1.5*len(rows_sbs_signatures)))
        # plot1, panel1 = plt.subplots(figsize=(5+1.4*len(xticklabels_list), 10+len(rows_sbs_signatures))) Title and mutation texts are not seen.
        plt.rc('axes', edgecolor='lightgray')
        # panel1 = plt.axes([0.04, 0.09, 0.95, 0.75])

        #make aspect ratio square
        panel1.set_aspect(1.0)
        #####################################################################


        ##################################################################################
        #set title
        if strand_bias==LAGGING_VERSUS_LEADING:
            title='Lagging versus Leading Strand Bias'
        elif strand_bias==TRANSCRIBED_VERSUS_UNTRANSCRIBED:
            title='Transcribed versus Untranscribed Strand Bias'
        elif strand_bias==GENIC_VERSUS_INTERGENIC:
            title='Genic versus Intergenic Strand Bias'
        panel1.text(len(percentage_strings)*3, len(rows_sbs_signatures)+2.5, title,  horizontalalignment='center', fontsize=60, fontweight='bold', fontname='Arial')
        ##################################################################################

        ##################################################################################
        #Colors from SigProfilerPlotting tool to be consistent
        colors = [[3 / 256, 189 / 256, 239 / 256],
                  [1 / 256, 1 / 256, 1 / 256],
                  [228 / 256, 41 / 256, 38 / 256],
                  [203 / 256, 202 / 256, 202 / 256],
                  [162 / 256, 207 / 256, 99 / 256],
                  [236 / 256, 199 / 256, 197 / 256]]

        #Put rectangles
        x = 0

        for i in range(0, len(mutation_types), 1):
            panel1.text((x+(len(percentage_strings)/2)-0.75), len(rows_sbs_signatures)+1.5, mutation_types[i], fontsize=55, fontweight='bold', fontname='Arial')
            panel1.add_patch(plt.Rectangle((x+.0415, len(rows_sbs_signatures)+0.75), len(percentage_strings)-(2*.0415), .5, facecolor=colors[i], clip_on=False))
            panel1.add_patch(plt.Rectangle((x, 0), len(percentage_strings), len(rows_sbs_signatures), facecolor=colors[i], zorder=0, alpha=0.25,edgecolor='grey'))
            x += len(percentage_strings)
        ##################################################################################

        ##################################################################################
        # CODE GOES HERE TO CENTER X-AXIS LABELS...
        panel1.set_xlim([0,len(mutation_types)*len(percentage_strings)])
        panel1.set_xticklabels([])

        panel1.tick_params(axis='x', which='minor', length=0, labelsize=35)

        #major ticks
        panel1.set_xticks(np.arange(0, len(mutation_types)*len(percentage_strings), 1))
        #minor ticks
        panel1.set_xticks(np.arange(0, len(mutation_types)*len(percentage_strings), 1)+0.5,minor=True)

        panel1.set_xticklabels(xticklabels_list,minor=True)

        panel1.xaxis.set_label_position('top')
        panel1.xaxis.set_ticks_position('top')

        plt.tick_params(
            axis='x',  # changes apply to the x-axis
            which='major',  # both major and minor ticks are affected
            bottom=False,  # ticks along the bottom edge are off
            top=False)  # labels along the bottom edge are off
        ##################################################################################


        ##################################################################################
        # CODE GOES HERE TO CENTER Y-AXIS LABELS...
        panel1.set_ylim([0,len(rows_sbs_signatures)])
        panel1.set_yticklabels([])

        panel1.tick_params(axis='y', which='minor', length=0, labelsize=40)

        #major ticks
        panel1.set_yticks(np.arange(0, len(rows_sbs_signatures), 1))
        #minor ticks
        panel1.set_yticks(np.arange(0, len(rows_sbs_signatures), 1)+0.5,minor=True)

        panel1.set_yticklabels(rows_sbs_signatures, minor=True)  # fontsize

        plt.tick_params(
            axis='y',  # changes apply to the x-axis
            which='major',  # both major and minor ticks are affected
            left=False)  # labels along the bottom edge are off
        ##################################################################################

        ##################################################################################
        # Gridlines based on major ticks
        panel1.grid(which='major', color='black', zorder=3)
        ##################################################################################

        ##################################################################################
        #Put the legend
        if strand_bias==TRANSCRIBED_VERSUS_UNTRANSCRIBED:
            legend_elements = [
                Line2D([0], [0], marker='o', color='white', label=TRANSCRIBED_STRAND, markerfacecolor='royalblue' ,markersize=40),
                Line2D([0], [0], marker='o', color='white', label=UNTRANSCRIBED_STRAND, markerfacecolor='yellowgreen',markersize=40)]
        elif strand_bias == GENIC_VERSUS_INTERGENIC:
                legend_elements = [
                    Line2D([0], [0], marker='o', color='white', label=GENIC, markerfacecolor='cyan',markersize=40),
                    Line2D([0], [0], marker='o', color='white', label=INTERGENIC, markerfacecolor='gray',markersize=40)]
        elif (strand_bias==LAGGING_VERSUS_LEADING):
            legend_elements = [
                Line2D([0], [0], marker='o', color='white', label=LAGGING, markerfacecolor='indianred', markersize=40),
                Line2D([0], [0], marker='o', color='white', label=LEADING, markerfacecolor='goldenrod', markersize=40)]

        panel1.legend(handles=legend_elements,ncol=len(legend_elements), bbox_to_anchor=(1, -0.1),loc='upper right', fontsize=40)
        ##################################################################################


        ######################################################################################################################################
        for percentage_diff_index, percentage_string in enumerate(percentage_strings):
             for mutation_type_index, mutation_type in enumerate(mutation_types):
                for row_sbs_signature_index, row_sbs_signature in enumerate(rows_sbs_signatures):
                    if (strand_bias==LAGGING_VERSUS_LEADING):
                        if row_sbs_signature in signature2mutation_type2strand2percentagedict:
                            if mutation_type in signature2mutation_type2strand2percentagedict[row_sbs_signature]:
                                lagging_percentage = None
                                leading_percentage = None

                                if (LAGGING in signature2mutation_type2strand2percentagedict[row_sbs_signature][mutation_type]) and (signature2mutation_type2strand2percentagedict[row_sbs_signature][mutation_type][LAGGING][percentage_string]==1):
                                    lagging_percentage = 100
                                if (LEADING in signature2mutation_type2strand2percentagedict[row_sbs_signature][mutation_type]) and (signature2mutation_type2strand2percentagedict[row_sbs_signature][mutation_type][LEADING][percentage_string]==1):
                                    leading_percentage = 100

                                if (lagging_percentage is not None) and (leading_percentage is None):
                                    radius = 0.49
                                    if (radius > 0):
                                        # print('Plot circle at x=%d y=%d for %s %s %s' % (mutation_type_index * len(percentage_strings) + percentage_diff_index, row_sbs_signature_index, row_sbs_signature,mutation_type, percentage_string))
                                        panel1.add_patch(plt.Circle((mutation_type_index * len(percentage_strings) + percentage_diff_index + 0.5, row_sbs_signature_index + 0.5),radius, color='indianred', fill=True))
                                elif (leading_percentage is not None) and (lagging_percentage is None):
                                    radius = 0.49
                                    if (radius > 0):
                                        # print('Plot circle at x=%d y=%d for %s %s %s' % (mutation_type_index * len(percentage_strings) + percentage_diff_index, row_sbs_signature_index, row_sbs_signature,mutation_type, percentage_string))
                                        panel1.add_patch(plt.Circle((mutation_type_index * len(percentage_strings) + percentage_diff_index + 0.5, row_sbs_signature_index + 0.5),radius, color='goldenrod', fill=True))

                                elif (lagging_percentage is not None) and (leading_percentage is not None):
                                    radius_lagging = 0.49
                                    radius_leading = 0.49
                                    if (radius_lagging > radius_leading):
                                        # First lagging
                                        panel1.add_patch(plt.Circle((mutation_type_index * len(percentage_strings) + percentage_diff_index + 0.5, row_sbs_signature_index + 0.5),radius_lagging, color='indianred', fill=True))
                                        # Second leading
                                        panel1.add_patch(plt.Circle((mutation_type_index * len(percentage_strings) + percentage_diff_index + 0.5, row_sbs_signature_index + 0.5),radius_leading, color='goldenrod', fill=True))
                                    else:
                                        # First leading
                                        panel1.add_patch(plt.Circle((mutation_type_index * len(percentage_strings) + percentage_diff_index + 0.5, row_sbs_signature_index + 0.5),radius_leading, color='goldenrod', fill=True))
                                        # Second lagging
                                        panel1.add_patch(plt.Circle((mutation_type_index * len(percentage_strings) + percentage_diff_index + 0.5, row_sbs_signature_index + 0.5),radius_lagging, color='indianred', fill=True))

                    elif (strand_bias == GENIC_VERSUS_INTERGENIC):
                        if row_sbs_signature in signature2mutation_type2strand2percentagedict:
                            if mutation_type in signature2mutation_type2strand2percentagedict[row_sbs_signature]:
                                genic_percentage = None
                                intergenic_percentage = None

                                if (GENIC in signature2mutation_type2strand2percentagedict[row_sbs_signature][mutation_type]) and (signature2mutation_type2strand2percentagedict[row_sbs_signature][mutation_type][GENIC][percentage_string]==1):
                                    genic_percentage = 100
                                if (INTERGENIC in signature2mutation_type2strand2percentagedict[row_sbs_signature][mutation_type]) and (signature2mutation_type2strand2percentagedict[row_sbs_signature][mutation_type][INTERGENIC][percentage_string]==1):
                                    intergenic_percentage = 100

                                if (genic_percentage is not None) and (intergenic_percentage is None):
                                    radius = 0.49
                                    if (radius > 0):
                                        # print('Plot circle at x=%d y=%d for %s %s %s' % (mutation_type_index * len(percentage_strings) + percentage_diff_index, row_sbs_signature_index,row_sbs_signature, mutation_type, percentage_string))
                                        panel1.add_patch(plt.Circle((mutation_type_index * len(percentage_strings) + percentage_diff_index + 0.5,row_sbs_signature_index + 0.5), radius, color='cyan',fill=True))

                                elif (intergenic_percentage is not None) and (genic_percentage is None):
                                    radius = 0.49
                                    if (radius > 0):
                                        # print('Plot circle at x=%d y=%d for %s %s %s' % (mutation_type_index * len(percentage_strings) + percentage_diff_index, row_sbs_signature_index,row_sbs_signature, mutation_type, percentage_string))
                                        panel1.add_patch(plt.Circle((mutation_type_index * len(percentage_strings) + percentage_diff_index + 0.5,row_sbs_signature_index + 0.5), radius, color='gray',fill=True))

                                elif (genic_percentage is not None) and (intergenic_percentage is not None):
                                    radius_genic = 0.49
                                    radius_intergenic = 0.49
                                    if (radius_genic > radius_intergenic):
                                        # First genic
                                        panel1.add_patch(plt.Circle((mutation_type_index * len(percentage_strings) + percentage_diff_index + 0.5,row_sbs_signature_index + 0.5), radius_genic,color='cyan', fill=True))
                                        # Second intergenic
                                        panel1.add_patch(plt.Circle((mutation_type_index * len(percentage_strings) + percentage_diff_index + 0.5,row_sbs_signature_index + 0.5), radius_intergenic,color='gray', fill=True))

                                    else:
                                        # First intergenic
                                        panel1.add_patch(plt.Circle((mutation_type_index * len(percentage_strings) + percentage_diff_index + 0.5,row_sbs_signature_index + 0.5), radius_intergenic, color='gray', fill=True))
                                        # Second genic
                                        panel1.add_patch(plt.Circle((mutation_type_index * len(percentage_strings) + percentage_diff_index + 0.5,row_sbs_signature_index + 0.5), radius_genic,color='cyan', fill=True))


                    elif (strand_bias == TRANSCRIBED_VERSUS_UNTRANSCRIBED):
                        if row_sbs_signature in signature2mutation_type2strand2percentagedict:
                            if mutation_type in signature2mutation_type2strand2percentagedict[row_sbs_signature]:
                                transcribed_percentage = None
                                untranscribed_percentage = None

                                if (TRANSCRIBED_STRAND in signature2mutation_type2strand2percentagedict[row_sbs_signature][mutation_type]) and (signature2mutation_type2strand2percentagedict[row_sbs_signature][mutation_type][TRANSCRIBED_STRAND][percentage_string]==1):
                                    transcribed_percentage = 100
                                if (UNTRANSCRIBED_STRAND in signature2mutation_type2strand2percentagedict[row_sbs_signature][mutation_type]) and (signature2mutation_type2strand2percentagedict[row_sbs_signature][mutation_type][UNTRANSCRIBED_STRAND][percentage_string]==1):
                                    untranscribed_percentage = 100

                                if (transcribed_percentage is not None) and (untranscribed_percentage is None):
                                    radius = 0.49
                                    if (radius > 0):
                                        # print('Plot circle at x=%d y=%d for %s %s %s' % (mutation_type_index * len(percentage_strings) + percentage_diff_index, row_sbs_signature_index,row_sbs_signature, mutation_type, percentage_string))
                                        panel1.add_patch(plt.Circle((mutation_type_index * len(percentage_strings) + percentage_diff_index + 0.5,row_sbs_signature_index + 0.5), radius, color='royalblue',fill=True))

                                elif (untranscribed_percentage is not None) and (transcribed_percentage is None):
                                    radius = 0.49
                                    if (radius > 0):
                                        # print('Plot circle at x=%d y=%d for %s %s %s' % (mutation_type_index * len(percentage_strings) + percentage_diff_index, row_sbs_signature_index,row_sbs_signature, mutation_type, percentage_string))
                                        panel1.add_patch(plt.Circle((mutation_type_index * len(percentage_strings) + percentage_diff_index + 0.5,row_sbs_signature_index + 0.5), radius, color='yellowgreen',fill=True))

                                elif (transcribed_percentage is not None) and (untranscribed_percentage is not None):
                                    radius_transcribed = 0.49
                                    radius_untranscribed = 0.49
                                    if (radius_transcribed > radius_untranscribed):
                                        # First transcribed
                                        panel1.add_patch(plt.Circle((mutation_type_index * len(percentage_strings) + percentage_diff_index + 0.5,row_sbs_signature_index + 0.5), radius_transcribed,color='royalblue', fill=True))
                                        # Second untranscribed
                                        panel1.add_patch(plt.Circle((mutation_type_index * len(percentage_strings) + percentage_diff_index + 0.5,row_sbs_signature_index + 0.5), radius_untranscribed,color='yellowgreen', fill=True))

                                    else:
                                        # First untranscribed
                                        panel1.add_patch(plt.Circle((mutation_type_index * len(percentage_strings) + percentage_diff_index + 0.5,row_sbs_signature_index + 0.5), radius_untranscribed,color='yellowgreen', fill=True))
                                        # Second transcribed
                                        panel1.add_patch(plt.Circle((mutation_type_index * len(percentage_strings) + percentage_diff_index + 0.5,row_sbs_signature_index + 0.5), radius_transcribed,color='royalblue', fill=True))
        ######################################################################################################################################


        ##################################################################################
        # create the directory if it does not exists
        filename = 'SBS_Signatures_%s_with_circle_plot_%s.png' % (strand_bias,str(significance_level).replace('.','_'))
        figFile = os.path.join(strandbias_figures_outputDir,CIRCLE_PLOTS, filename)
        plot1.savefig(figFile,bbox_inches='tight')
        plot1.tight_layout()

        plt.cla()
        plt.close(plot1)
        ##################################################################################

############################################################################################################################
