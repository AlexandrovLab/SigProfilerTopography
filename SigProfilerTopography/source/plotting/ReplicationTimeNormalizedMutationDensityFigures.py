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

import os
import sys
import numpy as np
import shutil
import pandas as pd

import matplotlib
BACKEND = 'Agg'
if matplotlib.get_backend().lower() != BACKEND.lower():
    # If backend is not set properly a call to describe will hang
    matplotlib.use(BACKEND)

from matplotlib import pyplot as plt

from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
rcParams.update({'figure.max_open_warning': 0})

from SigProfilerTopography.source.commons.TopographyCommons import SBS_CONTEXTS
from SigProfilerTopography.source.commons.TopographyCommons import ID
from SigProfilerTopography.source.commons.TopographyCommons import DBS

from SigProfilerTopography.source.commons.TopographyCommons import DATA
from SigProfilerTopography.source.commons.TopographyCommons import FIGURE
from SigProfilerTopography.source.commons.TopographyCommons import REPLICATIONTIME

from SigProfilerTopography.source.commons.TopographyCommons import SAMPLES

from SigProfilerTopography.source.commons.TopographyCommons import AGGREGATEDSUBSTITUTIONS
from SigProfilerTopography.source.commons.TopographyCommons import AGGREGATEDINDELS
from SigProfilerTopography.source.commons.TopographyCommons import AGGREGATEDDINUCS
from SigProfilerTopography.source.commons.TopographyCommons import SIGNATUREBASED
from SigProfilerTopography.source.commons.TopographyCommons import INDELBASED
from SigProfilerTopography.source.commons.TopographyCommons import MICROHOMOLOGY
from SigProfilerTopography.source.commons.TopographyCommons import REPEAT

from SigProfilerTopography.source.commons.TopographyCommons import takeAverage
from SigProfilerTopography.source.commons.TopographyCommons import getDictionary

from SigProfilerTopography.source.commons.TopographyCommons import Table_SBS_Signature_Cutoff_NumberofMutations_AverageProbability_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Table_DBS_Signature_Cutoff_NumberofMutations_AverageProbability_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Table_ID_Signature_Cutoff_NumberofMutations_AverageProbability_Filename

from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofSubsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofIndelsDict
from SigProfilerTopography.source.commons.TopographyCommons import Sample2NumberofDinucsDictFilename

from SigProfilerTopography.source.commons.TopographyCommons import getSample2SubsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2IndelsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import Sample2DinucsSignature2NumberofMutationsDictFilename

from SigProfilerTopography.source.commons.TopographyCommons import PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_TOOL
from SigProfilerTopography.source.commons.TopographyCommons import PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_MANUSCRIPT


# Same plot function for aggregated mutations and signature based mutations
def plotNormalizedMutationDensityFigureWithSimulations(title,
                                                       ylabel,
                                                       normalizedMutationDensityList,
                                                       sample,
                                                       signature,
                                                       analysesType,
                                                       indelType,
                                                       barcolor,
                                                       fillcolor,
                                                       outputDir,
                                                       jobname,
                                                       sample2NumberofMutationsDict,
                                                       signature_df,
                                                       sample2Signature2NumberofMutationsDict,
                                                       numberofSimulations,
                                                       plot_mode):

    #################################################################################
    ############################# For Simulations starts ############################
    #################################################################################
    listofSimulations = None

    # read the simulations
    if (numberofSimulations > 0):
        if (analysesType==SIGNATUREBASED):
            # If analysesType is SIGNATUREBASED  originalTitle holds the signature
            listofSimulations = readNormalizedMutationDataForSimulations(sample,signature,outputDir,jobname,numberofSimulations)
        elif (analysesType==INDELBASED):
            listofSimulations = readNormalizedMutationDataForSimulations(sample,indelType,outputDir,jobname,numberofSimulations)
        else:
            listofSimulations = readNormalizedMutationDataForSimulations(sample,analysesType,outputDir,jobname,numberofSimulations)

    simulationsLows, simulationsMeans, simulationsHighs = takeAverage(listofSimulations)

    #################################################################################
    ############################# For Simulations ends ##############################
    #################################################################################

    os.makedirs(os.path.join(outputDir, jobname, FIGURE, REPLICATIONTIME), exist_ok=True)

    # Note if you decrease the figure size decrease the fontsize accordingly
    # if plot_mode == PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_MANUSCRIPT and (title == 'Aggregated Substitutions' or
    #                                                                     title == 'Aggregated Dinucs' or
    #                                                                     title == 'Aggregated Indels'):
        # fig = plt.figure(figsize=(17, 17), dpi=300) # wider and higher for aggregated mutations
    # else:
    #     fig = plt.figure(figsize=(15, 15), dpi=300)

    # if plot_mode == PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_MANUSCRIPT:
    # else:
    #     ax = plt.gca()

    # Let's use same for aggregated mutations and signature based mutations
    fwidth = 15
    fheight = 15

    fig = plt.figure(figsize=(fwidth, fheight), facecolor=None)
    plt.style.use('ggplot')

    # define margins -> size in inches / figure dimension
    left_margin = 0.95 / fwidth
    right_margin = 0.2 / fwidth
    bottom_margin = 0.5 / fheight
    top_margin = 0.25 / fheight

    # create axes
    # dimensions are calculated relative to the figure size
    x = left_margin  # horiz. position of bottom-left corner
    y = bottom_margin  # vert. position of bottom-left corner
    w = 1 - (left_margin + right_margin)  # width of axes
    h = 1 - (bottom_margin + top_margin)  # height of axes
    ax = fig.add_axes([x, y, w, h])

    plt.style.use('ggplot')

    # This code makes the background white.
    ax.set_facecolor('white')

    # Note x get labels w.r.t. the order given here, 0 means get the 0th label from  xticks
    x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    width = 0.9  # the width of the bars
    plt.bar(x, normalizedMutationDensityList, width, label='Real Somatic Mutations', color=barcolor, edgecolor="black", linewidth=3, zorder=1)

    # plt.xticks(np.arange(10),('1st', '2nd', '3rd', '4th', '5th','6th','7th','8th','9th','10th'),rotation=20)
    # also works
    # plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])

    if plot_mode == PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_TOOL:
        plt.title(title, fontsize=40, fontweight='bold')
        plt.xlabel('Early <--- Replication Time ---> Late', fontsize=40, fontweight='semibold')
        plt.ylabel(ylabel, fontsize=40, fontweight='semibold')

        # Set label locations.
        plt.yticks(np.arange(0, 1.01, step=0.2))
        # This code puts some extra space below 0 and above 1.0
        plt.ylim(-0.01, 1.01)

        plt.tick_params(axis='y', which='major', labelsize=40, width=3, length=10)
        plt.tick_params(axis='y', which='minor', labelsize=40, width=3, length=10)

        for edge_i in ['left']:
            ax.spines[edge_i].set_edgecolor("black")
            ax.spines[edge_i].set_linewidth(3)
            # This code draws line only between [0,1]
            # This is not needed
            # ax.spines[edge_i].set_bounds(0, 1)

        plt.tick_params(
            axis='x',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            bottom=False,  # ticks along the bottom edge are off
            top=False,  # ticks along the top edge are off
            labelbottom=False)  # labels along the bottom edge are off

        if simulationsMeans is not None:
            sims_dashed_line = plt.plot(x, simulationsMeans, 'o--', color='black', label='Simulated Somatic Mutations', linewidth=5, zorder=4)
            if (simulationsLows is not None) and (simulationsHighs is not None):
                # if (len(simulationsLows)==len(simulationsHighs)):
                plt.fill_between(x, np.array(simulationsLows), np.array(simulationsHighs), facecolor=fillcolor, zorder=2)


    elif plot_mode == PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_MANUSCRIPT:
        plt.title(title, fontsize=90, pad=80, loc='center') # fontweight='bold'

        if title == 'Aggregated Substitutions' or title ==  'Aggregated Dinucs' or title == 'Aggregated Indels':
            plt.xlabel('Early <-------------> Late\nReplication Time', fontsize=80) # fontweight='semibold'
            plt.ylabel('Normalized\nMutation Density', fontsize=80,  labelpad=15) # fontweight='semibold'

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(True)

        for edge_i in ['left']:
            ax.spines[edge_i].set_edgecolor("black")
            ax.spines[edge_i].set_linewidth(3)
            # This code draws line only between [0,1]
            ax.spines[edge_i].set_bounds(0, 1)

        # set axis ticks
        # ax.tick_params(axis='both', which='both', length=0)
        ax.tick_params(axis='x', which='both', length=0)
        plt.tick_params(axis='y', which='major', labelsize=80, width=3, length=10)
        plt.tick_params(axis='y', which='minor', labelsize=80, width=3, length=10)

        # Set label locations.
        plt.yticks(np.arange(0, 1.01, step=0.2))
        # This code puts some extra space below 0 and above 1.0
        plt.ylim(-0.01, 1.01)

        # set axis labels
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax.get_yticklabels(), visible=True)
        ax.spines["bottom"].set_color('black')
        ax.spines["left"].set_color('black')
        ax.spines["top"].set_color('black')
        ax.spines["right"].set_color('black')

        # to put the legend's upper right-hand corner at (0.8,1) optimized for SBS6 in BreastCancer560 data for SigProfilerTopography Overview Figure
        # legend = ax.legend((bars[0], sims_dashed_line[0]),('Real', 'Simulated'),prop={'size': 50}, loc='upper right', bbox_to_anchor = (0.8, 1))
        # to put the legend's upper left corner at (x,y) optimized for SBS2 in BreastCancer560 data for Replication Overview Figure
        # legend = ax.legend((bars[0], sims_dashed_line[0]),('Real', 'Simulated'),prop={'size': 43}, loc='upper left', bbox_to_anchor = (0.02, 0.9))

        # if (legend is not None):
        #    frame = legend.get_frame()
        #    frame.set_facecolor('white')
        #    frame.set_edgecolor('black')

        if simulationsMeans is not None:
            sims_dashed_line = plt.plot(x, simulationsMeans, 'o--', color='black', label='Simulated Somatic Mutations', linewidth=6, zorder=4)
            if (simulationsLows is not None) and (simulationsHighs is not None):
                # if (len(simulationsLows)==len(simulationsHighs)):
                plt.fill_between(x, np.array(simulationsLows), np.array(simulationsHighs), facecolor=fillcolor, zorder=2)

    if sample is None:
        if (analysesType == INDELBASED):
            figureName = '%s_replication_time.png' % (indelType.replace(' ', ''))
        elif (analysesType == SIGNATUREBASED):
            # discreet mode [cancer_type     signature       cutoff  number_of_mutations     average_probability
            # samples_list    len(samples_list)       len(all_samples_list)   percentage_of_samples]
            #
            # probability mode [cancer_type     signature       number_of_mutations     number_of_all_mutations
            # average_probability     samples_list    len(samples_list)       len(all_samples_list)   percentage_of_samples]
            signature_based_num_of_mutations = int(signature_df[signature_df['signature'] == signature]['number_of_mutations'].values[0])
            figureName = '%s_%d_replication_time.png' % (signature.replace(' ', ''), signature_based_num_of_mutations)
        else:
            # AGGREGATEDSUBSTITUTIONS
            # AGGREGATEDINDELS
            # AGGREGATEDDINUCS
            figureName = '%s_replication_time.png' % (analysesType)

    else:
        if (analysesType == INDELBASED):
            figureName = '%s_%s_replication_timey.png' %(indelType.replace(' ', ''),sample)
        elif (analysesType == SIGNATUREBASED):
            numberofMutations = sample2Signature2NumberofMutationsDict[sample][signature]
            figureName = '%s_%s_%d_replication_time.png' % (signature.replace(' ', ''), sample, numberofMutations)
        else:
            # AGGREGATEDSUBSTITUTIONS
            # AGGREGATEDINDELS
            # AGGREGATEDDINUCS
            numberofMutations = sample2NumberofMutationsDict[sample]
            figureName = '%s_%d_replication_time.png' % (sample, numberofMutations)

    if (sample is None):
        figureFile = os.path.join(outputDir, jobname, FIGURE, REPLICATIONTIME, figureName)
    else:
        os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, REPLICATIONTIME), exist_ok=True)
        figureFile = os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, REPLICATIONTIME, figureName)

    plt.savefig(figureFile, dpi=100, bbox_inches="tight") # test plt instead of fig
    fig.clear()
    plt.close(fig)


# analysisType can be aggregated subs, aggregated indels, MICROHOMOLOGY, INDEL, subsSignature, indelsSignature, dbsSignature
def readNormalizedMutationData(sample, analysisType, outputDir, jobname):
    if sample is None:
        filename = '%s_NormalizedMutationDensity.txt' %(analysisType)
        if ((analysisType == AGGREGATEDSUBSTITUTIONS) or (analysisType == AGGREGATEDINDELS) or (analysisType == AGGREGATEDDINUCS) or (analysisType == MICROHOMOLOGY) or  (analysisType == REPEAT)):
            filepath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, analysisType, filename)
        else:
            filepath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, SIGNATUREBASED, filename)

    else:
        filename = '%s_%s_NormalizedMutationDensity.txt' %(sample,analysisType)
        if ((analysisType==AGGREGATEDSUBSTITUTIONS) or (analysisType==AGGREGATEDINDELS) or (analysisType==AGGREGATEDDINUCS) or (analysisType==MICROHOMOLOGY) or (analysisType==REPEAT)):
            filepath = os.path.join(outputDir, jobname, DATA, SAMPLES, sample, REPLICATIONTIME, analysisType, filename)
        else:
            filepath = os.path.join(outputDir, jobname, DATA, SAMPLES, sample, REPLICATIONTIME, SIGNATUREBASED, filename)

    #Check if filepath exists
    if os.path.exists(filepath):
        # normalizedMutationData_df = pd.read_csv(filepath, sep=" ", comment='#', header=None)
        # normalizedMutationData_df.dropna(axis=1, how='any', inplace=True)
        # return normalizedMutationData_df
        normalizedMutationData_array = np.loadtxt(filepath, dtype=float)
        return normalizedMutationData_array
    else:
        return None


#########################################################
def readNormalizedMutationDataForSimulations(sample, indelorSignatureorAnalysesType, outputDir, jobname,numberofSimulations):
    listofAverages = []

    ######################################################
    for simNum in range(1, numberofSimulations + 1):
        if sample is None:
            filename = '%s_sim%d_NormalizedMutationDensity.txt' % (indelorSignatureorAnalysesType,simNum)
            if ((indelorSignatureorAnalysesType == AGGREGATEDSUBSTITUTIONS) or (indelorSignatureorAnalysesType == AGGREGATEDINDELS) or (indelorSignatureorAnalysesType == AGGREGATEDDINUCS) or (indelorSignatureorAnalysesType == MICROHOMOLOGY) or (indelorSignatureorAnalysesType == REPEAT)):
                filepath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME,indelorSignatureorAnalysesType, filename)
            else:
                filepath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME,SIGNATUREBASED, filename)

        else:
            filename = '%s_%s_sim%d_NormalizedMutationDensity.txt' % (sample, indelorSignatureorAnalysesType,simNum)
            if ((indelorSignatureorAnalysesType == AGGREGATEDSUBSTITUTIONS) or (indelorSignatureorAnalysesType == AGGREGATEDINDELS) or (indelorSignatureorAnalysesType == AGGREGATEDDINUCS) or (indelorSignatureorAnalysesType == MICROHOMOLOGY) or (indelorSignatureorAnalysesType == REPEAT)):
                filepath = os.path.join(outputDir,jobname, DATA, SAMPLES, sample, REPLICATIONTIME, indelorSignatureorAnalysesType,filename)
            else:
                filepath = os.path.join(outputDir,jobname, DATA, SAMPLES, sample, REPLICATIONTIME, SIGNATUREBASED,filename)


        # Check if filepath exists
        if os.path.exists(filepath):
            normalizedMutationData = np.loadtxt(filepath, dtype=float)
            listofAverages.append(normalizedMutationData)
    ######################################################

    return listofAverages
#########################################################


def plotSignatureFigures(color,
                         fillcolor,
                         analysesType,
                         outputDir,
                         jobname,
                         numberofSimulations,
                         sample_based,
                         sample2NumberofMutationsDict,
                         signature_df,
                         sample2Signature2NumberofMutationsDict,
                         plot_mode):

    for signature in signature_df['signature'].unique():
        # We check such file exists or not
        normalizedMutationData = readNormalizedMutationData(None, signature, outputDir, jobname)

        if (normalizedMutationData is not None) and (not np.all(normalizedMutationData == 0)):
            # if not all([v == 0.0 for v in normalizedMutationData]):
            # use all generator for all true check
            plotNormalizedMutationDensityFigureWithSimulations(signature,
                                                               'Normalized\nsingle base substitution density',
                                                               normalizedMutationData,
                                                               None,
                                                               signature,
                                                               analysesType,
                                                               None,
                                                               color,
                                                               fillcolor,
                                                               outputDir,
                                                               jobname,
                                                               sample2NumberofMutationsDict,
                                                               signature_df,
                                                               sample2Signature2NumberofMutationsDict,
                                                               numberofSimulations,
                                                               plot_mode)

    if sample_based:
        for sample in sample2Signature2NumberofMutationsDict:
            for signature in sample2Signature2NumberofMutationsDict[sample]:
                normalizedMutationData = readNormalizedMutationData(sample, signature, outputDir, jobname)

                if (normalizedMutationData is not None) and (not np.all(normalizedMutationData == 0)):
                    # if not all([v == 0.0 for v in normalizedMutationData]):
                    # use all generator for all true check
                    plotNormalizedMutationDensityFigureWithSimulations('%s_%s' % (signature, sample),
                                                                       'Normalized\nsingle base substitution density',
                                                                       normalizedMutationData,
                                                                       sample,
                                                                       signature,
                                                                       analysesType,
                                                                       None,
                                                                       color,
                                                                       fillcolor,
                                                                       outputDir,
                                                                       jobname,
                                                                       sample2NumberofMutationsDict,
                                                                       signature_df,
                                                                       sample2Signature2NumberofMutationsDict,
                                                                       numberofSimulations,
                                                                       plot_mode)

def plotAllMutationTypesFigures(title,
                                color,
                                fillcolor,
                                analysesType,
                                indelType,
                                outputDir,
                                jobname,
                                numberofSimulations,
                                sample_based,
                                sample2NumberofMutationsDict,
                                signature_cutoff_numberofmutations_averageprobability_df,
                                sample2Signature2NumberofMutationsDict,
                                plot_mode):

    if (analysesType == INDELBASED):
        normalizedMutationData = readNormalizedMutationData(None, indelType, outputDir, jobname)
    else:
        normalizedMutationData = readNormalizedMutationData(None, analysesType, outputDir, jobname)

    if (normalizedMutationData is not None) and (not np.all(normalizedMutationData == 0)):
        plotNormalizedMutationDensityFigureWithSimulations(title,
                                                           '\nNormalized mutation density',
                                                           normalizedMutationData,
                                                           None,
                                                           None,
                                                           analysesType,
                                                           indelType,
                                                           color,
                                                           fillcolor,
                                                           outputDir,
                                                           jobname,
                                                           sample2NumberofMutationsDict,
                                                           signature_cutoff_numberofmutations_averageprobability_df,
                                                           sample2Signature2NumberofMutationsDict,
                                                           numberofSimulations,
                                                           plot_mode)

    if sample_based:
        for sample in sample2NumberofMutationsDict:
            if (analysesType == INDELBASED):
                normalizedMutationData = readNormalizedMutationData(sample, indelType, outputDir, jobname)
            else:
                normalizedMutationData = readNormalizedMutationData(sample, analysesType, outputDir, jobname)

            if (normalizedMutationData is not None) and (not np.all(normalizedMutationData == 0)):
                plotNormalizedMutationDensityFigureWithSimulations('%s %s' % (sample,title),
                                                                   '\nNormalized mutation density',
                                                                   normalizedMutationData,
                                                                   sample,
                                                                   None,
                                                                   analysesType,
                                                                   indelType,
                                                                   color,
                                                                   fillcolor,
                                                                   outputDir,
                                                                   jobname,
                                                                   sample2NumberofMutationsDict,
                                                                   signature_cutoff_numberofmutations_averageprobability_df,
                                                                   sample2Signature2NumberofMutationsDict,
                                                                   numberofSimulations,
                                                                   plot_mode)


# Plot legend only
def plot_replication_time_legend(output_dir, jobname):
    fig, ax = plt.subplots(figsize=(6, 3))

    # This code makes the background white.
    ax.set_facecolor('white')

    plt.gca().set_axis_off()

    real_subs_rectangle = mpatches.Patch(label='Real Subs', edgecolor='black', facecolor='royalblue', lw=3)
    real_dinucs_rectangle = mpatches.Patch(label='Real Dinucs', edgecolor='black', facecolor='crimson', lw=3)
    real_indels_rectangle = mpatches.Patch(label='Real Indels', edgecolor='black', facecolor='yellowgreen', lw=3)

    legend_elements = [
        real_subs_rectangle,
        real_dinucs_rectangle,
        real_indels_rectangle,
        Line2D([0], [2], linestyle="--", marker='.', lw=5, color='black', label='Simulations', markerfacecolor='black', markersize=30)]

    plt.legend(facecolor='white',handles=legend_elements, handlelength=5, ncol=1, loc="lower center", fontsize=30)  # bbox_to_anchor=(1, 0.5),
    plt.gca().set_axis_off()

    filename = 'Replication_Time_Legend.png'
    filepath = os.path.join(output_dir, jobname, FIGURE, REPLICATIONTIME, filename)
    fig.savefig(filepath, dpi=100, bbox_inches="tight")

    plt.cla()
    plt.close(fig)



def replicationTimeNormalizedMutationDensityFigures(outputDir,
                                                    jobname,
                                                    numberofSimulations,
                                                    sample_based,
                                                    plot_mode):

    # Initialize these dataframes as empty dataframe
    # We will read these dataframes if there is the corresponding data
    subsSignature_cutoff_numberofmutations_averageprobability_df = pd.DataFrame()
    dinucsSignature_cutoff_numberofmutations_averageprobability_df = pd.DataFrame()
    indelsSignature_cutoff_numberofmutations_averageprobability_df = pd.DataFrame()

    subsSignature_cutoff_numberofmutations_averageprobability_path = os.path.join(outputDir, jobname, DATA, Table_SBS_Signature_Cutoff_NumberofMutations_AverageProbability_Filename)
    dinucsSignature_cutoff_numberofmutations_averageprobability_path = os.path.join(outputDir, jobname, DATA, Table_DBS_Signature_Cutoff_NumberofMutations_AverageProbability_Filename)
    indelsSignature_cutoff_numberofmutations_averageprobability_path = os.path.join(outputDir, jobname, DATA, Table_ID_Signature_Cutoff_NumberofMutations_AverageProbability_Filename)

    if os.path.exists(subsSignature_cutoff_numberofmutations_averageprobability_path):
        subsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(subsSignature_cutoff_numberofmutations_averageprobability_path, sep='\t', header=0,dtype={'cutoff': np.float32, 'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})

    if os.path.exists(dinucsSignature_cutoff_numberofmutations_averageprobability_path):
        dinucsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(dinucsSignature_cutoff_numberofmutations_averageprobability_path, sep='\t', header=0,dtype={'cutoff': np.float32, 'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})

    if os.path.exists(indelsSignature_cutoff_numberofmutations_averageprobability_path):
        indelsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(indelsSignature_cutoff_numberofmutations_averageprobability_path,sep='\t', header=0,dtype={'cutoff': np.float32, 'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})

    if sample_based:
        sample2NumberofSubsDict = getSample2NumberofSubsDict(outputDir, jobname)
        sample2NumberofIndelsDict = getSample2NumberofIndelsDict(outputDir, jobname)
        sample2NumberofDinucsDict = getDictionary(outputDir, jobname, Sample2NumberofDinucsDictFilename)

        sample2SubsSignature2NumberofMutationsDict = getSample2SubsSignature2NumberofMutationsDict(outputDir,jobname)
        sample2IndelsSignature2NumberofMutationsDict = getSample2IndelsSignature2NumberofMutationsDict(outputDir,jobname)
        sample2DinucsSignature2NumberofMutationsDict = getDictionary(outputDir,jobname,Sample2DinucsSignature2NumberofMutationsDictFilename)

    else:
        sample2NumberofSubsDict = {}
        sample2NumberofIndelsDict = {}
        sample2NumberofDinucsDict = {}

        sample2SubsSignature2NumberofMutationsDict = {}
        sample2IndelsSignature2NumberofMutationsDict = {}
        sample2DinucsSignature2NumberofMutationsDict = {}


    ##########################################################################################
    ##########################  Plot figures starts  #########################################
    ##########################################################################################
    # Plot aggregated mutations figures
    plotAllMutationTypesFigures('Aggregated Substitutions', 'royalblue', 'lightblue', AGGREGATEDSUBSTITUTIONS, None,
                                outputDir, jobname, numberofSimulations, sample_based, sample2NumberofSubsDict,
                                subsSignature_cutoff_numberofmutations_averageprobability_df,
                                sample2SubsSignature2NumberofMutationsDict, plot_mode)

    plotAllMutationTypesFigures('Aggregated Dinucs', 'crimson', 'lightpink', AGGREGATEDDINUCS, None, outputDir, jobname,
                                numberofSimulations, sample_based, sample2NumberofDinucsDict,
                                dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                sample2DinucsSignature2NumberofMutationsDict, plot_mode)

    plotAllMutationTypesFigures('Aggregated Indels', 'yellowgreen', 'lightgreen', AGGREGATEDINDELS, None, outputDir,
                                jobname, numberofSimulations, sample_based, sample2NumberofIndelsDict,
                                indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                sample2IndelsSignature2NumberofMutationsDict, plot_mode)

    plotAllMutationTypesFigures(MICROHOMOLOGY, 'yellowgreen', 'lightgreen', INDELBASED, MICROHOMOLOGY, outputDir,
                                jobname, numberofSimulations, sample_based, sample2NumberofIndelsDict,
                                indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                sample2IndelsSignature2NumberofMutationsDict, plot_mode)

    plotAllMutationTypesFigures(REPEAT, 'yellowgreen', 'lightgreen', INDELBASED, REPEAT, outputDir, jobname,
                                numberofSimulations, sample_based, sample2NumberofIndelsDict,
                                indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                sample2IndelsSignature2NumberofMutationsDict, plot_mode)

    if plot_mode == PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_MANUSCRIPT :
        plot_replication_time_legend(outputDir, jobname)

    if (not subsSignature_cutoff_numberofmutations_averageprobability_df.empty):
        plotSignatureFigures('royalblue',
                             'lightblue',
                             SIGNATUREBASED,
                             outputDir,
                             jobname,
                             numberofSimulations,
                             sample_based,
                             sample2NumberofSubsDict,
                             subsSignature_cutoff_numberofmutations_averageprobability_df,
                             sample2SubsSignature2NumberofMutationsDict,
                             plot_mode)
    if (not dinucsSignature_cutoff_numberofmutations_averageprobability_df.empty):
        plotSignatureFigures('crimson',
                             'lightpink',
                             SIGNATUREBASED,
                             outputDir,
                             jobname,
                             numberofSimulations,
                             sample_based,
                             sample2NumberofDinucsDict,
                             dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                             sample2DinucsSignature2NumberofMutationsDict,
                             plot_mode)
    if (not indelsSignature_cutoff_numberofmutations_averageprobability_df.empty):
        plotSignatureFigures('yellowgreen',
                             'lightgreen',
                             SIGNATUREBASED,
                             outputDir,
                             jobname,
                             numberofSimulations,
                             sample_based,
                             sample2NumberofIndelsDict,
                             indelsSignature_cutoff_numberofmutations_averageprobability_df,
                             sample2IndelsSignature2NumberofMutationsDict,
                             plot_mode)
    ##########################################################################################
    ##########################  Plot figures ends  ###########################################
    ##########################################################################################