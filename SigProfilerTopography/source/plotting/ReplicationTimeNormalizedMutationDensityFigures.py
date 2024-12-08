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
import io
import pickle
import copy
import sys
import numpy as np
import shutil
import pandas as pd
import multiprocessing

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

from PIL import Image
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4

from sklearn.preprocessing import LabelEncoder

from matplotlib.backends.backend_pdf import PdfPages

from functools import reduce

from collections import OrderedDict

import sigProfilerPlotting as spplt

#from SigProfilerTopography.source.commons.TopographyCommons import SBS_CONTEXTS
from SigProfilerTopography.source.commons.TopographyCommons import SBS
from SigProfilerTopography.source.commons.TopographyCommons import ID
from SigProfilerTopography.source.commons.TopographyCommons import DBS

from SigProfilerTopography.source.commons.TopographyCommons import DATA
from SigProfilerTopography.source.commons.TopographyCommons import FIGURE
from SigProfilerTopography.source.commons.TopographyCommons import REPLICATIONTIME
from SigProfilerTopography.source.commons.TopographyCommons import MULTILEVEL

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

from SigProfilerTopography.source.commons.TopographyCommons import fill_signature_mutationtype_strand_replicationtime_arrays

from SigProfilerTopography.source.commons.TopographyCommons import Table_SBS_Signature_Cutoff_NumberofMutations_AverageProbability_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Table_DBS_Signature_Cutoff_NumberofMutations_AverageProbability_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Table_ID_Signature_Cutoff_NumberofMutations_AverageProbability_Filename

from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIPTIONSTRAND
from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIPTIONSTRANDBIAS
from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIBED
from SigProfilerTopography.source.commons.TopographyCommons import UNTRANSCRIBED
from SigProfilerTopography.source.commons.TopographyCommons import NONTRANSCRIBED
from SigProfilerTopography.source.commons.TopographyCommons import QUESTIONABLE

from SigProfilerTopography.source.commons.TopographyCommons import REPLICATIONSTRAND
from SigProfilerTopography.source.commons.TopographyCommons import REPLICATIONSTRANDBIAS
from SigProfilerTopography.source.commons.TopographyCommons import LAGGING
from SigProfilerTopography.source.commons.TopographyCommons import LEADING

from SigProfilerTopography.source.commons.TopographyCommons import GENICINTERGENICBIAS
from SigProfilerTopography.source.commons.TopographyCommons import GENIC
from SigProfilerTopography.source.commons.TopographyCommons import INTERGENIC

from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofSubsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofIndelsDict
from SigProfilerTopography.source.commons.TopographyCommons import Sample2NumberofDinucsDictFilename

from SigProfilerTopography.source.commons.TopographyCommons import getSample2SubsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2IndelsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import Sample2DinucsSignature2NumberofMutationsDictFilename

from SigProfilerTopography.source.commons.TopographyCommons import PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_TOOL
from SigProfilerTopography.source.commons.TopographyCommons import PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_MANUSCRIPT

ordered_SBS96_substitution_types = [
    "A[C>A]A",
    "A[C>A]C",
    "A[C>A]G",
    "A[C>A]T",
    "C[C>A]A",
    "C[C>A]C",
    "C[C>A]G",
    "C[C>A]T",
    "G[C>A]A",
    "G[C>A]C",
    "G[C>A]G",
    "G[C>A]T",
    "T[C>A]A",
    "T[C>A]C",
    "T[C>A]G",
    "T[C>A]T",
    "A[C>G]A",
    "A[C>G]C",
    "A[C>G]G",
    "A[C>G]T",
    "C[C>G]A",
    "C[C>G]C",
    "C[C>G]G",
    "C[C>G]T",
    "G[C>G]A",
    "G[C>G]C",
    "G[C>G]G",
    "G[C>G]T",
    "T[C>G]A",
    "T[C>G]C",
    "T[C>G]G",
    "T[C>G]T",
    "A[C>T]A",
    "A[C>T]C",
    "A[C>T]G",
    "A[C>T]T",
    "C[C>T]A",
    "C[C>T]C",
    "C[C>T]G",
    "C[C>T]T",
    "G[C>T]A",
    "G[C>T]C",
    "G[C>T]G",
    "G[C>T]T",
    "T[C>T]A",
    "T[C>T]C",
    "T[C>T]G",
    "T[C>T]T",
    "A[T>A]A",
    "A[T>A]C",
    "A[T>A]G",
    "A[T>A]T",
    "C[T>A]A",
    "C[T>A]C",
    "C[T>A]G",
    "C[T>A]T",
    "G[T>A]A",
    "G[T>A]C",
    "G[T>A]G",
    "G[T>A]T",
    "T[T>A]A",
    "T[T>A]C",
    "T[T>A]G",
    "T[T>A]T",
    "A[T>C]A",
    "A[T>C]C",
    "A[T>C]G",
    "A[T>C]T",
    "C[T>C]A",
    "C[T>C]C",
    "C[T>C]G",
    "C[T>C]T",
    "G[T>C]A",
    "G[T>C]C",
    "G[T>C]G",
    "G[T>C]T",
    "T[T>C]A",
    "T[T>C]C",
    "T[T>C]G",
    "T[T>C]T",
    "A[T>G]A",
    "A[T>G]C",
    "A[T>G]G",
    "A[T>G]T",
    "C[T>G]A",
    "C[T>G]C",
    "C[T>G]G",
    "C[T>G]T",
    "G[T>G]A",
    "G[T>G]C",
    "G[T>G]G",
    "G[T>G]T",
    "T[T>G]A",
    "T[T>G]C",
    "T[T>G]G",
    "T[T>G]T"
]

ordered_DBS78_doublet_types = ['AC>CA', 'AC>CG', 'AC>CT', 'AC>GA', 'AC>GG', 'AC>GT', 'AC>TA', 'AC>TG',
                               'AC>TT', 'AT>CA', 'AT>CC', 'AT>CG', 'AT>GA', 'AT>GC', 'AT>TA', 'CC>AA',
                               'CC>AG', 'CC>AT', 'CC>GA', 'CC>GG', 'CC>GT', 'CC>TA', 'CC>TG', 'CC>TT',
                               'CG>AT', 'CG>GC', 'CG>GT', 'CG>TA', 'CG>TC', 'CG>TT', 'CT>AA', 'CT>AC',
                               'CT>AG', 'CT>GA', 'CT>GC', 'CT>GG', 'CT>TA', 'CT>TC', 'CT>TG', 'GC>AA',
                               'GC>AG', 'GC>AT', 'GC>CA', 'GC>CG', 'GC>TA', 'TA>AT', 'TA>CG', 'TA>CT',
                               'TA>GC', 'TA>GG', 'TA>GT', 'TC>AA', 'TC>AG', 'TC>AT', 'TC>CA', 'TC>CG',
                               'TC>CT', 'TC>GA', 'TC>GG', 'TC>GT', 'TG>AA', 'TG>AC', 'TG>AT', 'TG>CA',
                               'TG>CC', 'TG>CT', 'TG>GA', 'TG>GC', 'TG>GT', 'TT>AA', 'TT>AC', 'TT>AG',
                               'TT>CA', 'TT>CC', 'TT>CG', 'TT>GA', 'TT>GC', 'TT>GG']

ordered_ID83_indel_types = [
    "1:Del:C:0",
    "1:Del:C:1",
    "1:Del:C:2",
    "1:Del:C:3",
    "1:Del:C:4",
    "1:Del:C:5",
    "1:Del:T:0",
    "1:Del:T:1",
    "1:Del:T:2",
    "1:Del:T:3",
    "1:Del:T:4",
    "1:Del:T:5",
    "1:Ins:C:0",
    "1:Ins:C:1",
    "1:Ins:C:2",
    "1:Ins:C:3",
    "1:Ins:C:4",
    "1:Ins:C:5",
    "1:Ins:T:0",
    "1:Ins:T:1",
    "1:Ins:T:2",
    "1:Ins:T:3",
    "1:Ins:T:4",
    "1:Ins:T:5",
    # >1bp INDELS
    "2:Del:R:0",
    "2:Del:R:1",
    "2:Del:R:2",
    "2:Del:R:3",
    "2:Del:R:4",
    "2:Del:R:5",
    "3:Del:R:0",
    "3:Del:R:1",
    "3:Del:R:2",
    "3:Del:R:3",
    "3:Del:R:4",
    "3:Del:R:5",
    "4:Del:R:0",
    "4:Del:R:1",
    "4:Del:R:2",
    "4:Del:R:3",
    "4:Del:R:4",
    "4:Del:R:5",
    "5:Del:R:0",
    "5:Del:R:1",
    "5:Del:R:2",
    "5:Del:R:3",
    "5:Del:R:4",
    "5:Del:R:5",
    "2:Ins:R:0",
    "2:Ins:R:1",
    "2:Ins:R:2",
    "2:Ins:R:3",
    "2:Ins:R:4",
    "2:Ins:R:5",
    "3:Ins:R:0",
    "3:Ins:R:1",
    "3:Ins:R:2",
    "3:Ins:R:3",
    "3:Ins:R:4",
    "3:Ins:R:5",
    "4:Ins:R:0",
    "4:Ins:R:1",
    "4:Ins:R:2",
    "4:Ins:R:3",
    "4:Ins:R:4",
    "4:Ins:R:5",
    "5:Ins:R:0",
    "5:Ins:R:1",
    "5:Ins:R:2",
    "5:Ins:R:3",
    "5:Ins:R:4",
    "5:Ins:R:5",
    # MicroHomology INDELS
    "2:Del:M:1",
    "3:Del:M:1",
    "3:Del:M:2",
    "4:Del:M:1",
    "4:Del:M:2",
    "4:Del:M:3",
    "5:Del:M:1",
    "5:Del:M:2",
    "5:Del:M:3",
    "5:Del:M:4",
    "5:Del:M:5",
]


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
        if (analysesType == SIGNATUREBASED):
            # If analysesType is SIGNATUREBASED  originalTitle holds the signature
            listofSimulations = readNormalizedMutationDataForSimulations(sample,signature,outputDir,jobname,numberofSimulations)
        elif (analysesType == INDELBASED):
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
        if ((analysisType == AGGREGATEDSUBSTITUTIONS) or (analysisType == AGGREGATEDINDELS) or (analysisType == AGGREGATEDDINUCS) or (analysisType == MICROHOMOLOGY) or (analysisType == REPEAT)):
            filepath = os.path.join(outputDir, jobname, DATA, SAMPLES, sample, REPLICATIONTIME, analysisType, filename)
        else:
            filepath = os.path.join(outputDir, jobname, DATA, SAMPLES, sample, REPLICATIONTIME, SIGNATUREBASED, filename)

    # Check if filepath exists
    if os.path.exists(filepath):
        # normalizedMutationData_df = pd.read_csv(filepath, sep=" ", comment='#', header=None)
        # normalizedMutationData_df.dropna(axis=1, how='any', inplace=True)
        # return normalizedMutationData_df
        normalizedMutationData_array = np.loadtxt(filepath, dtype=float)
        return normalizedMutationData_array
    else:
        return None


def readNormalizedMutationDataForSimulations(sample, indelorSignatureorAnalysesType, outputDir, jobname,numberofSimulations):
    listofAverages = []

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

    return listofAverages


def plotSignatureFigures(color,
                         fillcolor,
                         ylabel,
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
                                                               ylabel, # 'Normalized\nsingle base substitution density',
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
                                                                       ylabel, #'Normalized\nsingle base substitution density',
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


def plot_normalized_mutation_density_figure_with_simulations(title,
                                                             normalized_mutation_density,
                                                             simulationsLows,
                                                             simulationsMeans,
                                                             simulationsHighs,
                                                             color,
                                                             fillcolor,
                                                             ax):

    plt.style.use('ggplot')

    # This code makes the background white.
    ax.set_facecolor('white')

    # Note x get labels w.r.t. the order given here, 0 means get the 0th label from  xticks
    x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    width = 0.9  # the width of the bars
    ax.bar(x, normalized_mutation_density, width, label='Real Somatic Mutations', color=color, edgecolor="black", linewidth=3, zorder=1)

    ax.set_title(title, fontsize=49, fontweight='bold') #40
    ax.set_xlabel('Early <--- Replication Time ---> Late', fontsize=49, fontweight='semibold') #40
    ax.set_ylabel('Normalized\nMutation Density', fontsize=49, fontweight='semibold') #40

    # Set label locations.
    ax.set_yticks(np.arange(0, 1.01, step=0.2))
    # This code puts some extra space below 0 and above 1.0
    ax.set_ylim(-0.01, 1.01)

    ax.tick_params(axis='y', which='major', labelsize=40, width=3, length=10)
    ax.tick_params(axis='y', which='minor', labelsize=40, width=3, length=10)

    for edge_i in ['left']:
        ax.spines[edge_i].set_edgecolor("black")
        ax.spines[edge_i].set_linewidth(3)
        # This code draws line only between [0,1]
        # This is not needed
        # ax.spines[edge_i].set_bounds(0, 1)

    # Hide other spines
    for edge_i in ['top', 'right', 'bottom']:
        ax.spines[edge_i].set_visible(False)

    ax.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)  # labels along the bottom edge are off

    if simulationsMeans is not None:
        sims_dashed_line = plt.plot(x, simulationsMeans, 'o--', color='black', label='Simulated Somatic Mutations',
                                    linewidth=5, zorder=4)
        if (simulationsLows is not None) and (simulationsHighs is not None):
            # if (len(simulationsLows)==len(simulationsHighs)):
            plt.fill_between(x, np.array(simulationsLows), np.array(simulationsHighs), facecolor=fillcolor, zorder=2)



def get_normalized_mutation_density(arr1, arr2):
    # ValueError: operands could not be broadcast together with shapes (0,) (10,)
    if (arr1.size == 0) or (arr2.size == 0):
        return np.zeros(arr1.shape)
    arr = arr1/arr2
    return arr/max(arr)

def get_normalized_mutation_density_2D(arr1, arr2):
    # ValueError: operands could not be broadcast together with shapes (0,) (10,)
    if (arr1.size == 0) or (arr2.size == 0):
        return np.zeros(arr1.shape)

    arr = arr1/arr2
    return  arr / np.max(arr, axis=1)[:, None]


def get_combined3_df(mutation_main_type, mutation_sub_type, column_of_interest, combined2_df):
    if mutation_main_type == SBS:
        combined3_df = \
        combined2_df[combined2_df['Mutation'] == mutation_sub_type].groupby([column_of_interest, 'ReplicationTime'])[
            'MutationCount'].sum().reset_index()

    elif mutation_main_type == DBS:
        # mutation_sub_types = ['AC>NN', 'AT>NN', 'CC>NN', 'CG>NN', 'CT>NN', 'GC>NN', 'TA>NN', 'TC>NN', 'TG>NN', 'TT>NN']
        # MutationLong
        # N:C[TC>AA]C
        combined3_df = \
        combined2_df[combined2_df['Mutation'].str.contains(mutation_sub_type[0:3])].groupby([column_of_interest, 'ReplicationTime'])[
            'MutationCount'].sum().reset_index()

    elif mutation_main_type == ID:
        # mutation_sub_types = ['1bp Deletion', '1bp Insertion', '>1bp Deletion', '>1bp Insertion', 'Microhomology']
        # MutationLong
        # N:1:Del:T:5
        # T:1:Ins:T:5
        if mutation_sub_type == '1bp Deletion':
            combined3_df = \
            combined2_df[combined2_df['Mutation'].str.contains('1:Del')].groupby([column_of_interest, 'ReplicationTime'])[
                'MutationCount'].sum().reset_index()
        elif mutation_sub_type == '1bp Insertion':
            combined3_df = \
            combined2_df[combined2_df['Mutation'].str.contains('1:Ins')].groupby([column_of_interest, 'ReplicationTime'])[
                'MutationCount'].sum().reset_index()
        elif mutation_sub_type == '>1bp Deletion':
            combined3_df = \
            combined2_df[combined2_df['Mutation'].str.contains(':Del:R:')].groupby([column_of_interest, 'ReplicationTime'])[
                'MutationCount'].sum().reset_index()
        elif mutation_sub_type == '>1bp Insertion':
            combined3_df = \
            combined2_df[combined2_df['Mutation'].str.contains(':Ins:R:')].groupby([column_of_interest, 'ReplicationTime'])[
                'MutationCount'].sum().reset_index()
        elif mutation_sub_type == 'Microhomology':
            combined3_df = \
            combined2_df[combined2_df['Mutation'].str.contains(':Del:M:')].groupby([column_of_interest, 'ReplicationTime'])[
                'MutationCount'].sum().reset_index()

    return combined3_df

def get_specific_mutations_df(mutation_main_type, mutation_sub_type, combined2_df):

    if mutation_main_type == SBS:
        specific_mutations_df = combined2_df[combined2_df['Mutation'] == mutation_sub_type].groupby('ReplicationTime')[
            'MutationCount'].sum().reset_index()
    elif mutation_main_type == DBS:
        # mutation_sub_types = ['AC>NN', 'AT>NN', 'CC>NN', 'CG>NN', 'CT>NN', 'GC>NN', 'TA>NN', 'TC>NN', 'TG>NN', 'TT>NN']
        # MutationLong
        # N:C[TC>AA]C
        specific_mutations_df = \
        combined2_df[combined2_df['Mutation'].str.contains(mutation_sub_type[0:3])].groupby('ReplicationTime')[
            'MutationCount'].sum().reset_index()
    elif mutation_main_type == ID:
        # mutation_sub_types = ['1bp Deletion', '1bp Insertion', '>1bp Deletion', '>1bp Insertion', 'Microhomology']
        # MutationLong
        # N:1:Del:T:5
        # T:1:Ins:T:5
        if mutation_sub_type == '1bp Deletion':
            specific_mutations_df = \
                combined2_df[combined2_df['Mutation'].str.contains('1:Del')].groupby(
                    'ReplicationTime')['MutationCount'].sum().reset_index()
        elif mutation_sub_type == '1bp Insertion':
            specific_mutations_df = \
                combined2_df[combined2_df['Mutation'].str.contains('1:Ins')].groupby(
                    'ReplicationTime')['MutationCount'].sum().reset_index()
        elif mutation_sub_type == '>1bp Deletion':
            specific_mutations_df = \
                combined2_df[combined2_df['Mutation'].str.contains(':Del:R:')].groupby(
                    'ReplicationTime')['MutationCount'].sum().reset_index()
        elif mutation_sub_type == '>1bp Insertion':
            specific_mutations_df = \
                combined2_df[combined2_df['Mutation'].str.contains(':Ins:R:')].groupby(
                    'ReplicationTime')['MutationCount'].sum().reset_index()
        elif mutation_sub_type == 'Microhomology':
            specific_mutations_df = \
                combined2_df[combined2_df['Mutation'].str.contains(':Del:M:')].groupby(
                    'ReplicationTime')['MutationCount'].sum().reset_index()

    return specific_mutations_df


# Note that plt.close(), plt.clf(), and plt.cla() would not close memory
# Referenced the following post for the function below:
# https://stackoverflow.com/questions/28757348/how-to-clear-memory-completely-of-all-matplotlib-plots
def clear_plotting_memory():
    usedbackend = matplotlib.get_backend()
    matplotlib.use(usedbackend)
    allfignums = matplotlib.pyplot.get_fignums()
    for i in allfignums:
        fig = matplotlib.pyplot.figure(i)
        fig.clear()
        matplotlib.pyplot.close(fig)


def plot_in_pdf(outputDir,
                jobname,
                mutation_main_type,
                mutation_sub_type,
                asymmetry_type,
                strand1,
                strand2,
                strand4,
                all_normalized_mutation_density,
                list_all_sims_normalized_mutation_density,
                strand1_normalized_mutation_density,
                list_strand1_sims_normalized_mutation_density,
                strand2_normalized_mutation_density,
                list_strand2_sims_normalized_mutation_density,
                strand3_questionable_normalized_mutation_density,
                list_strand3_questionable_sims_normalized_mutation_density,
                strand4_nontranscribed_normalized_mutation_density,
                list_strand4_nontranscribed_sims_normalized_mutation_density,
                genic_normalized_mutation_density,
                list_genic_sims_normalized_mutation_density,
                pdf):

    # plotting
    if mutation_main_type == SBS:
        title = 'All Substitutions'
        color = 'royalblue'
        fillcolor = 'lightblue'
    elif mutation_main_type == DBS:
        title = 'All Doublets'
        color = 'crimson'
        fillcolor = 'lightpink'
    elif mutation_main_type == ID:
        title = 'All Indels'
        color = 'yellowgreen'
        fillcolor = 'lightgreen'


    if mutation_sub_type is not None:
        title = mutation_sub_type

    if asymmetry_type == REPLICATIONSTRANDBIAS:
        # SBS  	All Leading Lagging —> One row: 3
        # DBS 	All Leading Lagging Questionable —> Two rows: 1,3
        # ID    All Leading Lagging Questionable —> Two rows: 1,3

        if mutation_main_type == SBS:
            # SBS  	All Leading Lagging —> One row: 3
            fig_width, fig_height = 50, 15.3125  # Set a wider figure size for horizontal layout # 50, 15
            fig = plt.figure(figsize=(fig_width, fig_height))

            ax_left = fig.add_axes([0.0, 0.1, 0.3, 0.8])  # [left, middle, right]
            simulationsLows, simulationsMeans, simulationsHighs = takeAverage(list_all_sims_normalized_mutation_density)
            plot_normalized_mutation_density_figure_with_simulations(title,
                                                                     all_normalized_mutation_density,
                                                                     simulationsLows,
                                                                     simulationsMeans,
                                                                     simulationsHighs,
                                                                     color,
                                                                     fillcolor,
                                                                     ax_left)

            ax_middle = fig.add_axes([0.35, 0.1, 0.3, 0.8])  # [left, middle, right]
            simulationsLows, simulationsMeans, simulationsHighs = takeAverage(list_strand1_sims_normalized_mutation_density)
            plot_normalized_mutation_density_figure_with_simulations(title + ' on ' + strand1,
                                                                     strand1_normalized_mutation_density,
                                                                     simulationsLows,
                                                                     simulationsMeans,
                                                                     simulationsHighs,
                                                                     color,
                                                                     fillcolor,
                                                                     ax_middle)

            ax_right = fig.add_axes([0.7, 0.1, 0.3, 0.8])  # [left, middle, right]
            simulationsLows, simulationsMeans, simulationsHighs = takeAverage(list_strand2_sims_normalized_mutation_density)
            plot_normalized_mutation_density_figure_with_simulations(title + ' on ' + strand2,
                                                                     strand2_normalized_mutation_density,
                                                                     simulationsLows,
                                                                     simulationsMeans,
                                                                     simulationsHighs,
                                                                     color,
                                                                     fillcolor,
                                                                     ax_right)

        elif (mutation_main_type == DBS) or (mutation_main_type == ID):
            # DBS 	All Leading Lagging Questionable —> Two rows: 1,3
            # ID    All Leading Lagging Questionable —> Two rows: 1,3
            fig_width, fig_height = 50, 35  # Set a wider figure size for horizontal layout # 50, 35
            fig = plt.figure(figsize=(fig_width, fig_height))

            ax_top_middle = fig.add_axes([0.35, 0.5, 0.3, 0.35])
            simulationsLows, simulationsMeans, simulationsHighs = takeAverage(list_all_sims_normalized_mutation_density)
            plot_normalized_mutation_density_figure_with_simulations(title,
                                                                     all_normalized_mutation_density,
                                                                     simulationsLows,
                                                                     simulationsMeans,
                                                                     simulationsHighs,
                                                                     color,
                                                                     fillcolor,
                                                                     ax_top_middle)

            ax_bottom_left = fig.add_axes([0.0, 0.1, 0.3, 0.35])
            simulationsLows, simulationsMeans, simulationsHighs = takeAverage(list_strand1_sims_normalized_mutation_density)
            plot_normalized_mutation_density_figure_with_simulations(title + ' on ' + strand1,
                                                                     strand1_normalized_mutation_density,
                                                                     simulationsLows,
                                                                     simulationsMeans,
                                                                     simulationsHighs,
                                                                     color,
                                                                     fillcolor,
                                                                     ax_bottom_left)

            ax_bottom_middle = fig.add_axes([0.35, 0.1, 0.3, 0.35])
            simulationsLows, simulationsMeans, simulationsHighs = takeAverage(list_strand2_sims_normalized_mutation_density)
            plot_normalized_mutation_density_figure_with_simulations(title + ' on ' + strand2,
                                                                     strand2_normalized_mutation_density,
                                                                     simulationsLows,
                                                                     simulationsMeans,
                                                                     simulationsHighs,
                                                                     color,
                                                                     fillcolor,
                                                                     ax_bottom_middle)

            ax_bottom_right = fig.add_axes([0.7, 0.1, 0.3, 0.35])
            simulationsLows, simulationsMeans, simulationsHighs = takeAverage(list_strand3_questionable_sims_normalized_mutation_density)
            plot_normalized_mutation_density_figure_with_simulations(title + ' ' + QUESTIONABLE,
                                                                     strand3_questionable_normalized_mutation_density,
                                                                     simulationsLows,
                                                                     simulationsMeans,
                                                                     simulationsHighs,
                                                                     color,
                                                                     fillcolor,
                                                                     ax_bottom_right)

    elif asymmetry_type == TRANSCRIPTIONSTRANDBIAS:
        # SBS All Transcribed UnTranscribed NonTranscribed —> Two rows:  1,3
        # DBS All Transcribed UnTranscribed NonTranscribed Questionable —> Two rows: 3,2
        # ID  All Transcribed UnTranscribed NonTranscribed Questionable —> Two rows: 3,2
        if mutation_main_type == SBS:
            # SBS All Transcribed UnTranscribed NonTranscribed —> Two rows:  1,3
            fig_width, fig_height = 50, 35  # Set a wider figure size for horizontal layout
            fig = plt.figure(figsize=(fig_width, fig_height))

            ax_top_middle = fig.add_axes([0.35, 0.5, 0.3, 0.35])
            simulationsLows, simulationsMeans, simulationsHighs = takeAverage(list_all_sims_normalized_mutation_density)
            plot_normalized_mutation_density_figure_with_simulations(title,
                                                                     all_normalized_mutation_density,
                                                                     simulationsLows,
                                                                     simulationsMeans,
                                                                     simulationsHighs,
                                                                     color,
                                                                     fillcolor,
                                                                     ax_top_middle)

            ax_bottom_left = fig.add_axes([0.0, 0.1, 0.3, 0.35])
            simulationsLows, simulationsMeans, simulationsHighs = takeAverage(list_strand1_sims_normalized_mutation_density)
            plot_normalized_mutation_density_figure_with_simulations(title + ' on ' + strand1,
                                                                     strand1_normalized_mutation_density,
                                                                     simulationsLows,
                                                                     simulationsMeans,
                                                                     simulationsHighs,
                                                                     color,
                                                                     fillcolor,
                                                                     ax_bottom_left)

            ax_bottom_middle = fig.add_axes([0.35, 0.1, 0.3, 0.35])
            simulationsLows, simulationsMeans, simulationsHighs = takeAverage(list_strand2_sims_normalized_mutation_density)
            plot_normalized_mutation_density_figure_with_simulations(title + ' on ' + strand2,
                                                                     strand2_normalized_mutation_density,
                                                                     simulationsLows,
                                                                     simulationsMeans,
                                                                     simulationsHighs,
                                                                     color,
                                                                     fillcolor,
                                                                     ax_bottom_middle)

            ax_bottom_right = fig.add_axes([0.7, 0.1, 0.3, 0.35])
            simulationsLows, simulationsMeans, simulationsHighs = takeAverage(list_strand4_nontranscribed_sims_normalized_mutation_density)
            plot_normalized_mutation_density_figure_with_simulations(title + ' on ' + strand4,
                                                                     strand4_nontranscribed_normalized_mutation_density,
                                                                     simulationsLows,
                                                                     simulationsMeans,
                                                                     simulationsHighs,
                                                                     color,
                                                                     fillcolor,
                                                                     ax_bottom_right)

        elif (mutation_main_type == DBS) or (mutation_main_type == ID):
            # DBS All Transcribed UnTranscribed NonTranscribed Questionable —> Two rows: 3,2
            # ID  All Transcribed UnTranscribed NonTranscribed Questionable —> Two rows: 3,2
            fig_width, fig_height = 50, 35  # Set a wider figure size for horizontal layout
            fig = plt.figure(figsize=(fig_width, fig_height))

            ax_top_left = fig.add_axes([0.0, 0.3, 0.3, 0.35])
            simulationsLows, simulationsMeans, simulationsHighs = takeAverage(list_all_sims_normalized_mutation_density)
            plot_normalized_mutation_density_figure_with_simulations(title,
                                                                     all_normalized_mutation_density,
                                                                     simulationsLows,
                                                                     simulationsMeans,
                                                                     simulationsHighs,
                                                                     color,
                                                                     fillcolor,
                                                                     ax_top_left)

            ax_top_middle = fig.add_axes([0.35, 0.5, 0.3, 0.35])
            simulationsLows, simulationsMeans, simulationsHighs = takeAverage(list_strand1_sims_normalized_mutation_density)
            plot_normalized_mutation_density_figure_with_simulations(title + ' on ' + strand1,
                                                                     strand1_normalized_mutation_density,
                                                                     simulationsLows,
                                                                     simulationsMeans,
                                                                     simulationsHighs,
                                                                     color,
                                                                     fillcolor,
                                                                     ax_top_middle)

            ax_top_right = fig.add_axes([0.7, 0.5, 0.3, 0.35])
            simulationsLows, simulationsMeans, simulationsHighs = takeAverage(list_strand2_sims_normalized_mutation_density)
            plot_normalized_mutation_density_figure_with_simulations(title + ' on ' + strand2,
                                                                     strand2_normalized_mutation_density,
                                                                     simulationsLows,
                                                                     simulationsMeans,
                                                                     simulationsHighs,
                                                                     color,
                                                                     fillcolor,
                                                                     ax_top_right)

            ax_bottom_middle = fig.add_axes([0.35, 0.1, 0.3, 0.35])
            simulationsLows, simulationsMeans, simulationsHighs = takeAverage(list_strand4_nontranscribed_sims_normalized_mutation_density)
            plot_normalized_mutation_density_figure_with_simulations(title + ' on ' + strand4,
                                                                     strand4_nontranscribed_normalized_mutation_density,
                                                                     simulationsLows,
                                                                     simulationsMeans,
                                                                     simulationsHighs,
                                                                     color,
                                                                     fillcolor,
                                                                     ax_bottom_middle)

            ax_bottom_right = fig.add_axes([0.7, 0.1, 0.3, 0.35])
            simulationsLows, simulationsMeans, simulationsHighs = takeAverage(list_strand3_questionable_sims_normalized_mutation_density)
            plot_normalized_mutation_density_figure_with_simulations(title + ' ' + QUESTIONABLE,
                                                                     strand3_questionable_normalized_mutation_density,
                                                                     simulationsLows,
                                                                     simulationsMeans,
                                                                     simulationsHighs,
                                                                     color,
                                                                     fillcolor,
                                                                     ax_bottom_right)


    elif asymmetry_type == GENICINTERGENICBIAS:
        # SBS All Genic Intergenic —> One row: 3
        # DBS All Genic Intergenic Questionable —> Two rows: 1,3
        # ID  All Genic Intergenic Questionable —> Two rows: 1,3
        if mutation_main_type == SBS:
            # SBS  	All Leading Lagging —> One row: 3
            fig_width, fig_height = 50, 15.3125  # Set a wider figure size for horizontal layout
            fig = plt.figure(figsize=(fig_width, fig_height))

            ax_left = fig.add_axes([0.0, 0.1, 0.3, 0.8])  # [left, middle, right]
            simulationsLows, simulationsMeans, simulationsHighs = takeAverage(list_all_sims_normalized_mutation_density)
            plot_normalized_mutation_density_figure_with_simulations(title,
                                                                     all_normalized_mutation_density,
                                                                     simulationsLows,
                                                                     simulationsMeans,
                                                                     simulationsHighs,
                                                                     color,
                                                                     fillcolor,
                                                                     ax_left)

            ax_middle = fig.add_axes([0.35, 0.1, 0.3, 0.8])  # [left, middle, right]
            simulationsLows, simulationsMeans, simulationsHighs = takeAverage(list_genic_sims_normalized_mutation_density)
            plot_normalized_mutation_density_figure_with_simulations(title + ' on Genic Regions',
                                                                     genic_normalized_mutation_density,
                                                                     simulationsLows,
                                                                     simulationsMeans,
                                                                     simulationsHighs,
                                                                     color,
                                                                     fillcolor,
                                                                     ax_middle)

            ax_right = fig.add_axes([0.7, 0.1, 0.3, 0.8])  # [left, middle, right]
            simulationsLows, simulationsMeans, simulationsHighs = takeAverage(list_strand4_nontranscribed_sims_normalized_mutation_density)
            plot_normalized_mutation_density_figure_with_simulations(title + ' on Intergenic Regions',
                                                                     strand4_nontranscribed_normalized_mutation_density,
                                                                     simulationsLows,
                                                                     simulationsMeans,
                                                                     simulationsHighs,
                                                                     color,
                                                                     fillcolor,
                                                                     ax_right)

        elif (mutation_main_type == DBS) or (mutation_main_type == ID):
            # DBS All Genic Intergenic Questionable —> Two rows: 1,3
            # ID  All Genic Intergenic Questionable —> Two rows: 1,3
            fig_width, fig_height = 50, 35  # Set a wider figure size for horizontal layout
            fig = plt.figure(figsize=(fig_width, fig_height))

            ax_top_middle = fig.add_axes([0.35, 0.5, 0.3, 0.35])
            simulationsLows, simulationsMeans, simulationsHighs = takeAverage(list_all_sims_normalized_mutation_density)
            plot_normalized_mutation_density_figure_with_simulations(title,
                                                                     all_normalized_mutation_density,
                                                                     simulationsLows,
                                                                     simulationsMeans,
                                                                     simulationsHighs,
                                                                     color,
                                                                     fillcolor,
                                                                     ax_top_middle)

            ax_bottom_left = fig.add_axes([0.0, 0.1, 0.3, 0.35])
            simulationsLows, simulationsMeans, simulationsHighs = takeAverage(list_genic_sims_normalized_mutation_density)
            plot_normalized_mutation_density_figure_with_simulations(title + ' on Genic Regions',
                                                                     genic_normalized_mutation_density,
                                                                     simulationsLows,
                                                                     simulationsMeans,
                                                                     simulationsHighs,
                                                                     color,
                                                                     fillcolor,
                                                                     ax_bottom_left)

            ax_bottom_middle = fig.add_axes([0.35, 0.1, 0.3, 0.35])
            simulationsLows, simulationsMeans, simulationsHighs = takeAverage(list_strand4_nontranscribed_sims_normalized_mutation_density)
            plot_normalized_mutation_density_figure_with_simulations(title + ' on Intergenic Regions',
                                                                     strand4_nontranscribed_normalized_mutation_density,
                                                                     simulationsLows,
                                                                     simulationsMeans,
                                                                     simulationsHighs,
                                                                     color,
                                                                     fillcolor,
                                                                     ax_bottom_middle)

            ax_bottom_right = fig.add_axes([0.7, 0.1, 0.3, 0.35])
            simulationsLows, simulationsMeans, simulationsHighs = takeAverage(list_strand3_questionable_sims_normalized_mutation_density)
            plot_normalized_mutation_density_figure_with_simulations(title + ' ' + QUESTIONABLE,
                                                                     strand3_questionable_normalized_mutation_density,
                                                                     simulationsLows,
                                                                     simulationsMeans,
                                                                     simulationsHighs,
                                                                     color,
                                                                     fillcolor,
                                                                     ax_bottom_right)

    pdf.savefig(bbox_inches='tight')  # Save the current figure into the PDF

# Fill for missing replication timing bins with zero and return complete df
def get_strand_complete(main_df, column_of_interest, strand_value, complete_index):
    strand_df = main_df[main_df[column_of_interest] == strand_value][['ReplicationTime', 'MutationCount']]
    df_complete = pd.DataFrame({'ReplicationTime': complete_index})
    strand_df_complete = df_complete.merge(strand_df, on='ReplicationTime', how='left').fillna(0)
    return strand_df_complete

def get_sims_mutation_sub_type_specific_normalized_mutation_densities(asymmetry_type,
                                                                      mutation_main_type,
                                                                      mutation_sub_type,
                                                                      number_of_sims,
                                                                      sims_combined2_df,
                                                                      complete_index,
                                                                      number_of_attributable_bases_array,
                                                                      column_of_interest,
                                                                      strand1_value,
                                                                      strand2_value,
                                                                      strand3_value,
                                                                      strand4_value):

    list_all_sims_normalized_mutation_density = []
    list_strand1_sims_normalized_mutation_density = []
    list_strand2_sims_normalized_mutation_density = []
    list_strand3_questionable_sims_normalized_mutation_density = []
    list_strand4_nontranscribed_sims_normalized_mutation_density = []
    list_genic_sims_normalized_mutation_density = []

    # Group by 'Simulation' once
    grouped_sims = sims_combined2_df.groupby('Simulation')

    for sim in range(1, number_of_sims + 1):
        # Filter the DataFrame for the current simulation
        # sim_df = sims_combined2_df[sims_combined2_df['Simulation'] == sim]
        sim_df = grouped_sims.get_group(sim)

        # sim_df.columns.values: ['Simulation' 'Mutation' 'ReplicationStrand' 'ReplicationTime' 'MutationCount']
        specific_mutations_df = get_specific_mutations_df(mutation_main_type,
                                                          mutation_sub_type,
                                                          sim_df)

        df_complete = pd.DataFrame({'ReplicationTime': complete_index})
        specific_mutations_df_complete = df_complete.merge(specific_mutations_df, on='ReplicationTime', how='left').fillna(0)

        all_normalized_mutation_density = get_normalized_mutation_density(
            specific_mutations_df_complete['MutationCount'].values,
            number_of_attributable_bases_array)

        list_all_sims_normalized_mutation_density.append(all_normalized_mutation_density)

        combined3_df = get_combined3_df(mutation_main_type,
                                    mutation_sub_type,
                                    column_of_interest,
                                    sim_df)
        # combined3_df.columns.values: ['ReplicationStrand' 'ReplicationTime' 'MutationCount']

        # Process strand1
        # strand1 Leading for RSA, Transcribed for TSA, and Genic for GIGA
        strand1_df_complete = get_strand_complete(combined3_df, column_of_interest, strand1_value, complete_index)
        strand1_normalized_mutation_density = get_normalized_mutation_density(
            strand1_df_complete['MutationCount'].values, number_of_attributable_bases_array)
        list_strand1_sims_normalized_mutation_density.append(strand1_normalized_mutation_density)

        # Process strand2
        # strand2_df Lagging for RSA, UnTranscribed for TSA, and Intergenic for GIGA
        strand2_df_complete = get_strand_complete(combined3_df, column_of_interest, strand2_value, complete_index)
        strand2_normalized_mutation_density = get_normalized_mutation_density(
            strand2_df_complete['MutationCount'].values, number_of_attributable_bases_array)
        list_strand2_sims_normalized_mutation_density.append(strand2_normalized_mutation_density)

        # Process strand3
        # strand3 Questionable for RSA, TSA, and GIGA
        strand3_df_complete = get_strand_complete(combined3_df, column_of_interest, strand3_value, complete_index)
        strand3_questionable_normalized_mutation_density = get_normalized_mutation_density(
            strand3_df_complete['MutationCount'].values, number_of_attributable_bases_array)
        list_strand3_questionable_sims_normalized_mutation_density.append(strand3_questionable_normalized_mutation_density)

        if (asymmetry_type == GENICINTERGENICBIAS) or (asymmetry_type == TRANSCRIPTIONSTRANDBIAS):
            # Process strand4
            # strand4 NonTranscribed for TSA and GIGA
            strand4_df_complete = get_strand_complete(combined3_df, column_of_interest, strand4_value, complete_index)
            strand4_nontranscribed_normalized_mutation_density = get_normalized_mutation_density(
                strand4_df_complete['MutationCount'].values,
                number_of_attributable_bases_array)
            list_strand4_nontranscribed_sims_normalized_mutation_density.append(strand4_nontranscribed_normalized_mutation_density)

            # Process genic regions
            # strand1 + strand2 = Genic for GIGA
            genic_normalized_mutation_density = get_normalized_mutation_density(
                strand1_df_complete['MutationCount'].values + strand2_df_complete['MutationCount'].values,
                number_of_attributable_bases_array)
            list_genic_sims_normalized_mutation_density.append(genic_normalized_mutation_density)

    return list_all_sims_normalized_mutation_density, \
        list_strand1_sims_normalized_mutation_density,\
        list_strand2_sims_normalized_mutation_density, \
        list_strand3_questionable_sims_normalized_mutation_density, \
        list_strand4_nontranscribed_sims_normalized_mutation_density, \
        list_genic_sims_normalized_mutation_density


def get_sims_normalized_mutation_densities(asymmetry_type,
                                            number_of_sims,
                                            sims_combined2_df,
                                            complete_index,
                                            number_of_attributable_bases_array,
                                            column_of_interest,
                                            strand1_value,
                                            strand2_value,
                                            strand3_value,
                                            strand4_value):

    list_all_sims_normalized_mutation_density = []
    list_strand1_sims_normalized_mutation_density = []
    list_strand2_sims_normalized_mutation_density = []
    list_strand3_questionable_sims_normalized_mutation_density = []
    list_strand4_nontranscribed_sims_normalized_mutation_density = []
    list_genic_sims_normalized_mutation_density = []

    # sims_combined2_df.columns.values: ['Simulation' 'Mutation' 'ReplicationStrand' 'ReplicationTime' 'MutationCount']
    # sims_combined2_df.columns.values: ['Simulation' 'Mutation' 'TranscriptionStrand' 'ReplicationTime' 'MutationCount']

    # all_mutations_df.columns.values: ['Simulation' 'ReplicationTime' 'MutationCount']
    all_mutations_df = sims_combined2_df.groupby(['Simulation', 'ReplicationTime'])['MutationCount'].sum().reset_index()

    # e.g.: combined3_df.columns.values: ['Simulation' 'ReplicationStrand' 'ReplicationTime' 'MutationCount']
    # strand1, strand2 and strand3 mutations
    combined3_df = sims_combined2_df.groupby(['Simulation', column_of_interest, 'ReplicationTime'])['MutationCount'].sum().reset_index()
    # combined3_df.columns.values: ['Simulation' 'ReplicationStrand' 'ReplicationTime' 'MutationCount']

    # Create a dictionary to hold simulation data for quick access
    grouped_sims = all_mutations_df.groupby('Simulation')

    # Iterate through each unique simulation
    for sim in range(1, number_of_sims + 1):
        sim_df = grouped_sims.get_group(sim)
        df_complete = pd.DataFrame({'ReplicationTime': complete_index})
        all_mutations_df_complete = df_complete.merge(sim_df, on='ReplicationTime', how='left').fillna(0)

        all_normalized_mutation_density = get_normalized_mutation_density(
            all_mutations_df_complete['MutationCount'].values, number_of_attributable_bases_array)

        list_all_sims_normalized_mutation_density.append(all_normalized_mutation_density)

    # Group by 'Simulation' once
    grouped_sims = combined3_df.groupby('Simulation')

    for sim in range(1, number_of_sims + 1):
        # Filter the DataFrame for the current simulation
        sim_df = grouped_sims.get_group(sim)
        # sim_df.columns.values: ['Simulation' 'ReplicationStrand' 'ReplicationTime' 'MutationCount']

        # Process strand1
        # strand1 Leading for RSA, Transcribed for TSA, and Genic for GIGA
        strand1_df_complete = get_strand_complete(sim_df, column_of_interest, strand1_value, complete_index)
        strand1_normalized_mutation_density = get_normalized_mutation_density(
            strand1_df_complete['MutationCount'].values, number_of_attributable_bases_array)
        list_strand1_sims_normalized_mutation_density.append(strand1_normalized_mutation_density)

        # Process strand2
        # strand2 Lagging for RSA, UnTranscribed for TSA, and Intergenic for GIGA
        strand2_df_complete = get_strand_complete(sim_df, column_of_interest, strand2_value, complete_index)
        strand2_normalized_mutation_density = get_normalized_mutation_density(
            strand2_df_complete['MutationCount'].values, number_of_attributable_bases_array)
        list_strand2_sims_normalized_mutation_density.append(strand2_normalized_mutation_density)

        # Process strand3
        # strand3 Questionable for RSA, TSA, and GIGA
        strand3_questionable_df_complete = get_strand_complete(sim_df, column_of_interest, strand3_value, complete_index)
        strand3_questionable_normalized_mutation_density = get_normalized_mutation_density(
            strand3_questionable_df_complete['MutationCount'].values, number_of_attributable_bases_array)
        list_strand3_questionable_sims_normalized_mutation_density.append(strand3_questionable_normalized_mutation_density)

        if (asymmetry_type == GENICINTERGENICBIAS) or (asymmetry_type == TRANSCRIPTIONSTRANDBIAS):
            # Process strand4
            # strand4 NonTranscribed for TSA and GIGA
            strand4_nontranscribed_df_complete = get_strand_complete(sim_df, column_of_interest, strand4_value, complete_index)
            strand4_nontranscribed_normalized_mutation_density = get_normalized_mutation_density(
                strand4_nontranscribed_df_complete['MutationCount'].values, number_of_attributable_bases_array)
            list_strand4_nontranscribed_sims_normalized_mutation_density.append(strand4_nontranscribed_normalized_mutation_density)

            # Process Genic Regions
            genic_normalized_mutation_density = get_normalized_mutation_density(
                strand1_df_complete['MutationCount'].values + strand2_df_complete['MutationCount'].values,
                number_of_attributable_bases_array)
            list_genic_sims_normalized_mutation_density.append(genic_normalized_mutation_density)

    return list_all_sims_normalized_mutation_density, \
        list_strand1_sims_normalized_mutation_density,\
        list_strand2_sims_normalized_mutation_density, \
        list_strand3_questionable_sims_normalized_mutation_density, \
        list_strand4_nontranscribed_sims_normalized_mutation_density, \
        list_genic_sims_normalized_mutation_density


def get_mutation_sub_type_specific_normalized_mutation_densities(asymmetry_type,
                                                                 mutation_main_type,
                                                                 mutation_sub_type,
                                                                 combined2_df,
                                                                 complete_index,
                                                                 number_of_attributable_bases_array,
                                                                 column_of_interest,
                                                                 strand1_value,
                                                                 strand2_value,
                                                                 strand3_value,
                                                                 strand4_value):

    all_normalized_mutation_density = None
    strand1_normalized_mutation_density = None
    strand2_normalized_mutation_density = None
    strand3_questionable_normalized_mutation_density = None
    strand4_nontranscribed_normalized_mutation_density = None
    genic_normalized_mutation_density = None

    specific_mutations_df = get_specific_mutations_df(mutation_main_type,
                                                      mutation_sub_type,
                                                      combined2_df)

    df_complete = pd.DataFrame({'ReplicationTime': complete_index})
    specific_mutations_df_complete = df_complete.merge(specific_mutations_df, on='ReplicationTime', how='left').fillna(0)

    all_normalized_mutation_density = get_normalized_mutation_density(
        specific_mutations_df_complete['MutationCount'].values,
        number_of_attributable_bases_array)

    combined3_df = get_combined3_df(mutation_main_type,
                                    mutation_sub_type,
                                    column_of_interest,
                                    combined2_df)

    combined3_df = combined3_df[(combined3_df['ReplicationTime'] >= 1) & (combined3_df['ReplicationTime'] <= 10)]

    # strand1 Leading for RSA, Transcribed for TSA, and Genic for GIGA
    strand1_df_complete = get_strand_complete(combined3_df, column_of_interest, strand1_value, complete_index)
    strand1_normalized_mutation_density = get_normalized_mutation_density(
        strand1_df_complete['MutationCount'].values, number_of_attributable_bases_array)

    # strand2_df Lagging for RSA, UnTranscribed for TSA, and Intergenic for GIGA
    strand2_df_complete = get_strand_complete(combined3_df, column_of_interest, strand2_value, complete_index)
    strand2_normalized_mutation_density = get_normalized_mutation_density(
        strand2_df_complete['MutationCount'].values, number_of_attributable_bases_array)

    # strand3 Questionable for RSA, TSA, and GIGA
    strand3_df_complete = get_strand_complete(combined3_df, column_of_interest, strand3_value, complete_index)
    strand3_questionable_normalized_mutation_density = get_normalized_mutation_density(
        strand3_df_complete['MutationCount'].values, number_of_attributable_bases_array)

    if (asymmetry_type == GENICINTERGENICBIAS) or (asymmetry_type == TRANSCRIPTIONSTRANDBIAS):
        # strand4 NonTranscribed for TSA and GIGA
        strand4_df_complete = get_strand_complete(combined3_df, column_of_interest, strand4_value, complete_index)
        strand4_nontranscribed_normalized_mutation_density = get_normalized_mutation_density(
            strand4_df_complete['MutationCount'].values, number_of_attributable_bases_array)

        # strand1 + strand2 = Genic for GIGA
        genic_normalized_mutation_density = get_normalized_mutation_density(
            strand1_df_complete['MutationCount'].values + strand2_df_complete['MutationCount'].values,
            number_of_attributable_bases_array)

    return all_normalized_mutation_density, \
        strand1_normalized_mutation_density,\
        strand2_normalized_mutation_density, \
        strand3_questionable_normalized_mutation_density, \
        strand4_nontranscribed_normalized_mutation_density, \
        genic_normalized_mutation_density


def get_normalized_mutation_densities(asymmetry_type,
                                      combined2_df,
                                      complete_index,
                                      number_of_attributable_bases_array,
                                      column_of_interest,
                                      strand1_value,
                                      strand2_value,
                                      strand3_value,
                                      strand4_value):

    all_normalized_mutation_density = None
    strand1_normalized_mutation_density = None
    strand2_normalized_mutation_density = None
    strand3_questionable_normalized_mutation_density = None
    strand4_nontranscribed_normalized_mutation_density = None
    genic_normalized_mutation_density = None

    all_mutations_df = combined2_df.groupby('ReplicationTime')['MutationCount'].sum().reset_index()
    df_complete = pd.DataFrame({'ReplicationTime': complete_index})
    all_mutations_df_complete = df_complete.merge(all_mutations_df, on='ReplicationTime', how='left').fillna(0)

    all_normalized_mutation_density = get_normalized_mutation_density(
        all_mutations_df_complete['MutationCount'].values, number_of_attributable_bases_array)

    # strand1, strand2 and strand3 mutations
    combined3_df = combined2_df.groupby([column_of_interest, 'ReplicationTime'])['MutationCount'].sum().reset_index()
    combined3_df = combined3_df[(combined3_df['ReplicationTime'] >= 1) & (combined3_df['ReplicationTime'] <= 10)]

    # strand1 Leading for RSA, Transcribed for TSA, and Genic for GIGA
    strand1_df_complete = get_strand_complete(combined3_df, column_of_interest, strand1_value,
                                              complete_index)

    # strand1 Leading for RSA, Transcribed for TSA
    strand1_normalized_mutation_density = get_normalized_mutation_density(
        strand1_df_complete['MutationCount'].values, number_of_attributable_bases_array)

    # strand2 Lagging for RSA, UnTranscribed for TSA, and Intergenic for GIGA
    strand2_df_complete = get_strand_complete(combined3_df, column_of_interest, strand2_value,
                                              complete_index)

    # strand2 Lagging for RSA, UnTranscribed for TSA
    strand2_normalized_mutation_density = get_normalized_mutation_density(
        strand2_df_complete['MutationCount'].values, number_of_attributable_bases_array)

    # strand3 Questionable for RSA, TSA, and GIGA
    strand3_questionable_df_complete = get_strand_complete(combined3_df, column_of_interest, strand3_value,
                                                           complete_index)

    # strand3 Questionable for RSA and TSA
    strand3_questionable_normalized_mutation_density = get_normalized_mutation_density(
        strand3_questionable_df_complete['MutationCount'].values, number_of_attributable_bases_array)

    if (asymmetry_type == GENICINTERGENICBIAS) or (asymmetry_type == TRANSCRIPTIONSTRANDBIAS):
        # strand4 NonTranscribed for TSA and GIGA
        strand4_nontranscribed_df_complete = get_strand_complete(combined3_df, column_of_interest, strand4_value,
                                                                 complete_index)

        # strand4 for Nontranscribed for TSA
        strand4_nontranscribed_normalized_mutation_density = get_normalized_mutation_density(
            strand4_nontranscribed_df_complete['MutationCount'].values, number_of_attributable_bases_array)

        genic_normalized_mutation_density = get_normalized_mutation_density(
            strand1_df_complete['MutationCount'].values + strand2_df_complete['MutationCount'].values,
            number_of_attributable_bases_array)

    return all_normalized_mutation_density, \
        strand1_normalized_mutation_density,\
        strand2_normalized_mutation_density, \
        strand3_questionable_normalized_mutation_density, \
        strand4_nontranscribed_normalized_mutation_density, \
        genic_normalized_mutation_density

def get_signature_specific_normalized_mutation_density_from_np_array(asymmetry_type,
                                                  mutation_main_type,
                                                  signature_index,
                                                  number_of_attributable_bases_array,
                                                  simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                  simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                  simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                  simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                                                  simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                                                  simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array):

    all_normalized_mutation_density = None
    strand1_normalized_mutation_density = None
    strand2_normalized_mutation_density = None
    strand3_questionable_normalized_mutation_density = None
    strand4_nontranscribed_normalized_mutation_density = None
    genic_normalized_mutation_density = None

    if mutation_main_type == SBS:
        if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
            array_of_interest = simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array

        elif asymmetry_type == REPLICATIONSTRANDBIAS:
            array_of_interest = simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array

    elif mutation_main_type == DBS:
        if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
            array_of_interest = simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array

        elif asymmetry_type == REPLICATIONSTRANDBIAS:
            array_of_interest = simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array

    elif mutation_main_type == ID:
        if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
            array_of_interest = simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array

        elif asymmetry_type == REPLICATIONSTRANDBIAS:
            array_of_interest = simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array

    if asymmetry_type == TRANSCRIPTIONSTRANDBIAS:
        # strand1_value = 'T'
        # strand2_value = 'U'
        # strand3_value = 'Q'
        # strand4_value = 'N'
        # transcription_strand_mapping = {'T': 0, 'U': 1, 'N': 2, 'Q': 3}

        all_values = array_of_interest[0, signature_index, :, :, :]

        # strand1 = TRANSCRIBED
        strand1_values = array_of_interest[0, signature_index, :, 0, :]

        # strand2 = UNTRANSCRIBED
        strand2_values = array_of_interest[0, signature_index, :, 1, :]

        if mutation_main_type == DBS or mutation_main_type == ID:
            # strand3 = QUESTIONABLE
            strand3_values = array_of_interest[0, signature_index, :, 3, :]

        # strand4 = NONTRANSCRIBED
        strand4_values = array_of_interest[0, signature_index, :, 2, :]

    elif asymmetry_type == GENICINTERGENICBIAS:
        # strand1_value = 'T'
        # strand2_value = 'U'
        # strand3_value = 'Q'
        # strand4_value = 'N'
        # transcription_strand_mapping = {'T': 0, 'U': 1, 'N': 2, 'Q': 3}

        all_values = array_of_interest[0, signature_index, :, :, :]

        # strand1 = TRANSCRIBED
        strand1_values = array_of_interest[0, signature_index, :, 0, :]

        # strand2 = UNTRANSCRIBED
        strand2_values = array_of_interest[0, signature_index, :, 1, :]

        if mutation_main_type == DBS or mutation_main_type == ID:
            # strand3 = QUESTIONABLE
            strand3_values = array_of_interest[0, signature_index, :, 3, :]

        # strand4 = NONTRANSCRIBED
        strand4_values = array_of_interest[0, signature_index, :, 2, :]

        # TODO does questionable mean either transcribed or untranscribed?
        # TODO can we consider questionable values in genic values?
        genic_values =   strand1_values + strand2_values
        # intergenic_values =  strand4_values

    elif asymmetry_type == REPLICATIONSTRANDBIAS:
        # strand1_value = 'E'
        # strand2_value = 'A'
        # strand3_value = 'Q'
        # replication_strand_mapping = {'E': 0, 'A': 1, 'Q': 2}  # Map all possible strands

        all_values = array_of_interest[0, signature_index, :, :, :]

        # strand1 = LEADING
        strand1_values = array_of_interest[0, signature_index, :, 0, :]

        # strand2 = LAGGING
        strand2_values = array_of_interest[0, signature_index, :, 1, :]

        if (mutation_main_type == DBS) or (mutation_main_type == ID):
            # strand3 = QUESTIONABLE
            strand3_values = array_of_interest[0, signature_index, :, 2, :]

    # Accumulate (sum) values across all mutation types and strands
    # Here we sum across the mutation type axis (axis=0) and strand axis (axis=1)
    all_values = np.sum(all_values, axis=(0, 1))
    all_normalized_mutation_density = get_normalized_mutation_density(all_values, number_of_attributable_bases_array)

    strand1_replication_time_values = np.sum(strand1_values, axis=(0))
    strand1_normalized_mutation_density = get_normalized_mutation_density(strand1_replication_time_values, number_of_attributable_bases_array)

    strand2_replication_time_values = np.sum(strand2_values, axis=(0))
    strand2_normalized_mutation_density = get_normalized_mutation_density(strand2_replication_time_values, number_of_attributable_bases_array)

    if (mutation_main_type == DBS) or (mutation_main_type == ID):
        strand3_replication_time_values = np.sum(strand3_values, axis=(0))
        strand3_questionable_normalized_mutation_density = get_normalized_mutation_density(strand3_replication_time_values,
                                                                              number_of_attributable_bases_array)

    if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
        strand4_replication_time_values = np.sum(strand4_values, axis=(0))
        strand4_nontranscribed_normalized_mutation_density = get_normalized_mutation_density(strand4_replication_time_values,
                                                                              number_of_attributable_bases_array)

    if asymmetry_type == GENICINTERGENICBIAS:
        genic_replication_time_values = np.sum(genic_values, axis=(0))
        genic_normalized_mutation_density = get_normalized_mutation_density(genic_replication_time_values,
                                                                              number_of_attributable_bases_array)

    return all_normalized_mutation_density, \
        strand1_normalized_mutation_density, \
        strand2_normalized_mutation_density, \
        strand3_questionable_normalized_mutation_density, \
        strand4_nontranscribed_normalized_mutation_density, \
        genic_normalized_mutation_density


def get_normalized_mutation_density_from_np_array(asymmetry_type,
                                                  mutation_main_type,
                                                  number_of_attributable_bases_array,
                                                  simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                  simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                  simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                  simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                                                  simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                                                  simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array):

    all_normalized_mutation_density = None
    strand1_normalized_mutation_density = None
    strand2_normalized_mutation_density = None
    strand3_questionable_normalized_mutation_density = None
    strand4_nontranscribed_normalized_mutation_density = None
    genic_normalized_mutation_density = None

    if mutation_main_type == SBS:
        if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
            array_of_interest = simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array

        elif asymmetry_type == REPLICATIONSTRANDBIAS:
            array_of_interest = simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array

    elif mutation_main_type == DBS:
        if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
            array_of_interest = simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array

        elif asymmetry_type == REPLICATIONSTRANDBIAS:
            array_of_interest = simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array

    elif mutation_main_type == ID:
        if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
            array_of_interest = simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array

        elif asymmetry_type == REPLICATIONSTRANDBIAS:
            array_of_interest = simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array

    if asymmetry_type == TRANSCRIPTIONSTRANDBIAS:
        # strand1_value = 'T'
        # strand2_value = 'U'
        # strand3_value = 'Q'
        # strand4_value = 'N'
        # transcription_strand_mapping = {'T': 0, 'U': 1, 'N': 2, 'Q': 3}

        all_values = array_of_interest[0, -1, :, :, :]

        # strand1 = TRANSCRIBED
        strand1_values = array_of_interest[0, -1, :, 0, :]

        # strand2 = UNTRANSCRIBED
        strand2_values = array_of_interest[0, -1, :, 1, :]

        if mutation_main_type == DBS or mutation_main_type == ID:
            # strand3 = QUESTIONABLE
            strand3_values = array_of_interest[0, -1, :, 3, :]

        # strand4 = NONTRANSCRIBED
        strand4_values = array_of_interest[0, -1, :, 2, :]

    elif asymmetry_type == GENICINTERGENICBIAS:
        # strand1_value = 'T'
        # strand2_value = 'U'
        # strand3_value = 'Q'
        # strand4_value = 'N'
        # transcription_strand_mapping = {'T': 0, 'U': 1, 'N': 2, 'Q': 3}

        all_values = array_of_interest[0, -1, :, :, :]

        # strand1 = TRANSCRIBED
        strand1_values = array_of_interest[0, -1, :, 0, :]

        # strand2 = UNTRANSCRIBED
        strand2_values = array_of_interest[0, -1, :, 1, :]

        if mutation_main_type == DBS or mutation_main_type == ID:
            # strand3 = QUESTIONABLE
            strand3_values = array_of_interest[0, -1, :, 3, :]

        # strand4 = NONTRANSCRIBED
        strand4_values = array_of_interest[0, -1, :, 2, :]

        # TODO does questionable mean either transcribed or untranscribed?
        # TODO can we consider questionable values in genic values?
        genic_values =   strand1_values + strand2_values
        # intergenic_values =  strand4_values

    elif asymmetry_type == REPLICATIONSTRANDBIAS:
        # strand1_value = 'E'
        # strand2_value = 'A'
        # strand3_value = 'Q'
        # replication_strand_mapping = {'E': 0, 'A': 1, 'Q': 2}  # Map all possible strands

        all_values = array_of_interest[0, -1, :, :, :]

        # strand1 = LEADING
        strand1_values = array_of_interest[0, -1, :, 0, :]

        # strand2 = LAGGING
        strand2_values = array_of_interest[0, -1, :, 1, :]

        if (mutation_main_type == DBS) or (mutation_main_type == ID):
            # strand3 = QUESTIONABLE
            strand3_values = array_of_interest[0, -1, :, 2, :]

    # Accumulate (sum) values across all mutation types and strands
    # Here we sum across the mutation type axis (axis=0) and strand axis (axis=1)
    all_values = np.sum(all_values, axis=(0, 1))
    all_normalized_mutation_density = get_normalized_mutation_density(all_values, number_of_attributable_bases_array)

    strand1_replication_time_values = np.sum(strand1_values, axis=(0))
    strand1_normalized_mutation_density = get_normalized_mutation_density(strand1_replication_time_values, number_of_attributable_bases_array)

    strand2_replication_time_values = np.sum(strand2_values, axis=(0))
    strand2_normalized_mutation_density = get_normalized_mutation_density(strand2_replication_time_values, number_of_attributable_bases_array)

    if (mutation_main_type == DBS) or (mutation_main_type == ID):
        strand3_replication_time_values = np.sum(strand3_values, axis=(0))
        strand3_questionable_normalized_mutation_density = get_normalized_mutation_density(strand3_replication_time_values,
                                                                              number_of_attributable_bases_array)

    if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
        strand4_replication_time_values = np.sum(strand4_values, axis=(0))
        strand4_nontranscribed_normalized_mutation_density = get_normalized_mutation_density(strand4_replication_time_values,
                                                                              number_of_attributable_bases_array)

    if asymmetry_type == GENICINTERGENICBIAS:
        genic_replication_time_values = np.sum(genic_values, axis=(0))
        genic_normalized_mutation_density = get_normalized_mutation_density(genic_replication_time_values,
                                                                              number_of_attributable_bases_array)

    return all_normalized_mutation_density, \
        strand1_normalized_mutation_density, \
        strand2_normalized_mutation_density, \
        strand3_questionable_normalized_mutation_density, \
        strand4_nontranscribed_normalized_mutation_density, \
        genic_normalized_mutation_density


def get_mutation_specific_normalized_mutation_density_from_np_array(asymmetry_type,
                                                                    mutation_main_type,
                                                                    mutation_sub_type,
                                                                    ordered_SBS96_substitution_types,
                                                                    ordered_DBS78_doublet_types,
                                                                    ordered_ID83_indel_types,
                                                                    SBS96_mutation_mapping,
                                                                    DBS78_mutation_mapping,
                                                                    ID83_mutation_mapping,
                                                                    number_of_attributable_bases_array,
                                                                    simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                                    simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                                    simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                                    simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                                                                    simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                                                                    simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array):

    all_normalized_mutation_density = None
    strand1_normalized_mutation_density = None
    strand2_normalized_mutation_density = None
    strand3_questionable_normalized_mutation_density = None
    strand4_nontranscribed_normalized_mutation_density = None
    genic_normalized_mutation_density = None

    if mutation_main_type == SBS:
        mutation_indices = [SBS96_mutation_mapping[mutation_type] for mutation_type in ordered_SBS96_substitution_types if mutation_sub_type in mutation_type]

        if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
            array_of_interest = simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array

        elif asymmetry_type == REPLICATIONSTRANDBIAS:
            array_of_interest = simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array

    elif mutation_main_type == DBS:
        # mutation_sub_types = ['AC>NN', 'AT>NN', 'CC>NN', 'CG>NN', 'CT>NN', 'GC>NN', 'TA>NN', 'TC>NN', 'TG>NN', 'TT>NN']
        mutation_indices = [DBS78_mutation_mapping[mutation_type] for mutation_type in ordered_DBS78_doublet_types if mutation_sub_type[0:3] == mutation_type[0:3]]

        if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
            array_of_interest = simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array

        elif asymmetry_type == REPLICATIONSTRANDBIAS:
            array_of_interest = simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array

    elif mutation_main_type == ID:
        # mutation_sub_types = ['1bp Deletion', '1bp Insertion', '>1bp Deletion', '>1bp Insertion', 'Microhomology']

        if mutation_sub_type == '1bp Deletion':
            mutation_indices = [ID83_mutation_mapping[mutation_type] for mutation_type in ordered_ID83_indel_types if
                                '1:Del' in mutation_type]

        elif mutation_sub_type == '1bp Insertion':
            mutation_indices = [ID83_mutation_mapping[mutation_type] for mutation_type in ordered_ID83_indel_types if
                                '1:Ins' in mutation_type]

        elif mutation_sub_type == '>1bp Deletion':
            mutation_indices = [ID83_mutation_mapping[mutation_type] for mutation_type in ordered_ID83_indel_types if
                                ':Del:R:' in mutation_type]

        elif mutation_sub_type == '>1bp Insertion':
            mutation_indices = [ID83_mutation_mapping[mutation_type] for mutation_type in ordered_ID83_indel_types if
                                ':Ins:R:' in mutation_type]

        elif mutation_sub_type == 'Microhomology':
            mutation_indices = [ID83_mutation_mapping[mutation_type] for mutation_type in ordered_ID83_indel_types if
                                ':Del:M:' in mutation_type]

        if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
            array_of_interest = simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array

        elif asymmetry_type == REPLICATIONSTRANDBIAS:
            array_of_interest = simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array

    if asymmetry_type == TRANSCRIPTIONSTRANDBIAS:
        # strand1_value = 'T'
        # strand2_value = 'U'
        # strand3_value = 'Q'
        # strand4_value = 'N'
        # transcription_strand_mapping = {'T': 0, 'U': 1, 'N': 2, 'Q': 3, 'B': 4}  # Map all possible strands


        # all strands
        all_values = array_of_interest[0, -1, mutation_indices, :, :]

        # strand1 = TRANSCRIBED
        strand1_values = array_of_interest[0, -1, mutation_indices, 0, :]

        # strand2 = UNTRANSCRIBED
        strand2_values = array_of_interest[0, -1, mutation_indices, 1, :]

        if mutation_main_type == DBS or mutation_main_type == ID:
            # strand3 = QUESTIONABLE
            strand3_values = array_of_interest[0, -1, mutation_indices, 3, :]

        # strand4 = NONTRANSCRIBED
        strand4_values = array_of_interest[0, -1, mutation_indices, 2, :]

    elif asymmetry_type == GENICINTERGENICBIAS:
        # strand1_value = 'T'
        # strand2_value = 'U'
        # strand3_value = 'Q'
        # strand4_value = 'N'
        # transcription_strand_mapping = {'T': 0, 'U': 1, 'N': 2, 'Q': 3}

        all_values = array_of_interest[0, -1, mutation_indices, :, :]

        # strand1 = TRANSCRIBED
        strand1_values = array_of_interest[0, -1, mutation_indices, 0, :]

        # strand2 = UNTRANSCRIBED
        strand2_values = array_of_interest[0, -1, mutation_indices, 1, :]

        if mutation_main_type == DBS or mutation_main_type == ID:
            # strand3 = QUESTIONABLE
            strand3_values = array_of_interest[0, -1, mutation_indices, 3, :]

        # strand4 = NONTRANSCRIBED
        strand4_values = array_of_interest[0, -1, mutation_indices, 2, :]

        # TODO does questionable mean either transcribed or untranscribed?
        # TODO can we consider questionable values in genic values?
        genic_values =   strand1_values + strand2_values
        # intergenic_values =  strand4_values

    elif asymmetry_type == REPLICATIONSTRANDBIAS:
        # strand1_value = 'E'
        # strand2_value = 'A'
        # strand3_value = 'Q'
        # replication_strand_mapping = {'E': 0, 'A': 1, 'Q': 2}  # Map all possible strands

        all_values = array_of_interest[0, -1, mutation_indices, :, :]

        # strand1 = LEADING
        strand1_values = array_of_interest[0, -1, mutation_indices, 0, :]

        # strand2 = LAGGING
        strand2_values = array_of_interest[0, -1, mutation_indices, 1, :]

        if (mutation_main_type == DBS) or (mutation_main_type == ID):
            # strand3 = QUESTIONABLE
            strand3_values = array_of_interest[0, -1, mutation_indices, 2, :]

    # Accumulate (sum) values across all mutation types and strands
    # Here we sum across the mutation type axis (axis=0) and strand axis (axis=1)
    all_values = np.sum(all_values, axis=(0, 1))
    all_normalized_mutation_density = get_normalized_mutation_density(all_values, number_of_attributable_bases_array)

    strand1_replication_time_values = np.sum(strand1_values, axis=(0))
    strand1_normalized_mutation_density = get_normalized_mutation_density(strand1_replication_time_values, number_of_attributable_bases_array)

    strand2_replication_time_values = np.sum(strand2_values, axis=(0))
    strand2_normalized_mutation_density = get_normalized_mutation_density(strand2_replication_time_values, number_of_attributable_bases_array)

    if (mutation_main_type == DBS) or (mutation_main_type == ID):
        strand3_replication_time_values = np.sum(strand3_values, axis=(0))
        strand3_questionable_normalized_mutation_density = get_normalized_mutation_density(strand3_replication_time_values,
                                                                              number_of_attributable_bases_array)

    if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
        strand4_replication_time_values = np.sum(strand4_values, axis=(0))
        strand4_nontranscribed_normalized_mutation_density = get_normalized_mutation_density(strand4_replication_time_values,
                                                                              number_of_attributable_bases_array)

    if asymmetry_type == GENICINTERGENICBIAS:
        genic_replication_time_values = np.sum(genic_values, axis=(0))
        genic_normalized_mutation_density = get_normalized_mutation_density(genic_replication_time_values,
                                                                              number_of_attributable_bases_array)

    return all_normalized_mutation_density, \
        strand1_normalized_mutation_density, \
        strand2_normalized_mutation_density, \
        strand3_questionable_normalized_mutation_density, \
        strand4_nontranscribed_normalized_mutation_density, \
        genic_normalized_mutation_density


def get_mutation_specific_sims_normalized_mutation_density_from_np_array(asymmetry_type,
                        mutation_main_type,
                        mutation_sub_type,
                        ordered_SBS96_substitution_types,
                        ordered_DBS78_doublet_types,
                        ordered_ID83_indel_types,
                        SBS96_mutation_mapping,
                        DBS78_mutation_mapping,
                        ID83_mutation_mapping,
                        number_of_attributable_bases_array,
                        number_of_sims,
                        simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                        simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                        simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                        simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                        simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                        simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array):

    all_sims_normalized_mutation_density = None
    strand1_sims_normalized_mutation_density = None
    strand2_sims_normalized_mutation_density = None
    strand3_questionable_sims_normalized_mutation_density = None
    strand4_nontranscribed_sims_normalized_mutation_density = None
    genic_sims_normalized_mutation_density = None

    if mutation_main_type == SBS:

        mutation_indices = [SBS96_mutation_mapping[mutation_type] for mutation_type in ordered_SBS96_substitution_types if mutation_sub_type in mutation_type]

        if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
            array_of_interest = simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array

        elif asymmetry_type == REPLICATIONSTRANDBIAS:
            array_of_interest = simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array

    elif mutation_main_type == DBS:

        # mutation_sub_types = ['AC>NN', 'AT>NN', 'CC>NN', 'CG>NN', 'CT>NN', 'GC>NN', 'TA>NN', 'TC>NN', 'TG>NN', 'TT>NN']
        mutation_indices = [DBS78_mutation_mapping[mutation_type] for mutation_type in ordered_DBS78_doublet_types if mutation_sub_type[0:3] == mutation_type[0:3]]

        if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
            array_of_interest = simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array

        elif asymmetry_type == REPLICATIONSTRANDBIAS:
            array_of_interest = simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array

    elif mutation_main_type == ID:

        # mutation_sub_types = ['1bp Deletion', '1bp Insertion', '>1bp Deletion', '>1bp Insertion', 'Microhomology']

        if mutation_sub_type == '1bp Deletion':
            mutation_indices = [ID83_mutation_mapping[mutation_type] for mutation_type in ordered_ID83_indel_types if
                                '1:Del' in mutation_type]

        elif mutation_sub_type == '1bp Insertion':
            mutation_indices = [ID83_mutation_mapping[mutation_type] for mutation_type in ordered_ID83_indel_types if
                                '1:Ins' in mutation_type]

        elif mutation_sub_type == '>1bp Deletion':
            mutation_indices = [ID83_mutation_mapping[mutation_type] for mutation_type in ordered_ID83_indel_types if
                                ':Del:R:' in mutation_type]

        elif mutation_sub_type == '>1bp Insertion':
            mutation_indices = [ID83_mutation_mapping[mutation_type] for mutation_type in ordered_ID83_indel_types if
                                ':Ins:R:' in mutation_type]

        elif mutation_sub_type == 'Microhomology':
            mutation_indices = [ID83_mutation_mapping[mutation_type] for mutation_type in ordered_ID83_indel_types if
                                ':Del:M:' in mutation_type]

        if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
            array_of_interest = simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array

        elif asymmetry_type == REPLICATIONSTRANDBIAS:
            array_of_interest = simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array

    if asymmetry_type == TRANSCRIPTIONSTRANDBIAS:
        # strand1_value = 'T'
        # strand2_value = 'U'
        # strand3_value = 'Q'
        # strand4_value = 'N'
        # transcription_strand_mapping = {'T': 0, 'U': 1, 'N': 2, 'Q': 3}

        all_values = array_of_interest[1:number_of_sims+1, -1, mutation_indices, :, :]

        # strand1 = TRANSCRIBED
        strand1_values = array_of_interest[1:number_of_sims+1, -1, mutation_indices, 0, :]

        # strand2 = UNTRANSCRIBED
        strand2_values = array_of_interest[1:number_of_sims+1, -1, mutation_indices, 1, :]

        if mutation_main_type == DBS or mutation_main_type == ID:
            # strand3 = QUESTIONABLE
            strand3_values = array_of_interest[1:number_of_sims+1, -1, mutation_indices, 3, :]

        # strand4 = NONTRANSCRIBED
        strand4_values = array_of_interest[1:number_of_sims+1, -1, mutation_indices, 2, :]

    elif asymmetry_type == GENICINTERGENICBIAS:
        # strand1_value = 'T'
        # strand2_value = 'U'
        # strand3_value = 'Q'
        # strand4_value = 'N'
        # transcription_strand_mapping = {'T': 0, 'U': 1, 'N': 2, 'Q': 3}

        all_values = array_of_interest[1:number_of_sims+1, -1, mutation_indices, :, :]

        # strand1 = TRANSCRIBED
        strand1_values = array_of_interest[1:number_of_sims+1, -1, mutation_indices, 0, :]

        # strand2 = UNTRANSCRIBED
        strand2_values = array_of_interest[1:number_of_sims+1, -1, mutation_indices, 1, :]

        if mutation_main_type == DBS or mutation_main_type == ID:
            # strand3 = QUESTIONABLE
            strand3_values = array_of_interest[1:number_of_sims+1, -1, mutation_indices, 3, :]

        # strand4 = NONTRANSCRIBED
        strand4_values = array_of_interest[1:number_of_sims+1, -1, mutation_indices, 2, :]

        # TODO does questionable mean either transcribed or untranscribed?
        # TODO can we consider questionable values in genic values?
        genic_values =   strand1_values + strand2_values
        # intergenic_values =  strand4_values

    elif asymmetry_type == REPLICATIONSTRANDBIAS:
        # strand1_value = 'E'
        # strand2_value = 'A'
        # strand3_value = 'Q'
        # replication_strand_mapping = {'E': 0, 'A': 1, 'Q': 2}  # Map all possible strands

        all_values = array_of_interest[1:number_of_sims+1, -1, mutation_indices, :, :]

        # strand1 = LEADING
        strand1_values = array_of_interest[1:number_of_sims+1, -1, mutation_indices, 0, :]

        # strand2 = LAGGING
        strand2_values = array_of_interest[1:number_of_sims+1, -1, mutation_indices, 1, :]

        if mutation_main_type == DBS or mutation_main_type == ID:
            # strand3 = QUESTIONABLE
            strand3_values = array_of_interest[1:number_of_sims+1, -1, mutation_indices, 2, :]

    # Accumulate (sum) values across all mutation types and strands
    # Here we sum across the mutation type axis (axis=0) and strand axis (axis=1)

    # number_of_sims, mutation_types, num_of_strands, replication_times
    all_sims_replication_time_values = np.sum(all_values, axis=(1, 2))
    all_sims_normalized_mutation_density = get_normalized_mutation_density_2D(all_sims_replication_time_values, number_of_attributable_bases_array)

    # num_of_sims, mutation_types, replication_times
    strand1_sims_replication_time_values = np.sum(strand1_values, axis=(1))
    strand1_sims_normalized_mutation_density = get_normalized_mutation_density_2D(strand1_sims_replication_time_values, number_of_attributable_bases_array)

    # num_of_sims, mutation_types, replication_times
    strand2_sims_replication_time_values = np.sum(strand2_values, axis=(1))
    strand2_sims_normalized_mutation_density = get_normalized_mutation_density_2D(strand2_sims_replication_time_values, number_of_attributable_bases_array)

    if mutation_main_type == DBS or mutation_main_type == ID:
        strand3_sims_replication_time_values = np.sum(strand3_values, axis=(1))
        strand3_questionable_sims_normalized_mutation_density = get_normalized_mutation_density_2D(strand3_sims_replication_time_values, number_of_attributable_bases_array)

    if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
        strand4_sims_replication_time_values = np.sum(strand4_values, axis=(1))
        strand4_nontranscribed_sims_normalized_mutation_density = get_normalized_mutation_density_2D(strand4_sims_replication_time_values, number_of_attributable_bases_array)

    if asymmetry_type == GENICINTERGENICBIAS:
        genic_sims_replication_time_values = np.sum(genic_values, axis=(1))
        genic_sims_normalized_mutation_density = get_normalized_mutation_density_2D(genic_sims_replication_time_values, number_of_attributable_bases_array)

    return all_sims_normalized_mutation_density, \
        strand1_sims_normalized_mutation_density, \
        strand2_sims_normalized_mutation_density, \
        strand3_questionable_sims_normalized_mutation_density, \
        strand4_nontranscribed_sims_normalized_mutation_density, \
        genic_sims_normalized_mutation_density


def get_signature_specific_sims_normalized_mutation_density_from_np_array(asymmetry_type,
                                                       mutation_main_type,
                                                       signature_index,
                                                       number_of_attributable_bases_array,
                                                       number_of_sims,
                                                       simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                       simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                       simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                       simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                                                       simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                                                       simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array):

    all_sims_normalized_mutation_density = None
    strand1_sims_normalized_mutation_density = None
    strand2_sims_normalized_mutation_density = None
    strand3_questionable_sims_normalized_mutation_density = None
    strand4_nontranscribed_sims_normalized_mutation_density = None
    genic_sims_normalized_mutation_density = None

    if mutation_main_type == SBS:
        if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
            array_of_interest = simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array

        elif asymmetry_type == REPLICATIONSTRANDBIAS:
            array_of_interest = simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array

    elif mutation_main_type == DBS:
        if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
            array_of_interest = simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array

        elif asymmetry_type == REPLICATIONSTRANDBIAS:
            array_of_interest = simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array

    elif mutation_main_type == ID:
        if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
            array_of_interest = simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array

        elif asymmetry_type == REPLICATIONSTRANDBIAS:
            array_of_interest = simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array

    if asymmetry_type == TRANSCRIPTIONSTRANDBIAS:
        # strand1_value = 'T'
        # strand2_value = 'U'
        # strand3_value = 'Q'
        # strand4_value = 'N'
        # transcription_strand_mapping = {'T': 0, 'U': 1, 'N': 2, 'Q': 3}

        all_values = array_of_interest[1:number_of_sims+1, signature_index, :, :, :]

        # strand1 = TRANSCRIBED
        strand1_values = array_of_interest[1:number_of_sims+1, signature_index, :, 0, :]

        # strand2 = UNTRANSCRIBED
        strand2_values = array_of_interest[1:number_of_sims+1, signature_index, :, 1, :]

        if mutation_main_type == DBS or mutation_main_type == ID:
            # strand3 = QUESTIONABLE
            strand3_values = array_of_interest[1:number_of_sims+1, signature_index, :, 3, :]

        # strand4 = NONTRANSCRIBED
        strand4_values = array_of_interest[1:number_of_sims+1, signature_index, :, 2, :]

    elif asymmetry_type == GENICINTERGENICBIAS:
        # strand1_value = 'T'
        # strand2_value = 'U'
        # strand3_value = 'Q'
        # strand4_value = 'N'
        # transcription_strand_mapping = {'T': 0, 'U': 1, 'N': 2, 'Q': 3}

        all_values = array_of_interest[1:number_of_sims+1, signature_index, :, :, :]

        # strand1 = TRANSCRIBED
        strand1_values = array_of_interest[1:number_of_sims+1, signature_index, :, 0, :]

        # strand2 = UNTRANSCRIBED
        strand2_values = array_of_interest[1:number_of_sims+1, signature_index, :, 1, :]

        if mutation_main_type == DBS or mutation_main_type == ID:
            # strand3 = QUESTIONABLE
            strand3_values = array_of_interest[1:number_of_sims+1, signature_index, :, 3, :]

        # strand4 = NONTRANSCRIBED
        strand4_values = array_of_interest[1:number_of_sims+1, signature_index, :, 2, :]

        # TODO does questionable mean either transcribed or untranscribed?
        # TODO can we consider questionable values in genic values?
        genic_values =   strand1_values + strand2_values
        # intergenic_values =  strand4_values

    elif asymmetry_type == REPLICATIONSTRANDBIAS:
        # strand1_value = 'E'
        # strand2_value = 'A'
        # strand3_value = 'Q'
        # replication_strand_mapping = {'E': 0, 'A': 1, 'Q': 2}  # Map all possible strands

        all_values = array_of_interest[1:number_of_sims+1, signature_index, :, :, :]

        # strand1 = LEADING
        strand1_values = array_of_interest[1:number_of_sims+1, signature_index, :, 0, :]

        # strand2 = LAGGING
        strand2_values = array_of_interest[1:number_of_sims+1, signature_index, :, 1, :]

        if mutation_main_type == DBS or mutation_main_type == ID:
            # strand3 = QUESTIONABLE
            strand3_values = array_of_interest[1:number_of_sims+1, signature_index, :, 2, :]

    # Accumulate (sum) values across all mutation types and strands
    # Here we sum across the mutation type axis (axis=0) and strand axis (axis=1)

    # number_of_sims, mutation_types, num_of_strands, replication_times
    all_sims_replication_time_values = np.sum(all_values, axis=(1, 2))
    all_sims_normalized_mutation_density = get_normalized_mutation_density_2D(all_sims_replication_time_values, number_of_attributable_bases_array)

    # num_of_sims, mutation_types, replication_times
    strand1_sims_replication_time_values = np.sum(strand1_values, axis=(1))
    strand1_sims_normalized_mutation_density = get_normalized_mutation_density_2D(strand1_sims_replication_time_values, number_of_attributable_bases_array)

    # num_of_sims, mutation_types, replication_times
    strand2_sims_replication_time_values = np.sum(strand2_values, axis=(1))
    strand2_sims_normalized_mutation_density = get_normalized_mutation_density_2D(strand2_sims_replication_time_values, number_of_attributable_bases_array)

    if mutation_main_type == DBS or mutation_main_type == ID:
        strand3_sims_replication_time_values = np.sum(strand3_values, axis=(1))
        strand3_questionable_sims_normalized_mutation_density = get_normalized_mutation_density_2D(strand3_sims_replication_time_values, number_of_attributable_bases_array)

    if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
        strand4_sims_replication_time_values = np.sum(strand4_values, axis=(1))
        strand4_nontranscribed_sims_normalized_mutation_density = get_normalized_mutation_density_2D(strand4_sims_replication_time_values, number_of_attributable_bases_array)

    if asymmetry_type == GENICINTERGENICBIAS:
        genic_sims_replication_time_values = np.sum(genic_values, axis=(1))
        genic_sims_normalized_mutation_density = get_normalized_mutation_density_2D(genic_sims_replication_time_values, number_of_attributable_bases_array)

    return all_sims_normalized_mutation_density, \
        strand1_sims_normalized_mutation_density, \
        strand2_sims_normalized_mutation_density, \
        strand3_questionable_sims_normalized_mutation_density, \
        strand4_nontranscribed_sims_normalized_mutation_density, \
        genic_sims_normalized_mutation_density


def get_sims_normalized_mutation_density_from_np_array(asymmetry_type,
                                                       mutation_main_type,
                                                       number_of_attributable_bases_array,
                                                       number_of_sims,
                                                       simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                       simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                       simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                       simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                                                       simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                                                       simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array):

    all_sims_normalized_mutation_density = None
    strand1_sims_normalized_mutation_density = None
    strand2_sims_normalized_mutation_density = None
    strand3_questionable_sims_normalized_mutation_density = None
    strand4_nontranscribed_sims_normalized_mutation_density = None
    genic_sims_normalized_mutation_density = None

    if mutation_main_type == SBS:
        if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
            array_of_interest = simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array

        elif asymmetry_type == REPLICATIONSTRANDBIAS:
            array_of_interest = simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array

    elif mutation_main_type == DBS:
        if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
            array_of_interest = simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array

        elif asymmetry_type == REPLICATIONSTRANDBIAS:
            array_of_interest = simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array

    elif mutation_main_type == ID:
        if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
            array_of_interest = simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array

        elif asymmetry_type == REPLICATIONSTRANDBIAS:
            array_of_interest = simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array

    if asymmetry_type == TRANSCRIPTIONSTRANDBIAS:
        # strand1_value = 'T'
        # strand2_value = 'U'
        # strand3_value = 'Q'
        # strand4_value = 'N'
        # transcription_strand_mapping = {'T': 0, 'U': 1, 'N': 2, 'Q': 3}

        all_values = array_of_interest[1:number_of_sims+1, -1, :, :, :]

        # strand1 = TRANSCRIBED
        strand1_values = array_of_interest[1:number_of_sims+1, -1, :, 0, :]

        # strand2 = UNTRANSCRIBED
        strand2_values = array_of_interest[1:number_of_sims+1, -1, :, 1, :]

        if mutation_main_type == DBS or mutation_main_type == ID:
            # strand3 = QUESTIONABLE
            strand3_values = array_of_interest[1:number_of_sims+1, -1, :, 3, :]

        # strand4 = NONTRANSCRIBED
        strand4_values = array_of_interest[1:number_of_sims+1, -1, :, 2, :]

    elif asymmetry_type == GENICINTERGENICBIAS:
        # strand1_value = 'T'
        # strand2_value = 'U'
        # strand3_value = 'Q'
        # strand4_value = 'N'
        # transcription_strand_mapping = {'T': 0, 'U': 1, 'N': 2, 'Q': 3}

        all_values = array_of_interest[1:number_of_sims+1, -1, :, :, :]

        # strand1 = TRANSCRIBED
        strand1_values = array_of_interest[1:number_of_sims+1, -1, :, 0, :]

        # strand2 = UNTRANSCRIBED
        strand2_values = array_of_interest[1:number_of_sims+1, -1, :, 1, :]

        if mutation_main_type == DBS or mutation_main_type == ID:
            # strand3 = QUESTIONABLE
            strand3_values = array_of_interest[1:number_of_sims+1, -1, :, 3, :]

        # strand4 = NONTRANSCRIBED
        strand4_values = array_of_interest[1:number_of_sims+1, -1, :, 2, :]

        # TODO does questionable mean either transcribed or untranscribed?
        # TODO can we consider questionable values in genic values?
        genic_values =   strand1_values + strand2_values
        # intergenic_values =  strand4_values

    elif asymmetry_type == REPLICATIONSTRANDBIAS:
        # strand1_value = 'E'
        # strand2_value = 'A'
        # strand3_value = 'Q'
        # replication_strand_mapping = {'E': 0, 'A': 1, 'Q': 2}  # Map all possible strands

        all_values = array_of_interest[1:number_of_sims+1, -1, :, :, :]

        # strand1 = LEADING
        strand1_values = array_of_interest[1:number_of_sims+1, -1, :, 0, :]

        # strand2 = LAGGING
        strand2_values = array_of_interest[1:number_of_sims+1, -1, :, 1, :]

        if mutation_main_type == DBS or mutation_main_type == ID:
            # strand3 = QUESTIONABLE
            strand3_values = array_of_interest[1:number_of_sims+1, -1, :, 2, :]

    # Accumulate (sum) values across all mutation types and strands
    # Here we sum across the mutation type axis (axis=0) and strand axis (axis=1)

    # number_of_sims, mutation_types, num_of_strands, replication_times
    all_sims_replication_time_values = np.sum(all_values, axis=(1, 2))
    all_sims_normalized_mutation_density = get_normalized_mutation_density_2D(all_sims_replication_time_values, number_of_attributable_bases_array)

    # num_of_sims, mutation_types, replication_times
    strand1_sims_replication_time_values = np.sum(strand1_values, axis=(1))
    strand1_sims_normalized_mutation_density = get_normalized_mutation_density_2D(strand1_sims_replication_time_values, number_of_attributable_bases_array)

    # num_of_sims, mutation_types, replication_times
    strand2_sims_replication_time_values = np.sum(strand2_values, axis=(1))
    strand2_sims_normalized_mutation_density = get_normalized_mutation_density_2D(strand2_sims_replication_time_values, number_of_attributable_bases_array)

    if mutation_main_type == DBS or mutation_main_type == ID:
        strand3_sims_replication_time_values = np.sum(strand3_values, axis=(1))
        strand3_questionable_sims_normalized_mutation_density = get_normalized_mutation_density_2D(strand3_sims_replication_time_values, number_of_attributable_bases_array)

    if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
        strand4_sims_replication_time_values = np.sum(strand4_values, axis=(1))
        strand4_nontranscribed_sims_normalized_mutation_density = get_normalized_mutation_density_2D(strand4_sims_replication_time_values, number_of_attributable_bases_array)

    if asymmetry_type == GENICINTERGENICBIAS:
        genic_sims_replication_time_values = np.sum(genic_values, axis=(1))
        genic_sims_normalized_mutation_density = get_normalized_mutation_density_2D(genic_sims_replication_time_values, number_of_attributable_bases_array)

    return all_sims_normalized_mutation_density, \
        strand1_sims_normalized_mutation_density, \
        strand2_sims_normalized_mutation_density, \
        strand3_questionable_sims_normalized_mutation_density, \
        strand4_nontranscribed_sims_normalized_mutation_density, \
        genic_sims_normalized_mutation_density



def read_chrom_based_files(outputDir, jobname, number_of_sims, chrNamesList, mutation_main_type, column_of_interest):
    df_list = []
    sims_df_list = []

    for sim_num in range(0, number_of_sims + 1):
        for chrLong in chrNamesList:
            # e.g.: chr15_SBS_for_topography.txt
            filename = chrLong + '_' + mutation_main_type + '_for_topography.txt'

            if sim_num == 0:
                filepath = os.path.join(outputDir, jobname, DATA, 'chrbased', filename)
                if os.path.exists(filepath):
                    df = pd.read_csv(filepath, sep='\t', header=0)
                    df = df.groupby(['MutationLong', column_of_interest, 'ReplicationTime']).size().reset_index(
                        name='MutationCount')
                    df_list.append(df)

            else:
                filepath = os.path.join(outputDir, jobname, DATA, 'chrbased', 'sim' + str(sim_num), filename)
                if os.path.exists(filepath):
                    df = pd.read_csv(filepath, sep='\t', header=0)
                    df = df.groupby(['MutationLong', column_of_interest, 'ReplicationTime']).size().reset_index(
                        name='MutationCount')
                    df['Simulation'] = sim_num
                    sims_df_list.append(df)

    return df_list, sims_df_list

def get_strand_mutation_count_values(asymmetry_type,
                   mutation_main_type,
                   simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                   simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                   simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                   simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                   simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                   simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array):

    if mutation_main_type == SBS:
        if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
            array_of_interest = simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array

        elif asymmetry_type == REPLICATIONSTRANDBIAS:
            array_of_interest = simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array

    elif mutation_main_type == DBS:
        if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
            array_of_interest = simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array

        elif asymmetry_type == REPLICATIONSTRANDBIAS:
            array_of_interest = simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array

    elif mutation_main_type == ID:
        if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == GENICINTERGENICBIAS:
            array_of_interest = simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array

        elif asymmetry_type == REPLICATIONSTRANDBIAS:
            array_of_interest = simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array

    if asymmetry_type == TRANSCRIPTIONSTRANDBIAS:
        # strand1_value = 'T'
        # strand2_value = 'U'
        # strand3_value = 'Q'
        # strand4_value = 'N'
        # transcription_strand_mapping = {'T': 0, 'U': 1, 'N': 2, 'Q': 3, 'B':4}

        # strand1 = TRANSCRIBED
        strand1_values = array_of_interest[0, -1, :, 0, :]

        # strand2 = UNTRANSCRIBED
        strand2_values = array_of_interest[0, -1, :, 1, :]

        # strand3 = QUESTIONABLE
        strand3_values = array_of_interest[0, -1, :, 3, :]

        # strand4 = NONTRANSCRIBED
        strand4_values = array_of_interest[0, -1, :, 2, :]

    elif asymmetry_type == GENICINTERGENICBIAS:
        # strand1_value = 'T'
        # strand2_value = 'U'
        # strand3_value = 'Q'
        # strand4_value = 'N'
        # transcription_strand_mapping = {'T': 0, 'U': 1, 'N': 2, 'Q': 3, 'B':4}

        # strand1 = TRANSCRIBED
        strand1_values = array_of_interest[0, -1, :, 0, :]

        # strand2 = UNTRANSCRIBED
        strand2_values = array_of_interest[0, -1, :, 1, :]

        # strand3 = QUESTIONABLE
        strand3_values = array_of_interest[0, -1, :, 3, :]

        # strand4 = NONTRANSCRIBED
        strand4_values = array_of_interest[0, -1, :, 2, :]

        # TODO does questionable mean either transcribed or untranscribed?
        # TODO can we consider questionable values in genic values?
        genic_values = strand1_values + strand2_values
        intergenic_values =  strand4_values

    elif asymmetry_type == REPLICATIONSTRANDBIAS:
        # strand1_value = 'E'
        # strand2_value = 'A'
        # strand3_value = 'Q'
        # replication_strand_mapping = {'E': 0, 'A': 1, 'Q': 2, 'U': -1, 'B': 4}  # Map all possible strands

        # strand1 = LEADING
        strand1_values = array_of_interest[0, -1, :, 0, :]

        # strand2 = LAGGING
        strand2_values = array_of_interest[0, -1, :, 1, :]

        # strand3 = QUESTIONABLE
        strand3_values = array_of_interest[0, -1, :, 2, :]

    # Accumulate (sum) values across all mutation types and strands
    # Here we sum across the mutation type axis (axis=0) and strand axis (axis=1)
    strand1_mutationtypes_mutationcounts_np_array = np.sum(strand1_values, axis=(1))
    strand2_mutationtypes_mutationcounts_np_array = np.sum(strand2_values, axis=(1))
    strand3_questionable_mutationtypes_mutationcounts_np_array = np.sum(strand3_values, axis=(1))

    if asymmetry_type == GENICINTERGENICBIAS:
        genic_mutationtypes_mutationcounts_np_array = np.sum(genic_values, axis=(1))
        intergenic_mutationtypes_mutationcounts_np_array = np.sum(intergenic_values, axis=(1))

    if asymmetry_type == TRANSCRIPTIONSTRANDBIAS or asymmetry_type == REPLICATIONSTRANDBIAS:
        return strand1_mutationtypes_mutationcounts_np_array, \
            strand2_mutationtypes_mutationcounts_np_array, \
            strand3_questionable_mutationtypes_mutationcounts_np_array

    elif asymmetry_type == GENICINTERGENICBIAS:
        return genic_mutationtypes_mutationcounts_np_array, \
            intergenic_mutationtypes_mutationcounts_np_array, \
            strand3_questionable_mutationtypes_mutationcounts_np_array



def plot_strand_asymmetry_vs_replication_timing_figures(inputDir,
                                                        outputDir,
                                                        jobname,
                                                        number_of_sims,
                                                        asymmetry_types,
                                                        mutation_types,
                                                        ordered_SBS96_substitution_types,
                                                        ordered_DBS78_doublet_types,
                                                        ordered_ID83_indel_types,
                                                        SBS96_mutation_mapping,
                                                        DBS78_mutation_mapping,
                                                        ID83_mutation_mapping,
                                                        ordered_sbs_signatures_with_cutoffs,
                                                        ordered_dbs_signatures_with_cutoffs,
                                                        ordered_id_signatures_with_cutoffs,
                                                        ordered_sbs_signatures_cutoffs,
                                                        ordered_dbs_signatures_cutoffs,
                                                        ordered_id_signatures_cutoffs,
                                                        simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                        simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                        simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                        simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                                                        simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                                                        simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array):

    os.makedirs(os.path.join(outputDir, jobname, FIGURE, MULTILEVEL), exist_ok=True)

    strand1 = strand2 = strand3 = strand4 = None
    strand1_value = strand2_value = strand3_value = strand4_value = None

    # Get the number of attributable bases
    number_of_attributable_bases_file_name = 'NumberofAttributableBases.txt'
    number_of_attributable_bases_file_path = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME,
                                                          number_of_attributable_bases_file_name)
    number_of_attributable_bases_array = np.loadtxt(number_of_attributable_bases_file_path, dtype=int)

    # Create a complete index of ReplicationTime from 1 to 10
    complete_index = pd.Series(range(1, 11))

    for asymmetry_type in asymmetry_types:
        for mutation_main_type in mutation_types:

            # Set mutation subtypes
            # These mutation subtypes are shown in the mutation profile plots at the top
            if mutation_main_type == SBS:
                mutation_sub_types = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
                ordered_signatures_with_cutoffs = ordered_sbs_signatures_with_cutoffs

            elif mutation_main_type == DBS:
                mutation_sub_types = ['AC>NN', 'AT>NN', 'CC>NN', 'CG>NN', 'CT>NN', 'GC>NN', 'TA>NN', 'TC>NN', 'TG>NN',
                                      'TT>NN']
                ordered_signatures_with_cutoffs = ordered_dbs_signatures_with_cutoffs

            elif mutation_main_type == ID:
                mutation_sub_types = ['1bp Deletion', '1bp Insertion', '>1bp Deletion', '>1bp Insertion',
                                      'Microhomology']
                ordered_signatures_with_cutoffs = ordered_id_signatures_with_cutoffs

            if asymmetry_type == TRANSCRIPTIONSTRANDBIAS:
                column_of_interest = TRANSCRIPTIONSTRAND
                filename = mutation_main_type + '_TSA_versus_RT.pdf'
                strand1 = TRANSCRIBED
                strand2 = UNTRANSCRIBED
                strand3 = QUESTIONABLE
                strand4 = NONTRANSCRIBED
                strand1_value = 'T'
                strand2_value = 'U'
                strand3_value = 'Q'
                strand4_value = 'N'
                # Q: Questionable. This category is used to classify any mutations that are a mix of purines and pyrimidines and
                # thus can't be classified into one of the above 4 categories.

            elif asymmetry_type == GENICINTERGENICBIAS:
                column_of_interest = TRANSCRIPTIONSTRAND
                filename = mutation_main_type + '_GIGA_versus_RT.pdf'
                strand1 = TRANSCRIBED
                strand2 = UNTRANSCRIBED
                strand3 = QUESTIONABLE
                strand4 = NONTRANSCRIBED
                strand1_value = 'T'
                strand2_value = 'U'
                strand3_value = 'Q'
                strand4_value = 'N'
                # Q: Questionable. This category is used to classify any mutations that are a mix of purines and pyrimidines and
                # thus can't be classified into one of the above 4 categories.

            elif asymmetry_type == REPLICATIONSTRANDBIAS:
                column_of_interest = REPLICATIONSTRAND
                filename = mutation_main_type + '_RSA_versus_RT.pdf'
                strand1 = LEADING
                strand2 = LAGGING
                strand3 = QUESTIONABLE
                strand1_value = 'E'
                strand2_value = 'A'
                strand3_value = 'Q'
                # Q: Questionable. This category is used to classify any mutations that are a mix of purines and pyrimidines and
                # thus can't be classified into one of the above 4 categories.

            with PdfPages(os.path.join(outputDir, jobname, FIGURE, MULTILEVEL, filename)) as pdf:

                # Plot mutational profile and mutational profile with strand asymmetry for real mutations only using pivoted_df
                # pivoted_df columns e.g.: Mutation Leading  Lagging  Questionable
                # pivoted_df rows: SBS96, DBS78 or ID83 mutation types
                # pivoted_df values: MutationCount
                # pivoted_df is based on mutation_long_combined_df and mutation_long_df
                # mutation_long is converted into SBS96, DBS78 or ID83 mutation types

                strand1_values, \
                strand2_values, \
                strand3_values = get_strand_mutation_count_values(asymmetry_type,
                                mutation_main_type,
                                simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                                simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                                simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array)

                if mutation_main_type == SBS:
                    matrix_path = os.path.join(inputDir, 'output', mutation_main_type, jobname + '.SBS96.all')
                    alternative_matrix_path = os.path.join(inputDir, 'output', mutation_main_type,
                                                           jobname + '.all_samples.SBS96.all')
                    spmg_df = pd.read_csv(matrix_path, sep='\t')
                    spmg_df[jobname] = spmg_df.sum(axis=1)
                    spmg_df[['MutationType', jobname]].to_csv(alternative_matrix_path, sep='\t', index=False)
                    output_path = os.path.join(outputDir, jobname, FIGURE, MULTILEVEL, os.sep)
                    # original spp mutational profile plot
                    spplt.plotSBS(alternative_matrix_path,
                                  output_path,
                                  jobname,
                                  "96",
                                  percentage=True,
                                  savefig_format='png')
                    pdf.savefig(bbox_inches='tight')  # Save the current figure into the PDF

                    plotSBS_modified_from_spp(alternative_matrix_path,
                                              output_path,
                                              jobname,
                                              "96",
                                              # pivoted_df,
                                              strand1_values,
                                              strand2_values,
                                              asymmetry_type,
                                              percentage=True,
                                              savefig_format='png')
                    pdf.savefig(bbox_inches='tight')  # Save the current figure into the PDF

                elif mutation_main_type == DBS:
                    matrix_path = os.path.join(inputDir, 'output', mutation_main_type, jobname + '.DBS78.all')
                    alternative_matrix_path = os.path.join(inputDir, 'output', mutation_main_type,
                                                           jobname + '.all_samples.DBS78.all')
                    spmg_df = pd.read_csv(matrix_path, sep='\t')
                    spmg_df[jobname] = spmg_df.sum(axis=1)
                    spmg_df[['MutationType', jobname]].to_csv(alternative_matrix_path, sep='\t', index=False)
                    output_path = os.path.join(outputDir, jobname, FIGURE, MULTILEVEL, os.sep)
                    # original spp mutational profile plot
                    spplt.plotDBS(alternative_matrix_path, output_path, jobname, "78", percentage=True,
                                  savefig_format='png')
                    pdf.savefig(bbox_inches='tight')  # Save the current figure into the PDF

                    plotDBS_modified_from_spp(alternative_matrix_path,
                                              output_path,
                                              jobname,
                                              "78",
                                              # pivoted_df,
                                              strand1_values,
                                              strand2_values,
                                              strand3_values,
                                              asymmetry_type,
                                              percentage=True,
                                              savefig_format='png')
                    pdf.savefig(bbox_inches='tight')  # Save the current figure into the PDF

                elif mutation_main_type == ID:
                    matrix_path = os.path.join(inputDir, 'output', mutation_main_type, jobname + '.ID83.all')
                    alternative_matrix_path = os.path.join(inputDir, 'output', mutation_main_type,
                                                           jobname + '.all_samples.ID83.all')
                    spmg_df = pd.read_csv(matrix_path, sep='\t')

                    spmg_df[jobname] = spmg_df.sum(axis=1)
                    spmg_df[['MutationType', jobname]].to_csv(alternative_matrix_path, sep='\t', index=False)
                    output_path = os.path.join(outputDir, jobname, FIGURE, MULTILEVEL, os.sep)
                    # original spp mutational profile plot
                    spplt.plotID(alternative_matrix_path,
                                 output_path,
                                 jobname,
                                 "83",
                                 percentage=True,
                                 savefig_format='png')
                    pdf.savefig(bbox_inches='tight')  # Save the current figure into the PDF

                    plotID_modified_from_spp(alternative_matrix_path,
                                             output_path,
                                             jobname,
                                             "83",
                                             # pivoted_df,
                                             strand1_values,
                                             strand2_values,
                                             strand3_values,
                                             asymmetry_type,
                                             percentage=True,
                                             savefig_format='png')
                    pdf.savefig(bbox_inches='tight')  # Save the current figure into the PDF

                # # For real mutations
                all_normalized_mutation_density, \
                strand1_normalized_mutation_density,\
                strand2_normalized_mutation_density, \
                strand3_questionable_normalized_mutation_density, \
                strand4_nontranscribed_normalized_mutation_density, \
                genic_normalized_mutation_density = get_normalized_mutation_density_from_np_array(asymmetry_type,
                                                                                                mutation_main_type,
                                                                                                number_of_attributable_bases_array,
                                                                                                simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                                                                simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                                                                simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                                                                simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                                                                                                simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                                                                                                simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array)


                # For simulated mutations
                all_sims_normalized_mutation_density, \
                strand1_sims_normalized_mutation_density, \
                strand2_sims_normalized_mutation_density, \
                strand3_questionable_sims_normalized_mutation_density, \
                strand4_nontranscribed_sims_normalized_mutation_density, \
                genic_sims_normalized_mutation_density = get_sims_normalized_mutation_density_from_np_array(asymmetry_type,
                                                                                                                 mutation_main_type,
                                                                                                                 number_of_attributable_bases_array,
                                                                                                                 number_of_sims,
                                                                                                                 simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                                                                                 simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                                                                                 simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                                                                                 simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                                                                                                                 simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                                                                                                                 simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array)

                plot_in_pdf(outputDir,
                            jobname,
                            mutation_main_type,
                            None,
                            asymmetry_type,
                            strand1,
                            strand2,
                            strand4,
                            all_normalized_mutation_density,
                            all_sims_normalized_mutation_density,
                            strand1_normalized_mutation_density,
                            strand1_sims_normalized_mutation_density,
                            strand2_normalized_mutation_density,
                            strand2_sims_normalized_mutation_density,
                            strand3_questionable_normalized_mutation_density,
                            strand3_questionable_sims_normalized_mutation_density,
                            strand4_nontranscribed_normalized_mutation_density,
                            strand4_nontranscribed_sims_normalized_mutation_density,
                            genic_normalized_mutation_density,
                            genic_sims_normalized_mutation_density,
                            pdf)

                # Signature specific
                for signature_index, signature in  enumerate(ordered_signatures_with_cutoffs):

                    all_normalized_mutation_density, \
                        strand1_normalized_mutation_density, \
                        strand2_normalized_mutation_density, \
                        strand3_questionable_normalized_mutation_density, \
                        strand4_nontranscribed_normalized_mutation_density, \
                        genic_normalized_mutation_density = get_signature_specific_normalized_mutation_density_from_np_array(
                        asymmetry_type,
                        mutation_main_type,
                        signature_index,
                        number_of_attributable_bases_array,
                        simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                        simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                        simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                        simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                        simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                        simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array)

                    # For simulated mutations
                    all_sims_normalized_mutation_density, \
                        strand1_sims_normalized_mutation_density, \
                        strand2_sims_normalized_mutation_density, \
                        strand3_questionable_sims_normalized_mutation_density, \
                        strand4_nontranscribed_sims_normalized_mutation_density, \
                        genic_sims_normalized_mutation_density = get_signature_specific_sims_normalized_mutation_density_from_np_array(
                        asymmetry_type,
                        mutation_main_type,
                        signature_index,
                        number_of_attributable_bases_array,
                        number_of_sims,
                        simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                        simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                        simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                        simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                        simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                        simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array)

                    plot_in_pdf(outputDir,
                                jobname,
                                mutation_main_type,
                                signature,
                                asymmetry_type,
                                strand1,
                                strand2,
                                strand4,
                                all_normalized_mutation_density,
                                all_sims_normalized_mutation_density,
                                strand1_normalized_mutation_density,
                                strand1_sims_normalized_mutation_density,
                                strand2_normalized_mutation_density,
                                strand2_sims_normalized_mutation_density,
                                strand3_questionable_normalized_mutation_density,
                                strand3_questionable_sims_normalized_mutation_density,
                                strand4_nontranscribed_normalized_mutation_density,
                                strand4_nontranscribed_sims_normalized_mutation_density,
                                genic_normalized_mutation_density,
                                genic_sims_normalized_mutation_density,
                                pdf)

                # Mutation subtype specific
                for mutation_sub_type in mutation_sub_types:

                    # For real mutations
                    all_normalized_mutation_density, \
                    strand1_normalized_mutation_density, \
                    strand2_normalized_mutation_density, \
                    strand3_questionable_normalized_mutation_density, \
                    strand4_nontranscribed_normalized_mutation_density, \
                    genic_normalized_mutation_density = get_mutation_specific_normalized_mutation_density_from_np_array(
                        asymmetry_type,
                        mutation_main_type,
                        mutation_sub_type,
                        ordered_SBS96_substitution_types,
                        ordered_DBS78_doublet_types,
                        ordered_ID83_indel_types,
                        SBS96_mutation_mapping,
                        DBS78_mutation_mapping,
                        ID83_mutation_mapping,
                        number_of_attributable_bases_array,
                        simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                        simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                        simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                        simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                        simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                        simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array)

                    # For simulated mutations
                    all_sims_normalized_mutation_density, \
                    strand1_sims_normalized_mutation_density, \
                    strand2_sims_normalized_mutation_density, \
                    strand3_questionable_sims_normalized_mutation_density, \
                    strand4_nontranscribed_sims_normalized_mutation_density, \
                    genic_sims_normalized_mutation_density = get_mutation_specific_sims_normalized_mutation_density_from_np_array(
                        asymmetry_type,
                        mutation_main_type,
                        mutation_sub_type,
                        ordered_SBS96_substitution_types,
                        ordered_DBS78_doublet_types,
                        ordered_ID83_indel_types,
                        SBS96_mutation_mapping,
                        DBS78_mutation_mapping,
                        ID83_mutation_mapping,
                        number_of_attributable_bases_array,
                        number_of_sims,
                        simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                        simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                        simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                        simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                        simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                        simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array)

                    plot_in_pdf(outputDir,
                                jobname,
                                mutation_main_type,
                                mutation_sub_type,
                                asymmetry_type,
                                strand1,
                                strand2,
                                strand4,
                                all_normalized_mutation_density,
                                all_sims_normalized_mutation_density,
                                strand1_normalized_mutation_density,
                                strand1_sims_normalized_mutation_density,
                                strand2_normalized_mutation_density,
                                strand2_sims_normalized_mutation_density,
                                strand3_questionable_normalized_mutation_density,
                                strand3_questionable_sims_normalized_mutation_density,
                                strand4_nontranscribed_normalized_mutation_density,
                                strand4_nontranscribed_sims_normalized_mutation_density,
                                genic_normalized_mutation_density,
                                genic_sims_normalized_mutation_density,
                                pdf)


def nested_analyses_plot_strand_asymmetry_vs_replication_timing_figures_using_mp(inputDir,
                                                                                 outputDir,
                                                                                 jobname,
                                                                                 number_of_sims,
                                                                                 chromNamesList,
                                                                                 mutation_types,
                                                                                 asymmetry_types,
                                                                                 ordered_sbs_signatures_with_cutoffs,
                                                                                 ordered_dbs_signatures_with_cutoffs,
                                                                                 ordered_id_signatures_with_cutoffs,
                                                                                 ordered_sbs_signatures_cutoffs,
                                                                                 ordered_dbs_signatures_cutoffs,
                                                                                 ordered_id_signatures_cutoffs,
                                                                                 parallel_mode):

    transcription_strand_mapping = {'T': 0, 'U': 1, 'N': 2, 'Q': 3, 'B':4}  # Map all possible strands

    replication_strand_mapping = {'E': 0, 'A': 1, 'Q': 2, 'U': -1, 'B':4}  # Map all possible strands
    # mapping = {0: 'A', 1: 'E', 2: 'B', 3: 'X', -1: 'U'} From Replication Strand Asymmetry

    replication_time_mapping = {time: i for i, time in enumerate(range(1, 11))}

    SBS96_mutation_mapping = {mutation: i for i, mutation in
                              enumerate(ordered_SBS96_substitution_types)}  # Map all possible mutations

    DBS78_mutation_mapping = {mutation: i for i, mutation in
                              enumerate(ordered_DBS78_doublet_types)}  # Map all possible mutations

    ID83_mutation_mapping = {mutation: i for i, mutation in
                              enumerate(ordered_ID83_indel_types)}  # Map all possible mutations

    # Signature can be SBS signature or all
    # Mutation type will be SBS96 mutation types
    # Strand can be transcribed, untranscribed, nontranscribed, questionable
    # Replication time can be 1 to 10
    simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array = np.zeros((number_of_sims+1,
                                                                                                          np.size(ordered_sbs_signatures_with_cutoffs)+1,
                                                                                                          np.size(ordered_SBS96_substitution_types),
                                                                                                          5, # transcribed, untranscribed, nontranscribed, bidirectional, questionable
                                                                                                          10)) # Replication time

    simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array = np.zeros((number_of_sims+1,
                                                                                                          np.size(ordered_dbs_signatures_with_cutoffs)+1,
                                                                                                          np.size(ordered_DBS78_doublet_types),
                                                                                                          5, # transcribed, untranscribed, nontranscribed, bidirectional, questionable
                                                                                                          10)) # Replication time

    simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array = np.zeros((number_of_sims+1,
                                                                                                          np.size(ordered_id_signatures_with_cutoffs)+1,
                                                                                                          np.size(ordered_ID83_indel_types),
                                                                                                          5, # transcribed, untranscribed, nontranscribed, bidirectional, questionable
                                                                                                          10)) # Replication time

    simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array = np.zeros((number_of_sims+1,
                                                                                                          np.size(ordered_sbs_signatures_with_cutoffs)+1,
                                                                                                          np.size(ordered_SBS96_substitution_types),
                                                                                                          5, # Leading, Lagging
                                                                                                          10)) # Replication time

    simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array = np.zeros((number_of_sims+1,
                                                                                                          np.size(ordered_dbs_signatures_with_cutoffs)+1,
                                                                                                          np.size(ordered_DBS78_doublet_types),
                                                                                                          5, # Leading, Lagging, Questionable
                                                                                                          10)) # Replication time

    simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array = np.zeros((number_of_sims+1,
                                                                                                          np.size(ordered_id_signatures_with_cutoffs)+1,
                                                                                                          np.size(ordered_ID83_indel_types),
                                                                                                          5, # Leading, Lagging, Questionable
                                                                                                          10)) # Replication time



    def accumulate_np_arrays(result_tuple):
        chrLong = result_tuple[0]
        sim_num = result_tuple[1]

        sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array = result_tuple[2]
        dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array = result_tuple[3]
        id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array = result_tuple[4]

        sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array = result_tuple[5]
        dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array = result_tuple[6]
        id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array = result_tuple[7]

        print('MONITOR ACCUMULATE sim_num:', sim_num,  chrLong, flush=True)

        # accumulate
        simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array[sim_num] += sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array
        simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array[sim_num] += dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array
        simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array[sim_num] += id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array

        simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array[sim_num] += sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array
        simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array[sim_num] += dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array
        simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array[sim_num] += id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array


    sim_nums = range(0, number_of_sims + 1)
    sim_num_chr_tuples = ((sim_num, chrLong) for sim_num in sim_nums for chrLong in chromNamesList)

    if parallel_mode:

        jobs = []

        num_of_processes = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=num_of_processes)

        for simNum, chrLong in sim_num_chr_tuples:
            jobs.append(pool.apply_async(fill_signature_mutationtype_strand_replicationtime_arrays,
                args=(outputDir,
                      jobname,
                      chrLong,
                      simNum,
                      asymmetry_types,
                      ordered_sbs_signatures_with_cutoffs,
                      ordered_dbs_signatures_with_cutoffs,
                      ordered_id_signatures_with_cutoffs,
                      ordered_sbs_signatures_cutoffs,
                      ordered_dbs_signatures_cutoffs,
                      ordered_id_signatures_cutoffs,
                      ordered_SBS96_substitution_types,
                      ordered_DBS78_doublet_types,
                      ordered_ID83_indel_types,
                      transcription_strand_mapping,
                      replication_strand_mapping,
                      replication_time_mapping,
                      SBS96_mutation_mapping,
                      DBS78_mutation_mapping,
                      ID83_mutation_mapping,),
                callback=accumulate_np_arrays))

        pool.close()
        pool.join()

    else:
        for simNum, chrLong in sim_num_chr_tuples:
            result_tuple = fill_signature_mutationtype_strand_replicationtime_arrays(
                outputDir,
                jobname,
                chrLong,
                simNum,
                asymmetry_types,
                ordered_sbs_signatures_with_cutoffs,
                ordered_dbs_signatures_with_cutoffs,
                ordered_id_signatures_with_cutoffs,
                ordered_sbs_signatures_cutoffs,
                ordered_dbs_signatures_cutoffs,
                ordered_id_signatures_cutoffs,
                ordered_SBS96_substitution_types,
                ordered_DBS78_doublet_types,
                ordered_ID83_indel_types,
                transcription_strand_mapping,
                replication_strand_mapping,
                replication_time_mapping,
                SBS96_mutation_mapping,
                DBS78_mutation_mapping,
                ID83_mutation_mapping)

            accumulate_np_arrays(result_tuple)


    plot_strand_asymmetry_vs_replication_timing_figures(inputDir,
                                                        outputDir,
                                                        jobname,
                                                        number_of_sims,
                                                        asymmetry_types,
                                                        mutation_types,
                                                        ordered_SBS96_substitution_types,
                                                        ordered_DBS78_doublet_types,
                                                        ordered_ID83_indel_types,
                                                        SBS96_mutation_mapping,
                                                        DBS78_mutation_mapping,
                                                        ID83_mutation_mapping,
                                                        ordered_sbs_signatures_with_cutoffs,
                                                        ordered_dbs_signatures_with_cutoffs,
                                                        ordered_id_signatures_with_cutoffs,
                                                        ordered_sbs_signatures_cutoffs,
                                                        ordered_dbs_signatures_cutoffs,
                                                        ordered_id_signatures_cutoffs,
                                                        simnum_sbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                        simnum_dbs_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                        simnum_id_signature_mutationtype_transctiptionstrand_replicationtime_accumulated_np_array,
                                                        simnum_sbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                                                        simnum_dbs_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array,
                                                        simnum_id_signature_mutationtype_replicationstrand_replicationtime_accumulated_np_array)



def replicationTimeNormalizedMutationDensityFigures(outputDir,
                                                    jobname,
                                                    numberofSimulations,
                                                    mutation_types,
                                                    plot_sample_based,
                                                    plot_mode):

    # Initialize these dataframes as empty dataframe
    # We will read these dataframes if there is the corresponding data
    subsSignature_cutoff_numberofmutations_averageprobability_df = pd.DataFrame()
    dinucsSignature_cutoff_numberofmutations_averageprobability_df = pd.DataFrame()
    indelsSignature_cutoff_numberofmutations_averageprobability_df = pd.DataFrame()

    subsSignature_cutoff_numberofmutations_averageprobability_path = os.path.join(outputDir, jobname, DATA, Table_SBS_Signature_Cutoff_NumberofMutations_AverageProbability_Filename)
    dinucsSignature_cutoff_numberofmutations_averageprobability_path = os.path.join(outputDir, jobname, DATA, Table_DBS_Signature_Cutoff_NumberofMutations_AverageProbability_Filename)
    indelsSignature_cutoff_numberofmutations_averageprobability_path = os.path.join(outputDir, jobname, DATA, Table_ID_Signature_Cutoff_NumberofMutations_AverageProbability_Filename)

    if (SBS in mutation_types) and os.path.exists(subsSignature_cutoff_numberofmutations_averageprobability_path):
        subsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(subsSignature_cutoff_numberofmutations_averageprobability_path, sep='\t', header=0,dtype={'cutoff': np.float32, 'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})

    if (DBS in mutation_types) and os.path.exists(dinucsSignature_cutoff_numberofmutations_averageprobability_path):
        dinucsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(dinucsSignature_cutoff_numberofmutations_averageprobability_path, sep='\t', header=0,dtype={'cutoff': np.float32, 'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})

    if (ID in mutation_types) and os.path.exists(indelsSignature_cutoff_numberofmutations_averageprobability_path):
        indelsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(indelsSignature_cutoff_numberofmutations_averageprobability_path,sep='\t', header=0,dtype={'cutoff': np.float32, 'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})

    if plot_sample_based:
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
    if SBS in mutation_types:
        plotAllMutationTypesFigures('Aggregated Substitutions', 'royalblue', 'lightblue', AGGREGATEDSUBSTITUTIONS, None,
                                outputDir, jobname, numberofSimulations, plot_sample_based, sample2NumberofSubsDict,
                                subsSignature_cutoff_numberofmutations_averageprobability_df,
                                sample2SubsSignature2NumberofMutationsDict, plot_mode)

    if DBS in mutation_types:
        plotAllMutationTypesFigures('Aggregated Doublets', 'crimson', 'lightpink', AGGREGATEDDINUCS, None,
                                outputDir, jobname, numberofSimulations, plot_sample_based, sample2NumberofDinucsDict,
                                dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                sample2DinucsSignature2NumberofMutationsDict, plot_mode)
    if ID in mutation_types:
        plotAllMutationTypesFigures('Aggregated Indels', 'yellowgreen', 'lightgreen', AGGREGATEDINDELS, None,
                                outputDir, jobname, numberofSimulations, plot_sample_based, sample2NumberofIndelsDict,
                                indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                sample2IndelsSignature2NumberofMutationsDict, plot_mode)

    plotAllMutationTypesFigures(MICROHOMOLOGY, 'yellowgreen', 'lightgreen', INDELBASED, MICROHOMOLOGY,
                                outputDir, jobname, numberofSimulations, plot_sample_based, sample2NumberofIndelsDict,
                                indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                sample2IndelsSignature2NumberofMutationsDict, plot_mode)

    plotAllMutationTypesFigures(REPEAT, 'yellowgreen', 'lightgreen', INDELBASED, REPEAT,
                                outputDir, jobname, numberofSimulations, plot_sample_based, sample2NumberofIndelsDict,
                                indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                sample2IndelsSignature2NumberofMutationsDict, plot_mode)

    if plot_mode == PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_MANUSCRIPT :
        plot_replication_time_legend(outputDir, jobname)

    if (SBS in mutation_types) and (not subsSignature_cutoff_numberofmutations_averageprobability_df.empty):
        plotSignatureFigures('royalblue',
                             'lightblue',
                             'Normalized\nsingle base substitution density',
                             SIGNATUREBASED,
                             outputDir,
                             jobname,
                             numberofSimulations,
                             plot_sample_based,
                             sample2NumberofSubsDict,
                             subsSignature_cutoff_numberofmutations_averageprobability_df,
                             sample2SubsSignature2NumberofMutationsDict,
                             plot_mode)
    if (DBS in mutation_types) and (not dinucsSignature_cutoff_numberofmutations_averageprobability_df.empty):
        plotSignatureFigures('crimson',
                             'lightpink',
                             'Normalized\ndoublet base substitution density',
                             SIGNATUREBASED,
                             outputDir,
                             jobname,
                             numberofSimulations,
                             plot_sample_based,
                             sample2NumberofDinucsDict,
                             dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                             sample2DinucsSignature2NumberofMutationsDict,
                             plot_mode)
    if (ID in mutation_types) and (not indelsSignature_cutoff_numberofmutations_averageprobability_df.empty):
        plotSignatureFigures('yellowgreen',
                             'lightgreen',
                             'Normalized\nsmall insertion and deletion density',
                             SIGNATUREBASED,
                             outputDir,
                             jobname,
                             numberofSimulations,
                             plot_sample_based,
                             sample2NumberofIndelsDict,
                             indelsSignature_cutoff_numberofmutations_averageprobability_df,
                             sample2IndelsSignature2NumberofMutationsDict,
                             plot_mode)
    ##########################################################################################
    ##########################  Plot figures ends  ###########################################
    ##########################################################################################




def make_pickle_file_modified_from_spp(context="SBS96", return_plot_template=False, volume=None):

    if context == "SBS96":
        # newly added
        # Set default background colors
        plt.rcParams['figure.facecolor'] = 'white'  # For figure background
        plt.rcParams['axes.facecolor'] = 'white'  # For axes background

        plot_custom_text = False
        sig_probs = False
        pcawg = False

        # total_count = sum(sum(nuc.values()) for nuc in mutations[sample].values())
        # , extent=[-5, 80, -5, 30])
        plt.rcParams["axes.linewidth"] = 2
        plot1 = plt.figure(figsize=(43.93, 9.92))
        plt.rc("axes", edgecolor="lightgray")
        panel1 = plt.axes([0.04, 0.09, 0.95, 0.77])
        seq96 = [
            "A[C>A]A",
            "A[C>A]C",
            "A[C>A]G",
            "A[C>A]T",
            "C[C>A]A",
            "C[C>A]C",
            "C[C>A]G",
            "C[C>A]T",
            "G[C>A]A",
            "G[C>A]C",
            "G[C>A]G",
            "G[C>A]T",
            "T[C>A]A",
            "T[C>A]C",
            "T[C>A]G",
            "T[C>A]T",
            "A[C>G]A",
            "A[C>G]C",
            "A[C>G]G",
            "A[C>G]T",
            "C[C>G]A",
            "C[C>G]C",
            "C[C>G]G",
            "C[C>G]T",
            "G[C>G]A",
            "G[C>G]C",
            "G[C>G]G",
            "G[C>G]T",
            "T[C>G]A",
            "T[C>G]C",
            "T[C>G]G",
            "T[C>G]T",
            "A[C>T]A",
            "A[C>T]C",
            "A[C>T]G",
            "A[C>T]T",
            "C[C>T]A",
            "C[C>T]C",
            "C[C>T]G",
            "C[C>T]T",
            "G[C>T]A",
            "G[C>T]C",
            "G[C>T]G",
            "G[C>T]T",
            "T[C>T]A",
            "T[C>T]C",
            "T[C>T]G",
            "T[C>T]T",
            "A[T>A]A",
            "A[T>A]C",
            "A[T>A]G",
            "A[T>A]T",
            "C[T>A]A",
            "C[T>A]C",
            "C[T>A]G",
            "C[T>A]T",
            "G[T>A]A",
            "G[T>A]C",
            "G[T>A]G",
            "G[T>A]T",
            "T[T>A]A",
            "T[T>A]C",
            "T[T>A]G",
            "T[T>A]T",
            "A[T>C]A",
            "A[T>C]C",
            "A[T>C]G",
            "A[T>C]T",
            "C[T>C]A",
            "C[T>C]C",
            "C[T>C]G",
            "C[T>C]T",
            "G[T>C]A",
            "G[T>C]C",
            "G[T>C]G",
            "G[T>C]T",
            "T[T>C]A",
            "T[T>C]C",
            "T[T>C]G",
            "T[T>C]T",
            "A[T>G]A",
            "A[T>G]C",
            "A[T>G]G",
            "A[T>G]T",
            "C[T>G]A",
            "C[T>G]C",
            "C[T>G]G",
            "C[T>G]T",
            "G[T>G]A",
            "G[T>G]C",
            "G[T>G]G",
            "G[T>G]T",
            "T[T>G]A",
            "T[T>G]C",
            "T[T>G]G",
            "T[T>G]T",
        ]
        xlabels = []

        x = 0.4
        ymax = 0
        colors = [
            [3 / 256, 189 / 256, 239 / 256],
            [1 / 256, 1 / 256, 1 / 256],
            [228 / 256, 41 / 256, 38 / 256],
            [203 / 256, 202 / 256, 202 / 256],
            [162 / 256, 207 / 256, 99 / 256],
            [236 / 256, 199 / 256, 197 / 256],
        ]
        xlabels = [seq[0] + seq[2] + seq[6] for seq in seq96]
        i = 0

        x = 0.043
        y3 = 0.87
        y = int(ymax * 1.25)
        y2 = y + 2
        for i in range(0, 6, 1):
            panel1.add_patch(
                plt.Rectangle(
                    (x, y3),
                    0.15,
                    0.05,
                    facecolor=colors[i],
                    clip_on=False,
                    transform=plt.gcf().transFigure,
                )
            )
            x += 0.159

        yText = y3 + 0.06
        plt.text(
            0.1,
            yText,
            "C>A",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.255,
            yText,
            "C>G",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.415,
            yText,
            "C>T",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.575,
            yText,
            "T>A",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.735,
            yText,
            "T>C",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.89,
            yText,
            "T>G",
            fontsize=55,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )

        if y <= 4:
            y += 4

        while y % 4 != 0:
            y += 1
        # ytick_offest = int(y/4)
        y = ymax / 1.025
        ytick_offest = float(y / 3)
        labs = np.arange(0.375, 96.375, 1)

        panel1.set_xlim([0, 96])
        # panel1.set_ylim([0, y])
        panel1.set_xticks(labs)
        # panel1.set_yticks(ylabs)
        count = 0
        m = 0
        for i in range(0, 96, 1):
            plt.text(
                i / 101 + 0.0415,
                0.02,
                xlabels[i][0],
                fontsize=30,
                color="gray",
                rotation="vertical",
                verticalalignment="center",
                fontname="Courier New",
                transform=plt.gcf().transFigure,
            )
            plt.text(
                i / 101 + 0.0415,
                0.044,
                xlabels[i][1],
                fontsize=30,
                color=colors[m],
                rotation="vertical",
                verticalalignment="center",
                fontname="Courier New",
                fontweight="bold",
                transform=plt.gcf().transFigure,
            )
            plt.text(
                i / 101 + 0.0415,
                0.071,
                xlabels[i][2],
                fontsize=30,
                color="gray",
                rotation="vertical",
                verticalalignment="center",
                fontname="Courier New",
                transform=plt.gcf().transFigure,
            )
            count += 1
            if count == 16:
                count = 0
                m += 1

        plt.gca().yaxis.grid(True)
        plt.gca().grid(which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1)
        panel1.set_xlabel("")
        panel1.set_ylabel("")

        panel1.tick_params(
            axis="both",
            which="both",
            bottom=False,
            labelbottom=False,
            left=True,
            labelleft=True,
            right=True,
            labelright=False,
            top=False,
            labeltop=False,
            direction="in",
            length=25,
            colors="lightgray",
            width=2,
        )

        [i.set_color("black") for i in plt.gca().get_yticklabels()]

        # if return_plot_template == False:
        #     pickle.dump(plot1, open(path, "wb"))
        # else:
        #     pickle.dump(plot1, open(path, "wb"))
        #     return plot1

        return plot1

    elif context == "DBS78":
        # newly added
        # Set default background colors
        plt.rcParams['figure.facecolor'] = 'white'  # For figure background
        plt.rcParams['axes.facecolor'] = 'white'  # For axes background

        plot_custom_text = False
        pcawg = False
        sig_probs = False
        plt.rcParams["axes.linewidth"] = 4
        plot1 = plt.figure(figsize=(43.93, 9.92))
        plt.rc("axes", edgecolor="grey")
        panel1 = plt.axes([0.04, 0.09, 0.95, 0.77])
        xlabels = []

        x = 0.4
        ymax = 0
        colors = [
            [3 / 256, 189 / 256, 239 / 256],
            [3 / 256, 102 / 256, 204 / 256],
            [162 / 256, 207 / 256, 99 / 256],
            [1 / 256, 102 / 256, 1 / 256],
            [255 / 256, 153 / 256, 153 / 256],
            [228 / 256, 41 / 256, 38 / 256],
            [255 / 256, 178 / 256, 102 / 256],
            [255 / 256, 128 / 256, 1 / 256],
            [204 / 256, 153 / 256, 255 / 256],
            [76 / 256, 1 / 256, 153 / 256],
        ]

        x = 0.043
        y3 = 0.87
        y = int(ymax * 1.25)
        y2 = y + 2
        i = 0
        panel1.add_patch(
            plt.Rectangle(
                (0.043, y3),
                0.101,
                0.05,
                facecolor=colors[0],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (0.151, y3),
                0.067,
                0.05,
                facecolor=colors[1],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (0.225, y3),
                0.102,
                0.05,
                facecolor=colors[2],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (0.334, y3),
                0.067,
                0.05,
                facecolor=colors[3],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (0.408, y3),
                0.102,
                0.05,
                facecolor=colors[4],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (0.517, y3),
                0.067,
                0.05,
                facecolor=colors[5],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (0.591, y3),
                0.067,
                0.05,
                facecolor=colors[6],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (0.665, y3),
                0.102,
                0.05,
                facecolor=colors[7],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (0.774, y3),
                0.102,
                0.05,
                facecolor=colors[8],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (0.883, y3),
                0.102,
                0.05,
                facecolor=colors[9],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )

        yText = y3 + 0.06
        plt.text(
            0.07,
            yText,
            "AC>NN",
            fontsize=40,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.163,
            yText,
            "AT>NN",
            fontsize=40,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.255,
            yText,
            "CC>NN",
            fontsize=40,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.345,
            yText,
            "CG>NN",
            fontsize=40,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.435,
            yText,
            "CT>NN",
            fontsize=40,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.527,
            yText,
            "GC>NN",
            fontsize=40,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.6,
            yText,
            "TA>NN",
            fontsize=40,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.69,
            yText,
            "TC>NN",
            fontsize=40,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.8,
            yText,
            "TG>NN",
            fontsize=40,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.915,
            yText,
            "TT>NN",
            fontsize=40,
            fontweight="bold",
            fontname="Arial",
            transform=plt.gcf().transFigure,
        )

        if y <= 4:
            y += 4

        while y % 4 != 0:
            y += 1
        ytick_offest = int(y / 4)

        labs = np.arange(0.44, 78.44, 1)
        panel1.set_xlim([0, 78])
        panel1.set_ylim([0, y])
        panel1.set_xticks(labs)
        panel1.set_xticklabels(
            xlabels,
            rotation="vertical",
            fontsize=30,
            color="grey",
            fontname="Courier New",
            verticalalignment="top",
            fontweight="bold",
        )

        plt.gca().yaxis.grid(True)
        plt.gca().grid(which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1)
        panel1.set_xlabel("")
        panel1.set_ylabel("")
        panel1.tick_params(
            axis="both",
            which="both",
            bottom=False,
            labelbottom=True,
            left=True,
            labelleft=True,
            right=True,
            labelright=False,
            top=False,
            labeltop=False,
            direction="in",
            length=25,
            colors="lightgray",
            width=2,
        )

        [i.set_color("black") for i in plt.gca().get_yticklabels()]
        [i.set_color("grey") for i in plt.gca().get_xticklabels()]

        # if return_plot_template == False:
        #     pickle.dump(plot1, open(path, "wb"))
        # else:
        #     pickle.dump(plot1, open(path, "wb"))
        #     return plot1

        return plot1

    elif context == "ID83":
        # newly added
        # Set default background colors
        plt.rcParams['figure.facecolor'] = 'white'  # For figure background
        plt.rcParams['axes.facecolor'] = 'white'  # For axes background

        plt.rcParams["axes.linewidth"] = 2
        plot1 = plt.figure(figsize=(43.93, 12))
        plt.rc("axes", edgecolor="black")
        panel1 = plt.axes([0.045, 0.17, 0.92, 0.65])
        xlabels = []

        x = 0.4
        ymax = 0
        colors = [
            [253 / 256, 190 / 256, 111 / 256],
            [255 / 256, 128 / 256, 2 / 256],
            [176 / 256, 221 / 256, 139 / 256],
            [54 / 256, 161 / 256, 46 / 256],
            [253 / 256, 202 / 256, 181 / 256],
            [252 / 256, 138 / 256, 106 / 256],
            [241 / 256, 68 / 256, 50 / 256],
            [188 / 256, 25 / 256, 26 / 256],
            [208 / 256, 225 / 256, 242 / 256],
            [148 / 256, 196 / 256, 223 / 256],
            [74 / 256, 152 / 256, 201 / 256],
            [23 / 256, 100 / 256, 171 / 256],
            [226 / 256, 226 / 256, 239 / 256],
            [182 / 256, 182 / 256, 216 / 256],
            [134 / 256, 131 / 256, 189 / 256],
            [98 / 256, 64 / 256, 155 / 256],
        ]

        x = 0.0475
        y_top = 0.827
        y_bottom = 0.114
        y = int(ymax * 1.25)
        y2 = y + 2
        for i in range(0, 12, 1):
            panel1.add_patch(
                plt.Rectangle(
                    (x, y_top),
                    0.0595,
                    0.05,
                    facecolor=colors[i],
                    clip_on=False,
                    transform=plt.gcf().transFigure,
                )
            )
            panel1.add_patch(
                plt.Rectangle(
                    (x, y_bottom),
                    0.0595,
                    0.05,
                    facecolor=colors[i],
                    clip_on=False,
                    transform=plt.gcf().transFigure,
                )
            )
            x += 0.0665

        panel1.add_patch(
            plt.Rectangle(
                (x - 0.001, y_top),
                0.006,
                0.05,
                facecolor=colors[12],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (x - 0.001, y_bottom),
                0.006,
                0.05,
                facecolor=colors[12],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        x += 0.011
        panel1.add_patch(
            plt.Rectangle(
                (x, y_top),
                0.0155,
                0.05,
                facecolor=colors[13],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (x, y_bottom),
                0.0155,
                0.05,
                facecolor=colors[13],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        x += 0.022
        panel1.add_patch(
            plt.Rectangle(
                (x, y_top),
                0.027,
                0.05,
                facecolor=colors[14],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (x, y_bottom),
                0.027,
                0.05,
                facecolor=colors[14],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        x += 0.0335
        panel1.add_patch(
            plt.Rectangle(
                (x, y_top),
                0.049,
                0.05,
                facecolor=colors[15],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )
        panel1.add_patch(
            plt.Rectangle(
                (x, y_bottom),
                0.049,
                0.05,
                facecolor=colors[15],
                clip_on=False,
                transform=plt.gcf().transFigure,
            )
        )

        yText = y_top + 0.01
        plt.text(
            0.072,
            yText,
            "C",
            fontsize=40,
            fontname="Times New Roman",
            fontweight="bold",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.1385,
            yText,
            "T",
            fontsize=40,
            fontname="Times New Roman",
            fontweight="bold",
            color="white",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.205,
            yText,
            "C",
            fontsize=40,
            fontname="Times New Roman",
            fontweight="bold",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.2715,
            yText,
            "T",
            fontsize=40,
            fontname="Times New Roman",
            fontweight="bold",
            color="white",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.338,
            yText,
            "2",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.4045,
            yText,
            "3",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.471,
            yText,
            "4",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.5375,
            yText,
            "5+",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="white",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.604,
            yText,
            "2",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.6705,
            yText,
            "3",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.737,
            yText,
            "4",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.8035,
            yText,
            "5+",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="white",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.844,
            yText,
            "2",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.861,
            yText,
            "3",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.888,
            yText,
            "4",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.93,
            yText,
            "5+",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="white",
            transform=plt.gcf().transFigure,
        )

        yText_labels_top = yText + 0.075
        yText_labels_bottom = y_bottom - 0.03
        yText_labels_bottom_sec = yText_labels_bottom - 0.045

        plt.text(
            0.08,
            yText_labels_top,
            "1bp Deletion",
            fontsize=40,
            fontname="Times New Roman",
            weight="bold",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.21,
            yText_labels_top,
            "1bp Insertion",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.375,
            yText_labels_top,
            ">1bp Deletion at Repeats\n      (Deletion Length)",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.64,
            yText_labels_top,
            ">1bp Insertion at Repeats\n       (Insertion Length)",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.85,
            yText_labels_top,
            " Microhomology\n(Deletion Length)",
            fontsize=40,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )

        plt.text(
            0.058,
            yText_labels_bottom_sec,
            "Homopolymer Length",
            fontsize=35,
            fontname="Times New Roman",
            weight="bold",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.19,
            yText_labels_bottom_sec,
            "Homopolymer Length",
            fontsize=35,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.39,
            yText_labels_bottom_sec,
            "Number of Repeat Units",
            fontsize=35,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.65,
            yText_labels_bottom_sec,
            "Number of Repeat Units",
            fontsize=35,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        plt.text(
            0.85,
            yText_labels_bottom_sec,
            "Microhomology Length",
            fontsize=35,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )

        x = 0.0477
        for i in range(0, 8, 1):
            if i != 2 and i != 3:
                plt.text(
                    x,
                    yText_labels_bottom,
                    "1  2  3  4  5  6+",
                    fontsize=32,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )
            else:
                plt.text(
                    x,
                    yText_labels_bottom,
                    "0  1  2  3  4  5+",
                    fontsize=32,
                    fontweight="bold",
                    fontname="Times New Roman",
                    color="black",
                    transform=plt.gcf().transFigure,
                )

            x += 0.0665

        for i in range(0, 4, 1):
            plt.text(
                x,
                yText_labels_bottom,
                "0  1  2  3  4  5+",
                fontsize=32,
                fontweight="bold",
                fontname="Times New Roman",
                color="black",
                transform=plt.gcf().transFigure,
            )
            x += 0.0665

        plt.text(
            x,
            yText_labels_bottom,
            "1",
            fontsize=32,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        x += 0.011
        plt.text(
            x,
            yText_labels_bottom,
            "1  2",
            fontsize=32,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        x += 0.022
        plt.text(
            x,
            yText_labels_bottom,
            "1  2  3",
            fontsize=32,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )
        x += 0.0335
        plt.text(
            x,
            yText_labels_bottom,
            "1  2  3  4  5+",
            fontsize=32,
            fontweight="bold",
            fontname="Times New Roman",
            color="black",
            transform=plt.gcf().transFigure,
        )

        labs = np.arange(0.375, 83.375, 1)
        panel1.set_xlim([0, 83])
        panel1.set_ylim([0, y])
        panel1.set_xticks(labs)

        plt.gca().yaxis.grid(True)
        plt.gca().grid(which="major", axis="y", color=[0.6, 0.6, 0.6], zorder=1)
        panel1.set_xlabel("")
        panel1.set_ylabel("")

        panel1.tick_params(
            axis="both",
            which="both",
            bottom=False,
            labelbottom=False,
            left=False,
            labelleft=True,
            right=False,
            labelright=False,
            top=False,
            labeltop=False,
            direction="in",
            length=25,
            colors="gray",
            width=2,
        )

        [i.set_color("black") for i in plt.gca().get_yticklabels()]

        # if return_plot_template == False:
        #     pickle.dump(plot1, open(path, "wb"))
        # else:
        #     pickle.dump(plot1, open(path, "wb"))
        #     return plot1

        return plot1


def plotSBS_modified_from_spp(matrix_path,
                    output_path,
                    project,
                    plot_type,
                    # df,
                    strand1_values,
                    strand2_values,
                    asymmetry_type,
                    percentage=False,
                    custom_text_upper=None,
                    custom_text_middle=None,
                    custom_text_bottom=None,
                    savefig_format="pdf",
                    volume=None,
                    dpi=100):

    plot_custom_text = False
    sig_probs = False
    pcawg = False

    # load custom fonts for plotting
    spplt.load_custom_fonts()

    # create the output directory if it doesn't exist
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    if asymmetry_type == REPLICATIONSTRANDBIAS:
        strand1 = LEADING
        strand2 = LAGGING
        strand1_color = 'goldenrod'
        strand2_color = 'indianred'

    elif asymmetry_type == TRANSCRIPTIONSTRANDBIAS:
        strand1 = TRANSCRIBED
        strand2 = UNTRANSCRIBED
        strand1_color = 'royalblue'
        strand2_color = 'yellowgreen'

    elif asymmetry_type == GENICINTERGENICBIAS:
        # Mutation,B,Nontranscribed,Transcribed,UnTranscribed
        # A[C>A]A,0,189,57,44
        # A[C>A]C,0,96,32,31
        # df[GENIC] = df[TRANSCRIBED] + df[UNTRANSCRIBED]
        # df[INTERGENIC] = df[NONTRANSCRIBED]
        strand1 = GENIC
        strand2 = INTERGENIC
        strand1_color = 'cyan'
        strand2_color = 'gray'

    if plot_type == "96":
        try:
            data = spplt.process_input(matrix_path, plot_type)
            data = spplt.reindex_sbs96(data)
            sample_count = 0

            # buf = io.BytesIO()
            plot1 = make_pickle_file_modified_from_spp(context="SBS96",
                                                          return_plot_template=True,
                                                          volume=volume
            )

            # pickle.dump(fig_orig, buf)

            # figs = {}
            # buff_list = {}

            ctx = data.index  # [seq[0]+seq[2]+seq[6] for seq in data.index]
            colors = [
                [3 / 256, 189 / 256, 239 / 256],
                [1 / 256, 1 / 256, 1 / 256],
                [228 / 256, 41 / 256, 38 / 256],
                [203 / 256, 202 / 256, 202 / 256],
                [162 / 256, 207 / 256, 99 / 256],
                [236 / 256, 199 / 256, 197 / 256],
            ]
            colorsall = [
                [colors[j] for i in range(int(len(ctx) / 6))] for j in range(6)
            ]
            colors_flat_list = [item for sublist in colorsall for item in sublist]

            for sample in data.columns:
                # buf.seek(0)
                # figs[sample] = pickle.load(buf)
                # panel1 = figs[sample].axes[0]

                panel1 = plot1.axes[0]
                total_count = np.sum(data[sample].values)
                x = 0.4
                ymax = 0
                i = 0
                muts = data[sample].values

                # newly added
                strand1_muts = strand1_values # df[strand1].values
                strand2_muts = strand2_values # df[strand2].values
                total_muts = np.sum(strand1_muts) + np.sum(strand2_muts)

                if percentage:
                    if total_muts > 0:

                        # plt.bar(
                        #     np.arange(len(ctx)) + x,
                        #     muts / total_count * 100,
                        #     width=0.4,
                        #     color=colors_flat_list,
                        #     align="center",
                        #     zorder=1000,
                        # )
                        # ymax = np.max(muts / total_count * 100)

                        plt.bar(np.arange(len(ctx)) + x,
                                strand1_muts / total_muts * 100,
                                color=strand1_color,
                                label=strand1,
                                width=0.4,
                                align="center",
                                zorder=1000,)

                        plt.bar(np.arange(len(ctx)) + x,
                                strand2_muts / total_muts * 100,
                                color=strand2_color,
                                bottom=strand1_muts / total_muts * 100,
                                label=strand2,
                                width=0.4,
                                align="center",
                                zorder=1000,)

                        ymax = np.max((strand1_muts + strand2_muts) / total_muts * 100)

                        # newly added
                        plt.legend(loc='upper right', fontsize=30)

                    sig_probs = True
                else:
                    # plt.bar(
                    #     np.arange(len(ctx)) + x,
                    #     muts,
                    #     width=0.4,
                    #     color=colors_flat_list,
                    #     align="center",
                    #     zorder=1000,
                    # )
                    # ymax = np.max(muts)

                    plt.bar(np.arange(len(ctx)) + x,
                            strand1_muts,
                            color=strand1_color,
                            label=strand1,
                            width=0.4,
                            align="center",
                            zorder=1000, )

                    plt.bar(np.arange(len(ctx)) + x,
                            strand2_muts,
                            color=strand2_color,
                            bottom= strand1_values, # df[strand1],
                            label=strand2,
                            width=0.4,
                            align="center",
                            zorder=1000, )

                    ymax = np.max(strand1_muts + strand2_muts)

                    # newly added
                    plt.legend(loc='upper right', fontsize=30)

                x = 0.043
                y3 = 0.87
                y = int(ymax * 1.25)
                y2 = y + 2

                if y <= 4:
                    y += 4

                while y % 4 != 0:
                    y += 1

                y = ymax / 1.025
                ytick_offest = float(y / 3)

                if percentage:
                    ylabs = [
                        0,
                        round(ytick_offest, 1),
                        round(ytick_offest * 2, 1),
                        round(ytick_offest * 3, 1),
                        round(ytick_offest * 4, 1),
                    ]
                    ylabels = [
                        str(0),
                        str(round(ytick_offest, 1)) + "%",
                        str(round(ytick_offest * 2, 1)) + "%",
                        str(round(ytick_offest * 3, 1)) + "%",
                        str(round(ytick_offest * 4, 1)) + "%",
                    ]
                else:
                    ylabs = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]
                    ylabels = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]

                labs = np.arange(0.375, 96.375, 1)

                font_label_size = 30
                if not percentage:
                    if int(ylabels[3]) >= 1000:
                        font_label_size = 20

                if percentage:
                    if len(ylabels) > 2:
                        font_label_size = 20

                if not percentage:
                    ylabels = spplt.getylabels(ylabels)

                panel1.set_xlim([0, 96])
                panel1.set_ylim([0, y])
                panel1.set_yticks(ylabs)
                if sig_probs:
                    plt.text(
                        0.045,
                        0.75,
                        sample,
                        fontsize=60,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )
                else:
                    plt.text(
                        0.045,
                        0.75,
                        sample + ": " + "{:,}".format(int(total_count)) + " subs",
                        fontsize=60,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )

                panel1.set_yticklabels(ylabels, fontsize=font_label_size)
                plt.gca().yaxis.grid(True)
                plt.gca().grid(
                    which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1
                )
                panel1.set_xlabel("")
                panel1.set_ylabel("")

                custom_text_upper_plot = ""
                try:
                    custom_text_upper[sample_count]
                except:
                    custom_text_upper = False
                try:
                    custom_text_middle[sample_count]
                except:
                    custom_text_middle = False
                try:
                    custom_text_bottom[sample_count]
                except:
                    custom_text_bottom = False

                if custom_text_upper:
                    plot_custom_text = True
                    if len(custom_text_upper[sample_count]) > 40:
                        print(
                            "To add a custom text, please limit the string to <40 characters including spaces."
                        )
                        plot_custom_text = False
                if custom_text_middle:
                    if len(custom_text_middle[sample_count]) > 40:
                        print(
                            "To add a custom text, please limit the string to <40 characters including spaces."
                        )
                        plot_custom_text = False

                if plot_custom_text:
                    x_pos_custom = 0.98
                    if custom_text_upper and custom_text_middle:
                        custom_text_upper_plot = (
                            custom_text_upper[sample_count]
                            + "\n"
                            + custom_text_middle[sample_count]
                        )
                        if custom_text_bottom:
                            custom_text_upper_plot += (
                                "\n" + custom_text_bottom[sample_count]
                            )

                    if custom_text_upper and not custom_text_middle:
                        custom_text_upper_plot = custom_text_upper[sample_count]
                        panel1.text(
                            x_pos_custom,
                            0.78,
                            custom_text_upper_plot,
                            fontsize=40,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )

                    elif custom_text_upper and custom_text_middle:
                        if not custom_text_bottom:
                            panel1.text(
                                x_pos_custom,
                                0.72,
                                custom_text_upper_plot,
                                fontsize=40,
                                weight="bold",
                                color="black",
                                fontname="Arial",
                                transform=plt.gcf().transFigure,
                                ha="right",
                            )
                        else:
                            panel1.text(
                                x_pos_custom,
                                0.68,
                                custom_text_upper_plot,
                                fontsize=40,
                                weight="bold",
                                color="black",
                                fontname="Arial",
                                transform=plt.gcf().transFigure,
                                ha="right",
                            )

                    elif not custom_text_upper and custom_text_middle:
                        custom_text_upper_plot = custom_text_middle[sample_count]
                        panel1.text(
                            x_pos_custom,
                            0.78,
                            custom_text_upper_plot,
                            fontsize=40,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )

                if percentage:
                    plt.ylabel(
                        "Percentage of Single Base Substitutions",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )
                else:
                    plt.ylabel(
                        "Number of Single Base Substitutions",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )

                panel1.tick_params(
                    axis="both",
                    which="both",
                    bottom=False,
                    labelbottom=False,
                    left=True,
                    labelleft=True,
                    right=True,
                    labelright=False,
                    top=False,
                    labeltop=False,
                    direction="in",
                    length=25,
                    colors="lightgray",
                    width=2,
                )

                [i.set_color("black") for i in plt.gca().get_yticklabels()]
                sample_count += 1

            # return output_results(savefig_format, output_path, project, figs, "SBS_96", dpi=dpi)
        except:
            print("There may be an issue with the formatting of your matrix file.")
            pdf_path = output_path + "SBS_96_plots_" + project + ".pdf"
            if os.path.isfile(pdf_path):
                os.remove(pdf_path)


def plotDBS_modified_from_spp(matrix_path,
    output_path,
    project,
    plot_type,
    # df,
    strand1_values,
    strand2_values,
    strand3_values,
    asymmetry_type,
    percentage=False,
    custom_text_upper=None,
    custom_text_middle=None,
    custom_text_bottom=None,
    savefig_format="pdf",
    volume=None,
    dpi=100,):

    # create the output directory if it doesn't exist
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # load custom fonts for plotting
    spplt.load_custom_fonts()

    plot_custom_text = False
    pcawg = False
    sig_probs = False

    if asymmetry_type == REPLICATIONSTRANDBIAS:
        strand1 = LEADING
        strand2 = LAGGING
        strand3 = QUESTIONABLE
        strand1_color = 'goldenrod'
        strand2_color = 'indianred'
        strand3_color = 'darkviolet'

    elif asymmetry_type == TRANSCRIPTIONSTRANDBIAS:
        strand1 = TRANSCRIBED
        strand2 = UNTRANSCRIBED
        strand3 = QUESTIONABLE
        strand1_color = 'royalblue'
        strand2_color = 'yellowgreen'
        strand3_color = 'darkviolet'

    elif asymmetry_type == GENICINTERGENICBIAS:
        # Mutation,B,Nontranscribed,Transcribed,UnTranscribed
        # A[C>A]A,0,189,57,44
        # A[C>A]C,0,96,32,31
        # df[GENIC] = df[TRANSCRIBED] + df[UNTRANSCRIBED]
        # df[INTERGENIC] = df[NONTRANSCRIBED]
        strand1 = GENIC
        strand2 = INTERGENIC
        strand3 = QUESTIONABLE
        strand1_color = 'cyan'
        strand2_color = 'gray'
        strand3_color = 'darkviolet'

    if plot_type == "78" or plot_type == "78DBS" or plot_type == "DBS78":
        data = spplt.process_input(matrix_path, plot_type)

        dinucs = [
            "TT>GG",
            "TT>CG",
            "TT>AG",
            "TT>GC",
            "TT>CC",
            "TT>AC",
            "TT>GA",
            "TT>CA",
            "TT>AA",
            "AC>CA",
            "AC>CG",
            "AC>CT",
            "AC>GA",
            "AC>GG",
            "AC>GT",
            "AC>TA",
            "AC>TG",
            "AC>TT",
            "CT>AA",
            "CT>AC",
            "CT>AG",
            "CT>GA",
            "CT>GC",
            "CT>GG",
            "CT>TG",
            "CT>TC",
            "CT>TA",
            "AT>CA",
            "AT>CC",
            "AT>CG",
            "AT>GA",
            "AT>GC",
            "AT>TA",
            "TG>GT",
            "TG>CT",
            "TG>AT",
            "TG>GC",
            "TG>CC",
            "TG>AC",
            "TG>GA",
            "TG>CA",
            "TG>AA",
            "CC>AA",
            "CC>AG",
            "CC>AT",
            "CC>GA",
            "CC>GG",
            "CC>GT",
            "CC>TA",
            "CC>TG",
            "CC>TT",
            "CG>AT",
            "CG>GC",
            "CG>GT",
            "CG>TC",
            "CG>TA",
            "CG>TT",
            "TC>GT",
            "TC>CT",
            "TC>AT",
            "TC>GG",
            "TC>CG",
            "TC>AG",
            "TC>GA",
            "TC>CA",
            "TC>AA",
            "GC>AA",
            "GC>AG",
            "GC>AT",
            "GC>CA",
            "GC>CG",
            "GC>TA",
            "TA>GT",
            "TA>CT",
            "TA>AT",
            "TA>GG",
            "TA>CG",
            "TA>GC",
        ]

        revcompl = lambda x: "".join(
            [{"A": "T", "C": "G", "G": "C", "T": "A"}[B] for B in x][::-1]
        )
        mutations = OrderedDict()

        try:
            vals = set(list(data.index)) - set(dinucs)
            vals = sorted(list(vals))

            for ech in vals:
                ech_mod = ech.split(">")[0] + ">" + revcompl(ech.split(">")[1])
                data.rename(index={ech: ech_mod}, inplace=True)

            data = data.sort_index()
            ctx = data.index

            xlabels = [dn.split(">")[1] for dn in ctx]
            colors = [
                [3 / 256, 189 / 256, 239 / 256],
                [3 / 256, 102 / 256, 204 / 256],
                [162 / 256, 207 / 256, 99 / 256],
                [1 / 256, 102 / 256, 1 / 256],
                [255 / 256, 153 / 256, 153 / 256],
                [228 / 256, 41 / 256, 38 / 256],
                [255 / 256, 178 / 256, 102 / 256],
                [255 / 256, 128 / 256, 1 / 256],
                [204 / 256, 153 / 256, 255 / 256],
                [76 / 256, 1 / 256, 153 / 256],
            ]
            mainlist = [dn.split(">")[0] for dn in ctx]
            le = LabelEncoder()
            colors_idxs = le.fit_transform(mainlist)
            colors_flat_list = [colors[i] for i in colors_idxs]
            sample_count = 0

            # buf = io.BytesIO()
            plot1 = make_pickle_file_modified_from_spp(
                context="DBS78", return_plot_template=True, volume=volume
            )
            # pickle.dump(fig_orig, buf)
            # figs = {}

            for sample in data.columns:
                # buf.seek(0)
                # figs[sample] = pickle.load(buf)
                # panel1 = figs[sample].axes[0]
                panel1 = plot1.axes[0]

                total_count = np.sum(
                    data[sample].values
                )  # sum(sum(nuc.values()) for nuc in mutations[sample].values())

                x = 0.4
                muts = data[sample].values

                # newly added
                strand1_muts = strand1_values # df[strand1].values
                strand2_muts = strand2_values # df[strand2].values
                strand3_muts = strand3_values # df[strand3].values
                total_muts = np.sum(strand1_muts) + np.sum(strand2_muts) + np.sum(strand3_muts)

                if percentage:
                    if total_muts > 0:

                        # plt.bar(
                        #     np.asarray(range(len(ctx))) + x,
                        #     muts / total_count * 100,
                        #     width=0.4,
                        #     color=colors_flat_list,
                        #     align="center",
                        #     zorder=1000,
                        # )
                        # ymax = np.max(muts / total_count * 100)

                        plt.bar(np.arange(len(ctx)) + x,
                                strand1_muts / total_muts * 100,
                                color=strand1_color,
                                label=strand1 + ' Strand',
                                width=0.4,
                                align="center",
                                zorder=1000,)

                        plt.bar(np.arange(len(ctx)) + x,
                                strand2_muts / total_muts * 100,
                                color=strand2_color,
                                bottom=strand1_muts / total_muts * 100,
                                label=strand2 + ' Strand',
                                width=0.4,
                                align="center",
                                zorder=1000,)

                        if np.sum(strand3_muts) > 0:
                            plt.bar(np.arange(len(ctx)) + x,
                                    strand3_muts / total_muts * 100,
                                    color=strand3_color,
                                    label=strand3,
                                    width=0.4,
                                    align="center",
                                    zorder=1000,)

                        if np.sum(strand3_muts) > 0:
                            ymax = np.max((strand1_muts + strand2_muts + strand3_muts) / total_muts * 100)
                        else:
                            ymax = np.max((strand1_muts + strand2_muts) / total_muts * 100)

                        # newly added
                        plt.legend(loc='upper right', fontsize=30)

                    sig_probs = True
                else:
                    # plt.bar(
                    #     np.asarray(range(len(ctx))) + x,
                    #     muts,
                    #     width=0.4,
                    #     color=colors_flat_list,
                    #     align="center",
                    #     zorder=1000,
                    # )
                    # ymax = np.max(muts)

                    plt.bar(np.arange(len(ctx)) + x,
                            strand1_muts,
                            color=strand1_color,
                            label=strand1 + ' Strand',
                            width=0.4,
                            align="center",
                            zorder=1000, )

                    plt.bar(np.arange(len(ctx)) + x,
                            strand2_muts,
                            color=strand2_color,
                            bottom=strand1_values, # df[strand1],
                            label=strand2 + ' Strand',
                            width=0.4,
                            align="center",
                            zorder=1000, )

                    if np.sum(strand3_muts) > 0:
                        plt.bar(np.arange(len(ctx)) + x,
                                strand3_muts,
                                color=strand3_color,
                                label=strand3,
                                width=0.4,
                                align="center",
                                zorder=1000, )

                    if np.sum(strand3_muts) > 0:
                        ymax = np.max(strand1_muts + strand2_muts + strand3_muts)
                    else:
                        ymax = np.max(strand1_muts + strand2_muts)

                    # newly added
                    plt.legend(loc='upper right', fontsize=30)

                x = 0.043
                y3 = 0.87
                y = int(ymax * 1.25)
                y2 = y + 2
                i = 0

                if y <= 4:
                    y += 4

                while y % 4 != 0:
                    y += 1
                ytick_offest = int(y / 4)

                if percentage:
                    ylabs = [
                        0,
                        round(ytick_offest, 1),
                        round(ytick_offest * 2, 1),
                        round(ytick_offest * 3, 1),
                        round(ytick_offest * 4, 1),
                    ]
                    ylabels = [
                        str(0),
                        str(round(ytick_offest, 1)) + "%",
                        str(round(ytick_offest * 2, 1)) + "%",
                        str(round(ytick_offest * 3, 1)) + "%",
                        str(round(ytick_offest * 4, 1)) + "%",
                    ]
                else:
                    ylabs = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]
                    ylabels = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]

                if sig_probs:
                    plt.text(
                        0.045,
                        0.75,
                        sample,
                        fontsize=60,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )
                else:
                    plt.text(
                        0.045,
                        0.75,
                        sample
                        + ": "
                        + "{:,}".format(int(total_count))
                        + " double subs",
                        fontsize=60,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure,
                    )

                custom_text_upper_plot = ""
                try:
                    custom_text_upper[sample_count]
                except:
                    custom_text_upper = False
                try:
                    custom_text_middle[sample_count]
                except:
                    custom_text_middle = False
                try:
                    custom_text_bottom[sample_count]
                except:
                    custom_text_bottom = False

                if custom_text_upper:
                    plot_custom_text = True
                    if len(custom_text_upper[sample_count]) > 40:
                        print(
                            "To add a custom text, please limit the string to <40 characters including spaces."
                        )
                        plot_custom_text = False
                if custom_text_bottom:
                    if len(custom_text_bottom[sample_count]) > 40:
                        print(
                            "To add a custom text, please limit the string to <40 characters including spaces."
                        )
                        plot_custom_text = False

                if plot_custom_text:
                    x_pos_custom = 0.98
                    if custom_text_upper and custom_text_middle:
                        custom_text_upper_plot = (
                            custom_text_upper[sample_count]
                            + "\n"
                            + custom_text_middle[sample_count]
                        )
                        if custom_text_bottom:
                            custom_text_upper_plot += (
                                "\n" + custom_text_bottom[sample_count]
                            )

                    if custom_text_upper and not custom_text_middle:
                        custom_text_upper_plot = custom_text_upper[sample_count]
                        panel1.text(
                            x_pos_custom,
                            0.78,
                            custom_text_upper_plot,
                            fontsize=35,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )

                    elif custom_text_upper and custom_text_middle:
                        if not custom_text_bottom:
                            panel1.text(
                                x_pos_custom,
                                0.75,
                                custom_text_upper_plot,
                                fontsize=35,
                                weight="bold",
                                color="black",
                                fontname="Arial",
                                transform=plt.gcf().transFigure,
                                ha="right",
                            )
                        else:
                            panel1.text(
                                x_pos_custom,
                                0.7,
                                custom_text_upper_plot,
                                fontsize=35,
                                weight="bold",
                                color="black",
                                fontname="Arial",
                                transform=plt.gcf().transFigure,
                                ha="right",
                            )

                    elif not custom_text_upper and custom_text_middle:
                        custom_text_upper_plot = custom_text_middle[sample_count]
                        panel1.text(
                            x_pos_custom,
                            0.78,
                            custom_text_upper_plot,
                            fontsize=35,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure,
                            ha="right",
                        )

                if not percentage:
                    ylabels = spplt.getylabels(ylabels)

                labs = np.arange(0.44, 78.44, 1)
                panel1.set_xlim([0, 78])
                panel1.set_ylim([0, y])
                panel1.set_xticks(labs)
                panel1.set_yticks(ylabs)
                panel1.set_xticklabels(
                    xlabels,
                    rotation="vertical",
                    fontsize=30,
                    color="grey",
                    fontname="Courier New",
                    verticalalignment="top",
                    fontweight="bold",
                )

                panel1.set_yticklabels(ylabels, fontsize=25)
                plt.gca().yaxis.grid(True)
                plt.gca().grid(
                    which="major", axis="y", color=[0.93, 0.93, 0.93], zorder=1
                )
                panel1.set_xlabel("")
                panel1.set_ylabel("")

                if percentage:
                    plt.ylabel(
                        "Percentage of Double Base Substitutions",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )
                else:
                    plt.ylabel(
                        "Number of Double Base Substitutions",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )

                panel1.tick_params(
                    axis="both",
                    which="both",
                    bottom=False,
                    labelbottom=True,
                    left=True,
                    labelleft=True,
                    right=True,
                    labelright=False,
                    top=False,
                    labeltop=False,
                    direction="in",
                    length=25,
                    colors="lightgray",
                    width=2,
                )

                [i.set_color("black") for i in plt.gca().get_yticklabels()]
                [i.set_color("grey") for i in plt.gca().get_xticklabels()]
                sample_count += 1

            # return output_results(savefig_format, output_path, project, figs, "DBS_78", dpi=dpi)

        except:
            print("There may be an issue with the formatting of your matrix file.")
            pdf_path = output_path + "DBS_78_plots_" + project + ".pdf"
            if os.path.isfile(pdf_path):
                os.remove(pdf_path)



def plotID_modified_from_spp(matrix_path,
    output_path,
    project,
    plot_type,
    # df,
    strand1_values,
    strand2_values,
    strand3_values,
    asymmetry_type,
    percentage=False,
    custom_text_upper=None,
    custom_text_middle=None,
    custom_text_bottom=None,
    savefig_format="pdf",
    volume=None,
    dpi=100,):

    # create the output directory if it doesn't exist
    if not os.path.exists(output_path) and savefig_format.lower() != "pil_image":
        os.makedirs(output_path)

    # load custom fonts for plotting
    spplt.load_custom_fonts()

    plot_custom_text = False
    sig_probs = False
    pcawg = False

    if asymmetry_type == REPLICATIONSTRANDBIAS:
        strand1 = LEADING
        strand2 = LAGGING
        strand3 = QUESTIONABLE
        strand1_color = 'goldenrod'
        strand2_color = 'indianred'
        strand3_color = 'darkviolet'

    elif asymmetry_type == TRANSCRIPTIONSTRANDBIAS:
        strand1 = TRANSCRIBED
        strand2 = UNTRANSCRIBED
        strand3 = QUESTIONABLE
        strand1_color = 'royalblue'
        strand2_color = 'yellowgreen'
        strand3_color = 'darkviolet'

    elif asymmetry_type == GENICINTERGENICBIAS:
        # Mutation,B,Nontranscribed,Transcribed,UnTranscribed
        # A[C>A]A,0,189,57,44
        # A[C>A]C,0,96,32,31
        # df[GENIC] = df[TRANSCRIBED] + df[UNTRANSCRIBED]
        # df[INTERGENIC] = df[NONTRANSCRIBED]
        strand1 = GENIC
        strand2 = INTERGENIC
        strand3 = QUESTIONABLE
        strand1_color = 'cyan'
        strand2_color = 'gray'
        strand3_color = 'darkviolet'

    if (
        plot_type == "94"
        or plot_type == "ID94"
        or plot_type == "94ID"
        or plot_type == "83"
    ):
        data = spplt.process_input(matrix_path, plot_type)

        try:
            sample_count = 0
            # buf = io.BytesIO()
            # fig_orig = spplt.make_pickle_file(
            #     context="ID83", return_plot_template=True, volume=volume
            # )

            plot1 = make_pickle_file_modified_from_spp(context="ID83",
                                                       return_plot_template=True)

            # pickle.dump(fig_orig, buf)
            # figs = {}
            colors = [
                [253 / 256, 190 / 256, 111 / 256],
                [255 / 256, 128 / 256, 2 / 256],
                [176 / 256, 221 / 256, 139 / 256],
                [54 / 256, 161 / 256, 46 / 256],
                [253 / 256, 202 / 256, 181 / 256],
                [252 / 256, 138 / 256, 106 / 256],
                [241 / 256, 68 / 256, 50 / 256],
                [188 / 256, 25 / 256, 26 / 256],
                [208 / 256, 225 / 256, 242 / 256],
                [148 / 256, 196 / 256, 223 / 256],
                [74 / 256, 152 / 256, 201 / 256],
                [23 / 256, 100 / 256, 171 / 256],
                [226 / 256, 226 / 256, 239 / 256],
                [182 / 256, 182 / 256, 216 / 256],
                [134 / 256, 131 / 256, 189 / 256],
                [98 / 256, 64 / 256, 155 / 256],
            ]
            ctx = data.index
            xlabels = [
                i.split(":")[0] + i.split(":")[1] + i.split(":")[2]
                for i in ctx.to_list()
            ]
            xlables_set = [
                "1DelC",
                "1DelT",
                "1InsC",
                "1InsT",
                "2DelR",
                "3DelR",
                "4DelR",
                "5DelR",
                "2InsR",
                "3InsR",
                "4InsR",
                "5InsR",
                "2DelM",
                "3DelM",
                "4DelM",
                "5DelM",
            ]
            colors_idx = copy.deepcopy(xlabels)
            for ii in range(0, len(xlables_set)):
                colors_idx = [ii if x == xlables_set[ii] else x for x in colors_idx]

            colors_flat_list = [colors[i] for i in colors_idx]

            for sample in data.columns:  # mutations.keys():
                # buf.seek(0)
                # figs[sample] = pickle.load(buf)
                # panel1 = figs[sample].axes[0] # panel1 = figs[sample].axes[0]

                panel1 = plot1.axes[0]
                muts = data[sample].values
                total_count = np.sum(muts)
                x = 0.4

                # newly added
                strand1_muts = strand1_values # df[strand1].values
                strand2_muts = strand2_values # df[strand2].values
                strand3_muts = strand3_values # df[strand3].values
                total_muts = np.sum(strand1_muts) + np.sum(strand2_muts) + np.sum(strand3_muts)

                if percentage:
                    if total_muts > 0:
                        # plt.bar(
                        #     np.arange(len(ctx)) + x,
                        #     muts / total_count * 100,
                        #     width=0.4,
                        #     color=colors_flat_list,
                        #     align="center",
                        #     zorder=1000,
                        # )
                        # ymax = np.max(muts / total_count * 100)

                        plt.bar(np.arange(len(ctx)) + x,
                                strand1_muts / total_muts * 100,
                                color=strand1_color,
                                label=strand1 + ' Strand',
                                width=0.4,
                                align="center",
                                zorder=1000,)

                        plt.bar(np.arange(len(ctx)) + x,
                                strand2_muts / total_muts * 100,
                                color=strand2_color,
                                bottom=strand1_muts / total_muts * 100,
                                label=strand2 + ' Strand',
                                width=0.4,
                                align="center",
                                zorder=1000,)

                        if np.sum(strand3_muts) > 0:
                            plt.bar(np.arange(len(ctx)) + x,
                                    strand3_muts / total_muts * 100,
                                    color=strand3_color,
                                    label=strand3,
                                    width=0.4,
                                    align="center",
                                    zorder=1000,)

                        if np.sum(strand3_muts):
                            ymax = np.max((strand1_muts + strand2_muts + strand3_muts) / total_muts * 100)
                        else:
                            ymax = np.max((strand1_muts + strand2_muts) / total_muts * 100)

                        # newly added
                        plt.legend(loc='upper right', fontsize=30)

                    sig_probs = True
                else:
                    # plt.bar(
                    #     np.arange(len(ctx)) + x,
                    #     muts,
                    #     width=0.4,
                    #     color=colors_flat_list,
                    #     align="center",
                    #     zorder=1000,
                    # )
                    # ymax = np.max(muts)

                    plt.bar(np.arange(len(ctx)) + x,
                            strand1_muts,
                            color=strand1_color,
                            label=strand1 + ' Strand',
                            width=0.4,
                            align="center",
                            zorder=1000, )

                    plt.bar(np.arange(len(ctx)) + x,
                            strand2_muts,
                            color=strand2_color,
                            bottom=strand1_values, # df[strand1],
                            label=strand2 + ' Strand',
                            width=0.4,
                            align="center",
                            zorder=1000, )

                    if np.sum(strand3_muts) > 0:
                        plt.bar(np.arange(len(ctx)) + x,
                                strand3_muts,
                                color=strand3_color,
                                label=strand3,
                                width=0.4,
                                align="center",
                                zorder=1000, )

                    if np.sum(strand3_muts):
                        ymax = np.max(strand1_muts + strand2_muts + strand3_muts)
                    else:
                        ymax = np.max(strand1_muts + strand2_muts)

                    # newly added
                    plt.legend(loc='upper right', fontsize=30)

                x = 0.0475
                y_top = 0.827
                y_bottom = 0.114
                y = int(ymax * 1.25)
                y2 = y + 2

                if y <= 4:
                    y += 4

                while y % 4 != 0:
                    y += 1
                ytick_offest = int(y / 4)

                if percentage:
                    ylabs = [
                        0,
                        round(ytick_offest, 1),
                        round(ytick_offest * 2, 1),
                        round(ytick_offest * 3, 1),
                        round(ytick_offest * 4, 1),
                    ]
                    ylabels = [
                        str(0),
                        str(round(ytick_offest, 1)) + "%",
                        str(round(ytick_offest * 2, 1)) + "%",
                        str(round(ytick_offest * 3, 1)) + "%",
                        str(round(ytick_offest * 4, 1)) + "%",
                    ]
                else:
                    ylabs = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]
                    ylabels = [
                        0,
                        ytick_offest,
                        ytick_offest * 2,
                        ytick_offest * 3,
                        ytick_offest * 4,
                    ]

                labs = np.arange(0.375, 83.375, 1)

                if not percentage:
                    ylabels = spplt.getylabels(ylabels)

                panel1.set_xlim([0, 83]) # panel1.set_xlim([0, 83])
                panel1.set_ylim([0, y]) # panel1.set_ylim([0, y])
                panel1.set_xticks(labs) # panel1.set_xticks(labs)
                panel1.set_yticks(ylabs) # panel1.set_yticks(ylabs)
                if sig_probs:
                    plt.text( # plt.text(
                        0.0475,
                        0.75,
                        sample,
                        fontsize=60,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure, # plt.gcf().transFigure
                    )
                else:
                    plt.text( #  plt.text(
                        0.0475,
                        0.75,
                        sample + ": " + "{:,}".format(int(total_count)) + " indels",
                        fontsize=60,
                        weight="bold",
                        color="black",
                        fontname="Arial",
                        transform=plt.gcf().transFigure, # plt.gcf().transFigure,
                    )

                custom_text_upper_plot = ""
                try:
                    custom_text_upper[sample_count]
                except:
                    custom_text_upper = False
                try:
                    custom_text_middle[sample_count]
                except:
                    custom_text_middle = False
                try:
                    custom_text_bottom[sample_count]
                except:
                    custom_text_bottom = False

                if custom_text_upper:
                    plot_custom_text = True
                    if len(custom_text_upper[sample_count]) > 40:
                        print(
                            "To add a custom text, please limit the string to <40 characters including spaces."
                        )
                        plot_custom_text = False
                if custom_text_middle:
                    if len(custom_text_middle[sample_count]) > 40:
                        print(
                            "To add a custom text, please limit the string to <40 characters including spaces."
                        )
                        plot_custom_text = False

                if plot_custom_text:
                    x_pos_custom = 0.95
                    if custom_text_upper and custom_text_middle:
                        custom_text_upper_plot = (
                            custom_text_upper[sample_count]
                            + "\n"
                            + custom_text_middle[sample_count]
                        )
                        if custom_text_bottom:
                            custom_text_upper_plot += (
                                "\n" + custom_text_bottom[sample_count]
                            )

                    if custom_text_upper and not custom_text_middle:
                        custom_text_upper_plot = custom_text_upper[sample_count]
                        panel1.text( # panel1.text(
                            x_pos_custom,
                            0.78,
                            custom_text_upper_plot,
                            fontsize=35,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure, # plt.gcf().transFigure
                            ha="right",
                        )

                    elif custom_text_upper and custom_text_middle:
                        if not custom_text_bottom:
                            panel1.text( # panel1.text(
                                x_pos_custom,
                                0.72,
                                custom_text_upper_plot,
                                fontsize=35,
                                weight="bold",
                                color="black",
                                fontname="Arial",
                                transform=plt.gcf().transFigure, # plt.gcf().transFigure,
                                ha="right",
                            )
                        else:
                            panel1.text( # panel1.text(
                                x_pos_custom,
                                0.68,
                                custom_text_upper_plot,
                                fontsize=35,
                                weight="bold",
                                color="black",
                                fontname="Arial",
                                transform=plt.gcf().transFigure, # plt.gcf().transFigure,
                                ha="right",
                            )

                    elif not custom_text_upper and custom_text_middle:
                        custom_text_upper_plot = custom_text_middle[sample_count]
                        panel1.text( # panel1.text(
                            x_pos_custom,
                            0.78,
                            custom_text_upper_plot,
                            fontsize=35,
                            weight="bold",
                            color="black",
                            fontname="Arial",
                            transform=plt.gcf().transFigure, # plt.gcf().transFigure,
                            ha="right",
                        )

                panel1.set_yticklabels(ylabels, fontsize=30) # panel1.set_yticklabels(ylabels, fontsize=30)
                plt.gca().yaxis.grid(True) # plt.gca().yaxis.grid(True)
                plt.gca().grid(which="major", axis="y", color=[0.6, 0.6, 0.6], zorder=1) # plt.gca().
                panel1.set_xlabel("") # panel1.set_xlabel("")
                panel1.set_ylabel("") # panel1.set_ylabel("")

                if percentage:
                    plt.ylabel( # plt.ylabel
                        "Percentage of Indels",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )
                else:
                    plt.ylabel( # plt.ylabel(
                        "Number of Indels",
                        fontsize=35,
                        fontname="Times New Roman",
                        weight="bold",
                    )

                panel1.tick_params( # panel1.tick_params(
                    axis="both",
                    which="both",
                    bottom=False,
                    labelbottom=False,
                    left=False,
                    labelleft=True,
                    right=False,
                    labelright=False,
                    top=False,
                    labeltop=False,
                    direction="in",
                    length=25,
                    colors="gray",
                    width=2,
                )

                [i.set_color("black") for i in plt.gca().get_yticklabels()]

                sample_count += 1

            # return output_results(
            #     savefig_format, output_path, project, figs, pdf, "ID_83", dpi=dpi
            # )

        except:
            print("There may be an issue with the formatting of your matrix file.")
            pdf_path = output_path + "ID_83_plots_" + project + ".pdf"
            if os.path.isfile(pdf_path):
                os.remove(pdf_path)

    else:
        print(
            "Error: The function plotID does not support plot_type",
            plot_type,
            "so no plot has been generated."
        )