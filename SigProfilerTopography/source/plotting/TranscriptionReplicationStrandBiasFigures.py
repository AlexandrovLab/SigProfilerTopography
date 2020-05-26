# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018-2020 Burcak Otlu

import sys
import os
import shutil
import numpy as np
import scipy.stats as stats
import statsmodels.stats.multitest

import matplotlib
BACKEND = 'Agg'
if matplotlib.get_backend().lower() != BACKEND.lower():
    # If backend is not set properly a call to describe will hang
    matplotlib.use(BACKEND)

from matplotlib import pyplot as plt
import matplotlib as mpl
import pandas as pd

from SigProfilerTopography.source.commons.TopographyCommons import natural_key
from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIBED_STRAND
from SigProfilerTopography.source.commons.TopographyCommons import UNTRANSCRIBED_STRAND
from SigProfilerTopography.source.commons.TopographyCommons import LAGGING
from SigProfilerTopography.source.commons.TopographyCommons import LEADING

from SigProfilerTopography.source.commons.TopographyCommons import sixMutationTypes
from SigProfilerTopography.source.commons.TopographyCommons import STRANDBIAS
from SigProfilerTopography.source.commons.TopographyCommons import DATA
from SigProfilerTopography.source.commons.TopographyCommons import FIGURE
from SigProfilerTopography.source.commons.TopographyCommons import ALL
from SigProfilerTopography.source.commons.TopographyCommons import SCATTERPLOTS
from SigProfilerTopography.source.commons.TopographyCommons import SAMPLES

from SigProfilerTopography.source.commons.TopographyCommons import SUBS_STRAND_BIAS_NUMBER_OF_MUTATIONS_THRESHOLD

from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIPTIONSTRANDBIAS
from SigProfilerTopography.source.commons.TopographyCommons import REPLICATIONSTRANDBIAS

from SigProfilerTopography.source.commons.TopographyCommons import getDictionary

from SigProfilerTopography.source.commons.TopographyCommons import Table_SubsSignature_Cutoff_NumberofMutations_AverageProbability_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Table_IndelsSignature_Cutoff_NumberofMutations_AverageProbability_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Table_DinucsSignature_Cutoff_NumberofMutations_AverageProbability_Filename

from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofSubsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofIndelsDict
from SigProfilerTopography.source.commons.TopographyCommons import Sample2NumberofDinucsDictFilename

from SigProfilerTopography.source.commons.TopographyCommons import getSample2SubsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2IndelsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import Sample2DinucsSignature2NumberofMutationsDictFilename

from SigProfilerTopography.source.commons.TopographyCommons import readDictionary

from SigProfilerTopography.source.commons.TopographyCommons import Type2TranscriptionStrand2CountDict_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Signature2MutationType2TranscriptionStrand2CountDict_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Sample2Type2TranscriptionStrand2CountDict_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Type2Sample2TranscriptionStrand2CountDict_Filename

from SigProfilerTopography.source.commons.TopographyCommons import Type2ReplicationStrand2CountDict_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Signature2MutationType2ReplicationStrand2CountDict_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Sample2Type2ReplicationStrand2CountDict_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Type2Sample2ReplicationStrand2CountDict_Filename

from SigProfilerTopography.source.commons.TopographyCommons import writeDictionarySimple

transcriptionStrands = [TRANSCRIBED_STRAND, UNTRANSCRIBED_STRAND]
replicationStrands = [LAGGING, LEADING]

plt.rcParams.update({'figure.max_open_warning': 0})

########################################################################
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
    for mutationType in sixMutationTypes:
        if (mutationType in type2TranscriptionStrand2CountDict) and (mutationType in type2ReplicationStrand2CountDict):

            if ((TRANSCRIBED_STRAND in type2TranscriptionStrand2CountDict[mutationType]) and (UNTRANSCRIBED_STRAND in type2TranscriptionStrand2CountDict[mutationType])):
                transcriptionRatiosDict[mutationType]= np.log10(type2TranscriptionStrand2CountDict[mutationType][TRANSCRIBED_STRAND]/type2TranscriptionStrand2CountDict[mutationType][UNTRANSCRIBED_STRAND])

            if ((LAGGING in type2ReplicationStrand2CountDict[mutationType]) and (LEADING in type2ReplicationStrand2CountDict[mutationType])):
                replicationRatiosDict[mutationType] = np.log10(type2ReplicationStrand2CountDict[mutationType][LAGGING]/type2ReplicationStrand2CountDict[mutationType][LEADING])

            if (mutationType in replicationRatiosDict) and (mutationType in transcriptionRatiosDict):
                plt.scatter(replicationRatiosDict[mutationType],transcriptionRatiosDict[mutationType], label=mutationType)
    ########################################################################

    legend = plt.legend(loc='upper left', frameon=True, fancybox =False,labels=sixMutationTypes, bbox_to_anchor=(-0.0095, 1.0095))
    legend.get_frame().set_linewidth(1)

    frame = legend.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('black')

    plt.axvline(x=0.0, color='gray', linestyle='--')
    plt.axhline(y=0.0, color='gray', linestyle='--')

    if sample is None:
        figureName = 'all_mutationtypes_%s_scatterplots.png' %(STRANDBIAS)
        figureFile = os.path.join(outputDir, jobname, FIGURE, ALL, STRANDBIAS,SCATTERPLOTS,figureName)

    else:
        figureName = 'all_mutationtypes_%s_%d_%s_scatterplots.png' %(sample,numberofMutations,STRANDBIAS)
        os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, STRANDBIAS,SCATTERPLOTS), exist_ok=True)
        figureFile = os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, STRANDBIAS, SCATTERPLOTS, figureName)

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
            figureName = 'all_%s_signatures_%s_scatterplots.png' % (signatureType, STRANDBIAS)
            figureFile = os.path.join(outputDir, jobname, FIGURE, ALL, STRANDBIAS,SCATTERPLOTS,figureName)
        else:
            figureName = 'all_%s_signatures_%s_%d_%s_scatterplots.png' % (
            signatureType, sample, numberofMutations, STRANDBIAS)
            os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, STRANDBIAS,SCATTERPLOTS), exist_ok=True)
            figureFile = os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, STRANDBIAS, SCATTERPLOTS, figureName)

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
    for mutationType in sixMutationTypes:
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

    for mutationType in sixMutationTypes:
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
        figureFile = os.path.join(outputDir,jobname,FIGURE,ALL,STRANDBIAS,SCATTERPLOTS,figureName)
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
            figureFile = os.path.join(outputDir,jobname,FIGURE,ALL,STRANDBIAS,SCATTERPLOTS,figureName)
            fig.savefig(figureFile)
            plt.cla()
            plt.close(fig)
########################################################################


##################################################################
#Only this method supports simulations
#key can be a sample or a signature
def plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,key,isKeySample,numberofMutations,N,x_axis_labels,strand1_values,strand2_values,strand1_simulations_median_values ,strand2_simulations_median_values , fdr_bh_adjusted_pvalues, strand1Name, strand2Name, mutationsOrSignatures, color1, color2, figureName, width):
    # print('############# for debug starts Nov 12, 2018 ################')
    # print('Sample:%s --- Strand1:%s --- Strand2:%s --- mutationsOrSignatures:%s' %(sample,strand1Name,strand2Name,mutationsOrSignatures))
    # print('strand1_values: %s' %(strand1_values))
    # print('strand2_values: %s' %(strand2_values))
    # print('strand1_simulations_median_values: %s' %(strand1_simulations_median_values))
    # print('strand2_simulations_median_values: %s' %(strand2_simulations_median_values))
    # print('fdr_bh_adjusted_pvalues: %s' %(fdr_bh_adjusted_pvalues))
    # print('############# for debug ends Nov 12, 2018 ################')

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

    ##############################
    #To make the bar width not too wide
    if len(ind)<6:
        maxn = 6
        ax.set_xlim(-0.5, maxn - 0.5)
    ##############################

    ax.get_yaxis().set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

    rects1 = ax.bar(ind, strand1_values, width=width, color=color1)
    rects2 = ax.bar(ind + width, strand2_values, width=width, color=color2)

    if ((strand1_simulations_median_values is not None) and strand1_simulations_median_values):
        rects3 = ax.bar(ind+ 2*width, strand1_simulations_median_values, width=width, color=color1, hatch = '///')
    if ((strand2_simulations_median_values is not None) and strand2_simulations_median_values):
        rects4 = ax.bar(ind +3*width, strand2_simulations_median_values, width=width, color=color2, hatch = '///')

    # # Provide average or medians of the simulatons values
    # if (strand1_simulations_median_values is not None):
    #     ax.plot(ind,strand1_simulations_median_values,color=color1, marker='o', linestyle='dashed')
    # if (strand2_simulations_median_values is not None):
    #     ax.plot(ind + width, strand2_simulations_median_values, color=color2, marker='o', linestyle='dashed')

    ax.tick_params(axis='x', labelsize=35)
    ax.tick_params(axis='y', labelsize=35)

    if len(x_axis_labels)>6 :
        ax.set_xticklabels(x_axis_labels, fontsize=15,rotation=90)
    else:
        ax.set_xticklabels(x_axis_labels, fontsize=35)

    locs, labels = plt.yticks()
    ax.set_ylim(0, locs[-1] + 5000)

    # add some text for labels, title and axes ticks
    # ax.set_ylabel('Number of mutations', fontsize=30)
    if key is not None:
        ax.set_title('%s %s vs. %s %s' %(key,strand1Name,strand2Name,mutationsOrSignatures), fontsize=20,fontweight='bold')
    else:
        ax.set_title('%s vs. %s %s' %(strand1Name,strand2Name,mutationsOrSignatures), fontsize=20,fontweight='bold')


    #set the x axis tick locations
    if (numberofSimulations>0):
        ax.set_xticks(ind + (3 * width)/2)
        simulationsStrand1Name = 'Simulated %s' %(strand1Name)
        simulationsStrand2Name = 'Simulated %s' % (strand2Name)

        # legend = ax.legend((rects1[0], rects2[0], rects3[0], rects4[0]),
        #                    (strand1Name, strand2Name, simulationsStrand1Name, simulationsStrand2Name),
        #                    prop={'size': 25}, bbox_to_anchor = (0, 1.21), ncol = 2, loc = 'upper left')

        if ((rects1 is not None) and (rects2 is not None) and (rects3 is not None) and (rects4 is not None)):
            if ((len(rects1)>0) and (len(rects2)>0) and (len(rects3)>0) and (len(rects4)>0)):
                legend = ax.legend((rects1[0], rects2[0], rects3[0], rects4[0]),
                           (strand1Name, strand2Name, simulationsStrand1Name, simulationsStrand2Name),
                           prop={'size': 20}, ncol = 2, loc = 'best')

    else:
        #Old way with no simulations
        ax.set_xticks(ind + width/2)

        if ((rects1 is not None) and (rects2 is not None)):
            if ((len(rects1)>0) and (len(rects2)>0)):
                legend = ax.legend((rects1[0], rects2[0]), (strand1Name, strand2Name), prop={'size': 25}, loc='upper right')

    #Set the ylabel
    plt.ylabel('Number of single point mutations', fontsize=35, fontweight='normal')

    #Legend place is modified here.
    # ax.legend((rects1[0], rects2[0]), (strand1Name, strand2Name), prop={'size': 30}, loc='upper left')

    #To make the barplot background white
    ax.set_facecolor('white')
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
            if ((fdr_bh_adjusted_pvalue)<=0.0001):
                plt.annotate(
                    '***',  # Use `label` as label
                    (x_value, y_value),  # Place label at end of the bar
                    xytext=(0, space),  # Vertically shift label by `space`
                    textcoords="offset points",  # Interpret `xytext` as offset in points
                    ha='center',  # Horizontally center label
                    va=va,
                    fontsize=20)  # Vertically align label differently for

            elif ((fdr_bh_adjusted_pvalue)<=0.001):
                plt.annotate(
                    '**',  # Use `label` as label
                    (x_value, y_value),  # Place label at end of the bar
                    xytext=(0, space),  # Vertically shift label by `space`
                    textcoords="offset points",  # Interpret `xytext` as offset in points
                    ha='center',  # Horizontally center label
                    va=va,
                    fontsize=20)  # Vertically align label differently for

            elif ((fdr_bh_adjusted_pvalue)<=0.05):
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
        figureName = '%s.png' %(figureName)
        figureFile = os.path.join(outputDir, jobname, FIGURE, ALL, STRANDBIAS, figureName)
    elif (not isKeySample):
        figureName = '%s_%s.png' %(key,figureName)
        figureFile = os.path.join(outputDir, jobname, FIGURE, ALL, STRANDBIAS, figureName)
    else:
        figureName = '%s_%s_%d.png' %(figureName,key,numberofMutations)
        os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES, key, STRANDBIAS), exist_ok=True)
        figureFile = os.path.join(outputDir, jobname, FIGURE, SAMPLES, key, STRANDBIAS, figureName)

    fig.savefig(figureFile)
    plt.cla()
    plt.close(fig)
##################################################################


##################################################################
# Fill strand2TypeCountListDict using type2Strand2CountDict w.r.t. referenceTypes
def fillDictWRTReferenceTypes(types,strands,type2Strand2CountDict):
    strand2TypeCountListDict= {}

    for strand in strands:
        strand2TypeCountListDict[strand] = []
        for type in types:
            if type in type2Strand2CountDict:
                if strand in type2Strand2CountDict[type]:
                    count = type2Strand2CountDict[type][strand]
                    strand2TypeCountListDict[strand].append(count)
                else:
                    strand2TypeCountListDict[strand].append(0)
            else:
                strand2TypeCountListDict[strand].append(0)

    return strand2TypeCountListDict
##################################################################

##################################################################
def fillSimulationsType2StrandCountList(
        simNum2Type2Strand2CountDict,
        existingTypesList,
        numberofSimulations,
        strandbias):

    all_simulations_types_strand1_list = []
    all_simulations_types_strand2_list = []

    simulations_types_strand1_medians_list = []
    simulations_types_strand2_medians_list = []

    ##############################################################
    for simNum in range(1, numberofSimulations + 1):
        #new code
        type2Strand2CountDict = simNum2Type2Strand2CountDict[str(simNum)]

        if ((type2Strand2CountDict is not None) and (type2Strand2CountDict)):
            # Let's fill simulations w.r.t.  existingTypesList in the real data
            # if (simulationType2StrandCountDict is not None):
            if (strandbias == TRANSCRIPTIONSTRANDBIAS):
                strand2TypeStrandCountList_Dict = fillDictWRTReferenceTypes(existingTypesList,transcriptionStrands,type2Strand2CountDict)
                simulation_types_strand1_list = strand2TypeStrandCountList_Dict[TRANSCRIBED_STRAND]
                simulation_types_strand2_list = strand2TypeStrandCountList_Dict[UNTRANSCRIBED_STRAND]

            elif (strandbias == REPLICATIONSTRANDBIAS):
                strand2TypeStrandCountList_Dict = fillDictWRTReferenceTypes(existingTypesList,replicationStrands,type2Strand2CountDict)
                simulation_types_strand1_list = strand2TypeStrandCountList_Dict[LAGGING]
                simulation_types_strand2_list = strand2TypeStrandCountList_Dict[LEADING]

            if simulation_types_strand1_list:
                all_simulations_types_strand1_list.append(simulation_types_strand1_list)
            if simulation_types_strand2_list:
                all_simulations_types_strand2_list.append(simulation_types_strand2_list)
    ##############################################################

    ##################################################################
    if ((all_simulations_types_strand1_list is not None) and all_simulations_types_strand1_list):
        all_simulations_types_strand1_nparray = np.vstack(all_simulations_types_strand1_list)
        (rows, cols) = all_simulations_types_strand1_nparray.shape

        for col in range(cols):
            colwise_array = all_simulations_types_strand1_nparray[:, col]
            sorted_colwise_array = np.sort(colwise_array)
            simulations_types_strand1_medians_list.append(np.median(sorted_colwise_array))
    ##################################################################

    ##################################################################
    if ((all_simulations_types_strand2_list is not None) and all_simulations_types_strand2_list):
        all_simulations_types_strand2_nparray = np.vstack(all_simulations_types_strand2_list)
        (rows, cols) = all_simulations_types_strand2_nparray.shape

        for col in range(cols):
            colwise_array = all_simulations_types_strand2_nparray[:, col]
            sorted_colwise_array = np.sort(colwise_array)
            simulations_types_strand2_medians_list.append(np.median(sorted_colwise_array))
    ##################################################################

    return simulations_types_strand1_medians_list, simulations_types_strand2_medians_list
##################################################################


##################################################################
#TODO Update the function March 29, 2020
#TODO left here fill results_dict
def fillSimulationsSample2Type2StrandCountList(
        simNum2Sample2Type2Strand2CountDict,
        existingTypesList,
        numberofSimulations,
        signatures_or_samples,
        strandbias):

    sample2SimulationsTypesStrand1MediansListDict = {}
    sample2SimulationsTypesStrand2MediansListDict = {}

    sample2AllSimulationsTypesStrand1ListDict = {}
    sample2AllSimulationsTypesStrand2ListDict = {}

    # Fill sample2AllSimulationsMutationTypesTranscribedListDict
    # Fill sample2AllSimulationsMutationTypesUntranscribedListDict
    ##############################################################
    for sample in signatures_or_samples:
        sample2AllSimulationsTypesStrand1ListDict[sample] = []
        sample2AllSimulationsTypesStrand2ListDict[sample] = []

        #Do not include simNum=0 since it is the original data
        for simNum in range(1, numberofSimulations+1):
            #new code
            sample2TypeStrand2CountDict = simNum2Sample2Type2Strand2CountDict[str(simNum)]

            if ((sample2TypeStrand2CountDict is not None) and (sample2TypeStrand2CountDict)):

                if (strandbias == TRANSCRIPTIONSTRANDBIAS):
                    strand2TypeStrandCountList_Dict = fillDictWRTReferenceTypes(existingTypesList,transcriptionStrands,sample2TypeStrand2CountDict[sample])
                    sampleBased_types_strand1_list = strand2TypeStrandCountList_Dict[TRANSCRIBED_STRAND]
                    sampleBased_types_strand2_list = strand2TypeStrandCountList_Dict[UNTRANSCRIBED_STRAND]
                elif (strandbias == REPLICATIONSTRANDBIAS):
                    strand2TypeStrandCountList_Dict = fillDictWRTReferenceTypes(existingTypesList,replicationStrands,sample2TypeStrand2CountDict[sample])
                    sampleBased_types_strand1_list = strand2TypeStrandCountList_Dict[LAGGING]
                    sampleBased_types_strand2_list = strand2TypeStrandCountList_Dict[LEADING]

                if (sampleBased_types_strand1_list):
                    sample2AllSimulationsTypesStrand1ListDict[sample].append(sampleBased_types_strand1_list)
                if (sampleBased_types_strand2_list):
                    sample2AllSimulationsTypesStrand2ListDict[sample].append(sampleBased_types_strand2_list)
            ##############################################################
        #################################################################################

    ##############################################################
    #After filling the all simulations list calculate the medians for each sample
    # Fill sample2SimulationsMutationTypesTranscribedMediansListDict
    # Fill sample2SimulationsMutationTypesUntranscribedMediansListDict
    ##############################################################
    for sample in signatures_or_samples :

        sample2SimulationsTypesStrand1MediansListDict[sample] = []
        sample2SimulationsTypesStrand2MediansListDict[sample] = []

        all_simulations_types_strand1_list = sample2AllSimulationsTypesStrand1ListDict[sample]
        all_simulations_types_strand2_list = sample2AllSimulationsTypesStrand2ListDict[sample]

        ##################################################################
        if ((all_simulations_types_strand1_list is not None) and all_simulations_types_strand1_list):
            all_simulations_types_strand1_nparray = np.vstack(all_simulations_types_strand1_list)
            (rows, cols) = all_simulations_types_strand1_nparray.shape

            for col in range(cols):
                colwise_array = all_simulations_types_strand1_nparray[:, col]
                sorted_colwise_array = np.sort(colwise_array)
                sample2SimulationsTypesStrand1MediansListDict[sample].append(np.median(sorted_colwise_array))
        ##################################################################

        ##################################################################
        if ((all_simulations_types_strand2_list is not None) and all_simulations_types_strand2_list):
            all_simulations_types_strand2_nparray = np.vstack(all_simulations_types_strand2_list)
            (rows, cols) = all_simulations_types_strand2_nparray.shape

            for col in range(cols):
                colwise_array = all_simulations_types_strand2_nparray[:, col]
                sorted_colwise_array = np.sort(colwise_array)
                sample2SimulationsTypesStrand2MediansListDict[sample].append(np.median(sorted_colwise_array))
        ##################################################################
    ##############################################################

    return sample2SimulationsTypesStrand1MediansListDict, sample2SimulationsTypesStrand2MediansListDict
##################################################################

##################################################################
def calculate_p_values(types_strand1_list,simulations_types_strand1_list,types_strand2_list,simulations_types_strand2_list):
    types_strandbias_pvalues = []

    #If there are no simulations case
    if ((simulations_types_strand1_list is None) and  (simulations_types_strand2_list is None)):
        simulations_types_strand1_list = [(x + y) / 2 for x, y in zip(types_strand1_list, types_strand2_list)]
        simulations_types_strand2_list = simulations_types_strand1_list

    for count1, count1_simulations, count2, count2_simulations in zip(types_strand1_list,
                                                                      simulations_types_strand1_list,
                                                                      types_strand2_list,
                                                                      simulations_types_strand2_list):
        #Is this true? Yes, it is correct.
        # we compare whether there is  a significance difference between the counts
        # namely, counts coming from the original data and counts coming from the simulations
        #if pValue < 0.05 then the the lower or higher count that comes from original data is statistically significant from the count that comes from teh simulations
        contingency_table_array = [[count1, count1_simulations], [count2, count2_simulations]]

        if ((count1 < 3000000) and (count2 < 3000000) and (count1_simulations < 3000000) and (count2_simulations < 3000000)):
            oddsratio, pvalue_SBS = stats.fisher_exact(contingency_table_array)
        else:
            chi2, pvalue_SBS, dof, expected = stats.chi2_contingency(contingency_table_array)
        types_strandbias_pvalues.append(pvalue_SBS)

    return types_strandbias_pvalues
##################################################################


##################################################################
#Fill the dictionaries that are required for the bar plots
#key can be signature and sample
#type can  be mutation type or signature when key is sample
def fillPValuesDictionaries(strands,
                            types,
                            key2Type2Strand2CountDict,
                            key2SimulationsTypesStrand1MediansListDict,
                            key2SimulationsTypesStrand2MediansListDict):

    key2TypesStrand1CountListDict = {}
    key2TypesStrand2CountListDict = {}
    key2TypesStrandPValuesListDict = {}

    for key in key2Type2Strand2CountDict:
        strand2TypesStrandCountListDict = fillDictWRTReferenceTypes(types,strands,key2Type2Strand2CountDict[key])

        types_strand1_list = strand2TypesStrandCountListDict[strands[0]]
        types_strand2_list = strand2TypesStrandCountListDict[strands[1]]

        key2TypesStrand1CountListDict[key] = types_strand1_list
        key2TypesStrand2CountListDict[key] = types_strand2_list

        simulations_types_strand1_list = None
        simulations_types_strand2_list = None
        #Now we have
        if ((key2SimulationsTypesStrand1MediansListDict is not None) and (key in key2SimulationsTypesStrand1MediansListDict)):
            simulations_types_strand1_list = key2SimulationsTypesStrand1MediansListDict[key]

        if  ((key2SimulationsTypesStrand2MediansListDict is not None ) and (key in key2SimulationsTypesStrand2MediansListDict)):
            simulations_types_strand2_list = key2SimulationsTypesStrand2MediansListDict[key]


        # Please notice Fisher exact test does not return anything for large numbers therefore chi2 test is used.
        print(key)
        print(types)
        print('for %s types_strand1_list: %s' %(strands[0],types_strand1_list))
        print('for %s types_strand2_list: %s' %(strands[1],types_strand2_list))
        print('for %s simulations_types_strand1_list: %s' % (strands[0],simulations_types_strand1_list))
        print('for %s simulations_types_strand2_list: %s' % (strands[1],simulations_types_strand2_list))
        types_strand_pvalues = calculate_p_values(types_strand1_list,
                                                simulations_types_strand1_list,
                                                types_strand2_list,
                                                simulations_types_strand2_list)
        print('types_strand_pvalues:%s' %types_strand_pvalues)

        key2TypesStrandPValuesListDict[key] = types_strand_pvalues

    return key2TypesStrand1CountListDict, key2TypesStrand2CountListDict, key2TypesStrandPValuesListDict
##################################################################



##################################################################
#Key can be signature or sample
def plotBarPlots(outputDir,jobname,numberofSimulations,signature_cutoff_numberofmutations_averageprobability_df,isKeySample,
                 existingMutationTypesList,
                 key2SimulationsMutationTypes_Strand1_MediansListDict,key2SimulationsMutationTypes_Strand2_MediansListDict,
                 key2MutationType2Strand2CountDict,
                 key2MutationTypesStrand1CountListDict,key2MutationTypesStrand2CountListDict,key2MutationTypes_BH_FDR_Adjusted_PValues_Dict,width,strands,color1,color2,title,figureName):

    #######################################################
    for key in key2MutationType2Strand2CountDict:
        mutationtype_strand1_list = key2MutationTypesStrand1CountListDict[key]
        mutationtype_strand2_list = key2MutationTypesStrand2CountListDict[key]

        x_axis_labels = existingMutationTypesList
        N = len(x_axis_labels)

        if ((key in signature_cutoff_numberofmutations_averageprobability_df['signature'].unique()) and (key in key2MutationTypes_BH_FDR_Adjusted_PValues_Dict)):
            sample_mutationtypes_FDR_BH_adjusted_pvalues = key2MutationTypes_BH_FDR_Adjusted_PValues_Dict[key]
            #[signature cutoffProbability numberofMutations averageProbability]
            numberofMutations = int(signature_cutoff_numberofmutations_averageprobability_df[signature_cutoff_numberofmutations_averageprobability_df['signature']==key]['number_of_mutations'].values[0])

            samplebased_simulations_mutationtypes_strand1_medians_list= None
            samplebased_simulations_mutationtypes_strand2_medians_list = None

            if ((key2SimulationsMutationTypes_Strand1_MediansListDict is not None) and (key in key2SimulationsMutationTypes_Strand1_MediansListDict)):
                samplebased_simulations_mutationtypes_strand1_medians_list = key2SimulationsMutationTypes_Strand1_MediansListDict[key]

            if ((key2SimulationsMutationTypes_Strand2_MediansListDict is not None) and (key in key2SimulationsMutationTypes_Strand2_MediansListDict)):
                samplebased_simulations_mutationtypes_strand2_medians_list = key2SimulationsMutationTypes_Strand2_MediansListDict[key]

            plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,key,isKeySample,numberofMutations,N,x_axis_labels,mutationtype_strand1_list,mutationtype_strand2_list,samplebased_simulations_mutationtypes_strand1_medians_list,samplebased_simulations_mutationtypes_strand2_medians_list,sample_mutationtypes_FDR_BH_adjusted_pvalues,strands[0],strands[1],title,color1,color2,figureName,width)
    #######################################################

##################################################################

##################################################################
def fillSignificantStrandDict(type2Strand2CountDict,type2FDR_BH_Adjusted_PValue_Dict):
    type2SignificantStrandDict={}

    for type,fdr_bh_adjusted_pvalue in type2FDR_BH_Adjusted_PValue_Dict.items():
        strand_saved=None
        count_saved=0
        if (fdr_bh_adjusted_pvalue<=0.05):
            strand2CountDict = type2Strand2CountDict[type]
            for strand, count in strand2CountDict.items():
                #We do not consider NonTranscribed for transcription
                if ((strand!='NonTranscribed') and (count>count_saved)):
                    count_saved=count
                    strand_saved=strand
        type2SignificantStrandDict[type]=strand_saved

    return type2SignificantStrandDict
##################################################################


##################################################################
def fillMutationTypesSignificantStrandDict(subsSignature2MutationType2Strand2CountDict,subsSignature2MutationType_Strand_FDR_BH_Adjusted_PValues_Dict,sixMutationTypes):
    subsSignature2MutationTypesSignificantTranscriptionStrandsListDict={}

    for subSignature,FDR_BH_Adjusted_PValues_List in subsSignature2MutationType_Strand_FDR_BH_Adjusted_PValues_Dict.items():
        subsSignature2MutationTypesSignificantTranscriptionStrandsListDict[subSignature]=[]
        for i, fdr_bh_adjusted_pvalue in enumerate(FDR_BH_Adjusted_PValues_List,0):
            mutation_type=sixMutationTypes[i]
            if  fdr_bh_adjusted_pvalue<=0.05:
                strand_saved = None
                count_saved = 0
                strand2CountDict=subsSignature2MutationType2Strand2CountDict[subSignature][mutation_type]
                for strand, count in strand2CountDict.items():
                    if ((strand!='NonTranscribed') and (count>count_saved)):
                        count_saved=count
                        strand_saved=strand
                subsSignature2MutationTypesSignificantTranscriptionStrandsListDict[subSignature].append(strand_saved)
            else:
                subsSignature2MutationTypesSignificantTranscriptionStrandsListDict[subSignature].append(None)

    return subsSignature2MutationTypesSignificantTranscriptionStrandsListDict
##################################################################


##################################################################
def  fillPValuesDict(type2StrandBiasPValueDict, typesPValuesList,types):
    for type_index, my_type in enumerate(types):
        if (type_index<len(typesPValuesList)):
            type2StrandBiasPValueDict[my_type]=typesPValuesList[type_index]
##################################################################


###################################################################
# April 20, 2020
def togetherTranscriptionReplicationStrandBiasFigures(outputDir,jobname,figureAugmentation,numberofSimulations,sample_based):
    pass
###################################################################


###################################################################
# April 20, 2020
def transcriptionStrandBiasFigures(outputDir,jobname,figureAugmentation,numberofSimulations,sample_based):
    pass
###################################################################


###################################################################
# April 20, 2020
#Part 1 Transcription and Replication Together
#Part 2 Transcription and Replication Together
#Part 3 Transcription and Replication Together
#TODO how must it be ? Then do like that
# It can not read I guess.
def replicationStrandBiasFigures(outputDir,jobname,figureAugmentation,numberofSimulations,sample_based):

    jobnamePath = os.path.join(outputDir,jobname,FIGURE,ALL,STRANDBIAS)
    print('Topography.py jobnamePath:%s ' %jobnamePath)

    ############################################################
    if (os.path.exists(jobnamePath)):
        try:
            shutil.rmtree(jobnamePath)
        except OSError as e:
            print('Error: %s - %s.' % (e.filename, e.strerror))
    ############################################################

    isFigureAugmentation = False
    if (figureAugmentation == 'augmentation'):
        isFigureAugmentation = True

    #######################################################################################################################
    os.makedirs(os.path.join(outputDir, jobname, FIGURE, ALL, STRANDBIAS,SCATTERPLOTS), exist_ok=True)
    #######################################################################################################################

    ##########################################################################################
    #########################  Read dictionaries related with ################################
    #########################  signatures and samples starts  ################################
    ##########################################################################################
    subsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(os.path.join(outputDir, jobname, DATA,Table_SubsSignature_Cutoff_NumberofMutations_AverageProbability_Filename),sep='\t', header=0,dtype={'cutoff': np.float32,'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})
    indelsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(os.path.join(outputDir, jobname, DATA,Table_IndelsSignature_Cutoff_NumberofMutations_AverageProbability_Filename),sep='\t', header=0,dtype={'cutoff': np.float32,'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})
    dinucsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(os.path.join(outputDir, jobname, DATA,Table_DinucsSignature_Cutoff_NumberofMutations_AverageProbability_Filename),sep='\t', header=0,dtype={'cutoff': np.float32,'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})

    subsSignatures = subsSignature_cutoff_numberofmutations_averageprobability_df['signature'].unique()
    indelsSignatures = indelsSignature_cutoff_numberofmutations_averageprobability_df['signature'].unique()
    dinucsSignatures = dinucsSignature_cutoff_numberofmutations_averageprobability_df['signature'].unique()

    if sample_based:
        sample2NumberofSubsDict = getSample2NumberofSubsDict(outputDir, jobname)
        sample2NumberofIndelsDict = getSample2NumberofIndelsDict(outputDir, jobname)
        sample2NumberofDinucsDict = getDictionary(outputDir, jobname, Sample2NumberofDinucsDictFilename)
    else:
        sample2NumberofSubsDict = {}
        sample2NumberofIndelsDict = {}
        sample2NumberofDinucsDict = {}
    ##########################################################################################
    #########################  Read dictionaries related with ################################
    #########################  signatures and samples ends  ##################################
    ##########################################################################################

    ##########################################################################################
    #####################  Read dictionaries for the rest part starts ########################
    ##########################################################################################

    #########################   Replication starts   #########################################
    Type2ReplicationStrand2CountDict_Filepath = os.path.join(outputDir,jobname,DATA,REPLICATIONSTRANDBIAS,Type2ReplicationStrand2CountDict_Filename)
    simNum2Type2ReplicationStrand2CountDict = readDictionary(Type2ReplicationStrand2CountDict_Filepath)
    type2ReplicationStrand2CountDict = simNum2Type2ReplicationStrand2CountDict[str(0)]

    #These are subs signatures
    Signature2MutationType2ReplicationStrand2CountDict_Filepath = os.path.join(outputDir,jobname,DATA,REPLICATIONSTRANDBIAS,Signature2MutationType2ReplicationStrand2CountDict_Filename)
    simNum2SubsSignature2MutationType2ReplicationStrand2CountDict = readDictionary(Signature2MutationType2ReplicationStrand2CountDict_Filepath)
    subsSignature2MutationType2ReplicationStrand2CountDict = simNum2SubsSignature2MutationType2ReplicationStrand2CountDict[str(0)]

    if sample_based:
        Sample2Type2ReplicationStrand2CountDict_Filepath = os.path.join(outputDir,jobname,DATA,REPLICATIONSTRANDBIAS,Sample2Type2ReplicationStrand2CountDict_Filename)
        simNum2Sample2Type2ReplicationStrand2CountDict = readDictionary(Sample2Type2ReplicationStrand2CountDict_Filepath)
        sample2Type2ReplicationStrand2CountDict = simNum2Sample2Type2ReplicationStrand2CountDict[str(0)]

        Type2Sample2ReplicationStrand2CountDict_Filepath = os.path.join(outputDir,jobname,DATA,REPLICATIONSTRANDBIAS,Type2Sample2ReplicationStrand2CountDict_Filename)
        simNum2Type2Sample2ReplicationStrand2CountDict = readDictionary(Type2Sample2ReplicationStrand2CountDict_Filepath)
        type2Sample2ReplicationStrand2CountDict =  simNum2Type2Sample2ReplicationStrand2CountDict[str(0)]
    #########################   Replication ends   ###########################################

    ##########################################################################################
    #####################  Read dictionaries for the rest part ends ##########################
    ##########################################################################################


    ########################################################################
    ##########################  Part 4 starts ##############################
    ######## Bar plot starts includes sample based bar plots ###############
    ########################################################################
    # Step1: We need to calculate the p-values for each mutation type and signature using Fisher's exact test
    # Step2: We need to calculate FDR adjusted p-values for all mutation types and signatures.
    ########################################################################

    ########################################################################
    ############ MutationType based bar plot starts ########################
    ############ Signature based bar plot starts ###########################
    ########################################################################

    ##########################################################################################################
    ################################### Fill for simulations starts ##########################################
    ##########################################################################################################
    #######
    simulations_mutationtypes_lagging_medians_list = None
    simulations_mutationtypes_leading_medians_list = None
    sample2SimulationsMutationTypesLaggingMediansListDict = None
    sample2SimulationsMutationTypesLeadingMediansListDict = None
    #######

    #######
    simulations_subs_signatures_lagging_medians_list = None
    simulations_subs_signatures_leading_medians_list = None

    sample2SimulationsSubsSignatureTypesLaggingMediansListDict = None
    sample2SimulationsSubsSignatureTypesLeadingMediansListDict = None
    #######

    #######
    simulations_indels_signatures_lagging_medians_list = None
    simulations_indels_signatures_leading_medians_list = None

    sample2SimulationsIndelsSignatureTypesLaggingMediansListDict = None
    sample2SimulationsIndelsSignatureTypesLeadingMediansListDict = None
    #######

    #######
    simulations_dinucs_signatures_lagging_medians_list = None
    simulations_dinucs_signatures_leading_medians_list = None

    sample2SimulationsDinucsSignatureTypesLaggingMediansListDict = None
    sample2SimulationsDinucsSignatureTypesLeadingMediansListDict = None
    #######

    #######
    subsSignature2SimulationsMutationTypesLaggingMediansListDict = None
    subsSignature2SimulationsMutationTypesLeadingMediansListDict = None
    #######

    if (numberofSimulations > 0):

        #########################################################################################################################
        # simulations --- mutation type --- replication
        simulations_mutationtypes_lagging_medians_list, simulations_mutationtypes_leading_medians_list =\
            fillSimulationsType2StrandCountList(simNum2Type2ReplicationStrand2CountDict,
                                                sixMutationTypes,
                                                numberofSimulations,
                                                REPLICATIONSTRANDBIAS,)

        # simulations --- subs signatures --- replication
        simulations_subs_signatures_lagging_medians_list, simulations_subs_signatures_leading_medians_list = \
            fillSimulationsType2StrandCountList(simNum2Type2ReplicationStrand2CountDict,
                                                subsSignatures,
                                                numberofSimulations,
                                                REPLICATIONSTRANDBIAS)

        # simulations --- indels signatures --- replication
        simulations_indels_signatures_lagging_medians_list, simulations_indels_signatures_leading_medians_list = \
            fillSimulationsType2StrandCountList(simNum2Type2ReplicationStrand2CountDict,
                                                indelsSignatures,
                                                numberofSimulations,
                                                REPLICATIONSTRANDBIAS)


        # simulations --- dinucs signatures --- replication
        simulations_dinucs_signatures_lagging_medians_list, simulations_dinucs_signatures_leading_medians_list = \
            fillSimulationsType2StrandCountList(simNum2Type2ReplicationStrand2CountDict,
                                                dinucsSignatures,
                                                numberofSimulations,
                                                REPLICATIONSTRANDBIAS)
        #########################################################################################################################

        #########################################################################################################################
        # simulations signature --- mutation type --- replication

        print('DEBUG1 simNum2SubsSignature2MutationType2ReplicationStrand2CountDict: %s' %(simNum2SubsSignature2MutationType2ReplicationStrand2CountDict))

        subsSignature2SimulationsMutationTypesLaggingMediansListDict, subsSignature2SimulationsMutationTypesLeadingMediansListDict = \
            fillSimulationsSample2Type2StrandCountList(simNum2SubsSignature2MutationType2ReplicationStrand2CountDict,
                                                sixMutationTypes,
                                                numberofSimulations,
                                                subsSignature_cutoff_numberofmutations_averageprobability_df['signature'].unique(),
                                                REPLICATIONSTRANDBIAS)
        #########################################################################################################################

        #########################################################################################################################
        if sample_based:
            # samplebased --- simulations --- mutation type --- replication
            sample2SimulationsMutationTypesLaggingMediansListDict, sample2SimulationsMutationTypesLeadingMediansListDict =\
                fillSimulationsSample2Type2StrandCountList(simNum2Sample2Type2ReplicationStrand2CountDict,
                    sixMutationTypes,
                    numberofSimulations,
                    sample2NumberofSubsDict.keys(),
                    REPLICATIONSTRANDBIAS)

            # samplebased --- simulations --- signature --- replication
            sample2SimulationsSubsSignatureTypesLaggingMediansListDict, sample2SimulationsSubsSignatureTypesLeadingMediansListDict =\
                fillSimulationsSample2Type2StrandCountList(simNum2Sample2Type2ReplicationStrand2CountDict,
                    subsSignatures,
                    numberofSimulations,
                    sample2NumberofSubsDict.keys(),
                    REPLICATIONSTRANDBIAS)

            # samplebased --- simulations --- signature --- replication
            sample2SimulationsIndelsSignatureTypesLaggingMediansListDict, sample2SimulationsIndelsSignatureTypesLeadingMediansListDict =\
                fillSimulationsSample2Type2StrandCountList(simNum2Sample2Type2ReplicationStrand2CountDict,
                    indelsSignatures,
                    numberofSimulations,
                    sample2NumberofIndelsDict.keys(),
                    REPLICATIONSTRANDBIAS)

            # samplebased --- simulations --- signature --- replication
            sample2SimulationsDinucsSignatureTypesLaggingMediansListDict, sample2SimulationsDinucsSignatureTypesLeadingMediansListDict =\
                fillSimulationsSample2Type2StrandCountList(simNum2Sample2Type2ReplicationStrand2CountDict,
                    dinucsSignatures,
                    numberofSimulations,
                    sample2NumberofDinucsDict.keys(),
                    REPLICATIONSTRANDBIAS)
        #########################################################################################################################

    ##########################################################################################################
    ################################### Fill for simulations ends ############################################
    ##########################################################################################################

    #TODO These dictionaries will be written
    type2TranscriptionPValueDict={}
    type2ReplicationPValueDict={}

    #########################################################################################################
    #########################  Calculate p-values for SubsSignature Mutation Type starts  ###################
    #########################################################################################################

    print('DEBUG subsSignature2MutationType2ReplicationStrand2CountDict:%s' %(subsSignature2MutationType2ReplicationStrand2CountDict))
    print('DEBUG subsSignature2SimulationsMutationTypesLaggingMediansListDict:%s' %(subsSignature2SimulationsMutationTypesLaggingMediansListDict))
    print('DEBUG subsSignature2SimulationsMutationTypesLeadingMediansListDict:%s' %(subsSignature2SimulationsMutationTypesLeadingMediansListDict))
    subsSignature2MutationTypesLaggingCountListDict, \
    subsSignature2MutationTypesLeadingCountListDict, \
    subsSignature2MutationTypesReplicationPValuesListDict = fillPValuesDictionaries(
        replicationStrands,
        sixMutationTypes,
        subsSignature2MutationType2ReplicationStrand2CountDict,
        subsSignature2SimulationsMutationTypesLaggingMediansListDict,
        subsSignature2SimulationsMutationTypesLeadingMediansListDict)
    #########################################################################################################
    #########################  Calculate p-values for SubsSignature Mutation Type ends  #####################
    #########################################################################################################

    #########################################################################################################
    #############################  Calculate p-values for Mutation Type Starts    ###########################
    #########################################################################################################

    ########################################## Replication starts  ##########################################
    strandType_MutationType2CountList_Dict = fillDictWRTReferenceTypes(sixMutationTypes, replicationStrands,type2ReplicationStrand2CountDict)
    mutationtypes_lagging_list = strandType_MutationType2CountList_Dict[LAGGING]
    mutationtypes_leading_list = strandType_MutationType2CountList_Dict[LEADING]

    # Please notice Fisher exact test does not return anything for large numbers therefore chi2 test is used.
    mutationtype_replication_pvalues = calculate_p_values(mutationtypes_lagging_list,
                                                          simulations_mutationtypes_lagging_medians_list,
                                                          mutationtypes_leading_list,
                                                          simulations_mutationtypes_leading_medians_list)

    fillPValuesDict(type2ReplicationPValueDict, mutationtype_replication_pvalues,sixMutationTypes)
    ########################################## Replication ends  ############################################
    #########################################################################################################
    #############################  Calculate p-values for Mutation Type ends    #############################
    #########################################################################################################

    #########################################################################################################
    ########################  Calculate p-values for subs signatures starts    ##############################
    #########################################################################################################

    ########################################## Replication starts  ##########################################
    strandType_Signature2CountList_Dict = fillDictWRTReferenceTypes(subsSignatures,replicationStrands,type2ReplicationStrand2CountDict)
    subs_signatures_laggingStrandCount_list = strandType_Signature2CountList_Dict[LAGGING]
    subs_signatures_leadingStrandCount_list = strandType_Signature2CountList_Dict[LEADING]

    subs_signature_replication_pvalues = calculate_p_values(subs_signatures_laggingStrandCount_list,
                                                       simulations_subs_signatures_lagging_medians_list,
                                                       subs_signatures_leadingStrandCount_list,
                                                       simulations_subs_signatures_leading_medians_list)

    fillPValuesDict(type2ReplicationPValueDict, subs_signature_replication_pvalues,subsSignatures)
    ########################################## Replication ends  ############################################
    #########################################################################################################
    ########################  Calculate p-values for subs signatures ends    ################################
    #########################################################################################################


    #########################################################################################################
    ######################  Calculate p-values for indels signatures starts  ################################
    #########################################################################################################

    ########################################## Replication starts  ##########################################
    strandType_Signature2CountList_Dict = fillDictWRTReferenceTypes(indelsSignatures,replicationStrands,type2ReplicationStrand2CountDict)
    indels_signatures_laggingStrandCount_list = strandType_Signature2CountList_Dict[LAGGING]
    indels_signatures_leadingStrandCount_list = strandType_Signature2CountList_Dict[LEADING]

    indels_signature_replication_pvalues = calculate_p_values(indels_signatures_laggingStrandCount_list,
                                                       simulations_indels_signatures_lagging_medians_list,
                                                       indels_signatures_leadingStrandCount_list,
                                                       simulations_indels_signatures_leading_medians_list)

    fillPValuesDict(type2ReplicationPValueDict, indels_signature_replication_pvalues,indelsSignatures)
    ########################################## Replication ends  ############################################
    #########################################################################################################
    ######################  Calculate p-values for indels signatures ends    ################################
    #########################################################################################################


    #########################################################################################################
    ######################  Calculate p-values for dinucs signatures starts  ################################
    #########################################################################################################

    ########################################## Replication starts  ##########################################
    strandType_Signature2CountList_Dict = fillDictWRTReferenceTypes(dinucsSignatures,replicationStrands,type2ReplicationStrand2CountDict)
    dinucs_signatures_laggingStrandCount_list = strandType_Signature2CountList_Dict[LAGGING]
    dinucs_signatures_leadingStrandCount_list = strandType_Signature2CountList_Dict[LEADING]

    dinucs_signature_replication_pvalues = calculate_p_values(dinucs_signatures_laggingStrandCount_list,
                                                        simulations_dinucs_signatures_lagging_medians_list,
                                                        dinucs_signatures_leadingStrandCount_list,
                                                        simulations_dinucs_signatures_leading_medians_list)

    fillPValuesDict(type2ReplicationPValueDict, dinucs_signature_replication_pvalues,dinucsSignatures)
    ########################################## Replication ends  ############################################

    #########################################################################################################
    ######################  Calculate p-values for dinucs signatures ends    ################################
    #########################################################################################################

    if sample_based:
        #########################################################################################################
        ####################  Calculate p-values for sample based mutation types starts    ######################
        #########################################    Transcription    ###########################################
        ###########################################    Replication    ###########################################
        #########################################################################################################
        sample2MutationTypesLaggingCountListDict, sample2MutationTypesLeadingCountListDict, sample2MutationTypesReplicationPValuesListDict = fillPValuesDictionaries(replicationStrands,sixMutationTypes,sample2Type2ReplicationStrand2CountDict,sample2SimulationsMutationTypesLaggingMediansListDict,sample2SimulationsMutationTypesLeadingMediansListDict)
        #########################################################################################################
        ####################  Calculate p-values for sample based mutation types ends    ########################
        #########################################    Transcription    ###########################################
        ###########################################    Replication    ###########################################
        #########################################################################################################


        #########################################################################################################
        ####################  Calculate p-values for sample based subs signature types starts    ################
        #########################################    Transcription    ###########################################
        ###########################################    Replication    ###########################################
        #########################################################################################################
        sample2SubsSignatureTypesLaggingCountListDict, sample2SubsSignatureTypesLeadingCountListDict, sample2SubsSignature2ReplicationPValuesListDict = \
            fillPValuesDictionaries(replicationStrands,
                                        subsSignatures,
                                        sample2Type2ReplicationStrand2CountDict,
                                        sample2SimulationsSubsSignatureTypesLaggingMediansListDict,
                                        sample2SimulationsSubsSignatureTypesLeadingMediansListDict)
        #########################################################################################################
        ####################  Calculate p-values for sample based subs signature types ends    ##################
        #########################################    Transcription    ###########################################
        ###########################################    Replication    ###########################################
        #########################################################################################################


        #########################################################################################################
        ####################  Calculate p-values for sample based indels signature types starts    ##############
        #########################################    Transcription    ###########################################
        ###########################################    Replication    ###########################################
        #########################################################################################################
        sample2IndelsSignatureTypesLaggingCountListDict, sample2IndelsSignatureTypesLeadingCountListDict, sample2IndelsSignature2ReplicationPValuesListDict = \
            fillPValuesDictionaries(replicationStrands,
                                        indelsSignatures,
                                        sample2Type2ReplicationStrand2CountDict,
                                        sample2SimulationsIndelsSignatureTypesLaggingMediansListDict,
                                        sample2SimulationsIndelsSignatureTypesLeadingMediansListDict)
        #########################################################################################################
        ####################  Calculate p-values for sample based indels signature types ends    ################
        #########################################    Transcription    ###########################################
        ###########################################    Replication    ###########################################
        #########################################################################################################


        #########################################################################################################
        ####################  Calculate p-values for sample based dinucs signature types starts    ##############
        #########################################    Transcription    ###########################################
        ###########################################    Replication    ###########################################
        #########################################################################################################
        sample2DinucsSignatureTypesLaggingCountListDict, sample2DinucsSignatureTypesLeadingCountListDict, sample2DinucsSignature2ReplicationPValuesListDict = \
            fillPValuesDictionaries(replicationStrands,
                                    dinucsSignatures,
                                    sample2Type2ReplicationStrand2CountDict,
                                    sample2SimulationsDinucsSignatureTypesLaggingMediansListDict,
                                    sample2SimulationsDinucsSignatureTypesLeadingMediansListDict)
        #########################################################################################################
        ####################  Calculate p-values for sample based indels signature types ends    ################
        #########################################    Transcription    ###########################################
        ###########################################    Replication    ###########################################
        #########################################################################################################


    ######################################################
    ################# FDR BH starts ######################
    ######################################################
    #Order the samples
    #We must put the samples in the order they are put in all_p_values list
    all_p_values = []
    all_p_values.extend(mutationtype_replication_pvalues)
    all_p_values.extend(subs_signature_replication_pvalues)
    all_p_values.extend(indels_signature_replication_pvalues)
    all_p_values.extend(dinucs_signature_replication_pvalues)

    for subsSignature in subsSignatures:
        if subsSignature in subsSignature2MutationTypesReplicationPValuesListDict:
            all_p_values.extend(subsSignature2MutationTypesReplicationPValuesListDict[subsSignature])

    if sample_based:
        #############################################################################################
        ##################################  Sample Based starts  ####################################
        #############################################################################################
        sampleMutationTypesReplicationList = sample2MutationTypesReplicationPValuesListDict.keys()
        sampleSubsSignaturesReplicationList = sample2SubsSignature2ReplicationPValuesListDict.keys()
        sampleIndelsSignaturesReplicationList = sample2IndelsSignature2ReplicationPValuesListDict.keys()
        sampleDinucsSignaturesReplicationList = sample2DinucsSignature2ReplicationPValuesListDict.keys()

        for sample in sampleMutationTypesReplicationList:
            all_p_values.extend(sample2MutationTypesReplicationPValuesListDict[sample])

        for sample in sampleSubsSignaturesReplicationList:
            all_p_values.extend(sample2SubsSignature2ReplicationPValuesListDict[sample])

        for sample in sampleIndelsSignaturesReplicationList:
            all_p_values.extend(sample2IndelsSignature2ReplicationPValuesListDict[sample])

        for sample in sampleDinucsSignaturesReplicationList:
            all_p_values.extend(sample2DinucsSignature2ReplicationPValuesListDict[sample])
        #############################################################################################
        ##################################  Sample Based ends  ######################################
        #############################################################################################

    print('For information: Number of p values to be used in multiple testing: %d' %(len(all_p_values)))
    print('all_p_values')
    print(all_p_values)
    print('len(all_p_values)')
    print(len(all_p_values))
    print('######################################################')

    ####################################################################################
    mutationtype_replication_FDR_BH_adjusted_pvalues = []
    subsSignaturetype_replication_FDR_BH_adjusted_pvalues = []
    indelsSignaturetype_replication_FDR_BH_adjusted_pvalues = []
    dinucsSignaturetype_replication_FDR_BH_adjusted_pvalues = []

    #These dictionaries will be written
    type2Replication_FDR_BH_Adjusted_PValue_Dict={}
    subsSignature2MutationTypesReplication_FDR_BH_Adjusted_PValues_Dict = {}

    if ((all_p_values is  not None) and all_p_values):
        all_p_values_array = np.asarray(all_p_values)

        # rejected, all_FDR_BH_adjusted_p_values = statsmodels.stats.multitest.fdrcorrection(all_p_values_array, alpha=0.05, method='indep', is_sorted=False)
        rejected, all_FDR_BH_adjusted_p_values, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(all_p_values_array, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
        print('For information: Number of all_FDR_BH_adjusted_p_values: %d' % (len(all_FDR_BH_adjusted_p_values)))

        ################################################################
        #Replication MutationTypes
        start = 0
        end = len(sixMutationTypes)
        mutationtype_replication_FDR_BH_adjusted_pvalues = (all_FDR_BH_adjusted_p_values[start:end]).tolist()

        startIndex=start
        for mutationType in sixMutationTypes:
            type2Replication_FDR_BH_Adjusted_PValue_Dict[mutationType]=all_FDR_BH_adjusted_p_values[startIndex]
            startIndex+=1
        ################################################################

        ################################################################
        #Replication subsSignatures
        start = end
        end += len(subsSignatures)
        subsSignaturetype_replication_FDR_BH_adjusted_pvalues = all_FDR_BH_adjusted_p_values[start:end]

        startIndex=start
        for signature in subsSignatures:
            type2Replication_FDR_BH_Adjusted_PValue_Dict[signature] = all_FDR_BH_adjusted_p_values[startIndex]
            startIndex+=1
        ################################################################

        ################################################################
        #Replication indelsSignatures
        start = end
        end += len(indelsSignatures)
        indelsSignaturetype_replication_FDR_BH_adjusted_pvalues = all_FDR_BH_adjusted_p_values[start:end]

        startIndex=start
        for signature in indelsSignatures:
            type2Replication_FDR_BH_Adjusted_PValue_Dict[signature] = all_FDR_BH_adjusted_p_values[startIndex]
            startIndex+=1
        ################################################################

        ################################################################
        #Replication dinucsSignatures
        start = end
        end += len(dinucsSignatures)
        dinucsSignaturetype_replication_FDR_BH_adjusted_pvalues = all_FDR_BH_adjusted_p_values[start:end]

        startIndex=start
        for signature in dinucsSignatures:
            type2Replication_FDR_BH_Adjusted_PValue_Dict[signature] = all_FDR_BH_adjusted_p_values[startIndex]
            startIndex+=1
        ################################################################

        ################################################################
        for subsSignature in subsSignatures:
            if subsSignature in subsSignature2MutationTypesReplicationPValuesListDict:
                start = end
                end += len(sixMutationTypes)
                subsSignature2MutationTypesReplication_FDR_BH_Adjusted_PValues_Dict[subsSignature] = (all_FDR_BH_adjusted_p_values[start:end]).tolist()
        ################################################################

        print('######################################################')
        path = os.path.join(outputDir, jobname, FIGURE, ALL, STRANDBIAS)

        #Write Raw P Values dictionaries starts
        filename = 'type2ReplicationRawPValueDict.txt'
        writeDictionarySimple(type2ReplicationPValueDict, path, filename, None)

        filename = 'subsSignature2MutationTypesReplicationPValuesListDict.txt'
        writeDictionarySimple(subsSignature2MutationTypesReplicationPValuesListDict, path, filename, None)
        #Write Raw P Values dictionaries ends

        #type dictionaries
        filename = 'type2Replication_FDR_BH_Adjusted_PValue_Dict.txt'
        writeDictionarySimple(type2Replication_FDR_BH_Adjusted_PValue_Dict, path, filename, None)

        #subsSignatures2MutationTypes dictionaries
        filename = 'subsSignature2MutationTypesReplication_BH_FDR_Adjusted_PValues_Dict.txt'
        writeDictionarySimple(subsSignature2MutationTypesReplication_FDR_BH_Adjusted_PValues_Dict, path, filename, None)
        print('######################################################')

        type2SignificantReplicationStrandDict=fillSignificantStrandDict(type2ReplicationStrand2CountDict,type2Replication_FDR_BH_Adjusted_PValue_Dict)
        subsSignature2MutationTypesSignificantReplicationStrandsListDict=fillMutationTypesSignificantStrandDict(subsSignature2MutationType2ReplicationStrand2CountDict,subsSignature2MutationTypesReplication_FDR_BH_Adjusted_PValues_Dict,sixMutationTypes)

        print('######################################################')
        path = os.path.join(outputDir, jobname, FIGURE, ALL, STRANDBIAS)

        #type dictionaries
        filename = 'type2SignificantReplicationStrandDict.txt'
        writeDictionarySimple(type2SignificantReplicationStrandDict, path, filename, None)

        #subsSignatures2MutationTypes dictionaries
        filename = 'subsSignature2MutationTypesSignificantReplicationStrandsListDict.txt'
        writeDictionarySimple(subsSignature2MutationTypesSignificantReplicationStrandsListDict, path, filename, None)
        print('######################################################')

        print('######################################################')
        print('FDR BH adjusted p values starts')
        print('all_FDR_BH_adjusted_p_values')
        print(all_FDR_BH_adjusted_p_values)
        print('len(all_FDR_BH_adjusted_p_values)')
        print(len(all_FDR_BH_adjusted_p_values))
        print('FDR BH adjusted p values ends')
        print('######################################################')


        if sample_based:
            #We must get the BH FDR adjusted p values for each sample as in the order they are put in the all_p_values
            #Fill these dictionaries
            sample2MutationTypesReplication_BH_FDR_Adjusted_PValues_Dict = {}
            sample2SubsSignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict = {}
            sample2IndelsSignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict = {}
            sample2DinucsSignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict = {}


            # Sample MutationType Replication
            for sample in sampleMutationTypesReplicationList:
                start = end
                end += len(sixMutationTypes)
                sample2MutationTypesReplication_BH_FDR_Adjusted_PValues_Dict[sample] = all_FDR_BH_adjusted_p_values[start:end]

            # Sample Subs SignatureType Replication
            for sample in sampleSubsSignaturesReplicationList:
                start = end
                end += len(subsSignatures)
                sample2SubsSignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict[sample] = all_FDR_BH_adjusted_p_values[start:end]

            # Sample Indels SignatureType Replication
            for sample in sampleIndelsSignaturesReplicationList:
                start = end
                end += len(indelsSignatures)
                sample2IndelsSignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict[sample] = all_FDR_BH_adjusted_p_values[start:end]

            # Sample Dinucs SignatureType Replication
            for sample in sampleDinucsSignaturesReplicationList:
                start = end
                end += len(dinucsSignatures)
                sample2DinucsSignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict[sample] = all_FDR_BH_adjusted_p_values[start:end]
    ####################################################################################

    ######################################################
    ################# FDR BH ends ########################
    ######################################################

    isKeySample = False
    #######################################################
    ################# Plot mutations starts ###############
    #######################################################
    width = 0.20

    # Replication
    x_axis_labels = sixMutationTypes
    N = len(x_axis_labels)
    plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,None,isKeySample,None,N, x_axis_labels, mutationtypes_lagging_list, mutationtypes_leading_list,simulations_mutationtypes_lagging_medians_list, simulations_mutationtypes_leading_medians_list,mutationtype_replication_FDR_BH_adjusted_pvalues, 'Lagging', 'Leading','All Mutations', 'indianred', 'goldenrod', 'mutationtypes_replication_strand_bias',width)
    #######################################################
    ################# Plot mutations ends #################
    #######################################################

    #######################################################
    ########### Plot subs signatures starts ###############
    #######################################################
    # Replication
    x_axis_replication_labels = subsSignatures
    N_signatures_replication = len(x_axis_replication_labels)

    if (np.any(x_axis_replication_labels) and subs_signatures_laggingStrandCount_list and subs_signatures_leadingStrandCount_list and simulations_subs_signatures_lagging_medians_list and simulations_subs_signatures_leading_medians_list and (subsSignaturetype_replication_FDR_BH_adjusted_pvalues.size>0)):
        zipped = zip(x_axis_replication_labels, subs_signatures_laggingStrandCount_list, subs_signatures_leadingStrandCount_list, simulations_subs_signatures_lagging_medians_list, simulations_subs_signatures_leading_medians_list, subsSignaturetype_replication_FDR_BH_adjusted_pvalues)
        zipped = sorted(zipped, key=lambda x: natural_key(x[0]))
        x_axis_replication_labels, subs_signatures_laggingStrandCount_list, subs_signatures_leadingStrandCount_list, simulations_subs_signatures_lagging_medians_list, simulations_subs_signatures_leading_medians_list, subsSignaturetype_replication_FDR_BH_adjusted_pvalues = zip(*zipped)

        if (subs_signatures_laggingStrandCount_list and subs_signatures_leadingStrandCount_list):
            plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,None,isKeySample,None,N_signatures_replication,x_axis_replication_labels,subs_signatures_laggingStrandCount_list,subs_signatures_leadingStrandCount_list, simulations_subs_signatures_lagging_medians_list, simulations_subs_signatures_leading_medians_list,subsSignaturetype_replication_FDR_BH_adjusted_pvalues, 'Lagging','Leading', 'All Signatures', 'indianred', 'goldenrod','subs_signatures_replication_strand_bias', width)
    #######################################################
    ########### Plot subs signatures ends #################
    #######################################################

    #######################################################
    ########### Plot indels signatures starts #############
    #######################################################
    # Replication
    x_axis_replication_labels = indelsSignatures
    N_signatures_replication = len(x_axis_replication_labels)

    if (np.any(x_axis_replication_labels) and indels_signatures_laggingStrandCount_list and indels_signatures_leadingStrandCount_list and simulations_indels_signatures_lagging_medians_list and simulations_indels_signatures_leading_medians_list and (indelsSignaturetype_replication_FDR_BH_adjusted_pvalues.size>0)):
        zipped = zip(x_axis_replication_labels, indels_signatures_laggingStrandCount_list, indels_signatures_leadingStrandCount_list, simulations_indels_signatures_lagging_medians_list, simulations_indels_signatures_leading_medians_list, indelsSignaturetype_replication_FDR_BH_adjusted_pvalues)
        zipped = sorted(zipped, key=lambda x: natural_key(x[0]))
        x_axis_replication_labels, indels_signatures_laggingStrandCount_list, indels_signatures_leadingStrandCount_list, simulations_indels_signatures_lagging_medians_list, simulations_indels_signatures_leading_medians_list, indelsSignaturetype_replication_FDR_BH_adjusted_pvalues = zip(*zipped)

        if (indels_signatures_laggingStrandCount_list and indels_signatures_leadingStrandCount_list):
            plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,None,isKeySample,None,N_signatures_replication,x_axis_replication_labels,indels_signatures_laggingStrandCount_list,indels_signatures_leadingStrandCount_list, simulations_indels_signatures_lagging_medians_list, simulations_indels_signatures_leading_medians_list,indelsSignaturetype_replication_FDR_BH_adjusted_pvalues, 'Lagging','Leading', 'All Signatures', 'indianred', 'goldenrod','indels_signatures_replication_strand_bias', width)
    #######################################################
    ########### Plot indels signatures ends ###############
    #######################################################

    #######################################################
    ########### Plot dinucs signatures starts #############
    #######################################################
    # Replication
    x_axis_replication_labels = dinucsSignatures
    N_signatures_replication = len(x_axis_replication_labels)

    if (np.any(x_axis_replication_labels) and dinucs_signatures_laggingStrandCount_list and dinucs_signatures_leadingStrandCount_list and simulations_dinucs_signatures_lagging_medians_list and simulations_dinucs_signatures_leading_medians_list and (dinucsSignaturetype_replication_FDR_BH_adjusted_pvalues.size>0)):
        zipped = zip(x_axis_replication_labels, dinucs_signatures_laggingStrandCount_list, dinucs_signatures_leadingStrandCount_list, simulations_dinucs_signatures_lagging_medians_list, simulations_dinucs_signatures_leading_medians_list, dinucsSignaturetype_replication_FDR_BH_adjusted_pvalues)
        zipped = sorted(zipped, key=lambda x: natural_key(x[0]))
        x_axis_replication_labels, dinucs_signatures_laggingStrandCount_list, dinucs_signatures_leadingStrandCount_list, simulations_dinucs_signatures_lagging_medians_list, simulations_dinucs_signatures_leading_medians_list, dinucsSignaturetype_replication_FDR_BH_adjusted_pvalues = zip(*zipped)

        if (dinucs_signatures_laggingStrandCount_list and dinucs_signatures_leadingStrandCount_list):
            plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,None,isKeySample,None,N_signatures_replication,x_axis_replication_labels,dinucs_signatures_laggingStrandCount_list,dinucs_signatures_leadingStrandCount_list, simulations_dinucs_signatures_lagging_medians_list, simulations_dinucs_signatures_leading_medians_list,dinucsSignaturetype_replication_FDR_BH_adjusted_pvalues, 'Lagging','Leading', 'All Signatures', 'indianred', 'goldenrod','dinucs_signatures_replication_strand_bias', width)
    #######################################################
    ########### Plot dinucs signatures ends ###############
    #######################################################


    #################################################################
    ########### Plot sub signatures mutation types starts ###########
    #################################################################
    plotBarPlots(outputDir,jobname,numberofSimulations,subsSignature_cutoff_numberofmutations_averageprobability_df,isKeySample,sixMutationTypes,
                 subsSignature2SimulationsMutationTypesLaggingMediansListDict,
                 subsSignature2SimulationsMutationTypesLeadingMediansListDict,
                 subsSignature2MutationType2ReplicationStrand2CountDict,
                 subsSignature2MutationTypesLaggingCountListDict,
                 subsSignature2MutationTypesLeadingCountListDict,
                 subsSignature2MutationTypesReplication_FDR_BH_Adjusted_PValues_Dict,
                 width,replicationStrands,'indianred', 'goldenrod','All Mutations','mutationtypes_replication_strand_bias')
    #################################################################
    ########### Plot sub signatures mutation types ends #############
    #################################################################


    isKeySample=True
    if sample_based:
        #######################################################
        ######### Plot sample based mutations starts ##########
        #######################################################
        plotBarPlots(outputDir,jobname,numberofSimulations,sample2NumberofSubsDict,isKeySample,sixMutationTypes,
                     sample2SimulationsMutationTypesLaggingMediansListDict,
                     sample2SimulationsMutationTypesLeadingMediansListDict,
                     sample2Type2ReplicationStrand2CountDict,
                     sample2MutationTypesLaggingCountListDict,
                     sample2MutationTypesLeadingCountListDict,
                     sample2MutationTypesReplication_BH_FDR_Adjusted_PValues_Dict,
                     width,replicationStrands,'indianred', 'goldenrod','All Mutations','mutationtypes_replication_strand_bias')
        #######################################################
        ######### Plot sample based mutations ends ############
        #######################################################


        ############################################################
        #### Plot sample based subs signatures bar plot starts #####
        ############################################################
        plotBarPlots(outputDir,jobname,numberofSimulations,sample2NumberofSubsDict,isKeySample,subsSignatures,
                     sample2SimulationsSubsSignatureTypesLaggingMediansListDict,
                     sample2SimulationsSubsSignatureTypesLeadingMediansListDict,
                     sample2Type2ReplicationStrand2CountDict,
                     sample2SubsSignatureTypesLaggingCountListDict,
                     sample2SubsSignatureTypesLeadingCountListDict,
                     sample2SubsSignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict,
                     width,replicationStrands,'indianred', 'goldenrod','All Signatures','subs_signatures_replication_strand_bias')
        ############################################################
        #### Plot sample based subs signatures bar plot ends #######
        ############################################################

        ########################################################## New Part starts ##################################################################33
        ##############################################################
        #### Plot sample based indels signatures bar plot starts #####
        ##############################################################
        plotBarPlots(outputDir,jobname,numberofSimulations,sample2NumberofIndelsDict,isKeySample,indelsSignatures,
                     sample2SimulationsIndelsSignatureTypesLaggingMediansListDict,
                     sample2SimulationsIndelsSignatureTypesLeadingMediansListDict,
                     sample2Type2ReplicationStrand2CountDict,
                     sample2IndelsSignatureTypesLaggingCountListDict,
                     sample2IndelsSignatureTypesLeadingCountListDict,
                     sample2IndelsSignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict,
                     width,replicationStrands,'indianred', 'goldenrod','All Signatures', 'indels_signatures_replication_strand_bias')
        ############################################################
        #### Plot sample based indels signatures bar plot ends #####
        ############################################################
        ########################################################## New Part ends ####################################################################33

        ########################################################## New Part starts ##################################################################33
        ##############################################################
        #### Plot sample based dinucs signatures bar plot starts #####
        ##############################################################
        plotBarPlots(outputDir,jobname,numberofSimulations,sample2NumberofDinucsDict,isKeySample,dinucsSignatures,
                     sample2SimulationsDinucsSignatureTypesLaggingMediansListDict,
                     sample2SimulationsDinucsSignatureTypesLeadingMediansListDict,
                     sample2Type2ReplicationStrand2CountDict,
                     sample2DinucsSignatureTypesLaggingCountListDict,
                     sample2DinucsSignatureTypesLeadingCountListDict,
                     sample2DinucsSignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict,
                     width,replicationStrands,'indianred', 'goldenrod','All Signatures', 'dinucs_signatures_replication_strand_bias')
        ############################################################
        #### Plot sample based dinucs signatures bar plot ends #####
        ############################################################
        ########################################################## New Part ends ####################################################################33


    ########################################################################
    ######## Bar plot starts includes sample based bar plots ###############
    ##########################  Part 4 ends ################################
    ########################################################################

###################################################################




###################################################################
# This is used up to now.
# main function
def transcriptionReplicationStrandBiasFigures(outputDir,jobname,figureAugmentation,numberofSimulations,sample_based):

    # ####################################################
    # #TODO March 29, 2020
    # #TODO Append a dict for each result to results_list
    # #TODO Create a dataframe df = pd.DataFrame(results_list)
    # transcription_results_dict={}
    # replication_results_dict = {}
    #
    # #TODO Put values of results_dict into results_list
    # # for key in results_dict:
    # #     results_list.append(results_dict[key])
    # transcription_results_list=[]
    # replication_results_list = []
    #
    # #TODO create a dataframe using results_list
    # # df = pd.DataFrame(results_list)
    # # filePath = '/home/burcak/test.txt'
    # # df.to_csv(filePath, sep='\t', header=True, index=False)
    # # transcription_results_df
    # # replication_results_df
    #
    # # TODO fill transcription_results_dict and replication_results_dict in the following function fillSimulationsSample2Type2StrandCountList
    # ####################################################

    jobnamePath = os.path.join(outputDir,jobname,FIGURE,ALL,STRANDBIAS)
    print('Topography.py jobnamePath:%s ' %jobnamePath)

    # ############################################################
    # if (os.path.exists(jobnamePath)):
    #     try:
    #         shutil.rmtree(jobnamePath)
    #     except OSError as e:
    #         print('Error: %s - %s.' % (e.filename, e.strerror))
    # ############################################################

    isFigureAugmentation = False
    if (figureAugmentation == 'augmentation'):
        isFigureAugmentation = True

    #######################################################################################################################
    os.makedirs(os.path.join(outputDir, jobname, FIGURE, ALL, STRANDBIAS,SCATTERPLOTS), exist_ok=True)
    #######################################################################################################################

    ##########################################################################################
    #########################  Read dictionaries related with ################################
    #########################  signatures and samples starts  ################################
    ##########################################################################################
    subsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(os.path.join(outputDir, jobname, DATA,Table_SubsSignature_Cutoff_NumberofMutations_AverageProbability_Filename),sep='\t', header=0,dtype={'cutoff': np.float32,'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})
    indelsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(os.path.join(outputDir, jobname, DATA,Table_IndelsSignature_Cutoff_NumberofMutations_AverageProbability_Filename),sep='\t', header=0,dtype={'cutoff': np.float32,'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})
    dinucsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(os.path.join(outputDir, jobname, DATA,Table_DinucsSignature_Cutoff_NumberofMutations_AverageProbability_Filename),sep='\t', header=0,dtype={'cutoff': np.float32,'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})

    subsSignatures = subsSignature_cutoff_numberofmutations_averageprobability_df['signature'].unique()
    indelsSignatures = indelsSignature_cutoff_numberofmutations_averageprobability_df['signature'].unique()
    dinucsSignatures = dinucsSignature_cutoff_numberofmutations_averageprobability_df['signature'].unique()

    if sample_based:
        sample2NumberofSubsDict = getSample2NumberofSubsDict(outputDir, jobname)
        sample2NumberofIndelsDict = getSample2NumberofIndelsDict(outputDir, jobname)
        sample2NumberofDinucsDict = getDictionary(outputDir, jobname, Sample2NumberofDinucsDictFilename)
        sample2SubsSignature2NumberofMutationsDict = getSample2SubsSignature2NumberofMutationsDict(outputDir,jobname)
        sample2IndelsSignature2NumberofMutationsDict = getSample2IndelsSignature2NumberofMutationsDict(outputDir,jobname)
        sample2DinucsSignature2NumberofMutationsDict = getDictionary(outputDir,jobname, Sample2DinucsSignature2NumberofMutationsDictFilename)
    else:
        sample2NumberofSubsDict = {}
        sample2NumberofIndelsDict = {}
        sample2NumberofDinucsDict = {}
        sample2SubsSignature2NumberofMutationsDict = {}
        sample2IndelsSignature2NumberofMutationsDict = {}
        sample2DinucsSignature2NumberofMutationsDict = {}

    ##########################################################################################
    #########################  Read dictionaries related with ################################
    #########################  signatures and samples ends  ##################################
    ##########################################################################################

    ##########################################################################################
    #####################  Read dictionaries for the rest part starts ########################
    ##########################################################################################

    #########################   Transcription starts   #######################################
    Type2TranscriptionStrand2CountDict_Filepath = os.path.join(outputDir,jobname,DATA,TRANSCRIPTIONSTRANDBIAS,Type2TranscriptionStrand2CountDict_Filename)
    simNum2Type2TranscriptionStrand2CountDict = readDictionary(Type2TranscriptionStrand2CountDict_Filepath)
    type2TranscriptionStrand2CountDict = simNum2Type2TranscriptionStrand2CountDict[str(0)]

    #These are subs signatures
    Signature2MutationType2TranscriptionStrand2CountDict_Filepath = os.path.join(outputDir,jobname,DATA,TRANSCRIPTIONSTRANDBIAS,Signature2MutationType2TranscriptionStrand2CountDict_Filename)
    simNum2SubsSignature2MutationType2TranscriptionStrand2CountDict = readDictionary(Signature2MutationType2TranscriptionStrand2CountDict_Filepath)
    if str(0) in simNum2SubsSignature2MutationType2TranscriptionStrand2CountDict:
        subsSignature2MutationType2TranscriptionStrand2CountDict = simNum2SubsSignature2MutationType2TranscriptionStrand2CountDict[str(0)]
    else:
        subsSignature2MutationType2TranscriptionStrand2CountDict={}

    if sample_based:
        Sample2Type2TranscriptionStrand2CountDict_Filepath = os.path.join(outputDir,jobname,DATA,TRANSCRIPTIONSTRANDBIAS,Sample2Type2TranscriptionStrand2CountDict_Filename)
        simNum2Sample2Type2TranscriptionStrand2CountDict = readDictionary(Sample2Type2TranscriptionStrand2CountDict_Filepath)
        sample2Type2TranscriptionStrand2CountDict = simNum2Sample2Type2TranscriptionStrand2CountDict[str(0)]

        Type2Sample2TranscriptionStrand2CountDict_Filepath = os.path.join(outputDir,jobname,DATA,TRANSCRIPTIONSTRANDBIAS,Type2Sample2TranscriptionStrand2CountDict_Filename)
        simNum2Type2Sample2TranscriptionStrand2CountDict = readDictionary(Type2Sample2TranscriptionStrand2CountDict_Filepath)
        type2Sample2TranscriptionStrand2CountDict = simNum2Type2Sample2TranscriptionStrand2CountDict[str(0)]
    #########################   Transcription ends   #########################################

    #########################   Replication starts   #########################################
    Type2ReplicationStrand2CountDict_Filepath = os.path.join(outputDir,jobname,DATA,REPLICATIONSTRANDBIAS,Type2ReplicationStrand2CountDict_Filename)
    simNum2Type2ReplicationStrand2CountDict = readDictionary(Type2ReplicationStrand2CountDict_Filepath)
    type2ReplicationStrand2CountDict = simNum2Type2ReplicationStrand2CountDict[str(0)]

    #These are subs signatures
    Signature2MutationType2ReplicationStrand2CountDict_Filepath = os.path.join(outputDir,jobname,DATA,REPLICATIONSTRANDBIAS,Signature2MutationType2ReplicationStrand2CountDict_Filename)
    simNum2SubsSignature2MutationType2ReplicationStrand2CountDict = readDictionary(Signature2MutationType2ReplicationStrand2CountDict_Filepath)
    if str(0) in simNum2SubsSignature2MutationType2ReplicationStrand2CountDict:
        subsSignature2MutationType2ReplicationStrand2CountDict = simNum2SubsSignature2MutationType2ReplicationStrand2CountDict[str(0)]
    else:
        subsSignature2MutationType2ReplicationStrand2CountDict={}


    if sample_based:
        Sample2Type2ReplicationStrand2CountDict_Filepath = os.path.join(outputDir,jobname,DATA,REPLICATIONSTRANDBIAS,Sample2Type2ReplicationStrand2CountDict_Filename)
        simNum2Sample2Type2ReplicationStrand2CountDict = readDictionary(Sample2Type2ReplicationStrand2CountDict_Filepath)
        sample2Type2ReplicationStrand2CountDict = simNum2Sample2Type2ReplicationStrand2CountDict[str(0)]

        Type2Sample2ReplicationStrand2CountDict_Filepath = os.path.join(outputDir,jobname,DATA,REPLICATIONSTRANDBIAS,Type2Sample2ReplicationStrand2CountDict_Filename)
        simNum2Type2Sample2ReplicationStrand2CountDict = readDictionary(Type2Sample2ReplicationStrand2CountDict_Filepath)
        type2Sample2ReplicationStrand2CountDict =  simNum2Type2Sample2ReplicationStrand2CountDict[str(0)]
    #########################   Replication ends   ###########################################


    ##########################################################################################
    #####################  Read dictionaries for the rest part ends ##########################
    ##########################################################################################

    ########################################################################
    ##########################  Part 1 starts ##############################
    ########### MutationType All Samples Scatter Figures starts ############
    ######## SignatureType All Samples Scatter Figures starts ##############
    ########################################################################
    if sample_based:
        # Note: Since these are already in sample based no extra work is required.
        # MutationType SampleBased Figures --- No numberofMutations
        if ((type2Sample2TranscriptionStrand2CountDict is not None) and (type2Sample2ReplicationStrand2CountDict is not None)):
            plot_ncomms11383_Supp_FigE_MutationTypeBased_AllSamples_TranscriptionLog10Ratio_ReplicationLog10Ratio(type2Sample2TranscriptionStrand2CountDict,type2Sample2ReplicationStrand2CountDict,outputDir,jobname,isFigureAugmentation)

        #SubsSignature Sample Based Figures --- No numberofMutations
        if ((type2Sample2TranscriptionStrand2CountDict is not None) and (type2Sample2ReplicationStrand2CountDict is not None)):
            plot_ncomms11383_Supp_FigF_SignatureBased_AllSamples_TranscriptionLog10Ratio_ReplicationLog10Ratio(
                type2Sample2TranscriptionStrand2CountDict,
                type2Sample2ReplicationStrand2CountDict,
                subsSignatures,outputDir,jobname,isFigureAugmentation)

        #IndelsSignature Sample Based Figures --- No numberofMutations
        if ((type2Sample2TranscriptionStrand2CountDict is not None) and (type2Sample2ReplicationStrand2CountDict is not None)):
            plot_ncomms11383_Supp_FigF_SignatureBased_AllSamples_TranscriptionLog10Ratio_ReplicationLog10Ratio(
                type2Sample2TranscriptionStrand2CountDict,
                type2Sample2ReplicationStrand2CountDict,
                indelsSignatures,outputDir,jobname,isFigureAugmentation)

        #DinucsSignature Sample Based Figures
        if ((type2Sample2TranscriptionStrand2CountDict is not None) and (type2Sample2ReplicationStrand2CountDict is not None)):
            plot_ncomms11383_Supp_FigF_SignatureBased_AllSamples_TranscriptionLog10Ratio_ReplicationLog10Ratio(
                type2Sample2TranscriptionStrand2CountDict,
                type2Sample2ReplicationStrand2CountDict,
                dinucsSignatures,outputDir,jobname,isFigureAugmentation)
    ########################################################################
    ########### MutationType All Samples Scatter Figures ends ##############
    ######## SignatureType All Samples Scatter Figures ends ################
    ##########################  Part 1 ends ################################
    ########################################################################


    ########################################################################
    ##########################  Part 2 starts ##############################
    ############## Mutation Types Scatter Plots starts #####################
    ############## Signatures Scatter Plots starts #########################
    ########################################################################
    if ((type2TranscriptionStrand2CountDict is not None) and (type2ReplicationStrand2CountDict is not None)):
        plot_ncomms11383_Supp_FigG_AllMutationTypes_TranscriptionLog10Ratio_ReplicationLog10Ratio(None,None,type2TranscriptionStrand2CountDict, type2ReplicationStrand2CountDict,outputDir,jobname)

    if ((type2TranscriptionStrand2CountDict is not None) and (type2ReplicationStrand2CountDict is not None) and (np.any(subsSignatures))):
        plot_ncomms11383_Supp_FigH_AllSignatures_TranscriptionLog10Ratio_ReplicationLog10Ratio('subs',None,None,type2TranscriptionStrand2CountDict,type2ReplicationStrand2CountDict,subsSignature_cutoff_numberofmutations_averageprobability_df,outputDir,jobname)

    if ((type2TranscriptionStrand2CountDict is not None) and (type2ReplicationStrand2CountDict is not None) and (np.any(indelsSignatures))):
        plot_ncomms11383_Supp_FigH_AllSignatures_TranscriptionLog10Ratio_ReplicationLog10Ratio('indels',None,None,type2TranscriptionStrand2CountDict,type2ReplicationStrand2CountDict,indelsSignature_cutoff_numberofmutations_averageprobability_df,outputDir,jobname)

    if ((type2TranscriptionStrand2CountDict is not None) and (type2ReplicationStrand2CountDict is not None) and (np.any(dinucsSignatures))):
        plot_ncomms11383_Supp_FigH_AllSignatures_TranscriptionLog10Ratio_ReplicationLog10Ratio('dinucs',None,None,type2TranscriptionStrand2CountDict,type2ReplicationStrand2CountDict,dinucsSignature_cutoff_numberofmutations_averageprobability_df,outputDir,jobname)
    ########################################################################
    ############## Mutation Types Scatter Plots ends #######################
    ############## Signatures Scatter Plots ends ###########################
    ##########################  Part 2 ends ################################
    ########################################################################


    ########################################################################
    ##########################  Part 3 starts ##############################
    ###################     Sample Based    ################################
    ############## MutationType Scatter Plots starts #######################
    ############## Signature Based Scatter Plots starts ####################
    ########################################################################
    if sample_based:
        ###############################################################
        for sample in sample2NumberofSubsDict:
            if ((sample in sample2Type2TranscriptionStrand2CountDict) and (sample in sample2Type2ReplicationStrand2CountDict)):
                numberofMutations = sample2NumberofSubsDict[sample]
                plot_ncomms11383_Supp_FigG_AllMutationTypes_TranscriptionLog10Ratio_ReplicationLog10Ratio(sample,numberofMutations,sample2Type2TranscriptionStrand2CountDict[sample], sample2Type2ReplicationStrand2CountDict[sample],outputDir,jobname)
        ###############################################################

        ###############################################################
        for sample in sample2NumberofSubsDict:
            if ((sample in sample2Type2TranscriptionStrand2CountDict) and (sample in sample2Type2ReplicationStrand2CountDict)):
                numberofMutations = sample2NumberofSubsDict[sample]
                plot_ncomms11383_Supp_FigH_AllSignatures_TranscriptionLog10Ratio_ReplicationLog10Ratio('subs',sample,numberofMutations,sample2Type2TranscriptionStrand2CountDict[sample], sample2Type2ReplicationStrand2CountDict[sample], subsSignature_cutoff_numberofmutations_averageprobability_df, outputDir,jobname)
        ###############################################################

        ###############################################################
        for sample in sample2NumberofIndelsDict:
            if ((sample in sample2Type2TranscriptionStrand2CountDict) and (sample in sample2Type2ReplicationStrand2CountDict)):
                numberofMutations = sample2NumberofIndelsDict[sample]
                plot_ncomms11383_Supp_FigH_AllSignatures_TranscriptionLog10Ratio_ReplicationLog10Ratio('indels',sample,numberofMutations,sample2Type2TranscriptionStrand2CountDict[sample], sample2Type2ReplicationStrand2CountDict[sample], indelsSignature_cutoff_numberofmutations_averageprobability_df, outputDir,jobname)
        ###############################################################

        ###############################################################
        for sample in sample2NumberofDinucsDict:
            if ((sample in sample2Type2TranscriptionStrand2CountDict) and (sample in sample2Type2ReplicationStrand2CountDict)):
                numberofMutations = sample2NumberofDinucsDict[sample]
                plot_ncomms11383_Supp_FigH_AllSignatures_TranscriptionLog10Ratio_ReplicationLog10Ratio('dinucs',sample,numberofMutations,sample2Type2TranscriptionStrand2CountDict[sample], sample2Type2ReplicationStrand2CountDict[sample], dinucsSignature_cutoff_numberofmutations_averageprobability_df, outputDir,jobname)
        ###############################################################
    ########################################################################
    ###################     Sample Based    ################################
    ############## MutationType Scatter Plots ends #########################
    ############## Signature Based Scatter Plots ends ######################
    ##########################  Part 3 ends ################################
    ########################################################################

    ########################################################################
    ##########################  Part 4 starts ##############################
    ######## Bar plot starts includes sample based bar plots ###############
    ########################################################################
    # Step1: We need to calculate the p-values for each mutation type and signature using Fisher's exact test
    # Step2: We need to calculate FDR adjusted p-values for all mutation types and signatures.
    ########################################################################

    ########################################################################
    ############ MutationType based bar plot starts ########################
    ############ Signature based bar plot starts ###########################
    ########################################################################

    ##########################################################################################################
    ################################### Fill for simulations starts ##########################################
    ##########################################################################################################
    #######
    simulations_mutationtypes_transcribed_medians_list = None
    simulations_mutationtypes_untranscribed_medians_list = None
    simulations_mutationtypes_lagging_medians_list = None
    simulations_mutationtypes_leading_medians_list = None

    sample2SimulationsMutationTypesTranscribedMediansListDict = None
    sample2SimulationsMutationTypesUntranscribedMediansListDict = None
    sample2SimulationsMutationTypesLaggingMediansListDict = None
    sample2SimulationsMutationTypesLeadingMediansListDict = None
    #######

    #######
    simulations_subs_signatures_transcribed_medians_list = None
    simulations_subs_signatures_untranscribed_medians_list = None
    simulations_subs_signatures_lagging_medians_list = None
    simulations_subs_signatures_leading_medians_list = None

    sample2SimulationsSubsSignatureTypesTranscribedMediansListDict = None
    sample2SimulationsSubsSignatureTypesUntranscribedMediansListDict = None
    sample2SimulationsSubsSignatureTypesLaggingMediansListDict = None
    sample2SimulationsSubsSignatureTypesLeadingMediansListDict = None
    #######

    #######
    simulations_indels_signatures_transcribed_medians_list = None
    simulations_indels_signatures_untranscribed_medians_list = None
    simulations_indels_signatures_lagging_medians_list = None
    simulations_indels_signatures_leading_medians_list = None

    sample2SimulationsIndelsSignatureTypesTranscribedMediansListDict = None
    sample2SimulationsIndelsSignatureTypesUntranscribedMediansListDict = None
    sample2SimulationsIndelsSignatureTypesLaggingMediansListDict = None
    sample2SimulationsIndelsSignatureTypesLeadingMediansListDict = None
    #######

    #######
    simulations_dinucs_signatures_transcribed_medians_list = None
    simulations_dinucs_signatures_untranscribed_medians_list = None
    simulations_dinucs_signatures_lagging_medians_list = None
    simulations_dinucs_signatures_leading_medians_list = None

    sample2SimulationsDinucsSignatureTypesTranscribedMediansListDict = None
    sample2SimulationsDinucsSignatureTypesUntranscribedMediansListDict = None
    sample2SimulationsDinucsSignatureTypesLaggingMediansListDict = None
    sample2SimulationsDinucsSignatureTypesLeadingMediansListDict = None
    #######

    #######
    subsSignature2SimulationsMutationTypesTranscribedMediansListDict = None
    subsSignature2SimulationsMutationTypesUntranscribedMediansListDict = None
    subsSignature2SimulationsMutationTypesLaggingMediansListDict = None
    subsSignature2SimulationsMutationTypesLeadingMediansListDict = None
    #######

    if (numberofSimulations > 0):

        #########################################################################################################################
        # simulations --- mutation type --- transcription
        simulations_mutationtypes_transcribed_medians_list, simulations_mutationtypes_untranscribed_medians_list =\
            fillSimulationsType2StrandCountList(simNum2Type2TranscriptionStrand2CountDict,
                                                sixMutationTypes,
                                                numberofSimulations,
                                                TRANSCRIPTIONSTRANDBIAS)


        # simulations --- mutation type --- replication
        simulations_mutationtypes_lagging_medians_list, simulations_mutationtypes_leading_medians_list =\
            fillSimulationsType2StrandCountList(simNum2Type2ReplicationStrand2CountDict,
                                                sixMutationTypes,
                                                numberofSimulations,
                                                REPLICATIONSTRANDBIAS,)

        # simulations --- subs signatures --- transcription
        simulations_subs_signatures_transcribed_medians_list, simulations_subs_signatures_untranscribed_medians_list =\
            fillSimulationsType2StrandCountList(simNum2Type2TranscriptionStrand2CountDict,
                                                subsSignatures,
                                                numberofSimulations,
                                                TRANSCRIPTIONSTRANDBIAS)

        # simulations --- subs signatures --- replication
        simulations_subs_signatures_lagging_medians_list, simulations_subs_signatures_leading_medians_list = \
            fillSimulationsType2StrandCountList(simNum2Type2ReplicationStrand2CountDict,
                                                subsSignatures,
                                                numberofSimulations,
                                                REPLICATIONSTRANDBIAS)

        # simulations --- indels signatures --- transcription
        simulations_indels_signatures_transcribed_medians_list, simulations_indels_signatures_untranscribed_medians_list =\
            fillSimulationsType2StrandCountList(simNum2Type2TranscriptionStrand2CountDict,
                                                indelsSignatures,
                                                numberofSimulations,
                                                TRANSCRIPTIONSTRANDBIAS)

        # simulations --- indels signatures --- replication
        simulations_indels_signatures_lagging_medians_list, simulations_indels_signatures_leading_medians_list = \
            fillSimulationsType2StrandCountList(simNum2Type2ReplicationStrand2CountDict,
                                                indelsSignatures,
                                                numberofSimulations,
                                                REPLICATIONSTRANDBIAS)


        # simulations --- dinucs signatures --- transcription
        simulations_dinucs_signatures_transcribed_medians_list, simulations_dinucs_signatures_untranscribed_medians_list =\
            fillSimulationsType2StrandCountList(simNum2Type2TranscriptionStrand2CountDict,
                                                dinucsSignatures,
                                                numberofSimulations,
                                                TRANSCRIPTIONSTRANDBIAS)

        # simulations --- dinucs signatures --- replication
        simulations_dinucs_signatures_lagging_medians_list, simulations_dinucs_signatures_leading_medians_list = \
            fillSimulationsType2StrandCountList(simNum2Type2ReplicationStrand2CountDict,
                                                dinucsSignatures,
                                                numberofSimulations,
                                                REPLICATIONSTRANDBIAS)

        #########################################################################################################################

        #########################################################################################################################
        # simulations signature --- mutation type --- transcription
        subsSignature2SimulationsMutationTypesTranscribedMediansListDict, subsSignature2SimulationsMutationTypesUntranscribedMediansListDict = \
            fillSimulationsSample2Type2StrandCountList(simNum2SubsSignature2MutationType2TranscriptionStrand2CountDict,
                                                sixMutationTypes,
                                                numberofSimulations,
                                                subsSignature_cutoff_numberofmutations_averageprobability_df['signature'].unique(),
                                                TRANSCRIPTIONSTRANDBIAS)

        # simulations signature --- mutation type --- replication
        subsSignature2SimulationsMutationTypesLaggingMediansListDict, subsSignature2SimulationsMutationTypesLeadingMediansListDict = \
            fillSimulationsSample2Type2StrandCountList(simNum2SubsSignature2MutationType2ReplicationStrand2CountDict,
                                                sixMutationTypes,
                                                numberofSimulations,
                                                subsSignature_cutoff_numberofmutations_averageprobability_df['signature'].unique(),
                                                REPLICATIONSTRANDBIAS)
        #########################################################################################################################

        #########################################################################################################################
        if sample_based:
            #sample based starts
            # samplebased --- simulations --- mutation type --- transcription
            sample2SimulationsMutationTypesTranscribedMediansListDict, sample2SimulationsMutationTypesUntranscribedMediansListDict =\
                fillSimulationsSample2Type2StrandCountList(simNum2Sample2Type2TranscriptionStrand2CountDict,
                    sixMutationTypes,
                    numberofSimulations,
                    sample2NumberofSubsDict.keys(),
                    TRANSCRIPTIONSTRANDBIAS)

            # samplebased --- simulations --- mutation type --- replication
            sample2SimulationsMutationTypesLaggingMediansListDict, sample2SimulationsMutationTypesLeadingMediansListDict =\
                fillSimulationsSample2Type2StrandCountList(simNum2Sample2Type2ReplicationStrand2CountDict,
                    sixMutationTypes,
                    numberofSimulations,
                    sample2NumberofSubsDict.keys(),
                    REPLICATIONSTRANDBIAS)

            # samplebased --- simulations --- signature --- transcription
            sample2SimulationsSubsSignatureTypesTranscribedMediansListDict, sample2SimulationsSubsSignatureTypesUntranscribedMediansListDict =\
                fillSimulationsSample2Type2StrandCountList(simNum2Sample2Type2TranscriptionStrand2CountDict,
                    subsSignatures,
                    numberofSimulations,
                    sample2NumberofSubsDict.keys(),
                    TRANSCRIPTIONSTRANDBIAS)

            # samplebased --- simulations --- signature --- replication
            sample2SimulationsSubsSignatureTypesLaggingMediansListDict, sample2SimulationsSubsSignatureTypesLeadingMediansListDict =\
                fillSimulationsSample2Type2StrandCountList(simNum2Sample2Type2ReplicationStrand2CountDict,
                    subsSignatures,
                    numberofSimulations,
                    sample2NumberofSubsDict.keys(),
                    REPLICATIONSTRANDBIAS)

            # samplebased --- simulations --- signature --- transcription
            sample2SimulationsIndelsSignatureTypesTranscribedMediansListDict, sample2SimulationsIndelsSignatureTypesUntranscribedMediansListDict =\
                fillSimulationsSample2Type2StrandCountList(simNum2Sample2Type2TranscriptionStrand2CountDict,
                    indelsSignatures,
                    numberofSimulations,
                    sample2NumberofIndelsDict.keys(),
                    TRANSCRIPTIONSTRANDBIAS)

            # samplebased --- simulations --- signature --- replication
            sample2SimulationsIndelsSignatureTypesLaggingMediansListDict, sample2SimulationsIndelsSignatureTypesLeadingMediansListDict =\
                fillSimulationsSample2Type2StrandCountList(simNum2Sample2Type2ReplicationStrand2CountDict,
                    indelsSignatures,
                    numberofSimulations,
                    sample2NumberofIndelsDict.keys(),
                    REPLICATIONSTRANDBIAS)
            #sample based ends

            # samplebased --- simulations --- signature --- transcription
            sample2SimulationsDinucsSignatureTypesTranscribedMediansListDict, sample2SimulationsDinucsSignatureTypesUntranscribedMediansListDict =\
                fillSimulationsSample2Type2StrandCountList(simNum2Sample2Type2TranscriptionStrand2CountDict,
                    dinucsSignatures,
                    numberofSimulations,
                    sample2NumberofDinucsDict.keys(),
                    TRANSCRIPTIONSTRANDBIAS)

            # samplebased --- simulations --- signature --- replication
            sample2SimulationsDinucsSignatureTypesLaggingMediansListDict, sample2SimulationsDinucsSignatureTypesLeadingMediansListDict =\
                fillSimulationsSample2Type2StrandCountList(simNum2Sample2Type2ReplicationStrand2CountDict,
                    dinucsSignatures,
                    numberofSimulations,
                    sample2NumberofDinucsDict.keys(),
                    REPLICATIONSTRANDBIAS)
        #########################################################################################################################

    ##########################################################################################################
    ################################### Fill for simulations ends ############################################
    ##########################################################################################################

    #TODO These dictionaries will be written
    type2TranscriptionPValueDict={}
    type2ReplicationPValueDict={}

    #########################################################################################################
    #########################  Calculate p-values for SubsSignature Mutation Type starts  ###################
    #########################################################################################################
    subsSignature2MutationTypesTranscribedCountListDict, \
    subsSignature2MutationTypesUntranscribedCountListDict, \
    subsSignature2MutationTypesTranscriptionPValuesListDict = fillPValuesDictionaries(
        transcriptionStrands,
        sixMutationTypes,
        subsSignature2MutationType2TranscriptionStrand2CountDict,
        subsSignature2SimulationsMutationTypesTranscribedMediansListDict,
        subsSignature2SimulationsMutationTypesUntranscribedMediansListDict)

    subsSignature2MutationTypesLaggingCountListDict, \
    subsSignature2MutationTypesLeadingCountListDict, \
    subsSignature2MutationTypesReplicationPValuesListDict = fillPValuesDictionaries(
        replicationStrands,
        sixMutationTypes,
        subsSignature2MutationType2ReplicationStrand2CountDict,
        subsSignature2SimulationsMutationTypesLaggingMediansListDict,
        subsSignature2SimulationsMutationTypesLeadingMediansListDict)
    #########################################################################################################
    #########################  Calculate p-values for SubsSignature Mutation Type ends  #####################
    #########################################################################################################

    #########################################################################################################
    #############################  Calculate p-values for Mutation Type Starts    ###########################
    #########################################################################################################
    ########################################## Transcription starts  ########################################
    strandType_MutationType2CountList_Dict = fillDictWRTReferenceTypes(sixMutationTypes, transcriptionStrands,type2TranscriptionStrand2CountDict)
    mutationtypes_transcribed_list = strandType_MutationType2CountList_Dict[TRANSCRIBED_STRAND]
    mutationtypes_untranscribed_list = strandType_MutationType2CountList_Dict[UNTRANSCRIBED_STRAND]

    mutationtype_transcription_pvalues = calculate_p_values(mutationtypes_transcribed_list,
                                                            simulations_mutationtypes_transcribed_medians_list,
                                                            mutationtypes_untranscribed_list,
                                                            simulations_mutationtypes_untranscribed_medians_list)

    fillPValuesDict(type2TranscriptionPValueDict, mutationtype_transcription_pvalues,sixMutationTypes)
    ########################################## Transcription ends  ##########################################

    ########################################## Replication starts  ##########################################
    strandType_MutationType2CountList_Dict = fillDictWRTReferenceTypes(sixMutationTypes, replicationStrands,type2ReplicationStrand2CountDict)
    mutationtypes_lagging_list = strandType_MutationType2CountList_Dict[LAGGING]
    mutationtypes_leading_list = strandType_MutationType2CountList_Dict[LEADING]

    # Please notice Fisher exact test does not return anything for large numbers therefore chi2 test is used.
    mutationtype_replication_pvalues = calculate_p_values(mutationtypes_lagging_list,
                                                          simulations_mutationtypes_lagging_medians_list,
                                                          mutationtypes_leading_list,
                                                          simulations_mutationtypes_leading_medians_list)

    fillPValuesDict(type2ReplicationPValueDict, mutationtype_replication_pvalues,sixMutationTypes)
    ########################################## Replication ends  ############################################
    #########################################################################################################
    #############################  Calculate p-values for Mutation Type ends    #############################
    #########################################################################################################

    #########################################################################################################
    ########################  Calculate p-values for subs signatures starts    ##############################
    #########################################################################################################
    ########################################## Transcription starts  ########################################
    strandType_Signature2CountList_Dict = fillDictWRTReferenceTypes(subsSignatures,transcriptionStrands,type2TranscriptionStrand2CountDict)
    subs_signatures_transcribedStrandCount_list = strandType_Signature2CountList_Dict[TRANSCRIBED_STRAND]
    subs_signatures_untranscribedStrandCount_list = strandType_Signature2CountList_Dict[UNTRANSCRIBED_STRAND]

    subs_signature_transcription_pvalues = calculate_p_values(subs_signatures_transcribedStrandCount_list,
                                                         simulations_subs_signatures_transcribed_medians_list,
                                                         subs_signatures_untranscribedStrandCount_list,
                                                         simulations_subs_signatures_untranscribed_medians_list)

    fillPValuesDict(type2TranscriptionPValueDict, subs_signature_transcription_pvalues,subsSignatures)
    ########################################## Transcription ends  ##########################################

    ########################################## Replication starts  ##########################################
    strandType_Signature2CountList_Dict = fillDictWRTReferenceTypes(subsSignatures,replicationStrands,type2ReplicationStrand2CountDict)
    subs_signatures_laggingStrandCount_list = strandType_Signature2CountList_Dict[LAGGING]
    subs_signatures_leadingStrandCount_list = strandType_Signature2CountList_Dict[LEADING]

    subs_signature_replication_pvalues = calculate_p_values(subs_signatures_laggingStrandCount_list,
                                                       simulations_subs_signatures_lagging_medians_list,
                                                       subs_signatures_leadingStrandCount_list,
                                                       simulations_subs_signatures_leading_medians_list)

    fillPValuesDict(type2ReplicationPValueDict, subs_signature_replication_pvalues,subsSignatures)
    ########################################## Replication ends  ############################################
    #########################################################################################################
    ########################  Calculate p-values for subs signatures ends    ################################
    #########################################################################################################


    #########################################################################################################
    ######################  Calculate p-values for indels signatures starts  ################################
    #########################################################################################################
    ########################################## Transcription starts  ########################################
    strandType_Signature2CountList_Dict = fillDictWRTReferenceTypes(indelsSignatures,transcriptionStrands,type2TranscriptionStrand2CountDict)
    indels_signatures_transcribedStrandCount_list = strandType_Signature2CountList_Dict[TRANSCRIBED_STRAND]
    indels_signatures_untranscribedStrandCount_list = strandType_Signature2CountList_Dict[UNTRANSCRIBED_STRAND]


    indels_signature_transcription_pvalues = calculate_p_values(indels_signatures_transcribedStrandCount_list,
                                                         simulations_indels_signatures_transcribed_medians_list,
                                                         indels_signatures_untranscribedStrandCount_list,
                                                         simulations_indels_signatures_untranscribed_medians_list)

    fillPValuesDict(type2TranscriptionPValueDict, indels_signature_transcription_pvalues,indelsSignatures)
    ########################################## Transcription ends  #########################################

    ########################################## Replication starts  ##########################################
    strandType_Signature2CountList_Dict = fillDictWRTReferenceTypes(indelsSignatures,replicationStrands,type2ReplicationStrand2CountDict)
    indels_signatures_laggingStrandCount_list = strandType_Signature2CountList_Dict[LAGGING]
    indels_signatures_leadingStrandCount_list = strandType_Signature2CountList_Dict[LEADING]

    indels_signature_replication_pvalues = calculate_p_values(indels_signatures_laggingStrandCount_list,
                                                       simulations_indels_signatures_lagging_medians_list,
                                                       indels_signatures_leadingStrandCount_list,
                                                       simulations_indels_signatures_leading_medians_list)

    fillPValuesDict(type2ReplicationPValueDict, indels_signature_replication_pvalues,indelsSignatures)
    ########################################## Replication ends  ############################################
    #########################################################################################################
    ######################  Calculate p-values for indels signatures ends    ################################
    #########################################################################################################


    #########################################################################################################
    ######################  Calculate p-values for dinucs signatures starts  ################################
    #########################################################################################################
    ########################################## Transcription starts  ########################################
    strandType_Signature2CountList_Dict = fillDictWRTReferenceTypes(dinucsSignatures,transcriptionStrands,type2TranscriptionStrand2CountDict)
    dinucs_signatures_transcribedStrandCount_list = strandType_Signature2CountList_Dict[TRANSCRIBED_STRAND]
    dinucs_signatures_untranscribedStrandCount_list = strandType_Signature2CountList_Dict[UNTRANSCRIBED_STRAND]

    dinucs_signature_transcription_pvalues = calculate_p_values(dinucs_signatures_transcribedStrandCount_list,
                                                         simulations_dinucs_signatures_transcribed_medians_list,
                                                         dinucs_signatures_untranscribedStrandCount_list,
                                                         simulations_dinucs_signatures_untranscribed_medians_list)

    fillPValuesDict(type2TranscriptionPValueDict, dinucs_signature_transcription_pvalues,dinucsSignatures)
    ########################################## Transcription ends  ##########################################

    ########################################## Replication starts  ##########################################
    strandType_Signature2CountList_Dict = fillDictWRTReferenceTypes(dinucsSignatures,replicationStrands,type2ReplicationStrand2CountDict)
    dinucs_signatures_laggingStrandCount_list = strandType_Signature2CountList_Dict[LAGGING]
    dinucs_signatures_leadingStrandCount_list = strandType_Signature2CountList_Dict[LEADING]

    dinucs_signature_replication_pvalues = calculate_p_values(dinucs_signatures_laggingStrandCount_list,
                                                        simulations_dinucs_signatures_lagging_medians_list,
                                                        dinucs_signatures_leadingStrandCount_list,
                                                        simulations_dinucs_signatures_leading_medians_list)

    fillPValuesDict(type2ReplicationPValueDict, dinucs_signature_replication_pvalues,dinucsSignatures)
    ########################################## Replication ends  ############################################
    #########################################################################################################
    ######################  Calculate p-values for dinucs signatures ends    ################################
    #########################################################################################################

    if sample_based:
        #########################################################################################################
        ####################  Calculate p-values for sample based mutation types starts    ######################
        #########################################    Transcription    ###########################################
        ###########################################    Replication    ###########################################
        #########################################################################################################
        sample2MutationTypesTranscribedCountListDict, sample2MutationTypesUntranscribedCountListDict, sample2MutationTypesTranscriptionPValuesListDict = fillPValuesDictionaries(transcriptionStrands,sixMutationTypes,sample2Type2TranscriptionStrand2CountDict,sample2SimulationsMutationTypesTranscribedMediansListDict,sample2SimulationsMutationTypesUntranscribedMediansListDict)
        sample2MutationTypesLaggingCountListDict, sample2MutationTypesLeadingCountListDict, sample2MutationTypesReplicationPValuesListDict = fillPValuesDictionaries(replicationStrands,sixMutationTypes,sample2Type2ReplicationStrand2CountDict,sample2SimulationsMutationTypesLaggingMediansListDict,sample2SimulationsMutationTypesLeadingMediansListDict)
        #########################################################################################################
        ####################  Calculate p-values for sample based mutation types ends    ########################
        #########################################    Transcription    ###########################################
        ###########################################    Replication    ###########################################
        #########################################################################################################


        #########################################################################################################
        ####################  Calculate p-values for sample based subs signature types starts    ################
        #########################################    Transcription    ###########################################
        ###########################################    Replication    ###########################################
        #########################################################################################################
        sample2SubsSignatureTypesTranscribedCountListDict, sample2SubsSignatureTypesUntranscribedCountListDict, sample2SubsSignature2TranscriptionPValuesListDict = \
            fillPValuesDictionaries(transcriptionStrands,
                                        subsSignatures,
                                        sample2Type2TranscriptionStrand2CountDict,
                                        sample2SimulationsSubsSignatureTypesTranscribedMediansListDict,
                                        sample2SimulationsSubsSignatureTypesUntranscribedMediansListDict)

        sample2SubsSignatureTypesLaggingCountListDict, sample2SubsSignatureTypesLeadingCountListDict, sample2SubsSignature2ReplicationPValuesListDict = \
            fillPValuesDictionaries(replicationStrands,
                                        subsSignatures,
                                        sample2Type2ReplicationStrand2CountDict,
                                        sample2SimulationsSubsSignatureTypesLaggingMediansListDict,
                                        sample2SimulationsSubsSignatureTypesLeadingMediansListDict)
        #########################################################################################################
        ####################  Calculate p-values for sample based subs signature types ends    ##################
        #########################################    Transcription    ###########################################
        ###########################################    Replication    ###########################################
        #########################################################################################################


        #########################################################################################################
        ####################  Calculate p-values for sample based indels signature types starts    ##############
        #########################################    Transcription    ###########################################
        ###########################################    Replication    ###########################################
        #########################################################################################################
        sample2IndelsSignatureTypesTranscribedCountListDict, sample2IndelsSignatureTypesUntranscribedCountListDict, sample2IndelsSignature2TranscriptionPValuesListDict = \
            fillPValuesDictionaries(transcriptionStrands,
                                        indelsSignatures,
                                        sample2Type2TranscriptionStrand2CountDict,
                                        sample2SimulationsIndelsSignatureTypesTranscribedMediansListDict,
                                        sample2SimulationsIndelsSignatureTypesUntranscribedMediansListDict)

        sample2IndelsSignatureTypesLaggingCountListDict, sample2IndelsSignatureTypesLeadingCountListDict, sample2IndelsSignature2ReplicationPValuesListDict = \
            fillPValuesDictionaries(replicationStrands,
                                        indelsSignatures,
                                        sample2Type2ReplicationStrand2CountDict,
                                        sample2SimulationsIndelsSignatureTypesLaggingMediansListDict,
                                        sample2SimulationsIndelsSignatureTypesLeadingMediansListDict)
        #########################################################################################################
        ####################  Calculate p-values for sample based indels signature types ends    ################
        #########################################    Transcription    ###########################################
        ###########################################    Replication    ###########################################
        #########################################################################################################


        #########################################################################################################
        ####################  Calculate p-values for sample based dinucs signature types starts    ##############
        #########################################    Transcription    ###########################################
        ###########################################    Replication    ###########################################
        #########################################################################################################
        sample2DinucsSignatureTypesTranscribedCountListDict, sample2DinucsSignatureTypesUntranscribedCountListDict, sample2DinucsSignature2TranscriptionPValuesListDict = \
            fillPValuesDictionaries(transcriptionStrands,
                                    dinucsSignatures,
                                    sample2Type2TranscriptionStrand2CountDict,
                                    sample2SimulationsDinucsSignatureTypesTranscribedMediansListDict,
                                    sample2SimulationsDinucsSignatureTypesUntranscribedMediansListDict)

        sample2DinucsSignatureTypesLaggingCountListDict, sample2DinucsSignatureTypesLeadingCountListDict, sample2DinucsSignature2ReplicationPValuesListDict = \
            fillPValuesDictionaries(replicationStrands,
                                    dinucsSignatures,
                                    sample2Type2ReplicationStrand2CountDict,
                                    sample2SimulationsDinucsSignatureTypesLaggingMediansListDict,
                                    sample2SimulationsDinucsSignatureTypesLeadingMediansListDict)
        #########################################################################################################
        ####################  Calculate p-values for sample based indels signature types ends    ################
        #########################################    Transcription    ###########################################
        ###########################################    Replication    ###########################################
        #########################################################################################################


    ######################################################
    ################# FDR BH starts ######################
    ######################################################
    #Order the samples
    #We must put the samples in the order they are put in all_p_values list
    all_p_values = []
    all_p_values.extend(mutationtype_transcription_pvalues)
    all_p_values.extend(mutationtype_replication_pvalues)
    all_p_values.extend(subs_signature_transcription_pvalues)
    all_p_values.extend(subs_signature_replication_pvalues)
    all_p_values.extend(indels_signature_transcription_pvalues)
    all_p_values.extend(indels_signature_replication_pvalues)
    all_p_values.extend(dinucs_signature_transcription_pvalues)
    all_p_values.extend(dinucs_signature_replication_pvalues)

    for subsSignature in subsSignatures:
        if subsSignature in subsSignature2MutationTypesTranscriptionPValuesListDict:
            all_p_values.extend(subsSignature2MutationTypesTranscriptionPValuesListDict[subsSignature])

    for subsSignature in subsSignatures:
        if subsSignature in subsSignature2MutationTypesReplicationPValuesListDict:
            all_p_values.extend(subsSignature2MutationTypesReplicationPValuesListDict[subsSignature])

    if sample_based:
        #############################################################################################
        ##################################  Sample Based starts  ####################################
        #############################################################################################
        sampleMutationTypesTranscriptionList = sample2MutationTypesTranscriptionPValuesListDict.keys()
        sampleMutationTypesReplicationList = sample2MutationTypesReplicationPValuesListDict.keys()

        sampleSubsSignaturesTranscriptionList = sample2SubsSignature2TranscriptionPValuesListDict.keys()
        sampleSubsSignaturesReplicationList = sample2SubsSignature2ReplicationPValuesListDict.keys()

        sampleIndelsSignaturesTranscriptionList = sample2IndelsSignature2TranscriptionPValuesListDict.keys()
        sampleIndelsSignaturesReplicationList = sample2IndelsSignature2ReplicationPValuesListDict.keys()

        sampleDinucsSignaturesTranscriptionList = sample2DinucsSignature2TranscriptionPValuesListDict.keys()
        sampleDinucsSignaturesReplicationList = sample2DinucsSignature2ReplicationPValuesListDict.keys()

        for sample in sampleMutationTypesTranscriptionList:
            all_p_values.extend(sample2MutationTypesTranscriptionPValuesListDict[sample])

        for sample in sampleMutationTypesReplicationList:
            all_p_values.extend(sample2MutationTypesReplicationPValuesListDict[sample])

        for sample in sampleSubsSignaturesTranscriptionList:
            all_p_values.extend(sample2SubsSignature2TranscriptionPValuesListDict[sample])

        for sample in sampleSubsSignaturesReplicationList:
            all_p_values.extend(sample2SubsSignature2ReplicationPValuesListDict[sample])

        for sample in sampleIndelsSignaturesTranscriptionList:
            all_p_values.extend(sample2IndelsSignature2TranscriptionPValuesListDict[sample])

        for sample in sampleIndelsSignaturesReplicationList:
            all_p_values.extend(sample2IndelsSignature2ReplicationPValuesListDict[sample])

        for sample in sampleDinucsSignaturesTranscriptionList:
            all_p_values.extend(sample2DinucsSignature2TranscriptionPValuesListDict[sample])

        for sample in sampleDinucsSignaturesReplicationList:
            all_p_values.extend(sample2DinucsSignature2ReplicationPValuesListDict[sample])
        #############################################################################################
        ##################################  Sample Based ends  ######################################
        #############################################################################################

    print('For information: Number of p values to be used in multiple testing: %d' %(len(all_p_values)))
    print('######################################################')
    print('sixMutationTypes')
    print(sixMutationTypes)
    print('mutationtype_transcription_pvalues')
    print(mutationtype_transcription_pvalues)
    print('mutationtype_replication_pvalues')
    print(mutationtype_replication_pvalues)
    print('------------------------------------')

    print('subsSignatures')
    print(subsSignatures)
    print('number of subsSignatures: %d' %len(subsSignatures))
    print('subs_signature_transcription_pvalues')
    print(subs_signature_transcription_pvalues)
    print('subs_signature_replication_pvalues')
    print(subs_signature_replication_pvalues)
    print('------------------------------------')

    print('indelsSignatures')
    print(indelsSignatures)
    print('number of indelsSignatures: %d' %len(indelsSignatures))
    print('indels_signature_transcription_pvalues')
    print(indels_signature_transcription_pvalues)
    print('indels_signature_replication_pvalues')
    print(indels_signature_replication_pvalues)
    print('------------------------------------')

    print('dinucsSignatures')
    print(dinucsSignatures)
    print('number of dinucsSignatures: %d' %len(dinucsSignatures))
    print('dinucs_signature_transcription_pvalues')
    print(dinucs_signature_transcription_pvalues)
    print('dinucs_signature_replication_pvalues')
    print(dinucs_signature_replication_pvalues)
    print('------------------------------------')

    for subsSignature in subsSignatures:
        print(subsSignature)
        print('Transcription')
        print(subsSignature2MutationTypesTranscriptionPValuesListDict[subsSignature])
        print(len(subsSignature2MutationTypesTranscriptionPValuesListDict[subsSignature]))
        print('- - - - - - - - - - - - - - - - - - ')
    print('------------------------------------')

    for subsSignature in subsSignatures:
        print(subsSignature)
        print('Replication')
        print(subsSignature2MutationTypesReplicationPValuesListDict[subsSignature])
        print(len(subsSignature2MutationTypesReplicationPValuesListDict[subsSignature]))
        print('- - - - - - - - - - - - - - - - - - ')
    print('------------------------------------')

    print('all_p_values')
    print(all_p_values)
    print('len(all_p_values)')
    print(len(all_p_values))
    print('######################################################')

    ####################################################################################
    mutationtype_transcription_FDR_BH_adjusted_pvalues = []
    mutationtype_replication_FDR_BH_adjusted_pvalues = []
    subsSignaturetype_transcription_FDR_BH_adjusted_pvalues = []
    subsSignaturetype_replication_FDR_BH_adjusted_pvalues = []
    indelsSignaturetype_transcription_FDR_BH_adjusted_pvalues = []
    indelsSignaturetype_replication_FDR_BH_adjusted_pvalues = []
    dinucsSignaturetype_transcription_FDR_BH_adjusted_pvalues = []
    dinucsSignaturetype_replication_FDR_BH_adjusted_pvalues = []

    #These dictionaries will be written
    type2Transcription_FDR_BH_Adjusted_PValue_Dict={}
    type2Replication_FDR_BH_Adjusted_PValue_Dict={}
    subsSignature2MutationTypesTranscription_FDR_BH_Adjusted_PValues_Dict = {}
    subsSignature2MutationTypesReplication_FDR_BH_Adjusted_PValues_Dict = {}

    if ((all_p_values is  not None) and all_p_values):
        all_p_values_array = np.asarray(all_p_values)

        # rejected, all_FDR_BH_adjusted_p_values = statsmodels.stats.multitest.fdrcorrection(all_p_values_array, alpha=0.05, method='indep', is_sorted=False)
        rejected, all_FDR_BH_adjusted_p_values, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(all_p_values_array, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
        print('For information: Number of all_FDR_BH_adjusted_p_values: %d' % (len(all_FDR_BH_adjusted_p_values)))

        ################################################################
        #Transcription MutationTypes
        start = 0
        end = len(sixMutationTypes)
        mutationtype_transcription_FDR_BH_adjusted_pvalues = (all_FDR_BH_adjusted_p_values[start:end]).tolist()

        startIndex=start
        for mutationType in sixMutationTypes:
            type2Transcription_FDR_BH_Adjusted_PValue_Dict[mutationType]=all_FDR_BH_adjusted_p_values[startIndex]
            startIndex+=1

        #Replication MutationTypes
        start = end
        end += len(sixMutationTypes)
        mutationtype_replication_FDR_BH_adjusted_pvalues = (all_FDR_BH_adjusted_p_values[start:end]).tolist()

        startIndex=start
        for mutationType in sixMutationTypes:
            type2Replication_FDR_BH_Adjusted_PValue_Dict[mutationType]=all_FDR_BH_adjusted_p_values[startIndex]
            startIndex+=1
        ################################################################

        ################################################################
        #Trancription subsSignatures
        start = end
        end += len(subsSignatures)
        subsSignaturetype_transcription_FDR_BH_adjusted_pvalues = all_FDR_BH_adjusted_p_values[start:end]

        startIndex=start
        for signature in subsSignatures:
            type2Transcription_FDR_BH_Adjusted_PValue_Dict[signature] = all_FDR_BH_adjusted_p_values[startIndex]
            startIndex+=1

        #Replication subsSignatures
        start = end
        end += len(subsSignatures)
        subsSignaturetype_replication_FDR_BH_adjusted_pvalues = all_FDR_BH_adjusted_p_values[start:end]

        startIndex=start
        for signature in subsSignatures:
            type2Replication_FDR_BH_Adjusted_PValue_Dict[signature] = all_FDR_BH_adjusted_p_values[startIndex]
            startIndex+=1
        ################################################################

        ################################################################
        #Transcription indelsSignatures
        start = end
        end += len(indelsSignatures)
        indelsSignaturetype_transcription_FDR_BH_adjusted_pvalues = all_FDR_BH_adjusted_p_values[start:end]

        startIndex=start
        for signature in indelsSignatures:
            type2Transcription_FDR_BH_Adjusted_PValue_Dict[signature] = all_FDR_BH_adjusted_p_values[startIndex]
            startIndex+=1

        #Replication indelsSignatures
        start = end
        end += len(indelsSignatures)
        indelsSignaturetype_replication_FDR_BH_adjusted_pvalues = all_FDR_BH_adjusted_p_values[start:end]

        startIndex=start
        for signature in indelsSignatures:
            type2Replication_FDR_BH_Adjusted_PValue_Dict[signature] = all_FDR_BH_adjusted_p_values[startIndex]
            startIndex+=1
        ################################################################

        ################################################################
        #Transcription dinucsSignatures
        start = end
        end += len(dinucsSignatures)
        dinucsSignaturetype_transcription_FDR_BH_adjusted_pvalues = all_FDR_BH_adjusted_p_values[start:end]

        startIndex=start
        for signature in dinucsSignatures:
            type2Transcription_FDR_BH_Adjusted_PValue_Dict[signature] = all_FDR_BH_adjusted_p_values[startIndex]
            startIndex+=1

        #Replication dinucsSignatures
        start = end
        end += len(dinucsSignatures)
        dinucsSignaturetype_replication_FDR_BH_adjusted_pvalues = all_FDR_BH_adjusted_p_values[start:end]

        startIndex=start
        for signature in dinucsSignatures:
            type2Replication_FDR_BH_Adjusted_PValue_Dict[signature] = all_FDR_BH_adjusted_p_values[startIndex]
            startIndex+=1
        ################################################################

        ################################################################
        for subsSignature in subsSignatures:
            if subsSignature in subsSignature2MutationTypesTranscriptionPValuesListDict:
                start = end
                end += len(sixMutationTypes)
                subsSignature2MutationTypesTranscription_FDR_BH_Adjusted_PValues_Dict[subsSignature] = (all_FDR_BH_adjusted_p_values[start:end]).tolist()

        for subsSignature in subsSignatures:
            if subsSignature in subsSignature2MutationTypesReplicationPValuesListDict:
                start = end
                end += len(sixMutationTypes)
                subsSignature2MutationTypesReplication_FDR_BH_Adjusted_PValues_Dict[subsSignature] = (all_FDR_BH_adjusted_p_values[start:end]).tolist()
        ################################################################

        print('######################################################')
        path = os.path.join(outputDir, jobname, FIGURE, ALL, STRANDBIAS)

        #Write Raw P Values dictionaries starts
        filename='type2TranscriptionRawPValueDict.txt'
        writeDictionarySimple(type2TranscriptionPValueDict,path,filename,None)
        filename = 'type2ReplicationRawPValueDict.txt'
        writeDictionarySimple(type2ReplicationPValueDict, path, filename, None)

        filename = 'subsSignature2MutationTypesTranscriptionPValuesListDict.txt'
        writeDictionarySimple(subsSignature2MutationTypesTranscriptionPValuesListDict, path, filename, None)
        filename = 'subsSignature2MutationTypesReplicationPValuesListDict.txt'
        writeDictionarySimple(subsSignature2MutationTypesReplicationPValuesListDict, path, filename, None)
        #Write Raw P Values dictionaries ends

        #type dictionaries
        filename='type2Transcription_FDR_BH_Adjusted_PValue_Dict.txt'
        writeDictionarySimple(type2Transcription_FDR_BH_Adjusted_PValue_Dict,path,filename,None)
        filename = 'type2Replication_FDR_BH_Adjusted_PValue_Dict.txt'
        writeDictionarySimple(type2Replication_FDR_BH_Adjusted_PValue_Dict, path, filename, None)

        #subsSignatures2MutationTypes dictionaries
        filename = 'subsSignature2MutationTypesTranscription_BH_FDR_Adjusted_PValues_Dict.txt'
        writeDictionarySimple(subsSignature2MutationTypesTranscription_FDR_BH_Adjusted_PValues_Dict, path, filename, None)
        filename = 'subsSignature2MutationTypesReplication_BH_FDR_Adjusted_PValues_Dict.txt'
        writeDictionarySimple(subsSignature2MutationTypesReplication_FDR_BH_Adjusted_PValues_Dict, path, filename, None)
        print('######################################################')

        type2SignificantTranscriptionStrandDict=fillSignificantStrandDict(type2TranscriptionStrand2CountDict,type2Transcription_FDR_BH_Adjusted_PValue_Dict)
        type2SignificantReplicationStrandDict=fillSignificantStrandDict(type2ReplicationStrand2CountDict,type2Replication_FDR_BH_Adjusted_PValue_Dict)
        subsSignature2MutationTypesSignificantTranscriptionStrandsListDict=fillMutationTypesSignificantStrandDict(subsSignature2MutationType2TranscriptionStrand2CountDict,subsSignature2MutationTypesTranscription_FDR_BH_Adjusted_PValues_Dict,sixMutationTypes)
        subsSignature2MutationTypesSignificantReplicationStrandsListDict=fillMutationTypesSignificantStrandDict(subsSignature2MutationType2ReplicationStrand2CountDict,subsSignature2MutationTypesReplication_FDR_BH_Adjusted_PValues_Dict,sixMutationTypes)

        print('######################################################')
        path = os.path.join(outputDir, jobname, FIGURE, ALL, STRANDBIAS)

        #type dictionaries
        filename='type2SignificantTranscriptionStrandDict.txt'
        writeDictionarySimple(type2SignificantTranscriptionStrandDict,path,filename,None)
        filename = 'type2SignificantReplicationStrandDict.txt'
        writeDictionarySimple(type2SignificantReplicationStrandDict, path, filename, None)

        #subsSignatures2MutationTypes dictionaries
        filename = 'subsSignature2MutationTypesSignificantTranscriptionStrandsListDict.txt'
        writeDictionarySimple(subsSignature2MutationTypesSignificantTranscriptionStrandsListDict, path, filename, None)
        filename = 'subsSignature2MutationTypesSignificantReplicationStrandsListDict.txt'
        writeDictionarySimple(subsSignature2MutationTypesSignificantReplicationStrandsListDict, path, filename, None)
        print('######################################################')

        print('######################################################')
        print('FDR BH adjusted p values starts')
        print('sixMutationTypes')
        print(sixMutationTypes)
        print('mutationtype_transcription_FDR_BH_adjusted_pvalues')
        print(mutationtype_transcription_FDR_BH_adjusted_pvalues)
        print('mutationtype_replication_FDR_BH_adjusted_pvalues')
        print(mutationtype_replication_FDR_BH_adjusted_pvalues)
        print('------------------------------------')

        print('subsSignatures')
        print(subsSignatures)
        print('number of subsSignatures: %d' % len(subsSignatures))
        print('subsSignaturetype_transcription_FDR_BH_adjusted_pvalues')
        print(subsSignaturetype_transcription_FDR_BH_adjusted_pvalues)
        print('subsSignaturetype_replication_FDR_BH_adjusted_pvalues')
        print(subsSignaturetype_replication_FDR_BH_adjusted_pvalues)
        print('------------------------------------')

        print('indelsSignatures')
        print(indelsSignatures)
        print('number of indelsSignatures: %d' % len(indelsSignatures))
        print('indelsSignaturetype_transcription_FDR_BH_adjusted_pvalues')
        print(indelsSignaturetype_transcription_FDR_BH_adjusted_pvalues)
        print('indelsSignaturetype_replication_FDR_BH_adjusted_pvalues')
        print(indelsSignaturetype_replication_FDR_BH_adjusted_pvalues)
        print('------------------------------------')

        print('dinucsSignatures')
        print(dinucsSignatures)
        print('number of dinucsSignatures: %d' % len(dinucsSignatures))
        print('dinucsSignaturetype_transcription_FDR_BH_adjusted_pvalues')
        print(dinucsSignaturetype_transcription_FDR_BH_adjusted_pvalues)
        print('dinucsSignaturetype_replication_FDR_BH_adjusted_pvalues')
        print(dinucsSignaturetype_replication_FDR_BH_adjusted_pvalues)
        print('------------------------------------')

        for subsSignature in subsSignatures:
            print(subsSignature)
            print('Transcription')
            print(subsSignature2MutationTypesTranscription_FDR_BH_Adjusted_PValues_Dict[subsSignature])
            print(len(subsSignature2MutationTypesTranscription_FDR_BH_Adjusted_PValues_Dict[subsSignature]))
            print('- - - - - - - - - - - - - - - - - - ')
        print('------------------------------------')

        for subsSignature in subsSignatures:
            print(subsSignature)
            print('Replication')
            print(subsSignature2MutationTypesReplication_FDR_BH_Adjusted_PValues_Dict[subsSignature])
            print(len(subsSignature2MutationTypesReplication_FDR_BH_Adjusted_PValues_Dict[subsSignature]))
            print('- - - - - - - - - - - - - - - - - - ')
        print('------------------------------------')

        print('all_FDR_BH_adjusted_p_values')
        print(all_FDR_BH_adjusted_p_values)
        print('len(all_FDR_BH_adjusted_p_values)')
        print(len(all_FDR_BH_adjusted_p_values))
        print('FDR BH adjusted p values ends')
        print('######################################################')


        if sample_based:
            #We must get the BH FDR adjusted p values for each sample as in the order they are put in the all_p_values
            #Fill these dictionaries
            sample2MutationTypesTranscription_BH_FDR_Adjusted_PValues_Dict = {}
            sample2MutationTypesReplication_BH_FDR_Adjusted_PValues_Dict = {}
            sample2SubsSignatureTypesTranscription_BH_FDR_Adjusted_PValues_Dict = {}
            sample2SubsSignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict = {}
            sample2IndelsSignatureTypesTranscription_BH_FDR_Adjusted_PValues_Dict = {}
            sample2IndelsSignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict = {}
            sample2DinucsSignatureTypesTranscription_BH_FDR_Adjusted_PValues_Dict = {}
            sample2DinucsSignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict = {}


            # Sample MutationType Transcription
            for sample in sampleMutationTypesTranscriptionList:
                start = end
                end += len(sixMutationTypes)
                sample2MutationTypesTranscription_BH_FDR_Adjusted_PValues_Dict[sample] = all_FDR_BH_adjusted_p_values[start:end]

            # Sample MutationType Replication
            for sample in sampleMutationTypesReplicationList:
                start = end
                end += len(sixMutationTypes)
                sample2MutationTypesReplication_BH_FDR_Adjusted_PValues_Dict[sample] = all_FDR_BH_adjusted_p_values[start:end]

            # Sample Subs SignatureType Transcription
            for sample in sampleSubsSignaturesTranscriptionList:
                start = end
                end += len(subsSignatures)
                sample2SubsSignatureTypesTranscription_BH_FDR_Adjusted_PValues_Dict[sample] = all_FDR_BH_adjusted_p_values[start:end]

            # Sample Subs SignatureType Replication
            for sample in sampleSubsSignaturesReplicationList:
                start = end
                end += len(subsSignatures)
                sample2SubsSignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict[sample] = all_FDR_BH_adjusted_p_values[start:end]

            # Sample Indels SignatureType Transcription
            for sample in sampleIndelsSignaturesTranscriptionList:
                start = end
                end += len(indelsSignatures)
                sample2IndelsSignatureTypesTranscription_BH_FDR_Adjusted_PValues_Dict[sample] = all_FDR_BH_adjusted_p_values[start:end]

            # Sample Indels SignatureType Replication
            for sample in sampleIndelsSignaturesReplicationList:
                start = end
                end += len(indelsSignatures)
                sample2IndelsSignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict[sample] = all_FDR_BH_adjusted_p_values[start:end]

            # Sample Dinucs SignatureType Transcription
            for sample in sampleDinucsSignaturesTranscriptionList:
                start = end
                end += len(dinucsSignatures)
                sample2DinucsSignatureTypesTranscription_BH_FDR_Adjusted_PValues_Dict[sample] = all_FDR_BH_adjusted_p_values[start:end]

            # Sample Dinucs SignatureType Replication
            for sample in sampleDinucsSignaturesReplicationList:
                start = end
                end += len(dinucsSignatures)
                sample2DinucsSignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict[sample] = all_FDR_BH_adjusted_p_values[start:end]
    ####################################################################################

    ######################################################
    ################# FDR BH ends ########################
    ######################################################

    isKeySample = False
    #######################################################
    ################# Plot mutations starts ###############
    #######################################################
    width = 0.20

    # Transcription
    x_axis_labels = sixMutationTypes
    N = len(x_axis_labels)
    plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,None,isKeySample,None,N, x_axis_labels, mutationtypes_transcribed_list, mutationtypes_untranscribed_list,simulations_mutationtypes_transcribed_medians_list, simulations_mutationtypes_untranscribed_medians_list,mutationtype_transcription_FDR_BH_adjusted_pvalues, 'Transcribed','Untranscribed', 'All Mutations', 'royalblue', 'yellowgreen','mutationtypes_transcription_strand_bias', width)

    # Replication
    x_axis_labels = sixMutationTypes
    N = len(x_axis_labels)
    plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,None,isKeySample,None,N, x_axis_labels, mutationtypes_lagging_list, mutationtypes_leading_list,simulations_mutationtypes_lagging_medians_list, simulations_mutationtypes_leading_medians_list,mutationtype_replication_FDR_BH_adjusted_pvalues, 'Lagging', 'Leading','All Mutations', 'indianred', 'goldenrod', 'mutationtypes_replication_strand_bias',width)
    #######################################################
    ################# Plot mutations ends #################
    #######################################################

    #######################################################
    ########### Plot subs signatures starts ###############
    #######################################################
    # Transcription
    x_axis_transcription_labels = subsSignatures
    N_signatures_transcription = len(x_axis_transcription_labels)

    if (np.any(x_axis_transcription_labels) and subs_signatures_transcribedStrandCount_list and subs_signatures_untranscribedStrandCount_list and simulations_subs_signatures_transcribed_medians_list and simulations_subs_signatures_untranscribed_medians_list and (subsSignaturetype_transcription_FDR_BH_adjusted_pvalues.size>0)):
        zipped = zip(x_axis_transcription_labels, subs_signatures_transcribedStrandCount_list, subs_signatures_untranscribedStrandCount_list, simulations_subs_signatures_transcribed_medians_list, simulations_subs_signatures_untranscribed_medians_list, subsSignaturetype_transcription_FDR_BH_adjusted_pvalues)
        zipped = sorted(zipped, key=lambda x: natural_key(x[0]))
        x_axis_transcription_labels, subs_signatures_transcribedStrandCount_list, subs_signatures_untranscribedStrandCount_list, simulations_subs_signatures_transcribed_medians_list, simulations_subs_signatures_untranscribed_medians_list, subsSignaturetype_transcription_FDR_BH_adjusted_pvalues = zip(*zipped)

        if (subs_signatures_transcribedStrandCount_list and subs_signatures_untranscribedStrandCount_list):
            plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,None,isKeySample,None,N_signatures_transcription,x_axis_transcription_labels,subs_signatures_transcribedStrandCount_list,subs_signatures_untranscribedStrandCount_list,simulations_subs_signatures_transcribed_medians_list, simulations_subs_signatures_untranscribed_medians_list,subsSignaturetype_transcription_FDR_BH_adjusted_pvalues, 'Transcribed','Untranscribed', 'All Signatures', 'royalblue', 'yellowgreen','subs_signatures_transcription_strand_bias', width)

    # Replication
    x_axis_replication_labels = subsSignatures
    N_signatures_replication = len(x_axis_replication_labels)

    if (np.any(x_axis_replication_labels) and subs_signatures_laggingStrandCount_list and subs_signatures_leadingStrandCount_list and simulations_subs_signatures_lagging_medians_list and simulations_subs_signatures_leading_medians_list and (subsSignaturetype_replication_FDR_BH_adjusted_pvalues.size>0)):
        zipped = zip(x_axis_replication_labels, subs_signatures_laggingStrandCount_list, subs_signatures_leadingStrandCount_list, simulations_subs_signatures_lagging_medians_list, simulations_subs_signatures_leading_medians_list, subsSignaturetype_replication_FDR_BH_adjusted_pvalues)
        zipped = sorted(zipped, key=lambda x: natural_key(x[0]))
        x_axis_replication_labels, subs_signatures_laggingStrandCount_list, subs_signatures_leadingStrandCount_list, simulations_subs_signatures_lagging_medians_list, simulations_subs_signatures_leading_medians_list, subsSignaturetype_replication_FDR_BH_adjusted_pvalues = zip(*zipped)

        if (subs_signatures_laggingStrandCount_list and subs_signatures_leadingStrandCount_list):
            plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,None,isKeySample,None,N_signatures_replication,x_axis_replication_labels,subs_signatures_laggingStrandCount_list,subs_signatures_leadingStrandCount_list, simulations_subs_signatures_lagging_medians_list, simulations_subs_signatures_leading_medians_list,subsSignaturetype_replication_FDR_BH_adjusted_pvalues, 'Lagging','Leading', 'All Signatures', 'indianred', 'goldenrod','subs_signatures_replication_strand_bias', width)
    #######################################################
    ########### Plot subs signatures ends #################
    #######################################################

    #######################################################
    ########### Plot indels signatures starts #############
    #######################################################
    # Transcription
    x_axis_transcription_labels = indelsSignatures
    N_signatures_transcription = len(x_axis_transcription_labels)

    if (np.any(x_axis_transcription_labels) and indels_signatures_transcribedStrandCount_list and indels_signatures_untranscribedStrandCount_list and simulations_indels_signatures_transcribed_medians_list and simulations_indels_signatures_untranscribed_medians_list and (indelsSignaturetype_transcription_FDR_BH_adjusted_pvalues.size>0)):
        zipped = zip(x_axis_transcription_labels, indels_signatures_transcribedStrandCount_list, indels_signatures_untranscribedStrandCount_list, simulations_indels_signatures_transcribed_medians_list, simulations_indels_signatures_untranscribed_medians_list, indelsSignaturetype_transcription_FDR_BH_adjusted_pvalues)
        zipped = sorted(zipped, key=lambda x: natural_key(x[0]))
        x_axis_transcription_labels, indels_signatures_transcribedStrandCount_list, indels_signatures_untranscribedStrandCount_list, simulations_indels_signatures_transcribed_medians_list, simulations_indels_signatures_untranscribed_medians_list, indelsSignaturetype_transcription_FDR_BH_adjusted_pvalues = zip(*zipped)

        if (indels_signatures_transcribedStrandCount_list and indels_signatures_untranscribedStrandCount_list):
            plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,None,isKeySample,None,N_signatures_transcription,x_axis_transcription_labels, indels_signatures_transcribedStrandCount_list,indels_signatures_untranscribedStrandCount_list,simulations_indels_signatures_transcribed_medians_list, simulations_indels_signatures_untranscribed_medians_list,indelsSignaturetype_transcription_FDR_BH_adjusted_pvalues, 'Transcribed','Untranscribed', 'All Signatures', 'royalblue', 'yellowgreen','indels_signatures_transcription_strand_bias', width)

    # Replication
    x_axis_replication_labels = indelsSignatures
    N_signatures_replication = len(x_axis_replication_labels)

    if (np.any(x_axis_replication_labels) and indels_signatures_laggingStrandCount_list and indels_signatures_leadingStrandCount_list and simulations_indels_signatures_lagging_medians_list and simulations_indels_signatures_leading_medians_list and (indelsSignaturetype_replication_FDR_BH_adjusted_pvalues.size>0)):
        zipped = zip(x_axis_replication_labels, indels_signatures_laggingStrandCount_list, indels_signatures_leadingStrandCount_list, simulations_indels_signatures_lagging_medians_list, simulations_indels_signatures_leading_medians_list, indelsSignaturetype_replication_FDR_BH_adjusted_pvalues)
        zipped = sorted(zipped, key=lambda x: natural_key(x[0]))
        x_axis_replication_labels, indels_signatures_laggingStrandCount_list, indels_signatures_leadingStrandCount_list, simulations_indels_signatures_lagging_medians_list, simulations_indels_signatures_leading_medians_list, indelsSignaturetype_replication_FDR_BH_adjusted_pvalues = zip(*zipped)

        if (indels_signatures_laggingStrandCount_list and indels_signatures_leadingStrandCount_list):
            plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,None,isKeySample,None,N_signatures_replication,x_axis_replication_labels,indels_signatures_laggingStrandCount_list,indels_signatures_leadingStrandCount_list, simulations_indels_signatures_lagging_medians_list, simulations_indels_signatures_leading_medians_list,indelsSignaturetype_replication_FDR_BH_adjusted_pvalues, 'Lagging','Leading', 'All Signatures', 'indianred', 'goldenrod','indels_signatures_replication_strand_bias', width)
    #######################################################
    ########### Plot indels signatures ends ###############
    #######################################################

    #######################################################
    ########### Plot dinucs signatures starts #############
    #######################################################
    # Transcription
    x_axis_transcription_labels = dinucsSignatures
    N_signatures_transcription = len(x_axis_transcription_labels)

    if (np.any(x_axis_transcription_labels) and dinucs_signatures_transcribedStrandCount_list and dinucs_signatures_untranscribedStrandCount_list and simulations_dinucs_signatures_transcribed_medians_list and simulations_dinucs_signatures_untranscribed_medians_list and (dinucsSignaturetype_transcription_FDR_BH_adjusted_pvalues.size>0)):
        zipped = zip(x_axis_transcription_labels, dinucs_signatures_transcribedStrandCount_list, dinucs_signatures_untranscribedStrandCount_list, simulations_dinucs_signatures_transcribed_medians_list, simulations_dinucs_signatures_untranscribed_medians_list, dinucsSignaturetype_transcription_FDR_BH_adjusted_pvalues)
        zipped = sorted(zipped, key=lambda x: natural_key(x[0]))
        x_axis_transcription_labels, dinucs_signatures_transcribedStrandCount_list, dinucs_signatures_untranscribedStrandCount_list, simulations_dinucs_signatures_transcribed_medians_list, simulations_dinucs_signatures_untranscribed_medians_list, dinucsSignaturetype_transcription_FDR_BH_adjusted_pvalues = zip(*zipped)

        if (dinucs_signatures_transcribedStrandCount_list and dinucs_signatures_untranscribedStrandCount_list):
            plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,None,isKeySample,None,N_signatures_transcription,x_axis_transcription_labels, dinucs_signatures_transcribedStrandCount_list,dinucs_signatures_untranscribedStrandCount_list,simulations_dinucs_signatures_transcribed_medians_list, simulations_dinucs_signatures_untranscribed_medians_list,dinucsSignaturetype_transcription_FDR_BH_adjusted_pvalues, 'Transcribed','Untranscribed', 'All Signatures', 'royalblue', 'yellowgreen','dinucs_signatures_transcription_strand_bias', width)

    # Replication
    x_axis_replication_labels = dinucsSignatures
    N_signatures_replication = len(x_axis_replication_labels)

    if (np.any(x_axis_replication_labels) and dinucs_signatures_laggingStrandCount_list and dinucs_signatures_leadingStrandCount_list and simulations_dinucs_signatures_lagging_medians_list and simulations_dinucs_signatures_leading_medians_list and (dinucsSignaturetype_replication_FDR_BH_adjusted_pvalues.size>0)):
        zipped = zip(x_axis_replication_labels, dinucs_signatures_laggingStrandCount_list, dinucs_signatures_leadingStrandCount_list, simulations_dinucs_signatures_lagging_medians_list, simulations_dinucs_signatures_leading_medians_list, dinucsSignaturetype_replication_FDR_BH_adjusted_pvalues)
        zipped = sorted(zipped, key=lambda x: natural_key(x[0]))
        x_axis_replication_labels, dinucs_signatures_laggingStrandCount_list, dinucs_signatures_leadingStrandCount_list, simulations_dinucs_signatures_lagging_medians_list, simulations_dinucs_signatures_leading_medians_list, dinucsSignaturetype_replication_FDR_BH_adjusted_pvalues = zip(*zipped)

        if (dinucs_signatures_laggingStrandCount_list and dinucs_signatures_leadingStrandCount_list):
            plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,None,isKeySample,None,N_signatures_replication,x_axis_replication_labels,dinucs_signatures_laggingStrandCount_list,dinucs_signatures_leadingStrandCount_list, simulations_dinucs_signatures_lagging_medians_list, simulations_dinucs_signatures_leading_medians_list,dinucsSignaturetype_replication_FDR_BH_adjusted_pvalues, 'Lagging','Leading', 'All Signatures', 'indianred', 'goldenrod','dinucs_signatures_replication_strand_bias', width)
    #######################################################
    ########### Plot dinucs signatures ends ###############
    #######################################################


    #################################################################
    ########### Plot sub signatures mutation types starts ###########
    #################################################################
    plotBarPlots(outputDir,jobname,numberofSimulations,subsSignature_cutoff_numberofmutations_averageprobability_df,isKeySample,sixMutationTypes,
                 subsSignature2SimulationsMutationTypesTranscribedMediansListDict,
                 subsSignature2SimulationsMutationTypesUntranscribedMediansListDict,
                 subsSignature2MutationType2TranscriptionStrand2CountDict,
                 subsSignature2MutationTypesTranscribedCountListDict,
                 subsSignature2MutationTypesUntranscribedCountListDict,
                 subsSignature2MutationTypesTranscription_FDR_BH_Adjusted_PValues_Dict,
                 width,transcriptionStrands,'royalblue', 'yellowgreen', 'All Mutations','mutationtypes_transcription_strand_bias')

    plotBarPlots(outputDir,jobname,numberofSimulations,subsSignature_cutoff_numberofmutations_averageprobability_df,isKeySample,sixMutationTypes,
                 subsSignature2SimulationsMutationTypesLaggingMediansListDict,
                 subsSignature2SimulationsMutationTypesLeadingMediansListDict,
                 subsSignature2MutationType2ReplicationStrand2CountDict,
                 subsSignature2MutationTypesLaggingCountListDict,
                 subsSignature2MutationTypesLeadingCountListDict,
                 subsSignature2MutationTypesReplication_FDR_BH_Adjusted_PValues_Dict,
                 width,replicationStrands,'indianred', 'goldenrod','All Mutations','mutationtypes_replication_strand_bias')
    #################################################################
    ########### Plot sub signatures mutation types ends #############
    #################################################################


    isKeySample=True
    if sample_based:
        #######################################################
        ######### Plot sample based mutations starts ##########
        #######################################################
        plotBarPlots(outputDir,jobname,numberofSimulations,sample2NumberofSubsDict,isKeySample,sixMutationTypes,
                     sample2SimulationsMutationTypesTranscribedMediansListDict,
                     sample2SimulationsMutationTypesUntranscribedMediansListDict,
                     sample2Type2TranscriptionStrand2CountDict,
                     sample2MutationTypesTranscribedCountListDict,
                     sample2MutationTypesUntranscribedCountListDict,
                     sample2MutationTypesTranscription_BH_FDR_Adjusted_PValues_Dict,
                     width,transcriptionStrands,'royalblue', 'yellowgreen', 'All Mutations','mutationtypes_transcription_strand_bias')

        plotBarPlots(outputDir,jobname,numberofSimulations,sample2NumberofSubsDict,isKeySample,sixMutationTypes,
                     sample2SimulationsMutationTypesLaggingMediansListDict,
                     sample2SimulationsMutationTypesLeadingMediansListDict,
                     sample2Type2ReplicationStrand2CountDict,
                     sample2MutationTypesLaggingCountListDict,
                     sample2MutationTypesLeadingCountListDict,
                     sample2MutationTypesReplication_BH_FDR_Adjusted_PValues_Dict,
                     width,replicationStrands,'indianred', 'goldenrod','All Mutations','mutationtypes_replication_strand_bias')
        #######################################################
        ######### Plot sample based mutations ends ############
        #######################################################


        ############################################################
        #### Plot sample based subs signatures bar plot starts #####
        ############################################################
        plotBarPlots(outputDir,jobname,numberofSimulations,sample2NumberofSubsDict,isKeySample,subsSignatures,
                     sample2SimulationsSubsSignatureTypesTranscribedMediansListDict,
                     sample2SimulationsSubsSignatureTypesUntranscribedMediansListDict,
                     sample2Type2TranscriptionStrand2CountDict,
                     sample2SubsSignatureTypesTranscribedCountListDict,
                     sample2SubsSignatureTypesUntranscribedCountListDict,
                     sample2SubsSignatureTypesTranscription_BH_FDR_Adjusted_PValues_Dict,
                     width,transcriptionStrands,'royalblue', 'yellowgreen','All Signatures','subs_signatures_transcription_strand_bias')

        plotBarPlots(outputDir,jobname,numberofSimulations,sample2NumberofSubsDict,isKeySample,subsSignatures,
                     sample2SimulationsSubsSignatureTypesLaggingMediansListDict,
                     sample2SimulationsSubsSignatureTypesLeadingMediansListDict,
                     sample2Type2ReplicationStrand2CountDict,
                     sample2SubsSignatureTypesLaggingCountListDict,
                     sample2SubsSignatureTypesLeadingCountListDict,
                     sample2SubsSignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict,
                     width,replicationStrands,'indianred', 'goldenrod','All Signatures','subs_signatures_replication_strand_bias')
        ############################################################
        #### Plot sample based subs signatures bar plot ends #######
        ############################################################

        ########################################################## New Part starts ##################################################################33
        ##############################################################
        #### Plot sample based indels signatures bar plot starts #####
        ##############################################################
        plotBarPlots(outputDir,jobname,numberofSimulations,sample2NumberofIndelsDict,isKeySample,indelsSignatures,
                     sample2SimulationsIndelsSignatureTypesTranscribedMediansListDict,
                     sample2SimulationsIndelsSignatureTypesUntranscribedMediansListDict,
                     sample2Type2TranscriptionStrand2CountDict,
                     sample2IndelsSignatureTypesTranscribedCountListDict,
                     sample2IndelsSignatureTypesUntranscribedCountListDict,
                     sample2IndelsSignatureTypesTranscription_BH_FDR_Adjusted_PValues_Dict,
                     width,transcriptionStrands,'royalblue', 'yellowgreen','All Signatures', 'indels_signatures_transcription_strand_bias')

        plotBarPlots(outputDir,jobname,numberofSimulations,sample2NumberofIndelsDict,isKeySample,indelsSignatures,
                     sample2SimulationsIndelsSignatureTypesLaggingMediansListDict,
                     sample2SimulationsIndelsSignatureTypesLeadingMediansListDict,
                     sample2Type2ReplicationStrand2CountDict,
                     sample2IndelsSignatureTypesLaggingCountListDict,
                     sample2IndelsSignatureTypesLeadingCountListDict,
                     sample2IndelsSignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict,
                     width,replicationStrands,'indianred', 'goldenrod','All Signatures', 'indels_signatures_replication_strand_bias')
        ############################################################
        #### Plot sample based indels signatures bar plot ends #####
        ############################################################
        ########################################################## New Part ends ####################################################################33

        ########################################################## New Part starts ##################################################################33
        ##############################################################
        #### Plot sample based dinucs signatures bar plot starts #####
        ##############################################################
        plotBarPlots(outputDir,jobname,numberofSimulations,sample2NumberofDinucsDict,isKeySample,dinucsSignatures,
                     sample2SimulationsDinucsSignatureTypesTranscribedMediansListDict,
                     sample2SimulationsDinucsSignatureTypesUntranscribedMediansListDict,
                     sample2Type2TranscriptionStrand2CountDict,
                     sample2DinucsSignatureTypesTranscribedCountListDict,
                     sample2DinucsSignatureTypesUntranscribedCountListDict,
                     sample2DinucsSignatureTypesTranscription_BH_FDR_Adjusted_PValues_Dict,
                     width,transcriptionStrands,'royalblue', 'yellowgreen','All Signatures', 'dinucs_signatures_transcription_strand_bias')

        plotBarPlots(outputDir,jobname,numberofSimulations,sample2NumberofDinucsDict,isKeySample,dinucsSignatures,
                     sample2SimulationsDinucsSignatureTypesLaggingMediansListDict,
                     sample2SimulationsDinucsSignatureTypesLeadingMediansListDict,
                     sample2Type2ReplicationStrand2CountDict,
                     sample2DinucsSignatureTypesLaggingCountListDict,
                     sample2DinucsSignatureTypesLeadingCountListDict,
                     sample2DinucsSignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict,
                     width,replicationStrands,'indianred', 'goldenrod','All Signatures', 'dinucs_signatures_replication_strand_bias')
        ############################################################
        #### Plot sample based dinucs signatures bar plot ends #####
        ############################################################
        ########################################################## New Part ends ####################################################################33


    ########################################################################
    ######## Bar plot starts includes sample based bar plots ###############
    ##########################  Part 4 ends ################################
    ########################################################################