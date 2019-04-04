# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018 Burcak Otlu

import sys
import os
import scipy.stats as stats
import statsmodels.stats.multitest

import matplotlib
BACKEND = 'Agg'
if matplotlib.get_backend().lower() != BACKEND.lower():
    # If backend is not set properly a call to describe will hang
    matplotlib.use(BACKEND)

from matplotlib import pyplot as plt
import matplotlib as mpl

from SigProfilerTopography.source.commons.TopographyCommons import *

transcriptionStrands = [TRANSCRIBED_STRAND, UNTRANSCRIBED_STRAND]
replicationStrands = [LAGGING, LEADING]

########################################################################
#For Mutation Types
def plot_ncomms11383_Supp_FigG_AllMutationTypes_TranscriptionLog10Ratio_ReplicationLog10Ratio(sample,numberofMutations,mutationType2TranscriptionStrand2CountDict, mutationType2ReplicationStrand2CountDict,outputDir,jobname,isFigureAugmentation):

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

    if (isFigureAugmentation):
        plt.title(jobname, fontsize=20, fontweight='bold')

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
        if (mutationType in mutationType2TranscriptionStrand2CountDict) and (mutationType in mutationType2ReplicationStrand2CountDict):

            if ((TRANSCRIBED_STRAND in mutationType2TranscriptionStrand2CountDict[mutationType]) and (UNTRANSCRIBED_STRAND in mutationType2TranscriptionStrand2CountDict[mutationType])):
                transcriptionRatiosDict[mutationType]= np.log10(mutationType2TranscriptionStrand2CountDict[mutationType][TRANSCRIBED_STRAND]/mutationType2TranscriptionStrand2CountDict[mutationType][UNTRANSCRIBED_STRAND])

            if ((LAGGING in mutationType2ReplicationStrand2CountDict[mutationType]) and (LEADING in mutationType2ReplicationStrand2CountDict[mutationType])):
                replicationRatiosDict[mutationType] = np.log10(mutationType2ReplicationStrand2CountDict[mutationType][LAGGING]/mutationType2ReplicationStrand2CountDict[mutationType][LEADING])

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
        figureFile = os.path.join(outputDir, jobname, FIGURE, ALL, STRANDBIAS, figureName)

    else:
        figureName = 'all_mutationtypes_%s_%d_%s_scatterplots.png' %(sample,numberofMutations,STRANDBIAS)
        os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, STRANDBIAS), exist_ok=True)
        figureFile = os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, STRANDBIAS, figureName)

    fig.savefig(figureFile)
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
        signature2NumberofMutationsDict,
        outputDir,
        jobname,
        isFigureAugmentation):

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

    if (isFigureAugmentation):
        plt.title(jobname, fontsize=20, fontweight='bold')

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

    for signature in signature2NumberofMutationsDict:
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

        for signature in signature2NumberofMutationsDict:
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
            figureFile = os.path.join(outputDir, jobname, FIGURE, ALL, STRANDBIAS, figureName)
        else:
            figureName = 'all_%s_signatures_%s_%d_%s_scatterplots.png' % (
            signatureType, sample, numberofMutations, STRANDBIAS)
            os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, STRANDBIAS), exist_ok=True)
            figureFile = os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, STRANDBIAS, figureName)

        fig.savefig(figureFile)
        plt.close(fig)
########################################################################

########################################################################
#MutationTypeBased SampleBased Figures
def plot_ncomms11383_Supp_FigE_MutationTypeBased_AllSamples_TranscriptionLog10Ratio_ReplicationLog10Ratio(
        mutationType2Sample2TranscriptionStrand2CountDict,
        mutationType2Sample2ReplicationStrand2CountDict,
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
        if mutationType in mutationType2Sample2TranscriptionStrand2CountDict:
            for sample in mutationType2Sample2TranscriptionStrand2CountDict[mutationType].keys():
                if ((TRANSCRIBED_STRAND in mutationType2Sample2TranscriptionStrand2CountDict[mutationType][sample].keys()) and (UNTRANSCRIBED_STRAND in mutationType2Sample2TranscriptionStrand2CountDict[mutationType][sample].keys())):
                    transcriptionRatiosDict[mutationType][sample]= np.log10(mutationType2Sample2TranscriptionStrand2CountDict[mutationType][sample][TRANSCRIBED_STRAND]/mutationType2Sample2TranscriptionStrand2CountDict[mutationType][sample][UNTRANSCRIBED_STRAND])

        if mutationType in mutationType2Sample2ReplicationStrand2CountDict:
            for sample in mutationType2Sample2ReplicationStrand2CountDict[mutationType].keys():
                if ((LAGGING in mutationType2Sample2ReplicationStrand2CountDict[mutationType][sample].keys()) and (LEADING in mutationType2Sample2ReplicationStrand2CountDict[mutationType][sample].keys())):
                    replicationRatiosDict[mutationType][sample] = np.log10(mutationType2Sample2ReplicationStrand2CountDict[mutationType][sample][LAGGING]/mutationType2Sample2ReplicationStrand2CountDict[mutationType][sample][LEADING])

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

        if (mutationType in mutationType2Sample2TranscriptionStrand2CountDict):
            for sample in mutationType2Sample2TranscriptionStrand2CountDict[mutationType].keys():
                if ((sample in replicationRatiosDict[mutationType].keys()) and (sample in transcriptionRatiosDict[mutationType].keys())):
                    plt.scatter(replicationRatiosDict[mutationType][sample],transcriptionRatiosDict[mutationType][sample], facecolor='none', color=mutationType2ColorDict[mutationType])

        plt.axvline(x=0.0, color='gray', linestyle='--')
        plt.axhline(y=0.0, color='gray', linestyle='--')

        if (isFigureAugmentation):
            plt.title(jobname + ' ' + mutationType)

        newMutationType = mutationType.replace('>', '2')

        figureName = newMutationType + '_MutationType_' + STRANDBIAS + '.png'
        figureFile = os.path.join(outputDir,jobname,FIGURE,ALL,STRANDBIAS,figureName)
        fig.savefig(figureFile)
        plt.close(fig)
########################################################################


########################################################################
#SignatureBased SampleBased Figures
#Sig26 is very different
def plot_ncomms11383_Supp_FigF_SignatureBased_AllSamples_TranscriptionLog10Ratio_ReplicationLog10Ratio(mutationProbability2Signature2Sample2TranscriptionStrand2CountDict,mutationProbability2Signature2Sample2ReplicationStrand2CountDict,signatures,outputDir,jobname,isFigureAugmentation,mutationProbabilityThreshold):

    signature2Sample2TranscriptionStrand2CountDict= mutationProbability2Signature2Sample2TranscriptionStrand2CountDict[mutationProbabilityThreshold]
    signature2Sample2ReplicationStrand2CountDict = mutationProbability2Signature2Sample2ReplicationStrand2CountDict[mutationProbabilityThreshold]

    transcriptionRatiosDict = {}
    replicationRatiosDict = {}
    for signature in signatures:
        # initialization
        if signature not in transcriptionRatiosDict:
            transcriptionRatiosDict[signature] = {}
        if signature not in replicationRatiosDict:
            replicationRatiosDict[signature] = {}
        # Fill the dictionaries
        if signature in signature2Sample2TranscriptionStrand2CountDict:
            for sample in signature2Sample2TranscriptionStrand2CountDict[signature].keys():
                if (UNTRANSCRIBED_STRAND in signature2Sample2TranscriptionStrand2CountDict[signature][sample]) and (TRANSCRIBED_STRAND in signature2Sample2TranscriptionStrand2CountDict[signature][sample]):
                    transcriptionRatiosDict[signature][sample] = np.log10(signature2Sample2TranscriptionStrand2CountDict[signature][sample][TRANSCRIBED_STRAND] /signature2Sample2TranscriptionStrand2CountDict[signature][sample][UNTRANSCRIBED_STRAND])
                    # print(signature, sample)
                    # print(signature2Sample2TranscriptionStrand2CountDict[signature][sample][TRANSCRIBED_STRAND])
                    # print(signature2Sample2TranscriptionStrand2CountDict[signature][sample][UNTRANSCRIBED_STRAND])
                    # print(signature,sample,transcriptionRatiosDict[signature][sample])

        if signature in signature2Sample2ReplicationStrand2CountDict:
            for sample in signature2Sample2ReplicationStrand2CountDict[signature].keys():
                if (LAGGING in signature2Sample2ReplicationStrand2CountDict[signature][sample]) and (LEADING in signature2Sample2ReplicationStrand2CountDict[signature][sample]):
                    replicationRatiosDict[signature][sample] = np.log10(signature2Sample2ReplicationStrand2CountDict[signature][sample][LAGGING] /signature2Sample2ReplicationStrand2CountDict[signature][sample][LEADING])

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

            for sample in signature2Sample2TranscriptionStrand2CountDict[signature].keys():
                if (sample in replicationRatiosDict[signature]) and (sample in transcriptionRatiosDict[signature]):
                    plt.scatter(replicationRatiosDict[signature][sample], transcriptionRatiosDict[signature][sample],facecolor='none',color='green')

            plt.axvline(x=0.0, color='gray', linestyle='--')
            plt.axhline(y=0.0, color='gray', linestyle='--')

            if (isFigureAugmentation):
                plt.title(jobname + ' ' + signature)

            figureName = signature.replace(' ','') + '_Signature_' + STRANDBIAS + '.png'
            figureFile = os.path.join(outputDir,jobname,FIGURE,ALL,STRANDBIAS,figureName)
            fig.savefig(figureFile)
            plt.close(fig)
########################################################################


##################################################################
#Only this method supports simulations
#Only this method supports simulations
def plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,sample,numberofMutations,N,x_axis_labels,strand1_values,strand2_values,strand1_simulations_median_values ,strand2_simulations_median_values , fdr_bh_adjusted_pvalues, strand1Name, strand2Name, mutationsOrSignatures, color1, color2, figureName, width):
    # print('############# for debug starts Nov 12, 2018 ################')
    # print('Sample:%s --- Strand1:%s --- Strand2:%s --- mutationsOrSignatures:%s' %(sample,strand1Name,strand2Name,mutationsOrSignatures))
    # print('strand1_values: %s' %(strand1_values))
    # print('strand2_values: %s' %(strand2_values))
    # print('strand1_simulations_median_values: %s' %(strand1_simulations_median_values))
    # print('strand2_simulations_median_values: %s' %(strand2_simulations_median_values))
    # print('fdr_bh_adjusted_pvalues: %s' %(fdr_bh_adjusted_pvalues))
    # print('############# for debug ends Nov 12, 2018 ################')


    from matplotlib import rcParams
    rcParams.update({'figure.autolayout': True})

    ind = np.arange(N)  # the x locations for the groups

    fig, ax = plt.subplots(figsize=(16,10),dpi=300)

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
    ax.set_title('%s vs. %s %s' %(strand1Name,strand2Name,mutationsOrSignatures), fontsize=35,fontweight='bold')

    #set the x axis tick locations
    if (numberofSimulations>0):
        ax.set_xticks(ind + (3 * width)/2)
        simulationsStrand1Name = 'Simulated %s' %(strand1Name)
        simulationsStrand2Name = 'Simulated %s' % (strand2Name)

        # legend = ax.legend((rects1[0], rects2[0], rects3[0], rects4[0]),
        #                    (strand1Name, strand2Name, simulationsStrand1Name, simulationsStrand2Name),
        #                    prop={'size': 25}, bbox_to_anchor = (0, 1.21), ncol = 2, loc = 'upper left')

        legend = ax.legend((rects1[0], rects2[0], rects3[0], rects4[0]),
                           (strand1Name, strand2Name, simulationsStrand1Name, simulationsStrand2Name),
                           prop={'size': 25}, ncol = 2, loc = 'upper center')

    else:
        #Old way with no simulations
        ax.set_xticks(ind + width/2)
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

    if sample is None:
        figureName = '%s.png' %(figureName)
        figureFile = os.path.join(outputDir, jobname, FIGURE, ALL, STRANDBIAS, figureName)
    else:
        figureName = '%s_%s_%d.png' %(figureName,sample,numberofMutations)
        os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, STRANDBIAS), exist_ok=True)
        figureFile = os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, STRANDBIAS, figureName)

    fig.savefig(figureFile)
    plt.close(fig)
##################################################################


##################################################################
# Fill strand2TypeCountListDict using type2Strand2CountDict w.r.t. referenceTypes
def fillDictWRTReferenceTypes(referenceTypes,strands,type2Strand2CountDict):
    strand2TypeCountListDict= {}

    for strandType in strands:
        strand2TypeCountListDict[strandType] = []
        for type in referenceTypes:
            if type in type2Strand2CountDict:
                if strandType in type2Strand2CountDict[type]:
                    count = type2Strand2CountDict[type][strandType]
                    strand2TypeCountListDict[strandType].append(count)
                else:
                    strand2TypeCountListDict[strandType].append(0)
            else:
                strand2TypeCountListDict[strandType].append(0)

    return strand2TypeCountListDict
##################################################################

##################################################################
def fillSimulationsType2StrandCountList(
        outputDir,
        jobname,
        existingTypesList,
        numberofSimulations,
        strandbias,
        type2Strand2CountDict_Filename,
        type):

    all_simulations_types_strand1_list = []
    all_simulations_types_strand2_list = []

    simulations_types_strand1_medians_list = []
    simulations_types_strand2_medians_list = []

    ##############################################################
    for simNum in range(1, numberofSimulations + 1):
        simJobName = '%s_Sim%d' % (jobname, simNum)
        simulationType2StrandCountFilePath = os.path.join(outputDir,simJobName,DATA,strandbias,type2Strand2CountDict_Filename)

        simulationType2StrandCountDict = readDictionary(simulationType2StrandCountFilePath)

        if ((simulationType2StrandCountDict is not None) and (simulationType2StrandCountDict)):
            if (type == SIGNATURE):
                # Since it contains mutationprobability
                print('For debug, FEB 23, 2019 starts')
                print('simJobName')
                print(simJobName)
                print('simulationType2StrandCountDict')
                print(simulationType2StrandCountDict)
                print('For debug, FEB 23, 2019 ends')
                simulationType2StrandCountDict = simulationType2StrandCountDict[SUBSTITUTION_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD]

            # Let's fill simulations w.r.t.  existingTypesList in the real data
            # if (simulationType2StrandCountDict is not None):
            if (strandbias == TRANSCRIPTIONSTRANDBIAS):
                strand2TypeStrandCountList_Dict = fillDictWRTReferenceTypes(existingTypesList,transcriptionStrands,simulationType2StrandCountDict)
                simulation_types_strand1_list = strand2TypeStrandCountList_Dict[TRANSCRIBED_STRAND]
                simulation_types_strand2_list = strand2TypeStrandCountList_Dict[UNTRANSCRIBED_STRAND]

            elif (strandbias == REPLICATIONSTRANDBIAS):
                strand2TypeStrandCountList_Dict = fillDictWRTReferenceTypes(existingTypesList,replicationStrands,simulationType2StrandCountDict)
                simulation_types_strand1_list = strand2TypeStrandCountList_Dict[LAGGING]
                simulation_types_strand2_list = strand2TypeStrandCountList_Dict[LEADING]


            all_simulations_types_strand1_list.append(simulation_types_strand1_list)
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
def fillSimulationsSample2Type2StrandCountList(
        outputDir,
        jobname,
        existingTypesList,
        numberofSimulations,
        samplesWithAtLeast10KMutations2NumberofMutationsDict,
        strandbias,
        type2Sample2Strand2CountDict_Filename,
        type):

    sample2SimulationsTypesStrand1MediansListDict = {}
    sample2SimulationsTypesStrand2MediansListDict = {}

    sample2AllSimulationsTypesStrand1ListDict = {}
    sample2AllSimulationsTypesStrand2ListDict = {}

    # Fill sample2AllSimulationsMutationTypesTranscribedListDict
    # Fill sample2AllSimulationsMutationTypesUntranscribedListDict
    ##############################################################
    for sample in samplesWithAtLeast10KMutations2NumberofMutationsDict:
        sample2AllSimulationsTypesStrand1ListDict[sample] = []
        sample2AllSimulationsTypesStrand2ListDict[sample] = []

        for simNum in range(1, numberofSimulations+1):
            simJobName = '%s_Sim%d' % (jobname, simNum)

            simulationbased_type2Sample2StrandBias2CountFilePath = os.path.join(outputDir, simJobName, DATA,strandbias, type2Sample2Strand2CountDict_Filename)
            simulationBased_type2Sample2Strand2CountDict = readDictionary(simulationbased_type2Sample2StrandBias2CountFilePath)

            if ((simulationBased_type2Sample2Strand2CountDict is not None) and (simulationBased_type2Sample2Strand2CountDict)):

                if (type == SIGNATURE):
                    #Since it contains mutationprobability
                    simulationBased_type2Sample2Strand2CountDict = simulationBased_type2Sample2Strand2CountDict[SUBSTITUTION_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD]

                # if (simulationBased_type2Sample2Strand2CountDict is not None):
                simulationBased_sample2Type2Strand2CountDict = {}
                convert(simulationBased_type2Sample2Strand2CountDict,simulationBased_sample2Type2Strand2CountDict)

                if (strandbias == TRANSCRIPTIONSTRANDBIAS):
                    strand2TypeStrandCountList_Dict = fillDictWRTReferenceTypes(existingTypesList,transcriptionStrands,simulationBased_sample2Type2Strand2CountDict[sample])
                    sampleBased_types_strand1_list = strand2TypeStrandCountList_Dict[TRANSCRIBED_STRAND]
                    sampleBased_types_strand2_list = strand2TypeStrandCountList_Dict[UNTRANSCRIBED_STRAND]
                elif (strandbias == REPLICATIONSTRANDBIAS):
                    strand2TypeStrandCountList_Dict = fillDictWRTReferenceTypes(existingTypesList,replicationStrands,simulationBased_sample2Type2Strand2CountDict[sample])
                    sampleBased_types_strand1_list = strand2TypeStrandCountList_Dict[LAGGING]
                    sampleBased_types_strand2_list = strand2TypeStrandCountList_Dict[LEADING]

                sample2AllSimulationsTypesStrand1ListDict[sample].append(sampleBased_types_strand1_list)
                sample2AllSimulationsTypesStrand2ListDict[sample].append(sampleBased_types_strand2_list)
            ##############################################################
        #################################################################################

    ##############################################################
    #After filling the all simulations list calculate the medians for each sample
    # Fill sample2SimulationsMutationTypesTranscribedMediansListDict
    # Fill sample2SimulationsMutationTypesUntranscribedMediansListDict
    ##############################################################
    for sample in samplesWithAtLeast10KMutations2NumberofMutationsDict :

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
def convert(firstType2SecondType2Strand2CountDict, secondType2FirstType2Strand2CountDict):
    for firstType, secondType2Strand2CountDict in firstType2SecondType2Strand2CountDict.items():
        for secondType, strand2CountDict in secondType2Strand2CountDict.items():
            if secondType not in secondType2FirstType2Strand2CountDict:
                secondType2FirstType2Strand2CountDict[secondType] = {}
                secondType2FirstType2Strand2CountDict[secondType][firstType] = strand2CountDict
            else:
                secondType2FirstType2Strand2CountDict[secondType][firstType] = strand2CountDict
##################################################################


##################################################################
#Fill these dictionaries required for the bar plots
def fillSample2TypesPValuesDict(strands,
                                existingTypesList,
                                sample2Type2Strand2CountDict,
                                sample2SimulationsTypesStrand1MediansListDict,
                                sample2SimulationsTypesStrand2MediansListDict):

    sample2TypesStrand1CountListDict = {}
    sample2TypesStrand2CountListDict = {}
    sample2TypesStrandPValuesListDict = {}

    for sample in sample2Type2Strand2CountDict:
        strandType_MutationTypeStrandCountList_Dict = fillDictWRTReferenceTypes(existingTypesList,strands,sample2Type2Strand2CountDict[sample])

        sampleBased_mutationtypes_strand1_list = strandType_MutationTypeStrandCountList_Dict[strands[0]]
        sampleBased_mutationtypes_strand2_list = strandType_MutationTypeStrandCountList_Dict[strands[1]]

        sample2TypesStrand1CountListDict[sample] = sampleBased_mutationtypes_strand1_list
        sample2TypesStrand2CountListDict[sample] = sampleBased_mutationtypes_strand2_list

        simulations_sampleBased_mutationtypes_transcribed_list = None
        simulations_sampleBased_mutationtypes_untranscribed_list = None
        #Now we have
        if ((sample2SimulationsTypesStrand1MediansListDict is not None) and (sample in sample2SimulationsTypesStrand1MediansListDict)):
            simulations_sampleBased_mutationtypes_transcribed_list = sample2SimulationsTypesStrand1MediansListDict[sample]

        if  ((sample2SimulationsTypesStrand2MediansListDict is not None ) and (sample in sample2SimulationsTypesStrand2MediansListDict)):
            simulations_sampleBased_mutationtypes_untranscribed_list = sample2SimulationsTypesStrand2MediansListDict[sample]


        # Please notice Fisher exact test does not return anything for large numbers therefore chi2 test is used.
        samplebased_type_strand_pvalues = calculate_p_values(sampleBased_mutationtypes_strand1_list,
                                                                simulations_sampleBased_mutationtypes_transcribed_list,
                                                                sampleBased_mutationtypes_strand2_list,
                                                                simulations_sampleBased_mutationtypes_untranscribed_list)

        sample2TypesStrandPValuesListDict[sample] = samplebased_type_strand_pvalues

    return sample2TypesStrand1CountListDict, sample2TypesStrand2CountListDict, sample2TypesStrandPValuesListDict
##################################################################



##################################################################
def plotBarPlots(outputDir,jobname,numberofSimulations,sample2NumberofSubsDict,existingTranscriptionMutationTypesList,sample2SimulationsMutationTypesTranscribedMediansListDict,sample2SimulationsMutationTypesUntranscribedMediansListDict,sample2MutationType2TranscriptionStrand2CountDict,sample2MutationTypesTranscribedCountListDict,sample2MutationTypesUntranscribedCountListDict,sample2MutationTypesTranscription_BH_FDR_Adjusted_PValues_Dict,width,strands,color1,color2,title,figureName):

    #######################################################
    for sample in sample2MutationType2TranscriptionStrand2CountDict:
        mutationtype_2_transcribed_list = sample2MutationTypesTranscribedCountListDict[sample]
        mutationtype_2_untranscribed_list = sample2MutationTypesUntranscribedCountListDict[sample]

        x_axis_labels = existingTranscriptionMutationTypesList
        N = len(x_axis_labels)

        if ((sample in sample2NumberofSubsDict) and (sample in sample2MutationTypesTranscription_BH_FDR_Adjusted_PValues_Dict)):
            sample_mutationtypes_transcription_FDR_BH_adjusted_pvalues = sample2MutationTypesTranscription_BH_FDR_Adjusted_PValues_Dict[sample]
            numberofMutations = sample2NumberofSubsDict[sample]

            samplebased_simulations_mutationtypes_transcribed_medians_list= None
            samplebased_simulations_mutationtypes_untranscribed_medians_list = None

            if ((sample2SimulationsMutationTypesTranscribedMediansListDict is not None) and (sample in sample2SimulationsMutationTypesTranscribedMediansListDict)):
                samplebased_simulations_mutationtypes_transcribed_medians_list = sample2SimulationsMutationTypesTranscribedMediansListDict[sample]
            if ((sample2SimulationsMutationTypesUntranscribedMediansListDict is not None) and (sample in sample2SimulationsMutationTypesUntranscribedMediansListDict)):
                samplebased_simulations_mutationtypes_untranscribed_medians_list = sample2SimulationsMutationTypesUntranscribedMediansListDict[sample]

            plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,sample,numberofMutations,N,x_axis_labels,mutationtype_2_transcribed_list,mutationtype_2_untranscribed_list,samplebased_simulations_mutationtypes_transcribed_medians_list,samplebased_simulations_mutationtypes_untranscribed_medians_list,sample_mutationtypes_transcription_FDR_BH_adjusted_pvalues,strands[0],strands[1],title,color1,color2,figureName,width)
    #######################################################


##################################################################

# ##################################################################
#main function
def transcriptionReplicationStrandBiasFigures(outputDir,jobname,figureAugmentation,numberofSimulations):

    isFigureAugmentation = False
    if (figureAugmentation == 'augmentation'):
        isFigureAugmentation = True

    #######################################################################################################################
    os.makedirs(os.path.join(outputDir, jobname, FIGURE, ALL, STRANDBIAS), exist_ok=True)
    os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES), exist_ok=True)
    #######################################################################################################################

    ##########################################################################################
    #########################  Read dictionaries related with ################################
    #########################  signatures and samples starts  ################################
    ##########################################################################################
    sample2NumberofSubsDict = getSample2NumberofSubsDict(outputDir,jobname)
    sample2NumberofIndelsDict = getSample2NumberofIndelsDict(outputDir,jobname)
    subsSignature2NumberofMutationsDict = getSubsSignature2NumberofMutationsDict(outputDir,jobname)
    indelsSignature2NumberofMutationsDict = getIndelsSignature2NumberofMutationsDict(outputDir,jobname)
    sample2SubsSignature2NumberofMutationsDict = getSample2SubsSignature2NumberofMutationsDict(outputDir,jobname)
    sample2IndelsSignature2NumberofMutationsDict = getSample2IndelsSignature2NumberofMutationsDict(outputDir,jobname)

    subsSignatures= subsSignature2NumberofMutationsDict.keys()
    indelsSignatures = indelsSignature2NumberofMutationsDict.keys()
    ##########################################################################################
    #########################  Read dictionaries related with ################################
    #########################  signatures and samples ends  ##################################
    ##########################################################################################

    ##########################################################################################
    #####################  Read dictionaries for the rest part starts ########################
    ##########################################################################################

    ##########################################################################################
    mutationType2TranscriptionStrandCountFilePath = os.path.join(outputDir,jobname,DATA,TRANSCRIPTIONSTRANDBIAS,MutationType2TranscriptionStrand2CountDict_Filename)
    mutationType2TranscriptionStrandCountDict = readDictionary(mutationType2TranscriptionStrandCountFilePath)

    mutationType2ReplicationStrandCountFilePath = os.path.join(outputDir,jobname,DATA,REPLICATIONSTRANDBIAS,MutationType2ReplicationStrand2CountDict_Filename)
    mutationType2ReplicationStrandCountDict = readDictionary(mutationType2ReplicationStrandCountFilePath)

    mutationType2Sample2TranscriptionStrand2CountFilePath = os.path.join(outputDir, jobname, DATA,TRANSCRIPTIONSTRANDBIAS,MutationType2Sample2TranscriptionStrand2CountDict_Filename)
    mutationType2Sample2TranscriptionStrand2CountDict = readDictionary(mutationType2Sample2TranscriptionStrand2CountFilePath)

    mutationType2Sample2ReplicationStrand2CountFilePath = os.path.join(outputDir,jobname,DATA,REPLICATIONSTRANDBIAS,MutationType2Sample2ReplicationStrand2CountDict_Filename)
    mutationType2Sample2ReplicationStrand2CountDict = readDictionary(mutationType2Sample2ReplicationStrand2CountFilePath)
    ##########################################################################################

    ##########################################################################################
    mutationProbabilitySubs2Signature2TranscriptionStrand2CountFilePath = os.path.join(outputDir,jobname,DATA,TRANSCRIPTIONSTRANDBIAS,MutationProbability2SubsSignature2TranscriptionStrand2CountDict_Filename)
    mutationProbability2SubsSignature2TranscriptionStrand2CountDict = readDictionary(mutationProbabilitySubs2Signature2TranscriptionStrand2CountFilePath)

    mutationProbability2SubsSignature2ReplicationStrand2CountFilePath = os.path.join(outputDir, jobname,DATA, REPLICATIONSTRANDBIAS,MutationProbability2SubsSignature2ReplicationStrand2CountDict_Filename)
    mutationProbability2SubsSignature2ReplicationStrand2CountDict = readDictionary(mutationProbability2SubsSignature2ReplicationStrand2CountFilePath)

    mutationProbability2SubsSignature2Sample2TranscriptionStrand2CountFilePath = os.path.join(outputDir,jobname, DATA,TRANSCRIPTIONSTRANDBIAS,MutationProbability2SubsSignature2Sample2TranscriptionStrand2CountDict_Filename)
    mutationProbability2SubsSignature2Sample2TranscriptionStrand2CountDict = readDictionary(mutationProbability2SubsSignature2Sample2TranscriptionStrand2CountFilePath)

    mutationProbability2SubsSignature2Sample2ReplicationStrand2CountFilePath = os.path.join(outputDir,jobname,DATA,REPLICATIONSTRANDBIAS,MutationProbability2SubsSignature2Sample2ReplicationStrand2CountDict_Filename)
    mutationProbability2SubsSignature2Sample2ReplicationStrand2CountDict = readDictionary(mutationProbability2SubsSignature2Sample2ReplicationStrand2CountFilePath)
    ##########################################################################################

    ##########################################################################################
    #IndelsSignature starts
    mutationProbabilityIndels2Signature2TranscriptionStrand2CountFilePath = os.path.join(outputDir,jobname,DATA,TRANSCRIPTIONSTRANDBIAS,MutationProbability2IndelsSignature2TranscriptionStrand2CountDict_Filename)
    mutationProbability2IndelsSignature2TranscriptionStrand2CountDict = readDictionary(mutationProbabilityIndels2Signature2TranscriptionStrand2CountFilePath)

    mutationProbability2IndelsSignature2ReplicationStrand2CountFilePath = os.path.join(outputDir, jobname,DATA, REPLICATIONSTRANDBIAS,MutationProbability2IndelsSignature2ReplicationStrand2CountDict_Filename)
    mutationProbability2IndelsSignature2ReplicationStrand2CountDict = readDictionary(mutationProbability2IndelsSignature2ReplicationStrand2CountFilePath)

    mutationProbability2IndelsSignature2Sample2TranscriptionStrand2CountFilePath = os.path.join(outputDir,jobname, DATA,TRANSCRIPTIONSTRANDBIAS,MutationProbability2IndelsSignature2Sample2TranscriptionStrand2CountDict_Filename)
    mutationProbability2IndelsSignature2Sample2TranscriptionStrand2CountDict = readDictionary(mutationProbability2IndelsSignature2Sample2TranscriptionStrand2CountFilePath)

    mutationProbability2IndelsSignature2Sample2ReplicationStrand2CountFilePath = os.path.join(outputDir,jobname,DATA,REPLICATIONSTRANDBIAS,MutationProbability2IndelsSignature2Sample2ReplicationStrand2CountDict_Filename)
    mutationProbability2IndelsSignature2Sample2ReplicationStrand2CountDict = readDictionary(mutationProbability2IndelsSignature2Sample2ReplicationStrand2CountFilePath)
    #IndelsSignature ends
    ##########################################################################################

    #Initialize as empty lists
    subsSignature2Sample2TranscriptionStrand2CountDict = {}
    subsSignature2Sample2ReplicationStrand2CountDict = {}
    subsSignature2TranscriptionStrand2CountDict = {}
    subsSignature2ReplicationStrand2CountDict = {}

    #Initialize as empty lists
    indelsSignature2Sample2TranscriptionStrand2CountDict = {}
    indelsSignature2Sample2ReplicationStrand2CountDict = {}
    indelsSignature2TranscriptionStrand2CountDict = {}
    indelsSignature2ReplicationStrand2CountDict = {}

    if ((mutationProbability2SubsSignature2Sample2TranscriptionStrand2CountDict is not None) and (mutationProbability2SubsSignature2Sample2TranscriptionStrand2CountDict)):
        subsSignature2Sample2TranscriptionStrand2CountDict = mutationProbability2SubsSignature2Sample2TranscriptionStrand2CountDict[SUBSTITUTION_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING]

    if ((mutationProbability2SubsSignature2Sample2ReplicationStrand2CountDict) is not None and (mutationProbability2SubsSignature2Sample2ReplicationStrand2CountDict)):
        subsSignature2Sample2ReplicationStrand2CountDict = mutationProbability2SubsSignature2Sample2ReplicationStrand2CountDict[SUBSTITUTION_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING]

    if ((mutationProbability2SubsSignature2TranscriptionStrand2CountDict is not None) and (mutationProbability2SubsSignature2TranscriptionStrand2CountDict)):
        subsSignature2TranscriptionStrand2CountDict = mutationProbability2SubsSignature2TranscriptionStrand2CountDict[SUBSTITUTION_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING]

    if ((mutationProbability2SubsSignature2ReplicationStrand2CountDict is not None) and (mutationProbability2SubsSignature2ReplicationStrand2CountDict)):
        subsSignature2ReplicationStrand2CountDict = mutationProbability2SubsSignature2ReplicationStrand2CountDict[SUBSTITUTION_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING]

    if ((mutationProbability2IndelsSignature2Sample2TranscriptionStrand2CountDict is not None) and (mutationProbability2IndelsSignature2Sample2TranscriptionStrand2CountDict)):
        indelsSignature2Sample2TranscriptionStrand2CountDict = mutationProbability2IndelsSignature2Sample2TranscriptionStrand2CountDict[INDEL_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING]

    if ((mutationProbability2IndelsSignature2Sample2ReplicationStrand2CountDict) is not None and (mutationProbability2IndelsSignature2Sample2ReplicationStrand2CountDict)):
        indelsSignature2Sample2ReplicationStrand2CountDict = mutationProbability2IndelsSignature2Sample2ReplicationStrand2CountDict[INDEL_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING]

    if ((mutationProbability2IndelsSignature2TranscriptionStrand2CountDict is not None) and (mutationProbability2IndelsSignature2TranscriptionStrand2CountDict)):
        indelsSignature2TranscriptionStrand2CountDict = mutationProbability2IndelsSignature2TranscriptionStrand2CountDict[INDEL_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING]

    if ((mutationProbability2IndelsSignature2ReplicationStrand2CountDict is not None) and (mutationProbability2IndelsSignature2ReplicationStrand2CountDict)):
        indelsSignature2ReplicationStrand2CountDict = mutationProbability2IndelsSignature2ReplicationStrand2CountDict[INDEL_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING]


    existingTranscriptionMutationTypesList = sorted(list(mutationType2TranscriptionStrandCountDict.keys()))
    existingReplicationMutationTypesList = sorted(list(mutationType2ReplicationStrandCountDict.keys()))

    existingTranscriptionSubsSignatureTypesList = sorted(list(subsSignature2TranscriptionStrand2CountDict.keys()))
    existingReplicationSubsSignatureTypesList = sorted(list(subsSignature2ReplicationStrand2CountDict.keys()))

    existingTranscriptionIndelsSignatureTypesList = sorted(list(indelsSignature2TranscriptionStrand2CountDict.keys()))
    existingReplicationIndelsSignatureTypesList = sorted(list(indelsSignature2ReplicationStrand2CountDict.keys()))

    ##########################################################################################
    #####################  Read dictionaries for the rest part ends ##########################
    ##########################################################################################

    ########################################################################
    ##########################  Part 1 starts ##############################
    ########### MutationType All Samples Pooled Figures starts #############
    ######## SignatureType All Samples Pooled Figures starts ###############
    ########################################################################
    # Note: Since these are already in sample based no extra work is required.
    #MutationType SampleBased Figures --- No numberofMutations
    if ((mutationType2Sample2TranscriptionStrand2CountDict is not None) and (mutationType2Sample2ReplicationStrand2CountDict is not None)):
        plot_ncomms11383_Supp_FigE_MutationTypeBased_AllSamples_TranscriptionLog10Ratio_ReplicationLog10Ratio(mutationType2Sample2TranscriptionStrand2CountDict,mutationType2Sample2ReplicationStrand2CountDict,outputDir,jobname,isFigureAugmentation)

    #SubsSignature Sample Based Figures --- No numberofMutations
    if ((mutationProbability2SubsSignature2Sample2TranscriptionStrand2CountDict is not None) and (mutationProbability2SubsSignature2Sample2ReplicationStrand2CountDict is not None)):
        plot_ncomms11383_Supp_FigF_SignatureBased_AllSamples_TranscriptionLog10Ratio_ReplicationLog10Ratio(
            mutationProbability2SubsSignature2Sample2TranscriptionStrand2CountDict,
            mutationProbability2SubsSignature2Sample2ReplicationStrand2CountDict,
            subsSignatures,outputDir,jobname,isFigureAugmentation,SUBSTITUTION_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING)


    #IndelsSignature Sample Based Figures --- No numberofMutations
    if ((mutationProbability2IndelsSignature2Sample2TranscriptionStrand2CountDict is not None) and (mutationProbability2IndelsSignature2Sample2ReplicationStrand2CountDict is not None)):
        plot_ncomms11383_Supp_FigF_SignatureBased_AllSamples_TranscriptionLog10Ratio_ReplicationLog10Ratio(
            mutationProbability2IndelsSignature2Sample2TranscriptionStrand2CountDict,
            mutationProbability2IndelsSignature2Sample2ReplicationStrand2CountDict,
            indelsSignatures,outputDir,jobname,isFigureAugmentation,INDEL_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING)
    ########################################################################
    ########### MutationType All Samples Pooled Figures ends ###############
    ######## SignatureType All Samples Pooled Figures ends #################
    ##########################  Part 1 ends ################################
    ########################################################################


    ########################################################################
    ##########################  Part 2 starts ##############################
    ############## MutationType Scatter Plots starts #######################
    ############## Signature Based Scatter Plots starts ####################
    ########################################################################
    if ((mutationType2TranscriptionStrandCountDict is not None) and (mutationType2ReplicationStrandCountDict is not None)):
        plot_ncomms11383_Supp_FigG_AllMutationTypes_TranscriptionLog10Ratio_ReplicationLog10Ratio(None,None,mutationType2TranscriptionStrandCountDict, mutationType2ReplicationStrandCountDict,outputDir,jobname,isFigureAugmentation)

    if ((mutationProbability2SubsSignature2TranscriptionStrand2CountDict is not None) and (mutationProbability2SubsSignature2ReplicationStrand2CountDict is not None) and (subsSignatures)):
        subsSignature2TranscriptionStrand2CountDict = mutationProbability2SubsSignature2TranscriptionStrand2CountDict[SUBSTITUTION_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING]
        subsSignature2ReplicationStrand2CountDict = mutationProbability2SubsSignature2ReplicationStrand2CountDict[SUBSTITUTION_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING]
        plot_ncomms11383_Supp_FigH_AllSignatures_TranscriptionLog10Ratio_ReplicationLog10Ratio('subs',None,None,subsSignature2TranscriptionStrand2CountDict,subsSignature2ReplicationStrand2CountDict,subsSignature2NumberofMutationsDict,outputDir,jobname,isFigureAugmentation)

    if ((mutationProbability2IndelsSignature2TranscriptionStrand2CountDict is not None) and (mutationProbability2IndelsSignature2ReplicationStrand2CountDict is not None) and (indelsSignatures)):
        indelsSignature2TranscriptionStrand2CountDict = mutationProbability2IndelsSignature2TranscriptionStrand2CountDict[INDEL_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING]
        indelsSignature2ReplicationStrand2CountDict = mutationProbability2IndelsSignature2ReplicationStrand2CountDict[INDEL_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING]
        plot_ncomms11383_Supp_FigH_AllSignatures_TranscriptionLog10Ratio_ReplicationLog10Ratio('indels',None,None,indelsSignature2TranscriptionStrand2CountDict,indelsSignature2ReplicationStrand2CountDict,indelsSignature2NumberofMutationsDict,outputDir,jobname,isFigureAugmentation)
    ########################################################################
    ############## MutationType Scatter Plots ends #########################
    ############## Signature Based Scatter Plots ends ######################
    ##########################  Part 2 ends ################################
    ########################################################################


    ########################################################################
    ##########################  Part 3 starts ##############################
    ###################     Sample Based    ################################
    ############## MutationType Scatter Plots starts #######################
    ############## Signature Based Scatter Plots starts ####################
    ########################################################################

    ###############################################################
    sample2MutationType2TranscriptionStrand2CountDict = {}
    sample2MutationType2ReplicationStrand2CountDict = {}

    convert(mutationType2Sample2TranscriptionStrand2CountDict, sample2MutationType2TranscriptionStrand2CountDict)
    convert(mutationType2Sample2ReplicationStrand2CountDict, sample2MutationType2ReplicationStrand2CountDict)

    for sample in sample2NumberofSubsDict:
        if ((sample in sample2MutationType2TranscriptionStrand2CountDict) and (sample in sample2MutationType2ReplicationStrand2CountDict)):
            numberofMutations = sample2NumberofSubsDict[sample]
            plot_ncomms11383_Supp_FigG_AllMutationTypes_TranscriptionLog10Ratio_ReplicationLog10Ratio(sample,numberofMutations,sample2MutationType2TranscriptionStrand2CountDict[sample], sample2MutationType2ReplicationStrand2CountDict[sample],outputDir,jobname,isFigureAugmentation)
    ###############################################################

    ###############################################################
    sample2SubsSignature2TranscriptionStrand2CountDict = {}
    sample2SubsSignature2ReplicationStrand2CountDict = {}

    convert(subsSignature2Sample2TranscriptionStrand2CountDict, sample2SubsSignature2TranscriptionStrand2CountDict)
    convert(subsSignature2Sample2ReplicationStrand2CountDict, sample2SubsSignature2ReplicationStrand2CountDict)

    for sample in sample2NumberofSubsDict:
        if ((sample in sample2SubsSignature2TranscriptionStrand2CountDict) and (sample in sample2SubsSignature2ReplicationStrand2CountDict)):
            numberofMutations = sample2NumberofSubsDict[sample]
            plot_ncomms11383_Supp_FigH_AllSignatures_TranscriptionLog10Ratio_ReplicationLog10Ratio('subs',sample,numberofMutations,sample2SubsSignature2TranscriptionStrand2CountDict[sample], sample2SubsSignature2ReplicationStrand2CountDict[sample], subsSignature2NumberofMutationsDict, outputDir,jobname,isFigureAugmentation)
    ###############################################################


    ###############################################################
    sample2IndelsSignature2TranscriptionStrand2CountDict = {}
    sample2IndelsSignature2ReplicationStrand2CountDict = {}

    convert(indelsSignature2Sample2TranscriptionStrand2CountDict, sample2IndelsSignature2TranscriptionStrand2CountDict)
    convert(indelsSignature2Sample2ReplicationStrand2CountDict, sample2IndelsSignature2ReplicationStrand2CountDict)

    for sample in sample2NumberofIndelsDict:
        if ((sample in sample2IndelsSignature2TranscriptionStrand2CountDict) and (sample in sample2IndelsSignature2ReplicationStrand2CountDict)):
            numberofMutations = sample2NumberofIndelsDict[sample]
            plot_ncomms11383_Supp_FigH_AllSignatures_TranscriptionLog10Ratio_ReplicationLog10Ratio('indels',sample,numberofMutations,sample2IndelsSignature2TranscriptionStrand2CountDict[sample], sample2IndelsSignature2ReplicationStrand2CountDict[sample], indelsSignature2NumberofMutationsDict, outputDir,jobname,isFigureAugmentation)
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

    if (numberofSimulations > 0):
        #########################################################################################################################
        # simulations --- mutation type --- transcription
        simulations_mutationtypes_transcribed_medians_list, simulations_mutationtypes_untranscribed_medians_list =\
            fillSimulationsType2StrandCountList(outputDir,
                                                jobname,
                                                existingTranscriptionMutationTypesList,
                                                numberofSimulations,
                                                TRANSCRIPTIONSTRANDBIAS,
                                                MutationType2TranscriptionStrand2CountDict_Filename,
                                                MUTATION)

        # simulations --- mutation type --- replication
        simulations_mutationtypes_lagging_medians_list, simulations_mutationtypes_leading_medians_list =\
            fillSimulationsType2StrandCountList(outputDir,
                                                jobname,
                                                existingReplicationMutationTypesList,
                                                numberofSimulations,
                                                REPLICATIONSTRANDBIAS,
                                                MutationType2ReplicationStrand2CountDict_Filename,
                                                MUTATION)

        # simulations --- subs signatures --- transcription
        simulations_subs_signatures_transcribed_medians_list, simulations_subs_signatures_untranscribed_medians_list =\
            fillSimulationsType2StrandCountList(outputDir,
                                                jobname,
                                                existingTranscriptionSubsSignatureTypesList,
                                                numberofSimulations,
                                                TRANSCRIPTIONSTRANDBIAS,
                                                MutationProbability2SubsSignature2TranscriptionStrand2CountDict_Filename,
                                                SIGNATURE)

        # simulations --- subs signatures --- replication
        simulations_subs_signatures_lagging_medians_list, simulations_subs_signatures_leading_medians_list = \
            fillSimulationsType2StrandCountList(outputDir,
                                                jobname,
                                                existingReplicationSubsSignatureTypesList,
                                                numberofSimulations,
                                                REPLICATIONSTRANDBIAS,
                                                MutationProbability2SubsSignature2ReplicationStrand2CountDict_Filename,
                                                SIGNATURE)

        # simulations --- indels signatures --- transcription
        simulations_indels_signatures_transcribed_medians_list, simulations_indels_signatures_untranscribed_medians_list =\
            fillSimulationsType2StrandCountList(outputDir,
                                                jobname,
                                                existingTranscriptionIndelsSignatureTypesList,
                                                numberofSimulations,
                                                TRANSCRIPTIONSTRANDBIAS,
                                                MutationProbability2IndelsSignature2TranscriptionStrand2CountDict_Filename,
                                                SIGNATURE)

        # simulations --- indels signatures --- replication
        simulations_indels_signatures_lagging_medians_list, simulations_indels_signatures_leading_medians_list = \
            fillSimulationsType2StrandCountList(outputDir,
                                                jobname,
                                                existingReplicationIndelsSignatureTypesList,
                                                numberofSimulations,
                                                REPLICATIONSTRANDBIAS,
                                                MutationProbability2IndelsSignature2ReplicationStrand2CountDict_Filename,
                                                SIGNATURE)
        #########################################################################################################################


        #########################################################################################################################
        #sample based starts
        # samplebased --- simulations --- mutation type --- transcription
        sample2SimulationsMutationTypesTranscribedMediansListDict, sample2SimulationsMutationTypesUntranscribedMediansListDict =\
            fillSimulationsSample2Type2StrandCountList(outputDir,jobname,
                existingTranscriptionMutationTypesList,
                numberofSimulations,
                sample2NumberofSubsDict,
                TRANSCRIPTIONSTRANDBIAS,
                MutationType2Sample2TranscriptionStrand2CountDict_Filename,
                MUTATION)

        # samplebased --- simulations --- mutation type --- replication
        sample2SimulationsMutationTypesLaggingMediansListDict, sample2SimulationsMutationTypesLeadingMediansListDict =\
            fillSimulationsSample2Type2StrandCountList(outputDir,jobname,
                existingReplicationMutationTypesList,
                numberofSimulations,
                sample2NumberofSubsDict,
                REPLICATIONSTRANDBIAS,
                MutationType2Sample2ReplicationStrand2CountDict_Filename,
                MUTATION)

        # samplebased --- simulations --- signature --- transcription
        sample2SimulationsSubsSignatureTypesTranscribedMediansListDict, sample2SimulationsSubsSignatureTypesUntranscribedMediansListDict =\
            fillSimulationsSample2Type2StrandCountList(outputDir,jobname,
                existingTranscriptionSubsSignatureTypesList,
                numberofSimulations,
                sample2NumberofSubsDict,
                TRANSCRIPTIONSTRANDBIAS,
                MutationProbability2SubsSignature2Sample2TranscriptionStrand2CountDict_Filename,
                SIGNATURE)

        # samplebased --- simulations --- signature --- replication
        sample2SimulationsSubsSignatureTypesLaggingMediansListDict, sample2SimulationsSubsSignatureTypesLeadingMediansListDict =\
            fillSimulationsSample2Type2StrandCountList(outputDir,jobname,
                existingReplicationSubsSignatureTypesList,
                numberofSimulations,
                sample2NumberofSubsDict,
                REPLICATIONSTRANDBIAS,
                MutationProbability2SubsSignature2Sample2ReplicationStrand2CountDict_Filename,
                SIGNATURE)

        # samplebased --- simulations --- signature --- transcription
        sample2SimulationsIndelsSignatureTypesTranscribedMediansListDict, sample2SimulationsIndelsSignatureTypesUntranscribedMediansListDict =\
            fillSimulationsSample2Type2StrandCountList(outputDir,jobname,
                existingTranscriptionIndelsSignatureTypesList,
                numberofSimulations,
                sample2NumberofIndelsDict,
                TRANSCRIPTIONSTRANDBIAS,
                MutationProbability2IndelsSignature2Sample2TranscriptionStrand2CountDict_Filename,
                SIGNATURE)

        # samplebased --- simulations --- signature --- replication
        sample2SimulationsIndelsSignatureTypesLaggingMediansListDict, sample2SimulationsIndelsSignatureTypesLeadingMediansListDict =\
            fillSimulationsSample2Type2StrandCountList(outputDir,jobname,
                existingReplicationIndelsSignatureTypesList,
                numberofSimulations,
                sample2NumberofIndelsDict,
                REPLICATIONSTRANDBIAS,
                MutationProbability2IndelsSignature2Sample2ReplicationStrand2CountDict_Filename,
                SIGNATURE)
        #sample based ends
        #########################################################################################################################

    ##########################################################################################################
    ################################### Fill for simulations ends ############################################
    ##########################################################################################################


    #########################################################################################################
    #############################  Calculate p-values for Mutation Type Starts    ###########################
    #########################################################################################################

    ########################################## Transcription starts  ########################################
    strandType_MutationType2CountList_Dict = fillDictWRTReferenceTypes(existingTranscriptionMutationTypesList,transcriptionStrands,mutationType2TranscriptionStrandCountDict)
    mutationtypes_transcribed_list = strandType_MutationType2CountList_Dict[TRANSCRIBED_STRAND]
    mutationtypes_untranscribed_list = strandType_MutationType2CountList_Dict[UNTRANSCRIBED_STRAND]

    mutationtype_transcription_pvalues = calculate_p_values(mutationtypes_transcribed_list,
                                                            simulations_mutationtypes_transcribed_medians_list,
                                                            mutationtypes_untranscribed_list,
                                                            simulations_mutationtypes_untranscribed_medians_list)
    ########################################## Transcription ends  ##########################################


    ########################################## Replication starts  ##########################################
    strandType_MutationType2CountList_Dict = fillDictWRTReferenceTypes(existingReplicationMutationTypesList,replicationStrands,mutationType2ReplicationStrandCountDict)
    mutationtypes_lagging_list = strandType_MutationType2CountList_Dict[LAGGING]
    mutationtypes_leading_list = strandType_MutationType2CountList_Dict[LEADING]

    #Please notice Fisher exact test does not return anything for large numbers therefore chi2 test is used.
    mutationtype_replication_pvalues = calculate_p_values(mutationtypes_lagging_list,
                                                          simulations_mutationtypes_lagging_medians_list,
                                                          mutationtypes_leading_list,
                                                          simulations_mutationtypes_leading_medians_list)
    ########################################## Replication ends  ############################################

    #########################################################################################################
    #############################  Calculate p-values for Mutation Type ends    #############################
    #########################################################################################################


    #########################################################################################################
    ########################  Calculate p-values for subs signatures starts    ##############################
    #########################################################################################################

    ########################################## Transcription starts  ########################################
    strandType_Signature2CountList_Dict = fillDictWRTReferenceTypes(existingTranscriptionSubsSignatureTypesList,transcriptionStrands,subsSignature2TranscriptionStrand2CountDict)
    subs_signatures_transcribedStrandCount_list = strandType_Signature2CountList_Dict[TRANSCRIBED_STRAND]
    subs_signatures_untranscribedStrandCount_list = strandType_Signature2CountList_Dict[UNTRANSCRIBED_STRAND]

    subs_signature_transcription_pvalues = calculate_p_values(subs_signatures_transcribedStrandCount_list,
                                                         simulations_subs_signatures_transcribed_medians_list,
                                                         subs_signatures_untranscribedStrandCount_list,
                                                         simulations_subs_signatures_untranscribed_medians_list)
    ########################################## Transcription ends  ##########################################

    ########################################## Replication starts  ##########################################
    strandType_Signature2CountList_Dict = fillDictWRTReferenceTypes(existingReplicationSubsSignatureTypesList,replicationStrands,subsSignature2ReplicationStrand2CountDict)
    subs_signatures_laggingStrandCount_list = strandType_Signature2CountList_Dict[LAGGING]
    subs_signatures_leadingStrandCount_list = strandType_Signature2CountList_Dict[LEADING]

    subs_signature_replication_pvalues = calculate_p_values(subs_signatures_laggingStrandCount_list,
                                                       simulations_subs_signatures_lagging_medians_list,
                                                       subs_signatures_leadingStrandCount_list,
                                                       simulations_subs_signatures_leading_medians_list)
    ########################################## Replication ends  ############################################

    #########################################################################################################
    ########################  Calculate p-values for subs signatures ends    ################################
    #########################################################################################################


    #########################################################################################################
    ######################  Calculate p-values for indels signatures ends    ################################
    #########################################################################################################

    ########################################## Transcription starts  ########################################
    strandType_Signature2CountList_Dict = fillDictWRTReferenceTypes(existingTranscriptionIndelsSignatureTypesList,transcriptionStrands,indelsSignature2TranscriptionStrand2CountDict)
    indels_signatures_transcribedStrandCount_list = strandType_Signature2CountList_Dict[TRANSCRIBED_STRAND]
    indels_signatures_untranscribedStrandCount_list = strandType_Signature2CountList_Dict[UNTRANSCRIBED_STRAND]


    indels_signature_transcription_pvalues = calculate_p_values(indels_signatures_transcribedStrandCount_list,
                                                         simulations_indels_signatures_transcribed_medians_list,
                                                         indels_signatures_untranscribedStrandCount_list,
                                                         simulations_indels_signatures_untranscribed_medians_list)
    ########################################## Transcription ends  ##########################################

    ########################################## Replication starts  ##########################################
    strandType_Signature2CountList_Dict = fillDictWRTReferenceTypes(existingReplicationIndelsSignatureTypesList,replicationStrands,indelsSignature2ReplicationStrand2CountDict)
    indels_signatures_laggingStrandCount_list = strandType_Signature2CountList_Dict[LAGGING]
    indels_signatures_leadingStrandCount_list = strandType_Signature2CountList_Dict[LEADING]

    indels_signature_replication_pvalues = calculate_p_values(indels_signatures_laggingStrandCount_list,
                                                       simulations_indels_signatures_lagging_medians_list,
                                                       indels_signatures_leadingStrandCount_list,
                                                       simulations_indels_signatures_leading_medians_list)
    ########################################## Replication ends  ############################################

    #########################################################################################################
    ######################  Calculate p-values for indels signatures ends    ################################
    #########################################################################################################


    #########################################################################################################
    ####################  Calculate p-values for sample based mutation types starts    ######################
    #########################################    Transcription    ###########################################
    ###########################################    Replication    ###########################################
    #########################################################################################################
    sample2MutationTypesTranscribedCountListDict, sample2MutationTypesUntranscribedCountListDict, sample2MutationTypesTranscriptionPValuesListDict = fillSample2TypesPValuesDict(transcriptionStrands,existingTranscriptionMutationTypesList,sample2MutationType2TranscriptionStrand2CountDict,sample2SimulationsMutationTypesTranscribedMediansListDict,sample2SimulationsMutationTypesUntranscribedMediansListDict)
    sample2MutationTypesLaggingCountListDict, sample2MutationTypesLeadingCountListDict, sample2MutationTypesReplicationPValuesListDict = fillSample2TypesPValuesDict(replicationStrands,existingReplicationMutationTypesList,sample2MutationType2ReplicationStrand2CountDict,sample2SimulationsMutationTypesLaggingMediansListDict,sample2SimulationsMutationTypesLeadingMediansListDict)
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

    ##################################################################################
    sample2SubsSignatureTypesTranscribedCountListDict, sample2SubsSignatureTypesUntranscribedCountListDict, sample2SubsSignature2TranscriptionPValuesListDict = \
        fillSample2TypesPValuesDict(transcriptionStrands,
                                    existingTranscriptionSubsSignatureTypesList,
                                    sample2SubsSignature2TranscriptionStrand2CountDict,
                                    sample2SimulationsSubsSignatureTypesTranscribedMediansListDict,
                                    sample2SimulationsSubsSignatureTypesUntranscribedMediansListDict)
    ##################################################################################

    ##################################################################################
    sample2SubsSignatureTypesLaggingCountListDict, sample2SubsSignatureTypesLeadingCountListDict, sample2SubsSignature2ReplicationPValuesListDict = \
        fillSample2TypesPValuesDict(replicationStrands,
                                    existingReplicationSubsSignatureTypesList,
                                    sample2SubsSignature2ReplicationStrand2CountDict,
                                    sample2SimulationsSubsSignatureTypesLaggingMediansListDict,
                                    sample2SimulationsSubsSignatureTypesLeadingMediansListDict)
    ##################################################################################

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

    ##################################################################################
    sample2IndelsSignatureTypesTranscribedCountListDict, sample2IndelsSignatureTypesUntranscribedCountListDict, sample2IndelsSignature2TranscriptionPValuesListDict = \
        fillSample2TypesPValuesDict(transcriptionStrands,
                                    existingTranscriptionIndelsSignatureTypesList,
                                    sample2IndelsSignature2TranscriptionStrand2CountDict,
                                    sample2SimulationsIndelsSignatureTypesTranscribedMediansListDict,
                                    sample2SimulationsIndelsSignatureTypesUntranscribedMediansListDict)
    ##################################################################################

    ##################################################################################
    sample2IndelsSignatureTypesLaggingCountListDict, sample2IndelsSignatureTypesLeadingCountListDict, sample2IndelsSignature2ReplicationPValuesListDict = \
        fillSample2TypesPValuesDict(replicationStrands,
                                    existingReplicationIndelsSignatureTypesList,
                                    sample2IndelsSignature2ReplicationStrand2CountDict,
                                    sample2SimulationsIndelsSignatureTypesLaggingMediansListDict,
                                    sample2SimulationsIndelsSignatureTypesLeadingMediansListDict)
    ##################################################################################

    #########################################################################################################
    ####################  Calculate p-values for sample based indels signature types ends    ################
    #########################################    Transcription    ###########################################
    ###########################################    Replication    ###########################################
    #########################################################################################################


    ######################################################
    ################# FDR BH starts ######################
    ######################################################
    all_p_values = []
    all_p_values.extend(mutationtype_transcription_pvalues)
    all_p_values.extend(mutationtype_replication_pvalues)
    all_p_values.extend(subs_signature_transcription_pvalues)
    all_p_values.extend(subs_signature_replication_pvalues)
    all_p_values.extend(indels_signature_transcription_pvalues)
    all_p_values.extend(indels_signature_replication_pvalues)

    #Order the samples
    #We must put the samples in the order they are put in all_p_values list

    # sortedSamplesList = list(sorted(samplesWithAtLeast10KMutations2NumberofMutationsDict.keys()))
    sortedSampleMutationTypesTranscriptionList = list(sorted(sample2MutationTypesTranscriptionPValuesListDict.keys()))
    sortedSampleMutationTypesReplicationList = list(sorted(sample2MutationTypesReplicationPValuesListDict.keys()))

    sortedSampleSubsSignaturesTranscriptionList = list(sorted(sample2SubsSignature2TranscriptionPValuesListDict.keys()))
    sortedSampleSubsSignaturesReplicationList = list(sorted(sample2SubsSignature2ReplicationPValuesListDict.keys()))

    sortedSampleIndelsSignaturesTranscriptionList = list(sorted(sample2IndelsSignature2TranscriptionPValuesListDict.keys()))
    sortedSampleIndelsSignaturesReplicationList = list(sorted(sample2IndelsSignature2ReplicationPValuesListDict.keys()))


    for sample in sortedSampleMutationTypesTranscriptionList:
        all_p_values.extend(sample2MutationTypesTranscriptionPValuesListDict[sample])

    for sample in sortedSampleMutationTypesReplicationList:
        all_p_values.extend(sample2MutationTypesReplicationPValuesListDict[sample])

    for sample in sortedSampleSubsSignaturesTranscriptionList:
        all_p_values.extend(sample2SubsSignature2TranscriptionPValuesListDict[sample])

    for sample in sortedSampleSubsSignaturesReplicationList:
        all_p_values.extend(sample2SubsSignature2ReplicationPValuesListDict[sample])

    for sample in sortedSampleIndelsSignaturesTranscriptionList:
        all_p_values.extend(sample2IndelsSignature2TranscriptionPValuesListDict[sample])

    for sample in sortedSampleIndelsSignaturesReplicationList:
        all_p_values.extend(sample2IndelsSignature2ReplicationPValuesListDict[sample])

    print('For information: Number of p values to be used in multiple testing: %d' %(len(all_p_values)))
    ####################################################################################
    mutationtype_transcription_FDR_BH_adjusted_pvalues = []
    mutationtype_replication_FDR_BH_adjusted_pvalues = []
    subsSignaturetype_transcription_FDR_BH_adjusted_pvalues = []
    subsSignaturetype_replication_FDR_BH_adjusted_pvalues = []
    indelsSignaturetype_transcription_FDR_BH_adjusted_pvalues = []
    indelsSignaturetype_replication_FDR_BH_adjusted_pvalues = []

    if ((all_p_values is  not None) and all_p_values):
        all_p_values_array = np.asarray(all_p_values)

        # rejected, all_FDR_BH_adjusted_p_values = statsmodels.stats.multitest.fdrcorrection(all_p_values_array, alpha=0.05, method='indep', is_sorted=False)
        rejected, all_FDR_BH_adjusted_p_values, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(all_p_values_array, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

        start = 0
        end = len(existingTranscriptionMutationTypesList)
        mutationtype_transcription_FDR_BH_adjusted_pvalues = all_FDR_BH_adjusted_p_values[start:end]

        start = end
        end += len(existingReplicationMutationTypesList)
        mutationtype_replication_FDR_BH_adjusted_pvalues = all_FDR_BH_adjusted_p_values[start:end]

        start = end
        end += len(existingTranscriptionSubsSignatureTypesList)
        subsSignaturetype_transcription_FDR_BH_adjusted_pvalues = all_FDR_BH_adjusted_p_values[start:end]

        start = end
        end += len(existingReplicationSubsSignatureTypesList)
        subsSignaturetype_replication_FDR_BH_adjusted_pvalues = all_FDR_BH_adjusted_p_values[start:end]

        start = end
        end += len(existingTranscriptionIndelsSignatureTypesList)
        indelsSignaturetype_transcription_FDR_BH_adjusted_pvalues = all_FDR_BH_adjusted_p_values[start:end]

        start = end
        end += len(existingReplicationIndelsSignatureTypesList)
        indelsSignaturetype_replication_FDR_BH_adjusted_pvalues = all_FDR_BH_adjusted_p_values[start:end]

        #We must get the BH FDR adjusted p values for each sample as in the order they are put in the all_p_values
        #Fill these dictionaries
        sample2MutationTypesTranscription_BH_FDR_Adjusted_PValues_Dict = {}
        sample2MutationTypesReplication_BH_FDR_Adjusted_PValues_Dict = {}
        sample2SubsSignatureTypesTranscription_BH_FDR_Adjusted_PValues_Dict = {}
        sample2SubsSignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict = {}
        sample2IndelsSignatureTypesTranscription_BH_FDR_Adjusted_PValues_Dict = {}
        sample2IndelsSignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict = {}

        # Sample MutationType Transcription
        for sample in sortedSampleMutationTypesTranscriptionList:
            start = end
            end += len(existingTranscriptionMutationTypesList)
            sample2MutationTypesTranscription_BH_FDR_Adjusted_PValues_Dict[sample] = all_FDR_BH_adjusted_p_values[start:end]

        # Sample MutationType Replication
        for sample in sortedSampleMutationTypesReplicationList:
            start = end
            end += len(existingReplicationMutationTypesList)
            sample2MutationTypesReplication_BH_FDR_Adjusted_PValues_Dict[sample] = all_FDR_BH_adjusted_p_values[start:end]

        # Sample Subs SignatureType Transcription
        for sample in sortedSampleSubsSignaturesTranscriptionList:
            start = end
            end += len(existingTranscriptionSubsSignatureTypesList)
            sample2SubsSignatureTypesTranscription_BH_FDR_Adjusted_PValues_Dict[sample] = all_FDR_BH_adjusted_p_values[start:end]

        # Sample Subs SignatureType Replication
        for sample in sortedSampleSubsSignaturesReplicationList:
            start = end
            end += len(existingReplicationSubsSignatureTypesList)
            sample2SubsSignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict[sample] = all_FDR_BH_adjusted_p_values[start:end]

        # Sample Indels SignatureType Transcription
        for sample in sortedSampleIndelsSignaturesTranscriptionList:
            start = end
            end += len(existingTranscriptionIndelsSignatureTypesList)
            sample2IndelsSignatureTypesTranscription_BH_FDR_Adjusted_PValues_Dict[sample] = all_FDR_BH_adjusted_p_values[start:end]

        # Sample Indels SignatureType Replication
        for sample in sortedSampleIndelsSignaturesReplicationList:
            start = end
            end += len(existingReplicationIndelsSignatureTypesList)
            sample2IndelsSignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict[sample] = all_FDR_BH_adjusted_p_values[start:end]
    ####################################################################################

    ######################################################
    ################# FDR BH ends ########################
    ######################################################


    #######################################################
    ################# Plot mutations starts ###############
    #######################################################
    # Transcription
    x_axis_labels = existingTranscriptionMutationTypesList
    N = len(x_axis_labels)

    width = 0.20
    plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,None,None,N, x_axis_labels, mutationtypes_transcribed_list, mutationtypes_untranscribed_list,simulations_mutationtypes_transcribed_medians_list, simulations_mutationtypes_untranscribed_medians_list,mutationtype_transcription_FDR_BH_adjusted_pvalues, 'Transcribed','Untranscribed', 'All Mutations', 'royalblue', 'yellowgreen','mutationtypes_transcription_strand_bias', width)

    # Replication
    x_axis_labels = existingReplicationMutationTypesList
    N = len(x_axis_labels)

    plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,None,None,N, x_axis_labels, mutationtypes_lagging_list, mutationtypes_leading_list,simulations_mutationtypes_lagging_medians_list, simulations_mutationtypes_leading_medians_list,mutationtype_replication_FDR_BH_adjusted_pvalues, 'Lagging', 'Leading','All Mutations', 'indianred', 'goldenrod', 'mutationtypes_replication_strand_bias',width)
    #######################################################
    ################# Plot mutations ends #################
    #######################################################

    #######################################################
    ########### Plot subs signatures starts ###############
    #######################################################
    x_axis_transcription_labels = existingTranscriptionSubsSignatureTypesList
    N_signatures_transcription = len(x_axis_transcription_labels)

    x_axis_replication_labels = existingReplicationSubsSignatureTypesList
    N_signatures_replication = len(x_axis_replication_labels)

    # Transcription
    plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,None,None,N_signatures_transcription,x_axis_transcription_labels, subs_signatures_transcribedStrandCount_list,subs_signatures_untranscribedStrandCount_list,simulations_subs_signatures_transcribed_medians_list, simulations_subs_signatures_untranscribed_medians_list,subsSignaturetype_transcription_FDR_BH_adjusted_pvalues, 'Transcribed','Untranscribed', 'All Signatures', 'royalblue', 'yellowgreen','subs_signatures_transcription_strand_bias', width)

    # Replication
    plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,None,None,N_signatures_replication,x_axis_replication_labels,subs_signatures_laggingStrandCount_list,subs_signatures_leadingStrandCount_list, simulations_subs_signatures_lagging_medians_list, simulations_subs_signatures_leading_medians_list,subsSignaturetype_replication_FDR_BH_adjusted_pvalues, 'Lagging','Leading', 'All Signatures', 'indianred', 'goldenrod','subs_signatures_replication_strand_bias', width)
    #######################################################
    ########### Plot subs signatures ends #################
    #######################################################

    #######################################################
    ########### Plot indels signatures starts #############
    #######################################################
    x_axis_transcription_labels = existingTranscriptionIndelsSignatureTypesList
    N_signatures_transcription = len(x_axis_transcription_labels)

    x_axis_replication_labels = existingReplicationIndelsSignatureTypesList
    N_signatures_replication = len(x_axis_replication_labels)

    # Transcription
    plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,None,None,N_signatures_transcription,x_axis_transcription_labels, indels_signatures_transcribedStrandCount_list,indels_signatures_untranscribedStrandCount_list,simulations_indels_signatures_transcribed_medians_list, simulations_indels_signatures_untranscribed_medians_list,indelsSignaturetype_transcription_FDR_BH_adjusted_pvalues, 'Transcribed','Untranscribed', 'All Signatures', 'royalblue', 'yellowgreen','indels_signatures_transcription_strand_bias', width)

    # Replication
    plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,None,None,N_signatures_replication,x_axis_replication_labels,indels_signatures_laggingStrandCount_list,indels_signatures_leadingStrandCount_list, simulations_indels_signatures_lagging_medians_list, simulations_indels_signatures_leading_medians_list,indelsSignaturetype_replication_FDR_BH_adjusted_pvalues, 'Lagging','Leading', 'All Signatures', 'indianred', 'goldenrod','indels_signatures_replication_strand_bias', width)
    #######################################################
    ########### Plot indels signatures ends ###############
    #######################################################

    #######################################################
    ######### Plot sample based mutations starts ##########
    #######################################################
    plotBarPlots(outputDir,jobname,numberofSimulations,sample2NumberofSubsDict,existingTranscriptionMutationTypesList,
                 sample2SimulationsMutationTypesTranscribedMediansListDict,
                 sample2SimulationsMutationTypesUntranscribedMediansListDict,
                 sample2MutationType2TranscriptionStrand2CountDict,
                 sample2MutationTypesTranscribedCountListDict,
                 sample2MutationTypesUntranscribedCountListDict,
                 sample2MutationTypesTranscription_BH_FDR_Adjusted_PValues_Dict,
                 width,transcriptionStrands,'royalblue', 'yellowgreen', 'All Mutations','mutationtypes_transcription_strand_bias')

    plotBarPlots(outputDir,jobname,numberofSimulations,sample2NumberofSubsDict,existingReplicationMutationTypesList,
                 sample2SimulationsMutationTypesLaggingMediansListDict,
                 sample2SimulationsMutationTypesLeadingMediansListDict,
                 sample2MutationType2ReplicationStrand2CountDict,
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
    plotBarPlots(outputDir,jobname,numberofSimulations,sample2NumberofSubsDict,existingTranscriptionSubsSignatureTypesList,
                 sample2SimulationsSubsSignatureTypesTranscribedMediansListDict,
                 sample2SimulationsSubsSignatureTypesUntranscribedMediansListDict,
                 sample2SubsSignature2TranscriptionStrand2CountDict,
                 sample2SubsSignatureTypesTranscribedCountListDict,
                 sample2SubsSignatureTypesUntranscribedCountListDict,
                 sample2SubsSignatureTypesTranscription_BH_FDR_Adjusted_PValues_Dict,
                 width,transcriptionStrands,'royalblue', 'yellowgreen','All Signatures','subs_signatures_transcription_strand_bias')

    plotBarPlots(outputDir,jobname,numberofSimulations,sample2NumberofSubsDict,existingReplicationSubsSignatureTypesList,
                 sample2SimulationsSubsSignatureTypesLaggingMediansListDict,
                 sample2SimulationsSubsSignatureTypesLeadingMediansListDict,
                 sample2SubsSignature2ReplicationStrand2CountDict,
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
    plotBarPlots(outputDir,jobname,numberofSimulations,sample2NumberofIndelsDict,existingTranscriptionIndelsSignatureTypesList,
                 sample2SimulationsIndelsSignatureTypesTranscribedMediansListDict,
                 sample2SimulationsIndelsSignatureTypesUntranscribedMediansListDict,
                 sample2IndelsSignature2TranscriptionStrand2CountDict,
                 sample2IndelsSignatureTypesTranscribedCountListDict,
                 sample2IndelsSignatureTypesUntranscribedCountListDict,
                 sample2IndelsSignatureTypesTranscription_BH_FDR_Adjusted_PValues_Dict,
                 width,transcriptionStrands,'royalblue', 'yellowgreen','All Signatures', 'indels_signatures_transcription_strand_bias')

    plotBarPlots(outputDir,jobname,numberofSimulations,sample2NumberofIndelsDict,existingReplicationIndelsSignatureTypesList,
                 sample2SimulationsIndelsSignatureTypesLaggingMediansListDict,
                 sample2SimulationsIndelsSignatureTypesLeadingMediansListDict,
                 sample2IndelsSignature2ReplicationStrand2CountDict,
                 sample2IndelsSignatureTypesLaggingCountListDict,
                 sample2IndelsSignatureTypesLeadingCountListDict,
                 sample2IndelsSignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict,
                 width,replicationStrands,'indianred', 'goldenrod','All Signatures', 'indels_signatures_replication_strand_bias')
    ############################################################
    #### Plot sample based indels signatures bar plot ends #####
    ############################################################
    ########################################################## New Part ends ####################################################################33



    ########################################################################
    ######## Bar plot starts includes sample based bar plots ###############
    ##########################  Part 4 ends ################################
    ########################################################################