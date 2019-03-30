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
        figureName = 'allmutationtypes_%s_scatterplots.png' %(STRANDBIAS)
        figureFile = os.path.join(outputDir, jobname, FIGURE, ALL, STRANDBIAS, figureName)

    else:
        figureName = 'allmutationtypes_%s_%d_%s_scatterplots.png' %(sample,numberofMutations,STRANDBIAS)
        os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, STRANDBIAS), exist_ok=True)
        figureFile = os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, STRANDBIAS, figureName)

    fig.savefig(figureFile)
    plt.close(fig)
########################################################################

########################################################################
#May 9, 2018 starts
#For Signatures
def plot_ncomms11383_Supp_FigH_AllSignatures_TranscriptionLog10Ratio_ReplicationLog10Ratio(
        sample,
        numberofMutations,
        signature2TranscriptionStrand2CountDict,
        signature2ReplicationStrand2CountDict,
        signatures,
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

    for signature in signatures:
        #################################################################################################
        #First check whether we have this signature or not
        if ((signature in signature2TranscriptionStrand2CountDict) and (TRANSCRIBED_STRAND in (signature2TranscriptionStrand2CountDict[signature])) and
             (UNTRANSCRIBED_STRAND in (signature2TranscriptionStrand2CountDict[signature])) ):

            if ((signature2TranscriptionStrand2CountDict[signature][TRANSCRIBED_STRAND]+signature2TranscriptionStrand2CountDict[signature][UNTRANSCRIBED_STRAND]) >= ONE_THOUSAND):
                transcriptionRatiosDict[signature]= np.log10(signature2TranscriptionStrand2CountDict[signature][TRANSCRIBED_STRAND]/signature2TranscriptionStrand2CountDict[signature][UNTRANSCRIBED_STRAND])
        #################################################################################################

        #################################################################################################
        # First check whether we have this signature or not
        if ((signature in signature2ReplicationStrand2CountDict) and (LAGGING in (signature2ReplicationStrand2CountDict[signature])) and
                (LEADING in (signature2ReplicationStrand2CountDict[signature]))):

            if ((signature2ReplicationStrand2CountDict[signature][LAGGING]+signature2ReplicationStrand2CountDict[signature][LEADING])>= ONE_THOUSAND):
                replicationRatiosDict[signature] = np.log10(signature2ReplicationStrand2CountDict[signature][LAGGING]/signature2ReplicationStrand2CountDict[signature][LEADING])
        #################################################################################################

    signaturesShownInLegend= []

    for signature in signatures:
        if ((signature in replicationRatiosDict.keys()) and (signature in transcriptionRatiosDict.keys())):
            signaturesShownInLegend.append(signature)
            plt.scatter(replicationRatiosDict[signature],transcriptionRatiosDict[signature],label = signature)

    legend = plt.legend(loc='upper left',frameon=True,fancybox= False,labels=signaturesShownInLegend,bbox_to_anchor=(-0.0095, 1.0095))
    legend.get_frame().set_linewidth(1)

    frame = legend.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('black')

    plt.axvline(x=0.0, color='gray', linestyle='--')
    plt.axhline(y=0.0, color='gray', linestyle='--')

    if sample is None:
        figureName = 'allsignatures_%s_scatterplots.png' %(STRANDBIAS)
        figureFile = os.path.join(outputDir, jobname, FIGURE, ALL, STRANDBIAS, figureName)

    else:
        figureName = 'allsignatures_%s_%d_%s_scatterplots.png' %(sample,numberofMutations,STRANDBIAS)
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

        figureName = newMutationType + '_' + STRANDBIAS + '_MutationTypeBased_SampleBased.png'
        figureFile = os.path.join(outputDir,jobname,FIGURE,ALL,STRANDBIAS,figureName)
        fig.savefig(figureFile)
        plt.close(fig)
########################################################################


########################################################################
#SignatureBased SampleBased Figures
#Sig26 is very different
def plot_ncomms11383_Supp_FigF_SignatureBased_AllSamples_TranscriptionLog10Ratio_ReplicationLog10Ratio(mutationProbability2Signature2Sample2TranscriptionStrand2CountDict,mutationProbability2Signature2Sample2ReplicationStrand2CountDict,signatures,outputDir,jobname,isFigureAugmentation):

    #TODO Check it. Dow we need string instead of float?
    signature2Sample2TranscriptionStrand2CountDict= mutationProbability2Signature2Sample2TranscriptionStrand2CountDict[SUBSTITUTION_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD]
    signature2Sample2ReplicationStrand2CountDict = mutationProbability2Signature2Sample2ReplicationStrand2CountDict[SUBSTITUTION_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD]

    # old code
    # signature2Sample2ReplicationStrand2CountDict = mutationProbability2Signature2Sample2ReplicationStrand2CountDict[MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING]

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

            figureName = signature.replace(' ','') + '_' + STRANDBIAS + '_SignatureBased_SampleBased.png'
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

    #AttributeError: 'function' object has no attribute 'get_frame'
    # frame = ax.legend.get_frame()
    # frame.set_facecolor('white')
    # frame.set_edgecolor('black')

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
                simulationType2StrandCountDict = simulationType2StrandCountDict[MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING]

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
                    simulationBased_type2Sample2Strand2CountDict = simulationBased_type2Sample2Strand2CountDict[MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING]

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

# ##################################################################
def transcriptionReplicationStrandBiasFigures(outputDir,jobname,figureAugmentation,numberofSimulations):

    isFigureAugmentation = False
    if (figureAugmentation == 'augmentation'):
        isFigureAugmentation = True

    #######################################################################################################################
    # new code
    os.makedirs(os.path.join(outputDir, jobname, FIGURE, ALL, STRANDBIAS), exist_ok=True)
    os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES), exist_ok=True)
    #######################################################################################################################


    ##########################################################################################
    #########################  Read dictionaries related with ################################
    #########################  signatures and samples starts  ################################
    ##########################################################################################
    # Get signatures
    signaturesFilePath = os.path.join(outputDir, jobname, DATA, SignatureFilename)
    signatures = readSignatureList(signaturesFilePath)

    ##########################################################################################
    #Please note that this data structure is being used
    signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = {}
    signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath = os.path.join(outputDir,jobname, DATA,
                                                                                           SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilename)

    if (os.path.exists(signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath)):
        signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = readDictionary(signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath)
    ##########################################################################################


    ##########################################################################################
    #Please note that this data structure is being used
    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = {}
    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath = os.path.join(outputDir, jobname, DATA,Sample2SubsSignature2NumberofMutationsDictFilename)

    if (os.path.exists(sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath)):
        sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = readDictionary(sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath)
    ##########################################################################################



    ##########################################################################################
    #Load samplesWithAtLeast10KMutations2NumberofMutationsDict
    samplesWithAtLeast10KMutations2NumberofMutationsDict = {}
    samplesWithAtLeast10KMutations2NumberofMutationsDictFilePath = os.path.join(outputDir,jobname,DATA,Sample2NumberofSubsDictFilename)

    if (os.path.exists(samplesWithAtLeast10KMutations2NumberofMutationsDictFilePath)):
        samplesWithAtLeast10KMutations2NumberofMutationsDict = readDictionary(samplesWithAtLeast10KMutations2NumberofMutationsDictFilePath)
    ##########################################################################################

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
    mutationProbability2Signature2TranscriptionStrand2CountFilePath = os.path.join(outputDir,jobname,DATA,TRANSCRIPTIONSTRANDBIAS,MutationProbability2Signature2TranscriptionStrand2CountDict_Filename)
    mutationProbability2Signature2TranscriptionStrand2CountDict = readDictionary(mutationProbability2Signature2TranscriptionStrand2CountFilePath)

    mutationProbability2Signature2ReplicationStrand2CountFilePath = os.path.join(outputDir, jobname,DATA, REPLICATIONSTRANDBIAS,MutationProbability2Signature2ReplicationStrand2CountDict_Filename)
    mutationProbability2Signature2ReplicationStrand2CountDict = readDictionary(mutationProbability2Signature2ReplicationStrand2CountFilePath)

    mutationProbability2Signature2Sample2TranscriptionStrand2CountFilePath = os.path.join(outputDir,jobname, DATA,TRANSCRIPTIONSTRANDBIAS,MutationProbability2Signature2Sample2TranscriptionStrand2CountDict_Filename)
    mutationProbability2Signature2Sample2TranscriptionStrand2CountDict = readDictionary(mutationProbability2Signature2Sample2TranscriptionStrand2CountFilePath)

    mutationProbability2Signature2Sample2ReplicationStrand2CountFilePath = os.path.join(outputDir,jobname,DATA,REPLICATIONSTRANDBIAS,MutationProbability2Signature2Sample2ReplicationStrand2CountDict_Filename)
    mutationProbability2Signature2Sample2ReplicationStrand2CountDict = readDictionary(mutationProbability2Signature2Sample2ReplicationStrand2CountFilePath)
    ##########################################################################################

    #Initialize as empty lists
    signature2Sample2TranscriptionStrand2CountDict = {}
    signature2Sample2ReplicationStrand2CountDict = {}
    signature2TranscriptionStrand2CountDict = {}
    signature2ReplicationStrand2CountDict = {}

    #get signature based for mutation probability = 0.5
    if ((mutationProbability2Signature2Sample2TranscriptionStrand2CountDict is not None) and (mutationProbability2Signature2Sample2TranscriptionStrand2CountDict)):
        signature2Sample2TranscriptionStrand2CountDict = mutationProbability2Signature2Sample2TranscriptionStrand2CountDict[MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING]

    if (mutationProbability2Signature2Sample2ReplicationStrand2CountDict is not None and (mutationProbability2Signature2Sample2ReplicationStrand2CountDict)):
        signature2Sample2ReplicationStrand2CountDict = mutationProbability2Signature2Sample2ReplicationStrand2CountDict[MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING]

    if ((mutationProbability2Signature2TranscriptionStrand2CountDict is not None) and (mutationProbability2Signature2TranscriptionStrand2CountDict)):
        signature2TranscriptionStrand2CountDict = mutationProbability2Signature2TranscriptionStrand2CountDict[MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING]

    if ((mutationProbability2Signature2ReplicationStrand2CountDict is not None) and (mutationProbability2Signature2ReplicationStrand2CountDict)):
        signature2ReplicationStrand2CountDict = mutationProbability2Signature2ReplicationStrand2CountDict[MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING]

    existingTranscriptionMutationTypesList = sorted(list(mutationType2TranscriptionStrandCountDict.keys()))
    existingReplicationMutationTypesList = sorted(list(mutationType2ReplicationStrandCountDict.keys()))
    existingTranscriptionSignatureTypesList = sorted(list(signature2TranscriptionStrand2CountDict.keys()))
    existingReplicationSignatureTypesList = sorted(list(signature2ReplicationStrand2CountDict.keys()))

    ##########################################################################################
    #####################  Read dictionaries for the rest part ends ##########################
    ##########################################################################################

    ########################################################################
    ##########################  Part 1 starts ##############################
    ################## MutationType All Samples Figures starts #############
    ############### SignatureType All Samples Figures starts ###############
    ########################################################################
    # Note: Since these are already in sample based no extra work is required.
    #MutationType SampleBased Figures --- No numberofMutations
    if ((mutationType2Sample2TranscriptionStrand2CountDict is not None) and (mutationType2Sample2ReplicationStrand2CountDict is not None)):
        plot_ncomms11383_Supp_FigE_MutationTypeBased_AllSamples_TranscriptionLog10Ratio_ReplicationLog10Ratio(mutationType2Sample2TranscriptionStrand2CountDict,mutationType2Sample2ReplicationStrand2CountDict,outputDir,jobname,isFigureAugmentation)

    #SignatureType Sample Based Figures --- No numberofMutations
    if ((mutationProbability2Signature2Sample2TranscriptionStrand2CountDict is not None) and (mutationProbability2Signature2Sample2ReplicationStrand2CountDict is not None)):
        plot_ncomms11383_Supp_FigF_SignatureBased_AllSamples_TranscriptionLog10Ratio_ReplicationLog10Ratio(mutationProbability2Signature2Sample2TranscriptionStrand2CountDict,mutationProbability2Signature2Sample2ReplicationStrand2CountDict,signatures,outputDir,jobname,isFigureAugmentation)
    ########################################################################
    ################## MutationType All Samples Figures ends ###############
    ############### SignatureType All Samples Figures ends #################
    ##########################  Part 1 ends ################################
    ########################################################################


    ########################################################################
    ##########################  Part 2 starts ##############################
    ############## MutationType Scatter Plots starts #######################
    ############## Signature Based Scatter Plots starts ####################
    ########################################################################
    if ((mutationType2TranscriptionStrandCountDict is not None) and (mutationType2ReplicationStrandCountDict is not None)):
        plot_ncomms11383_Supp_FigG_AllMutationTypes_TranscriptionLog10Ratio_ReplicationLog10Ratio(None,None,mutationType2TranscriptionStrandCountDict, mutationType2ReplicationStrandCountDict,outputDir,jobname,isFigureAugmentation)

    if ((mutationProbability2Signature2TranscriptionStrand2CountDict is not None) and (mutationProbability2Signature2ReplicationStrand2CountDict is not None) and (signatures)):
        signature2TranscriptionStrand2CountDict = mutationProbability2Signature2TranscriptionStrand2CountDict[MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING]
        signature2ReplicationStrand2CountDict = mutationProbability2Signature2ReplicationStrand2CountDict[MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING]
        plot_ncomms11383_Supp_FigH_AllSignatures_TranscriptionLog10Ratio_ReplicationLog10Ratio(None,None,signature2TranscriptionStrand2CountDict,signature2ReplicationStrand2CountDict,signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,outputDir,jobname,isFigureAugmentation)
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

    for sample in samplesWithAtLeast10KMutations2NumberofMutationsDict:
        if ((sample in sample2MutationType2TranscriptionStrand2CountDict) and (sample in sample2MutationType2ReplicationStrand2CountDict)):
            numberofMutations = samplesWithAtLeast10KMutations2NumberofMutationsDict[sample]
            plot_ncomms11383_Supp_FigG_AllMutationTypes_TranscriptionLog10Ratio_ReplicationLog10Ratio(sample,numberofMutations,sample2MutationType2TranscriptionStrand2CountDict[sample], sample2MutationType2ReplicationStrand2CountDict[sample],outputDir,jobname,isFigureAugmentation)
    ###############################################################

    ###############################################################
    sample2Signature2TranscriptionStrand2CountDict = {}
    sample2Signature2ReplicationStrand2CountDict = {}

    convert(signature2Sample2TranscriptionStrand2CountDict, sample2Signature2TranscriptionStrand2CountDict)
    convert(signature2Sample2ReplicationStrand2CountDict, sample2Signature2ReplicationStrand2CountDict)

    for sample in samplesWithAtLeast10KMutations2NumberofMutationsDict:
        if ((sample in sample2Signature2TranscriptionStrand2CountDict) and (sample in sample2Signature2ReplicationStrand2CountDict)):
            numberofMutations = samplesWithAtLeast10KMutations2NumberofMutationsDict[sample]
            plot_ncomms11383_Supp_FigH_AllSignatures_TranscriptionLog10Ratio_ReplicationLog10Ratio(sample,numberofMutations,sample2Signature2TranscriptionStrand2CountDict[sample], sample2Signature2ReplicationStrand2CountDict[sample], signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict, outputDir,jobname,isFigureAugmentation)
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
    simulations_mutationtypes_transcribed_medians_list = None
    simulations_mutationtypes_untranscribed_medians_list = None

    samplebased_simulations_mutationtypes_transcribed_medians_list = None
    samplebased_simulations_mutationtypes_untranscribed_medians_list = None

    #######################################################

    # This is required for both mutationtypes transcription and
    # This is required for both simulations mutationtypes transcription

    ##########################################################################################################
    ####################### Fill simulations_mutationtypes_transcribed_medians_list starts ###################
    ####################### Fill simulations_mutationtypes_untranscribed_medians_list starts #################
    ########################### uses existingTranscriptionMutationTypesList ##################################
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
    simulations_signatures_transcribed_medians_list = None
    simulations_signatures_untranscribed_medians_list = None
    simulations_signatures_lagging_medians_list = None
    simulations_signatures_leading_medians_list = None

    sample2SimulationsSignatureTypesTranscribedMediansListDict = None
    sample2SimulationsSignatureTypesUntranscribedMediansListDict = None
    sample2SimulationsSignatureTypesLaggingMediansListDict = None
    sample2SimulationsSignatureTypesLeadingMediansListDict = None
    #######

    if (numberofSimulations > 0):
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
        # simulations --- signatures --- transcription
        simulations_signatures_transcribed_medians_list, simulations_signatures_untranscribed_medians_list =\
            fillSimulationsType2StrandCountList(outputDir,
                                                jobname,
                                                existingTranscriptionSignatureTypesList,
                                                numberofSimulations,
                                                TRANSCRIPTIONSTRANDBIAS,
                                                MutationProbability2Signature2TranscriptionStrand2CountDict_Filename,
                                                SIGNATURE)

        # simulations --- signatures --- replication
        simulations_signatures_lagging_medians_list, simulations_signatures_leading_medians_list = \
            fillSimulationsType2StrandCountList(outputDir,
                                                jobname,
                                                existingReplicationSignatureTypesList,
                                                numberofSimulations,
                                                REPLICATIONSTRANDBIAS,
                                                MutationProbability2Signature2ReplicationStrand2CountDict_Filename,
                                                SIGNATURE)


        #sample based starts
        # samplebased --- simulations --- mutation type --- transcription
        sample2SimulationsMutationTypesTranscribedMediansListDict, sample2SimulationsMutationTypesUntranscribedMediansListDict =\
            fillSimulationsSample2Type2StrandCountList(outputDir,jobname,
                existingTranscriptionMutationTypesList,
                numberofSimulations,
                samplesWithAtLeast10KMutations2NumberofMutationsDict,
                TRANSCRIPTIONSTRANDBIAS,
                MutationType2Sample2TranscriptionStrand2CountDict_Filename,
                MUTATION)

        # samplebased --- simulations --- mutation type --- replication
        sample2SimulationsMutationTypesLaggingMediansListDict, sample2SimulationsMutationTypesLeadingMediansListDict =\
            fillSimulationsSample2Type2StrandCountList(outputDir,jobname,
                existingReplicationMutationTypesList,
                numberofSimulations,
                samplesWithAtLeast10KMutations2NumberofMutationsDict,
                REPLICATIONSTRANDBIAS,
                MutationType2Sample2ReplicationStrand2CountDict_Filename,
                MUTATION)

        # samplebased --- simulations --- signature --- transcription
        sample2SimulationsSignatureTypesTranscribedMediansListDict, sample2SimulationsSignatureTypesUntranscribedMediansListDict =\
            fillSimulationsSample2Type2StrandCountList(outputDir,jobname,
                existingTranscriptionSignatureTypesList,
                numberofSimulations,
                samplesWithAtLeast10KMutations2NumberofMutationsDict,
                TRANSCRIPTIONSTRANDBIAS,
                MutationProbability2Signature2Sample2TranscriptionStrand2CountDict_Filename,
                SIGNATURE)

        # samplebased --- simulations --- signature --- replication
        sample2SimulationsSignatureTypesLaggingMediansListDict, sample2SimulationsSignatureTypesLeadingMediansListDict =\
            fillSimulationsSample2Type2StrandCountList(outputDir,jobname,
                existingReplicationSignatureTypesList,
                numberofSimulations,
                samplesWithAtLeast10KMutations2NumberofMutationsDict,
                REPLICATIONSTRANDBIAS,
                MutationProbability2Signature2Sample2ReplicationStrand2CountDict_Filename,
                SIGNATURE)
        #sample based ends
    ##########################################################################################################
    ####################### Fill simulations_mutationtypes_transcribed_medians_list ends #####################
    ####################### Fill simulations_mutationtypes_untranscribed_medians_list ends ###################
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
    #############################  Calculate p-values for signatures starts    ##############################
    #########################################################################################################

    ########################################## Transcription starts  ########################################
    strandType_Signature2CountList_Dict = fillDictWRTReferenceTypes(existingTranscriptionSignatureTypesList,transcriptionStrands,signature2TranscriptionStrand2CountDict)
    signatures_transcribedStrandCount_list = strandType_Signature2CountList_Dict[TRANSCRIBED_STRAND]
    signatures_untranscribedStrandCount_list = strandType_Signature2CountList_Dict[UNTRANSCRIBED_STRAND]


    signature_transcription_pvalues = calculate_p_values(signatures_transcribedStrandCount_list,
                                                         simulations_signatures_transcribed_medians_list,
                                                         signatures_untranscribedStrandCount_list,
                                                         simulations_signatures_untranscribed_medians_list)
    ########################################## Transcription ends  ##########################################

    ########################################## Replication starts  ##########################################
    strandType_Signature2CountList_Dict = fillDictWRTReferenceTypes(existingReplicationSignatureTypesList,replicationStrands,signature2ReplicationStrand2CountDict)
    signatures_laggingStrandCount_list = strandType_Signature2CountList_Dict[LAGGING]
    signatures_leadingStrandCount_list = strandType_Signature2CountList_Dict[LEADING]

    signature_replication_pvalues = calculate_p_values(signatures_laggingStrandCount_list,
                                                       simulations_signatures_lagging_medians_list,
                                                       signatures_leadingStrandCount_list,
                                                       simulations_signatures_leading_medians_list)
    ########################################## Replication ends  ############################################

    #########################################################################################################
    #############################  Calculate p-values for signatures ends    ################################
    #########################################################################################################


    #########################################################################################################
    ####################  Calculate p-values for sample based mutation types starts    ######################
    #########################################    Transcription    ###########################################
    ###########################################    Replication    ###########################################
    #########################################################################################################
    #Fill these dictionaries required for the bar plots
    sample2MutationTypesTranscribedCountListDict = {}
    sample2MutationTypesUntranscribedCountListDict = {}

    sample2MutationTypesLaggingCountListDict = {}
    sample2MutationTypesLeadingCountListDict = {}

    #Fill these dictionaries
    sample2MutationTypesTranscriptionPValuesListDict= {}
    sample2MutationTypesReplicationPValuesListDict= {}

    ##################################################################################
    for sample in sample2MutationType2TranscriptionStrand2CountDict:
        strandType_MutationTypeTranscriptionCountList_Dict = fillDictWRTReferenceTypes(existingTranscriptionMutationTypesList,transcriptionStrands,sample2MutationType2TranscriptionStrand2CountDict[sample])

        sampleBased_mutationtypes_transcribed_list = strandType_MutationTypeTranscriptionCountList_Dict[TRANSCRIBED_STRAND]
        sampleBased_mutationtypes_untranscribed_list = strandType_MutationTypeTranscriptionCountList_Dict[UNTRANSCRIBED_STRAND]

        sample2MutationTypesTranscribedCountListDict[sample] = sampleBased_mutationtypes_transcribed_list
        sample2MutationTypesUntranscribedCountListDict[sample] = sampleBased_mutationtypes_untranscribed_list

        simulations_sampleBased_mutationtypes_transcribed_list = None
        simulations_sampleBased_mutationtypes_untranscribed_list = None
        #Now we have
        if ((sample2SimulationsMutationTypesTranscribedMediansListDict is not None) and (sample in sample2SimulationsMutationTypesTranscribedMediansListDict)):
            simulations_sampleBased_mutationtypes_transcribed_list = sample2SimulationsMutationTypesTranscribedMediansListDict[sample]
        if  ((sample2SimulationsMutationTypesUntranscribedMediansListDict is not None ) and (sample in sample2SimulationsMutationTypesUntranscribedMediansListDict)):
            simulations_sampleBased_mutationtypes_untranscribed_list = sample2SimulationsMutationTypesUntranscribedMediansListDict[sample]


        if ((simulations_sampleBased_mutationtypes_transcribed_list is not None) and (simulations_sampleBased_mutationtypes_untranscribed_list is not None)):

            # Please notice Fisher exact test does not return anything for large numbers therefore chi2 test is used.
            samplebased_mutationtype_transcription_pvalues = calculate_p_values(sampleBased_mutationtypes_transcribed_list,
                                                                    simulations_sampleBased_mutationtypes_transcribed_list,
                                                                    sampleBased_mutationtypes_untranscribed_list,
                                                                    simulations_sampleBased_mutationtypes_untranscribed_list)

            sample2MutationTypesTranscriptionPValuesListDict[sample] = samplebased_mutationtype_transcription_pvalues
    ##################################################################################

    ##################################################################################
    for sample in sample2MutationType2ReplicationStrand2CountDict:
        strandType_MutationTypeReplicationCountList_Dict = fillDictWRTReferenceTypes(existingReplicationMutationTypesList,replicationStrands,sample2MutationType2ReplicationStrand2CountDict[sample])

        sampleBased_mutationtypes_lagging_list = strandType_MutationTypeReplicationCountList_Dict[LAGGING]
        sampleBased_mutationtypes_leading_list = strandType_MutationTypeReplicationCountList_Dict[LEADING]

        sample2MutationTypesLaggingCountListDict[sample] = sampleBased_mutationtypes_lagging_list
        sample2MutationTypesLeadingCountListDict[sample] = sampleBased_mutationtypes_leading_list

        # Please notice Fisher exact test does not return anything for large numbers therefore chi2 test is used.

        simulations_sampleBased_mutationtypes_lagging_list = None
        simulations_sampleBased_mutationtypes_leading_list = None

        if ((sample2SimulationsMutationTypesLaggingMediansListDict is not None) and (sample in sample2SimulationsMutationTypesLaggingMediansListDict)):
            simulations_sampleBased_mutationtypes_lagging_list = sample2SimulationsMutationTypesLaggingMediansListDict[sample]

        if ((sample2SimulationsMutationTypesLeadingMediansListDict is not None) and (sample in sample2SimulationsMutationTypesLeadingMediansListDict)):
            simulations_sampleBased_mutationtypes_leading_list = sample2SimulationsMutationTypesLeadingMediansListDict[sample]

        if ((simulations_sampleBased_mutationtypes_lagging_list is not None) and (simulations_sampleBased_mutationtypes_leading_list is not None)):
            samplebased_mutationtype_replication_pvalues = calculate_p_values(sampleBased_mutationtypes_lagging_list,
                                                                  simulations_sampleBased_mutationtypes_lagging_list,
                                                                  sampleBased_mutationtypes_leading_list,
                                                                  simulations_sampleBased_mutationtypes_leading_list)
            sample2MutationTypesReplicationPValuesListDict[sample] = samplebased_mutationtype_replication_pvalues
    ##################################################################################

    #########################################################################################################
    ####################  Calculate p-values for sample based mutation types ends    ########################
    #########################################    Transcription    ###########################################
    ###########################################    Replication    ###########################################
    #########################################################################################################



    #########################################################################################################
    ####################  Calculate p-values for sample based signature types starts    #####################
    #########################################    Transcription    ###########################################
    ###########################################    Replication    ###########################################
    #########################################################################################################
    #Fill these dictionaries required for the bar plots
    sample2SignatureTypesTranscribedCountListDict = {}
    sample2SignatureTypesUntranscribedCountListDict = {}

    sample2SignatureTypesLaggingCountListDict = {}
    sample2SignatureTypesLeadingCountListDict = {}

    #Fill these dictionaries
    sample2Signature2TranscriptionPValuesListDict = {}
    sample2Signature2ReplicationPValuesListDict = {}

    ##################################################################################
    for sample in sample2Signature2TranscriptionStrand2CountDict:
        strandType_Signature2CountList_Dict = fillDictWRTReferenceTypes(existingTranscriptionSignatureTypesList,transcriptionStrands,sample2Signature2TranscriptionStrand2CountDict[sample])

        sampleBased_signatures_transcribedStrandCount_list = strandType_Signature2CountList_Dict[TRANSCRIBED_STRAND]
        sampleBased_signatures_untranscribedStrandCount_list = strandType_Signature2CountList_Dict[UNTRANSCRIBED_STRAND]

        sample2SignatureTypesTranscribedCountListDict[sample] = sampleBased_signatures_transcribedStrandCount_list
        sample2SignatureTypesUntranscribedCountListDict[sample] = sampleBased_signatures_untranscribedStrandCount_list

        simulations_sampleBased_signatures_transcribedStrandCount_list = None
        simulations_sampleBased_signatures_untranscribedStrandCount_list = None

        # Transcription: Transcribed/(Transcribed+Untranscribed)
        if ((sample2SimulationsSignatureTypesTranscribedMediansListDict is not None) and (sample in sample2SimulationsSignatureTypesTranscribedMediansListDict)):
            simulations_sampleBased_signatures_transcribedStrandCount_list = sample2SimulationsSignatureTypesTranscribedMediansListDict[sample]

        if ((sample2SimulationsSignatureTypesUntranscribedMediansListDict is not None) and (sample in sample2SimulationsSignatureTypesUntranscribedMediansListDict)):
            simulations_sampleBased_signatures_untranscribedStrandCount_list = sample2SimulationsSignatureTypesUntranscribedMediansListDict[sample]

        if ((simulations_sampleBased_signatures_transcribedStrandCount_list is not None) and (simulations_sampleBased_signatures_untranscribedStrandCount_list is not None)):
            samplebased_signature_transcription_pvalues = calculate_p_values(sampleBased_signatures_transcribedStrandCount_list,
                                                                 simulations_sampleBased_signatures_transcribedStrandCount_list,
                                                                 sampleBased_signatures_untranscribedStrandCount_list,
                                                                 simulations_sampleBased_signatures_untranscribedStrandCount_list)

            sample2Signature2TranscriptionPValuesListDict[sample] = samplebased_signature_transcription_pvalues
    ##################################################################################

    ##################################################################################
    for sample in sample2Signature2ReplicationStrand2CountDict:
        strandType_Signature2CountList_Dict = fillDictWRTReferenceTypes(existingReplicationSignatureTypesList,replicationStrands,sample2Signature2ReplicationStrand2CountDict[sample])

        sampleBased_signatures_laggingStrandCount_list = strandType_Signature2CountList_Dict[LAGGING]
        sampleBased_signatures_leadingStrandCount_list = strandType_Signature2CountList_Dict[LEADING]

        sample2SignatureTypesLaggingCountListDict[sample] = sampleBased_signatures_laggingStrandCount_list
        sample2SignatureTypesLeadingCountListDict[sample] = sampleBased_signatures_leadingStrandCount_list

        simulations_sampleBased_signatures_laggingStrandCount_list = None
        simulations_sampleBased_signatures_leadingStrandCount_list = None

        # Replication: Lagging/(Lagging+Leading)
        if ((sample2SimulationsSignatureTypesLaggingMediansListDict is not None) and (sample in sample2SimulationsSignatureTypesLaggingMediansListDict)):
            simulations_sampleBased_signatures_laggingStrandCount_list = sample2SimulationsSignatureTypesLaggingMediansListDict[sample]

        if ((sample2SimulationsSignatureTypesLeadingMediansListDict is not None) and (sample in sample2SimulationsSignatureTypesLeadingMediansListDict)):
            simulations_sampleBased_signatures_leadingStrandCount_list = sample2SimulationsSignatureTypesLeadingMediansListDict[sample]

        if ((simulations_sampleBased_signatures_laggingStrandCount_list is not None) and (simulations_sampleBased_signatures_leadingStrandCount_list is not None)):
            samplebased_signature_replication_pvalues = calculate_p_values(sampleBased_signatures_laggingStrandCount_list,
                                                               simulations_sampleBased_signatures_laggingStrandCount_list,
                                                               sampleBased_signatures_leadingStrandCount_list,
                                                               simulations_sampleBased_signatures_leadingStrandCount_list)

            sample2Signature2ReplicationPValuesListDict[sample] = samplebased_signature_replication_pvalues
    ##################################################################################

    #########################################################################################################
    ####################  Calculate p-values for sample based signature types ends    #######################
    #########################################    Transcription    ###########################################
    ###########################################    Replication    ###########################################
    #########################################################################################################


    ######################################################
    ################# FDR BH starts ######################
    ######################################################
    all_p_values = []
    all_p_values.extend(mutationtype_transcription_pvalues)
    all_p_values.extend(mutationtype_replication_pvalues)
    all_p_values.extend(signature_transcription_pvalues)
    all_p_values.extend(signature_replication_pvalues)

    #Order the samples
    #We must put the samples in the order they are put in all_p_values list

    # sortedSamplesList = list(sorted(samplesWithAtLeast10KMutations2NumberofMutationsDict.keys()))
    sortedSampleMutationTypesTranscriptionList = list(sorted(sample2MutationTypesTranscriptionPValuesListDict.keys()))
    sortedSampleMutationTypesReplicationList = list(sorted(sample2MutationTypesReplicationPValuesListDict.keys()))

    sortedSampleSignaturesTranscriptionList = list(sorted(sample2Signature2TranscriptionPValuesListDict.keys()))
    sortedSampleSignaturesReplicationList = list(sorted(sample2Signature2ReplicationPValuesListDict.keys()))

    for sample in sortedSampleMutationTypesTranscriptionList:
        all_p_values.extend(sample2MutationTypesTranscriptionPValuesListDict[sample])

    for sample in sortedSampleMutationTypesReplicationList:
        all_p_values.extend(sample2MutationTypesReplicationPValuesListDict[sample])

    for sample in sortedSampleSignaturesTranscriptionList:
        all_p_values.extend(sample2Signature2TranscriptionPValuesListDict[sample])

    for sample in sortedSampleSignaturesReplicationList:
        all_p_values.extend(sample2Signature2ReplicationPValuesListDict[sample])



    ####################################################################################
    mutationtype_transcription_FDR_BH_adjusted_pvalues = []
    mutationtype_replication_FDR_BH_adjusted_pvalues = []
    signaturetype_transcription_FDR_BH_adjusted_pvalues = []
    signaturetype_replication_FDR_BH_adjusted_pvalues = []

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
        end += len(existingTranscriptionSignatureTypesList)
        signaturetype_transcription_FDR_BH_adjusted_pvalues = all_FDR_BH_adjusted_p_values[start:end]

        start = end
        end += len(existingReplicationSignatureTypesList)
        signaturetype_replication_FDR_BH_adjusted_pvalues = all_FDR_BH_adjusted_p_values[start:end]

        #We must get the BH FDR adjusted p values for each sample as in the order they are put in the all_p_values
        #Fill these dictionaries
        sample2MutationTypesTranscription_BH_FDR_Adjusted_PValues_Dict = {}
        sample2MutationTypesReplication_BH_FDR_Adjusted_PValues_Dict = {}
        sample2SignatureTypesTranscription_BH_FDR_Adjusted_PValues_Dict = {}
        sample2SignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict = {}

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

        # Sample SignatureType Transcription
        for sample in sortedSampleSignaturesTranscriptionList:
            start = end
            end += len(existingTranscriptionSignatureTypesList)
            sample2SignatureTypesTranscription_BH_FDR_Adjusted_PValues_Dict[sample] = all_FDR_BH_adjusted_p_values[start:end]

        # Sample SignatureType Replication
        for sample in sortedSampleSignaturesReplicationList:
            start = end
            end += len(existingReplicationSignatureTypesList)
            sample2SignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict[sample] = all_FDR_BH_adjusted_p_values[start:end]
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

    ########################################################## New Part starts ##################################################################33
    #######################################################
    ######### Plot sample based mutations starts ##########
    #######################################################

    #######################################################
    for sample in sample2MutationType2TranscriptionStrand2CountDict:
        mutationtype_2_transcribed_list = sample2MutationTypesTranscribedCountListDict[sample]
        mutationtype_2_untranscribed_list = sample2MutationTypesUntranscribedCountListDict[sample]

        x_axis_labels = existingTranscriptionMutationTypesList
        N = len(x_axis_labels)

        if ((sample in samplesWithAtLeast10KMutations2NumberofMutationsDict) and (sample in sample2MutationTypesTranscription_BH_FDR_Adjusted_PValues_Dict)):
            sample_mutationtypes_transcription_FDR_BH_adjusted_pvalues = sample2MutationTypesTranscription_BH_FDR_Adjusted_PValues_Dict[sample]
            numberofMutations = samplesWithAtLeast10KMutations2NumberofMutationsDict[sample]
            samplebased_simulations_mutationtypes_transcribed_medians_list = sample2SimulationsMutationTypesTranscribedMediansListDict[sample]
            samplebased_simulations_mutationtypes_untranscribed_medians_list = sample2SimulationsMutationTypesUntranscribedMediansListDict[sample]
            plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,sample,numberofMutations,N,x_axis_labels,mutationtype_2_transcribed_list,mutationtype_2_untranscribed_list,samplebased_simulations_mutationtypes_transcribed_medians_list,samplebased_simulations_mutationtypes_untranscribed_medians_list,sample_mutationtypes_transcription_FDR_BH_adjusted_pvalues,'Transcribed','Untranscribed','All Mutations','royalblue','yellowgreen','mutationtypes_transcription_strand_bias',width)
    #######################################################

    #######################################################
    for sample in sample2MutationType2ReplicationStrand2CountDict:
        mutationtype_2_lagging_list = sample2MutationTypesLaggingCountListDict[sample]
        mutationtype_2_leading_list = sample2MutationTypesLeadingCountListDict[sample]

        x_axis_labels = existingReplicationMutationTypesList
        N = len(x_axis_labels)

        # Replication
        if ((sample in samplesWithAtLeast10KMutations2NumberofMutationsDict) and (sample in sample2MutationTypesReplication_BH_FDR_Adjusted_PValues_Dict)):
            sample_mutationtypes_replication_FDR_BH_adjusted_pvalues = sample2MutationTypesReplication_BH_FDR_Adjusted_PValues_Dict[sample]
            numberofMutations = samplesWithAtLeast10KMutations2NumberofMutationsDict[sample]
            samplebased_simulations_mutationtypes_lagging_medians_list = sample2SimulationsMutationTypesLaggingMediansListDict[sample]
            samplebased_simulations_mutationtypes_leading_medians_list = sample2SimulationsMutationTypesLeadingMediansListDict[sample]
            plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,sample,numberofMutations,N,x_axis_labels,mutationtype_2_lagging_list,mutationtype_2_leading_list,samplebased_simulations_mutationtypes_lagging_medians_list,samplebased_simulations_mutationtypes_leading_medians_list,sample_mutationtypes_replication_FDR_BH_adjusted_pvalues,'Lagging', 'Leading','All Mutations', 'indianred', 'goldenrod', 'mutationtypes_replication_strand_bias',width)
    #######################################################

    #######################################################
    ######### Plot sample based mutations ends ############
    #######################################################
    ########################################################## New Part ends ####################################################################33


    #######################################################
    ################ Plot signatures starts ###############
    #######################################################
    x_axis_transcription_labels = existingTranscriptionSignatureTypesList
    N_signatures_transcription = len(x_axis_transcription_labels)

    x_axis_replication_labels = existingReplicationSignatureTypesList
    N_signatures_replication = len(x_axis_replication_labels)

    # Transcription
    plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,None,None,N_signatures_transcription,x_axis_transcription_labels, signatures_transcribedStrandCount_list,signatures_untranscribedStrandCount_list,simulations_signatures_transcribed_medians_list, simulations_signatures_untranscribed_medians_list,signaturetype_transcription_FDR_BH_adjusted_pvalues, 'Transcribed','Untranscribed', 'All Signatures', 'royalblue', 'yellowgreen','signatures_transcription_strand_bias', width)

    # Replication
    plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,None,None,N_signatures_replication,x_axis_replication_labels,signatures_laggingStrandCount_list,signatures_leadingStrandCount_list, simulations_signatures_lagging_medians_list, simulations_signatures_leading_medians_list,signaturetype_replication_FDR_BH_adjusted_pvalues, 'Lagging','Leading', 'All Signatures', 'indianred', 'goldenrod','signatures_replication_strand_bias', width)
    #######################################################
    ################ Plot signatures ends #################
    #######################################################

    ########################################################## New Part starts ##################################################################33
    #######################################################
    #### Plot sample based signatures bar plot starts #####
    #######################################################
    for sample in sample2Signature2TranscriptionStrand2CountDict:
        if sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
            signature_2_transcribedStrandCount_list = sample2SignatureTypesTranscribedCountListDict[sample]
            signature_2_untranscribedStrandCount_list = sample2SignatureTypesUntranscribedCountListDict[sample]

            x_axis_transcription_labels = existingTranscriptionSignatureTypesList
            N_signatures_transcription = len(x_axis_transcription_labels)

            # Transcription
            if ((sample in samplesWithAtLeast10KMutations2NumberofMutationsDict) and (sample in sample2SignatureTypesTranscription_BH_FDR_Adjusted_PValues_Dict)):
                sample_signaturetypes_transcription_FDR_BH_adjusted_pvalues = sample2SignatureTypesTranscription_BH_FDR_Adjusted_PValues_Dict[sample]
                numberofMutations = samplesWithAtLeast10KMutations2NumberofMutationsDict[sample]
                samplebased_simulations_signaturetypes_transcribed_medians_list = sample2SimulationsSignatureTypesTranscribedMediansListDict[sample]
                samplebased_simulations_signaturetypes_untranscribed_medians_list = sample2SimulationsSignatureTypesUntranscribedMediansListDict[sample]

                plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,sample,numberofMutations, N_signatures_transcription, x_axis_transcription_labels,
                                             signature_2_transcribedStrandCount_list,signature_2_untranscribedStrandCount_list,
                                            samplebased_simulations_signaturetypes_transcribed_medians_list, samplebased_simulations_signaturetypes_untranscribed_medians_list,
                                             sample_signaturetypes_transcription_FDR_BH_adjusted_pvalues,
                                             'Transcribed','Untranscribed', 'All Signatures', 'royalblue', 'yellowgreen',
                                             'signatures_transcription_strand_bias', width)

    for sample in sample2Signature2ReplicationStrand2CountDict:
        if sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
            signature_2_laggingStrandCount_list = sample2SignatureTypesLaggingCountListDict[sample]
            signature_2_leadingStrandCount_list = sample2SignatureTypesLeadingCountListDict[sample]

            x_axis_replication_labels = existingReplicationSignatureTypesList
            N_signatures_replication = len(x_axis_replication_labels)

            # Replication
            if ((sample in samplesWithAtLeast10KMutations2NumberofMutationsDict) and (sample in sample2SignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict)):
                sample_signaturetypes_replication_FDR_BH_adjusted_pvalues = sample2SignatureTypesReplication_BH_FDR_Adjusted_PValues_Dict[sample]
                numberofMutations= samplesWithAtLeast10KMutations2NumberofMutationsDict[sample]
                samplebased_simulations_signaturetypes_lagging_medians_list = sample2SimulationsSignatureTypesLaggingMediansListDict[sample]
                samplebased_simulations_signaturetypes_leading_medians_list = sample2SimulationsSignatureTypesLeadingMediansListDict[sample]

                plotStrandBiasFigureWithBarPlots(outputDir,jobname,numberofSimulations,sample,numberofMutations,N_signatures_replication,x_axis_replication_labels,
                                                 signature_2_laggingStrandCount_list,signature_2_leadingStrandCount_list,samplebased_simulations_signaturetypes_lagging_medians_list,samplebased_simulations_signaturetypes_leading_medians_list,
                                                 sample_signaturetypes_replication_FDR_BH_adjusted_pvalues,
                                                 'Lagging','Leading', 'All Signatures', 'indianred', 'goldenrod',
                                             'signatures_replication_strand_bias', width)
    #######################################################
    #### Plot sample based signatures bar plot ends #######
    #######################################################
    ########################################################## New Part ends ####################################################################33


    ########################################################################
    ############ MutationType based bar plot ends ##########################
    ############ Signature based bar plot ends ## ##########################
    ########################################################################

    ########################################################################
    ######## Bar plot starts includes sample based bar plots ###############
    ##########################  Part 4 ends ################################
    ########################################################################