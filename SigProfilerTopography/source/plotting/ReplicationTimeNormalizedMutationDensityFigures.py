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

from SigProfilerTopography.source.commons.TopographyCommons import Table_SubsSignature_Cutoff_NumberofMutations_AverageProbability_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Table_IndelsSignature_Cutoff_NumberofMutations_AverageProbability_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Table_DinucsSignature_Cutoff_NumberofMutations_AverageProbability_Filename

from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofSubsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofIndelsDict
from SigProfilerTopography.source.commons.TopographyCommons import Sample2NumberofDinucsDictFilename

from SigProfilerTopography.source.commons.TopographyCommons import getSample2SubsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2IndelsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import Sample2DinucsSignature2NumberofMutationsDictFilename

from SigProfilerTopography.source.commons.TopographyCommons import PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_TOOL
from SigProfilerTopography.source.commons.TopographyCommons import PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_MANUSCRIPT

plt.rcParams.update({'figure.max_open_warning': 0})

########################################################
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
                                                       signature_cutoff_numberofmutations_averageprobability_df,
                                                       sample2Signature2NumberofMutationsDict,
                                                       numberofSimulations,
                                                       plot_mode):

    #################################################################################
    ############################# For Simulations starts ############################
    #################################################################################
    listofSimulations = None

    #read the simulations
    if (numberofSimulations > 0):
        if (analysesType==SIGNATUREBASED):
            #If analysesType is SIGNATUREBASED  originalTitle holds the signature
            listofSimulations = readNormalizedMutationDataForSimulations(sample,signature,outputDir,jobname,numberofSimulations)
        elif (analysesType==INDELBASED):
            listofSimulations = readNormalizedMutationDataForSimulations(sample,indelType,outputDir,jobname,numberofSimulations)
        else:
            listofSimulations = readNormalizedMutationDataForSimulations(sample,analysesType,outputDir,jobname,numberofSimulations)

    simulationsLows, simulationsMeans, simulationsHighs = takeAverage(listofSimulations)

    #################################################################################
    ############################# For Simulations ends ##############################
    #################################################################################

    ##################### legacy code starts ##########################
    os.makedirs(os.path.join(outputDir, jobname, FIGURE, REPLICATIONTIME), exist_ok=True)

    from matplotlib import rcParams
    rcParams.update({'figure.autolayout': True})

    # Note if you decrease the figure size decrease the fontsize accordingly
    fig = plt.figure(figsize=(15, 15), dpi=300)

    plt.style.use('ggplot')

    ax = plt.gca()
    # This code makes the background white.
    ax.set_facecolor('white')

    ####################################################
    # Note x get labels w.r.t. the order given here, 0 means get the 0th label from  xticks
    x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    width = 0.9  # the width of the bars
    bars=plt.bar(x, normalizedMutationDensityList, width, label='Real Somatic Mutations', color=barcolor, edgecolor="black", linewidth=3, zorder=1)

    # plt.xticks(np.arange(10),('1st', '2nd', '3rd', '4th', '5th','6th','7th','8th','9th','10th'),rotation=20)
    # also works
    # plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])

    if simulationsMeans is not None:
        sims_dashed_line=plt.plot(x, simulationsMeans, 'o--', color='black', label='Simulated Somatic Mutations', linewidth=5, zorder =2)
        if (simulationsLows is not None) and (simulationsHighs is not None):
            # if (len(simulationsLows)==len(simulationsHighs)):
            plt.fill_between(x, np.array(simulationsLows), np.array(simulationsHighs),facecolor=fillcolor, zorder =2)
    ####################################################


    ####################################################
    if plot_mode==PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_TOOL:
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


    elif plot_mode==PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_MANUSCRIPT:
        # set axis ticks
        # ax.tick_params(axis='both', which='both', length=0)
        ax.tick_params(axis='x', which='both', length=0)
        ax.tick_params(axis='y', which='both', length=0)
        # set axis labels
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax.get_yticklabels(), visible=False)
        ax.spines["bottom"].set_color('black')
        ax.spines["left"].set_color('black')
        ax.spines["top"].set_color('black')
        ax.spines["right"].set_color('black')

        # to put the legend's upper right-hand corner at (0.8,1) optimized for SBS6 in BreastCancer560 data for SigProfilerTopography Overview Figure
        # legend = ax.legend((bars[0], sims_dashed_line[0]),('Real', 'Simulated'),prop={'size': 50}, loc='upper right', bbox_to_anchor = (0.8, 1))
        # to put the legend's upper left corner at (x,y) optimized for SBS2 in BreastCancer560 data for Replication Overview Figure
        legend = ax.legend((bars[0], sims_dashed_line[0]),('Real', 'Simulated'),prop={'size': 43}, loc='upper left', bbox_to_anchor = (0.02, 0.9))

        if (legend is not None):
            frame = legend.get_frame()
            frame.set_facecolor('white')
            frame.set_edgecolor('black')
    ####################################################


    ########################################################################
    if sample is None:
        if (analysesType == INDELBASED):
            figureName = '%s_replication_time.png' % (indelType.replace(' ', ''))
        elif (analysesType == SIGNATUREBASED):
            #[signature cutoff numberofMutations averageProbability]
            numberofMutations=int(signature_cutoff_numberofmutations_averageprobability_df[signature_cutoff_numberofmutations_averageprobability_df['signature']==signature]['number_of_mutations'].values[0])
            figureName = '%s_%d_replication_time.png' % (signature.replace(' ', ''), numberofMutations)
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
    ########################################################################


    ########################################################################
    if (sample is None):
        figureFile = os.path.join(outputDir, jobname, FIGURE, REPLICATIONTIME, figureName)
    else:
        os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, REPLICATIONTIME), exist_ok=True)
        figureFile = os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, REPLICATIONTIME, figureName)
    ########################################################################

    fig.savefig(figureFile)
    plt.cla()
    plt.close(fig)
    ##################### legacy code ends ############################
########################################################


#########################################################
#analysisType can be aggregated subs, aggregated indels, MICROHOMOLOGY, INDEL, subsSignature, indelsSignature, dbsSignature
def readNormalizedMutationData(sample,analysisType,outputDir,jobname):
    if sample is None:
        filename = '%s_NormalizedMutationDensity.txt' %(analysisType)
        if ((analysisType==AGGREGATEDSUBSTITUTIONS) or (analysisType==AGGREGATEDINDELS) or (analysisType==AGGREGATEDDINUCS) or (analysisType==MICROHOMOLOGY) or (analysisType==REPEAT)):
            filepath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, analysisType, filename)
        else:
            filepath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, SIGNATUREBASED, filename)

    else:
        filename = '%s_%s_NormalizedMutationDensity.txt' %(sample,analysisType)
        if ((analysisType==AGGREGATEDSUBSTITUTIONS) or (analysisType==AGGREGATEDINDELS) or (analysisType==AGGREGATEDDINUCS) or (analysisType==MICROHOMOLOGY) or (analysisType==REPEAT)):
            filepath = os.path.join(outputDir,jobname,DATA,SAMPLES,sample,REPLICATIONTIME,analysisType,filename)
        else:
            filepath = os.path.join(outputDir,jobname,DATA,SAMPLES,sample,REPLICATIONTIME,SIGNATUREBASED,filename)

    #Check if filepath exists
    if os.path.exists(filepath):
        # normalizedMutationData_df = pd.read_table(filepath, sep=" ", comment='#', header=None)
        normalizedMutationData_df = pd.read_csv(filepath, sep=" ", comment='#', header=None)
        normalizedMutationData_df.dropna(axis=1, how='any', inplace=True)
        return normalizedMutationData_df
    else:
        return None
#########################################################


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


#########################################################
def plotSignatureFigures(color,fillcolor,analysesType,outputDir, jobname,numberofSimulations,sample_based,sample2NumberofMutationsDict,signature_cutoff_numberofmutations_averageprobability_df,sample2Signature2NumberofMutationsDict,plot_mode):
    for signature in signature_cutoff_numberofmutations_averageprobability_df['signature'].unique():
        # We check such file exists or not
        normalizedMutationData = readNormalizedMutationData(None, signature, outputDir, jobname)

        if (normalizedMutationData is not None):
            normalizedMutationData = normalizedMutationData.iloc[0].tolist()
            # if not all([v == 0.0 for v in normalizedMutationData]):
            # use all generator for all true check
            if not all(v == 0.0 for v in normalizedMutationData):
                # print('For %s: plot signature based replication time figure' % signature)
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
                                                                   signature_cutoff_numberofmutations_averageprobability_df,
                                                                   sample2Signature2NumberofMutationsDict,
                                                                   numberofSimulations,
                                                                   plot_mode)

    if sample_based:
        for sample in sample2Signature2NumberofMutationsDict:
            for signature in sample2Signature2NumberofMutationsDict[sample]:
                normalizedMutationData = readNormalizedMutationData(sample, signature, outputDir, jobname)

                if (normalizedMutationData is not None):
                    normalizedMutationData = normalizedMutationData.iloc[0].tolist()
                    # if not all([v == 0.0 for v in normalizedMutationData]):
                    # use all generator for all true check
                    if not all(v == 0.0 for v in normalizedMutationData):
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
                                                                           signature_cutoff_numberofmutations_averageprobability_df,
                                                                           sample2Signature2NumberofMutationsDict,
                                                                           numberofSimulations,
                                                                           plot_mode)
#########################################################

#########################################################
def plotAllMutationTypesFigures(title,color,fillcolor,analysesType,indelType,outputDir,jobname,numberofSimulations,sample_based,sample2NumberofMutationsDict,signature_cutoff_numberofmutations_averageprobability_df,sample2Signature2NumberofMutationsDict,plot_mode):
    if (analysesType == INDELBASED):
        normalizedMutationData = readNormalizedMutationData(None,indelType,outputDir,jobname)
    else:
        normalizedMutationData = readNormalizedMutationData(None, analysesType, outputDir, jobname)

    if (normalizedMutationData is not None):
        normalizedMutationData = normalizedMutationData.iloc[0].tolist()
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

    ######## Sample Based starts ########
    if sample_based:
        for sample in sample2NumberofMutationsDict:
            if (analysesType == INDELBASED):
                normalizedMutationData = readNormalizedMutationData(sample, indelType, outputDir, jobname)
            else:
                normalizedMutationData = readNormalizedMutationData(sample, analysesType, outputDir, jobname)
            if (normalizedMutationData is not None):
                normalizedMutationData = normalizedMutationData.iloc[0].tolist()
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
    ######## Sample Based ends #########
#########################################################

##################################################################
def replicationTimeNormalizedMutationDensityFigures(outputDir,jobname,numberofSimulations,sample_based,mutationTypes,plot_mode):

    # Initialize these dataframes as empty dataframe
    # We will read these dataframes if there is the corresponding data
    subsSignature_cutoff_numberofmutations_averageprobability_df = pd.DataFrame()
    dinucsSignature_cutoff_numberofmutations_averageprobability_df = pd.DataFrame()
    indelsSignature_cutoff_numberofmutations_averageprobability_df = pd.DataFrame()

    ##########################################################################################
    for mutation_type_context in mutationTypes:
        if (mutation_type_context in SBS_CONTEXTS):
            subsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(os.path.join(outputDir, jobname, DATA,Table_SubsSignature_Cutoff_NumberofMutations_AverageProbability_Filename),sep='\t', header=0,dtype={'cutoff': np.float32,'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})
    if (DBS in mutationTypes):
        dinucsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(os.path.join(outputDir, jobname, DATA,Table_DinucsSignature_Cutoff_NumberofMutations_AverageProbability_Filename),sep='\t', header=0,dtype={'cutoff': np.float32,'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})
    if (ID in mutationTypes):
        indelsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(os.path.join(outputDir, jobname, DATA,Table_IndelsSignature_Cutoff_NumberofMutations_AverageProbability_Filename),sep='\t', header=0,dtype={'cutoff': np.float32,'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})

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

    ##########################################################################################
    ##########################  Plot figures starts  #########################################
    ##########################################################################################
    if (not subsSignature_cutoff_numberofmutations_averageprobability_df.empty):
        #Formerly it was yellowgreen
        plotAllMutationTypesFigures('Aggregated Substitutions','royalblue','lightblue',AGGREGATEDSUBSTITUTIONS ,None,outputDir, jobname,numberofSimulations,sample_based,sample2NumberofSubsDict,subsSignature_cutoff_numberofmutations_averageprobability_df,sample2SubsSignature2NumberofMutationsDict,plot_mode)
        plotSignatureFigures('royalblue','lightblue',SIGNATUREBASED, outputDir, jobname, numberofSimulations,sample_based, sample2NumberofSubsDict,subsSignature_cutoff_numberofmutations_averageprobability_df,sample2SubsSignature2NumberofMutationsDict,plot_mode)
    if (not dinucsSignature_cutoff_numberofmutations_averageprobability_df.empty):
        plotAllMutationTypesFigures('Aggregated Dinucs','crimson','lightpink',AGGREGATEDDINUCS ,None,outputDir, jobname,numberofSimulations,sample_based,sample2NumberofDinucsDict,dinucsSignature_cutoff_numberofmutations_averageprobability_df,sample2DinucsSignature2NumberofMutationsDict,plot_mode)
        plotSignatureFigures('crimson','lightpink', SIGNATUREBASED, outputDir, jobname, numberofSimulations,sample_based,sample2NumberofDinucsDict,dinucsSignature_cutoff_numberofmutations_averageprobability_df,sample2DinucsSignature2NumberofMutationsDict,plot_mode)
    if (not indelsSignature_cutoff_numberofmutations_averageprobability_df.empty):
        plotAllMutationTypesFigures('Aggregated Indels','yellowgreen','lightgreen',AGGREGATEDINDELS,None,outputDir, jobname,numberofSimulations,sample_based,sample2NumberofIndelsDict,indelsSignature_cutoff_numberofmutations_averageprobability_df,sample2IndelsSignature2NumberofMutationsDict,plot_mode)
        plotAllMutationTypesFigures(MICROHOMOLOGY,'yellowgreen','lightgreen',INDELBASED,MICROHOMOLOGY,outputDir, jobname, numberofSimulations,sample_based, sample2NumberofIndelsDict, indelsSignature_cutoff_numberofmutations_averageprobability_df,sample2IndelsSignature2NumberofMutationsDict,plot_mode)
        plotAllMutationTypesFigures(REPEAT, 'yellowgreen','lightgreen',INDELBASED, REPEAT, outputDir, jobname, numberofSimulations,sample_based, sample2NumberofIndelsDict, indelsSignature_cutoff_numberofmutations_averageprobability_df,sample2IndelsSignature2NumberofMutationsDict,plot_mode)
        plotSignatureFigures('yellowgreen','lightgreen', SIGNATUREBASED, outputDir, jobname, numberofSimulations, sample_based,sample2NumberofIndelsDict, indelsSignature_cutoff_numberofmutations_averageprobability_df,sample2IndelsSignature2NumberofMutationsDict,plot_mode)
    ##########################################################################################
    ##########################  Plot figures ends  ###########################################
    ##########################################################################################

##################################################################
