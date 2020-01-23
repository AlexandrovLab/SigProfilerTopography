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

from SigProfilerTopography.source.commons.TopographyCommons import SBS96
from SigProfilerTopography.source.commons.TopographyCommons import ID
from SigProfilerTopography.source.commons.TopographyCommons import DBS

from SigProfilerTopography.source.commons.TopographyCommons import DATA
from SigProfilerTopography.source.commons.TopographyCommons import FIGURE
from SigProfilerTopography.source.commons.TopographyCommons import ALL
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

from SigProfilerTopography.source.commons.TopographyCommons import SubsSignature2PropertiesListDictFilename
from SigProfilerTopography.source.commons.TopographyCommons import IndelsSignature2PropertiesListDictFilename
from SigProfilerTopography.source.commons.TopographyCommons import DinucsSignature2PropertiesListDictFilename

from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofSubsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofIndelsDict
from SigProfilerTopography.source.commons.TopographyCommons import Sample2NumberofDinucsDictFilename

from SigProfilerTopography.source.commons.TopographyCommons import getSample2SubsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2IndelsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import Sample2DinucsSignature2NumberofMutationsDictFilename

plt.rcParams.update({'figure.max_open_warning': 0})

########################################################
def plotNormalizedMutationDensityFigureWithSimulations(title, ylabel, normalizedMutationDensityList, sample, signature, analysesType,indelType,barcolor,
                                                    outputDir,jobname,isFigureAugmentation,
                                                    sample2NumberofMutationsDict,
                                                    signature2PropertiesListDict,
                                                    sample2Signature2NumberofMutationsDict,
                                                    numberofSimulations):

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


    # new way
    simulationsLows, simulationsMeans, simulationsHighs = takeAverage(listofSimulations)

    #old way
    #Find Lows Medians Highs
    # if ((listofSimulations is not None) and listofSimulations):
    #
    #     stackedSimulations = np.vstack(listofSimulations)
    #
    #     (rows, cols) = stackedSimulations.shape
    #
    #     CI_indexLow = int(round(.05 * rows, 0))+1
    #     CI_indexHigh = int(round(.95 * rows, 0))-1
    #
    #     simulationsLows = []
    #     simulationsMedians = []
    #     simulationsHighs = []
    #
    #     for col in range(cols):
    #         colwise_array = stackedSimulations[:, col]
    #         sorted_colwise_array = np.sort(colwise_array)
    #
    #         if (CI_indexLow < sorted_colwise_array.shape[0]):
    #             simulationsLows.append(sorted_colwise_array[CI_indexLow])
    #         simulationsMedians.append(np.median(sorted_colwise_array))
    #         if (CI_indexHigh < sorted_colwise_array.shape[0]):
    #             simulationsHighs.append(sorted_colwise_array[CI_indexHigh])
    #################################################################################
    ############################# For Simulations ends ##############################
    #################################################################################

    ##################### legacy code starts ##########################
    os.makedirs(os.path.join(outputDir, jobname, FIGURE, ALL, REPLICATIONTIME), exist_ok=True)

    from matplotlib import rcParams
    rcParams.update({'figure.autolayout': True})

    # Note if you decrease the figure size decrease the fontsize accordingly
    fig = plt.figure(figsize=(15, 15), dpi=300)

    plt.style.use('ggplot')

    ax = plt.gca()
    # This code makes the background white.
    ax.set_facecolor('white')
    for edge_i in ['left']:
        ax.spines[edge_i].set_edgecolor("black")
        ax.spines[edge_i].set_linewidth(3)
        # This code draws line only between [0,1]
        ax.spines[edge_i].set_bounds(0, 1)

    width = 0.9  # the width of the bars

    # Note x get labels w.r.t. the order given here, 0 means get the 0th label from  xticks
    x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    # plt.xticks(np.arange(10),('1st', '2nd', '3rd', '4th', '5th','6th','7th','8th','9th','10th'),rotation=20)
    plt.yticks(np.arange(0, 1.01, step=0.2))

    # also works
    # plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])

    plt.bar(x, normalizedMutationDensityList, width, color=barcolor, edgecolor="black", linewidth=3, zorder=1)

    if simulationsMeans is not None:
        plt.plot(x, simulationsMeans, 'o--', color='navy', label='Average Simulations Aggregated indels', linewidth=1.0, zorder =2)
        if (simulationsLows is not None) and (simulationsHighs is not None):
            # if (len(simulationsLows)==len(simulationsHighs)):
            plt.fill_between(x, np.array(simulationsLows), np.array(simulationsHighs),facecolor='lightblue', zorder =2)

    # This code puts some extra space below 0 and above 1.0
    plt.ylim(-0.01, 1.01)

    plt.tick_params(axis='y', which='major', labelsize=40, width=3, length=10)
    plt.tick_params(axis='y', which='minor', labelsize=40, width=3, length=10)

    plt.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)  # labels along the bottom edge are off

    ####################################################
    if (isFigureAugmentation):
        plt.title(jobname + ' ' + title, fontsize=40, fontweight='bold')
    else:
        plt.title(title, fontsize=40, fontweight='bold')
    ####################################################

    plt.xlabel('Early <--- Replication Time ---> Late', fontsize=40, fontweight='semibold')
    plt.ylabel(ylabel, fontsize=40, fontweight='semibold')

    ########################################################################
    if sample is None:
        if (analysesType == INDELBASED):
            figureName = '%s_replication_time.png' % (indelType.replace(' ', ''))
        elif (analysesType == SIGNATUREBASED):
            #[cutoff numberofMutations averageProbability]
            numberofMutations = signature2PropertiesListDict[signature][1]
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
        # figureFile = os.path.join(outputDir, jobname, FIGURE, ALL, REPLICATIONTIME, analysesType, figureName)
        figureFile = os.path.join(outputDir, jobname, FIGURE, ALL, REPLICATIONTIME, figureName)

    else:
        # os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, REPLICATIONTIME, analysesType), exist_ok=True)
        os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, REPLICATIONTIME), exist_ok=True)
        # figureFile = os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, REPLICATIONTIME, analysesType, figureName)
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
        #TODO Change it to read_csv
        normalizedMutationData_df = pd.read_table(filepath, sep=" ", comment='#', header=None)
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
def plotSignatureFigures(color,analysesType,outputDir, jobname,numberofSimulations,sample_based,isFigureAugmentation,sample2NumberofMutationsDict,signature2PropertiesListDict,sample2Signature2NumberofMutationsDict):
    for signature in signature2PropertiesListDict:
        # We check such file exists or not
        normalizedMutationData = readNormalizedMutationData(None, signature, outputDir, jobname)

        if (normalizedMutationData is not None):
            normalizedMutationData = normalizedMutationData.iloc[0].tolist()
            # if not all([v == 0.0 for v in normalizedMutationData]):
            # use all generator for all true check
            if not all(v == 0.0 for v in normalizedMutationData):
                # print('For %s: plot signature based replication time figure' % signature)
                plotNormalizedMutationDensityFigureWithSimulations(signature,
                                                                   'Normalized\nsingle point mutation density',
                                                                   normalizedMutationData, None, signature,
                                                                   analysesType, None,
                                                                   color, outputDir, jobname,
                                                                   isFigureAugmentation,
                                                                   sample2NumberofMutationsDict,
                                                                   signature2PropertiesListDict,
                                                                   sample2Signature2NumberofMutationsDict,
                                                                   numberofSimulations)

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
                                                                           'Normalized\nsingle point mutation density',
                                                                           normalizedMutationData, sample, signature,
                                                                           analysesType, None,
                                                                           color, outputDir, jobname,
                                                                           isFigureAugmentation,
                                                                           sample2NumberofMutationsDict,
                                                                           signature2PropertiesListDict,
                                                                           sample2Signature2NumberofMutationsDict,
                                                                           numberofSimulations)
#########################################################

#########################################################
def plotAllMutationTypesFigures(title,color,analysesType,indelType,outputDir,jobname,numberofSimulations,sample_based,isFigureAugmentation,sample2NumberofMutationsDict,signature2PropertiesListDict,sample2Signature2NumberofMutationsDict):
    if (analysesType == INDELBASED):
        normalizedMutationData = readNormalizedMutationData(None,indelType,outputDir,jobname)
    else:
        normalizedMutationData = readNormalizedMutationData(None, analysesType, outputDir, jobname)

    if (normalizedMutationData is not None):
        normalizedMutationData = normalizedMutationData.iloc[0].tolist()
        plotNormalizedMutationDensityFigureWithSimulations(title, '\nNormalized mutation density',
                                                           normalizedMutationData, None, None, analysesType, indelType,
                                                           color,outputDir, jobname, isFigureAugmentation,
                                                           sample2NumberofMutationsDict,
                                                           signature2PropertiesListDict,
                                                           sample2Signature2NumberofMutationsDict,
                                                           numberofSimulations)

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
                                                                   normalizedMutationData, sample, None, analysesType, indelType,
                                                                   color, outputDir, jobname, isFigureAugmentation,
                                                                   sample2NumberofMutationsDict,
                                                                   signature2PropertiesListDict,
                                                                   sample2Signature2NumberofMutationsDict,
                                                                   numberofSimulations)
    ######## Sample Based ends #########
#########################################################

##################################################################
def replicationTimeNormalizedMutationDensityFigures(outputDir,jobname,figureAugmentation,numberofSimulations,sample_based,mutationTypes):

    jobnamePath = os.path.join(outputDir,jobname,FIGURE,ALL,REPLICATIONTIME)
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


    ##########################################################################################
    subsSignature2PropertiesListDict = getDictionary(outputDir, jobname,SubsSignature2PropertiesListDictFilename)
    indelsSignature2PropertiesListDict = getDictionary(outputDir, jobname,IndelsSignature2PropertiesListDictFilename)
    dinucsSignature2PropertiesListDict = getDictionary(outputDir, jobname,DinucsSignature2PropertiesListDictFilename)

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
    if (SBS96 in mutationTypes):
        plotAllMutationTypesFigures('Aggregated Substitutions','yellowgreen',AGGREGATEDSUBSTITUTIONS ,None,outputDir, jobname,numberofSimulations,sample_based,isFigureAugmentation,sample2NumberofSubsDict,subsSignature2PropertiesListDict,sample2SubsSignature2NumberofMutationsDict)
        plotSignatureFigures('yellowgreen',SIGNATUREBASED, outputDir, jobname, numberofSimulations,sample_based,isFigureAugmentation, sample2NumberofSubsDict,subsSignature2PropertiesListDict,sample2SubsSignature2NumberofMutationsDict)
    if (ID in mutationTypes):
        plotAllMutationTypesFigures('Aggregated Indels','indianred',AGGREGATEDINDELS,None,outputDir, jobname,numberofSimulations,sample_based,isFigureAugmentation,sample2NumberofIndelsDict,indelsSignature2PropertiesListDict,sample2IndelsSignature2NumberofMutationsDict)
        plotAllMutationTypesFigures(MICROHOMOLOGY,'indianred',INDELBASED,MICROHOMOLOGY,outputDir, jobname, numberofSimulations,sample_based, isFigureAugmentation, sample2NumberofIndelsDict, indelsSignature2PropertiesListDict,sample2IndelsSignature2NumberofMutationsDict)
        plotAllMutationTypesFigures(REPEAT, 'indianred',INDELBASED, REPEAT, outputDir, jobname, numberofSimulations,sample_based,isFigureAugmentation, sample2NumberofIndelsDict, indelsSignature2PropertiesListDict,sample2IndelsSignature2NumberofMutationsDict)
        plotSignatureFigures('indianred', SIGNATUREBASED, outputDir, jobname, numberofSimulations, sample_based, isFigureAugmentation,sample2NumberofIndelsDict, indelsSignature2PropertiesListDict,sample2IndelsSignature2NumberofMutationsDict)
    if (DBS in mutationTypes):
        plotAllMutationTypesFigures('Aggregated Dinucs','crimson',AGGREGATEDDINUCS ,None,outputDir, jobname,numberofSimulations,sample_based,isFigureAugmentation,sample2NumberofDinucsDict,dinucsSignature2PropertiesListDict,sample2DinucsSignature2NumberofMutationsDict)
        plotSignatureFigures('crimson', SIGNATUREBASED, outputDir, jobname, numberofSimulations,sample_based, isFigureAugmentation,sample2NumberofDinucsDict,dinucsSignature2PropertiesListDict,sample2DinucsSignature2NumberofMutationsDict)
    ##########################################################################################
    ##########################  Plot figures ends  ###########################################
    ##########################################################################################

##################################################################
