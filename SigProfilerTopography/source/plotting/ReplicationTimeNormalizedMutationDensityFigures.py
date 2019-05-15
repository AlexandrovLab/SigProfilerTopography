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

import matplotlib
BACKEND = 'Agg'
if matplotlib.get_backend().lower() != BACKEND.lower():
    # If backend is not set properly a call to describe will hang
    matplotlib.use(BACKEND)

from matplotlib import pyplot as plt
from SigProfilerTopography.source.commons.TopographyCommons import *


########################################################
def plotNormalizedMutationDensityFigureWithSimulations(title, ylabel, normalizedMutationDensityList, sample, signature, analysesType,indelType,barcolor,
                                                    outputDir,jobname,isFigureAugmentation,
                                                    sample2NumberofMutationsDict,
                                                    signature2NumberofMutationsDict,
                                                    sample2Signature2NumberofMutationsDict,
                                                    numberofSimulations):

    #################################################################################
    ############################# For Simulations starts ############################
    #################################################################################
    simulationsLows = None
    simulationsMedians = None
    simulationsHighs = None

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


    #Find Lows Medians Highs
    if ((listofSimulations is not None) and listofSimulations):

        stackedSimulations = np.vstack(listofSimulations)

        (rows, cols) = stackedSimulations.shape

        CI_indexLow = int(round(.05 * rows, 0))+1
        CI_indexHigh = int(round(.95 * rows, 0))-1

        simulationsLows = []
        simulationsMedians = []
        simulationsHighs = []

        for col in range(cols):
            colwise_array = stackedSimulations[:, col]
            sorted_colwise_array = np.sort(colwise_array)
            simulationsLows.append(sorted_colwise_array[CI_indexLow])
            simulationsMedians.append(np.median(sorted_colwise_array))
            simulationsHighs.append(sorted_colwise_array[CI_indexHigh])
        #################################################################################
        ############################# For Simulations ends ##############################
        #################################################################################

    ##################### legacy code starts ##########################
    # os.makedirs(os.path.join(outputDir, jobname, FIGURE, ALL, REPLICATIONTIME, analysesType), exist_ok=True)
    os.makedirs(os.path.join(outputDir, jobname, FIGURE, ALL, REPLICATIONTIME), exist_ok=True)
    os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES), exist_ok=True)

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
    if simulationsMedians is not None:
        plt.plot(x, simulationsMedians, 'o--', color='navy', label='Average Simulations Aggregated indels', linewidth=1.0, zorder =2)
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
            numberofMutations = signature2NumberofMutationsDict[signature]
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
    plt.close(fig)
    ##################### legacy code ends ############################
########################################################


#########################################################
#type can be aggregated subs, aggregated indels, MICROHOMOLOGY, INDEL, subsSignature, indelsSignature
def readNormalizedMutationData(sample,type,outputDir,jobname):
    if sample is None:
        filename = '%s_NormalizedMutationDensity.txt' %(type)
    else:
        filename = '%s_%s_NormalizedMutationDensity.txt' %(sample,type)

    filepath = os.path.join(outputDir,jobname,DATA,REPLICATIONTIME,filename)

    #Check if filepath exists
    if os.path.exists(filepath):
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
        else:
            filename = '%s_%s_sim%d_NormalizedMutationDensity.txt' % (sample, indelorSignatureorAnalysesType,simNum)
        filepath = os.path.join(outputDir,jobname, DATA,REPLICATIONTIME, filename)

        # Check if filepath exists
        if os.path.exists(filepath):
            normalizedMutationData = np.loadtxt(filepath, dtype=float)
            listofAverages.append(normalizedMutationData)
    ######################################################

    return listofAverages
#########################################################


#########################################################
def plotSignatureFigures(color,analysesType,outputDir, jobname,numberofSimulations,isFigureAugmentation,sample2NumberofMutationsDict,signature2NumberofMutationsDict,sample2Signature2NumberofMutationsDict):
    for signature in signature2NumberofMutationsDict:
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
                                                                   signature2NumberofMutationsDict,
                                                                   sample2Signature2NumberofMutationsDict,
                                                                   numberofSimulations)

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
                                                                       signature2NumberofMutationsDict,
                                                                       sample2Signature2NumberofMutationsDict,
                                                                       numberofSimulations)
#########################################################

#########################################################
def plotAllMutationTypesFigures(title,color,analysesType,indelType,outputDir,jobname,numberofSimulations,isFigureAugmentation,sample2NumberofMutationsDict,signature2NumberofMutationsDict,sample2Signature2NumberofMutationsDict):
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
                                                           signature2NumberofMutationsDict,
                                                           sample2Signature2NumberofMutationsDict,
                                                           numberofSimulations)

    ######## Sample Based AGGREGATEDINDELS starts ########
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
                                                               signature2NumberofMutationsDict,
                                                               sample2Signature2NumberofMutationsDict,
                                                               numberofSimulations)
    ######## Sample Based AGGREGATEDINDELS ends #########
#########################################################

##################################################################
def replicationTimeNormalizedMutationDensityFigures(outputDir,jobname,figureAugmentation,numberofSimulations,mutationTypes):

    isFigureAugmentation = False
    if (figureAugmentation == 'augmentation'):
        isFigureAugmentation = True

    # We know the indels type we are interested in.
    indeltypes = [MICROHOMOLOGY, REPEAT]
    analysesTypes = [AGGREGATEDINDELS, AGGREGATEDSUBSTITUTIONS, INDELBASED, SIGNATUREBASED]

    ##########################################################################################
    sample2NumberofSubsDict = getSample2NumberofSubsDict(outputDir,jobname)
    sample2NumberofIndelsDict = getSample2NumberofIndelsDict(outputDir,jobname)
    sample2NumberofDinucsDict = getDictionary(outputDir, jobname,Sample2NumberofDinucsDictFilename)

    subsSignature2NumberofMutationsDict = getSubsSignature2NumberofMutationsDict(outputDir,jobname)
    indelsSignature2NumberofMutationsDict = getIndelsSignature2NumberofMutationsDict(outputDir,jobname)
    dinucsSignature2NumberofMutationsDict = getDictionary(outputDir,jobname,DinucsSignature2NumberofMutationsDictFilename)

    sample2SubsSignature2NumberofMutationsDict = getSample2SubsSignature2NumberofMutationsDict(outputDir,jobname)
    sample2IndelsSignature2NumberofMutationsDict = getSample2IndelsSignature2NumberofMutationsDict(outputDir,jobname)
    sample2DinucsSignature2NumberofMutationsDict = getDictionary(outputDir,jobname,Sample2DinucsSignature2NumberofMutationsDictFilename)
    ##########################################################################################

    ##########################################################################################
    ##########################  Plot figures starts  #########################################
    ##########################################################################################
    if (SUBS in mutationTypes):
        plotAllMutationTypesFigures('Aggregated Substitutions','yellowgreen',AGGREGATEDSUBSTITUTIONS ,None,outputDir, jobname,numberofSimulations,isFigureAugmentation,sample2NumberofSubsDict,subsSignature2NumberofMutationsDict,sample2SubsSignature2NumberofMutationsDict)
        plotSignatureFigures('yellowgreen',SIGNATUREBASED, outputDir, jobname, numberofSimulations,isFigureAugmentation, sample2NumberofSubsDict,subsSignature2NumberofMutationsDict,sample2SubsSignature2NumberofMutationsDict)
    if (INDELS in mutationTypes):
        plotAllMutationTypesFigures('Aggregated Indels','indianred',AGGREGATEDINDELS,None,outputDir, jobname,numberofSimulations,isFigureAugmentation,sample2NumberofIndelsDict,indelsSignature2NumberofMutationsDict,sample2IndelsSignature2NumberofMutationsDict)
        plotAllMutationTypesFigures(MICROHOMOLOGY,'indianred',INDELBASED,MICROHOMOLOGY,outputDir, jobname, numberofSimulations, isFigureAugmentation, sample2NumberofIndelsDict, indelsSignature2NumberofMutationsDict,sample2IndelsSignature2NumberofMutationsDict)
        plotAllMutationTypesFigures(REPEAT, 'indianred',INDELBASED, REPEAT, outputDir, jobname, numberofSimulations,isFigureAugmentation, sample2NumberofIndelsDict, indelsSignature2NumberofMutationsDict,sample2IndelsSignature2NumberofMutationsDict)
        plotSignatureFigures('indianred', SIGNATUREBASED, outputDir, jobname, numberofSimulations, isFigureAugmentation,sample2NumberofIndelsDict, indelsSignature2NumberofMutationsDict,sample2IndelsSignature2NumberofMutationsDict)
    if (DINUCS in mutationTypes):
        plotAllMutationTypesFigures('Aggregated Dinucs','crimson',AGGREGATEDDINUCS ,None,outputDir, jobname,numberofSimulations,isFigureAugmentation,sample2NumberofDinucsDict,dinucsSignature2NumberofMutationsDict,sample2DinucsSignature2NumberofMutationsDict)
        plotSignatureFigures('crimson', SIGNATUREBASED, outputDir, jobname, numberofSimulations, isFigureAugmentation,sample2NumberofDinucsDict,dinucsSignature2NumberofMutationsDict,sample2DinucsSignature2NumberofMutationsDict)
    ##########################################################################################
    ##########################  Plot figures ends  ###########################################
    ##########################################################################################

##################################################################
