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

#############################################################
current_abs_path = os.path.abspath(os.path.dirname(__file__))
print('ReplicationTimeNormalizedMutationDensityFigures.py current_abs_path:%s' %(current_abs_path))
#############################################################

commonsPath = os.path.join(current_abs_path, '..','commons')
sys.path.append(commonsPath)

from SigProfilerTopography.source.commons.TopographyCommons import *

#########################################################
def autolabel(rects):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        #plt.text(rect.get_x() + rect.get_width()/2., 1.01*height,'%d' % int(height),ha='center', va='bottom')
        plt.text(rect.get_x() + rect.get_width()/2., 1.01*height,'%f' % height,ha='center', va='bottom')

#########################################################


########################################################
#Nov 7, 2018 with Simulations
def plotNormalizedMutationDensityFigureWithSimulations(originalTitle, ylabel, normalizedMutationDensityList, sample, signature, analysesType,barcolor,jobname,isFigureAugmentation,
                                        sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
                                        samplesWithAtLeast10KMutations2NumberofMutationsDict,
                                        sample2NumberofIndelsDict,
                                        signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
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
            listofSimulations = readNormalizedMutationDataForSimulations(sample,signature,jobname,numberofSimulations)
        elif (analysesType==INDELBASED):
            listofSimulations = readNormalizedMutationDataForSimulations(sample,originalTitle,jobname,numberofSimulations)
        else:
            listofSimulations = readNormalizedMutationDataForSimulations(sample,analysesType,jobname,numberofSimulations)


    #Find Lows Medians Highs
    if ((listofSimulations is not None) and listofSimulations):

        stackedSimulations = np.vstack(listofSimulations)
        # print('stackedSimulations.shape')
        # print(stackedSimulations.shape)

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
    os.makedirs(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, FIGURE, ALL, REPLICATIONTIME, analysesType), exist_ok=True)
    os.makedirs(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, FIGURE, SAMPLES), exist_ok=True)

    from matplotlib import rcParams
    rcParams.update({'figure.autolayout': True})

    # Note if you decrease teh figure size decrease the fontsize accordingly
    fig = plt.figure(figsize=(15, 15), dpi=300)
    # fig = plt.figure(dpi=300)

    plt.style.use('ggplot')
    # print(normalizedMutationDensityList)

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
        # plt.plot(x, simulationsMedians, color='navy', linestyle='--',label='Average Simulations Aggregated indels', linewidth=1.0)
        # plt.plot(x, simulationsMedians, 'o-', color='navy', label='Average Simulations Aggregated indels', linewidth=1.0)
        plt.plot(x, simulationsMedians, 'o--', color='navy', label='Average Simulations Aggregated indels', linewidth=1.0, zorder =2)

        # '.-' i too small
        # plt.plot(x, simulationsMedians, '.-',color='navy', label='Average Simulations Aggregated indels', linewidth=1.0)
        #interpolate is meaningful when two curves cross each other
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
    if (analysesType == INDELBASED):
        if originalTitle == MICROHOMOLOGY:
            title = 'Microhomology-mediated \nindels'
        elif originalTitle == REPEAT:
            title = 'Repeat-mediated \nindels'
    else:
        title = originalTitle
    ####################################################

    ####################################################
    if (isFigureAugmentation):
        plt.title(jobname + ' ' + title, fontsize=50, fontweight='bold')
    else:
        plt.title(title, fontsize=50, fontweight='bold')
    ####################################################

    plt.xlabel('Early <--- Replication Time ---> Late', fontsize=40, fontweight='semibold')
    plt.ylabel(ylabel, fontsize=40, fontweight='semibold')

    ########################################################################
    if sample is None:
        if (analysesType == INDELBASED):
            figureName = '%s_%s_NormalizedMutationDensity.png' % (originalTitle.replace(' ', ''), analysesType)
        elif (analysesType == SIGNATUREBASED):
            numberofMutations = signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict[originalTitle]
            figureName = '%s_%d_%s_NormalizedMutationDensity.png' % (
            originalTitle.replace(' ', ''), numberofMutations, analysesType)
        else:
            # AGGREGATEDSUBSTITUTIONS
            # AGGREGATEDINDELS
            figureName = '%s_NormalizedMutationDensity.png' % (analysesType)

    else:
        if (analysesType == INDELBASED):
            figureName = '%s_%s_%s_NormalizedMutationDensity.png' % (
            originalTitle.replace(' ', ''), sample, analysesType)
        elif (analysesType == SIGNATUREBASED):
            numberofMutations = sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict[sample][
                originalTitle]
            figureName = '%s_%s_%d_%s_NormalizedMutationDensity.png' % (
            originalTitle.replace(' ', ''), sample, numberofMutations, analysesType)
        elif (analysesType == AGGREGATEDSUBSTITUTIONS):
            numberofMutations = samplesWithAtLeast10KMutations2NumberofMutationsDict[sample]
            figureName = '%s_%d_%s_NormalizedMutationDensity.png' % (sample, numberofMutations, analysesType)
        elif (analysesType == AGGREGATEDINDELS):
            numberofIndels = sample2NumberofIndelsDict[sample]
            figureName = '%s_%d_%s_NormalizedMutationDensity.png' % (sample, numberofIndels, analysesType)
    ########################################################################


    ########################################################################
    if (sample is None):
        figureFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, FIGURE, ALL, REPLICATIONTIME, analysesType, figureName)
    else:
        os.makedirs(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, FIGURE, SAMPLES, sample, REPLICATIONTIME, analysesType), exist_ok=True)
        figureFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, FIGURE, SAMPLES, sample, REPLICATIONTIME, analysesType, figureName)
    ########################################################################

    fig.savefig(figureFile)
    plt.close(fig)
    ##################### legacy code ends ############################
########################################################


#########################################################
def readNormalizedMutationData(sample,indelorSignatureorAnalysesType,jobname):
    if sample is None:
        filename = '%s_NormalizedMutationDensity.txt' %(indelorSignatureorAnalysesType)
    else:
        filename = '%s_%s_NormalizedMutationDensity.txt' %(sample,indelorSignatureorAnalysesType)

    filepath = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,OUTPUT,jobname,DATA,REPLICATIONTIME,NORMALIZED_MUTATION_DENSITY,filename)

    #Check if filepath exists
    if os.path.exists(filepath):
        normalizedMutationData_df = pd.read_table(filepath, sep=" ", comment='#', header=None)
        normalizedMutationData_df.dropna(axis=1, how='any', inplace=True)
        return normalizedMutationData_df
    else:
        return None
#########################################################


#########################################################
def readNormalizedMutationDataForSimulations(sample, indelorSignatureorAnalysesType, jobname,numberofSimulations):
    listofAverages = []

    if sample is None:
        filename = '%s_NormalizedMutationDensity.txt' % (indelorSignatureorAnalysesType)
    else:
        filename = '%s_%s_NormalizedMutationDensity.txt' % (sample,indelorSignatureorAnalysesType)

    ######################################################
    for i in range(1, numberofSimulations + 1):
        simjobname = '%s_Sim%d' %(jobname,i)
        filepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, simjobname, DATA,REPLICATIONTIME, NORMALIZED_MUTATION_DENSITY, filename)

        # Check if filepath exists
        if os.path.exists(filepath):
            #new way
            normalizedMutationData = np.loadtxt(filepath, dtype=float)
            listofAverages.append(normalizedMutationData)
            # print('for debug filepath: %s' %(filepath))
            # print('for debug normalizedMutationData')
            # print(normalizedMutationData)
    ######################################################

    return listofAverages
#########################################################





##################################################################
def replicationTimeNormalizedMutationDensityFigures(jobname,figureAugmentation,numberofSimulations):

    isFigureAugmentation = False
    if (figureAugmentation == 'augmentation'):
        isFigureAugmentation = True


    # We know the indels type we are interested in.
    indeltypes = [MICROHOMOLOGY, REPEAT]

    #Current one
    #SIGNATUREBASED_SAMPLEBASED is removed
    analysesTypes = [AGGREGATEDINDELS, AGGREGATEDSUBSTITUTIONS, INDELBASED, SIGNATUREBASED]

    ##########################################################################################
    ################### Read the necessary dictionaries starts ###############################
    ##########################################################################################

    ##########################################################################################
    signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = {}
    signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath = os.path.join(current_abs_path,
                                                                                           ONE_DIRECTORY_UP,
                                                                                           ONE_DIRECTORY_UP, OUTPUT,
                                                                                           jobname, DATA,
                                                                                           SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilename)

    if (os.path.exists(signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath)):
        signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = readDictionary(signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath)
    ##########################################################################################

    ##########################################################################################
    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = {}
    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath = os.path.join(current_abs_path,
                                                                                                  ONE_DIRECTORY_UP,
                                                                                                  ONE_DIRECTORY_UP,
                                                                                                  OUTPUT, jobname, DATA,
                                                                                                  Sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilename)

    if (os.path.exists(sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath)):
        sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = readDictionary(sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilePath)
    ##########################################################################################

    ##########################################################################################
    samplesWithAtLeast10KMutations2NumberofMutationsDict = {}
    samplesWithAtLeast10KMutations2NumberofMutationsDictFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP,
                                                                                ONE_DIRECTORY_UP, OUTPUT, jobname, DATA,
                                                                                SamplesWithAtLeast10KMutations2NumberofMutationsDictFilename)

    if (os.path.exists(samplesWithAtLeast10KMutations2NumberofMutationsDictFilePath)):
        samplesWithAtLeast10KMutations2NumberofMutationsDict = readDictionary(samplesWithAtLeast10KMutations2NumberofMutationsDictFilePath)
    ##########################################################################################


    ##########################################################################################
    sample2NumberofIndelsDict = {}
    sample2NumberofIndelsDictFilePath = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,OUTPUT,jobname,DATA,Samples2NumberofIndelsDictFilename)

    if (os.path.exists(sample2NumberofIndelsDictFilePath)):
        sample2NumberofIndelsDict = readDictionary(sample2NumberofIndelsDictFilePath)
    ##########################################################################################

    ##########################################################################################
    ################### Read the necessary dictionaries ends #################################
    ##########################################################################################


    ##########################################################################################
    ##########################  Plot figures starts  #########################################
    ##########################################################################################


    ##########################################################################################
    ##########################  Simulations starts  ##########################################
    ##########################################################################################
    for analysesType in analysesTypes:
        if (analysesType == AGGREGATEDINDELS):
            normalizedMutationData = readNormalizedMutationData(None, analysesType, jobname)
            if (normalizedMutationData is not None):
                normalizedMutationData = normalizedMutationData.iloc[0].tolist()
                plotNormalizedMutationDensityFigureWithSimulations('Aggregated Indels', '\nNormalized mutation density',
                                                    normalizedMutationData, None, None,analysesType, 'indianred',
                                                    jobname, isFigureAugmentation,
                                                    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
                                                    samplesWithAtLeast10KMutations2NumberofMutationsDict,
                                                    sample2NumberofIndelsDict,
                                                    signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,numberofSimulations)

            ######## Sample Based AGGREGATEDINDELS starts ########
            for sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
                normalizedMutationData = readNormalizedMutationData(sample, analysesType, jobname)
                if (normalizedMutationData is not None):
                    normalizedMutationData = normalizedMutationData.iloc[0].tolist()
                    # TODO
                    plotNormalizedMutationDensityFigureWithSimulations('Aggregated Indels', '\nNormalized mutation density',
                                                        normalizedMutationData, sample,None, analysesType,
                                                        'indianred', jobname, isFigureAugmentation,
                                                        sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
                                                        samplesWithAtLeast10KMutations2NumberofMutationsDict,
                                                        sample2NumberofIndelsDict,
                                                        signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,numberofSimulations)
            ######## Sample Based AGGREGATEDINDELS ends #########



        # Data is generated in ReplicationTimeAnalysis.py
        # If indels based figures not needed/wanted  following lines can be commented.
        elif (analysesType == INDELBASED):
            for indelType in indeltypes:
                # if (checkValidnessForIndelsBased(analysesType,jobname,indelType)):
                normalizedMutationData = readNormalizedMutationData(None, indelType, jobname)
                if (normalizedMutationData is not None):
                    normalizedMutationData = normalizedMutationData.iloc[0].tolist()
                    plotNormalizedMutationDensityFigureWithSimulations(indelType, 'Normalized mutation density',
                                                        normalizedMutationData, None, None,analysesType,
                                                        'indianred', jobname, isFigureAugmentation,
                                                        sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
                                                        samplesWithAtLeast10KMutations2NumberofMutationsDict,
                                                        sample2NumberofIndelsDict,
                                                        signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,numberofSimulations)

                ######## Sample Based INDELBASED starts ########
                for sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
                    normalizedMutationData = readNormalizedMutationData(sample, indelType, jobname)
                    if (normalizedMutationData is not None):
                        normalizedMutationData = normalizedMutationData.iloc[0].tolist()
                        # TODO
                        plotNormalizedMutationDensityFigureWithSimulations(indelType, 'Normalized mutation density',
                                                            normalizedMutationData, sample, None,analysesType,
                                                            'indianred', jobname, isFigureAugmentation,
                                                            sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
                                                            samplesWithAtLeast10KMutations2NumberofMutationsDict,
                                                            sample2NumberofIndelsDict,
                                                            signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,numberofSimulations)
                ######## Sample Based INDELBASED starts ########

        elif (analysesType == AGGREGATEDSUBSTITUTIONS):
            normalizedMutationData = readNormalizedMutationData(None, analysesType, jobname)

            if (normalizedMutationData is not None):
                normalizedMutationData = normalizedMutationData.iloc[0].tolist()
                plotNormalizedMutationDensityFigureWithSimulations('Aggregated Substitutions',
                                                    'Normalized\nsingle point mutation density',
                                                    normalizedMutationData, None, None, analysesType, 'yellowgreen',
                                                    jobname, isFigureAugmentation,
                                                    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
                                                    samplesWithAtLeast10KMutations2NumberofMutationsDict,
                                                    sample2NumberofIndelsDict,
                                                    signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,numberofSimulations)

                # #Do we need this anymore?
                # #get the signatures, they may have spaces before, after or in between
                # signaturesFilePath = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,OUTPUT,jobname,DATA,SignatureFilename)
                # signatures = readSignatureList(signaturesFilePath)

            for sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
                normalizedMutationData = readNormalizedMutationData(sample, analysesType, jobname)

                if (normalizedMutationData is not None):
                    normalizedMutationData = normalizedMutationData.iloc[0].tolist()
                    plotNormalizedMutationDensityFigureWithSimulations('Aggregated Substitutions',
                                                        'Normalized\nsingle point mutation density',
                                                        normalizedMutationData, sample,None, analysesType,
                                                        'yellowgreen', jobname, isFigureAugmentation,
                                                        sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
                                                        samplesWithAtLeast10KMutations2NumberofMutationsDict,
                                                        sample2NumberofIndelsDict,
                                                        signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,numberofSimulations)


        elif (analysesType == SIGNATUREBASED):
            for signature in signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
                # We check such file exists or not
                normalizedMutationData = readNormalizedMutationData(None, signature, jobname)

                if (normalizedMutationData is not None):
                    normalizedMutationData = normalizedMutationData.iloc[0].tolist()
                    # if not all([v == 0.0 for v in normalizedMutationData]):
                    # use all generator for all true check
                    if not all(v == 0.0 for v in normalizedMutationData):
                        # print('For %s: plot signature based replication time figure' % signature)
                        plotNormalizedMutationDensityFigureWithSimulations(signature,
                                                            'Normalized\nsingle point mutation density',
                                                            normalizedMutationData, None, signature, analysesType,
                                                            'yellowgreen', jobname, isFigureAugmentation,
                                                            sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
                                                            samplesWithAtLeast10KMutations2NumberofMutationsDict,
                                                            sample2NumberofIndelsDict,
                                                            signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,numberofSimulations)

                for sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
                    if signature in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict[sample]:
                        normalizedMutationData = readNormalizedMutationData(sample, signature, jobname)

                        if (normalizedMutationData is not None):
                            normalizedMutationData = normalizedMutationData.iloc[0].tolist()
                            # if not all([v == 0.0 for v in normalizedMutationData]):
                            # use all generator for all true check
                            if not all(v == 0.0 for v in normalizedMutationData):
                                # print('For %s: plot signature based replication time figure' % signature)
                                plotNormalizedMutationDensityFigureWithSimulations(signature,
                                                                    'Normalized\nsingle point mutation density',
                                                                    normalizedMutationData, sample,signature,
                                                                    analysesType, 'yellowgreen', jobname,
                                                                    isFigureAugmentation,
                                                                    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
                                                                    samplesWithAtLeast10KMutations2NumberofMutationsDict,
                                                                    sample2NumberofIndelsDict,
                                                                    signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,numberofSimulations)

    ##########################################################################################
    ##########################  Simulations ends  ############################################
    ##########################################################################################


    ##########################################################################################
    ##########################  Plot figures ends  ###########################################
    ##########################################################################################

##################################################################
