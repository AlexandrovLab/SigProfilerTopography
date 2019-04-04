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

plusOrMinus = 1000
windowSize = plusOrMinus*2+1


#############################################################################
##################### Read Average as Pandas Series #########################
#############################################################################
#Jobname has to be only jobname given in the argument
def readAverage(sample,signatureName,analyseType,outputDir,jobname):

    #####################################################
    if (analyseType == SIGNATUREBASED):
        filename = '%s_AverageNucleosomeSignalArray.txt' %(signatureName)
        averageFilePath = os.path.join(outputDir, jobname, DATA,NUCLEOSOMEOCCUPANCY, analyseType, filename)
    elif (analyseType== AGGREGATEDSUBSTITUTIONS or analyseType == AGGREGATEDINDELS):
        # new way
        filename = '%s_AverageNucleosomeSignalArray.txt' %(jobname)
        averageFilePath = os.path.join(outputDir, jobname, DATA,NUCLEOSOMEOCCUPANCY, analyseType, filename)
    #####################################################

    #####################################################
    if (analyseType == SAMPLEBASED_SIGNATUREBASED):
        filename = '%s_%s_AverageNucleosomeSignalArray.txt' % (signatureName, sample)
        averageFilePath = os.path.join(outputDir, jobname, DATA,NUCLEOSOMEOCCUPANCY, SIGNATUREBASED, filename)

    elif (analyseType == SAMPLEBASED_AGGREGATEDSUBSTITUTIONS):
        filename = '%s_%s_AverageNucleosomeSignalArray.txt' % (sample, jobname)
        averageFilePath = os.path.join(outputDir, jobname, DATA,NUCLEOSOMEOCCUPANCY, AGGREGATEDSUBSTITUTIONS, filename)


    elif (analyseType == SAMPLEBASED_AGGREGATEDINDELS):
        filename = '%s_%s_AverageNucleosomeSignalArray.txt' % (sample, jobname)
        averageFilePath = os.path.join(outputDir, jobname, DATA,NUCLEOSOMEOCCUPANCY, AGGREGATEDINDELS, filename)
    #####################################################

    return readAsNumpyArray(averageFilePath)
#############################################################################
##################### Read Average as Pandas Series #########################
#############################################################################




#############################################################################
def readAsNumpyArray(averageFilePath):
    if os.path.exists(averageFilePath):
        average = np.loadtxt(averageFilePath, dtype=float, delimiter='\t')
        average = average[plusOrMinus:(plusOrMinus*3)+1]
        return average
    else:
        return None
#############################################################################


#############################################################################
##################### Read Average for Simulations start ####################
#############################################################################
def readAverageForSimulations(sample,signature,analyseType,outputDir,jobname,numberofSimulations):

    listofAverages = []

    # Read the file w.r.t. the current folder

    if (analyseType == SIGNATUREBASED):
        # for signature based name may contain empty spaces
        #new way
        for i in range(1,numberofSimulations+1):
            simulationJobName = '%s_Sim%d' %(jobname,i)
            filename = '%s_AverageNucleosomeSignalArray.txt' %(signature)
            averageFilePath = os.path.join(outputDir, simulationJobName, DATA,NUCLEOSOMEOCCUPANCY, analyseType, filename)
            averageNucleosomeOccupancy = readAsNumpyArray(averageFilePath)
            if (averageNucleosomeOccupancy is not None):
                listofAverages.append(averageNucleosomeOccupancy)

    elif (analyseType== AGGREGATEDSUBSTITUTIONS or analyseType == AGGREGATEDINDELS):
        # i=0 for real/original data
        # i=1 for Sim1
        # ...
        # i=numberofSimulations for the last Simulations numberofSimulations
        for i in range(1,numberofSimulations+1):
            simulationJobName = '%s_Sim%d' %(jobname,i)
            filename = '%s_AverageNucleosomeSignalArray.txt' %(simulationJobName)
            averageFilePath = os.path.join(outputDir, simulationJobName, DATA,NUCLEOSOMEOCCUPANCY, analyseType, filename)
            averageNucleosomeOccupancy = readAsNumpyArray(averageFilePath)
            if (averageNucleosomeOccupancy is not None):
                listofAverages.append(averageNucleosomeOccupancy)

    #########################################################################
    if (analyseType == SAMPLEBASED_SIGNATUREBASED):
        for i in range(1,numberofSimulations+1):
            simulationJobName = '%s_Sim%d' %(jobname,i)
            filename = '%s_%s_AverageNucleosomeSignalArray.txt' %(signature,sample)
            averageFilePath = os.path.join(outputDir, simulationJobName, DATA,NUCLEOSOMEOCCUPANCY, SIGNATUREBASED, filename)
            averageNucleosomeOccupancy = readAsNumpyArray(averageFilePath)
            if (averageNucleosomeOccupancy is not None):
                listofAverages.append(averageNucleosomeOccupancy)

    elif (analyseType == SAMPLEBASED_AGGREGATEDSUBSTITUTIONS):
        for i in range(1, numberofSimulations + 1):
            simulationJobName = '%s_Sim%d' % (jobname, i)
            filename = '%s_%s_AverageNucleosomeSignalArray.txt' % (sample,simulationJobName)
            averageFilePath = os.path.join(outputDir,simulationJobName, DATA, NUCLEOSOMEOCCUPANCY,AGGREGATEDSUBSTITUTIONS, filename)
            averageNucleosomeOccupancy = readAsNumpyArray(averageFilePath)
            if (averageNucleosomeOccupancy is not None):
                listofAverages.append(averageNucleosomeOccupancy)

    elif (analyseType == SAMPLEBASED_AGGREGATEDINDELS):
        for i in range(1, numberofSimulations + 1):
            simulationJobName = '%s_Sim%d' % (jobname, i)
            filename = '%s_%s_AverageNucleosomeSignalArray.txt' % (sample,simulationJobName)
            averageFilePath = os.path.join(outputDir,simulationJobName, DATA, NUCLEOSOMEOCCUPANCY,AGGREGATEDINDELS, filename)
            averageNucleosomeOccupancy = readAsNumpyArray(averageFilePath)
            if (averageNucleosomeOccupancy is not None):
                listofAverages.append(averageNucleosomeOccupancy)
    #########################################################################


    return listofAverages
#############################################################################
##################### Read Average for Simulations end ######################
#############################################################################


#############################################################################
########################## Plot Figure starts  ##############################
#############################################################################
def plotSignatureBasedAverageNucleosomeOccupancyFigureWithSimulations(sample,signature,numberofMutations,xlabel,ylabel,outputDir,jobname, isFigureAugmentation,numberofSimulations,color):

    simulationsSignatureBasedMedians = None
    listofSimulationsSignatureBased = None

    if ((sample is not None) and (signature is not None)):
        realAverage = readAverage(sample, signature, SAMPLEBASED_SIGNATUREBASED, outputDir, jobname)
        figurename = '%s_%s_%d' % (signature, sample, numberofMutations)
        title = '%s_%s' % (signature, sample)
        if (numberofSimulations>0):
            listofSimulationsSignatureBased = readAverageForSimulations(sample, signature, SAMPLEBASED_SIGNATUREBASED, outputDir, jobname,numberofSimulations)
    else:
        realAverage = readAverage(None, signature, SIGNATUREBASED, outputDir, jobname)
        figurename = '%s_%d' % (signature, numberofMutations)
        title = '%s' % (signature)
        if (numberofSimulations>0):
            listofSimulationsSignatureBased = readAverageForSimulations(sample, signature, SIGNATUREBASED, outputDir, jobname,numberofSimulations)

    if ((realAverage is not None) and (pd.notna(realAverage).any(axis=0)) and (np.any(realAverage))):

        # print('sample:%s signature:%s' %(sample,signature))

        ###########################################################################
        #Get rid of divide by zero error.
        if ((listofSimulationsSignatureBased is not None) and listofSimulationsSignatureBased):
            # averageofSimulationsSignatureBased = sumofSimulationsSignatureBased / len(listofSimulationsSignatureBased)
            stackedSimulationsSignatureBased = np.vstack(listofSimulationsSignatureBased)
            (rows, cols) = stackedSimulationsSignatureBased.shape

            #Note that CI_indexLow makes the lower part of the nucleosome occupancy figures with simulations.
            CI_indexLow = int(round(.05 * rows, 0))+1
            CI_indexHigh = int(round(.95 * rows, 0))-1

            simulationsSignatureBasedLows = []
            simulationsSignatureBasedMedians = []
            simulationsSignatureBasedHighs = []

            for col in range(cols):
                colwise_array = stackedSimulationsSignatureBased[:, col]
                sorted_colwise_array = np.sort(colwise_array)
                simulationsSignatureBasedLows.append(sorted_colwise_array[CI_indexLow])
                # simulationsSignatureBasedMedians.append(np.median(sorted_colwise_array))
                simulationsSignatureBasedMedians.append(np.mean(sorted_colwise_array))
                simulationsSignatureBasedHighs.append(sorted_colwise_array[CI_indexHigh])
        ###########################################################################

        from matplotlib import rcParams
        rcParams.update({'figure.autolayout': True})

        # print('Plot signature based figure for %s starts' %signature)
        fig = plt.figure(figsize=(20,10),facecolor=None,dpi=300)
        plt.style.use('ggplot')

        # This code makes the background white.
        ax = plt.gca()
        ax.set_facecolor('white')
        # This code puts the edge line
        for edge_i in ['left', 'bottom','right', 'top']:
            ax.spines[edge_i].set_edgecolor("black")
            ax.spines[edge_i].set_linewidth(3)

        x = np.arange(-plusOrMinus ,plusOrMinus+1 ,1)
        # print('x.shape:%s' %x.shape)
        # plt.plot(x, average,'b-',label='test')


        listofLegends = []

        if (realAverage is not None):
            original = plt.plot(x, realAverage, color=color, label='Aggregated substitutions',linewidth=3)
            listofLegends.append(original[0])

        if (simulationsSignatureBasedMedians is not None):
            simulations = plt.plot(x, simulationsSignatureBasedMedians, color='gray', linestyle='--',  label='Average Simulations Aggregated substitutions', linewidth=3)
            listofLegends.append(simulations[0])
            plt.fill_between(x, np.array(simulationsSignatureBasedLows), np.array(simulationsSignatureBasedHighs),facecolor='lightblue')

        plt.legend(loc= 'lower left', handles=listofLegends, prop={'size': 24}, shadow=False, edgecolor='white', facecolor='white')

        #put the number of snps
        text = '%d subs' %(numberofMutations)
        plt.text(0.99, 0.99, text, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, fontsize=24)

        #Put vertical line at x=0
        # plt.axvline(x=0, ymin=0, ymax=1, color='gray', linestyle='--')
        plt.axvline(x=0, color='gray', linestyle='--')

        # This code provides the x and y tick marks and labels
        plt.xticks(np.arange(-1000, +1001, step=500), fontsize=30)

        #July 27, 2018
        plt.xlim((-1000, 1000))

        #Comment these lines to see the values out of range
        plt.yticks(np.arange(0.55, 1.4, step=0.1), fontsize=30)
        # This code overwrites the y tick labels
        ax.set_yticklabels(['0.55','','0.75','','0.95','','1.15','','1.35'])
        # This code provides some extra space
        plt.ylim((0.54, 1.4))

        # This code puts the tick marks
        plt.tick_params(axis='both', which='major', labelsize=30,width=3,length=10)
        plt.tick_params(axis='both', which='minor', labelsize=30,width=3,length=10)

        if (isFigureAugmentation):
            plt.title(jobname + ' ' + title, fontsize=40,fontweight='bold')
        else:
            plt.title(title, fontsize=40,fontweight='bold')

        # plt.xlabel(xlabel,fontsize=30,fontweight='bold')
        # plt.ylabel(ylabel,fontsize=30,fontweight='bold')

        plt.xlabel(xlabel,fontsize=32,fontweight='semibold')
        plt.ylabel(ylabel,fontsize=32,fontweight='semibold')

        filename = figurename.replace(' ', '') + '_NucleosomeOccupancy.png'

        #######################################################################
        # new code
        if (sample is None):
            figureFile = os.path.join(outputDir, jobname, FIGURE, ALL, NUCLEOSOMEOCCUPANCY, filename)
        else:
            os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, NUCLEOSOMEOCCUPANCY), exist_ok=True)
            figureFile = os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, NUCLEOSOMEOCCUPANCY, filename)
        #######################################################################

        fig.savefig(figureFile)
        plt.close(fig)
#############################################################################
########################## Plot Figure starts  ##############################
#############################################################################


#############################################################################
############## Plot AggregatedIndels For Simulations starts #################
#############################################################################
def plotAggregatedIndelsWithSimulations(xlabel,ylabel,sample,signature,analyseType,outputDir,jobname,isFigureAugmentation,numberofIndels,numberofSimulations):

    realAggregatedIndels = readAverage(sample, signature, analyseType, outputDir, jobname)
    simulationsAggregatedIndelsMedians = None
    listofSimulationsAggregatedIndels = None

    if (sample is None):
        figurename = '%s_AggregatedIndels.png' % (jobname)
        if (numberofSimulations>0):
            listofSimulationsAggregatedIndels =  readAverageForSimulations(sample,None, AGGREGATEDINDELS, outputDir, jobname,numberofSimulations)
    else:
        figurename = '%s_%d_AggregatedIndels_%s.png' % (sample,numberofIndels,jobname)
        if (numberofSimulations>0):
            listofSimulationsAggregatedIndels = readAverageForSimulations(sample, None, SAMPLEBASED_AGGREGATEDINDELS, outputDir, jobname, numberofSimulations)

    #####################################################################
    #####################################################################
    #####################################################################

    #Find Lows Medians Highs
    if (listofSimulationsAggregatedIndels is not None):
        stackedSimulationAggregatedIndels = np.vstack(listofSimulationsAggregatedIndels)
        (rows, cols) = stackedSimulationAggregatedIndels.shape

        # Note that CI_indexLow makes the lower part of the nucleosome occupancy figures with simulations.
        CI_indexLow = int(round(.05 * rows, 0))+1
        CI_indexHigh = int(round(.95 * rows, 0))-1

        simulationsAggregatedIndelsLows = []
        simulationsAggregatedIndelsMedians = []
        simulationsAggregatedIndelsHighs = []

        for col in range(cols):
            print('-------------    ')
            colwise_array = stackedSimulationAggregatedIndels[:, col]
            sorted_colwise_array = np.sort(colwise_array)
            simulationsAggregatedIndelsLows.append(sorted_colwise_array[CI_indexLow])
            # simulationsAggregatedIndelsMedians.append(np.median(sorted_colwise_array))
            simulationsAggregatedIndelsMedians.append(np.mean(sorted_colwise_array))
            simulationsAggregatedIndelsHighs.append(sorted_colwise_array[CI_indexHigh])

    #####################################################################
    #####################################################################
    #####################################################################

    #####################################################################
    #####################################################################
    #####################################################################
    print('Plot AggregatedIndels With Simulations figure starts')
    fig = plt.figure(figsize=(30, 10), facecolor=None)
    plt.style.use('ggplot')

    # This code makes the background white.
    ax = plt.gca()
    ax.set_facecolor('white')

    # This code puts the edge line
    for edge_i in ['left', 'bottom','right', 'top']:
        ax.spines[edge_i].set_edgecolor("black")
        ax.spines[edge_i].set_linewidth(3)

    x = np.arange(-plusOrMinus, plusOrMinus + 1, 1)
    # print('x.shape:%s' % x.shape)

    listofLegends = []

    if (realAggregatedIndels is not None):
        original = plt.plot(x, realAggregatedIndels, 'darkgreen', label='Aggregated indels',linewidth=3)
        listofLegends.append(original[0])
    if simulationsAggregatedIndelsMedians is not None:
        simulations = plt.plot(x, simulationsAggregatedIndelsMedians, color='gray', linestyle=':',label='Average Simulations Aggregated indels', linewidth=3)
        listofLegends.append(simulations[0])
        plt.fill_between(x,np.array(simulationsAggregatedIndelsLows),np.array(simulationsAggregatedIndelsHighs),facecolor='lightgreen')

    #test it
    # plt.legend(loc='lower left', prop={'size': 24},  shadow=False, edgecolor='white', facecolor ='white')
    plt.legend(loc= 'lower left',handles=listofLegends, prop={'size': 24}, shadow=False, edgecolor='white', facecolor='white')

    # put the number of snps
    text = '%d indels' % (numberofIndels)
    plt.text(0.99, 0.99, text, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, fontsize=24)

    #Put vertical line at x=0
    plt.axvline(x=0, color='gray', linestyle='--')

    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.tick_params(axis='both', which='minor', labelsize=24)

    # This code puts the tick marks
    plt.tick_params(axis='both', which='major', labelsize=24,width=3,length=10)
    plt.tick_params(axis='both', which='minor', labelsize=24,width=3,length=10)

    # This code provides the x and y tick marks and labels
    plt.xticks(np.arange(-1000, +1001, step=500), fontsize=30)
    plt.yticks(np.arange(0.7, 1.1, step=0.05), fontsize=30)

    # This code overwrites the y tick labels
    ax.set_yticklabels(['0.7','','0.8','','0.9','','1.0','','1.1'])

    #July 27, 2018
    plt.xlim((-1000, 1000))
    # This code provides some extra space
    plt.ylim((0.65,1.15))

    if (isFigureAugmentation):
        plt.title(jobname, fontsize=40)

    plt.xlabel(xlabel, fontsize=30)
    plt.ylabel(ylabel, fontsize=30)

    ######################################################################
    #new code
    if (sample is None):
        figureFile = os.path.join(outputDir,jobname, FIGURE,ALL,NUCLEOSOMEOCCUPANCY, figurename)
    else:
        os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample,NUCLEOSOMEOCCUPANCY), exist_ok=True)
        figureFile = os.path.join(outputDir, jobname, FIGURE,SAMPLES, sample,NUCLEOSOMEOCCUPANCY, figurename)
    ######################################################################

    fig.savefig(figureFile)
    print('Plot AggregatedIndels With Simulations figure ends')
    plt.close(fig)
    #####################################################################
    #####################################################################
    #####################################################################

#############################################################################
############## Plot AggregatedIndels For Simulations ends ###################
#############################################################################



#############################################################################
########### Plot AggregatedSubstitutions For Simulations starts #############
#############################################################################
def plotAggregatedSubstitutionsWithSimulations(xlabel,ylabel,sample,signature,analyseType,outputDir,jobname,isFigureAugmentation,sampleBasedNumberofMutations,numberofSimulations):

    simulationsAggregatedSubstitutionsMedians = None
    listofSimulationsAggregatedSubstitutions = None

    realAggregatedSubstitutions = readAverage(sample, signature, analyseType, outputDir, jobname)

    if (sample is None):
        filename = '%s_Aggregated_Substitutions.png' %(jobname)
        if (numberofSimulations>0):
            listofSimulationsAggregatedSubstitutions = readAverageForSimulations(None, None, AGGREGATEDSUBSTITUTIONS,outputDir,jobname, numberofSimulations)

    else:
        filename = '%s_%d_Aggregated_Substitutions_%s.png' %(sample,sampleBasedNumberofMutations,outputDir,jobname)
        if (numberofSimulations>0):
            listofSimulationsAggregatedSubstitutions = readAverageForSimulations(sample, None, SAMPLEBASED_AGGREGATEDSUBSTITUTIONS, outputDir,jobname,numberofSimulations)

    #####################################################################
    #####################################################################
    #####################################################################
    # print('sample: %s' %(sample))
    # print('len(listofSimulationsAggregatedSubstitutions): %d' %(len(listofSimulationsAggregatedSubstitutions)))

    # Take Average
    if (listofSimulationsAggregatedSubstitutions is not None):

        stackedSimulationAggregatedSubstitutions = np.vstack(listofSimulationsAggregatedSubstitutions)
        (rows, cols) = stackedSimulationAggregatedSubstitutions.shape

        # Note that CI_indexLow makes the lower part of the nucleosome occupancy figures with simulations.
        CI_indexLow = int(round(.05 * rows, 0))+1
        CI_indexHigh = int(round(.95 * rows, 0))-1

        simulationsAggregatedSubstitutionsLows = []
        simulationsAggregatedSubstitutionsMedians = []
        simulationsAggregatedSubstitutionsHighs = []

        for col in range(cols):
            colwise_array = stackedSimulationAggregatedSubstitutions[:, col]
            sorted_colwise_array = np.sort(colwise_array)
            simulationsAggregatedSubstitutionsLows.append(sorted_colwise_array[CI_indexLow])
            # simulationsAggregatedSubstitutionsMedians.append(np.median(sorted_colwise_array))
            simulationsAggregatedSubstitutionsMedians.append(np.mean(sorted_colwise_array))
            simulationsAggregatedSubstitutionsHighs.append(sorted_colwise_array[CI_indexHigh])
    #####################################################################
    #####################################################################
    #####################################################################


    #####################################################################
    #####################################################################
    #####################################################################
    fig = plt.figure(figsize=(30, 10), facecolor=None)
    plt.style.use('ggplot')

    # This code makes the background white.
    ax = plt.gca()
    ax.set_facecolor('white')

    # This code puts the edge line
    for edge_i in ['left', 'bottom','right', 'top']:
        ax.spines[edge_i].set_edgecolor("black")
        ax.spines[edge_i].set_linewidth(3)

    x = np.arange(-plusOrMinus, plusOrMinus + 1, 1)
    # print('x.shape:%s' % x.shape)

    listofLegends = []

    if (realAggregatedSubstitutions is not None):
        original = plt.plot(x, realAggregatedSubstitutions, 'royalblue', label='Aggregated substitutions',linewidth=3)
        listofLegends.append(original[0])
    if (simulationsAggregatedSubstitutionsMedians is not None):
        simulations = plt.plot(x, simulationsAggregatedSubstitutionsMedians, color='gray',linestyle='--', label='Average Simulations Aggregated substitutions',linewidth=3)
        listofLegends.append(simulations[0])
        plt.fill_between(x,np.array(simulationsAggregatedSubstitutionsLows),np.array(simulationsAggregatedSubstitutionsHighs),facecolor='lightblue')

    # plt.legend(loc='lower left', prop={'size': 24},  shadow=False, edgecolor='white', facecolor ='white')
    plt.legend(loc= 'lower left',handles=listofLegends, prop={'size': 24}, shadow=False, edgecolor='white', facecolor='white')

    # put the number of snps
    text = '%d subs' % (sampleBasedNumberofMutations)
    plt.text(0.99, 0.99, text, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, fontsize=24)

    #Put vertical line at x=0
    plt.axvline(x=0, color='gray', linestyle='--')

    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.tick_params(axis='both', which='minor', labelsize=24)

    # This code puts the tick marks
    plt.tick_params(axis='both', which='major', labelsize=24,width=3,length=10)
    plt.tick_params(axis='both', which='minor', labelsize=24,width=3,length=10)

    # This code provides the x and y tick marks and labels
    plt.xticks(np.arange(-1000, +1001, step=500), fontsize=30)
    plt.yticks(np.arange(0.7, 1.1, step=0.05), fontsize=30)

    # This code overwrites the y tick labels
    ax.set_yticklabels(['0.7','','0.8','','0.9','','1.0','','1.1'])

    #July 27, 2018
    plt.xlim((-1000, 1000))
    # This code provides some extra space
    plt.ylim((0.65,1.15))

    if (isFigureAugmentation):
        plt.title(jobname, fontsize=40)

    plt.xlabel(xlabel, fontsize=30)
    plt.ylabel(ylabel, fontsize=30)

    #######################################################################
    # new code
    if (sample is None):
        figureFile = os.path.join(outputDir, jobname, FIGURE, ALL, NUCLEOSOMEOCCUPANCY, filename)
    else:
        os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample,NUCLEOSOMEOCCUPANCY), exist_ok=True)
        figureFile = os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, NUCLEOSOMEOCCUPANCY, filename)
    #######################################################################

    fig.savefig(figureFile)
    plt.close(fig)
    #####################################################################
    #####################################################################
    #####################################################################

#############################################################################
########### Plot AggregatedSubstitutions For Simulations ends ###############
#############################################################################


#############################################################################
############################ Plot Figure ####################################
#############################################################################
def plotAggregatedSubstitutionsandAggregatedIndelsWithSimulations(xlabel,ylabel,sample,outputDir,jobname,isFigureAugmentation,numberofSPMs,numberofIndels,numberofSimulations):

    simulationsAggregatedSubstitutionsMedians = None
    simulationsAggregatedIndelsMedians = None

    listofSimulationsAggregatedIndels = None
    listofSimulationsAggregatedSubstitutions = None

    #######################################################################################################################
    if (sample is None):
        filename = '%s_Aggregated_Substitutions_Indels.png' % (jobname)
        realAggregatedIndels = readAverage(None,None, AGGREGATEDINDELS, outputDir,jobname)
        realAggregatedSubstitutions = readAverage(None,None,AGGREGATEDSUBSTITUTIONS,outputDir,jobname)

        if (numberofSimulations>0):
            listofSimulationsAggregatedIndels = readAverageForSimulations(None, None, AGGREGATEDINDELS, outputDir, jobname,numberofSimulations)
            listofSimulationsAggregatedSubstitutions =  readAverageForSimulations(None,None,AGGREGATEDSUBSTITUTIONS,outputDir,jobname,numberofSimulations)
    #######################################################################################################################


    #######################################################################################################################
    else:
        filename = '%s_%s_Aggregated_Substitutions_%d_Indels_%d.png' % (sample, jobname, numberofSPMs, numberofIndels)
        realAggregatedIndels = readAverage(sample,None, SAMPLEBASED_AGGREGATEDINDELS, outputDir, jobname)
        realAggregatedSubstitutions = readAverage(sample,None,SAMPLEBASED_AGGREGATEDSUBSTITUTIONS, outputDir, jobname)

        if (numberofSimulations>0):
            listofSimulationsAggregatedIndels = readAverageForSimulations(sample, None, SAMPLEBASED_AGGREGATEDINDELS, outputDir, jobname, numberofSimulations)
            listofSimulationsAggregatedSubstitutions = readAverageForSimulations(sample, None, SAMPLEBASED_AGGREGATEDSUBSTITUTIONS, outputDir, jobname,numberofSimulations)
    #######################################################################################################################


    #####################################################################
    #####################################################################
    #####################################################################

    #Take Median
    if ((listofSimulationsAggregatedIndels is not None) and listofSimulationsAggregatedIndels):

        stackedSimulationAggregatedIndels = np.vstack(listofSimulationsAggregatedIndels)
        (rows, cols) = stackedSimulationAggregatedIndels.shape

        # Note that CI_indexLow makes the lower part of the nucleosome occupancy figures with simulations.
        CI_indexLow = int(round(.05 * rows, 0))+1
        CI_indexHigh = int(round(.95 * rows, 0))-1

        simulationsAggregatedIndelsLows = []
        simulationsAggregatedIndelsHighs = []
        simulationsAggregatedIndelsMedians = []

        for col in range(cols):
            colwise_array = stackedSimulationAggregatedIndels[:, col]
            sorted_colwise_array = np.sort(colwise_array)

            simulationsAggregatedIndelsLows.append(sorted_colwise_array[CI_indexLow])
            #Take median
            # simulationsAggregatedIndelsMedians.append(np.median(sorted_colwise_array))
            #Take mean
            simulationsAggregatedIndelsMedians.append(np.mean(sorted_colwise_array))
            simulationsAggregatedIndelsHighs.append(sorted_colwise_array[CI_indexHigh])
    #####################################################################
    #####################################################################
    #####################################################################

    #####################################################################
    #####################################################################
    #####################################################################
    # Take Average
    if ((listofSimulationsAggregatedSubstitutions is not None) and listofSimulationsAggregatedSubstitutions):

        stackedSimulationAggregatedSubstitutions = np.vstack(listofSimulationsAggregatedSubstitutions)

        (rows, cols) = stackedSimulationAggregatedSubstitutions.shape

        CI_indexLow = int(round(.05 * rows, 0))+1
        CI_indexHigh = int(round(.95 * rows, 0))-1

        simulationsAggregatedSubstitutionsLows = []
        simulationsAggregatedSubstitutionsMedians = []
        simulationsAggregatedSubstitutionsHighs = []

        for col in range(cols):
            colwise_array = stackedSimulationAggregatedSubstitutions[:, col]
            sorted_colwise_array = np.sort(colwise_array)
            simulationsAggregatedSubstitutionsLows.append(sorted_colwise_array[CI_indexLow])
            # simulationsAggregatedSubstitutionsMedians.append(np.median(sorted_colwise_array))
            simulationsAggregatedSubstitutionsMedians.append(np.mean(sorted_colwise_array))
            simulationsAggregatedSubstitutionsHighs.append(sorted_colwise_array[CI_indexHigh])
    #####################################################################
    #####################################################################
    #####################################################################


    #####################################################################
    #####################################################################
    #####################################################################
    fig = plt.figure(figsize=(30, 10), facecolor=None)
    plt.style.use('ggplot')

    # This code makes the background white.
    ax = plt.gca()
    ax.set_facecolor('white')

    # This code puts the edge line
    for edge_i in ['left', 'bottom','right', 'top']:
        ax.spines[edge_i].set_edgecolor("black")
        ax.spines[edge_i].set_linewidth(3)

    x = np.arange(-plusOrMinus, plusOrMinus + 1, 1)
    # print('x.shape:%s' % x.shape)

    listofLegends = []

    if (realAggregatedSubstitutions is not None):
        aggSubs = plt.plot(x, realAggregatedSubstitutions, 'royalblue', label='Aggregated substitutions',linewidth=3,zorder=10)
        listofLegends.append(aggSubs[0])
    if (simulationsAggregatedSubstitutionsMedians is not None):
        simsAggSubs = plt.plot(x, simulationsAggregatedSubstitutionsMedians, color='gray',linestyle='--', label='Average Simulations Aggregated substitutions',linewidth=3,zorder=10)
        listofLegends.append(simsAggSubs[0])
        plt.fill_between(x,np.array(simulationsAggregatedSubstitutionsLows),np.array(simulationsAggregatedSubstitutionsHighs),facecolor='lightblue',zorder=10)

    if (realAggregatedIndels is not None):
        aggIndels = plt.plot(x, realAggregatedIndels, 'darkgreen', label='Aggregated indels',linewidth=3,zorder=10)
        listofLegends.append(aggIndels[0])
    if simulationsAggregatedIndelsMedians is not None:
        simsAggIndels = plt.plot(x, simulationsAggregatedIndelsMedians, color='gray', linestyle=':',label='Average Simulations Aggregated indels', linewidth=3,zorder=10)
        listofLegends.append(simsAggIndels[0])
        plt.fill_between(x,np.array(simulationsAggregatedIndelsLows),np.array(simulationsAggregatedIndelsHighs),facecolor='lightgreen',zorder=5)


    # old code
    # plt.legend(loc='lower left', prop={'size': 24},  shadow=False, edgecolor='white', facecolor ='white')
    plt.legend(loc= 'lower left',handles = listofLegends, prop={'size': 24}, shadow=False, edgecolor='white', facecolor ='white')

    # put the number of snps and indels
    if sample is not None:
        text = '%d subs, %d indels' %(numberofSPMs,numberofIndels)
        plt.text(0.99, 0.99, text, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, fontsize=24)

    #Put vertical line at x=0
    plt.axvline(x=0, color='gray', linestyle='--')

    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.tick_params(axis='both', which='minor', labelsize=24)

    # This code puts the tick marks
    plt.tick_params(axis='both', which='major', labelsize=24,width=3,length=10)
    plt.tick_params(axis='both', which='minor', labelsize=24,width=3,length=10)

    # This code provides the x and y tick marks and labels
    plt.xticks(np.arange(-1000, +1001, step=500), fontsize=30)
    plt.yticks(np.arange(0.7, 1.1, step=0.05), fontsize=30)

    # This code overwrites the y tick labels
    ax.set_yticklabels(['0.7','','0.8','','0.9','','1.0','','1.1'])

    #July 27, 2018
    plt.xlim((-1000, 1000))
    # This code provides some extra space
    plt.ylim((0.65,1.15))

    if (sample is not None):
        plt.title(sample, fontsize=40, fontweight='bold')
    else:
        plt.title(jobname, fontsize=40, fontweight='bold')


    plt.xlabel(xlabel, fontsize=30)
    plt.ylabel(ylabel, fontsize=30)

    ######################################################################################
    if (sample is None):
        figureFile = os.path.join(outputDir,jobname,FIGURE,ALL,NUCLEOSOMEOCCUPANCY,filename)
    else:
        os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample,NUCLEOSOMEOCCUPANCY), exist_ok=True)
        figureFile = os.path.join(outputDir,jobname,FIGURE,SAMPLES,sample,NUCLEOSOMEOCCUPANCY,filename)
    ######################################################################################

    fig.savefig(figureFile)
    plt.close(fig)
    #####################################################################
    #####################################################################
    #####################################################################

#############################################################################
############################ Plot Figure ####################################
#############################################################################


#########################################################
def checkValidness(analsesType,outputDir,jobname):
    #Check whether directory exists and there are files under it.

    data_file_path = os.path.join(outputDir,jobname,DATA,NUCLEOSOMEOCCUPANCY,analsesType)

    #If directory exists and if there are files under it.
    if (os.path.exists(data_file_path)):
        filenames = [f for f in os.listdir(data_file_path) if os.path.isfile(os.path.join(data_file_path, f))]
        if (len(filenames)>0):
            return True
        else:
            return False
    else:
        return False
#########################################################


#########################################################
def nucleosomeOccupancyAverageSignalFigures(outputDir,jobname,figureAugmentation,numberofSimulations):

    #######################################################################################################################
    os.makedirs(os.path.join(outputDir, jobname, FIGURE, ALL, NUCLEOSOMEOCCUPANCY), exist_ok=True)
    os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES), exist_ok=True)
    #######################################################################################################################

    isFigureAugmentation = False
    if (figureAugmentation=='augmentation'):
        isFigureAugmentation = True

    ############## Read necessary dictionaries starts ########################################
    sample2NumberofSubsDict = getSample2NumberofSubsDict(outputDir,jobname)
    sample2NumberofIndelsDict = getSample2NumberofIndelsDict(outputDir,jobname)

    sample2SubsSignature2NumberofMutationsDict = getSample2SubsSignature2NumberofMutationsDict(outputDir,jobname)
    sample2IndelsSignature2NumberofMutationsDict = getSample2IndelsSignature2NumberofMutationsDict(outputDir, jobname)

    subsSignature2NumberofMutationsDict = getSubsSignature2NumberofMutationsDict(outputDir,jobname)
    indelsSignature2NumberofMutationsDict = getIndelsSignature2NumberofMutationsDict(outputDir,jobname)
    ############## Read necessary dictionaries ends ##########################################


    # Plot All Together : Aggregated Substitutions and Aggregated Indels
    # Or plot one of them:  Aggregated Substitutions or Aggregated Indels
    if checkValidness(AGGREGATEDSUBSTITUTIONS,outputDir,jobname) and checkValidness(AGGREGATEDINDELS,outputDir,jobname):

        ##############################################################
        plotAggregatedSubstitutionsandAggregatedIndelsWithSimulations('Interval around variant (bp)','Average nucleosome signal',
                                                                          None,outputDir,jobname,isFigureAugmentation,
                                                                          0,0,numberofSimulations)
        ##############################################################


        ##############################################################
        #Arrays are filled w.r.t. sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict
        samplesfromSubs  = sample2NumberofSubsDict.keys()
        samplesfromIndels = sample2NumberofIndelsDict.keys()

        for sample in (samplesfromSubs | samplesfromIndels):
            numberofSubs = 0
            numberofIndels = 0
            if sample in sample2NumberofSubsDict:
                numberofSubs = sample2NumberofSubsDict[sample]
            if sample in sample2NumberofIndelsDict:
                numberofIndels = sample2NumberofIndelsDict[sample]

            plotAggregatedSubstitutionsandAggregatedIndelsWithSimulations('Interval around variant (bp)','Average nucleosome signal',
                                                                        sample,outputDir,jobname, isFigureAugmentation,numberofSubs, numberofIndels,numberofSimulations)

        ##############################################################

    ######################################################################################
    elif (checkValidness(AGGREGATEDSUBSTITUTIONS,outputDir,jobname)):
        sampleBasedNumberofMutations = 0
        plotAggregatedSubstitutionsWithSimulations('Interval around single point mutation (bp)', 'Average nucleosome signal',
                                                    None, None, AGGREGATEDSUBSTITUTIONS,
                                                    outputDir,jobname,isFigureAugmentation,sampleBasedNumberofMutations,numberofSimulations)

        for sample in sample2SubsSignature2NumberofMutationsDict:
        # for sample in samplesWithAtLeast10KMutations2NumberofMutationsDict:
            sampleBasedNumberofMutations = sample2NumberofSubsDict[sample]
            plotAggregatedSubstitutionsWithSimulations('Interval around single point mutation (bp)', 'Average nucleosome signal',
                                                        sample, None, SAMPLEBASED_AGGREGATEDSUBSTITUTIONS,
                                                       outputDir,jobname, isFigureAugmentation, sampleBasedNumberofMutations,numberofSimulations)
    ######################################################################################



    ######################################################################################
    elif (checkValidness(AGGREGATEDINDELS,outputDir,jobname)):
        numberofIndels = 0
        plotAggregatedIndelsWithSimulations('Interval around variant (bp)', 'Average nucleosome signal',
                                            None, None,
                                            AGGREGATEDINDELS, outputDir, jobname, isFigureAugmentation,numberofIndels,numberofSimulations)

        for sample in sample2SubsSignature2NumberofMutationsDict:
        # for sample in samplesWithAtLeast10KMutations2NumberofMutationsDict:
            if sample in sample2NumberofIndelsDict:
                numberofIndels = sample2NumberofIndelsDict[sample]
                plotAggregatedIndelsWithSimulations('Interval around variant (bp)', 'Average nucleosome signal',
                                                sample, None,
                                                SAMPLEBASED_AGGREGATEDINDELS, outputDir, jobname, isFigureAugmentation,numberofIndels,numberofSimulations)
    ######################################################################################


    #############################################################################################################################################
    #Plot Signature Based
    #ncomms11383 Fig3b signature based average nucleosome occupancy figures
    if checkValidness(SIGNATUREBASED,outputDir,jobname):
        #Subs Signatures
        for subsSignature in subsSignature2NumberofMutationsDict:
            subsSignatureBasedNumberofMutations = subsSignature2NumberofMutationsDict[subsSignature]
            plotSignatureBasedAverageNucleosomeOccupancyFigureWithSimulations(None,subsSignature,subsSignatureBasedNumberofMutations,
                                                                                  'Interval around single point mutation (bp)','Average nucleosome signal',
                                                                              outputDir, jobname,isFigureAugmentation,numberofSimulations,'royalblue')


        #Indels Signatures
        for indelSignature in indelsSignature2NumberofMutationsDict:
            indelsSignatureBasedNumberofMutations = indelsSignature2NumberofMutationsDict[indelSignature]
            plotSignatureBasedAverageNucleosomeOccupancyFigureWithSimulations(None,indelSignature,indelsSignatureBasedNumberofMutations,'Interval around variant (bp)','Average nucleosome signal',
                                                                              outputDir, jobname,isFigureAugmentation,numberofSimulations,'darkgreen')


        # SampleBased Subs SignatureBased Nucleosome Occupancy Figures
        for sample in sample2SubsSignature2NumberofMutationsDict:
            for subsSignature in sample2SubsSignature2NumberofMutationsDict[sample]:
                sampleBasedSignatureBasedNumberofMutations = sample2SubsSignature2NumberofMutationsDict[sample][subsSignature]
                plotSignatureBasedAverageNucleosomeOccupancyFigureWithSimulations(sample,subsSignature,sampleBasedSignatureBasedNumberofMutations,
                                                                                      'Interval around single point mutation (bp)','Average nucleosome signal',
                                                                                  outputDir, jobname, isFigureAugmentation,numberofSimulations,'royalblue')


        # SampleBased Indels SignatureBased Nucleosome Occupancy Figures
        for sample in sample2IndelsSignature2NumberofMutationsDict:
            for indelsSignature in sample2IndelsSignature2NumberofMutationsDict[sample]:
                sampleBasedSignatureBasedNumberofMutations = sample2IndelsSignature2NumberofMutationsDict[sample][indelsSignature]
                plotSignatureBasedAverageNucleosomeOccupancyFigureWithSimulations(sample,indelsSignature,sampleBasedSignatureBasedNumberofMutations,
                                                                                      'Interval around variant (bp)','Average nucleosome signal',
                                                                                  outputDir, jobname, isFigureAugmentation,numberofSimulations,'darkgreen')
    #############################################################################################################################################