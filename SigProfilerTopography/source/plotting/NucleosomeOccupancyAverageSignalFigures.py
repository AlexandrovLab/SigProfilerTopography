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
    elif ((analyseType== AGGREGATEDSUBSTITUTIONS) or (analyseType == AGGREGATEDINDELS) or (analyseType == AGGREGATEDDINUCS)):
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

    elif (analyseType == SAMPLEBASED_AGGREGATEDDINUCS):
        filename = '%s_%s_AverageNucleosomeSignalArray.txt' % (sample, jobname)
        averageFilePath = os.path.join(outputDir, jobname, DATA,NUCLEOSOMEOCCUPANCY, AGGREGATEDDINUCS, filename)

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
        for simNum in range(1,numberofSimulations+1):
            filename = '%s_sim%d_AverageNucleosomeSignalArray.txt' %(signature,simNum)
            averageFilePath = os.path.join(outputDir,jobname,DATA,NUCLEOSOMEOCCUPANCY,analyseType,filename)
            averageNucleosomeOccupancy = readAsNumpyArray(averageFilePath)
            if (averageNucleosomeOccupancy is not None):
                listofAverages.append(averageNucleosomeOccupancy)

    elif (analyseType== AGGREGATEDSUBSTITUTIONS or analyseType == AGGREGATEDINDELS or analyseType == AGGREGATEDDINUCS):
        # i=0 for real/original data
        # i=1 for Sim1
        # ...
        # i=numberofSimulations for the last Simulations numberofSimulations
        for simNum in range(1,numberofSimulations+1):
            filename = '%s_sim%d_AverageNucleosomeSignalArray.txt' %(jobname,simNum)
            averageFilePath = os.path.join(outputDir, jobname, DATA,NUCLEOSOMEOCCUPANCY, analyseType, filename)
            averageNucleosomeOccupancy = readAsNumpyArray(averageFilePath)
            if (averageNucleosomeOccupancy is not None):
                listofAverages.append(averageNucleosomeOccupancy)

    #########################################################################
    if (analyseType == SAMPLEBASED_SIGNATUREBASED):
        for simNum in range(1,numberofSimulations+1):
            filename = '%s_%s_sim%d_AverageNucleosomeSignalArray.txt' %(signature,sample,simNum)
            averageFilePath = os.path.join(outputDir, jobname, DATA,NUCLEOSOMEOCCUPANCY, SIGNATUREBASED, filename)
            averageNucleosomeOccupancy = readAsNumpyArray(averageFilePath)
            if (averageNucleosomeOccupancy is not None):
                listofAverages.append(averageNucleosomeOccupancy)

    elif (analyseType == SAMPLEBASED_AGGREGATEDSUBSTITUTIONS):
        for simNum in range(1, numberofSimulations + 1):
            filename = '%s_%s_sim%d_AverageNucleosomeSignalArray.txt' % (sample,jobname,simNum)
            averageFilePath = os.path.join(outputDir,jobname, DATA, NUCLEOSOMEOCCUPANCY,AGGREGATEDSUBSTITUTIONS, filename)
            averageNucleosomeOccupancy = readAsNumpyArray(averageFilePath)
            if (averageNucleosomeOccupancy is not None):
                listofAverages.append(averageNucleosomeOccupancy)

    elif (analyseType == SAMPLEBASED_AGGREGATEDINDELS):
        for simNum in range(1, numberofSimulations + 1):
            filename = '%s_%s_sim%d_AverageNucleosomeSignalArray.txt' % (sample,jobname,simNum)
            averageFilePath = os.path.join(outputDir,jobname, DATA, NUCLEOSOMEOCCUPANCY,AGGREGATEDINDELS, filename)
            averageNucleosomeOccupancy = readAsNumpyArray(averageFilePath)
            if (averageNucleosomeOccupancy is not None):
                listofAverages.append(averageNucleosomeOccupancy)

    elif (analyseType == SAMPLEBASED_AGGREGATEDDINUCS):
        for simNum in range(1, numberofSimulations + 1):
            filename = '%s_%s_sim%d_AverageNucleosomeSignalArray.txt' % (sample,jobname,simNum)
            averageFilePath = os.path.join(outputDir,jobname, DATA, NUCLEOSOMEOCCUPANCY,AGGREGATEDDINUCS, filename)
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
def plotAllSamplesPooledAndSampleBasedSignaturesFiguresInOneFigure(signature2NumberofMutationsDict,sample2Signature2NumberofMutationsDict,outputDir,jobname,color,xlabel,ylabel):

    for signature in signature2NumberofMutationsDict:
        label2NumpyArrayDict = {}
        signatureBasedNumberofMutations = signature2NumberofMutationsDict[signature]
        realAverage = readAverage(None, signature, SIGNATUREBASED, outputDir, jobname)
        label2NumpyArrayDict[signature] = realAverage

        for sample in sample2Signature2NumberofMutationsDict:
            if signature in sample2Signature2NumberofMutationsDict[sample]:
                realAverage = readAverage(sample,signature,SAMPLEBASED_SIGNATUREBASED,outputDir,jobname)
                label = '%s_%s' %(sample,signature)
                label2NumpyArrayDict[label] = realAverage

        ####################################### plottoing starts #######################################
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
        for label in label2NumpyArrayDict:
            realAverage = label2NumpyArrayDict[label]
            if (realAverage is not None):
                if (label==signature):
                    original = plt.plot(x,realAverage,color=color,label=label,linewidth=3,zorder=10)
                else:
                    original = plt.plot(x,realAverage,color='gray',label=label,linewidth=3,zorder=5)
                listofLegends.append(original[0])


        # plt.legend(loc= 'lower left', handles=listofLegends, prop={'size': 12}, shadow=False, edgecolor='white', facecolor='white')
        plt.legend(loc='lower left', handles=listofLegends, prop={'size': 12}, shadow=False, edgecolor='white',facecolor='white', ncol=5, fancybox=True)


        # text = '%d subs' %(numberofMutations)
        # plt.text(0.99, 0.99, text, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, fontsize=24)

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

        title = signature
        plt.title(title, fontsize=40,fontweight='bold')

        plt.xlabel(xlabel,fontsize=32,fontweight='semibold')
        plt.ylabel(ylabel,fontsize=32,fontweight='semibold')

        filename = '%s_AllInOne_NucleosomeOccupancy.png' %(signature)
        figureFile = os.path.join(outputDir, jobname, FIGURE, ALL, NUCLEOSOMEOCCUPANCY, filename)

        fig.savefig(figureFile)
        plt.close(fig)
        ####################################### plottoing ends #######################################


#############################################################################
########################## Plot Figure ends  ################################
#############################################################################


#############################################################################
########################## Plot Figure starts  ##############################
#############################################################################
def plotSignatureBasedAverageNucleosomeOccupancyFigureWithSimulations(sample,signature,numberofMutations,xlabel,ylabel,label,text,outputDir,jobname, isFigureAugmentation,numberofSimulations,color,fillcolor):

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
            original = plt.plot(x, realAverage, color=color, label=label,linewidth=3)
            listofLegends.append(original[0])

        if (simulationsSignatureBasedMedians is not None):
            label = 'Average Simulations %s' %(label)
            simulations = plt.plot(x, simulationsSignatureBasedMedians, color='gray', linestyle='--',  label=label, linewidth=3)
            listofLegends.append(simulations[0])
            plt.fill_between(x, np.array(simulationsSignatureBasedLows), np.array(simulationsSignatureBasedHighs),facecolor=fillcolor)

        plt.legend(loc= 'lower left', handles=listofLegends, prop={'size': 24}, shadow=False, edgecolor='white', facecolor='white')

        #put the number of snps
        text = '%d %s' %(numberofMutations,text)
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
def takeAverage(listofSimulationsAggregatedMutations):
    simulationsAggregatedMutationsLows = None
    simulationsAggregatedMutationsHighs = None
    simulationsAggregatedMutationsMedians = None

    if ((listofSimulationsAggregatedMutations is not None) and listofSimulationsAggregatedMutations):
        stackedSimulationAggregatedMutations = np.vstack(listofSimulationsAggregatedMutations)
        (rows, cols) = stackedSimulationAggregatedMutations.shape

        # Note that CI_indexLow makes the lower part of the nucleosome occupancy figures with simulations.
        CI_indexLow = int(round(.05 * rows, 0))+1
        CI_indexHigh = int(round(.95 * rows, 0))-1

        simulationsAggregatedMutationsLows = []
        simulationsAggregatedMutationsHighs = []
        simulationsAggregatedMutationsMedians = []

        for col in range(cols):
            colwise_array = stackedSimulationAggregatedMutations[:, col]
            sorted_colwise_array = np.sort(colwise_array)

            simulationsAggregatedMutationsLows.append(sorted_colwise_array[CI_indexLow])
            #Take median
            # simulationsAggregatedIndelsMedians.append(np.median(sorted_colwise_array))
            #Take mean
            simulationsAggregatedMutationsMedians.append(np.mean(sorted_colwise_array))
            simulationsAggregatedMutationsHighs.append(sorted_colwise_array[CI_indexHigh])

    return  simulationsAggregatedMutationsLows,simulationsAggregatedMutationsMedians,simulationsAggregatedMutationsHighs
#############################################################################


#############################################################################
############################ Plot Figure ####################################
#############################################################################
def plotAllMutationsPooledWithSimulations(xlabel,ylabel,sample,outputDir,jobname,numberofSPMs,numberofIndels,numberofDinucs,numberofSimulations,mutationTypes):
    realAggregatedSubstitutions = None
    realAggregatedIndels = None
    realAggregatedDinucs = None

    listofSimulationsAggregatedSubstitutions = None
    listofSimulationsAggregatedIndels = None
    listofSimulationsAggregatedDinucs = None

    #######################################################################################################################
    if (sample is None):
        filename = 'Aggregated_All_Mutations_NucleosomeOccupancy.png'
        if (SUBS in mutationTypes):
            realAggregatedSubstitutions = readAverage(None,None,AGGREGATEDSUBSTITUTIONS,outputDir,jobname)
        if (INDELS in mutationTypes):
            realAggregatedIndels = readAverage(None,None, AGGREGATEDINDELS, outputDir,jobname)
        if (DINUCS in mutationTypes):
            realAggregatedDinucs = readAverage(None,None,AGGREGATEDDINUCS,outputDir,jobname)

        if (numberofSimulations>0):
            if (SUBS in mutationTypes):
                listofSimulationsAggregatedSubstitutions = readAverageForSimulations(None,None,AGGREGATEDSUBSTITUTIONS,outputDir,jobname,numberofSimulations)
            if (INDELS in mutationTypes):
                listofSimulationsAggregatedIndels = readAverageForSimulations(None,None,AGGREGATEDINDELS,outputDir,jobname,numberofSimulations)
            if (DINUCS in mutationTypes):
                listofSimulationsAggregatedDinucs = readAverageForSimulations(None,None,AGGREGATEDDINUCS,outputDir,jobname,numberofSimulations)
    #######################################################################################################################


    #######################################################################################################################
    else:
        # filename = '%s_%s_Aggregated_Substitutions_%d_Indels_%d.png' % (sample, jobname, numberofSPMs, numberofIndels)
        filename = '%s_Aggregated_All_Mutations_NucleosomeOccupancy.png' % (sample)
        if (SUBS in mutationTypes):
            realAggregatedSubstitutions = readAverage(sample,None,SAMPLEBASED_AGGREGATEDSUBSTITUTIONS,outputDir,jobname)
        if (INDELS in mutationTypes):
            realAggregatedIndels = readAverage(sample,None,SAMPLEBASED_AGGREGATEDINDELS,outputDir,jobname)
        if (DINUCS in mutationTypes):
            realAggregatedDinucs = readAverage(sample,None,SAMPLEBASED_AGGREGATEDDINUCS,outputDir,jobname)

        if (numberofSimulations>0):
            if (SUBS in mutationTypes):
                listofSimulationsAggregatedSubstitutions = readAverageForSimulations(sample, None, SAMPLEBASED_AGGREGATEDSUBSTITUTIONS, outputDir,jobname,numberofSimulations)
            if (INDELS in mutationTypes):
                listofSimulationsAggregatedIndels = readAverageForSimulations(sample, None,SAMPLEBASED_AGGREGATEDINDELS, outputDir,jobname,numberofSimulations)
            if (DINUCS in mutationTypes):
                listofSimulationsAggregatedDinucs = readAverageForSimulations(sample, None, SAMPLEBASED_AGGREGATEDDINUCS, outputDir,jobname,numberofSimulations)
    #######################################################################################################################


    #####################################################################
    #Take Median
    simulationsAggregatedSubstitutionsLows,simulationsAggregatedSubstitutionsMedians,simulationsAggregatedSubstitutionsHighs = takeAverage(listofSimulationsAggregatedSubstitutions)
    simulationsAggregatedIndelsLows,simulationsAggregatedIndelsMedians,simulationsAggregatedIndelsHighs = takeAverage(listofSimulationsAggregatedIndels)
    simulationsAggregatedDinucsLows,simulationsAggregatedDinucsMedians,simulationsAggregatedDinucsHighs = takeAverage(listofSimulationsAggregatedDinucs)
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
        simsAggSubs = plt.plot(x, simulationsAggregatedSubstitutionsMedians, color='gray',linestyle='--', label='Average Simulations Aggregated Substitutions',linewidth=3,zorder=10)
        listofLegends.append(simsAggSubs[0])
        plt.fill_between(x,np.array(simulationsAggregatedSubstitutionsLows),np.array(simulationsAggregatedSubstitutionsHighs),facecolor='lightblue',zorder=10)

    if (realAggregatedIndels is not None):
        aggIndels = plt.plot(x, realAggregatedIndels, 'darkgreen', label='Aggregated indels',linewidth=3,zorder=10)
        listofLegends.append(aggIndels[0])
    if simulationsAggregatedIndelsMedians is not None:
        simsAggIndels = plt.plot(x, simulationsAggregatedIndelsMedians, color='gray', linestyle=':',label='Average Simulations Aggregated Indels', linewidth=3,zorder=10)
        listofLegends.append(simsAggIndels[0])
        plt.fill_between(x,np.array(simulationsAggregatedIndelsLows),np.array(simulationsAggregatedIndelsHighs),facecolor='lightgreen',zorder=5)


    if (realAggregatedDinucs is not None):
        aggDinucs = plt.plot(x, realAggregatedDinucs, 'crimson', label='Aggregated dinucs',linewidth=3,zorder=10)
        listofLegends.append(aggDinucs[0])
    if simulationsAggregatedDinucsMedians is not None:
        simsAggDinucs = plt.plot(x, simulationsAggregatedDinucsMedians, color='gray', linestyle=':',label='Average Simulations Aggregated Dinucs', linewidth=3,zorder=10)
        listofLegends.append(simsAggDinucs[0])
        plt.fill_between(x,np.array(simulationsAggregatedDinucsLows),np.array(simulationsAggregatedDinucsHighs),facecolor='lightpink',zorder=5)


    # old code
    # plt.legend(loc='lower left', prop={'size': 24},  shadow=False, edgecolor='white', facecolor ='white')
    plt.legend(loc= 'lower left',handles = listofLegends, prop={'size': 24}, shadow=False, edgecolor='white', facecolor ='white')

    # put the number of snps and indels
    if sample is not None:
        text = '%d subs, %d indels %d dinucs' %(numberofSPMs,numberofIndels,numberofDinucs)
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
def plotSignatureBasedFigures(mutationType,signature2NumberofMutationsDict,sample2Signature2NumberofMutationsDict,outputDir,jobname,isFigureAugmentation,numberofSimulations):
    if (mutationType==SUBS):
        xlabel = 'Interval around single point mutation (bp)'
        ylabel = 'Average nucleosome signal'
        label = 'Aggregated Substitutions'
        text = 'subs'
        color = 'royalblue'
        fillcolor = 'lightblue'
    elif (mutationType==INDELS):
        xlabel = 'Interval around indel (bp)'
        ylabel = 'Average nucleosome signal'
        label = 'Aggregated Indels'
        text = 'indels'
        color = 'darkgreen'
        fillcolor = 'lightgreen'
    elif (mutationType==DINUCS):
        xlabel = 'Interval around dinuc (bp)'
        ylabel = 'Average nucleosome signal'
        label = 'Aggregated Dinucs'
        text = 'dinucs'
        color = 'crimson'
        fillcolor = 'lightpink'

    for signature in signature2NumberofMutationsDict:
        signatureBasedNumberofMutations = signature2NumberofMutationsDict[signature]
        plotSignatureBasedAverageNucleosomeOccupancyFigureWithSimulations(None, signature,
                                                                          signatureBasedNumberofMutations,
                                                                          xlabel,ylabel,label,text,
                                                                          outputDir, jobname, isFigureAugmentation,
                                                                          numberofSimulations, color,fillcolor)

    # SampleBased Subs SignatureBased Nucleosome Occupancy Figures
    for sample in sample2Signature2NumberofMutationsDict:
        for signature in sample2Signature2NumberofMutationsDict[sample]:
            sampleBasedSignatureBasedNumberofMutations = sample2Signature2NumberofMutationsDict[sample][signature]
            plotSignatureBasedAverageNucleosomeOccupancyFigureWithSimulations(sample, signature,
                                                                              sampleBasedSignatureBasedNumberofMutations,
                                                                              xlabel,ylabel,label,text,
                                                                              outputDir, jobname, isFigureAugmentation,
                                                                              numberofSimulations,color,fillcolor)


#########################################################

#########################################################
def nucleosomeOccupancyAverageSignalFigures(outputDir,jobname,figureAugmentation,numberofSimulations,mutationTypes):

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
    sample2NumberofDinucsDict = getDictionary(outputDir,jobname,Sample2NumberofDinucsDictFilename)

    subsSignature2NumberofMutationsDict = getSubsSignature2NumberofMutationsDict(outputDir,jobname)
    indelsSignature2NumberofMutationsDict = getIndelsSignature2NumberofMutationsDict(outputDir,jobname)
    dinucsSignature2NumberofMutationsDict = getDictionary(outputDir, jobname,DinucsSignature2NumberofMutationsDictFilename)

    sample2SubsSignature2NumberofMutationsDict = getSample2SubsSignature2NumberofMutationsDict(outputDir,jobname)
    sample2IndelsSignature2NumberofMutationsDict = getSample2IndelsSignature2NumberofMutationsDict(outputDir, jobname)
    sample2DinucsSignature2NumberofMutationsDict = getDictionary(outputDir, jobname,Sample2DinucsSignature2NumberofMutationsDictFilename)
    ############## Read necessary dictionaries ends ##########################################

    ##############################################################
    #Plot "all samples pooled" and "sample based" signature based in one figure
    plotAllSamplesPooledAndSampleBasedSignaturesFiguresInOneFigure(subsSignature2NumberofMutationsDict,sample2SubsSignature2NumberofMutationsDict,outputDir,jobname,'royalblue','Interval around single point mutation (bp)','Average nucleosome signal')
    plotAllSamplesPooledAndSampleBasedSignaturesFiguresInOneFigure(indelsSignature2NumberofMutationsDict,sample2IndelsSignature2NumberofMutationsDict,outputDir,jobname,'darkgreen','Interval around indel (bp)','Average nucleosome signal')
    plotAllSamplesPooledAndSampleBasedSignaturesFiguresInOneFigure(dinucsSignature2NumberofMutationsDict,sample2DinucsSignature2NumberofMutationsDict,outputDir,jobname,'crimson','Interval around indel (bp)','Average nucleosome signal')
    ##############################################################

    ##############################################################
    plotAllMutationsPooledWithSimulations('Interval around variant (bp)','Average nucleosome signal',None,outputDir,jobname,0,0,0,numberofSimulations,mutationTypes)
    ##############################################################

    ##############################################################
    #Arrays are filled w.r.t. sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict
    samplesfromSubs  = sample2NumberofSubsDict.keys()
    samplesfromIndels = sample2NumberofIndelsDict.keys()
    samplesfromDinucs = sample2NumberofDinucsDict.keys()

    for sample in (samplesfromSubs | samplesfromIndels | samplesfromDinucs):
        numberofSubs = 0
        numberofIndels = 0
        numberofDinucs = 0
        if sample in sample2NumberofSubsDict:
            numberofSubs = sample2NumberofSubsDict[sample]
        if sample in sample2NumberofIndelsDict:
            numberofIndels = sample2NumberofIndelsDict[sample]
        if sample in sample2NumberofDinucsDict:
            numberofDinucs = sample2NumberofDinucsDict[sample]

            plotAllMutationsPooledWithSimulations('Interval around variant (bp)','Average nucleosome signal',sample,outputDir,jobname,numberofSubs, numberofIndels,numberofDinucs,numberofSimulations,mutationTypes)
    ##############################################################



    #############################################################################################################################################
    #Plot Signature Based
    #ncomms11383 Fig3b signature based average nucleosome occupancy figures
    if checkValidness(SIGNATUREBASED,outputDir,jobname):
        if (SUBS in mutationTypes):
            #Subs Signatures
            plotSignatureBasedFigures(SUBS,subsSignature2NumberofMutationsDict,sample2SubsSignature2NumberofMutationsDict,outputDir,jobname,isFigureAugmentation,numberofSimulations)
        if (INDELS in mutationTypes):
            #Indels Signatures
            plotSignatureBasedFigures(INDELS,indelsSignature2NumberofMutationsDict,sample2IndelsSignature2NumberofMutationsDict,outputDir,jobname,isFigureAugmentation,numberofSimulations)
        if (DINUCS in mutationTypes):
            # Dinucs Signatures
            plotSignatureBasedFigures(DINUCS,dinucsSignature2NumberofMutationsDict,sample2DinucsSignature2NumberofMutationsDict,outputDir,jobname,isFigureAugmentation,numberofSimulations)
    #############################################################################################################################################
