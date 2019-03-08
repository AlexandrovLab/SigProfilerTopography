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

import matplotlib.patches as mpatches


#############################################################
current_abs_path = os.path.abspath(os.path.dirname(__file__))
print('NucleosomeOccupancyAverageSignalFigures.py current_abs_path:%s' %(current_abs_path))
#############################################################


commonsPath = os.path.join(current_abs_path, '..','commons')
sys.path.append(commonsPath)

from SigProfilerTopography.source.commons.TopographyCommons import *

plusOrMinus = 1000
windowSize = plusOrMinus*2+1



# #############################################################################
# Not used
# class AnyObject1(object):
#     pass
#
# class AnyObject2(object):
#     pass
#
# class AnyObject3(object):
#     pass
#
# class AnyObject4(object):
#     pass
#
#
# class AnyObjectHandler1(object):
#     def legend_artist(self, legend, orig_handle, fontsize, handlebox):
#         x0, y0 = handlebox.xdescent, handlebox.ydescent
#         width, height = handlebox.width, handlebox.height
#         patch = mpatches.Rectangle([x0, y0], width, height, facecolor='royalblue',
#                                    edgecolor='royalblue', lw=1,
#                                    transform=handlebox.get_transform())
#         handlebox.add_artist(patch)
#         return patch
#
# class AnyObjectHandler2(object):
#     def legend_artist(self, legend, orig_handle, fontsize, handlebox):
#         x0, y0 = handlebox.xdescent, handlebox.ydescent
#         width, height = handlebox.width, handlebox.height
#         patch = mpatches.Rectangle([x0, y0], width, height, facecolor='lightblue',
#                                    edgecolor='gray', lw=1, linestyle='--',
#                                    transform=handlebox.get_transform())
#         handlebox.add_artist(patch)
#         return patch
#
#
# class AnyObjectHandler3(object):
#     def legend_artist(self, legend, orig_handle, fontsize, handlebox):
#         x0, y0 = handlebox.xdescent, handlebox.ydescent
#         width, height = handlebox.width, handlebox.height
#         patch = mpatches.Rectangle([x0, y0], width, height, facecolor='darkgreen',
#                                    edgecolor='darkgreen', lw=1,
#                                    transform=handlebox.get_transform())
#         handlebox.add_artist(patch)
#         return patch
#
#
# class AnyObjectHandler4(object):
#     def legend_artist(self, legend, orig_handle, fontsize, handlebox):
#         x0, y0 = handlebox.xdescent, handlebox.ydescent
#         width, height = handlebox.width, handlebox.height
#         patch = mpatches.Rectangle([x0, y0], width, height, facecolor='lightgreen',
#                                    edgecolor='gray', lw=1, linestyle='--',
#                                    transform=handlebox.get_transform())
#         handlebox.add_artist(patch)
#         return patch
# #############################################################################




#############################################################################
##################### Read Average as Pandas Series #########################
#############################################################################
#Jobname has to be only jobname given in the argument
def readAverage(sample,signatureName,analyseType,jobname):
    # Read the file w.r.t. the current folder

    if (analyseType == SIGNATUREBASED):
        # for signature based name may contain empty spaces
        #new way
        filename = '%s_AverageNucleosomeSignalArray.txt' %(signatureName)
        averageFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, DATA,NUCLEOSOMEOCCUPANCY, analyseType, filename)

    elif (analyseType== AGGREGATEDSUBSTITUTIONS or analyseType == AGGREGATEDINDELS):
        # new way
        filename = '%s_AverageNucleosomeSignalArray.txt' %(jobname)
        averageFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, DATA,NUCLEOSOMEOCCUPANCY, analyseType, filename)

    #####################################################
    if (analyseType == SAMPLEBASED_SIGNATUREBASED):
        # for signature based name may contain empty spaces
        # new way
        filename = '%s_%s_AverageNucleosomeSignalArray.txt' % (signatureName, sample)

        averageFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, DATA,NUCLEOSOMEOCCUPANCY, SIGNATUREBASED, filename)

    elif (analyseType == SAMPLEBASED_AGGREGATEDSUBSTITUTIONS):
        # new way
        filename = '%s_%s_AverageNucleosomeSignalArray.txt' % (sample, jobname)
        averageFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, DATA,NUCLEOSOMEOCCUPANCY, AGGREGATEDSUBSTITUTIONS, filename)

    elif (analyseType == SAMPLEBASED_AGGREGATEDINDELS):
        # new way
        filename = '%s_%s_AverageNucleosomeSignalArray.txt' % (sample, jobname)
        averageFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, DATA,NUCLEOSOMEOCCUPANCY, AGGREGATEDINDELS, filename)

    #####################################################

    return readAsNumpyArray(averageFilePath)
#############################################################################
##################### Read Average as Pandas Series #########################
#############################################################################


#############################################################################
def readAsPandasSeries(averageFilePath):
    if os.path.exists(averageFilePath):
        average = pd.read_table(averageFilePath, sep="\t", header=None, comment='#', dtype=float)
        print(average.shape)
        print('number of rows:%s' % average.shape[0])
        print('number of columns:%s' % average.shape[1])
        print(average.head())

        #former
        #average = average.loc[0,0:(windowSize-1)]

        #In the analyses, I provide data for 0:4000 for windowSize of 2000
        #However in the plots we are interested in [1000:3000]
        average = average.loc[0, plusOrMinus:(plusOrMinus*3)]
        # print('Read average.shape')
        # print(average.shape)
        return average
    else:
        return None
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
def readAverageForSimulations(sample,signature,analyseType,jobname,numberofSimulations):

    listofAverages = []

    # Read the file w.r.t. the current folder

    if (analyseType == SIGNATUREBASED):
        # for signature based name may contain empty spaces
        #new way
        for i in range(1,numberofSimulations+1):
            simulationJobName = '%s_Sim%d' %(jobname,i)
            filename = '%s_AverageNucleosomeSignalArray.txt' %(signature)
            averageFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, simulationJobName, DATA,NUCLEOSOMEOCCUPANCY, analyseType, filename)
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
            averageFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, simulationJobName, DATA,NUCLEOSOMEOCCUPANCY, analyseType, filename)
            averageNucleosomeOccupancy = readAsNumpyArray(averageFilePath)
            if (averageNucleosomeOccupancy is not None):
                listofAverages.append(averageNucleosomeOccupancy)

    #########################################################################
    if (analyseType == SAMPLEBASED_SIGNATUREBASED):
        for i in range(1,numberofSimulations+1):
            simulationJobName = '%s_Sim%d' %(jobname,i)
            filename = '%s_%s_AverageNucleosomeSignalArray.txt' %(signature,sample)
            averageFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, simulationJobName, DATA,NUCLEOSOMEOCCUPANCY, SIGNATUREBASED, filename)
            averageNucleosomeOccupancy = readAsNumpyArray(averageFilePath)
            if (averageNucleosomeOccupancy is not None):
                listofAverages.append(averageNucleosomeOccupancy)

    elif (analyseType == SAMPLEBASED_AGGREGATEDSUBSTITUTIONS):
        for i in range(1, numberofSimulations + 1):
            simulationJobName = '%s_Sim%d' % (jobname, i)
            filename = '%s_%s_AverageNucleosomeSignalArray.txt' % (sample,simulationJobName)
            averageFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT,simulationJobName, DATA, NUCLEOSOMEOCCUPANCY,AGGREGATEDSUBSTITUTIONS, filename)
            averageNucleosomeOccupancy = readAsNumpyArray(averageFilePath)
            if (averageNucleosomeOccupancy is not None):
                listofAverages.append(averageNucleosomeOccupancy)

    elif (analyseType == SAMPLEBASED_AGGREGATEDINDELS):
        for i in range(1, numberofSimulations + 1):
            simulationJobName = '%s_Sim%d' % (jobname, i)
            filename = '%s_%s_AverageNucleosomeSignalArray.txt' % (sample,simulationJobName)
            averageFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT,simulationJobName, DATA, NUCLEOSOMEOCCUPANCY,AGGREGATEDINDELS, filename)
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
def plotSignatureBasedAverageNucleosomeOccupancyFigureWithSimulations(sample,signature,numberofMutations,xlabel,ylabel,jobname, isFigureAugmentation,numberofSimulations):

    simulationsSignatureBasedMedians = None
    listofSimulationsSignatureBased = None

    if ((sample is not None) and (signature is not None)):
        realAverage = readAverage(sample, signature, SAMPLEBASED_SIGNATUREBASED, jobname)
        figurename = '%s_%s_%d' % (signature, sample, numberofMutations)
        title = '%s_%s' % (signature, sample)
        if (numberofSimulations>0):
            listofSimulationsSignatureBased = readAverageForSimulations(sample, signature, SAMPLEBASED_SIGNATUREBASED, jobname,numberofSimulations)
    else:
        realAverage = readAverage(None, signature, SIGNATUREBASED, jobname)
        figurename = '%s_%d' % (signature, numberofMutations)
        title = '%s' % (signature)
        if (numberofSimulations>0):
            listofSimulationsSignatureBased = readAverageForSimulations(sample, signature, SIGNATUREBASED, jobname,numberofSimulations)

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

        # color = signatureBasedColors[name]
        # color = 'navy'
        color = 'royalblue'

        listofLegends = []

        if (realAverage is not None):
            original = plt.plot(x, realAverage, color='royalblue', label='Aggregated substitutions',linewidth=3)
            listofLegends.append(original[0])

        if (simulationsSignatureBasedMedians is not None):
            simulations = plt.plot(x, simulationsSignatureBasedMedians, color='gray', linestyle='--',  label='Average Simulations Aggregated substitutions', linewidth=3)
            listofLegends.append(simulations[0])
            plt.fill_between(x, np.array(simulationsSignatureBasedLows), np.array(simulationsSignatureBasedHighs),facecolor='lightblue')

        plt.legend(loc= 'lower left', handles=listofLegends, prop={'size': 24}, shadow=False, edgecolor='white', facecolor='white')

        #put the number of snps
        text = '%d subs' %(numberofMutations)
        plt.text(0.85,0.95, text, ha='center', va='center', transform=ax.transAxes, fontsize=24)

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

        filename = figurename.replace(' ', '') + '_AverageNucleosomeOccupancySignal.png'

        #######################################################################
        # new code
        if (sample is None):
            figureFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, FIGURE, ALL, NUCLEOSOMEOCCUPANCY, filename)
        else:
            os.makedirs(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, FIGURE, SAMPLES, sample, NUCLEOSOMEOCCUPANCY), exist_ok=True)
            figureFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, FIGURE, SAMPLES, sample, NUCLEOSOMEOCCUPANCY, filename)
        #######################################################################

        fig.savefig(figureFile)
        plt.close(fig)
#############################################################################
########################## Plot Figure starts  ##############################
#############################################################################


#############################################################################
############## Plot AggregatedIndels For Simulations starts #################
#############################################################################
def plotAggregatedIndelsWithSimulations(xlabel,ylabel,sample,signature,analyseType,jobname,isFigureAugmentation,numberofIndels,numberofSimulations):

    realAggregatedIndels = readAverage(sample, signature, analyseType, jobname)
    simulationsAggregatedIndelsMedians = None
    listofSimulationsAggregatedIndels = None

    if (sample is None):
        figurename = '%s_AggregatedIndels.png' % (jobname)
        if (numberofSimulations>0):
            listofSimulationsAggregatedIndels =  readAverageForSimulations(sample,None, AGGREGATEDINDELS, jobname,numberofSimulations)
    else:
        figurename = '%s_%d_AggregatedIndels_%s.png' % (sample,numberofIndels,jobname)
        if (numberofSimulations>0):
            listofSimulationsAggregatedIndels = readAverageForSimulations(sample, None, SAMPLEBASED_AGGREGATEDINDELS,jobname, numberofSimulations)

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
    plt.text(0.85, 0.95, text, ha='center', va='center', transform=ax.transAxes, fontsize=24)

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
        figureFile = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,OUTPUT,jobname, FIGURE,ALL,NUCLEOSOMEOCCUPANCY, figurename)
    else:
        os.makedirs(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, FIGURE, SAMPLES, sample,NUCLEOSOMEOCCUPANCY), exist_ok=True)
        figureFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, FIGURE,SAMPLES, sample,NUCLEOSOMEOCCUPANCY, figurename)
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
def plotAggregatedSubstitutionsWithSimulations(xlabel,ylabel,sample,signature,analyseType,jobname,isFigureAugmentation,sampleBasedNumberofMutations,numberofSimulations):

    simulationsAggregatedSubstitutionsMedians = None
    listofSimulationsAggregatedSubstitutions = None

    realAggregatedSubstitutions = readAverage(sample, signature, analyseType, jobname)

    if (sample is None):
        filename = '%s_Aggregated_Substitutions.png' %(jobname)
        if (numberofSimulations>0):
            listofSimulationsAggregatedSubstitutions = readAverageForSimulations(None, None, AGGREGATEDSUBSTITUTIONS,jobname, numberofSimulations)

    else:
        filename = '%s_%d_Aggregated_Substitutions_%s.png' %(sample,sampleBasedNumberofMutations,jobname)
        if (numberofSimulations>0):
            listofSimulationsAggregatedSubstitutions = readAverageForSimulations(sample, None, SAMPLEBASED_AGGREGATEDSUBSTITUTIONS, jobname,numberofSimulations)

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
    plt.text(0.85, 0.95, text, ha='center', va='center', transform=ax.transAxes, fontsize=24)

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
        figureFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, FIGURE, ALL, NUCLEOSOMEOCCUPANCY, filename)
    else:
        os.makedirs(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, FIGURE, SAMPLES, sample,NUCLEOSOMEOCCUPANCY), exist_ok=True)
        figureFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, FIGURE, SAMPLES, sample, NUCLEOSOMEOCCUPANCY, filename)
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
#For Debugging starts FEB 26, 2019
#############################################################################
def plotSignalsandCountsForDebug(sample,jobname,numberofSimulations):
    if sample is None:
        filename = '%s_Signals_Counts.png' % (jobname)
    else:
        filename = '%s_%s_Signals_Counts.png' % (sample,jobname)

    listofLegends = []

    listofSimulationsAggregatedSubstitutions_Signal = None
    listofSimulationsAggregatedSubstitutions_Count = None

    original_Signal = readSignalorCount(sample,jobname,'AccumulatedSignalArray.txt',AGGREGATEDSUBSTITUTIONS)
    original_Count = readSignalorCount(sample,jobname, 'AccumulatedCountArray.txt',AGGREGATEDSUBSTITUTIONS)

    if (numberofSimulations > 0):
        listofSimulationsAggregatedSubstitutions_Signal = readSignalorCountSimulations(sample,jobname,'AccumulatedSignalArray.txt', AGGREGATEDSUBSTITUTIONS, numberofSimulations)
        listofSimulationsAggregatedSubstitutions_Count = readSignalorCountSimulations(sample,jobname,'AccumulatedCountArray.txt', AGGREGATEDSUBSTITUTIONS, numberofSimulations)

    stackedSimulations_Signal = np.vstack(listofSimulationsAggregatedSubstitutions_Signal)
    stackedSimulations_Count = np.vstack(listofSimulationsAggregatedSubstitutions_Count)

    (rowsSignal, colsSignal) = stackedSimulations_Signal.shape
    print('For Debug FEB26, 2019 rowsSignal:%d colsSignal:%d' %(rowsSignal,colsSignal))
    (rowsCount, colsCount) = stackedSimulations_Count.shape
    print('For Debug FEB26, 2019 rowsCount:%d colsCount:%d' % (rowsCount,colsCount))

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

    if (stackedSimulations_Signal  is not None):
        for row in range(rowsSignal):
            print('Simulation Signal %d' %(row))
            simulation_signal_array = stackedSimulations_Signal[row,:]
            simulation_signal = plt.plot(x, simulation_signal_array, 'blue', label='Simulations Signal',linewidth=3,zorder=5)
        listofLegends.append(simulation_signal[0])

    if (stackedSimulations_Count  is not None):
        for row in range(rowsCount):
            print('Simulation Count %d' %(row))
            simulation_count_array = stackedSimulations_Count[row,:]
            # simulation_count = plt.plot(x, simulation_count_array, 'green', label='Simulations Count', linestyle='--',linewidth=3,zorder=5)
            simulation_count = plt.plot(x, simulation_count_array, 'green', label='Simulations Count',linewidth=3, zorder=5)
        listofLegends.append(simulation_count[0])

    if (original_Signal is not None):
        aggSubs = plt.plot(x, original_Signal, 'black', label='Original Signal',linewidth=5,zorder=10)
        listofLegends.append(aggSubs[0])

    if (original_Count is not None):
        aggSubs = plt.plot(x, original_Count, 'red', label='Original Count',linewidth=5,zorder=10)
        listofLegends.append(aggSubs[0])


    plt.legend(loc= 'upper left',handles = listofLegends, prop={'size': 24}, shadow=False, edgecolor='white', facecolor ='white')

    #Put vertical line at x=0
    plt.axvline(x=0, color='gray', linestyle='--')

    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.tick_params(axis='both', which='minor', labelsize=24)

    # This code puts the tick marks
    plt.tick_params(axis='both', which='major', labelsize=24,width=3,length=10)
    plt.tick_params(axis='both', which='minor', labelsize=24,width=3,length=10)

    # This code provides the x and y tick marks and labels
    plt.xticks(np.arange(-1000, +1001, step=500), fontsize=30)

    plt.xlim((-1000, 1000))
    # This code provides some extra space
    # plt.ylim((0.65,1.15))

    if sample is not None:
        title = '%s_%s' %(sample,jobname)
    else:
        title = jobname

    plt.title(title, fontsize=40, fontweight='bold')

    plt.xlabel('Single Point Substitutions', fontsize=30)
    plt.ylabel('Signal & Count', fontsize=30)

    figureFile = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,OUTPUT,jobname,FIGURE,ALL,NUCLEOSOMEOCCUPANCY,filename)

    fig.savefig(figureFile)
    plt.close(fig)
#############################################################################
#For Debugging ends FEB 26, 2019
#############################################################################


#############################################################################
#For Debugging starts FEB 26, 2019
#############################################################################
def readSignalorCount(sample,jobname,filename,analyseType):
    if sample is not None:
        filename = '%s_%s_%s' % (sample,jobname,filename)
    else:
        filename = '%s_%s' % (jobname,filename)

    filepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, DATA, NUCLEOSOMEOCCUPANCY, analyseType, filename)
    return readAsNumpyArray(filepath)
#############################################################################
#For Debugging ends FEB 26, 2019
#############################################################################


#############################################################################
#For Debugging starts FEB 26, 2019
#############################################################################
def readSignalorCountSimulations(sample,jobname,filename,analyseType,numberofSimulations):
    listofArrays = []

    for i in range(1, numberofSimulations + 1):
        simulationJobName = '%s_Sim%d' % (jobname, i)
        if sample is None:
            newfilename = '%s_%s' % (simulationJobName,filename)
        else:
            newfilename = '%s_%s_%s' % (sample,simulationJobName,filename)

        filepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, simulationJobName,DATA, NUCLEOSOMEOCCUPANCY, analyseType, newfilename)
        print('For Debug %s' %(filepath))
        signalorCountArray = readAsNumpyArray(filepath)
        if (signalorCountArray is not None):
            listofArrays.append(signalorCountArray)

    print('len(listofArrays): %d' %len(listofArrays))
    return listofArrays
#############################################################################
#For Debugging ends FEB 26, 2019
#############################################################################


#############################################################################
############################ Plot Figure ####################################
#############################################################################
def plotAggregatedSubstitutionsandAggregatedIndelsWithSimulations(xlabel,ylabel,sample,jobname,isFigureAugmentation,numberofSPMs,numberofIndels,numberofSimulations):

    simulationsAggregatedSubstitutionsMedians = None
    simulationsAggregatedIndelsMedians = None

    listofSimulationsAggregatedIndels = None
    listofSimulationsAggregatedSubstitutions = None


    #######################################################################################################################
    if (sample is None):
        filename = '%s_Aggregated_Substitutions_Indels.png' % (jobname)
        realAggregatedIndels = readAverage(None,None, AGGREGATEDINDELS, jobname)
        realAggregatedSubstitutions = readAverage(None,None,AGGREGATEDSUBSTITUTIONS,jobname)

        if (numberofSimulations>0):
            listofSimulationsAggregatedIndels = readAverageForSimulations(None, None, AGGREGATEDINDELS, jobname,numberofSimulations)
            listofSimulationsAggregatedSubstitutions =  readAverageForSimulations(None,None,AGGREGATEDSUBSTITUTIONS,jobname,numberofSimulations)
    #######################################################################################################################


    #######################################################################################################################
    else:
        filename = '%s_%s_Aggregated_Substitutions_%d_Indels_%d.png' % (sample, jobname, numberofSPMs, numberofIndels)
        realAggregatedIndels = readAverage(sample,None, SAMPLEBASED_AGGREGATEDINDELS, jobname)
        realAggregatedSubstitutions = readAverage(sample,None,SAMPLEBASED_AGGREGATEDSUBSTITUTIONS,jobname)

        if (numberofSimulations>0):
            listofSimulationsAggregatedIndels = readAverageForSimulations(sample, None, SAMPLEBASED_AGGREGATEDINDELS,jobname, numberofSimulations)
            listofSimulationsAggregatedSubstitutions = readAverageForSimulations(sample, None, SAMPLEBASED_AGGREGATEDSUBSTITUTIONS, jobname,numberofSimulations)
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

    # plt.legend([AnyObject1(), AnyObject2(), AnyObject3(), AnyObject4()],
    #            ['Aggregated substitutions',
    #             'Average simulations aggregated substitutions',
    #             'Aggregated indels',
    #             'Average simulations aggregated indels'],
    #            handler_map={AnyObject1: AnyObjectHandler1(), AnyObject2: AnyObjectHandler2(),
    #                         AnyObject3: AnyObjectHandler3(), AnyObject4: AnyObjectHandler4()},
    #            loc='lower left', shadow=False, edgecolor='white', facecolor='white',prop={'size': 24})

    # put the number of snps and indels
    if sample is not None:
        text = '%d subs, %d indels' %(numberofSPMs,numberofIndels)
        plt.text(0.85, 0.95, text, ha='center', va='center', transform=ax.transAxes, fontsize=24)

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

    # if (isFigureAugmentation):
    #     plt.title(jobname, fontsize=36)

    if (sample is not None):
        plt.title(sample, fontsize=40, fontweight='bold')
    else:
        plt.title(jobname, fontsize=40, fontweight='bold')


    plt.xlabel(xlabel, fontsize=30)
    plt.ylabel(ylabel, fontsize=30)

    ######################################################################################
    #new code
    if (sample is None):
        figureFile = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,OUTPUT,jobname,FIGURE,ALL,NUCLEOSOMEOCCUPANCY,filename)
    else:
        os.makedirs(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, FIGURE, SAMPLES, sample,NUCLEOSOMEOCCUPANCY), exist_ok=True)
        figureFile = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,OUTPUT,jobname,FIGURE,SAMPLES,sample,NUCLEOSOMEOCCUPANCY,filename)
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
def checkValidness(analsesType,jobname):
    #Check whether directory exists and there are files under it.

    data_file_path = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,OUTPUT,jobname,DATA,NUCLEOSOMEOCCUPANCY,analsesType)

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
def nucleosomeOccupancyAverageSignalFigures(jobname,figureAugmentation,numberofSimulations):

    current_abs_path = os.path.dirname(os.path.realpath(__file__))

    #######################################################################################################################
    os.makedirs(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, FIGURE, ALL, NUCLEOSOMEOCCUPANCY), exist_ok=True)
    os.makedirs(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, FIGURE, SAMPLES), exist_ok=True)
    #######################################################################################################################

    isFigureAugmentation = False
    if (figureAugmentation=='augmentation'):
        isFigureAugmentation = True

    #Signature Based starts
    analyseType = SIGNATUREBASED
    #get the signatures, they may have spaces before, after or in between

    ##########################################################################################
    ############## Read necessary dictionaries starts ########################################
    ##########################################################################################

    ##########################################################################################
    sample2NumberofIndelsDict = {}
    sample2NumberofIndelsDictFilePath = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,OUTPUT,jobname,DATA,Samples2NumberofIndelsDictFilename)

    if (os.path.exists(sample2NumberofIndelsDictFilePath)):
        sample2NumberofIndelsDict = readDictionary(sample2NumberofIndelsDictFilePath)
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
    ############## Read necessary dictionaries ends ##########################################
    ##########################################################################################


    # Plot All Together : Aggregated Substitutions and Aggregated Indels
    # Or plot one of them:  Aggregated Substitutions or Aggregated Indels
    if checkValidness(AGGREGATEDSUBSTITUTIONS,jobname) and checkValidness(AGGREGATEDINDELS,jobname):

        ##############################################################
        plotAggregatedSubstitutionsandAggregatedIndelsWithSimulations('Interval around variant (bp)','Average nucleosome signal',
                                                                          None,jobname,isFigureAugmentation,
                                                                          0,0,numberofSimulations)
        ##############################################################


        ##############################################################
        #Arrays are filled w.r.t. sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict
        for sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
        # for sample in samplesWithAtLeast10KMutations2NumberofMutationsDict:
            numberofEligibleSPMs = 0
            numberofIndels = 0
            if sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
                numberofEligibleSPMs = sum(sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict[sample].values())
            if sample in sample2NumberofIndelsDict:
                numberofIndels = sample2NumberofIndelsDict[sample]

            plotAggregatedSubstitutionsandAggregatedIndelsWithSimulations('Interval around variant (bp)','Average nucleosome signal',
                                                                        sample,jobname, isFigureAugmentation,
                                                                        numberofEligibleSPMs, numberofIndels,numberofSimulations)

            #For debug FEB 28, 2019
            plotSignalsandCountsForDebug(sample,jobname, numberofSimulations)
        ##############################################################

    # ######################################################################################
    # elif (checkValidness(AGGREGATEDSUBSTITUTIONS,jobname)):
    #     sampleBasedNumberofMutations = 0
    #     plotAggregatedSubstitutionsWithSimulations('Interval around single point mutation (bp)', 'Average nucleosome signal',
    #                                                 None, None, AGGREGATEDSUBSTITUTIONS,
    #                                                 jobname,isFigureAugmentation,sampleBasedNumberofMutations,numberofSimulations)
    #
    #     for sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
    #     # for sample in samplesWithAtLeast10KMutations2NumberofMutationsDict:
    #         sampleBasedNumberofMutations = samplesWithAtLeast10KMutations2NumberofMutationsDict[sample]
    #         plotAggregatedSubstitutionsWithSimulations('Interval around single point mutation (bp)', 'Average nucleosome signal',
    #                                     sample, None, SAMPLEBASED_AGGREGATEDSUBSTITUTIONS,
    #                                     jobname, isFigureAugmentation, sampleBasedNumberofMutations,numberofSimulations)
    # ######################################################################################

    # ######################################################################################
    # elif (checkValidness(AGGREGATEDINDELS,jobname)):
    #     numberofIndels = 0
    #     plotAggregatedIndelsWithSimulations('Interval around variant (bp)', 'Average nucleosome signal',
    #                                         None, None,
    #                                         AGGREGATEDINDELS, jobname, isFigureAugmentation,numberofIndels,numberofSimulations)
    #
    #     for sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
    #     # for sample in samplesWithAtLeast10KMutations2NumberofMutationsDict:
    #         if sample in sample2NumberofIndelsDict:
    #             numberofIndels = sample2NumberofIndelsDict[sample]
    #             plotAggregatedIndelsWithSimulations('Interval around variant (bp)', 'Average nucleosome signal',
    #                                             sample, None,
    #                                             SAMPLEBASED_AGGREGATEDINDELS, jobname, isFigureAugmentation,numberofIndels,numberofSimulations)
    # ######################################################################################



    # #Plot Signature Based
    # #Plot ncomms11383 Fig3b signature based average nucleosome occupancy figures
    # if checkValidness(SIGNATUREBASED,jobname):
    #     #SignatureBased
    #     for signature in signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
    #         signatureBasedNumberofMutations = signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict[signature]
    #         plotSignatureBasedAverageNucleosomeOccupancyFigureWithSimulations(None,signature,signatureBasedNumberofMutations,
    #                                                                               'Interval around single point mutation (bp)','Average nucleosome signal',
    #                                                                               jobname,isFigureAugmentation,numberofSimulations)
    #
    #     # SampleBased SignatureBased Nucleosome Occupancy Figures
    #     for sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
    #         for signature in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict[sample]:
    #             sampleBasedSignatureBasedNumberofMutations = sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict[sample][signature]
    #             plotSignatureBasedAverageNucleosomeOccupancyFigureWithSimulations(sample, signature,sampleBasedSignatureBasedNumberofMutations,
    #                                                                                   'Interval around single point mutation (bp)','Average nucleosome signal',
    #                                                                                   jobname, isFigureAugmentation,numberofSimulations)

    #########################################################

    #For debugging FEB 26, 2019
    plotSignalsandCountsForDebug(None,jobname,numberofSimulations)
