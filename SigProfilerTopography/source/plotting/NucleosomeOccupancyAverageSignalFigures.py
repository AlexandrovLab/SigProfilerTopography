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

# plusOrMinus = 1000
# windowSize = plusOrMinus*2+1

plt.rcParams.update({'figure.max_open_warning': 0})

#############################################################################
##################### Read Average as Pandas Series #########################
#############################################################################
#Jobname has to be only jobname given in the argument
def readData(sample,signatureName,analyseType,outputDir,jobname,occupancy_type,libraryFilenameMemo):

    partial_file_name='AverageSignalArray'

    #####################################################
    if (analyseType == SIGNATUREBASED):
        if libraryFilenameMemo is None:
            filename = '%s_%s.txt' %(signatureName,partial_file_name)
        else:
            filename = '%s_%s_%s.txt' %(signatureName,libraryFilenameMemo,partial_file_name)
        averageFilePath = os.path.join(outputDir, jobname, DATA, occupancy_type, analyseType, filename)
    elif ((analyseType== AGGREGATEDSUBSTITUTIONS) or (analyseType == AGGREGATEDINDELS) or (analyseType == AGGREGATEDDINUCS)):
        # new way
        if libraryFilenameMemo is None:
            filename = '%s_%s.txt' %(jobname,partial_file_name)
        else:
            filename = '%s_%s_%s.txt' %(jobname,libraryFilenameMemo,partial_file_name)
        averageFilePath = os.path.join(outputDir, jobname, DATA, occupancy_type, analyseType, filename)
    #####################################################

    #####################################################
    if (analyseType == SAMPLEBASED_SIGNATUREBASED):
        if libraryFilenameMemo is None:
            filename = '%s_%s_%s.txt' % (signatureName, sample,partial_file_name)
        else:
            filename = '%s_%s_%s_%s.txt' % (signatureName, sample, libraryFilenameMemo,partial_file_name)
        averageFilePath = os.path.join(outputDir, jobname, DATA, SAMPLES, sample, occupancy_type, SIGNATUREBASED, filename)

    elif (analyseType == SAMPLEBASED_AGGREGATEDSUBSTITUTIONS):
        if libraryFilenameMemo is None:
            filename = '%s_%s_%s.txt' % (sample, jobname,partial_file_name)
        else:
            filename = '%s_%s_%s_%s.txt' % (sample, jobname,libraryFilenameMemo,partial_file_name)
        averageFilePath = os.path.join(outputDir, jobname, DATA, SAMPLES, sample, occupancy_type, AGGREGATEDSUBSTITUTIONS, filename)

    elif (analyseType == SAMPLEBASED_AGGREGATEDINDELS):
        if libraryFilenameMemo is None:
            filename = '%s_%s_%s.txt' % (sample, jobname,partial_file_name)
        else:
            filename = '%s_%s_%s_%s.txt' % (sample, jobname,libraryFilenameMemo,partial_file_name)
        averageFilePath = os.path.join(outputDir, jobname, DATA, SAMPLES, sample, occupancy_type, AGGREGATEDINDELS, filename)

    elif (analyseType == SAMPLEBASED_AGGREGATEDDINUCS):
        if libraryFilenameMemo is None:
            filename = '%s_%s_%s.txt' % (sample, jobname, partial_file_name)
        else:
            filename = '%s_%s_%s_%s.txt' % (sample, jobname,libraryFilenameMemo,partial_file_name)
        averageFilePath = os.path.join(outputDir, jobname, DATA, SAMPLES, sample, occupancy_type, AGGREGATEDDINUCS, filename)
    #####################################################

    return readAsNumpyArray(averageFilePath)
#############################################################################
##################### Read Average as Pandas Series #########################
#############################################################################




#############################################################################
def readAsNumpyArray(averageFilePath):
    if os.path.exists(averageFilePath):
        average = np.loadtxt(averageFilePath, dtype=float, delimiter='\t')
        # average = average[(plusOrMinus/2):((plusOrMinus/2)*3)+1]
        return average
    else:
        return None
#############################################################################


#############################################################################
##################### Read Average for Simulations start ####################
#############################################################################
def readDataForSimulations(sample,signature,analyseType,outputDir,jobname,numberofSimulations,occupancy_type,libraryFilenameMemo):
    partial_file_name = 'AverageSignalArray'

    listofAverages = []
    # Read the file w.r.t. the current folder
    if (analyseType == SIGNATUREBASED):
        # for signature based name may contain empty spaces
        #new way
        for simNum in range(1,numberofSimulations+1):
            if (libraryFilenameMemo is None):
                filename = '%s_sim%d_%s.txt' %(signature,simNum,partial_file_name)
            else:
                filename = '%s_sim%d_%s_%s.txt' %(signature,simNum,libraryFilenameMemo,partial_file_name)
            averageFilePath = os.path.join(outputDir,jobname,DATA, occupancy_type,analyseType,filename)
            averageNucleosomeOccupancy = readAsNumpyArray(averageFilePath)
            if (averageNucleosomeOccupancy is not None):
                listofAverages.append(averageNucleosomeOccupancy)

    elif (analyseType== AGGREGATEDSUBSTITUTIONS or analyseType == AGGREGATEDINDELS or analyseType == AGGREGATEDDINUCS):
        # i=0 for real/original data
        # i=1 for Sim1
        # ...
        # i=numberofSimulations for the last Simulations numberofSimulations
        for simNum in range(1,numberofSimulations+1):
            if (libraryFilenameMemo is None):
                filename = '%s_sim%d_%s.txt' %(jobname,simNum,partial_file_name)
            else:
                filename = '%s_sim%d_%s_%s.txt' %(jobname,simNum,libraryFilenameMemo,partial_file_name)
            averageFilePath = os.path.join(outputDir, jobname, DATA, occupancy_type, analyseType, filename)
            averageNucleosomeOccupancy = readAsNumpyArray(averageFilePath)
            if (averageNucleosomeOccupancy is not None):
                listofAverages.append(averageNucleosomeOccupancy)

    #########################################################################
    if (analyseType == SAMPLEBASED_SIGNATUREBASED):
        for simNum in range(1,numberofSimulations+1):
            if (libraryFilenameMemo is None):
                filename = '%s_%s_sim%d_%s.txt' %(signature,sample,simNum,partial_file_name)
            else:
                filename = '%s_%s_sim%d_%s_%s.txt' %(signature,sample,simNum,libraryFilenameMemo,partial_file_name)
            averageFilePath = os.path.join(outputDir, jobname, DATA, SAMPLES, sample, occupancy_type, SIGNATUREBASED, filename)
            averageNucleosomeOccupancy = readAsNumpyArray(averageFilePath)
            if (averageNucleosomeOccupancy is not None):
                listofAverages.append(averageNucleosomeOccupancy)

    elif (analyseType == SAMPLEBASED_AGGREGATEDSUBSTITUTIONS):
        for simNum in range(1, numberofSimulations + 1):
            if (libraryFilenameMemo is None):
                filename = '%s_%s_sim%d_%s.txt' % (sample,jobname,simNum,partial_file_name)
            else:
                filename = '%s_%s_sim%d_%s_%s.txt' % (sample,jobname,simNum,libraryFilenameMemo,partial_file_name)
            averageFilePath = os.path.join(outputDir,jobname, DATA, SAMPLES, sample, occupancy_type,AGGREGATEDSUBSTITUTIONS, filename)
            averageNucleosomeOccupancy = readAsNumpyArray(averageFilePath)
            if (averageNucleosomeOccupancy is not None):
                listofAverages.append(averageNucleosomeOccupancy)

    elif (analyseType == SAMPLEBASED_AGGREGATEDINDELS):
        for simNum in range(1, numberofSimulations + 1):
            if (libraryFilenameMemo is None):
                filename = '%s_%s_sim%d_%s.txt' % (sample,jobname,simNum,partial_file_name)
            else:
                filename = '%s_%s_sim%d_%s_%s.txt' % (sample,jobname,simNum,libraryFilenameMemo,partial_file_name)
            averageFilePath = os.path.join(outputDir,jobname, DATA, SAMPLES, sample, occupancy_type,AGGREGATEDINDELS, filename)
            averageNucleosomeOccupancy = readAsNumpyArray(averageFilePath)
            if (averageNucleosomeOccupancy is not None):
                listofAverages.append(averageNucleosomeOccupancy)

    elif (analyseType == SAMPLEBASED_AGGREGATEDDINUCS):
        for simNum in range(1, numberofSimulations + 1):
            if (libraryFilenameMemo is None):
                filename = '%s_%s_sim%d_%s.txt' % (sample,jobname,simNum,partial_file_name)
            else:
                filename = '%s_%s_sim%d_%s_%s.txt' % (sample,jobname,simNum,libraryFilenameMemo,partial_file_name)
            averageFilePath = os.path.join(outputDir,jobname, DATA, SAMPLES, sample, occupancy_type,AGGREGATEDDINUCS, filename)
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
def plotAllSamplesPooledAndSampleBasedSignaturesFiguresInOneFigure(signature2NumberofMutationsDict,sample2Signature2NumberofMutationsDict,outputDir,jobname,color,xlabel,ylabel,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus):

    if (occupancy_type==NUCLEOSOMEOCCUPANCY):
        filenameEnd = 'NucleosomeOccupancy'
    elif (occupancy_type==EPIGENOMICSOCCUPANCY):
        filenameEnd = 'EpigenomicsOccupancy'

    for signature in signature2NumberofMutationsDict:
        min_list = []
        max_list = []

        label2NumpyArrayDict = {}
        signatureBasedNumberofMutations = signature2NumberofMutationsDict[signature]
        realAverage = readData(None, signature, SIGNATUREBASED, outputDir, jobname,occupancy_type,libraryFilenameMemo)
        label2NumpyArrayDict[signature] = realAverage

        for sample in sample2Signature2NumberofMutationsDict:
            if signature in sample2Signature2NumberofMutationsDict[sample]:
                realAverage = readData(sample,signature,SAMPLEBASED_SIGNATUREBASED,outputDir,jobname,occupancy_type,libraryFilenameMemo)
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
                max_list.append(np.amax(realAverage))
                min_list.append(np.amin(realAverage))

                if (label==signature):
                    original = plt.plot(x,realAverage,color=color,label=label,linewidth=3,zorder=10)
                else:
                    original = plt.plot(x,realAverage,color='gray',label=label,linewidth=3,zorder=5)
                listofLegends.append(original[0])


        # plt.legend(loc= 'lower left', handles=listofLegends, prop={'size': 12}, shadow=False, edgecolor='white', facecolor='white')
        plt.legend(loc='lower left', handles=listofLegends, prop={'size': 10}, shadow=False, edgecolor='white',facecolor='white', ncol=8, framealpha=0)

        # text = '%d subs' %(numberofMutations)
        # plt.text(0.99, 0.99, text, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, fontsize=24)

        #put the library filename
        plt.text(0.01, 0.99, libraryFilename, verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, fontsize=24)

        #Put vertical line at x=0
        # plt.axvline(x=0, ymin=0, ymax=1, color='gray', linestyle='--')
        plt.axvline(x=0, color='gray', linestyle='--')

        # This code provides the x and y tick marks and labels
        # plt.xticks(np.arange(-1000, +1001, step=500), fontsize=30)
        plt.xticks(np.arange(-plusOrMinus, plusOrMinus+1, step=500), fontsize=30)

        #July 27, 2018
        # plt.xlim((-1000, 1000))
        plt.xlim((-plusOrMinus, plusOrMinus))

        #Let's comment for HM figures
        # ###################################################################################
        # if (len(min_list)>0 and len(max_list)>0):
        #     ###################################################################################
        #     min_average_nucleosome_signal = np.amin(min_list)
        #     max_average_nucleosome_signal = np.amax(max_list)
        #
        #     if (min_average_nucleosome_signal >= 0.55 and max_average_nucleosome_signal <= 1.35):
        #         # Comment these lines to see the values out of range
        #         plt.yticks(np.arange(0.55, 1.4, step=0.1), fontsize=30)
        #         # This code overwrites the y tick labels
        #         ax.set_yticklabels(['0.55', '', '0.75', '', '0.95', '', '1.15', '', '1.35'])
        #         # This code provides some extra space
        #         plt.ylim((0.54, 1.4))
        #     else:
        #         yticklabels = []
        #
        #         ymin = 0.55
        #         while (ymin > min_average_nucleosome_signal):
        #             ymin -= 0.1
        #
        #         ymax = 1.35
        #         while (ymax < max_average_nucleosome_signal):
        #             ymax += 0.1
        #
        #         ymin = round(ymin, 2)
        #         ymax = round(ymax, 2)
        #
        #         ytick = ymin
        #         while (ytick <= ymax):
        #             yticklabels.append(round(ytick, 2))
        #             if (ytick < ymax):
        #                 yticklabels.append('')
        #             ytick = round(ytick + 0.1, 2)
        #
        #         if ymax not in yticklabels:
        #             yticklabels.append(ymax)
        #
        #         plt.yticks(np.arange(ymin, ymax, step=0.05), fontsize=30)
        #         ax.set_yticklabels(yticklabels)
        #         plt.ylim((ymin - 0.01, ymax + 0.01))
        #     ###################################################################################
        # else:
        #     # Comment these lines to see the values out of range
        #     plt.yticks(np.arange(0.55, 1.4, step=0.1), fontsize=30)
        #     # This code overwrites the y tick labels
        #     ax.set_yticklabels(['0.55', '', '0.75', '', '0.95', '', '1.15', '', '1.35'])
        #     # This code provides some extra space
        #     plt.ylim((0.54, 1.4))
        # ###################################################################################

        # This code puts the tick marks
        plt.tick_params(axis='both', which='major', labelsize=30,width=3,length=10)
        plt.tick_params(axis='both', which='minor', labelsize=30,width=3,length=10)

        title = signature
        plt.title(title, fontsize=40,fontweight='bold')

        plt.xlabel(xlabel,fontsize=32,fontweight='semibold')
        plt.ylabel(ylabel,fontsize=32,fontweight='semibold')

        if (libraryFilenameMemo is None):
            filename = '%s_AllInOne_%s.png' %(signature,filenameEnd)
        else:
            filename = '%s_AllInOne_%s_%s.png' %(signature,libraryFilenameMemo,filenameEnd)

        figureFile = os.path.join(outputDir, jobname, FIGURE, ALL, occupancy_type, filename)

        fig.savefig(figureFile)
        plt.cla()
        plt.close(fig)
        ####################################### plottoing ends #######################################

#############################################################################
########################## Plot Figure ends  ################################
#############################################################################


#############################################################################
########################## Plot Figure starts  ##############################
#############################################################################
#Called by plotSignatureBasedFigures
def plotSignatureBasedAverageNucleosomeOccupancyFigureWithSimulations(sample,signature,numberofMutations,xlabel,ylabel,label,text,outputDir,jobname,isFigureAugmentation,numberofSimulations,color,linestyle,fillcolor,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus):

    if (occupancy_type==NUCLEOSOMEOCCUPANCY):
        figurenameEnd='_NucleosomeOccupancy.png'
    elif (occupancy_type==EPIGENOMICSOCCUPANCY):
        figurenameEnd='_EpigenomicsOccupancy.png'

    min_list=[]
    max_list=[]
    min_average_nucleosome_signal=0
    max_average_nucleosome_signal=0

    simulationsSignatureBasedLows = None
    simulationsSignatureBasedMedians = None
    simulationsSignatureBasedHighs = None
    listofSimulationsSignatureBased = None

    if ((sample is not None) and (signature is not None)):
        realAverage = readData(sample, signature, SAMPLEBASED_SIGNATUREBASED, outputDir, jobname,occupancy_type,libraryFilenameMemo)
        if (libraryFilenameMemo is None):
            figurename = '%s_%s_%d' % (signature, sample, numberofMutations)
        else:
            figurename = '%s_%s_%d_%s' % (signature, sample, numberofMutations,libraryFilenameMemo)

        title = '%s_%s' % (signature, sample)
        if (numberofSimulations>0):
            listofSimulationsSignatureBased = readDataForSimulations(sample, signature, SAMPLEBASED_SIGNATUREBASED, outputDir, jobname,numberofSimulations,occupancy_type,libraryFilenameMemo)
    else:
        realAverage = readData(None, signature, SIGNATUREBASED, outputDir, jobname,occupancy_type,libraryFilenameMemo)
        if (libraryFilenameMemo is None):
            figurename = '%s_%d' % (signature, numberofMutations)
        else:
            figurename = '%s_%d_%s' % (signature, numberofMutations,libraryFilenameMemo)

        title = '%s' % (signature)
        if (numberofSimulations>0):
            listofSimulationsSignatureBased = readDataForSimulations(sample, signature, SIGNATUREBASED, outputDir, jobname,numberofSimulations,occupancy_type,libraryFilenameMemo)

    #For information
    # Is there any nan in realAverage?
    # (np.argwhere(np.isnan(realAverage))).size>0 can give the following type error
    # TypeError: ufunc 'isnan' not supported for the input types, and the inputs could not be safely coerced to any supported types according to the casting rule ''safe''
    # Therefore commented
    #use pd.isnull instead of np.isnan
    if (realAverage is not None):
        # if (np.argwhere(np.isnan(realAverage))).size>0:
        if (np.argwhere(pd.isnull(realAverage))).size > 0:
            print('Attention: There are nan in realAverage')
            print(np.argwhere(pd.isnull(realAverage)))

    if ((realAverage is not None) and (pd.notna(realAverage).any(axis=0)) and (np.any(realAverage))):
        min_list.append(np.amin(realAverage))
        max_list.append(np.amax(realAverage))

        ###########################################################################
        #Get rid of divide by zero error.
        if ((listofSimulationsSignatureBased is not None) and listofSimulationsSignatureBased):
            # averageofSimulationsSignatureBased = sumofSimulationsSignatureBased / len(listofSimulationsSignatureBased)
            stackedSimulationsSignatureBased = np.vstack(listofSimulationsSignatureBased)
            (rows, cols) = stackedSimulationsSignatureBased.shape
            # print('rows:%d cols:%d' %(rows,cols))

            #For information
            if (stackedSimulationsSignatureBased is not None):
                # if ((np.argwhere(np.isnan(stackedSimulationsSignatureBased))).size>0):
                if ((np.argwhere(pd.isnull(stackedSimulationsSignatureBased))).size > 0):
                    print('Sample:%s signature:%s' %(sample,signature))
                    print('Attention: There are nan in stackedSimulationsSignatureBased')
                    print('np.argwhere(pd.isnull(stackedSimulationsSignatureBased))')
                    print(np.argwhere(pd.isnull(stackedSimulationsSignatureBased)))

            #Note that CI_indexLow makes the lower part of the nucleosome occupancy figures with simulations.
            CI_indexLow = int(round(.05 * rows, 0))+1
            CI_indexHigh = int(round(.95 * rows, 0))-1
            # print('CI_indexLow:%d CI_indexHigh:%d' %(CI_indexLow,CI_indexHigh))

            simulationsSignatureBasedLows = []
            simulationsSignatureBasedMedians = []
            simulationsSignatureBasedHighs = []

            for col in range(cols):
                colwise_array = stackedSimulationsSignatureBased[:, col]
                sorted_colwise_array = np.sort(colwise_array)
                # print('sorted_colwise_array')
                # print(sorted_colwise_array)
                if (CI_indexLow<sorted_colwise_array.size):
                    simulationsSignatureBasedLows.append(sorted_colwise_array[CI_indexLow])
                simulationsSignatureBasedMedians.append(np.mean(sorted_colwise_array))
                # simulationsSignatureBasedMedians.append(np.median(sorted_colwise_array))
                if (CI_indexHigh<sorted_colwise_array.size):
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

        original = plt.plot(x, realAverage, color=color, label=label,linewidth=3)
        listofLegends.append(original[0])

        if (simulationsSignatureBasedMedians is not None):
            label = 'Average Simulations %s' %(label)
            simulations = plt.plot(x, simulationsSignatureBasedMedians, color='gray', linestyle=linestyle,  label=label, linewidth=3)
            listofLegends.append(simulations[0])
            if (len(simulationsSignatureBasedLows) == len(simulationsSignatureBasedHighs)):
                plt.fill_between(x, np.array(simulationsSignatureBasedLows), np.array(simulationsSignatureBasedHighs),facecolor=fillcolor)

        if (simulationsSignatureBasedLows is not None and simulationsSignatureBasedLows):
            min_list.append(np.nanmin(simulationsSignatureBasedLows))
        if (simulationsSignatureBasedHighs is not None and simulationsSignatureBasedHighs):
            max_list.append(np.nanmax(simulationsSignatureBasedHighs))

        plt.legend(loc= 'lower left', handles=listofLegends, prop={'size': 24}, shadow=False, edgecolor='white', facecolor='white',framealpha=0)

        #put the number of snps
        tobeWrittenText = "{:,}".format(numberofMutations)
        tobeWrittenText=tobeWrittenText + " " + text
        plt.text(0.99, 0.99, tobeWrittenText, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, fontsize=24)

        #put the library filename
        plt.text(0.01, 0.99, libraryFilename, verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, fontsize=24)

        #Put vertical line at x=0
        # plt.axvline(x=0, ymin=0, ymax=1, color='gray', linestyle='--')
        plt.axvline(x=0, color='gray', linestyle='--')

        # This code provides the x and y tick marks and labels
        # plt.xticks(np.arange(-1000, +1001, step=500), fontsize=30)
        plt.xticks(np.arange(-plusOrMinus, plusOrMinus+1, step=500), fontsize=30)

        #July 27, 2018
        # plt.xlim((-1000, 1000))
        plt.xlim((-plusOrMinus, plusOrMinus))

        #Let's not set ylim for HMs
        # ###################################################################################
        # min_average_nucleosome_signal= np.amin(min_list)
        # max_average_nucleosome_signal= np.amax(max_list)
        #
        # if (min_average_nucleosome_signal>= 0.55 and max_average_nucleosome_signal<=1.35):
        #     #Comment these lines to see the values out of range
        #     plt.yticks(np.arange(0.55, 1.4, step=0.1), fontsize=30)
        #     # This code overwrites the y tick labels
        #     ax.set_yticklabels(['0.55','','0.75','','0.95','','1.15','','1.35'])
        #     # This code provides some extra space
        #     plt.ylim((0.54, 1.4))
        # else:
        #     yticklabels = []
        #
        #     ymin = 0.55
        #     while (ymin>min_average_nucleosome_signal):
        #         ymin -= 0.1
        #
        #     ymax = 1.35
        #     while (ymax<max_average_nucleosome_signal):
        #         ymax += 0.1
        #
        #     ymin=round(ymin,2)
        #     ymax=round(ymax,2)
        #
        #     ytick=ymin
        #     while (ytick<=ymax):
        #         yticklabels.append(round(ytick, 2))
        #         if (ytick<ymax):
        #             yticklabels.append('')
        #         ytick = round(ytick + 0.2, 2)
        #
        #     if ymax not in yticklabels:
        #         yticklabels.append(ymax)
        #
        #     plt.yticks(np.arange(ymin,ymax, step=0.1), fontsize=30)
        #     ax.set_yticklabels(yticklabels)
        #     plt.ylim((ymin-0.01,ymax+0.01))
        # ###################################################################################

        # This code puts the tick marks
        plt.tick_params(axis='both', which='major', labelsize=30,width=3,length=10)
        plt.tick_params(axis='both', which='minor', labelsize=30,width=3,length=10)

        if (isFigureAugmentation):
            plt.title(jobname + ' ' + title, fontsize=40,fontweight='bold')
        else:
            plt.title(title, fontsize=40,fontweight='bold')

        plt.xlabel(xlabel,fontsize=32,fontweight='semibold')
        plt.ylabel(ylabel,fontsize=32,fontweight='semibold')

        filename = figurename.replace(' ', '') + figurenameEnd

        #######################################################################
        # new code
        if (sample is None):
            figureFile = os.path.join(outputDir, jobname, FIGURE, ALL, occupancy_type, filename)
        else:
            os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, occupancy_type), exist_ok=True)
            figureFile = os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, occupancy_type, filename)
        #######################################################################

        fig.savefig(figureFile)
        #Clears the axis without removing the axis itself
        plt.cla()
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
def plotAllMutationsPooledWithSimulations(xlabel,ylabel,sample,outputDir,jobname,numberofSubs,numberofIndels,numberofDinucs,numberofSimulations,mutationTypes,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus):

    if (occupancy_type==NUCLEOSOMEOCCUPANCY):
        filenameEnd='NucleosomeOccupancy'
    elif (occupancy_type==EPIGENOMICSOCCUPANCY):
        filenameEnd='EpigenomicsOccupancy'

    realAggregatedSubstitutions = None
    realAggregatedIndels = None
    realAggregatedDinucs = None

    listofSimulationsAggregatedSubstitutions = None
    listofSimulationsAggregatedIndels = None
    listofSimulationsAggregatedDinucs = None

    min_list=[]
    max_list=[]
    min_average_nucleosome_signal=0
    max_average_nucleosome_signal=0

    #######################################################################################################################
    if (sample is None):
        if libraryFilenameMemo is None:
            filename = 'Aggregated_All_Mutations_%s.png' %(filenameEnd)
        else:
            filename = 'Aggregated_All_Mutations_%s_%s.png' %(libraryFilenameMemo,filenameEnd)
        if (SBS96 in mutationTypes):
            realAggregatedSubstitutions = readData(None,None,AGGREGATEDSUBSTITUTIONS,outputDir,jobname,occupancy_type,libraryFilenameMemo)
        if (ID in mutationTypes):
            realAggregatedIndels = readData(None,None, AGGREGATEDINDELS, outputDir,jobname,occupancy_type,libraryFilenameMemo)
        if (DBS in mutationTypes):
            realAggregatedDinucs = readData(None,None,AGGREGATEDDINUCS,outputDir,jobname,occupancy_type,libraryFilenameMemo)

        if (numberofSimulations>0):
            if (SBS96 in mutationTypes):
                listofSimulationsAggregatedSubstitutions = readDataForSimulations(None,None,AGGREGATEDSUBSTITUTIONS,outputDir,jobname,numberofSimulations,occupancy_type,libraryFilenameMemo)
            if (ID in mutationTypes):
                listofSimulationsAggregatedIndels = readDataForSimulations(None,None,AGGREGATEDINDELS,outputDir,jobname,numberofSimulations,occupancy_type,libraryFilenameMemo)
            if (DBS in mutationTypes):
                listofSimulationsAggregatedDinucs = readDataForSimulations(None,None,AGGREGATEDDINUCS,outputDir,jobname,numberofSimulations,occupancy_type,libraryFilenameMemo)
    #######################################################################################################################


    #######################################################################################################################
    else:
        # filename = '%s_%s_Aggregated_Substitutions_%d_Indels_%d.png' % (sample, jobname, numberofSPMs, numberofIndels)
        if (libraryFilenameMemo is None):
            filename = '%s_Aggregated_All_Mutations_%s.png' % (sample,filenameEnd)
        else:
            filename = '%s_Aggregated_All_Mutations_%s_%s.png' % (sample,libraryFilenameMemo,filenameEnd)
        if (SBS96 in mutationTypes):
            realAggregatedSubstitutions = readData(sample,None,SAMPLEBASED_AGGREGATEDSUBSTITUTIONS,outputDir,jobname,occupancy_type,libraryFilenameMemo)
        if (ID in mutationTypes):
            realAggregatedIndels = readData(sample,None,SAMPLEBASED_AGGREGATEDINDELS,outputDir,jobname,occupancy_type,libraryFilenameMemo)
        if (DBS in mutationTypes):
            realAggregatedDinucs = readData(sample,None,SAMPLEBASED_AGGREGATEDDINUCS,outputDir,jobname,occupancy_type,libraryFilenameMemo)

        if (numberofSimulations>0):
            if (SBS96 in mutationTypes):
                listofSimulationsAggregatedSubstitutions = readDataForSimulations(sample, None, SAMPLEBASED_AGGREGATEDSUBSTITUTIONS, outputDir,jobname,numberofSimulations,occupancy_type,libraryFilenameMemo)
            if (ID in mutationTypes):
                listofSimulationsAggregatedIndels = readDataForSimulations(sample, None,SAMPLEBASED_AGGREGATEDINDELS, outputDir,jobname,numberofSimulations,occupancy_type,libraryFilenameMemo)
            if (DBS in mutationTypes):
                listofSimulationsAggregatedDinucs = readDataForSimulations(sample, None, SAMPLEBASED_AGGREGATEDDINUCS, outputDir,jobname,numberofSimulations,occupancy_type,libraryFilenameMemo)
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

    x = np.arange(-plusOrMinus, plusOrMinus+1, 1)
    # print('x.shape:%s' % x.shape)

    listofLegends = []

    if (realAggregatedSubstitutions is not None):
        min_list.append(np.amin(realAggregatedSubstitutions))
        max_list.append(np.amax(realAggregatedSubstitutions))
        aggSubs = plt.plot(x, realAggregatedSubstitutions, 'royalblue', label='Aggregated substitutions',linewidth=3,zorder=10)
        listofLegends.append(aggSubs[0])
    if (simulationsAggregatedSubstitutionsLows is not None):
        min_list.append(np.amin(simulationsAggregatedSubstitutionsLows))
    if (simulationsAggregatedSubstitutionsHighs is not None):
        max_list.append(np.amax(simulationsAggregatedSubstitutionsHighs))
    if (simulationsAggregatedSubstitutionsMedians is not None):
        simsAggSubs = plt.plot(x, simulationsAggregatedSubstitutionsMedians, color='gray',linestyle='dashed', label='Average Simulations Aggregated Substitutions',linewidth=3,zorder=10)
        listofLegends.append(simsAggSubs[0])
        plt.fill_between(x,np.array(simulationsAggregatedSubstitutionsLows),np.array(simulationsAggregatedSubstitutionsHighs),facecolor='lightblue',zorder=10)

    if (realAggregatedIndels is not None):
        min_list.append(np.amin(realAggregatedIndels))
        max_list.append(np.amax(realAggregatedIndels))
        aggIndels = plt.plot(x, realAggregatedIndels, 'darkgreen', label='Aggregated indels',linewidth=3,zorder=10)
        listofLegends.append(aggIndels[0])
    if (simulationsAggregatedIndelsLows is not None):
        min_list.append(np.amin(simulationsAggregatedIndelsLows))
    if (simulationsAggregatedIndelsHighs is not None):
        max_list.append(np.amax(simulationsAggregatedIndelsHighs))
    if simulationsAggregatedIndelsMedians is not None:
        simsAggIndels = plt.plot(x, simulationsAggregatedIndelsMedians, color='gray', linestyle='dotted',label='Average Simulations Aggregated Indels', linewidth=3,zorder=10)
        listofLegends.append(simsAggIndels[0])
        plt.fill_between(x,np.array(simulationsAggregatedIndelsLows),np.array(simulationsAggregatedIndelsHighs),facecolor='lightgreen',zorder=5)

    if (realAggregatedDinucs is not None):
        min_list.append(np.amin(realAggregatedDinucs))
        max_list.append(np.amax(realAggregatedDinucs))
        aggDinucs = plt.plot(x, realAggregatedDinucs, 'crimson', label='Aggregated dinucs',linewidth=3,zorder=10)
        listofLegends.append(aggDinucs[0])
    if (simulationsAggregatedDinucsLows is not None):
        min_list.append(np.amin(simulationsAggregatedDinucsLows))
    if (simulationsAggregatedDinucsHighs is not None):
        max_list.append(np.amax(simulationsAggregatedDinucsHighs))
    if simulationsAggregatedDinucsMedians is not None:
        simsAggDinucs = plt.plot(x, simulationsAggregatedDinucsMedians, color='gray', linestyle='dashdot',label='Average Simulations Aggregated Dinucs', linewidth=3,zorder=10)
        listofLegends.append(simsAggDinucs[0])
        plt.fill_between(x,np.array(simulationsAggregatedDinucsLows),np.array(simulationsAggregatedDinucsHighs),facecolor='lightpink',zorder=5)

    # old code
    # plt.legend(loc='lower left', prop={'size': 24},  shadow=False, edgecolor='white', facecolor ='white')
    plt.legend(loc= 'lower left',handles = listofLegends, prop={'size': 24}, shadow=False, edgecolor='white', facecolor ='white',framealpha=0)

    # put the number of subs, indels and dinucs
    text=""

    #Subs
    if numberofSubs>0:
        subs_text="{:,} subs".format(numberofSubs)
        text=subs_text

    #Indels
    if numberofIndels>0:
        indels_text = "{:,} indels".format(numberofIndels)
        if len(text)>0:
            text= text + ', ' + indels_text
        else:
            text= indels_text
    #Dinucs
    if numberofDinucs>0:
        dinucs_text = "{:,} dinucs".format(numberofDinucs)
        if len(text)>0:
            text= text + ', ' + dinucs_text
        else:
            text= dinucs_text

    plt.text(0.99, 0.99, text, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, fontsize=24)

    # put the library filename
    plt.text(0.01, 0.99, libraryFilename, verticalalignment='top', horizontalalignment='left', transform=ax.transAxes,fontsize=24)

    #Put vertical line at x=0
    plt.axvline(x=0, color='gray', linestyle='--')

    plt.tick_params(axis='both', which='major', labelsize=24)
    plt.tick_params(axis='both', which='minor', labelsize=24)

    # This code puts the tick marks
    plt.tick_params(axis='both', which='major', labelsize=24,width=3,length=10)
    plt.tick_params(axis='both', which='minor', labelsize=24,width=3,length=10)

    # This code provides the x and y tick marks and labels
    # plt.xticks(np.arange(-1000, +1001, step=500), fontsize=30)
    plt.xticks(np.arange(-plusOrMinus, plusOrMinus+1, step=500), fontsize=30)

    #July 27, 2018
    # plt.xlim((-1000, 1000))
    plt.xlim((-plusOrMinus, plusOrMinus))

    # Let's comment this for HM figures.
    # ###################################################################################
    # if (len(min_list)>0 and len(max_list)>0):
    #     ###################################################################################
    #     min_average_nucleosome_signal=np.amin(min_list)
    #     max_average_nucleosome_signal=np.amax(max_list)
    #
    #     if (min_average_nucleosome_signal >= 0.7 and max_average_nucleosome_signal <= 1.1):
    #         plt.yticks(np.arange(0.7, 1.1, step=0.05), fontsize=30)
    #         # This code overwrites the y tick labels
    #         ax.set_yticklabels(['0.7', '', '0.8', '', '0.9', '', '1.0', '', '1.1'])
    #         # This code provides some extra space
    #         plt.ylim((0.65, 1.15))
    #     else:
    #         yticklabels = []
    #
    #         ymin = 0.7
    #         while (ymin > min_average_nucleosome_signal):
    #             ymin -= 0.1
    #
    #         ymax = 1.1
    #         while (ymax < max_average_nucleosome_signal):
    #             ymax += 0.1
    #
    #         ymin = round(ymin, 2)
    #         ymax = round(ymax, 2)
    #
    #         ytick = ymin
    #         while (ytick <= ymax):
    #             yticklabels.append(round(ytick, 2))
    #             if (ytick < ymax):
    #                 yticklabels.append('')
    #             ytick = round(ytick + 0.1, 2)
    #
    #         if ymax not in yticklabels:
    #             yticklabels.append(ymax)
    #
    #         plt.yticks(np.arange(ymin, ymax, step=0.05), fontsize=30)
    #         ax.set_yticklabels(yticklabels)
    #         plt.ylim((ymin-0.01, ymax+0.01))
    #     ###################################################################################
    # else:
    #     plt.yticks(np.arange(0.7, 1.1, step=0.05), fontsize=30)
    #     # This code overwrites the y tick labels
    #     ax.set_yticklabels(['0.7', '', '0.8', '', '0.9', '', '1.0', '', '1.1'])
    #     # This code provides some extra space
    #     plt.ylim((0.65, 1.15))
    # ###################################################################################

    if (sample is not None):
        plt.title(sample, fontsize=40, fontweight='bold')
    else:
        plt.title(jobname, fontsize=40, fontweight='bold')

    plt.xlabel(xlabel, fontsize=30)
    plt.ylabel(ylabel, fontsize=30)

    ######################################################################################
    if (sample is None):
        figureFile = os.path.join(outputDir,jobname,FIGURE,ALL,occupancy_type,filename)
    else:
        os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample,occupancy_type), exist_ok=True)
        figureFile = os.path.join(outputDir,jobname,FIGURE,SAMPLES,sample,occupancy_type,filename)
    ######################################################################################

    fig.savefig(figureFile)
    plt.cla()
    plt.close(fig)
    #####################################################################
    #####################################################################
    #####################################################################

#############################################################################
############################ Plot Figure ####################################
#############################################################################


#########################################################
def checkValidness(analsesType,outputDir,jobname,occupancy_type):
    #Check whether directory exists and there are files under it.
    data_file_path = os.path.join(outputDir,jobname,DATA,occupancy_type,analsesType)

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
def plotSignatureBasedFigures(mutationType,signature2NumberofMutationsDict,sample2Signature2NumberofMutationsDict,outputDir,jobname,isFigureAugmentation,numberofSimulations,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus):
    if (occupancy_type==NUCLEOSOMEOCCUPANCY):
        ylabel = 'Average nucleosome signal'
    elif (occupancy_type==EPIGENOMICSOCCUPANCY):
        ylabel = 'Average epigenomics signal'

    if (mutationType==SBS96):
        xlabel = 'Interval around single point mutation (bp)'
        label = 'Aggregated Substitutions'
        text = 'subs'
        color = 'royalblue'
        fillcolor = 'lightblue'
        linestyle='dashed'
    elif (mutationType==ID):
        xlabel = 'Interval around indel (bp)'
        label = 'Aggregated Indels'
        text = 'indels'
        color = 'darkgreen'
        fillcolor = 'lightgreen'
        linestyle='dotted'
    elif (mutationType==DBS):
        xlabel = 'Interval around dinuc (bp)'
        label = 'Aggregated Dinucs'
        text = 'dinucs'
        color = 'crimson'
        fillcolor = 'lightpink'
        linestyle='dashdot'

    for signature in signature2NumberofMutationsDict:
        signatureBasedNumberofMutations = signature2NumberofMutationsDict[signature]
        plotSignatureBasedAverageNucleosomeOccupancyFigureWithSimulations(None, signature,
                                                                          signatureBasedNumberofMutations,
                                                                          xlabel,ylabel,label,text,
                                                                          outputDir, jobname, isFigureAugmentation,
                                                                          numberofSimulations, color,linestyle,fillcolor,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus)

    # SampleBased Subs SignatureBased Nucleosome Occupancy Figures
    for sample in sample2Signature2NumberofMutationsDict:
        for signature in sample2Signature2NumberofMutationsDict[sample]:
            sampleBasedSignatureBasedNumberofMutations = sample2Signature2NumberofMutationsDict[sample][signature]
            plotSignatureBasedAverageNucleosomeOccupancyFigureWithSimulations(sample, signature,
                                                                              sampleBasedSignatureBasedNumberofMutations,
                                                                              xlabel,ylabel,label,text,
                                                                              outputDir, jobname, isFigureAugmentation,
                                                                              numberofSimulations,color,linestyle,fillcolor,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus)
#########################################################

#########################################################
def occupancyAverageSignalFigures(outputDir,jobname,figureAugmentation,numberofSimulations,sample_based,mutationTypes,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus):
    if (occupancy_type==NUCLEOSOMEOCCUPANCY):
        ylabel='Average nucleosome signal'
    elif (occupancy_type==EPIGENOMICSOCCUPANCY):
        ylabel='Average epigenomics signal'

    #######################################################################################################################
    os.makedirs(os.path.join(outputDir, jobname, FIGURE, ALL, occupancy_type), exist_ok=True)
    #######################################################################################################################

    isFigureAugmentation = False
    if (figureAugmentation=='augmentation'):
        isFigureAugmentation = True

    ############## Read necessary dictionaries starts ########################################
    mutationType2NumberofMutationsDict = getDictionary(outputDir, jobname, MutationType2NumberofMutatiosDictFilename)

    subsSignature2NumberofMutationsDict = getSubsSignature2NumberofMutationsDict(outputDir,jobname)
    indelsSignature2NumberofMutationsDict = getIndelsSignature2NumberofMutationsDict(outputDir,jobname)
    dinucsSignature2NumberofMutationsDict = getDictionary(outputDir, jobname,DinucsSignature2NumberofMutationsDictFilename)

    if sample_based:
        sample2NumberofSubsDict = getSample2NumberofSubsDict(outputDir, jobname)
        sample2NumberofIndelsDict = getSample2NumberofIndelsDict(outputDir, jobname)
        sample2NumberofDinucsDict = getDictionary(outputDir, jobname, Sample2NumberofDinucsDictFilename)

        sample2SubsSignature2NumberofMutationsDict = getSample2SubsSignature2NumberofMutationsDict(outputDir,jobname)
        sample2IndelsSignature2NumberofMutationsDict = getSample2IndelsSignature2NumberofMutationsDict(outputDir, jobname)
        sample2DinucsSignature2NumberofMutationsDict = getDictionary(outputDir, jobname,Sample2DinucsSignature2NumberofMutationsDictFilename)
    else:
        sample2NumberofSubsDict = {}
        sample2NumberofIndelsDict = {}
        sample2NumberofDinucsDict = {}

        sample2SubsSignature2NumberofMutationsDict = {}
        sample2IndelsSignature2NumberofMutationsDict = {}
        sample2DinucsSignature2NumberofMutationsDict = {}
    ############## Read necessary dictionaries ends ##########################################


    ##############################################################
    numberofSubs=0
    numberofIndels=0
    numberofDinucs=0

    if (SUBS in mutationType2NumberofMutationsDict):
        numberofSubs=mutationType2NumberofMutationsDict[SUBS]
    if (INDELS in mutationType2NumberofMutationsDict):
        numberofIndels = mutationType2NumberofMutationsDict[INDELS]
    if (DINUCS in mutationType2NumberofMutationsDict):
        numberofDinucs = mutationType2NumberofMutationsDict[DINUCS]

    #tissue based
    plotAllMutationsPooledWithSimulations('Interval around variant (bp)',ylabel,None,outputDir,jobname,numberofSubs,numberofIndels,numberofDinucs,numberofSimulations,mutationTypes,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus)
    ##############################################################

    #############################################################################################################################################
    #Plot Signature Based
    #ncomms11383 Fig3b signature based average nucleosome occupancy figures
    if checkValidness(SIGNATUREBASED,outputDir,jobname,occupancy_type):
        if (SBS96 in mutationTypes):
            #Subs Signatures
            plotSignatureBasedFigures(SBS96,subsSignature2NumberofMutationsDict,sample2SubsSignature2NumberofMutationsDict,outputDir,jobname,isFigureAugmentation,numberofSimulations,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus)
        if (ID in mutationTypes):
            #Indels Signatures
            plotSignatureBasedFigures(ID,indelsSignature2NumberofMutationsDict,sample2IndelsSignature2NumberofMutationsDict,outputDir,jobname,isFigureAugmentation,numberofSimulations,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus)
        if (DBS in mutationTypes):
            # Dinucs Signatures
            plotSignatureBasedFigures(DBS,dinucsSignature2NumberofMutationsDict,sample2DinucsSignature2NumberofMutationsDict,outputDir,jobname,isFigureAugmentation,numberofSimulations,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus)
    #############################################################################################################################################


    ##############################################################
    if sample_based:
        #ALL SAMPLES IN ONE
        #Plot "all samples pooled" and "sample based" signature based in one figure
        plotAllSamplesPooledAndSampleBasedSignaturesFiguresInOneFigure(subsSignature2NumberofMutationsDict,sample2SubsSignature2NumberofMutationsDict,outputDir,jobname,'royalblue','Interval around single point mutation (bp)',ylabel,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus)
        plotAllSamplesPooledAndSampleBasedSignaturesFiguresInOneFigure(indelsSignature2NumberofMutationsDict,sample2IndelsSignature2NumberofMutationsDict,outputDir,jobname,'darkgreen','Interval around indel (bp)',ylabel,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus)
        plotAllSamplesPooledAndSampleBasedSignaturesFiguresInOneFigure(dinucsSignature2NumberofMutationsDict,sample2DinucsSignature2NumberofMutationsDict,outputDir,jobname,'crimson','Interval around indel (bp)',ylabel,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus)

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

            #sample based
            plotAllMutationsPooledWithSimulations('Interval around variant (bp)',ylabel,sample,outputDir,jobname,numberofSubs, numberofIndels,numberofDinucs,numberofSimulations,mutationTypes,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus)
    ##############################################################
