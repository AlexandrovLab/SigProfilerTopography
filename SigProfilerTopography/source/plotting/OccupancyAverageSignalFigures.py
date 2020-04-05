# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018 Burcak Otlu

import os
import numpy as np
import pandas as pd
import multiprocessing

import matplotlib
BACKEND = 'Agg'
if matplotlib.get_backend().lower() != BACKEND.lower():
    # If backend is not set properly a call to describe will hang
    matplotlib.use(BACKEND)

from matplotlib import pyplot as plt

from SigProfilerTopography.source.commons.TopographyCommons import SBS96

from SigProfilerTopography.source.commons.TopographyCommons import SBS
from SigProfilerTopography.source.commons.TopographyCommons import ID
from SigProfilerTopography.source.commons.TopographyCommons import DBS

from SigProfilerTopography.source.commons.TopographyCommons import SUBS
from SigProfilerTopography.source.commons.TopographyCommons import INDELS
from SigProfilerTopography.source.commons.TopographyCommons import DINUCS

from SigProfilerTopography.source.commons.TopographyCommons import DATA
from SigProfilerTopography.source.commons.TopographyCommons import FIGURE
from SigProfilerTopography.source.commons.TopographyCommons import ALL
from SigProfilerTopography.source.commons.TopographyCommons import HEATMAPS

from SigProfilerTopography.source.commons.TopographyCommons import SIGNATUREBASED
from SigProfilerTopography.source.commons.TopographyCommons import AGGREGATEDSUBSTITUTIONS
from SigProfilerTopography.source.commons.TopographyCommons import AGGREGATEDINDELS
from SigProfilerTopography.source.commons.TopographyCommons import AGGREGATEDDINUCS

from SigProfilerTopography.source.commons.TopographyCommons import SAMPLEBASED_SIGNATUREBASED
from SigProfilerTopography.source.commons.TopographyCommons import SAMPLEBASED_AGGREGATEDSUBSTITUTIONS
from SigProfilerTopography.source.commons.TopographyCommons import SAMPLEBASED_AGGREGATEDINDELS
from SigProfilerTopography.source.commons.TopographyCommons import SAMPLEBASED_AGGREGATEDDINUCS
from SigProfilerTopography.source.commons.TopographyCommons import SAMPLES

from SigProfilerTopography.source.commons.TopographyCommons import NUCLEOSOMEOCCUPANCY
from SigProfilerTopography.source.commons.TopographyCommons import EPIGENOMICSOCCUPANCY

from SigProfilerTopography.source.commons.TopographyCommons import takeAverage
from SigProfilerTopography.source.commons.TopographyCommons import getDictionary

from SigProfilerTopography.source.commons.TopographyCommons import Table_MutationType_NumberofMutations_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Table_SubsSignature_Cutoff_NumberofMutations_AverageProbability_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Table_IndelsSignature_Cutoff_NumberofMutations_AverageProbability_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Table_DinucsSignature_Cutoff_NumberofMutations_AverageProbability_Filename

from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofSubsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofIndelsDict
from SigProfilerTopography.source.commons.TopographyCommons import Sample2NumberofDinucsDictFilename

from SigProfilerTopography.source.commons.TopographyCommons import getSample2SubsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2IndelsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import Sample2DinucsSignature2NumberofMutationsDictFilename

from SigProfilerTopography.source.commons.TopographyCommons import BIOSAMPLE_UNDECLARED
from SigProfilerTopography.source.commons.TopographyCommons import natural_key

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
def plotAllSamplesPooledAndSampleBasedSignaturesFiguresInOneFigure(signature_cutoff_numberofmutations_averageprobability_df,sample2Signature2NumberofMutationsDict,outputDir,jobname,color,xlabel,ylabel,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus):

    if (occupancy_type==NUCLEOSOMEOCCUPANCY):
        filenameEnd = 'NucleosomeOccupancy'
    elif (occupancy_type==EPIGENOMICSOCCUPANCY):
        filenameEnd = 'EpigenomicsOccupancy'

    for signature in signature_cutoff_numberofmutations_averageprobability_df['signature'].unique():
        min_list = []
        max_list = []

        label2NumpyArrayDict = {}
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
        plt.clf()
        plt.cla()
        plt.close(fig)
        ####################################### plotting ends #######################################

#############################################################################
########################## Plot Figure ends  ################################
#############################################################################


#############################################################################
########################## Plot Figure starts  ##############################
#############################################################################
#Called by plotSignatureBasedFigures
def plotSignatureBasedAverageNucleosomeOccupancyFigureWithSimulations(sample,signature,numberofMutations,xlabel,ylabel,label,text,outputDir,jobname,isFigureAugmentation,numberofSimulations,color,linestyle,fillcolor,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus,verbose):

    if (occupancy_type==NUCLEOSOMEOCCUPANCY):
        figurenameEnd='_NucleosomeOccupancy.png'
    elif (occupancy_type==EPIGENOMICSOCCUPANCY):
        figurenameEnd='_EpigenomicsOccupancy.png'

    min_list=[]
    max_list=[]
    min_average_nucleosome_signal=0
    max_average_nucleosome_signal=0

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
            print('Attention: There are %d nans in realAverage in %s for %s' %(len(np.argwhere(np.isnan(realAverage))),libraryFilenameMemo,signature))
            if verbose: print('\tVerbose %s' %(np.argwhere(pd.isnull(realAverage))))

    if ((realAverage is not None) and (pd.notna(realAverage).any(axis=0)) and (np.any(realAverage))):
        min_list.append(np.amin(realAverage))
        max_list.append(np.amax(realAverage))

        #95%CI
        simulationsSignatureBasedLows, simulationsSignatureBasedMeans, simulationsSignatureBasedHighs = takeAverage(listofSimulationsSignatureBased)

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

        if (simulationsSignatureBasedMeans is not None):
            label = 'Average Simulations %s' %(label)
            simulations = plt.plot(x, simulationsSignatureBasedMeans, color='gray', linestyle=linestyle,  label=label, linewidth=3)
            listofLegends.append(simulations[0])
            if (simulationsSignatureBasedLows is not None) and (simulationsSignatureBasedHighs is not None):
                plt.fill_between(x, np.array(simulationsSignatureBasedLows), np.array(simulationsSignatureBasedHighs),facecolor=fillcolor)

        if ((simulationsSignatureBasedLows is not None) and (not np.all(np.isnan(simulationsSignatureBasedLows)))):
            min_list.append(np.nanmin(simulationsSignatureBasedLows))
        if ((simulationsSignatureBasedHighs is not None) and (not np.all(np.isnan(simulationsSignatureBasedHighs)))):
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
        plt.clf()
        #Clears the axis without removing the axis itself
        plt.cla()
        plt.close()
#############################################################################
########################## Plot Figure starts  ##############################
#############################################################################




#############################################################################
############################ Plot Figure ####################################
#############################################################################
def plotAllMutationsPooledWithSimulations(xlabel,ylabel,sample,outputDir,jobname,numberofSubs,numberofIndels,numberofDinucs,numberofSimulations,mutationType,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus):
    if mutationType==SBS96:
        to_be_added_to_the_filename='SBS%s' %(mutationType)
    else:
        to_be_added_to_the_filename=mutationType

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

    min_average_nucleosome_signal=0
    max_average_nucleosome_signal=0

    #######################################################################################################################
    if (sample is None):
        if libraryFilenameMemo is None:
            filename = 'Aggregated_All_Mutations_%s_%s.png' %(to_be_added_to_the_filename,filenameEnd)
        else:
            filename = 'Aggregated_All_Mutations_%s_%s_%s.png' %(to_be_added_to_the_filename,libraryFilenameMemo,filenameEnd)

        if (SBS96==mutationType):
            realAggregatedSubstitutions = readData(None,None,AGGREGATEDSUBSTITUTIONS,outputDir,jobname,occupancy_type,libraryFilenameMemo)
        if (ID==mutationType):
            realAggregatedIndels = readData(None,None, AGGREGATEDINDELS, outputDir,jobname,occupancy_type,libraryFilenameMemo)
        if (DBS==mutationType):
            realAggregatedDinucs = readData(None,None,AGGREGATEDDINUCS,outputDir,jobname,occupancy_type,libraryFilenameMemo)

        if (numberofSimulations>0):
            if (SBS96==mutationType):
                listofSimulationsAggregatedSubstitutions = readDataForSimulations(None,None,AGGREGATEDSUBSTITUTIONS,outputDir,jobname,numberofSimulations,occupancy_type,libraryFilenameMemo)
            if (ID==mutationType):
                listofSimulationsAggregatedIndels = readDataForSimulations(None,None,AGGREGATEDINDELS,outputDir,jobname,numberofSimulations,occupancy_type,libraryFilenameMemo)
            if (DBS==mutationType):
                listofSimulationsAggregatedDinucs = readDataForSimulations(None,None,AGGREGATEDDINUCS,outputDir,jobname,numberofSimulations,occupancy_type,libraryFilenameMemo)
    #######################################################################################################################


    #######################################################################################################################
    else:
        # filename = '%s_%s_Aggregated_Substitutions_%d_Indels_%d.png' % (sample, jobname, numberofSPMs, numberofIndels)
        if (libraryFilenameMemo is None):
            filename = '%s_Aggregated_All_Mutations_%s_%s.png' % (sample,to_be_added_to_the_filename,filenameEnd)
        else:
            filename = '%s_Aggregated_All_Mutations_%s_%s_%s.png' % (sample,to_be_added_to_the_filename,libraryFilenameMemo,filenameEnd)

        if (SBS96 == mutationType):
            realAggregatedSubstitutions = readData(sample,None,SAMPLEBASED_AGGREGATEDSUBSTITUTIONS,outputDir,jobname,occupancy_type,libraryFilenameMemo)
        if (ID == mutationType):
            realAggregatedIndels = readData(sample,None,SAMPLEBASED_AGGREGATEDINDELS,outputDir,jobname,occupancy_type,libraryFilenameMemo)
        if (DBS == mutationType):
            realAggregatedDinucs = readData(sample,None,SAMPLEBASED_AGGREGATEDDINUCS,outputDir,jobname,occupancy_type,libraryFilenameMemo)

        if (numberofSimulations>0):
            if (SBS96 == mutationType):
                listofSimulationsAggregatedSubstitutions = readDataForSimulations(sample, None, SAMPLEBASED_AGGREGATEDSUBSTITUTIONS, outputDir,jobname,numberofSimulations,occupancy_type,libraryFilenameMemo)
            if (ID == mutationType):
                listofSimulationsAggregatedIndels = readDataForSimulations(sample, None,SAMPLEBASED_AGGREGATEDINDELS, outputDir,jobname,numberofSimulations,occupancy_type,libraryFilenameMemo)
            if (DBS == mutationType):
                listofSimulationsAggregatedDinucs = readDataForSimulations(sample, None, SAMPLEBASED_AGGREGATEDDINUCS, outputDir,jobname,numberofSimulations,occupancy_type,libraryFilenameMemo)
    #######################################################################################################################

    #####################################################################
    #95%CI
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

    ##############################################################
    if (realAggregatedSubstitutions is not None):
        aggSubs = plt.plot(x, realAggregatedSubstitutions, 'royalblue', label='Aggregated substitutions',linewidth=3,zorder=10)
        listofLegends.append(aggSubs[0])

    if (simulationsAggregatedSubstitutionsMedians is not None):
        simsAggSubs = plt.plot(x, simulationsAggregatedSubstitutionsMedians, color='gray',linestyle='dashed', label='Average Simulations Aggregated Substitutions',linewidth=3,zorder=10)
        listofLegends.append(simsAggSubs[0])
        if (simulationsAggregatedSubstitutionsLows is not None) and (simulationsAggregatedSubstitutionsHighs is not None):
            plt.fill_between(x,np.array(simulationsAggregatedSubstitutionsLows),np.array(simulationsAggregatedSubstitutionsHighs),facecolor='lightblue',zorder=10)
    ##############################################################

    ##############################################################
    if (realAggregatedIndels is not None):
        aggIndels = plt.plot(x, realAggregatedIndels, 'darkgreen', label='Aggregated indels',linewidth=3,zorder=10)
        listofLegends.append(aggIndels[0])

    if simulationsAggregatedIndelsMedians is not None:
        simsAggIndels = plt.plot(x, simulationsAggregatedIndelsMedians, color='gray', linestyle='dashed',label='Average Simulations Aggregated Indels', linewidth=3,zorder=10)
        listofLegends.append(simsAggIndels[0])
        if (simulationsAggregatedIndelsLows is not None) and (simulationsAggregatedIndelsHighs is not None):
            plt.fill_between(x,np.array(simulationsAggregatedIndelsLows),np.array(simulationsAggregatedIndelsHighs),facecolor='lightgreen',zorder=5)
    ##############################################################

    ##############################################################
    if (realAggregatedDinucs is not None):
        aggDinucs = plt.plot(x, realAggregatedDinucs, 'crimson', label='Aggregated dinucs',linewidth=3,zorder=10)
        listofLegends.append(aggDinucs[0])

    if simulationsAggregatedDinucsMedians is not None:
        simsAggDinucs = plt.plot(x, simulationsAggregatedDinucsMedians, color='gray', linestyle='dashed',label='Average Simulations Aggregated Dinucs', linewidth=3,zorder=10)
        listofLegends.append(simsAggDinucs[0])
        if (simulationsAggregatedDinucsLows is not None) and (simulationsAggregatedDinucsHighs is not None):
            plt.fill_between(x,np.array(simulationsAggregatedDinucsLows),np.array(simulationsAggregatedDinucsHighs),facecolor='lightpink',zorder=5)
    ##############################################################

    plt.legend(loc= 'lower left',handles = listofLegends, prop={'size': 24}, shadow=False, edgecolor='white', facecolor ='white',framealpha=0)

    #######################################################################
    # put the number of subs, indels and dinucs
    text=""

    #Subs
    if (mutationType==SBS96) and (numberofSubs>0):
        subs_text="{:,} subs".format(numberofSubs)
        text=subs_text

    #Indels
    if (mutationType==ID) and (numberofIndels>0):
        indels_text = "{:,} indels".format(numberofIndels)
        if len(text)>0:
            text= text + ', ' + indels_text
        else:
            text= indels_text
    #Dinucs
    if (mutationType==DBS) and (numberofDinucs>0):
        dinucs_text = "{:,} dinucs".format(numberofDinucs)
        if len(text)>0:
            text= text + ', ' + dinucs_text
        else:
            text= dinucs_text

    #put number of mutations
    plt.text(0.99, 0.99, text, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, fontsize=24)
    #######################################################################

    #put the library filename
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
        os.makedirs(os.path.join(outputDir, jobname,FIGURE,SAMPLES,sample,occupancy_type), exist_ok=True)
        figureFile = os.path.join(outputDir,jobname,FIGURE,SAMPLES,sample,occupancy_type,filename)
    ######################################################################################

    fig.savefig(figureFile)
    plt.clf()
    plt.cla()
    plt.close()
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
        return True
        #Takes very long time especially in epigenomics
        # filenames = [f for f in os.listdir(data_file_path) if os.path.isfile(os.path.join(data_file_path, f))]
        # print('DEBUG len(filenames):%d' %(len(filenames)))
        # if (len(filenames)>0):
        #     return True
        # else:
        #     return False
    else:
        return False
#########################################################


#########################################################
def plotSignatureBasedFigures(mutationType,signature_cutoff_numberofmutations_averageprobability_df,sample2Signature2NumberofMutationsDict,outputDir,jobname,isFigureAugmentation,numberofSimulations,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus,verbose):
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
        linestyle='dashed'
        # linestyle='dotted'
    elif (mutationType==DBS):
        xlabel = 'Interval around dinuc (bp)'
        label = 'Aggregated Dinucs'
        text = 'dinucs'
        color = 'crimson'
        fillcolor = 'lightpink'
        linestyle='dashed'
        # linestyle='dashdot'

    for signature in signature_cutoff_numberofmutations_averageprobability_df['signature'].unique():
        #[signature cutoff numberofMutations  averageProbability]
        signatureBasedNumberofMutations = int(signature_cutoff_numberofmutations_averageprobability_df[signature_cutoff_numberofmutations_averageprobability_df['signature']==signature]['number_of_mutations'].values[0])
        plotSignatureBasedAverageNucleosomeOccupancyFigureWithSimulations(None, signature,
                                                                          signatureBasedNumberofMutations,
                                                                          xlabel,ylabel,label,text,
                                                                          outputDir, jobname, isFigureAugmentation,
                                                                          numberofSimulations, color,linestyle,fillcolor,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus,verbose)

    # SampleBased Subs SignatureBased Nucleosome Occupancy Figures
    for sample in sample2Signature2NumberofMutationsDict:
        for signature in sample2Signature2NumberofMutationsDict[sample]:
            sampleBasedSignatureBasedNumberofMutations = sample2Signature2NumberofMutationsDict[sample][signature]
            plotSignatureBasedAverageNucleosomeOccupancyFigureWithSimulations(sample, signature,
                                                                              sampleBasedSignatureBasedNumberofMutations,
                                                                              xlabel,ylabel,label,text,
                                                                              outputDir, jobname, isFigureAugmentation,
                                                                              numberofSimulations,color,linestyle,fillcolor,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus,verbose)
#########################################################



########################################################
def updateDictWithTwoLevels(signature2ENCODEHM2FoldChangeDict, signature, hm, fold_change):
    if signature in signature2ENCODEHM2FoldChangeDict:
        if hm in signature2ENCODEHM2FoldChangeDict[signature]:
            print('There is a problem')
        else:
            signature2ENCODEHM2FoldChangeDict[signature][hm]=fold_change
    else:
        signature2ENCODEHM2FoldChangeDict[signature] = {}
        signature2ENCODEHM2FoldChangeDict[signature][hm] = fold_change
########################################################



########################################################
def updateDictWithThreeLevels(signature2Biosample2ENCODEHM2FoldChangeDict, signature, biosample, hm, fold_change):
    if signature in signature2Biosample2ENCODEHM2FoldChangeDict:
        if biosample in signature2Biosample2ENCODEHM2FoldChangeDict[signature]:
            if hm in signature2Biosample2ENCODEHM2FoldChangeDict[signature][biosample]:
                print('There is a problem')
            else:
                signature2Biosample2ENCODEHM2FoldChangeDict[signature][biosample][hm]=fold_change
        else:
            signature2Biosample2ENCODEHM2FoldChangeDict[signature][biosample]={}
            signature2Biosample2ENCODEHM2FoldChangeDict[signature][biosample][hm]=fold_change
    else:
        signature2Biosample2ENCODEHM2FoldChangeDict[signature] = {}
        signature2Biosample2ENCODEHM2FoldChangeDict[signature][biosample] = {}
        signature2Biosample2ENCODEHM2FoldChangeDict[signature][biosample][hm] = fold_change
#######################################################



########################################################
def calculate_fold_change(output_dir,numberofSimulations,signature,cancer_type,dna_element,plusOrMinus_epigenomics,verbose):
    occupancy_type=EPIGENOMICSOCCUPANCY
    center=plusOrMinus_epigenomics
    plusorMinus=250
    start=center-plusorMinus
    end= center+plusorMinus+1

    if verbose: print('\n\tVerbose ----------> %s %s %s' % (signature, cancer_type, dna_element))
    avg_real_signal = None

    # SBS1_sim1_ENCFF330CCJ_osteoblast_H3K79me2-human_AverageSignalArray.txt
    avg_real_data_signal_array = readData(None, signature, SIGNATUREBASED, output_dir, cancer_type, occupancy_type,dna_element)

    if avg_real_data_signal_array is not None:
        #If there is nan in the list np.mean returns nan.
        # 1st way
        # avg_real_data_signal_array[np.isnan(avg_real_data_signal_array)] = 0
        # avg_real_signal = np.mean(avg_real_data_signal_array[1750:2251])
        # 2nd way
        if not np.all(np.isnan(avg_real_data_signal_array[start:end])):
            avg_real_signal = np.nanmean(avg_real_data_signal_array[start:end])

    avg_simulated_signal = None
    if (numberofSimulations > 0):
        listofSimulationsSignatureBased = readDataForSimulations(None, signature, SIGNATUREBASED, output_dir,cancer_type, numberofSimulations, occupancy_type,dna_element)

        if ((listofSimulationsSignatureBased is not None) and listofSimulationsSignatureBased):
            stackedSimulationsSignatureBased = np.vstack(listofSimulationsSignatureBased)
            (rows, cols) = stackedSimulationsSignatureBased.shape
            if verbose: print('\tVerbose After np.vstack --- stackedSimulationsSignatureBased rows:%d cols:%d' % (rows, cols))

            #Get the part that we are interested
            stackedSimulationsSignatureBased_of_interest= stackedSimulationsSignatureBased[:,start:end]
            (rows, cols) = stackedSimulationsSignatureBased_of_interest.shape
            if verbose: print(
                '\tVerbose After get part interested in --- stackedSimulationsSignatureBased_of_interest rows:%d cols:%d' % (rows, cols))

            # Get rid of all the rows with all nans
            stackedSimulationsSignatureBased_of_interest = stackedSimulationsSignatureBased_of_interest[~np.isnan(stackedSimulationsSignatureBased_of_interest).all(axis=1)]
            (rows, cols) = stackedSimulationsSignatureBased_of_interest.shape
            if verbose: print(
                '\tVerbose After removing all nan rows --- stackedSimulationsSignatureBased_of_interest rows:%d cols:%d' % (rows, cols))

            simulationsSignatureBasedMeans = np.nanmean(stackedSimulationsSignatureBased_of_interest, axis=1)
            if not np.all(np.isnan(simulationsSignatureBasedMeans)):
                avg_simulated_signal = np.nanmean(simulationsSignatureBasedMeans)

    if (avg_real_signal is not None) and (avg_simulated_signal is not None):
        fold_change = avg_real_signal / avg_simulated_signal
        if verbose: print('\tVerbose ----------> %s %s %s avg_real_signal:%f\tavg_simulated_signal:%f\tfold change:%f' % (signature, cancer_type, dna_element,avg_real_signal, avg_simulated_signal, fold_change))
        return fold_change
    else:
        if verbose: print('\tVerbose ----------> %s %s %s avg_real_signal:%s\tavg_simulated_signal:%s\tfold change: nan' % (signature,cancer_type,dna_element,str(avg_real_signal),str(avg_simulated_signal)))
        return None
########################################################



########################################################
#Constraint epigenomic biosample name must be contained in epigenomic file memo
def fill_signature2CancerType2Biosample2HM2FoldChangeDict(output_dir,numberofSimulations,signature,signature_type,cancer_type,epigenomics_files_memos,epigenomics_biosamples,epigenomics_dna_elements,occupancy_type,plusOrMinus_epigenomics,verbose):

    list=[]

    #Step1
    signature2EpigenomicsFileMemo2FoldChangeDict = {}

    #Step2
    signature2Biosample2EpigenomicsFileMemo2FoldChangeDict = {}

    #Step3
    signature2Biosample2DNAElementPooled2AverageFoldChangeDict = {}

    #Step4
    signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict = {}

    ########################################################################
    #Step1
    if verbose: print('\tVerbose ===================================================================================================')
    if verbose: print('\tVerbose signature:%s' %(signature))

    # epigenomics_files_memos have been used in naming
    # ID8_sim56_H3K27me3_Breast_Epithelium_AccumulatedCountArray.txt
    # ID8_sim56_H3K27me3_Breast_Epithelium_AccumulatedSignalArray.txt

    if verbose: print('\tVerbose #####################################################')
    if verbose: print('\tVerbose DNA Elements %s' %(epigenomics_files_memos))

    for epigenomics_files_memo in epigenomics_files_memos:
        fold_change=calculate_fold_change(output_dir,numberofSimulations,signature,cancer_type,epigenomics_files_memo,plusOrMinus_epigenomics,verbose)
        if fold_change is not None:
            updateDictWithTwoLevels(signature2EpigenomicsFileMemo2FoldChangeDict,signature,epigenomics_files_memo,fold_change)

    if verbose: print('\tVerbose #######################################################')
    if verbose: print('\tVerbose Step1 After %s' % (signature))
    if verbose: print('\tVerbose signature2EpigenomicsFileMemo2FoldChangeDict')
    if verbose: print('\tVerbose %s' %signature2EpigenomicsFileMemo2FoldChangeDict)
    ########################################################################

    ########################################################################
    #Step2
    if epigenomics_biosamples is not None:
        for biosample in epigenomics_biosamples:
            for signature in signature2EpigenomicsFileMemo2FoldChangeDict:
                for epigenomics_file_memo in signature2EpigenomicsFileMemo2FoldChangeDict[signature]:
                    if (biosample==BIOSAMPLE_UNDECLARED):
                        fold_change=signature2EpigenomicsFileMemo2FoldChangeDict[signature][epigenomics_file_memo]
                        updateDictWithThreeLevels(signature2Biosample2EpigenomicsFileMemo2FoldChangeDict, signature, BIOSAMPLE_UNDECLARED, epigenomics_file_memo, fold_change)
                    else:
                        if (biosample.lower() in epigenomics_file_memo.lower()):
                            fold_change=signature2EpigenomicsFileMemo2FoldChangeDict[signature][epigenomics_file_memo]
                            updateDictWithThreeLevels(signature2Biosample2EpigenomicsFileMemo2FoldChangeDict, signature, biosample, epigenomics_file_memo, fold_change)

    if verbose: print('\tVerbose #######################################################')
    if verbose: print('\tVerbose Step2 After %s' % (signature))
    if verbose: print('\tVerbose signature2Biosample2EpigenomicsFileMemo2FoldChangeDict')
    if verbose: print('\tVerbose %s' %signature2Biosample2EpigenomicsFileMemo2FoldChangeDict)
    ########################################################################

    ########################################################
    #Step3 Pool for DNA elements
    #Will be filled in Step3
    if epigenomics_dna_elements is not None:
        for signature in signature2Biosample2EpigenomicsFileMemo2FoldChangeDict:
            for biosample in signature2Biosample2EpigenomicsFileMemo2FoldChangeDict[signature]:
                #Pool for each dna_element in epigenomics_dna_elements
                for epigenomics_dna_element in epigenomics_dna_elements:
                    signature_biosample_dna_element_fold_change_list = []
                    for signature_biosample_dna_element in signature2Biosample2EpigenomicsFileMemo2FoldChangeDict[signature][biosample]:
                        if (epigenomics_dna_element.lower() in signature_biosample_dna_element.lower()):
                            signature_biosample_dna_element_fold_change_list.append(signature2Biosample2EpigenomicsFileMemo2FoldChangeDict[signature][biosample][signature_biosample_dna_element])
                    if len(signature_biosample_dna_element_fold_change_list)>0:
                        average_fold_change=np.nanmean(signature_biosample_dna_element_fold_change_list)
                    else:
                        average_fold_change=np.nan
                    #update
                    # signature2Biosample2DNAElementPooled2AverageFoldChangeDict[signature][biosample][epigenomics_dna_element]=average_fold_change
                    updateDictWithThreeLevels(signature2Biosample2DNAElementPooled2AverageFoldChangeDict, signature, biosample,epigenomics_dna_element, average_fold_change)
    ########################################################


    if verbose: print('\tVerbose #######################################################')
    if verbose: print('\tVerbose Step3 After %s' % (signature))
    if verbose: print('\tVerbose signature2Biosample2DNAElementPooled2AverageFoldChangeDict')
    if verbose: print('\tVerbose %s' %signature2Biosample2DNAElementPooled2AverageFoldChangeDict)
    ########################################################################

    ########################################################
    #Step4 Pool for Biosamples
    #Will be filled in Step4
    if epigenomics_biosamples is not None:
        for signature in signature2Biosample2DNAElementPooled2AverageFoldChangeDict:
            signature_based_biosamples=signature2Biosample2DNAElementPooled2AverageFoldChangeDict[signature]
            for epigenomics_dna_element in epigenomics_dna_elements:
                signature_dna_element_fold_change_list=[]
                for biosample in signature_based_biosamples:
                    signature_dna_element_fold_change_list.append(signature2Biosample2DNAElementPooled2AverageFoldChangeDict[signature][biosample][epigenomics_dna_element])
                if (len(signature_dna_element_fold_change_list)>0 and not np.all(np.isnan(signature_dna_element_fold_change_list))):
                    average_fold_change=np.nanmean(signature_dna_element_fold_change_list)
                else:
                    average_fold_change=np.nan
                #update
                # signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict[signature][epigenomics_dna_element]=average_fold_change
                updateDictWithTwoLevels(signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict, signature, epigenomics_dna_element,average_fold_change)
    ########################################################

    if verbose: print('\tVerbose #######################################################')
    if verbose: print('\tVerbose Step4 After %s' % (signature))
    if verbose: print('\tVerbose signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict')
    if verbose: print('\tVerbose %s' %signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict)
    ########################################################################

    list.append(output_dir)
    list.append(signature)
    list.append(signature_type)
    list.append(cancer_type)
    list.append(occupancy_type)
    list.append(signature2Biosample2EpigenomicsFileMemo2FoldChangeDict)
    list.append(signature2Biosample2DNAElementPooled2AverageFoldChangeDict)
    list.append(signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict)
    list.append(verbose)

    return list
########################################################


########################################################
#Plot heatmaps
# For Step4
def fill_average_fold_change_array_rows_signatures_columns_dna_elements(signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict):
    signatures=[]
    dna_elements=[]

    for signature in signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict:
        if signature not in signatures:
            signatures.append(signature)
        for dna_element in signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict[signature]:
            if dna_element not in dna_elements:
                dna_elements.append(dna_element)

    #sort the dna_elements and signatures
    dna_elements=sorted(dna_elements,key=natural_key)
    signatures=sorted(signatures,key=natural_key)

    #Initialize
    average_fold_change_array = np.zeros((len(signatures), len(dna_elements)))

    #fill the average_fold_change_array
    for signature_index, signature in enumerate(signatures,0):
        for dna_element_index, dna_element in enumerate(dna_elements,0):
            if signature in signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict:
                if dna_element in signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict[signature]:
                    average_fold_change_array[signature_index,dna_element_index]=signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict[signature][dna_element]
                else:
                    average_fold_change_array[signature_index,dna_element_index]=np.nan
            else:
                average_fold_change_array[signature_index, dna_element_index] = np.nan

    return signatures, dna_elements, average_fold_change_array
########################################################


########################################################
#Plot heatmaps
# For Step3
def fill_average_fold_change_array_rows_biosamples_columns_dna_elements(biosample2DNAElementPooled2AverageFoldChangeDict):
    biosamples=[]
    dna_elements=[]

    for biosample in biosample2DNAElementPooled2AverageFoldChangeDict:
        biosamples.append(biosample)
        for dna_element in biosample2DNAElementPooled2AverageFoldChangeDict[biosample]:
            if dna_element not in dna_elements:
                dna_elements.append(dna_element)

    # sort the biosamples
    biosamples = sorted(biosamples,key=natural_key)

    #sort the hms
    dna_elements=sorted(dna_elements,key=natural_key)

    #Initialize
    average_fold_change_array = np.zeros((len(biosamples), len(dna_elements)))

    for biosample_index,biosample in enumerate(biosamples,0):
        for dna_element_index,dna_element in enumerate(dna_elements,0):
            if biosample in biosample2DNAElementPooled2AverageFoldChangeDict:
                if dna_element in biosample2DNAElementPooled2AverageFoldChangeDict[biosample]:
                    average_fold_change = biosample2DNAElementPooled2AverageFoldChangeDict[biosample][dna_element]
                    average_fold_change_array[biosample_index,dna_element_index]=average_fold_change
                else:
                    average_fold_change_array[biosample_index,dna_element_index]=np.nan
            else:
                average_fold_change_array[biosample_index, dna_element_index] = np.nan

    return biosamples, dna_elements, average_fold_change_array
########################################################

########################################################
#For Step2
def fill_average_fold_change_array_row_each_biosample(ENCODEHM2FoldChangeDict):
    encode_hms=[]
    for encode_hm in ENCODEHM2FoldChangeDict:
        if encode_hm not in encode_hms:
            encode_hms.append(encode_hm)

    # sort the hms
    encode_hms = sorted(encode_hms, key=natural_key)

    #Initialize
    average_fold_change_array = np.zeros((1, len(encode_hms)))

    for hm_index,hm in enumerate(encode_hms,0):
        if hm in ENCODEHM2FoldChangeDict:
            average_fold_change = ENCODEHM2FoldChangeDict[hm]
            average_fold_change_array[0,hm_index]=average_fold_change
        else:
            average_fold_change_array[0,hm_index]=np.nan

    return encode_hms, average_fold_change_array
########################################################


###################################################################
def heatmap(data, row_labels, col_labels,ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    #nans are handled here
    data = np.ma.masked_invalid(data)
    ax.patch.set(hatch='x', edgecolor='black')

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom",fontsize=50,labelpad=25)
    cbar.ax.tick_params(labelsize=50)

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels,fontsize=50)
    ax.set_yticklabels(row_labels,fontsize=50)

    # Let the horizontal axes labeling appear on top.
    # ax.tick_params(top=True, bottom=False,labeltop=True, labelbottom=False)
    ax.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True)

    # Rotate the tick labels and set their alignment.
    # plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",rotation_mode="anchor")
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right",rotation_mode="anchor")

    # Turn spines off and create white grid.
    # for edge, spine in ax.spines.items():
        # spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)

    # ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.grid(which="minor", color="black", linestyle='-', linewidth=3)

    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar
###################################################################


###################################################################
def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=["black", "white"],
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A list or array of two color specifications.  The first is used for
        values below a threshold, the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            # kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            kw.update(color='white')
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts
###################################################################



########################################################
#March 30, 2020
#For Step4
def plot_heatmap_rows_signatures(signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict,cancer_type,signatureType,heatmap_output_path,verbose):

    signatures, dna_elements, average_fold_change_array= fill_average_fold_change_array_rows_signatures_columns_dna_elements(signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict)

    ##########################################################################
    # fig, ax = plt.subplots() #does not plot and blocks multiprocessing
    if len(signatures) > 20:
        # if signature in ['SBS1', 'SBS5','SBS40','DBS2','DBS4','ID1','ID2','ID5']:
        fig, ax = plt.subplots(figsize=(60, 100))
        text_fontsize=60
    else:
        fig, ax = plt.subplots(figsize=(50, 40))
        text_fontsize=70

    # im, cbar = heatmap(fold_change_array, signatures, hms, ax=ax,cmap="seismic", cbarlabel="Fold Change [real/simulated]")
    try:
        print('min:%f max:%f' % (np.min(average_fold_change_array), np.max(average_fold_change_array)))
    except ValueError:
        print('average_fold_change_array')
        print(average_fold_change_array)
        print('average_fold_change_array.size')
        print(average_fold_change_array.size)
        print('average_fold_change_array.shape')
        print(average_fold_change_array.shape)

    # Blue White Re
    im, cbar = heatmap(average_fold_change_array, signatures, dna_elements, ax=ax, cmap='seismic',cbarlabel="Fold Change [Real mutations/Simulated Mutations]", vmin=0.25, vmax=1.75)

    #Here we set the fontsize
    texts = annotate_heatmap(im, valfmt="{x:.2f} ", fontsize=text_fontsize)
    print('texts:%s' % texts)

    plt.title(cancer_type, fontsize=90, y=1.01)
    # Results in big squares when array is small
    # plt.tight_layout()

    #############################################################################################################
    #Add signature type
    if len(signatures)==1:
        filename = 'Step4_%s_rows_biosamples_pooled_columns_dna_elements_heatmap.png' %(signatures[0])
    else:
        filename = 'Step5_rows_%s_signatures_columns_dna_elements_heatmap.png' %(signatureType)
    #############################################################################################################

    figureFile = os.path.join(heatmap_output_path, filename)

    fig.savefig(figureFile)

    plt.close()
    ##########################################################################

########################################################



########################################################
#March 30, 2020
#For Step3
def plot_heatmap_rows_biosamples_given_signature(signature2Biosample2DNAElementPooled2AverageFoldChangeDict,heatmap_output_path):

    for signature in signature2Biosample2DNAElementPooled2AverageFoldChangeDict:
        biosamples, dna_elements, average_fold_change_array = fill_average_fold_change_array_rows_biosamples_columns_dna_elements(signature2Biosample2DNAElementPooled2AverageFoldChangeDict[signature])
        print('%s %s %s %s' %(signature,biosamples,dna_elements,average_fold_change_array))

        ##########################################################################
        # fig, ax = plt.subplots() #does not plot and blocks multiprocessing
        if len(biosamples)>20:
            fig, ax = plt.subplots(figsize=(60, 100))
        else:
            fig, ax = plt.subplots(figsize=(50, 40))

        print('min:%f max:%f' %(np.min(average_fold_change_array),np.max(average_fold_change_array)))

        #Blue White Red
        print('average_fold_change_array')
        print(average_fold_change_array)
        print('biosamples:%s' %(biosamples))
        print('dna_elements:%s' % (dna_elements))
        im, cbar = heatmap(average_fold_change_array, biosamples, dna_elements, ax=ax, cmap='seismic', cbarlabel="Fold Change [Real mutations/Simulated Mutations]",vmin=0.25, vmax=1.75)

        texts = annotate_heatmap(im, valfmt="{x:.2f} ", fontsize=70)
        print('texts:%s' %texts)

        plt.title('%s' % (signature), fontsize=90, y=1.01)
        # Results in big squares when array is small
        plt.tight_layout()
        filename = 'Step3_%s_rows_biosamples_columns_dna_elements_heatmap.png' % (signature)
        figureFile = os.path.join(heatmap_output_path,filename)

        fig.savefig(figureFile)

        plt.close()
        ##########################################################################

########################################################


########################################################
def plot_heatmap_row_each_biosample_given_signature(signature2Biosample2DNAElement2FoldChangeDict,cancer_type,heatmap_output_path,verbose):
    if verbose: print('\n\tVerbose Plotting starts for Step2 using signature2Biosample2DNAElement2FoldChangeDict')
    for signature in signature2Biosample2DNAElement2FoldChangeDict:
        for biosample in signature2Biosample2DNAElement2FoldChangeDict[signature]:
            encode_hms, average_fold_change_array = fill_average_fold_change_array_row_each_biosample(signature2Biosample2DNAElement2FoldChangeDict[signature][biosample])

            if verbose: print('\tVerbose %s %s len(encode_hms):%d %s' %(signature, biosample,len(encode_hms),encode_hms))
            cancer_type_encode_biosamples_list=[]
            cancer_type_encode_biosamples_list.append('%s, %s' %(cancer_type,biosample))

            ##########################################################################
            # fig, ax = plt.subplots() #does not plot and blocks multiprocessing
            fig, ax = plt.subplots(figsize=(100, 40))

            # im, cbar = heatmap(fold_change_array, signatures, hms, ax=ax,cmap="seismic", cbarlabel="Fold Change [real/simulated]")
            if verbose: print('\tVerbose min:%f max:%f' %(np.min(average_fold_change_array),np.max(average_fold_change_array)))

            #Blue White Red
            im, cbar = heatmap(average_fold_change_array, cancer_type_encode_biosamples_list, encode_hms, ax=ax, cmap='seismic', cbarlabel="Fold Change [Real mutations/Simulated Mutations]",vmin=0, vmax=2)

            texts = annotate_heatmap(im, valfmt="{x:.2f} ", fontsize=50)
            if verbose: print('\tVerbose texts:%s' %texts)

            plt.title('%s' % (signature), fontsize=50, y=1.01)
            #Results in big squares when array is small
            # plt.tight_layout()
            filename = 'Step2_%s_rows_%s_%s_columns_encode_hms_heatmap.png' % (signature,cancer_type,biosample)
            figureFile = os.path.join(heatmap_output_path,filename)

            # Results in big squares when array is small
            plt.tight_layout()

            fig.savefig(figureFile)
            plt.clf()
            plt.cla()
            plt.close(fig)
            ##########################################################################
    if verbose: print('\tVerbose Plotting ends for Step2 using signature2Biosample2DNAElement2FoldChangeDict')

########################################################


########################################################
def accumulate(signature,
    signature_type,
    signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict,
    accumulated_subsSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict,
    accumulated_dbsSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict,
    accumulated_idSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict):

    if (signature_type==SBS):
        if (signature not in accumulated_subsSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict):
            if signature in signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict:
                accumulated_subsSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict[signature]=signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict[signature]
        else:
            print('There is a situation: %s' %accumulated_subsSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict[signature])

    if (signature_type==DBS):
        if (signature not in accumulated_dbsSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict):
            if signature in signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict:
                accumulated_dbsSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict[signature]=signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict[signature]
        else:
            print('There is a situation: %s' %accumulated_dbsSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict[signature])

    if (signature_type==ID):
        if (signature not in accumulated_idSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict):
            if signature in signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict:
                accumulated_idSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict[signature]=signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict[signature]
        else:
            print('There is a situation: %s' %accumulated_idSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict[signature])
########################################################


#########################################################
def plot_heatmaps(outputDir,jobname,numberofSimulations,epigenomics_files_memos,epigenomics_biosamples,epigenomics_dna_elements,occupancy_type,plusOrMinus_epigenomics,verbose):

    if verbose: print('\tVerbose epigenomics_files_memos:%s' % (epigenomics_files_memos))
    if verbose: print('\tVerbose epigenomics_biosamples:%s' %(epigenomics_biosamples))
    if verbose: print('\tVerbose epigenomics_dna_elements:%s' %(epigenomics_dna_elements))

    accumulated_subsSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict={}
    accumulated_dbsSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict={}
    accumulated_idSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict={}

    #We can plot heatmaps if there are at least one simulation
    if (numberofSimulations>=1):

        cancer_type = jobname
        subsSignature_cutoff_numberofmutations_averageprobability_df=pd.read_csv(os.path.join(outputDir,jobname,DATA,Table_SubsSignature_Cutoff_NumberofMutations_AverageProbability_Filename),sep='\t', header=0, dtype={'cutoff':np.float32,'signature':str, 'number_of_mutations':np.int32,'average_probability':np.float32})
        indelsSignature_cutoff_numberofmutations_averageprobability_df= pd.read_csv(os.path.join(outputDir,jobname,DATA,Table_IndelsSignature_Cutoff_NumberofMutations_AverageProbability_Filename),sep='\t', header=0, dtype={'cutoff':np.float32,'signature':str, 'number_of_mutations':np.int32,'average_probability':np.float32})
        dinucsSignature_cutoff_numberofmutations_averageprobability_df=pd.read_csv(os.path.join(outputDir,jobname,DATA,Table_DinucsSignature_Cutoff_NumberofMutations_AverageProbability_Filename),sep='\t',header=0, dtype={'cutoff':np.float32,'signature':str, 'number_of_mutations':np.int32,'average_probability':np.float32})

        ###########################################
        signatures = []

        #For real run
        for signature in subsSignature_cutoff_numberofmutations_averageprobability_df['signature'].unique():
            signatures.append((signature,SBS))
        for signature in indelsSignature_cutoff_numberofmutations_averageprobability_df['signature'].unique():
            signatures.append((signature,ID))
        for signature in dinucsSignature_cutoff_numberofmutations_averageprobability_df['signature'].unique():
            signatures.append((signature,DBS))

        #For tests
        # signatures.append(('DBS6',DBS))
        # signatures.append(('SBS17b',SBS))
        # signatures.append(('SBS18',SBS))
        # signatures.append(('SBS33',SBS))
        ###########################################

        ########################################################
        def plot_heatmap(list):

            output_dir = list[0]
            signature = list[1]
            signature_type = list[2]
            cancer_type = list[3]
            occupancy_type = list[4]
            signature2Biosample2DNAElement2FoldChangeDict = list[5]
            signature2Biosample2DNAElementPooled2AverageFoldChangeDict = list[6]
            signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict = list[7]
            verbose = list[8]

            os.makedirs(os.path.join(output_dir, cancer_type, FIGURE, ALL, occupancy_type, HEATMAPS), exist_ok=True)
            heatmap_output_path = os.path.join(output_dir, cancer_type, FIGURE, ALL, occupancy_type, HEATMAPS)

            # Step2 Plot a heatmap for Step2 using signature2Biosample2DNAElement2FoldChangeDict
            if signature2Biosample2DNAElement2FoldChangeDict:
                plot_heatmap_row_each_biosample_given_signature(signature2Biosample2DNAElement2FoldChangeDict, cancer_type,heatmap_output_path, verbose)

            # Step3
            if signature2Biosample2DNAElementPooled2AverageFoldChangeDict:
                plot_heatmap_rows_biosamples_given_signature(signature2Biosample2DNAElementPooled2AverageFoldChangeDict,heatmap_output_path)

            # Step4
            if signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict:
                plot_heatmap_rows_signatures(signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict, cancer_type,signature_type,heatmap_output_path, verbose)

            accumulate(signature,
                       signature_type,
                       signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict,
                       accumulated_subsSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict,
                       accumulated_dbsSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict,
                       accumulated_idSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict)
            ##########################################################################
        ########################################################

        #######################################################################
        numofProcesses = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(numofProcesses)
        jobs=[]

        for signature,signature_type in signatures:
            jobs.append(pool.apply_async(fill_signature2CancerType2Biosample2HM2FoldChangeDict, args=(outputDir,numberofSimulations,signature,signature_type,cancer_type,epigenomics_files_memos,epigenomics_biosamples,epigenomics_dna_elements,occupancy_type,plusOrMinus_epigenomics,verbose,), callback=plot_heatmap))
            #Sequential Run
            # list=fill_signature2CancerType2Biosample2HM2FoldChangeDict(outputDir,numberofSimulations,signature,signature_type,cancer_type,epigenomics_files_memos,epigenomics_biosamples,epigenomics_dna_elements,occupancy_type,plusOrMinus_epigenomics,verbose)
            # plot_heatmap(list)

        for job in jobs:
            if verbose: print('\tVerbose %s Worker pid %s job.get():%s ' % (occupancy_type, str(os.getpid()), job.get()))

        heatmap_output_path = os.path.join(outputDir, cancer_type, FIGURE, ALL, occupancy_type, HEATMAPS)

        if accumulated_subsSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict:
            plot_heatmap_rows_signatures(accumulated_subsSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict,'Mutographs ESCC SBS Signatures',SBS,heatmap_output_path, verbose)
        if accumulated_dbsSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict:
            plot_heatmap_rows_signatures(accumulated_dbsSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict, 'Mutographs ESCC DBS Signatures',DBS,heatmap_output_path, verbose)
        if accumulated_idSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict:
            plot_heatmap_rows_signatures(accumulated_idSignature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict, 'Mutographs ESCC ID Signatures',ID,heatmap_output_path, verbose)

        pool.close()
        pool.join()
        ########################################################################

#########################################################



#########################################################
#Original old call
def occupancyAverageSignalFigures(outputDir,jobname,figureAugmentation,numberofSimulations,sample_based,mutationTypes,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus,verbose):
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
    mutationtype_numberofmutations_df = pd.read_csv(os.path.join(outputDir, jobname, DATA,Table_MutationType_NumberofMutations_Filename),sep='\t', header=0,dtype={'mutation_type': str,'number_of_mutations': np.int32})
    subsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(os.path.join(outputDir, jobname, DATA,Table_SubsSignature_Cutoff_NumberofMutations_AverageProbability_Filename),sep='\t', header=0,dtype={'cutoff': np.float32,'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})
    indelsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(os.path.join(outputDir, jobname, DATA,Table_IndelsSignature_Cutoff_NumberofMutations_AverageProbability_Filename),sep='\t', header=0,dtype={'cutoff': np.float32,'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})
    dinucsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(os.path.join(outputDir, jobname, DATA,Table_DinucsSignature_Cutoff_NumberofMutations_AverageProbability_Filename),sep='\t', header=0,dtype={'cutoff': np.float32,'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})

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

    if (SUBS in mutationtype_numberofmutations_df['mutation_type'].unique()):
        numberofSubs=mutationtype_numberofmutations_df.loc[ mutationtype_numberofmutations_df['mutation_type']==SUBS,'number_of_mutations'].values[0]
    if (INDELS in mutationtype_numberofmutations_df['mutation_type'].unique()):
        numberofIndels = mutationtype_numberofmutations_df.loc[ mutationtype_numberofmutations_df['mutation_type']==INDELS,'number_of_mutations'].values[0]
    if (DINUCS in mutationtype_numberofmutations_df['mutation_type'].unique()):
        numberofDinucs = mutationtype_numberofmutations_df.loc[ mutationtype_numberofmutations_df['mutation_type']==DINUCS,'number_of_mutations'].values[0]

    #tissue based
    for mutationType in mutationTypes:
        if verbose: print('\tVerbose Worker pid %s Plot all mutations pooled %s\t%s' %(str(os.getpid()),str(mutationType),libraryFilenameMemo))
        plotAllMutationsPooledWithSimulations('Interval around variant (bp)',ylabel,None,outputDir,jobname,numberofSubs,numberofIndels,numberofDinucs,numberofSimulations,mutationType,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus)
    ##############################################################

    #############################################################################################################################################
    #Plot Signature Based
    #ncomms11383 Fig3b signature based average nucleosome occupancy figures
    if checkValidness(SIGNATUREBASED,outputDir,jobname,occupancy_type):
        if (SBS96 in mutationTypes):
            #Subs Signatures
            if verbose: print('\tVerbose Worker pid %s Plot signature based SBS96 %s' % (str(os.getpid()),libraryFilenameMemo))
            plotSignatureBasedFigures(SBS96,subsSignature_cutoff_numberofmutations_averageprobability_df,sample2SubsSignature2NumberofMutationsDict,outputDir,jobname,isFigureAugmentation,numberofSimulations,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus,verbose)
        if (ID in mutationTypes):
            #Indels Signatures
            if verbose: print('\tVerbose Worker pid %s Plot signature based ID %s' % (str(os.getpid()),libraryFilenameMemo))
            plotSignatureBasedFigures(ID,indelsSignature_cutoff_numberofmutations_averageprobability_df,sample2IndelsSignature2NumberofMutationsDict,outputDir,jobname,isFigureAugmentation,numberofSimulations,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus,verbose)
        if (DBS in mutationTypes):
            # Dinucs Signatures
            if verbose: print('\tVerbose Worker pid %s Plot signature based DBS %s' % (str(os.getpid()),libraryFilenameMemo))
            plotSignatureBasedFigures(DBS,dinucsSignature_cutoff_numberofmutations_averageprobability_df,sample2DinucsSignature2NumberofMutationsDict,outputDir,jobname,isFigureAugmentation,numberofSimulations,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus,verbose)
    #############################################################################################################################################

    ##############################################################
    if sample_based:
        #ALL SAMPLES IN ONE
        #Plot "all samples pooled" and "sample based" signature based in one figure
        plotAllSamplesPooledAndSampleBasedSignaturesFiguresInOneFigure(subsSignature_cutoff_numberofmutations_averageprobability_df,sample2SubsSignature2NumberofMutationsDict,outputDir,jobname,'royalblue','Interval around single point mutation (bp)',ylabel,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus)
        plotAllSamplesPooledAndSampleBasedSignaturesFiguresInOneFigure(indelsSignature_cutoff_numberofmutations_averageprobability_df,sample2IndelsSignature2NumberofMutationsDict,outputDir,jobname,'darkgreen','Interval around indel (bp)',ylabel,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus)
        plotAllSamplesPooledAndSampleBasedSignaturesFiguresInOneFigure(dinucsSignature_cutoff_numberofmutations_averageprobability_df,sample2DinucsSignature2NumberofMutationsDict,outputDir,jobname,'crimson','Interval around indel (bp)',ylabel,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus)

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
            for mutationType in mutationTypes:
                plotAllMutationsPooledWithSimulations('Interval around variant (bp)',ylabel,sample,outputDir,jobname,numberofSubs, numberofIndels,numberofDinucs,numberofSimulations,mutationType,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus)
    ##############################################################
