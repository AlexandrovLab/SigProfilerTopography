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

import matplotlib as mpl
# BACKEND = 'Agg'
# if mpl.get_backend().lower() != BACKEND.lower():
#     # If backend is not set properly a call to describe will hang
#     mpl.use(BACKEND)

from matplotlib import pyplot as plt
from mpl_toolkits import axes_grid1

from statsmodels.stats.weightstats import ztest
import statsmodels.stats.multitest

from decimal import Decimal
import scipy

from SigProfilerTopography.source.commons.TopographyCommons import SBS96
from SigProfilerTopography.source.commons.TopographyCommons import AVERAGE_SIGNAL_ARRAY
from SigProfilerTopography.source.commons.TopographyCommons import ACCUMULATED_COUNT_ARRAY

from SigProfilerTopography.source.commons.TopographyCommons import SBS_CONTEXTS
from SigProfilerTopography.source.commons.TopographyCommons import SBS
from SigProfilerTopography.source.commons.TopographyCommons import ID
from SigProfilerTopography.source.commons.TopographyCommons import DBS

from SigProfilerTopography.source.commons.TopographyCommons import SUBS
from SigProfilerTopography.source.commons.TopographyCommons import INDELS
from SigProfilerTopography.source.commons.TopographyCommons import DINUCS

from SigProfilerTopography.source.commons.TopographyCommons import DATA
from SigProfilerTopography.source.commons.TopographyCommons import FIGURE
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
from SigProfilerTopography.source.commons.TopographyCommons import NUCLEOSOME_DNA_ELEMENT
from SigProfilerTopography.source.commons.TopographyCommons import ATAC_DNA_ELEMENT
from SigProfilerTopography.source.commons.TopographyCommons import OPEN_CHROMATIN

from SigProfilerTopography.source.commons.TopographyCommons import PLOTS
from SigProfilerTopography.source.commons.TopographyCommons import TABLES
from SigProfilerTopography.source.commons.TopographyCommons import DETAILED
from SigProfilerTopography.source.commons.TopographyCommons import EXCEL_FILES

from SigProfilerTopography.source.commons.TopographyCommons import takeAverage
from SigProfilerTopography.source.commons.TopographyCommons import getDictionary

from SigProfilerTopography.source.commons.TopographyCommons import Table_MutationType_NumberofMutations_NumberofSamples_SamplesList_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Table_SubsSignature_Cutoff_NumberofMutations_AverageProbability_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Table_IndelsSignature_Cutoff_NumberofMutations_AverageProbability_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Table_DinucsSignature_Cutoff_NumberofMutations_AverageProbability_Filename

from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofSubsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofIndelsDict
from SigProfilerTopography.source.commons.TopographyCommons import Sample2NumberofDinucsDictFilename

from SigProfilerTopography.source.commons.TopographyCommons import getSample2SubsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2IndelsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import Sample2DinucsSignature2NumberofMutationsDictFilename

from SigProfilerTopography.source.commons.TopographyCommons import natural_key

from SigProfilerTopography.source.commons.TopographyCommons import COLORBAR_SEISMIC
from SigProfilerTopography.source.commons.TopographyCommons import COLORBAR_DISCREET

from SigProfilerTopography.source.commons.TopographyCommons import PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_TOOL
from SigProfilerTopography.source.commons.TopographyCommons import PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_MANUSCRIPT
from SigProfilerTopography.source.commons.TopographyCommons import PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_MANUSCRIPT_OCCUPANCY_ANALYSIS_FIGURE

from SigProfilerTopography.source.commons.TopographyCommons import write_excel_file


plt.rcParams.update({'figure.max_open_warning': 0})

INDIVIDUAL_SIGNATURE='INDIVIDUAL_SIGNATURE'
ALL_SIGNATURES='ALL_SIGNATURES'

ALL_MUTATIONS = [AGGREGATEDSUBSTITUTIONS, AGGREGATEDDINUCS, AGGREGATEDINDELS]

from SigProfilerTopography.source.commons.TopographyCommons import UNDECLARED


#############################################################################
##################### Read Average as Pandas Series #########################
#############################################################################
#Jobname has to be only jobname given in the argument
def readData(sample,signatureName,analyseType,outputDir,jobname,occupancy_type,libraryFilenameMemo,partial_file_name):

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
def readDataForSimulations(sample,signature,analyseType,outputDir,jobname,numberofSimulations,occupancy_type,libraryFilenameMemo,partial_file_name):
    # partial_file_name = 'AverageSignalArray'

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
    else:
        filenameEnd = 'EpigenomicsOccupancy'

    for signature in signature_cutoff_numberofmutations_averageprobability_df['signature'].unique():
        min_list = []
        max_list = []

        label2NumpyArrayDict = {}
        realAverage = readData(None, signature, SIGNATUREBASED, outputDir, jobname,occupancy_type,libraryFilenameMemo,AVERAGE_SIGNAL_ARRAY)
        label2NumpyArrayDict[signature] = realAverage

        for sample in sample2Signature2NumberofMutationsDict:
            if signature in sample2Signature2NumberofMutationsDict[sample]:
                realAverage = readData(sample,signature,SAMPLEBASED_SIGNATUREBASED,outputDir,jobname,occupancy_type,libraryFilenameMemo,AVERAGE_SIGNAL_ARRAY)
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

        figureFile = os.path.join(outputDir, jobname, FIGURE, occupancy_type, filename)

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
def plotSignatureBasedAverageOccupancyFigureWithSimulations(sample,signature,numberofMutations,xlabel,ylabel,label,text,outputDir,jobname,numberofSimulations,color,linestyle,fillcolor,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus,verbose,plot_mode):

    if (occupancy_type==NUCLEOSOMEOCCUPANCY):
        figurenameEnd='_NucleosomeOccupancy.png'
    elif (occupancy_type==EPIGENOMICSOCCUPANCY):
        figurenameEnd='_EpigenomicsOccupancy.png'
    else:
        figurenameEnd='_EpigenomicsOccupancy.png'

    if plot_mode== PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_TOOL:
        real_label='Real %s' %(label)
        simulated_label = 'Simulated %s' % (label)
    elif (plot_mode==PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_MANUSCRIPT) or (plot_mode==PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_MANUSCRIPT_OCCUPANCY_ANALYSIS_FIGURE):
        real_label='Real'
        simulated_label = 'Simulated'

    min_list=[]
    max_list=[]
    min_average_nucleosome_signal=0
    max_average_nucleosome_signal=0

    listofSimulationsSignatureBased = None

    if ((sample is not None) and (signature is not None)):
        realAverage = readData(sample, signature, SAMPLEBASED_SIGNATUREBASED, outputDir, jobname,occupancy_type,libraryFilenameMemo,AVERAGE_SIGNAL_ARRAY)
        if (libraryFilenameMemo is None):
            figurename = '%s_%s_%d' % (signature, sample, numberofMutations)
        else:
            figurename = '%s_%s_%d_%s' % (signature, sample, numberofMutations,libraryFilenameMemo)

        title = '%s_%s' % (signature, sample)
        if (numberofSimulations>0):
            listofSimulationsSignatureBased = readDataForSimulations(sample, signature, SAMPLEBASED_SIGNATUREBASED, outputDir, jobname,numberofSimulations,occupancy_type,libraryFilenameMemo,AVERAGE_SIGNAL_ARRAY)
    else:
        realAverage = readData(None, signature, SIGNATUREBASED, outputDir, jobname,occupancy_type,libraryFilenameMemo,AVERAGE_SIGNAL_ARRAY)
        if (libraryFilenameMemo is None):
            figurename = '%s_%d' % (signature, numberofMutations)
        else:
            figurename = '%s_%d_%s' % (signature, numberofMutations,libraryFilenameMemo)

        title = '%s' % (signature)
        if (numberofSimulations>0):
            listofSimulationsSignatureBased = readDataForSimulations(sample, signature, SIGNATUREBASED, outputDir, jobname,numberofSimulations,occupancy_type,libraryFilenameMemo,AVERAGE_SIGNAL_ARRAY)

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
        fig = plt.figure(figsize=(20,10),facecolor=None,dpi=100)
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

        original = plt.plot(x, realAverage, color=color, label=real_label,linewidth=5)
        listofLegends.append(original[0])

        if (simulationsSignatureBasedMeans is not None):
            simulations = plt.plot(x, simulationsSignatureBasedMeans, color='gray', linestyle=linestyle,  label=simulated_label, linewidth=5)
            listofLegends.append(simulations[0])
            if (simulationsSignatureBasedLows is not None) and (simulationsSignatureBasedHighs is not None):
                plt.fill_between(x, np.array(simulationsSignatureBasedLows), np.array(simulationsSignatureBasedHighs),facecolor=fillcolor)


        if ((simulationsSignatureBasedLows is not None) and (not np.all(np.isnan(simulationsSignatureBasedLows)))):
            min_list.append(np.nanmin(simulationsSignatureBasedLows))
        if ((simulationsSignatureBasedHighs is not None) and (not np.all(np.isnan(simulationsSignatureBasedHighs)))):
            max_list.append(np.nanmax(simulationsSignatureBasedHighs))

        ##############################################################################
        if plot_mode==PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_TOOL:
            # plt.legend(loc='best', handles=listofLegends, prop={'size': 24}, shadow=False, edgecolor='white',facecolor='white', framealpha=0)
            plt.legend(loc='lower left', handles=listofLegends, prop={'size': 24}, shadow=False, edgecolor='white',facecolor='white', framealpha=0)

            #put the number of snps
            tobeWrittenText = "{:,}".format(numberofMutations)
            tobeWrittenText=tobeWrittenText + " " + text
            plt.text(0.99, 0.99, tobeWrittenText, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, fontsize=24)

            #put the library filename
            plt.text(0.01, 0.99, libraryFilename, verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, fontsize=24)

            # This code provides the x and y tick marks and labels
            # plt.xticks(np.arange(-1000, +1001, step=500), fontsize=30)
            plt.xticks(np.arange(-plusOrMinus, plusOrMinus + 1, step=500), fontsize=30)

            # July 27, 2018
            # plt.xlim((-1000, 1000))
            plt.xlim((-plusOrMinus, plusOrMinus))

            # This code puts the tick marks
            plt.tick_params(axis='both', which='major', labelsize=30, width=3, length=10)
            plt.tick_params(axis='both', which='minor', labelsize=30, width=3, length=10)

            plt.title(title, fontsize=40, fontweight='bold')

            plt.xlabel(xlabel, fontsize=32, fontweight='semibold')
            plt.ylabel(ylabel, fontsize=32, fontweight='semibold')

        elif plot_mode==PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_MANUSCRIPT_OCCUPANCY_ANALYSIS_FIGURE:
            plt.legend(loc='upper left', handles=listofLegends, prop={'size': 24}, shadow=False, edgecolor='white',facecolor='white', framealpha=0)

            #put the number of mutations
            tobeWrittenText = "{:,}".format(numberofMutations)
            tobeWrittenText=tobeWrittenText + " " + text
            plt.text(0.99, 0.99, tobeWrittenText, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, fontsize=24)

            #put the library filename
            # plt.text(0.01, 0.99, libraryFilename, verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, fontsize=24)

            # July 27, 2018
            # plt.xlim((-1000, 1000))
            plt.xlim((-plusOrMinus, plusOrMinus))

            # put/set axis ticks
            plt.tick_params(axis='both', which='major', labelsize=30, width=3, length=10)
            plt.tick_params(axis='both', which='minor', labelsize=30, width=3, length=10)

            # put/set axis labels
            plt.setp(ax.get_xticklabels(), visible=True)
            plt.setp(ax.get_yticklabels(), visible=True)

            plt.title(title, fontsize=40, fontweight='bold')
            plt.xlabel(xlabel, fontsize=32, fontweight='semibold')
            plt.ylabel(ylabel, fontsize=32, fontweight='semibold')

        elif plot_mode==PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_MANUSCRIPT:
            # legend=ax.legend(loc='upper left', handles=listofLegends, prop={'size': 35}, shadow=False, edgecolor='black',facecolor='white', framealpha=0)
            legend = ax.legend(handles=listofLegends, prop={'size': 35}, loc='upper left')
            if (legend is not None):
                frame = legend.get_frame()
                frame.set_facecolor('white')
                frame.set_edgecolor('black')

            # set axis ticks
            # ax.tick_params(axis='both', which='both', length=0)
            ax.tick_params(axis='x', which='both', length=0)
            ax.tick_params(axis='y', which='both', length=0)
            # set axis labels
            plt.setp(ax.get_xticklabels(), visible=False)
            plt.setp(ax.get_yticklabels(), visible=False)
        ##############################################################################

        #Put vertical line at x=0
        # plt.axvline(x=0, ymin=0, ymax=1, color='gray', linestyle='--')
        plt.axvline(x=0, color='gray', linestyle='--')

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

        filename = figurename.replace(' ', '') + figurenameEnd

        #######################################################################
        # new code
        if (sample is None):
            if occupancy_type==NUCLEOSOMEOCCUPANCY:
                figureFile = os.path.join(outputDir, jobname, FIGURE, occupancy_type, filename)
            elif occupancy_type==EPIGENOMICSOCCUPANCY:
                figureFile = os.path.join(outputDir, jobname, FIGURE, occupancy_type, PLOTS, filename)
        else:
            os.makedirs(os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, occupancy_type), exist_ok=True)
            figureFile = os.path.join(outputDir, jobname, FIGURE, SAMPLES, sample, occupancy_type, filename)
        #######################################################################

        fig.savefig(figureFile, dpi=100, bbox_inches="tight")
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
def plotAllMutationsPooledWithSimulations(xlabel,ylabel,sample,outputDir,jobname,numberofSubs,numberofIndels,numberofDinucs,numberofSimulations,mutationType,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus,plot_mode):
    if mutationType==SBS96:
        to_be_added_to_the_filename='SBS%s' %(mutationType)
    else:
        to_be_added_to_the_filename=mutationType

    if (occupancy_type==NUCLEOSOMEOCCUPANCY):
        filenameEnd='NucleosomeOccupancy'
    elif (occupancy_type==EPIGENOMICSOCCUPANCY):
        filenameEnd='EpigenomicsOccupancy'
    else:
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
            realAggregatedSubstitutions = readData(None,None,AGGREGATEDSUBSTITUTIONS,outputDir,jobname,occupancy_type,libraryFilenameMemo,AVERAGE_SIGNAL_ARRAY)
        if (ID==mutationType):
            realAggregatedIndels = readData(None,None, AGGREGATEDINDELS, outputDir,jobname,occupancy_type,libraryFilenameMemo,AVERAGE_SIGNAL_ARRAY)
        if (DBS==mutationType):
            realAggregatedDinucs = readData(None,None,AGGREGATEDDINUCS,outputDir,jobname,occupancy_type,libraryFilenameMemo,AVERAGE_SIGNAL_ARRAY)

        if (numberofSimulations>0):
            if (SBS96==mutationType):
                listofSimulationsAggregatedSubstitutions = readDataForSimulations(None,None,AGGREGATEDSUBSTITUTIONS,outputDir,jobname,numberofSimulations,occupancy_type,libraryFilenameMemo,AVERAGE_SIGNAL_ARRAY)
            if (ID==mutationType):
                listofSimulationsAggregatedIndels = readDataForSimulations(None,None,AGGREGATEDINDELS,outputDir,jobname,numberofSimulations,occupancy_type,libraryFilenameMemo,AVERAGE_SIGNAL_ARRAY)
            if (DBS==mutationType):
                listofSimulationsAggregatedDinucs = readDataForSimulations(None,None,AGGREGATEDDINUCS,outputDir,jobname,numberofSimulations,occupancy_type,libraryFilenameMemo,AVERAGE_SIGNAL_ARRAY)
    #######################################################################################################################


    #######################################################################################################################
    else:
        # filename = '%s_%s_Aggregated_Substitutions_%d_Indels_%d.png' % (sample, jobname, numberofSPMs, numberofIndels)
        if (libraryFilenameMemo is None):
            filename = '%s_Aggregated_All_Mutations_%s_%s.png' % (sample,to_be_added_to_the_filename,filenameEnd)
        else:
            filename = '%s_Aggregated_All_Mutations_%s_%s_%s.png' % (sample,to_be_added_to_the_filename,libraryFilenameMemo,filenameEnd)

        if (SBS96 == mutationType):
            realAggregatedSubstitutions = readData(sample,None,SAMPLEBASED_AGGREGATEDSUBSTITUTIONS,outputDir,jobname,occupancy_type,libraryFilenameMemo,AVERAGE_SIGNAL_ARRAY)
        if (ID == mutationType):
            realAggregatedIndels = readData(sample,None,SAMPLEBASED_AGGREGATEDINDELS,outputDir,jobname,occupancy_type,libraryFilenameMemo,AVERAGE_SIGNAL_ARRAY)
        if (DBS == mutationType):
            realAggregatedDinucs = readData(sample,None,SAMPLEBASED_AGGREGATEDDINUCS,outputDir,jobname,occupancy_type,libraryFilenameMemo,AVERAGE_SIGNAL_ARRAY)

        if (numberofSimulations>0):
            if (SBS96 == mutationType):
                listofSimulationsAggregatedSubstitutions = readDataForSimulations(sample, None, SAMPLEBASED_AGGREGATEDSUBSTITUTIONS, outputDir,jobname,numberofSimulations,occupancy_type,libraryFilenameMemo,AVERAGE_SIGNAL_ARRAY)
            if (ID == mutationType):
                listofSimulationsAggregatedIndels = readDataForSimulations(sample, None,SAMPLEBASED_AGGREGATEDINDELS, outputDir,jobname,numberofSimulations,occupancy_type,libraryFilenameMemo,AVERAGE_SIGNAL_ARRAY)
            if (DBS == mutationType):
                listofSimulationsAggregatedDinucs = readDataForSimulations(sample, None, SAMPLEBASED_AGGREGATEDDINUCS, outputDir,jobname,numberofSimulations,occupancy_type,libraryFilenameMemo,AVERAGE_SIGNAL_ARRAY)
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
    fig = plt.figure(figsize=(20, 10), facecolor=None, dpi=100)
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

    if (plot_mode == PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_TOOL):
        if (SBS96 == mutationType):
            real_label = 'Real Aggregated SBSs'
            simulated_label = 'Simulated Aggregated SBSs'
        elif (DBS == mutationType):
            real_label = 'Real Aggregated DBSs'
            simulated_label = 'Simulated Aggregated DBSs'
        elif (ID == mutationType):
            real_label = 'Real Aggregated IDs'
            simulated_label = 'Simulated Aggregated IDs'
    elif ((plot_mode == PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_MANUSCRIPT) or (plot_mode==PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_MANUSCRIPT_OCCUPANCY_ANALYSIS_FIGURE)):
        real_label = 'Real'
        simulated_label = 'Simulated'

    ##############################################################
    if (realAggregatedSubstitutions is not None):
        aggSubs = plt.plot(x, realAggregatedSubstitutions, 'royalblue', label=real_label,linewidth=3, zorder=10)
        listofLegends.append(aggSubs[0])
    if (simulationsAggregatedSubstitutionsMedians is not None):
        simsAggSubs = plt.plot(x, simulationsAggregatedSubstitutionsMedians, color='gray',linestyle='dashed', label=simulated_label, linewidth=3,zorder=10)
        listofLegends.append(simsAggSubs[0])
        if (simulationsAggregatedSubstitutionsLows is not None) and (simulationsAggregatedSubstitutionsHighs is not None):
            plt.fill_between(x,np.array(simulationsAggregatedSubstitutionsLows),np.array(simulationsAggregatedSubstitutionsHighs),facecolor='lightblue',zorder=10)
    ##############################################################

    ##############################################################
    if (realAggregatedIndels is not None):
        aggIndels = plt.plot(x, realAggregatedIndels, 'darkgreen', label=real_label,linewidth=3,zorder=10)
        listofLegends.append(aggIndels[0])
    if simulationsAggregatedIndelsMedians is not None:
        simsAggIndels = plt.plot(x, simulationsAggregatedIndelsMedians, color='gray', linestyle='dashed',label=simulated_label, linewidth=3,zorder=10)
        listofLegends.append(simsAggIndels[0])
        if (simulationsAggregatedIndelsLows is not None) and (simulationsAggregatedIndelsHighs is not None):
            plt.fill_between(x,np.array(simulationsAggregatedIndelsLows),np.array(simulationsAggregatedIndelsHighs),facecolor='lightgreen',zorder=5)
    ##############################################################

    ##############################################################
    if (realAggregatedDinucs is not None):
        aggDinucs = plt.plot(x, realAggregatedDinucs, 'crimson', label=real_label,linewidth=3,zorder=10)
        listofLegends.append(aggDinucs[0])
    if simulationsAggregatedDinucsMedians is not None:
        simsAggDinucs = plt.plot(x, simulationsAggregatedDinucsMedians, color='gray', linestyle='dashed',label=simulated_label, linewidth=3,zorder=10)
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

    plt.xlabel(xlabel, fontsize=32, fontweight='semibold')
    plt.ylabel(ylabel, fontsize=32, fontweight='semibold')

    ######################################################################################
    if (sample is None):
        if occupancy_type==EPIGENOMICSOCCUPANCY:
            figureFile = os.path.join(outputDir,jobname,FIGURE,occupancy_type,PLOTS,filename)
        elif occupancy_type==NUCLEOSOMEOCCUPANCY:
            figureFile = os.path.join(outputDir,jobname,FIGURE,occupancy_type,filename)
    else:
        os.makedirs(os.path.join(outputDir, jobname,FIGURE,SAMPLES,sample,occupancy_type), exist_ok=True)
        figureFile = os.path.join(outputDir,jobname,FIGURE,SAMPLES,sample,occupancy_type,filename)
    ######################################################################################

    fig.savefig(figureFile,dpi=100, bbox_inches="tight")
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
        # if (len(filenames)>0):
        #     return True
        # else:
        #     return False
    else:
        return False
#########################################################


#########################################################
def plotSignatureBasedFigures(mutationType,signature_cutoff_numberofmutations_averageprobability_df,sample2Signature2NumberofMutationsDict,outputDir,jobname,numberofSimulations,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus,verbose,plot_mode):
    if (occupancy_type==NUCLEOSOMEOCCUPANCY):
        ylabel = 'Average nucleosome signal'
    elif (occupancy_type==EPIGENOMICSOCCUPANCY):
        ylabel = 'Average epigenomics signal'
    else:
        #For epigenomics, epigenomics_dir_name can be different from EPIGENOMICSOCCUPANCY
        ylabel = 'Average epigenomics signal'

    if (mutationType==SBS96):
        xlabel = 'Interval around single base substitution (bp)'
        # label = 'Single Base Substitutions'
        label = 'SBSs'
        text = 'subs'
        color = 'royalblue'
        fillcolor = 'lightblue'
        linestyle='dashed'
    elif (mutationType==ID):
        xlabel = 'Interval around insertion and deletion (bp)'
        # label = 'Insertions and Deletions'
        label = 'IDs'
        text = 'indels'
        color = 'darkgreen'
        fillcolor = 'lightgreen'
        linestyle='dashed'
    elif (mutationType==DBS):
        xlabel = 'Interval around doublet base substitution (bp)'
        # label = 'Doublet Base Substitutions'
        label = 'DBSs'
        text = 'dinucs'
        color = 'crimson'
        fillcolor = 'lightpink'
        linestyle='dashed'

    for signature in signature_cutoff_numberofmutations_averageprobability_df['signature'].unique():
        #[signature cutoff numberofMutations  averageProbability]
        signatureBasedNumberofMutations = int(signature_cutoff_numberofmutations_averageprobability_df[signature_cutoff_numberofmutations_averageprobability_df['signature']==signature]['number_of_mutations'].values[0])
        plotSignatureBasedAverageOccupancyFigureWithSimulations(None, signature,
                                                                          signatureBasedNumberofMutations,
                                                                          xlabel,ylabel,label,text,
                                                                          outputDir, jobname,
                                                                          numberofSimulations, color,linestyle,fillcolor,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus,verbose,plot_mode)

    # SampleBased Subs SignatureBased Nucleosome Occupancy Figures
    for sample in sample2Signature2NumberofMutationsDict:
        for signature in sample2Signature2NumberofMutationsDict[sample]:
            sampleBasedSignatureBasedNumberofMutations = sample2Signature2NumberofMutationsDict[sample][signature]
            plotSignatureBasedAverageOccupancyFigureWithSimulations(sample, signature,
                                                                              sampleBasedSignatureBasedNumberofMutations,
                                                                              xlabel,ylabel,label,text,
                                                                              outputDir, jobname,
                                                                              numberofSimulations,color,linestyle,fillcolor,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus,verbose,plot_mode)
#########################################################


# ########################################################
# def calculate_fold_change(output_dir,numberofSimulations,signature,cancer_type,dna_element,occupancy_type,plusOrMinus_epigenomics,verbose):
#     center=plusOrMinus_epigenomics
#     plusorMinus=250
#     start=center-plusorMinus
#     end= center+plusorMinus+1
#
#     if verbose: print('\n\tVerbose ----------> %s %s %s' % (signature, cancer_type, dna_element))
#     avg_real_signal = None
#
#     # SBS1_sim1_ENCFF330CCJ_osteoblast_H3K79me2-human_AverageSignalArray.txt
#     avg_real_data_signal_array = readData(None, signature, SIGNATUREBASED, output_dir, cancer_type, occupancy_type,dna_element,AVERAGE_SIGNAL_ARRAY)
#
#     if avg_real_data_signal_array is not None:
#         #If there is nan in the list np.mean returns nan.
#         # 1st way
#         # avg_real_data_signal_array[np.isnan(avg_real_data_signal_array)] = 0
#         # avg_real_signal = np.mean(avg_real_data_signal_array[1750:2251])
#         # 2nd way
#         if not np.all(np.isnan(avg_real_data_signal_array[start:end])):
#             avg_real_signal = np.nanmean(avg_real_data_signal_array[start:end])
#
#     avg_simulated_signal = None
#     if (numberofSimulations > 0):
#         listofSimulationsSignatureBased = readDataForSimulations(None, signature, SIGNATUREBASED, output_dir,cancer_type, numberofSimulations, occupancy_type,dna_element,AVERAGE_SIGNAL_ARRAY)
#
#         if ((listofSimulationsSignatureBased is not None) and listofSimulationsSignatureBased):
#             stackedSimulationsSignatureBased = np.vstack(listofSimulationsSignatureBased)
#             (rows, cols) = stackedSimulationsSignatureBased.shape
#             if verbose: print('\tVerbose After np.vstack --- stackedSimulationsSignatureBased rows:%d cols:%d' % (rows, cols))
#
#             #Get the part that we are interested
#             stackedSimulationsSignatureBased_of_interest= stackedSimulationsSignatureBased[:,start:end]
#             (rows, cols) = stackedSimulationsSignatureBased_of_interest.shape
#             if verbose: print(
#                 '\tVerbose After get part interested in --- stackedSimulationsSignatureBased_of_interest rows:%d cols:%d' % (rows, cols))
#
#             # Get rid of all the rows with all nans
#             stackedSimulationsSignatureBased_of_interest = stackedSimulationsSignatureBased_of_interest[~np.isnan(stackedSimulationsSignatureBased_of_interest).all(axis=1)]
#             (rows, cols) = stackedSimulationsSignatureBased_of_interest.shape
#             if verbose: print(
#                 '\tVerbose After removing all nan rows --- stackedSimulationsSignatureBased_of_interest rows:%d cols:%d' % (rows, cols))
#
#             simulationsSignatureBasedMeans = np.nanmean(stackedSimulationsSignatureBased_of_interest, axis=1)
#             if not np.all(np.isnan(simulationsSignatureBasedMeans)):
#                 avg_simulated_signal = np.nanmean(simulationsSignatureBasedMeans)
#
#     if (avg_real_signal is not None) and (avg_simulated_signal is not None):
#         fold_change = avg_real_signal / avg_simulated_signal
#         if verbose: print('\tVerbose ----------> %s %s %s avg_real_signal:%f\tavg_simulated_signal:%f\tfold change:%f' % (signature, cancer_type, dna_element,avg_real_signal, avg_simulated_signal, fold_change))
#         return fold_change
#     else:
#         if verbose: print('\tVerbose ----------> %s %s %s avg_real_signal:%s\tavg_simulated_signal:%s\tfold change: nan' % (signature,cancer_type,dna_element,str(avg_real_signal),str(avg_simulated_signal)))
#         return None
# ########################################################



########################################################
#Plot heatmaps
# Used for after Step2
# For Step4 and Step5
def fill_average_fold_change_array_rows_signatures_columns_dna_elements(signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict,epigenomics_dna_elements):
    signatures=[]
    dna_elements=[]

    for signature in signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict:
        if signature not in signatures:
            signatures.append(signature)
        for dna_element in signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict[signature]:
            if dna_element not in dna_elements:
                dna_elements.append(dna_element)

    # Enlarge dna_elements with epigenomics_dna_elements
    for epigenomics_dna_element in epigenomics_dna_elements:
        if epigenomics_dna_element not in dna_elements:
            dna_elements.append(epigenomics_dna_element)

    # Sort the dna_elements and signatures
    dna_elements=sorted(dna_elements,key=natural_key)
    signatures=sorted(signatures,key=natural_key)

    # Initialize
    average_fold_change_array = np.zeros((len(signatures), len(dna_elements)))

    # Fill the average_fold_change_array
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
            # before
            # average_fold_change = ENCODEHM2FoldChangeDict[hm]
            #updated for combined p value method
            average_fold_change = ENCODEHM2FoldChangeDict[hm][6]

            average_fold_change_array[0,hm_index]=average_fold_change
        else:
            average_fold_change_array[0,hm_index]=np.nan

    return encode_hms, average_fold_change_array
########################################################


def add_colorbar(im, aspect=20, pad_fraction=0.5, **kwargs):
    """Add a vertical color bar to an image plot."""
    divider = axes_grid1.make_axes_locatable(im.axes)
    width = axes_grid1.axes_size.AxesY(im.axes, aspect=1./aspect)
    pad = axes_grid1.axes_size.Fraction(pad_fraction, width)
    current_ax = plt.gca()
    cax = divider.append_axes("right", size=width, pad=pad)
    plt.sca(current_ax)
    return im.axes.figure.colorbar(im, cax=cax, **kwargs)

###################################################################
def heatmap(data, row_labels, col_labels,ax=None, fontsize=None,
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

    # # Create colorbar
    # cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    # cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom",fontsize=fontsize,labelpad=25)
    # cbar.ax.tick_params(labelsize=fontsize)

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.

    ax.set_xticklabels(col_labels,fontsize=fontsize)
    ax.set_yticklabels(row_labels,fontsize=fontsize)

    # Let the x axes labeling appear on bottom.
    ax.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True)
    ax.tick_params(which="minor", bottom=False, left=False)

    # Rotate the tick labels and set their alignment.
    # plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",rotation_mode="anchor")
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right",rotation_mode="anchor")

    # Turn spines off and create white grid.
    # for edge, spine in ax.spines.items():
    # spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)

    ax.grid(which="minor", color="black", linestyle='-', linewidth=3)
    ax.grid(b=False, which="major")

    # return im, cbar
    return im

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
        valfmt = mpl.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            # kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            kw.update(color='black')
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts
###################################################################

###################################################################
#November 10, 2020
def plot_heatmap_rows_biosamples_columns_pooled_DNA_elements(step2_signature2biosample2dna_element2avg_fold_change_dict,
                                                             epigenomics_dna_elements,
                                                             cancer_type,
                                                             heatmaps_output_dir,
                                                             verbose):


    for signature in step2_signature2biosample2dna_element2avg_fold_change_dict:

        biosamples, dna_elements, average_fold_change_array= fill_average_fold_change_array_rows_signatures_columns_dna_elements(step2_signature2biosample2dna_element2avg_fold_change_dict[signature],epigenomics_dna_elements)

        #Update ATAC-Seq to Chromatin
        dna_elements = ['Chromatin' if ATAC_DNA_ELEMENT in dna_element else dna_element for dna_element in dna_elements]

        if verbose:
            print('biosamples: %s' %(biosamples))
            print('dna_elements: %s' %(dna_elements))
            print('average_fold_change_array: %s' %(average_fold_change_array))

        ##########################################################################
        fig, ax = plt.subplots(figsize=(len(dna_elements),len(biosamples)))

        try:
            if verbose: print('\tVerbose min:%f max:%f' % (np.min(average_fold_change_array), np.max(average_fold_change_array)))
        except ValueError:
            if verbose: print('\tVerbose average_fold_change_array:%s' %(average_fold_change_array))
            if verbose: print('\tVerbose average_fold_change_array.size')
            if verbose: print(average_fold_change_array.size)
            if verbose: print('\tVerbose average_fold_change_array.shape')
            if verbose: print(average_fold_change_array.shape)

        # Blue White Red
        im = heatmap(average_fold_change_array, biosamples, dna_elements, ax=ax, cmap='seismic',cbarlabel="Fold Change [Real mutations/Simulated Mutations]", vmin=0.25, vmax=1.75)
        texts = annotate_heatmap(im, valfmt="{x:.2f} ")
        title="%s %s" %(cancer_type,signature)
        plt.title(title, y=1.01)

        # Results in big squares when array is small
        plt.tight_layout()

        #############################################################################################################
        filename = 'Step2_%s_rows_biosamples_columns_dna_elements_heatmap.png' %(signature)
        figureFile = os.path.join(heatmaps_output_dir, filename)
        fig.savefig(figureFile,bbox_inches='tight')
        plt.close()
        ##########################################################################

###################################################################


########################################################
#March 30, 2020
#After Filtering, After Filtering and Significant Ones
def plot_heatmap_rows_signatures_columns_pooled_DNA_elements(signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict,
                                                             signature2dna_element2significancedict,
                                                             epigenomics_dna_elements,
                                                             cancer_type,
                                                             signatureType,
                                                             heatmaps_output_dir,
                                                             filename_text,
                                                             verbose):

    signatures, dna_elements, average_fold_change_array = fill_average_fold_change_array_rows_signatures_columns_dna_elements(signature2BiosamplePooledDNAElementPooled2AverageFoldChangeDict,epigenomics_dna_elements)

    #Update ATAC-Seq to Chromatin
    dna_elements = [OPEN_CHROMATIN if ATAC_DNA_ELEMENT in dna_element else dna_element for dna_element in dna_elements]

    if signatureType==SBS:
        signatures = ['ALL SUBSTITUTIONS' if AGGREGATEDSUBSTITUTIONS in signature else signature for signature in signatures]
    elif signatureType==DBS:
        signatures = ['ALL DINUCLEOTIDES' if AGGREGATEDDINUCS in signature else signature for signature in signatures]
    elif signatureType==ID:
        signatures = ['ALL INDELS' if AGGREGATEDINDELS in signature else signature for signature in signatures]

    if verbose:
        print('signatures: %s' %(signatures))
        print('dna_elements: %s' %(dna_elements))
        print('average_fold_change_array: %s' %(average_fold_change_array))

    ##########################################################################
    fig, ax = plt.subplots(figsize=(len(dna_elements),len(signatures)))

    try:
        if verbose: print('\tVerbose min:%f max:%f' % (np.min(average_fold_change_array), np.max(average_fold_change_array)))
    except ValueError:
        if verbose: print('\tVerbose average_fold_change_array:%s' %(average_fold_change_array))
        if verbose: print('\tVerbose average_fold_change_array.size')
        if verbose: print(average_fold_change_array.size)
        if verbose: print('\tVerbose average_fold_change_array.shape')
        if verbose: print(average_fold_change_array.shape)

    # If dna_element is not statistically significant makes ir cell color white starts
    # backup_average_fold_change_array=average_fold_change_array.copy()
    #
    # for signature_index, signature in enumerate(signatures,0):
    #     for dna_element_index, dna_element in enumerate(dna_elements,0):
    #         if signature2dna_element2significancedict is not  None:
    #             if signature=='ALL SUBSTITUTIONS':
    #                 signature=AGGREGATEDSUBSTITUTIONS
    #             elif signature=='ALL DINUCLEOTIDES':
    #                 signature=AGGREGATEDDINUCS
    #             elif signature=='ALL INDELS':
    #                 signature=AGGREGATEDINDELS
    #             if signature in signature2dna_element2significancedict:
    #                 if dna_element not in signature2dna_element2significancedict[signature]:
    #                     backup_average_fold_change_array[signature_index, dna_element_index]=1.0
    #
    # im, cbar = heatmap(backup_average_fold_change_array, signatures, dna_elements, ax=ax, cmap='seismic',cbarlabel="Fold Change [Real mutations/Simulated Mutations]", vmin=0.25, vmax=1.75)
    # If dna_element is not statistically significant makes ir cell color white ends

    # Color the heatmap
    im = heatmap(average_fold_change_array,
                       signatures,
                       dna_elements,
                       ax=ax,
                       cmap='seismic',
                       cbarlabel="Fold Change [Real mutations/Simulated Mutations]",
                       vmin=0.25,
                       vmax=1.75)

    cbar = add_colorbar(im)
    cbar.ax.set_ylabel('Fold Change [Real mutations/Simulated Mutations]', rotation=-90, va="bottom", fontsize=None, labelpad=25)
    cbar.ax.tick_params(labelsize=None)

    # Write average fold change w/wo star in each heatmap cell
    for signature_index, signature in enumerate(signatures,0):
        for dna_element_index, dna_element in enumerate(dna_elements,0):
            if signature2dna_element2significancedict is not None:
                if signature=='ALL SUBSTITUTIONS':
                    signature=AGGREGATEDSUBSTITUTIONS
                elif signature=='ALL DINUCLEOTIDES':
                    signature=AGGREGATEDDINUCS
                elif signature=='ALL INDELS':
                    signature=AGGREGATEDINDELS
                if signature in signature2dna_element2significancedict:
                    if dna_element in signature2dna_element2significancedict[signature]:
                        text = ax.text(dna_element_index, signature_index, "%.2f*" %(average_fold_change_array[signature_index, dna_element_index]), ha="center", va="center", color="k")
                    elif (dna_element==OPEN_CHROMATIN) and any([True for key in signature2dna_element2significancedict[signature].keys() if ATAC_DNA_ELEMENT in key]):
                        text = ax.text(dna_element_index, signature_index, "%.2f*" %(average_fold_change_array[signature_index, dna_element_index]), ha="center", va="center", color="k")
                    else:
                        text = ax.text(dna_element_index, signature_index, "%.2f" %(average_fold_change_array[signature_index, dna_element_index]), ha="center", va="center", color="k")
                else:
                    text = ax.text(dna_element_index, signature_index,"%.2f" % (average_fold_change_array[signature_index, dna_element_index]), ha="center",va="center", color="k")
            else:
                text = ax.text(dna_element_index, signature_index,"%.2f" % (average_fold_change_array[signature_index, dna_element_index]), ha="center",va="center", color="k")


    # texts = annotate_heatmap(im, valfmt="{x:.2f} ")
    plt.title(cancer_type, y=1.01)

    # Results in big squares when array is small
    plt.tight_layout()

    #############################################################################################################
    #Add signature type
    if filename_text:
       filename = '%s_rows_%s_signatures_columns_dna_elements_heatmap.png' %(filename_text,signatureType)
    else:
        filename = 'rows_%s_signatures_columns_dna_elements_heatmap.png' % (signatureType)
    #############################################################################################################

    figureFile = os.path.join(heatmaps_output_dir, filename)
    fig.savefig(figureFile,bbox_inches='tight')
    plt.close()
    ##########################################################################

########################################################



########################################################
# Detailed heatmap before filtering: row biosample, column DNA elements (unpooled and uncombined)
# Under heatmaps/detailed
def plot_heatmap_one_row_only_for_each_biosample_given_signature(signature2Biosample2DNAElement2FoldChangeDict,cancer_type,heatmap_output_path,verbose):
    if verbose: print('\n\tVerbose Plotting starts for Step2 using signature2Biosample2DNAElement2FoldChangeDict')
    for signature in signature2Biosample2DNAElement2FoldChangeDict:
        for biosample in signature2Biosample2DNAElement2FoldChangeDict[signature]:
            encode_hms, average_fold_change_array = fill_average_fold_change_array_row_each_biosample(signature2Biosample2DNAElement2FoldChangeDict[signature][biosample])

            if verbose: print('\tVerbose %s %s len(encode_hms):%d %s' %(signature, biosample,len(encode_hms),encode_hms))
            cancer_type_encode_biosamples_list=[]
            cancer_type_encode_biosamples_list.append('%s, %s' %(cancer_type,biosample))

            ##########################################################################
            fig, ax = plt.subplots(figsize=(len(encode_hms), len(cancer_type_encode_biosamples_list)))

            if verbose: print('\tVerbose min:%f max:%f' %(np.min(average_fold_change_array),np.max(average_fold_change_array)))

            #Blue White Red
            im = heatmap(average_fold_change_array, cancer_type_encode_biosamples_list, encode_hms, ax=ax, cmap='seismic', cbarlabel="Fold Change [Real mutations/Simulated Mutations]",vmin=0, vmax=2)
            texts = annotate_heatmap(im, valfmt="{x:.2f} ")
            if verbose: print('\tVerbose texts:%s' %texts)

            plt.title('%s' % (signature), y=1.01)
            filename = 'Step1_%s_%s_rows_%s_columns_dna_elements_heatmap.png' % (cancer_type,signature,biosample)
            figureFile = os.path.join(heatmap_output_path,filename)

            plt.tight_layout()
            fig.savefig(figureFile,bbox_inches='tight')
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



########################################################
#Sep24, 2020
def compute_fold_change_with_p_values_plot_heatmaps(combine_p_values_method,
                                          fold_change_window_size,
                                          num_of_real_data_avg_overlap,
                                          outputDir,
                                          jobname,
                                          numberofSimulations,
                                          mutation_types_contexts,
                                          nucleosome_file,
                                          nucleosome_biosample,
                                          epigenomics_files_memos,
                                          epigenomics_biosamples,
                                          epigenomics_dna_elements,
                                          plusOrMinus_epigenomics,
                                          plusOrMinus_nucleosome,
                                          epigenomics_heatmap_significance_level,
                                          verbose):

    if verbose: print('\tVerbose epigenomics_files_memos:%s' % (epigenomics_files_memos))
    if verbose: print('\tVerbose epigenomics_biosamples:%s' %(epigenomics_biosamples))
    if verbose: print('\tVerbose epigenomics_dna_elements:%s' %(epigenomics_dna_elements))

    #Initialize
    sbs_signatures=[]
    dbs_signatures=[]
    id_signatures=[]

    #We can plot heatmaps if there is at least one simulation
    if (numberofSimulations>=1):
        for mutation_type_context in mutation_types_contexts:
            if (mutation_type_context in SBS_CONTEXTS):
                subsSignature_cutoff_numberofmutations_averageprobability_df=pd.read_csv(os.path.join(outputDir,jobname,DATA,Table_SubsSignature_Cutoff_NumberofMutations_AverageProbability_Filename),sep='\t', header=0, dtype={'cutoff':np.float32,'signature':str, 'number_of_mutations':np.int32,'average_probability':np.float32})
                sbs_signatures = subsSignature_cutoff_numberofmutations_averageprobability_df['signature'].unique().tolist()

        if (DBS in mutation_types_contexts):
            dinucsSignature_cutoff_numberofmutations_averageprobability_df=pd.read_csv(os.path.join(outputDir,jobname,DATA,Table_DinucsSignature_Cutoff_NumberofMutations_AverageProbability_Filename),sep='\t',header=0, dtype={'cutoff':np.float32,'signature':str, 'number_of_mutations':np.int32,'average_probability':np.float32})
            dbs_signatures = dinucsSignature_cutoff_numberofmutations_averageprobability_df['signature'].unique().tolist()

        if (ID in mutation_types_contexts):
            indelsSignature_cutoff_numberofmutations_averageprobability_df= pd.read_csv(os.path.join(outputDir,jobname,DATA,Table_IndelsSignature_Cutoff_NumberofMutations_AverageProbability_Filename),sep='\t', header=0, dtype={'cutoff':np.float32,'signature':str, 'number_of_mutations':np.int32,'average_probability':np.float32})
            id_signatures=indelsSignature_cutoff_numberofmutations_averageprobability_df['signature'].unique().tolist()


        sbs_signatures.append(AGGREGATEDSUBSTITUTIONS)
        dbs_signatures.append(AGGREGATEDDINUCS)
        id_signatures.append(AGGREGATEDINDELS)

        #For tests
        # signatures.append(('DBS6',DBS))
        # signatures.append(('SBS17b',SBS))
        # signatures.append(('SBS18',SBS))
        # signatures.append(('SBS33',SBS))

        # Combined p value and plot heatmaps
        compute_fold_change_with_combined_p_values_plot_heatmaps(combine_p_values_method,
                                                    fold_change_window_size,
                                                    num_of_real_data_avg_overlap,
                                                    plusOrMinus_epigenomics,
                                                    plusOrMinus_nucleosome,
                                                    epigenomics_heatmap_significance_level,
                                                    outputDir,
                                                    jobname,
                                                    numberofSimulations,
                                                    nucleosome_file,
                                                    nucleosome_biosample,
                                                    epigenomics_files_memos,
                                                    epigenomics_biosamples,
                                                    epigenomics_dna_elements,
                                                    sbs_signatures,
                                                    dbs_signatures,
                                                    id_signatures,
                                                    verbose)


########################################################

########################################################
#Sep24, 2020
def compute_fold_change_with_combined_p_values_plot_heatmaps(combine_p_values_method,
                                                    fold_change_window_size,
                                                    num_of_real_data_avg_overlap,
                                                    epigenomics_center,
                                                    nucleosome_center,
                                                    epigenomics_heatmap_significance_level,
                                                    outputDir,
                                                    jobname,
                                                    numberofSimulations,
                                                    nucleosome_file,
                                                    nucleosome_biosample,
                                                    epigenomics_files_memos,
                                                    epigenomics_biosamples,
                                                    epigenomics_dna_elements,
                                                    sbs_signatures,
                                                    dbs_signatures,
                                                    id_signatures,
                                                    verbose):


    signatures= []
    signatures.extend(sbs_signatures)
    signatures.extend(dbs_signatures)
    signatures.extend(id_signatures)

    os.makedirs(os.path.join(outputDir, jobname, FIGURE, EPIGENOMICSOCCUPANCY, HEATMAPS, TABLES), exist_ok=True)
    os.makedirs(os.path.join(outputDir, jobname, FIGURE, EPIGENOMICSOCCUPANCY, HEATMAPS, DETAILED), exist_ok=True)
    os.makedirs(os.path.join(outputDir, jobname, FIGURE, EPIGENOMICSOCCUPANCY, HEATMAPS, EXCEL_FILES), exist_ok=True)

    heatmaps_output_dir = os.path.join(outputDir, jobname, FIGURE, EPIGENOMICSOCCUPANCY, HEATMAPS)
    heatmaps_detailed_output_dir = os.path.join(outputDir, jobname, FIGURE, EPIGENOMICSOCCUPANCY,HEATMAPS,DETAILED)

    # Step1 Calculate p value using z-test
    # Epigenomics Signatures
    # Epigenomics All Mutations (SUBS, INDELS, DINUCS)
    # Nucleosome Signatures
    # Nucleosome All Mutations (SUBS, INDELS, DINUCS)
    # complete_list
    #[jobname, signature, biosample, dna_element, avg_real_signal, avg_sim_signal, fold_change, min_sim_signal,
    #  max_sim_signal, pvalue, num_of_sims, num_of_sims_with_not_nan_avgs, real_data_avg_count, sim_avg_count,
    #  list(simulationsHorizontalMeans)]
    step1_p_value_df, step1_signature2Biosample2DNAElement2PValueDict= step1_calculate_p_value(fold_change_window_size,
                        epigenomics_center,
                        nucleosome_center,
                        outputDir,
                        jobname,
                        numberofSimulations,
                        nucleosome_file,
                        nucleosome_biosample,
                        epigenomics_files_memos,
                        epigenomics_biosamples,
                        signatures,
                        heatmaps_output_dir)

    #Plot heatmaps using step1 under epigenomics_occupancy/heatmaps/detailed/
    #Plot heatmap unfiltered --- row biosample, columns DNA elements (unpooled, uncombined) under detailed/
    plot_heatmap_one_row_only_for_each_biosample_given_signature(step1_signature2Biosample2DNAElement2PValueDict, jobname, heatmaps_detailed_output_dir, verbose)

    # Step2 uses Step1
    # Combine p values using Fisher's method
    # Pool for epigenomics_dna_elements only
    # Filter Step1 such that fold_change is not (nan,None), p_value is not (nan,None), real_data_avg_count>=100 and get Step2
    step2_combined_p_value_df, step2_signature2biosample2dna_element2combined_p_value_list_dict, step2_signature2biosample2dna_element2avg_fold_change_dict = step2_combine_p_value(step1_signature2Biosample2DNAElement2PValueDict,
                                                                                                                                                         heatmaps_output_dir,
                                                                                                                                                         combine_p_values_method,
                                                                                                                                                         num_of_real_data_avg_overlap,
                                                                                                                                                         nucleosome_file,
                                                                                                                                                         epigenomics_dna_elements)

    #Plot heatmaps using step2 under epigenomics_occupancy/heatmaps/detailed
    #Rows biosamples
    #Columns dna_elements
    plot_heatmap_rows_biosamples_columns_pooled_DNA_elements(step2_signature2biosample2dna_element2avg_fold_change_dict,
                                                             epigenomics_dna_elements,
                                                             jobname,
                                                             heatmaps_detailed_output_dir,
                                                             verbose)

    # Step3 uses Step1
    # Combine p values using Fisher's method
    # Pool for biosample and epigenomics_dna_elements
    # Filter Step1 such that fold_change is not (nan,None), p_value is not (nan,None), real_data_avg_count>=num_of_real_data_avg_overlap and get Step3
    step3_combined_p_value_df, step3_signature2dna_element2combined_p_value_list_dict, step3_signature2dna_element2avg_fold_change_dict = step3_combine_p_value(step1_signature2Biosample2DNAElement2PValueDict,
                          heatmaps_output_dir,
                          combine_p_values_method,
                          num_of_real_data_avg_overlap,
                          nucleosome_file,
                          epigenomics_dna_elements)

    # Plot Step3 heatmaps: Step3 have combined p-values pooled for  biosample and epigenomics_dna_elements but multiple testing correction
    step3_sbs_signature2dna_element2avg_fold_change_dict, \
    step3_dbs_signature2dna_element2avg_fold_change_dict, \
    step3_id_signature2dna_element2avg_fold_change_dict = breakdown_signatures(step3_signature2dna_element2avg_fold_change_dict,sbs_signatures,dbs_signatures,id_signatures)

    if any(step3_sbs_signature2dna_element2avg_fold_change_dict):
        plot_heatmap_rows_signatures_columns_pooled_DNA_elements(step3_sbs_signature2dna_element2avg_fold_change_dict,
                                                                 None,
                                                                 epigenomics_dna_elements,
                                                                 jobname,
                                                                 SBS,
                                                                 heatmaps_detailed_output_dir,
                                                                 "Step3",
                                                                 verbose)

    if any(step3_dbs_signature2dna_element2avg_fold_change_dict):
        plot_heatmap_rows_signatures_columns_pooled_DNA_elements(step3_dbs_signature2dna_element2avg_fold_change_dict,
                                                                 None,
                                                                 epigenomics_dna_elements,
                                                                 jobname,
                                                                 DBS,
                                                                 heatmaps_detailed_output_dir,
                                                                 "Step3",
                                                                 verbose)

    if any(step3_id_signature2dna_element2avg_fold_change_dict):
        plot_heatmap_rows_signatures_columns_pooled_DNA_elements(step3_id_signature2dna_element2avg_fold_change_dict,
                                                                 None,
                                                                 epigenomics_dna_elements,
                                                                 jobname,
                                                                 ID,
                                                                 heatmaps_detailed_output_dir,
                                                                 "Step3",
                                                                 verbose)

    # Step4 Correct combined p values
    # combined p value list
    # [fold_change_list,avg_fold_change,p_value_list,combined_p_value]
    step4_q_value_df,step4_signature2dna_element2q_value_list_dict = step4_apply_multiple_tests_correction(step3_signature2dna_element2combined_p_value_list_dict,heatmaps_output_dir)

    #Step5
    #Filter using q values (combined_q_value<=significance_level and (avg_fold_change>=enriched_fold_change or avg_fold_change<=depleted_fold_change))
    # (signature, cancer_type, dna_element) with combined q_value <= 0.01 and (avg_fold_change >= 1.1 or <=0.9)
    #[fold_change_list, avg_fold_change, q_value_list, combined_q_value]
    step5_filtered_q_value_df,step5_signature2dna_element2average_fold_changedict,signature2dna_element2significancedict = step5_filter_signature_dna_element(step4_signature2dna_element2q_value_list_dict,
                                                                                                                                                            heatmaps_output_dir,
                                                                                                                                                            epigenomics_heatmap_significance_level)

    # Plot heatmaps
    step5_sbs_signature2dna_element2avg_fold_change_dict, \
    step5_dbs_signature2dna_element2avg_fold_change_dict, \
    step5_id_signature2dna_element2avg_fold_change_dict = breakdown_signatures(step5_signature2dna_element2average_fold_changedict,sbs_signatures,dbs_signatures,id_signatures)

    if any(step5_sbs_signature2dna_element2avg_fold_change_dict):
        plot_heatmap_rows_signatures_columns_pooled_DNA_elements(step5_sbs_signature2dna_element2avg_fold_change_dict,
                                                                 signature2dna_element2significancedict,
                                                                 epigenomics_dna_elements,
                                                                 jobname,
                                                                 SBS,
                                                                 heatmaps_output_dir,
                                                                 'Final',
                                                                 verbose)

    if any(step5_dbs_signature2dna_element2avg_fold_change_dict):
        plot_heatmap_rows_signatures_columns_pooled_DNA_elements(step5_dbs_signature2dna_element2avg_fold_change_dict,
                                                                 signature2dna_element2significancedict,
                                                                 epigenomics_dna_elements,
                                                                 jobname,
                                                                 DBS,
                                                                 heatmaps_output_dir,
                                                                 'Final',
                                                                 verbose)

    if any(step5_id_signature2dna_element2avg_fold_change_dict):
        plot_heatmap_rows_signatures_columns_pooled_DNA_elements(step5_id_signature2dna_element2avg_fold_change_dict,
                                                                 signature2dna_element2significancedict,
                                                                 epigenomics_dna_elements,
                                                                 jobname,
                                                                 ID,
                                                                 heatmaps_output_dir,
                                                                 'Final',
                                                                 verbose)

    #write excel files
    excel_file_path = os.path.join(heatmaps_output_dir,EXCEL_FILES,'Epigenomics_Occupancy.xlsx')
    df_list = [step1_p_value_df,step2_combined_p_value_df, step3_combined_p_value_df, step4_q_value_df, step5_filtered_q_value_df]
    sheet_list = ['step1_p_value', 'step2_combined_p_value', 'step3_combined_p_value', 'step4_q_value', 'step5_filtered_q_value']
    write_excel_file(df_list, sheet_list, excel_file_path)
########################################################


########################################################
def breakdown_signatures(signature2dna_element2avg_fold_change_dict,sbs_signatures,dbs_signatures,id_signatures):
    sbs_signature2dna_element2avg_fold_change_dict={}
    dbs_signature2dna_element2avg_fold_change_dict={}
    id_signature2dna_element2avg_fold_change_dict={}

    for signature in signature2dna_element2avg_fold_change_dict:
        if signature in sbs_signatures:
            sbs_signature2dna_element2avg_fold_change_dict[signature]=signature2dna_element2avg_fold_change_dict[signature]
        elif signature in dbs_signatures:
            dbs_signature2dna_element2avg_fold_change_dict[signature]=signature2dna_element2avg_fold_change_dict[signature]
        elif signature in id_signatures:
            id_signature2dna_element2avg_fold_change_dict[signature]=signature2dna_element2avg_fold_change_dict[signature]

    return sbs_signature2dna_element2avg_fold_change_dict, dbs_signature2dna_element2avg_fold_change_dict, id_signature2dna_element2avg_fold_change_dict
########################################################

########################################################
#Used for search for dna_elements
def get_dna_element(dna_element_long,epigenomics_dna_elements,nucleosome_file):
    if epigenomics_dna_elements:
        for dna_element in epigenomics_dna_elements:
            dna_element_with_underscore= "_%s" %dna_element
            if dna_element_with_underscore in dna_element_long:
                return dna_element
    if (nucleosome_file is not None) and (dna_element_long in nucleosome_file):
        return NUCLEOSOME_DNA_ELEMENT
    print('Put ---- %s --- into epigenomics_dna_elements or epigenomics_biosamples' %(dna_element_long))
    return UNDECLARED
########################################################


########################################################
#Used for search for biosamples
def get_biosample(file_memo,biosample_list,nucleosome_file):
    if biosample_list:
        for biosample in biosample_list:
            biosample_with_underscore= "_%s_" %biosample
            if biosample_with_underscore in file_memo:
                return biosample
    if (nucleosome_file is not None) and (file_memo in nucleosome_file):
        return NUCLEOSOME_DNA_ELEMENT
    print('Put ---- %s --- into epigenomics_dna_elements or epigenomics_biosamples' %(file_memo))
    return UNDECLARED
########################################################


########################################################
#Sep24 2020 #Enrichment is done in this function.
#Always ztest
#one sample or two_sample?
#I decided to use one sample because for simulations I will get vertical vector and average of that vertical vector  must be equal to avg_simulated_signal, there is a way to self verification
#Comparing one mean with means of n simulations gives a more realistic p-value.
#In case of comparison of two samples, ztest and ttest gives either 0 or very low p-values.
def calculate_fold_change_real_over_sim(center,
                                        plusorMinus,
                                        output_dir,
                                        jobname,
                                        numberofSimulations,
                                        signature,
                                        nucleosome_file,
                                        nucleosome_biosample,
                                        epigenomics_file_memo,
                                        epigenomics_biosamples,
                                        occupancy_type):

    avg_real_signal=None
    avg_sim_signal=None
    fold_change = None
    min_sim_signal=None
    max_sim_signal=None
    pvalue=None
    num_of_sims=None
    num_of_sims_with_not_nan_avgs=None
    real_data_avg_count=None
    sim_avg_count=None
    simulationsHorizontalMeans = None

    biosample=None
    dna_element=None
    dna_element_to_be_read=None

    if occupancy_type==EPIGENOMICSOCCUPANCY:
        biosample=get_biosample(epigenomics_file_memo,epigenomics_biosamples,nucleosome_file)
        dna_element=epigenomics_file_memo
        dna_element_to_be_read = epigenomics_file_memo
    elif occupancy_type==NUCLEOSOMEOCCUPANCY:
        biosample = nucleosome_biosample
        dna_element = os.path.basename(nucleosome_file)
        dna_element_to_be_read = None

    start=center-plusorMinus
    end=center+plusorMinus+1

    ####################################################################################################################
    # SBS1_sim1_ENCFF330CCJ_osteoblast_H3K79me2-human_AverageSignalArray.txt
    if (signature==AGGREGATEDSUBSTITUTIONS) or (signature==AGGREGATEDDINUCS) or (signature==AGGREGATEDINDELS):
        real_data_avg_signal_array = readData(None, None, signature, output_dir, jobname, occupancy_type, dna_element_to_be_read,AVERAGE_SIGNAL_ARRAY)
    else:
        real_data_avg_signal_array = readData(None, signature, SIGNATUREBASED, output_dir, jobname, occupancy_type, dna_element_to_be_read,AVERAGE_SIGNAL_ARRAY)

    if real_data_avg_signal_array is not None:
        #If there is nan in the list np.mean returns nan.
        avg_real_signal = np.nanmean(real_data_avg_signal_array[start:end])
    ####################################################################################################################

    ####################################################################################################################
    #Read accumulated_count_array
    if (signature==AGGREGATEDSUBSTITUTIONS) or (signature==AGGREGATEDDINUCS) or (signature==AGGREGATEDINDELS):
        real_data_accumulated_count_array = readData(None, None, signature, output_dir, jobname, occupancy_type, dna_element_to_be_read,ACCUMULATED_COUNT_ARRAY)
    else:
        real_data_accumulated_count_array = readData(None, signature, SIGNATUREBASED, output_dir, jobname, occupancy_type, dna_element_to_be_read,ACCUMULATED_COUNT_ARRAY)

    if real_data_accumulated_count_array is not None:
        #If there is nan in the list np.mean returns nan.
        real_data_avg_count = np.nanmean(real_data_accumulated_count_array[start:end])
    ####################################################################################################################


    ####################################################################################################################
    if (numberofSimulations > 0):
        if (signature==AGGREGATEDSUBSTITUTIONS) or (signature==AGGREGATEDDINUCS) or (signature==AGGREGATEDINDELS):
            listofSimulationsSignatureBased = readDataForSimulations(None, None, signature, output_dir,jobname, numberofSimulations, occupancy_type,dna_element_to_be_read,AVERAGE_SIGNAL_ARRAY)
        else:
            listofSimulationsSignatureBased = readDataForSimulations(None, signature, SIGNATUREBASED, output_dir,jobname, numberofSimulations, occupancy_type,dna_element_to_be_read,AVERAGE_SIGNAL_ARRAY)

        if ((listofSimulationsSignatureBased is not None) and listofSimulationsSignatureBased):
            #This is the simulations data
            stackedSimulationsSignatureBased = np.vstack(listofSimulationsSignatureBased)
            (rows, cols) = stackedSimulationsSignatureBased.shape
            num_of_sims = rows

            #One sample way
            # print('stackedSimulationsSignatureBased.shape')
            # print(stackedSimulationsSignatureBased.shape)
            stackedSimulationsSignatureBased_of_interest=stackedSimulationsSignatureBased[:,start:end]
            # print('stackedSimulationsSignatureBased_of_interest.shape')
            # print(stackedSimulationsSignatureBased_of_interest.shape)

            #Get rid of rows with all nans
            stackedSimulationsSignatureBased_of_interest=stackedSimulationsSignatureBased_of_interest[~np.isnan(stackedSimulationsSignatureBased_of_interest).all(axis=1)]

            #Take mean row-wise
            simulationsHorizontalMeans = np.nanmean(stackedSimulationsSignatureBased_of_interest, axis=1)
            avg_sim_signal = np.nanmean(simulationsHorizontalMeans)
            min_sim_signal = np.nanmin(simulationsHorizontalMeans)
            max_sim_signal = np.nanmax(simulationsHorizontalMeans)
            # print('avg_sim_signal:%f' %(avg_sim_signal))

            # print('%s %s %s Number of nans in simulationsHorizontalMeans: %d' %(signature, jobname, dna_element,len(np.argwhere(np.isnan(simulationsHorizontalMeans)))))
            # print('Before')
            # print('simulationsHorizontalMeans.shape')
            # print(simulationsHorizontalMeans.shape)

            #Get rid of nans in simulationsHorizontalMeans
            #simulationsHorizontalMeans is used in p-value calculation
            simulationsHorizontalMeans = simulationsHorizontalMeans[~np.isnan(simulationsHorizontalMeans)]
            # print('After')
            # print('simulationsHorizontalMeans.shape')
            # print(simulationsHorizontalMeans.shape)
            num_of_sims_with_not_nan_avgs=simulationsHorizontalMeans.shape[0]
            # print('number of not nan simulations:%d' %num_of_sims_with_not_nan_avgs)
            # print('%s %s %s Number of nans in simulationsHorizontalMeans: %d' %(signature, jobname, dna_element,len(np.argwhere(np.isnan(simulationsHorizontalMeans)))))
            # print(np.argwhere(np.isnan(simulationsHorizontalMeans)))
            # print(simulationsHorizontalMeans)
    ####################################################################################################################

    ####################################################################################################################
    if (numberofSimulations > 0):
        if (signature==AGGREGATEDSUBSTITUTIONS) or (signature==AGGREGATEDDINUCS) or (signature==AGGREGATEDINDELS):
            listofSimulationsSignatureBased = readDataForSimulations(None, None, signature, output_dir,jobname, numberofSimulations, occupancy_type,dna_element_to_be_read,ACCUMULATED_COUNT_ARRAY)
        else:
            listofSimulationsSignatureBased = readDataForSimulations(None, signature, SIGNATUREBASED, output_dir,jobname, numberofSimulations, occupancy_type,dna_element_to_be_read,ACCUMULATED_COUNT_ARRAY)

        if ((listofSimulationsSignatureBased is not None) and listofSimulationsSignatureBased):
            #This is the simulations data
            stackedSimulationsSignatureBased = np.vstack(listofSimulationsSignatureBased)
            (rows, cols) = stackedSimulationsSignatureBased.shape

            #One sample way
            stackedSimulationsSignatureBased_of_interest=stackedSimulationsSignatureBased[:,start:end]

            #Get rid of rows with all nans
            stackedSimulationsSignatureBased_of_interest=stackedSimulationsSignatureBased_of_interest[~np.isnan(stackedSimulationsSignatureBased_of_interest).all(axis=1)]

            #Take mean row-wise
            simulationsHorizontalCountMeans = np.nanmean(stackedSimulationsSignatureBased_of_interest, axis=1)
            sim_avg_count = np.nanmean(simulationsHorizontalCountMeans)
    ####################################################################################################################


    if (avg_real_signal is not None) and (avg_sim_signal is not None):
        fold_change = avg_real_signal / avg_sim_signal
        # print('avg_real_signal:%f' % (avg_real_signal))
        # print('avg_sim_signal:%f' % (avg_sim_signal))
        # print('fold change:%f' % (fold_change))

        if (simulationsHorizontalMeans is not None):
            # zstat, pvalue_ztest_1sample = ztest(simulationsHorizontalMeans, value=avg_real_signal)
            # if there is only one simulation mean in simulationsHorizontalMeans, then pvalue is nan
            zstat, pvalue = ztest(simulationsHorizontalMeans,[avg_real_signal])
            # print('%s %s %s  avg_real_signal:%f avg_sim_signal:%f min_sim_signal:%f max_sim_signal:%f fold_change:%f p_value: %.2E' %(signature, jobname, dna_element,avg_real_signal,avg_sim_signal,min_sim_signal,max_sim_signal, fold_change, Decimal(pvalue) ))
            # print('###############################################################################################################################')

        return [jobname, signature, biosample, dna_element, avg_real_signal, avg_sim_signal, fold_change, min_sim_signal, max_sim_signal, pvalue, num_of_sims, num_of_sims_with_not_nan_avgs, real_data_avg_count, sim_avg_count, list(simulationsHorizontalMeans)]
    else:
        if (simulationsHorizontalMeans is not None):
            return [jobname, signature, biosample, dna_element, avg_real_signal, avg_sim_signal, fold_change, min_sim_signal, max_sim_signal, pvalue, num_of_sims, num_of_sims_with_not_nan_avgs, real_data_avg_count, sim_avg_count, list(simulationsHorizontalMeans)]
        else:
            return [jobname, signature, biosample, dna_element, avg_real_signal, avg_sim_signal, fold_change, min_sim_signal, max_sim_signal, pvalue, num_of_sims, num_of_sims_with_not_nan_avgs, real_data_avg_count, sim_avg_count, None]
########################################################


########################################################
# complete list with p value
#[jobname,
# signature,
# biosample,
# dna_element,
# avg_real_signal,
# avg_sim_signal,
# fold_change,
# min_sim_signal,
# max_sim_signal,
# pvalue,
# num_of_sims,
# num_of_sims_with_not_nan_avgs,
# real_data_avg_count,
# sim_avg_count,
# list(simulationsHorizontalMeans)]
def write_dictionary_as_dataframe_step1_p_value(step1_signature2Biosample2DNAElement2PValueDict,filepath):
    L = sorted([(complete_list[0], signature, biosample, dna_element, complete_list[4], complete_list[5], complete_list[6], complete_list[7], complete_list[8], complete_list[9], complete_list[10],complete_list[11], complete_list[12], complete_list[13], complete_list[14])
                for signature, a in step1_signature2Biosample2DNAElement2PValueDict.items()
                  for biosample, b in a.items()
                   for dna_element, complete_list in b.items()])
    df = pd.DataFrame(L, columns=['cancer_type', 'signature', 'biosample', 'dna_element', 'avg_real_signal','avg_simulated_signal', 'fold_change', 'min_sim_signal', 'max_sim_signal','p_value', 'num_of_sims', 'num_of_sims_with_not_nan_avgs', 'real_data_avg_count', 'sim_avg_count', 'sim_signals'])
    df.to_csv(filepath, sep='\t', header=True, index=False)

    return df
########################################################


########################################################
#Combined p value
#[avg_real_signal_list,avg_sim_signal_list,fold_change_list,avg_fold_change,p_value_list,combined_p_value]
def write_dictionary_as_dataframe_step2_combined_p_value(signature2biosample2pooled_dna_element2combined_p_value_list_dict,filepath):
    L = sorted([(signature, biosample, dna_element, combined_p_value_list[0], combined_p_value_list[1], combined_p_value_list[2], combined_p_value_list[3], combined_p_value_list[4], combined_p_value_list[5], combined_p_value_list[6])
                for signature, a in signature2biosample2pooled_dna_element2combined_p_value_list_dict.items()
                for biosample, b in a.items()
                for dna_element, combined_p_value_list in b.items()])
    df = pd.DataFrame(L, columns=['signature', 'biosample','dna_element', 'dna_element_long_list' ,'avg_real_signal_list','avg_sim_signal_list','fold_change_list', 'avg_fold_change', 'p_value_list', 'combined_p_value'])
    df.to_csv(filepath, sep='\t', header=True, index=False)

    return df
########################################################


########################################################
#Combined p value
#[avg_real_signal_list,avg_sim_signal_list,fold_change_list,avg_fold_change,p_value_list,combined_p_value]
def write_dictionary_as_dataframe_step3_combined_p_value(step2_signature2dna_element2combined_p_value_list_dict,filepath):
    L = sorted([(signature, dna_element, combined_p_value_list[0], combined_p_value_list[1],combined_p_value_list[2], combined_p_value_list[3], combined_p_value_list[4], combined_p_value_list[5],combined_p_value_list[6])
                for signature, a in step2_signature2dna_element2combined_p_value_list_dict.items()
                  for dna_element, combined_p_value_list in a.items()])
    df = pd.DataFrame(L, columns=['signature', 'dna_element', 'dna_element_long_list', 'avg_real_signal_list', 'avg_sim_signal_list', 'fold_change_list', 'avg_fold_change' , 'p_value_list', 'combined_p_value'])
    df.to_csv(filepath, sep='\t', header=True, index=False)

    return df
########################################################


########################################################
# Q Value List
#[avg_real_signal_list, avg_sim_signal_list, fold_change_list, avg_fold_change, p_value_list, combined_p_value]
def write_dictionary_as_dataframe_step4_q_value(step3_signature2dna_element2q_value_list_dict,filepath):
    L = sorted([(signature, dna_element, q_value_list[0], q_value_list[1], q_value_list[2], q_value_list[3], q_value_list[4], q_value_list[5], q_value_list[6], q_value_list[7])
                for signature, a in step3_signature2dna_element2q_value_list_dict.items()
                  for dna_element, q_value_list in a.items()])
    df = pd.DataFrame(L, columns=['signature', 'dna_element',
                                  'dna_element_long_list', 'avg_real_signal_list', 'avg_sim_signal_list', 'fold_change_list', 'avg_fold_change', 'p_value_list', 'combined_p_value','q_value'])
    df.to_csv(filepath, sep='\t', header=True, index=False)

    return df
########################################################


########################################################
#[avg_real_signal_list, avg_sim_signal_list, fold_change_list, avg_fold_change, p_value_list, combined_p_value,q_value]
def write_dictionary_as_dataframe_step5_filtered_q_value(step4_signature2dna_element2filtered_q_list_dict,epigenomics_heatmap_significance_level,filepath):
    L = sorted([(signature, dna_element, filtered_q_value_list[0], filtered_q_value_list[1], filtered_q_value_list[2], filtered_q_value_list[3], filtered_q_value_list[4], filtered_q_value_list[5], filtered_q_value_list[6], filtered_q_value_list[7])
                for signature, a in step4_signature2dna_element2filtered_q_list_dict.items()
                  for dna_element, filtered_q_value_list in a.items()])
    df = pd.DataFrame(L, columns=['signature', 'dna_element',
                                  'dna_element_long_list','avg_real_signal_list','avg_sim_signal_list','fold_change_list', 'avg_fold_change', 'p_value_list', 'combined_p_value','filtered_q_value'])

    df=df[df['filtered_q_value']<=epigenomics_heatmap_significance_level]

    df.to_csv(filepath, sep='\t', header=True, index=False)

    return df
########################################################


########################################################
#Sep24, 2020
# Step1
# Epigenomics Signatures
# Epigenomics AGGREGATEDSUBSTITUTIONS AGGREGATEDDINUCS AGGREGATEDINDELS
# Nucleosome Signatures
# Nucleosome AGGREGATEDSUBSTITUTIONS AGGREGATEDDINUCS AGGREGATEDINDELS
def step1_calculate_p_value(fold_change_window_size,
                        epigenomics_center,
                        nucleosome_center,
                        output_dir,
                        jobname,
                        numberofSimulations,
                        nucleosome_file,
                        nucleosome_biosample,
                        epigenomics_files_memos,
                        epigenomics_biosamples,
                        signatures,
                        heatmaps_output_dir):

    #######################################################################
    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)
    #######################################################################

    ########################################################
    signature2Biosample2DNAElement2PValueDict = {}
    ########################################################

    ########################################################
    plusorMinus = fold_change_window_size//2
    ########################################################

    ########################################################
    def update_dictionary(complete_list):
        # complete_list:[jobname, signature, biosample, dna_element, avg_real_signal, avg_sim_signal, fold_change, min_sim_signal,
        #  max_sim_signal, pvalue, num_of_sims, num_of_sims_with_not_nan_avgs, real_data_avg_count, sim_avg_count,
        #  list(simulationsHorizontalMeans)]
        signature=complete_list[1]
        biosample=complete_list[2]
        dna_element=complete_list[3]

        if signature in signature2Biosample2DNAElement2PValueDict:
            if biosample in signature2Biosample2DNAElement2PValueDict[signature]:
                if dna_element in signature2Biosample2DNAElement2PValueDict[signature][biosample]:
                    print('There is a problem: %s %s %s' %(signature,biosample,dna_element))
                else:
                    signature2Biosample2DNAElement2PValueDict[signature][biosample][dna_element] = complete_list
            else:
                signature2Biosample2DNAElement2PValueDict[signature][biosample] = {}
                signature2Biosample2DNAElement2PValueDict[signature][biosample][dna_element] = complete_list
        else:
            signature2Biosample2DNAElement2PValueDict[signature] = {}
            signature2Biosample2DNAElement2PValueDict[signature][biosample] = {}
            signature2Biosample2DNAElement2PValueDict[signature][biosample][dna_element] = complete_list
    ########################################################

    ########################################################################
    for signature in signatures:
        #Epigenomics
        occupancy_type=EPIGENOMICSOCCUPANCY
        for epigenomics_file_memo in epigenomics_files_memos:
            pool.apply_async(calculate_fold_change_real_over_sim,
                             args=(epigenomics_center,plusorMinus,output_dir,jobname,numberofSimulations,signature,nucleosome_file,nucleosome_biosample,epigenomics_file_memo,epigenomics_biosamples,occupancy_type,),
                             callback=update_dictionary)
        #Nucleosome
        # if data files are ready it returns otherwise it returns None
        occupancy_type=NUCLEOSOMEOCCUPANCY
        epigenomics_file_memo=None
        pool.apply_async(calculate_fold_change_real_over_sim,
                 args=(nucleosome_center,plusorMinus,output_dir,jobname,numberofSimulations,signature,nucleosome_file,nucleosome_biosample,epigenomics_file_memo,epigenomics_biosamples,occupancy_type,),
                 callback=update_dictionary)
    ########################################################################

    #######################################################################
    pool.close()
    pool.join()
    #######################################################################

    # print('##############################################################')
    # print('Step1 Getting p-values')
    # Write dictionary as a dataframe
    df_filename = 'Step1_Signature_Biosample_DNAElement_PValue.txt'
    filepath = os.path.join(heatmaps_output_dir,TABLES, df_filename)
    step1_p_value_df = write_dictionary_as_dataframe_step1_p_value(signature2Biosample2DNAElement2PValueDict,filepath)
    # print('##############################################################')

    return step1_p_value_df, signature2Biosample2DNAElement2PValueDict
########################################################


########################################################
#November 10, 2020
def step2_combine_p_value(signature2Biosample2DNAElement2PValueDict,
                          heatmaps_output_dir,
                          combine_p_values_method,
                          num_of_real_data_avg_overlap,
                          nucleosome_file,
                          epigenomics_dna_elements):

    #Fill and return this dictionary
    signature2biosample2pooled_dna_element2combined_p_value_list_dict = {}
    signature2biosample2pooled_dna_element2avg_fold_change_dict = {}

    #Pooling for biosample and dna_element combine q_values
    for signature in signature2Biosample2DNAElement2PValueDict:
        for biosample in signature2Biosample2DNAElement2PValueDict[signature]:
            # dna_element_long <-- epigenomics_file_memo
            # dna_element_long <-- os.path.basename(nucleosome_file)
            for dna_element_long in signature2Biosample2DNAElement2PValueDict[signature][biosample]:
                dna_element = get_dna_element(dna_element_long,epigenomics_dna_elements,nucleosome_file)
                complete_list = signature2Biosample2DNAElement2PValueDict[signature][biosample][dna_element_long]

                # p value complete list has [signature, cancer_type, biosample, dna_element, avg_real_signal, avg_sim_signal, fold_change, min_sim_signal, max_sim_signal, pvalue, num_of_sims, num_of_sims_with_not_nan_avgs, real_data_avg_count, list(simulationsHorizontalMeans)]
                # complete_list:[
                # 0 jobname,
                # 1 signature,
                # 2 biosample,
                # 3 dna_element,
                # 4 avg_real_signal,
                # 5 avg_sim_signal,
                # 6 fold_change,
                # 7 min_sim_signal,
                # 8 max_sim_signal,
                # 9 pvalue,
                # 10 num_of_sims,
                # 11 num_of_sims_with_not_nan_avgs,
                # 12 real_data_avg_count,
                # 13 sim_avg_count,
                # 14 list(simulationsHorizontalMeans)]
                avg_real_signal=complete_list[4]
                avg_sim_signal = complete_list[5]
                fold_change=complete_list[6]
                p_value=complete_list[9]
                real_data_avg_count=complete_list[12]
                sim_data_avg_count=complete_list[13]

                if ((fold_change is not None) and (not np.isnan(np.array([fold_change], dtype=np.float)).any()) and (str(fold_change) != 'nan')) and \
                        ((p_value is not None) and (not np.isnan(np.array([p_value], dtype=np.float)).any()) and (str(p_value) != 'nan')) and \
                        ((real_data_avg_count>=num_of_real_data_avg_overlap) or (sim_data_avg_count>=num_of_real_data_avg_overlap)):

                    if signature in signature2biosample2pooled_dna_element2combined_p_value_list_dict:
                        if biosample in signature2biosample2pooled_dna_element2combined_p_value_list_dict[signature]:
                            if dna_element in signature2biosample2pooled_dna_element2combined_p_value_list_dict[signature][biosample]:

                                #Add dna_element_long to the list
                                signature2biosample2pooled_dna_element2combined_p_value_list_dict[signature][biosample][dna_element][0].append(dna_element_long)
                                #Add avg_real_signal to the list
                                signature2biosample2pooled_dna_element2combined_p_value_list_dict[signature][biosample][dna_element][1].append(avg_real_signal)
                                # Add avg_sim_signal to the list
                                signature2biosample2pooled_dna_element2combined_p_value_list_dict[signature][biosample][dna_element][2].append(avg_sim_signal)
                                #Add to the fold change list
                                signature2biosample2pooled_dna_element2combined_p_value_list_dict[signature][biosample][dna_element][3].append(fold_change)
                                #Add to the p value list
                                signature2biosample2pooled_dna_element2combined_p_value_list_dict[signature][biosample][dna_element][4].append(p_value)
                            else:
                                signature2biosample2pooled_dna_element2combined_p_value_list_dict[signature][biosample][dna_element]=[[dna_element_long], [avg_real_signal], [avg_sim_signal],[fold_change],[p_value]]
                        else:
                            signature2biosample2pooled_dna_element2combined_p_value_list_dict[signature][biosample] = {}
                            signature2biosample2pooled_dna_element2combined_p_value_list_dict[signature][biosample][dna_element] = [[dna_element_long], [avg_real_signal], [avg_sim_signal],[fold_change], [p_value]]
                    else:
                        signature2biosample2pooled_dna_element2combined_p_value_list_dict[signature] = {}
                        signature2biosample2pooled_dna_element2combined_p_value_list_dict[signature][biosample] = {}
                        signature2biosample2pooled_dna_element2combined_p_value_list_dict[signature][biosample][dna_element] = [[dna_element_long], [avg_real_signal], [avg_sim_signal], [fold_change], [p_value]]

    for signature in signature2biosample2pooled_dna_element2combined_p_value_list_dict:
        for biosample in signature2biosample2pooled_dna_element2combined_p_value_list_dict[signature]:
            for dna_element in signature2biosample2pooled_dna_element2combined_p_value_list_dict[signature][biosample]:

                dna_element_long_list=signature2biosample2pooled_dna_element2combined_p_value_list_dict[signature][biosample][dna_element][0]
                avg_real_signal_list=signature2biosample2pooled_dna_element2combined_p_value_list_dict[signature][biosample][dna_element][1]
                avg_sim_signal_list=signature2biosample2pooled_dna_element2combined_p_value_list_dict[signature][biosample][dna_element][2]
                fold_change_list=signature2biosample2pooled_dna_element2combined_p_value_list_dict[signature][biosample][dna_element][3]
                p_value_list=signature2biosample2pooled_dna_element2combined_p_value_list_dict[signature][biosample][dna_element][4]

                avg_fold_change= np.nanmean(fold_change_list)
                p_values_array=np.asarray(p_value_list)

                try:
                    test_statistic,combined_p_value=scipy.stats.combine_pvalues(p_values_array, method=combine_p_values_method, weights=None)
                except FloatingPointError:
                    print('signature:%s dna_element:%s fold_change_list:%s p_value_list:%s' %(signature,dna_element,fold_change_list,p_value_list))
                    if len(p_value_list)>0:
                        combined_p_value=p_value_list[0]

                signature2biosample2pooled_dna_element2combined_p_value_list_dict[signature][biosample][dna_element]=[dna_element_long_list,
                                                                                                                      avg_real_signal_list,
                                                                                                                      avg_sim_signal_list,
                                                                                                                      fold_change_list,
                                                                                                                      avg_fold_change,
                                                                                                                      p_value_list,
                                                                                                                      combined_p_value]

                if signature in signature2biosample2pooled_dna_element2avg_fold_change_dict:
                    if biosample in signature2biosample2pooled_dna_element2avg_fold_change_dict[signature]:
                        if dna_element in signature2biosample2pooled_dna_element2avg_fold_change_dict[signature][biosample]:
                            print('There is a problem')
                        else:
                            signature2biosample2pooled_dna_element2avg_fold_change_dict[signature][biosample][dna_element]=avg_fold_change
                    else:
                        signature2biosample2pooled_dna_element2avg_fold_change_dict[signature][biosample] = {}
                        signature2biosample2pooled_dna_element2avg_fold_change_dict[signature][biosample][dna_element] = avg_fold_change
                else:
                    signature2biosample2pooled_dna_element2avg_fold_change_dict[signature] = {}
                    signature2biosample2pooled_dna_element2avg_fold_change_dict[signature][biosample] = {}
                    signature2biosample2pooled_dna_element2avg_fold_change_dict[signature][biosample][dna_element] = avg_fold_change

    # Write dictionary as a pandas dataframe
    df_filename = 'Step2_Signature_Biosample_DNAElement_CombinedPValue.txt'
    filepath = os.path.join(heatmaps_output_dir,TABLES, df_filename)
    step2_combined_p_value_df=write_dictionary_as_dataframe_step2_combined_p_value(signature2biosample2pooled_dna_element2combined_p_value_list_dict,filepath)

    return step2_combined_p_value_df, signature2biosample2pooled_dna_element2combined_p_value_list_dict, signature2biosample2pooled_dna_element2avg_fold_change_dict
########################################################


########################################################
def step3_combine_p_value(signature2Biosample2DNAElement2PValueDict,
                          heatmaps_output_dir,
                          combine_p_values_method,
                          num_of_real_data_avg_overlap,
                          nucleosome_file,
                          epigenomics_dna_elements):

    #Fill and return this dictionary
    signature2dna_element2combined_p_value_list_dict = {}
    signature2dna_element2avg_fold_change_dict = {}

    #Pooling for biosample and dna_element combine q_values
    for signature in signature2Biosample2DNAElement2PValueDict:
        for biosample in signature2Biosample2DNAElement2PValueDict[signature]:
            # dna_element_long <-- epigenomics_file_memo
            # dna_element_long <-- os.path.basename(nucleosome_file)
            for dna_element_long in signature2Biosample2DNAElement2PValueDict[signature][biosample]:
                dna_element = get_dna_element(dna_element_long,epigenomics_dna_elements,nucleosome_file)
                complete_list = signature2Biosample2DNAElement2PValueDict[signature][biosample][dna_element_long]

                # p value complete list has [signature, cancer_type, biosample, dna_element, avg_real_signal, avg_sim_signal, fold_change, min_sim_signal, max_sim_signal, pvalue, num_of_sims, num_of_sims_with_not_nan_avgs, real_data_avg_count, list(simulationsHorizontalMeans)]
                # complete_list:[
                # 0 jobname,
                # 1 signature,
                # 2 biosample,
                # 3 dna_element,
                # 4 avg_real_signal,
                # 5 avg_sim_signal,
                # 6 fold_change,
                # 7 min_sim_signal,
                # 8 max_sim_signal,
                # 9 pvalue,
                # 10 num_of_sims,
                # 11 num_of_sims_with_not_nan_avgs,
                # 12 real_data_avg_count,
                # 13 sim_avg_count,
                # 14 list(simulationsHorizontalMeans)]
                avg_real_signal=complete_list[4]
                avg_sim_signal=complete_list[5]
                fold_change=complete_list[6]
                p_value=complete_list[9]
                real_data_avg_count=complete_list[12]
                sim_data_avg_count=complete_list[13]

                if ((fold_change is not None) and (not np.isnan(np.array([fold_change], dtype=np.float)).any()) and (str(fold_change) != 'nan')) and \
                        ((p_value is not None) and (not np.isnan(np.array([p_value], dtype=np.float)).any()) and (str(p_value) != 'nan')) and \
                        (real_data_avg_count>=num_of_real_data_avg_overlap or sim_data_avg_count>=num_of_real_data_avg_overlap):

                    if signature in signature2dna_element2combined_p_value_list_dict:
                        if dna_element in signature2dna_element2combined_p_value_list_dict[signature]:
                            #Add dna_element_long
                            signature2dna_element2combined_p_value_list_dict[signature][dna_element][0].append(dna_element_long)
                            #Add to avg_real_signal_list
                            signature2dna_element2combined_p_value_list_dict[signature][dna_element][1].append(avg_real_signal)
                            #Add to avg_sim_signal_list
                            signature2dna_element2combined_p_value_list_dict[signature][dna_element][2].append(avg_sim_signal)
                            #Add to the fold change list
                            signature2dna_element2combined_p_value_list_dict[signature][dna_element][3].append(fold_change)
                            #Add to the p value list
                            signature2dna_element2combined_p_value_list_dict[signature][dna_element][4].append(p_value)
                        else:
                            signature2dna_element2combined_p_value_list_dict[signature][dna_element]=[[dna_element_long], [avg_real_signal], [avg_sim_signal],[fold_change],[p_value]]
                    else:
                        signature2dna_element2combined_p_value_list_dict[signature] = {}
                        signature2dna_element2combined_p_value_list_dict[signature][dna_element] = [[dna_element_long],[avg_real_signal], [avg_sim_signal], [fold_change], [p_value]]

    for signature in signature2dna_element2combined_p_value_list_dict:
        for dna_element in signature2dna_element2combined_p_value_list_dict[signature]:

            dna_element_long_list=signature2dna_element2combined_p_value_list_dict[signature][dna_element][0]
            avg_real_signal_list=signature2dna_element2combined_p_value_list_dict[signature][dna_element][1]
            avg_sim_signal_list=signature2dna_element2combined_p_value_list_dict[signature][dna_element][2]
            fold_change_list=signature2dna_element2combined_p_value_list_dict[signature][dna_element][3]
            p_value_list=signature2dna_element2combined_p_value_list_dict[signature][dna_element][4]

            avg_fold_change= np.nanmean(fold_change_list)
            p_values_array=np.asarray(p_value_list)

            try:
                test_statistic,combined_p_value=scipy.stats.combine_pvalues(p_values_array, method=combine_p_values_method, weights=None)
            except FloatingPointError:
                print('signature:%s dna_element:%s fold_change_list:%s p_value_list:%s' %(signature,dna_element,fold_change_list,p_value_list))
                if len(p_value_list)>0:
                    combined_p_value=p_value_list[0]

            signature2dna_element2combined_p_value_list_dict[signature][dna_element]=[dna_element_long_list, avg_real_signal_list,avg_sim_signal_list,fold_change_list,avg_fold_change,p_value_list,combined_p_value]

            if signature in signature2dna_element2avg_fold_change_dict:
                if dna_element in signature2dna_element2avg_fold_change_dict[signature]:
                    print('There is a problem')
                else:
                    signature2dna_element2avg_fold_change_dict[signature][dna_element]=avg_fold_change
            else:
                signature2dna_element2avg_fold_change_dict[signature] = {}
                signature2dna_element2avg_fold_change_dict[signature][dna_element] = avg_fold_change

    # Write dictionary as a pandas dataframe
    df_filename = 'Step3_Signature_DNAElement_CombinedPValue.txt'
    filepath = os.path.join(heatmaps_output_dir,TABLES, df_filename)
    step3_combined_p_value_df = write_dictionary_as_dataframe_step3_combined_p_value(signature2dna_element2combined_p_value_list_dict,filepath)

    return step3_combined_p_value_df, signature2dna_element2combined_p_value_list_dict, signature2dna_element2avg_fold_change_dict
########################################################

########################################################
#[dna_element_long_list, avg_real_signal_list, avg_sim_signal_list, fold_change_list, avg_fold_change, p_value_list, combined_p_value]
def step4_apply_multiple_tests_correction(signature2dna_element2combined_p_value_list_dict,heatmaps_output_dir):
    signature2dna_element2q_value_list_dict={}

    all_p_values = []
    all_p_values_element_names = []
    all_FDR_BH_adjusted_p_values = None

    for signature in signature2dna_element2combined_p_value_list_dict:
        for dna_element in signature2dna_element2combined_p_value_list_dict[signature]:
            combined_p_value = signature2dna_element2combined_p_value_list_dict[signature][dna_element][6]
            if (combined_p_value is not None) and (not np.isnan(np.array([combined_p_value], dtype=np.float)).any()) and (str(combined_p_value)!='nan'):
                element_name = (signature, dna_element)
                all_p_values.append(combined_p_value)
                all_p_values_element_names.append(element_name)
            else:
                print('combined_p_value is None or nan: %s %s %s' % (signature, dna_element, signature2dna_element2combined_p_value_list_dict[signature][dna_element]))

    all_p_values_array = np.asarray(all_p_values)

    #If there a p_values in the array
    if(all_p_values_array.size>0):
        rejected, all_FDR_BH_adjusted_p_values, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(all_p_values_array, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

    # print('#######################################')
    # print('len(all_p_values):%d' %len(all_p_values))

    # if ((all_FDR_BH_adjusted_p_values is not None) and (all_FDR_BH_adjusted_p_values.size>0)):
        # print('len(all_FDR_BH_adjusted_p_values):%d' %(len(all_FDR_BH_adjusted_p_values)))
    # print('#######################################')

    for element_index, element_name in enumerate(all_p_values_element_names,0):
        signature, dna_element = element_name
        q_value=all_FDR_BH_adjusted_p_values[element_index]

        # combined_p_value_list
        # [dna_element_long_list,avg_real_signal_list,avg_sim_signal_list,fold_change_list,avg_fold_change,p_value_list,combined_p_value]
        dna_element_long_list=signature2dna_element2combined_p_value_list_dict[signature][dna_element][0]
        avg_real_signal_list=signature2dna_element2combined_p_value_list_dict[signature][dna_element][1]
        avg_sim_signal_list=signature2dna_element2combined_p_value_list_dict[signature][dna_element][2]
        fold_change_list=signature2dna_element2combined_p_value_list_dict[signature][dna_element][3]
        avg_fold_change=signature2dna_element2combined_p_value_list_dict[signature][dna_element][4]
        p_value_list=signature2dna_element2combined_p_value_list_dict[signature][dna_element][5]
        combined_p_value=signature2dna_element2combined_p_value_list_dict[signature][dna_element][6]

        if signature in signature2dna_element2q_value_list_dict:
            if dna_element in signature2dna_element2q_value_list_dict[signature]:
                print('There is a situation')
            else:
                signature2dna_element2q_value_list_dict[signature][dna_element]=[dna_element_long_list,avg_real_signal_list,avg_sim_signal_list,fold_change_list,avg_fold_change,p_value_list,combined_p_value,q_value]
        else:
            signature2dna_element2q_value_list_dict[signature] = {}
            signature2dna_element2q_value_list_dict[signature][dna_element] = [dna_element_long_list,avg_real_signal_list,avg_sim_signal_list,fold_change_list,avg_fold_change,p_value_list,combined_p_value,q_value]

    # Write dictionary as a dataframe
    df_filename = 'Step4_Signature_CancerType_DNAElement_QValue.txt'
    filepath = os.path.join(heatmaps_output_dir, TABLES, df_filename)
    step4_q_value_df=write_dictionary_as_dataframe_step4_q_value(signature2dna_element2q_value_list_dict,filepath)

    return step4_q_value_df,signature2dna_element2q_value_list_dict
########################################################


########################################################
def step5_filter_signature_dna_element(signature2dna_element2q_value_list_dict,heatmaps_output_dir,epigenomics_heatmap_significance_level):
    signature2dna_element2filtered_q_list_dict={}
    signature2dna_element2average_fold_changedict={}
    signature2dna_element2significancedict = {}

    for signature in signature2dna_element2q_value_list_dict:
        for dna_element in signature2dna_element2q_value_list_dict[signature]:
            q_value_list=signature2dna_element2q_value_list_dict[signature][dna_element]

            dna_element_long_list=q_value_list[0]
            avg_real_signal_list=q_value_list[1]
            avg_sim_signal_list=q_value_list[2]
            fold_change_list=q_value_list[3]
            avg_fold_change=q_value_list[4]
            p_value_list=q_value_list[5]
            combined_p_value=q_value_list[6]
            q_value=q_value_list[7]

            #Filter here
            # if (q_value<=significance_level and (avg_fold_change>=enriched_fold_change or avg_fold_change<=depleted_fold_change)):
            # Let's plot all
            if (q_value <= epigenomics_heatmap_significance_level):
                if signature in signature2dna_element2significancedict:
                    if dna_element in signature2dna_element2significancedict[signature]:
                        print('There is a problem')
                    else:
                        signature2dna_element2significancedict[signature][dna_element]="*"
                else:
                    signature2dna_element2significancedict[signature] = {}
                    signature2dna_element2significancedict[signature][dna_element] = "*"


            if signature in signature2dna_element2filtered_q_list_dict:
                    if dna_element in signature2dna_element2filtered_q_list_dict[signature]:
                        print('There is a problem')
                    else:
                        signature2dna_element2filtered_q_list_dict[signature][dna_element]=[dna_element_long_list,avg_real_signal_list, avg_sim_signal_list, fold_change_list, avg_fold_change, p_value_list, combined_p_value, q_value]
                        signature2dna_element2average_fold_changedict[signature][dna_element]=avg_fold_change
            else:
                signature2dna_element2filtered_q_list_dict[signature] = {}
                signature2dna_element2filtered_q_list_dict[signature][dna_element] = [dna_element_long_list,avg_real_signal_list, avg_sim_signal_list, fold_change_list, avg_fold_change, p_value_list, combined_p_value, q_value]
                signature2dna_element2average_fold_changedict[signature] = {}
                signature2dna_element2average_fold_changedict[signature][dna_element] = avg_fold_change

    # Write dictionary as a dataframe
    df_filename = 'Step5_Signature_CancerType_DNAElement_FilteredQValue.txt'
    filepath = os.path.join(heatmaps_output_dir, TABLES, df_filename)
    #Filter rows in write_dictionary_as_dataframe_step5_filtered_q_value
    step5_filtered_q_value_df=write_dictionary_as_dataframe_step5_filtered_q_value(signature2dna_element2filtered_q_list_dict,epigenomics_heatmap_significance_level,filepath)

    return step5_filtered_q_value_df,signature2dna_element2average_fold_changedict,signature2dna_element2significancedict
########################################################



#########################################################
def occupancyAverageSignalFigures(outputDir,jobname,numberofSimulations,sample_based,mutation_types_contexts,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus,verbose,plot_mode):
    if (occupancy_type==NUCLEOSOMEOCCUPANCY):
        ylabel='Average nucleosome signal'
    else:
        #For epigenomics, epigenomics_dir_name can be different than EPIGENOMICSOCCUPANCY
        ylabel='Average epigenomics signal'

    # Initialize these dataframes as empty dataframe
    # We will read these dataframes if there is the corresponding data
    subsSignature_cutoff_numberofmutations_averageprobability_df = pd.DataFrame()
    dinucsSignature_cutoff_numberofmutations_averageprobability_df = pd.DataFrame()
    indelsSignature_cutoff_numberofmutations_averageprobability_df = pd.DataFrame()

    if occupancy_type==NUCLEOSOMEOCCUPANCY:
        os.makedirs(os.path.join(outputDir, jobname, FIGURE, occupancy_type), exist_ok=True)
    elif occupancy_type==EPIGENOMICSOCCUPANCY:
        os.makedirs(os.path.join(outputDir, jobname, FIGURE, occupancy_type, PLOTS), exist_ok=True)

    # Read necessary dictionaries
    mutationtype_numberofmutations_numberofsamples_sampleslist_df = pd.read_csv(os.path.join(outputDir, jobname, DATA,Table_MutationType_NumberofMutations_NumberofSamples_SamplesList_Filename),sep='\t', header=0,dtype={'mutation_type': str,'number_of_mutations': np.int32})

    for mutation_type_context in mutation_types_contexts:
        if (mutation_type_context in SBS_CONTEXTS):
            subsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(os.path.join(outputDir, jobname, DATA,Table_SubsSignature_Cutoff_NumberofMutations_AverageProbability_Filename),sep='\t', header=0,dtype={'cutoff': np.float32,'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})
    if (DBS in mutation_types_contexts):
        dinucsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(os.path.join(outputDir, jobname, DATA,Table_DinucsSignature_Cutoff_NumberofMutations_AverageProbability_Filename),sep='\t', header=0,dtype={'cutoff': np.float32,'signature': str,'number_of_mutations': np.int32,'average_probability': np.float32})
    if (ID in mutation_types_contexts):
        indelsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(os.path.join(outputDir, jobname, DATA,Table_IndelsSignature_Cutoff_NumberofMutations_AverageProbability_Filename), sep='\t',header=0, dtype={'cutoff': np.float32, 'signature': str, 'number_of_mutations': np.int32,'average_probability': np.float32})

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

    numberofSubs=0
    numberofIndels=0
    numberofDinucs=0

    if (SUBS in mutationtype_numberofmutations_numberofsamples_sampleslist_df['mutation_type'].unique()):
        numberofSubs=mutationtype_numberofmutations_numberofsamples_sampleslist_df.loc[ mutationtype_numberofmutations_numberofsamples_sampleslist_df['mutation_type']==SUBS,'number_of_mutations'].values[0]
    if (INDELS in mutationtype_numberofmutations_numberofsamples_sampleslist_df['mutation_type'].unique()):
        numberofIndels = mutationtype_numberofmutations_numberofsamples_sampleslist_df.loc[ mutationtype_numberofmutations_numberofsamples_sampleslist_df['mutation_type']==INDELS,'number_of_mutations'].values[0]
    if (DINUCS in mutationtype_numberofmutations_numberofsamples_sampleslist_df['mutation_type'].unique()):
        numberofDinucs = mutationtype_numberofmutations_numberofsamples_sampleslist_df.loc[ mutationtype_numberofmutations_numberofsamples_sampleslist_df['mutation_type']==DINUCS,'number_of_mutations'].values[0]

    # Tissue based
    for mutationType in mutation_types_contexts:
        if verbose: print('\tVerbose Worker pid %s Plot all mutations pooled %s\t%s' %(str(os.getpid()),str(mutationType),libraryFilenameMemo))
        plotAllMutationsPooledWithSimulations('Interval around variant (bp)',ylabel,None,outputDir,jobname,numberofSubs,numberofIndels,numberofDinucs,numberofSimulations,mutationType,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus,plot_mode)

    #Plot Signature Based
    #ncomms11383 Fig3b signature based average nucleosome occupancy figures
    if checkValidness(SIGNATUREBASED,outputDir,jobname,occupancy_type):
        if (not subsSignature_cutoff_numberofmutations_averageprobability_df.empty):
            #Subs Signatures
            if verbose: print('\tVerbose Worker pid %s Plot signature based SBS96 %s' % (str(os.getpid()),libraryFilenameMemo))
            plotSignatureBasedFigures(SBS96,subsSignature_cutoff_numberofmutations_averageprobability_df,sample2SubsSignature2NumberofMutationsDict,outputDir,jobname,numberofSimulations,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus,verbose,plot_mode)
        if (not dinucsSignature_cutoff_numberofmutations_averageprobability_df.empty):
            # Dinucs Signatures
            if verbose: print('\tVerbose Worker pid %s Plot signature based DBS %s' % (str(os.getpid()),libraryFilenameMemo))
            plotSignatureBasedFigures(DBS,dinucsSignature_cutoff_numberofmutations_averageprobability_df,sample2DinucsSignature2NumberofMutationsDict,outputDir,jobname,numberofSimulations,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus,verbose,plot_mode)
        if (not indelsSignature_cutoff_numberofmutations_averageprobability_df.empty):
            # Indels Signatures
            if verbose: print('\tVerbose Worker pid %s Plot signature based ID %s' % (str(os.getpid()), libraryFilenameMemo))
            plotSignatureBasedFigures(ID, indelsSignature_cutoff_numberofmutations_averageprobability_df,sample2IndelsSignature2NumberofMutationsDict, outputDir, jobname, numberofSimulations,libraryFilename, libraryFilenameMemo, occupancy_type, plusOrMinus, verbose, plot_mode)

    if sample_based:
        #ALL SAMPLES IN ONE
        #Plot "all samples pooled" and "sample based" signature based in one figure
        if (not subsSignature_cutoff_numberofmutations_averageprobability_df.empty):
            plotAllSamplesPooledAndSampleBasedSignaturesFiguresInOneFigure(subsSignature_cutoff_numberofmutations_averageprobability_df,sample2SubsSignature2NumberofMutationsDict,outputDir,jobname,'royalblue','Interval around single point mutation (bp)',ylabel,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus)
        if (not dinucsSignature_cutoff_numberofmutations_averageprobability_df.empty):
            plotAllSamplesPooledAndSampleBasedSignaturesFiguresInOneFigure(dinucsSignature_cutoff_numberofmutations_averageprobability_df,sample2DinucsSignature2NumberofMutationsDict,outputDir,jobname,'crimson','Interval around indel (bp)',ylabel,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus)
        if (not indelsSignature_cutoff_numberofmutations_averageprobability_df.empty):
            plotAllSamplesPooledAndSampleBasedSignaturesFiguresInOneFigure(indelsSignature_cutoff_numberofmutations_averageprobability_df,sample2IndelsSignature2NumberofMutationsDict,outputDir,jobname,'darkgreen','Interval around indel (bp)',ylabel,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus)

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
            for mutationType in mutation_types_contexts:
                plotAllMutationsPooledWithSimulations('Interval around variant (bp)',ylabel,sample,outputDir,jobname,numberofSubs, numberofIndels,numberofDinucs,numberofSimulations,mutationType,libraryFilename,libraryFilenameMemo,occupancy_type,plusOrMinus)
