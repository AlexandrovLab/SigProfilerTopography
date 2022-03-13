# !/usr/bin/env python3

# Author: burcakotlu

# Contact: burcakotlu@eng.ucsd.edu

import os
import matplotlib
import multiprocessing
import pandas as pd
import numpy as np

BACKEND = 'Agg'
if matplotlib.get_backend().lower() != BACKEND.lower():
    # If backend is not set properly a call to describe will hang
    matplotlib.use(BACKEND)

from matplotlib import pyplot as plt
#To plot for bi range of x axis such as 250M
plt.rcParams['agg.path.chunksize'] = 10000

from SigProfilerTopography.source.commons.TopographyCommons import  getChromSizesDict
from SigProfilerTopography.source.commons.TopographyCommons import  ONE_DIRECTORY_UP
from SigProfilerTopography.source.commons.TopographyCommons import  LIB
from SigProfilerTopography.source.commons.TopographyCommons import  NUCLEOSOME
from SigProfilerTopography.source.commons.TopographyCommons import  CHRBASED

from SigProfilerTopography.source.commons.TopographyCommons import  CHROM
from SigProfilerTopography.source.commons.TopographyCommons import  START
from SigProfilerTopography.source.commons.TopographyCommons import  END
from SigProfilerTopography.source.commons.TopographyCommons import  SIGNAL

current_abs_path = os.path.dirname(os.path.realpath(__file__))

######################################################################
#Used by plotting
def readChromBasedNucleosomeDF(chrLong,nucleosomeFilename):
    chrBasedNucleosmeFilename = '%s_%s' %(chrLong,nucleosomeFilename)
    chrBasedNucleosmeFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME, CHRBASED, chrBasedNucleosmeFilename)

    if (os.path.exists(chrBasedNucleosmeFilePath)):
        column_names = [CHROM, START, END, SIGNAL]
        chrbased_nucleosome_df = pd.read_csv(chrBasedNucleosmeFilePath, sep='\t', header=None, comment='#',names=column_names, dtype={CHROM: 'category', START: np.int32, END: np.int32, SIGNAL: np.float32})
        return chrbased_nucleosome_df
    else:
        return None
######################################################################


######################################################################
# For plotting
# main function parallel
# parallel it does not end for big chromosomes
def plotChrBasedNucleosomeOccupancyFigures(genome,nucleosomeFilename):
    #read chromnames for this nucleosome data

    nucleosomeFilenameWoExtension = nucleosomeFilename[0:-4]

    chromSizesDict = getChromSizesDict(genome)
    chromNamesList = list(chromSizesDict.keys())

    #Start the pool
    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)

    poolInputList = []
    for chrLong in chromNamesList:
        chromBasedNucleosomeDF = readChromBasedNucleosomeDF(chrLong,nucleosomeFilename)
        if chromBasedNucleosomeDF is not None:
            inputList = []
            inputList.append(chrLong)
            inputList.append(nucleosomeFilenameWoExtension)
            poolInputList.append(inputList)

    pool.map(plotNucleosomeOccupancySignalCountFiguresInParallel,poolInputList)

    ################################
    pool.close()
    pool.join()
    ################################

######################################################################

######################################################################
#main function sequential
def plotChrBasedNucleosomeOccupancyFiguresSequential(genome,nucleosomeFilename):
    #read chromnames for this nucleosome data
    nucleosomeFilenameWoExtension = nucleosomeFilename[0:-4]

    chromSizesDict = getChromSizesDict(genome)
    chromNamesList = list(chromSizesDict.keys())

    for chrLong in chromNamesList:
        #Actually no need for such a check
        chromBasedNucleosomeDF = readChromBasedNucleosomeDF(chrLong,nucleosomeFilename)
        if chromBasedNucleosomeDF is not None:
            inputList = []
            inputList.append(chrLong)
            inputList.append(nucleosomeFilenameWoExtension)
            plotNucleosomeOccupancySignalCountFiguresInParallel(inputList)
######################################################################



######################################################################
#main function
def plotChrBasedNucleosomeOccupancyFiguresFromText(genome,nucleosomeFilename):
    #read chromnames for this nucleosome data

    nucleosomeFilenameWoExtension = nucleosomeFilename[0:-4]

    chromSizesDict = getChromSizesDict(genome)
    chromNamesList = list(chromSizesDict.keys())

    #Start the pool
    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)

    poolInputList = []
    for chrLong in chromNamesList:
        chromBasedNucleosomeDF = readChromBasedNucleosomeDF(chrLong,nucleosomeFilename)
        if chromBasedNucleosomeDF is not None:
            inputList = []
            inputList.append(chrLong)
            inputList.append(nucleosomeFilenameWoExtension)
            poolInputList.append(inputList)

    pool.map(plotNucleosomeOccupancySignalCountFiguresInParallelFromTextFiles,poolInputList)

    ################################
    pool.close()
    pool.join()
    ################################

######################################################################


#####################################################################
# For plotting
def plotNucleosomeOccupancySignalCountFiguresInParallel(inputList):
    chrLong = inputList[0]
    nucleosomeFilenameWoExtension = inputList[1]

    ##############################################################
    signalArrayFilename = '%s_signal_%s.npy' %(chrLong,nucleosomeFilenameWoExtension)
    countArrayFilename = '%s_count_%s.npy' % (chrLong, nucleosomeFilenameWoExtension)

    chrBasedSignalNucleosmeFile = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,NUCLEOSOME,CHRBASED,signalArrayFilename)
    chrBasedCountNucleosmeFile = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,NUCLEOSOME,CHRBASED,countArrayFilename)

    #################################################################################################################
    if (os.path.exists(chrBasedSignalNucleosmeFile) and os.path.exists(chrBasedCountNucleosmeFile)):

        signal_array_npy = np.load(chrBasedSignalNucleosmeFile)
        count_array_npy = np.load(chrBasedCountNucleosmeFile)

        fig = plt.figure(figsize=(30, 10), facecolor=None)
        plt.style.use('ggplot')

        figureFilename = '%s_NucleosomeOccupancy_Signal_Count_from_npy_scatter_allChr.png' %(chrLong)
        figureFilepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB,NUCLEOSOME,CHRBASED,figureFilename)

        # This code makes the background white.
        ax = plt.gca()
        ax.set_facecolor('white')

        # This code puts the edge line
        for edge_i in ['left', 'bottom', 'right', 'top']:
            ax.spines[edge_i].set_edgecolor("black")
            ax.spines[edge_i].set_linewidth(3)

        #All Chrom
        chromSize = signal_array_npy.size
        x = np.arange(0, chromSize, 1)
        plt.scatter(x, signal_array_npy, color='black', label='Signal Array', zorder=1)
        plt.scatter(x, count_array_npy, color='red', label='Count Array',zorder=1)

        # #Small Portion of Chrom
        # if chrLong == 'chrM':
        #     chromSize = signal_array_npy.size
        #     x = np.arange(0, chromSize, 1)
        #     # signalPlot = plt.plot(x, signal_array_npy, 'black', label='Signal Array', linewidth=1,zorder=1)
        #     # countPlot = plt.plot(x, count_array_npy, 'red', label='Count Array', linewidth=1,zorder=1)
        #     signalPlot = plt.scatter(x, signal_array_npy, color='black', label='Signal Array',zorder=1)
        #     countPlot = plt.scatter(x, count_array_npy, color ='red', label='Count Array',zorder=1)
        # else:
        #     x=np.arange(41100000,41200000,1)
        #     # signalPlot = plt.plot(x, signal_array_npy[41100000:41200000], 'black', label='Signal Array', linewidth=1,zorder=1)
        #     # countPlot = plt.plot(x, count_array_npy[41100000:41200000], 'red', label='Count Array', linewidth=1,zorder=1)
        #     signalPlot = plt.scatter(x, signal_array_npy[41100000:41200000], color='black', label='Signal Array',zorder=1)
        #     countPlot = plt.scatter(x, count_array_npy[41100000:41200000], color='red', label='Count Array',zorder=1)

        print('chr %s' %chrLong)
        print('x shape and size')
        print(x.size)
        print(x.shape)

        print('signal array shape and size')
        print(signal_array_npy.size)
        print(signal_array_npy.shape)

        print('count array shape and size')
        print(count_array_npy.size)
        print(count_array_npy.shape)

        # listofLegends = []
        # listofLegends.append(signalPlot[0])
        # listofLegends.append(countPlot[0])
        # plt.legend(loc='upper left', handles=listofLegends, prop={'size': 24}, shadow=False, edgecolor='white',facecolor='white')

        # Put vertical line at x=0
        # plt.axvline(x=0, color='gray', linestyle='--')

        #Put ylim
        plt.ylim((0,10))

        plt.title('For %s' %(chrLong), fontsize=40, fontweight='bold')
        # plt.show()
        fig.savefig(figureFilepath)
        plt.close(fig)
    #################################################################################################################

#####################################################################


#####################################################################
#For plotting
def plotNucleosomeOccupancySignalCountFiguresInParallelFromTextFiles(inputList):
    chrLong = inputList[0]
    nucleosomeFilenameWoExtension = inputList[1]

    ##############################################################
    signalArrayFilename = '%s_signal_%s.txt' %(chrLong,nucleosomeFilenameWoExtension)
    countArrayFilename = '%s_count_%s.txt' % (chrLong, nucleosomeFilenameWoExtension)

    chrBasedSignalNucleosmeFile = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,NUCLEOSOME,CHRBASED,signalArrayFilename)
    chrBasedCountNucleosmeFile = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,NUCLEOSOME,CHRBASED,countArrayFilename)

    #################################################################################################################
    if (os.path.exists(chrBasedSignalNucleosmeFile) and os.path.exists(chrBasedCountNucleosmeFile)):

        signal_array_txt = np.loadtxt(chrBasedSignalNucleosmeFile)
        count_array_txt = np.loadtxt(chrBasedCountNucleosmeFile)

        fig = plt.figure(figsize=(30, 10), facecolor=None)
        plt.style.use('ggplot')

        figureFilename = '%s_NucleosomeOccupancy_Signal_Count_from_text.png' %(chrLong)
        figureFilepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB,NUCLEOSOME,CHRBASED,figureFilename)

        # This code makes the background white.
        ax = plt.gca()
        ax.set_facecolor('white')

        # This code puts the edge line
        for edge_i in ['left', 'bottom', 'right', 'top']:
            ax.spines[edge_i].set_edgecolor("black")
            ax.spines[edge_i].set_linewidth(3)


        #All Chrom
        chromSize = signal_array_txt.size
        x = np.arange(0, chromSize, 1)
        signalPlot = plt.plot(x, signal_array_txt, 'black', label='Signal Array', marker='.', zorder=1)
        countPlot = plt.plot(x, count_array_txt, 'red', label='Count Array', marker='.',zorder=1)

        # # Small Portion of Chrom
        # if chrLong == 'chrM':
        #     chromSize = signal_array_txt.size
        #     x = np.arange(0, chromSize, 1)
        #     signalPlot = plt.plot(x, signal_array_txt, 'black', label='Signal Array', linewidth=1,zorder=1)
        #     countPlot = plt.plot(x, count_array_txt, 'red', label='Count Array', linewidth=1,zorder=1)
        # else:
        #     x=np.arange(41100000,41200000,1)
        #     signalPlot = plt.plot(x, signal_array_txt[41100000:41200000], 'black', label='Signal Array', linewidth=1,zorder=1)
        #     countPlot = plt.plot(x, count_array_txt[41100000:41200000], 'red', label='Count Array', linewidth=1,zorder=1)

        print('chr %s' %chrLong)
        print('signal array shape and size')
        print(signal_array_txt.size)
        print(signal_array_txt.shape)
        print('count array shape and size')
        print(count_array_txt.size)
        print(count_array_txt.shape)

        listofLegends = []
        listofLegends.append(signalPlot[0])
        listofLegends.append(countPlot[0])
        plt.legend(loc='upper left', handles=listofLegends, prop={'size': 24}, shadow=False, edgecolor='white',facecolor='white')

        # Put vertical line at x=0
        # plt.axvline(x=0, color='gray', linestyle='--')

        #Put ylim
        plt.ylim((0,10))

        plt.title('For %s' %(chrLong), fontsize=40, fontweight='bold')
        # plt.show()
        fig.savefig(figureFilepath)
        plt.close(fig)
#####################################################################
