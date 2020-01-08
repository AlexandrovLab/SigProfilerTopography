# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time and strand bias.
# Copyright (C) 2018 Burcak Otlu

###################################################################################
################ This python code read the nucleosome occupancy filename ##########
######################## Remove the outliers or not ###############################
###### Write the chromosome based average nucleosome occupancy signal arrays ######
###################################################################################

#This code must be working for bed or wig files.
#Converting bed and wig files into chr based dataframes and then chr based arrays
#Please notice that unnesting requires too much memory because of that jobs sleep and then become zombie

import os
import matplotlib
import sys
from sys import getsizeof
import pandas as pd
import numpy as np
import multiprocessing
import shutil


BACKEND = 'Agg'
if matplotlib.get_backend().lower() != BACKEND.lower():
    # If backend is not set properly a call to describe will hang
    matplotlib.use(BACKEND)

from matplotlib import pyplot as plt
#To plot for bi range of x axis such as 250M
plt.rcParams['agg.path.chunksize'] = 10000

#############################################################
current_abs_path = os.path.dirname(os.path.realpath(__file__))
#############################################################

commonsPath = os.path.join(current_abs_path, '..','commons')
sys.path.append(commonsPath)

from SigProfilerTopography.source.commons.TopographyCommons import EPIGENOMICSOCCUPANCY
from SigProfilerTopography.source.commons.TopographyCommons import NUCLEOSOMEOCCUPANCY

from SigProfilerTopography.source.commons.TopographyCommons import ONE_DIRECTORY_UP
from SigProfilerTopography.source.commons.TopographyCommons import LIB
from SigProfilerTopography.source.commons.TopographyCommons import DATA
from SigProfilerTopography.source.commons.TopographyCommons import EPIGENOMICS
from SigProfilerTopography.source.commons.TopographyCommons import NUCLEOSOME
from SigProfilerTopography.source.commons.TopographyCommons import CHRBASED

from SigProfilerTopography.source.commons.TopographyCommons import chrom
from SigProfilerTopography.source.commons.TopographyCommons import start
from SigProfilerTopography.source.commons.TopographyCommons import end
from SigProfilerTopography.source.commons.TopographyCommons import signal

from SigProfilerTopography.source.commons.TopographyCommons import getChromSizesDict
from SigProfilerTopography.source.commons.TopographyCommons import memory_usage
from SigProfilerTopography.source.commons.TopographyCommons import readFileInBEDFormat

#########################################################################
def unnesting(df, explode):
    idx=df.index.repeat(df[explode[0]].str.len())
    df1=pd.concat([pd.DataFrame({x:np.concatenate(df[x].values)} )for x in explode],axis=1)
    df1.index=idx
    return df1.join(df.drop(explode,1),how='left')
#########################################################################


######################################################################
def updateSignalArrays(row,signalArray):
    signalArray[row[start]:row[end]] += row[signal]
######################################################################

######################################################################
def updateSignalArraysForListComprehension(row,signalArray):
    #row [chrom start end signal]
    signalArray[row[1]:row[2]] += row[3]
######################################################################

######################################################################
# This is used right now.
def writeChrBasedOccupancySignalArray(inputList):
    outputDir = inputList[0]
    jobname = inputList[1]
    chrLong = inputList[2]
    chromSize = inputList[3]
    chrBasedFileDF = inputList[4]
    filename = inputList[5]
    occupancy_type= inputList[6]
    max_signal=inputList[7]
    min_signal = inputList[8]

    if (occupancy_type==EPIGENOMICSOCCUPANCY):
        os.makedirs(os.path.join(outputDir,jobname,DATA,occupancy_type,LIB,CHRBASED), exist_ok=True)
    elif (occupancy_type==NUCLEOSOMEOCCUPANCY):
        os.makedirs(os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,NUCLEOSOME,CHRBASED), exist_ok=True)

    filenameWoExtension = os.path.splitext(os.path.basename(filename))[0]

    if (np.finfo(np.float16).min<=min_signal) and (max_signal<=np.finfo(np.float16).max):
        signalArray = np.zeros(chromSize, dtype=np.float16)
    else:
        signalArray = np.zeros(chromSize,dtype=np.float32)

    # Using apply
    # chrBasedFileDF.apply(updateSignalArrays, signalArray=signalArray, axis=1)

    # Use list comprehension
    [updateSignalArraysForListComprehension(row,signalArray) for row in chrBasedFileDF.values]

    #############################  Save as npy starts ################################
    signalArrayFilename = '%s_signal_%s' %(chrLong,filenameWoExtension)
    # signalArrayFilenameText = '%s_signal_%s.txt' %(chrLong,filenameWoExtension)
    if occupancy_type==EPIGENOMICSOCCUPANCY:
        #When we submit multiple jobs then they try to write on the same location on disk and it causes error.
        #Therefore save location is changed.
        # chrBasedSignalFile = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,EPIGENOMICS,CHRBASED,signalArrayFilename)
        chrBasedSignalFile = os.path.join(outputDir,jobname,DATA,occupancy_type,LIB,CHRBASED,signalArrayFilename)
        np.save(chrBasedSignalFile,signalArray)
    if occupancy_type==NUCLEOSOMEOCCUPANCY:
        chrBasedSignalFile = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,NUCLEOSOME,CHRBASED,signalArrayFilename)
        np.save(chrBasedSignalFile, signalArray)
    #############################  Save as npy ends ##################################

######################################################################




######################################################################
#Used by nuclesome occupancy read all file write chrom based npy files.
def readNucleosomeOccupancyData(quantileValue,nucleosomeFilename):

    column_names = [chrom, start, end, signal]
    nucleosome_df = pd.read_table(nucleosomeFilename, sep="\t", header=None, comment='#', names=column_names, dtype={chrom: np.category, start: np.int32, end: np.int32, signal: np.float32})

    print('After nucleosome occupancy is loaded into memory')
    print('Memory usage in %s MB' % memory_usage())

    #########################################################
    if (quantileValue < 1.0):
        # remove the outliers
        q = nucleosome_df[signal].quantile(quantileValue)
        print('q:%f' % q)
        print('before %d' % (nucleosome_df.shape[0]))
        nucleosome_df = nucleosome_df[nucleosome_df[signal] < q]

        max_signal=nucleosome_df[signal].max()
        min_signal=nucleosome_df[signal].min()

        print('After nucleosome_df is subset')
        print('Memory usage in %s MB' % memory_usage())
    #########################################################

    nucleosome_df_grouped = nucleosome_df.groupby(chrom)

    print('After nucleosome occupancy grouped by')
    print('Memory usage in %s MB' % memory_usage())

    return nucleosome_df_grouped, max_signal, min_signal
######################################################################


######################################################################
def readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArraysSequentially(genome, quantileValue, nucleosomeFilename):
    chromSizesDict = getChromSizesDict(genome)

    print('Before nucleosome occupancy is loaded into memory')
    print('Memory usage in %s MB' %memory_usage())

    if os.path.exists(nucleosomeFilename):

        nucleosome_df_grouped, max_signal, min_signal = readNucleosomeOccupancyData(quantileValue,nucleosomeFilename)

        print('After nucleosome occupancy grouped by')
        print('Memory usage in %s MB' % memory_usage())

        for chrLong, chromBasedNucleosomeDF in nucleosome_df_grouped:
            print('Debug June 13, 2019 For %s write nucleosome signal and count array' %(chrLong))
            chromSize = chromSizesDict[chrLong]
            inputList = []
            inputList.append(chrLong)
            inputList.append(chromSize)
            inputList.append(chromBasedNucleosomeDF)
            inputList.append(nucleosomeFilename)
            inputList.append(NUCLEOSOMEOCCUPANCY)
            inputList.append(max_signal)
            inputList.append(min_signal)
            writeChrBasedOccupancySignalArray(inputList)

        print('After all chr based files are written')
        print('Memory usage in %s MB' % memory_usage())

######################################################################


#######################################################
def deleteChrBasedNpyFiles(chrBasedNpyFilesPath):
    #############################################
    # Delete the chrBasedNpyFilesPath if exists

    ################################################
    if (os.path.exists(chrBasedNpyFilesPath)):
        try:
            shutil.rmtree(chrBasedNpyFilesPath,ignore_errors=True)
        except OSError as e:
            print('Error: %s - %s.' % (e.filename, e.strerror))
    ################################################
#######################################################



######################################################################
#Dec 2, 2019
def readBEDandWriteChromBasedSignalArrays(outputDir, jobname, genome,BEDFileWithPath,occupancy_type):
    chromSizesDict = getChromSizesDict(genome)
    numofProcesses = multiprocessing.cpu_count()

    #To reduce memory footprint
    # Delete old chr based epigenomics files
    if occupancy_type==EPIGENOMICSOCCUPANCY:
        chrBasedNpyFilesPath=os.path.join(outputDir,jobname,DATA,occupancy_type,LIB,CHRBASED)
        deleteChrBasedNpyFiles(chrBasedNpyFilesPath)

    if os.path.exists(BEDFileWithPath):
        bedFile=os.path.basename(BEDFileWithPath)

        bedfilename, bedfile_extention= os.path.splitext(bedFile)

        if (bedfile_extention.lower()=='.bed' or bedfile_extention.lower()=='.np' or bedfile_extention.lower()=='.narrowpeak'):
            file_df, max_signal, min_signal=readFileInBEDFormat(BEDFileWithPath)
            file_df_grouped = file_df.groupby(chrom)
            pool = multiprocessing.Pool(numofProcesses)

            poolInputList = []

            for chrLong, chromBasedFileDF in file_df_grouped:
                chromSize = chromSizesDict[chrLong]
                inputList = []
                inputList.append(outputDir)
                inputList.append(jobname)
                inputList.append(chrLong)
                inputList.append(chromSize)
                inputList.append(chromBasedFileDF)
                inputList.append(bedfilename)
                inputList.append(occupancy_type)
                inputList.append(max_signal)
                inputList.append(min_signal)
                poolInputList.append(inputList)

            pool.map(writeChrBasedOccupancySignalArray, poolInputList)

            ################################
            pool.close()
            pool.join()
            ################################

######################################################################



######################################################################
def readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArraysInParallel(genome, quantileValue, nucleosomeFilename):
    chromSizesDict = getChromSizesDict(genome)

    # Start the pool
    numofProcesses = multiprocessing.cpu_count()
    print('Number of processors:%d' %(numofProcesses))

    print('Before nucleosome occupancy is loaded into memory')
    print('Memory usage in %s MB' %memory_usage())

    column_names = [chrom, start, end, signal]
    if os.path.exists(nucleosomeFilename):

        nucleosome_df_grouped, max_signal, min_signal = readNucleosomeOccupancyData(quantileValue, nucleosomeFilename)
        pool = multiprocessing.Pool(numofProcesses)

        print('After pool is initialized')
        print('Memory usage in %s MB' % memory_usage())

        poolInputList = []

        for chrLong, chromBasedNucleosomeDF in nucleosome_df_grouped:
            print('for %s write nucleosome signal and count array' %(chrLong))
            chromSize = chromSizesDict[chrLong]
            inputList = []
            inputList.append(chrLong)
            inputList.append(chromSize)
            inputList.append(chromBasedNucleosomeDF)
            inputList.append(nucleosomeFilename)
            inputList.append(NUCLEOSOMEOCCUPANCY)
            inputList.append(max_signal)
            inputList.append(min_signal)
            poolInputList.append(inputList)

        pool.map(writeChrBasedOccupancySignalArray, poolInputList)

        ################################
        pool.close()
        pool.join()
        ################################

        print('After pool is closed and joined')
        print('Memory usage in %s MB' % memory_usage())

######################################################################





