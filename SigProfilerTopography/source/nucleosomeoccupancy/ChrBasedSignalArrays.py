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

from SigProfilerTopography.source.commons.TopographyCommons import generateIntervalVersion

from SigProfilerTopography.source.commons.TopographyCommons import CHROM 
from SigProfilerTopography.source.commons.TopographyCommons import START 
from SigProfilerTopography.source.commons.TopographyCommons import END
from SigProfilerTopography.source.commons.TopographyCommons import SIGNAL

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

from SigProfilerTopography.source.commons.TopographyCommons import CHROM
from SigProfilerTopography.source.commons.TopographyCommons import START
from SigProfilerTopography.source.commons.TopographyCommons import END
from SigProfilerTopography.source.commons.TopographyCommons import SIGNAL

from SigProfilerTopography.source.commons.TopographyCommons import getChromSizesDict
from SigProfilerTopography.source.commons.TopographyCommons import memory_usage
from SigProfilerTopography.source.commons.TopographyCommons import readFileInBEDFormat

######################################################################
def updateSignalArray(row,signalArray):
    # signalArray[row[START]:row[END]] += row[SIGNAL]
    signalArray[row[1]:row[2]] += row[3]

######################################################################

######################################################################
def updateSignalArraysForListComprehension(row,signalArrayDict):
    #row [chrom start end signal]
    signalArray = signalArrayDict[row[0]]
    signalArray[row[1]:row[2]] += row[3]
######################################################################

######################################################################
# This write chrBasedSignalArray is used right now.
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
    [updateSignalArray(row,signalArray) for row in chrBasedFileDF.values]

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
#Used by nuclesome occupancy read all file at once
#Filter w.r.t. a quantileValue
def readNucleosomeOccupancyData(quantileValue,nucleosomeFilename):

    column_names = [CHROM, START, END, SIGNAL]
    nucleosome_df = pd.read_table(nucleosomeFilename, sep="\t", header=None, comment='#', names=column_names, dtype={CHROM: np.category, START: np.int32, END: np.int32, SIGNAL: np.float32})

    print('After nucleosome occupancy is loaded into memory')
    print('Memory usage in %s MB' % memory_usage())

    #########################################################
    if (quantileValue < 1.0):
        # remove the outliers
        q = nucleosome_df[SIGNAL].quantile(quantileValue)
        print('q:%f' % q)
        print('before %d' % (nucleosome_df.shape[0]))
        nucleosome_df = nucleosome_df[nucleosome_df[SIGNAL] < q]

        max_signal=nucleosome_df[SIGNAL].max()
        min_signal=nucleosome_df[SIGNAL].min()

        print('After nucleosome_df is subset')
        print('Memory usage in %s MB' % memory_usage())
    #########################################################

    nucleosome_df_grouped = nucleosome_df.groupby(CHROM)

    print('After nucleosome occupancy grouped by')
    print('Memory usage in %s MB' % memory_usage())

    return nucleosome_df_grouped, max_signal, min_signal
######################################################################


######################################################################
#Read all data at once, filter w.r.t. a quantileValue and write chrom based signal npy files sequentially
#This function is used for GM12878 and K562 nucleosome occupancy chrom based signal npy files creation sequentially
#For in house use, also has parallel write version
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


##################################################################
#JAN 10, 2020
#If file is originally a  wig file use this function
#Read at once, write in parallel
#
#If original file is bigWig and converted from bigWig to wig and has lines starts with #bedGraph it means that it is now bedGraph format
#In fact this is bedgraph format
# There can be wig files prepared from bedGraph which includes many lines start with #bedGraph
#This code is intended for bigWig nucleosome occupancy files
# bedGraph section chr1:10111-60880
# chr1    10111   10112   4.9
# chr1    10112   10114   5
# chr1    10114   10116   5.1
# bedGraph section chrX:155161349-155260605
# chrX    155161349       155161419       0
# chrX    155161470       155161530       0
# chrX    155161597       155161634       0
#Please note that end is exclusive
def readWig_with_fixedStep_variableStep_writeChrBasedSignalArrays(outputDir, jobname, genome,wig_file_path,occupancy_type,quantileValue,remove_outliers):
    wig_file_name=os.path.splitext(os.path.basename(wig_file_path))[0]

    chromSizesDict = getChromSizesDict(genome)
    numofProcesses = multiprocessing.cpu_count()

    #To reduce memory footprint
    # Delete old chr based epigenomics files
    if occupancy_type==EPIGENOMICSOCCUPANCY:
        chrBasedNpyFilesPath=os.path.join(outputDir,jobname,DATA,occupancy_type,LIB,CHRBASED)
        deleteChrBasedNpyFiles(chrBasedNpyFilesPath)

    if os.path.exists(wig_file_path):
        # Read the wavelet signal
        wig_unprocessed_df = pd.read_table(wig_file_path, sep="\t", comment='#', header=None)

        # Process the wavelet signal, convert into interval version
        # Add column names
        wigfile_interval_version_df = generateIntervalVersion(wig_unprocessed_df)
        max_signal = wigfile_interval_version_df[SIGNAL].max()
        min_signal = wigfile_interval_version_df[SIGNAL].min()

        # Outlier elimination starts
        if ((remove_outliers==True) and (quantileValue < 1.0)):
            # remove the outliers
            q = wigfile_interval_version_df[SIGNAL].quantile(quantileValue)
            print('Signal greater than %f is removed for quantile: %f' % (q,quantileValue))
            print('before outlier removal number of rows: %d' % (wigfile_interval_version_df.shape[0]))
            file_df = wigfile_interval_version_df[wigfile_interval_version_df[SIGNAL] < q]
            print('after outlier removal number of rows: %d' % (file_df.shape[0]))
        # Outlier elimination ends

        file_df_grouped = file_df.groupby(CHROM)
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
            inputList.append(wig_file_name)
            inputList.append(occupancy_type)
            inputList.append(max_signal)
            inputList.append(min_signal)
            poolInputList.append(inputList)

        pool.map(writeChrBasedOccupancySignalArray, poolInputList)

        ################################
        pool.close()
        pool.join()
        ################################

##################################################################

######################################################################
#Read in chunks, at the end have all signal arrays in memory and write sequentially
def readWig_write_derived_from_bedgraph(outputDir, jobname, genome, library_file_with_path, occupancy_type):
    wig_file_name_wo_extension = os.path.splitext(os.path.basename(library_file_with_path))[0]
    chromSizesDict = getChromSizesDict(genome)

    #To reduce memory footprint
    # Delete old chr based epigenomics files
    if occupancy_type==EPIGENOMICSOCCUPANCY:
        chrBasedNpyFilesPath=os.path.join(outputDir,jobname,DATA,occupancy_type,LIB,CHRBASED)
        deleteChrBasedNpyFiles(chrBasedNpyFilesPath)

    if os.path.exists(library_file_with_path):
        column_names = [CHROM, START, END, SIGNAL]

        #chromSize
        signalArrayDict={}
        max_signal=np.finfo(np.float32).min
        min_signal=np.finfo(np.float32).max

        for chunk_df in  pd.read_table(library_file_with_path,chunksize=5000000, sep="\t", header=None, comment='#', names=column_names, dtype={CHROM: str, START: np.int32, END: np.int32, SIGNAL: np.float32}):
            print(chunk_df[CHROM].unique())

            chrom_list=chunk_df[CHROM].unique()
            if (chunk_df[SIGNAL].max()>max_signal):
                max_signal=chunk_df[SIGNAL].max()
            if (chunk_df[SIGNAL].min()<min_signal):
                min_signal=chunk_df[SIGNAL].min()

            for chrom in chrom_list:
                if chrom not in signalArrayDict:
                    signalArray=np.zeros(chromSizesDict[chrom], dtype=np.float32)
                    signalArrayDict[chrom]=signalArray

            # Use list comprehension
            [updateSignalArraysForListComprehension(row, signalArrayDict) for row in chunk_df.values]
            print('#################')

        print('\n#####################################################')
        print('min_signal:%f max_signal:%f' %(min_signal,max_signal))
        print('#####################################################\n')

        if (occupancy_type == EPIGENOMICSOCCUPANCY):
            os.makedirs(os.path.join(outputDir, jobname, DATA, occupancy_type, LIB, CHRBASED), exist_ok=True)
        elif (occupancy_type == NUCLEOSOMEOCCUPANCY):
            os.makedirs(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME, CHRBASED),exist_ok=True)

        for key in signalArrayDict.keys():
            signalArray=signalArrayDict[key]
            print('%s np.sum(signalArray,dtype=np.float32):%f' %(key,np.sum(signalArray,dtype=np.float32)))
            signalArrayFilename = '%s_signal_%s' %(key,wig_file_name_wo_extension)
            if (occupancy_type == EPIGENOMICSOCCUPANCY):
                chrBasedSignalFile = os.path.join(outputDir, jobname, DATA, occupancy_type, LIB, CHRBASED, signalArrayFilename)
            elif (occupancy_type == NUCLEOSOMEOCCUPANCY):
                chrBasedSignalFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME,CHRBASED, signalArrayFilename)
            np.save(chrBasedSignalFile,signalArray)
######################################################################


######################################################################
#Dec 2, 2019
#This function is used when user provides its own bed/np/narrowpeak files for occupancy analysis
#Read at once writes in parallel
def readBEDandWriteChromBasedSignalArrays(outputDir, jobname, genome,BEDFileWithPath,occupancy_type,quantileValue,remove_outliers):
    chromSizesDict = getChromSizesDict(genome)
    numofProcesses = multiprocessing.cpu_count()

    #To reduce memory footprint
    # Delete old chr based epigenomics files
    if occupancy_type==EPIGENOMICSOCCUPANCY:
        chrBasedNpyFilesPath=os.path.join(outputDir,jobname,DATA,occupancy_type,LIB,CHRBASED)
        deleteChrBasedNpyFiles(chrBasedNpyFilesPath)

    if os.path.exists(BEDFileWithPath):
        bedfilename, bedfile_extention= os.path.splitext(os.path.basename(BEDFileWithPath))

        if (bedfile_extention.lower()=='.bed' or bedfile_extention.lower()=='.np' or bedfile_extention.lower()=='.narrowpeak'):
            discard_signal=False
            file_df, max_signal, min_signal=readFileInBEDFormat(BEDFileWithPath,discard_signal)

            #Outlier elimination starts
            if ((remove_outliers==True) and (quantileValue < 1.0)):
                # remove the outliers
                q = file_df[SIGNAL].quantile(quantileValue)
                print('Signal greater than %f is removed for quantile: %f' % (q,quantileValue))
                print('before outlier removal number of rows: %d' % (file_df.shape[0]))
                file_df = file_df[file_df[SIGNAL] < q]
                print('after outlier removal number of rows: %d' % (file_df.shape[0]))
            #Outlier elimination ends

            file_df_grouped = file_df.groupby(CHROM)
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
#Read all data at once, filter w.r.t. a quantileValue and write chrom based in parallel
#This function is used for GM12878 and K562 nucleosome occupancy chrom based signal npy files creation in parallel
#For in house use, also has sequential write version
def readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArraysInParallel(genome, quantileValue, nucleosomeFilename):
    chromSizesDict = getChromSizesDict(genome)

    # Start the pool
    numofProcesses = multiprocessing.cpu_count()
    print('Number of processors:%d' %(numofProcesses))

    print('Before nucleosome occupancy is loaded into memory')
    print('Memory usage in %s MB' %memory_usage())

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