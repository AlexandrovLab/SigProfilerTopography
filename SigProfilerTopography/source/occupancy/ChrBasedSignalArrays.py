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

from SigProfilerTopography.source.commons.TopographyCommons import NAME
from SigProfilerTopography.source.commons.TopographyCommons import SCORE
from SigProfilerTopography.source.commons.TopographyCommons import STRAND

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

######################################################################
def updateSignalArray(row,signalArray):
    # signalArray[row[START]:row[END]] += row[SIGNAL]
    signalArray[row[1]:row[2]] += row[3]
######################################################################

######################################################################
def updateSignalArraysForListComprehension(row,signalArray):
    #row [chrom start end signal]
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
    file_name_with_path = inputList[5]
    occupancy_type= inputList[6]
    max_signal=inputList[7]
    min_signal = inputList[8]

    if (occupancy_type==NUCLEOSOMEOCCUPANCY):
        os.makedirs(os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,NUCLEOSOME,CHRBASED), exist_ok=True)
    elif (occupancy_type==EPIGENOMICSOCCUPANCY):
        os.makedirs(os.path.join(outputDir,jobname,DATA,occupancy_type,LIB,CHRBASED), exist_ok=True)
    else:
        os.makedirs(os.path.join(outputDir,jobname,DATA,occupancy_type,LIB,CHRBASED), exist_ok=True)

    filenameWoExtension = os.path.splitext(os.path.basename(file_name_with_path))[0]

    # if (np.finfo(np.float16).min<=min_signal) and (max_signal<=np.finfo(np.float16).max):
    #     print("Signal array dtype is set as np.float16")
    #     signalArray = np.zeros(chromSize, dtype=np.float16)
    # else:
    #     print("Signal array dtype is set as np.float32")
    #     signalArray = np.zeros(chromSize,dtype=np.float32)

    #To avoid overflow
    # print("Signal array dtype is set as np.float32")
    signalArray = np.zeros(chromSize,dtype=np.float32)

    # Using apply
    # chrBasedFileDF.apply(updateSignalArrays, signalArray=signalArray, axis=1)

    # Use list comprehension
    [updateSignalArray(row,signalArray) for row in chrBasedFileDF.values]

    #############################  Save as npy starts ################################
    signalArrayFilename = '%s_signal_%s' %(chrLong,filenameWoExtension)
    # signalArrayFilenameText = '%s_signal_%s.txt' %(chrLong,filenameWoExtension)
    if occupancy_type==NUCLEOSOMEOCCUPANCY:
        chrBasedSignalFile = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,NUCLEOSOME,CHRBASED,signalArrayFilename)
        np.save(chrBasedSignalFile, signalArray)
    else:
        #When we submit multiple jobs then they try to write on the same location on disk and it causes error.
        #Therefore save location is changed.
        # chrBasedSignalFile = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,EPIGENOMICS,CHRBASED,signalArrayFilename)
        chrBasedSignalFile = os.path.join(outputDir,jobname,DATA,occupancy_type,LIB,CHRBASED,signalArrayFilename)
        np.save(chrBasedSignalFile,signalArray)
    #############################  Save as npy ends ##################################

######################################################################


######################################################################
#Used by nuclesome occupancy read all file at once
#Filter w.r.t. a quantileValue
def readNucleosomeOccupancyData(quantileValue,nucleosomeFilename):

    column_names = [CHROM, START, END, SIGNAL]
    nucleosome_df = pd.read_csv(nucleosomeFilename, sep='\t', header=None, comment='#', names=column_names, dtype={CHROM: 'category', START: np.int32, END: np.int32, SIGNAL: np.float32})

    max_signal = nucleosome_df[SIGNAL].max()
    min_signal = nucleosome_df[SIGNAL].min()
    mean_signal = nucleosome_df[SIGNAL].mean()
    std_signal = nucleosome_df[SIGNAL].std()
    print('\n#########################################')
    print('Before outlier elimination')
    print('Max Signal: %f' % max_signal)
    print('Min Signal: %f' % min_signal)
    print('Mean Signal: %f' % mean_signal)
    print('Std Signal: %f' % std_signal)
    print('Memory usage in %s MB' % memory_usage())

    #########################################################
    if (quantileValue < 1.0):
        # remove the outliers
        q = nucleosome_df[SIGNAL].quantile(quantileValue)
        print('\n#########################################')
        print('q:%f' % q)
        print('before outlier elimination number of rows: %d' % (nucleosome_df.shape[0]))
        nucleosome_df = nucleosome_df[nucleosome_df[SIGNAL] < q]
        print('after outlier elimination number of rows: %d' % (nucleosome_df.shape[0]))

        max_signal = nucleosome_df[SIGNAL].max()
        min_signal = nucleosome_df[SIGNAL].min()
        mean_signal = nucleosome_df[SIGNAL].mean()
        std_signal = nucleosome_df[SIGNAL].std()
        print('\n#########################################')
        print('After outlier elimination')
        print('Max Signal: %f' % max_signal)
        print('Min Signal: %f' % min_signal)
        print('Mean Signal: %f' % mean_signal)
        print('Std Signal: %f' % std_signal)
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
        print('max_signal:%d min_signal:%d' %(max_signal,min_signal))
        print('np.finfo(np.float16).min:%f np.finfo(np.float16).max:%f' % (np.finfo(np.float16).min,np.finfo(np.float16).max))

        print('After nucleosome occupancy grouped by')
        print('Memory usage in %s MB' % memory_usage())

        for chrLong, chromBasedNucleosomeDF in nucleosome_df_grouped:
            print('For %s write nucleosome signal and count array' %(chrLong))
            chromSize = chromSizesDict[chrLong]
            inputList = []
            inputList.append(None)
            inputList.append(None)
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
    chrBasedNpyFilesPath=os.path.join(outputDir,jobname,DATA,occupancy_type,LIB,CHRBASED)
    deleteChrBasedNpyFiles(chrBasedNpyFilesPath)

    if os.path.exists(wig_file_path):
        # Read the wavelet signal
        wig_unprocessed_df = pd.read_csv(wig_file_path, sep='\t', comment='#', header=None)

        # Process the wavelet signal, convert into interval version
        # Add column names
        wigfile_interval_version_df = generateIntervalVersion(wig_unprocessed_df)

        max_signal = wigfile_interval_version_df[SIGNAL].max()
        min_signal = wigfile_interval_version_df[SIGNAL].min()
        mean_signal = wigfile_interval_version_df[SIGNAL].mean()
        std_signal = wigfile_interval_version_df[SIGNAL].std()
        print('Max Signal: %f' % max_signal)
        print('Min Signal: %f' % min_signal)
        print('Mean Signal: %f' % mean_signal)
        print('Std Signal: %f' % std_signal)


        # Outlier elimination starts
        if ((remove_outliers==True) and (quantileValue < 1.0)):
            # remove the outliers
            q = wigfile_interval_version_df[SIGNAL].quantile(quantileValue)
            print('Signal greater than %f is removed for quantile: %f' % (q,quantileValue))
            print('before outlier removal number of rows: %d' % (wigfile_interval_version_df.shape[0]))
            file_df = wigfile_interval_version_df[wigfile_interval_version_df[SIGNAL] < q]
            print('after outlier removal number of rows: %d' % (file_df.shape[0]))
        else:
            file_df = wigfile_interval_version_df

        file_df_grouped = file_df.groupby(CHROM)
        pool = multiprocessing.Pool(numofProcesses)

        poolInputList = []

        for chrLong, chromBasedFileDF in file_df_grouped:
            if chrLong in chromSizesDict:
                chromSize = chromSizesDict[chrLong]
                inputList = []
                inputList.append(outputDir)
                inputList.append(jobname)
                inputList.append(chrLong)
                inputList.append(chromSize)
                inputList.append(chromBasedFileDF)
                inputList.append(wig_file_path)
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

##################################################################
def fill_signal_array_dict(chrname,chr_based_df_list,chromSizesDict,occupancy_type,verbose):

    if verbose: print('\tVerbose %s Worker pid %s Fill Signal Array Dict starts:  current_mem_usage %.2f (mb)' % (occupancy_type, str(os.getpid()), memory_usage()))
    process_signal_array_dict={}

    process_signal_array = np.zeros(chromSizesDict[chrname], dtype=np.float32)
    max_signal = np.finfo(np.float32).min
    min_signal = np.finfo(np.float32).max

    if verbose: print('\tVerbose %s Worker pid %s Fill Signal Array Dict --- len(chr_based_df_list):%d current_mem_usage %.2f (mb)' % (occupancy_type, str(os.getpid()), len(chr_based_df_list),memory_usage()))
    if verbose: print('\tVerbose %s Worker pid %s Fill Signal Array Dict --- chrname:%s current_mem_usage %.2f (mb)' % (occupancy_type, str(os.getpid()), chrname,memory_usage()))

    for chr_based_df in chr_based_df_list:
        if (chr_based_df[SIGNAL].max() > max_signal):
            max_signal = chr_based_df[SIGNAL].max()
        if (chr_based_df[SIGNAL].min() < min_signal):
            min_signal = chr_based_df[SIGNAL].min()
        # Use list comprehension
        [updateSignalArraysForListComprehension(row, process_signal_array) for row in chr_based_df.values]

    # Initialzie the list, you will return this list
    process_signal_array_dict['min_signal'] = min_signal
    process_signal_array_dict['max_signal'] = max_signal
    process_signal_array_dict[chrname]=process_signal_array

    if verbose: print('\tVerbose %s Worker pid %s Fill Signal Array Dict: min_signal: %f max_signal: %f current_mem_usage %.2f (mb)' % (occupancy_type, str(os.getpid()), process_signal_array_dict['min_signal'], process_signal_array_dict['max_signal'], memory_usage()))
    if verbose: print('\tVerbose %s Worker pid %s Fill Signal Array Dict: signal_array_dict.keys():%s current_mem_usage %.2f (mb)' % (occupancy_type, str(os.getpid()), process_signal_array_dict.keys(), memory_usage()))
    if verbose: print('\tVerbose %s Worker pid %s Fill Signal Array Dict ends: current_mem_usage %.2f (mb)\n' % (occupancy_type, str(os.getpid()), memory_usage()))

    return process_signal_array_dict
##################################################################

##################################################################
def update_min_max_array(process_signal_array_dict, signal_array_dict):
    process_min_signal = process_signal_array_dict['min_signal']
    process_max_signal = process_signal_array_dict['max_signal']

    if (process_min_signal < signal_array_dict['min_signal']):
        signal_array_dict['min_signal']=process_min_signal
    if (process_max_signal > signal_array_dict['max_signal']):
        signal_array_dict['max_signal']=process_max_signal

    process_all_keys=set(process_signal_array_dict.keys())
    process_chrom_keys=process_all_keys.difference(set(['min_signal','max_signal']))

    for chrom_key in process_chrom_keys:
        if chrom_key not in signal_array_dict:
            signal_array_dict[chrom_key]=process_signal_array_dict[chrom_key]
        else:
            signal_array_dict[chrom_key]+=process_signal_array_dict[chrom_key]
##################################################################



##################################################################
def append_chr_based_df_list_dict(chrname,chrom_based_chunk_df,chr2DataframeListDict):
    if chrname in chr2DataframeListDict:
        chr2DataframeListDict[chrname].append(chrom_based_chunk_df)
    else:
        chr2DataframeListDict[chrname]=[chrom_based_chunk_df]
##################################################################


######################################################################
#No pool not used anymore
def readWig_write_derived_from_bedgraph(outputDir, jobname, genome, library_file_with_path, occupancy_type,verbose,chunksize = 10 ** 7):
    chr2DataframeList={}
    signal_array_dict = {}
    signal_array_dict['max_signal'] = np.finfo(np.float32).min
    signal_array_dict['min_signal'] = np.finfo(np.float32).max

    wig_file_name_wo_extension = os.path.splitext(os.path.basename(library_file_with_path))[0]
    chromSizesDict = getChromSizesDict(genome)

    possible_chrom_list=list(chromSizesDict.keys())

    #To reduce memory footprint
    # Delete old chr based epigenomics files
    if occupancy_type==EPIGENOMICSOCCUPANCY:
        chrBasedNpyFilesPath=os.path.join(outputDir,jobname,DATA,occupancy_type,LIB,CHRBASED)
        deleteChrBasedNpyFiles(chrBasedNpyFilesPath)

    if os.path.exists(library_file_with_path):

        data=pd.read_csv(library_file_with_path, delimiter='\t', iterator=True, chunksize=chunksize)
        print(type(data))
        print(sys.getsizeof(data))
        print('################################')

        for chunk_df in data:
            print(type(chunk_df))
            print(sys.getsizeof(chunk_df))
            chunk_df.columns=[CHROM, START, END, SIGNAL]
            chunk_df[CHROM] = chunk_df[CHROM].astype('category')
            chunk_df[START] = chunk_df[START].astype(np.int32)
            chunk_df[END] = chunk_df[END].astype(np.int32)
            chunk_df[SIGNAL] = chunk_df[SIGNAL].astype(np.float32)
            print(sys.getsizeof(chunk_df))
            print('################################')

            chunk_df = chunk_df[chunk_df[CHROM].isin(possible_chrom_list)]
            chunk_df_grouped = chunk_df.groupby(CHROM)
            for chrname, chrom_based_chunk_df in chunk_df_grouped:
                if chrname in possible_chrom_list:
                    append_chr_based_df_list_dict(chrname,chrom_based_chunk_df,chr2DataframeList)

        for chrname in chr2DataframeList:
            chr_based_df_list=chr2DataframeList[chrname]
            process_signal_array_dict = {}

            process_signal_array = np.zeros(chromSizesDict[chrname], dtype=np.float32)
            max_signal = np.finfo(np.float32).min
            min_signal = np.finfo(np.float32).max

            for chr_based_df in chr_based_df_list:
                if (chr_based_df[SIGNAL].max() > max_signal):
                    max_signal = chr_based_df[SIGNAL].max()
                if (chr_based_df[SIGNAL].min() < min_signal):
                    min_signal = chr_based_df[SIGNAL].min()
                    # Use list comprehension
                [updateSignalArraysForListComprehension(row, process_signal_array) for row in chr_based_df.values]

            # Initialzie the list, you will return this list
            process_signal_array_dict['min_signal'] = min_signal
            process_signal_array_dict['max_signal'] = max_signal
            process_signal_array_dict[chrname]=process_signal_array

            update_min_max_array(process_signal_array_dict, signal_array_dict)

        if (occupancy_type == EPIGENOMICSOCCUPANCY):
            os.makedirs(os.path.join(outputDir, jobname, DATA, occupancy_type, LIB, CHRBASED), exist_ok=True)
        elif (occupancy_type == NUCLEOSOMEOCCUPANCY):
            os.makedirs(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME, CHRBASED),exist_ok=True)

        all_keys=set(signal_array_dict.keys())
        chrom_keys=all_keys.difference(set(['min_signal','max_signal']))

        for key in chrom_keys:
            signalArray=signal_array_dict[key]
            print('%s np.sum(signalArray,dtype=np.float32):%f' %(key,np.sum(signalArray,dtype=np.float32)))
            signalArrayFilename = '%s_signal_%s' %(key,wig_file_name_wo_extension)
            if (occupancy_type == EPIGENOMICSOCCUPANCY):
                chrBasedSignalFile = os.path.join(outputDir, jobname, DATA, occupancy_type, LIB, CHRBASED, signalArrayFilename)
            elif (occupancy_type == NUCLEOSOMEOCCUPANCY):
                chrBasedSignalFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME,CHRBASED, signalArrayFilename)
            np.save(chrBasedSignalFile,signalArray)
######################################################################

######################################################################
#For bedgraph and wig derived from bedgraph files
#Pool is used.
#Read at once.
#Write in parallel
#Use quantile for outliers
def readWig_write_derived_from_bedgraph_using_pool_read_all(outputDir, jobname, genome, library_file_with_path, occupancy_type,remove_outliers,verbose,quantileValue):
    chromSizesDict = getChromSizesDict(genome)
    possible_chrom_list = list(chromSizesDict.keys())
    if verbose: print('\tVerbose possible_chrom_list:%s' % (possible_chrom_list))

    # To reduce memory footprint
    # Delete old chr based epigenomics signal/library files
    chrBasedNpyFilesPath=os.path.join(outputDir,jobname,DATA,occupancy_type,LIB,CHRBASED)
    deleteChrBasedNpyFiles(chrBasedNpyFilesPath)

    if os.path.exists(library_file_with_path):
        numofProcesses = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(numofProcesses)
        jobs = []

        column_names = [CHROM, START, END, SIGNAL]
        data_df=pd.read_csv(library_file_with_path, sep='\t', header=None, comment='#' ,names=column_names, dtype={CHROM: 'category', START: np.int32, END: np.int32, SIGNAL: np.float32})

        if verbose: print('\tVerbose Before data_df[CHROM].unique():%s' % (data_df[CHROM].unique()))
        data_df = data_df[data_df[CHROM].isin(possible_chrom_list)]
        if verbose: print('\tVerbose After data_df[CHROM].unique():%s' % (data_df[CHROM].unique()))

        max_signal = data_df[SIGNAL].max()
        min_signal = data_df[SIGNAL].min()
        mean_skipna_True_signal = data_df[SIGNAL].mean(skipna=True)
        mean_skipna_False_signal = data_df[SIGNAL].mean(skipna=False)
        std_signal = data_df[SIGNAL].std()

        print('\nfile_df.shape:(%d,%d)' % (data_df.shape[0], data_df.shape[1]))
        print('Max Signal: %f' % max_signal)
        print('Min Signal: %f' % min_signal)
        print('Mean Signal skipna=True: %f' % mean_skipna_True_signal)
        print('Mean Signal skipna=False: %f' % mean_skipna_False_signal)
        print('Std Signal: %f' % std_signal)

        if verbose: print('\tVerbose %s Worker pid %s type(data):%s' %(occupancy_type, str(os.getpid()),type(data_df)))
        if verbose: print('\tVerbose %s Worker pid %s sys.getsizeof(data):%s' %(occupancy_type, str(os.getpid()),sys.getsizeof(data_df)))
        if verbose: print('\tVerbose ################################')

        if (remove_outliers and (quantileValue < 1.0)):
            # remove the outliers
            q = data_df[SIGNAL].quantile(quantileValue)
            print('\tq:%f' % q)
            print('\tbefore outlier elimination number of rows: %d' % (data_df.shape[0]))
            data_df = data_df[data_df[SIGNAL] < q]
            print('\tafter outlier elimination number of rows: %d' % (data_df.shape[0]))

        data_df_grouped = data_df.groupby(CHROM)

        for chrname,chrBased_data_df in data_df_grouped:
            if chrname in possible_chrom_list:
                inputList=[]
                inputList.append(outputDir)
                inputList.append(jobname)
                inputList.append(chrname)
                inputList.append(chromSizesDict[chrname])
                inputList.append(chrBased_data_df)
                inputList.append(library_file_with_path)
                inputList.append(occupancy_type)
                inputList.append(max_signal)
                inputList.append(min_signal)
                jobs.append(pool.apply_async(writeChrBasedOccupancySignalArray,args=(inputList,)))
                if verbose: print('\tVerbose JOB chrName:%s chrBased_data_df(%d,%d)' % (chrname, chrBased_data_df.shape[0], chrBased_data_df.shape[1]))

        # wait for all jobs to finish
        for job in jobs:
            if verbose: print('\tVerbose %s Worker pid %s job.get():%s ' % (occupancy_type, str(os.getpid()), job.get()))

        pool.close()
        pool.join()
######################################################################

######################################################################
#Pool is used
#Read in chunks, at the end have all signal arrays in memory and write sequentially
def readWig_write_derived_from_bedgraph_using_pool_chunks(outputDir, jobname, genome, library_file_with_path, occupancy_type,remove_outliers,verbose,chunksize = 10 ** 7,number_of_stds_for_outlier=10):
    chr2DataframeListDict={}
    signal_array_dict = {}
    signal_array_dict['max_signal'] = np.finfo(np.float32).min
    signal_array_dict['min_signal'] = np.finfo(np.float32).max

    wig_file_name_wo_extension = os.path.splitext(os.path.basename(library_file_with_path))[0]
    chromSizesDict = getChromSizesDict(genome)

    possible_chrom_list=list(chromSizesDict.keys())
    if verbose: print('\tVerbose %s Worker pid %s chrom_list:%s' % (occupancy_type, str(os.getpid()), possible_chrom_list))

    #To reduce memory footprint
    # Delete old chr based epigenomics files
    if occupancy_type==EPIGENOMICSOCCUPANCY:
        chrBasedNpyFilesPath=os.path.join(outputDir,jobname,DATA,occupancy_type,LIB,CHRBASED)
        deleteChrBasedNpyFiles(chrBasedNpyFilesPath)

    if os.path.exists(library_file_with_path):
        numofProcesses = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(numofProcesses)
        jobs = []

        #########################################################################################
        def accumulate_apply_async_result(process_signal_array_dict):
            if verbose: print('\tVerbose %s Worker pid %s Accumulate Apply Sync Result starts --- process_signal_array_dict.keys():%s' % (occupancy_type, str(os.getpid()), process_signal_array_dict.keys()))
            update_min_max_array(process_signal_array_dict,signal_array_dict)
            if verbose: print('\tVerbose %s Worker pid %s Accumulate Apply Sync Result ends --- process_signal_array_dict.keys():%s' % (occupancy_type, str(os.getpid()), process_signal_array_dict.keys()))
        #########################################################################################

        data=pd.read_csv(library_file_with_path, delimiter='\t', iterator=True, chunksize=chunksize)
        if verbose: print('\tVerbose %s Worker pid %s type(data):%s' %(occupancy_type, str(os.getpid()),type(data)))
        if verbose: print('\tVerbose %s Worker pid %s sys.getsizeof(data):%s' %(occupancy_type, str(os.getpid()),sys.getsizeof(data)))
        if verbose: print('\tVerbose ################################')

        for chunk_df in data:
            print(type(chunk_df))
            print(sys.getsizeof(chunk_df))
            chunk_df.columns=[CHROM, START, END, SIGNAL]
            chunk_df[CHROM] = chunk_df[CHROM].astype('category')
            chunk_df[START] = chunk_df[START].astype(np.int32)
            chunk_df[END] = chunk_df[END].astype(np.int32)
            chunk_df[SIGNAL] = chunk_df[SIGNAL].astype(np.float32)
            if verbose: print('\tVerbose %s Worker pid %s sys.getsizeof(chunk_df):%s' % (occupancy_type, str(os.getpid()), sys.getsizeof(chunk_df)))
            if verbose: print('\tVerbose ################################')

            chunk_df = chunk_df[chunk_df[CHROM].isin(possible_chrom_list)]
            chunk_df_grouped = chunk_df.groupby(CHROM)
            for chrname, chrom_based_chunk_df in chunk_df_grouped:
                if chrname in possible_chrom_list:
                    append_chr_based_df_list_dict(chrname,chrom_based_chunk_df,chr2DataframeListDict)

        for chrname in chr2DataframeListDict:
            chr_based_df_list=chr2DataframeListDict[chrname]
            jobs.append(pool.apply_async(fill_signal_array_dict,args=(chrname,chr_based_df_list, chromSizesDict, occupancy_type, verbose,),callback=accumulate_apply_async_result))
            if verbose: print('\tVerbose JOB chrName:%s chrom_based_chunk_df(%d,%d)' % (chrname, chrom_based_chunk_df.shape[0], chrom_based_chunk_df.shape[1]))

        # wait for all jobs to finish
        for job in jobs:
            if verbose: print('\tVerbose %s Worker pid %s job.get():%s ' % (occupancy_type, str(os.getpid()), job.get()))

        pool.close()
        pool.join()

        if verbose: print('\n\tVerbose %s Worker pid %s library_file:%s min_signal:%f max_signal:%f' % (occupancy_type, str(os.getpid()), library_file_with_path, signal_array_dict['min_signal'], signal_array_dict['max_signal']))
        if (occupancy_type == EPIGENOMICSOCCUPANCY):
            os.makedirs(os.path.join(outputDir, jobname, DATA, occupancy_type, LIB, CHRBASED), exist_ok=True)
        elif (occupancy_type == NUCLEOSOMEOCCUPANCY):
            os.makedirs(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME, CHRBASED),exist_ok=True)

        all_keys=set(signal_array_dict.keys())
        chrom_keys=all_keys.difference(set(['min_signal','max_signal']))

        if verbose: print('\tVerbose %s Worker pid %s  library_file:%s chrom_keys:%s\n' % (occupancy_type, str(os.getpid()), library_file_with_path,chrom_keys))

        for key in chrom_keys:
            signalArray=signal_array_dict[key]
            print('%s np.sum(signalArray,dtype=np.float32):%f' %(key,np.sum(signalArray,dtype=np.float32)))
            signalArrayFilename = '%s_signal_%s' %(key,wig_file_name_wo_extension)
            if (occupancy_type == EPIGENOMICSOCCUPANCY):
                chrBasedSignalFile = os.path.join(outputDir, jobname, DATA, occupancy_type, LIB, CHRBASED, signalArrayFilename)
            elif (occupancy_type == NUCLEOSOMEOCCUPANCY):
                chrBasedSignalFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME,CHRBASED, signalArrayFilename)

            if remove_outliers:
                mean=signalArray.mean()
                std=signalArray.std()
                # number_of_elements_as_outliers=np.count_nonzero(signalArray>(mean+4*std))
                signalArray[signalArray>(mean+number_of_stds_for_outlier*std)]=0
            np.save(chrBasedSignalFile,signalArray)
######################################################################


##################################################################
#SIGNAL must have type np.float32
#No Outlier Elimination is done
def readFileInBEDFormat(file_with_path,discard_signal):
    file_df=None

    print('############################################')
    if os.path.exists(file_with_path):
        file_df = pd.read_csv(file_with_path, header=None, nrows=1,sep='\t')  # 2.25 GB
        ncols=file_df.shape[1]

        if (ncols<=3):
            print('There is no enough columns in this bed file')
        elif (ncols==4):
            print('SigProfilerTopography assumes that score column is in the 4th column of this bed file and there is no header')
            file_df=pd.read_csv(file_with_path,
                                header=None,
                                usecols=[0, 1, 2, 3],
                                names = [CHROM,START,END,SIGNAL],
                                dtype={0: 'category', 1: np.int32, 2: np.int32, 3: np.float32},
                                sep='\t')

        elif ((ncols==10) or (ncols==9)):
            # ENCODE narrowpeak BED6+4 ncols=10
            # ENCODE broadpeak BED6+3 ncols=9

            if (ncols==10):
                print('ENCODE narrowpeak BED6+4')
            elif (ncols==9):
                print('ENCODE narrowpeak BED6+3')

            if discard_signal==True:
                file_df = pd.read_csv(file_with_path, header=None, usecols=[0,1,2],
                                        names=[CHROM,START,END],
                                        dtype={0: 'category', 1: np.int32, 2: np.int32},sep='\t')

            else:
                print('SigProfilerTopography assumes that signal column is in the 7th column of this bed file and there is no header')
                file_df = pd.read_csv(file_with_path, header=None, usecols=[0, 1, 2, 3, 4, 5, 6],
                                        names=[CHROM, START, END, NAME, SCORE, STRAND, SIGNAL],
                                        dtype={0: 'category', 1: np.int32, 2: np.int32, 3: str, 4: np.int32,
                                               5: 'category', 6: np.float32},sep='\t')

                # file_df.drop([3,4,5], inplace=True, axis=1)
                file_df.drop([NAME, SCORE, STRAND], inplace=True, axis=1)


        elif (ncols>=5):
            print('SigProfilerTopography assumes that score column is in the 5th column of this bed file and there is no header')
            file_df=pd.read_csv(file_with_path,header=None, usecols=[0, 1, 2, 4], names = [CHROM,START,END,SIGNAL], dtype={0: 'category', 1: np.int32, 2: np.int32, 4: np.float32},sep='\t')

        print("file_df.dtypes")
        print(file_df.dtypes)
        print('\nfile_df.columns.values')
        print(file_df.columns.values)

        print('\nfile_df.shape:(%d,%d)' %(file_df.shape[0],file_df.shape[1]))
        # print(file_df.head())

        if SIGNAL in file_df.columns.values:
            max_signal=file_df[SIGNAL].max()
            min_signal=file_df[SIGNAL].min()
            mean_skipna_True_signal=file_df[SIGNAL].mean(skipna=True)
            mean_skipna_False_signal=file_df[SIGNAL].mean(skipna=False)
            std_signal=file_df[SIGNAL].std()

            print('Max Signal: %f' %max_signal)
            print('Min Signal: %f' %min_signal)
            print('Mean Signal skipna=True: %f' %mean_skipna_True_signal)
            print('Mean Signal skipna=False: %f' %mean_skipna_False_signal)
            print('Std Signal: %f' %std_signal)
            print('file_df[SIGNAL].isnull().sum():%d' %file_df[SIGNAL].isnull().sum())

            return file_df, max_signal, min_signal
        else:
            return file_df
        print('############################################')

    return file_df
##################################################################



######################################################################
#Dec 2, 2019
#This function is used when user provides its own bed/np/narrowpeak files for occupancy analysis
#Read at once writes in parallel
def readBEDandWriteChromBasedSignalArrays(outputDir, jobname, genome, BEDFileWithPath, occupancy_type, quantileValue, remove_outliers):
    chromSizesDict = getChromSizesDict(genome)
    numofProcesses = multiprocessing.cpu_count()

    #To reduce memory footprint
    # Delete old chr based epigenomics files
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
                if chrLong in chromSizesDict:
                    chromSize = chromSizesDict[chrLong]
                    inputList = []
                    inputList.append(outputDir)
                    inputList.append(jobname)
                    inputList.append(chrLong)
                    inputList.append(chromSize)
                    inputList.append(chromBasedFileDF)
                    inputList.append(BEDFileWithPath)
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
#quantileValue set as 0.97 For  GM12878 and K562 MNase-seq data
def readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArraysInParallel(genome, quantileValue, nucleosomeFilename):
    chromSizesDict = getChromSizesDict(genome)

    # Start the pool
    numofProcesses = multiprocessing.cpu_count()
    print('Number of processors:%d' %(numofProcesses))

    print('Before nucleosome occupancy is loaded into memory')
    print('Memory usage in %s MB' %memory_usage())

    if os.path.exists(nucleosomeFilename):
        nucleosome_df_grouped, max_signal, min_signal = readNucleosomeOccupancyData(quantileValue, nucleosomeFilename)
        print('max_signal:%d min_signal:%d' %(max_signal,min_signal))
        print('np.finfo(np.float16).min:%f np.finfo(np.float16).max:%f' % (np.finfo(np.float16).min,np.finfo(np.float16).max))

        pool = multiprocessing.Pool(numofProcesses)

        print('After pool is initialized')
        print('Memory usage in %s MB' % memory_usage())
        poolInputList = []

        for chrLong, chromBasedNucleosomeDF in nucleosome_df_grouped:
            print('for %s write nucleosome signal and count array' %(chrLong))
            chromSize = chromSizesDict[chrLong]
            inputList = []
            inputList.append(None)
            inputList.append(None)
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
