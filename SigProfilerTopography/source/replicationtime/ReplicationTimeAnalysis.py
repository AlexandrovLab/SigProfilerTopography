# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018-2020 Burcak Otlu

#############################################################
# This python code generates normalized mutation density data for different analyses types which are
# Using all mutations (all chromosomes, all samples, signatures) which leads to Aggregated Substitutions
# Using all mutations but in signature based manner leads to Signature Based
# Using all mutations in signature based and sample based manner leads to Signature based and Sample Based
# Using all indels (all chromosomes, all samples) which leads to Aggregated Indels
# Using all indels (diving indels of length >=3 as microhomology-mediated indels and indels of length <3 as repeat-mediated  indels)

## In this python code
## I read the wavelet-transformed signal wig across whole genome
## Sort the signal in descending order
## Divide the data into 10 equal deciles
## Return 10 deciles as data frames
## Read the mutation file
## Group by mutation file by chromosome
## Group by decile file by chromosome
## For each chromosome (Done in the loop sequentially)
## Make a list of input
## chromosome based mutation file mutation file
## chromosome based decile files coming from 10 different deciles (Done in parallel for all of the 10 deciles)
## In the function I will generate interval tree from decile file
## Overlap the mutations with the interval tree
## Return the total number of mutations that overlap between mutation file and decile file
## Get chromosome based mutation file
## Get chromosome based decile file
## Accumulate the total number of mutations that overlap for each decile file
#############################################################

#############################################################
# Constraints, Thresholds
# Please note that for sample based and signature based replication time analysis
# We consider samples with at least 3000 mutations at total in all deciles.
#############################################################

import multiprocessing
import sys
import os
import math
import numpy as np
import pandas as pd


from SigProfilerTopography.source.commons.TopographyCommons import DATA
from SigProfilerTopography.source.commons.TopographyCommons import REPLICATIONTIME
from SigProfilerTopography.source.commons.TopographyCommons import SAMPLES
from SigProfilerTopography.source.commons.TopographyCommons import SIGNATUREBASED

from SigProfilerTopography.source.commons.TopographyCommons import NUMOFBASES

from SigProfilerTopography.source.commons.TopographyCommons import CHROM
from SigProfilerTopography.source.commons.TopographyCommons import START
from SigProfilerTopography.source.commons.TopographyCommons import END
from SigProfilerTopography.source.commons.TopographyCommons import SIGNAL

from SigProfilerTopography.source.commons.TopographyCommons import SAMPLE
from SigProfilerTopography.source.commons.TopographyCommons import LENGTH
from SigProfilerTopography.source.commons.TopographyCommons import SIMULATION_NUMBER

from SigProfilerTopography.source.commons.TopographyCommons import TYPE
from SigProfilerTopography.source.commons.TopographyCommons import SUBS
from SigProfilerTopography.source.commons.TopographyCommons import INDELS
from SigProfilerTopography.source.commons.TopographyCommons import DINUCS

from SigProfilerTopography.source.commons.TopographyCommons import AGGREGATEDSUBSTITUTIONS
from SigProfilerTopography.source.commons.TopographyCommons import AGGREGATEDINDELS
from SigProfilerTopography.source.commons.TopographyCommons import AGGREGATEDDINUCS

from SigProfilerTopography.source.commons.TopographyCommons import MICROHOMOLOGY
from SigProfilerTopography.source.commons.TopographyCommons import REPEAT

from SigProfilerTopography.source.commons.TopographyCommons import readWig_with_fixedStep_variableStep
from SigProfilerTopography.source.commons.TopographyCommons import memory_usage

from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM
from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT

from SigProfilerTopography.source.commons.TopographyCommons import Sample2NumberofDinucsDictFilename
from SigProfilerTopography.source.commons.TopographyCommons import Sample2DinucsSignature2NumberofMutationsDictFilename

from SigProfilerTopography.source.commons.TopographyCommons import getDictionary
from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofSubsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2NumberofIndelsDict

from SigProfilerTopography.source.commons.TopographyCommons import getSample2SubsSignature2NumberofMutationsDict
from SigProfilerTopography.source.commons.TopographyCommons import getSample2IndelsSignature2NumberofMutationsDict

from SigProfilerTopography.source.commons.TopographyCommons import getChrShort

from SigProfilerTopography.source.commons.TopographyCommons import MAXIMUM_NUMBER_JOBS_IN_THE_POOL_AT_ONCE
from SigProfilerTopography.source.commons.TopographyCommons import get_chrBased_simBased_combined_df_split

from SigProfilerTopography.source.commons.TopographyCommons import  MEGABYTE_IN_BYTES
from SigProfilerTopography.source.commons.TopographyCommons import  decideFileType

from SigProfilerTopography.source.commons.TopographyCommons import  get_chrBased_simBased_combined_df

##################################################################
# Please note that this dictionary is copied from .../SigProfilerMatrixGenerator/SigProfilerMatrixGenerator/scripts/SigProfilerMatrixGeneratorFunc.py line 101
# If this dictionary is updated in SigProfilerMatrixGeneratorFunc.py, it has to be updated in ReplicationTimeAnalysis.py
# Provides the reference file conversion from binary to base information
tsb_ref = {0: ['N', 'A'], 1: ['N', 'C'], 2: ['N', 'G'], 3: ['N', 'T'],
           4: ['T', 'A'], 5: ['T', 'C'], 6: ['T', 'G'], 7: ['T', 'T'],
           8: ['U', 'A'], 9: ['U', 'C'], 10: ['U', 'G'], 11: ['U', 'T'],
           12: ['B', 'A'], 13: ['B', 'C'], 14: ['B', 'G'], 15: ['B', 'T'],
           16: ['N', 'N'], 17: ['T', 'N'], 18: ['U', 'N'], 19: ['B', 'N']}
##################################################################


##################################################################
#Higher the replication time signal earlier the replication is
#Regions with high values indicate domains of early replication where initiaion occurs earlier in S-phase or early in a higher proportion of cells.
def readRepliSeqTimeData(genome,chromNamesList,repliseqDataFilename,matrix_generator_path,verbose):

    ###################################################################
    ############### Read RepliSeq Time data starts ####################
    ###################################################################

    isFileTypeBEDGRAPH = decideFileType(repliseqDataFilename)

    if isFileTypeBEDGRAPH:
        column_names = [CHROM, START, END, SIGNAL]
        replication_time_interval_version_df = pd.read_csv(repliseqDataFilename, sep='\t', header=None, comment='#', names=column_names,dtype={CHROM: 'category', START: np.int32, END: np.int32, SIGNAL: np.float32})
    else:
        #JAN 7, 2020
        replication_time_interval_version_df = readWig_with_fixedStep_variableStep(repliseqDataFilename)

    chrNamesInReplicationTimeDataArray = replication_time_interval_version_df[CHROM].unique()
    print('Before --- Chromosome names in replication time signal data: %s replication_time_interval_version_df.shape(%d,%d)\n' %(chrNamesInReplicationTimeDataArray,replication_time_interval_version_df.shape[0],replication_time_interval_version_df.shape[1]))

    #Remove rows with chromosomes that are not in chromNamesList
    replication_time_interval_version_df=replication_time_interval_version_df[replication_time_interval_version_df[CHROM].isin(chromNamesList)]

    chrNamesInReplicationTimeDataArray = replication_time_interval_version_df[CHROM].unique()
    print('After considering only chromosomes in chromNamesList --- Chromosome names in replication time signal data: %s replication_time_interval_version_df.shape(%d,%d)\n' %(chrNamesInReplicationTimeDataArray,replication_time_interval_version_df.shape[0],replication_time_interval_version_df.shape[1]))

    #Augment wavelet_processed_df with numberofAttributableBases
    wavelet_processed_augmented_df = augment(genome,replication_time_interval_version_df,matrix_generator_path,verbose)

    #Return 10 deciles: the first decile is the earliest one and the tenth decile is the latest one

    #Sort in descending order
    #Higher the replication time signal earlier the replication is
    wavelet_processed_augmented_df.sort_values(SIGNAL, ascending=False, inplace=True)


    # print('############ after sort wavelet_processed_augmented_df ###################')
    # print(wavelet_processed_augmented_df.head())
    # print('############ after sort wavelet_processed_augmented_df ###################')

    #Split wavelet_processed_augmented_df into 10 deciles
    deciles_df_list = np.array_split(wavelet_processed_augmented_df,10)
    # print('Number of decile:%d' %len(deciles))
    #deciles is a list and each decile is a dataframe <class 'pandas.core.frame.DataFrame'>
    #The first decile is the earliest one
    #The last decile is the latest one
    # print('type(deciles_df_list):%s' %type(deciles_df_list))
    # type(deciles_df_list) --> <class 'list'>

    # ############################################################
    # #For information
    # totalNumberofIntervals = 0
    # for decileIndex, decile_df in enumerate(deciles_df_list,1):
    #     # print('decileIndex: %d' %decileIndex)
    #     # print('type(decile_df)')
    #     # print(type(decile_df))
    #     # type(decile_df) --> <class 'pandas.core.frame.DataFrame'>
    #     totalNumberofIntervals += len(decile_df)
    # print('totalNumberofIntervals in all deciles : %d' %totalNumberofIntervals)
    # ############################################################

    return chrNamesInReplicationTimeDataArray, deciles_df_list
    ###################################################################
    ############### Read RepliSeq Time data ends ######################
    ###################################################################

##################################################################

##################################################################
#Jan 13, 2020
def getNumberofAttributableBasesUsingMatrixGeneratorGenome(wavelet_row,chrom_string):
    start =wavelet_row[1]
    end = wavelet_row[2]

    #old way
    #In my code ends are inclusive
    #twobitreader uses ends exlusive
    # seq_old_way = chrBasedGenome.get_slice(start, end+1)

    seq = ''
    for i in range(start,end+1,1):
        seq += tsb_ref[chrom_string[i - 1]][1]

    numofAttributableBases = seq.count('A') + seq.count('T') + seq.count('G') + seq.count('C') + seq.count('a') + seq.count('t') + seq.count('g') + seq.count('c')
    # print('######### debug starts ##############')
    # print(wavelet_row)
    # print('len(seq):%d' %len(seq))
    # print('numofAttributableBases:%d' %numofAttributableBases)
    # print('##########  debug ends #############')
    return numofAttributableBases
##################################################################


##################################################################
def addNumofAttributableBasesColumnForApplyAsync(chrLong,chrBased_wavelet_processed_df_group,chrbased_file_path,verbose):
    if verbose: print('\tVerbose Worker pid %s %s Before Augment with number of attributable bases memory_usage %.2f MB' % (str(os.getpid()), chrLong, memory_usage()))

    # 1st way Slower. Not tested in terms fo final outcome/ results.
    # chrom_string = np.memmap(chrbased_file_path, dtype=np.byte, mode='r')

    # 2nd way Faster than 1st way
    with open(chrbased_file_path, "rb") as f2:
        chrom_string = f2.read()

    resulting_df = chrBased_wavelet_processed_df_group.apply(getNumberofAttributableBasesUsingMatrixGeneratorGenome,chrom_string=chrom_string, axis= 1)

    if (len(chrBased_wavelet_processed_df_group)!=len(resulting_df)):
        print('There is a situation: len(chrBased_wavelet_processed_df_group) is not equal  to len(resulting_df)')

    chrBased_wavelet_processed_df_group[NUMOFBASES] = resulting_df

    if verbose: print('\tVerbose Worker pid %s %s After Augment with number of attributable bases memory_usage %.2f MB' % (str(os.getpid()),chrLong, memory_usage()))

    return (chrLong,chrBased_wavelet_processed_df_group)
##################################################################

##################################################################
def addNumofAttributableBasesColumn(inputList):
    chrLong = inputList[0]
    chrBased_wavelet_processed_df_group =inputList[1]

    #new way
    chrom_string=inputList[2]

    # old way
    # wholeGenome = inputList[3]
    # chrBasedGenome = wholeGenome[chrLong]
    # resulting_df = chrBased_wavelet_processed_df_group.apply(getNumberofAttributableBases, chrBasedGenome = wholeGenome[chrLong], axis= 1)

    resulting_df = chrBased_wavelet_processed_df_group.apply(getNumberofAttributableBasesUsingMatrixGeneratorGenome,chrom_string=chrom_string, axis= 1)

    # print('######## debug numberofAttributableBases starts ########')
    # print('for %s starts' %chrLong)
    # print('len(chrBased_wavelet_processed_df_group): %d ' %len(chrBased_wavelet_processed_df_group))
    # print('len(resulting_df): %d' %len(resulting_df))
    # print('type(chrBased_wavelet_processed_df_group): %s' % type(chrBased_wavelet_processed_df_group))
    # print('type(resulting_df): %s' %type(resulting_df))

    if (len(chrBased_wavelet_processed_df_group)!=len(resulting_df)):
        print('There is a situation: len(chrBased_wavelet_processed_df_group) is not equal  to len(resulting_df)')

    # print('debug starts')
    # print(resulting_df)
    # print(chrBased_wavelet_processed_df_group.columns.values.tolist())
    # print('debug ends')

    #Add the resulting_df as a new column to chrBased_wavelet_processed_df_group
    chrBased_wavelet_processed_df_group[NUMOFBASES] = resulting_df

    # print('debug starts')
    # print(chrBased_wavelet_processed_df_group.columns.values.tolist())
    # print('debug ends')

    # print('for %s ends' % chrLong)
    # print('######## debug numberofAttributableBases ends ########')

    return (chrLong,chrBased_wavelet_processed_df_group)
##################################################################

##################################################################
# April 6, 2020
def search_for_each_mutation_using_list_comprehension(
        mutation_row,
        sample2NumberofSubsDict,
        sample2NumberofIndelsDict,
        sample2NumberofDinucsDict,
        sample2SubsSignature2NumberofMutationsDict,
        sample2IndelsSignature2NumberofMutationsDict,
        sample2DinucsSignature2NumberofMutationsDict,
        chrBasedReplicationTimeDataArrayWithDecileIndex,
        simNum2Type2DecileIndex2NumberofMutationsDict,
        simNum2Sample2Type2DecileIndex2NumberofMutationsDict,
        subsSignature_cutoff_numberofmutations_averageprobability_df,
        indelsSignature_cutoff_numberofmutations_averageprobability_df,
        dinucsSignature_cutoff_numberofmutations_averageprobability_df,
        sample_based,
        df_columns):


    ###########################################
    indexofType = df_columns.index(TYPE)
    mutation_row_type = mutation_row[indexofType]

    indexofSimulationNumber = df_columns.index(SIMULATION_NUMBER)
    mutation_row_sim_num = mutation_row[indexofSimulationNumber]

    if mutation_row_type==SUBS:
        my_type=AGGREGATEDSUBSTITUTIONS
        sample2NumberofMutationsDict=sample2NumberofSubsDict
        sample2Signature2NumberofMutationsDict=sample2SubsSignature2NumberofMutationsDict
        signature_cutoff_numberofmutations_averageprobability_df=subsSignature_cutoff_numberofmutations_averageprobability_df
    elif mutation_row_type == INDELS:
        my_type=AGGREGATEDINDELS
        sample2NumberofMutationsDict=sample2NumberofIndelsDict
        sample2Signature2NumberofMutationsDict=sample2IndelsSignature2NumberofMutationsDict
        signature_cutoff_numberofmutations_averageprobability_df=indelsSignature_cutoff_numberofmutations_averageprobability_df
    elif mutation_row_type == DINUCS:
        my_type=AGGREGATEDDINUCS
        sample2NumberofMutationsDict=sample2NumberofDinucsDict
        sample2Signature2NumberofMutationsDict=sample2DinucsSignature2NumberofMutationsDict
        signature_cutoff_numberofmutations_averageprobability_df=dinucsSignature_cutoff_numberofmutations_averageprobability_df
    ###########################################


    # For single point mutations start and end are the same, therefore we need to add 1 to mutation_row[END]
    indexofSample = df_columns.index(SAMPLE)
    sample = mutation_row[indexofSample]

    indexofStart = df_columns.index(START)
    start=mutation_row[indexofStart]

    #end is exclusive for subs, indels and dincuc provided by readChrBased methods
    # start= mutation_row[START]

    if (my_type==AGGREGATEDINDELS):
        # end = start + mutation_row[LENGTH]
        indexofLength = df_columns.index(LENGTH)
        length=mutation_row[indexofLength]
        end = start + int(length)
    elif (my_type==AGGREGATEDSUBSTITUTIONS):
        end = start + 1
    elif (my_type==AGGREGATEDDINUCS):
        end = start + 2

    slicedArray = chrBasedReplicationTimeDataArrayWithDecileIndex[start:end]

    # np.nonzero returns the indices of the elements that are non-zero.
    # np.unique finds the unique elements of an array returns ndarray the sorted unique values.
    # np.nditer efficient multi-dimensional iterator object to iterate over arrays.
    uniqueIndexesArray = np.unique(slicedArray[np.nonzero(slicedArray)])

    if (uniqueIndexesArray.size>0):
        for decileIndex in np.nditer(uniqueIndexesArray):
            # type(decileIndex) is numpy.ndarray
            decileIndexNum = int(decileIndex)

            ################### AGGREGATED starts #####################
            if mutation_row_sim_num in simNum2Type2DecileIndex2NumberofMutationsDict:
                if my_type in simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num]:
                    if decileIndexNum in simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][my_type]:
                        simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][my_type][decileIndexNum] += 1
                    else:
                        simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][my_type][decileIndexNum] = 1
                else:
                    simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][my_type]={}
                    simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][my_type][decileIndexNum] = 1
            else:
                simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num] = {}
                simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][my_type] = {}
                simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][my_type][decileIndexNum] = 1

            if (my_type==AGGREGATEDINDELS):
                if length >= 3:
                    if mutation_row_sim_num in simNum2Type2DecileIndex2NumberofMutationsDict:
                        if MICROHOMOLOGY in simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num]:
                            if decileIndexNum in simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][MICROHOMOLOGY]:
                                simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][MICROHOMOLOGY][decileIndexNum] += 1
                            else:
                                simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][MICROHOMOLOGY][decileIndexNum] = 1
                        else:
                            simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][MICROHOMOLOGY] = {}
                            simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][MICROHOMOLOGY][decileIndexNum] = 1
                    else:
                        simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num]={}
                        simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][MICROHOMOLOGY] = {}
                        simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][MICROHOMOLOGY][decileIndexNum] = 1

                else:
                    if mutation_row_sim_num in simNum2Type2DecileIndex2NumberofMutationsDict:
                        if REPEAT in simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num]:
                            if decileIndexNum in simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][REPEAT]:
                                simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][REPEAT][decileIndexNum] += 1
                            else:
                                simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][REPEAT][decileIndexNum] = 1
                        else:
                            simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][REPEAT] = {}
                            simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][REPEAT][decileIndexNum] = 1
                    else:
                        simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num] = {}
                        simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][REPEAT] = {}
                        simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][REPEAT][decileIndexNum] = 1
            ################### AGGREGATED ends #######################

            ########################### Signatures start ###########################
            for signature in signature_cutoff_numberofmutations_averageprobability_df['signature'].unique():
                indexofSignature=df_columns.index(signature)
                cutoff=float(signature_cutoff_numberofmutations_averageprobability_df[signature_cutoff_numberofmutations_averageprobability_df['signature']==signature]['cutoff'].values[0])
                if mutation_row[indexofSignature] >= cutoff:
                    if mutation_row_sim_num in simNum2Type2DecileIndex2NumberofMutationsDict:
                        if signature in simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num]:
                            if decileIndexNum in simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][signature]:
                                simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][signature][decileIndexNum] += 1
                            else:
                                simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][signature][decileIndexNum] = 1
                        else:
                            simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][signature] = {}
                            simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][signature][decileIndexNum] = 1
                    else:
                        simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num] = {}
                        simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][signature] = {}
                        simNum2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][signature][decileIndexNum] = 1

            ########################### Signatures end #############################

            # ############################# Sample Based starts ###########################
            # if sample_based:
            #
            #     if sample in sample2NumberofMutationsDict:
            #         ############## Sample Based Aggregated Substitutions starts ############
            #         if mutation_row_sim_num in simNum2Sample2Type2DecileIndex2NumberofMutationsDict:
            #             if sample in simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num]:
            #                 if my_type in simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample]:
            #
            #                     if decileIndexNum in simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][my_type]:
            #                         simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][my_type][decileIndexNum] += 1
            #                     else:
            #                         simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][my_type][decileIndexNum] = 1
            #
            #                     if (my_type == AGGREGATEDINDELS):
            #                         if length >= 3:
            #                             if MICROHOMOLOGY in simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample]:
            #                                 if decileIndexNum in simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][MICROHOMOLOGY]:
            #                                     simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][MICROHOMOLOGY][decileIndexNum] += 1
            #                                 else:
            #                                     simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][MICROHOMOLOGY][decileIndexNum] = 1
            #                             else:
            #                                 simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][MICROHOMOLOGY] = {}
            #                                 simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][MICROHOMOLOGY][decileIndexNum] = 1
            #
            #                         else:
            #                             if REPEAT in simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample]:
            #                                 if decileIndexNum in simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][REPEAT]:
            #                                     simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][REPEAT][decileIndexNum] += 1
            #                                 else:
            #                                     simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][REPEAT][decileIndexNum] = 1
            #                             else:
            #                                 simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][REPEAT] = {}
            #                                 simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][REPEAT][decileIndexNum] = 1
            #
            #                 else:
            #                     simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][my_type] = {}
            #                     simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][my_type][decileIndexNum] = 1
            #                     if my_type == AGGREGATEDINDELS:
            #                         if length >= 3:
            #                             simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][MICROHOMOLOGY] = {}
            #                             simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][MICROHOMOLOGY][decileIndexNum] = 1
            #                         else:
            #                             simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][REPEAT] = {}
            #                             simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][REPEAT][decileIndexNum] = 1
            #
            #             else:
            #                 simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample] = {}
            #                 simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][my_type] = {}
            #                 simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][my_type][decileIndexNum] = 1
            #                 if my_type == AGGREGATEDINDELS:
            #                     if length >= 3:
            #                         simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][MICROHOMOLOGY] = {}
            #                         simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][MICROHOMOLOGY][decileIndexNum] = 1
            #                     else:
            #                         simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][REPEAT] = {}
            #                         simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][REPEAT][decileIndexNum] = 1
            #
            #         else:
            #             simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num]={}
            #             simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample] = {}
            #             simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][my_type] = {}
            #             simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][my_type][decileIndexNum] = 1
            #             if my_type == AGGREGATEDINDELS:
            #                 if length >= 3:
            #                     simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][MICROHOMOLOGY] = {}
            #                     simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][MICROHOMOLOGY][decileIndexNum] = 1
            #                 else:
            #                     simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][REPEAT] = {}
            #                     simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][REPEAT][decileIndexNum] = 1
            #         ############## Sample Based Aggregated Substitutions ends ##############
            #
            #
            # ############## Sample Based Signatures ends ############################
            # if sample in sample2Signature2NumberofMutationsDict:
            #     for signature in sample2Signature2NumberofMutationsDict[sample]:
            #         indexofSignature = df_columns.index(signature)
            #         cutoff=float(signature_cutoff_numberofmutations_averageprobability_df[signature_cutoff_numberofmutations_averageprobability_df['signature']==signature]['cutoff'].values[0])
            #         if mutation_row[indexofSignature] >= cutoff:
            #             if mutation_row_sim_num in simNum2Sample2Type2DecileIndex2NumberofMutationsDict:
            #                 if decileIndexNum in simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][signature]:
            #                     simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][signature][decileIndexNum] += 1
            #                 else:
            #                     simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][signature][decileIndexNum] = 1
            #
            #             else:
            #                 simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num]={}
            #                 simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample] = {}
            #                 simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][signature] = {}
            #                 simNum2Sample2Type2DecileIndex2NumberofMutationsDict[mutation_row_sim_num][sample][signature][decileIndexNum] = 1
            # ############## Sample Based Signatures ends ############################
            #
            # ############################# Sample Based ends ############################

##################################################################



##################################################################
#April 6, 2020
def searchforAllMutations(
        sample2NumberofSubsDict,
        sample2NumberofIndelsDict,
        sample2NumberofDinucsDict,
        sample2SubsSignature2NumberofMutationsDict,
        sample2IndelsSignature2NumberofMutationsDict,
        sample2DinucsSignature2NumberofMutationsDict,
        chrBased_simBased_combined_df_split,
        chrBasedReplicationTimeDataArrayWithDecileIndex,
       subsSignature_cutoff_numberofmutations_averageprobability_df,
       indelsSignature_cutoff_numberofmutations_averageprobability_df,
       dinucsSignature_cutoff_numberofmutations_averageprobability_df,
       sample_based):

    ######################################################################
    simNum2Type2DecileIndex2NumberofMutationsDict = {}
    simNum2Sample2Type2DecileIndex2NumberofMutationsDict = {}

    df_columns=list(chrBased_simBased_combined_df_split.columns.values)
    [search_for_each_mutation_using_list_comprehension(mutation_row,
                                                   sample2NumberofSubsDict,
                                                   sample2NumberofIndelsDict,
                                                   sample2NumberofDinucsDict,
                                                   sample2SubsSignature2NumberofMutationsDict,
                                                   sample2IndelsSignature2NumberofMutationsDict,
                                                   sample2DinucsSignature2NumberofMutationsDict,
                                                   chrBasedReplicationTimeDataArrayWithDecileIndex,
                                                    simNum2Type2DecileIndex2NumberofMutationsDict,
                                                    simNum2Sample2Type2DecileIndex2NumberofMutationsDict,
                                                    subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                    indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                                    dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                                    sample_based,
                                                    df_columns) for mutation_row in chrBased_simBased_combined_df_split.values]

    return simNum2Type2DecileIndex2NumberofMutationsDict, simNum2Sample2Type2DecileIndex2NumberofMutationsDict

##################################################################



##################################################################
#Please notice that replication time data are not overlappig data therefore only setting one decileIndex will be correct.
# e.g.:
# start end
#   10 1009
# 1010 2009
def fillArray(chrBased_replicationtimedata_row,chrBasedDecileIndexArray,decileIndex):
    start= chrBased_replicationtimedata_row[START]
    end = chrBased_replicationtimedata_row[END] + 1
    chrBasedDecileIndexArray[start:end] = decileIndex
##################################################################

##################################################################
#Explanation of this function
#This function is called for each chromosome
#For each chromosome, a numpy array is filled.
#If there is a chromosome locus with a decile index 8 let's say, in the array that locus is filled with 8.
#Decile index can be between 1-10.
#In this function, each available interval with an index is filled in the corresponding numpy array with chrom size
def fillChrBasedReplicationTimeNPArray(chrLong,chromSize,chrBased_grouped_decile_df_list):
    #We can set the starting index as 1 in builtin function enumerate
    #First chrBased_grouped_decile has index of 1
    #Last chrBased_grouped_decile has index of 10

    # int8	Byte (-128 to 127)
    chrBasedReplicationTimeDataArrayWithDecileIndex = np.zeros(chromSize, dtype=np.int8)

    #First decileIndex is 1, last decile index is 10.
    for decileIndex, chrBased_grouped_decile_df in enumerate(chrBased_grouped_decile_df_list,1):
        # Solution to keyError
        for name, chrBased_replicationtimedata_df in chrBased_grouped_decile_df:
            if (chrLong==name) and (chrBased_replicationtimedata_df is not None) and  (not chrBased_replicationtimedata_df.empty):
                #what is chrBased_decile's type? DataFrame
                chrBased_replicationtimedata_df.apply(fillArray,chrBasedDecileIndexArray=chrBasedReplicationTimeDataArrayWithDecileIndex,decileIndex=decileIndex, axis=1)

        # # Please note that although key exists if dataframe to be returned is empty dataframe it gives KeyError
        # if chrLong in chrBased_grouped_decile_df.groups.keys():
        #     if chrLong in chrBased_grouped_decile_df.groups:
        #         chrBased_replicationtimedata_df = chrBased_grouped_decile_df.get_group(chrLong)
        #         print('DEBUG %s decileIndex:%d chrBased_replicationtimedata_df.shape(%d,%d)' %(chrLong,decileIndex,chrBased_replicationtimedata_df.shape[0],chrBased_replicationtimedata_df.shape[1]))
        #     else:
        #         chrBased_replicationtimedata_df=None
        #     #what is chrBased_decile's type? DataFrame
        #     if ((chrBased_replicationtimedata_df is not None) and (not chrBased_replicationtimedata_df.empty)):
        #         chrBased_replicationtimedata_df.apply(fillArray,chrBasedDecileIndexArray=chrBasedReplicationTimeDataArrayWithDecileIndex,decileIndex=decileIndex, axis=1)

    return chrBasedReplicationTimeDataArrayWithDecileIndex
##################################################################

##################################################################
#April 6, 2020
def fillDictionaries_chromBased_simBased_for_all_mutations(
        sample2NumberofSubsDict,
        sample2NumberofIndelsDict,
        sample2NumberofDinucsDict,
        sample2SubsSignature2NumberofMutationsDict,
        sample2IndelsSignature2NumberofMutationsDict,
        sample2DinucsSignature2NumberofMutationsDict,
        chrBasedReplicationTimeDataArrayWithDecileIndex,
        chrBased_simBased_combined_df_split,
        sample_based,
        simNum,
        chrLong,
        subsSignature_cutoff_numberofmutations_averageprobability_df,
        indelsSignature_cutoff_numberofmutations_averageprobability_df,
        dinucsSignature_cutoff_numberofmutations_averageprobability_df,
        verbose):


    simNum2Type2DecileIndex2NumberofMutationsDict = {}
    simNum2Sample2Type2DecileIndex2NumberofMutationsDict = {}

    ##################################################################################################################
    if ((chrBased_simBased_combined_df_split is not None) and (not chrBased_simBased_combined_df_split.empty)):
        if verbose: print('\tVerbose Worker %s Search for all mutations STARTS %s simNum:%d %s MB' % (str(os.getpid()), chrLong, simNum, memory_usage()),flush=True)
        simNum2Type2DecileIndex2NumberofMutationsDict, simNum2Sample2Type2DecileIndex2NumberofMutationsDict = searchforAllMutations(
            sample2NumberofSubsDict,
            sample2NumberofIndelsDict,
            sample2NumberofDinucsDict,
            sample2SubsSignature2NumberofMutationsDict,
            sample2IndelsSignature2NumberofMutationsDict,
            sample2DinucsSignature2NumberofMutationsDict,
            chrBased_simBased_combined_df_split,
            chrBasedReplicationTimeDataArrayWithDecileIndex,
            subsSignature_cutoff_numberofmutations_averageprobability_df,
            indelsSignature_cutoff_numberofmutations_averageprobability_df,
            dinucsSignature_cutoff_numberofmutations_averageprobability_df,
            sample_based)
        if verbose: print('\tVerbose Worker %s Search for all mutations ENDS %s simNum:%d %s MB' % (str(os.getpid()), chrLong, simNum, memory_usage()),flush=True)
    else:
        print('%s simNum:%d None' %(chrLong,simNum),flush=True)
    ##################################################################################################################

    return  simNum2Type2DecileIndex2NumberofMutationsDict, simNum2Sample2Type2DecileIndex2NumberofMutationsDict
##################################################################

##################################################################
def combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray_simbased_chrbased_splitbased(outputDir,
                                                                                        jobname,
                                                                                        chrLong,
                                                                                        chromSize,
                                                                                        simNum,
                                                                                        splitIndex,
                                                                                        chrBased_grouped_decile_df_list,
                                                                                        subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                                        indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                                        dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                                        sample_based,
                                                                                        sample2NumberofSubsDict,
                                                                                        sample2NumberofIndelsDict,
                                                                                        sample2NumberofDinucsDict,
                                                                                        sample2SubsSignature2NumberofMutationsDict,
                                                                                        sample2IndelsSignature2NumberofMutationsDict,
                                                                                        sample2DinucsSignature2NumberofMutationsDict,
                                                                                        verbose):

    chrBased_simBased_combined_df_split = get_chrBased_simBased_combined_df_split(outputDir, jobname, chrLong, simNum, splitIndex)

    return combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray(chrLong,
                                                                        chromSize,
                                                                        simNum,
                                                                        chrBased_grouped_decile_df_list,
                                                                        chrBased_simBased_combined_df_split,
                                                                        subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                        indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                        dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                        sample_based,
                                                                        sample2NumberofSubsDict,
                                                                        sample2NumberofIndelsDict,
                                                                        sample2NumberofDinucsDict,
                                                                        sample2SubsSignature2NumberofMutationsDict,
                                                                        sample2IndelsSignature2NumberofMutationsDict,
                                                                        sample2DinucsSignature2NumberofMutationsDict,
                                                                        verbose)
##################################################################


##################################################################
# May 10, 2020
def combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray_simbased_chrbased(outputDir, jobname, chrLong, chromSize,simNum,chrBased_grouped_decile_df_list,
                                                                                        subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                                        indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                                        dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                                        sample_based,
                                                                                        sample2NumberofSubsDict,
                                                                                        sample2NumberofIndelsDict,
                                                                                        sample2NumberofDinucsDict,
                                                                                        sample2SubsSignature2NumberofMutationsDict,
                                                                                        sample2IndelsSignature2NumberofMutationsDict,
                                                                                        sample2DinucsSignature2NumberofMutationsDict,
                                                                                        verbose):

    chrBased_simBased_combined_df = get_chrBased_simBased_combined_df(outputDir, jobname, chrLong, simNum)

    return combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray(chrLong,
                                                                       chromSize,
                                                                       simNum,
                                                                       chrBased_grouped_decile_df_list,
                                                                       chrBased_simBased_combined_df,
                                                                        subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                        indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                        dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                        sample_based,
                                                                        sample2NumberofSubsDict,
                                                                        sample2NumberofIndelsDict,
                                                                        sample2NumberofDinucsDict,
                                                                        sample2SubsSignature2NumberofMutationsDict,
                                                                        sample2IndelsSignature2NumberofMutationsDict,
                                                                        sample2DinucsSignature2NumberofMutationsDict,
                                                                        verbose)

##################################################################



##################################################################
#April 9, 2020
def combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray(chrLong,chromSize,simNum,chrBased_grouped_decile_df_list,
                                                                       chrBased_simBased_combined_df_split,
                                                                        subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                        indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                        dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                        sample_based,
                                                                        sample2NumberofSubsDict,
                                                                        sample2NumberofIndelsDict,
                                                                        sample2NumberofDinucsDict,
                                                                        sample2SubsSignature2NumberofMutationsDict,
                                                                        sample2IndelsSignature2NumberofMutationsDict,
                                                                        sample2DinucsSignature2NumberofMutationsDict,
                                                                        verbose):

    #Fill replication time numpy array
    chrBasedReplicationTimeDataArrayWithDecileIndex = fillChrBasedReplicationTimeNPArray(chrLong, chromSize,chrBased_grouped_decile_df_list)

    #Search for (chrLong,simNum,splitIndex) tuple
    simNum2Type2DecileIndex2NumberofMutationsDict, simNum2Sample2Type2DecileIndex2NumberofMutationsDict = fillDictionaries_chromBased_simBased_for_all_mutations(
        sample2NumberofSubsDict,
        sample2NumberofIndelsDict,
        sample2NumberofDinucsDict,
        sample2SubsSignature2NumberofMutationsDict,
        sample2IndelsSignature2NumberofMutationsDict,
        sample2DinucsSignature2NumberofMutationsDict,
        chrBasedReplicationTimeDataArrayWithDecileIndex,
        chrBased_simBased_combined_df_split,
        sample_based,
        simNum,
        chrLong,
        subsSignature_cutoff_numberofmutations_averageprobability_df,
        indelsSignature_cutoff_numberofmutations_averageprobability_df,
        dinucsSignature_cutoff_numberofmutations_averageprobability_df,
        verbose)

    return simNum2Type2DecileIndex2NumberofMutationsDict, simNum2Sample2Type2DecileIndex2NumberofMutationsDict
##################################################################


##################################################################
def getNormalizedMutationDensityList(mutationDensityDict):
    densities = mutationDensityDict.values()
    maxDensity = max(densities)
    if maxDensity>0:
        normalizedMutationDensities = [x/maxDensity for x in densities]
    else:
        normalizedMutationDensities = densities

    return normalizedMutationDensities
##################################################################



##################################################################
#March 22, 2019 starts
def getMutationDensityDictNewVersion(decileIndex2NumberofAttributableBasesDict,decileBasedAllChrAccumulatedCountDict):
    decileBasedMutationDensityDict = {}
    numberofMutations = 0

    numberofMutationsList = []
    numberofAttributableBasesList = []

    for decileIndex in decileIndex2NumberofAttributableBasesDict:
        if (decileIndex in decileBasedAllChrAccumulatedCountDict):
            count = decileBasedAllChrAccumulatedCountDict[decileIndex]
            numberofMutations += count
            numofAttBases = decileIndex2NumberofAttributableBasesDict[decileIndex]
            mutationDensity = float(count) / numofAttBases
            decileBasedMutationDensityDict[decileIndex] = mutationDensity
            numberofMutationsList.append(count)
            numberofAttributableBasesList.append(numofAttBases)

        else:
            decileBasedMutationDensityDict[decileIndex] = 0
            numberofMutationsList.append(0)
            numberofAttributableBasesList.append(0)


    return numberofMutations, decileBasedMutationDensityDict, numberofMutationsList, numberofAttributableBasesList
#March 22, 2019 starts
##################################################################

##################################################################
def getMutationDensityDict(decile_df_list,decileBasedAllChrAccumulatedCountDict):
    decileBasedMutationDensityDict = {}
    numberofMutations = 0

    numberofMutationsList = []
    numberofAttributableBasesList = []

    #Modifiled as enumerate(deciles,1) formerly it was enumerate(deciles,0)
    for i,decile_df in enumerate(decile_df_list,1):
        if (i in decileBasedAllChrAccumulatedCountDict):
            count = decileBasedAllChrAccumulatedCountDict[i]
            numberofMutations += count
            numofAttBases = decile_df[NUMOFBASES].sum()
            mutationDensity=float(count)/numofAttBases
            decileBasedMutationDensityDict[i] = mutationDensity
            numberofMutationsList.append(count)
            numberofAttributableBasesList.append(numofAttBases)
        else:
            decileBasedMutationDensityDict[i] = 0
            numberofMutationsList.append(0)
            numberofAttributableBasesList.append(0)

        # print('decile: %d numofAttBases: %d' %(i,numofAttBases))

    return numberofMutations, decileBasedMutationDensityDict, numberofMutationsList, numberofAttributableBasesList
##################################################################


#########################################################################################
def accumulate_result(simNum2Type2DecileIndex2NumberofMutationsDict,simNum2Sample2Type2DecileIndex2NumberofMutationsDict,simNum2Type2DecileIndex2AccumulatedNumberofMutationsDict,simNum2Sample2Type2DecileIndex2AccumulatedNumberofMutationsDict,verbose):

    if verbose: print('\tVerbose Worker pid %s memory_usage %.2f MB type2DecileIndex2NumberofMutationsDict: %.2f MB sample2Type2DecileIndex2NumberofMutationsDict: %.2f MB  simNum2Type2DecileIndex2AccumulatedNumberofMutationsDict: %.2f MB simNum2Sample2Type2DecileIndex2AccumulatedNumberofMutationsDict: %.2f MB BEFORE ACCUMULATE' % (str(os.getpid()), memory_usage(),sys.getsizeof(simNum2Type2DecileIndex2NumberofMutationsDict)/MEGABYTE_IN_BYTES,sys.getsizeof(simNum2Sample2Type2DecileIndex2NumberofMutationsDict)/MEGABYTE_IN_BYTES,sys.getsizeof(simNum2Type2DecileIndex2AccumulatedNumberofMutationsDict)/MEGABYTE_IN_BYTES,sys.getsizeof(simNum2Sample2Type2DecileIndex2AccumulatedNumberofMutationsDict)/MEGABYTE_IN_BYTES), flush=True)

    ###########################################################################
    for simNum in simNum2Type2DecileIndex2NumberofMutationsDict:
        for my_type in simNum2Type2DecileIndex2NumberofMutationsDict[simNum]:
            for decileIndex in simNum2Type2DecileIndex2NumberofMutationsDict[simNum][my_type]:
                if simNum in simNum2Type2DecileIndex2AccumulatedNumberofMutationsDict:
                    if my_type in simNum2Type2DecileIndex2AccumulatedNumberofMutationsDict[simNum]:
                        if decileIndex in simNum2Type2DecileIndex2AccumulatedNumberofMutationsDict[simNum][my_type]:
                            simNum2Type2DecileIndex2AccumulatedNumberofMutationsDict[simNum][my_type][decileIndex] += simNum2Type2DecileIndex2NumberofMutationsDict[simNum][my_type][decileIndex]
                        else:
                            simNum2Type2DecileIndex2AccumulatedNumberofMutationsDict[simNum][my_type][decileIndex] = simNum2Type2DecileIndex2NumberofMutationsDict[simNum][my_type][decileIndex]
                    else:
                        simNum2Type2DecileIndex2AccumulatedNumberofMutationsDict[simNum][my_type]={}
                        simNum2Type2DecileIndex2AccumulatedNumberofMutationsDict[simNum][my_type][decileIndex] = simNum2Type2DecileIndex2NumberofMutationsDict[simNum][my_type][decileIndex]
                else:
                    simNum2Type2DecileIndex2AccumulatedNumberofMutationsDict[simNum] = {}
                    simNum2Type2DecileIndex2AccumulatedNumberofMutationsDict[simNum][my_type] = {}
                    simNum2Type2DecileIndex2AccumulatedNumberofMutationsDict[simNum][my_type][decileIndex] = simNum2Type2DecileIndex2NumberofMutationsDict[simNum][my_type][decileIndex]
    ###########################################################################


    ######################## Sample Based starts #######################################
    for simNum in simNum2Sample2Type2DecileIndex2NumberofMutationsDict:
        for sample in simNum2Sample2Type2DecileIndex2NumberofMutationsDict[simNum]:
            for my_type in simNum2Sample2Type2DecileIndex2NumberofMutationsDict[simNum][sample]:
                for decileIndex in simNum2Sample2Type2DecileIndex2NumberofMutationsDict[simNum][sample][my_type]:
                    if simNum in simNum2Sample2Type2DecileIndex2AccumulatedNumberofMutationsDict:
                        if sample in simNum2Sample2Type2DecileIndex2AccumulatedNumberofMutationsDict[simNum]:
                            if my_type in simNum2Sample2Type2DecileIndex2AccumulatedNumberofMutationsDict[simNum][sample]:
                                if decileIndex in simNum2Sample2Type2DecileIndex2AccumulatedNumberofMutationsDict[simNum][sample][my_type]:
                                    simNum2Sample2Type2DecileIndex2AccumulatedNumberofMutationsDict[simNum][sample][my_type][decileIndex] += simNum2Sample2Type2DecileIndex2NumberofMutationsDict[simNum][sample][my_type][decileIndex]
                                else:
                                    simNum2Sample2Type2DecileIndex2AccumulatedNumberofMutationsDict[simNum][sample][my_type][decileIndex] = simNum2Sample2Type2DecileIndex2NumberofMutationsDict[simNum][sample][my_type][decileIndex]
                            else:
                                simNum2Sample2Type2DecileIndex2AccumulatedNumberofMutationsDict[simNum][sample][my_type] = {}
                                simNum2Sample2Type2DecileIndex2AccumulatedNumberofMutationsDict[simNum][sample][my_type][decileIndex] = simNum2Sample2Type2DecileIndex2NumberofMutationsDict[simNum][sample][my_type][decileIndex]
                        else:
                            simNum2Sample2Type2DecileIndex2AccumulatedNumberofMutationsDict[simNum][sample] = {}
                            simNum2Sample2Type2DecileIndex2AccumulatedNumberofMutationsDict[simNum][sample][my_type] = {}
                            simNum2Sample2Type2DecileIndex2AccumulatedNumberofMutationsDict[simNum][sample][my_type][decileIndex] = simNum2Sample2Type2DecileIndex2NumberofMutationsDict[simNum][sample][my_type][decileIndex]
                    else:
                        simNum2Sample2Type2DecileIndex2AccumulatedNumberofMutationsDict[simNum] = {}
                        simNum2Sample2Type2DecileIndex2AccumulatedNumberofMutationsDict[simNum][sample] = {}
                        simNum2Sample2Type2DecileIndex2AccumulatedNumberofMutationsDict[simNum][sample][my_type] = {}
                        simNum2Sample2Type2DecileIndex2AccumulatedNumberofMutationsDict[simNum][sample][my_type][decileIndex] = simNum2Sample2Type2DecileIndex2NumberofMutationsDict[simNum][sample][my_type][decileIndex]
    ######################## Sample Based ends #########################################

    if verbose: print('\tVerbose Worker pid %s memory_usage %.2f MB type2DecileIndex2NumberofMutationsDict: %.2f MB sample2Type2DecileIndex2NumberofMutationsDict: %.2f MB  simNum2Type2DecileIndex2AccumulatedNumberofMutationsDict: %.2f MB simNum2Sample2Type2DecileIndex2AccumulatedNumberofMutationsDict: %.2f MB AFTER ACCUMULATE' % (str(os.getpid()), memory_usage(),sys.getsizeof(simNum2Type2DecileIndex2NumberofMutationsDict)/MEGABYTE_IN_BYTES,sys.getsizeof(simNum2Sample2Type2DecileIndex2NumberofMutationsDict)/MEGABYTE_IN_BYTES,sys.getsizeof(simNum2Type2DecileIndex2AccumulatedNumberofMutationsDict)/MEGABYTE_IN_BYTES,sys.getsizeof(simNum2Sample2Type2DecileIndex2AccumulatedNumberofMutationsDict)/MEGABYTE_IN_BYTES), flush=True)
#########################################################################################



##################################################################
def calculateCountsForMutationsFillingReplicationTimeNPArrayRuntime(computationType,
        outputDir,
        jobname,
        numofSimulations,
        job_tuples,
        sample2NumberofSubsDict,
        sample2NumberofIndelsDict,
        sample2NumberofDinucsDict,
        sample2SubsSignature2NumberofMutationsDict,
        sample2IndelsSignature2NumberofMutationsDict,
        sample2DinucsSignature2NumberofMutationsDict,
        chromSizesDict,
        chromNamesList,
        chrBased_grouped_decile_df_list,
        subsSignature_cutoff_numberofmutations_averageprobability_df,
        indelsSignature_cutoff_numberofmutations_averageprobability_df,
        dinucsSignature_cutoff_numberofmutations_averageprobability_df,
        sample_based,
        verbose):

    ################################
    simNum2Type2DecileIndex2AccumulatedNumberofMutationsDict = {}
    simNum2Sample2Type2DecileIndex2AccumulatedNumberofMutationsDict = {}
    ################################

    #########################################################################################
    def accumulate_apply_async_result(result_tuple):
        simNum2Type2DecileIndex2NumberofMutationsDict = result_tuple[0]
        simNum2Sample2Type2DecileIndex2NumberofMutationsDict = result_tuple[1]

        print('MONITOR ACCUMULATE', flush=True)

        accumulate_result(simNum2Type2DecileIndex2NumberofMutationsDict,
                          simNum2Sample2Type2DecileIndex2NumberofMutationsDict,
                          simNum2Type2DecileIndex2AccumulatedNumberofMutationsDict,
                          simNum2Sample2Type2DecileIndex2AccumulatedNumberofMutationsDict,
                          verbose)
    #########################################################################################


    #######################################################################################################################
    if (computationType == USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM):

        sim_nums = range(0, numofSimulations + 1)
        sim_num_chr_tuples = ((sim_num, chrLong) for sim_num in sim_nums for chrLong in chromNamesList)

        ################################
        numofProcesses = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=numofProcesses)
        ################################

        ################################
        jobs = []
        ################################

        for simNum, chrLong in sim_num_chr_tuples:
            chromSize = chromSizesDict[chrLong]
            jobs.append(pool.apply_async(combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray_simbased_chrbased,
                                 args=(outputDir,
                                       jobname,
                                       chrLong,
                                       chromSize,
                                       simNum,
                                       chrBased_grouped_decile_df_list,
                                       subsSignature_cutoff_numberofmutations_averageprobability_df,
                                       indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                       dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                       sample_based,
                                       sample2NumberofSubsDict,
                                       sample2NumberofIndelsDict,
                                       sample2NumberofDinucsDict,
                                       sample2SubsSignature2NumberofMutationsDict,
                                       sample2IndelsSignature2NumberofMutationsDict,
                                       sample2DinucsSignature2NumberofMutationsDict,
                                       verbose,),
                                 callback=accumulate_apply_async_result))
            print('MONITOR %s simNum:%d len(jobs):%d' % (chrLong, simNum, len(jobs)), flush=True)
        ################################################################################

        ##############################################################################
        # wait for all jobs to finish
        for job in jobs:
            if verbose: print('\tVerbose Replication Strand Bias Worker pid %s job.get():%s ' % (str(os.getpid()), job.get()))
        ##############################################################################

        ################################
        pool.close()
        pool.join()
        ################################

    #######################################################################################################################

    ######################## starts July 22 2020 ###############################
    elif (computationType==USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT):

        ################################
        numofProcesses = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=numofProcesses)
        ################################

        ################################
        jobs=[]
        ################################

        ################################
        for chrLong, simNum, splitIndex in job_tuples:
            chromSize = chromSizesDict[chrLong]
            jobs.append(pool.apply_async(combined_generateReplicationTimeNPArrayAndSearchMutationsOnNPArray_simbased_chrbased_splitbased,
                                 args=(outputDir,
                                       jobname,
                                       chrLong,
                                       chromSize,
                                       simNum,
                                       splitIndex,
                                       chrBased_grouped_decile_df_list,
                                       subsSignature_cutoff_numberofmutations_averageprobability_df,
                                       indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                       dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                       sample_based,
                                       sample2NumberofSubsDict,
                                       sample2NumberofIndelsDict,
                                       sample2NumberofDinucsDict,
                                       sample2SubsSignature2NumberofMutationsDict,
                                       sample2IndelsSignature2NumberofMutationsDict,
                                       sample2DinucsSignature2NumberofMutationsDict,
                                       verbose,),
                                 callback=accumulate_apply_async_result))

            print('MONITOR %s %d len(jobs):%d' % (chrLong, simNum, len(jobs)), flush=True)
        ################################

        ##############################################################################
        # wait for all jobs to finish
        for job in jobs:
            if verbose: print('\tVerbose Replication Strand Bias Worker pid %s job.get():%s ' % (str(os.getpid()), job.get()))
        ##############################################################################

        ################################
        pool.close()
        pool.join()
        ################################

    ####################### ends July 22 2020   ###############################


    return simNum2Type2DecileIndex2AccumulatedNumberofMutationsDict, simNum2Sample2Type2DecileIndex2AccumulatedNumberofMutationsDict
##################################################################




##################################################################
# Augment wavelet_processed_df with numberofAttributableBases
def augment(genome,wavelet_processed_df,matrix_generator_path,verbose):

    ################################
    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)
    ################################

    #Augment for each chromosome
    chrBased_wavelet_processed_df_groups = wavelet_processed_df.groupby(CHROM)

    frames=[]

    ####################################################################################
    # tuple contains (chrLong,chrBased_wavelet_processed_df_group)
    def accumulate_apply_async_result(tuple):
        chrBased_augmented_df = tuple[1]
        frames.append(chrBased_augmented_df)
    ####################################################################################

    for chrLong, chrBased_wavelet_processed_df_group in chrBased_wavelet_processed_df_groups:
        chrShort=getChrShort(chrLong)
        chrbased_filename = chrShort + ".txt"
        chrbased_file_path = os.path.join(matrix_generator_path, 'references', 'chromosomes', 'tsb', genome, chrbased_filename)
        if os.path.exists(chrbased_file_path):
            pool.apply_async(addNumofAttributableBasesColumnForApplyAsync, (chrLong,chrBased_wavelet_processed_df_group,chrbased_file_path,verbose), callback=accumulate_apply_async_result)
    #JAN 21, 2020 ends

    ################################
    pool.close()
    pool.join()
    ################################

    augment_df = pd.concat(frames, ignore_index=True)

    return augment_df
##################################################################


##################################################################
def writeReplicationTimeData(outputDir,jobname,sample_based,decile_df_list,decileIndex2NumberofAttributableBasesDict,simNum2Type2DecileBasedAllChrAccumulatedCountDict,simNum2Sample2Type2DecileBasedAllChrAccumulatedCountDict):
    os.makedirs(os.path.join(outputDir, jobname, DATA, REPLICATIONTIME), exist_ok=True)

    ##############################################################
    for simNum in simNum2Type2DecileBasedAllChrAccumulatedCountDict:
        type2DecileBasedAllChrAccumulatedCountDict = simNum2Type2DecileBasedAllChrAccumulatedCountDict[simNum]
        for type in type2DecileBasedAllChrAccumulatedCountDict:
            decileBasedAllChrAccumulatedCountDict = type2DecileBasedAllChrAccumulatedCountDict[type]

            if (simNum==0):
                normalizedMutationDensityFilename = '%s_NormalizedMutationDensity.txt' %(type)
                numberofMutationsFilename = '%s_NumberofMutations.txt' %(type)
                numberofAttributabelBasesFilename = '%s_NumberofAttributableBases.txt' %(type)
            else:
                normalizedMutationDensityFilename = '%s_sim%d_NormalizedMutationDensity.txt' %(type,simNum)
                numberofMutationsFilename = '%s_sim%d_NumberofMutations.txt' %(type,simNum)
                numberofAttributabelBasesFilename = '%s_sim%d_NumberofAttributableBases.txt' %(type,simNum)

            if (type==AGGREGATEDSUBSTITUTIONS) or (type==AGGREGATEDINDELS) or (type==AGGREGATEDDINUCS) or (type==MICROHOMOLOGY) or (type==REPEAT):
                os.makedirs(os.path.join(outputDir, jobname, DATA, REPLICATIONTIME,type), exist_ok=True)
                normalizedMutationDensityFilePath = os.path.join(outputDir,jobname, DATA, REPLICATIONTIME,type,normalizedMutationDensityFilename)
                numberofMutationsFilePath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME,type,numberofMutationsFilename)
                numberofAttributabelBasesFilePath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME,type,numberofAttributabelBasesFilename)
            else:
                os.makedirs(os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, SIGNATUREBASED), exist_ok=True)
                normalizedMutationDensityFilePath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, SIGNATUREBASED,normalizedMutationDensityFilename)
                numberofMutationsFilePath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, SIGNATUREBASED,numberofMutationsFilename)
                numberofAttributabelBasesFilePath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME, SIGNATUREBASED,numberofAttributabelBasesFilename)

            # If decileBasedAllChrAccumulatedCountDict is not empty
            if (decileBasedAllChrAccumulatedCountDict):

                if (decile_df_list is not None):
                    numberofMutations, mutationDensityDict,numberofMutationsList, numberofAttributableBasesList = getMutationDensityDict(decile_df_list,decileBasedAllChrAccumulatedCountDict)
                elif (decileIndex2NumberofAttributableBasesDict is not None):
                    numberofMutations, mutationDensityDict,numberofMutationsList, numberofAttributableBasesList = getMutationDensityDictNewVersion(decileIndex2NumberofAttributableBasesDict, decileBasedAllChrAccumulatedCountDict)

                normalizedMutationDensityList = getNormalizedMutationDensityList(mutationDensityDict)

                #Normalized Mutation Density
                with open(normalizedMutationDensityFilePath, 'w') as file:
                    for normalizedMutationDensity in normalizedMutationDensityList:
                        file.write(str(normalizedMutationDensity) + ' ')
                    file.write('\n')

                #Number of Mutations
                with open(numberofMutationsFilePath, 'w') as file:
                    for numberofMutations in numberofMutationsList:
                        file.write(str(numberofMutations) + ' ')
                    file.write('\n')

                #Number of Attributable Bases
                with open(numberofAttributabelBasesFilePath, 'w') as file:
                    for numberofAttributabelBases in numberofAttributableBasesList:
                        file.write(str(numberofAttributabelBases) + ' ')
                    file.write('\n')
    ##############################################################


    ######################### Sample Based starts #####################
    if sample_based:
        for simNum in simNum2Sample2Type2DecileBasedAllChrAccumulatedCountDict:
            sample2Type2DecileBasedAllChrAccumulatedCountDict= simNum2Sample2Type2DecileBasedAllChrAccumulatedCountDict[simNum]
            for sample in sample2Type2DecileBasedAllChrAccumulatedCountDict:
                os.makedirs(os.path.join(outputDir, jobname, DATA, SAMPLES, sample, REPLICATIONTIME), exist_ok=True)

                for type in sample2Type2DecileBasedAllChrAccumulatedCountDict[sample]:
                    decileBasedAllChrAccumulatedCountDict = sample2Type2DecileBasedAllChrAccumulatedCountDict[sample][type]

                    if (simNum==0):
                        normalizedMutationDensityFilename = '%s_%s_NormalizedMutationDensity.txt' %(sample,type)
                    else:
                        normalizedMutationDensityFilename = '%s_%s_sim%d_NormalizedMutationDensity.txt' %(sample,type,simNum)

                    if (type == AGGREGATEDSUBSTITUTIONS) or (type == AGGREGATEDINDELS) or (type == AGGREGATEDDINUCS) or (type == MICROHOMOLOGY) or (type == REPEAT):
                        os.makedirs(os.path.join(outputDir, jobname, DATA,SAMPLES,sample, REPLICATIONTIME, type), exist_ok=True)
                        normalizedMutationDensityFilePath = os.path.join(outputDir, jobname, DATA, SAMPLES, sample,REPLICATIONTIME,type, normalizedMutationDensityFilename)
                    else:
                        os.makedirs(os.path.join(outputDir, jobname, DATA,SAMPLES,sample, REPLICATIONTIME, SIGNATUREBASED), exist_ok=True)
                        normalizedMutationDensityFilePath = os.path.join(outputDir, jobname, DATA, SAMPLES, sample,REPLICATIONTIME,SIGNATUREBASED, normalizedMutationDensityFilename)

                    #If decileBasedAllChrAccumulatedCountDict is not empty
                    if (decileBasedAllChrAccumulatedCountDict):
                        if (decile_df_list is not None):
                            numberofMutations, mutationDensityDict, numberofMutationsList, numberofAttributableBasesList = getMutationDensityDict(decile_df_list, decileBasedAllChrAccumulatedCountDict)
                        elif (decileIndex2NumberofAttributableBasesDict is not None):
                            numberofMutations, mutationDensityDict, numberofMutationsList, numberofAttributableBasesList = getMutationDensityDictNewVersion(decileIndex2NumberofAttributableBasesDict, decileBasedAllChrAccumulatedCountDict)
                        normalizedMutationDensityList = getNormalizedMutationDensityList(mutationDensityDict)

                        with open(normalizedMutationDensityFilePath, 'w') as file:
                            for normalizedMutationDensity in normalizedMutationDensityList:
                                file.write(str(normalizedMutationDensity) + ' ')
                            file.write('\n')
    ######################### Sample Based ends #######################

##################################################################


##################################################################
#main function
def replicationTimeAnalysis(computationType,sample_based,genome,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,job_tuples,repliseqDataFilename,subsSignature_cutoff_numberofmutations_averageprobability_df,indelsSignature_cutoff_numberofmutations_averageprobability_df,dinucsSignature_cutoff_numberofmutations_averageprobability_df,verbose,matrix_generator_path):

    print('\n#################################################################################')
    print('--- ReplicationTimeAnalysis starts')
    print('--- Replication Time Analyis Computation Type:%s' % (computationType))

    #########################################################################
    # Analysis Type can be
    # AggregatedSubstitutions: All in one
    # AggregatedIndels : All in one
    # AggregatedDinucs : All in one
    # IndelsBased : Microhomology, Repeat
    # SignatureBased: Subs Signatures Sig1, Sig2, ... and Indels Signatures  ID1, ID2, ..., ... and Dinucs Signatures  DBS1, DBS2, ...

    # We know the indels type we are interested in.
    # Microhomology indels --- len(indels) >= 3
    # Repeat indels --- len(indels) < 3
    #########################################################################

    ##########################################################################################
    if (sample_based):
        sample2NumberofSubsDict = getSample2NumberofSubsDict(outputDir,jobname)
        sample2NumberofIndelsDict = getSample2NumberofIndelsDict(outputDir,jobname)
        sample2NumberofDinucsDict = getDictionary(outputDir,jobname,Sample2NumberofDinucsDictFilename)

        sample2SubsSignature2NumberofMutationsDict = getSample2SubsSignature2NumberofMutationsDict(outputDir,jobname)
        sample2IndelsSignature2NumberofMutationsDict = getSample2IndelsSignature2NumberofMutationsDict(outputDir,jobname)
        sample2DinucsSignature2NumberofMutationsDict = getDictionary(outputDir,jobname,Sample2DinucsSignature2NumberofMutationsDictFilename)
    else:
        sample2NumberofSubsDict = {}
        sample2NumberofIndelsDict = {}
        sample2NumberofDinucsDict = {}

        sample2SubsSignature2NumberofMutationsDict = {}
        sample2IndelsSignature2NumberofMutationsDict = {}
        sample2DinucsSignature2NumberofMutationsDict = {}
    ##########################################################################################

    # Fill replication np arrays during runtime

    ###################################################################
    ############### Read RepliSeq Time data starts ####################
    ###################################################################
    # Fist decile_df in decile_df_list contains the intervals that are replicated the earliest.
    # Last decile_df in decile_df_list contains the intervals that are replicated the latest.
    # Please note that each decile_df contains intervals from all chroms (mixed chroms)
    #What is the type of deciles? Deciles is a list of dataframes.
    # if verbose: print('\tVerbose Worker pid %s READ Repliseq DATA STARTS  %s MB' % (str(os.getpid()), memory_usage()))
    if verbose: print('\tVerbose READ Repliseq DATA STARTS')
    #Whole genome is needed here
    #Formerly I was reading 2bit files using twobitreader
    #Now, formerly downloaded matrix generator reference genome is being used.
    chrNamesInReplicationTimeDataArray, decile_df_list = readRepliSeqTimeData(genome,chromNamesList,repliseqDataFilename,matrix_generator_path,verbose)
    # if verbose: print('\tVerbose Worker pid %s READ Repliseq DATA ENDS  %s MB' % (str(os.getpid()), memory_usage()))
    if verbose: print('\tVerbose READ Repliseq DATA ENDS')

    #Get chrBased grouped deciles
    chrBased_grouped_decile_df_list = []
    #Get this chrLong of each decile
    #The first decile is the earliest one with index 1
    #The last decile is the latest one with index 10
    #We provide these indexes later in the fillChrBasedDecileIndexArray function by enumerate
    for decile_df in decile_df_list:
        chrBased_grouped_decile_df = decile_df.groupby(CHROM)
        chrBased_grouped_decile_df_list.append(chrBased_grouped_decile_df)
    ###################################################################
    ############### Read RepliSeq Time data ends ######################
    ###################################################################


    #######################################################################################################
    ################################### Replication Time Data Analysis starts #############################
    #######################################################################################################
    simNum2Type2DecileBasedAllChrAccumulatedCountDict, simNum2Sample2Type2DecileBasedAllChrAccumulatedCountDict = calculateCountsForMutationsFillingReplicationTimeNPArrayRuntime(
                                                                                                                            computationType,
                                                                                                                            outputDir,
                                                                                                                            jobname,
                                                                                                                            numofSimulations,
                                                                                                                            job_tuples,
                                                                                                                            sample2NumberofSubsDict,
                                                                                                                            sample2NumberofIndelsDict,
                                                                                                                            sample2NumberofDinucsDict,
                                                                                                                            sample2SubsSignature2NumberofMutationsDict,
                                                                                                                            sample2IndelsSignature2NumberofMutationsDict,
                                                                                                                            sample2DinucsSignature2NumberofMutationsDict,
                                                                                                                            chromSizesDict,
                                                                                                                            chromNamesList,
                                                                                                                            chrBased_grouped_decile_df_list,
                                                                                                                            subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                                                                            indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                                                                            dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                                                                            sample_based,
                                                                                                                            verbose)

    writeReplicationTimeData(outputDir,jobname,sample_based,decile_df_list,None,simNum2Type2DecileBasedAllChrAccumulatedCountDict,simNum2Sample2Type2DecileBasedAllChrAccumulatedCountDict)
    #######################################################################################################
    ################################### Replication Time Data Analysis ends ###############################
    #######################################################################################################

    print('--- ReplicationTimeAnalysis ends')
    print('#################################################################################\n')

##################################################################