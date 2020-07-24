# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018-2020 Burcak Otlu

#############################################################
# This version use np.arrays
# Right now transcription strand bias analysis works for single point mutations and signatures.
# Constraints, Thresholds
# Please note that for sample based transcription strand bias analysis
# We consider samples with at least 1000 mutations both on transcribed and non-transcribed strands.
#############################################################

#############################################################
# What is transcription strand bias?
# It is the ratio of = (number of mutations on transcribed strand) / (number of mutations on un-transcribed strand)
#############################################################

import multiprocessing
import numpy as np
import pandas as pd
import os
import math

from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIPTIONSTRAND
from SigProfilerTopography.source.commons.TopographyCommons import SAMPLE
from SigProfilerTopography.source.commons.TopographyCommons import MUTATION

from SigProfilerTopography.source.commons.TopographyCommons import TYPE

from SigProfilerTopography.source.commons.TopographyCommons import SUBS
from SigProfilerTopography.source.commons.TopographyCommons import INDELS
from SigProfilerTopography.source.commons.TopographyCommons import DINUCS

from SigProfilerTopography.source.commons.TopographyCommons import UNTRANSCRIBED_STRAND
from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIBED_STRAND
from SigProfilerTopography.source.commons.TopographyCommons import NONTRANSCRIBED_STRAND
from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIPTIONSTRANDBIAS

from SigProfilerTopography.source.commons.TopographyCommons import updateDictionaries_simulations_integrated
from SigProfilerTopography.source.commons.TopographyCommons import updateDictionaries_simulations_integrated_for_list_comprehension

from SigProfilerTopography.source.commons.TopographyCommons import readChrBasedMutationsDF
from SigProfilerTopography.source.commons.TopographyCommons import accumulate_simulations_integrated_for_each_tuple

from SigProfilerTopography.source.commons.TopographyCommons import writeDictionary
from SigProfilerTopography.source.commons.TopographyCommons import write_signature_mutation_type_strand_bias_dictionary_as_dataframe
from SigProfilerTopography.source.commons.TopographyCommons import write_type_strand_bias_dictionary_as_dataframe


from SigProfilerTopography.source.commons.TopographyCommons import get_chrBased_simBased_combined_df_split
from SigProfilerTopography.source.commons.TopographyCommons import get_chrBased_simBased_combined_df

from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM
from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT

from SigProfilerTopography.source.commons.TopographyCommons import Type2TranscriptionStrand2CountDict_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Signature2MutationType2TranscriptionStrand2CountDict_Filename

from SigProfilerTopography.source.commons.TopographyCommons import Sample2Type2TranscriptionStrand2CountDict_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Type2Sample2TranscriptionStrand2CountDict_Filename
from SigProfilerTopography.source.commons.TopographyCommons import memory_usage
from SigProfilerTopography.source.commons.TopographyCommons import NUMBER_OF_MUTATIONS_IN_EACH_SPLIT
from SigProfilerTopography.source.commons.TopographyCommons import MAXIMUM_NUMBER_JOBS_IN_THE_POOL_AT_ONCE


########################################################################
#April 24, 2020
def searchAllMutationUsingTranscriptionStrandColumn_using_list_comprehension(mutation_row,
                                                                                  simNum2Type2TranscriptionStrand2CountDict,
                                                                                  simNum2Sample2Type2TranscriptionStrand2CountDict,
                                                                                  simNum2Type2Sample2TranscriptionStrand2CountDict,
                                                                                  simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                                                                  subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                                  indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                                  dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                                  sample_based,df_columns):

    mutationType = None

    indexofTranscriptionStrand = df_columns.index(TRANSCRIPTIONSTRAND)
    mutationTranscriptionStrand = mutation_row[indexofTranscriptionStrand]

    indexofSample = df_columns.index(SAMPLE)
    mutationSample = mutation_row[indexofSample]

    ##########################################
    indexofType = df_columns.index(TYPE)
    my_type=mutation_row[indexofType]

    if (my_type==SUBS):
        signature_cutoff_numberofmutations_averageprobability_df=subsSignature_cutoff_numberofmutations_averageprobability_df
    elif (my_type==INDELS):
        signature_cutoff_numberofmutations_averageprobability_df=indelsSignature_cutoff_numberofmutations_averageprobability_df
    elif (my_type==DINUCS):
        signature_cutoff_numberofmutations_averageprobability_df=dinucsSignature_cutoff_numberofmutations_averageprobability_df

    if (my_type==SUBS):
        #e.g.: C>A
        indexofMutation = df_columns.index(MUTATION)
        mutationType = mutation_row[indexofMutation]
    ##########################################

    #Values on TranscriptionStrand column
    # N --> Non-transcribed
    # T --> Transcribed
    # U --> Untranscribed
    # Q --> Question Not known

    if (mutationTranscriptionStrand == 'U'):
        updateDictionaries_simulations_integrated_for_list_comprehension(mutation_row,
                                mutationType,
                                mutationSample,
                                sample_based,
                                simNum2Type2TranscriptionStrand2CountDict,
                                simNum2Sample2Type2TranscriptionStrand2CountDict,
                                simNum2Type2Sample2TranscriptionStrand2CountDict,
                                simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                UNTRANSCRIBED_STRAND,
                                signature_cutoff_numberofmutations_averageprobability_df,
                                df_columns)

    elif (mutationTranscriptionStrand == 'T'):
        updateDictionaries_simulations_integrated_for_list_comprehension(mutation_row,
                                mutationType,
                                mutationSample,
                                sample_based,
                                simNum2Type2TranscriptionStrand2CountDict,
                                simNum2Sample2Type2TranscriptionStrand2CountDict,
                                simNum2Type2Sample2TranscriptionStrand2CountDict,
                                simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                TRANSCRIBED_STRAND,
                                signature_cutoff_numberofmutations_averageprobability_df,
                                df_columns)
    elif (mutationTranscriptionStrand == 'B'):
        updateDictionaries_simulations_integrated_for_list_comprehension(mutation_row,
                                mutationType,
                                mutationSample,
                                sample_based,
                               simNum2Type2TranscriptionStrand2CountDict,
                               simNum2Sample2Type2TranscriptionStrand2CountDict,
                               simNum2Type2Sample2TranscriptionStrand2CountDict,
                               simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                               UNTRANSCRIBED_STRAND,
                               signature_cutoff_numberofmutations_averageprobability_df,
                               df_columns)
        updateDictionaries_simulations_integrated_for_list_comprehension(mutation_row,
                                mutationType,
                                mutationSample,
                                sample_based,
                                simNum2Type2TranscriptionStrand2CountDict,
                                simNum2Sample2Type2TranscriptionStrand2CountDict,
                                simNum2Type2Sample2TranscriptionStrand2CountDict,
                                simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                TRANSCRIBED_STRAND,
                                signature_cutoff_numberofmutations_averageprobability_df,
                                df_columns)
    elif (mutationTranscriptionStrand == 'N'):
        updateDictionaries_simulations_integrated_for_list_comprehension(mutation_row,
                                mutationType,
                                mutationSample,
                                sample_based,
                                simNum2Type2TranscriptionStrand2CountDict,
                                simNum2Sample2Type2TranscriptionStrand2CountDict,
                                simNum2Type2Sample2TranscriptionStrand2CountDict,
                                simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                NONTRANSCRIBED_STRAND,
                                signature_cutoff_numberofmutations_averageprobability_df,
                                df_columns)

########################################################################

########################################################################
#April 5, 2020
def searchAllMutationUsingTranscriptionStrandColumn_using_apply(
        mutation_row,
        simNum2Type2TranscriptionStrand2CountDict,
        simNum2Sample2Type2TranscriptionStrand2CountDict,
        simNum2Type2Sample2TranscriptionStrand2CountDict,
        simNum2Signature2MutationType2TranscriptionStrand2CountDict,
        subsSignature_cutoff_numberofmutations_averageprobability_df,
        indelsSignature_cutoff_numberofmutations_averageprobability_df,
        dinucsSignature_cutoff_numberofmutations_averageprobability_df,
        sample_based):

    mutationType = None
    mutationTranscriptionStrand = mutation_row[TRANSCRIPTIONSTRAND]
    mutationSample = mutation_row[SAMPLE]

    ##########################################
    type=mutation_row[TYPE]

    if (type==SUBS):
        signature_cutoff_numberofmutations_averageprobability_df=subsSignature_cutoff_numberofmutations_averageprobability_df
    elif (type==INDELS):
        signature_cutoff_numberofmutations_averageprobability_df=indelsSignature_cutoff_numberofmutations_averageprobability_df
    elif (type==DINUCS):
        signature_cutoff_numberofmutations_averageprobability_df=dinucsSignature_cutoff_numberofmutations_averageprobability_df

    if (type==SUBS):
        #e.g.: C>A
        mutationType = mutation_row[MUTATION]
    ##########################################

    #Values on TranscriptionStrand column
    # N --> Non-transcribed
    # T --> Transcribed
    # U --> Untranscribed
    # Q --> Question Not known

    if (mutationTranscriptionStrand == 'U'):
        updateDictionaries_simulations_integrated(mutation_row,
                                mutationType,
                                mutationSample,
                                sample_based,
                                simNum2Type2TranscriptionStrand2CountDict,
                                simNum2Sample2Type2TranscriptionStrand2CountDict,
                                simNum2Type2Sample2TranscriptionStrand2CountDict,
                                simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                UNTRANSCRIBED_STRAND,
                                signature_cutoff_numberofmutations_averageprobability_df)

    elif (mutationTranscriptionStrand == 'T'):
        updateDictionaries_simulations_integrated(mutation_row,
                                mutationType,
                                mutationSample,
                                sample_based,
                                simNum2Type2TranscriptionStrand2CountDict,
                                simNum2Sample2Type2TranscriptionStrand2CountDict,
                                simNum2Type2Sample2TranscriptionStrand2CountDict,
                                simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                TRANSCRIBED_STRAND,
                                signature_cutoff_numberofmutations_averageprobability_df)
    elif (mutationTranscriptionStrand == 'B'):
        updateDictionaries_simulations_integrated(mutation_row,
                                mutationType,
                                mutationSample,
                                sample_based,
                               simNum2Type2TranscriptionStrand2CountDict,
                               simNum2Sample2Type2TranscriptionStrand2CountDict,
                               simNum2Type2Sample2TranscriptionStrand2CountDict,
                               simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                               UNTRANSCRIBED_STRAND,
                               signature_cutoff_numberofmutations_averageprobability_df)
        updateDictionaries_simulations_integrated(mutation_row,
                                mutationType,
                                mutationSample,
                                sample_based,
                                simNum2Type2TranscriptionStrand2CountDict,
                                simNum2Sample2Type2TranscriptionStrand2CountDict,
                                simNum2Type2Sample2TranscriptionStrand2CountDict,
                                simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                TRANSCRIBED_STRAND,
                                signature_cutoff_numberofmutations_averageprobability_df)
    elif (mutationTranscriptionStrand == 'N'):
        updateDictionaries_simulations_integrated(mutation_row,
                                mutationType,
                                mutationSample,
                                sample_based,
                                simNum2Type2TranscriptionStrand2CountDict,
                                simNum2Sample2Type2TranscriptionStrand2CountDict,
                                simNum2Type2Sample2TranscriptionStrand2CountDict,
                                simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                NONTRANSCRIBED_STRAND,
                                signature_cutoff_numberofmutations_averageprobability_df)
########################################################################


########################################################################
# April 30, 2020
def searchAllMutations(chrBased_simBased_combined_df,sample_based,subsSignature_cutoff_numberofmutations_averageprobability_df,indelsSignature_cutoff_numberofmutations_averageprobability_df,dinucsSignature_cutoff_numberofmutations_averageprobability_df,verbose):

    ################################################################################
    simNum2Type2TranscriptionStrand2CountDict = {}
    simNum2Sample2Type2TranscriptionStrand2CountDict = {}
    simNum2Type2Sample2TranscriptionStrand2CountDict = {}
    simNum2Signature2MutationType2TranscriptionStrand2CountDict = {}
    ################################################################################

    ################################################################################
    if ((chrBased_simBased_combined_df is not None) and (not chrBased_simBased_combined_df.empty)):
        if verbose: print('Worker pid %s searchMutationUsingTranscriptionStrandColumn_simulations_integrated starts %s MB' % (str(os.getpid()), memory_usage()))

        # #####################################################################################
        # #Using Apply
        # chrBased_simBased_combined_df_split.apply(searchAllMutationUsingTranscriptionStrandColumn_using_apply,
        #                              simNum2Type2TranscriptionStrand2CountDict=simNum2Type2TranscriptionStrand2CountDict,
        #                              simNum2Sample2Type2TranscriptionStrand2CountDict=simNum2Sample2Type2TranscriptionStrand2CountDict,
        #                              simNum2Type2Sample2TranscriptionStrand2CountDict=simNum2Type2Sample2TranscriptionStrand2CountDict,
        #                              simNum2Signature2MutationType2TranscriptionStrand2CountDict=simNum2Signature2MutationType2TranscriptionStrand2CountDict,
        #                             subsSignature_cutoff_numberofmutations_averageprobability_df=subsSignature_cutoff_numberofmutations_averageprobability_df,
        #                             indelsSignature_cutoff_numberofmutations_averageprobability_df=indelsSignature_cutoff_numberofmutations_averageprobability_df,
        #                             dinucsSignature_cutoff_numberofmutations_averageprobability_df=dinucsSignature_cutoff_numberofmutations_averageprobability_df,
        #                              sample_based=sample_based,
        #                              axis=1)
        # #####################################################################################

        #####################################################################################
        # Using list comprehension
        df_columns = list(chrBased_simBased_combined_df.columns.values)
        [searchAllMutationUsingTranscriptionStrandColumn_using_list_comprehension(mutation_row,
                                                                                  simNum2Type2TranscriptionStrand2CountDict,
                                                                                  simNum2Sample2Type2TranscriptionStrand2CountDict,
                                                                                  simNum2Type2Sample2TranscriptionStrand2CountDict,
                                                                                  simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                                                                  subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                                  indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                                  dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                                  sample_based,
                                                                                  df_columns) for mutation_row in chrBased_simBased_combined_df.values]
        #####################################################################################

        if verbose: print('Worker pid %s searchMutationUsingTranscriptionStrandColumn_simulations_integrated ends %s MB' % (str(os.getpid()), memory_usage()))
    ################################################################################

    return (simNum2Type2TranscriptionStrand2CountDict,
            simNum2Sample2Type2TranscriptionStrand2CountDict,
            simNum2Type2Sample2TranscriptionStrand2CountDict,
            simNum2Signature2MutationType2TranscriptionStrand2CountDict)
########################################################################


########################################################################
def searchAllMutations_simbased_chrombased_splitbased(outputDir,
                                           jobname,
                                           chrLong,
                                           simNum,
                                           splitIndex,
                                           sample_based,
                                           subsSignature_cutoff_numberofmutations_averageprobability_df,
                                           indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                           dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                           verbose):

    chrBased_simBased_combined_df_split = get_chrBased_simBased_combined_df_split(outputDir, jobname, chrLong, simNum, splitIndex)

    return searchAllMutations(chrBased_simBased_combined_df_split,
                              sample_based,
                              subsSignature_cutoff_numberofmutations_averageprobability_df,
                              indelsSignature_cutoff_numberofmutations_averageprobability_df,
                              dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                              verbose)
########################################################################

########################################################################
def searchAllMutations_simbased_chrombased(outputDir,
                                       jobname,
                                       chrLong,
                                       simNum,
                                       sample_based,
                                       subsSignature_cutoff_numberofmutations_averageprobability_df,
                                       indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                       dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                       verbose):

    chrBased_simBased_combined_df = get_chrBased_simBased_combined_df(outputDir, jobname, chrLong, simNum)

    return  searchAllMutations(chrBased_simBased_combined_df,
                               sample_based,
                               subsSignature_cutoff_numberofmutations_averageprobability_df,
                               indelsSignature_cutoff_numberofmutations_averageprobability_df,
                               dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                               verbose)
########################################################################


########################################################################
#main function
def transcriptionStrandBiasAnalysis(computationType,sample_based,chromNamesList,outputDir,jobname,numofSimulations,job_tuples,subsSignature_cutoff_numberofmutations_averageprobability_df,indelsSignature_cutoff_numberofmutations_averageprobability_df,dinucsSignature_cutoff_numberofmutations_averageprobability_df,verbose):

    print('\n#################################################################################')
    print('--- TranscriptionStrandBias Analysis starts')

    strandBias = TRANSCRIPTIONSTRANDBIAS

    ###############################################################################
    #Accumulate Results
    simNum2Type2TranscriptionStrand2AccumulatedCountDict = {}
    simNum2Sample2Type2TranscriptionStrand2AccumulatedCountDict = {}
    simNum2Type2Sample2TranscriptionStrand2AccumulatedCountDict = {}
    simNum2Signature2MutationType2TranscriptionStrand2AccumulatedCountDict = {}
    ###############################################################################

    #########################################################################################
    def accumulate_apply_async_result(result_tuple):
        simNum2Type2Strand2CountDict = result_tuple[0]
        simNum2Sample2Type2Strand2CountDict = result_tuple[1]
        simNum2Type2Sample2Strand2CountDict = result_tuple[2]
        simNum2Signature2MutationType2Strand2CountDict = result_tuple[3]

        print('MONITOR ACCUMULATE', flush=True)

        accumulate_simulations_integrated_for_each_tuple(
            simNum2Type2Strand2CountDict,
            simNum2Sample2Type2Strand2CountDict,
            simNum2Type2Sample2Strand2CountDict,
            simNum2Signature2MutationType2Strand2CountDict,
            simNum2Type2TranscriptionStrand2AccumulatedCountDict,
            simNum2Sample2Type2TranscriptionStrand2AccumulatedCountDict,
            simNum2Type2Sample2TranscriptionStrand2AccumulatedCountDict,
            simNum2Signature2MutationType2TranscriptionStrand2AccumulatedCountDict)
    #########################################################################################

    ###############################################################################
    # April 30, 2020
    # Read the chrom based sim based mutations data in the worker process
    if (computationType==USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM):

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
            jobs.append(pool.apply_async(searchAllMutations_simbased_chrombased,
                                    args=(outputDir,
                                          jobname,
                                          chrLong,
                                          simNum,
                                          sample_based,
                                          subsSignature_cutoff_numberofmutations_averageprobability_df,
                                          indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                          dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                          verbose,),
                                    callback=accumulate_apply_async_result))
            print('MONITOR %s simNum:%d len(jobs):%d' % (chrLong, simNum, len(jobs)), flush=True)
        ################################################################################

        ##############################################################################
        # wait for all jobs to finish
        for job in jobs:
            if verbose: print('\tVerbose Transcription Strand Bias Worker pid %s job.get():%s ' % (str(os.getpid()), job.get()))
        ##############################################################################

        ################################
        pool.close()
        pool.join()
        ################################

    ###############################################################################

    ###############################################################################
    if (computationType==USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT):

        ################################
        numofProcesses = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=numofProcesses)
        ################################

        ################################
        jobs = []
        ################################

        for chrLong, simNum, splitIndex in job_tuples:
            jobs.append(pool.apply_async(searchAllMutations_simbased_chrombased_splitbased,
                                    args=(outputDir,
                                          jobname,
                                          chrLong,
                                          simNum,
                                          splitIndex,
                                          sample_based,
                                          subsSignature_cutoff_numberofmutations_averageprobability_df,
                                          indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                          dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                          verbose,),
                                    callback=accumulate_apply_async_result))
            print('MONITOR %s simNum:%d len(jobs):%d' % (chrLong, simNum, len(jobs)), flush=True)
        ################################################################################

        ##############################################################################
        # wait for all jobs to finish
        for job in jobs:
            if verbose: print('\tVerbose Transcription Strand Bias Worker pid %s job.get():%s ' % (str(os.getpid()), job.get()))
        ##############################################################################

        ################################
        pool.close()
        pool.join()
        ################################

    ###############################################################################

    #################################################################################################################
    ##########################################      Output starts      ##############################################
    #################################################################################################################
    transcription_atrands = [TRANSCRIBED_STRAND, UNTRANSCRIBED_STRAND, NONTRANSCRIBED_STRAND]

    write_signature_mutation_type_strand_bias_dictionary_as_dataframe(simNum2Signature2MutationType2TranscriptionStrand2AccumulatedCountDict,
                                                                      strandBias,
                                                                      transcription_atrands,
                                                                      outputDir,
                                                                      jobname)


    write_type_strand_bias_dictionary_as_dataframe(simNum2Type2TranscriptionStrand2AccumulatedCountDict,
                                                   strandBias,
                                                   transcription_atrands,
                                                   outputDir,
                                                   jobname)


    if sample_based:
        writeDictionary(simNum2Sample2Type2TranscriptionStrand2AccumulatedCountDict,outputDir,jobname,Sample2Type2TranscriptionStrand2CountDict_Filename,strandBias,None)
        writeDictionary(simNum2Type2Sample2TranscriptionStrand2AccumulatedCountDict, outputDir, jobname,Type2Sample2TranscriptionStrand2CountDict_Filename, strandBias, None)
        #TODO write as dataframe if sample_based will be maintained
    #################################################################################################################
    ##########################################      Output ends      ################################################
    #################################################################################################################

    print('--- TranscriptionStrandBias Analysis ends')
    print('#################################################################################\n')

########################################################################