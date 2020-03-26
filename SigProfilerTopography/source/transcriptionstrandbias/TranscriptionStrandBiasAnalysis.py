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
import os

from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIPTIONSTRAND
from SigProfilerTopography.source.commons.TopographyCommons import SAMPLE
from SigProfilerTopography.source.commons.TopographyCommons import MUTATION

from SigProfilerTopography.source.commons.TopographyCommons import SUBS
from SigProfilerTopography.source.commons.TopographyCommons import INDELS
from SigProfilerTopography.source.commons.TopographyCommons import DINUCS

from SigProfilerTopography.source.commons.TopographyCommons import UNTRANSCRIBED_STRAND
from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIBED_STRAND
from SigProfilerTopography.source.commons.TopographyCommons import NONTRANSCRIBED_STRAND
from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIPTIONSTRANDBIAS

from SigProfilerTopography.source.commons.TopographyCommons import updateDictionaries_simulations_integrated
from SigProfilerTopography.source.commons.TopographyCommons import readChrBasedMutationsDF
from SigProfilerTopography.source.commons.TopographyCommons import accumulate_simulations_integrated
from SigProfilerTopography.source.commons.TopographyCommons import accumulate_simulations_integrated_for_each_tuple
from SigProfilerTopography.source.commons.TopographyCommons import writeDictionary

from SigProfilerTopography.source.commons.TopographyCommons import COMPUTATION_CHROMOSOMES_SEQUENTIAL_ALL_SIMULATIONS_PARALLEL
from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC

from SigProfilerTopography.source.commons.TopographyCommons import Type2TranscriptionStrand2CountDict_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Signature2MutationType2TranscriptionStrand2CountDict_Filename

from SigProfilerTopography.source.commons.TopographyCommons import Sample2Type2TranscriptionStrand2CountDict_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Type2Sample2TranscriptionStrand2CountDict_Filename
from SigProfilerTopography.source.commons.TopographyCommons import memory_usage


########################################################################
def searchMutationUsingTranscriptionStrandColumn_simulations_integrated(
        mutation_row,
        simNum2Type2TranscriptionStrand2CountDict,
        simNum2Sample2Type2TranscriptionStrand2CountDict,
        simNum2Type2Sample2TranscriptionStrand2CountDict,
        simNum2Signature2MutationType2TranscriptionStrand2CountDict,
        signature_cutoff_numberofmutations_averageprobability_df,
        type,
        sample_based):

    mutationType = None
    mutationTranscriptionStrand = mutation_row[TRANSCRIPTIONSTRAND]
    mutationSample = mutation_row[SAMPLE]

    if (type==SUBS):
        #e.g.: C>A
        mutationType = mutation_row[MUTATION]

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
#DEC 24, 2019
#Called by USING_APPLY_ASYNC
def searchMutationsForApplySync(chrBased_simBased_subs_df,
                                chrBased_simBased_indels_df,
                                chrBased_simBased_dinucs_df,
                                numofSimulations,
                                sample_based,
                                subsSignature_cutoff_numberofmutations_averageprobability_df,
                                indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                verbose):

    simNum2Type2TranscriptionStrand2CountDict = {}
    simNum2Sample2Type2TranscriptionStrand2CountDict = {}
    simNum2Type2Sample2TranscriptionStrand2CountDict = {}
    simNum2Signature2MutationType2TranscriptionStrand2CountDict = {}

    for simNum in range(0,numofSimulations+1):
        simNum2Type2TranscriptionStrand2CountDict[simNum]={}
        simNum2Sample2Type2TranscriptionStrand2CountDict[simNum]={}
        simNum2Type2Sample2TranscriptionStrand2CountDict[simNum]={}
        simNum2Signature2MutationType2TranscriptionStrand2CountDict[simNum]={}

    if ((chrBased_simBased_subs_df is not None) and (not chrBased_simBased_subs_df.empty)):
        if verbose: print('Worker pid %s SBS searchMutationUsingTranscriptionStrandColumn_simulations_integrated starts %s MB' % (str(os.getpid()), memory_usage()))
        chrBased_simBased_subs_df.apply(searchMutationUsingTranscriptionStrandColumn_simulations_integrated,
                                     simNum2Type2TranscriptionStrand2CountDict=simNum2Type2TranscriptionStrand2CountDict,
                                     simNum2Sample2Type2TranscriptionStrand2CountDict=simNum2Sample2Type2TranscriptionStrand2CountDict,
                                     simNum2Type2Sample2TranscriptionStrand2CountDict=simNum2Type2Sample2TranscriptionStrand2CountDict,
                                     simNum2Signature2MutationType2TranscriptionStrand2CountDict=simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                     signature_cutoff_numberofmutations_averageprobability_df=subsSignature_cutoff_numberofmutations_averageprobability_df,
                                     type=SUBS,
                                     sample_based=sample_based,
                                     axis=1)
        if verbose: print('Worker pid %s SBS searchMutationUsingTranscriptionStrandColumn_simulations_integrated ends %s MB' % (str(os.getpid()), memory_usage()))

    if ((chrBased_simBased_indels_df is not None) and (not chrBased_simBased_indels_df.empty)):
        if verbose: print('Worker pid %s ID searchMutationUsingTranscriptionStrandColumn_simulations_integrated starts %s MB' % (str(os.getpid()), memory_usage()))
        chrBased_simBased_indels_df.apply(searchMutationUsingTranscriptionStrandColumn_simulations_integrated,
                                       simNum2Type2TranscriptionStrand2CountDict=simNum2Type2TranscriptionStrand2CountDict,
                                       simNum2Sample2Type2TranscriptionStrand2CountDict=simNum2Sample2Type2TranscriptionStrand2CountDict,
                                       simNum2Type2Sample2TranscriptionStrand2CountDict=simNum2Type2Sample2TranscriptionStrand2CountDict,
                                       simNum2Signature2MutationType2TranscriptionStrand2CountDict=simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                       signature_cutoff_numberofmutations_averageprobability_df=indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                       type=INDELS,
                                       sample_based=sample_based,
                                       axis=1)
        if verbose: print('Worker pid %s ID searchMutationUsingTranscriptionStrandColumn_simulations_integrated ends %s MB' % (str(os.getpid()), memory_usage()))

    if ((chrBased_simBased_dinucs_df is not None) and (not chrBased_simBased_dinucs_df.empty)):
        if verbose: print('Worker pid %s DBS searchMutationUsingTranscriptionStrandColumn_simulations_integrated starts %s MB' % (str(os.getpid()), memory_usage()))
        chrBased_simBased_dinucs_df.apply(searchMutationUsingTranscriptionStrandColumn_simulations_integrated,
                                       simNum2Type2TranscriptionStrand2CountDict=simNum2Type2TranscriptionStrand2CountDict,
                                       simNum2Sample2Type2TranscriptionStrand2CountDict=simNum2Sample2Type2TranscriptionStrand2CountDict,
                                       simNum2Type2Sample2TranscriptionStrand2CountDict=simNum2Type2Sample2TranscriptionStrand2CountDict,
                                       simNum2Signature2MutationType2TranscriptionStrand2CountDict=simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                       signature_cutoff_numberofmutations_averageprobability_df=dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                       type=DINUCS,
                                       sample_based=sample_based,
                                       axis=1)
        if verbose: print('Worker pid %s DBS searchMutationUsingTranscriptionStrandColumn_simulations_integrated ends %s MB' % (str(os.getpid()), memory_usage()))

    return (simNum2Type2TranscriptionStrand2CountDict,
            simNum2Sample2Type2TranscriptionStrand2CountDict,
            simNum2Type2Sample2TranscriptionStrand2CountDict,
            simNum2Signature2MutationType2TranscriptionStrand2CountDict)
########################################################################


########################################################################
#Called from COMPUTATION_CHROMOSOMES_SEQUENTIAL_ALL_SIMULATIONS_PARALLEL
def searchMutations(inputList):
    chrBased_subs_split_df = inputList[0]
    chrBased_indels_split_df = inputList[1]
    chrBased_dinucs_split_df = inputList[2]
    numofSimulations = inputList[3]
    sample_based = inputList[4]
    subsSignature_cutoff_numberofmutations_averageprobability_df=inputList[5]
    indelsSignature_cutoff_numberofmutations_averageprobability_df = inputList[6]
    dinucsSignature_cutoff_numberofmutations_averageprobability_df = inputList[7]

    simNum2Type2TranscriptionStrand2CountDict = {}
    simNum2Sample2Type2TranscriptionStrand2CountDict = {}
    simNum2Type2Sample2TranscriptionStrand2CountDict = {}
    simNum2Signature2MutationType2TranscriptionStrand2CountDict = {}

    #Initialization
    for simNum in range(0,numofSimulations+1):
        simNum2Type2TranscriptionStrand2CountDict[simNum]={}
        simNum2Sample2Type2TranscriptionStrand2CountDict[simNum]={}
        simNum2Type2Sample2TranscriptionStrand2CountDict[simNum]={}
        simNum2Signature2MutationType2TranscriptionStrand2CountDict[simNum]={}

    if ((chrBased_subs_split_df is not None) and (not chrBased_subs_split_df.empty)):
        chrBased_subs_split_df.apply(searchMutationUsingTranscriptionStrandColumn_simulations_integrated,
                                     simNum2Type2TranscriptionStrand2CountDict=simNum2Type2TranscriptionStrand2CountDict,
                                     simNum2Sample2Type2TranscriptionStrand2CountDict=simNum2Sample2Type2TranscriptionStrand2CountDict,
                                     simNum2Type2Sample2TranscriptionStrand2CountDict=simNum2Type2Sample2TranscriptionStrand2CountDict,
                                     simNum2Signature2MutationType2TranscriptionStrand2CountDict=simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                     signature_cutoff_numberofmutations_averageprobability_df=subsSignature_cutoff_numberofmutations_averageprobability_df,
                                     type=SUBS,
                                     sample_based=sample_based,
                                     axis=1)

    if ((chrBased_indels_split_df is not None) and (not chrBased_indels_split_df.empty)):
        chrBased_indels_split_df.apply(searchMutationUsingTranscriptionStrandColumn_simulations_integrated,
                                       simNum2Type2TranscriptionStrand2CountDict=simNum2Type2TranscriptionStrand2CountDict,
                                       simNum2Sample2Type2TranscriptionStrand2CountDict=simNum2Sample2Type2TranscriptionStrand2CountDict,
                                       simNum2Type2Sample2TranscriptionStrand2CountDict=simNum2Type2Sample2TranscriptionStrand2CountDict,
                                       simNum2Signature2MutationType2TranscriptionStrand2CountDict=simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                       signature_cutoff_numberofmutations_averageprobability_df=indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                       type=INDELS,
                                       sample_based=sample_based,
                                       axis=1)

    if ((chrBased_dinucs_split_df is not None) and (not chrBased_dinucs_split_df.empty)):
        chrBased_dinucs_split_df.apply(searchMutationUsingTranscriptionStrandColumn_simulations_integrated,
                                       simNum2Type2TranscriptionStrand2CountDict=simNum2Type2TranscriptionStrand2CountDict,
                                       simNum2Sample2Type2TranscriptionStrand2CountDict=simNum2Sample2Type2TranscriptionStrand2CountDict,
                                       simNum2Type2Sample2TranscriptionStrand2CountDict=simNum2Type2Sample2TranscriptionStrand2CountDict,
                                       simNum2Signature2MutationType2TranscriptionStrand2CountDict=simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                       signature_cutoff_numberofmutations_averageprobability_df=dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                       type=DINUCS,
                                       sample_based=sample_based,
                                       axis=1)

    return (simNum2Type2TranscriptionStrand2CountDict,
            simNum2Sample2Type2TranscriptionStrand2CountDict,
            simNum2Type2Sample2TranscriptionStrand2CountDict,
            simNum2Signature2MutationType2TranscriptionStrand2CountDict)
########################################################################


########################################################################
#main function
def transcriptionStrandBiasAnalysis(computationType,sample_based,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,subsSignature_cutoff_numberofmutations_averageprobability_df,indelsSignature_cutoff_numberofmutations_averageprobability_df,dinucsSignature_cutoff_numberofmutations_averageprobability_df,verbose):

    print('\n#################################################################################')
    print('--- TranscriptionStrandBias Analysis starts')

    #############################################
    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)
    #############################################

    ##################### Read Transcripts starts ######################
    #NCBI has the long chromosome names such as: chr1, chr2, chr3, chr4, ... , chr21, chr22, chrX, chrY, chrMT
    # transcriptsSource = NCBI
    # GRCh37_hg19_Transcripts_df = readTranscriptsNCBI()

    # Let's make SigProfiler use the same transcripts file
    #Ensembl has the short chromosome names such as: 1,2,3,4, ... ,21, 22, X, Y, MT
    # transcriptsSource = ENSEMBL
    # transcripts_df = readTrancriptsENSEMBL(genome)

    # print('transcripts_df.shape')
    # print(transcripts_df.shape)
    ##################### Read Transcripts ends ########################

    strandBias = TRANSCRIPTIONSTRANDBIAS

    #Accumulate chrBased Results
    accumulatedAllChromosomesType2TranscriptionStrand2CountDict = {}
    accumulatedAllChromosomesSample2Type2TranscriptionStrand2CountDict = {}
    accumulatedAllChromosomesType2Sample2TranscriptionStrand2CountDict = {}
    accumulatedAllChromosomesSignature2MutationType2TranscriptionStrand2CountDict = {}

    if (computationType == COMPUTATION_CHROMOSOMES_SEQUENTIAL_ALL_SIMULATIONS_PARALLEL):

        ####################################################################################################
        for chrLong in chromNamesList:
            chromSize = chromSizesDict[chrLong]
            poolInputList = []

            #You need to initialize to None so that you don't use former for loop values accidentally

            ####################################################################
            for simNum in range(0,numofSimulations+1):
                inputList = []
                chrBased_subs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, SUBS, simNum)
                chrBased_indels_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, INDELS, simNum)
                chrBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, DINUCS, simNum)
                inputList.append(chrBased_subs_df)  # each time different split
                inputList.append(chrBased_indels_df)
                inputList.append(chrBased_dinucs_df)
                inputList.append(numofSimulations)
                inputList.append(sample_based)
                inputList.append(subsSignature_cutoff_numberofmutations_averageprobability_df)
                inputList.append(indelsSignature_cutoff_numberofmutations_averageprobability_df)
                inputList.append(dinucsSignature_cutoff_numberofmutations_averageprobability_df)
                poolInputList.append(inputList)
            ####################################################################

            listofTuples = pool.map(searchMutations,poolInputList)

            accumulate_simulations_integrated(listofTuples,
                                              accumulatedAllChromosomesType2TranscriptionStrand2CountDict,
                                              accumulatedAllChromosomesSample2Type2TranscriptionStrand2CountDict,
                                              accumulatedAllChromosomesType2Sample2TranscriptionStrand2CountDict,
                                              accumulatedAllChromosomesSignature2MutationType2TranscriptionStrand2CountDict)
        ####################################################################################################

    #DEC 24, 2019
    elif (computationType==USING_APPLY_ASYNC):

        ####################################################################################################
        for chrLong in chromNamesList:
            chromSize = chromSizesDict[chrLong]

            ####################################################################
            def accumulate_apply_async_result(result_tuple):
                chrBased_SimNum2Type2Strand2CountDict=result_tuple[0]
                chrBased_SimNum2Sample2Type2Strand2CountDict=result_tuple[1]
                chrBased_SimNum2Type2Sample2Strand2CountDict=result_tuple[2]
                chrBased_SimNum2Signature2MutationType2Strand2CountDict=result_tuple[3]

                if verbose: print('Worker pid %s Accumulate Transcription Strand Bias %s MB' % (str(os.getpid()), memory_usage()))

                accumulate_simulations_integrated_for_each_tuple(
                            chrBased_SimNum2Type2Strand2CountDict,
                            chrBased_SimNum2Sample2Type2Strand2CountDict,
                            chrBased_SimNum2Type2Sample2Strand2CountDict,
                            chrBased_SimNum2Signature2MutationType2Strand2CountDict,
                            accumulatedAllChromosomesType2TranscriptionStrand2CountDict,
                            accumulatedAllChromosomesSample2Type2TranscriptionStrand2CountDict,
                            accumulatedAllChromosomesType2Sample2TranscriptionStrand2CountDict,
                            accumulatedAllChromosomesSignature2MutationType2TranscriptionStrand2CountDict)
            ####################################################################

            ####################################################################
            for simNum in range(0,numofSimulations+1):
                chrBased_simBased_subs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, SUBS, simNum)
                chrBased_simBased_indels_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, INDELS, simNum)
                chrBased_simBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, DINUCS, simNum)
                pool.apply_async(searchMutationsForApplySync, (chrBased_simBased_subs_df,chrBased_simBased_indels_df,chrBased_simBased_dinucs_df,numofSimulations,sample_based,subsSignature_cutoff_numberofmutations_averageprobability_df,indelsSignature_cutoff_numberofmutations_averageprobability_df,dinucsSignature_cutoff_numberofmutations_averageprobability_df,verbose),callback=accumulate_apply_async_result)
            ####################################################################


        ####################################################################################################

    ################################
    pool.close()
    pool.join()
    ################################

    #################################################################################################################
    ##########################################      Output starts      ##############################################
    #################################################################################################################
    #############################################################################
    writeDictionary(accumulatedAllChromosomesType2TranscriptionStrand2CountDict,outputDir,jobname,Type2TranscriptionStrand2CountDict_Filename,strandBias,None)
    writeDictionary(accumulatedAllChromosomesSignature2MutationType2TranscriptionStrand2CountDict,outputDir,jobname,Signature2MutationType2TranscriptionStrand2CountDict_Filename,strandBias,None)

    if sample_based:
        writeDictionary(accumulatedAllChromosomesSample2Type2TranscriptionStrand2CountDict,outputDir,jobname,Sample2Type2TranscriptionStrand2CountDict_Filename,strandBias,None)
        writeDictionary(accumulatedAllChromosomesType2Sample2TranscriptionStrand2CountDict, outputDir, jobname,Type2Sample2TranscriptionStrand2CountDict_Filename, strandBias, None)
    #################################################################################################################
    ##########################################      Output ends      ################################################
    #################################################################################################################

    print('--- TranscriptionStrandBias Analysis ends')
    print('#################################################################################\n')

########################################################################