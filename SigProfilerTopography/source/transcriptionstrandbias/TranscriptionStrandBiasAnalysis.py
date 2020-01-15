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

from SigProfilerTopography.source.commons.TopographyCommons import START
from SigProfilerTopography.source.commons.TopographyCommons import PYRAMIDINESTRAND
from SigProfilerTopography.source.commons.TopographyCommons import LENGTH
from SigProfilerTopography.source.commons.TopographyCommons import REF
from SigProfilerTopography.source.commons.TopographyCommons import ALT

from SigProfilerTopography.source.commons.TopographyCommons import allPyrimidine
from SigProfilerTopography.source.commons.TopographyCommons import allPurine

from SigProfilerTopography.source.commons.TopographyCommons import SUBS
from SigProfilerTopography.source.commons.TopographyCommons import INDELS
from SigProfilerTopography.source.commons.TopographyCommons import DINUCS

from SigProfilerTopography.source.commons.TopographyCommons import UNTRANSCRIBED_STRAND
from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIBED_STRAND
from SigProfilerTopography.source.commons.TopographyCommons import NONTRANSCRIBED_STRAND
from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIPTIONSTRANDBIAS

from SigProfilerTopography.source.commons.TopographyCommons import updateDictionaries_simulations_integrated
from SigProfilerTopography.source.commons.TopographyCommons import readTrancriptsENSEMBL
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
# Fill chrBased array once for each chromosome
# int8	Byte (-128 to 127)
# We use two arrays; one for positive strand and one for negative strand.
# 0 --> no gene
# 1 --> there is a gene on positive strand
# -2 --> there is a gene on negative strand

# Legacy Comments
# Each row is a pandas series in fact
# labels = ['#bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
# using cdsStart and cdsEnd give error since there are intervals of 0 length in that case
# such as cdsStart 33623833 cdsEnd 33623833
#Please notice that same genomic loci on a certain strand can be transcribed and untranscribed at the same time.
#However same genomic loci on a certain strand can be leading or lagging but not both
def fillTranscriptionArray(transcription_row,chrBased_genes_on_positive_strand,chrBased_genes_on_negative_strand):
    #gene on positive strand
    if (transcription_row['strand']==1):
        chrBased_genes_on_positive_strand[transcription_row['txStart']:(transcription_row['txEnd']+1)] = 1
    #gene on negative strand
    elif (transcription_row['strand']==-1):
        chrBased_genes_on_negative_strand[transcription_row['txStart']:(transcription_row['txEnd']+1)] = -2
########################################################################

########################################################################
def searchMutationUsingTranscriptionStrandColumn_simulations_integrated(
        mutation_row,
        simNum2Type2TranscriptionStrand2CountDict,
        simNum2Sample2Type2TranscriptionStrand2CountDict,
        simNum2Type2Sample2TranscriptionStrand2CountDict,
        simNum2Signature2MutationType2TranscriptionStrand2CountDict,
        signature2PropertiesListDict,
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
                                signature2PropertiesListDict)

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
                                signature2PropertiesListDict)
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
                               signature2PropertiesListDict)
        updateDictionaries_simulations_integrated(mutation_row,
                                mutationType,
                                mutationSample,
                                sample_based,
                                simNum2Type2TranscriptionStrand2CountDict,
                                simNum2Sample2Type2TranscriptionStrand2CountDict,
                                simNum2Type2Sample2TranscriptionStrand2CountDict,
                                simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                TRANSCRIBED_STRAND,
                                signature2PropertiesListDict)
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
                                signature2PropertiesListDict)

########################################################################


########################################################################
#This code uses transcription strand array
# TODO Consider NONTRANSCRIBED_STRAND
def searchMutationOnTranscriptionStrandArray(
        mutation_row,
        chrBased_gene_array,
        simNum2Type2TranscriptionStrand2CountDict,
        simNum2Sample2Type2TranscriptionStrand2CountDict,
        simNum2Type2Sample2TranscriptionStrand2CountDict,
        simNum2Signature2MutationType2TranscriptionStrand2CountDict,
        signature2PropertiesListDict,
        type,
        sample_based):

    mutationStart = mutation_row[START]
    mutationType = None
    mutationPyramidineStrand = mutation_row[PYRAMIDINESTRAND]
    mutationSample = mutation_row[SAMPLE]

    if (type==SUBS):
        mutationEnd = mutationStart+1
        #e.g.: C>A
        mutationType = mutation_row[MUTATION]
    elif (type==INDELS):
        mutationEnd = mutationStart+mutation_row[LENGTH]
        ref = mutation_row[REF]
        alt = mutation_row[ALT]

        if (len(ref)>len(alt)):
            mutation = ref[len(alt):]
        elif (len(alt)>len(ref)):
            mutation = alt[len(ref):]

        if (allPyrimidine(mutation)):
            mutationPyramidineStrand = 1
        elif (allPurine(mutation)):
            mutationPyramidineStrand = -1
        else:
            mutationPyramidineStrand = 0

    elif(type==DINUCS):
        mutationEnd = mutationStart+2

    # mutationPyramidineStrand= 1 --> pyrimidine mutation is on the + strand
    # mutationPyramidineStrand=-1 --> pyrimidine mutation is on the - strand
    # mutationPyramidineStrand= 0 --> mutation is not all pyrimidine or purine, therefore we can not decide the mutationPyramidineStrand (This is a case for DINUCS and INDELS)

    #Values on chrBased_transcription_array and their meanings
    # 0 --> no-transcription
    # 1 --> transcription on positive strand
    # -2 --> transcription on negative strand
    # -1 --> transcription on both positive and negative strands

    # uniqueIndexesArray = np.unique(chrBased_transcription_array[mutationStart:mutationEnd + 1])
    # mutationEnd is exclusive
    uniqueIndexesArray = np.unique(chrBased_gene_array[mutationStart:mutationEnd])

    if ((-1 in uniqueIndexesArray) or ((1 in uniqueIndexesArray) and (-2 in uniqueIndexesArray))):
        if (mutationPyramidineStrand != 0):
            updateDictionaries_simulations_integrated(mutation_row,
                                    mutationType,
                                    mutationSample,
                                    sample_based,
                                  simNum2Type2TranscriptionStrand2CountDict,
                                  simNum2Sample2Type2TranscriptionStrand2CountDict,
                                  simNum2Type2Sample2TranscriptionStrand2CountDict,
                                  simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                  TRANSCRIBED_STRAND,
                                  signature2PropertiesListDict)

            updateDictionaries_simulations_integrated(mutation_row,
                                    mutationType,
                                    mutationSample,
                                    sample_based,
                                  simNum2Type2TranscriptionStrand2CountDict,
                                  simNum2Sample2Type2TranscriptionStrand2CountDict,
                                  simNum2Type2Sample2TranscriptionStrand2CountDict,
                                  simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                  UNTRANSCRIBED_STRAND,
                                  signature2PropertiesListDict)

    elif (1 in uniqueIndexesArray):
        if (mutationPyramidineStrand == 1):
            # Transcription is on positive strand, if mutation pyramidine strand is + then increment untranscribed
            updateDictionaries_simulations_integrated(mutation_row,
                                    mutationType,
                                    mutationSample,
                                    sample_based,
                                    simNum2Type2TranscriptionStrand2CountDict,
                                    simNum2Sample2Type2TranscriptionStrand2CountDict,
                                    simNum2Type2Sample2TranscriptionStrand2CountDict,
                                    simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                    UNTRANSCRIBED_STRAND,
                                    signature2PropertiesListDict)

        elif (mutationPyramidineStrand == -1):
            # Transcription is on positive strand, if mutation pyramidine strand is - then increment transcribed
            updateDictionaries_simulations_integrated(mutation_row,
                                    mutationType,
                                    mutationSample,
                                    sample_based,
                                    simNum2Type2TranscriptionStrand2CountDict,
                                    simNum2Sample2Type2TranscriptionStrand2CountDict,
                                    simNum2Type2Sample2TranscriptionStrand2CountDict,
                                    simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                    TRANSCRIBED_STRAND,
                                    signature2PropertiesListDict)


    elif (-2 in uniqueIndexesArray):
        if (mutationPyramidineStrand == 1):
            # Transcription is on negative strand, if mutation pyramidine strand is + then increment transcribed
            updateDictionaries_simulations_integrated(mutation_row,
                                    mutationType,
                                    mutationSample,
                                    sample_based,
                                    simNum2Type2TranscriptionStrand2CountDict,
                                    simNum2Sample2Type2TranscriptionStrand2CountDict,
                                    simNum2Type2Sample2TranscriptionStrand2CountDict,
                                    simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                    TRANSCRIBED_STRAND,
                                    signature2PropertiesListDict)

        if (mutationPyramidineStrand == -1):
            # Transcription is on negative strand, if mutation pyramidine strand is - then increment untranscribed
            updateDictionaries_simulations_integrated(mutation_row,
                                    mutationType,
                                    mutationSample,
                                    sample_based,
                                    simNum2Type2TranscriptionStrand2CountDict,
                                    simNum2Sample2Type2TranscriptionStrand2CountDict,
                                    simNum2Type2Sample2TranscriptionStrand2CountDict,
                                    simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                    UNTRANSCRIBED_STRAND,
                                    signature2PropertiesListDict)
########################################################################



########################################################################
#DEC 24, 2019
#Called by USING_APPLY_ASYNC
def searchMutationsForApplySync(chrBased_simBased_subs_df,
                                chrBased_simBased_indels_df,
                                chrBased_simBased_dinucs_df,
                                chrBased_gene_array,
                                numofSimulations,
                                sample_based,
                                subsSignature2PropertiesListDict,
                                indelsSignature2PropertiesListDict,
                                dinucsSignature2PropertiesListDict,
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

    if (chrBased_gene_array is not None):
        if ((chrBased_simBased_subs_df is not None) and (not chrBased_simBased_subs_df.empty)):
            chrBased_simBased_subs_df.apply(searchMutationOnTranscriptionStrandArray,
                                    chrBased_gene_array=chrBased_gene_array,
                                    simNum2Type2TranscriptionStrand2CountDict = simNum2Type2TranscriptionStrand2CountDict,
                                    simNum2Sample2Type2TranscriptionStrand2CountDict = simNum2Sample2Type2TranscriptionStrand2CountDict,
                                    simNum2Type2Sample2TranscriptionStrand2CountDict = simNum2Type2Sample2TranscriptionStrand2CountDict,
                                    simNum2Signature2MutationType2TranscriptionStrand2CountDict = simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                    signature2PropertiesListDict=subsSignature2PropertiesListDict,
                                    type=SUBS,
                                    sample_based=sample_based,
                                    axis=1)

        if ((chrBased_simBased_indels_df is not None) and (not chrBased_simBased_indels_df.empty)):
            chrBased_simBased_indels_df.apply(searchMutationOnTranscriptionStrandArray,
                                    chrBased_gene_array=chrBased_gene_array,
                                    simNum2Type2TranscriptionStrand2CountDict = simNum2Type2TranscriptionStrand2CountDict,
                                    simNum2Sample2Type2TranscriptionStrand2CountDict = simNum2Sample2Type2TranscriptionStrand2CountDict,
                                    simNum2Type2Sample2TranscriptionStrand2CountDict = simNum2Type2Sample2TranscriptionStrand2CountDict,
                                    simNum2Signature2MutationType2TranscriptionStrand2CountDict =simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                    signature2PropertiesListDict=indelsSignature2PropertiesListDict,
                                    type= INDELS,
                                    sample_based=sample_based,
                                    axis=1)

        if ((chrBased_simBased_dinucs_df is not None) and (not chrBased_simBased_dinucs_df.empty)):
            chrBased_simBased_dinucs_df.apply(searchMutationOnTranscriptionStrandArray,
                                    chrBased_gene_array=chrBased_gene_array,
                                    simNum2Type2TranscriptionStrand2CountDict = simNum2Type2TranscriptionStrand2CountDict,
                                    simNum2Sample2Type2TranscriptionStrand2CountDict = simNum2Sample2Type2TranscriptionStrand2CountDict,
                                    simNum2Type2Sample2TranscriptionStrand2CountDict= simNum2Type2Sample2TranscriptionStrand2CountDict,
                                    simNum2Signature2MutationType2TranscriptionStrand2CountDict = simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                    signature2PropertiesListDict=dinucsSignature2PropertiesListDict,
                                    type=DINUCS,
                                    sample_based=sample_based,
                                    axis=1)

    elif (chrBased_gene_array is None):
        if ((chrBased_simBased_subs_df is not None) and (not chrBased_simBased_subs_df.empty)):
            if verbose: print('Worker pid %s SBS searchMutationUsingTranscriptionStrandColumn_simulations_integrated starts %s MB' % (str(os.getpid()), memory_usage()))
            chrBased_simBased_subs_df.apply(searchMutationUsingTranscriptionStrandColumn_simulations_integrated,
                                         simNum2Type2TranscriptionStrand2CountDict=simNum2Type2TranscriptionStrand2CountDict,
                                         simNum2Sample2Type2TranscriptionStrand2CountDict=simNum2Sample2Type2TranscriptionStrand2CountDict,
                                         simNum2Type2Sample2TranscriptionStrand2CountDict=simNum2Type2Sample2TranscriptionStrand2CountDict,
                                         simNum2Signature2MutationType2TranscriptionStrand2CountDict=simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                         signature2PropertiesListDict=subsSignature2PropertiesListDict,
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
                                           signature2PropertiesListDict=indelsSignature2PropertiesListDict,
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
                                           signature2PropertiesListDict=dinucsSignature2PropertiesListDict,
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
#Called by COMPUTATION_CHROMOSOMES_SEQUENTIAL_ALL_SIMULATIONS_PARALLEL
def searchMutations(inputList):
    chrBased_subs_split_df = inputList[0]
    chrBased_indels_split_df = inputList[1]
    chrBased_dinucs_split_df = inputList[2]
    chrBased_gene_array = inputList[3]
    numofSimulations = inputList[4]
    sample_based = inputList[5]
    subsSignature2PropertiesListDict=inputList[6]
    indelsSignature2PropertiesListDict = inputList[7]
    dinucsSignature2PropertiesListDict = inputList[8]

    simNum2Type2TranscriptionStrand2CountDict = {}
    simNum2Sample2Type2TranscriptionStrand2CountDict = {}
    simNum2Type2Sample2TranscriptionStrand2CountDict = {}
    simNum2Signature2MutationType2TranscriptionStrand2CountDict = {}

    for simNum in range(0,numofSimulations+1):
        simNum2Type2TranscriptionStrand2CountDict[simNum]={}
        simNum2Sample2Type2TranscriptionStrand2CountDict[simNum]={}
        simNum2Type2Sample2TranscriptionStrand2CountDict[simNum]={}
        simNum2Signature2MutationType2TranscriptionStrand2CountDict[simNum]={}

    if (chrBased_gene_array is not None):
        if ((chrBased_subs_split_df is not None) and (not chrBased_subs_split_df.empty)):
            chrBased_subs_split_df.apply(searchMutationOnTranscriptionStrandArray,
                                    chrBased_gene_array=chrBased_gene_array,
                                    simNum2Type2TranscriptionStrand2CountDict = simNum2Type2TranscriptionStrand2CountDict,
                                    simNum2Sample2Type2TranscriptionStrand2CountDict = simNum2Sample2Type2TranscriptionStrand2CountDict,
                                    simNum2Type2Sample2TranscriptionStrand2CountDict = simNum2Type2Sample2TranscriptionStrand2CountDict,
                                    simNum2Signature2MutationType2TranscriptionStrand2CountDict = simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                    signature2PropertiesListDict=subsSignature2PropertiesListDict,
                                    type=SUBS,
                                    sample_based=sample_based,
                                    axis=1)

        if ((chrBased_indels_split_df is not None) and (not chrBased_indels_split_df.empty)):
            chrBased_indels_split_df.apply(searchMutationOnTranscriptionStrandArray,
                                    chrBased_gene_array=chrBased_gene_array,
                                    simNum2Type2TranscriptionStrand2CountDict = simNum2Type2TranscriptionStrand2CountDict,
                                    simNum2Sample2Type2TranscriptionStrand2CountDict = simNum2Sample2Type2TranscriptionStrand2CountDict,
                                    simNum2Type2Sample2TranscriptionStrand2CountDict = simNum2Type2Sample2TranscriptionStrand2CountDict,
                                    simNum2Signature2MutationType2TranscriptionStrand2CountDict =simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                    signature2PropertiesListDict=indelsSignature2PropertiesListDict,
                                    type= INDELS,
                                    sample_based=sample_based,
                                    axis=1)


        if ((chrBased_dinucs_split_df is not None) and (not chrBased_dinucs_split_df.empty)):
            chrBased_dinucs_split_df.apply(searchMutationOnTranscriptionStrandArray,
                                    chrBased_gene_array=chrBased_gene_array,
                                    simNum2Type2TranscriptionStrand2CountDict = simNum2Type2TranscriptionStrand2CountDict,
                                    simNum2Sample2Type2TranscriptionStrand2CountDict = simNum2Sample2Type2TranscriptionStrand2CountDict,
                                    simNum2Type2Sample2TranscriptionStrand2CountDict= simNum2Type2Sample2TranscriptionStrand2CountDict,
                                    simNum2Signature2MutationType2TranscriptionStrand2CountDict = simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                    signature2PropertiesListDict=dinucsSignature2PropertiesListDict,
                                    type=DINUCS,
                                    sample_based=sample_based,
                                    axis=1)

    elif (chrBased_gene_array is None):
        if ((chrBased_subs_split_df is not None) and (not chrBased_subs_split_df.empty)):
            chrBased_subs_split_df.apply(searchMutationUsingTranscriptionStrandColumn_simulations_integrated,
                                         simNum2Type2TranscriptionStrand2CountDict=simNum2Type2TranscriptionStrand2CountDict,
                                         simNum2Sample2Type2TranscriptionStrand2CountDict=simNum2Sample2Type2TranscriptionStrand2CountDict,
                                         simNum2Type2Sample2TranscriptionStrand2CountDict=simNum2Type2Sample2TranscriptionStrand2CountDict,
                                         simNum2Signature2MutationType2TranscriptionStrand2CountDict=simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                         signature2PropertiesListDict=subsSignature2PropertiesListDict,
                                         type=SUBS,
                                         sample_based=sample_based,
                                         axis=1)

        if ((chrBased_indels_split_df is not None) and (not chrBased_indels_split_df.empty)):
            chrBased_indels_split_df.apply(searchMutationUsingTranscriptionStrandColumn_simulations_integrated,
                                           simNum2Type2TranscriptionStrand2CountDict=simNum2Type2TranscriptionStrand2CountDict,
                                           simNum2Sample2Type2TranscriptionStrand2CountDict=simNum2Sample2Type2TranscriptionStrand2CountDict,
                                           simNum2Type2Sample2TranscriptionStrand2CountDict=simNum2Type2Sample2TranscriptionStrand2CountDict,
                                           simNum2Signature2MutationType2TranscriptionStrand2CountDict=simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                           signature2PropertiesListDict=indelsSignature2PropertiesListDict,
                                           type=INDELS,
                                           sample_based=sample_based,
                                           axis=1)

        if ((chrBased_dinucs_split_df is not None) and (not chrBased_dinucs_split_df.empty)):
            chrBased_dinucs_split_df.apply(searchMutationUsingTranscriptionStrandColumn_simulations_integrated,
                                           simNum2Type2TranscriptionStrand2CountDict=simNum2Type2TranscriptionStrand2CountDict,
                                           simNum2Sample2Type2TranscriptionStrand2CountDict=simNum2Sample2Type2TranscriptionStrand2CountDict,
                                           simNum2Type2Sample2TranscriptionStrand2CountDict=simNum2Type2Sample2TranscriptionStrand2CountDict,
                                           simNum2Signature2MutationType2TranscriptionStrand2CountDict=simNum2Signature2MutationType2TranscriptionStrand2CountDict,
                                           signature2PropertiesListDict=dinucsSignature2PropertiesListDict,
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
def transcriptionStrandBiasAnalysis(computationType,sample_based,useTranscriptionStrandColumn,genome,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,verbose):

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
    if (not useTranscriptionStrandColumn):
        transcripts_df = readTrancriptsENSEMBL(genome)

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
            if (useTranscriptionStrandColumn):
                chrBased_transcription_array = None
            else:
                ################################################################################
                # Get chrBased ensembl transcripts
                chrBased_transcripts_df = transcripts_df[transcripts_df['chrom'] == chrLong]

                chrBased_transcription_array = np.zeros(chromSize, dtype=np.int8)
                chrBased_transcripts_df.apply(fillTranscriptionArray,chrBased_transcription_array=chrBased_transcription_array,axis=1)
                ################################################################################

            ####################################################################
            for simNum in range(0,numofSimulations+1):
                inputList = []
                chrBased_subs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, SUBS, simNum)
                chrBased_indels_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, INDELS, simNum)
                chrBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, DINUCS, simNum)
                inputList.append(chrBased_subs_df)  # each time different split
                inputList.append(chrBased_indels_df)
                inputList.append(chrBased_dinucs_df)
                inputList.append(chrBased_transcription_array) # same for all
                inputList.append(numofSimulations)
                inputList.append(sample_based)
                inputList.append(subsSignature2PropertiesListDict)
                inputList.append(indelsSignature2PropertiesListDict)
                inputList.append(dinucsSignature2PropertiesListDict)
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

            #You need to initialize to None so that you don't use former for loop values accidentally
            if (useTranscriptionStrandColumn):
                chrBased_transcription_array = None
            else:
                ################################################################################
                # Get chrBased ensembl transcripts
                chrBased_transcripts_df = transcripts_df[transcripts_df['chrom'] == chrLong]

                chrBased_transcription_array = np.zeros(chromSize, dtype=np.int8)
                chrBased_transcripts_df.apply(fillTranscriptionArray,chrBased_transcription_array=chrBased_transcription_array,axis=1)
                ################################################################################

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
                pool.apply_async(searchMutationsForApplySync, (chrBased_simBased_subs_df,chrBased_simBased_indels_df,chrBased_simBased_dinucs_df,chrBased_transcription_array,numofSimulations,sample_based,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,verbose),callback=accumulate_apply_async_result)
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