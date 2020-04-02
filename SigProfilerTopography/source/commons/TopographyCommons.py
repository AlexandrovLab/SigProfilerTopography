# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time and strand bias.
# Copyright (C) 2018-2020 Burcak Otlu

import os
import pandas as pd
import numpy as np
import math
import json
import pickle
import re
import shutil
import psutil

import scipy
from scipy.stats import sem, t
from scipy import mean

#To handle warnings as errors
# import warnings
# warnings.filterwarnings("error")

LINUX='linux'
UNIX='unix'
WINDOWS='windows'

current_abs_path = os.path.dirname(os.path.realpath(__file__))

LEADING= 'Leading'
LAGGING = 'Lagging'

UNTRANSCRIBED_STRAND = 'UnTranscribed'
TRANSCRIBED_STRAND = 'Transcribed'
NONTRANSCRIBED_STRAND = 'NonTranscribed'

PLUS = '+'
MINUS = '-'

MAXIMUM_CHROMOSOME_LENGTH = 250000000
# 1024*1024*1024 = 1073741824
# 1024*1024 = 1048576
GIGABYTE_IN_BYTES = 1073741824
MEGABYTE_IN_BYTES = 1048576

BIGWIG='BIGWIG'
BIGBED='BIGBED'
WIG='WIG'
BED='BED'
NARROWPEAK='narrowpeak'
LIBRARY_FILE_TYPE_OTHER='LIBRARY_FILE_TYPE_OTHER'

BED_6PLUS4='BED6+4'
BED_9PLUS2='BED9+2'

SAMPLE_MMR_DEFICIENT_THRESHOLD = 10000
SUBS_STRAND_BIAS_NUMBER_OF_MUTATIONS_THRESHOLD = 1000

DISPLAY = 'display'
NODISPLAY = 'nodisplay'

INPUT = 'input'
OUTPUT = 'output'
SIMULATION = 'simulation'
SIMULATIONS = 'simulations'
SIMULATIONS_FOR_TOPOGRAPHY = 'simulations_for_topography'
DATA = 'data'
CHRBASED = 'chrbased'
FIGURE = 'figure'

LIB = 'lib'
TRANSCRIPTS = 'transcripts'
NUCLEOSOME = 'nucleosome'
REPLICATION = 'replication'
UCSCGENOME = 'ucscgenome'

HISTONE_MODIFICATION= "histone modification"
TRANSCRIPTION_FACTOR= "transcription factor"

NCBI = 'ncbi'
ENSEMBL = 'ensembl'

GM12878 = 'GM12878'
K562 = 'K562'
MCF7='MCF7'
HEPG2='HEPG2'
HELAS3='HELAS3'
SKNSH='SKNSH'
IMR90='IMR90'
NHEK='NHEK'
BJ='BJ'
HUVEC='HUVEC'
BG02ES='BG02ES'
GM06990='GM06990'
GM12801='GM12801'
GM12812='GM12812'
GM12813='GM12813'

#NUCLEOSOME OCCUPANCY FILES
GM12878_NUCLEOSOME_OCCUPANCY_FILE = 'wgEncodeSydhNsomeGm12878Sig.bigWig'
K562_NUCLEOSOME_OCCUPANCY_FILE = 'wgEncodeSydhNsomeK562Sig.bigWig'

#REPLICATION  TIME FILES
MCF7_REPLICATION_TIME_SIGNAL_FILE = 'GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.wig'
MCF7_REPLICATION_TIME_VALLEY_FILE = 'GSM923442_hg19_wgEncodeUwRepliSeqMcf7ValleysRep1.bed'
MCF7_REPLICATION_TIME_PEAK_FILE = 'GSM923442_hg19_wgEncodeUwRepliSeqMcf7PkRep1.bed'

HEPG2_REPLICATION_TIME_SIGNAL_FILE = 'wgEncodeUwRepliSeqHepg2WaveSignalRep1.wig'
HEPG2_REPLICATION_TIME_VALLEY_FILE = 'wgEncodeUwRepliSeqHepg2ValleysRep1.bed'
HEPG2_REPLICATION_TIME_PEAK_FILE = 'wgEncodeUwRepliSeqHepg2PkRep1.bed'

HELAS3_REPLICATION_TIME_SIGNAL_FILE = 'wgEncodeUwRepliSeqHelas3WaveSignalRep1.wig'
HELAS3_REPLICATION_TIME_VALLEY_FILE = 'wgEncodeUwRepliSeqHelas3ValleysRep1.bed'
HELAS3_REPLICATION_TIME_PEAK_FILE = 'wgEncodeUwRepliSeqHelas3PkRep1.bed'

SKNSH_REPLICATION_TIME_SIGNAL_FILE = 'wgEncodeUwRepliSeqSknshWaveSignalRep1.wig'
SKNSH_REPLICATION_TIME_VALLEY_FILE = 'wgEncodeUwRepliSeqSknshValleysRep1.bed'
SKNSH_REPLICATION_TIME_PEAK_FILE = 'wgEncodeUwRepliSeqSknshPkRep1.bed'

K562_REPLICATION_TIME_SIGNAL_FILE = 'wgEncodeUwRepliSeqK562WaveSignalRep1.wig'
K562_REPLICATION_TIME_VALLEY_FILE = 'wgEncodeUwRepliSeqK562ValleysRep1.bed'
K562_REPLICATION_TIME_PEAK_FILE = 'wgEncodeUwRepliSeqK562PkRep1.bed'

IMR90_REPLICATION_TIME_SIGNAL_FILE = 'wgEncodeUwRepliSeqImr90WaveSignalRep1.wig'
IMR90_REPLICATION_TIME_VALLEY_FILE = 'wgEncodeUwRepliSeqImr90ValleysRep1.bed'
IMR90_REPLICATION_TIME_PEAK_FILE = 'wgEncodeUwRepliSeqImr90PkRep1.bed'

NHEK_REPLICATION_TIME_SIGNAL_FILE = 'wgEncodeUwRepliSeqNhekWaveSignalRep1.wig'
NHEK_REPLICATION_TIME_VALLEY_FILE = 'wgEncodeUwRepliSeqNhekValleysRep1.bed'
NHEK_REPLICATION_TIME_PEAK_FILE = 'wgEncodeUwRepliSeqNhekPkRep1.bed'

BJ_REPLICATION_TIME_SIGNAL_FILE = 'wgEncodeUwRepliSeqBjWaveSignalRep1.wig'
BJ_REPLICATION_TIME_VALLEY_FILE = 'wgEncodeUwRepliSeqBjValleysRep1.bed'
BJ_REPLICATION_TIME_PEAK_FILE = 'wgEncodeUwRepliSeqBjPkRep1.bed'

HUVEC_REPLICATION_TIME_SIGNAL_FILE = 'wgEncodeUwRepliSeqHuvecWaveSignalRep1.wig'
HUVEC_REPLICATION_TIME_VALLEY_FILE = 'wgEncodeUwRepliSeqHuvecValleysRep1.bed'
HUVEC_REPLICATION_TIME_PEAK_FILE = 'wgEncodeUwRepliSeqHuvecPkRep1.bed'

BG02ES_REPLICATION_TIME_SIGNAL_FILE = 'wgEncodeUwRepliSeqBg02esWaveSignalRep1.wig'
BG02ES_REPLICATION_TIME_VALLEY_FILE = 'wgEncodeUwRepliSeqBg02esValleysRep1.bed'
BG02ES_REPLICATION_TIME_PEAK_FILE = 'wgEncodeUwRepliSeqBg02esPkRep1.bed'

GM12878_REPLICATION_TIME_SIGNAL_FILE = 'wgEncodeUwRepliSeqGm12878WaveSignalRep1.wig'
GM12878_REPLICATION_TIME_VALLEY_FILE = 'wgEncodeUwRepliSeqGm12878ValleysRep1.bed'
GM12878_REPLICATION_TIME_PEAK_FILE = 'wgEncodeUwRepliSeqGm12878PkRep1.bed'

GM06990_REPLICATION_TIME_SIGNAL_FILE = 'wgEncodeUwRepliSeqGm06990WaveSignalRep1.wig'
GM06990_REPLICATION_TIME_VALLEY_FILE = 'wgEncodeUwRepliSeqGm06990ValleysRep1.bed'
GM06990_REPLICATION_TIME_PEAK_FILE = 'wgEncodeUwRepliSeqGm06990PkRep1.bed'

GM12801_REPLICATION_TIME_SIGNAL_FILE = 'wgEncodeUwRepliSeqGm12801WaveSignalRep1.wig'
GM12801_REPLICATION_TIME_VALLEY_FILE = 'wgEncodeUwRepliSeqGm12801ValleysRep1.bed'
GM12801_REPLICATION_TIME_PEAK_FILE = 'wgEncodeUwRepliSeqGm12801PkRep1.bed'

GM12812_REPLICATION_TIME_SIGNAL_FILE = 'wgEncodeUwRepliSeqGm12812WaveSignalRep1.wig'
GM12812_REPLICATION_TIME_VALLEY_FILE = 'wgEncodeUwRepliSeqGm12812ValleysRep1.bed'
GM12812_REPLICATION_TIME_PEAK_FILE = 'wgEncodeUwRepliSeqGm12812PkRep1.bed'

GM12813_REPLICATION_TIME_SIGNAL_FILE = 'wgEncodeUwRepliSeqGm12813WaveSignalRep1.wig'
GM12813_REPLICATION_TIME_VALLEY_FILE = 'wgEncodeUwRepliSeqGm12813ValleysRep1.bed'
GM12813_REPLICATION_TIME_PEAK_FILE = 'wgEncodeUwRepliSeqGm12813PkRep1.bed'

available_nucleosome_biosamples=[GM12878,K562]
available_replication_time_biosamples=[GM12878,K562,MCF7,HEPG2,HELAS3,SKNSH,IMR90,NHEK,BJ,HUVEC,BG02ES,GM06990,GM12801,GM12812,GM12813]

GRCh37_hg19_NCBIREFSEQCURATED = 'GRCh37_hg19_NCBIRefSeqCurated'
GRCh37_ENSEMBL = 'GRCh37_transcripts.txt'

HG19_CHROM_SIZES = 'hg19.chrom.sizes.txt'
HG38_CHROM_SIZES = 'hg38.chrom.sizes.txt'

ENCODE_NUCLEOSOME_GM12878_BIGWIG_URL = 'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeSydhNsome/wgEncodeSydhNsomeGm12878Sig.bigWig'
ENCODE_NUCLEOSOME_K562_BIGWIG_URL = 'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeSydhNsome/wgEncodeSydhNsomeK562Sig.bigWig'
# BIGWIG_TO_WIG_EXECUTABLE_LINUX_X86_64_URL = 'http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig'
BIGWIG_TO_WIG_EXECUTABLE_LINUX_X86_64_URL = 'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v369/bigWigToWig'

HG19 = 'hg19'
HG38 = 'hg38'
MM9 = 'mm9'
MM10 = 'mm10'
GRCh37 = 'GRCh37'
GRCh38 = 'GRCh38'

FIXED_STEP='fixedStep'
VARIABLE_STEP='variableStep'

###################################################################################################
#Tables
Table_AllCutoff_SubsSignature_NumberofMutations_AverageProbability_Filename = "Table_AllCutoffs_SBS_Signature_NumberofMutations_AverageProbability.txt"
Table_AllCutoff_IndelsSignature_NumberofMutations_AverageProbability_Filename = "Table_AllCutoffs_ID_Signature_NumberofMutations_AverageProbability.txt"
Table_AllCutoff_DinucsSignature_NumberofMutations_AverageProbability_Filename = "Table_AllCutoffs_DBS_Signature_NumberofMutations_AverageProbability.txt"

#Tables
Table_SubsSignature_Cutoff_NumberofMutations_AverageProbability_Filename = "Table_SBS_Signature_Cutoff_NumberofMutations_AverageProbability.txt"
Table_IndelsSignature_Cutoff_NumberofMutations_AverageProbability_Filename = "Table_ID_Signature_Cutoff_NumberofMutations_AverageProbability.txt"
Table_DinucsSignature_Cutoff_NumberofMutations_AverageProbability_Filename = "Table_DBS_Signature_Cutoff_NumberofMutations_AverageProbability.txt"

#Table
Table_MutationType_NumberofMutations_Filename='Table_MutationType_NumberofMutations.txt'

#For Subs
Sample2NumberofSubsDictFilename = 'Sample2NumberofSubsDict.txt'
SubsSignature2NumberofMutationsDictFilename = 'SubsSignature2NumberofMutationsDict.txt'
Sample2SubsSignature2NumberofMutationsDictFilename = 'Sample2SubsSignature2NumberofMutationsDict.txt'

#For Indels
Sample2NumberofIndelsDictFilename = 'Sample2NumberofIndelsDict.txt'
IndelsSignature2NumberofMutationsDictFilename = 'IndelsSignature2NumberofMutationsDict.txt'
Sample2IndelsSignature2NumberofMutationsDictFilename = 'Sample2IndelsSignature2NumberofMutationsDict.txt'

#For Dinucs
Sample2NumberofDinucsDictFilename = 'Sample2NumberofDinucsDict.txt'
DinucsSignature2NumberofMutationsDictFilename = 'DinucsSignature2NumberofMutationsDict.txt'
Sample2DinucsSignature2NumberofMutationsDictFilename = 'Sample2DinucsSignature2NumberofMutationsDict.txt'
###################################################################################################

#For Replication
DecileIndex2NumfAttributableBasesDictFilename = 'DecileIndex2NumfAttributableBasesDict.txt'

ONE_DIRECTORY_UP = '..'

GRCh37ChromSizesDictFilename = 'hg19ChromSizesDict.txt'
GRCh38ChromSizesDictFilename = 'hg38ChromSizesDict.txt'
INDELBASED = 'indelbased'
SIGNATUREBASED = 'signaturebased'
SAMPLEBASED = 'samplebased'

SAMPLEBASED_SIGNATUREBASED = 'samplebased_signaturebased'

SAMPLEBASED_AGGREGATEDSUBSTITUTIONS = 'samplebased_aggregatedsubstitutions'
SAMPLEBASED_AGGREGATEDINDELS = 'samplebased_aggregatedindels'
SAMPLEBASED_AGGREGATEDDINUCS = 'samplebased_aggregateddinucs'

AGGREGATEDDINUCS = 'aggregateddinucs'
AGGREGATEDINDELS = 'aggregatedindels'
AGGREGATEDSUBSTITUTIONS = 'aggregatedsubstitutions'

SOURCE = 'source'
COMMONS = 'commons'
EPIGENOMICS = 'epigenomics'
NUCLEOSOMEOCCUPANCY = 'nucleosome_occupancy'
EPIGENOMICSOCCUPANCY = 'epigenomics_occupancy'
TRANSCRIPTIONFACTOROCCUPANCY='tf_occupancy'
REPLICATIONTIME = 'replication_time'
PROCESSIVITY = 'processivity'
STRANDBIAS = 'strand_bias'
SCATTERPLOTS= 'scatter_plots'
HEATMAPS='heatmaps'
PLOTTING = 'plotting'
TRANSCRIPTIONSTRANDBIAS = 'transcription_strand_bias'
REPLICATIONSTRANDBIAS = 'replication_strand_bias'

ALL = 'all'
SAMPLES = 'samples'

TRANSCRIPTION_LABEL = 'Transcription'
REPLICATION_LABEL = 'Replication'

NUMBER_OF_MUTATIONS = 'number_of_mutations'
MUTATION_DENSITY = 'mutation_density'
NORMALIZED_MUTATION_DENSITY = 'normalized_mutation_density'

MICROHOMOLOGY = 'Microhomology-mediated indels'
REPEAT = 'Repeat-mediated indels'

DEFICIENT = 'deficient'
PROFICIENT = 'proficient'

#Mutation Types
C2A = 'C>A'
C2G = 'C>G'
C2T = 'C>T'
T2A = 'T>A'
T2C = 'T>C'
T2G = 'T>G'

sixMutationTypes = [C2A,C2G,C2T,T2A,T2C,T2G]

SIGNATURE = 'Signature'
MUTATION = 'Mutation'
MUTATIONLONG = 'MutationLong'

COMPUTATION_CHROMOSOMES_SEQUENTIAL_ALL_SIMULATIONS_PARALLEL = 'COMPUTATION_CHROMOSOMES_SEQUENTIAL_ALL_SIMULATIONS_PARALLEL'
USING_IMAP_UNORDERED='USING_IMAP_UNORDERED'
USING_APPLY_ASYNC='USING_APPLY_ASYNC'

############################################
#Column Names
PROJECT = 'Project'
SAMPLE = 'Sample'
GENOME = 'Genome'

CHROM = 'Chrom'
START = 'Start'
END = 'End'
SIGNAL = 'Signal'
NAME='Name'
SCORE='Score'
STRAND='Strand'

SIZE='Size'

REF = 'Ref'
ALT = 'Alt'
PYRAMIDINESTRAND = 'PyramidineStrand'
TRANSCRIPTIONSTRAND = 'TranscriptionStrand'
MUTATION = 'Mutation'
MUTATIONS = 'Mutations'
CONTEXT = 'Context'

#For Replication Time Data
NUMOFBASES = 'numofBases'

COUNT = 'Count'
MMR = 'MMR'

TYPE = 'Type'
LENGTH = 'Length'
CATEGORY = 'Category'

INDEL = 'INDEL'
INDEL_GT_3BP = 'INDEL(>3bp)'
INDEL_LTE_3BP = 'INDEL(<=3bp)'

DOT= 'Dot'
DATASOURCE = 'DataSource'

SIMULATION_NUMBER= 'Simulation_Number'
############################################

Type2ReplicationStrand2CountDict_Filename = 'Type2ReplicationStrand2CountDict.txt'
Sample2Type2ReplicationStrand2CountDict_Filename = 'Sample2Type2ReplicationStrand2CountDict.txt'
Type2Sample2ReplicationStrand2CountDict_Filename = 'Type2Sample2ReplicationStrand2CountDict.txt'
Signature2MutationType2ReplicationStrand2CountDict_Filename  = 'Signature2MutationType2ReplicationStrand2CountDict.txt'

Type2TranscriptionStrand2CountDict_Filename = 'Type2TranscriptionStrand2CountDict.txt'
Sample2Type2TranscriptionStrand2CountDict_Filename = 'Sample2Type2TranscriptionStrand2CountDict.txt'
Type2Sample2TranscriptionStrand2CountDict_Filename = 'Type2Sample2TranscriptionStrand2CountDict.txt'
Signature2MutationType2TranscriptionStrand2CountDict_Filename  = 'Signature2MutationType2TranscriptionStrand2CountDict.txt'

DATA_Folder = 'data'

BONFERRONI_CORRECTION = 'BONFERRONI_CORRECTION'
FDR_BH_CORRECTION = 'FDR_BH_CORRECTION'

USING_ONE_SAMPLE_TTEST = 'USING_ONE_SAMPLE_TTEST'
USING_NULL_DISTRIBUTION = 'USING_NULL_DISTRIBUTION'
USING_GAUSSIAN_KDE = 'USING_GAUSSIAN_KDE'
USING_ZSCORE = 'USING_ZSCORE'

SBS96   = '96'
SBS192  = '192'
SBS384  = '384'
SBS1536 = '1536'
SBS3072 = '3072'

SNV='SNV'
ID = 'ID'
DBS= 'DBS'
SBS='SBS'
SBS_CONTEXTS = [SBS96,SBS192,SBS384,SBS1536,SBS3072]

# Used for dictionaries
# MutationType2NumberofMutationsDict keys
SUBS = 'SUBS'
INDELS = 'INDELS'
DINUCS = 'DINUCS'

DEFAULT_HISTONE_OCCUPANCY_FILE1 = 'ENCFF291WFP_H3K27me3_breast_epithelium.bed'
DEFAULT_HISTONE_OCCUPANCY_FILE2 = 'ENCFF906MJM_H3K36me3_breast_epithelium.bed'
DEFAULT_HISTONE_OCCUPANCY_FILE3 = 'ENCFF065FJK_H3K9me3_breast_epithelium.bed'
DEFAULT_HISTONE_OCCUPANCY_FILE4 = 'ENCFF154XFN_H3K27ac_breast_epithelium.bed'
DEFAULT_HISTONE_OCCUPANCY_FILE5 = 'ENCFF336DDM_H3K4me1_breast_epithelium.bed'
DEFAULT_HISTONE_OCCUPANCY_FILE6 = 'ENCFF065TIH_H3K4me3_breast_epithelium.bed'

BIOSAMPLE_UNDECLARED='Biosample_Undeclared'


########################################################
def get_splits(outputDir, jobname, simNum,chrLong):
    #################################################################################################################
    # If library file does not exists there is no library file to use and fill the signal and count arrays
    # Nucleosomes have chrM
    # SinglePointMutations and Indels have chrMT

    chrLong_for_mutations_data = chrLong
    if (chrLong == 'chrM'):
        chrLong_for_mutations_data = 'chrMT'

    chrBased_subs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong_for_mutations_data, SUBS, simNum)
    chrBased_indels_df = readChrBasedMutationsDF(outputDir, jobname, chrLong_for_mutations_data, INDELS, simNum)
    chrBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong_for_mutations_data, DINUCS, simNum)

    #default behaviour
    chrBased_subs_df_split_list = [chrBased_subs_df]
    chrBased_indels_df_split_list = [chrBased_indels_df]
    chrBased_dinucs_df_split_list = [chrBased_dinucs_df]

    number_of_mutations_list = []
    order_in_number_of_mutations_list = []

    #There is an order in the way we put the elements in the list
    if (chrBased_subs_df is not None):
        number_of_subs_mutations = chrBased_subs_df.shape[0]
        number_of_mutations_list.append(number_of_subs_mutations)
        order_in_number_of_mutations_list.append(AGGREGATEDSUBSTITUTIONS)
    if (chrBased_indels_df is not None):
        number_of_indels_mutations = chrBased_indels_df.shape[0]
        number_of_mutations_list.append(number_of_indels_mutations)
        order_in_number_of_mutations_list.append(AGGREGATEDINDELS)
    if (chrBased_dinucs_df is not None):
        number_of_dinucs_mutations = chrBased_dinucs_df.shape[0]
        number_of_mutations_list.append(number_of_dinucs_mutations)
        order_in_number_of_mutations_list.append(AGGREGATEDDINUCS)

    ##############################################################################
    if (len(number_of_mutations_list) == 2):
        #There are 3 possible content for the list

        #Subs, Indels
        if order_in_number_of_mutations_list[0]==AGGREGATEDSUBSTITUTIONS and order_in_number_of_mutations_list[1]==AGGREGATEDINDELS:
            if number_of_indels_mutations <= number_of_subs_mutations:
                number_of_splits = number_of_subs_mutations / number_of_indels_mutations
                chrBased_subs_df_split_list = np.array_split(chrBased_subs_df, number_of_splits)
            elif number_of_subs_mutations <= number_of_indels_mutations:
                number_of_splits = number_of_indels_mutations / number_of_subs_mutations
                chrBased_indels_df_split_list = np.array_split(chrBased_indels_df, number_of_splits)

        #Subs, Dinucs
        elif order_in_number_of_mutations_list[0]==AGGREGATEDSUBSTITUTIONS and order_in_number_of_mutations_list[1]==AGGREGATEDDINUCS:
            if number_of_dinucs_mutations <= number_of_subs_mutations:
                number_of_splits = number_of_subs_mutations / number_of_dinucs_mutations
                chrBased_subs_df_split_list = np.array_split(chrBased_subs_df, number_of_splits)
            elif number_of_subs_mutations <= number_of_dinucs_mutations:
                number_of_splits = number_of_dinucs_mutations / number_of_subs_mutations
                chrBased_dinucs_df_split_list = np.array_split(chrBased_dinucs_df, number_of_splits)

        #Indels, Dinucs
        elif order_in_number_of_mutations_list[0]==AGGREGATEDINDELS and order_in_number_of_mutations_list[1]==AGGREGATEDDINUCS:
            if number_of_dinucs_mutations <= number_of_indels_mutations:
                number_of_splits = number_of_indels_mutations / number_of_dinucs_mutations
                chrBased_indels_df_split_list = np.array_split(chrBased_indels_df, number_of_splits)
            elif number_of_indels_mutations <= number_of_dinucs_mutations:
                number_of_splits = number_of_dinucs_mutations / number_of_indels_mutations
                chrBased_dinucs_df_split_list = np.array_split(chrBased_dinucs_df, number_of_splits)
    ##############################################################################


    ##############################################################################
    elif len(number_of_mutations_list) == 3:

        #########################################################
        #  Dinucs <= Indels <= Subs  Most probable case
        if ((number_of_dinucs_mutations <= number_of_indels_mutations) and (number_of_indels_mutations <= number_of_subs_mutations)):
            number_of_splits = number_of_subs_mutations / number_of_indels_mutations
            chrBased_subs_df_split_list = np.array_split(chrBased_subs_df, number_of_splits)

        # Indels <= Dinucs <= Subs
        elif ((number_of_indels_mutations <= number_of_dinucs_mutations) and (number_of_dinucs_mutations <= number_of_subs_mutations)):
            # split number_of_mutations_list[0]
            number_of_splits = number_of_subs_mutations / number_of_dinucs_mutations
            chrBased_subs_df_split_list = np.array_split(chrBased_subs_df, number_of_splits)
        #########################################################

        #########################################################
        # Dinucs <= Subs <= Indels
        # number_of_mutations_list[0] is between
        # order  2-0-1
        elif ((number_of_dinucs_mutations <= number_of_subs_mutations) and (number_of_subs_mutations <= number_of_indels_mutations)):
            # split number_of_mutations_list[1]
            number_of_splits = number_of_indels_mutations / number_of_subs_mutations
            chrBased_indels_df_split_list = np.array_split(chrBased_indels_df, number_of_splits)

        # Subs <= Dinucs <= Indels
        # order  0-2-1
        # number_of_mutations_list[2] is between
        elif ((number_of_subs_mutations <= number_of_dinucs_mutations) and (number_of_dinucs_mutations <= number_of_indels_mutations)):
            # split number_of_mutations_list[1]
            number_of_splits = number_of_indels_mutations / number_of_dinucs_mutations
            chrBased_indels_df_split_list = np.array_split(chrBased_indels_df, number_of_splits)
        #########################################################

        #########################################################
        # Subs <= Indels <= Dinucs  Less Likely
        # number_of_mutations_list[1] is between
        # order  0-1-2
        elif ((number_of_subs_mutations <= number_of_indels_mutations) and (number_of_indels_mutations <= number_of_dinucs_mutations)):
            # split number_of_mutations_list[2]
            number_of_splits = number_of_dinucs_mutations / number_of_indels_mutations
            chrBased_dinucs_df_split_list = np.array_split(chrBased_dinucs_df, number_of_splits)

        # Indels <= Subs <= Dinucs
        # number_of_mutations_list[0] is between
        # order  1-0-2
        elif ((number_of_indels_mutations <= number_of_subs_mutations) and (number_of_subs_mutations <= number_of_dinucs_mutations)):
            # split number_of_mutations_list[2]
            number_of_splits = number_of_dinucs_mutations / number_of_subs_mutations
            chrBased_dinucs_df_split_list = np.array_split(chrBased_dinucs_df, number_of_splits)
        #########################################################
    ##############################################################################

    return chrBased_subs_df_split_list, chrBased_indels_df_split_list, chrBased_dinucs_df_split_list
########################################################

########################################################
# e.g: chrLong='chr3' chrShort='3'
def getChrShort(chrLong):
    chrShort = chrLong[3:]
    return chrShort
########################################################

########################################################
def natural_key(string_):
    """See http://www.codinghorror.com/blog/archives/001018.html"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]
########################################################


#######################################################
def getNucleosomeFile(nucleosome_biosample):
    nucleosome_file=None

    if (nucleosome_biosample is not None):
        if (nucleosome_biosample==GM12878):
            nucleosome_file = GM12878_NUCLEOSOME_OCCUPANCY_FILE
        elif (nucleosome_biosample==K562):
            nucleosome_file = K562_NUCLEOSOME_OCCUPANCY_FILE

    return nucleosome_file
#######################################################

#######################################################
def getReplicationTimeFiles(replication_time_biosample):
    replication_time_signal_file=None
    replication_time_valley_file=None
    replication_time_peak_file=None

    if replication_time_biosample is not None:
        if (replication_time_biosample==GM12878):
            replication_time_signal_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB, REPLICATION,GM12878_REPLICATION_TIME_SIGNAL_FILE)
            replication_time_valley_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB, REPLICATION,GM12878_REPLICATION_TIME_VALLEY_FILE)
            replication_time_peak_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB, REPLICATION,GM12878_REPLICATION_TIME_PEAK_FILE)
        elif (replication_time_biosample == K562):
            replication_time_signal_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB, REPLICATION,K562_REPLICATION_TIME_SIGNAL_FILE)
            replication_time_valley_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB, REPLICATION,K562_REPLICATION_TIME_VALLEY_FILE)
            replication_time_peak_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB, REPLICATION,K562_REPLICATION_TIME_PEAK_FILE)
        elif (replication_time_biosample == MCF7):
            replication_time_signal_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB, REPLICATION,MCF7_REPLICATION_TIME_SIGNAL_FILE)
            replication_time_valley_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,MCF7_REPLICATION_TIME_VALLEY_FILE)
            replication_time_peak_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,MCF7_REPLICATION_TIME_PEAK_FILE)
        elif (replication_time_biosample == HEPG2):
            replication_time_signal_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,HEPG2_REPLICATION_TIME_SIGNAL_FILE)
            replication_time_valley_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,HEPG2_REPLICATION_TIME_VALLEY_FILE)
            replication_time_peak_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,HEPG2_REPLICATION_TIME_PEAK_FILE)
        elif (replication_time_biosample == HELAS3):
            replication_time_signal_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,HELAS3_REPLICATION_TIME_SIGNAL_FILE)
            replication_time_valley_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,HELAS3_REPLICATION_TIME_VALLEY_FILE)
            replication_time_peak_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,HELAS3_REPLICATION_TIME_PEAK_FILE)
        elif (replication_time_biosample == SKNSH):
            replication_time_signal_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,SKNSH_REPLICATION_TIME_SIGNAL_FILE)
            replication_time_valley_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,SKNSH_REPLICATION_TIME_VALLEY_FILE)
            replication_time_peak_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,SKNSH_REPLICATION_TIME_PEAK_FILE)
        elif (replication_time_biosample == IMR90):
            replication_time_signal_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,IMR90_REPLICATION_TIME_SIGNAL_FILE)
            replication_time_valley_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,IMR90_REPLICATION_TIME_VALLEY_FILE)
            replication_time_peak_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,IMR90_REPLICATION_TIME_PEAK_FILE)
        elif (replication_time_biosample == NHEK):
            replication_time_signal_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,NHEK_REPLICATION_TIME_SIGNAL_FILE)
            replication_time_valley_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,NHEK_REPLICATION_TIME_VALLEY_FILE)
            replication_time_peak_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,NHEK_REPLICATION_TIME_PEAK_FILE)
        elif (replication_time_biosample == BJ):
            replication_time_signal_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,BJ_REPLICATION_TIME_SIGNAL_FILE)
            replication_time_valley_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,BJ_REPLICATION_TIME_VALLEY_FILE)
            replication_time_peak_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,BJ_REPLICATION_TIME_PEAK_FILE)
        elif (replication_time_biosample == HUVEC):
            replication_time_signal_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,HUVEC_REPLICATION_TIME_SIGNAL_FILE)
            replication_time_valley_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,HUVEC_REPLICATION_TIME_VALLEY_FILE)
            replication_time_peak_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION, HUVEC_REPLICATION_TIME_PEAK_FILE)
        elif (replication_time_biosample == BG02ES):
            replication_time_signal_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,BG02ES_REPLICATION_TIME_SIGNAL_FILE)
            replication_time_valley_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,BG02ES_REPLICATION_TIME_VALLEY_FILE)
            replication_time_peak_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION, BG02ES_REPLICATION_TIME_PEAK_FILE)
        elif (replication_time_biosample == GM06990):
            replication_time_signal_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,GM06990_REPLICATION_TIME_SIGNAL_FILE)
            replication_time_valley_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,GM06990_REPLICATION_TIME_VALLEY_FILE)
            replication_time_peak_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION, GM06990_REPLICATION_TIME_PEAK_FILE)
        elif (replication_time_biosample == GM12801):
            replication_time_signal_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,GM12801_REPLICATION_TIME_SIGNAL_FILE)
            replication_time_valley_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,GM12801_REPLICATION_TIME_VALLEY_FILE)
            replication_time_peak_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,GM12801_REPLICATION_TIME_PEAK_FILE)
        elif (replication_time_biosample == GM12812):
            replication_time_signal_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,GM12812_REPLICATION_TIME_SIGNAL_FILE)
            replication_time_valley_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,GM12812_REPLICATION_TIME_VALLEY_FILE)
            replication_time_peak_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION, GM12812_REPLICATION_TIME_PEAK_FILE)
        elif (replication_time_biosample == GM12813):
            replication_time_signal_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,GM12813_REPLICATION_TIME_SIGNAL_FILE)
            replication_time_valley_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION,GM12813_REPLICATION_TIME_VALLEY_FILE)
            replication_time_peak_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, LIB, REPLICATION, GM12813_REPLICATION_TIME_PEAK_FILE)

    return replication_time_signal_file, replication_time_valley_file, replication_time_peak_file
#######################################################


#############################################################################
#Modified DEC 9, 2019
def takeAverage(listofSimulationsAggregatedMutations):
    simulationsAggregatedMutationsLows = None
    simulationsAggregatedMutationsHighs = None
    simulationsAggregatedMutationsMeans = None

    #Number of simulations >= 2
    if ((listofSimulationsAggregatedMutations is not None) and len(listofSimulationsAggregatedMutations)>=2):
        stackedSimulationAggregatedMutations = np.vstack(listofSimulationsAggregatedMutations)
        (rows, cols) = stackedSimulationAggregatedMutations.shape

        simulationsAggregatedMutationsLows = []
        simulationsAggregatedMutationsHighs = []
        simulationsAggregatedMutationsMeans = []

        for col in range(cols):
            colwise_array = stackedSimulationAggregatedMutations[:, col]

            # ########################
            # #1st way : When there are 2 (small number of) simulations it provides very distant start and end
            # confidence = 0.95
            # n = len(colwise_array)
            # m = mean(colwise_array)
            # std_err = sem(colwise_array)
            # h = std_err * t.ppf((1 + confidence) / 2, n - 1)
            # start = m - h
            # end = m + h
            # ########################

            ########################
            #2nd way
            mu = mean(colwise_array)
            sigma = sem(colwise_array)
            start, end = scipy.stats.norm.interval(0.95, loc=mu, scale=sigma)
            ########################

            simulationsAggregatedMutationsLows.append(start)
            simulationsAggregatedMutationsMeans.append(np.mean(colwise_array))
            simulationsAggregatedMutationsHighs.append(end)

    # Number of simulations == 1
    elif ((listofSimulationsAggregatedMutations is not None) and len(listofSimulationsAggregatedMutations)==1):
        stackedSimulationAggregatedMutations = np.vstack(listofSimulationsAggregatedMutations)
        (rows, cols) = stackedSimulationAggregatedMutations.shape

        simulationsAggregatedMutationsMeans = []

        for col in range(cols):
            colwise_array = stackedSimulationAggregatedMutations[:, col]
            simulationsAggregatedMutationsMeans.append(np.mean(colwise_array))


    return  simulationsAggregatedMutationsLows,simulationsAggregatedMutationsMeans,simulationsAggregatedMutationsHighs
#############################################################################

###################################################################
def calc_chunksize(n_workers, len_iterable, factor=4):
    """Calculate chunksize argument for Pool-methods.

    Resembles source-code within `multiprocessing.pool.Pool._map_async`.
    """
    chunksize, extra = divmod(len_iterable, n_workers * factor)
    if extra:
        chunksize += 1
    return chunksize
###################################################################

###########################################################
def memory_usage():
    pid = os.getpid()
    process = psutil.Process(pid)
    # memoryUseInGB = process.memory_info()[0]/2.**30  #memory use in GB
    memoryUseInMB = process.memory_info()[0]/2.**20  # memory use in MB
    return memoryUseInMB
    # print('************** Current Memory Use: '+ str(round(memoryUse1,2))+" GB *****************\n")
###########################################################


########################################################################################
def getShortNames(chromNamesList):
    return [chrName[3:] for chrName in chromNamesList]
########################################################################################

########################################################################################
#Uses readDictionary
def getDictionary(outputDir,jobname, DictFilename):
    dictionary = {}
    DictFilePath = os.path.join(outputDir,jobname,DATA,DictFilename)

    if (os.path.exists(DictFilePath)):
        dictionary = readDictionary(DictFilePath)

    return dictionary
########################################################################################

########################################################################################
def getSample2NumberofIndelsDict(outputDir,jobname):
    sample2NumberofIndelsDict = {}
    sample2NumberofIndelsDictFilePath = os.path.join(outputDir,jobname,DATA,Sample2NumberofIndelsDictFilename)

    if (os.path.exists(sample2NumberofIndelsDictFilePath)):
        sample2NumberofIndelsDict = readDictionary(sample2NumberofIndelsDictFilePath)

    return sample2NumberofIndelsDict
########################################################################################


########################################################################################
def getSample2SubsSignature2NumberofMutationsDict(outputDir,jobname):
    sample2SubsSignature2NumberofMutationsDict= {}

    sample2SubsSignature2NumberofMutationsDictFilePath = os.path.join(outputDir,jobname,DATA,Sample2SubsSignature2NumberofMutationsDictFilename)

    if (os.path.exists(sample2SubsSignature2NumberofMutationsDictFilePath)):
        sample2SubsSignature2NumberofMutationsDict = readDictionary(sample2SubsSignature2NumberofMutationsDictFilePath)

    return sample2SubsSignature2NumberofMutationsDict
########################################################################################


########################################################################################
def getSample2NumberofSubsDict(outputDir,jobname):
    #Load samplesWithAtLeast10KMutations2NumberofMutationsDict
    sample2NumberofSubsDict = {}

    Sample2NumberofSubsDictFilePath = os.path.join(outputDir,jobname,DATA,Sample2NumberofSubsDictFilename)

    if (os.path.exists(Sample2NumberofSubsDictFilePath)):
        sample2NumberofSubsDict = readDictionary(Sample2NumberofSubsDictFilePath)

    return sample2NumberofSubsDict
########################################################################################

########################################################################################
def getSample2IndelsSignature2NumberofMutationsDict(outputDir,jobname):
    sample2IndelsSignature2NumberofMutationsDict = {}

    sample2IndelsSignature2NumberofMutationsDictFilePath = os.path.join(outputDir,jobname,DATA,Sample2IndelsSignature2NumberofMutationsDictFilename)

    if (os.path.exists(sample2IndelsSignature2NumberofMutationsDictFilePath)):
        sample2IndelsSignature2NumberofMutationsDict = readDictionary(sample2IndelsSignature2NumberofMutationsDictFilePath)

    return sample2IndelsSignature2NumberofMutationsDict
########################################################################################


###################################################################
def getChromSizesDict(genome):
    chromSizesDict = {}

    #TODO Do for mouse genomes
    if (genome==GRCh37):
        chromSizesDictPath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB,UCSCGENOME,GRCh37ChromSizesDictFilename)
    elif (genome==GRCh38):
        chromSizesDictPath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB,UCSCGENOME,GRCh38ChromSizesDictFilename)

    if (os.path.exists(chromSizesDictPath)):
        chromSizesDict = readDictionary(chromSizesDictPath)

    return chromSizesDict
###################################################################


###################################################################
def readTrancriptsENSEMBL(genome):

    if (genome==GRCh37):
        transcriptsFilenamePath = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,TRANSCRIPTS,GRCh37_ENSEMBL)

    if (os.path.exists(transcriptsFilenamePath)):
        ensembl_transcripts_df = pd.read_csv(transcriptsFilenamePath, header=0,sep='\t')

        #Gene stable ID  Transcript stable ID    Chromosome/scaffold name        Strand  Transcript start (bp)   Transcript end (bp)     Transcript type

        #Change the column name     Chromosome/scaffold name    -->     chrom
        #                           Strand                      -->     strand
        #                           Transcript start (bp)       -->     txStart
        #                           Transcript end (bp)         -->     txEnd
        ensembl_transcripts_df.rename(columns={'Chromosome/scaffold name': 'chrom', 'Strand': 'strand', 'Transcript start (bp)':'txStart', 'Transcript end (bp)':'txEnd'}, inplace=True)
        print('Chromosome names in transcripts data: %s' % (ensembl_transcripts_df['chrom'].unique()))

        # print(GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df.shape)
        # print(GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df.head())
        # print(GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df.tail())
        # print(GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df.info())

        return ensembl_transcripts_df
    else:
        return None
# Gene stable ID	Transcript stable ID	Chromosome/scaffold name	Strand	Transcript start (bp)	Transcript end (bp)	Transcript type
###################################################################

###################################################################
#Not used anymore
def readTranscriptsNCBI():
    transcriptsFilenamePath = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,TRANSCRIPTS,GRCh37_hg19_NCBIREFSEQCURATED)

    if (os.path.exists(transcriptsFilenamePath)):
        GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df = pd.read_csv(transcriptsFilenamePath, header = 0, sep='\t')
        print('debug GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df starts')
        print(GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df.columns.values)
        # print(GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df.shape)
        # print(GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df.head())
        # print(GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df.tail())
        # print(GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df.info())
        print('debug GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df ends')
        return GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df
    else:
        return None
##################################################################


##################################################################
def getSignatures(chrBased_mutation_df):
    # We assume that after the column named 'Mutation' there are the signature columns in tab separated way.
    columnNamesList = list(chrBased_mutation_df.columns.values)
    mutationIndex = columnNamesList.index(MUTATION)
    #Exclude the last one because it contains simulation_number column at the end.
    if (SIMULATION_NUMBER==columnNamesList[-1]):
        signatures = columnNamesList[(mutationIndex+1):-1]
    else:
        signatures = columnNamesList[(mutationIndex+1):]

    return signatures
##################################################################


##################################################################
class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NpEncoder, self).default(obj)
##################################################################

##################################################################
def writeUserFriendlySignatureBasedCutoff(outputDir,jobname,DATA,signature2PropertiesListDict,table_signature_cutoff_numberofmutations_averageprobability_filename):

    allcutoffs_signature_numberofmutations_averageprobability = open(os.path.join(outputDir,jobname,DATA, table_signature_cutoff_numberofmutations_averageprobability_filename), 'w')

    # header line
    allcutoffs_signature_numberofmutations_averageprobability.write('cancer_type\tsignature\tcutoff\tnumber_of_mutations\taverage_probability\n')

    sorted_signatures=sorted(signature2PropertiesListDict.keys(), key=natural_key)
    for signature in sorted_signatures:
        cutoff_numbeofmutations_averageprobability = signature2PropertiesListDict[signature]
        cutoff=cutoff_numbeofmutations_averageprobability[0]
        number_of_mutations = cutoff_numbeofmutations_averageprobability[1]
        average_probability = cutoff_numbeofmutations_averageprobability[2]
        allcutoffs_signature_numberofmutations_averageprobability.write('%s\t%s\t%s\t%s\t%s\n' %(jobname,signature,cutoff,number_of_mutations,average_probability))

    allcutoffs_signature_numberofmutations_averageprobability.close()
##################################################################


##################################################################
#two header lines
#rows will be cutoff
#columns will be for each signature number_of_mutations and average_probability
def writeUserFriendlyAllCutoffs(outputDir, jobname, DATA,cutoff2Signature2NumberofMutationsAverageProbabilityListDict,number_of_samples,table_allcutoffs_signature_numberofmutations_averageprobability_filename):
    signature_list=[]
    cutoff_list=sorted(cutoff2Signature2NumberofMutationsAverageProbabilityListDict.keys())

    for cutoff in cutoff2Signature2NumberofMutationsAverageProbabilityListDict:
        for signature in cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff]:
            if signature not in signature_list:
                signature_list.append(signature)

    signature_list=sorted(signature_list,key=natural_key)

    allcutoffs_signature_numberofmutations_averageprobability = open(os.path.join(outputDir,jobname,DATA, table_allcutoffs_signature_numberofmutations_averageprobability_filename), 'w')

    # 1st header line
    allcutoffs_signature_numberofmutations_averageprobability.write('number_of_samples\t%d\t' %(number_of_samples))
    for signature in signature_list:
        allcutoffs_signature_numberofmutations_averageprobability.write('%s\t%s\t' %(signature,signature))
    allcutoffs_signature_numberofmutations_averageprobability.write('\n')

    # 2nd header line
    allcutoffs_signature_numberofmutations_averageprobability.write('cancer_type\tcutoff\t')
    for signature in signature_list:
        allcutoffs_signature_numberofmutations_averageprobability.write('number_of_mutations\taverage_probability\t')
    allcutoffs_signature_numberofmutations_averageprobability.write('\n')

    for cutoff in cutoff_list:
        allcutoffs_signature_numberofmutations_averageprobability.write('%s\t%f\t' %(jobname,float(cutoff)))
        for signature in signature_list:
            if signature in cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff]:
                allcutoffs_signature_numberofmutations_averageprobability.write('%d\t%f\t' % (cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature][0],cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature][1]))
            else:
                allcutoffs_signature_numberofmutations_averageprobability.write('\t\t' )
        allcutoffs_signature_numberofmutations_averageprobability.write('\n')

    allcutoffs_signature_numberofmutations_averageprobability.close()
##################################################################



##################################################################
def fillCutoff2Signature2PropertiesListDictionary(outputDir,jobname,chromNamesList,type,cutoffs,average_probability,num_of_sbs_required,num_of_id_required,num_of_dbs_required,mutationType2NumberofMutationsDict):

    #Filled in the first part
    #PropertiesList consists of[sum_of_number of mutations, sum of probabilities]
    cutoff2Signature2PropertiesListDict={}

    #Filled in the second part
    # [number of mutations, average probability]
    cutoff2Signature2NumberofMutationsAverageProbabilityListDict={}

    #Filled in the third part
    #PropertiesList=[CufoffProbability NumberofMutations AverageMutationProbability]
    signature2PropertiesListDict={}

    all_samples=set()

    for chrLong in chromNamesList:
        chrbased_samples, chrBased_mutation_df = readChrBasedMutationsDF(outputDir,jobname,chrLong,type,0,return_number_of_samples=True)
        all_samples = all_samples.union(chrbased_samples)

        if ((chrBased_mutation_df is not None) and (not chrBased_mutation_df.empty)):
            # Sample  Chrom   Start   PyrimidineStrand        Mutation        DBS2    DBS4    DBS6    DBS7    DBS11
            # PD10011a        10      24033661        1       TC>AA   0.0     0.7656325053758131      0.15420390829468886     0.07918943063517644     0.000974155694321615
            signatures = getSignatures(chrBased_mutation_df)

            if type in mutationType2NumberofMutationsDict:
                mutationType2NumberofMutationsDict[type] += chrBased_mutation_df.shape[0]
            else:
                mutationType2NumberofMutationsDict[type] = chrBased_mutation_df.shape[0]

            # First part starts
            # First accumulate number of mutations and sum of probabilities with mutations with probability >= cutoff probability
            for cutoff in cutoffs:
                for signature in signatures:
                    # chrBased_mutation_df[signature]=chrBased_mutation_df[signature].astype(np.float64)
                    number_of_mutations = len(chrBased_mutation_df[chrBased_mutation_df[signature]>=float(cutoff)])
                    # This results in infinity
                    # sum_of_probabilities = (chrBased_mutation_df[chrBased_mutation_df[signature]>=float(cutoff)])[signature].sum()
                    sum_of_probabilities = np.sum(((chrBased_mutation_df[chrBased_mutation_df[signature]>=float(cutoff)])[signature]).values,dtype=np.float64)

                    if cutoff not in cutoff2Signature2PropertiesListDict:
                        cutoff2Signature2PropertiesListDict[cutoff]={}
                        cutoff2Signature2PropertiesListDict[cutoff][signature]=[]
                        cutoff2Signature2PropertiesListDict[cutoff][signature].append(np.int64(0))
                        cutoff2Signature2PropertiesListDict[cutoff][signature].append(np.float64(0.0))
                        cutoff2Signature2PropertiesListDict[cutoff][signature][0]+=number_of_mutations
                        cutoff2Signature2PropertiesListDict[cutoff][signature][1]+=sum_of_probabilities
                    elif signature not in cutoff2Signature2PropertiesListDict[cutoff]:
                        cutoff2Signature2PropertiesListDict[cutoff][signature]=[]
                        cutoff2Signature2PropertiesListDict[cutoff][signature].append(np.int64(0))
                        cutoff2Signature2PropertiesListDict[cutoff][signature].append(np.float64(0.0))
                        cutoff2Signature2PropertiesListDict[cutoff][signature][0]+=number_of_mutations
                        cutoff2Signature2PropertiesListDict[cutoff][signature][1]+=sum_of_probabilities
                    else:
                        cutoff2Signature2PropertiesListDict[cutoff][signature][0]+=number_of_mutations
                        cutoff2Signature2PropertiesListDict[cutoff][signature][1]+=sum_of_probabilities
        #First part ends
        #Accumulation ended for each chromosome


    #Second part starts
    #Second find average probability
    for cutoff in cutoff2Signature2PropertiesListDict:
        for signature in cutoff2Signature2PropertiesListDict[cutoff]:
            if (cutoff2Signature2PropertiesListDict[cutoff][signature][0]>0):
                if cutoff not in cutoff2Signature2NumberofMutationsAverageProbabilityListDict:
                    cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff]={}
                    cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature] = []
                    cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature].append(np.int64(cutoff2Signature2PropertiesListDict[cutoff][signature][0]))
                    cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature].append(np.float64(cutoff2Signature2PropertiesListDict[cutoff][signature][1]/cutoff2Signature2PropertiesListDict[cutoff][signature][0]))
                elif signature not in cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff]:
                    cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature] = []
                    cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature].append(np.int64(cutoff2Signature2PropertiesListDict[cutoff][signature][0]))
                    cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature].append(np.float64(cutoff2Signature2PropertiesListDict[cutoff][signature][1]/cutoff2Signature2PropertiesListDict[cutoff][signature][0]))
                else:
                    print('There is a situation in fillCutoff2Signature2PropertiesListDictionary method.')
                    cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature].append(np.int64(cutoff2Signature2PropertiesListDict[cutoff][signature][0]))
                    cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature].append(np.float64(cutoff2Signature2PropertiesListDict[cutoff][signature][1]/cutoff2Signature2PropertiesListDict[cutoff][signature][0]))
    #Second part ends


    #Set the filenames and number of required mutations
    if (type==SUBS):
        number_of_required_mutations=num_of_sbs_required
        table_allcutoffs_signature_numberofmutations_averageprobability_filename = Table_AllCutoff_SubsSignature_NumberofMutations_AverageProbability_Filename
        table_signature_cutoff_numberofmutations_averageprobability_filename=Table_SubsSignature_Cutoff_NumberofMutations_AverageProbability_Filename
    elif (type == INDELS):
        number_of_required_mutations=num_of_id_required
        table_allcutoffs_signature_numberofmutations_averageprobability_filename = Table_AllCutoff_IndelsSignature_NumberofMutations_AverageProbability_Filename
        table_signature_cutoff_numberofmutations_averageprobability_filename=Table_IndelsSignature_Cutoff_NumberofMutations_AverageProbability_Filename
    elif (type== DINUCS):
        number_of_required_mutations=num_of_dbs_required
        table_allcutoffs_signature_numberofmutations_averageprobability_filename = Table_AllCutoff_DinucsSignature_NumberofMutations_AverageProbability_Filename
        table_signature_cutoff_numberofmutations_averageprobability_filename=Table_DinucsSignature_Cutoff_NumberofMutations_AverageProbability_Filename


    ####################################################################
    number_of_samples=len(all_samples)
    writeUserFriendlyAllCutoffs(outputDir,jobname,DATA,cutoff2Signature2NumberofMutationsAverageProbabilityListDict,number_of_samples,table_allcutoffs_signature_numberofmutations_averageprobability_filename)
    ####################################################################

    #Third find the signature based cufoff probability with number of mutations >= required number of mutations  and averega mutation probability >=0.9
    sorted_cutoffs=sorted(cutoff2Signature2NumberofMutationsAverageProbabilityListDict.keys())
    for cutoff in sorted_cutoffs:
        for signature in cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff]:
            if signature not in signature2PropertiesListDict:
                if (cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature][0]>=number_of_required_mutations) and (cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature][1]>=average_probability):
                    signature2PropertiesListDict[signature]=[cutoff,np.int(cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature][0]),np.float(cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature][1])]
    #Third part ends


    ####################################################################
    #Update signature2PropertiesListDict with signatures not satisfying our constraints
    #There are two options

    #Option1: Use the lowest existing cutoff with the highest number of mutations
    # sorted_cutoffs=sorted(cutoff2Signature2NumberofMutationsAverageProbabilityListDict.keys())
    # for cutoff in sorted_cutoffs:
    #     for signature in cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff]:
    #         if signature not in signature2PropertiesListDict:
    #             signature2PropertiesListDict[signature]=[cutoff,np.int(cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature][0]),np.float(cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature][1])]

    #Option2: Use the highest existing cutoff with the required number of mutations
    # sorted_cutoffs=sorted(cutoff2Signature2NumberofMutationsAverageProbabilityListDict.keys(),reverse=True)
    # if (len(sorted_cutoffs)>0):
    #     min_cutoff=sorted_cutoffs[-1]
    #     signatures=cutoff2Signature2NumberofMutationsAverageProbabilityListDict[min_cutoff].keys()
    #
    #     for signature in signatures:
    #         if signature not in signature2PropertiesListDict:
    #             for cutoff in sorted_cutoffs:
    #                 if (cutoff in cutoff2Signature2NumberofMutationsAverageProbabilityListDict) and (signature in cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff]):
    #                     if (cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature][0]>=number_of_required_mutations):
    #                         signature2PropertiesListDict[signature]=[cutoff,np.int(cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature][0]),np.float(cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature][1])]
    #                         break
    #         #If not found apply option1
    #         if signature not in signature2PropertiesListDict:
    #             signature2PropertiesListDict[signature] = [min_cutoff, np.int(cutoff2Signature2NumberofMutationsAverageProbabilityListDict[min_cutoff][signature][0]),np.float(cutoff2Signature2NumberofMutationsAverageProbabilityListDict[min_cutoff][signature][1])]
    ####################################################################

    #Write the dictionaries
    os.makedirs(os.path.join(outputDir,jobname,DATA), exist_ok=True)

    ####################################################################
    writeUserFriendlySignatureBasedCutoff(outputDir,jobname,DATA,signature2PropertiesListDict,table_signature_cutoff_numberofmutations_averageprobability_filename)
    ####################################################################


    # ####################################################################
    # #Write dictionary as a dataframe starts
    # filePath = os.path.join(outputDir, jobname, DATA, table_allcutoffs_signature_numberofmutations_averageprobability_filename)
    #
    # L = sorted([(cutoff, signature, numberofmutations_averageprobability_list[0],numberofmutations_averageprobability_list[1])
    #             for cutoff, a in cutoff2Signature2NumberofMutationsAverageProbabilityListDict.items()
    #              for signature, numberofmutations_averageprobability_list in a.items()])
    # df = pd.DataFrame(L,columns=['cutoff','signature', 'number_of_mutations','average_probability'])
    #
    # #write this dataframe
    # df.to_csv(filePath, sep='\t', header=True, index=False)
    # #Write dictionary as a dataframe ends
    # ####################################################################


    ####################################################################
    L = sorted([(signature, cutoff_numberofmutations_averageprobability[0],cutoff_numberofmutations_averageprobability[1],cutoff_numberofmutations_averageprobability[2])
                 for signature, cutoff_numberofmutations_averageprobability in signature2PropertiesListDict.items()])

    signature_cutoff_numberofmutations_averageprobability_df = pd.DataFrame(L,columns=['signature','cutoff', 'number_of_mutations','average_probability'])

    signature_cutoff_numberofmutations_averageprobability_df['signature']=signature_cutoff_numberofmutations_averageprobability_df['signature'].astype(str)
    signature_cutoff_numberofmutations_averageprobability_df['cutoff']=signature_cutoff_numberofmutations_averageprobability_df['cutoff'].astype(np.float32)
    signature_cutoff_numberofmutations_averageprobability_df['number_of_mutations']=signature_cutoff_numberofmutations_averageprobability_df['number_of_mutations'].astype(np.int32)
    signature_cutoff_numberofmutations_averageprobability_df['average_probability'] = signature_cutoff_numberofmutations_averageprobability_df['average_probability'].astype(np.float32)
    ####################################################################

    return signature_cutoff_numberofmutations_averageprobability_df
##################################################################


##################################################################
# We are filling and writing
# Sample2NumberofMutationsDictFilename
# Signature2NumberofMutationsDictFilename
# Sample2Signature2NumberofMutationsDictFilename

# We are writing Signature2NumberofMutationsDictFilename for the signatures in signature2PropertiesListDict using chrBased_mutation_df
# At the end, Signature2NumberofMutationsDictFilename and signature2PropertiesListDict must match
# It is like double check
def fill_mutations_dictionaries_write(outputDir, jobname, chromNamesList, type, signature_cutoff_numberofmutations_averageprobability_df,num_of_sbs_required,num_of_id_required,num_of_dbs_required):
    sample2NumberofMutationsDict = {}
    signature2NumberofMutationsDict = {}
    sample2Signature2NumberofMutationsDict = {}

    #Fill dictionaries with conditions satisfied
    if (type==SUBS):
        Sample2NumberofMutationsDictFilename = Sample2NumberofSubsDictFilename
        Signature2NumberofMutationsDictFilename = SubsSignature2NumberofMutationsDictFilename
        Sample2Signature2NumberofMutationsDictFilename = Sample2SubsSignature2NumberofMutationsDictFilename
        minimum_number_of_mutations_required = num_of_sbs_required
    elif (type==INDELS):
        Sample2NumberofMutationsDictFilename = Sample2NumberofIndelsDictFilename
        Signature2NumberofMutationsDictFilename = IndelsSignature2NumberofMutationsDictFilename
        Sample2Signature2NumberofMutationsDictFilename = Sample2IndelsSignature2NumberofMutationsDictFilename
        minimum_number_of_mutations_required = num_of_id_required
    elif (type==DINUCS):
        Sample2NumberofMutationsDictFilename = Sample2NumberofDinucsDictFilename
        Signature2NumberofMutationsDictFilename = DinucsSignature2NumberofMutationsDictFilename
        Sample2Signature2NumberofMutationsDictFilename = Sample2DinucsSignature2NumberofMutationsDictFilename
        minimum_number_of_mutations_required = num_of_dbs_required

    #######################################################################
    for chrLong in chromNamesList:
        chrBased_mutation_df = readChrBasedMutationsDF(outputDir,jobname,chrLong,type,0)

        if ((chrBased_mutation_df is not None) and (not chrBased_mutation_df.empty)):
            # Sample  Chrom   Start   PyrimidineStrand        Mutation        DBS2    DBS4    DBS6    DBS7    DBS11
            # PD10011a        10      24033661        1       TC>AA   0.0     0.7656325053758131      0.15420390829468886     0.07918943063517644     0.000974155694321615
            #old way
            # signatures = getSignatures(chrBased_mutation_df)
            chrBased_mutation_df_sample_grouped = chrBased_mutation_df.groupby('Sample')

            for sample, chrBased_mutation_df_sample_group_df in chrBased_mutation_df_sample_grouped:
                number_of_mutations = chrBased_mutation_df_sample_group_df.shape[0]
                if sample in sample2NumberofMutationsDict:
                    sample2NumberofMutationsDict[sample] += number_of_mutations
                else:
                    sample2NumberofMutationsDict[sample] = number_of_mutations

                #new way
                for signature in signature_cutoff_numberofmutations_averageprobability_df['signature'].unique():
                    cutoff=float(signature_cutoff_numberofmutations_averageprobability_df[signature_cutoff_numberofmutations_averageprobability_df['signature']==signature]['cutoff'].values[0])
                    number_of_mutations= len(chrBased_mutation_df_sample_group_df[chrBased_mutation_df_sample_group_df[signature]>=cutoff])
                    if signature in signature2NumberofMutationsDict:
                        signature2NumberofMutationsDict[signature] += number_of_mutations
                    else:
                        signature2NumberofMutationsDict[signature] = number_of_mutations
                    if sample in sample2Signature2NumberofMutationsDict:
                        if signature in sample2Signature2NumberofMutationsDict[sample]:
                            sample2Signature2NumberofMutationsDict[sample][signature] += number_of_mutations
                        else:
                            sample2Signature2NumberofMutationsDict[sample][signature] = number_of_mutations
                    else:
                        sample2Signature2NumberofMutationsDict[sample] = {}
                        sample2Signature2NumberofMutationsDict[sample][signature] = number_of_mutations

    #######################################################################

    new_sample2NumberofMutatiosDict = {k: v for k, v in sample2NumberofMutationsDict.items() if sample2NumberofMutationsDict[k] >= minimum_number_of_mutations_required}
    new_signature2NumberofMutationsDict = {k: v for k, v in signature2NumberofMutationsDict.items() if signature2NumberofMutationsDict[k] >= minimum_number_of_mutations_required}

    new_sample2Signature2NumberofMutationsDict = {}
    for k1, v1 in sample2Signature2NumberofMutationsDict.items():
        for k2, v2 in sample2Signature2NumberofMutationsDict[k1].items():
            if (sample2Signature2NumberofMutationsDict[k1][k2] >= minimum_number_of_mutations_required):
                if k1 in new_sample2Signature2NumberofMutationsDict:
                    new_sample2Signature2NumberofMutationsDict[k1][k2] = v2
                else:
                    new_sample2Signature2NumberofMutationsDict[k1] = {}
                    new_sample2Signature2NumberofMutationsDict[k1][k2] = v2

    #Write Conditions satisfied Dictionaries
    writeDictionaryUnderDataDirectory(new_sample2NumberofMutatiosDict,outputDir,jobname,Sample2NumberofMutationsDictFilename)
    writeDictionaryUnderDataDirectory(new_signature2NumberofMutationsDict,outputDir,jobname,Signature2NumberofMutationsDictFilename)
    writeDictionaryUnderDataDirectory(new_sample2Signature2NumberofMutationsDict,outputDir,jobname,Sample2Signature2NumberofMutationsDictFilename)
##################################################################

##################################################################
#Used for all kinds of mutations SUBs and Dinucs and Indels
def readChrBasedMutationsDF(outputDir,jobname,chrLong,type,simulationNumber,return_number_of_samples=False):
    filename = '%s_%s_for_topography.txt' %(chrLong,type)

    if (simulationNumber==0):
        chrBasedMutationDFFilePath = os.path.join(outputDir,jobname,DATA,CHRBASED,filename)
    else:
        simulation = 'sim%s' % (simulationNumber)
        chrBasedMutationDFFilePath = os.path.join(outputDir,jobname,DATA,CHRBASED,simulation,filename)

    chrBased_mutation_df = None
    #############################################
    if (os.path.exists(chrBasedMutationDFFilePath)):

        #############################################
        only_header_chrBased_mutation_df = pd.read_csv(chrBasedMutationDFFilePath, sep='\t', comment='#', nrows=1)
        columnNamesList = list(only_header_chrBased_mutation_df.columns.values)

        # Please note that we assume that after the column named 'Mutation' there are the signature columns in tab separated way.
        mutationtIndex = columnNamesList.index(MUTATION)
        signatures = columnNamesList[(mutationtIndex + 1):]
        #############################################

        #################################################
        mydtypes = {}

        # Subs
        # Sample  Chrom   Start   MutationLong    PyramidineStrand        TranscriptionStrand     Mutation        SBS1    SBS2    SBS3    SBS4    SBS5    SBS6    SBS7a   SBS7b   SBS7c   SBS7d   SBS8    SBS9    SBS10a  SBS10b  SBS11   SBS12   SBS13   SBS14   SBS15   SBS16   SBS17a  SBS17b  SBS18   SBS19   SBS20   SBS21   SBS22
        # Indels
        # Sample  Chrom   Start   MutationLong    Ref     Alt     Length  PyramidineStrand        TranscriptionStrand     Mutation        ID1     ID2     ID3     ID4     ID5     ID6     ID7     ID8     ID9     ID10    ID11    ID12    ID13    ID14    ID15    ID16    ID17
        # Dinucs
        # Sample  Chrom   Start   MutationLong    PyramidineStrand        TranscriptionStrand     Mutation        DBS1    DBS2    DBS3    DBS4    DBS5    DBS6    DBS7    DBS8    DBS9    DBS10   DBS11

        #np.float16 Half precision float: sign bit, 5 bits exponent, 10 bits mantissa
        #np.int8 Byte (-128 to 127)
        #np.int16 Integer (-32768 to 32767)
        #np.int32   Integer (-2147483648 to 2147483647)

        #To lower the dataframe size
        for signature in signatures:
            mydtypes[signature] = np.float16

        if ((type==SUBS) or (type==DINUCS)):
            mydtypes[SAMPLE] = 'category'
            mydtypes[CHROM] = 'category'
            mydtypes[START] = np.int32
            mydtypes[MUTATIONLONG] = 'category'
            mydtypes[PYRAMIDINESTRAND] = np.int8
            mydtypes[TRANSCRIPTIONSTRAND] = 'category'
            mydtypes[MUTATION] = 'category'

        if (type==INDELS):
            mydtypes[SAMPLE] = 'category'
            mydtypes[CHROM] = 'category'
            mydtypes[START] = np.int32
            mydtypes[MUTATIONLONG] = 'category'
            mydtypes[REF]=str
            mydtypes[ALT]=str
            mydtypes[LENGTH]=np.int16
            mydtypes[PYRAMIDINESTRAND] = np.int8
            mydtypes[TRANSCRIPTIONSTRAND] = 'category'
            mydtypes[MUTATION] = 'category'
        #################################################

        chrBased_mutation_df = pd.read_csv(chrBasedMutationDFFilePath,sep='\t', header=0, dtype=mydtypes)
        chrBased_mutation_df[SIMULATION_NUMBER]=simulationNumber
    #############################################

    if return_number_of_samples:
        if (chrBased_mutation_df is not None):
            chrbased_samples = chrBased_mutation_df[SAMPLE].unique()
            return chrbased_samples, chrBased_mutation_df
        else:
            return set(), chrBased_mutation_df
    else:
        return chrBased_mutation_df
##################################################################

##################################################################
def doesSimulationsAlreadyExists(outputDir,jobname,numofSimulations):
    if (numofSimulations>0):
        for simNum in range(1,numofSimulations+1):
            sim='sim%d' %(simNum)
            simDir=os.path.join(outputDir,jobname,DATA,CHRBASED,sim)
            if os.path.exists(simDir):
                numberoffiles=len([name for name in os.listdir(simDir) if os.path.isfile(os.path.join(simDir,name))])
                if numberoffiles==0:
                    return False
            else:
                return False
        return True
    else:
        return False
##################################################################


##########################################################################################
############### Common functions for Nucleosome Occupancy Analysis starts ################
##########################################################################################

########################################################################################
#We will accumulate signature2SplitArrayDict in signature2AccumulatedSplitsChrBasedArrayDict
def accumulateSimulationBasedTypeBasedArrays(simNum2Type2AccumulatedSplitsChrBasedArrayDict, simNum2Type2SplitArrayDict):
    for simNum in simNum2Type2SplitArrayDict:
        type2SplitArrayDict = simNum2Type2SplitArrayDict[simNum]

        if simNum in simNum2Type2AccumulatedSplitsChrBasedArrayDict:
            type2AccumulatedSplitsChrBasedArrayDict = simNum2Type2AccumulatedSplitsChrBasedArrayDict[simNum]
        else:
            type2AccumulatedSplitsChrBasedArrayDict = {}
            simNum2Type2AccumulatedSplitsChrBasedArrayDict[simNum] = type2AccumulatedSplitsChrBasedArrayDict

        for type in type2SplitArrayDict.keys():
            if type in type2AccumulatedSplitsChrBasedArrayDict:
                type2AccumulatedSplitsChrBasedArrayDict[type] += type2SplitArrayDict[type]
            else:
                type2AccumulatedSplitsChrBasedArrayDict[type] = type2SplitArrayDict[type]
########################################################################################

########################################################################################
def accumulateSimulationBasedSampleBasedTypeBasedArrays(simNum2Sample2Type2AccumulatedSplitsChrBasedArrayDict,simNum2Sample2Type2SplitArrayDict):
    for simNum in simNum2Sample2Type2SplitArrayDict:
        sample2Type2SplitArrayDict = simNum2Sample2Type2SplitArrayDict[simNum]
        if simNum in simNum2Sample2Type2AccumulatedSplitsChrBasedArrayDict:
            sample2Type2AccumulatedSplitsChrBasedArrayDict = simNum2Sample2Type2AccumulatedSplitsChrBasedArrayDict[simNum]
        else:
            sample2Type2AccumulatedSplitsChrBasedArrayDict= {}
            simNum2Sample2Type2AccumulatedSplitsChrBasedArrayDict[simNum] = sample2Type2AccumulatedSplitsChrBasedArrayDict

        for sample in sample2Type2SplitArrayDict.keys():
            for type in sample2Type2SplitArrayDict[sample].keys():
                if (sample in sample2Type2AccumulatedSplitsChrBasedArrayDict):
                    if (type in sample2Type2AccumulatedSplitsChrBasedArrayDict[sample]):
                        sample2Type2AccumulatedSplitsChrBasedArrayDict[sample][type] += sample2Type2SplitArrayDict[sample][type]
                    else:
                        sample2Type2AccumulatedSplitsChrBasedArrayDict[sample][type] = sample2Type2SplitArrayDict[sample][type]
                else:
                    sample2Type2AccumulatedSplitsChrBasedArrayDict[sample]= {}
                    sample2Type2AccumulatedSplitsChrBasedArrayDict[sample][type] = sample2Type2SplitArrayDict[sample][type]
########################################################################################


########################################################################################
#window_array is of size 2*plusOrMinus
# mutation_row_start will be at position=plusOrMinus of the window_array
def func_addSignal(window_array, entry_start, entry_end, entry_signal, mutation_row_start,plusOrMinus):
    max_start=max(entry_start,mutation_row_start-plusOrMinus)
    min_end=min(entry_end,mutation_row_start+plusOrMinus)
    window_array[max_start-(mutation_row_start-plusOrMinus):min_end-(mutation_row_start-plusOrMinus)+1]+=entry_signal
    # print('window_array[%d:%d]+=%f' %(max_start-(mutation_row_start-plusOrMinus),min_end-(mutation_row_start-plusOrMinus),entry_signal))
########################################################################################


########################################################################################
def computeAverageNucleosomeOccupancyArray(plusorMinus,signalArray,countArray):
    windowSize=plusorMinus*2+1
    averageArray =  np.zeros(windowSize)

    np.seterr(divide='ignore', invalid='ignore')
    if (np.any(countArray)):
        averageArray = np.divide(signalArray,countArray)
    np.seterr(divide='raise', invalid='ignore')

    #October 27, 2019
    #Assume that there is no signal 0
    #Assume that there is no count 0
    #Then average signal must be 0
    #If we want to know that there is no signal there we can still keep nan values.
    # We can discard nans by np.nanmean()
    #Since later on, conversion of nans into zeros may lead to unrealistic fold changes between original data and simulations
    # averageArray[np.isnan(averageArray)] = 0

    return averageArray
########################################################################################

########################################################################################
def writeSimulationBasedAverageNucleosomeOccupancy(
        occupancy_type,
        sample_based,
        plusorMinus,
        simNum2Type2AccumulatedSignalArrayDict,
        simNum2Type2AccumulatedCountArrayDict,
        simNum2Sample2Type2AccumulatedSignalArrayDict,
        simNum2Sample2Type2AccumulatedCountArrayDict,
        outputDir,
        jobname,
        library_file_memo,verbose):

    os.makedirs(os.path.join(outputDir, jobname, DATA, occupancy_type),exist_ok=True)

    ####################################################################################
    for simNum in simNum2Type2AccumulatedSignalArrayDict:
        type2AccumulatedSignalArrayDict= simNum2Type2AccumulatedSignalArrayDict[simNum]
        type2AccumulatedCountArrayDict = simNum2Type2AccumulatedCountArrayDict[simNum]
        if (AGGREGATEDSUBSTITUTIONS in type2AccumulatedSignalArrayDict):
            writeAverageNucleosomeOccupancyFiles(occupancy_type,plusorMinus,type2AccumulatedSignalArrayDict[AGGREGATEDSUBSTITUTIONS],type2AccumulatedCountArrayDict[AGGREGATEDSUBSTITUTIONS],outputDir, jobname,library_file_memo, AGGREGATEDSUBSTITUTIONS,simNum)
        if (AGGREGATEDINDELS in type2AccumulatedSignalArrayDict):
            writeAverageNucleosomeOccupancyFiles(occupancy_type,plusorMinus,type2AccumulatedSignalArrayDict[AGGREGATEDINDELS],type2AccumulatedCountArrayDict[AGGREGATEDINDELS], outputDir,jobname,library_file_memo, AGGREGATEDINDELS,simNum)
        if (AGGREGATEDDINUCS in type2AccumulatedSignalArrayDict):
            writeAverageNucleosomeOccupancyFiles(occupancy_type,plusorMinus,type2AccumulatedSignalArrayDict[AGGREGATEDDINUCS],type2AccumulatedCountArrayDict[AGGREGATEDDINUCS], outputDir,jobname,library_file_memo, AGGREGATEDDINUCS,simNum)
        writeSignatureBasedAverageNucleosomeOccupancyFiles(occupancy_type,plusorMinus,type2AccumulatedSignalArrayDict,type2AccumulatedCountArrayDict,outputDir,jobname,library_file_memo,simNum)
    ####################################################################################


    ####################################################################################
    if sample_based:
        for simNum in simNum2Sample2Type2AccumulatedSignalArrayDict:
            sample2Type2AccumulatedSignalArrayDict = simNum2Sample2Type2AccumulatedSignalArrayDict[simNum]
            sample2Type2AccumulatedCountArrayDict = simNum2Sample2Type2AccumulatedCountArrayDict[simNum]
            ####################################################################################
            writeSampleBasedAverageNucleosomeOccupancyFiles(occupancy_type,plusorMinus,sample2Type2AccumulatedSignalArrayDict,sample2Type2AccumulatedCountArrayDict,outputDir,jobname,library_file_memo,AGGREGATEDSUBSTITUTIONS,simNum)
            writeSampleBasedAverageNucleosomeOccupancyFiles(occupancy_type,plusorMinus,sample2Type2AccumulatedSignalArrayDict,sample2Type2AccumulatedCountArrayDict,outputDir,jobname,library_file_memo,AGGREGATEDINDELS,simNum)
            writeSampleBasedAverageNucleosomeOccupancyFiles(occupancy_type,plusorMinus,sample2Type2AccumulatedSignalArrayDict,sample2Type2AccumulatedCountArrayDict, outputDir,jobname,library_file_memo, AGGREGATEDDINUCS,simNum)
            writeSampleBasedSignatureBasedAverageNucleosomeOccupancyFiles(occupancy_type,plusorMinus,sample2Type2AccumulatedSignalArrayDict,sample2Type2AccumulatedCountArrayDict,outputDir,jobname,library_file_memo,simNum)
    ####################################################################################

########################################################################################



########################################################################################
#Both "all single point mutations" and "all indels" use this function
#simulationNumber == 0 means original data
#simulationNumber > 0 means simulation data
def writeAverageNucleosomeOccupancyFiles(occupancy_type,plusorMinus,allMutationsAccumulatedAllChromsSignalArray,allMutationsAccumulatedAllChromsCountArray,outputDir,jobname,library_file_memo,nucleosomeOccupancyAnalysisType,simulationNumber):

    os.makedirs(os.path.join(outputDir, jobname, DATA, occupancy_type, nucleosomeOccupancyAnalysisType),exist_ok=True)
    averageNucleosomeSignalArray = computeAverageNucleosomeOccupancyArray(plusorMinus,allMutationsAccumulatedAllChromsSignalArray, allMutationsAccumulatedAllChromsCountArray)


    if (simulationNumber==0):
        if library_file_memo is not None:
            accumulatedSignalFilename = '%s_%s_AccumulatedSignalArray.txt' % (jobname,library_file_memo)
            accumulatedCountFilename = '%s_%s_AccumulatedCountArray.txt' % (jobname,library_file_memo)
            averageNucleosomeSignalFilename = '%s_%s_AverageSignalArray.txt' % (jobname,library_file_memo)
        else:
            accumulatedSignalFilename = '%s_AccumulatedSignalArray.txt' %(jobname)
            accumulatedCountFilename = '%s_AccumulatedCountArray.txt' %(jobname)
            averageNucleosomeSignalFilename = '%s_AverageSignalArray.txt' %(jobname)
    else:
        if library_file_memo is not None:
            accumulatedSignalFilename = '%s_sim%d_%s_AccumulatedSignalArray.txt' %(jobname,simulationNumber,library_file_memo)
            accumulatedCountFilename = '%s_sim%d_%s_AccumulatedCountArray.txt' %(jobname,simulationNumber,library_file_memo)
            averageNucleosomeSignalFilename = '%s_sim%d_%s_AverageSignalArray.txt' %(jobname,simulationNumber,library_file_memo)

        else:
            accumulatedSignalFilename = '%s_sim%d_AccumulatedSignalArray.txt' %(jobname,simulationNumber)
            accumulatedCountFilename = '%s_sim%d_AccumulatedCountArray.txt' %(jobname,simulationNumber)
            averageNucleosomeSignalFilename = '%s_sim%d_AverageSignalArray.txt' %(jobname,simulationNumber)

    accumulatedSignalFilePath = os.path.join(outputDir, jobname,DATA, occupancy_type, nucleosomeOccupancyAnalysisType, accumulatedSignalFilename)
    accumulatedCountFilePath = os.path.join(outputDir, jobname, DATA,occupancy_type, nucleosomeOccupancyAnalysisType, accumulatedCountFilename)
    averageNucleosomeSignalFilePath = os.path.join(outputDir,jobname, DATA, occupancy_type, nucleosomeOccupancyAnalysisType,averageNucleosomeSignalFilename)

    allMutationsAccumulatedAllChromsSignalArray.tofile(file=accumulatedSignalFilePath, sep="\t", format="%s")
    allMutationsAccumulatedAllChromsCountArray.tofile(file=accumulatedCountFilePath, sep="\t", format="%s")
    averageNucleosomeSignalArray.tofile(file=averageNucleosomeSignalFilePath, sep="\t", format="%s")
########################################################################################

########################################################################################
def writeSignatureBasedAverageNucleosomeOccupancyFiles(occupancy_type,plusorMinus,type2AccumulatedAllChromsSignalArrayDict,type2AccumulatedAllChromsCountArrayDict,outputDir,jobname,library_file_memo,simulationNumber):

    os.makedirs(os.path.join(outputDir,jobname,DATA,occupancy_type,SIGNATUREBASED), exist_ok=True)

    for type in type2AccumulatedAllChromsSignalArrayDict.keys():
        if ((type!=AGGREGATEDSUBSTITUTIONS) and (type!=AGGREGATEDINDELS) and (type!=AGGREGATEDDINUCS)):
            signalArray = type2AccumulatedAllChromsSignalArrayDict[type]
            countArray = type2AccumulatedAllChromsCountArrayDict[type]
            averageNucleosomeSignalArray = computeAverageNucleosomeOccupancyArray(plusorMinus,signalArray,countArray)

            #To provide filename with no space in signature name
            #signatureWithNoSpace = signature.replace(' ','')

            if (simulationNumber==0):
                if library_file_memo is not None:
                    accumulatedSignalFilename = '%s_%s_AccumulatedSignalArray.txt' % (type,library_file_memo)
                    accumulatedCountFilename = '%s_%s_AccumulatedCountArray.txt' % (type,library_file_memo)
                    averageNucleosomeSignalFilename = '%s_%s_AverageSignalArray.txt' % (type,library_file_memo)
                else:
                    accumulatedSignalFilename = '%s_AccumulatedSignalArray.txt' %(type)
                    accumulatedCountFilename = '%s_AccumulatedCountArray.txt' %(type)
                    averageNucleosomeSignalFilename = '%s_AverageSignalArray.txt' %(type)
            else:
                if library_file_memo is not None:
                    accumulatedSignalFilename = '%s_sim%d_%s_AccumulatedSignalArray.txt' % (type, simulationNumber,library_file_memo)
                    accumulatedCountFilename = '%s_sim%d_%s_AccumulatedCountArray.txt' % (type, simulationNumber,library_file_memo)
                    averageNucleosomeSignalFilename = '%s_sim%d_%s_AverageSignalArray.txt' % (type, simulationNumber,library_file_memo)
                else:
                    accumulatedSignalFilename = '%s_sim%d_AccumulatedSignalArray.txt' %(type,simulationNumber)
                    accumulatedCountFilename = '%s_sim%d_AccumulatedCountArray.txt' %(type,simulationNumber)
                    averageNucleosomeSignalFilename = '%s_sim%d_AverageSignalArray.txt' %(type,simulationNumber)

            accumulatedSignalFilePath = os.path.join(outputDir,jobname,DATA,occupancy_type,SIGNATUREBASED,accumulatedSignalFilename)
            accumulatedCountFilePath = os.path.join(outputDir,jobname,DATA,occupancy_type,SIGNATUREBASED,accumulatedCountFilename)
            averageNucleosomeSignalFilePath = os.path.join(outputDir,jobname,DATA,occupancy_type,SIGNATUREBASED,averageNucleosomeSignalFilename)

            signalArray.tofile(file=accumulatedSignalFilePath, sep="\t",format="%s")
            countArray.tofile(file=accumulatedCountFilePath, sep="\t", format="%s")
            averageNucleosomeSignalArray.tofile(file=averageNucleosomeSignalFilePath, sep="\t", format="%s")
########################################################################################


########################################################################################
def writeSampleBasedSignatureBasedAverageNucleosomeOccupancyFiles(occupancy_type,plusorMinus,sample2Type2AccumulatedAllChromsSignalArrayDict,
                                                                sample2Type2AccumulatedAllChromsCountArrayDict,outputDir,jobname,library_file_memo,simulationNumber):


    for sample in sample2Type2AccumulatedAllChromsSignalArrayDict:
        for type in sample2Type2AccumulatedAllChromsSignalArrayDict[sample]:
            if ((type!=AGGREGATEDSUBSTITUTIONS) and (type!=AGGREGATEDINDELS) and (type!=AGGREGATEDDINUCS)):
                signalArray = sample2Type2AccumulatedAllChromsSignalArrayDict[sample][type]
                countArray = sample2Type2AccumulatedAllChromsCountArrayDict[sample][type]
                averageNucleosomeSignalArray = computeAverageNucleosomeOccupancyArray(plusorMinus,signalArray, countArray)

                if (simulationNumber==0):
                    if (library_file_memo is not None):
                        accumulatedSignalFilename = '%s_%s_%s_AccumulatedSignalArray.txt' %(type,sample,library_file_memo)
                        accumulatedCountFilename = '%s_%s_%s_AccumulatedCountArray.txt' %(type,sample,library_file_memo)
                        averageNucleosomeSignalFilename = '%s_%s_%s_AverageSignalArray.txt' %(type,sample,library_file_memo)
                    else:
                        accumulatedSignalFilename = '%s_%s_AccumulatedSignalArray.txt' % (type, sample)
                        accumulatedCountFilename = '%s_%s_AccumulatedCountArray.txt' % (type, sample)
                        averageNucleosomeSignalFilename = '%s_%s_AverageSignalArray.txt' % (type, sample)
                else:
                    if (library_file_memo is not None):
                        accumulatedSignalFilename = '%s_%s_sim%d_%s_AccumulatedSignalArray.txt' %(type,sample,simulationNumber,library_file_memo)
                        accumulatedCountFilename = '%s_%s_sim%d_%s_AccumulatedCountArray.txt' %(type,sample,simulationNumber,library_file_memo)
                        averageNucleosomeSignalFilename = '%s_%s_sim%d_%s_AverageSignalArray.txt' %(type,sample,simulationNumber,library_file_memo)

                    else:
                        accumulatedSignalFilename = '%s_%s_sim%d_AccumulatedSignalArray.txt' %(type,sample,simulationNumber)
                        accumulatedCountFilename = '%s_%s_sim%d_AccumulatedCountArray.txt' %(type,sample,simulationNumber)
                        averageNucleosomeSignalFilename = '%s_%s_sim%d_AverageSignalArray.txt' %(type,sample,simulationNumber)

                os.makedirs(os.path.join(outputDir, jobname, DATA, SAMPLES, sample, occupancy_type,SIGNATUREBASED), exist_ok=True)

                accumulatedSignalFilePath = os.path.join(outputDir,jobname, DATA, SAMPLES, sample, occupancy_type, SIGNATUREBASED,accumulatedSignalFilename)
                accumulatedCountFilePath = os.path.join(outputDir,jobname, DATA, SAMPLES, sample, occupancy_type, SIGNATUREBASED,accumulatedCountFilename)
                averageNucleosomeSignalFilePath = os.path.join(outputDir,jobname, DATA, SAMPLES, sample,  occupancy_type, SIGNATUREBASED,averageNucleosomeSignalFilename)

                signalArray.tofile(file=accumulatedSignalFilePath, sep="\t", format="%s")
                countArray.tofile(file=accumulatedCountFilePath, sep="\t", format="%s")
                averageNucleosomeSignalArray.tofile(file=averageNucleosomeSignalFilePath, sep="\t", format="%s")
########################################################################################


########################################################################################
#All Mutations can be single point mutations or indels
def writeSampleBasedAverageNucleosomeOccupancyFiles(occupancy_type,
                                                    plusorMinus,
                                                    sample2Type2SignalArrayDict,
                                                    sample2Types2CountArrayDict,
                                                    outputDir,
                                                    jobname,
                                                    library_file_memo,
                                                    nucleosomeOccupancyAnalysisType,
                                                    simulationNumber):


    for sample in sample2Type2SignalArrayDict:
        type2SignalArray = sample2Type2SignalArrayDict[sample]
        type2CountArray = sample2Types2CountArrayDict[sample]

        if nucleosomeOccupancyAnalysisType in type2SignalArray:
            signalArray = type2SignalArray[nucleosomeOccupancyAnalysisType]
            countArray = type2CountArray[nucleosomeOccupancyAnalysisType]

            averageNucleosomeSignalArray = computeAverageNucleosomeOccupancyArray(plusorMinus,signalArray,countArray)

            if (simulationNumber==0):
                if library_file_memo is not None:
                    accumulatedSignalFilename = '%s_%s_%s_ccumulatedSignalArray.txt' % (sample, jobname,library_file_memo)
                    accumulatedCountFilename = '%s_%s_%s_AccumulatedCountArray.txt' % (sample, jobname,library_file_memo)
                    averageNucleosomeSignalFilename = '%s_%s_%s_AverageSignalArray.txt' % (sample, jobname,library_file_memo)
                else:
                    accumulatedSignalFilename = '%s_%s_AccumulatedSignalArray.txt' %(sample,jobname)
                    accumulatedCountFilename = '%s_%s_AccumulatedCountArray.txt' %(sample,jobname)
                    averageNucleosomeSignalFilename = '%s_%s_AverageSignalArray.txt' %(sample,jobname)
            else:
                if (library_file_memo is not None):
                    accumulatedSignalFilename = '%s_%s_sim%d_%s_AccumulatedSignalArray.txt' %(sample,jobname,simulationNumber,library_file_memo)
                    accumulatedCountFilename = '%s_%s_sim%d_%s_AccumulatedCountArray.txt' %(sample,jobname,simulationNumber,library_file_memo)
                    averageNucleosomeSignalFilename = '%s_%s_sim%d_%s_AverageSignalArray.txt' %(sample,jobname,simulationNumber,library_file_memo)
                else:
                    accumulatedSignalFilename = '%s_%s_sim%d_AccumulatedSignalArray.txt' %(sample,jobname,simulationNumber)
                    accumulatedCountFilename = '%s_%s_sim%d_AccumulatedCountArray.txt' %(sample,jobname,simulationNumber)
                    averageNucleosomeSignalFilename = '%s_%s_sim%d_AverageSignalArray.txt' %(sample,jobname,simulationNumber)

            os.makedirs(os.path.join(outputDir, jobname, DATA, SAMPLES, sample, occupancy_type, nucleosomeOccupancyAnalysisType), exist_ok=True)

            accumulatedSignalFilePath = os.path.join(outputDir, jobname, DATA, SAMPLES, sample, occupancy_type, nucleosomeOccupancyAnalysisType, accumulatedSignalFilename)
            accumulatedCountFilePath = os.path.join(outputDir, jobname, DATA, SAMPLES, sample, occupancy_type, nucleosomeOccupancyAnalysisType, accumulatedCountFilename)
            averageNucleosomeSignalFilePath = os.path.join(outputDir,jobname, DATA, SAMPLES, sample, occupancy_type, nucleosomeOccupancyAnalysisType, averageNucleosomeSignalFilename)

            signalArray.tofile(file=accumulatedSignalFilePath, sep="\t", format="%s")
            countArray.tofile(file=accumulatedCountFilePath, sep="\t", format="%s")
            averageNucleosomeSignalArray.tofile(file=averageNucleosomeSignalFilePath, sep="\t", format="%s")
########################################################################################


##########################################################################################
############### Common functions for Nucleosome Occupancy Analysis ends ##################
##########################################################################################

############################################################
#Used by DataPreparationCommons.py
#Kept here for further possible usage
# Notice that [::-1] provides visiting x from the last base to the first base
revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N'}[B] for B in x][::-1])
############################################################

############################################################
#Tested works correctly.
def getNucleotides(chromosomeShort,start,end,humanGenome):
    if chromosomeShort == 'MT':
        chromosomeShort = 'M'
    elif chromosomeShort == '23':
        chromosomeShort = 'X'
    elif chromosomeShort == '24':
        chromosomeShort = 'Y'

    chromosomeLong = 'chr' + chromosomeShort
    chrBased_humanGenome = humanGenome[chromosomeLong]
    seq = chrBased_humanGenome.get_slice(start,end)
    seq = seq.upper()
    return seq
############################################################


######################################################################
#BED files: end is not inclusive
#BED files: score - A score between 0 and 1000. See track lines, below, for ways to configure the display style of scored data.
def updateChrBasedSignalArray(data_row,chrBasedSignalArray):
    chrBasedSignalArray[data_row[START]:data_row[END]] += data_row[SIGNAL]
######################################################################


######################################################################
def fillNumpyArray(start,end,signal,chrBasedSignalArray):
    chrBasedSignalArray[start:end]+=signal
######################################################################


######################################################################
def decideFileType(library_file_with_path):
    with open(library_file_with_path, "r") as f:
        for line in f:
            if 'bedGraph' in line:
                return True
            elif ('fixedStep' in line) or ('variableStep' in line):
                return False
            else:
                return True
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
            file_df=pd.read_csv(file_with_path,header=None, usecols=[0, 1, 2, 3],names = [CHROM,START,END,SIGNAL],dtype={0: 'category', 1: np.int32, 2: np.int32, 3: np.float32},sep='\t')

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
            mean_signal=file_df[SIGNAL].mean()
            std_signal=file_df[SIGNAL].std()
            print('Max Signal: %f' %max_signal)
            print('Min Signal: %f' %min_signal)
            print('Mean Signal: %f' % mean_signal)
            print('Std Signal: %f' % std_signal)
            return file_df, max_signal, min_signal
        else:
            return file_df
        print('############################################')

    return file_df
##################################################################

##################################################################
#JAN 7, 2020
def generateIntervalVersion(wig_unprocessed_df):
    #Read the wig file and generate chr start end signal
    #default initialization
    step = 1
    span=1

    rows_list = []

    # e.g. rows
    # fixedStep chrom=chr1 start=24500 step=1000 span=1000
    # 57.4679
    # 57.467
    # 57.4651
    # 57.4623
    # 57.4586
    # 57.454
    # 57.4484
    # 57.442
    # 57.4347
    # 57.4266
    # 57.4176

    # variableStep chrom=chr2
    # 300701 12.5
    # 300702 12.5
    # 300703 12.5
    # 300704 12.5
    # 300705 12.5

    for row in wig_unprocessed_df.itertuples(index=True, name='Pandas'):
        # row's type is <class 'pandas.core.frame.Pandas'>
        # row[0] is the index
        # row[1] is fixedStep chrom=chr1 start=24500 step=1000 span=1000 or 57.4679
        if (row[1].startswith('fixedStep')):
            #This is the information line
            # chrom = row[1].split()[1].split('=')[1]
            # start = int(row[1].split()[2].split('=')[1])
            # step = int(row[1].split()[3].split('=')[1])
            step_type=FIXED_STEP
            for i in range(1,len(row[1].split())):
                element=row[1].split()[i]
                key=element.split('=')[0]
                value=element.split('=')[1]
                if key=='chrom':
                    chrom = value
                elif key=='start':
                    start=int(value)
                elif key=='step':
                    step=int(value)
                elif key=='span':
                    span=int(value)

        #Please notice that we do  not expect step in variableStep
        elif (row[1].startswith('variableStep')):
            step_type=VARIABLE_STEP
            for i in range(1,len(row[1].split())):
                element=row[1].split()[i]
                key=element.split('=')[0]
                value=element.split('=')[1]
                if key=='chrom':
                    chrom = value
                elif key=='start':
                    start=int(value)
                elif key=='step':
                    step=int(value)
                elif key=='span':
                    span=int(value)
        elif not (row[1].startswith('track') or row[1].startswith('browser')):
            if step_type==FIXED_STEP:
                #We read only signal
                signal = float(row[1])
                ##################
                end = start + span-1
                list = [chrom, start, end, signal]
                rows_list.append(list)
                start += step
            elif step_type==VARIABLE_STEP:
                #We read start and signal
                #49304701 10.0
                start = int(row[1].split()[0])
                signal = float(row[1].split()[1])
                ##################
                end = start + span-1
                list = [chrom, start, end, signal]
                rows_list.append(list)

    # print('Number of intervals to be inserted in wavelet_processed_df: %d' %len(rows_list))
    #rows_list contain the list of row where each row is a dictionary

    wig_chrom_start_end_signal_version_df = pd.DataFrame(rows_list, columns=[CHROM,START,END,SIGNAL])
    # print('wig_chrom_start_end_signal_version_df.dtypes')
    # print(wig_chrom_start_end_signal_version_df.dtypes)

    wig_chrom_start_end_signal_version_df[CHROM] = wig_chrom_start_end_signal_version_df[CHROM].astype(str)
    wig_chrom_start_end_signal_version_df[START] = wig_chrom_start_end_signal_version_df[START].astype(np.int32)
    wig_chrom_start_end_signal_version_df[END] = wig_chrom_start_end_signal_version_df[END].astype(np.int32)
    wig_chrom_start_end_signal_version_df[SIGNAL] = wig_chrom_start_end_signal_version_df[SIGNAL].astype(np.float32)

    # print('wig_chrom_start_end_signal_version_df.dtypes')
    # print(wig_chrom_start_end_signal_version_df.dtypes)

    return wig_chrom_start_end_signal_version_df
##################################################################


##################################################################
#JAN 7, 2020
#If file is originally a  wig file use this function
#No Outlier Elimination is done
def readWig_with_fixedStep_variableStep(wig_file_path):

    #Read the wavelet signal
    wig_unprocessed_df = pd.read_csv(wig_file_path, sep='\t', comment='#', header=None)

    #Process the wavelet signal, convert into interval version
    #Add column names
    wigfile_interval_version_df = generateIntervalVersion(wig_unprocessed_df)

    return wigfile_interval_version_df
##################################################################


######################################################################
def updateSignalArraysForListComprehension(row,signalArrayDict):
    #row [CHROM START END SIGNAL]
    signalArray=signalArrayDict[row[0]]
    signalArray[row[1]:row[2]] += row[3]
######################################################################


##################################################################
def updateDictionary(type2Strand2CountDict,mutationType,strand):
    if (mutationType in type2Strand2CountDict):
        if strand in type2Strand2CountDict[mutationType]:
            type2Strand2CountDict[mutationType][strand] += 1
        else:
            type2Strand2CountDict[mutationType][strand] = 1
    else:
        type2Strand2CountDict[mutationType] = {}
        type2Strand2CountDict[mutationType][strand] = 1
##################################################################

##################################################################
def updateSampleBasedDictionary(sample2Type2Strand2CountDict,signature_or_mutationType,mutationSample,strand):
    if (mutationSample in sample2Type2Strand2CountDict):
        if signature_or_mutationType in sample2Type2Strand2CountDict[mutationSample]:
            if strand in sample2Type2Strand2CountDict[mutationSample][signature_or_mutationType]:
                sample2Type2Strand2CountDict[mutationSample][signature_or_mutationType][strand] += 1
            else:
                sample2Type2Strand2CountDict[mutationSample][signature_or_mutationType][strand] = 1
        else:
            sample2Type2Strand2CountDict[mutationSample][signature_or_mutationType] = {}
            sample2Type2Strand2CountDict[mutationSample][signature_or_mutationType][strand] = 1
    else:
        sample2Type2Strand2CountDict[mutationSample] = {}
        sample2Type2Strand2CountDict[mutationSample][signature_or_mutationType] = {}
        sample2Type2Strand2CountDict[mutationSample][signature_or_mutationType][strand] = 1
##################################################################

##################################################################
#type can be signature or mutationType
def updateTypeBasedDictionary(type2Sample2Strand2CountDict,type,mutationSample,strand):
    if type in type2Sample2Strand2CountDict:
        if mutationSample in type2Sample2Strand2CountDict[type]:
            if strand in type2Sample2Strand2CountDict[type][mutationSample]:
                type2Sample2Strand2CountDict[type][mutationSample][strand] += 1
            else:
                type2Sample2Strand2CountDict[type][mutationSample][strand] = 1
        else:
            type2Sample2Strand2CountDict[type][mutationSample]= {}
            type2Sample2Strand2CountDict[type][mutationSample][strand] = 1
    else:
        type2Sample2Strand2CountDict[type]={}
        type2Sample2Strand2CountDict[type][mutationSample]={}
        type2Sample2Strand2CountDict[type][mutationSample][strand] =1
##################################################################

########################################################################
def updateDictionaries_simulations_integrated(mutation_row,
                        mutationType,
                        mutationSample,
                        sample_based,
                        simNum2Type2Strand2CountDict,
                        simNum2Sample2Type2Strand2CountDict,
                        simNum2Type2Sample2Strand2CountDict,
                        simNum2Signature2MutationType2Strand2CountDict,
                        strand,
                        signature_cutoff_numberofmutations_averageprobability_df):

    simNum = mutation_row[SIMULATION_NUMBER]
    type2Strand2CountDict = simNum2Type2Strand2CountDict[simNum]
    sample2Type2Strand2CountDict = simNum2Sample2Type2Strand2CountDict[simNum]
    type2Sample2Strand2CountDict = simNum2Type2Sample2Strand2CountDict[simNum]
    signature2MutationType2Strand2CountDict = simNum2Signature2MutationType2Strand2CountDict[simNum]

    #################################################################################################
    # Update1: update mutationType in type2Strand2CountDict
    if (mutationType is not None):
        updateDictionary(type2Strand2CountDict,mutationType,strand)
    #################################################################################################

    #################################################################################################
    # Update2: signature in type2Strand2CountDict
    for signature in signature_cutoff_numberofmutations_averageprobability_df['signature'].unique():
        #[cutoff numberofMutations averageProbability]
        cutoff=float(signature_cutoff_numberofmutations_averageprobability_df[signature_cutoff_numberofmutations_averageprobability_df['signature']==signature]['cutoff'].values[0])
        if (mutation_row[signature] >= cutoff):
            if signature in type2Strand2CountDict:
                if strand in type2Strand2CountDict[signature]:
                    type2Strand2CountDict[signature][strand] += 1
                else:
                    type2Strand2CountDict[signature][strand] = 1
            else:
                type2Strand2CountDict[signature] = {}
                type2Strand2CountDict[signature][strand] = 1
    #################################################################################################


    #################################################################################################
    #Update3 signature2MutationType2Strand2CountDict
    if (mutationType is not None):
        for signature in signature_cutoff_numberofmutations_averageprobability_df['signature'].unique():
            cutoff = float(signature_cutoff_numberofmutations_averageprobability_df[signature_cutoff_numberofmutations_averageprobability_df['signature'] == signature]['cutoff'].values[0])
            if (mutation_row[signature] >= cutoff):
                if signature in signature2MutationType2Strand2CountDict:
                    if mutationType in signature2MutationType2Strand2CountDict[signature]:
                        if strand in signature2MutationType2Strand2CountDict[signature][mutationType]:
                            signature2MutationType2Strand2CountDict[signature][mutationType][strand] +=1
                        else:
                            signature2MutationType2Strand2CountDict[signature][mutationType][strand] = 1

                    else:
                        signature2MutationType2Strand2CountDict[signature][mutationType] = {}
                        signature2MutationType2Strand2CountDict[signature][mutationType][strand] = 1
                else:
                    signature2MutationType2Strand2CountDict[signature] = {}
                    signature2MutationType2Strand2CountDict[signature][mutationType]={}
                    signature2MutationType2Strand2CountDict[signature][mutationType][strand] = 1
    #################################################################################################

    #################################################################################################
    # Update4: sample and mutationType in sample2Type2Strand2CountDict
    if sample_based:
        if (mutationType is not None):
            updateSampleBasedDictionary(sample2Type2Strand2CountDict,mutationType,mutationSample,strand)
    #################################################################################################

    #################################################################################################
    # Update5: sample and signature in sample2Type2Strand2CountDict
    if sample_based:
        for signature in signature_cutoff_numberofmutations_averageprobability_df['signature'].unique():
            #[cutoff numberofMutations averageProbability]
            cutoff = float(signature_cutoff_numberofmutations_averageprobability_df[signature_cutoff_numberofmutations_averageprobability_df['signature'] == signature]['cutoff'].values[0])
            if (mutation_row[signature] >= cutoff):
                updateSampleBasedDictionary(sample2Type2Strand2CountDict, signature, mutationSample, strand)
    #################################################################################################

    #################################################################################################
    #Update6 sample and mutationType type2Sample2TranscriptionStrand2CountDict
    if sample_based:
        if (type2Sample2Strand2CountDict is not None) and (mutationType is not None):
            updateTypeBasedDictionary(type2Sample2Strand2CountDict,mutationType,mutationSample,strand)
    #################################################################################################

    #################################################################################################
    #Update7: sample and signature in type2Sample2TranscriptionStrand2CountDict
    if sample_based:
        for signature in signature_cutoff_numberofmutations_averageprobability_df['signature'].unique():
            #[cutoff numberofMutations averageProbability]
            cutoff = float(signature_cutoff_numberofmutations_averageprobability_df[signature_cutoff_numberofmutations_averageprobability_df['signature'] == signature]['cutoff'].values[0])
            if (mutation_row[signature] >= cutoff):
                if (type2Sample2Strand2CountDict is not None):
                    updateTypeBasedDictionary(type2Sample2Strand2CountDict,signature,mutationSample,strand)
    #################################################################################################

########################################################################


########################################################################
def accumulate_simulations_integrated_for_each_tuple(
    chrBased_SimNum2Type2Strand2CountDict,
    chrBased_SimNum2Sample2Type2Strand2CountDict,
    chrBased_SimNum2Type2Sample2Strand2CountDict,
    chrBased_SimNum2Signature2MutationType2Strand2CountDict,
    accumulatedAllChromosomes_SimNum2Type2Strand2CountDict,
    accumulatedAllChromosomes_SimNum2Sample2Type2Strand2CountDict,
    accumulatedAllChromosomes_SimNum2Type2Sample2Strand2CountDict,
    accumulatedAllChromosomes_SimNum2Signature2MutationType2Strand2CountDict):

    # Accumulate1 chrBasedType2Strand2CountDict in accumulatedAllChromosomesType2Strand2CountDict
    if (chrBased_SimNum2Type2Strand2CountDict is not None):
        for simNum, chrBasedType2Strand2CountDict in chrBased_SimNum2Type2Strand2CountDict.items():
            if simNum in accumulatedAllChromosomes_SimNum2Type2Strand2CountDict:
                for type, strand2CountDict in chrBasedType2Strand2CountDict.items():
                    if type in accumulatedAllChromosomes_SimNum2Type2Strand2CountDict[simNum]:
                        for strand, count in strand2CountDict.items():
                            if strand in accumulatedAllChromosomes_SimNum2Type2Strand2CountDict[simNum][type]:
                                accumulatedAllChromosomes_SimNum2Type2Strand2CountDict[simNum][type][strand] += count
                            else:
                                accumulatedAllChromosomes_SimNum2Type2Strand2CountDict[simNum][type][strand] = count
                    else:
                        accumulatedAllChromosomes_SimNum2Type2Strand2CountDict[simNum][type] = strand2CountDict
            else:
                accumulatedAllChromosomes_SimNum2Type2Strand2CountDict[simNum] = chrBasedType2Strand2CountDict

    # Accumulate2 chrBasedSample2Type2ReplicationStrand2CountDict in accumulatedAllChromosomesSample2Type2ReplicationStrand2CountDict
    if (chrBased_SimNum2Sample2Type2Strand2CountDict is not None):
        for simNum, chrBasedSample2Type2Strand2CountDict in chrBased_SimNum2Sample2Type2Strand2CountDict.items():
            if simNum in accumulatedAllChromosomes_SimNum2Sample2Type2Strand2CountDict:
                for sample, type2Strand2CountDict in chrBasedSample2Type2Strand2CountDict.items():
                    if sample in accumulatedAllChromosomes_SimNum2Sample2Type2Strand2CountDict[simNum]:
                        for type, strand2CountDict in type2Strand2CountDict.items():
                            if type in accumulatedAllChromosomes_SimNum2Sample2Type2Strand2CountDict[simNum][sample]:
                                for strand, count in strand2CountDict.items():
                                    if strand in accumulatedAllChromosomes_SimNum2Sample2Type2Strand2CountDict[simNum][
                                        sample][type]:
                                        accumulatedAllChromosomes_SimNum2Sample2Type2Strand2CountDict[simNum][sample][
                                            type][strand] += count
                                    else:
                                        accumulatedAllChromosomes_SimNum2Sample2Type2Strand2CountDict[simNum][sample][
                                            type][strand] = count
                            else:
                                accumulatedAllChromosomes_SimNum2Sample2Type2Strand2CountDict[simNum][sample][
                                    type] = strand2CountDict
                    else:
                        accumulatedAllChromosomes_SimNum2Sample2Type2Strand2CountDict[simNum][
                            sample] = type2Strand2CountDict
            else:
                accumulatedAllChromosomes_SimNum2Sample2Type2Strand2CountDict[
                    simNum] = chrBasedSample2Type2Strand2CountDict

    # Accumulate3 chrBasedType2Sample2Strand2CountDict in accumulatedAllChromosomesType2Sample2Strand2CountDict
    if (chrBased_SimNum2Type2Sample2Strand2CountDict is not None):
        for simNum, chrBasedType2Sample2Strand2CountDict in chrBased_SimNum2Type2Sample2Strand2CountDict.items():
            if simNum in accumulatedAllChromosomes_SimNum2Type2Sample2Strand2CountDict:
                for type, sample2Strand2CountDict in chrBasedType2Sample2Strand2CountDict.items():
                    if type in accumulatedAllChromosomes_SimNum2Type2Sample2Strand2CountDict[simNum]:
                        for sample, strand2CountDict in sample2Strand2CountDict.items():
                            if sample in accumulatedAllChromosomes_SimNum2Type2Sample2Strand2CountDict[simNum][type]:
                                for strand, count in strand2CountDict.items():
                                    if strand in \
                                            accumulatedAllChromosomes_SimNum2Type2Sample2Strand2CountDict[simNum][type][
                                                sample]:
                                        accumulatedAllChromosomes_SimNum2Type2Sample2Strand2CountDict[simNum][type][
                                            sample][strand] += count
                                    else:
                                        accumulatedAllChromosomes_SimNum2Type2Sample2Strand2CountDict[simNum][type][
                                            sample][strand] = count
                            else:
                                accumulatedAllChromosomes_SimNum2Type2Sample2Strand2CountDict[simNum][type][
                                    sample] = strand2CountDict
                    else:
                        accumulatedAllChromosomes_SimNum2Type2Sample2Strand2CountDict[simNum][
                            type] = sample2Strand2CountDict
            else:
                accumulatedAllChromosomes_SimNum2Type2Sample2Strand2CountDict[
                    simNum] = chrBasedType2Sample2Strand2CountDict

    # Accumulate4 chrBasedSignature2MutationType2Strand2CountDict in accumulatedAllChromosomesSignature2MutationType2Strand2CountDict
    if (chrBased_SimNum2Signature2MutationType2Strand2CountDict is not None):
        for simNum, chrBasedSignature2MutationType2Strand2CountDict in chrBased_SimNum2Signature2MutationType2Strand2CountDict.items():
            if simNum in accumulatedAllChromosomes_SimNum2Signature2MutationType2Strand2CountDict:
                for signature, mutationType2Strand2CountDict in chrBasedSignature2MutationType2Strand2CountDict.items():
                    if signature in accumulatedAllChromosomes_SimNum2Signature2MutationType2Strand2CountDict[simNum]:
                        for mutationType, strand2CountDict in mutationType2Strand2CountDict.items():
                            if mutationType in \
                                    accumulatedAllChromosomes_SimNum2Signature2MutationType2Strand2CountDict[simNum][
                                        signature]:
                                for strand, count in strand2CountDict.items():
                                    if strand in \
                                            accumulatedAllChromosomes_SimNum2Signature2MutationType2Strand2CountDict[
                                                simNum][signature][mutationType]:
                                        accumulatedAllChromosomes_SimNum2Signature2MutationType2Strand2CountDict[
                                            simNum][signature][mutationType][strand] += count
                                    else:
                                        accumulatedAllChromosomes_SimNum2Signature2MutationType2Strand2CountDict[
                                            simNum][signature][mutationType][strand] = count
                            else:
                                accumulatedAllChromosomes_SimNum2Signature2MutationType2Strand2CountDict[simNum][
                                    signature][mutationType] = strand2CountDict
                    else:
                        accumulatedAllChromosomes_SimNum2Signature2MutationType2Strand2CountDict[simNum][
                            signature] = mutationType2Strand2CountDict
            else:
                accumulatedAllChromosomes_SimNum2Signature2MutationType2Strand2CountDict[
                    simNum] = chrBasedSignature2MutationType2Strand2CountDict
########################################################################


########################################################################
def accumulate_simulations_integrated(listofTuples,
                                      accumulatedAllChromosomes_SimNum2Type2Strand2CountDict,
                                      accumulatedAllChromosomes_SimNum2Sample2Type2Strand2CountDict,
                                      accumulatedAllChromosomes_SimNum2Type2Sample2Strand2CountDict,
                                      accumulatedAllChromosomes_SimNum2Signature2MutationType2Strand2CountDict):

    for mytuple in listofTuples:
        chrBased_SimNum2Type2Strand2CountDict = mytuple[0]
        chrBased_SimNum2Sample2Type2Strand2CountDict = mytuple[1]
        chrBased_SimNum2Type2Sample2Strand2CountDict = mytuple[2]
        chrBased_SimNum2Signature2MutationType2Strand2CountDict = mytuple[3]

        accumulate_simulations_integrated_for_each_tuple(
            chrBased_SimNum2Type2Strand2CountDict,
            chrBased_SimNum2Sample2Type2Strand2CountDict,
            chrBased_SimNum2Type2Sample2Strand2CountDict,
            chrBased_SimNum2Signature2MutationType2Strand2CountDict,
            accumulatedAllChromosomes_SimNum2Type2Strand2CountDict,
            accumulatedAllChromosomes_SimNum2Sample2Type2Strand2CountDict,
            accumulatedAllChromosomes_SimNum2Type2Sample2Strand2CountDict,
            accumulatedAllChromosomes_SimNum2Signature2MutationType2Strand2CountDict)
########################################################################

########################################################################
#To write samples with signatures with at least 10K eligible mutations
def writeDictionaryUnderDataDirectory(dictionary,outputDir,jobname,filename):

    os.makedirs(os.path.join(outputDir,jobname,DATA), exist_ok=True)
    filePath = os.path.join(outputDir,jobname,DATA,filename)

    with open(filePath, 'w') as file:
        file.write(json.dumps(dictionary))
########################################################################


########################################################################
class MyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(MyEncoder, self).default(obj)
########################################################################

########################################################################
def writeDictionaryUsingPickle(dictionary,filepath):
    fileDir = os.path.dirname(filepath)
    os.makedirs(fileDir, exist_ok=True)

    with open(filepath, "wb") as file:
        pickle.dump(dictionary, file)
########################################################################

########################################################################
def writeDictionarySimple(dictionary,path,filename,customJSONEncoder):
    os.makedirs(os.path.join(path), exist_ok=True)

    filePath = os.path.join(path,filename)
    with open(filePath, 'w') as file:
        file.write(json.dumps(dictionary, cls=customJSONEncoder))
########################################################################


########################################################################
def writeDictionary(dictionary,outputDir,jobname,filename,subDirectory,customJSONEncoder):
    os.makedirs(os.path.join(outputDir,jobname,DATA,subDirectory), exist_ok=True)

    filePath = os.path.join(outputDir,jobname,DATA,subDirectory,filename)
    with open(filePath, 'w') as file:
        file.write(json.dumps(dictionary, cls=customJSONEncoder))
########################################################################

########################################################################
def readDictionaryUsingPickle(filePath):
    with open(filePath, "rb") as file:
        dictionary  = pickle.load(file)
    return  dictionary
########################################################################

########################################################################
def readDictionary(filePath):
    if (os.path.exists(filePath) and (os.path.getsize(filePath) > 0)):
        with open(filePath,'r') as json_data:
            dictionary = json.load(json_data)
        return dictionary
    else:
        # return None
        # Provide empty dictionary for not to fail for loops on None type dictionary
        return {}
########################################################################


#########################################################################
# copy file 1.maf from inputDir/output/simulations/jobname_simulations_genome_96/1.maf to inputDir/output/simulations/sim1/96/1.maf
# copy file 2.maf from inputDir/output/simulations/jobname_simulations_genome_96/2.maf to inputDir/output/simulations/sim2/96/2.maf
# copy file N.maf from inputDir/output/simulations/jobname_simulations_genome_96/N.maf to inputDir/output/simulations/simN/96/N.maf

# copy file 1.maf from inputDir/output/simulations/jobname_simulations_genome_ID/1.maf to inputDir/output/simulations/sim1/ID/1.maf
# copy file 2.maf from inputDir/output/simulations/jobname_simulations_genome_ID/2.maf to inputDir/output/simulations/sim2/ID/2.maf
# copy file N.maf from inputDir/output/simulations/jobname_simulations_genome_ID/N.maf to inputDir/output/simulations/simN/ID/N.maf

# copy file 1.maf from inputDir/output/simulations/jobname_simulations_genome_DBS/1.maf to inputDir/output/simulations/sim1/DBS/1.maf
# copy file 2.maf from inputDir/output/simulations/jobname_simulations_genome_DBS/2.maf to inputDir/output/simulations/sim2/DBS/2.maf
# copy file N.maf from inputDir/output/simulations/jobname_simulations_genome_DBS/N.maf to inputDir/output/simulations/simN/DBS/N.maf
def copyMafFiles(copyFromDir,copyToMainDir,mutation_type_context,numberofSimulations):
    # copyFromDir = os.path.join(inputDir, 'output', 'simulations', dirName)
    # copyToMainDir = os.path.join(inputDir, 'output', 'simulations')

    # No need to  go to this directory
    # os.chdir(copyFromDir)

    for simNum in range(1,numberofSimulations+1):
        simDir= 'sim%d' %(simNum)
        fname='%d.maf' %(simNum)
        fileDir=os.path.join(copyFromDir,fname)
        copyToDir = os.path.join(copyToMainDir, simDir, mutation_type_context)
        shutil.copy(fileDir, copyToDir)
#########################################################################


########################################################################
#for control purposes
def checkForSumInMutationType2Strand2CountDict(mutationType2Strand2CountDict,mutationType2Sample2Strand2CountDict):

    controlMutationType2Strand2CountDict = {}

    #Accumulate controlMutationType2Strand2CountDict
    for mutationType, sample2Strand2CountDict in mutationType2Sample2Strand2CountDict.items():
        if mutationType not in controlMutationType2Strand2CountDict:
            controlMutationType2Strand2CountDict[mutationType] = {}
        for sample, strand2CountDict in sample2Strand2CountDict.items():
            for strand, count in strand2CountDict.items():
                if strand in controlMutationType2Strand2CountDict[mutationType]:
                    controlMutationType2Strand2CountDict[mutationType][strand] += count
                else:
                    controlMutationType2Strand2CountDict[mutationType][strand] = count


    #Now check whether they havethe same counts or not
    for mutationType, strand2CountDict in mutationType2Strand2CountDict.items():
        for strand, count in strand2CountDict.items():
            if controlMutationType2Strand2CountDict[mutationType][strand] != count:
                return False;

    return True
########################################################################


########################################################################
#for control purposes
def checkForSumInMutationProbability2Signature2Strand2CountDict(mutationProbability2Signature2Strand2CountDict,mutationProbability2Signature2Sample2Strand2CountDict):

    controlMutationProbability2Signature2Strand2CountDict = {}

    #Accumulate controlMutationProbability2Signature2Strand2CountDict
    for mutationProbability, signature2Sample2Strand2CountDict in mutationProbability2Signature2Sample2Strand2CountDict.items():
        if mutationProbability not in controlMutationProbability2Signature2Strand2CountDict:
            controlMutationProbability2Signature2Strand2CountDict[mutationProbability] = {}
        for signature, sample2Strand2CountDict in signature2Sample2Strand2CountDict.items():
            if signature not in controlMutationProbability2Signature2Strand2CountDict[mutationProbability]:
                controlMutationProbability2Signature2Strand2CountDict[mutationProbability][signature] = {}
            for sample, strand2CountDict in sample2Strand2CountDict.items():
                for strand, count in strand2CountDict.items():
                     if strand not in controlMutationProbability2Signature2Strand2CountDict[mutationProbability][signature]:
                        controlMutationProbability2Signature2Strand2CountDict[mutationProbability][signature][strand] = count
                     else:
                        controlMutationProbability2Signature2Strand2CountDict[mutationProbability][signature][strand] += count

    #Now check whether they have the same counts or not
    for mutationProbability, signature2Strand2CountDict in mutationProbability2Signature2Strand2CountDict.items():
        for signature, strand2CountDict in signature2Strand2CountDict.items():
            for strand, count in strand2CountDict.items():
                if controlMutationProbability2Signature2Strand2CountDict[mutationProbability][signature][strand] != count:
                    return False;

    return True
########################################################################

########################################################################
# A helper function
# Maybe can be used later
#accumulatedAllChromosomesMutationProbability2Signature2StrandCountDict is full
#accumulatedAllChromosomesMutationProbability2Signature2RatioDict is empty, filled in this function
def calculateRatio(accumulatedAllChromosomesMutationProbability2Signature2StrandCountDict,accumulatedAllChromosomesMutationProbability2Signature2RatioDict,strandList):
    numeratorStrand = strandList[0]
    denominatorStrand = strandList[1]

    for mutationProbability, signature2StrandCountDict in accumulatedAllChromosomesMutationProbability2Signature2StrandCountDict.items():
        accumulatedAllChromosomesMutationProbability2Signature2RatioDict[mutationProbability] = {}
        for signature, strandCountDict in signature2StrandCountDict.items():
            if ((numeratorStrand in strandCountDict) and (denominatorStrand in strandCountDict)):
                if((strandCountDict[numeratorStrand] + strandCountDict[denominatorStrand])>=SUBS_STRAND_BIAS_NUMBER_OF_MUTATIONS_THRESHOLD):
                    accumulatedAllChromosomesMutationProbability2Signature2RatioDict[mutationProbability][signature] = (strandCountDict[numeratorStrand]) / (strandCountDict[numeratorStrand] + strandCountDict[denominatorStrand])
########################################################################


########################################################################
# A helper function
# Maybe can be used later
#mutationProbability2Signature2ReplicationStrand2CountDict contains all the accumulated counts for all chromosomes
#Ratio is laggingCount/(laggingCount + leadingCount)
#Ratio is the transcribed/ (transcribed + non_transcribed)
# strandNameList = [LAGGING,LEADING] for replicationStrandAnalysis
# strandNameList = [TRANSCRIBED_STRAND, NON_TRANSCRIBED_STRAND]    for transcribedStrandAnalysis
def convert(mutationProbability2Signature2Strand2CountDict, strandNameList):
    signature2SumofMutationProbabilitiesDict = {}
    signature2MutationProbability2RatioDict = {}
    signature2WeightedAverageRatioDict = {}
    signature2StdErrorDict = {}

    numeratorStrand = strandNameList[0]
    denominatorStrand = strandNameList[1]

    #Fill signature2WeightedAverageRatioDict
    for mutationProbability, signature2Strand2CountDict in mutationProbability2Signature2Strand2CountDict.items():
        for signature, strand2CountDict in signature2Strand2CountDict.items():
            #Only once
            if signature not in signature2MutationProbability2RatioDict:
                signature2MutationProbability2RatioDict[signature] = {}

            if ((numeratorStrand in strand2CountDict) and (denominatorStrand in strand2CountDict)):
                #In order to consider there must be at least ONE_THOUSAND mutations on the leading and lagging strands
                # In order to consider there must be at least ONE_THOUSAND mutations on the transcribed and non-transcribed strands
                if (strand2CountDict[numeratorStrand]+strand2CountDict[denominatorStrand] >= SUBS_STRAND_BIAS_NUMBER_OF_MUTATIONS_THRESHOLD):
                    if signature in signature2SumofMutationProbabilitiesDict:
                        signature2SumofMutationProbabilitiesDict[signature] += mutationProbability
                    else:
                        signature2SumofMutationProbabilitiesDict[signature] = mutationProbability

                    ratio = (strand2CountDict[numeratorStrand])/(strand2CountDict[numeratorStrand] + strand2CountDict[denominatorStrand])
                    signature2MutationProbability2RatioDict[signature][mutationProbability]= ratio
                    #First Accumulate
                    if signature in signature2WeightedAverageRatioDict:
                        signature2WeightedAverageRatioDict[signature] += mutationProbability * ratio
                    else:
                        signature2WeightedAverageRatioDict[signature] = mutationProbability * ratio

    #Then divide by sumofMutationProbabilities
    for signature in signature2WeightedAverageRatioDict.keys():
        if (signature2SumofMutationProbabilitiesDict[signature]!=0):
            signature2WeightedAverageRatioDict[signature] /= signature2SumofMutationProbabilitiesDict[signature]
        else:
            # print('debug starts')
            # print('For signature: %s signature2SumofMutationProbabilitiesDict[%s] is zero' %(signature,signature))
            # print('debug ends')
            signature2WeightedAverageRatioDict[signature] = 0

    #Calculate the signature2StdErrorDict
    for signature, weightedAverageRatio in signature2WeightedAverageRatioDict.items():
        variance= 0
        sampleSize = 0
        if signature in signature2MutationProbability2RatioDict:
            mutationProbability2RatioDict = signature2MutationProbability2RatioDict[signature]
            for mutationProbability,ratio in mutationProbability2RatioDict.items():
                difference = weightedAverageRatio-ratio
                variance += difference**2
                sampleSize +=1

            stddev = math.sqrt(variance/sampleSize)
            stderr = stddev/math.sqrt(sampleSize)
            signature2StdErrorDict[signature] = stderr

    return signature2WeightedAverageRatioDict,signature2StdErrorDict, signature2SumofMutationProbabilitiesDict
########################################################################


#################################################################
#JAN 20, 2020
# From DataPreparationCommons.py


############################################################
#Same for Release and old_PCAWG Matlab Probabilities
# example line for original data
#UCEC-US_SP89389 10      2017540 N:AT[T>A]CA     1
# example line for simulated data
#UCEC-US_SP89389_1       10      1575080 T:AT[C>T]TG     1
def readChrBasedMutations(chr_based_mutation_filepath,mutation_type_context):

    if (os.path.exists(chr_based_mutation_filepath)):
        try:
            mutations_with_genomic_positions_df = pd.read_csv(chr_based_mutation_filepath, sep='\t', header=None)
        except pd.errors.EmptyDataError:
            mutations_with_genomic_positions_df = pd.DataFrame()

        if (not mutations_with_genomic_positions_df.empty):
            if (mutation_type_context==DBS):
                # For DBS MatrixGenerator provides
                # UAD-US_SP50263 10      110099884       Q:T[GC>AG]C     0
                # For DBS Extractor has
                # Sample Names    MutationTypes   DBS2    DBS4    DBS5    DBS6    DBS9    DBS11
                # LUAD-US_SP50518 AC>CA   0.004819278307958045    0.09604751880153346     0.0     0.07540775450119858     0.8237254483893098      0.0
                mutations_with_genomic_positions_df.columns = [SAMPLE,CHROM,START,MUTATIONLONG,PYRAMIDINESTRAND]
                mutations_with_genomic_positions_df[SAMPLE] = mutations_with_genomic_positions_df[SAMPLE].astype(str)
                mutations_with_genomic_positions_df[CHROM] = mutations_with_genomic_positions_df[CHROM].astype(str)
                mutations_with_genomic_positions_df[START] = mutations_with_genomic_positions_df[START].astype(int)
                mutations_with_genomic_positions_df[MUTATIONLONG] = mutations_with_genomic_positions_df[MUTATIONLONG].astype(str)
                mutations_with_genomic_positions_df[PYRAMIDINESTRAND] = mutations_with_genomic_positions_df[PYRAMIDINESTRAND].astype(int)
                #Add new columns
                # MatrixGenerator generates Q:A[AC>TT]A
                # PCAWG_Matlab dbs probabilities has  AT>GC
                mutations_with_genomic_positions_df[TRANSCRIPTIONSTRAND] = mutations_with_genomic_positions_df[MUTATIONLONG].str[0]
                mutations_with_genomic_positions_df[MUTATION] = mutations_with_genomic_positions_df[MUTATIONLONG].str[4:9]
            elif (mutation_type_context==ID):
                # For ID MatrixGenerator provides
                # LUAD-US_SP50263 10      8045169 U:2:Ins:R:5     T       TTC     1
                # For ID Extractor has
                # Sample Names    MutationTypes   ID1     ID2     ID3     ID4     ID5     ID6     ID8     ID9     ID13
                # LUAD-US_SP50518 1:Del:C:0       0.002114485363152394    0.0     0.4891560241408412      0.01586903032961574     0.2531175852230711      0.0     0.23974287494331953     0.0     0.0
                mutations_with_genomic_positions_df.columns = [SAMPLE,CHROM,START,MUTATIONLONG,REF,ALT,PYRAMIDINESTRAND]
                mutations_with_genomic_positions_df[SAMPLE] = mutations_with_genomic_positions_df[SAMPLE].astype(str)
                mutations_with_genomic_positions_df[CHROM] = mutations_with_genomic_positions_df[CHROM].astype(str)
                mutations_with_genomic_positions_df[START] = mutations_with_genomic_positions_df[START].astype(int)
                mutations_with_genomic_positions_df[MUTATIONLONG] = mutations_with_genomic_positions_df[MUTATIONLONG].astype(str)
                mutations_with_genomic_positions_df[PYRAMIDINESTRAND] = mutations_with_genomic_positions_df[PYRAMIDINESTRAND].astype(int)
                #Add new column
                # MatrixGenerator generates N:1:Ins:T:5
                # PCAWG_Matlab id probabilities has 1:Ins:T:1
                mutations_with_genomic_positions_df[LENGTH] = mutations_with_genomic_positions_df[REF].apply(len)
                mutations_with_genomic_positions_df[TRANSCRIPTIONSTRAND] = mutations_with_genomic_positions_df[MUTATIONLONG].str[0]
                mutations_with_genomic_positions_df[MUTATION] = mutations_with_genomic_positions_df[MUTATIONLONG].str[2:]
                #order the columns make CONTEXT at the end.
                ordered_column_names = [SAMPLE,CHROM,START,MUTATIONLONG,REF,ALT,LENGTH,PYRAMIDINESTRAND,TRANSCRIPTIONSTRAND,MUTATION]
                mutations_with_genomic_positions_df = mutations_with_genomic_positions_df[ordered_column_names]
            elif(mutation_type_context in SBS_CONTEXTS):
                # For SNV MatrixGenerator provides
                # LUAD-US_SP50263 10      440625  U:GG[C>T]AG     -1
                # For SNV Extractor has
                # Sample    Mutation   SBS1    SBS2    SBS3    SBS4    SBS5    SBS13   SBS17a  SBS17b  SBS18   SBS28   SBS40
                # LUAD-US_SP50518 A[C>A]A 0.005281537126491598    1.091660854697097e-06   0.0     0.0     0.12212513236310162     0.0022549935716281717   0.0     0.0     0.0     0.0     0.8703372452779239
                mutations_with_genomic_positions_df.columns = [SAMPLE,CHROM,START,MUTATIONLONG,PYRAMIDINESTRAND]
                mutations_with_genomic_positions_df[SAMPLE] = mutations_with_genomic_positions_df[SAMPLE].astype(str)
                mutations_with_genomic_positions_df[CHROM] = mutations_with_genomic_positions_df[CHROM].astype(str)
                mutations_with_genomic_positions_df[START] = mutations_with_genomic_positions_df[START].astype(int)
                mutations_with_genomic_positions_df[MUTATIONLONG] = mutations_with_genomic_positions_df[MUTATIONLONG].astype(str)
                mutations_with_genomic_positions_df[PYRAMIDINESTRAND] = mutations_with_genomic_positions_df[PYRAMIDINESTRAND].astype(int)
                #Add new column
                # Add Context Column from T:TG[C>T]GC to  G[C>T]G
                # MatrixGenerator generates T:TG[C>T]GC
                # Extractor and PCAWG_Matlab sbs probabilities has G[C>T]G
                mutations_with_genomic_positions_df[TRANSCRIPTIONSTRAND] = mutations_with_genomic_positions_df[MUTATIONLONG].str[0]
                mutations_with_genomic_positions_df[MUTATION] = mutations_with_genomic_positions_df[MUTATIONLONG].str[3:10]
            return mutations_with_genomic_positions_df

    return None
############################################################

############################################################
def readChrBasedMutationsMergeWithProbabilitiesAndWrite(inputList):
    chrShort = inputList[0]
    outputDir = inputList[1]
    jobname = inputList[2]
    chr_based_mutation_filepath = inputList[3]
    mutations_probabilities_df = inputList[4]
    mutation_type_context = inputList[5]
    simNum = inputList[6]
    PCAWG=inputList[7]

    ###############################################################################################
    chr_based_mutation_df = readChrBasedMutations(chr_based_mutation_filepath,mutation_type_context)
    ###############################################################################################

    if ((chr_based_mutation_df is not None) and (mutations_probabilities_df is not None)):

        ############################################################################
        #Step2 SigProfilerTopography Python Package

        #For PCAWG_Matlab
        #Convert CMDI-UK_SP116871_1 --> SP116871 # if(simNum>0): simNum=1 Simulation1
        #Convert CMDI-UK_SP116871 --> SP116871 # simNum=0 Original Data
        if PCAWG:
            chr_based_mutation_df[SAMPLE] = chr_based_mutation_df[SAMPLE].str.split('_', expand=True)[1]

        #For Release SigProfilerTopography Python Package
        # For SNV
        # LUAD-US_SP50263 10      440625  U:GG[C>T]AG     -1
        # For ID
        # LUAD-US_SP50263 10      8045169 U:2:Ins:R:5     T       TTC     1
        # For DBS
        # UAD-US_SP50263 10      110099884       Q:T[GC>AG]C     0
        elif simNum>=1:
            # Get rid of simulation number at the end
            chr_based_mutation_df[SAMPLE] = chr_based_mutation_df[SAMPLE].str.rsplit('_', 1, expand=True)[0]
        ############################################################################

        ############################################################################
        merged_df = pd.merge(chr_based_mutation_df,mutations_probabilities_df, how='inner', left_on=[SAMPLE, MUTATION],right_on=[SAMPLE, MUTATION])
        ############################################################################

        if ((merged_df is not None) and (chr_based_mutation_df.shape[0]!=merged_df.shape[0])):
            print('##############################')
            print('There is a situation/problem. For simNum:%s chr:%s All %s mutation context type mutations are not merged with signature probabilities'  %(simNum,chrShort,mutation_type_context))
            print('For simNum:%s chr:%s mutation context type:%s chr_based_mutation_df.shape(%d,%d)-- merged_df.shape(%d,%d) ' % (simNum, chrShort, mutation_type_context,chr_based_mutation_df.shape[0],chr_based_mutation_df.shape[1],merged_df.shape[0],merged_df.shape[1]))

        if ((merged_df is not None) and (not merged_df.empty)):
            if (mutation_type_context in SBS_CONTEXTS):
                chrBasedMergedMutationsFileName = 'chr%s_%s_for_topography.txt' %(chrShort,SUBS)
            elif (mutation_type_context == DBS):
                chrBasedMergedMutationsFileName = 'chr%s_%s_for_topography.txt' %(chrShort,DINUCS)
            elif (mutation_type_context==ID):
                chrBasedMergedMutationsFileName = 'chr%s_%s_for_topography.txt' %(chrShort,INDELS)

            if (simNum ==0):
                chr_based_merged_mutations_file_path = os.path.join(outputDir, jobname, DATA, CHRBASED,chrBasedMergedMutationsFileName)
            else:
                simDir = 'sim%d' %(simNum)
                chr_based_merged_mutations_file_path = os.path.join(outputDir, jobname, DATA, CHRBASED,simDir,chrBasedMergedMutationsFileName)

            # #After test uncomment
            # if ('MutationLong' in merged_df.columns.values):
            #     merged_df.drop(['MutationLong'], inplace=True, axis=1)

            #After merge
            if (mutation_type_context in SBS_CONTEXTS):
                merged_df[MUTATION] = merged_df[MUTATION].str[2:5]

            merged_df.to_csv(chr_based_merged_mutations_file_path, sep='\t', header=True, index=False)
        else:
            print('-------------No merge file for sim%d mutation_type_context:%s for chr%s' %(simNum,mutation_type_context,chrShort))
############################################################


############################################################
def readProbabilities(probabilitiesFile):

    #For Release and PCAWG_Matlab
    #This is same for Release and PCAWG_Matlab
    # There is header in the first column
    #Sample names can be composed of numbers

    probabilities_df = pd.read_csv(probabilitiesFile, sep='\t', header=0,nrows=1)

    #No need
    # probabilities_df[SAMPLE]=probabilities_df[SAMPLE].astype('category')
    # probabilities_df[MUTATION] = probabilities_df[MUTATION].astype('category')

    if ('Sample Names' in probabilities_df.columns.values) and ('MutationTypes' in probabilities_df.columns.values):
        probabilities_df = pd.read_csv(probabilitiesFile, sep='\t', header=0,dtype={'Sample Names': 'category', 'MutationTypes': 'category'})
        probabilities_df.rename(columns={'Sample Names': SAMPLE, 'MutationTypes': MUTATION}, inplace=True)
    else:
        probabilities_df = pd.read_csv(probabilitiesFile, sep='\t', header=0,dtype={SAMPLE: 'category', MUTATION: 'category'})

    # Mutation_Probabilities.txt for SBS96
    # Sample Names    MutationTypes   SBS1    SBS2    SBS3    SBS4    SBS5    SBS13   SBS17a  SBS17b  SBS18   SBS28   SBS40
    # LUAD-US_SP50518 A[C>A]A 0.005281537126491598    1.091660854697097e-06   0.0     0.0     0.12212513236310162     0.0022549935716281717   0.0     0.0     0.0     0.0     0.8703372452779239

    # Mutation_Probabilities.txt for ID83
    # Sample Names    MutationTypes   ID1     ID2     ID3     ID4     ID5     ID6     ID8     ID9     ID13
    # LUAD-US_SP50518 1:Del:C:0       0.002114485363152394    0.0     0.4891560241408412      0.01586903032961574     0.2531175852230711      0.0     0.23974287494331953     0.0     0.0

    # Mutation_Probabilities.txt for DBS78
    # Sample Names    MutationTypes   DBS2    DBS4    DBS5    DBS6    DBS9    DBS11
    # LUAD-US_SP50518 AC>CA   0.004819278307958045    0.09604751880153346     0.0     0.07540775450119858     0.8237254483893098      0.0

    print('Probabilities information starts')
    print('probabilities_df.shape')
    print(probabilities_df.shape)
    print('probabilities_df.head()')
    print(probabilities_df.head())
    print('probabilities_df.dtypes')
    print(probabilities_df.dtypes)

    if (SAMPLE in probabilities_df.columns.values):
        print('Unique samples in probabilities_df')
        print(probabilities_df[SAMPLE].unique())
        print('# of unique samples in probabilities_df: %d\n' %(len(probabilities_df[SAMPLE].unique())))

    if (MUTATION in probabilities_df.columns.values):
        print('Unique MutationTypes in probabilities_df')
        print(probabilities_df[MUTATION].unique())
        print('# of unique mutation types in probabilities_df: %d' %(len(probabilities_df[MUTATION].unique())))
    print('Probabilities information ends')
    print('##############################')

    return probabilities_df
############################################################