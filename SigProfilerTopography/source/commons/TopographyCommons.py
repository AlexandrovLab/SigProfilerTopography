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

import scipy.stats as stats

#To handle warnings as errors
# import warnings
# warnings.filterwarnings("error")

LINUX='linux'
UNIX='unix'
WINDOWS='windows'

current_abs_path = os.path.dirname(os.path.realpath(__file__))

LEADING= 'Leading'
LAGGING = 'Lagging'

GENIC='Genic'
INTERGENIC='Intergenic'

UNTRANSCRIBED_STRAND = 'UnTranscribed'
TRANSCRIBED_STRAND = 'Transcribed'
NONTRANSCRIBED_STRAND = 'NonTranscribed'

LAGGING_VERSUS_LEADING='Lagging_Versus_Leading'
TRANSCRIBED_VERSUS_UNTRANSCRIBED='Transcribed_Versus_Untranscribed'
GENIC_VERSUS_INTERGENIC='Genic_Versus_Intergenic'

LAGGING_VERSUS_LEADING_P_VALUE='lagging_versus_leading_p_value'
TRANSCRIBED_VERSUS_UNTRANSCRIBED_P_VALUE='transcribed_versus_untranscribed_p_value'
GENIC_VERSUS_INTERGENIC_P_VALUE='genic_versus_intergenic_p_value'

LAGGING_VERSUS_LEADING_Q_VALUE='lagging_versus_leading_q_value'
TRANSCRIBED_VERSUS_UNTRANSCRIBED_Q_VALUE='transcribed_versus_untranscribed_q_value'
GENIC_VERSUS_INTERGENIC_Q_VALUE='genic_versus_intergenic_q_value'

TRANSCRIBED_REAL_COUNT='Transcribed_real_count'
UNTRANSCRIBED_REAL_COUNT='UnTranscribed_real_count'

GENIC_REAL_COUNT='genic_real_count'
INTERGENIC_REAL_COUNT='intergenic_real_count'

LAGGING_REAL_COUNT ='Lagging_real_count'
LEADING_REAL_COUNT ='Leading_real_count'

TRANSCRIBED_SIMULATIONS_MEAN_COUNT='Transcribed_mean_sims_count'
UNTRANSCRIBED_SIMULATIONS_MEAN_COUNT='UnTranscribed_mean_sims_count'

GENIC_SIMULATIONS_MEAN_COUNT ='genic_mean_sims_count'
INTERGENIC_SIMULATIONS_MEAN_COUNT ='intergenic_mean_sims_count'

LAGGING_SIMULATIONS_MEAN_COUNT ='Lagging_mean_sims_count'
LEADING_SIMULATIONS_MEAN_COUNT ='Leading_mean_sims_count'

AT_LEAST_5_PERCENT_DIFF='5%'
AT_LEAST_10_PERCENT_DIFF='10%'
AT_LEAST_20_PERCENT_DIFF='20%'
AT_LEAST_30_PERCENT_DIFF='30%'
AT_LEAST_50_PERCENT_DIFF='50%'
AT_LEAST_75_PERCENT_DIFF='75%'
AT_LEAST_100_PERCENT_DIFF='100%'

percentage_numbers = [10, 20, 30, 50, 75, 100]
percentage_strings = [AT_LEAST_10_PERCENT_DIFF, AT_LEAST_20_PERCENT_DIFF, AT_LEAST_30_PERCENT_DIFF, AT_LEAST_50_PERCENT_DIFF, AT_LEAST_75_PERCENT_DIFF, AT_LEAST_100_PERCENT_DIFF]

PLUS = '+'
MINUS = '-'

MAXIMUM_CHROMOSOME_LENGTH = 250000000
# 1024*1024*1024 = 1073741824
# 1024*1024 = 1048576
GIGABYTE_IN_BYTES = 1073741824
MEGABYTE_IN_BYTES = 1048576

NUMBER_OF_MUTATIONS_IN_EACH_SPLIT=100000
MAXIMUM_NUMBER_JOBS_IN_THE_POOL_AT_ONCE=50

BIGWIG='BIGWIG'
BIGBED='BIGBED'
WIG='WIG'
BED='BED'
BEDGRAPH='BEDGRAPH'
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

MEF="MEF"

#Epigenomics File for GRCh37
DEFAULT_H3K27ME3_OCCUPANCY_FILE = 'ENCFF291WFP_breast_epithelium_H3K27me3.bed'
DEFAULT_H3K36ME3_OCCUPANCY_FILE = 'ENCFF906MJM_breast_epithelium_H3K36me3.bed'
DEFAULT_H3K9ME3_OCCUPANCY_FILE = 'ENCFF065FJK_breast_epithelium_H3K9me3.bed'
DEFAULT_H3K27AC_OCCUPANCY_FILE = 'ENCFF154XFN_breast_epithelium_H3K27ac.bed'
DEFAULT_H3K4ME1_OCCUPANCY_FILE = 'ENCFF336DDM_breast_epithelium_H3K4me1.bed'
DEFAULT_H3K4ME3_OCCUPANCY_FILE = 'ENCFF065TIH_breast_epithelium_H3K4me3.bed'
DEFAULT_CTCF_OCCUPANCY_FILE = 'ENCFF782GCQ_breast_epithelium_Normal_CTCF-human.bed'
DEFAULT_ATAC_SEQ_OCCUPANCY_FILE = 'ENCFF035ICJ_breast_epithelium_Normal_ATAC-seq.wig'

#Epigenomics File for mm10
ENCFF575PMI_mm10_embryonic_facial_prominence_ATAC_seq="ENCFF575PMI_mm10_embryonic_facial_prominence_ATAC_seq.wig"
ENCFF993SRY_mm10_embryonic_fibroblast_H3K4me1="ENCFF993SRY_mm10_embryonic_fibroblast_H3K4me1.bed"
ENCFF912DNP_mm10_embryonic_fibroblast_H3K4me3="ENCFF912DNP_mm10_embryonic_fibroblast_H3K4me3.bed"
ENCFF611HDQ_mm10_embryonic_fibroblast_CTCF="ENCFF611HDQ_mm10_embryonic_fibroblast_CTCF.bed"
ENCFF152DUV_mm10_embryonic_fibroblast_POLR2A="ENCFF152DUV_mm10_embryonic_fibroblast_POLR2A.bed"
ENCFF114VLZ_mm10_embryonic_fibroblast_H3K27ac="ENCFF114VLZ_mm10_embryonic_fibroblast_H3K27ac.bed"

#NUCLEOSOME OCCUPANCY FILES
MM10_MEF_NUCLEOSOME_FILE='Mus_musculus.mmNuc0030101.nucleosome.shift.bw'
GM12878_NUCLEOSOME_OCCUPANCY_FILE = 'wgEncodeSydhNsomeGm12878Sig.bigWig'
K562_NUCLEOSOME_OCCUPANCY_FILE = 'wgEncodeSydhNsomeK562Sig.bigWig'

SIGPROFILERTOPOGRAPHY_DEFAULT_FILES=[GM12878_NUCLEOSOME_OCCUPANCY_FILE,K562_NUCLEOSOME_OCCUPANCY_FILE,DEFAULT_ATAC_SEQ_OCCUPANCY_FILE]

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

MEF_REPLICATION_TIME_SIGNAL_FILE = 'ENCFF001JVQ_mm10_embryonic_fibroblast_wavelet_smoothed_signal.wig'


available_nucleosome_biosamples=[GM12878,K562,MEF]
available_replication_time_biosamples=[GM12878,K562,MCF7,HEPG2,HELAS3,SKNSH,IMR90,NHEK,BJ,HUVEC,BG02ES,GM06990,GM12801,GM12812,GM12813]

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
COMBINE_P_VALUES_METHOD_FISHER='fisher'
COMBINE_P_VALUES_METHOD_STOUFFER='stouffer'
WEIGHTED_AVERAGE_METHOD='WEIGHTED_AVERAGE_METHOD'
COLORBAR_SEISMIC='seismic'
COLORBAR_DISCREET='discreet'
NUCLEOSOME_BIOSAMPLE="K562"
PLOTS='plots'
TABLES='tables'
DETAILED='detailed'
EXCEL_FILES='excel_files'
NUCLEOSOME_DNA_ELEMENT='Nucleosome'
ATAC_DNA_ELEMENT='ATAC'
OPEN_CHROMATIN='Open\nChromatin'
###################################################################################################

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
Table_MutationType_NumberofMutations_NumberofSamples_SamplesList_Filename='Table_MutationType_NumberofMutations_NumberofSamples_SamplesList.txt'
Table_ChrLong_NumberofMutations_Filename='Table_ChrLong_NumberofMutations.txt'

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

#These are used for getting chromosome names
GRCh37ChromSizesDictFilename = 'hg19ChromSizesDict.txt'
GRCh38ChromSizesDictFilename = 'hg38ChromSizesDict.txt'
MM9ChromSizesDictFilename = 'mm9ChromSizesDict.txt'
MM10ChromSizesDictFilename = 'mm10ChromSizesDict.txt'

INDELBASED = 'indelbased'
SIGNATUREBASED = 'signaturebased'
SAMPLEBASED = 'samplebased'

SAMPLE_BASED = 'sample_based'

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
SCATTER_PLOTS= 'scatter_plots'
BAR_PLOTS= 'bar_plots'
CIRCLE_PLOTS= 'circle_plots'
CIRCLE_BAR_PLOTS = 'circle_bar_plots'
HEATMAPS='heatmaps'
PLOTTING = 'plotting'
TRANSCRIPTIONSTRANDBIAS = 'transcription_strand_bias'
REPLICATIONSTRANDBIAS = 'replication_strand_bias'

SAMPLES = 'samples'

TRANSCRIPTION_LABEL = 'Transcription'
REPLICATION_LABEL = 'Replication'

NUMBER_OF_MUTATIONS = 'number_of_mutations'
MUTATION_DENSITY = 'mutation_density'
NORMALIZED_MUTATION_DENSITY = 'normalized_mutation_density'

MICROHOMOLOGY = 'Microhomology_mediated_indels'
REPEAT = 'Repeat_mediated_indels'

DEFICIENT = 'deficient'
PROFICIENT = 'proficient'

#Mutation Types
C2A = 'C>A'
C2G = 'C>G'
C2T = 'C>T'
T2A = 'T>A'
T2C = 'T>C'
T2G = 'T>G'

six_mutation_types = [C2A,C2G,C2T,T2A,T2C,T2G]

SIGNATURE = 'Signature'
MUTATION = 'Mutation'
MUTATIONLONG = 'MutationLong'

COMPUTATION_CHROMOSOMES_SEQUENTIAL_ALL_SIMULATIONS_PARALLEL = 'COMPUTATION_CHROMOSOMES_SEQUENTIAL_ALL_SIMULATIONS_PARALLEL'
USING_IMAP_UNORDERED='USING_IMAP_UNORDERED'
USING_APPLY_ASYNC='USING_APPLY_ASYNC'
USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM='USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM'
USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT='USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT'
USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT_USING_POOL_INPUT_LIST='USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT_USING_POOL_INPUT_LIST'

CONSIDER_DISTANCE='CONSIDER_DISTANCE' #default
CONSIDER_COUNT='CONSIDER_COUNT'
CONSIDER_DISTANCE_ALL_SAMPLES_TOGETHER='CONSIDER_DISTANCE_ALL_SAMPLES_TOGETHER'

MISSING_SIGNAL='MISSING_SIGNAL'  #default
NO_SIGNAL='NO_SIGNAL'

RELAXED='relaxed'
STRINGENT='stringent'

DEFAULT_AVERAGE_PROBABILITY=0.75

DEFAULT_NUM_OF_SBS_REQUIRED=2000
DEFAULT_NUM_OF_DBS_REQUIRED=200
DEFAULT_NUM_OF_ID_REQUIRED=1000

DEFAULT_NUM_OF_REAL_DATA_OVERLAP_REQUIRED=100
NUMBER_OF_REQUIRED_MUTATIONS_FOR_STRAND_BIAS_BAR_PLOT=1

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
MUTATION_LONG = 'MutationLong'
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
SBS288 = '288'
SBS384  = '384'
SBS1536 = '1536'
SBS3072 = '3072'

#These 3 (SNV,ID,DBS) are used for matrixgenerator creared directories
SNV='SNV'
ID = 'ID'
DBS= 'DBS'

SBS='SBS'
SBS_CONTEXTS = [SBS96,SBS192,SBS288,SBS384,SBS1536,SBS3072]

# Used for dictionaries
# MutationType2NumberofMutationsDict keys
SUBS = 'SUBS'
INDELS = 'INDELS'
DINUCS = 'DINUCS'



UNDECLARED = 'Undeclared'

AVERAGE_SIGNAL_ARRAY = 'AverageSignalArray'
ACCUMULATED_COUNT_ARRAY = 'AccumulatedCountArray'
ACCUMULATED_SIGNAL_ARRAY = 'AccumulatedSignalArray'

PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_TOOL='PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_TOOL'
PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_MANUSCRIPT='PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_MANUSCRIPT'
PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_MANUSCRIPT_OCCUPANCY_ANALYSIS_FIGURE='PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_MANUSCRIPT_OCCUPANCY_ANALYSIS_FIGURE'

############################################################################################################################
#sheet name must be less than 31 characters
def write_excel_file(df_list, sheet_list, file_name):
    writer = pd.ExcelWriter(file_name,engine='xlsxwriter')
    for dataframe, sheet in zip(df_list, sheet_list):
        dataframe.to_excel(writer, sheet_name=sheet, startrow=0 , startcol=0, index=False)
    writer.save()
############################################################################################################################



# ############################################################
def get_mutation_type_context_for_probabilities_file(mutation_types_contexts_for_signature_probabilities,mutation_type):
    if (mutation_type==SUBS) and (SBS96 in mutation_types_contexts_for_signature_probabilities):
        return SBS96
    elif (mutation_type==SUBS) and (SBS192 in mutation_types_contexts_for_signature_probabilities):
        return SBS192
    elif (mutation_type==SUBS) and (SBS288 in mutation_types_contexts_for_signature_probabilities):
        return SBS288
    elif  (mutation_type==SUBS) and (SBS384 in mutation_types_contexts_for_signature_probabilities):
        return SBS384
    elif  (mutation_type==SUBS) and (SBS1536 in mutation_types_contexts_for_signature_probabilities):
        return SBS1536
    elif  (mutation_type==SUBS) and (SBS3072 in mutation_types_contexts_for_signature_probabilities):
        return SBS3072
    elif (mutation_type==INDELS) and (ID in mutation_types_contexts_for_signature_probabilities):
        return ID
    elif (mutation_type==DINUCS) and (DBS in mutation_types_contexts_for_signature_probabilities):
        return DBS
    else:
        return None
# ############################################################

##################################################################
#Dec 17, 2020
def get_chrBased_simBased_dfs(outputDir,jobname,chrLong,simNum):
    chrBased_simBased_subs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, SUBS, simNum)
    chrBased_simBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, DINUCS, simNum)
    chrBased_simBased_indels_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, INDELS, simNum)

    return chrBased_simBased_subs_df, chrBased_simBased_dinucs_df, chrBased_simBased_indels_df
##################################################################

##################################################################
# April 27, 2020
def get_chrBased_simBased_combined_df(outputDir,jobname,chrLong,simNum):

    chrBased_simBased_subs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, SUBS, simNum)
    chrBased_simBased_indels_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, INDELS, simNum)
    chrBased_simBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, DINUCS, simNum)

    if (chrBased_simBased_subs_df is not None):
        chrBased_simBased_subs_df[TYPE] = SUBS

    if (chrBased_simBased_indels_df is not None):
        chrBased_simBased_indels_df[TYPE] = INDELS

    if (chrBased_simBased_dinucs_df is not None):
        chrBased_simBased_dinucs_df[TYPE] = DINUCS

    if (chrBased_simBased_subs_df is not None) or (chrBased_simBased_indels_df is not None) or (chrBased_simBased_dinucs_df is not None):
        chrBased_simBased_combined_df = pd.concat([chrBased_simBased_subs_df, chrBased_simBased_indels_df, chrBased_simBased_dinucs_df], ignore_index=True,axis=0)
        return chrBased_simBased_combined_df
    else:
        return None
##################################################################


##################################################################
#April 26, 2020
def index_marks(nrows, chunk_size):
    return range(chunk_size, math.ceil(nrows / chunk_size) * chunk_size, chunk_size)
##################################################################


##################################################################
#April 26, 2020
#Returns split df list
def get_chrBased_simBased_combined_chunks_df(outputDir,jobname,chrLong,simNum):

    chrBased_simBased_subs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, SUBS, simNum)
    chrBased_simBased_indels_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, INDELS, simNum)
    chrBased_simBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, DINUCS, simNum)

    if (chrBased_simBased_subs_df is not None):
        chrBased_simBased_subs_df[TYPE] = SUBS

    if (chrBased_simBased_indels_df is not None):
        chrBased_simBased_indels_df[TYPE] = INDELS

    if (chrBased_simBased_dinucs_df is not None):
        chrBased_simBased_dinucs_df[TYPE] = DINUCS

    if (chrBased_simBased_subs_df is not None) or (chrBased_simBased_indels_df is not None) or (chrBased_simBased_dinucs_df is not None):
        chrBased_simBased_combined_df = pd.concat([chrBased_simBased_subs_df, chrBased_simBased_indels_df, chrBased_simBased_dinucs_df], ignore_index=True,axis=0)

        chrBased_simBased_number_of_mutations = chrBased_simBased_combined_df.shape[0]
        indices=index_marks(chrBased_simBased_number_of_mutations,NUMBER_OF_MUTATIONS_IN_EACH_SPLIT)
        # print(type(indices),flush=True)
        # print(indices,flush=True)
        #np.split returns a list of dataframes
        split_df_list=np.split(chrBased_simBased_combined_df, indices)
        return split_df_list
    else:
        return None
##################################################################

##################################################################
#April 9, 2020
#Returns split df
def get_chrBased_simBased_combined_df_split(outputDir,jobname,chrLong,simNum,splitIndex):
    chrBased_simBased_subs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, SUBS, simNum)
    chrBased_simBased_indels_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, INDELS, simNum)
    chrBased_simBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, DINUCS, simNum)

    if (chrBased_simBased_subs_df is not None):
        chrBased_simBased_subs_df[TYPE] = SUBS

    if (chrBased_simBased_indels_df is not None):
        chrBased_simBased_indels_df[TYPE] = INDELS

    if (chrBased_simBased_dinucs_df is not None):
        chrBased_simBased_dinucs_df[TYPE] = DINUCS

    if (chrBased_simBased_subs_df is not None) or (chrBased_simBased_indels_df is not None) or (chrBased_simBased_dinucs_df is not None):
        chrBased_simBased_combined_df = pd.concat([chrBased_simBased_subs_df, chrBased_simBased_indels_df, chrBased_simBased_dinucs_df], ignore_index=True,axis=0)

        ###########################################################
        chrBased_simBased_number_of_mutations = chrBased_simBased_combined_df.shape[0]
        number_of_splits = math.ceil(chrBased_simBased_number_of_mutations / NUMBER_OF_MUTATIONS_IN_EACH_SPLIT)

        #################################
        split_start_end_tuples = []
        start = 0
        for split in range(1, number_of_splits + 1):
            end = start + NUMBER_OF_MUTATIONS_IN_EACH_SPLIT
            if end > chrBased_simBased_combined_df.shape[0]:
                end = chrBased_simBased_combined_df.shape[0]
            split_start_end_tuples.append((start, end))
            start = end
        #################################

        #################################
        if (splitIndex<len(split_start_end_tuples)):
            split_start, split_end = split_start_end_tuples[splitIndex]
            # print('MONITOR %s simNum:%d splitIndex:%d split_start:%d split_end:%d len(split_start_end_tuples):%d split_start_end_tuples:%s' % (chrLong, simNum, splitIndex, split_start, split_end, len(split_start_end_tuples),split_start_end_tuples), flush=True)
            chrBased_simBased_combined_df_split = chrBased_simBased_combined_df.iloc[split_start:split_end, :]
            return chrBased_simBased_combined_df_split
        else:
            # print('MONITOR %s simNum:%d splitIndex:%d Does not exists ' % (chrLong, simNum, splitIndex), flush=True)
            return None
        #################################
    else:
        return None
##################################################################


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
        if (nucleosome_biosample==MEF):
            nucleosome_file=MM10_MEF_NUCLEOSOME_FILE
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
        if (replication_time_biosample==MEF):
            replication_time_signal_file = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB, REPLICATION,MEF_REPLICATION_TIME_SIGNAL_FILE)
        elif (replication_time_biosample==GM12878):
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
            mu = np.mean(colwise_array)
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
#Bookkeeping
#http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.chrom.sizes
#http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
#http://hgdownload.cse.ucsc.edu/goldenpath/mm9/bigZips/mm9.chrom.sizes
#http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.chrom.sizes
def getChromSizesDict(genome):
    chromSizesDict = {}

    if (genome==GRCh37):
        chromSizesDictPath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB,UCSCGENOME,GRCh37ChromSizesDictFilename)
    elif (genome==GRCh38):
        chromSizesDictPath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB,UCSCGENOME,GRCh38ChromSizesDictFilename)
    elif (genome==MM9):
        chromSizesDictPath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB,UCSCGENOME,MM9ChromSizesDictFilename)
    elif (genome==MM10):
        chromSizesDictPath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB,UCSCGENOME,MM10ChromSizesDictFilename)

    if (os.path.exists(chromSizesDictPath)):
        chromSizesDict = readDictionary(chromSizesDictPath)

    return chromSizesDict
###################################################################


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
#two header lines
#rows will be cutoff
#columns will be for each signature number_of_mutations and average_probability
def writeAllCutoffs(outputDir, jobname, DATA,cutoff2Signature2NumberofMutationsAverageProbabilityListDict,table_allcutoffs_signature_numberofmutations_averageprobability_filename):
    signature_list=[]
    cutoff_list=sorted(cutoff2Signature2NumberofMutationsAverageProbabilityListDict.keys())

    for cutoff in cutoff2Signature2NumberofMutationsAverageProbabilityListDict:
        for signature in cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff]:
            if signature not in signature_list:
                signature_list.append(signature)

    signature_list=sorted(signature_list,key=natural_key)

    allcutoffs_signature_numberofmutations_averageprobability = open(os.path.join(outputDir,jobname,DATA, table_allcutoffs_signature_numberofmutations_averageprobability_filename), 'w')

    # 1st header line
    allcutoffs_signature_numberofmutations_averageprobability.write('\t\t' )
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
def fillCutoff2Signature2PropertiesListDictionary(outputDir,
                                                  jobname,
                                                  chromNamesList,
                                                  mutation_type,
                                                  cutoffs,
                                                  average_probability,
                                                  num_of_sbs_required,
                                                  num_of_id_required,
                                                  num_of_dbs_required,
                                                  mutationType2PropertiesDict,
                                                  chrLong2NumberofMutationsDict):

    #Filled in the first part
    #PropertiesList consists of [sum_of_number of mutations, sum of probabilities, samples_list]
    cutoff2Signature2PropertiesListDict={}

    #Filled in the second part
    # [number of mutations, average probability, samples_list]
    cutoff2Signature2NumberofMutationsAverageProbabilityListDict={}

    #Filled in the third part
    #PropertiesList=[ cutoff number_of_mutations average_mutation_probability samples_list]
    signature2PropertiesListDict={}

    #This samples are for this mutation type
    all_samples=set()

    for chrLong in chromNamesList:
        chrbased_samples, chrBased_mutation_df = readChrBasedMutationsDF(outputDir,jobname,chrLong,mutation_type,0,return_number_of_samples=True)

        all_samples = all_samples.union(chrbased_samples)

        if ((chrBased_mutation_df is not None) and (not chrBased_mutation_df.empty)):
            # Sample  Chrom   Start   PyrimidineStrand        Mutation        DBS2    DBS4    DBS6    DBS7    DBS11
            # PD10011a        10      24033661        1       TC>AA   0.0     0.7656325053758131      0.15420390829468886     0.07918943063517644     0.000974155694321615
            signatures = getSignatures(chrBased_mutation_df)

            if mutation_type in mutationType2PropertiesDict:
                mutationType2PropertiesDict[mutation_type]['number_of_mutations'] += chrBased_mutation_df.shape[0]
                mutationType2PropertiesDict[mutation_type]['number_of_samples'] = len(all_samples)
                mutationType2PropertiesDict[mutation_type]['samples_list'] = list(all_samples)
            else:
                mutationType2PropertiesDict[mutation_type] = {}
                mutationType2PropertiesDict[mutation_type]['number_of_mutations']=chrBased_mutation_df.shape[0]
                mutationType2PropertiesDict[mutation_type]['number_of_samples'] = len(all_samples)
                mutationType2PropertiesDict[mutation_type]['samples_list'] = list(all_samples)

            if chrLong in chrLong2NumberofMutationsDict:
                chrLong2NumberofMutationsDict[chrLong] += chrBased_mutation_df.shape[0]
            else:
                chrLong2NumberofMutationsDict[chrLong] = chrBased_mutation_df.shape[0]

            # First part starts
            # First accumulate number of mutations and sum of probabilities with mutations with probability >= cutoff probability
            for cutoff in cutoffs:
                for signature in signatures:
                    # chrBased_mutation_df[signature]=chrBased_mutation_df[signature].astype(np.float64)
                    number_of_mutations = len(chrBased_mutation_df[chrBased_mutation_df[signature]>=float(cutoff)])
                    # Samples having mutations with prob ge cutoff for given cutoff and signature
                    samples_array=chrBased_mutation_df[chrBased_mutation_df[signature]>=float(cutoff)]['Sample'].unique()

                    # This results in infinity
                    # sum_of_probabilities = (chrBased_mutation_df[chrBased_mutation_df[signature]>=float(cutoff)])[signature].sum()
                    sum_of_probabilities = np.sum(((chrBased_mutation_df[chrBased_mutation_df[signature]>=float(cutoff)])[signature]).values,dtype=np.float64)

                    if cutoff not in cutoff2Signature2PropertiesListDict:
                        cutoff2Signature2PropertiesListDict[cutoff]={}
                        cutoff2Signature2PropertiesListDict[cutoff][signature]=[]
                        cutoff2Signature2PropertiesListDict[cutoff][signature].append(np.int64(0))
                        cutoff2Signature2PropertiesListDict[cutoff][signature].append(np.float64(0.0))
                        cutoff2Signature2PropertiesListDict[cutoff][signature].append(samples_array.tolist())
                        cutoff2Signature2PropertiesListDict[cutoff][signature][0]+=number_of_mutations
                        cutoff2Signature2PropertiesListDict[cutoff][signature][1]+=sum_of_probabilities
                    elif signature not in cutoff2Signature2PropertiesListDict[cutoff]:
                        cutoff2Signature2PropertiesListDict[cutoff][signature]=[]
                        cutoff2Signature2PropertiesListDict[cutoff][signature].append(np.int64(0))
                        cutoff2Signature2PropertiesListDict[cutoff][signature].append(np.float64(0.0))
                        cutoff2Signature2PropertiesListDict[cutoff][signature].append(samples_array.tolist())
                        cutoff2Signature2PropertiesListDict[cutoff][signature][0]+=number_of_mutations
                        cutoff2Signature2PropertiesListDict[cutoff][signature][1]+=sum_of_probabilities
                    else:
                        cutoff2Signature2PropertiesListDict[cutoff][signature][0]+=number_of_mutations
                        cutoff2Signature2PropertiesListDict[cutoff][signature][1]+=sum_of_probabilities
                        existing_samples_list = cutoff2Signature2PropertiesListDict[cutoff][signature][2]
                        #Update existing_samples with new samples
                        cutoff2Signature2PropertiesListDict[cutoff][signature][2]=list(set(existing_samples_list).union(set(samples_array)))
        #First part ends
        #Accumulation ended for each chromosome


    #Second part starts
    #Second: Find average probability
    # Fill another dictionary: cutoff2Signature2NumberofMutationsAverageProbabilityListDict
    # In fact, in this part, we can fill a dataframe
    for cutoff in cutoff2Signature2PropertiesListDict:
        for signature in cutoff2Signature2PropertiesListDict[cutoff]:
            if (cutoff2Signature2PropertiesListDict[cutoff][signature][0]>0):
                if cutoff not in cutoff2Signature2NumberofMutationsAverageProbabilityListDict:
                    cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff]={}
                    cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature] = []
                    cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature].append(np.int64(cutoff2Signature2PropertiesListDict[cutoff][signature][0]))
                    cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature].append(np.float64(cutoff2Signature2PropertiesListDict[cutoff][signature][1]/cutoff2Signature2PropertiesListDict[cutoff][signature][0]))
                    cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature].append(cutoff2Signature2PropertiesListDict[cutoff][signature][2])
                elif signature not in cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff]:
                    cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature] = []
                    cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature].append(np.int64(cutoff2Signature2PropertiesListDict[cutoff][signature][0]))
                    cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature].append(np.float64(cutoff2Signature2PropertiesListDict[cutoff][signature][1]/cutoff2Signature2PropertiesListDict[cutoff][signature][0]))
                    cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature].append(cutoff2Signature2PropertiesListDict[cutoff][signature][2])
                else:
                    print('There is a situation in fillCutoff2Signature2PropertiesListDictionary method.')
                    cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature].append(np.int64(cutoff2Signature2PropertiesListDict[cutoff][signature][0]))
                    cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature].append(np.float64(cutoff2Signature2PropertiesListDict[cutoff][signature][1]/cutoff2Signature2PropertiesListDict[cutoff][signature][0]))
                    cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature].append(cutoff2Signature2PropertiesListDict[cutoff][signature][2])
    #Second part ends


    #Set the filenames and number of required mutations
    if (mutation_type==SUBS):
        number_of_required_mutations=num_of_sbs_required
        table_allcutoffs_signature_numberofmutations_averageprobability_filename = Table_AllCutoff_SubsSignature_NumberofMutations_AverageProbability_Filename
        table_signature_cutoff_numberofmutations_averageprobability_filename = Table_SubsSignature_Cutoff_NumberofMutations_AverageProbability_Filename
    elif (mutation_type == INDELS):
        number_of_required_mutations=num_of_id_required
        table_allcutoffs_signature_numberofmutations_averageprobability_filename = Table_AllCutoff_IndelsSignature_NumberofMutations_AverageProbability_Filename
        table_signature_cutoff_numberofmutations_averageprobability_filename = Table_IndelsSignature_Cutoff_NumberofMutations_AverageProbability_Filename
    elif (mutation_type== DINUCS):
        number_of_required_mutations=num_of_dbs_required
        table_allcutoffs_signature_numberofmutations_averageprobability_filename = Table_AllCutoff_DinucsSignature_NumberofMutations_AverageProbability_Filename
        table_signature_cutoff_numberofmutations_averageprobability_filename = Table_DinucsSignature_Cutoff_NumberofMutations_AverageProbability_Filename

    #Third find the signature based cufoff probability with number of mutations >= required number of mutations  and averega mutation probability >=required_probability
    sorted_cutoffs=sorted(cutoff2Signature2NumberofMutationsAverageProbabilityListDict.keys())
    for cutoff in sorted_cutoffs:
        for signature in cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff]:
            if signature not in signature2PropertiesListDict:
                if (cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature][0]>=number_of_required_mutations) and (cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature][1]>=average_probability):
                    signature2PropertiesListDict[signature]=[cutoff,
                                                             np.int(cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature][0]),
                                                             np.float(cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature][1]),
                                                             cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature][2]]
    #Third part ends


    #Write the dictionaries
    os.makedirs(os.path.join(outputDir,jobname,DATA), exist_ok=True)

    ####################################################################
    # Jobname is added here as the first column
    # AverageProbabilityList contains  number_of_mutations average_probability samples_list
    # Let's not write the samples_list to avoid repetions
    writeAllCutoffs(outputDir,jobname,DATA,cutoff2Signature2NumberofMutationsAverageProbabilityListDict,table_allcutoffs_signature_numberofmutations_averageprobability_filename)
    ####################################################################


    ####################################################################
    L = sorted([(jobname,
                 signature,
                 cutoff_numberofmutations_averageprobability[0],
                 cutoff_numberofmutations_averageprobability[1],
                 cutoff_numberofmutations_averageprobability[2],
                 cutoff_numberofmutations_averageprobability[3],
                 len(cutoff_numberofmutations_averageprobability[3]))
                 for signature, cutoff_numberofmutations_averageprobability in signature2PropertiesListDict.items()])

    if L:
        signature_cutoff_numberofmutations_averageprobability_df = pd.DataFrame(L,columns=['cancer_type',
                                                                                           'signature',
                                                                                           'cutoff',
                                                                                           'number_of_mutations',
                                                                                           'average_probability',
                                                                                           'samples_list',
                                                                                           'len(samples_list)'])

        signature_cutoff_numberofmutations_averageprobability_df['cancer_type']=signature_cutoff_numberofmutations_averageprobability_df['cancer_type'].astype(str)
        signature_cutoff_numberofmutations_averageprobability_df['signature']=signature_cutoff_numberofmutations_averageprobability_df['signature'].astype(str)
        signature_cutoff_numberofmutations_averageprobability_df['cutoff']=signature_cutoff_numberofmutations_averageprobability_df['cutoff'].astype(np.float32)
        signature_cutoff_numberofmutations_averageprobability_df['number_of_mutations']=signature_cutoff_numberofmutations_averageprobability_df['number_of_mutations'].astype(np.int32)
        signature_cutoff_numberofmutations_averageprobability_df['average_probability'] = signature_cutoff_numberofmutations_averageprobability_df['average_probability'].astype(np.float32)
    else:
        #Create empty dataframe
        signature_cutoff_numberofmutations_averageprobability_df = pd.DataFrame(columns=['cancer_type',
                                                                                           'signature',
                                                                                           'cutoff',
                                                                                           'number_of_mutations',
                                                                                           'average_probability',
                                                                                           'samples_list',
                                                                                           'len(samples_list)'])

    signature_cutoff_numberofmutations_averageprobability_df.to_csv(os.path.join(outputDir,jobname,DATA, table_signature_cutoff_numberofmutations_averageprobability_filename), sep='\t', index=False)

    signature_cutoff_numberofmutations_averageprobability_df.drop(['cancer_type', 'samples_list', 'len(samples_list)'], inplace=True, axis=1)

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
def readChrBasedMutationsDF(outputDir,jobname,chrLong,mutation_type,simulationNumber,return_number_of_samples=False):
    filename = '%s_%s_for_topography.txt' %(chrLong,mutation_type)

    if (simulationNumber==0):
        chrBasedMutationDFFilePath = os.path.join(outputDir,jobname,DATA,CHRBASED,filename)
    else:
        simulation = 'sim%s' % (simulationNumber)
        chrBasedMutationDFFilePath = os.path.join(outputDir,jobname,DATA,CHRBASED,simulation,filename)

    chrBased_mutation_df = None
    #############################################
    if (os.path.exists(chrBasedMutationDFFilePath)):
        #############################################
        try:
         only_header_chrBased_mutation_df = pd.read_csv(chrBasedMutationDFFilePath, sep='\t', comment='#', nrows=1)
         # Please note thar encoding and engine slow done and increase memory usage
         # only_header_chrBased_mutation_df = pd.read_csv(chrBasedMutationDFFilePath, sep='\t', comment='#', nrows=1,encoding='utf8',engine='python')
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

         # np.float16 Half precision float: sign bit, 5 bits exponent, 10 bits mantissa
         # np.int8 Byte (-128 to 127)
         # np.int16 Integer (-32768 to 32767)
         # np.int32   Integer (-2147483648 to 2147483647)

         # To lower the dataframe size
         for signature in signatures:
             mydtypes[signature] = np.float16

         if ((mutation_type == SUBS) or (mutation_type == DINUCS)):
             mydtypes[SAMPLE] = 'category'
             mydtypes[CHROM] = 'category'
             mydtypes[START] = np.int32
             mydtypes[MUTATIONLONG] = 'category'
             mydtypes[PYRAMIDINESTRAND] = np.int8
             mydtypes[TRANSCRIPTIONSTRAND] = 'category'
             mydtypes[MUTATION] = 'category'

         if (mutation_type == INDELS):
             mydtypes[SAMPLE] = 'category'
             mydtypes[CHROM] = 'category'
             mydtypes[START] = np.int32
             mydtypes[MUTATIONLONG] = 'category'
             mydtypes[REF] = str
             mydtypes[ALT] = str
             mydtypes[LENGTH] = np.int16
             mydtypes[PYRAMIDINESTRAND] = np.int8
             mydtypes[TRANSCRIPTIONSTRAND] = 'category'
             mydtypes[MUTATION] = 'category'
         #################################################

         chrBased_mutation_df = pd.read_csv(chrBasedMutationDFFilePath, sep='\t', header=0, dtype=mydtypes)
         # chrBased_mutation_df = pd.read_csv(chrBasedMutationDFFilePath,sep='\t', header=0, dtype=mydtypes,encoding='utf8',engine='python')
         chrBased_mutation_df[SIMULATION_NUMBER] = simulationNumber

        except pd.errors.EmptyDataError:
         chrBased_mutation_df = pd.DataFrame()
    #############################################

    if return_number_of_samples:
        if ((chrBased_mutation_df is not None) and (not chrBased_mutation_df.empty)):
            chrbased_samples = chrBased_mutation_df[SAMPLE].unique()
            return chrbased_samples, chrBased_mutation_df
        else:
            return set(), chrBased_mutation_df
    else:
        return chrBased_mutation_df
##################################################################


##########################################################################################
############### Common functions for Nucleosome Occupancy Analysis starts ################
##########################################################################################


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
#July 26, 2020 Using Numpy Arrays
def writeSimulationBasedAverageNucleosomeOccupancyUsingNumpyArray(occupancy_type,
                                                   sample_based,
                                                   plusorMinus,
                                                   subsSignatures,
                                                   dinucsSignatures,
                                                   indelsSignatures,
                                                   allSims_subsSignature_accumulated_signal_np_array,
                                                   allSims_dinucsSignature_accumulated_signal_np_array,
                                                   allSims_indelsSignature_accumulated_signal_np_array,
                                                   allSims_subsSignature_accumulated_count_np_array,
                                                   allSims_dinucsSignature_accumulated_count_np_array,
                                                   allSims_indelsSignature_accumulated_count_np_array,
                                                   outputDir,
                                                   jobname,
                                                   numofSimulations,
                                                   library_file_memo):

    os.makedirs(os.path.join(outputDir, jobname, DATA, occupancy_type),exist_ok=True)

    ####################################################################################
    for simNum in range(0, numofSimulations + 1):

        subsSignature_accumulated_signal_np_array = allSims_subsSignature_accumulated_signal_np_array[simNum]
        dinucsSignature_accumulated_signal_np_array =  allSims_dinucsSignature_accumulated_signal_np_array[simNum]
        indelsSignature_accumulated_signal_np_array = allSims_indelsSignature_accumulated_signal_np_array[simNum]

        subsSignature_accumulated_count_np_array = allSims_subsSignature_accumulated_count_np_array[simNum]
        dinucsSignature_accumulated_count_np_array = allSims_dinucsSignature_accumulated_count_np_array[simNum]
        indelsSignature_accumulated_count_np_array = allSims_indelsSignature_accumulated_count_np_array[simNum]

        # Last row contains AGGREGATEDSUBSTITUTIONS
        writeAverageNucleosomeOccupancyFiles(occupancy_type,plusorMinus,subsSignature_accumulated_signal_np_array[-1],subsSignature_accumulated_count_np_array[-1],outputDir, jobname,library_file_memo, AGGREGATEDSUBSTITUTIONS,simNum)

        # Last row contains AGGREGATEDDINUCS
        writeAverageNucleosomeOccupancyFiles(occupancy_type,plusorMinus,dinucsSignature_accumulated_signal_np_array[-1],dinucsSignature_accumulated_count_np_array[-1], outputDir,jobname,library_file_memo, AGGREGATEDDINUCS,simNum)

        # Last row contains AGGREGATEDINDELS
        writeAverageNucleosomeOccupancyFiles(occupancy_type,plusorMinus,indelsSignature_accumulated_signal_np_array[-1],indelsSignature_accumulated_count_np_array[-1], outputDir,jobname,library_file_memo, AGGREGATEDINDELS,simNum)

        #Signatures
        writeSignatureBasedAverageNucleosomeOccupancyFilesUsingNumpyArray(occupancy_type,
                                                           plusorMinus,
                                                           subsSignatures,
                                                           dinucsSignatures,
                                                           indelsSignatures,
                                                           subsSignature_accumulated_signal_np_array,
                                                           dinucsSignature_accumulated_signal_np_array,
                                                           indelsSignature_accumulated_signal_np_array,
                                                           subsSignature_accumulated_count_np_array,
                                                           dinucsSignature_accumulated_count_np_array,
                                                           indelsSignature_accumulated_count_np_array,
                                                           outputDir,
                                                           jobname,
                                                           library_file_memo,
                                                           simNum)
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
#July 26, 2020 Using Numpy Array
def writeSignatureBasedAverageNucleosomeOccupancyFilesUsingNumpyArray(occupancy_type,
                                                           plusorMinus,
                                                           subsSignatures,
                                                           dinucsSignatures,
                                                           indelsSignatures,
                                                           subsSignature_accumulated_signal_np_array,
                                                           dinucsSignature_accumulated_signal_np_array,
                                                           indelsSignature_accumulated_signal_np_array,
                                                           subsSignature_accumulated_count_np_array,
                                                           dinucsSignature_accumulated_count_np_array,
                                                           indelsSignature_accumulated_count_np_array,
                                                           outputDir,
                                                           jobname,
                                                           library_file_memo,
                                                           simNum):

    os.makedirs(os.path.join(outputDir,jobname,DATA,occupancy_type,SIGNATUREBASED), exist_ok=True)

    all_signatures=[subsSignatures,dinucsSignatures,indelsSignatures]
    all_signal_arrays=[subsSignature_accumulated_signal_np_array,dinucsSignature_accumulated_signal_np_array,indelsSignature_accumulated_signal_np_array]
    all_count_arrays= [subsSignature_accumulated_count_np_array,dinucsSignature_accumulated_count_np_array,indelsSignature_accumulated_count_np_array]

    for signatures_index,signatures in enumerate(all_signatures,0):
        number_of_signatures = signatures.size
        signal_arrays=all_signal_arrays[signatures_index]
        count_arrays=all_count_arrays[signatures_index]

        #Why -1, because last one is aggregated mutations and it is not written here.
        for signature_index in range(0,number_of_signatures-1):
            signature = signatures[signature_index]
            signalArray = signal_arrays[signature_index]
            countArray = count_arrays[signature_index]
            averageNucleosomeSignalArray = computeAverageNucleosomeOccupancyArray(plusorMinus,signalArray,countArray)

            #To provide filename with no space in signature name
            #signatureWithNoSpace = signature.replace(' ','')

            if (simNum==0):
                if library_file_memo is not None:
                    accumulatedSignalFilename = '%s_%s_AccumulatedSignalArray.txt' % (signature,library_file_memo)
                    accumulatedCountFilename = '%s_%s_AccumulatedCountArray.txt' % (signature,library_file_memo)
                    averageNucleosomeSignalFilename = '%s_%s_AverageSignalArray.txt' % (signature,library_file_memo)
                else:
                    accumulatedSignalFilename = '%s_AccumulatedSignalArray.txt' %(signature)
                    accumulatedCountFilename = '%s_AccumulatedCountArray.txt' %(signature)
                    averageNucleosomeSignalFilename = '%s_AverageSignalArray.txt' %(signature)
            else:
                if library_file_memo is not None:
                    accumulatedSignalFilename = '%s_sim%d_%s_AccumulatedSignalArray.txt' % (signature, simNum,library_file_memo)
                    accumulatedCountFilename = '%s_sim%d_%s_AccumulatedCountArray.txt' % (signature, simNum,library_file_memo)
                    averageNucleosomeSignalFilename = '%s_sim%d_%s_AverageSignalArray.txt' % (signature, simNum,library_file_memo)
                else:
                    accumulatedSignalFilename = '%s_sim%d_AccumulatedSignalArray.txt' %(signature,simNum)
                    accumulatedCountFilename = '%s_sim%d_AccumulatedCountArray.txt' %(signature,simNum)
                    averageNucleosomeSignalFilename = '%s_sim%d_AverageSignalArray.txt' %(signature,simNum)

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

    wig_chrom_start_end_signal_version_df[CHROM] = wig_chrom_start_end_signal_version_df[CHROM].astype('category')
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
#Used by Replication Time and Replication Strand Bias
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
#Lagging_Count Leading_Count
#Transcribed_Count UnTranscribed_Count
def write_sample_based_strand1_strand2_as_dataframe(output_dir,
                                                    jobname,
                                                    num_of_simulations,
                                                    strand_bias,
                                                    all_samples_np_array,
                                                    all_types_np_array,
                                                    all_sims_all_samples_all_types_strand1_np_array,
                                                    all_sims_all_samples_all_types_strand2_np_array):

    if strand_bias==REPLICATIONSTRANDBIAS:
        sample_type_strand1_strand2_ratio_file_name = 'Sample_Type_%s_Strand_Table.txt' %(LAGGING_VERSUS_LEADING)
        strand1_column="lagging_count"
        strand2_column="leading_count"
    elif strand_bias==TRANSCRIPTIONSTRANDBIAS:
        sample_type_strand1_strand2_ratio_file_name = 'Sample_Type_%s_Strand_Table.txt' %(TRANSCRIBED_VERSUS_UNTRANSCRIBED)
        strand1_column="transcribed_count"
        strand2_column="untranscribed_count"

    os.makedirs(os.path.join(output_dir, jobname, DATA, strand_bias, SAMPLE_BASED), exist_ok=True)
    sample_type_strand1_strand2_ratio_file_path = os.path.join(output_dir, jobname, DATA, strand_bias, SAMPLE_BASED,sample_type_strand1_strand2_ratio_file_name)

    with open(sample_type_strand1_strand2_ratio_file_path, 'w') as f:
        i, j, k = all_sims_all_samples_all_types_strand1_np_array.shape
        f.write('sim_num\tsample\ttype\t%s\t%s\n' %(strand1_column,strand2_column))

        for sim_index in range(num_of_simulations+1):
            for sample_index in range(j):
                for type_index in range(k):
                    f.write('%d\t%s\t%s\t%d\t%d\n'% (sim_index,
                                                     all_samples_np_array[sample_index],
                                                     all_types_np_array[type_index],
                                                     all_sims_all_samples_all_types_strand1_np_array[sim_index][sample_index][type_index],
                                                     all_sims_all_samples_all_types_strand2_np_array[sim_index][sample_index][type_index]))
########################################################################

########################################################################
#Main function for type
#Fills a dictionary and writes it as a dataframe
def write_type_strand_bias_np_array_as_dataframe(all_sims_all_types_strand_np_arrays_list,
                                                all_types_np_array,
                                                strand_bias,
                                                strands,
                                                outputDir,
                                                jobname):

    #Fill type2Strand2ListDict using all_sims_all_types_strand_np_arrays_list
    type2Strand2ListDict = {}

    ##################################################################################################
    for strand_index, strand in enumerate(strands,0):
        all_sims_all_types_strand_np_array = all_sims_all_types_strand_np_arrays_list[strand_index]
        num_of_sims, num_of_types = all_sims_all_types_strand_np_array.shape

        for sim_index in range(0,num_of_sims):
            for type_index in range(0,num_of_types):
                my_type = all_types_np_array[type_index]

                if my_type in type2Strand2ListDict:
                    if strand in type2Strand2ListDict[my_type]:
                        strand_list=type2Strand2ListDict[my_type][strand]
                        if (sim_index==0):
                            #Set real_data
                            strand_list[0]=all_sims_all_types_strand_np_array[sim_index, type_index]
                        else:
                            #Append to sims_data_list
                            strand_list[1].append(all_sims_all_types_strand_np_array[sim_index, type_index])
                    else:
                        type2Strand2ListDict[my_type][strand]=[0,[]]
                        if (sim_index==0):
                            type2Strand2ListDict[my_type][strand][0] = all_sims_all_types_strand_np_array[sim_index, type_index]
                        else:
                            type2Strand2ListDict[my_type][strand][1].append(all_sims_all_types_strand_np_array[sim_index, type_index])
                else:
                    type2Strand2ListDict[my_type]={}
                    type2Strand2ListDict[my_type][strand]=[0,[]]
                    if (sim_index==0):
                        type2Strand2ListDict[my_type][strand][0] = all_sims_all_types_strand_np_array[sim_index, type_index]
                    else:
                        type2Strand2ListDict[my_type][strand][1].append(all_sims_all_types_strand_np_array[sim_index, type_index])
    ##################################################################################################

    ##################################################################################################
    # In strand_list we have
    # real_data
    # sims_data_list
    #Add these to information strand_list
    # mean_sims_data
    # min_sims_data
    # max_sims_data
    #strand_list=[real_data, sims_data_list, mean_sims_data, min_sims_data, max_sims_data]
    for my_type in type2Strand2ListDict:
        for strand in type2Strand2ListDict[my_type]:
            sims_data_list=type2Strand2ListDict[my_type][strand][1]
            mean_sims=np.nanmean(sims_data_list)
            min_sims=np.min(sims_data_list)
            max_sims=np.max(sims_data_list)
            type2Strand2ListDict[my_type][strand].append(mean_sims)
            type2Strand2ListDict[my_type][strand].append(min_sims)
            type2Strand2ListDict[my_type][strand].append(max_sims)
    ##################################################################################################


    ##################################################################################################
    #Calculate p-value only
    for my_type in type2Strand2ListDict:
        if strand_bias==REPLICATIONSTRANDBIAS:
            lagging_strand=strands[0]
            leading_strand=strands[1]

            lagging_real_count=type2Strand2ListDict[my_type][lagging_strand][0]
            # lagging_sims_list=type2Strand2ListDict[my_type][lagging_strand][1]
            lagging_sims_mean_count = type2Strand2ListDict[my_type][lagging_strand][2]

            leading_real_count=type2Strand2ListDict[my_type][leading_strand][0]
            # leading_sims_list=type2Strand2ListDict[my_type][leading_strand][1]
            leading_sims_mean_count = type2Strand2ListDict[my_type][leading_strand][2]

            #Calculate p value
            contingency_table_array = [[lagging_real_count, lagging_sims_mean_count], [leading_real_count, leading_sims_mean_count]]
            oddsratio, lagging_versus_leading_p_value = stats.fisher_exact(contingency_table_array)

            #Set p_value
            type2Strand2ListDict[my_type][LAGGING_VERSUS_LEADING_P_VALUE]=lagging_versus_leading_p_value

        elif strand_bias==TRANSCRIPTIONSTRANDBIAS:
            transcribed_strand = strands[0]
            untranscribed_strand = strands[1]
            nontranscribed_strand = strands[2]

            transcribed_real_count = type2Strand2ListDict[my_type][transcribed_strand][0]
            # transcribed_sims_list = type2Strand2ListDict[my_type][transcribed_strand][1]
            transcribed_sims_mean_count = type2Strand2ListDict[my_type][transcribed_strand][2]

            untranscribed_real_count = type2Strand2ListDict[my_type][untranscribed_strand][0]
            # untranscribed_sims_list = type2Strand2ListDict[my_type][untranscribed_strand][1]
            untranscribed_sims_mean_count = type2Strand2ListDict[my_type][untranscribed_strand][2]

            nontranscribed_real_count = type2Strand2ListDict[my_type][nontranscribed_strand][0]
            # nontranscribed_sims_list = type2Strand2ListDict[my_type][nontranscribed_strand][1]
            nontranscribed_sims_mean_count = type2Strand2ListDict[my_type][nontranscribed_strand][2]

            # Calculate p value
            contingency_table_array = [[transcribed_real_count, transcribed_sims_mean_count],[untranscribed_real_count, untranscribed_sims_mean_count]]
            oddsratio, transcribed_versus_untranscribed_p_value = stats.fisher_exact(contingency_table_array)

            genic_real_count = transcribed_real_count + untranscribed_real_count
            genic_sims_mean_count = transcribed_sims_mean_count + untranscribed_sims_mean_count

            # Calculate p value (transcribed + untranscribed) versus nontranscribed
            contingency_table_array = [[genic_real_count, genic_sims_mean_count],[nontranscribed_real_count, nontranscribed_sims_mean_count]]
            oddsratio, genic_versus_intergenic_p_value = stats.fisher_exact(contingency_table_array)

            # Set p_values
            type2Strand2ListDict[my_type][TRANSCRIBED_VERSUS_UNTRANSCRIBED_P_VALUE]=transcribed_versus_untranscribed_p_value
            type2Strand2ListDict[my_type][GENIC_VERSUS_INTERGENIC_P_VALUE]= genic_versus_intergenic_p_value

    #Calculate q-value and significant_strand will be done during plotting figures
    ##################################################################################################

    ##################################################################################################
    if strand_bias == TRANSCRIPTIONSTRANDBIAS:
        type_strand_count_table_file_name = 'Type_%s_Strand_Table.txt' %(TRANSCRIBED_VERSUS_UNTRANSCRIBED)
        type_strand_table_filepath = os.path.join(outputDir, jobname, DATA, strand_bias,type_strand_count_table_file_name)
        write_type_transcription_dataframe(strands, type2Strand2ListDict, jobname , TRANSCRIBED_VERSUS_UNTRANSCRIBED, type_strand_table_filepath)

        type_strand_count_table_file_name = 'Type_%s_Strand_Table.txt' %(GENIC_VERSUS_INTERGENIC)
        type_strand_table_filepath = os.path.join(outputDir, jobname, DATA, strand_bias,type_strand_count_table_file_name)
        write_type_transcription_dataframe(strands, type2Strand2ListDict, jobname, GENIC_VERSUS_INTERGENIC, type_strand_table_filepath)

    elif strand_bias == REPLICATIONSTRANDBIAS:
        type_strand_count_table_file_name = 'Type_%s_Strand_Table.txt' %(LAGGING_VERSUS_LEADING)
        type_strand_table_filepath = os.path.join(outputDir, jobname, DATA, strand_bias,type_strand_count_table_file_name)
        write_type_replication_dataframe(strands, type2Strand2ListDict, jobname, type_strand_table_filepath)
    ##################################################################################################

########################################################################

##############################################
#subfunction for type
def write_type_replication_dataframe(strands, type2Strand2ListDict, cancer_type, filepath):
    # strand_list=[0:real_data, 1:sims_data_list, 2:mean_sims_data, 3:min_sims_data, 4:max_sims_data]

    # strand_list contains
    # [ real_data,
    # sims_data_list,
    # mean_sims_data,
    # min_sims_data,
    # max_sims_data ]

    strand1=strands[0]
    strand2=strands[1]

    strand1_real_data="%s_real_count" %(strand1)
    strand1_sims_data_list = "%s_sims_count_list" % (strand1)
    strand1_mean_sims_data = "%s_mean_sims_count" % (strand1)
    strand1_min_sims_data = "%s_min_sims_count" % (strand1)
    strand1_max_sims_data = "%s_max_sims_count" % (strand1)

    strand2_real_data="%s_real_count" %(strand2)
    strand2_sims_data_list = "%s_sims_count_list" % (strand2)
    strand2_mean_sims_data = "%s_mean_sims_count" % (strand2)
    strand2_min_sims_data = "%s_min_sims_count" % (strand2)
    strand2_max_sims_data = "%s_max_sims_count" % (strand2)

    L = sorted([(cancer_type, my_type,
                 a[strand1][0], a[strand2][0],
                 a[strand1][2], a[strand2][2],
                 a[LAGGING_VERSUS_LEADING_P_VALUE],
                 a[strand1][0], a[strand1][2], a[strand1][3], a[strand1][4], a[strand1][1],
                 a[strand2][0], a[strand2][2], a[strand2][3], a[strand2][4], a[strand2][1])
                for my_type, a in type2Strand2ListDict.items()])
    df = pd.DataFrame(L, columns=['cancer_type', 'type',
                                  strand1_real_data, strand2_real_data,
                                  strand1_mean_sims_data, strand2_mean_sims_data,
                                  LAGGING_VERSUS_LEADING_P_VALUE,
                                  strand1_real_data, strand1_mean_sims_data, strand1_min_sims_data, strand1_max_sims_data, strand1_sims_data_list,
                                  strand2_real_data, strand2_mean_sims_data, strand2_min_sims_data, strand2_max_sims_data, strand2_sims_data_list])
    df.to_csv(filepath, sep='\t', header=True, index=False)
##############################################


########################################################################
#subfunction for type
def write_type_transcription_dataframe(strands, type2Strand2ListDict, cancer_type, strand_bias_subtype, filepath):
    # strand_list=[0:real_data, 1:sims_data_list, 2:mean_sims_data, 3:min_sims_data, 4:max_sims_data]

    # strand_list contains
    # [ real_data,
    # sims_data_list,
    # mean_sims_data,
    # min_sims_data,
    # max_sims_data ]

    strand1=strands[0]
    strand2=strands[1]
    strand3=strands[2]

    strand1_real_data="%s_real_count" %(strand1)
    strand1_sims_data_list = "%s_sims_count_list" % (strand1)
    strand1_mean_sims_data = "%s_mean_sims_count" % (strand1)
    strand1_min_sims_data = "%s_min_sims_count" % (strand1)
    strand1_max_sims_data = "%s_max_sims_count" % (strand1)

    strand2_real_data="%s_real_count" %(strand2)
    strand2_sims_data_list = "%s_sims_count_list" % (strand2)
    strand2_mean_sims_data = "%s_mean_sims_count" % (strand2)
    strand2_min_sims_data = "%s_min_sims_count" % (strand2)
    strand2_max_sims_data = "%s_max_sims_count" % (strand2)

    strand3_real_data = "%s_real_count" %(strand3)
    strand3_sims_data_list = "%s_sims_count_list" %(strand3)
    strand3_mean_sims_data = "%s_mean_sims_count" %(strand3)
    strand3_min_sims_data = "%s_min_sims_count" %(strand3)
    strand3_max_sims_data = "%s_max_sims_count" %(strand3)

    if (strand_bias_subtype==TRANSCRIBED_VERSUS_UNTRANSCRIBED):
        # # Formerly writes Transcribed Untranscribed
        # L = sorted([(cancer_type, my_type,
        #              a[strand1][0], a[strand2][0],
        #              a[strand1][2], a[strand2][2],
        #              a[TRANSCRIBED_VERSUS_UNTRANSCRIBED_P_VALUE],
        #              a[strand1][0], a[strand1][2], a[strand1][3], a[strand1][4], a[strand1][1],
        #              a[strand2][0], a[strand2][2], a[strand2][3], a[strand2][4], a[strand2][1])
        #             for my_type, a in type2Strand2ListDict.items()])
        # df = pd.DataFrame(L, columns=['cancer_type', 'type',
        #                               strand1_real_data, strand2_real_data,
        #                               strand1_mean_sims_data, strand2_mean_sims_data,
        #                               TRANSCRIBED_VERSUS_UNTRANSCRIBED_P_VALUE,
        #                               strand1_real_data, strand1_mean_sims_data, strand1_min_sims_data, strand1_max_sims_data, strand1_sims_data_list,
        #                               strand2_real_data, strand2_mean_sims_data, strand2_min_sims_data, strand2_max_sims_data, strand2_sims_data_list])

        # Now writes Transcribed Untranscribed Nontranscribed
        L = sorted([(cancer_type, my_type,
                     a[strand1][0], a[strand2][0], a[strand3][0],
                     a[strand1][2], a[strand2][2], a[strand3][2],
                     a[TRANSCRIBED_VERSUS_UNTRANSCRIBED_P_VALUE],
                     a[strand1][0], a[strand1][2], a[strand1][3], a[strand1][4], a[strand1][1],
                     a[strand2][0], a[strand2][2], a[strand2][3], a[strand2][4], a[strand2][1],
                     a[strand3][0], a[strand3][2], a[strand3][3], a[strand3][4], a[strand3][1])
                    for my_type, a in type2Strand2ListDict.items()])
        df = pd.DataFrame(L, columns=['cancer_type', 'type',
                                      strand1_real_data, strand2_real_data, strand3_real_data,
                                      strand1_mean_sims_data, strand2_mean_sims_data, strand3_mean_sims_data,
                                      TRANSCRIBED_VERSUS_UNTRANSCRIBED_P_VALUE,
                                      strand1_real_data, strand1_mean_sims_data, strand1_min_sims_data, strand1_max_sims_data, strand1_sims_data_list,
                                      strand2_real_data, strand2_mean_sims_data, strand2_min_sims_data, strand2_max_sims_data, strand2_sims_data_list,
                                      strand3_real_data, strand3_mean_sims_data, strand3_min_sims_data, strand3_max_sims_data, strand3_sims_data_list])

        df.to_csv(filepath, sep='\t', header=True, index=False)

    elif strand_bias_subtype==GENIC_VERSUS_INTERGENIC:
        L = sorted([(cancer_type, my_type,
                     (a[strand1][0] + a[strand2][0]), a[strand3][0],
                     (a[strand1][2] + a[strand2][2]), a[strand3][2],
                     a[GENIC_VERSUS_INTERGENIC_P_VALUE],
                     a[strand1][0], a[strand1][2], a[strand1][3], a[strand1][4], a[strand1][1],
                     a[strand2][0], a[strand2][2], a[strand2][3], a[strand2][4], a[strand2][1],
                     a[strand3][0], a[strand3][2], a[strand3][3], a[strand3][4], a[strand3][1])
                    for my_type, a in type2Strand2ListDict.items()])
        df = pd.DataFrame(L, columns=['cancer_type', 'type',
                                      'genic_real_count', 'intergenic_real_count',
                                      'genic_mean_sims_count', 'intergenic_mean_sims_count',
                                      GENIC_VERSUS_INTERGENIC_P_VALUE,
                                      strand1_real_data, strand1_mean_sims_data, strand1_min_sims_data, strand1_max_sims_data, strand1_sims_data_list,
                                      strand2_real_data, strand2_mean_sims_data, strand2_min_sims_data, strand2_max_sims_data, strand2_sims_data_list,
                                      strand3_real_data, strand3_mean_sims_data, strand3_min_sims_data, strand3_max_sims_data, strand3_sims_data_list])
        df.to_csv(filepath, sep='\t', header=True, index=False)
########################################################################


# Write replication strand bias for real data for each signature
def write_sbs_signature_sbs96_mutation_type_replication_strand_bias(subs_signature_SBS96_mutation_type_lagging_np_array,
                                                                subs_signature_SBS96_mutation_type_leading_np_array,
                                                                SBS96_mutation_types_np_array,
                                                                ordered_sbs_signatures,
                                                                strand_bias,
                                                                outputDir,
                                                                jobname):

    # /restricted/alexandrov-group/burcak/SigProfilerTopographyRuns/Combined_PCAWG_nonPCAWG_4th_iteration/Liver-HCC/data/replication_strand_bias
    for sbs_signature_idx, sbs_signature in enumerate(ordered_sbs_signatures):
        file_name = "%s_%s_real_data.txt" %(sbs_signature, strand_bias)
        with open(os.path.join(outputDir, jobname, DATA, strand_bias, file_name),'w') as writer:
            writer.write('MutationType\tNumber_of_Mutations\n')
            for SBS96_mutation_type_idx, SBS96_mutation_type in enumerate(SBS96_mutation_types_np_array):
                lagging_count = subs_signature_SBS96_mutation_type_lagging_np_array[sbs_signature_idx][SBS96_mutation_type_idx]
                leading_count = subs_signature_SBS96_mutation_type_leading_np_array[sbs_signature_idx][SBS96_mutation_type_idx]
                writer.write('A:' + SBS96_mutation_type + '\t' + str(lagging_count) + '\n')
                writer.write('E:' + SBS96_mutation_type + '\t' + str(leading_count) + '\n')
        writer.close()



########################################################################
# Main function for signature -- mutation type
# Fills a dictionary and writes it as a dataframe
def write_signature_mutation_type_strand_bias_np_array_as_dataframe(all_sims_subs_signature_mutation_type_strand_np_arrays_list,
                                                                    SBS6_mutation_types_np_array,
                                                                    subs_signatures_np_array,
                                                                    strand_bias,
                                                                    strands,
                                                                    outputDir,
                                                                    jobname):

    #Fill signature2MutationType2Strand2ListDict using np_arrays_list
    signature2MutationType2Strand2ListDict = {}

    ##################################################################################################
    for strand_index, strand in enumerate(strands,0):
        all_sims_subs_signature_mutation_type_strand_np_array = all_sims_subs_signature_mutation_type_strand_np_arrays_list[strand_index]
        num_of_sims, num_of_subs_signatures, num_of_mutation_types = all_sims_subs_signature_mutation_type_strand_np_array.shape

        for sim_index in range(0,num_of_sims):
            for subs_signature_index in range(0,num_of_subs_signatures):
                signature = subs_signatures_np_array[subs_signature_index]
                for mutation_type_index in range(0,num_of_mutation_types):
                    mutation_type = SBS6_mutation_types_np_array[mutation_type_index]

                    if signature in signature2MutationType2Strand2ListDict:
                        if mutation_type in signature2MutationType2Strand2ListDict[signature]:
                            if strand in signature2MutationType2Strand2ListDict[signature][mutation_type]:
                                strand_list = signature2MutationType2Strand2ListDict[signature][mutation_type][strand]
                                if sim_index==0:
                                    strand_list[0] = all_sims_subs_signature_mutation_type_strand_np_array[sim_index,subs_signature_index,mutation_type_index]
                                else:
                                    strand_list[1].append(all_sims_subs_signature_mutation_type_strand_np_array[sim_index,subs_signature_index,mutation_type_index])
                            else:
                                signature2MutationType2Strand2ListDict[signature][mutation_type][strand]=[0,[]]
                                if (sim_index==0):
                                    signature2MutationType2Strand2ListDict[signature][mutation_type][strand][0] = all_sims_subs_signature_mutation_type_strand_np_array[sim_index,subs_signature_index,mutation_type_index]
                                else:
                                    signature2MutationType2Strand2ListDict[signature][mutation_type][strand][1].append(all_sims_subs_signature_mutation_type_strand_np_array[sim_index,subs_signature_index,mutation_type_index])
                        else:
                            signature2MutationType2Strand2ListDict[signature][mutation_type]={}
                            signature2MutationType2Strand2ListDict[signature][mutation_type][strand] = [0, []]
                            if (sim_index==0):
                                signature2MutationType2Strand2ListDict[signature][mutation_type][strand][0] = all_sims_subs_signature_mutation_type_strand_np_array[sim_index,subs_signature_index,mutation_type_index]
                            else:
                                signature2MutationType2Strand2ListDict[signature][mutation_type][strand][1].append(all_sims_subs_signature_mutation_type_strand_np_array[sim_index,subs_signature_index,mutation_type_index])

                    else:
                        signature2MutationType2Strand2ListDict[signature]={}
                        signature2MutationType2Strand2ListDict[signature][mutation_type]={}
                        signature2MutationType2Strand2ListDict[signature][mutation_type][strand] = [0, []]
                        if (sim_index==0):
                            signature2MutationType2Strand2ListDict[signature][mutation_type][strand][0] = all_sims_subs_signature_mutation_type_strand_np_array[sim_index,subs_signature_index,mutation_type_index]
                        else:
                            signature2MutationType2Strand2ListDict[signature][mutation_type][strand][1].append(all_sims_subs_signature_mutation_type_strand_np_array[sim_index,subs_signature_index,mutation_type_index])
    ##################################################################################################

    ##################################################################################################
    # In strand_list we have
    # real_data
    # sims_data_list
    #Add these to information strand_list
    # mean_sims_data
    # min_sims_data
    # max_sims_data
    #strand_list=[real_data, sims_data_list, mean_sims_data, min_sims_data, max_sims_data]
    for signature in signature2MutationType2Strand2ListDict:
        for mutation_type in signature2MutationType2Strand2ListDict[signature]:
            for strand in signature2MutationType2Strand2ListDict[signature][mutation_type]:
                sims_data_list=signature2MutationType2Strand2ListDict[signature][mutation_type][strand][1]
                mean_sims=np.nanmean(sims_data_list)
                min_sims=np.min(sims_data_list)
                max_sims=np.max(sims_data_list)
                signature2MutationType2Strand2ListDict[signature][mutation_type][strand].append(mean_sims)
                signature2MutationType2Strand2ListDict[signature][mutation_type][strand].append(min_sims)
                signature2MutationType2Strand2ListDict[signature][mutation_type][strand].append(max_sims)
    ##################################################################################################

    ##################################################################################################
    #Calculate p-value only
    for signature in signature2MutationType2Strand2ListDict:
        for mutation_type in signature2MutationType2Strand2ListDict[signature]:
            if (strand_bias==REPLICATIONSTRANDBIAS):
                lagging_strand=strands[0]
                leading_strand=strands[1]

                lagging_real_count=signature2MutationType2Strand2ListDict[signature][mutation_type][lagging_strand][0]
                lagging_sims_list = signature2MutationType2Strand2ListDict[signature][mutation_type][lagging_strand][1]
                lagging_sims_mean_count = signature2MutationType2Strand2ListDict[signature][mutation_type][lagging_strand][2]

                leading_real_count=signature2MutationType2Strand2ListDict[signature][mutation_type][leading_strand][0]
                leading_sims_list=signature2MutationType2Strand2ListDict[signature][mutation_type][leading_strand][1]
                leading_sims_mean_count = signature2MutationType2Strand2ListDict[signature][mutation_type][leading_strand][2]

                #Calculate p value using Fisher's exact test
                contingency_table_array = [[lagging_real_count, lagging_sims_mean_count], [leading_real_count, leading_sims_mean_count]]
                oddsratio, lagging_versus_leading_p_value = stats.fisher_exact(contingency_table_array)

                #Set p_value
                signature2MutationType2Strand2ListDict[signature][mutation_type][LAGGING_VERSUS_LEADING_P_VALUE]=lagging_versus_leading_p_value

            elif (strand_bias==TRANSCRIPTIONSTRANDBIAS):
                transcribed_strand = strands[0]
                untranscribed_strand = strands[1]
                nontranscribed_strand = strands[2]

                transcribed_real_count = signature2MutationType2Strand2ListDict[signature][mutation_type][transcribed_strand][0]
                transcribed_sims_list = signature2MutationType2Strand2ListDict[signature][mutation_type][transcribed_strand][1]
                transcribed_sims_mean_count = signature2MutationType2Strand2ListDict[signature][mutation_type][transcribed_strand][2]

                untranscribed_real_count = signature2MutationType2Strand2ListDict[signature][mutation_type][untranscribed_strand][0]
                untranscribed_sims_list = signature2MutationType2Strand2ListDict[signature][mutation_type][untranscribed_strand][1]
                untranscribed_sims_mean_count = signature2MutationType2Strand2ListDict[signature][mutation_type][untranscribed_strand][2]

                nontranscribed_real_count = signature2MutationType2Strand2ListDict[signature][mutation_type][nontranscribed_strand][0]
                nontranscribed_sims_list = signature2MutationType2Strand2ListDict[signature][mutation_type][nontranscribed_strand][1]
                nontranscribed_sims_mean_count = signature2MutationType2Strand2ListDict[signature][mutation_type][nontranscribed_strand][2]

                # Calculate p value
                contingency_table_array = [[transcribed_real_count, transcribed_sims_mean_count],[untranscribed_real_count, untranscribed_sims_mean_count]]
                oddsratio, transcribed_versus_untranscribed_p_value = stats.fisher_exact(contingency_table_array)

                genic_real_count = transcribed_real_count + untranscribed_real_count
                genic_sims_mean_count = transcribed_sims_mean_count + untranscribed_sims_mean_count

                #Calculate p value (transcribed + untranscribed) versus nontranscribed
                contingency_table_array = [[genic_real_count, genic_sims_mean_count],[nontranscribed_real_count, nontranscribed_sims_mean_count]]
                oddsratio, genic_versus_intergenic_p_value = stats.fisher_exact(contingency_table_array)

                #Set p_values
                signature2MutationType2Strand2ListDict[signature][mutation_type][TRANSCRIBED_VERSUS_UNTRANSCRIBED_P_VALUE]=transcribed_versus_untranscribed_p_value
                signature2MutationType2Strand2ListDict[signature][mutation_type][GENIC_VERSUS_INTERGENIC_P_VALUE]=genic_versus_intergenic_p_value
    ##################################################################################################

    ##################################################################################################
    if (strand_bias==REPLICATIONSTRANDBIAS):
        signature_mutation_type_strand_count_table_file_name = 'Signature_Mutation_Type_%s_Strand_Table.txt' %(LAGGING_VERSUS_LEADING)
        signature_mutation_type_strand_table_filepath = os.path.join(outputDir, jobname, DATA, strand_bias,signature_mutation_type_strand_count_table_file_name)
        write_signature_mutation_type_replication_dataframe(strands, signature2MutationType2Strand2ListDict,jobname, signature_mutation_type_strand_table_filepath)

    elif (strand_bias==TRANSCRIPTIONSTRANDBIAS):
        signature_mutation_type_strand_count_table_file_name = 'Signature_Mutation_Type_%s_Strand_Table.txt' %(TRANSCRIBED_VERSUS_UNTRANSCRIBED)
        signature_mutation_type_strand_table_filepath = os.path.join(outputDir, jobname, DATA, strand_bias,signature_mutation_type_strand_count_table_file_name)
        write_signature_mutation_type_transcription_dataframe(strands,signature2MutationType2Strand2ListDict,jobname,TRANSCRIBED_VERSUS_UNTRANSCRIBED, signature_mutation_type_strand_table_filepath)

        signature_mutation_type_strand_count_table_file_name = 'Signature_Mutation_Type_%s_Strand_Table.txt' %(GENIC_VERSUS_INTERGENIC)
        signature_mutation_type_strand_table_filepath = os.path.join(outputDir, jobname, DATA, strand_bias,signature_mutation_type_strand_count_table_file_name)
        write_signature_mutation_type_transcription_dataframe(strands,signature2MutationType2Strand2ListDict,jobname, GENIC_VERSUS_INTERGENIC,signature_mutation_type_strand_table_filepath)
    ##################################################################################################

########################################################################



########################################################################
#sub function for signature -- mutation type
def write_signature_mutation_type_transcription_dataframe(strands,signature2MutationType2Strand2ListDict,cancer_type,strand_bias_subtype, filepath):
    # strand_list=[0:real_data, 1:sims_data_list, 2:mean_sims_data, 3:min_sims_data, 4:max_sims_data, 5:p_value]

    strand1=strands[0]
    strand2=strands[1]
    strand3=strands[2]

    # strand_list contains
    # [ real_data,
    # sims_data_list,
    # mean_sims_data,
    # min_sims_data,
    # max_sims_data ]

    strand1_real_data="%s_real_count" %(strand1)
    strand1_sims_data_list = "%s_sims_count_list" % (strand1)
    strand1_mean_sims_data = "%s_mean_sims_count" % (strand1)
    strand1_min_sims_data = "%s_min_sims_count" % (strand1)
    strand1_max_sims_data = "%s_max_sims_count" % (strand1)

    strand2_real_data="%s_real_count" %(strand2)
    strand2_sims_data_list = "%s_sims_count_list" % (strand2)
    strand2_mean_sims_data = "%s_mean_sims_count" % (strand2)
    strand2_min_sims_data = "%s_min_sims_count" % (strand2)
    strand2_max_sims_data = "%s_max_sims_count" % (strand2)

    strand3_real_data = "%s_real_count" %(strand3)
    strand3_sims_data_list = "%s_sims_count_list" %(strand3)
    strand3_mean_sims_data = "%s_mean_sims_count" %(strand3)
    strand3_min_sims_data = "%s_min_sims_count" %(strand3)
    strand3_max_sims_data = "%s_max_sims_count" %(strand3)

    if (strand_bias_subtype==TRANSCRIBED_VERSUS_UNTRANSCRIBED):
        # Before writes Transcribed UnTranscribed
        # L = sorted([(cancer_type, signature, mutation_type,
        #              b[strand1][0], b[strand2][0],
        #              b[strand1][2], b[strand2][2],
        #              b[TRANSCRIBED_VERSUS_UNTRANSCRIBED_P_VALUE],
        #              b[strand1][0], b[strand1][2], b[strand1][3], b[strand1][4], b[strand1][1],
        #              b[strand2][0], b[strand2][2], b[strand2][3], b[strand2][4], b[strand2][1])
        #             for signature, a in signature2MutationType2Strand2ListDict.items()
        #             for mutation_type, b in a.items()])
        # df = pd.DataFrame(L, columns=['cancer_type', 'signature', 'mutation_type',
        #                               strand1_real_data, strand2_real_data,
        #                               strand1_mean_sims_data, strand2_mean_sims_data,
        #                               TRANSCRIBED_VERSUS_UNTRANSCRIBED_P_VALUE,
        #                               strand1_real_data, strand1_mean_sims_data, strand1_min_sims_data, strand1_max_sims_data, strand1_sims_data_list,
        #                               strand2_real_data, strand2_mean_sims_data, strand2_min_sims_data, strand2_max_sims_data, strand2_sims_data_list])
        # df.to_csv(filepath, sep='\t', header=True, index=False)

        # Now writes Transcribed UnTranscribed Nontranscribed
        L = sorted([(cancer_type, signature, mutation_type,
                     b[strand1][0], b[strand2][0], b[strand3][0],
                     b[strand1][2], b[strand2][2], b[strand3][2],
                     b[TRANSCRIBED_VERSUS_UNTRANSCRIBED_P_VALUE],
                     b[strand1][0], b[strand1][2], b[strand1][3], b[strand1][4], b[strand1][1],
                     b[strand2][0], b[strand2][2], b[strand2][3], b[strand2][4], b[strand2][1],
                     b[strand3][0], b[strand3][2], b[strand3][3], b[strand3][4], b[strand3][1])
                    for signature, a in signature2MutationType2Strand2ListDict.items()
                     for mutation_type, b in a.items()])
        df = pd.DataFrame(L, columns=['cancer_type', 'signature', 'mutation_type',
                                      strand1_real_data, strand2_real_data, strand3_real_data,
                                      strand1_mean_sims_data, strand2_mean_sims_data, strand3_mean_sims_data,
                                      TRANSCRIBED_VERSUS_UNTRANSCRIBED_P_VALUE,
                                      strand1_real_data, strand1_mean_sims_data, strand1_min_sims_data, strand1_max_sims_data, strand1_sims_data_list,
                                      strand2_real_data, strand2_mean_sims_data, strand2_min_sims_data, strand2_max_sims_data, strand2_sims_data_list,
                                      strand3_real_data, strand3_mean_sims_data, strand3_min_sims_data, strand3_max_sims_data, strand3_sims_data_list])
        df.to_csv(filepath, sep='\t', header=True, index=False)

    elif (strand_bias_subtype==GENIC_VERSUS_INTERGENIC):
        L = sorted([(cancer_type, signature, mutation_type,
                     (b[strand1][0] + b[strand2][0]), b[strand3][0],
                     (b[strand1][2] + b[strand2][2]), b[strand3][2],
                     b[GENIC_VERSUS_INTERGENIC_P_VALUE],
                     b[strand1][0], b[strand1][2], b[strand1][3], b[strand1][4], b[strand1][1],
                     b[strand2][0], b[strand2][2], b[strand2][3], b[strand2][4], b[strand2][1],
                     b[strand3][0], b[strand3][2], b[strand3][3], b[strand3][4], b[strand3][1])
                    for signature, a in signature2MutationType2Strand2ListDict.items()
                        for mutation_type, b in a.items()])
        df = pd.DataFrame(L, columns=['cancer_type', 'signature', 'mutation_type',
                                      'genic_real_count', 'intergenic_real_count',
                                      'genic_mean_sims_count', 'intergenic_mean_sims_count',
                                      GENIC_VERSUS_INTERGENIC_P_VALUE,
                                      strand1_real_data, strand1_mean_sims_data, strand1_min_sims_data, strand1_max_sims_data, strand1_sims_data_list,
                                      strand2_real_data, strand2_mean_sims_data, strand2_min_sims_data, strand2_max_sims_data, strand2_sims_data_list,
                                      strand3_real_data, strand3_mean_sims_data, strand3_min_sims_data, strand3_max_sims_data, strand3_sims_data_list])
        df.to_csv(filepath, sep='\t', header=True, index=False)
########################################################################


########################################################################
#July 24, 2020 written using write_signature_mutation_type_replication_dataframe
def get_replication_signature_df(strands, signature_of_interest, signature2MutationType2Strand2ListDict,cancer_type):
    strand1 = strands[0]
    strand2 = strands[1]

    strand1_real_data = "%s_real_count" % (strand1)
    strand1_sims_data_list = "%s_sims_count_list" % (strand1)
    strand1_mean_sims_data = "%s_mean_sims_count" % (strand1)
    strand1_min_sims_data = "%s_min_sims_count" % (strand1)
    strand1_max_sims_data = "%s_max_sims_count" % (strand1)

    strand2_real_data = "%s_real_count" % (strand2)
    strand2_sims_data_list = "%s_sims_count_list" % (strand2)
    strand2_mean_sims_data = "%s_mean_sims_count" % (strand2)
    strand2_min_sims_data = "%s_min_sims_count" % (strand2)
    strand2_max_sims_data = "%s_max_sims_count" % (strand2)

    column_names=['cancer_type', 'signature', 'mutation_type',
                                  strand1_real_data, strand2_real_data,
                                  strand1_mean_sims_data, strand2_mean_sims_data,
                                  LAGGING_VERSUS_LEADING_P_VALUE,
                                  strand1_real_data, strand1_mean_sims_data, strand1_min_sims_data,
                                  strand1_max_sims_data, strand1_sims_data_list,
                                  strand2_real_data, strand2_mean_sims_data, strand2_min_sims_data,
                                  strand2_max_sims_data, strand2_sims_data_list]

    L = sorted([(cancer_type, signature, mutation_type,
                 b[strand1][0], b[strand2][0],
                 b[strand1][2], b[strand2][2],
                 b[LAGGING_VERSUS_LEADING_P_VALUE],
                 b[strand1][0], b[strand1][2], b[strand1][3], b[strand1][4], b[strand1][1],
                 b[strand2][0], b[strand2][2], b[strand2][3], b[strand2][4], b[strand2][1])
                for signature, a in signature2MutationType2Strand2ListDict.items()
                 if signature == signature_of_interest
                  for mutation_type, b in a.items()])
    df = pd.DataFrame(L, columns=column_names)

    return df,column_names
########################################################################




########################################################################
#sub function for signature -- mutation type
def write_signature_mutation_type_replication_dataframe(strands,signature2MutationType2Strand2ListDict,cancer_type,filepath):
    # strand_list=[0:real_data, 1:sims_data_list, 2:mean_sims_data, 3:min_sims_data, 4:max_sims_data, 5:p_value]

    # strand_list contains
    # [ real_data,
    # sims_data_list,
    # mean_sims_data,
    # min_sims_data,
    # max_sims_data ]

    strand1=strands[0]
    strand2=strands[1]

    strand1_real_data="%s_real_count" %(strand1)
    strand1_sims_data_list = "%s_sims_count_list" % (strand1)
    strand1_mean_sims_data = "%s_mean_sims_count" % (strand1)
    strand1_min_sims_data = "%s_min_sims_count" % (strand1)
    strand1_max_sims_data = "%s_max_sims_count" % (strand1)

    strand2_real_data="%s_real_count" %(strand2)
    strand2_sims_data_list = "%s_sims_count_list" % (strand2)
    strand2_mean_sims_data = "%s_mean_sims_count" % (strand2)
    strand2_min_sims_data = "%s_min_sims_count" % (strand2)
    strand2_max_sims_data = "%s_max_sims_count" % (strand2)

    L = sorted([(cancer_type, signature, mutation_type,
                 b[strand1][0], b[strand2][0],
                 b[strand1][2], b[strand2][2],
                 b[LAGGING_VERSUS_LEADING_P_VALUE],
                 b[strand1][0], b[strand1][2], b[strand1][3], b[strand1][4], b[strand1][1],
                 b[strand2][0], b[strand2][2], b[strand2][3], b[strand2][4], b[strand2][1])
                for signature, a in signature2MutationType2Strand2ListDict.items()
                    for mutation_type, b in a.items()])
    df = pd.DataFrame(L, columns=['cancer_type', 'signature', 'mutation_type',
                                  strand1_real_data, strand2_real_data,
                                  strand1_mean_sims_data, strand2_mean_sims_data,
                                  LAGGING_VERSUS_LEADING_P_VALUE,
                                  strand1_real_data, strand1_mean_sims_data, strand1_min_sims_data, strand1_max_sims_data, strand1_sims_data_list,
                                  strand2_real_data, strand2_mean_sims_data, strand2_min_sims_data, strand2_max_sims_data, strand2_sims_data_list ])
    df.to_csv(filepath, sep='\t', header=True, index=False)
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


############################################################
#JAN 20, 2020
# From DataPreparationCommons.py
############################################################

############################################################
#Same for Release and old_PCAWG Matlab Probabilities
# example line for original data
#UCEC-US_SP89389 10      2017540 N:AT[T>A]CA     1
# example line for simulated data
#UCEC-US_SP89389_1       10      1575080 T:AT[C>T]TG     1
def readChrBasedMutations(chr_based_mutation_filepath,mutation_type_context,mutation_type_context_for_probabilities):

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
                mutations_with_genomic_positions_df[SAMPLE] = mutations_with_genomic_positions_df[SAMPLE].astype('category')
                mutations_with_genomic_positions_df[CHROM] = mutations_with_genomic_positions_df[CHROM].astype('category')
                mutations_with_genomic_positions_df[START] = mutations_with_genomic_positions_df[START].astype(int)
                mutations_with_genomic_positions_df[MUTATIONLONG] = mutations_with_genomic_positions_df[MUTATIONLONG].astype('category')
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
                mutations_with_genomic_positions_df[SAMPLE] = mutations_with_genomic_positions_df[SAMPLE].astype('category')
                mutations_with_genomic_positions_df[CHROM] = mutations_with_genomic_positions_df[CHROM].astype('category')
                mutations_with_genomic_positions_df[START] = mutations_with_genomic_positions_df[START].astype(int)
                mutations_with_genomic_positions_df[MUTATIONLONG] = mutations_with_genomic_positions_df[MUTATIONLONG].astype('category')
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
            elif(mutation_type_context in SBS_CONTEXTS) and (mutation_type_context_for_probabilities in SBS_CONTEXTS):
                # For SNV MatrixGenerator provides
                # LUAD-US_SP50263 10      440625  U:GG[C>T]AG     -1
                # For SNV Extractor has
                # Sample    Mutation   SBS1    SBS2    SBS3    SBS4    SBS5    SBS13   SBS17a  SBS17b  SBS18   SBS28   SBS40
                # LUAD-US_SP50518 A[C>A]A 0.005281537126491598    1.091660854697097e-06   0.0     0.0     0.12212513236310162     0.0022549935716281717   0.0     0.0     0.0     0.0     0.8703372452779239
                mutations_with_genomic_positions_df.columns = [SAMPLE,CHROM,START,MUTATIONLONG,PYRAMIDINESTRAND]
                mutations_with_genomic_positions_df[SAMPLE] = mutations_with_genomic_positions_df[SAMPLE].astype('category')
                mutations_with_genomic_positions_df[CHROM] = mutations_with_genomic_positions_df[CHROM].astype('category')
                mutations_with_genomic_positions_df[START] = mutations_with_genomic_positions_df[START].astype(int)
                mutations_with_genomic_positions_df[MUTATIONLONG] = mutations_with_genomic_positions_df[MUTATIONLONG].astype('category')
                mutations_with_genomic_positions_df[PYRAMIDINESTRAND] = mutations_with_genomic_positions_df[PYRAMIDINESTRAND].astype(int)
                #Add new column
                # Add Context Column from T:TG[C>T]GC to  G[C>T]G
                # MatrixGenerator generates T:TG[C>T]GC
                # Extractor and PCAWG_Matlab sbs probabilities has G[C>T]G
                # SBS288 has T:A[C>G]A
                mutations_with_genomic_positions_df[TRANSCRIPTIONSTRAND] = mutations_with_genomic_positions_df[MUTATIONLONG].str[0]
                #TODO Make conversion for each possible mutation context that can be in probabilities file.
                if (mutation_type_context_for_probabilities==SBS288):
                    # Since SBS288 probabilities mutation context does not contain B we are assigning B to N
                    # TODO More correct way is to assign half of B as T and other half of B as U
                    mutations_with_genomic_positions_df.loc[mutations_with_genomic_positions_df[MUTATIONLONG].str[0] == 'B', MUTATION] = 'N:' + mutations_with_genomic_positions_df[MUTATIONLONG].str[3:10]
                    mutations_with_genomic_positions_df.loc[mutations_with_genomic_positions_df[MUTATIONLONG].str[0] != 'B', MUTATION] = mutations_with_genomic_positions_df[MUTATIONLONG].str[0:2] + mutations_with_genomic_positions_df[MUTATIONLONG].str[3:10]
                elif (mutation_type_context_for_probabilities==SBS96):
                    mutations_with_genomic_positions_df[MUTATION] = mutations_with_genomic_positions_df[MUTATIONLONG].str[3:10]
                else:
                    #TODO to be updated as needed
                    mutations_with_genomic_positions_df[MUTATION] = mutations_with_genomic_positions_df[MUTATIONLONG].str[3:10]

            #Set dtype as 'category'
            mutations_with_genomic_positions_df[MUTATION]=mutations_with_genomic_positions_df[MUTATION].astype('category')
            mutations_with_genomic_positions_df[TRANSCRIPTIONSTRAND]=mutations_with_genomic_positions_df[TRANSCRIPTIONSTRAND].astype('category')

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
    mutation_type_context_for_probabilities=inputList[5]
    mutation_type_context = inputList[6]
    simNum = inputList[7]
    PCAWG=inputList[8]

    ###############################################################################################
    chr_based_mutation_df = readChrBasedMutations(chr_based_mutation_filepath,mutation_type_context,mutation_type_context_for_probabilities)
    ###############################################################################################

    #######################################################################################
    if ((chr_based_mutation_df is not None) and (mutations_probabilities_df is not None)):

        ############################################################################
        #Step2 SigProfilerTopography Python Package

        #For PCAWG_Matlab
        #Convert CMDI-UK_SP116871_1 --> SP116871 # if(simNum>0): simNum=1 Simulation1
        #Convert CMDI-UK_SP116871 --> SP116871 # simNum=0 Original Data
        if PCAWG:
            chr_based_mutation_df[SAMPLE] = chr_based_mutation_df[SAMPLE].str.split('_', expand=True)[1]
            chr_based_mutation_df[SAMPLE]= chr_based_mutation_df[SAMPLE].astype('category')

        #For Release SigProfilerTopography Python Package
        # For SNV
        # LUAD-US_SP50263 10      440625  U:GG[C>T]AG     -1
        # For ID
        # LUAD-US_SP50263 10      8045169 U:2:Ins:R:5     T       TTC     1
        # For DBS
        # UAD-US_SP50263 10      110099884       Q:T[GC>AG]C     0
        if simNum>=1:
            # Get rid of simulation number at the end
            chr_based_mutation_df[SAMPLE] = chr_based_mutation_df[SAMPLE].str.rsplit('_', 1, expand=True)[0]
            chr_based_mutation_df[SAMPLE]= chr_based_mutation_df[SAMPLE].astype('category')
        ############################################################################


        ############################################################################
        if SAMPLE not in mutations_probabilities_df.columns.values:
            merged_df = pd.merge(chr_based_mutation_df,mutations_probabilities_df, how='inner', left_on=[MUTATION],right_on=[MUTATION])
        else:
            merged_df = pd.merge(chr_based_mutation_df, mutations_probabilities_df, how='inner',left_on=[SAMPLE, MUTATION], right_on=[SAMPLE, MUTATION])
        ############################################################################

        if ((merged_df is not None) and (chr_based_mutation_df.shape[0]!=merged_df.shape[0])):
            print('##############################')
            print('There is a situation/problem. For simNum:%s chr:%s All mutation context type: %s mutations are not merged with signature probabilities'  %(simNum,chrShort,mutation_type_context))
            print('For simNum:%s chr:%s mutation context type:%s chr_based_mutation_df.shape(%d,%d)-- merged_df.shape(%d,%d) ' % (simNum, chrShort, mutation_type_context,chr_based_mutation_df.shape[0],chr_based_mutation_df.shape[1],merged_df.shape[0],merged_df.shape[1]))
            print('Which samples are not merged?: %s' %(set(chr_based_mutation_df[SAMPLE].unique()).difference(set(merged_df[SAMPLE].unique()))))
            temp_df = pd.merge(chr_based_mutation_df, mutations_probabilities_df, how='outer',left_on=[SAMPLE, MUTATION], right_on=[SAMPLE, MUTATION], indicator=True)
            print('which rows of chr_based_mutation_df are not merged?')
            print(temp_df[temp_df['_merge']=='left_only'])
            print("chr_based_mutation_df[MUTATION].unique()")
            print(chr_based_mutation_df[MUTATION].unique())
            print("chr_based_mutation_df[SAMPLE].unique()")
            print(chr_based_mutation_df[SAMPLE].unique())
            print("mutations_probabilities_df[MUTATION].unique()")
            print(mutations_probabilities_df[MUTATION].unique())
            print("mutations_probabilities_df[SAMPLE].unique()")
            print(mutations_probabilities_df[SAMPLE].unique())

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
            # Make MUTATION column one of six mutation types
            # TODO Make conversion for every possible probabilities mutation type context
            if (mutation_type_context_for_probabilities == SBS288):
                # SBS288 has T:A[C>G]A
                merged_df[MUTATION] = merged_df[MUTATION].str[4:7]
            else:
                #SBS96 has A[C>G]A
                merged_df[MUTATION] = merged_df[MUTATION].str[2:5]

            merged_df.to_csv(chr_based_merged_mutations_file_path, sep='\t', header=True, index=False)
        else:
            print('-------------No merge file for sim%d mutation_type_context:%s for chr%s' %(simNum,mutation_type_context,chrShort))
    #######################################################################################


    #######################################################################################
    elif ((chr_based_mutation_df is not None) and (mutations_probabilities_df is None)):

        ############################################################################
        #Step2 SigProfilerTopography Python Package

        #For PCAWG_Matlab
        #Convert CMDI-UK_SP116871_1 --> SP116871 # if(simNum>0): simNum=1 Simulation1
        #Convert CMDI-UK_SP116871 --> SP116871 # simNum=0 Original Data
        if PCAWG:
            chr_based_mutation_df[SAMPLE] = chr_based_mutation_df[SAMPLE].str.split('_', expand=True)[1]

        if simNum>=1:
            # Get rid of simulation number at the end
            chr_based_mutation_df[SAMPLE] = chr_based_mutation_df[SAMPLE].str.rsplit('_', 1, expand=True)[0]
        ############################################################################


        if (mutation_type_context in SBS_CONTEXTS):
            chrBasedMergedMutationsFileName = 'chr%s_%s_for_topography.txt' %(chrShort,SUBS)
        elif (mutation_type_context == DBS):
            chrBasedMergedMutationsFileName = 'chr%s_%s_for_topography.txt' %(chrShort,DINUCS)
        elif (mutation_type_context==ID):
            chrBasedMergedMutationsFileName = 'chr%s_%s_for_topography.txt' %(chrShort,INDELS)

        if (simNum==0):
            chr_based_merged_mutations_file_path = os.path.join(outputDir, jobname, DATA, CHRBASED,chrBasedMergedMutationsFileName)
        else:
            simDir = 'sim%d' %(simNum)
            chr_based_merged_mutations_file_path = os.path.join(outputDir, jobname, DATA, CHRBASED,simDir,chrBasedMergedMutationsFileName)

        # #After test uncomment
        # if ('MutationLong' in merged_df.columns.values):
        #     merged_df.drop(['MutationLong'], inplace=True, axis=1)

        #After merge
        # Make MUTATION column one of six mutation types
        # TODO Make conversion for every possible probabilities mutation type context
        if (mutation_type_context_for_probabilities == SBS288):
            # SBS288 has T:A[C>G]A
            chr_based_mutation_df[MUTATION] = chr_based_mutation_df[MUTATION].str[4:7]
        else:
            #SBS96 has A[C>G]A
            chr_based_mutation_df[MUTATION] = chr_based_mutation_df[MUTATION].str[2:5]

        #write
        chr_based_mutation_df.to_csv(chr_based_merged_mutations_file_path, sep='\t', header=True, index=False)
    #######################################################################################

############################################################


############################################################
def readProbabilities(probabilitiesFile,verbose):

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
    elif ('Sample Names' not in probabilities_df.columns.values) and ('MutationTypes' in probabilities_df.columns.values):
        probabilities_df = pd.read_csv(probabilitiesFile, sep='\t', header=0,dtype={'MutationTypes': 'category'})
        probabilities_df.rename(columns={'MutationTypes': MUTATION}, inplace=True)
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

    if verbose:
        print('\tVerbose Probabilities information starts')
        print('\tVerbose probabilities_df.shape')
        print(probabilities_df.shape)
        print('\tVerbose probabilities_df.head()')
        print(probabilities_df.head())
        print('\tVerbose probabilities_df.dtypes')
        print(probabilities_df.dtypes)

        if (SAMPLE in probabilities_df.columns.values):
            print('\tVerbose Unique samples in probabilities_df')
            print(probabilities_df[SAMPLE].unique())
            print('\tVerbose # of unique samples in probabilities_df: %d\n' %(len(probabilities_df[SAMPLE].unique())))

        if (MUTATION in probabilities_df.columns.values):
            print('\tVerbose Unique MutationTypes in probabilities_df')
            print(probabilities_df[MUTATION].unique())
            print('\tVerbose # of unique mutation types in probabilities_df: %d' %(len(probabilities_df[MUTATION].unique())))
        print('\tVerbose Probabilities information ends')
        print('\tVerbose ##############################')

    return probabilities_df
############################################################

