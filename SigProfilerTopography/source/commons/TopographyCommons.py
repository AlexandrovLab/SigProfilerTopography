# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time and strand bias.
# Copyright (C) 2018 Burcak Otlu

import os
import pandas as pd
import numpy as np
import math
import json
import sys
import multiprocessing
import pickle
import twobitreader
import urllib.request
import shutil

#To handle warnings as errors
# import warnings
# warnings.filterwarnings("error")

current_abs_path = os.path.dirname(os.path.realpath(__file__))

LEADING= 'Leading'
LAGGING = 'Lagging'

UNTRANSCRIBED_STRAND = 'UnTranscribed'
TRANSCRIBED_STRAND = 'Transcribed'
NONTRANSCRIBED_STRAND = 'NonTranscribed'

PLUS = '+'
MINUS = '-'

MAXIMUM_CHROMOSOME_LENGTH = 250000000
GIGABYTE_IN_BYTES = 1000000000

BIGWIG='BIGWIG'
BIGBED='BIGBED'
WIG='WIG'
BED='BED'
NARROWPEAK='narrowpeak'
LIBRARY_FILE_TYPE_OTHER='LIBRARY_FILE_TYPE_OTHER'

BED_6PLUS4='BED6+4'
BED_9PLUS2='BED9+2'

############################################################
#Constraints , Thresholds
SUBSTITUTION_NUMBER_OF_MINIMUM_REQUIRED_MUTATIONS = 5000
INDEL_NUMBER_OF_MINIMUM_REQUIRED_MUTATIONS = 1000
DINUC_NUMBER_OF_MINIMUM_REQUIRED_MUTATIONS = 200

PROCESSIVITY_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD = round(0.9,2)
SUBSTITUTION_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD = round(0.9,2)
INDEL_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD = round(0.5,2)
DINUC_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD = round(0.5,2)
############################################################

SUBSTITUTION_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING = str(SUBSTITUTION_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD)
INDEL_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING = str(INDEL_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD)
DINUC_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING = str(DINUC_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD)

SAMPLE_MMR_DEFICIENT_THRESHOLD = 10000
SUBS_STRAND_BIAS_NUMBER_OF_MUTATIONS_THRESHOLD = 1000

#Global Constants for Nucleosome Occupancy Analysis
# plusOrMinus = 2000
# windowSize = plusOrMinus*2+1

NOTSET = 'notSet'

DISPLAY = 'display'
NODISPLAY = 'nodisplay'

INPUT = 'input'
OUTPUT = 'output'
SIMULATION = 'simulation'
SIMULATIONS = 'simulations'
SIMULATIONS_FOR_TOPOGRAPHY = 'simulations_for_topography'
DATA = 'data'
SIGNAL = 'signal'
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

ENCODE_NUCLEOSOME_GM12878_BIGWIG = 'wgEncodeSydhNsomeGm12878Sig.bigWig'
ENCODE_NUCLEOSOME_K562_BIGWIG = 'wgEncodeSydhNsomeK562Sig.bigWig'

ENCODE_NUCLEOSOME_GM12878_WIG = 'wgEncodeSydhNsomeGm12878Sig.wig'
ENCODE_NUCLEOSOME_K562_WIG = 'wgEncodeSydhNsomeK562Sig.wig'

BIGWIG2WIG = 'bigWigToWig'

GSM923442_HG19_ENCODE_REPLISEQ_MCF7_WAVELET_SIGNAL_WIG = 'GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.wig'
GSM923442_HG19_ENCODE_REPLISEQ_MCF7_VALLEY_BED = 'GSM923442_hg19_wgEncodeUwRepliSeqMcf7ValleysRep1.bed'
GSM923442_HG19_ENCODE_REPLISEQ_MCF7_PEAK_BED = 'GSM923442_hg19_wgEncodeUwRepliSeqMcf7PkRep1.bed'
GRCh37_hg19_NCBIREFSEQCURATED = 'GRCh37_hg19_NCBIRefSeqCurated'
GRCh37_ENSEMBL = 'GRCh37_transcripts.txt'

HG19_CHROM_SIZES = 'hg19.chrom.sizes.txt'
HG38_CHROM_SIZES = 'hg38.chrom.sizes.txt'

HG19_2BIT = 'hg19.2bit'
HG38_2BIT = 'hg38.2bit'

HG19_URL = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit'
HG38_URL = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit'

MM9_2BIT = 'mm9.2bit'
MM10_2BIT = 'mm10.2bit'

MM9_URL = 'http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit'
MM10_URL = 'http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit'

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

MutationType2NumberofMutatiosDictFilename='MutationType2NumberofMutatiosDict.txt'

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

#For Replication
DecileIndex2NumfAttributableBasesDictFilename = 'DecileIndex2NumfAttributableBasesDict.txt'

ONE_DIRECTORY_UP = '..'
AVAILABLE_LIBRARY_FILENAMES = 'AvailableLibraryFiles.txt'
AVAILABLE_LIBRARY_FILENAMES_PATH = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,AVAILABLE_LIBRARY_FILENAMES)
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
PLOTTING = 'plotting'
TRANSCRIPTIONSTRANDBIAS = 'transcriptionstrandbias'
REPLICATIONSTRANDBIAS = 'replicationstrandbias'

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
CHR = 'chr'


COMPUTATION_CHROMOSOMES_SEQUENTIAL = 'COMPUTATION_CHROMOSOMES_SEQUENTIAL'

#For small number of samples
COMPUTATION_ALL_CHROMOSOMES_PARALLEL = 'COMPUTATION_ALL_CHROMOSOMES_PARALLEL'

#For big number of samples
COMPUTATION_CHROMOSOMES_SEQUENTIAL_ALL_SIMULATIONS_PARALLEL = 'COMPUTATION_CHROMOSOMES_SEQUENTIAL_ALL_SIMULATIONS_PARALLEL'
COMPUTATION_CHROMOSOMES_SEQUENTIAL_SIMULATIONS_SEQUENTIAL = 'COMPUTATION_CHROMOSOMES_SEQUENTIAL_SIMULATIONS_SEQUENTIAL'

COMPUTATION_CHROMOSOMES_SEQUENTIAL_CHROMOSOME_SPLITS_PARALLEL = 'COMPUTATION_CHROMOSOMES_SEQUENTIAL_CHROMOSOME_SPLITS_PARALLEL'

############################################
#Column Names
PROJECT = 'Project'
SAMPLE = 'Sample'
GENOME = 'Genome'
CHROM = 'Chrom'
CHROMOSOME = 'Chromosome' # to be depreceated
START = 'Start'
END = 'End'
STRAND = 'Strand'
REF = 'Ref'
ALT = 'Alt'
PYRAMIDINESTRAND = 'PyramidineStrand'
TRANSCRIPTIONSTRAND = 'TranscriptionStrand'
MUTATION = 'Mutation'
MUTATIONS = 'Mutations'
CONTEXT = 'Context'

TRANSCRIPTIONSTRAND = 'TranscriptionStrand'

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

############################################
#Nucleosome Occupancy Data  Columns
chrom= 'chrom'
start= 'start'
end= 'end'
signal= 'signal'

#Roadmap Epigenomics Data Columns
name='name'
score='score'
strand='strand'
############################################

# ############################################
# #data_type
# AVERAGE_SIGNAL='AVERAGE_SIGNAL'
# ACCUMULATED_SIGNAL='ACCUMULATED_SIGNAL'
# ACCUMULATED_COUNT='ACCUMULATED_COUNT'
# ############################################

#For double check purposes
OriginalPaperSignature2Sample2ReplicationStrand2CountDict_Filename = 'OriginalPaperSignature2Sample2ReplicationStrand2CountDict.txt'
OriginalPaperSignature2Sample2TranscriptionStrand2CountDict_Filename = 'OriginalPaperSignature2Sample2TranscriptionStrand2CountDict.txt'

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
ID = 'ID'
DBS= 'DBS'
SBS_CONTEXTS = [SBS96,SBS192,SBS384,SBS1536,SBS3072]

# Used for dictionaries
# MutationType2NumberofMutationsDict keys
SUBS = 'SUBS'
INDELS = 'INDELS'
DINUCS = 'DINUCS'

DEFAULT_HISTONE_OCCUPANCY_FILE = 'ENCFF530PJF_liver_H3K9me3-human.bigBed'
DEFAULT_NUCLEOSOME_OCCUPANCY_FILE = 'wgEncodeSydhNsomeGm12878Sig.bigWig'
DEFAULT_REPLICATION_TIME_SIGNAL_FILE = 'GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.wig'
DEFAULT_REPLICATION_TIME_VALLEY_FILE = 'GSM923442_hg19_wgEncodeUwRepliSeqMcf7ValleysRep1.bed'
DEFAULT_REPLICATION_TIME_PEAK_FILE = 'GSM923442_hg19_wgEncodeUwRepliSeqMcf7PkRep1.bed'

###################################################################
#Works on tscc on python but not from pbs file
#urllib.error.URLError: <urlopen error [Errno 101] Network is unreachable>
def downloadFromWeb(url,filepath_to_be_saved):
    # Download the file from `url` and save it locally under `file_name`:
    with urllib.request.urlopen(url) as response, open(filepath_to_be_saved, 'wb') as out_file:
        shutil.copyfileobj(response, out_file)
###################################################################

# ###################################################################
# def downloadFromWebUsingWGET(url,filepath_to_be_saved):
#     wget.download(url,filepath_to_be_saved)
# ##################################################################

###########################################################
import psutil

def memory_usage():
    pid = os.getpid()
    process = psutil.Process(pid)
    memoryUse1 = process.memory_info()[0]/2.**30  #memory use in GB
    print('************** Current Memory Use: '+ str(round(memoryUse1,2))+" GB *****************\n")
###########################################################

########################################################################################
def getShortNames(chromNamesList):
    return [chrName[3:] for chrName in chromNamesList]
########################################################################################

########################################################################################
def getAvailableLibraryFilenamesList():
    availableLibraryFilenamesPath = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,AVAILABLE_LIBRARY_FILENAMES)

    if (os.path.exists(availableLibraryFilenamesPath)):
        availableLibraryFilenamesList = readAsAList(availableLibraryFilenamesPath)

    return availableLibraryFilenamesList
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
def getDecileIndex2NumberofAttributableBasesDict(replicationTimeFilename_wo_extension):
    decileIndex2NumberofAttributableBasesDict = {}
    decileIndex2NumberofAttributableBasesDictFilePath = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,REPLICATION,replicationTimeFilename_wo_extension,DecileIndex2NumfAttributableBasesDictFilename)

    if (os.path.exists(decileIndex2NumberofAttributableBasesDictFilePath)):
        decileIndex2NumberofAttributableBasesDict = readDictionaryUsingPickle(decileIndex2NumberofAttributableBasesDictFilePath)

    return decileIndex2NumberofAttributableBasesDict
########################################################################################


########################################################################################
#Comment this function
#Get signature2PropertiesListDict
def getSubsSignature2NumberofMutationsDict(outputDir,jobname):
    #Load signaturesWithAtLeast10KEligibleMutations
    subsSignature2NumberofMutationsDict = {}

    SubsSignature2NumberofMutationsDictFilePath = os.path.join(outputDir,jobname,DATA,SubsSignature2NumberofMutationsDictFilename)

    if (os.path.exists(SubsSignature2NumberofMutationsDictFilePath)):
        subsSignature2NumberofMutationsDict = readDictionary(SubsSignature2NumberofMutationsDictFilePath)

    return subsSignature2NumberofMutationsDict
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
def getIndelsSignature2NumberofMutationsDict(outputDir,jobname):
    indelsSignature2NumberofMutationsDict = {}
    indelsSignature2NumberofMutationsDictFilePath = os.path.join(outputDir,jobname,DATA,IndelsSignature2NumberofMutationsDictFilename)

    if (os.path.exists(indelsSignature2NumberofMutationsDictFilePath)):
        indelsSignature2NumberofMutationsDict = readDictionary(indelsSignature2NumberofMutationsDictFilePath)

    return indelsSignature2NumberofMutationsDict
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
def readAsAList(filename):
    with open(filename) as f:
        content = f.readlines()
    # you may also want to remove whitespace characters like `\n` at the end of each line
    lineList = [x.strip() for x in content]

    return lineList
###################################################################

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
        ensembl_transcripts_df = pd.read_table(transcriptsFilenamePath, header=0,sep="\t")

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
def readTranscriptsNCBI():
    transcriptsFilenamePath = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,TRANSCRIPTS,GRCh37_hg19_NCBIREFSEQCURATED)

    if (os.path.exists(transcriptsFilenamePath)):
        GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df = pd.read_table(transcriptsFilenamePath, header = 0, sep="\t")
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
def fillCutoff2Signature2PropertiesListDictionary(outputDir,jobname,chromNamesList,type,cutoffs):
    #Filled in the first part
    #PrpertiesList consists of [number of mutations, sum of probabilities]
    cutoff2Signature2PropertiesListDict={}

    #Filled in the second part
    cutoff2Signature2NumberofMutationsAverageProbabilityListDict={}

    #Filled in the third part
    #PropertiesList=[CufoffProbability NumberofMutations AverageMutationProbability]
    signature2PropertiesListDict={}

    for chrLong in chromNamesList:
        if (type==SUBS):
            chrBased_mutation_df = readChrBasedSubsDF(outputDir,jobname,chrLong,type,0)
        elif (type == INDELS):
            chrBased_mutation_df = readChrBasedMutationsDF(outputDir,jobname,chrLong,type,0)
        elif (type== DINUCS):
            chrBased_mutation_df = readChrBasedMutationsDF(outputDir, jobname, chrLong,type,0)

        if ((chrBased_mutation_df is not None) and (not chrBased_mutation_df.empty)):
            # Sample  Chrom   Start   PyrimidineStrand        Mutation        DBS2    DBS4    DBS6    DBS7    DBS11
            # PD10011a        10      24033661        1       TC>AA   0.0     0.7656325053758131      0.15420390829468886     0.07918943063517644     0.000974155694321615
            signatures = getSignatures(chrBased_mutation_df)


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
        #Fisrt part ends
        #Accumulation ended for each chromosome

    # print('Part1 Results cutoff2Signature2PropertiesListDict')
    # print(cutoff2Signature2PropertiesListDict)

    #Second part starts
    #Second find average cutoff
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


    # print('Part2 Results cutoff2Signature2NumberofMutationsAverageProbabilityListDict')
    # print(cutoff2Signature2NumberofMutationsAverageProbabilityListDict)

    #Set the filenames and number of required mutations
    if (type==SUBS):
        number_of_required_mutations=SUBSTITUTION_NUMBER_OF_MINIMUM_REQUIRED_MUTATIONS
        cutoff2Signature2NumberofMutationsAverageProbabilityListDict_filename = "Cutoff2SubsSignature2NumberofMutationsAverageProbabilityListDict.txt"
        signature2PropertiesList_filename = "SubsSignature2PropertiesListDict.txt"
    elif (type == INDELS):
        number_of_required_mutations=INDEL_NUMBER_OF_MINIMUM_REQUIRED_MUTATIONS
        cutoff2Signature2NumberofMutationsAverageProbabilityListDict_filename = "Cutoff2IndelsSignature2NumberofMutationsAverageProbabilityListDict.txt"
        signature2PropertiesList_filename = "IndelsSignature2PropertiesListDict.txt"
    elif (type== DINUCS):
        number_of_required_mutations=DINUC_NUMBER_OF_MINIMUM_REQUIRED_MUTATIONS
        cutoff2Signature2NumberofMutationsAverageProbabilityListDict_filename = "Cutoff2DinucsSignature2NumberofMutationsAverageProbabilityListDict.txt"
        signature2PropertiesList_filename = "DinucsSignature2PropertiesListDict.txt"

    #Third find the signature based cufoff probability with number of mutations >= required number of mutations  and averega mutation probability >=0.9
    sorted_cutoffs=sorted(cutoff2Signature2NumberofMutationsAverageProbabilityListDict.keys())
    for cutoff in sorted_cutoffs:
        for signature in cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff]:
            if signature not in signature2PropertiesListDict:
                if (cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature][0]>=number_of_required_mutations) and (cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature][1]>=0.9):
                    signature2PropertiesListDict[signature]=[cutoff,np.int(cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature][0]),np.float(cutoff2Signature2NumberofMutationsAverageProbabilityListDict[cutoff][signature][1])]
    #Third part ends

    # print('Part3 Results signature2PropertiesListDict')
    # print(signature2PropertiesListDict)

    #Write the dictionaries
    os.makedirs(os.path.join(outputDir,jobname,DATA), exist_ok=True)

    # TypeError: Object of type int64 is not JSON serializable
    # filePath = os.path.join(outputDir,jobname,DATA,cutoff2Signature2PropertiesListDict_filename)
    # json.dump(cutoff2Signature2PropertiesListDict,open(filePath,'w'))

    # filePath = os.path.join(outputDir,jobname,DATA,cutoff2Signature2NumberofMutationsAverageProbabilityListDict_filename)
    # json.dump(cutoff2Signature2NumberofMutationsAverageProbabilityListDict,open(filePath,'w'))

    filePath = os.path.join(outputDir,jobname,DATA,signature2PropertiesList_filename)
    json.dump(signature2PropertiesListDict,open(filePath,'w'))

    return signature2PropertiesListDict
##################################################################


##################################################################
# We are filling and writing
# Sample2NumberofMutationsDictFilename
# Signature2NumberofMutationsDictFilename
# Sample2Signature2NumberofMutationsDictFilename

# We are writing Signature2NumberofMutationsDictFilename for the signatures in signature2PropertiesListDict using chrBased_mutation_df
# At the end, Signature2NumberofMutationsDictFilename and signature2PropertiesListDict must match
# It is like double check
def fill_mutations_dictionaries_write(outputDir, jobname, chromNamesList, type, mutationType2NumberofMutationsDict,signature2PropertiesListDict):
    sample2NumberofMutationsDict = {}
    signature2NumberofMutationsDict = {}
    sample2Signature2NumberofMutationsDict = {}

    #Fill dictionaries with conditions satisfied
    if (type==SUBS):
        Sample2NumberofMutationsDictFilename = Sample2NumberofSubsDictFilename
        Signature2NumberofMutationsDictFilename = SubsSignature2NumberofMutationsDictFilename
        Sample2Signature2NumberofMutationsDictFilename = Sample2SubsSignature2NumberofMutationsDictFilename
        minimum_number_of_mutations_required = SUBSTITUTION_NUMBER_OF_MINIMUM_REQUIRED_MUTATIONS
    elif (type==INDELS):
        Sample2NumberofMutationsDictFilename = Sample2NumberofIndelsDictFilename
        Signature2NumberofMutationsDictFilename = IndelsSignature2NumberofMutationsDictFilename
        Sample2Signature2NumberofMutationsDictFilename = Sample2IndelsSignature2NumberofMutationsDictFilename
        minimum_number_of_mutations_required = INDEL_NUMBER_OF_MINIMUM_REQUIRED_MUTATIONS
    elif (type==DINUCS):
        Sample2NumberofMutationsDictFilename = Sample2NumberofDinucsDictFilename
        Signature2NumberofMutationsDictFilename = DinucsSignature2NumberofMutationsDictFilename
        Sample2Signature2NumberofMutationsDictFilename = Sample2DinucsSignature2NumberofMutationsDictFilename
        minimum_number_of_mutations_required = DINUC_NUMBER_OF_MINIMUM_REQUIRED_MUTATIONS

    #######################################################################
    for chrLong in chromNamesList:
        if (type==SUBS):
            chrBased_mutation_df = readChrBasedSubsDF(outputDir,jobname,chrLong,type,0)
        elif (type == INDELS):
            chrBased_mutation_df = readChrBasedMutationsDF(outputDir,jobname,chrLong,type,0)
        elif (type== DINUCS):
            chrBased_mutation_df = readChrBasedMutationsDF(outputDir, jobname, chrLong,type,0)

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

                if type in mutationType2NumberofMutationsDict:
                    mutationType2NumberofMutationsDict[type] += number_of_mutations
                else:
                    mutationType2NumberofMutationsDict[type] = number_of_mutations

                #new way
                for signature in signature2PropertiesListDict:
                    number_of_mutations= len(chrBased_mutation_df_sample_group_df[chrBased_mutation_df_sample_group_df[signature]>=float(signature2PropertiesListDict[signature][0])])
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
#Used for Dinucs and Indels
def readChrBasedMutationsDF(outputDir,jobname,chrLong,type,simulationNumber):
    filename = '%s_%s_for_topography.txt' %(chrLong,type)

    if (simulationNumber==0):
        chrBasedMutationDFFilePath = os.path.join(outputDir,jobname,DATA,CHRBASED,filename)
    else:
        simulation = 'sim%s' % (simulationNumber)
        chrBasedMutationDFFilePath = os.path.join(outputDir,jobname,DATA,CHRBASED,simulation,filename)

    chrBased_mutation_df = None
    #############################################
    if (os.path.exists(chrBasedMutationDFFilePath)):
        chrBased_mutation_df = pd.read_table(chrBasedMutationDFFilePath,sep="\t", header=0)
        chrBased_mutation_df[SIMULATION_NUMBER]=simulationNumber
    else:
        print('%s does not exist' %(chrBasedMutationDFFilePath))
    #############################################

    return chrBased_mutation_df
##################################################################

##################################################################
def doesSimulationsAlreadyExits(outputDir,jobname,numofSimulations):
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

##################################################################
#TODO To be deleted
#TODO This was making big dataframe and will be depreceated
def getCombinedChrBasedDF(outputDir, jobname, chrLong,original_chrBased_mutations_df,mutationType,numofSimulations):
    # Simulation case
    # Read also simulated data and append them vertically
    if (original_chrBased_mutations_df is not None):
        frames = []
        frames.append(original_chrBased_mutations_df)
        for simNumber in range(1, numofSimulations + 1):
            # Read simulation data and append
            if (mutationType==SUBS):
                chrBased_mutation_df = readChrBasedSubsDF(outputDir, jobname, chrLong, mutationType, simNumber)
            else:
                chrBased_mutation_df= readChrBasedMutationsDF(outputDir,jobname,chrLong,mutationType,simNumber)
            frames.append(chrBased_mutation_df)

        combined_chrBased_mutations_df = pd.concat(frames, ignore_index=True)
        return combined_chrBased_mutations_df
    else:
        return None
##################################################################

##################################################################
def readChrBasedSubsDF(outputDir,jobname,chrLong,type,simulationNumber):

    filename = '%s_%s_for_topography.txt' % (chrLong, type)

    if (simulationNumber==0):
        chrBasedMutationDFFilePath = os.path.join(outputDir,jobname,DATA,CHRBASED,filename)
    else:
        simDir = 'sim%s' % (simulationNumber)
        chrBasedMutationDFFilePath = os.path.join(outputDir,jobname,DATA,CHRBASED,simDir,filename)

    ##########################################################################################
    if (os.path.exists(chrBasedMutationDFFilePath)):

        #############################################
        only_header_chrBased_mutation_df = pd.read_table(chrBasedMutationDFFilePath, sep="\t", comment='#', nrows=1)
        columnNamesList = list(only_header_chrBased_mutation_df.columns.values)

        # We assume that after the column named 'Mutation' there are the signature columns in tab separated way.
        mutationtIndex = columnNamesList.index(MUTATION)
        signatures = columnNamesList[(mutationtIndex + 1):]
        #############################################

        #################################################
        mydtypes = {}

        #np.float16 Half precision float: sign bit, 5 bits exponent, 10 bits mantissa
        #np.int32   Integer (-2147483648 to 2147483647)
        #np.int8 Byte (-128 to 127)

        #To lower the dataframe size
        for signature in signatures:
            mydtypes[signature] = np.float16

        mydtypes[SAMPLE] = str
        mydtypes[CHROM] = str
        mydtypes[START] = np.int32
        mydtypes[MUTATIONLONG] = str
        mydtypes[PYRAMIDINESTRAND] = np.int8
        mydtypes[TRANSCRIPTIONSTRAND] = str
        mydtypes[MUTATION] = str
        #################################################
    else:
        print('%s does not exist' %(chrBasedMutationDFFilePath))
    ##########################################################################################

    #############################################
    if os.path.exists(chrBasedMutationDFFilePath):
        chrBased_mutation_df = pd.read_table(chrBasedMutationDFFilePath,sep="\t", comment='#',dtype=mydtypes)
        #Add original or simulation column
        #DataType 0 means original data
        #DataType other than 0 means simulation data
        chrBased_mutation_df[SIMULATION_NUMBER]=simulationNumber
        return chrBased_mutation_df
    else:
        return None
##################################################################



##################################################################
#mutationsWithSignatureBasedProbabilitiesFileName 'breast_cancer_mutation_probabilities_final.txt'
def readSubs(subsWithSignatureBasedProbabilitiesFileName):

    subsFilePath = os.path.join(subsWithSignatureBasedProbabilitiesFileName)

    #################################################
    #First read only first row
    subs_df = pd.read_table(subsFilePath, sep="\t", comment='#', dtype={SAMPLE: str, CHROM: str},nrows=1)
    columnNamesList = list(subs_df.columns.values)

    contextIndex = columnNamesList.index(CONTEXT)

    # We assume that after the column named 'Context' there are the signature columns in tab separated way.
    signatures = columnNamesList[(contextIndex + 1):]
    #################################################

    sample2NumberofSubsDict = {}
    subsSignature2NumberofMutationsDict = {}
    sample2SubsSignature2NumberofMutationsDict = {}

    #################################################
    mydtypes = {}
    #np.float16 Half precision float: sign bit, 5 bits exponent, 10 bits mantissa
    #np.int32   Integer (-2147483648 to 2147483647)
    #np.int8 Byte (-128 to 127)

    for signature in signatures:
        mydtypes[signature] = np.float32

    mydtypes[SAMPLE] = str
    mydtypes[CHROM] = str
    mydtypes[START] = np.int32
    mydtypes[END] = np.int32
    mydtypes[PYRAMIDINESTRAND] = np.int8
    mydtypes[MUTATION] = str
    mydtypes[CONTEXT] = str
    #################################################

    #################################################
    # mutation_df = pd.read_table(mutationFilePath, sep="\t", comment='#',dtype={'Sample':str,'Chromosome': str, 'Start': int, 'End':int, 'PyramidineStrand': int, 'Mutation':str, 'Context':str})
    subs_df = pd.read_table(subsFilePath, sep="\t", comment='#',dtype=mydtypes)
    #################################################

    listofSamples = subs_df[SAMPLE].unique()
    print('Number of samples in single point mutations file: %d' %(len(listofSamples)))

    ##############################################################
    for sample in listofSamples:
        numberofSubs =  len(subs_df[subs_df[SAMPLE] == sample])
        if (numberofSubs>=SUBSTITUTION_NUMBER_OF_MINIMUM_REQUIRED_MUTATIONS):
            sample2NumberofSubsDict[sample] = numberofSubs
    ##############################################################

    ##############################################################
    for signature in signatures:
        signaturebased_df = subs_df[subs_df[signature] >= SUBSTITUTION_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD]
        numberofSubs = len(signaturebased_df)
        if (numberofSubs>= SUBSTITUTION_NUMBER_OF_MINIMUM_REQUIRED_MUTATIONS):
            subsSignature2NumberofMutationsDict[signature] = numberofSubs
    ##############################################################


    ##############################################################

    for sample in sample2NumberofSubsDict:
        for signature in subsSignature2NumberofMutationsDict:
            #check if there are at least 10K mutations with probability >= 0.5 for this (sample,signature) pair
            numberofMutations = len(subs_df[ ((subs_df[SAMPLE]==sample) & (subs_df[signature]>= SUBSTITUTION_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD)) ])
            if (numberofMutations>= SUBSTITUTION_NUMBER_OF_MINIMUM_REQUIRED_MUTATIONS):
                if sample in sample2SubsSignature2NumberofMutationsDict:
                    sample2SubsSignature2NumberofMutationsDict[sample][signature] = numberofMutations
                else:
                    sample2SubsSignature2NumberofMutationsDict[sample]={}
                    sample2SubsSignature2NumberofMutationsDict[sample][signature] = numberofMutations
    ##############################################################

    # print('###########################')
    # print('size of %s in %d Bytes -- %f in GB ' %(subsWithSignatureBasedProbabilitiesFileName, sys.getsizeof(subs_df),sys.getsizeof(subs_df)/GIGABYTE_IN_BYTES))
    # print('###########################')

    return signatures, \
           sample2NumberofSubsDict, \
           subsSignature2NumberofMutationsDict, \
           sample2SubsSignature2NumberofMutationsDict, \
           subs_df
##################################################################

##################################################################
def readMutationsAndWriteChrBased(mutationsWithSignatureBasedProbabilitiesFileName):
    signatures, samplesWithAtLeast10KMutationsList, mutation_df = readSubs(mutationsWithSignatureBasedProbabilitiesFileName)

    mutation_df_grouped= mutation_df.groupby(CHROM)

    os.makedirs(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, INPUT, CHRBASED),exist_ok=True)

    #########################################################
    ############### Write signatures starts #################
    #########################################################
    signatures_array = np.array(signatures)

    # Write signatures_array to a file
    filename = 'Signatures_%s' %(mutationsWithSignatureBasedProbabilitiesFileName)
    SignaturesFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, INPUT, filename)

    np.savetxt(SignaturesFile, signatures_array, delimiter='\t', fmt='%s')
    #########################################################
    ############### Write signatures ends ###################
    #########################################################


    #########################################################
    ############### Write Unique Chrnames starts ############
    #########################################################
    # Get the unique chromsome names in mutation_df
    uniqueChrNames = mutation_df[CHROM].unique()
    # The unique values are returned as a NumPy array

    # Write uniqueChrNames to a file
    filename = 'ChrNames_%s' %(mutationsWithSignatureBasedProbabilitiesFileName)
    ChrNamesFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, INPUT, filename)

    # chrnames_df = pd.DataFrame(uniqueChrNames)
    # chrnames_df.to_csv(ChrNamesFile, sep='\t', header=None, index=None)
    np.savetxt(ChrNamesFile, uniqueChrNames, delimiter='\t', fmt='%s')
    #########################################################
    ############### Write Unique Chrnames ends ##############
    #########################################################

    for chr, chrBased_mutation_df in mutation_df_grouped:
        chrBasedMutationFileName = 'chr%s_%s' %(chr,mutationsWithSignatureBasedProbabilitiesFileName)
        chrBasedMutationFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, INPUT, CHRBASED, chrBasedMutationFileName)
        chrBased_mutation_df.to_csv(chrBasedMutationFile, index=None, sep='\t', mode='w')
        print('writeChrBasedNucleosome:%s for %s ends' % (type(chrBased_mutation_df), chr))
##################################################################

##################################################################
#For Parallel Writing indels and single point mutations
def writeChrBasedMutationDF(inputList):
    chr = inputList[0]
    mutationsFileName = inputList[1]
    chrBased_mutation_df = inputList[2]
    outputDir = inputList[3]
    jobname =  inputList[4]

    chrBasedMutationFileName = 'chr%s_%s' % (chr, mutationsFileName)
    chrBasedMutationFile = os.path.join(outputDir,jobname,DATA, CHRBASED,chrBasedMutationFileName)
    # lock.acquire()

    if (chrBased_mutation_df is not None):
        chrBased_mutation_df.to_csv(chrBasedMutationFile, index=None, sep='\t', mode='w')

    # lock.release()
    # print('write for %s ends' %(chr))
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
def accumulateSampleBasedArrays(sample2AllSinglePointMutationsAccumulatedSplitsChrBasedArrayDict,sample2AllSinglePointMutationsSplitArrayDict):
    for sample in sample2AllSinglePointMutationsSplitArrayDict.keys():
        if sample in sample2AllSinglePointMutationsAccumulatedSplitsChrBasedArrayDict:
            sample2AllSinglePointMutationsAccumulatedSplitsChrBasedArrayDict[sample] += sample2AllSinglePointMutationsSplitArrayDict[sample]
        else:
            sample2AllSinglePointMutationsAccumulatedSplitsChrBasedArrayDict[sample] = sample2AllSinglePointMutationsSplitArrayDict[sample]
########################################################################################


########################################################################################
def func_addSignal(window_array, entry_start, entry_end, entry_signal, mutation_row_start,plusOrMinus):
    max_start=max(entry_start,mutation_row_start-plusOrMinus)
    min_end=min(entry_end,mutation_row_start+plusOrMinus)
    window_array[max_start-(mutation_row_start-plusOrMinus):min_end-(mutation_row_start-plusOrMinus)+1]+=entry_signal
    # print('window_array[%d:%d]+=%f' %(max_start-(mutation_row_start-plusOrMinus),min_end-(mutation_row_start-plusOrMinus),entry_signal))
########################################################################################


########################################################################################
#Using list comprehension
#September 18, 2019
#You need to send mutation_row[START], mutation_row[SAMPLE], mutation_row[SIMULATION_NUMBER], and mutation_row[signature]
def fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated_using_pyBigWig_using_list_comp(
        row,
        chrLong,
        library_file,
        chrBasedSignalArray,
        library_file_type,
        signal_index,
        my_upperBound,
        maximum_chrom_size,
        sample2NumberofMutationsDict,
        sample2Signature2NumberofMutationsDict,
        simNum2Type2SignalArrayDict,
        simNum2Type2CountArrayDict,
        simNum2Sample2Type2SignalArrayDict,
        simNum2Sample2Type2CountArrayDict,
        signature2PropertiesListDict,
        my_type,
        plusOrMinus,
        mycolumns):

    window_array=None
    windowSize=plusOrMinus*2+1

    # mycolumns = [SAMPLE, START, SIMULATION_NUMBER]
    mutation_row_sample=row[0]
    mutation_row_start=row[1]
    mutation_row_simulation_number=row[2]

    #Get or fill window_array using Case1, Case2, and Case3
    # Case 1: start is very close to the chromosome start
    if (mutation_row_start<plusOrMinus):
        print('Case 1: start is very close to the chromosome start --- mutation[Start]:%d' %(mutation_row_start))
        if ((library_file_type==BED) or (library_file_type==NARROWPEAK)):
            window_array = chrBasedSignalArray[0:(mutation_row_start + plusOrMinus + 1)]
            window_array = np.pad(window_array, (plusOrMinus - mutation_row_start, 0), 'constant', constant_values=(0, 0))

        elif (library_file_type==BIGWIG):
            #Important: The bigWig format does not support overlapping intervals.
            window_array=library_file.values(chrLong,0,(mutation_row_start+plusOrMinus+1),numpy=True)
            # How do you handle outliers?
            window_array[np.isnan(window_array)] = 0
            window_array[window_array>my_upperBound]=my_upperBound
            window_array = np.pad(window_array, (plusOrMinus - mutation_row_start, 0), 'constant',constant_values=(0, 0))

        elif (library_file_type==BIGBED):
            #We assume that in the 7th column there is signal data
            list_of_entries=library_file.entries(chrLong,0,(mutation_row_start+plusOrMinus+1))
            if list_of_entries is not None:
                window_array = np.zeros((windowSize,),dtype=np.float32)
                # We did not handle outliers for BigBed files.

                #From DNA methylation get the 7th
                # library_file_bed_format==BED_6PLUS4):
                # (713235, 713435, 'Peak_40281\t15\t.\t3.48949\t5.67543\t3.79089\t158')
                #signal_index=3
                #library_file_bed_format==BED_9PLUS2):
                #[(10810, 10811, 'MCF7_NoStarve_B1__GC_\t3\t+\t10810\t10811\t255,0,0\t3\t100'), (10812, 10813, 'MCF7_NoStarve_B1__GC_\t3\t+\t10812\t10813\t255,0,0\t3\t100'), (10815, 10816, 'MCF7_NoStarve_B1__GC_\t3\t+\t10815\t10816\t0,255,0\t3\t0')]
                #signal_index=7
                [(func_addSignal(window_array, entry[0], entry[1], np.float32(entry[2].split()[signal_index]),mutation_row_start, plusOrMinus) if len(entry) >= 3 else (func_addSignal(window_array, entry[0], entry[1], 1, mutation_row_start, plusOrMinus))) for entry in list_of_entries]

    # Case 2: start is very close to the chromosome end
    elif (mutation_row_start+plusOrMinus > maximum_chrom_size):
        print('Case2: start is very close to the chromosome end ---  mutation[Start]:%d' %(mutation_row_start))

        if ((library_file_type==BED) or (library_file_type==NARROWPEAK)):
            window_array = chrBasedSignalArray[(mutation_row_start-plusOrMinus):maximum_chrom_size]
            window_array = np.pad(window_array, (0,mutation_row_start+plusOrMinus-maximum_chrom_size+1),'constant',constant_values=(0,0))

        elif (library_file_type==BIGWIG):
            #Important: The bigWig format does not support overlapping intervals.
            window_array = library_file.values(chrLong,(mutation_row_start-plusOrMinus),maximum_chrom_size,numpy=True)
            # How do you handle outliers?
            window_array[np.isnan(window_array)] = 0
            window_array[window_array>my_upperBound]=my_upperBound
            window_array = np.pad(window_array, (0,mutation_row_start+plusOrMinus-maximum_chrom_size+1),'constant',constant_values=(0,0))

        elif (library_file_type==BIGBED):
            # print('Case2 Debug Sep 5, 2019 %s mutation_row[START]:%d mutation_row[START]-plusOrMinus:%d maximum_chrom_size:%d' %(chrLong,mutation_row[START],mutation_row[START]-plusOrMinus,maximum_chrom_size))
            if ((mutation_row_start-plusOrMinus)<maximum_chrom_size):
                list_of_entries=library_file.entries(chrLong,(mutation_row_start-plusOrMinus),maximum_chrom_size)
                if list_of_entries is not None:
                    window_array = np.zeros((windowSize,),dtype=np.float32)
                    # We did not handle outliers for BigBed files.
                    [(func_addSignal(window_array, entry[0], entry[1], np.float32(entry[2].split()[signal_index]),mutation_row_start,plusOrMinus) if len(entry) >= 3 else (func_addSignal(window_array, entry[0], entry[1],1, mutation_row_start,plusOrMinus))) for entry in list_of_entries]


    #Case 3: No problem
    else:
        if ((library_file_type==BED) or (library_file_type==NARROWPEAK)):
            window_array = chrBasedSignalArray[(mutation_row_start-plusOrMinus):(mutation_row_start+plusOrMinus+1)]

        elif (library_file_type==BIGWIG):
            #Important: You have to go over intervals if there are overlapping intervals.
            window_array = library_file.values(chrLong, (mutation_row_start-plusOrMinus), (mutation_row_start+plusOrMinus+1),numpy=True)
            #How do you handle outliers?
            window_array[np.isnan(window_array)] = 0
            window_array[window_array>my_upperBound]=my_upperBound

        elif (library_file_type==BIGBED):
            # print('Case3 Debug Sep 5, 2019 %s mutation_row[START]:%d mutation_row[START]-plusOrMinus:%d mutation_row[START]+plusOrMinus+1:%d' %(chrLong,mutation_row[START],mutation_row[START]-plusOrMinus,mutation_row[START]+plusOrMinus+1))
            if ((mutation_row_start+plusOrMinus+1)<=maximum_chrom_size):
                list_of_entries=library_file.entries(chrLong, (mutation_row_start-plusOrMinus), (mutation_row_start+plusOrMinus+1))
                if list_of_entries is not None:
                    window_array = np.zeros((windowSize,),dtype=np.float32)
                    # We did not handle outliers for BigBed files.
                    [(func_addSignal(window_array, entry[0], entry[1], np.float32(entry[2].split()[signal_index]),mutation_row_start,plusOrMinus) if len(entry) >= 3 else (func_addSignal(window_array, entry[0], entry[1],1, mutation_row_start,plusOrMinus))) for entry in list_of_entries]

    ##########################################################

    #Get the sample at this mutation_row
    sample = mutation_row_sample
    simulationNumber= mutation_row_simulation_number

    #####################################################
    if simulationNumber not in simNum2Type2SignalArrayDict:
        simNum2Type2SignalArrayDict[simulationNumber] = {}
        simNum2Type2CountArrayDict[simulationNumber] = {}

    type2SignalArrayDict = simNum2Type2SignalArrayDict[simulationNumber]
    type2CountArrayDict =  simNum2Type2CountArrayDict[simulationNumber]
    #####################################################

    #####################################################
    if simulationNumber not in simNum2Sample2Type2SignalArrayDict:
        simNum2Sample2Type2SignalArrayDict[simulationNumber] = {}
        simNum2Sample2Type2CountArrayDict[simulationNumber] = {}
    sample2Type2SignalArrayDict = simNum2Sample2Type2SignalArrayDict[simulationNumber]
    sample2Type2CountArrayDict = simNum2Sample2Type2CountArrayDict[simulationNumber]
    #####################################################

    #Fill dictionaries uisng window_array
    if (window_array is not None) and (np.any(window_array)):
        #TODO: Is there a faster way than using for loop?
        ################# Signatures starts #######################
        #mutation_row[signature] mutation probability for that signature
        #signature2PropertiesListDict[signature][0] cutoff probability for that signature
        for signature in signature2PropertiesListDict:
            indexofSignature = mycolumns.index(signature)
            mutation_row_signature = row[indexofSignature]
            if (mutation_row_signature >= float(signature2PropertiesListDict[signature][0])):
                if (signature in type2SignalArrayDict):
                    type2SignalArrayDict[signature] += window_array
                    type2CountArrayDict[signature] += (window_array>0)
                else:
                    type2SignalArrayDict[signature] = np.zeros(windowSize)
                    type2CountArrayDict[signature] = np.zeros(windowSize, dtype=int)
                    type2SignalArrayDict[signature] += window_array
                    type2CountArrayDict[signature] += (window_array>0)

                ####################################################
                if (sample in sample2Signature2NumberofMutationsDict) and (signature in sample2Signature2NumberofMutationsDict[sample]):
                    if sample in sample2Type2SignalArrayDict:
                        if signature in sample2Type2SignalArrayDict[sample]:
                            sample2Type2SignalArrayDict[sample][signature] += window_array
                            sample2Type2CountArrayDict[sample][signature] += (window_array>0)
                        else:
                            sample2Type2SignalArrayDict[sample][signature] = np.zeros(windowSize)
                            sample2Type2CountArrayDict[sample][signature] = np.zeros(windowSize, dtype=int)
                            sample2Type2SignalArrayDict[sample][signature] += window_array
                            sample2Type2CountArrayDict[sample][signature] += (window_array > 0)

                    else:
                        sample2Type2SignalArrayDict[sample] = {}
                        sample2Type2CountArrayDict[sample] = {}
                        sample2Type2SignalArrayDict[sample][signature] = np.zeros(windowSize)
                        sample2Type2CountArrayDict[sample][signature] = np.zeros(windowSize, dtype=int)
                        sample2Type2SignalArrayDict[sample][signature] += window_array
                        sample2Type2CountArrayDict[sample][signature] += (window_array > 0)
                ####################################################
        ################# Signatures ends #########################

        ######################################################################
        if my_type in type2SignalArrayDict:
            type2SignalArrayDict[my_type] += window_array
            type2CountArrayDict[my_type] += (window_array > 0)
        else:
            type2SignalArrayDict[my_type] = np.zeros(windowSize)
            type2CountArrayDict[my_type] = np.zeros(windowSize, dtype=int)
            type2SignalArrayDict[my_type] += window_array
            type2CountArrayDict[my_type] += (window_array > 0)

        if (sample in sample2NumberofMutationsDict):
            if sample in sample2Type2SignalArrayDict:
                if my_type in sample2Type2SignalArrayDict[sample]:
                    sample2Type2SignalArrayDict[sample][my_type] += window_array
                    sample2Type2CountArrayDict[sample][my_type] += (window_array > 0)
                else:
                    sample2Type2SignalArrayDict[sample][my_type] = np.zeros(windowSize)
                    sample2Type2CountArrayDict[sample][my_type] = np.zeros(windowSize, dtype=int)
                    sample2Type2SignalArrayDict[sample][my_type] += window_array
                    sample2Type2CountArrayDict[sample][my_type] += (window_array > 0)
            else:
                sample2Type2SignalArrayDict[sample] = {}
                sample2Type2CountArrayDict[sample] = {}
                sample2Type2SignalArrayDict[sample][my_type] = np.zeros(windowSize)
                sample2Type2CountArrayDict[sample][my_type] = np.zeros(windowSize, dtype=int)
                sample2Type2SignalArrayDict[sample][my_type] += window_array
                sample2Type2CountArrayDict[sample][my_type] += (window_array > 0)
        ######################################################################


########################################################################################



########################################################################################
#Using pyBigWig for bigBed and bigWig files starts
#Using bed files prepared on the fly starts
#Original data and simulations data
#SimNum 0 means original
#SimNum 1 means sim1
#my_type can be AGGREGATEDSUBSTITUTIONS,AGGREGATEDINDELS and AGGREGATEDDINUCS
#Signatures are also treated as a type
def fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated_using_pyBigWig(mutation_row,
                        chrLong,
                        library_file,
                        chrBasedSignalArray,
                        library_file_type,
                        my_upperBound,
                        maximum_chrom_size,
                        sample2NumberofMutationsDict,
                        sample2Signature2NumberofMutationsDict,
                        simNum2Type2SignalArrayDict,
                        simNum2Type2CountArrayDict,
                        simNum2Sample2Type2SignalArrayDict,
                        simNum2Sample2Type2CountArrayDict,
                        signature2PropertiesListDict,
                        my_type,
                        plusOrMinus):

    #signature2PropertiesListDict PropertiesList[ cutoff_probability, number_of_mutations_with_probability_ge_cutoff_probability, average_probability]
    window_array=None
    windowSize=plusOrMinus*2+1

    #Get or fill window_array using Case1, Case2, and Case3
    # Case 1: start is very close to the chromosome start
    if (mutation_row[START]<plusOrMinus):
        print('Case 1: start is very close to the chromosome start --- mutation[Start]:%d' %(mutation_row[START]))

        if (library_file_type==BED):
            window_array = chrBasedSignalArray[0:(mutation_row[START] + plusOrMinus + 1)]
            window_array = np.pad(window_array, (plusOrMinus - mutation_row[START], 0), 'constant', constant_values=(0, 0))

        elif (library_file_type==BIGWIG):
            #Important: The bigWig format does not support overlapping intervals.
            window_array=library_file.values(chrLong,0,(mutation_row[START]+plusOrMinus+1),numpy=True)
            # How do you handle outliers?
            window_array[np.isnan(window_array)] = 0
            window_array[window_array>my_upperBound]=my_upperBound
            window_array = np.pad(window_array, (plusOrMinus - mutation_row[START], 0), 'constant',constant_values=(0, 0))

        elif (library_file_type==BIGBED):
            # print('Case1 Debug Sep 5, 2019 %s mutation_row[START]:%d mutation_row[START]+plusOrMinus+1:%d' %(chrLong,mutation_row[START],mutation_row[START]+plusOrMinus+1))
            list_of_tuples=library_file.entries(chrLong,0,(mutation_row[START]+plusOrMinus+1))
            if list_of_tuples is not None:
                window_array = np.zeros((windowSize,),dtype=np.float32)
                window_array_start = mutation_row[START] - plusOrMinus

                for i in list_of_tuples:
                    #In each tuple there is start end  and space separated string column
                    if len(i)>=3:
                        window_array[i[0]-window_array_start:i[1]-window_array_start]+=np.float32(i[2].split()[1])
                    else:
                        window_array[i[0]-window_array_start:i[1]-window_array_start]+=1
                # We did not handle outliers for BigBed files.

    # Case 2: start is very close to the chromosome end
    elif (mutation_row[START]+plusOrMinus > maximum_chrom_size):
        print('Case2: start is very close to the chromosome end ---  mutation[Start]:%d' %(mutation_row[START]))

        if (library_file_type==BED):
            window_array = chrBasedSignalArray[(mutation_row[START]-plusOrMinus):maximum_chrom_size]
            window_array = np.pad(window_array, (0,mutation_row[START]+plusOrMinus-maximum_chrom_size+1),'constant',constant_values=(0,0))

        elif (library_file_type==BIGWIG):
            #Important: The bigWig format does not support overlapping intervals.
            window_array = library_file.values(chrLong,(mutation_row[START]-plusOrMinus),maximum_chrom_size,numpy=True)
            # How do you handle outliers?
            window_array[np.isnan(window_array)] = 0
            window_array[window_array>my_upperBound]=my_upperBound
            window_array = np.pad(window_array, (0,mutation_row[START]+plusOrMinus-maximum_chrom_size+1),'constant',constant_values=(0,0))

        elif (library_file_type==BIGBED):
            # print('Case2 Debug Sep 5, 2019 %s mutation_row[START]:%d mutation_row[START]-plusOrMinus:%d maximum_chrom_size:%d' %(chrLong,mutation_row[START],mutation_row[START]-plusOrMinus,maximum_chrom_size))
            if ((mutation_row[START]-plusOrMinus)<maximum_chrom_size):
                list_of_tuples=library_file.entries(chrLong,(mutation_row[START]-plusOrMinus),maximum_chrom_size)
                if list_of_tuples is not None:
                    window_array = np.zeros((windowSize,),dtype=np.float32)
                    window_array_start = mutation_row[START] - plusOrMinus
                    for i in list_of_tuples:
                        #In each tuple there is start end  and space separated string column
                        if len(i)>=3:
                            window_array[i[0]-window_array_start:i[1]-window_array_start]+=np.float32(i[2].split()[1])
                        else:
                            window_array[i[0]-window_array_start:i[1]-window_array_start]+=1
                    # We did not handle outliers for BigBed files.

    #Case 3: No problem
    else:
        if (library_file_type==BED):
            window_array = chrBasedSignalArray[(mutation_row[START]-plusOrMinus):(mutation_row[START]+plusOrMinus+1)]

        elif (library_file_type==BIGWIG):
            #Important: You have to go over intervals if there are overlapping intervals.
            window_array = library_file.values(chrLong, (mutation_row[START]-plusOrMinus), (mutation_row[START]+plusOrMinus+1),numpy=True)
            #How do you handle outliers?
            window_array[np.isnan(window_array)] = 0
            window_array[window_array>my_upperBound]=my_upperBound

        elif (library_file_type==BIGBED):
            # print('Case3 Debug Sep 5, 2019 %s mutation_row[START]:%d mutation_row[START]-plusOrMinus:%d mutation_row[START]+plusOrMinus+1:%d' %(chrLong,mutation_row[START],mutation_row[START]-plusOrMinus,mutation_row[START]+plusOrMinus+1))
            if ((mutation_row[START]+plusOrMinus+1)<=maximum_chrom_size):
                list_of_tuples=library_file.entries(chrLong, (mutation_row[START]-plusOrMinus), (mutation_row[START]+plusOrMinus+1))
                if list_of_tuples is not None:
                    window_array = np.zeros((windowSize,),dtype=np.float32)
                    window_array_start=mutation_row[START] - plusOrMinus
                    for i in list_of_tuples:
                        #In each tuple there is start end  and space separated string column
                        if len(i)>=3:
                            window_array[i[0]-window_array_start:i[1]-window_array_start]+=np.float32(i[2].split()[1])
                        else:
                            window_array[i[0]-window_array_start:i[1]-window_array_start]+=1
                    # We did not handle outliers for BigBed files.

    ##########################################################

    #Get the sample at this mutation_row
    sample = mutation_row[SAMPLE]
    simulationNumber= mutation_row[SIMULATION_NUMBER]

    #####################################################
    if simulationNumber not in simNum2Type2SignalArrayDict:
        simNum2Type2SignalArrayDict[simulationNumber] = {}
        simNum2Type2CountArrayDict[simulationNumber] = {}

    type2SignalArrayDict = simNum2Type2SignalArrayDict[simulationNumber]
    type2CountArrayDict =  simNum2Type2CountArrayDict[simulationNumber]
    #####################################################

    #####################################################
    if simulationNumber not in simNum2Sample2Type2SignalArrayDict:
        simNum2Sample2Type2SignalArrayDict[simulationNumber] = {}
        simNum2Sample2Type2CountArrayDict[simulationNumber] = {}
    sample2Type2SignalArrayDict = simNum2Sample2Type2SignalArrayDict[simulationNumber]
    sample2Type2CountArrayDict = simNum2Sample2Type2CountArrayDict[simulationNumber]
    #####################################################

    #Fill dictionaries uisng window_array
    if (window_array is not None) and (np.any(window_array)):
        #TODO: Is there a faster way than using for loop?
        ################# Signatures starts #######################
        #mutation_row[signature] mutation probability for that signature
        #signature2PropertiesListDict[signature][0] cutoff probability for that signature
        for signature in signature2PropertiesListDict:
            if (mutation_row[signature] >= float(signature2PropertiesListDict[signature][0])):
                if (signature in type2SignalArrayDict):
                    type2SignalArrayDict[signature] += window_array
                    type2CountArrayDict[signature] += (window_array>0)
                else:
                    type2SignalArrayDict[signature] = np.zeros(windowSize)
                    type2CountArrayDict[signature] = np.zeros(windowSize, dtype=int)
                    type2SignalArrayDict[signature] += window_array
                    type2CountArrayDict[signature] += (window_array>0)

                ####################################################
                if (sample in sample2Signature2NumberofMutationsDict) and (signature in sample2Signature2NumberofMutationsDict[sample]):
                    if sample in sample2Type2SignalArrayDict:
                        if signature in sample2Type2SignalArrayDict[sample]:
                            sample2Type2SignalArrayDict[sample][signature] += window_array
                            sample2Type2CountArrayDict[sample][signature] += (window_array>0)
                        else:
                            sample2Type2SignalArrayDict[sample][signature] = np.zeros(windowSize)
                            sample2Type2CountArrayDict[sample][signature] = np.zeros(windowSize, dtype=int)
                            sample2Type2SignalArrayDict[sample][signature] += window_array
                            sample2Type2CountArrayDict[sample][signature] += (window_array > 0)

                    else:
                        sample2Type2SignalArrayDict[sample] = {}
                        sample2Type2CountArrayDict[sample] = {}
                        sample2Type2SignalArrayDict[sample][signature] = np.zeros(windowSize)
                        sample2Type2CountArrayDict[sample][signature] = np.zeros(windowSize, dtype=int)
                        sample2Type2SignalArrayDict[sample][signature] += window_array
                        sample2Type2CountArrayDict[sample][signature] += (window_array > 0)
                ####################################################
        ################# Signatures ends #########################

        ######################################################################
        if my_type in type2SignalArrayDict:
            type2SignalArrayDict[my_type] += window_array
            type2CountArrayDict[my_type] += (window_array > 0)
        else:
            type2SignalArrayDict[my_type] = np.zeros(windowSize)
            type2CountArrayDict[my_type] = np.zeros(windowSize, dtype=int)
            type2SignalArrayDict[my_type] += window_array
            type2CountArrayDict[my_type] += (window_array > 0)

        if (sample in sample2NumberofMutationsDict):
            if sample in sample2Type2SignalArrayDict:
                if my_type in sample2Type2SignalArrayDict[sample]:
                    sample2Type2SignalArrayDict[sample][my_type] += window_array
                    sample2Type2CountArrayDict[sample][my_type] += (window_array > 0)
                else:
                    sample2Type2SignalArrayDict[sample][my_type] = np.zeros(windowSize)
                    sample2Type2CountArrayDict[sample][my_type] = np.zeros(windowSize, dtype=int)
                    sample2Type2SignalArrayDict[sample][my_type] += window_array
                    sample2Type2CountArrayDict[sample][my_type] += (window_array > 0)
            else:
                sample2Type2SignalArrayDict[sample] = {}
                sample2Type2CountArrayDict[sample] = {}
                sample2Type2SignalArrayDict[sample][my_type] = np.zeros(windowSize)
                sample2Type2CountArrayDict[sample][my_type] = np.zeros(windowSize, dtype=int)
                sample2Type2SignalArrayDict[sample][my_type] += window_array
                sample2Type2CountArrayDict[sample][my_type] += (window_array > 0)
        ######################################################################

#Using pyBigWig for bigBed and bigWig files ends
#Using bed files prepared on the fly ends
########################################################################################


########################################################################################
#Depreceated. Now using fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated_using_pyBigWig
#Original data and simulations data
#SimNum 0 means original
#SimNum 1 means sim1
#my_type can be AGGREGATEDSUBSTITUTIONS,AGGREGATEDINDELS and AGGREGATEDDINUCS
#Signatures are also treated as a type
def fillSignalArrayAndCountArrayForMutationsSimulationsIntegrated(mutation_row,
                        chrBasedSignalArray,
                        maximum_chrom_size,
                        sample2NumberofMutationsDict,
                        sample2Signature2NumberofMutationsDict,
                        simNum2Type2SignalArrayDict,
                        simNum2Type2CountArrayDict,
                        simNum2Sample2Type2SignalArrayDict,
                        simNum2Sample2Type2CountArrayDict,
                        signature2PropertiesListDict,
                        my_type,
                        plusOrMinus):


    windowSize=plusOrMinus*2+1
    #signature2PropertiesListDict PropertiesList[ cutoff_probability, number_of_mutations_with_probability_ge_cutoff_probability, average_probability]

    # Case 1: start is very close to the chromosome start
    if (mutation_row[START]<plusOrMinus):
        print('Case 1: start is very close to the chromosome start --- mutation[Start]:%d' %(mutation_row[START]))
        window_array = chrBasedSignalArray[0:(mutation_row[START]+plusOrMinus+1)]
        window_array = np.pad(window_array, (plusOrMinus-mutation_row[START],0),'constant',constant_values=(0,0))
    # Case 2: start is very close to the chromosome end
    elif (mutation_row[START]+plusOrMinus > maximum_chrom_size):
        print('Case2: start is very close to the chromosome end ---  mutation[Start]:%d' %(mutation_row[START]))
        window_array = chrBasedSignalArray[(mutation_row[START]-plusOrMinus):maximum_chrom_size]
        window_array = np.pad(window_array, (0,mutation_row[START]+plusOrMinus-maximum_chrom_size+1),'constant',constant_values=(0,0))
    #Case 3: No problem
    else:
        window_array = chrBasedSignalArray[(mutation_row[START]-plusOrMinus):(mutation_row[START]+plusOrMinus+1)]
    ##########################################################

    #Get the sample at this mutation_row
    sample = mutation_row[SAMPLE]
    simulationNumber= mutation_row[SIMULATION_NUMBER]

    #####################################################
    if simulationNumber not in simNum2Type2SignalArrayDict:
        simNum2Type2SignalArrayDict[simulationNumber] = {}
        simNum2Type2CountArrayDict[simulationNumber] = {}

    type2SignalArrayDict = simNum2Type2SignalArrayDict[simulationNumber]
    type2CountArrayDict =  simNum2Type2CountArrayDict[simulationNumber]
    #####################################################

    #####################################################
    if simulationNumber not in simNum2Sample2Type2SignalArrayDict:
        simNum2Sample2Type2SignalArrayDict[simulationNumber] = {}
        simNum2Sample2Type2CountArrayDict[simulationNumber] = {}
    sample2Type2SignalArrayDict = simNum2Sample2Type2SignalArrayDict[simulationNumber]
    sample2Type2CountArrayDict = simNum2Sample2Type2CountArrayDict[simulationNumber]
    #####################################################

    #TODO: Is there a faster way than using for loop?
    ################# Signatures starts #######################
    #mutation_row[signature] mutation probability for that signature
    #signature2PropertiesListDict[signature][0] cutoff probability for that signature
    for signature in signature2PropertiesListDict:
        if (mutation_row[signature] >= float(signature2PropertiesListDict[signature][0])):
            if (signature in type2SignalArrayDict):
                type2SignalArrayDict[signature] += window_array
                type2CountArrayDict[signature] += (window_array>0)
            else:
                type2SignalArrayDict[signature] = np.zeros(windowSize)
                type2CountArrayDict[signature] = np.zeros(windowSize, dtype=int)
                type2SignalArrayDict[signature] += window_array
                type2CountArrayDict[signature] += (window_array>0)

            ####################################################
            if (sample in sample2Signature2NumberofMutationsDict) and (signature in sample2Signature2NumberofMutationsDict[sample]):
                if sample in sample2Type2SignalArrayDict:
                    if signature in sample2Type2SignalArrayDict[sample]:
                        sample2Type2SignalArrayDict[sample][signature] += window_array
                        sample2Type2CountArrayDict[sample][signature] += (window_array>0)
                    else:
                        sample2Type2SignalArrayDict[sample][signature] = np.zeros(windowSize)
                        sample2Type2CountArrayDict[sample][signature] = np.zeros(windowSize, dtype=int)
                        sample2Type2SignalArrayDict[sample][signature] += window_array
                        sample2Type2CountArrayDict[sample][signature] += (window_array > 0)

                else:
                    sample2Type2SignalArrayDict[sample] = {}
                    sample2Type2CountArrayDict[sample] = {}
                    sample2Type2SignalArrayDict[sample][signature] = np.zeros(windowSize)
                    sample2Type2CountArrayDict[sample][signature] = np.zeros(windowSize, dtype=int)
                    sample2Type2SignalArrayDict[sample][signature] += window_array
                    sample2Type2CountArrayDict[sample][signature] += (window_array > 0)
            ####################################################
    ################# Signatures ends #########################

    ######################################################################
    if my_type in type2SignalArrayDict:
        type2SignalArrayDict[my_type] += window_array
        type2CountArrayDict[my_type] += (window_array > 0)
    else:
        type2SignalArrayDict[my_type] = np.zeros(windowSize)
        type2CountArrayDict[my_type] = np.zeros(windowSize, dtype=int)
        type2SignalArrayDict[my_type] += window_array
        type2CountArrayDict[my_type] += (window_array > 0)

    if (sample in sample2NumberofMutationsDict):
        if sample in sample2Type2SignalArrayDict:
            if my_type in sample2Type2SignalArrayDict[sample]:
                sample2Type2SignalArrayDict[sample][my_type] += window_array
                sample2Type2CountArrayDict[sample][my_type] += (window_array > 0)
            else:
                sample2Type2SignalArrayDict[sample][my_type] = np.zeros(windowSize)
                sample2Type2CountArrayDict[sample][my_type] = np.zeros(windowSize, dtype=int)
                sample2Type2SignalArrayDict[sample][my_type] += window_array
                sample2Type2CountArrayDict[sample][my_type] += (window_array > 0)
        else:
            sample2Type2SignalArrayDict[sample] = {}
            sample2Type2CountArrayDict[sample] = {}
            sample2Type2SignalArrayDict[sample][my_type] = np.zeros(windowSize)
            sample2Type2CountArrayDict[sample][my_type] = np.zeros(windowSize, dtype=int)
            sample2Type2SignalArrayDict[sample][my_type] += window_array
            sample2Type2CountArrayDict[sample][my_type] += (window_array > 0)
    ######################################################################


########################################################################################


########################################################################################
# Use only one fillSignalArrayAndCountArrayForMutations method
def fillSignalArrayAndCountArrayForMutations(mutation_row,
                        nucleosome_array,
                        maximum_chrom_size,
                        sample2NumberofMutationsDict,
                        signature2NumberofMutationsDict,
                        sample2Signature2NumberofMutationsDict,
                        type2SignalArrayDict,
                        type2CountArrayDict,
                        sample2Type2SignalArrayDict,
                        sample2Type2CountArrayDict,
                        MUTATION_SIGNATURE_PROBABILITY_THRESHOLD,
                        type,
                        plusOrMinus):

    windowSize=plusOrMinus*2+1

    # Case 1: start is very close to the chromosome start
    if (mutation_row[START]<plusOrMinus):
        print('Case 1: start is very close to the chromosome start --- mutation[Start]:%d' %(mutation_row[START]))
        window_array = nucleosome_array[0:(mutation_row[START]+plusOrMinus+1)]
        window_array = np.pad(window_array, (plusOrMinus-mutation_row[START],0),'constant',constant_values=(0,0))
    # Case 2: start is very close to the chromosome end
    elif (mutation_row[START]+plusOrMinus > maximum_chrom_size):
        print('Case2: start is very close to the chromosome end ---  mutation[Start]:%d' %(mutation_row[START]))
        window_array = nucleosome_array[(mutation_row[START]-plusOrMinus):maximum_chrom_size]
        window_array = np.pad(window_array, (0,mutation_row[START]+plusOrMinus-maximum_chrom_size+1),'constant',constant_values=(0,0))
    #Case 3: No problem
    else:
        window_array = nucleosome_array[(mutation_row[START]-plusOrMinus):(mutation_row[START]+plusOrMinus+1)]
    ##########################################################

    #Get the sample at this mutation_row
    sample = mutation_row[SAMPLE]

    #TODO: Is there a faster way than using for loop?
    ################# Signatures starts #######################
    for signature in signature2NumberofMutationsDict:
        if (mutation_row[signature] >= MUTATION_SIGNATURE_PROBABILITY_THRESHOLD):
            if (signature in type2SignalArrayDict):
                type2SignalArrayDict[signature] += window_array
                type2CountArrayDict[signature] += (window_array>0)
            else:
                type2SignalArrayDict[signature] = np.zeros(windowSize)
                type2CountArrayDict[signature] = np.zeros(windowSize, dtype=int)
                type2SignalArrayDict[signature] += window_array
                type2CountArrayDict[signature] += (window_array>0)

            ####################################################
            if (sample in sample2Signature2NumberofMutationsDict) and (signature in sample2Signature2NumberofMutationsDict[sample]):
                if sample in sample2Type2SignalArrayDict:
                    if signature in sample2Type2SignalArrayDict[sample]:
                        sample2Type2SignalArrayDict[sample][signature] += window_array
                        sample2Type2CountArrayDict[sample][signature] += (window_array>0)
                    else:
                        sample2Type2SignalArrayDict[sample][signature] = np.zeros(windowSize)
                        sample2Type2CountArrayDict[sample][signature] = np.zeros(windowSize, dtype=int)
                        sample2Type2SignalArrayDict[sample][signature] += window_array
                        sample2Type2CountArrayDict[sample][signature] += (window_array > 0)

                else:
                    sample2Type2SignalArrayDict[sample] = {}
                    sample2Type2CountArrayDict[sample] = {}
                    sample2Type2SignalArrayDict[sample][signature] = np.zeros(windowSize)
                    sample2Type2CountArrayDict[sample][signature] = np.zeros(windowSize, dtype=int)
                    sample2Type2SignalArrayDict[sample][signature] += window_array
                    sample2Type2CountArrayDict[sample][signature] += (window_array > 0)
            ####################################################
    ################# Signatures ends #########################

    ######################################################################
    if type in type2SignalArrayDict:
        type2SignalArrayDict[type] += window_array
        type2CountArrayDict[type] += (window_array > 0)
    else:
        type2SignalArrayDict[type] = np.zeros(windowSize)
        type2CountArrayDict[type] = np.zeros(windowSize, dtype=int)
        type2SignalArrayDict[type] += window_array
        type2CountArrayDict[type] += (window_array > 0)


    if (sample in sample2NumberofMutationsDict):
        if sample in sample2Type2SignalArrayDict:
            if type in sample2Type2SignalArrayDict[sample]:
                sample2Type2SignalArrayDict[sample][type] += window_array
                sample2Type2CountArrayDict[sample][type] += (window_array > 0)
            else:
                sample2Type2SignalArrayDict[sample][type] = np.zeros(windowSize)
                sample2Type2CountArrayDict[sample][type] = np.zeros(windowSize, dtype=int)
                sample2Type2SignalArrayDict[sample][type] += window_array
                sample2Type2CountArrayDict[sample][type] += (window_array > 0)
        else:
            sample2Type2SignalArrayDict[sample] = {}
            sample2Type2CountArrayDict[sample] = {}
            sample2Type2SignalArrayDict[sample][type] = np.zeros(windowSize)
            sample2Type2CountArrayDict[sample][type] = np.zeros(windowSize, dtype=int)
            sample2Type2SignalArrayDict[sample][type] += window_array
            sample2Type2CountArrayDict[sample][type] += (window_array > 0)
    ######################################################################

########################################################################################


########################################################################################
def computeAverageNucleosomeOccupancyArray(plusorMinus,signalArray,countArray):
    windowSize=plusorMinus*2+1
    averageArray =  np.zeros(windowSize)

    np.seterr(divide='ignore', invalid='ignore')
    if (np.any(countArray)):
        averageArray = np.divide(signalArray,countArray)
    np.seterr(divide='raise', invalid='ignore')

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
        library_file_memo):


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

##################################################################
def init(l):
    global lock
    lock = l
##################################################################


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

############################################################
def addPyramidineStrandColumnToDF(inputList):
    chrShort= inputList[0]
    chr_based_indels_df = inputList[1]
    genome = inputList[2]

    chr_based_indels_df = chr_based_indels_df.apply(addPyramidineStrandColumn,genome=genome,axis=1)

    return chr_based_indels_df
############################################################

############################################################
def allPyrimidine(ref):
    if (('A' in ref) or ('G' in ref)):
        return False
    else:
        return True
############################################################


############################################################
def allPurine(ref):
    if (('T' in ref) or ('C' in ref)):
        return False
    else:
        return True
############################################################

##################################################################
def addPyramidineStrandColumn(mutation_row, genome):
    chromosomeShort = mutation_row[CHROM]
    start = mutation_row[START]
    end = mutation_row[END]
    ref = mutation_row[REF]

    #Make  zero based
    start = start-1
    end = end-1

    referenceGenomeDNASequence = getNucleotides(chromosomeShort, start, end, genome)

    if (ref == referenceGenomeDNASequence) and (allPyrimidine(ref)):
        mutation_row[PYRAMIDINESTRAND] = +1
    elif (ref == referenceGenomeDNASequence) and (allPurine(ref)):
        mutation_row[PYRAMIDINESTRAND] = -1
    else:
        mutation_row[PYRAMIDINESTRAND] = 0

    return mutation_row
##################################################################


##################################################################
def readIndelsAndWriteChrBasedParallel(genome,outputDir,jobname,indelsFileName):
    sample2NumberofIndelsDict, \
    indelsSignatures2NumberofMutationsDict, \
    sample2IndelsSignatures2NumberofMutationsDict, \
    indels_df = readIndels(indelsFileName)

    # print('For debugging purposes indels_df.columns.values')
    # print(indels_df.columns.values)

    ##############################################################################################
    #Do we have PYRAMIDINESTRAND column? If not add
    if (PYRAMIDINESTRAND not in indels_df.columns.values):
        if (genome == GRCh37):
            genome = twobitreader.TwoBitFile(os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB, UCSCGENOME, HG19_2BIT))
        elif (genome == GRCh38):
            genome = twobitreader.TwoBitFile(os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB, UCSCGENOME, HG38_2BIT))

        ###########################################################
        numofProcesses = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(numofProcesses)

        poolInputList = []
        indels_df_grouped = indels_df.groupby(CHROM)
        for chrom, chr_based_indels_df in indels_df_grouped:
            inputList = []
            inputList.append(chrom)
            inputList.append(chr_based_indels_df)
            inputList.append(genome)
            poolInputList.append(inputList)

        list_of_chrBased_indels_df = pool.map(addPyramidineStrandColumnToDF,poolInputList)

        pool.close()
        pool.join()
        ###########################################################

        indels_df = pd.concat(list_of_chrBased_indels_df, axis=0)

    ##############################################################################################

    indels_df_grouped= indels_df.groupby(CHROM)

    # print('len(indels_df_grouped)')
    # print(len(indels_df_grouped))

    # print('Debug March 27, 2019 indels_df.head()')
    # print(indels_df.iloc[0])
    # print(indels_df.iloc[1])
    # print(indels_df.iloc[2])
    # print(indels_df.iloc[3])
    # print(indels_df.iloc[4])
    # print(indels_df.iloc[5])

    os.makedirs(os.path.join(outputDir,jobname,DATA, CHRBASED),exist_ok=True)

    #########################################################
    ############### Write Unique Chrnames starts ############
    #########################################################
    uniqueChrNames = indels_df[CHROM].unique()
    print('Chromosome names in indels data: %s' %(uniqueChrNames))
    # The unique values are returned as a NumPy array

    # # Write uniqueChrNames to a file
    # filename = ChrNamesInIndelsFilename
    # ChrNamesFile = os.path.join(outputDir,jobname,DATA,filename)
    # np.savetxt(ChrNamesFile,uniqueChrNames, delimiter='\t', fmt='%s')
    #########################################################
    ############### Write Unique Chrnames ends ##############
    #########################################################

    #############################################################################################################################
    ####################################################### Write starts ########################################################
    #############################################################################################################################
    writeDictionaryUnderDataDirectory(sample2NumberofIndelsDict,outputDir,jobname,Sample2NumberofIndelsDictFilename)
    writeDictionaryUnderDataDirectory(indelsSignatures2NumberofMutationsDict,outputDir,jobname,IndelsSignature2NumberofMutationsDictFilename)
    writeDictionaryUnderDataDirectory(sample2IndelsSignatures2NumberofMutationsDict,outputDir,jobname,Sample2IndelsSignature2NumberofMutationsDictFilename)
    #############################################################################################################################
    ####################################################### Write ends ##########################################################
    #############################################################################################################################

    #########################################################
    # l = multiprocessing.Lock()
    numofProcesses = multiprocessing.cpu_count()
    # pool = multiprocessing.Pool(numofProcesses,initializer=init, initargs=(l,))
    pool = multiprocessing.Pool(numofProcesses)
    #########################################################

    print('Number of cores: %d' %(numofProcesses))
    poolInputList =[]

    #Get the filename at the end of the full path
    indelsFileName = os.path.basename(indelsFileName)

    #########################################################
    for chr, chrBased_indels_df in indels_df_grouped:
        inputList = []
        inputList.append(chr)
        inputList.append(indelsFileName)
        inputList.append(chrBased_indels_df)
        inputList.append(outputDir)
        inputList.append(jobname)
        poolInputList.append(inputList)
    #########################################################

    pool.map(writeChrBasedMutationDF,poolInputList)

    #########################################################
    pool.close()
    pool.join()
    #########################################################

##################################################################



##################################################################
def readSubsAndWriteChrBasedParallel(outputDir,jobname,subsWithSignatureBasedProbabilitiesFileName):
    signatures, \
    sample2NumberofSubsDict, \
    subsSignature2NumberofMutationsDict, \
    sample2SubsSignature2NumberofMutationsDict, \
    subs_df = readSubs(subsWithSignatureBasedProbabilitiesFileName)

    print('Substitutions Signatures')
    print(signatures)

    mutation_df_grouped= subs_df.groupby(CHROM)
    os.makedirs(os.path.join(outputDir,jobname,DATA,CHRBASED),exist_ok=True)

    #############################################################################################################################
    ################################################# Write Dictionaries starts #################################################
    #############################################################################################################################
    writeDictionaryUnderDataDirectory(sample2NumberofSubsDict,outputDir,jobname,Sample2NumberofSubsDictFilename)
    writeDictionaryUnderDataDirectory(subsSignature2NumberofMutationsDict,outputDir,jobname,SubsSignature2NumberofMutationsDictFilename)
    writeDictionaryUnderDataDirectory(sample2SubsSignature2NumberofMutationsDict,outputDir,jobname,Sample2SubsSignature2NumberofMutationsDictFilename)
    #############################################################################################################################
    ################################################# Write Dictionaries ends ###################################################
    #############################################################################################################################


    #########################################################
    ############### Write Unique Chrnames starts ############
    #########################################################
    uniqueChrNames = subs_df[CHROM].unique()
    print('Chromosome names in single point mutations data: %s' %(uniqueChrNames))

    # ChrNamesFile = os.path.join(outputDir,jobname,DATA,ChrNamesInSPMsFilename)
    # np.savetxt(ChrNamesFile,uniqueChrNames, delimiter='\t', fmt='%s')
    #########################################################
    ############### Write Unique Chrnames ends ##############
    #########################################################

    #########################################################
    # l = multiprocessing.Lock()
    numofProcesses = multiprocessing.cpu_count()
    # pool = multiprocessing.Pool(numofProcesses,initializer=init, initargs=(l,))
    pool = multiprocessing.Pool(numofProcesses)
    #########################################################

    poolInputList =[]

    # Get the filename at the end of the full path
    mutationsWithSignatureBasedProbabilitiesFileName = os.path.basename(subsWithSignatureBasedProbabilitiesFileName)

    #########################################################
    for chr, chrBased_mutation_df in mutation_df_grouped:

        inputList = []
        inputList.append(chr)
        inputList.append(mutationsWithSignatureBasedProbabilitiesFileName)
        inputList.append(chrBased_mutation_df)
        inputList.append(outputDir)
        inputList.append(jobname)
        poolInputList.append(inputList)
    #########################################################

    pool.map(writeChrBasedMutationDF,poolInputList)

    #########################################################
    pool.close()
    pool.join()
    #########################################################

##################################################################


##################################################################
#Please note that we add Count and MMR colum in this method
def readIndels(allIndelsFileName):
    # Columns: ['Sample', 'Chrom', 'Start', 'End', 'Ref','Alt', 'Type', 'Length', 'Category']
    # Sample  Chrom      Start   End     Ref    Alt Type    Length  Category
    allIndelsFilePath = os.path.join(allIndelsFileName)

    #################################################
    #First read only first row
    indels_df = pd.read_table(allIndelsFilePath, sep="\t", comment='#', dtype={SAMPLE: str, CHROM: str},nrows=1)
    columnNamesList = list(indels_df.columns.values)

    contextIndex = columnNamesList.index(CONTEXT)

    # We assume that after the column named 'Context' there are the signature columns in tab separated way.
    signatures = columnNamesList[(contextIndex + 1):]
    print('Indels Signatures')
    print(signatures)
    #################################################

    # indels_df = pd.read_table(allIndelsFilePath,sep="\t",dtype={'Sample': str, 'Chromosome': str}, header=0)
    indels_df = pd.read_table(allIndelsFilePath,sep="\t", header=0)

    indels_df[SAMPLE] = indels_df[SAMPLE].astype(str)
    indels_df[CHROM] = indels_df[CHROM].astype(str)
    indels_df[START] = indels_df[START].astype(int)
    indels_df[END] = indels_df[END].astype(int)

    ###########################################################################
    listofSamples = indels_df[SAMPLE].unique()

    sample2NumberofIndelsDict = {}
    for sample in listofSamples:
        numberofIndels =  len(indels_df[indels_df[SAMPLE] == sample])
        if (numberofIndels>=INDEL_NUMBER_OF_MINIMUM_REQUIRED_MUTATIONS):
            sample2NumberofIndelsDict[sample] = numberofIndels
    ###########################################################################

    print('Number of samples in indels file: %d' %(len(indels_df[SAMPLE].unique())))
    print('Number of indels: %d' %(indels_df.shape[0]))
    # print('Indels file shape')
    # print(indels_df.shape)

    # drop columns 'Reference','Mutation','Type','Length','Category'
    # Do not drop Length column. If length >= 3bp indel type is Microhomology otherwise < 3bp indel type is repeat-med indel
    # indels_df.drop([REF,ALT,TYPE,CATEGORY], axis=1, inplace=True, errors='ignore')
    indels_df.drop([TYPE, CATEGORY], axis=1, inplace=True, errors='ignore')

    #Add a new column called Count
    grouped_mutation_df_by_sample = indels_df.groupby(SAMPLE)
    indels_df[COUNT] = indels_df.groupby(SAMPLE)[SAMPLE].transform('count')

    # print('Number of samples in indels file:%s' % len(grouped_mutation_df_by_sample))

    # Add a new column called MMR (Mis Match Repair)
    if 'Count' in indels_df.columns:
        indels_df[MMR] = np.where(indels_df[COUNT] >= SAMPLE_MMR_DEFICIENT_THRESHOLD ,DEFICIENT,PROFICIENT)
    #################################################################################################
    indelsSignature2NumberofMutationsDict = {}
    sample2IndelsSignature2NumberofMutationsDict = {}

    for signature in signatures:
        signaturebased_df = indels_df[indels_df[signature] >= INDEL_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD]
        numberofIndels = len(signaturebased_df)
        if (numberofIndels>INDEL_NUMBER_OF_MINIMUM_REQUIRED_MUTATIONS):
            indelsSignature2NumberofMutationsDict[signature] = numberofIndels

    for sample in sample2NumberofIndelsDict:
        for signature in indelsSignature2NumberofMutationsDict:
            numberofMutations = len(indels_df[ ((indels_df[SAMPLE]==sample) & (indels_df[signature]>= INDEL_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD)) ])
            if (numberofMutations> INDEL_NUMBER_OF_MINIMUM_REQUIRED_MUTATIONS):
                if sample in sample2IndelsSignature2NumberofMutationsDict:
                    sample2IndelsSignature2NumberofMutationsDict[sample][signature] = numberofMutations
                else:
                    sample2IndelsSignature2NumberofMutationsDict[sample]={}
                    sample2IndelsSignature2NumberofMutationsDict[sample][signature] = numberofMutations
    #################################################################################################

    return sample2NumberofIndelsDict,indelsSignature2NumberofMutationsDict,sample2IndelsSignature2NumberofMutationsDict,indels_df
##################################################################


######################################################################
#BED files: end is not inclusive
#BED files: score - A score between 0 and 1000. See track lines, below, for ways to configure the display style of scored data.
def updateChrBasedSignalArray(data_row,chrBasedSignalArray):
    chrBasedSignalArray[data_row[start]:data_row[end]] += data_row[signal]
######################################################################


######################################################################
def fillNumpyArray(start,end,signal,chrBasedSignalArray):
    chrBasedSignalArray[start:end]+=signal
######################################################################


##################################################################
#Make signal from bed files and signal from wig files have the same type np.float32
#No outlier elimination in bed files
def readFileInBEDFormat(file_with_path):
    file_df=None

    if os.path.exists(file_with_path):
        file_df = pd.read_table(file_with_path, nrows=1)  # 2.25 GB
        ncols=file_df.shape[1]

        if (ncols<=3):
            print('There is no enough columns in this bed file')
        elif (ncols==4):
            print('SigProfilerTopogarphy assumes that score column is in the 4th column of this bed file and there is no header')
            file_df=pd.read_table(file_with_path,header=None, usecols=[0, 1, 2, 3],names = [chrom,start,end,signal],dtype={0: 'category', 1: np.int32, 2: np.int32, 3: np.float32})
        elif (ncols>=5):
            print('SigProfilerTopogarphy assumes that score column is in the 5th column of this bed file and there is no header')
            file_df=pd.read_table(file_with_path,header=None, usecols=[0, 1, 2, 4], names = [chrom,start,end,signal], dtype={0: 'category', 1: np.int32, 2: np.int32, 4: np.float32})

        # print("file_df.dtypes")
        # print(file_df.dtypes)

        # print('file_df.columns.values')
        # print(file_df.columns.values)

    return file_df
##################################################################


##################################################################
#Make signal from narrowPeak files type np.float32
#No outlier elimination in bed files
def readFileInNarrowPeakFormat(file_with_path):
#ENCODE narrowPeak: Narrow (or Point-Source) Peaks format
#This format is used to provide called peaks of signal enrichment based on pooled, normalized (interpreted) data. It is a BED6+4 format

#Discrad any line that starts with
# track type=narrowPeak visibility=3 db=hg19 name="nPk" description="ENCODE narrowPeak Example"
# browser position chr1:9356000-9365000

    file_df = None

    if os.path.exists(file_with_path):
        file_df = pd.read_table(file_with_path, nrows=1)  # 2.25 GB
        ncols = file_df.shape[1]

        if (ncols == 10):
            print('SigProfilerTopogarphy assumes that signal column is in the 7th column of this bed file and there is no header')
            file_df = pd.read_table(file_with_path, header=None, usecols=[0,1,2,3,4,5,6],
                                    names=[chrom, start, end, name,score,strand,signal],
                                    dtype={0: 'category', 1: np.int32, 2: np.int32, 3:str,4:np.int32, 5:'category', 6:np.float32})

            file_df.drop([3,4,5], inplace=True, axis=1)
    return file_df
##################################################################


##################################################################
def readBED(bedFilename):
    columns = ['chr', 'start', 'end', 'column4','column5','column6','column7','column8','column9']

    #BED files are seperated by tab
    bed_df = pd.read_table(bedFilename, sep='\t', comment='#', header=None, names=columns)

    #BED files are end exclusive by default.
    #Make end inclusive
    bed_df['end'] =bed_df['end']-1

    # print('debug starts')
    # print('############## bed_df starts ##############')
    # print(bed_df.shape)
    # print(bed_df.head())
    # print(bed_df.info())
    # print('############## bed_df ends ##############')
    # print('debug ends')

    return bed_df
##################################################################


#todo which method to use?
#readRepliSeqSignal uses * as separator
#readWaveletSmoothedRepliSeqSignal

##################################################################
# WaveletSmoothedSignal Filename: 'GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.wig'
# SumSignal: GSM923442_hg19_wgEncodeUwRepliSeqMcf7SumSignalRep1.wig
def readRepliSeqSignal(repliseqDataFilename):
    #TODO is there a better/more clean way to handle this?
    #This sep='*' is used for reading line by line without different number of columns error for different lines
    repliSeq_unprocessed_df = pd.read_table(repliseqDataFilename, sep='*', comment='#', header=None)
    # print('debug starts')
    # print('############## wavelet_unprocessed_df ##############')
    # print(repliSeq_unprocessed_df.shape)
    # print(repliSeq_unprocessed_df.head())
    # print('############## wavelet_unprocessed_df ##############')
    # print('debug ends')

    return repliSeq_unprocessed_df
##################################################################



##################################################################
# Original Filename: GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.bigWig
# WaveletSmoothedSignal Filename: 'GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.wig'
def readWaveletSmoothedSignalReplicationTime(repliseqDataFilename):
    replication_time_wavelet_signal_unprocessed_df = pd.read_table(repliseqDataFilename, sep="\t", comment='#', header=None)

    return replication_time_wavelet_signal_unprocessed_df
##################################################################


##################################################################
def processSmoothedWaveletSignal(wavelet_unprocessed_df):
    #Read the file and generate chr start end signal wavelet_smoothed_signal
    #Then sort the data w.r.t. signal in descending order
    #Divide the data into 10 equal deciles
    #Return 10 deciles: the first decile is the earliest one and the tenth decile is the latest one

    columns = ['chr','start','end','signal']

    #Create an empty dataframe
    wavelet_processed_df = pd.DataFrame(columns=columns)

    # print('############### empty wavelet_processed_df starts #################')
    # print(wavelet_processed_df.shape)
    # print(wavelet_processed_df.head())
    # print('############### empty wavelet_processed_df ends #################')

    #dummy initialization
    i = 0
    chrom = 'chr1'
    start = 0
    step = 1000

    rows_list = []

    for row in wavelet_unprocessed_df.itertuples(index=True, name='Pandas'):
        # row's type is <class 'pandas.core.frame.Pandas'>
        # row[0] is the index
        # row[1] is fixedStep chrom=chr1 start=24500 step=1000 span=1000 or 57.4679
        if (row[1].startswith('fixedStep')):
            #This is the information line
            chrom = row[1].split()[1].split('=')[1]
            start = int(row[1].split()[2].split('=')[1])
            step = int(row[1].split()[3].split('=')[1])
        else:
            signal = float(row[1])
            chr = chrom
            start = start
            end = start + step-1
            dict = {'chr':chr, 'start':start, 'end':end, 'signal':signal}
            rows_list.append(dict)
            start += step

    # print('Number of intervals to be inserted in wavelet_processed_df: %d' %len(rows_list))

    #rows_list contain the list of row where each row is a dictionary
    wavelet_processed_df = pd.DataFrame(rows_list, columns=['chr','start','end','signal'])

    # print('############### wavelet_processed_df is filled starts #################')
    # print(wavelet_processed_df.shape)
    # print(wavelet_processed_df.head())
    # print('############### wavelet_processed_df is filled ends #################')

    return wavelet_processed_df
##################################################################

##################################################################
def processSumSignal(sum_signal_unprocessed_df):
    columns = ['chr', 'start', 'end','signal']

    #Create an empty dataframe
    sum_signal_processed_df = pd.DataFrame(columns=columns)

    # print('############### empty sum_signal_processed_df starts #################')
    # print(sum_signal_processed_df.shape)
    # print(sum_signal_processed_df.head())
    # print('############### empty sum_signal_processed_df ends #################')

    #dummy initialization
    i = 0
    chrom = 'chr1'
    start = 0
    step = 1000

    rows_list = []

    for row in sum_signal_unprocessed_df.itertuples(index=True, name='Pandas'):
        # row's type is <class 'pandas.core.frame.Pandas'>
        # row[0] is the index
        # row[1] is fixedStep chrom=chr1 start=24500 step=1000 span=1000 or 57.4679
        if (row[1].startswith('variableStep')):
            # e.g. row[1] variableStep chrom=chr1 span=1000
            #This is the information line
            chrom = row[1].split()[1].split('=')[1]
            step = int(row[1].split()[2].split('=')[1])
        else:
            # e.g. row[1] '24500\t35'
            start= int(row[1].split('\t')[0])
            signal = float(row[1].split('\t')[1])
            chr = chrom
            end = start + step-1

            dict = {'chr':chr, 'start':start, 'end':end, 'signal':signal}
            rows_list.append(dict)

    # print('Number of intervals to be inserted in wavelet_processed_df: %d' %len(rows_list))

    #rows_list contain the list of row where each row is a dictionary
    sum_signal_processed_df = pd.DataFrame(rows_list, columns=['chr','start','end','signal'])

    # print('############### wavelet_processed_df is filled starts #################')
    # print(sum_signal_processed_df.shape)
    # print(sum_signal_processed_df.head())
    # print('############### wavelet_processed_df is filled ends #################')

    return sum_signal_processed_df
##################################################################


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
                        signature2PropertiesListDict):

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
    for signature in signature2PropertiesListDict:
        if (mutation_row[signature] >= float(signature2PropertiesListDict[signature][0])):
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
        for signature in signature2PropertiesListDict:
            if (mutation_row[signature] >= float(signature2PropertiesListDict[signature][0])):
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
        for signature in signature2PropertiesListDict:
            if (mutation_row[signature] >= float(signature2PropertiesListDict[signature][0])):
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
        for signature in signature2PropertiesListDict:
            if (mutation_row[signature] >= float(signature2PropertiesListDict[signature][0])):
                if (type2Sample2Strand2CountDict is not None):
                    updateTypeBasedDictionary(type2Sample2Strand2CountDict,signature,mutationSample,strand)
    #################################################################################################


########################################################################

########################################################################
#legacy code
def updateDictionaries(mutation_row,
                        mutationType,
                        mutationSample,
                        type2Strand2CountDict,
                        sample2Type2Strand2CountDict,
                        type2Sample2Strand2CountDict,
                        signature2MutationType2Strand2CountDict,
                        strand,
                        signature2NumberofMutationsDict,
                        mutationProbability):

    #################################################################################################
    # Update1: update mutationType in type2Strand2CountDict
    if ((type2Strand2CountDict is not None) and (mutationType is not None)):
        updateDictionary(type2Strand2CountDict,mutationType,strand)
    #################################################################################################

    #################################################################################################
    # Update2: signature in type2Strand2CountDict
    for signature in signature2NumberofMutationsDict:
        if (mutation_row[signature] >= mutationProbability):
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
    if ((signature2MutationType2Strand2CountDict is not None) and (mutationType is not None)):
        for signature in signature2NumberofMutationsDict:
            if (mutation_row[signature] >= mutationProbability):
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
    if (sample2Type2Strand2CountDict is not None) and (mutationType is not None):
        updateSampleBasedDictionary(sample2Type2Strand2CountDict,mutationType,mutationSample,strand)
    #################################################################################################

    #################################################################################################
    # Update5: sample and signature in sample2Type2Strand2CountDict
    for signature in signature2NumberofMutationsDict:
        if (mutation_row[signature] >= mutationProbability):
            updateSampleBasedDictionary(sample2Type2Strand2CountDict, signature, mutationSample, strand)
    #################################################################################################

    #################################################################################################
    #Update6 sample and mutationType type2Sample2TranscriptionStrand2CountDict
    if (type2Sample2Strand2CountDict is not None) and (mutationType is not None):
        updateTypeBasedDictionary(type2Sample2Strand2CountDict,mutationType,mutationSample,strand)
    #################################################################################################

    #################################################################################################
    #Update7: sample and signature in type2Sample2TranscriptionStrand2CountDict
    for signature in signature2NumberofMutationsDict:
        if (mutation_row[signature] >= mutationProbability):
            if (type2Sample2Strand2CountDict is not None):
                updateTypeBasedDictionary(type2Sample2Strand2CountDict,signature,mutationSample,strand)
    #################################################################################################

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
                                        if strand in accumulatedAllChromosomes_SimNum2Sample2Type2Strand2CountDict[simNum][sample][type]:
                                            accumulatedAllChromosomes_SimNum2Sample2Type2Strand2CountDict[simNum][sample][type][strand] += count
                                        else:
                                            accumulatedAllChromosomes_SimNum2Sample2Type2Strand2CountDict[simNum][sample][type][strand] = count
                                else:
                                    accumulatedAllChromosomes_SimNum2Sample2Type2Strand2CountDict[simNum][sample][type] = strand2CountDict
                        else:
                            accumulatedAllChromosomes_SimNum2Sample2Type2Strand2CountDict[simNum][sample] = type2Strand2CountDict
                else:
                    accumulatedAllChromosomes_SimNum2Sample2Type2Strand2CountDict[simNum] = chrBasedSample2Type2Strand2CountDict


        # Accumulate3 chrBasedType2Sample2Strand2CountDict in accumulatedAllChromosomesType2Sample2Strand2CountDict
        if (chrBased_SimNum2Type2Sample2Strand2CountDict is not None):
            for simNum, chrBasedType2Sample2Strand2CountDict in chrBased_SimNum2Type2Sample2Strand2CountDict.items():
                if simNum in accumulatedAllChromosomes_SimNum2Type2Sample2Strand2CountDict:
                    for type, sample2Strand2CountDict in chrBasedType2Sample2Strand2CountDict.items():
                        if type in accumulatedAllChromosomes_SimNum2Type2Sample2Strand2CountDict[simNum]:
                            for sample, strand2CountDict in sample2Strand2CountDict.items():
                                if sample in accumulatedAllChromosomes_SimNum2Type2Sample2Strand2CountDict[simNum][type]:
                                    for strand, count in strand2CountDict.items():
                                        if strand in accumulatedAllChromosomes_SimNum2Type2Sample2Strand2CountDict[simNum][type][sample]:
                                            accumulatedAllChromosomes_SimNum2Type2Sample2Strand2CountDict[simNum][type][sample][strand] += count
                                        else:
                                            accumulatedAllChromosomes_SimNum2Type2Sample2Strand2CountDict[simNum][type][sample][strand] = count
                                else:
                                    accumulatedAllChromosomes_SimNum2Type2Sample2Strand2CountDict[simNum][type][sample] = strand2CountDict
                        else:
                            accumulatedAllChromosomes_SimNum2Type2Sample2Strand2CountDict[simNum][type] = sample2Strand2CountDict
                else:
                    accumulatedAllChromosomes_SimNum2Type2Sample2Strand2CountDict[simNum] = chrBasedType2Sample2Strand2CountDict


        # Accumulate4 chrBasedSignature2MutationType2Strand2CountDict in accumulatedAllChromosomesSignature2MutationType2Strand2CountDict
        if (chrBased_SimNum2Signature2MutationType2Strand2CountDict is not None):
            for simNum, chrBasedSignature2MutationType2Strand2CountDict in chrBased_SimNum2Signature2MutationType2Strand2CountDict.items():
                if simNum in accumulatedAllChromosomes_SimNum2Signature2MutationType2Strand2CountDict:
                    for signature, mutationType2Strand2CountDict in chrBasedSignature2MutationType2Strand2CountDict.items():
                        if signature in accumulatedAllChromosomes_SimNum2Signature2MutationType2Strand2CountDict[simNum]:
                            for mutationType, strand2CountDict in mutationType2Strand2CountDict.items():
                                if mutationType in accumulatedAllChromosomes_SimNum2Signature2MutationType2Strand2CountDict[simNum][signature]:
                                    for strand, count in strand2CountDict.items():
                                        if strand in accumulatedAllChromosomes_SimNum2Signature2MutationType2Strand2CountDict[simNum][signature][mutationType]:
                                            accumulatedAllChromosomes_SimNum2Signature2MutationType2Strand2CountDict[simNum][signature][mutationType][strand] += count
                                        else:
                                            accumulatedAllChromosomes_SimNum2Signature2MutationType2Strand2CountDict[simNum][signature][mutationType][strand] = count
                                else:
                                    accumulatedAllChromosomes_SimNum2Signature2MutationType2Strand2CountDict[simNum][signature][mutationType] = strand2CountDict
                        else:
                            accumulatedAllChromosomes_SimNum2Signature2MutationType2Strand2CountDict[simNum][signature] = mutationType2Strand2CountDict
                else:
                    accumulatedAllChromosomes_SimNum2Signature2MutationType2Strand2CountDict[simNum] = chrBasedSignature2MutationType2Strand2CountDict

########################################################################

########################################################################
def accumulate(listofTuples,
               accumulatedAllChromosomesType2Strand2CountDict,
               accumulatedAllChromosomesSample2Type2Strand2CountDict,
               accumulatedAllChromosomesType2Sample2Strand2CountDict,
               accumulatedAllChromosomesSignature2MutationType2Strand2CountDict):

    for mytuple in listofTuples:
        chrBasedType2Strand2CountDict = mytuple[0]
        chrBasedSample2Type2Strand2CountDict = mytuple[1]
        chrBasedType2Sample2Strand2CountDict = mytuple[2]
        chrBasedSignature2MutationType2Strand2CountDict = mytuple[3]

        #Accumulate1 chrBasedType2Strand2CountDict in accumulatedAllChromosomesType2Strand2CountDict
        if (chrBasedType2Strand2CountDict is not None):
            for type,strand2CountDict in chrBasedType2Strand2CountDict.items():
                if type in accumulatedAllChromosomesType2Strand2CountDict:
                    for strand, count in strand2CountDict.items():
                        if strand in accumulatedAllChromosomesType2Strand2CountDict[type]:
                            accumulatedAllChromosomesType2Strand2CountDict[type][strand] += count
                        else:
                            accumulatedAllChromosomesType2Strand2CountDict[type][strand] = count
                else:
                    accumulatedAllChromosomesType2Strand2CountDict[type] = strand2CountDict

        #Accumulate2 chrBasedSample2Type2ReplicationStrand2CountDict in accumulatedAllChromosomesSample2Type2ReplicationStrand2CountDict
        if (chrBasedSample2Type2Strand2CountDict is not None):
            for sample, type2Strand2CountDict in chrBasedSample2Type2Strand2CountDict.items():
                if sample in accumulatedAllChromosomesSample2Type2Strand2CountDict:
                    for type,strand2CountDict in chrBasedSample2Type2Strand2CountDict[sample].items():
                        if type in accumulatedAllChromosomesSample2Type2Strand2CountDict[sample]:
                            for strand, count in chrBasedSample2Type2Strand2CountDict[sample][type].items():
                                if strand in accumulatedAllChromosomesSample2Type2Strand2CountDict[sample][type]:
                                    accumulatedAllChromosomesSample2Type2Strand2CountDict[sample][type][strand] += count
                                else:
                                    accumulatedAllChromosomesSample2Type2Strand2CountDict[sample][type][strand] = count
                        else:
                            accumulatedAllChromosomesSample2Type2Strand2CountDict[sample][type] = strand2CountDict
                else:
                    accumulatedAllChromosomesSample2Type2Strand2CountDict[sample] = type2Strand2CountDict


        #Accumulate3 chrBasedType2Sample2Strand2CountDict in accumulatedAllChromosomesType2Sample2Strand2CountDict
        if (chrBasedType2Sample2Strand2CountDict is not None):
            for type, sample2Strand2CountDict in chrBasedType2Sample2Strand2CountDict.items():
                if type in accumulatedAllChromosomesType2Sample2Strand2CountDict:
                    for sample, strand2CountDict in chrBasedType2Sample2Strand2CountDict[type].items():
                        if sample in accumulatedAllChromosomesType2Sample2Strand2CountDict[type]:
                            for strand,count in chrBasedType2Sample2Strand2CountDict[type][sample].items():
                                if strand in accumulatedAllChromosomesType2Sample2Strand2CountDict[type][sample]:
                                    accumulatedAllChromosomesType2Sample2Strand2CountDict[type][sample][strand] += count
                                else:
                                    accumulatedAllChromosomesType2Sample2Strand2CountDict[type][sample][strand] = count
                        else:
                            accumulatedAllChromosomesType2Sample2Strand2CountDict[type][sample] = strand2CountDict
                else:
                    accumulatedAllChromosomesType2Sample2Strand2CountDict[type] = sample2Strand2CountDict


        #Accumulate4 chrBasedSignature2MutationType2Strand2CountDict in accumulatedAllChromosomesSignature2MutationType2Strand2CountDict
        if (chrBasedSignature2MutationType2Strand2CountDict is not None):
            for signature,mutationType2Strand2CountDict in chrBasedSignature2MutationType2Strand2CountDict.items():
                if signature in accumulatedAllChromosomesSignature2MutationType2Strand2CountDict:
                    for mutationType, strand2CountDict in chrBasedSignature2MutationType2Strand2CountDict[signature].items():
                        if mutationType in accumulatedAllChromosomesSignature2MutationType2Strand2CountDict[signature]:
                            for strand, count in chrBasedSignature2MutationType2Strand2CountDict[signature][mutationType].items():
                                if strand in accumulatedAllChromosomesSignature2MutationType2Strand2CountDict[signature][mutationType]:
                                    accumulatedAllChromosomesSignature2MutationType2Strand2CountDict[signature][mutationType][strand] += count
                                else:
                                    accumulatedAllChromosomesSignature2MutationType2Strand2CountDict[signature][mutationType][strand] = count

                        else:
                            accumulatedAllChromosomesSignature2MutationType2Strand2CountDict[signature][mutationType] = strand2CountDict
                else:
                    accumulatedAllChromosomesSignature2MutationType2Strand2CountDict[signature] = mutationType2Strand2CountDict

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

########################################################################
def createFiles(outputDir,jobname,filename):
    os.makedirs(os.path.join(outputDir, jobname, DATA), exist_ok=True)
    filePath = os.path.join(outputDir, jobname, DATA, filename)
    open(filePath,'w+')
########################################################################

########################################################################
#To write samples with signatures with at least 10K eligible mutations
def appendDictionaryUnderDataDirectory(dictionary,outputDir,jobname,filename):
    os.makedirs(os.path.join(outputDir,jobname,DATA), exist_ok=True)
    filePath = os.path.join(outputDir,jobname,DATA,filename)

    with open(filePath, 'a') as file:
        file.write(json.dumps(dictionary))
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
def writeList2File(list,filePath):
    with open(filePath, "w") as f:
        for item in list:
            f.write("%s\n" % item)
########################################################################

########################################################################
def append2File(item,filePath):
    if (os.path.exists(filePath)):
        with open(filePath, "a+") as f:
            f.write("%s\n" % item)
########################################################################



########################################################################
#Will be depreceated
def readSignatureList(filePath):
    list = []
    if (os.path.exists(filePath)):
        with open(filePath, "r") as f:
            for line in f:
                list.append(line.strip())
        return list
    else:
        return list
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
def readChromSizes(genomeAssemby):

    if genomeAssemby==HG19:
        chromSizesFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,UCSCGENOME, HG19_CHROM_SIZES)
    elif genomeAssemby==HG38:
        chromSizesFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,UCSCGENOME, HG38_CHROM_SIZES)

    chrom2size_df = pd.read_table(chromSizesFilePath, sep='\t', header = None, names = {'chrom','size'})

    chrom2size_dict = dict(zip(chrom2size_df['chrom'], chrom2size_df['size']))

    return chrom2size_dict
########################################################################