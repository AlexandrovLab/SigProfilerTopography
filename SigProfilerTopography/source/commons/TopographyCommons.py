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

from multiprocessing import Lock

#############################################################
current_abs_path = os.path.dirname(os.path.realpath(__file__))
#############################################################

LEADING= 'Leading'
LAGGING = 'Lagging'

UNTRANSCRIBED_STRAND = 'UnTranscribed Strand'
TRANSCRIBED_STRAND = 'Transcribed Strand'

PLUS = '+'
MINUS = '-'

MAXIMUM_CHROMOSOME_LENGTH = 250000000
GIGABYTE_IN_BYTES = 1000000000


#Constraints , Thresholds
ONE_THOUSAND = 1000
TEN_THOUSAND = 10000

# #TODO This will be used
# SUBSTITUTION_SIGNATURE_PROBABILITY_THRESHOLD = round(0.5,2)
# DINUC_SIGNATURE_PROBABILITY_THRESHOLD = round(0.5,2)
# INDEL_SIGNATURE_PROBABILITY_THRESHOLD = round(0.5,2)

MUTATION_SIGNATURE_PROBABILITY_THRESHOLD = round(0.9,2)
MUTATION_SIGNATURE_PROBABILITY_THRESHOLD_STRING = str(MUTATION_SIGNATURE_PROBABILITY_THRESHOLD)

#Global Constants for Nucleosome Occupancy Analysis
plusOrMinus = 2000
windowSize = plusOrMinus*2+1

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

NCBI = 'ncbi'
ENSEMBL = 'ensembl'

ENCODE_NUCLEOSOME_K562_WIG = 'wgEncodeSydhNsomeK562Sig.wig'
GSM923442_HG19_ENCODE_REPLISEQ_MCF7_WAVELET_SIGNAL_WIG = 'GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.wig'
GSM923442_HG19_ENCODE_REPLISEQ_MCF7_VALLEY_BED = 'GSM923442_hg19_wgEncodeUwRepliSeqMcf7ValleysRep1.bed'
GSM923442_HG19_ENCODE_REPLISEQ_MCF7_PEAK_BED = 'GSM923442_hg19_wgEncodeUwRepliSeqMcf7PkRep1.bed'
GRCh37_hg19_NCBIREFSEQCURATED = 'GRCh37_hg19_NCBIRefSeqCurated'
GRCh37_hg19_ENSEMBL = 'GRCh37_transcripts.txt'

HG19_CHROM_SIZES = 'hg19.chrom.sizes.txt'
HG38_CHROM_SIZES = 'hg38.chrom.sizes.txt'

HG19_2BIT = 'hg19.2bit'
HG38_2BIT = 'hg38.2bit'

HG19 = 'hg19'
HG38 = 'hg38'
MM10 = 'mm10'

GRCh37 = 'GRCh37'
GRCh38 = 'GRCh38'

#For each job we will write signatures, chrnames in nucleosomes, single point mutations, and indels
SignatureFilename = 'Signatures.txt'
SamplesWithAtLeast10KMutations2NumberofMutationsDictFilename = 'SamplesWithAtLeast10KMutations2NumberofMutationsDict.txt'
Samples2NumberofIndelsDictFilename = 'Samples2NumberofIndelsDict.txt'

SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilename = 'SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict.txt'
Sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilename = 'Sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict.txt'
ChrNamesInNucleosomesFilename = 'ChrNamesInNucleosomes.txt'
ChrNamesInIndelsFilename = 'ChrNamesInIndels.txt'
ChrNamesInSPMsFilename = 'ChrNamesInSPMs.txt'

AVAILABLENUCLEOSOMEOCCUPANCYFILESNAME = 'AvailableNucleosomeOccupancyFiles.txt'
GRCh37ChromSizesDictFilename = 'hg19ChromSizesDict.txt'
GRCh38ChromSizesDictFilename = 'hg38ChromSizesDict.txt'

ONE_DIRECTORY_UP = '..'

INDELBASED = 'indelbased'
SIGNATUREBASED = 'signaturebased'
SAMPLEBASED = 'samplebased'

SAMPLEBASED_SIGNATUREBASED = 'samplebased_signaturebased'
SAMPLEBASED_AGGREGATEDINDELS = 'samplebased_aggregatedindels'
SAMPLEBASED_AGGREGATEDSUBSTITUTIONS = 'samplebased_aggregatedsubstitutions'

AGGREGATEDINDELS = 'aggregatedindels'
AGGREGATEDSUBSTITUTIONS = 'aggregatedsubstitutions'

SOURCE = 'source'
COMMONS = 'commons'
NUCLEOSOMEOCCUPANCY = 'nucleosome_occupancy'
REPLICATIONTIME = 'replication_time'
PROCESSIVITY = 'processivity'
STRANDBIAS = 'strand_bias'
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

MICROHOMOLOGY = 'Microhomology'
REPEAT = 'Repeat'

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

CHR = 'chr'

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
MUTATION = 'Mutation'
MUTATIONS = 'Mutations'
CONTEXT = 'Context'

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
############################################

############################################
#Nucleosome Occupancy Data  Columns
chrom= 'chrom'
start= 'start'
end= 'end'
signal= 'signal'
############################################

#For double check purposes
OriginalPaperSignature2Sample2ReplicationStrand2CountDict_Filename = 'OriginalPaperSignature2Sample2ReplicationStrand2CountDict.txt'
OriginalPaperSignature2Sample2TranscriptionStrand2CountDict_Filename = 'OriginalPaperSignature2Sample2TranscriptionStrand2CountDict.txt'

MutationType2ReplicationStrand2CountDict_Filename           = 'MutationType2ReplicationStrand2CountDict.txt'
MutationType2Sample2ReplicationStrand2CountDict_Filename    = 'MutationType2Sample2ReplicationStrand2CountDict.txt'
MutationProbability2Signature2ReplicationStrand2CountDict_Filename = 'MutationProbability2Signature2ReplicationStrand2CountDict.txt'
MutationProbability2Signature2Sample2ReplicationStrand2CountDict_Filename = 'MutationProbability2Signature2Sample2ReplicationStrand2CountDict.txt'
Signature2ReplicationWeightedAverageRatioDict_Filename = 'Signature2ReplicationWeightedAverageRatioDict.txt'
Signature2ReplicationStdErrorDict_Filename = 'Signature2ReplicationStdErrorDict.txt'
Signature2ReplicationSumofMutationProbabilitiesDict_Filename = 'Signature2ReplicationSumofMutationProbabilitiesDict.txt'

MutationType2TranscriptionStrand2CountDict_Filename = 'MutationType2TranscriptionStrand2CountDict.txt'
MutationType2Sample2TranscriptionStrand2CountDict_Filename = 'MutationType2Sample2TranscriptionStrand2CountDict.txt'
MutationProbability2Signature2TranscriptionStrand2CountDict_Filename = 'MutationProbability2Signature2TranscriptionStrand2CountDict.txt'
MutationProbability2Signature2Sample2TranscriptionStrand2CountDict_Filename = 'MutationProbability2Signature2Sample2TranscriptionStrand2CountDict.txt'
Signature2TranscriptionWeightedAverageRatioDict_Filename = 'Signature2TranscriptionWeightedAverageRatioDict.txt'
Signature2TranscriptionStdErrorDict_Filename = 'Signature2TranscriptionStdErrorDict.txt'
Signature2TranscriptionSumofMutationProbabilitiesDict_Filename = 'Signature2TranscriptionSumofMutationProbabilitiesDict.txt'

DATA_Folder = 'data'

BONFERRONI_CORRECTION = 'BONFERRONI_CORRECTION'
FDR_BH_CORRECTION = 'FDR_BH_CORRECTION'

USING_POISSON_DISTRIBUTION = 'USING_POISSON_DISTRIBUTION'
USING_NULL_DISTRIBUTION = 'USING_NULL_DISTRIBUTION'
USING_GAUSSIAN_KDE = 'USING_GAUSSIAN_KDE'


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
def readTrancriptsENSEMBL():
    transcriptsFilenamePath = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,TRANSCRIPTS,GRCh37_hg19_ENSEMBL)

    if (os.path.exists(transcriptsFilenamePath)):
        GRCh37_hg19_ensembl_df = pd.read_table(transcriptsFilenamePath, header=0,sep="\t")
        print('debug GRCh37_hg19_ensembl_df starts')
        print('before columns renaming GRCh37_hg19_ensembl_df.columns.values')
        print(GRCh37_hg19_ensembl_df.columns.values)
        GRCh37_hg19_ensembl_df.rename(columns={'Chromosome/scaffold name': 'chrom', 'Strand': 'strand', 'Transcript start (bp)':'txStart', 'Transcript end (bp)':'txEnd'}, inplace=True)
        print('after columns renaming GRCh37_hg19_ensembl_df.columns.values')
        print(GRCh37_hg19_ensembl_df.columns.values)
        # print(GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df.shape)
        # print(GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df.head())
        # print(GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df.tail())
        # print(GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df.info())
        print('debug GRCh37_hg19_ensembl_df ends')
        #Change the column name     Chromosome/scaffold name    -->     chrom
        #                           Strand                      -->     strand
        #                           Transcript start (bp)       -->     txStart
        #                           Transcript end (bp)         -->     txEnd
        return GRCh37_hg19_ensembl_df
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
def readChrBasedIndelsDF(outputDir,jobname,chrLong,filename):
    filename = '%s_%s' %(chrLong,filename)

    chrBasedMutationDFFilePath = os.path.join(outputDir,jobname,DATA,CHRBASED,filename)

    #############################################
    if os.path.exists(chrBasedMutationDFFilePath):
        chrBased_mutation_df = pd.read_table(chrBasedMutationDFFilePath,sep="\t", comment='#')
        return chrBased_mutation_df
    else:
        return None
##################################################################


##################################################################
def readChrBasedMutationDF(outputDir,jobname,chrLong,filename):

    filename = chrLong + '_' + filename

    chrBasedMutationDFFilePath = os.path.join(outputDir,jobname,DATA,CHRBASED,filename)

    ##########################################################################################
    if (os.path.exists(chrBasedMutationDFFilePath)):

        #############################################
        only_header_chrBased_mutation_df = pd.read_table(chrBasedMutationDFFilePath, sep="\t", comment='#', nrows=1)
        columnNamesList = list(only_header_chrBased_mutation_df.columns.values)

        contextIndex = columnNamesList.index(CONTEXT)

        # We assume that after the column named 'Context' there are the signature columns in tab separated way.
        signatures = columnNamesList[(contextIndex + 1):]
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
        mydtypes[END] = np.int32
        mydtypes[PYRAMIDINESTRAND] = np.int8
        mydtypes[MUTATION] = str
        mydtypes[CONTEXT] = str
        #################################################
    ##########################################################################################

    #############################################
    if os.path.exists(chrBasedMutationDFFilePath):
        chrBased_mutation_df = pd.read_table(chrBasedMutationDFFilePath,sep="\t", comment='#',dtype=mydtypes)
        return chrBased_mutation_df
    else:
        return None
##################################################################



##################################################################
#mutationsWithSignatureBasedProbabilitiesFileName 'breast_cancer_mutation_probabilities_final.txt'
def readMutations(mutationsWithSignatureBasedProbabilitiesFileName):

    # Read the file w.r.t. the current folder
    mutationFilePath = os.path.join(mutationsWithSignatureBasedProbabilitiesFileName)
    # mutationFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, INPUT, jobname,mutationsWithSignatureBasedProbabilitiesFileName)

    #################################################
    #First read only first row
    mutation_df = pd.read_table(mutationFilePath, sep="\t", comment='#', dtype={SAMPLE: str, CHROM: str},nrows=1)
    columnNamesList = list(mutation_df.columns.values)

    contextIndex = columnNamesList.index(CONTEXT)

    # We assume that after the column named 'Context' there are the signature columns in tab separated way.
    signatures = columnNamesList[(contextIndex + 1):]
    #################################################

    signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = {}

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
    mydtypes[END] = np.int32
    mydtypes[PYRAMIDINESTRAND] = np.int8
    mydtypes[MUTATION] = str
    mydtypes[CONTEXT] = str
    #################################################

    #################################################
    #Read All rows
    #old way
    # mutation_df = pd.read_table(mutationFilePath, sep="\t", comment='#',dtype={'Sample':str,'Chromosome': str, 'Start': int, 'End':int, 'PyramidineStrand': int, 'Mutation':str, 'Context':str})
    #new way
    mutation_df = pd.read_table(mutationFilePath, sep="\t", comment='#',dtype=mydtypes)
    #################################################

    print('mutation_df[Chrom].unique()')
    print(mutation_df[CHROM].unique())

    listofSamples = mutation_df[SAMPLE].unique()
    print('Number of samples in single point mutations file: %d' %(len(listofSamples)))

    ##############################################################
    samplesWithAtLeast10KMutations2NumberofMutationsDict = {}

    for sample in listofSamples:
        numberofMutations =  len(mutation_df[mutation_df[SAMPLE] == sample])
        if (numberofMutations>=TEN_THOUSAND):
            samplesWithAtLeast10KMutations2NumberofMutationsDict[sample] = numberofMutations

    ########## debug starts #############
    print('samplesWithAtLeast10KMutations2NumberofMutationsDict')
    print(samplesWithAtLeast10KMutations2NumberofMutationsDict)

    print('len(samplesWithAtLeast10KMutations2NumberofMutationsDict)')
    print(len(samplesWithAtLeast10KMutations2NumberofMutationsDict))
    ########## debug ends ###############


    ##############################################################

    ##############################################################
    # first write down the header
    print('Signature\tNumberofMutationsWithProbability>=%f\tNumberofSamples' %MUTATION_SIGNATURE_PROBABILITY_THRESHOLD)
    for signature in signatures:
        signaturebased_df = mutation_df[mutation_df[signature] >= MUTATION_SIGNATURE_PROBABILITY_THRESHOLD]
        numberofMutations = len(signaturebased_df)
        if (numberofMutations>= TEN_THOUSAND):
            signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict[signature] = numberofMutations

        signatureBased_numberofSamples = len(signaturebased_df[SAMPLE].unique())
        print('%s\t%d\t%d' % (signature, signaturebased_df.shape[0], signatureBased_numberofSamples))
    ##############################################################


    ##############################################################
    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict = {}

    #We can not get rid of this part since we need the signatures with at least 10K mutations  for that sample
    # samplesWithAtLeast10KMutationsList and signaturesWithAtLeast10KEligibleMutationsList do not provide this information
    #Get the list of signatures for each sample that has at least 10K mutations with probability >= 0.5
    for sample in samplesWithAtLeast10KMutations2NumberofMutationsDict:
        for signature in signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
            #check if there are at least 10K mutations with probability >= 0.5 for this (sample,signature) pair
            numberofMutations = len(mutation_df[ ((mutation_df[SAMPLE]==sample) & (mutation_df[signature]>= MUTATION_SIGNATURE_PROBABILITY_THRESHOLD)) ])
            if (numberofMutations>= TEN_THOUSAND):
                if sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
                    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict[sample][signature] = numberofMutations
                else:
                    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict[sample]={}
                    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict[sample][signature] = numberofMutations
    ##############################################################

    print('###########################')
    print('size of %s in %d Bytes -- %f in GB ' %(mutationsWithSignatureBasedProbabilitiesFileName, sys.getsizeof(mutation_df),sys.getsizeof(mutation_df)/GIGABYTE_IN_BYTES))
    print('###########################')

    return signatures, samplesWithAtLeast10KMutations2NumberofMutationsDict, signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict, sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict, mutation_df
##################################################################

##################################################################
def readMutationsAndWriteChrBased(mutationsWithSignatureBasedProbabilitiesFileName):
    signatures, samplesWithAtLeast10KMutationsList, mutation_df = readMutations(mutationsWithSignatureBasedProbabilitiesFileName)

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
    print('write for %s ends' %(chr))
##################################################################



##########################################################################################
############### Common functions for Nucleosome Occupancy Analysis starts ################
##########################################################################################

########################################################################################
#We will accumulate signature2SplitArrayDict in signature2AccumulatedSplitsChrBasedArrayDict
def accumulateSignatureBasedArrays(signature2AccumulatedSplitsChrBasedArrayDict, signature2SplitArrayDict):
    for signatureKey in signature2SplitArrayDict.keys():
        signature2AccumulatedSplitsChrBasedArrayDict[signatureKey] += signature2SplitArrayDict[signatureKey]
########################################################################################

########################################################################################
def accumulateSampleBasedSignatureBasedArrays(sample2Signature2AccumulatedSplitsChrBasedArrayDict,sample2Signature2SplitArrayDict):
    for sample in sample2Signature2SplitArrayDict.keys():
        for signature in sample2Signature2SplitArrayDict[sample].keys():
            if sample in sample2Signature2AccumulatedSplitsChrBasedArrayDict:
                if signature in sample2Signature2AccumulatedSplitsChrBasedArrayDict[sample]:
                    sample2Signature2AccumulatedSplitsChrBasedArrayDict[sample][signature] += sample2Signature2SplitArrayDict[sample][signature]
                else:
                    sample2Signature2AccumulatedSplitsChrBasedArrayDict[sample][signature] = sample2Signature2SplitArrayDict[sample][signature]
            else:
                sample2Signature2AccumulatedSplitsChrBasedArrayDict[sample]= {}
                sample2Signature2AccumulatedSplitsChrBasedArrayDict[sample][signature] = sample2Signature2SplitArrayDict[sample][signature]
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
# Common for nucleosome occupancy analyis
# head chr3_wgEncodeSydhNsomeK562Sig.wig
#   chr3    60039   60043   1.2
#   chr3    60043   60049   1.3
#It seems that ends are exclusive therefore we do not need to add 1 to the end during slicing
def fillNucleosomeSignalArray(nucleosome_row,nucleosome_array):
    nucleosome_array[nucleosome_row['start']:nucleosome_row['end']] += nucleosome_row['signal']
########################################################################################

########################################################################################
# For nucleosome occupancy analyis
#This function is general can be move to TopographyCommons.py
#This is for single point mutations and signatures
def fillSplitBasedSignalArrayAndCountArrayForSPMs(mutation_row,nucleosome_array,signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,signature2SignalArrayDict,signature2CountArrayDict,allSinglePointMutationsSignalArray,allSinglePointMutationsCountArray):

    #Do this only once.
    ##########################################################
    # Case 1: start is very close to the chromosome start
    if (mutation_row[START]<plusOrMinus):
        print('Case 1: start is very close to the chromosome start --- mutation[Start]:%d mutation[End]:%d' %(mutation_row[START],mutation_row[END]))
        window_array = nucleosome_array[0:(mutation_row[END]+plusOrMinus+1)]
        window_array = np.pad(window_array, (plusOrMinus-mutation_row[START],0),'constant')
    # Case 2: start is very close to the chromosome end
    elif (mutation_row[END]+plusOrMinus > MAXIMUM_CHROMOSOME_LENGTH):
        print('mutation[Start]:%d mutation[End]:%d' %(mutation_row[START],mutation_row[END]))
        window_array = nucleosome_array[(mutation_row[START]-plusOrMinus):MAXIMUM_CHROMOSOME_LENGTH]
        window_array = np.pad(window_array, (0,mutation_row[END]+plusOrMinus-MAXIMUM_CHROMOSOME_LENGTH),'constant')
    #Case 3: No problem
    else:
        window_array = nucleosome_array[(mutation_row[START]-plusOrMinus):(mutation_row[END]+plusOrMinus+1)]
    ##########################################################

    #TODO: Is there a faster way than using for loop?
    ################# Signatures starts #######################
    for signature in signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
        if (mutation_row[signature] >= MUTATION_SIGNATURE_PROBABILITY_THRESHOLD):
            signature2SignalArrayDict[signature] += window_array
            signature2CountArrayDict[signature] += (window_array>0)
    ################# Signatures ends #########################

    ################ All single mutations starts ##############
    allSinglePointMutationsSignalArray += window_array
    allSinglePointMutationsCountArray += (window_array>0)
    ################ All single mutations ends ################
########################################################################################


########################################################################################
# For nucleosome occupancy analyis
#This function is general can be move to TopographyCommons.py
#This is for single point mutations and signatures
def  fillSplitBasedSignalArrayAndCountArrayForSPMsWithExtraSampleBased(mutation_row,
                                                  nucleosome_array,
                                                  maximum_chrom_size,
                                                  signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
                                                  samplesWithAtLeast10KMutations2NumberofMutationsDict,
                                                  sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,
                                                  signature2SignalArrayDict,
                                                  signature2CountArrayDict,
                                                  allSinglePointMutationsSignalArray,
                                                  allSinglePointMutationsCountArray,
                                                  sample2Signature2SignalArrayDict,
                                                  sample2Signature2CountArrayDict,
                                                  sample2AllSinglePointMutationsSignalArrayDict,
                                                  sample2AllSinglePointMutationsCountArrayDict):

    #Do this only once.
    #Subs have START and END same
    ##########################################################


    # Case 1: start is very close to the chromosome start
    if (mutation_row[START]<plusOrMinus):
        print('Case 1: start is very close to the chromosome start --- mutation[Start]:%d mutation[End]:%d' %(mutation_row[START],mutation_row[END]))
        window_array = nucleosome_array[0:(mutation_row[START]+plusOrMinus+1)]
        window_array = np.pad(window_array, (plusOrMinus-mutation_row[START],0),'constant',constant_values=(0,0))
    # Case 2: start is very close to the chromosome end
    elif (mutation_row[START]+plusOrMinus > maximum_chrom_size):
        print('Case2: start is very close to the chromosome end ---  mutation[Start]:%d mutation[End]:%d' %(mutation_row[START],mutation_row[END]))
        window_array = nucleosome_array[(mutation_row[START]-plusOrMinus):maximum_chrom_size]
        window_array = np.pad(window_array, (0,mutation_row[START]+plusOrMinus-maximum_chrom_size+1),'constant',constant_values=(0,0))
    #Case 3: No problem
    else:
        window_array = nucleosome_array[(mutation_row[START]-plusOrMinus):(mutation_row[END]+plusOrMinus+1)]
    ##########################################################

    #Get the sample at this mutation_row
    sample = mutation_row[SAMPLE]

    #TODO: Is there a faster way than using for loop?
    ################# Signatures starts #######################
    #We do not want to consider signatures with eligible mutations less than 10K
    for signature in signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict:
        if (mutation_row[signature] >= MUTATION_SIGNATURE_PROBABILITY_THRESHOLD):

            ####################################################
            signature2SignalArrayDict[signature] += window_array
            signature2CountArrayDict[signature] += (window_array>0)
            ####################################################

            ####################################################
            if (sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict) and (signature in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict[sample]):
                sample2Signature2SignalArrayDict[sample][signature] += window_array
                sample2Signature2CountArrayDict[sample][signature] += (window_array>0)
            ####################################################

    ################# Signatures ends #########################


    ################ All single mutations starts ##############
    allSinglePointMutationsSignalArray += window_array
    allSinglePointMutationsCountArray += (window_array>0)
    ################ All single mutations ends ################

    ################ Sample Based All single mutations starts ##############
    if (sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict):
        sample2AllSinglePointMutationsSignalArrayDict[sample] += window_array
        sample2AllSinglePointMutationsCountArrayDict[sample] += (window_array > 0)
    ################ Sample Based All single mutations ends ##############

    # #For debug starts FEB 28, 2019
    # #for each mutation
    # if (sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict and sample=='PD24197a'):
    #     count_array = np.zeros(windowSize)
    #     count_array += (window_array > 0)
    #
    #     signal_array_sum = np.sum(window_array)
    #     signal_array_max = np.amax(window_array)
    #     signal_array_min = np.amin(window_array)
    #
    #     count_array_sum = np.sum(count_array)
    #     count_array_max = np.amax(count_array)
    #     count_array_min = np.amin(count_array)
    #
    #     if (signal_array_sum>0):
    #         signalFilename = '%s_chr%s_%d_signal_array_sum_%s_max_%s_min_%s.txt' % (sample, mutation_row[CHROM], mutation_row[START], str(signal_array_sum).replace('.','_point'), str(signal_array_max).replace('.','_point'),str(signal_array_min).replace('.','_point'))
    #         signalFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, '560_BRCA_WGS',DATA, NUCLEOSOMEOCCUPANCY, signalFilename)
    #         window_array.tofile(file=signalFilePath, sep="\t", format="%s"),
    #
    #     if (count_array_sum>0):
    #         countFilename = '%s_chr%s_%d_count_array_sum_%s_max_%s_min_%s.txt' % (sample, mutation_row[CHROM], mutation_row[START], str(count_array_sum).replace('.','_point'), str(count_array_max).replace('.','_point'),str(count_array_min).replace('.','_point'))
    #         countFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, '560_BRCA_WGS',DATA, NUCLEOSOMEOCCUPANCY, countFilename)
    #         count_array.tofile(file=countFilePath, sep="\t", format="%s")
    # #For debug ends FEB 28, 2019
########################################################################################


########################################################################################
def fillSplitBasedSignalArrayAndCountArrayForIndelsWithExtraSampleBased(indel_row,nucleosome_array,maximum_chrom_size,sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,allIndelsSignalArray,allIndelsCountArray,sample2AllIndelsSignalArrayDict,sample2AllIndelsCountArrayDict):

    ##########################################################
    # Case 1: start is very close to the chromosome start
    if (indel_row[START]<plusOrMinus):
        print('Case 1: start is very close to the chromosome start --- mutation[Start]:%d mutation[Start]:%d' %(indel_row[START],indel_row[START]))
        window_array = nucleosome_array[0:(indel_row[START]+plusOrMinus+1)]
        window_array = np.pad(window_array, (plusOrMinus-indel_row[START],0),'constant')
    # Case 2: start is very close to the chromosome end
    elif (indel_row[START]+plusOrMinus > maximum_chrom_size):
        print('mutation[Start]:%d mutation[End]:%d' %(indel_row[START],indel_row[START]))
        window_array = nucleosome_array[(indel_row[START]-plusOrMinus):maximum_chrom_size]
        window_array = np.pad(window_array, (0,indel_row[START]+plusOrMinus-maximum_chrom_size+1),'constant')
    #Case 3: No problem
    else:
        window_array = nucleosome_array[(indel_row[START]-plusOrMinus):(indel_row[START]+plusOrMinus+1)]
    ##########################################################

    #Get the sample at this mutation_row
    sample = indel_row[SAMPLE]


    ################ All single mutations starts ##############
    allIndelsSignalArray += window_array
    allIndelsCountArray += (window_array>0)
    ################ All single mutations ends ################

    if (sample in sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict):
        sample2AllIndelsSignalArrayDict[sample] += window_array
        sample2AllIndelsCountArrayDict[sample] += (window_array>0)
########################################################################################


########################################################################################
# For nucleosome occupancy analyis
#This is for indels
#Please notice that I have used indel_row['Start'] only not indel_row['End']
def fillSplitBasedSignalArrayAndCountArrayForIndels(indel_row,nucleosome_array,allIndelsSignalArray,allIndelsCountArray):

    #Do this only once.
    ##########################################################
    # Case 1: start is very close to the chromosome start
    if (indel_row[START]<plusOrMinus):
        print('Case 1: start is very close to the chromosome start --- mutation[Start]:%d mutation[Start]:%d' %(indel_row[START],indel_row[START]))
        window_array = nucleosome_array[0:(indel_row[START]+plusOrMinus+1)]
        window_array = np.pad(window_array, (plusOrMinus-indel_row[START],0),'constant')
    # Case 2: start is very close to the chromosome end
    elif (indel_row[START]+plusOrMinus > MAXIMUM_CHROMOSOME_LENGTH):
        print('mutation[Start]:%d mutation[End]:%d' %(indel_row[START],indel_row[START]))
        window_array = nucleosome_array[(indel_row[START]-plusOrMinus):MAXIMUM_CHROMOSOME_LENGTH]
        window_array = np.pad(window_array, (0,indel_row[START]+plusOrMinus-MAXIMUM_CHROMOSOME_LENGTH),'constant')
    #Case 3: No problem
    else:
        window_array = nucleosome_array[(indel_row[START]-plusOrMinus):(indel_row[START]+plusOrMinus+1)]
    ##########################################################

    ################ All single mutations starts ##############
    allIndelsSignalArray += window_array
    allIndelsCountArray += (window_array>0)
    ################ All single mutations ends ################

########################################################################################


########################################################################################
def computeAverageNucleosomeOccupancyArray(signalArray,countArray):
    averageArray =  np.zeros(windowSize)

    np.seterr(divide='ignore', invalid='ignore')
    if (np.any(countArray)):
        averageArray = np.divide(signalArray,countArray)
    np.seterr(divide='raise', invalid='ignore')

    return averageArray
########################################################################################

########################################################################################
#Both "all single point mutations" and "all indels" use this function
def writeAverageNucleosomeOccupancyFiles(allMutationsAccumulatedAllChromsSignalArray,allMutationsAccumulatedAllChromsCountArray,outputDir,jobname,nucleosomeOccupancyAnalysisType):
    os.makedirs(os.path.join(outputDir,jobname,DATA,NUCLEOSOMEOCCUPANCY,nucleosomeOccupancyAnalysisType), exist_ok=True)

    averageNucleosomeSignalArray = computeAverageNucleosomeOccupancyArray(allMutationsAccumulatedAllChromsSignalArray, allMutationsAccumulatedAllChromsCountArray)

    accumulatedSignalFilename = '%s_AccumulatedSignalArray.txt' %(jobname)
    accumulatedCountFilename = '%s_AccumulatedCountArray.txt' %(jobname)
    averageNucleosomeSignalFilename = '%s_AverageNucleosomeSignalArray.txt' %(jobname)

    accumulatedSignalFilePath = os.path.join(outputDir, jobname,DATA, NUCLEOSOMEOCCUPANCY, nucleosomeOccupancyAnalysisType, accumulatedSignalFilename)
    accumulatedCountFilePath = os.path.join(outputDir, jobname, DATA,NUCLEOSOMEOCCUPANCY, nucleosomeOccupancyAnalysisType, accumulatedCountFilename)
    averageNucleosomeSignalFilePath = os.path.join(outputDir,jobname, DATA, NUCLEOSOMEOCCUPANCY, nucleosomeOccupancyAnalysisType,averageNucleosomeSignalFilename)

    allMutationsAccumulatedAllChromsSignalArray.tofile(file=accumulatedSignalFilePath, sep="\t", format="%s")
    allMutationsAccumulatedAllChromsCountArray.tofile(file=accumulatedCountFilePath, sep="\t", format="%s")
    averageNucleosomeSignalArray.tofile(file=averageNucleosomeSignalFilePath, sep="\t", format="%s")
########################################################################################

########################################################################################
def writeSignatureBasedAverageNucleosomeOccupancyFiles(signature2AccumulatedAllChromsSignalArrayDict,signature2AccumulatedAllChromsCountArrayDict,outputDir,jobname):
    os.makedirs(os.path.join(outputDir,jobname,DATA,NUCLEOSOMEOCCUPANCY,SIGNATUREBASED), exist_ok=True)

    for signature in signature2AccumulatedAllChromsSignalArrayDict.keys():
        signalArray = signature2AccumulatedAllChromsSignalArrayDict[signature]
        countArray = signature2AccumulatedAllChromsCountArrayDict[signature]
        averageNucleosomeSignalArray = computeAverageNucleosomeOccupancyArray(signalArray,countArray)

        #To provide filename with no space in signature name
        #signatureWithNoSpace = signature.replace(' ','')

        accumulatedSignalFilename = '%s_AccumulatedSignalArray.txt' %(signature)
        accumulatedCountFilename = '%s_AccumulatedCountArray.txt' %(signature)
        averageNucleosomeSignalFilename = '%s_AverageNucleosomeSignalArray.txt' %(signature)

        accumulatedSignalFilePath = os.path.join(outputDir,jobname,DATA,NUCLEOSOMEOCCUPANCY,SIGNATUREBASED,accumulatedSignalFilename)
        accumulatedCountFilePath = os.path.join(outputDir,jobname,DATA,NUCLEOSOMEOCCUPANCY,SIGNATUREBASED,accumulatedCountFilename)
        averageNucleosomeSignalFilePath = os.path.join(outputDir,jobname,DATA,NUCLEOSOMEOCCUPANCY,SIGNATUREBASED,averageNucleosomeSignalFilename)

        signalArray.tofile(file=accumulatedSignalFilePath, sep="\t",format="%s")
        countArray.tofile(file=accumulatedCountFilePath, sep="\t", format="%s")
        averageNucleosomeSignalArray.tofile(file=averageNucleosomeSignalFilePath, sep="\t", format="%s")
########################################################################################


########################################################################################
def writeSampleBasedSignatureBasedAverageNucleosomeOccupancyFiles(sample2Signature2AccumulatedAllChromsSignalArrayDict,
                                                                sample2Signature2AccumulatedAllChromsCountArrayDict,outputDir,jobname):


    os.makedirs(os.path.join(outputDir,jobname,DATA,NUCLEOSOMEOCCUPANCY,SIGNATUREBASED), exist_ok=True)

    for sample in sample2Signature2AccumulatedAllChromsSignalArrayDict:
        for signature in sample2Signature2AccumulatedAllChromsSignalArrayDict[sample]:
            signalArray = sample2Signature2AccumulatedAllChromsSignalArrayDict[sample][signature]
            countArray = sample2Signature2AccumulatedAllChromsCountArrayDict[sample][signature]
            averageNucleosomeSignalArray = computeAverageNucleosomeOccupancyArray(signalArray, countArray)

            # To provide filename with no space in signature name
            # signatureWithNoSpace = signature.replace(' ','')

            accumulatedSignalFilename = '%s_%s_AccumulatedSignalArray.txt' %(signature,sample)
            accumulatedCountFilename = '%s_%s_AccumulatedCountArray.txt' %(signature,sample)
            averageNucleosomeSignalFilename = '%s_%s_AverageNucleosomeSignalArray.txt' %(signature,sample)

            accumulatedSignalFilePath = os.path.join(outputDir,jobname, DATA, NUCLEOSOMEOCCUPANCY, SIGNATUREBASED,accumulatedSignalFilename)
            accumulatedCountFilePath = os.path.join(outputDir,jobname, DATA, NUCLEOSOMEOCCUPANCY, SIGNATUREBASED,accumulatedCountFilename)
            averageNucleosomeSignalFilePath = os.path.join(outputDir,jobname, DATA, NUCLEOSOMEOCCUPANCY, SIGNATUREBASED,averageNucleosomeSignalFilename)

            signalArray.tofile(file=accumulatedSignalFilePath, sep="\t", format="%s")
            countArray.tofile(file=accumulatedCountFilePath, sep="\t", format="%s")
            averageNucleosomeSignalArray.tofile(file=averageNucleosomeSignalFilePath, sep="\t", format="%s")
########################################################################################


########################################################################################
#All Mutations can be single point mutations or indels
def writeSampleBasedAverageNucleosomeOccupancyFiles(sample2AllMutationsAccumulatedAllChromsSignalArrayDict,
                                                    sample2AllMutationsAccumulatedAllChromsCountArrayDict,
                                                    outputDir,
                                                    jobname,
                                                    nucleosomeOccupancyAnalysisType):

    os.makedirs(os.path.join(outputDir,jobname,DATA,NUCLEOSOMEOCCUPANCY,nucleosomeOccupancyAnalysisType), exist_ok=True)


    for sample in sample2AllMutationsAccumulatedAllChromsSignalArrayDict:
        allMutationsAccumulatedAllChromsSignalArray = sample2AllMutationsAccumulatedAllChromsSignalArrayDict[sample]
        allMutationsAccumulatedAllChromsCountArray = sample2AllMutationsAccumulatedAllChromsCountArrayDict[sample]

        averageNucleosomeSignalArray = computeAverageNucleosomeOccupancyArray(allMutationsAccumulatedAllChromsSignalArray, allMutationsAccumulatedAllChromsCountArray)

        accumulatedSignalFilename = '%s_%s_AccumulatedSignalArray.txt' %(sample,jobname)
        accumulatedCountFilename = '%s_%s_AccumulatedCountArray.txt' %(sample,jobname)
        averageNucleosomeSignalFilename = '%s_%s_AverageNucleosomeSignalArray.txt' %(sample,jobname)

        accumulatedSignalFilePath = os.path.join(outputDir, jobname,DATA, NUCLEOSOMEOCCUPANCY, nucleosomeOccupancyAnalysisType, accumulatedSignalFilename)
        accumulatedCountFilePath = os.path.join(outputDir, jobname, DATA,NUCLEOSOMEOCCUPANCY, nucleosomeOccupancyAnalysisType, accumulatedCountFilename)
        averageNucleosomeSignalFilePath = os.path.join(outputDir,jobname, DATA, NUCLEOSOMEOCCUPANCY, nucleosomeOccupancyAnalysisType,averageNucleosomeSignalFilename)

        allMutationsAccumulatedAllChromsSignalArray.tofile(file=accumulatedSignalFilePath, sep="\t", format="%s")
        allMutationsAccumulatedAllChromsCountArray.tofile(file=accumulatedCountFilePath, sep="\t", format="%s")
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

##################################################################
def readIndelsAndWriteChrBasedParallel(outputDir,jobname,indelsFileName):
    sample2NumberofIndelsDict,indels_df = readIndels(indelsFileName)

    print('For debugging purposes indels_df.columns')
    print(indels_df.columns)

    indels_df_grouped= indels_df.groupby(CHROM)

    print('len(indels_df_grouped)')
    print(len(indels_df_grouped))

    os.makedirs(os.path.join(outputDir,jobname,DATA, CHRBASED),exist_ok=True)

    #########################################################
    ############### Write Unique Chrnames starts ############
    #########################################################
    uniqueChrNames = indels_df[CHROM].unique()
    # The unique values are returned as a NumPy array

    # Write uniqueChrNames to a file
    filename = ChrNamesInIndelsFilename

    ChrNamesFile = os.path.join(outputDir,jobname,DATA,filename)

    np.savetxt(ChrNamesFile,uniqueChrNames, delimiter='\t', fmt='%s')
    #########################################################
    ############### Write Unique Chrnames ends ##############
    #########################################################

    #############################################################################################################################
    ######################## Write samplesWithAtLeast10KMutations2NumberofMutationsDict starts ##################################
    #############################################################################################################################
    writeDictionaryUnderDataDirectory(sample2NumberofIndelsDict,outputDir,jobname,Samples2NumberofIndelsDictFilename)
    #############################################################################################################################
    ######################## Write samplesWithAtLeast10KMutations2NumberofMutationsDict starts ##################################
    #############################################################################################################################

    #########################################################
    # l = multiprocessing.Lock()
    numofProcesses = multiprocessing.cpu_count()
    # pool = multiprocessing.Pool(numofProcesses,initializer=init, initargs=(l,))
    pool = multiprocessing.Pool(numofProcesses)
    #########################################################

    print('For debugging purposes number of cores: %d' %(numofProcesses))
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
def readMutationsAndWriteChrBasedParallel(outputDir,jobname,mutationsWithSignatureBasedProbabilitiesFileName):
    signatures, \
    samplesWithAtLeast10KMutations2NumberofMutationsDict, \
    signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,\
    sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,\
    mutation_df = readMutations(mutationsWithSignatureBasedProbabilitiesFileName)

    mutation_df_grouped= mutation_df.groupby(CHROM)
    os.makedirs(os.path.join(outputDir,jobname,DATA,CHRBASED),exist_ok=True)

    #########################################################
    ############### Write signatures starts #################
    #########################################################
    signatures_array = np.array(signatures)

    # Write signatures_array to a file
    SignaturesFile = os.path.join(outputDir,jobname,DATA, SignatureFilename)

    np.savetxt(SignaturesFile,signatures_array,delimiter='\t', fmt='%s')
    #########################################################
    ############### Write signatures ends ###################
    #########################################################

    #############################################################################################################################
    ######################## Write samplesWithAtLeast10KMutations2NumberofMutationsDict starts ##################################
    #############################################################################################################################
    writeDictionaryUnderDataDirectory(samplesWithAtLeast10KMutations2NumberofMutationsDict,outputDir,jobname,SamplesWithAtLeast10KMutations2NumberofMutationsDictFilename)
    #############################################################################################################################
    ######################## Write samplesWithAtLeast10KMutations2NumberofMutationsDict starts ##################################
    #############################################################################################################################

    #############################################################################################################################
    ######################## Write signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict starts ##################################
    #############################################################################################################################
    writeDictionaryUnderDataDirectory(signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,outputDir,jobname,SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilename)
    #############################################################################################################################
    ######################## Write signaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict starts ##################################
    #############################################################################################################################

    #############################################################################################################################
    ############### Write sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict starts #########################
    #############################################################################################################################
    writeDictionaryUnderDataDirectory(sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict,outputDir,jobname,Sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDictFilename)
    #############################################################################################################################
    ############### Write sample2SignaturesWithAtLeast10KEligibleMutations2NumberofMutationsDict ends ###########################
    #############################################################################################################################


    #########################################################
    ############### Write Unique Chrnames starts ############
    #########################################################
    uniqueChrNames = mutation_df[CHROM].unique()

    ChrNamesFile = os.path.join(outputDir,jobname,DATA,ChrNamesInSPMsFilename)
    np.savetxt(ChrNamesFile,uniqueChrNames, delimiter='\t', fmt='%s')
    #########################################################
    ############### Write Unique Chrnames ends ##############
    #########################################################

    #########################################################
    l = multiprocessing.Lock()
    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses,initializer=init, initargs=(l,))
    #########################################################

    poolInputList =[]

    # Get the filename at the end of the full path
    mutationsWithSignatureBasedProbabilitiesFileName = os.path.basename(mutationsWithSignatureBasedProbabilitiesFileName)

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
        numberofMutations =  len(indels_df[indels_df[SAMPLE] == sample])
        sample2NumberofIndelsDict[sample] = numberofMutations
    ###########################################################################

    print('Number of samples in indels file: %d' %(len(indels_df[SAMPLE].unique())))
    print('Number of indels: %d' %(indels_df.shape[0]))
    print('Indels file shape')
    print(indels_df.shape)

    # drop columns 'Reference','Mutation','Type','Length','Category'
    # Do not drop Length column. If length >= 3bp indel type is Microhomology otherwise < 3bp indel type is repeat-med indel
    indels_df.drop([REF,ALT,TYPE,CATEGORY], axis=1, inplace=True, errors='ignore')

    #Add a new column called Count
    grouped_mutation_df_by_sample = indels_df.groupby(SAMPLE)
    indels_df[COUNT] = indels_df.groupby(SAMPLE)[SAMPLE].transform('count')

    # grouped_mutation_df.groups
    print('Number of samples in indels file:%s' % len(grouped_mutation_df_by_sample))

    # Add a new column called MMR (Mis Match Repair)
    if 'Count' in indels_df.columns:
        indels_df[MMR] = np.where(indels_df[COUNT] >= TEN_THOUSAND ,DEFICIENT,PROFICIENT)

    return sample2NumberofIndelsDict,indels_df
##################################################################




##################################################################
def readBED(bedFilename):
    # Read the file w.r.t. the current folder
    bedFilePath = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,REPLICATION,bedFilename)

    columns = ['chr', 'start', 'end', 'column4','column5','column6','column7','column8','column9']

    #BED files are seperated by tab
    bed_df = pd.read_table(bedFilePath, sep='\t', comment='#', header=None, names=columns)

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
    # Read the file w.r.t. the current folder

    repliSeqFilePath = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,REPLICATION,repliseqDataFilename)
    # print(repliSeqFilePath)
    #todo is there a better/more clean way to handle this?
    #This sep='*' is used for reading line by line without different number of columns error for different lines
    repliSeq_unprocessed_df = pd.read_table(repliSeqFilePath, sep='*', comment='#', header=None)
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
def readWaveletSmoothedRepliSeqSignal(repliseqDataFilename):
    # Read the file w.r.t. the current folder

    waveletFilePath = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,REPLICATION,repliseqDataFilename)
    wavelet_unprocessed_df = pd.read_table(waveletFilePath, sep="\t", comment='#', header=None)
    # print('debug starts')
    # print('############## wavelet_unprocessed_df ##############')
    # print(wavelet_unprocessed_df.shape)
    # print(wavelet_unprocessed_df.head())
    # print('############## wavelet_unprocessed_df ##############')
    # print('debug ends')

    return wavelet_unprocessed_df
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
def prepareMutationProbabilityList(startMutationProbability, endMutationProbability, step):
    # print('type(startMutationProbability):%s type(endMutationProbability):%s type(step):%s', (type(startMutationProbability), type(endMutationProbability), type(step)))
    startMutationProbability = float("{0:.2f}".format(startMutationProbability))

    mutationProbabilityList = []
    while startMutationProbability <= endMutationProbability:
        #Put the rounded mutation probability to the list
        mutationProbabilityList.append(round(startMutationProbability,2))
        startMutationProbability += step
        startMutationProbability = float("{0:.2f}".format(startMutationProbability))

    # endMutationProbability += step
    # mutationProbabilityList = np.arange(startMutationProbability,endMutationProbability,step)

    return mutationProbabilityList
##################################################################

########################################################################
#this is called from TranscriptionStrandBiasAnalysis.py and ReplicationStrandBiasAnalysis.py
def updateDictionaries(mutation_row,
                        mutationType,
                        mutationSample,
                        mutationType2Strand2CountDict,
                        mutationType2Sample2Strand2CountDict,
                        mutationProbability2Signature2Strand2CountDict,
                        mutationProbability2Signature2Sample2Strand2CountDict,
                        strand,
                        signatureList,
                        mutationProbabilityList):

    #################################################################################################
    # Update1: mutationType2TranscriptionStrand2CountDict
    if mutationType in mutationType2Strand2CountDict:
        if strand in mutationType2Strand2CountDict[mutationType]:
            mutationType2Strand2CountDict[mutationType][strand] += 1
        else:
            mutationType2Strand2CountDict[mutationType][strand] = 1
    else:
        mutationType2Strand2CountDict[mutationType] = {}
        mutationType2Strand2CountDict[mutationType][strand] = 1
    #################################################################################################

    #################################################################################################
    # Update2: mutationType2Sample2TranscriptionStrand2CountDict
    if mutationType in mutationType2Sample2Strand2CountDict:
        if mutationSample in mutationType2Sample2Strand2CountDict[mutationType]:
            if strand in mutationType2Sample2Strand2CountDict[mutationType][mutationSample]:
                mutationType2Sample2Strand2CountDict[mutationType][mutationSample][strand] += 1
            else:
                mutationType2Sample2Strand2CountDict[mutationType][mutationSample][strand] = 1
        else:
            mutationType2Sample2Strand2CountDict[mutationType][mutationSample] = {}
            mutationType2Sample2Strand2CountDict[mutationType][mutationSample][strand] = 1

    else:
        mutationType2Sample2Strand2CountDict[mutationType] = {}
        mutationType2Sample2Strand2CountDict[mutationType][mutationSample] = {}
        mutationType2Sample2Strand2CountDict[mutationType][mutationSample][strand] = 1
    #################################################################################################

    #################################################################################################
    # Update3: mutationProbability2Signature2TranscriptionStrand2CountDict
    for mutationProbability in mutationProbabilityList:
        for signature in signatureList:
            if (mutation_row[signature] >= mutationProbability):
                if mutationProbability in mutationProbability2Signature2Strand2CountDict:
                    if signature in mutationProbability2Signature2Strand2CountDict[mutationProbability]:
                        if strand in mutationProbability2Signature2Strand2CountDict[mutationProbability][signature]:
                            mutationProbability2Signature2Strand2CountDict[mutationProbability][signature][strand] += 1
                        else:
                            mutationProbability2Signature2Strand2CountDict[mutationProbability][signature][strand] = 1
                    else:
                        mutationProbability2Signature2Strand2CountDict[mutationProbability][signature] = {}
                        mutationProbability2Signature2Strand2CountDict[mutationProbability][signature][strand] = 1
                else:
                    mutationProbability2Signature2Strand2CountDict[mutationProbability] = {}
                    mutationProbability2Signature2Strand2CountDict[mutationProbability][signature] = {}
                    mutationProbability2Signature2Strand2CountDict[mutationProbability][signature][strand] = 1
    #################################################################################################

    #################################################################################################
    # Update4: mutationProbability2Signature2Sample2TranscriptionStrand2CountDict
    for mutationProbability in mutationProbabilityList:
        for signature in signatureList:
            if (mutation_row[signature] >= mutationProbability):
                if mutationProbability in mutationProbability2Signature2Sample2Strand2CountDict:
                    if signature in mutationProbability2Signature2Sample2Strand2CountDict[mutationProbability]:
                        if mutationSample in mutationProbability2Signature2Sample2Strand2CountDict[mutationProbability][signature]:
                            if strand in mutationProbability2Signature2Sample2Strand2CountDict[mutationProbability][signature][mutationSample]:
                                mutationProbability2Signature2Sample2Strand2CountDict[mutationProbability][signature][mutationSample][strand] += 1
                            else:
                                mutationProbability2Signature2Sample2Strand2CountDict[mutationProbability][signature][mutationSample][strand] = 1
                        else:
                            mutationProbability2Signature2Sample2Strand2CountDict[mutationProbability][signature][mutationSample] = {}
                            mutationProbability2Signature2Sample2Strand2CountDict[mutationProbability][signature][mutationSample][strand] = 1
                    else:
                        mutationProbability2Signature2Sample2Strand2CountDict[mutationProbability][signature] = {}
                        mutationProbability2Signature2Sample2Strand2CountDict[mutationProbability][signature][mutationSample] = {}
                        mutationProbability2Signature2Sample2Strand2CountDict[mutationProbability][signature][mutationSample][strand] = 1
                else:
                    mutationProbability2Signature2Sample2Strand2CountDict[mutationProbability] = {}
                    mutationProbability2Signature2Sample2Strand2CountDict[mutationProbability][signature] = {}
                    mutationProbability2Signature2Sample2Strand2CountDict[mutationProbability][signature][mutationSample] = {}
                    mutationProbability2Signature2Sample2Strand2CountDict[mutationProbability][signature][mutationSample][strand] = 1
    #################################################################################################

########################################################################



########################################################################
def accumulate(listofTuples,
               accumulatedAllChromosomesMutationType2Strand2CountDict,
               accumulatedAllChromosomesMutationType2Sample2Strand2CountDict,
               accumulatedAllChromosomesMutationProbability2Signature2Strand2CountDict,
               accumulatedAllChromosomesMutationProbability2Signature2Sample2Strand2CountDict):

    for mytuple in listofTuples:
        chrBasedMutationType2Strand2CountDict = mytuple[0]
        chrBasedMutationType2Sample2Strand2CountDict = mytuple[1]
        chrBasedMutationProbability2Signature2Strand2CountDict = mytuple[2]
        chrBasedMutationProbability2Signature2Sample2Strand2CountDict = mytuple[3]

        #Accumulate mutationType2Strand2CountDict
        if (chrBasedMutationType2Strand2CountDict is not None):
            for mutationType,strand2CountDict in chrBasedMutationType2Strand2CountDict.items():
                if mutationType in accumulatedAllChromosomesMutationType2Strand2CountDict:
                    for strand, count in strand2CountDict.items():
                        if strand in accumulatedAllChromosomesMutationType2Strand2CountDict[mutationType]:
                            accumulatedAllChromosomesMutationType2Strand2CountDict[mutationType][strand] += count
                        else:
                            accumulatedAllChromosomesMutationType2Strand2CountDict[mutationType][strand] = count
                else:
                    accumulatedAllChromosomesMutationType2Strand2CountDict[mutationType] = strand2CountDict

        #Accumulate mutationType2Sample2Strand2CountDict
        if (chrBasedMutationType2Sample2Strand2CountDict is not None):
            for mutationType, sample2Strand2CountDict in chrBasedMutationType2Sample2Strand2CountDict.items():
                if mutationType in accumulatedAllChromosomesMutationType2Sample2Strand2CountDict:
                    for sample, strand2CountDict in sample2Strand2CountDict.items():
                        if sample in accumulatedAllChromosomesMutationType2Sample2Strand2CountDict[mutationType]:
                            for strand, count in strand2CountDict.items():
                                if strand in accumulatedAllChromosomesMutationType2Sample2Strand2CountDict[mutationType][sample]:
                                    accumulatedAllChromosomesMutationType2Sample2Strand2CountDict[mutationType][sample][strand] += count
                                else:
                                    accumulatedAllChromosomesMutationType2Sample2Strand2CountDict[mutationType][sample][strand] = count
                        else:
                            accumulatedAllChromosomesMutationType2Sample2Strand2CountDict[mutationType][sample] = strand2CountDict
                else:
                    accumulatedAllChromosomesMutationType2Sample2Strand2CountDict[mutationType] = sample2Strand2CountDict

        #Accumulate mutationProbability2Signature2Strand2CountDict
        if (chrBasedMutationProbability2Signature2Strand2CountDict is not None):
            for mutationProbability, signature2StrandCountDict in chrBasedMutationProbability2Signature2Strand2CountDict.items():
                if mutationProbability in accumulatedAllChromosomesMutationProbability2Signature2Strand2CountDict:
                    for signature, strand2CountDict in signature2StrandCountDict.items():
                        if signature in accumulatedAllChromosomesMutationProbability2Signature2Strand2CountDict[mutationProbability]:
                            for strand, count in strand2CountDict.items():
                                if strand in accumulatedAllChromosomesMutationProbability2Signature2Strand2CountDict[mutationProbability][signature]:
                                    accumulatedAllChromosomesMutationProbability2Signature2Strand2CountDict[mutationProbability][signature][strand] += count
                                else:
                                    accumulatedAllChromosomesMutationProbability2Signature2Strand2CountDict[mutationProbability][signature][strand] = count
                        else:
                            accumulatedAllChromosomesMutationProbability2Signature2Strand2CountDict[mutationProbability][signature] = strand2CountDict
                else:
                    accumulatedAllChromosomesMutationProbability2Signature2Strand2CountDict[mutationProbability] = signature2StrandCountDict

        #Accumulate mutationProbability2Signature2Sample2Strand2CountDict
        if (chrBasedMutationProbability2Signature2Sample2Strand2CountDict is not None):
            for mutationProbability,signature2Sample2Strand2CountDict in chrBasedMutationProbability2Signature2Sample2Strand2CountDict.items():
                if mutationProbability in accumulatedAllChromosomesMutationProbability2Signature2Sample2Strand2CountDict:
                    for signature, sample2Strand2CountDict in signature2Sample2Strand2CountDict.items():
                        if signature in accumulatedAllChromosomesMutationProbability2Signature2Sample2Strand2CountDict[mutationProbability]:
                            for sample, strand2CountDict in sample2Strand2CountDict.items():
                                if sample in accumulatedAllChromosomesMutationProbability2Signature2Sample2Strand2CountDict[mutationProbability][signature]:
                                    for strand, count in strand2CountDict.items():
                                        if strand in accumulatedAllChromosomesMutationProbability2Signature2Sample2Strand2CountDict[mutationProbability][signature][sample]:
                                            accumulatedAllChromosomesMutationProbability2Signature2Sample2Strand2CountDict[mutationProbability][signature][sample][strand] += count
                                        else:
                                            accumulatedAllChromosomesMutationProbability2Signature2Sample2Strand2CountDict[mutationProbability][signature][sample][strand] = count
                                else:
                                    accumulatedAllChromosomesMutationProbability2Signature2Sample2Strand2CountDict[mutationProbability][signature][sample] = strand2CountDict
                        else:
                            accumulatedAllChromosomesMutationProbability2Signature2Sample2Strand2CountDict[mutationProbability][signature] = sample2Strand2CountDict
                else:
                    accumulatedAllChromosomesMutationProbability2Signature2Sample2Strand2CountDict[mutationProbability] = signature2Sample2Strand2CountDict
########################################################################



########################################################################
def calculateRatio(accumulatedAllChromosomesMutationProbability2Signature2StrandCountDict,accumulatedAllChromosomesMutationProbability2Signature2RatioDict,strandList):
    numeratorStrand = strandList[0]
    denominatorStrand = strandList[1]

    for mutationProbability, signature2StrandCountDict in accumulatedAllChromosomesMutationProbability2Signature2StrandCountDict.items():
        accumulatedAllChromosomesMutationProbability2Signature2RatioDict[mutationProbability] = {}
        for signature, strandCountDict in signature2StrandCountDict.items():
            if ((numeratorStrand in strandCountDict) and (denominatorStrand in strandCountDict)):
                if((strandCountDict[numeratorStrand] + strandCountDict[denominatorStrand])>=ONE_THOUSAND):
                    accumulatedAllChromosomesMutationProbability2Signature2RatioDict[mutationProbability][signature] = (strandCountDict[numeratorStrand]) / (strandCountDict[numeratorStrand] + strandCountDict[denominatorStrand])
########################################################################


########################################################################
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
                if (strand2CountDict[numeratorStrand]+strand2CountDict[denominatorStrand] >= ONE_THOUSAND):
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
#To write samples with signatures with at least 10K eligible mutations
def writeDictionaryUnderDataDirectory(dictionary,outputDir,jobname,filename):

    os.makedirs(os.path.join(outputDir,jobname,DATA), exist_ok=True)
    filePath = os.path.join(outputDir,jobname,DATA,filename)

    with open(filePath, 'w') as file:
        file.write(json.dumps(dictionary))
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
    os.makedirs(os.path.join(filePath),exist_ok=True)

    with open(filePath, "w") as f:
        for item in list:
            f.write("%s\n" % item)
########################################################################

########################################################################
def append2File(item,filePath):
    os.makedirs(os.path.join(filePath),exist_ok=True)

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


