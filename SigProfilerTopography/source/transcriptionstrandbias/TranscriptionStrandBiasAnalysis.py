# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018 Burcak Otlu



# Version2
# This version use np.arrays
# Right now transcription strand bias analysis works for single point mutations and signatures.


#############################################################
# Constraints, Thresholds
# Please note that for sample based transcription strand bias analysis
# We consider samples with at least 1000 mutations both on transcribed and non-transcribed strands.
#############################################################

#############################################################
# What is transcription strand bias?
# It is the ratio of = (number of mutations on transcribed strand) / (number of mutations on un-transcribed strand)
#############################################################

import sys
import os

#############################################################
current_abs_path = os.path.abspath(os.path.dirname(__file__))
print('TranscriptionStrandBiasAnalysis.py current_abs_path:%s' %(current_abs_path))
#############################################################

commonsPath = os.path.join(current_abs_path, '..','commons')
sys.path.append(commonsPath)

from SigProfilerTopography.source.commons.TopographyCommons import *

########################################################################
# Each row is a pandas series in fact
# labels = ['#bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
# using cdsStart and cdsEnd give error since there are intervals of 0 length in that case
# such as cdsStart 33623833 cdsEnd 33623833
#Please notice that same genomic loci on a certain strand can be transcribed and nontranscribed at the same time.
#However same genomic loci on a certain strand can be leading or lagging but not both
#TODO Is (transcription_row['txEnd']+1) allright?
def fillTranscriptionArrays(transcription_row,chrBased_transcription_plus_array,chrBased_transcription_minus_array,transcriptsSource):

    if (transcriptsSource==NCBI):
        if (transcription_row['strand']==PLUS):
            chrBased_transcription_plus_array[transcription_row['txStart']:(transcription_row['txEnd']+1)] = 1
        elif (transcription_row['strand']==MINUS):
            chrBased_transcription_minus_array[transcription_row['txStart']:(transcription_row['txEnd']+1)] = -1

    elif (transcriptsSource==ENSEMBL):
        if (transcription_row['strand']==1):
            chrBased_transcription_plus_array[transcription_row['txStart']:(transcription_row['txEnd']+1)] = 1
        elif (transcription_row['strand']==-1):
            chrBased_transcription_minus_array[transcription_row['txStart']:(transcription_row['txEnd']+1)] = -1

########################################################################



########################################################################
#Summary
#They are the same then increment UNTRANSCRIBED_STRAND count
#They are the opposite then increment TRANSCRIBED_STRAND count
def searchMutationOnTranscriptionArrays(
        mutation_row,
        chrBased_transcription_plus_array,
        chrBased_transcription_minus_array,
        mutationType2TranscriptionStrand2CountDict,
        mutationType2Sample2TranscriptionStrand2CountDict,
        mutationProbability2Signature2TranscriptionStrand2CountDict,
        mutationProbability2Signature2Sample2TranscriptionStrand2CountDict,
        signatureList,
        mutationProbabilityList):

    mutationStart = mutation_row['Start']
    mutationEnd = mutation_row['End']
    mutationPyramidineStrand = mutation_row['PyramidineStrand']
    mutationType = mutation_row['Mutation']
    mutationSample = mutation_row['Sample']

    #############################################################################################################
    #if there is overlap with chrBased_transcription_plus_array
    if (np.any(chrBased_transcription_plus_array[mutationStart:mutationEnd+1])):

        #########################################################################################################
        # Case1: mutation and transcription strands are on the same strand: increment untranscribed
        #########################################################################################################
        if  (mutationPyramidineStrand==1):
            updateDictionaries(mutation_row,
                                mutationType,
                                mutationSample,
                                mutationType2TranscriptionStrand2CountDict,
                                mutationType2Sample2TranscriptionStrand2CountDict,
                                mutationProbability2Signature2TranscriptionStrand2CountDict,
                                mutationProbability2Signature2Sample2TranscriptionStrand2CountDict,
                                UNTRANSCRIBED_STRAND,
                                signatureList,
                                mutationProbabilityList)

        #########################################################################################################
        # Case2: mutation and transcription strands are the opposite increment transcribed
        #########################################################################################################
        elif (mutationPyramidineStrand==-1):
            updateDictionaries(mutation_row,
                                mutationType,
                                mutationSample,
                                mutationType2TranscriptionStrand2CountDict,
                                mutationType2Sample2TranscriptionStrand2CountDict,
                                mutationProbability2Signature2TranscriptionStrand2CountDict,
                                mutationProbability2Signature2Sample2TranscriptionStrand2CountDict,
                                TRANSCRIBED_STRAND,
                                signatureList,
                                mutationProbabilityList)

    #############################################################################################################


    #############################################################################################################
    # if there is overlap with chrBased_transcription_minus_array
    if (np.any(chrBased_transcription_minus_array[mutationStart:mutationEnd+1])):

        #########################################################################################################
        # Case1: mutation and transcription strands are on the opposite strand: increment transcribed
        #########################################################################################################
        if (mutationPyramidineStrand == 1):
            updateDictionaries(mutation_row,
                               mutationType,
                               mutationSample,
                               mutationType2TranscriptionStrand2CountDict,
                               mutationType2Sample2TranscriptionStrand2CountDict,
                               mutationProbability2Signature2TranscriptionStrand2CountDict,
                               mutationProbability2Signature2Sample2TranscriptionStrand2CountDict,
                               TRANSCRIBED_STRAND,
                               signatureList,
                               mutationProbabilityList)

        #########################################################################################################
        # Case2: mutation and transcription strands are on the same strands increment untranscribed
        #########################################################################################################
        elif (mutationPyramidineStrand == -1):
            updateDictionaries(mutation_row,
                               mutationType,
                               mutationSample,
                               mutationType2TranscriptionStrand2CountDict,
                               mutationType2Sample2TranscriptionStrand2CountDict,
                               mutationProbability2Signature2TranscriptionStrand2CountDict,
                               mutationProbability2Signature2Sample2TranscriptionStrand2CountDict,
                               UNTRANSCRIBED_STRAND,
                               signatureList,
                               mutationProbabilityList)
    #############################################################################################################


########################################################################

########################################################################
def searchMutationsOnTranscriptionArrays(inputList):
    chrBased_spms_split_df = inputList[0]
    chrBased_transcription_on_positive_strand_array = inputList[1]
    chrBased_transcription_on_negative_strand_array = inputList[2]
    signatureList = inputList[3]
    mutationProbabilityList = inputList[4]

    #Initialize empty dictionaries
    mutationType2TranscriptionStrand2CountDict = {}
    mutationType2Sample2TranscriptionStrand2CountDict = {}
    mutationProbability2Signature2TranscriptionStrand2CountDict = {}
    mutationProbability2Signature2Sample2TranscriptionStrand2CountDict = {}

    chrBased_spms_split_df.apply(searchMutationOnTranscriptionArrays,
                            chrBased_transcription_plus_array=chrBased_transcription_on_positive_strand_array,
                            chrBased_transcription_minus_array=chrBased_transcription_on_negative_strand_array,
                            mutationType2TranscriptionStrand2CountDict=mutationType2TranscriptionStrand2CountDict,
                            mutationType2Sample2TranscriptionStrand2CountDict=mutationType2Sample2TranscriptionStrand2CountDict,
                            mutationProbability2Signature2TranscriptionStrand2CountDict=mutationProbability2Signature2TranscriptionStrand2CountDict,
                            mutationProbability2Signature2Sample2TranscriptionStrand2CountDict=mutationProbability2Signature2Sample2TranscriptionStrand2CountDict,
                            signatureList=signatureList,
                            mutationProbabilityList=mutationProbabilityList,
                            axis=1)


    return (mutationType2TranscriptionStrand2CountDict,
            mutationType2Sample2TranscriptionStrand2CountDict,
            mutationProbability2Signature2TranscriptionStrand2CountDict,
            mutationProbability2Signature2Sample2TranscriptionStrand2CountDict)
########################################################################


########################################################################
def fillTranscriptionNPArrayAndSearchMutationsOnThisArray(inputList):

    chrLong = inputList[0]
    chrBased_spms_df = inputList[1]
    chrBased_GRCh37_hg19_NCBI_Curated_RefSeq_Transcripts_df = inputList[2]
    signatureList = inputList[3]
    mutationProbabilityList = inputList[4]

    # int8	Byte (-128 to 127)
    #0 means non-transcribed
    #1 means transcription on positive strand
    #-1 means  transcription on negative strand
    chrBased_transcription_on_positive_strand_array = np.zeros(MAXIMUM_CHROMOSOME_LENGTH, dtype=np.int8)
    chrBased_transcription_on_negative_strand_array = np.zeros(MAXIMUM_CHROMOSOME_LENGTH, dtype=np.int8)

    chrBased_GRCh37_hg19_NCBI_Curated_RefSeq_Transcripts_df.apply(fillTranscriptionArrays,
                                                                    chrBased_transcription_plus_array=chrBased_transcription_on_positive_strand_array,
                                                                    chrBased_transcription_minus_array=chrBased_transcription_on_negative_strand_array,
                                                                    axis=1)
    #Initialize empty dictionaries
    mutationType2TranscriptionStrand2CountDict = {}
    mutationType2Sample2TranscriptionStrand2CountDict = {}
    mutationProbability2Signature2TranscriptionStrand2CountDict = {}
    mutationProbability2Signature2Sample2TranscriptionStrand2CountDict = {}


    chrBased_spms_df.apply(searchMutationOnTranscriptionArrays,
                            chrBased_transcription_plus_array=chrBased_transcription_on_positive_strand_array,
                            chrBased_transcription_minus_array=chrBased_transcription_on_negative_strand_array,
                            mutationType2TranscriptionStrand2CountDict=mutationType2TranscriptionStrand2CountDict,
                            mutationType2Sample2TranscriptionStrand2CountDict=mutationType2Sample2TranscriptionStrand2CountDict,
                            mutationProbability2Signature2TranscriptionStrand2CountDict=mutationProbability2Signature2TranscriptionStrand2CountDict,
                            mutationProbability2Signature2Sample2TranscriptionStrand2CountDict=mutationProbability2Signature2Sample2TranscriptionStrand2CountDict,
                            signatureList=signatureList,
                            mutationProbabilityList=mutationProbabilityList,
                            axis=1)


    return (mutationType2TranscriptionStrand2CountDict,
            mutationType2Sample2TranscriptionStrand2CountDict,
            mutationProbability2Signature2TranscriptionStrand2CountDict,
            mutationProbability2Signature2Sample2TranscriptionStrand2CountDict)
########################################################################



########################################################################
def transcriptionStrandBiasAnalysis(jobname,singlePointMutationsFilename,mutationProbabilityStart,mutationProbabilityEnd,mutationProbabilityStep):
# if __name__ == '__main__':

    print('########################## TranscriptionStrandBias Analysis starts ##########################')
    print('#################### TranscriptionStrandBias Analysis system arguments: #####################')
    print(sys.argv)

    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)

    # jobname = sys.argv[1]
    # singlePointMutationsFilename = sys.argv[2]
    # mutationProbabilityStart  = sys.argv[3]
    # mutationProbabilityEnd = sys.argv[4]
    # mutationProbabilityStep = sys.argv[5]

    ##################### Read GRCh37 HG19 NCBI RefSeq Curated Transcripts starts ######################
    #NCBI has the long chromosome names such as: chr1, chr2, chr3, chr4, ... , chr21, chr22, chrX, chrY, chrMT
    # transcriptsSource = NCBI
    # GRCh37_hg19_Transcripts_df = readTranscriptsNCBI()

    # Let's make SigProfiler use the same transcripts file
    #Ensembl has the short chromosome names such as: 1,2,3,4, ... ,21, 22, X, Y, MT
    transcriptsSource = ENSEMBL
    GRCh37_hg19_Transcripts_df = readTrancriptsENSEMBL()

    # print('GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df[chrom].unique()')
    # print(GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df['chrom'].unique())
    # GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df[chrom].unique()
    # ['chr1' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chrX'
    #  'chrY' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17'
    #  'chr18' 'chr19' 'chr20' 'chr21' 'chr22' 'chr6_apd_hap1' 'chr6_cox_hap2'
    #  'chr6_dbb_hap3' 'chr6_mcf_hap5' 'chr6_qbl_hap6' 'chr4_ctg9_hap1'
    #  'chr6_mann_hap4' 'chr6_ssto_hap7' 'chrUn_gl000211' 'chrUn_gl000212'
    #  'chrUn_gl000213' 'chrUn_gl000218' 'chrUn_gl000219' 'chrUn_gl000220'
    #  'chrUn_gl000224' 'chrUn_gl000228' 'chrUn_gl000241' 'chr17_ctg5_hap1'
    #  'chr1_gl000192_random' 'chr4_gl000193_random' 'chr4_gl000194_random'
    #  'chr7_gl000195_random' 'chr17_gl000205_random']

    # chrNamesInGRCh37 = GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df['chrom'].unique()
    ##################### Read GRCh37 HG19 NCBI RefSeq Curated Transcripts ends ########################


    #################### Prepare mutation Probability List starts ######################################
    mutationProbabilityList = prepareMutationProbabilityList(float(mutationProbabilityStart),float(mutationProbabilityEnd),float(mutationProbabilityStep))
    #################### Prepare mutation Probability List ends ########################################

    transcriptionStrands = [TRANSCRIBED_STRAND, UNTRANSCRIBED_STRAND]
    strandBias = TRANSCRIPTIONSTRANDBIAS

    #Load the signatures
    signatures = []
    SignaturesFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT,jobname,DATA,SignatureFilename)
    if (os.path.exists(SignaturesFilePath)):
        signaturesArray = np.loadtxt(SignaturesFilePath,dtype=str, delimiter='\t')
        signatures = list(signaturesArray)

    #Load the chrnames in single point mutations
    filename = ChrNamesInSPMsFilename
    ChrNamesFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname,DATA,filename)
    if (os.path.exists(ChrNamesFile)):
        chrNamesArray = np.loadtxt(ChrNamesFile,dtype=str, delimiter='\t')
        chrNamesInSPMs = chrNamesArray.tolist()

    # For more information
    accumulatedAllChromosomesMutationProbability2Signature2RatioDict = {}

    #Accumulate chrBased Results
    accumulatedAllChromosomesMutationType2TranscriptStrand2CountDict = {}
    accumulatedAllChromosomesMutationType2Sample2TranscriptStrand2CountDict = {}
    accumulatedAllChromosomesMutationProbability2Signature2TranscriptStrand2CountDict = {}
    accumulatedAllChromosomesMutationProbability2Signature2Sample2TranscriptStrand2CountDict = {}

    # ############################################################################################################
    # #####################################      Version 1 starts      ###########################################
    # ###############################      All Chromosomes in parallel      ######################################
    # ############################################################################################################
    # poolInputList = []
    # ####################################################################################################
    # for chrName in chrNamesInSPMs:
    #     # THEN READ CHRBASED MUTATION
    #     chrLong = 'chr%s' %chrName
    #     chrBased_spms_df = readChrBasedMutationDF(jobname, chrLong, singlePointMutationsFilename)
    #     chrBased_GRCh37_hg19_NCBI_Curated_RefSeq_Transcripts_df = GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df[GRCh37_hg19_NCBI_Curated_RefSeq_Curated_Transcripts_df['chrom'] == chrLong]
    #
    #     if ((chrBased_spms_df is  not None) and (chrBased_GRCh37_hg19_NCBI_Curated_RefSeq_Transcripts_df is  not None)):
    #         inputList = []
    #         inputList.append(chrLong)
    #         inputList.append(chrBased_spms_df)
    #         inputList.append(chrBased_GRCh37_hg19_NCBI_Curated_RefSeq_Transcripts_df)
    #         inputList.append(signatures)  # same for all
    #         inputList.append(mutationProbabilityList)
    #         poolInputList.append(inputList)
    #
    # listofTuples = pool.map(fillTranscriptionNPArrayAndSearchMutationsOnThisArray, poolInputList)
    #
    # accumulate(listofTuples,
    #            accumulatedAllChromosomesMutationType2TranscriptStrand2CountDict,
    #            accumulatedAllChromosomesMutationType2Sample2TranscriptStrand2CountDict,
    #            accumulatedAllChromosomesMutationProbability2Signature2TranscriptStrand2CountDict,
    #            accumulatedAllChromosomesMutationProbability2Signature2Sample2TranscriptStrand2CountDict)
    # ############################################################################################################
    # #####################################      Version 1 ends      #############################################
    # ###############################      All Chromosomes in parallel      ######################################
    # ############################################################################################################



    ############################################################################################################
    #####################################      Version2  starts      ###########################################
    ###############################       Chromosomes sequentially      ########################################
    ###############################      All ChrBased Splits sequentially     ##################################
    ############################################################################################################

    ####################################################################################################
    for chrName in chrNamesInSPMs:
        # THEN READ CHRBASED MUTATION
        chrLong = 'chr%s' %chrName

        #Read chrBased Single Point Mutations
        chrBased_spms_df = readChrBasedMutationDF(jobname, chrLong, singlePointMutationsFilename)

        #Get chrBased ncbi refeq genes
        if (transcriptsSource==NCBI):
            chrBased_GRCh37_hg19_Transcripts_df = GRCh37_hg19_Transcripts_df[GRCh37_hg19_Transcripts_df['chrom'] == chrLong]
        elif (transcriptsSource==ENSEMBL):
            chrBased_GRCh37_hg19_Transcripts_df = GRCh37_hg19_Transcripts_df[GRCh37_hg19_Transcripts_df['chrom'] == chrName]

            if (chrBased_GRCh37_hg19_Transcripts_df is None):
                print('debug August 31 2018 starts')
                print('chrBased_GRCh37_hg19_Transcripts_df is None')
                print('debug August 31 2018 ends')
            elif (chrBased_GRCh37_hg19_Transcripts_df.empty):
                print('debug August 31 2018 starts')
                print('chrBased_GRCh37_hg19_Transcripts_df is empty')
                print('debug August 31 2018 ends')


        #Fill these chrBased arrays once for each chromosome
        # int8	Byte (-128 to 127)
        # 0 means non-transcribed
        # 1 means transcription on positive strand
        # -1 means  transcription on negative strand
        chrBased_transcription_on_positive_strand_array = np.zeros(MAXIMUM_CHROMOSOME_LENGTH, dtype=np.int8)
        chrBased_transcription_on_negative_strand_array = np.zeros(MAXIMUM_CHROMOSOME_LENGTH, dtype=np.int8)



        chrBased_GRCh37_hg19_Transcripts_df.apply(fillTranscriptionArrays,
                                                    chrBased_transcription_plus_array=chrBased_transcription_on_positive_strand_array,
                                                    chrBased_transcription_minus_array=chrBased_transcription_on_negative_strand_array,
                                                    transcriptsSource=transcriptsSource,
                                                    axis=1)



        poolInputList = []
        #split chrBased_spms_df into splits
        if ((chrBased_spms_df is  not None) and (not chrBased_spms_df.empty)):
            chrBasedSPMsDFSplits = np.array_split(chrBased_spms_df, numofProcesses)

            for chrBasedSPMsDFSplit in chrBasedSPMsDFSplits:
                inputList = []
                inputList.append(chrBasedSPMsDFSplit)  # each time different split
                inputList.append(chrBased_transcription_on_positive_strand_array) # same for all
                inputList.append(chrBased_transcription_on_negative_strand_array)  # same for all
                inputList.append(signatures)  # same for all
                inputList.append(mutationProbabilityList)  # same for all
                poolInputList.append(inputList)

            listofTuples = pool.map(searchMutationsOnTranscriptionArrays, poolInputList)

            accumulate(listofTuples,
                       accumulatedAllChromosomesMutationType2TranscriptStrand2CountDict,
                       accumulatedAllChromosomesMutationType2Sample2TranscriptStrand2CountDict,
                       accumulatedAllChromosomesMutationProbability2Signature2TranscriptStrand2CountDict,
                       accumulatedAllChromosomesMutationProbability2Signature2Sample2TranscriptStrand2CountDict)

    ############################################################################################################
    #####################################      Version2  ends      #############################################
    ###############################       Chromosomes sequentially      ########################################
    ###############################      All ChrBased Splits sequentially     ##################################
    ############################################################################################################


    ####################################################################################################


    #debug starts
    # print('Keys are not sorted.')
    print('accumulatedAllChromosomesMutationType2TranscriptStrand2CountDict')
    print(accumulatedAllChromosomesMutationType2TranscriptStrand2CountDict)

    print('accumulatedAllChromosomesMutationType2Sample2TranscriptStrand2CountDict')
    print(accumulatedAllChromosomesMutationType2Sample2TranscriptStrand2CountDict)

    print('accumulatedAllChromosomesMutationProbability2Signature2TranscriptStrand2CountDict')
    print(accumulatedAllChromosomesMutationProbability2Signature2TranscriptStrand2CountDict)

    print('accumulatedAllChromosomesMutationProbability2Signature2Sample2TranscriptStrand2CountDict')
    print(accumulatedAllChromosomesMutationProbability2Signature2Sample2TranscriptStrand2CountDict)
    #debug ends



    #################################################################################################################
    ##########################################      Output starts      ##############################################
    #################################################################################################################
    calculateRatio(accumulatedAllChromosomesMutationProbability2Signature2TranscriptStrand2CountDict,
                   accumulatedAllChromosomesMutationProbability2Signature2RatioDict,transcriptionStrands)


    #############################################################################
    #Highlight the results for mutation probability threshold
    print('###############################################################')
    print('For mutation probability threshold: %f' %MUTATION_SIGNATURE_PROBABILITY_THRESHOLD)
    #header line
    print('Signature\t%s\t%s' %(TRANSCRIBED_STRAND,UNTRANSCRIBED_STRAND))
    signature2TranscriptStrand2CountDict=  accumulatedAllChromosomesMutationProbability2Signature2TranscriptStrand2CountDict[MUTATION_SIGNATURE_PROBABILITY_THRESHOLD]

    for signature in signatures:
        if (signature in signature2TranscriptStrand2CountDict.keys()):
            transcriptStrand2CountDict = signature2TranscriptStrand2CountDict[signature]
            if  ((TRANSCRIBED_STRAND in transcriptStrand2CountDict.keys()) and (UNTRANSCRIBED_STRAND in transcriptStrand2CountDict.keys())):
                print('%s\t%d\t%d' %(signature,transcriptStrand2CountDict[TRANSCRIBED_STRAND],transcriptStrand2CountDict[UNTRANSCRIBED_STRAND]))
    print('TranscriptionStrandBias ratio: number of mutations on transcribed/(transcribed + nontranscribed)')
    print('###############################################################')
    #############################################################################

    #convert
    signature2WeightedAverageRatioDict, signature2StdErrorDict, signature2SumofMutationProbabilitiesDict = convert(accumulatedAllChromosomesMutationProbability2Signature2TranscriptStrand2CountDict,transcriptionStrands)

    #To be used for plotting starts
    writeDictionary(accumulatedAllChromosomesMutationType2TranscriptStrand2CountDict,jobname,MutationType2TranscriptionStrand2CountDict_Filename,strandBias,None)
    writeDictionary(accumulatedAllChromosomesMutationType2Sample2TranscriptStrand2CountDict,jobname,MutationType2Sample2TranscriptionStrand2CountDict_Filename,strandBias,None)
    writeDictionary(accumulatedAllChromosomesMutationProbability2Signature2TranscriptStrand2CountDict,jobname,MutationProbability2Signature2TranscriptionStrand2CountDict_Filename,strandBias,None)
    writeDictionary(accumulatedAllChromosomesMutationProbability2Signature2Sample2TranscriptStrand2CountDict,jobname,MutationProbability2Signature2Sample2TranscriptionStrand2CountDict_Filename,strandBias,None)

    writeDictionary(signature2WeightedAverageRatioDict,jobname,Signature2TranscriptionWeightedAverageRatioDict_Filename,strandBias,None)
    writeDictionary(signature2StdErrorDict, jobname,Signature2TranscriptionStdErrorDict_Filename,strandBias,None)
    writeDictionary(signature2SumofMutationProbabilitiesDict,jobname,Signature2TranscriptionSumofMutationProbabilitiesDict_Filename,strandBias,None)
    #To be used for plotting ends

    #################################################################################################################
    ##########################################      Output ends      ################################################
    #################################################################################################################

    print('########################## TranscriptionStrandBias Analysis ends ############################')

########################################################################
