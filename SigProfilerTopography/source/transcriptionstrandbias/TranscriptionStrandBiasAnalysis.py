# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018 Burcak Otlu

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


from SigProfilerTopography.source.commons.TopographyCommons import *


########################################################################
# Fill chrBased array once for each chromosome
# int8	Byte (-128 to 127)
# We use only one array.
# 0 --> no-transcription
# 1 --> transcription on positive strand
# 2 --> transcription on negative strand
# 3 --> transcription on both positive and negative strands

# Legacy Comments
# Each row is a pandas series in fact
# labels = ['#bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
# using cdsStart and cdsEnd give error since there are intervals of 0 length in that case
# such as cdsStart 33623833 cdsEnd 33623833
#Please notice that same genomic loci on a certain strand can be transcribed and nontranscribed at the same time.
#However same genomic loci on a certain strand can be leading or lagging but not both

# Legacy code
# if (transcriptsSource == NCBI):
#     if (transcription_row['strand'] == PLUS):
#         chrBased_transcription_plus_array[transcription_row['txStart']:(transcription_row['txEnd'] + 1)] = 1
#     elif (transcription_row['strand'] == MINUS):
#         chrBased_transcription_minus_array[transcription_row['txStart']:(transcription_row['txEnd'] + 1)] = -1
def fillTranscriptionArray(transcription_row,chrBased_transcription_array):
    if (transcription_row['strand']==1):
        chrBased_transcription_array[transcription_row['txStart']:(transcription_row['txEnd']+1)] += 1
    elif (transcription_row['strand']==-1):
        chrBased_transcription_array[transcription_row['txStart']:(transcription_row['txEnd']+1)] += 2
########################################################################


########################################################################
# TODO Consider NONTRANSCRIBED_STRAND
def searchMutationUsingTranscriptionStrandColumn(
        mutation_row,
        type2TranscriptionStrand2CountDict,
        sample2Type2TranscriptionStrand2CountDict,
        type2Sample2TranscriptionStrand2CountDict,
        signature2MutationType2TranscriptionStrand2CountDict,
        signature2NumberofMutationsDict,
        mutationProbabilityThreshold,
        type):

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
        updateDictionaries(mutation_row,
                                mutationType,
                                mutationSample,
                                type2TranscriptionStrand2CountDict,
                                sample2Type2TranscriptionStrand2CountDict,
                                type2Sample2TranscriptionStrand2CountDict,
                                signature2MutationType2TranscriptionStrand2CountDict,
                                UNTRANSCRIBED_STRAND,
                                signature2NumberofMutationsDict,
                                mutationProbabilityThreshold)
    elif (mutationTranscriptionStrand == 'T'):
        updateDictionaries(mutation_row,
                                mutationType,
                                mutationSample,
                                type2TranscriptionStrand2CountDict,
                                sample2Type2TranscriptionStrand2CountDict,
                                type2Sample2TranscriptionStrand2CountDict,
                                signature2MutationType2TranscriptionStrand2CountDict,
                                TRANSCRIBED_STRAND,
                                signature2NumberofMutationsDict,
                                mutationProbabilityThreshold)
    elif (mutationTranscriptionStrand == 'B'):
        updateDictionaries(mutation_row,
                                mutationType,
                                mutationSample,
                                type2TranscriptionStrand2CountDict,
                                sample2Type2TranscriptionStrand2CountDict,
                                type2Sample2TranscriptionStrand2CountDict,
                                signature2MutationType2TranscriptionStrand2CountDict,
                                UNTRANSCRIBED_STRAND,
                                signature2NumberofMutationsDict,
                                mutationProbabilityThreshold)
        updateDictionaries(mutation_row,
                                mutationType,
                                mutationSample,
                                type2TranscriptionStrand2CountDict,
                                sample2Type2TranscriptionStrand2CountDict,
                                type2Sample2TranscriptionStrand2CountDict,
                                signature2MutationType2TranscriptionStrand2CountDict,
                                TRANSCRIBED_STRAND,
                                signature2NumberofMutationsDict,
                                mutationProbabilityThreshold)
########################################################################

########################################################################
# TODO Consider NONTRANSCRIBED_STRAND
def searchMutationOnTranscriptionStrandArray(
        mutation_row,
        chrBased_transcription_array,
        type2TranscriptionStrand2CountDict,
        sample2Type2TranscriptionStrand2CountDict,
        type2Sample2TranscriptionStrand2CountDict,
        signature2MutationType2TranscriptionStrand2CountDict,
        signature2NumberofMutationsDict,
        mutationProbabilityThreshold,
        type):

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
    elif(type==DINUCS):
        mutationEnd = mutationStart+2

    # mutationPyramidineStrand= 1 --> pyrimidine mutation is on the + strand
    # mutationPyramidineStrand=-1 --> pyrimidine mutation is on the - strand
    # mutationPyramidineStrand= 0 --> mutation is not all pyrimidine or purine, therefore we can not decide the mutationPyramidineStrand (This is a case for DINUCS and INDELS)

    #Values on chrBased_transcription_array and their meanings
    # 0 --> no-transcription
    # 1 --> transcription on positive strand
    # 2 --> transcription on negative strand
    # 3 --> transcription on both positive and negative strands

    # uniqueIndexesArray = np.unique(chrBased_transcription_array[mutationStart:mutationEnd + 1])
    # mutationEnd is exclusive
    uniqueIndexesArray = np.unique(chrBased_transcription_array[mutationStart:mutationEnd])

    if ((3 in uniqueIndexesArray) or ((1 in uniqueIndexesArray) and (2 in uniqueIndexesArray))):
        if (mutationPyramidineStrand != 0):
            updateDictionaries(mutation_row,
                                    mutationType,
                                    mutationSample,
                                    type2TranscriptionStrand2CountDict,
                                    sample2Type2TranscriptionStrand2CountDict,
                                    type2Sample2TranscriptionStrand2CountDict,
                                    signature2MutationType2TranscriptionStrand2CountDict,
                                    TRANSCRIBED_STRAND,
                                    signature2NumberofMutationsDict,
                                    mutationProbabilityThreshold)

            updateDictionaries(mutation_row,
                                    mutationType,
                                    mutationSample,
                                    type2TranscriptionStrand2CountDict,
                                    sample2Type2TranscriptionStrand2CountDict,
                                    type2Sample2TranscriptionStrand2CountDict,
                                    signature2MutationType2TranscriptionStrand2CountDict,
                                    UNTRANSCRIBED_STRAND,
                                    signature2NumberofMutationsDict,
                                    mutationProbabilityThreshold)

    elif (1 in uniqueIndexesArray):
        if (mutationPyramidineStrand == 1):
            # Transcription is on positive strand, if mutation pyramidine strand is + then increment untranscribed
            updateDictionaries(mutation_row,
                                    mutationType,
                                    mutationSample,
                                    type2TranscriptionStrand2CountDict,
                                    sample2Type2TranscriptionStrand2CountDict,
                                    type2Sample2TranscriptionStrand2CountDict,
                                    signature2MutationType2TranscriptionStrand2CountDict,
                                    UNTRANSCRIBED_STRAND,
                                    signature2NumberofMutationsDict,
                                    mutationProbabilityThreshold)

        elif (mutationPyramidineStrand == -1):
            # Transcription is on positive strand, if mutation pyramidine strand is - then increment transcribed
            updateDictionaries(mutation_row,
                                    mutationType,
                                    mutationSample,
                                    type2TranscriptionStrand2CountDict,
                                    sample2Type2TranscriptionStrand2CountDict,
                                    type2Sample2TranscriptionStrand2CountDict,
                                    signature2MutationType2TranscriptionStrand2CountDict,
                                    TRANSCRIBED_STRAND,
                                    signature2NumberofMutationsDict,
                                    mutationProbabilityThreshold)


    elif (2 in uniqueIndexesArray):
        if (mutationPyramidineStrand == 1):
            # Transcription is on negative strand, if mutation pyramidine strand is + then increment transcribed
            updateDictionaries(mutation_row,
                                    mutationType,
                                    mutationSample,
                                    type2TranscriptionStrand2CountDict,
                                    sample2Type2TranscriptionStrand2CountDict,
                                    type2Sample2TranscriptionStrand2CountDict,
                                    signature2MutationType2TranscriptionStrand2CountDict,
                                    TRANSCRIBED_STRAND,
                                    signature2NumberofMutationsDict,
                                    mutationProbabilityThreshold)

        if (mutationPyramidineStrand == -1):
            # Transcription is on negative strand, if mutation pyramidine strand is - then increment untranscribed
            updateDictionaries(mutation_row,
                                    mutationType,
                                    mutationSample,
                                    type2TranscriptionStrand2CountDict,
                                    sample2Type2TranscriptionStrand2CountDict,
                                    type2Sample2TranscriptionStrand2CountDict,
                                    signature2MutationType2TranscriptionStrand2CountDict,
                                    UNTRANSCRIBED_STRAND,
                                    signature2NumberofMutationsDict,
                                    mutationProbabilityThreshold)
########################################################################


########################################################################
def searchMutations(inputList):
    chrBased_subs_split_df = inputList[0]
    chrBased_indels_split_df = inputList[1]
    chrBased_dinucs_split_df = inputList[2]
    chrBased_transcription_array = inputList[3]
    subsSignature2NumberofMutationsDict = inputList[4]
    indelsSignature2NumberofMutationsDict = inputList[5]
    dinucsSignature2NumberofMutationsDict = inputList[6]

    type2TranscriptionStrand2CountDict= {}
    sample2Type2TranscriptionStrand2CountDict= {}
    type2Sample2TranscriptionStrand2CountDict= {}
    signature2MutationType2TranscriptionStrand2CountDict = {}

    if (chrBased_transcription_array is not None):
        if ((chrBased_subs_split_df is not None) and (not chrBased_subs_split_df.empty)):
            chrBased_subs_split_df.apply(searchMutationOnTranscriptionStrandArray,
                                    chrBased_transcription_array=chrBased_transcription_array,
                                    type2TranscriptionStrand2CountDict = type2TranscriptionStrand2CountDict,
                                    sample2Type2TranscriptionStrand2CountDict =sample2Type2TranscriptionStrand2CountDict,
                                    type2Sample2TranscriptionStrand2CountDict =type2Sample2TranscriptionStrand2CountDict,
                                    signature2MutationType2TranscriptionStrand2CountDict =signature2MutationType2TranscriptionStrand2CountDict,
                                    signature2NumberofMutationsDict=subsSignature2NumberofMutationsDict,
                                    mutationProbabilityThreshold=SUBSTITUTION_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD,
                                    type=SUBS,
                                    axis=1)

        if ((chrBased_indels_split_df is not None) and (not chrBased_indels_split_df.empty)):
            chrBased_indels_split_df.apply(searchMutationOnTranscriptionStrandArray,
                                    chrBased_transcription_array=chrBased_transcription_array,
                                    type2TranscriptionStrand2CountDict = type2TranscriptionStrand2CountDict,
                                    sample2Type2TranscriptionStrand2CountDict =sample2Type2TranscriptionStrand2CountDict,
                                    type2Sample2TranscriptionStrand2CountDict =type2Sample2TranscriptionStrand2CountDict,
                                    signature2MutationType2TranscriptionStrand2CountDict =signature2MutationType2TranscriptionStrand2CountDict,
                                    signature2NumberofMutationsDict=indelsSignature2NumberofMutationsDict,
                                    mutationProbabilityThreshold=INDEL_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD,
                                    type= INDELS,
                                    axis=1)


        if ((chrBased_dinucs_split_df is not None) and (not chrBased_dinucs_split_df.empty)):
            chrBased_dinucs_split_df.apply(searchMutationOnTranscriptionStrandArray,
                                    chrBased_transcription_array=chrBased_transcription_array,
                                    type2TranscriptionStrand2CountDict = type2TranscriptionStrand2CountDict,
                                    sample2Type2TranscriptionStrand2CountDict =sample2Type2TranscriptionStrand2CountDict,
                                    type2Sample2TranscriptionStrand2CountDict=type2Sample2TranscriptionStrand2CountDict,
                                    signature2MutationType2TranscriptionStrand2CountDict =signature2MutationType2TranscriptionStrand2CountDict,
                                    signature2NumberofMutationsDict=dinucsSignature2NumberofMutationsDict,
                                    mutationProbabilityThreshold=DINUC_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD,
                                    type=DINUCS,
                                    axis=1)

    elif (chrBased_transcription_array is None):
        if ((chrBased_subs_split_df is not None) and (not chrBased_subs_split_df.empty)):
            chrBased_subs_split_df.apply(searchMutationUsingTranscriptionStrandColumn,
                                         type2TranscriptionStrand2CountDict=type2TranscriptionStrand2CountDict,
                                         sample2Type2TranscriptionStrand2CountDict=sample2Type2TranscriptionStrand2CountDict,
                                         type2Sample2TranscriptionStrand2CountDict=type2Sample2TranscriptionStrand2CountDict,
                                         signature2MutationType2TranscriptionStrand2CountDict=signature2MutationType2TranscriptionStrand2CountDict,
                                         signature2NumberofMutationsDict=subsSignature2NumberofMutationsDict,
                                         mutationProbabilityThreshold=SUBSTITUTION_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD,
                                         type=SUBS,
                                         axis=1)

        if ((chrBased_indels_split_df is not None) and (not chrBased_indels_split_df.empty)):
            chrBased_indels_split_df.apply(searchMutationUsingTranscriptionStrandColumn,
                                           type2TranscriptionStrand2CountDict=type2TranscriptionStrand2CountDict,
                                           sample2Type2TranscriptionStrand2CountDict=sample2Type2TranscriptionStrand2CountDict,
                                           type2Sample2TranscriptionStrand2CountDict=type2Sample2TranscriptionStrand2CountDict,
                                           signature2MutationType2TranscriptionStrand2CountDict=signature2MutationType2TranscriptionStrand2CountDict,
                                           signature2NumberofMutationsDict=indelsSignature2NumberofMutationsDict,
                                           mutationProbabilityThreshold=INDEL_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD,
                                           type=INDELS,
                                           axis=1)

        if ((chrBased_dinucs_split_df is not None) and (not chrBased_dinucs_split_df.empty)):
            chrBased_dinucs_split_df.apply(searchMutationUsingTranscriptionStrandColumn,
                                           type2TranscriptionStrand2CountDict=type2TranscriptionStrand2CountDict,
                                           sample2Type2TranscriptionStrand2CountDict=sample2Type2TranscriptionStrand2CountDict,
                                           type2Sample2TranscriptionStrand2CountDict=type2Sample2TranscriptionStrand2CountDict,
                                           signature2MutationType2TranscriptionStrand2CountDict=signature2MutationType2TranscriptionStrand2CountDict,
                                           signature2NumberofMutationsDict=dinucsSignature2NumberofMutationsDict,
                                           mutationProbabilityThreshold=DINUC_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD,
                                           type=DINUCS,
                                           axis=1)

    return (type2TranscriptionStrand2CountDict,
            sample2Type2TranscriptionStrand2CountDict,
            type2Sample2TranscriptionStrand2CountDict,
            signature2MutationType2TranscriptionStrand2CountDict)
########################################################################


########################################################################
#main function
def transcriptionStrandBiasAnalysis(mutationTypes,computationType,useTranscriptionStrandColumn,genome,chromSizesDict,chromNamesList,outputDir,jobname):

    print('########################## TranscriptionStrandBias Analysis starts ##########################')
    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)

    ##################### Read Transcripts starts ######################
    #NCBI has the long chromosome names such as: chr1, chr2, chr3, chr4, ... , chr21, chr22, chrX, chrY, chrMT
    # transcriptsSource = NCBI
    # GRCh37_hg19_Transcripts_df = readTranscriptsNCBI()

    # Let's make SigProfiler use the same transcripts file
    #Ensembl has the short chromosome names such as: 1,2,3,4, ... ,21, 22, X, Y, MT
    # transcriptsSource = ENSEMBL
    transcripts_df = readTrancriptsENSEMBL(genome)

    print('transcripts_df.shape')
    print(transcripts_df.shape)
    ##################### Read Transcripts ends ########################

    strandBias = TRANSCRIPTIONSTRANDBIAS

    subsSignature2NumberofMutationsDict = getSubsSignature2NumberofMutationsDict(outputDir,jobname)
    indelsSignature2NumberofMutationsDict = getIndelsSignature2NumberofMutationsDict(outputDir,jobname)
    dinucsSignature2NumberofMutationsDict = getDictionary(outputDir,jobname,DinucsSignature2NumberofMutationsDictFilename)

    #Accumulate chrBased Results
    accumulatedAllChromosomesType2TranscriptionStrand2CountDict = {}
    accumulatedAllChromosomesSample2Type2TranscriptionStrand2CountDict = {}
    accumulatedAllChromosomesType2Sample2TranscriptionStrand2CountDict = {}
    accumulatedAllChromosomesSignature2MutationType2TranscriptionStrand2CountDict = {}

    if (computationType==COMPUTATION_ALL_CHROMOSOMES_PARALLEL):
        ############################################################################################################
        #####################################      Version 1 starts      ###########################################
        ###############################      All Chromosomes in parallel      ######################################
        ############################################################################################################
        poolInputList = []
        ####################################################################################################
        for chrLong in chromNamesList:
            # THEN READ CHRBASED MUTATION
            chromSize = chromSizesDict[chrLong]

            if (SUBS in mutationTypes):
                chrBased_subs_df = readChrBasedSubsDF(outputDir,jobname, chrLong, SUBS)
            if (INDELS in mutationTypes):
                chrBased_indels_df = readChrBasedMutationsDF(outputDir,jobname, chrLong, INDELS)
            if (DINUCS in mutationTypes):
                chrBased_dinucs_df = readChrBasedMutationsDF(outputDir,jobname, chrLong, DINUCS)

            # Get chrBased ensembl transcripts
            #For transcripts_df chrShort is needed
            if (chrLong!='chrM'):
                chrBased_transcripts_df = transcripts_df[transcripts_df['chrom'] == chrLong[3:]]
            elif (chrLong=='chrM'):
                chrBased_transcripts_df = transcripts_df[transcripts_df['chrom'] == 'MT']

            if (useTranscriptionStrandColumn):
                chrBased_transcription_array = None
            else:
                ################################################################################
                chrBased_transcription_array = np.zeros(chromSize, dtype=np.int8)
                chrBased_transcripts_df.apply(fillTranscriptionArray,
                                            chrBased_transcription_array=chrBased_transcription_array,
                                            axis=1)
                ################################################################################

            inputList = []
            inputList.append(chrBased_subs_df)  # each time different split
            inputList.append(chrBased_indels_df)
            inputList.append(chrBased_dinucs_df)
            inputList.append(chrBased_transcription_array)  # same for all
            inputList.append(subsSignature2NumberofMutationsDict)  # same for all
            inputList.append(indelsSignature2NumberofMutationsDict)  # same for all
            inputList.append(dinucsSignature2NumberofMutationsDict)  # same for all
            poolInputList.append(inputList)

        listofTuples = pool.map(searchMutations, poolInputList)

        accumulate(listofTuples,
                   accumulatedAllChromosomesType2TranscriptionStrand2CountDict,
                   accumulatedAllChromosomesSample2Type2TranscriptionStrand2CountDict,
                   accumulatedAllChromosomesType2Sample2TranscriptionStrand2CountDict,
                   accumulatedAllChromosomesSignature2MutationType2TranscriptionStrand2CountDict)

        ############################################################################################################
        #####################################      Version 1 ends      #############################################
        ###############################      All Chromosomes in parallel      ######################################
        ############################################################################################################

    elif (computationType==COMPUTATION_CHROMOSOMES_SEQUENTIAL_CHROMOSOME_SPLITS_PARALLEL):

        ############################################################################################################
        #####################################      Version2  starts      ###########################################
        ###############################       Chromosomes sequentially      ########################################
        ###############################      All ChrBased Splits sequentially     ##################################
        ############################################################################################################

        ####################################################################################################
        for chrLong in chromNamesList:

            chromSize = chromSizesDict[chrLong]

            if (SUBS in mutationTypes):
                chrBased_subs_df = readChrBasedSubsDF(outputDir, jobname, chrLong, SUBS)
            if (INDELS in mutationTypes):
                chrBased_indels_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, INDELS)
            if (DINUCS in mutationTypes):
                chrBased_dinucs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong,DINUCS)

            # Get chrBased ensembl transcripts
            chrBased_transcripts_df = transcripts_df[transcripts_df['chrom'] == chrLong]

            #You need to initialize to None so that you don't use former for loop values accidentally
            chrBased_subs_df_splits_list = None
            chrBased_indels_df_splits_list = None
            chrBased_dinucs_df_splits_list = None

            if ((chrBased_subs_df is not None) and (not chrBased_subs_df.empty)):
                chrBased_subs_df_splits_list = np.array_split(chrBased_subs_df, numofProcesses)

            if ((chrBased_indels_df is not None) and (not chrBased_indels_df.empty)):
                chrBased_indels_df_splits_list = np.array_split(chrBased_indels_df, numofProcesses)

            if ((chrBased_dinucs_df is not None) and (not chrBased_dinucs_df.empty)):
                chrBased_dinucs_df_splits_list = np.array_split(chrBased_dinucs_df, numofProcesses)


            if (useTranscriptionStrandColumn):
                chrBased_transcription_array = None
            else:
                ################################################################################
                chrBased_transcription_array = np.zeros(chromSize, dtype=np.int8)
                chrBased_transcripts_df.apply(fillTranscriptionArray,
                                            chrBased_transcription_array=chrBased_transcription_array,
                                            axis=1)
                ################################################################################


            ####################################################################
            poolInputList = []
            for split_index in range(numofProcesses):
                chrBased_subs_split_df = None
                chrBased_indels_split_df = None
                if ((chrBased_subs_df_splits_list is not None) and (len(chrBased_subs_df_splits_list))):
                    chrBased_subs_split_df = chrBased_subs_df_splits_list[split_index]
                if ((chrBased_indels_df_splits_list is not None) and (len(chrBased_indels_df_splits_list))):
                    chrBased_indels_split_df = chrBased_indels_df_splits_list[split_index]
                if ((chrBased_dinucs_df_splits_list is not None) and (len(chrBased_dinucs_df_splits_list))):
                    chrBased_dinucs_split_df = chrBased_dinucs_df_splits_list[split_index]


                inputList = []
                inputList.append(chrBased_subs_split_df)  # each time different split
                inputList.append(chrBased_indels_split_df)
                inputList.append(chrBased_dinucs_split_df)
                inputList.append(chrBased_transcription_array) # same for all
                inputList.append(subsSignature2NumberofMutationsDict)  # same for all
                inputList.append(indelsSignature2NumberofMutationsDict)  # same for all
                inputList.append(dinucsSignature2NumberofMutationsDict)  # same for all
                poolInputList.append(inputList)
            ####################################################################

            listofTuples = pool.map(searchMutations,poolInputList)

            accumulate(listofTuples,
                           accumulatedAllChromosomesType2TranscriptionStrand2CountDict,
                           accumulatedAllChromosomesSample2Type2TranscriptionStrand2CountDict,
                           accumulatedAllChromosomesType2Sample2TranscriptionStrand2CountDict,
                           accumulatedAllChromosomesSignature2MutationType2TranscriptionStrand2CountDict)
        ############################################################################################################
        #####################################      Version2  ends      #############################################
        ###############################       Chromosomes sequentially      ########################################
        ###############################      All ChrBased Splits sequentially     ##################################
        ############################################################################################################

    #################################################################################################################
    ##########################################      Output starts      ##############################################
    #################################################################################################################
    #############################################################################
    writeDictionary(accumulatedAllChromosomesType2TranscriptionStrand2CountDict,outputDir,jobname,Type2TranscriptionStrand2CountDict_Filename,strandBias,None)
    writeDictionary(accumulatedAllChromosomesSample2Type2TranscriptionStrand2CountDict,outputDir,jobname,Sample2Type2TranscriptionStrand2CountDict_Filename,strandBias,None)
    writeDictionary(accumulatedAllChromosomesType2Sample2TranscriptionStrand2CountDict, outputDir, jobname,Type2Sample2TranscriptionStrand2CountDict_Filename, strandBias, None)
    writeDictionary(accumulatedAllChromosomesSignature2MutationType2TranscriptionStrand2CountDict,outputDir,jobname,Signature2MutationType2TranscriptionStrand2CountDict_Filename,strandBias,None)
    #################################################################################################################
    ##########################################      Output ends      ################################################
    #################################################################################################################
    print('########################## TranscriptionStrandBias Analysis ends ############################')

########################################################################
