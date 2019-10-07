# This code contains the common functions that are used for Data Preparation

import os
import sys
import twobitreader

complementDict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

ONE_DIRECTORY_UP = '..'
INPUT = 'input'
OUTPUT = 'output'
LIB = 'lib'
UCSCGENOME = 'ucscgenome'
NOTSET= 'notset'

current_abs_path = os.path.abspath(os.path.dirname(__file__))
sys.path.append(current_abs_path)

from TopographyCommons import *

############################################################
def readProbabilities(probabilitiesFile,logger):

    probabilities_df = pd.read_table(probabilitiesFile,sep="\t", header=0)

    logger.info('Probabilities information starts')
    logger.info('probabilities_df.shape')
    logger.info(probabilities_df.shape)
    logger.info('probabilities_df.head()')
    logger.info(probabilities_df.head())

    if ('Sample' in probabilities_df.columns.values):
        logger.info('Unique samples in probabilities_df')
        logger.info(probabilities_df['Sample'].unique())
        logger.info('# of unique samples in probabilities_df: %d\n' %(len(probabilities_df['Sample'].unique())))

    if ('Mutations' in probabilities_df.columns.values):
        logger.info('Unique mutations in probabilities_df')
        logger.info(probabilities_df['Mutations'].unique())
        logger.info('# of unique mutations in probabilities_df: %d\n' %(len(probabilities_df['Mutations'].unique())))

    if ('MutationTypes' in probabilities_df.columns.values):
        logger.info('Unique MutationTypes in probabilities_df')
        logger.info(probabilities_df['MutationTypes'].unique())
        logger.info('# of unique MutationTypes in probabilities_df: %d' %(len(probabilities_df['MutationTypes'].unique())))
    logger.info('Probabilities information ends')
    logger.info('##############################')

    return probabilities_df
############################################################



############################################################
# Notice that [::-1] provides visiting x from the last base to the first base
revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N'}[B] for B in x][::-1])
############################################################


# ############################################################
# def complement(upperNucleotides):
#     complemented = ''
#     for upperNucleotide in upperNucleotides:
#         complemented += complementDict[upperNucleotide]
#     return complemented
# ############################################################


############################################################
def getMutationPyramidineStrandInformation(row,hg19_genome,hg38_genome):
    originalNucleotide = row[REF]
    mutatedNucleotide = row[ALT]

    if (row[GENOME] == GRCh37 or row[GENOME]==HG19):
        threeNucleotidesUpper = getNucleotides(row[CHROM], row[START]-2, row[START]+1, hg19_genome)
        oneNucleotide = getNucleotides(row[CHROM], row[START]-1, row[START], hg19_genome)

    elif (row[GENOME] == GRCh38 or row[GENOME]==HG38):
        threeNucleotidesUpper = getNucleotides(row[CHROM], row[START]-2, row[START]+1, hg38_genome)
        oneNucleotide = getNucleotides(row[CHROM], row[START]-1, row[START], hg38_genome)

    if (oneNucleotide==originalNucleotide):
        strand = +1
    else:
        strand = -1

    if (originalNucleotide=='G' or originalNucleotide=='A' or originalNucleotide=='g' or originalNucleotide=='a'):
        originalNucleotide = revcompl(originalNucleotide.upper())
        mutatedNucleotide = revcompl(mutatedNucleotide.upper())
        threeNucleotidesUpper = revcompl(threeNucleotidesUpper)
        pyramidineStrand = -1*strand
    else:
        pyramidineStrand = strand

    mutations = '%s[%s>%s]%s' %(threeNucleotidesUpper[0],originalNucleotide,mutatedNucleotide,threeNucleotidesUpper[2])

    row[MUTATIONS] = mutations
    row[PYRAMIDINESTRAND] = pyramidineStrand
    row[MUTATION] = '%s>%s' %(originalNucleotide,mutatedNucleotide)
    row[CONTEXT] = threeNucleotidesUpper

    return row
############################################################

############################################################
def getMutationInformation(row,hg19_genome, hg38_genome):
    originalNucleotide = row[REF]
    mutatedNucleotide = row[ALT]
    strand = convertStrandIntoNum(row[STRAND])

    if row[GENOME] == 'GRCh37':
        threeNucleotidesUpper = getNucleotides(row[CHROM], row[START]-2, row[START]+1, hg19_genome)
    elif row[GENOME] == 'GRCh38':
        threeNucleotidesUpper = getNucleotides(row[CHROM], row[START]-2, row[START]+1, hg38_genome)

    if (originalNucleotide=='G' or originalNucleotide=='A' or originalNucleotide=='g' or originalNucleotide=='a'):
        originalNucleotide = revcompl(originalNucleotide.upper())
        mutatedNucleotide = revcompl(mutatedNucleotide.upper())
        threeNucleotidesUpper = revcompl(threeNucleotidesUpper)
        pyramidineStrand = -1* strand
    else:
        pyramidineStrand = strand

    mutations = '%s[%s>%s]%s' %(threeNucleotidesUpper[0],originalNucleotide,mutatedNucleotide,threeNucleotidesUpper[2])

    row[MUTATIONS] = mutations
    row[PYRAMIDINESTRAND] = pyramidineStrand
    row[MUTATION] = '%s>%s' %(originalNucleotide,mutatedNucleotide)
    row[CONTEXT] = threeNucleotidesUpper

    return row
############################################################

############################################################
#This code is customized w.r.t. PCAWG Matlab Probabilities
#TODO However during release we must be able to use MUTATIONLONG
# example line for original data
#UCEC-US_SP89389 10      2017540 N:AT[T>A]CA     1
# example line for simulated data
#UCEC-US_SP89389_1       10      1575080 T:AT[C>T]TG     1
def readChrBasedMutations(chr_based_mutation_filepath,mutation_type_context):

    if (os.path.exists(chr_based_mutation_filepath)):
        try:
            mutations_with_genomic_positions_df = pd.read_table(chr_based_mutation_filepath, sep="\t", header=None)
        except pd.errors.EmptyDataError:
            mutations_with_genomic_positions_df = pd.DataFrame()

        if (not mutations_with_genomic_positions_df.empty):
            if (mutation_type_context==DBS):
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
                mutations_with_genomic_positions_df.columns = [SAMPLE,CHROM,START,MUTATIONLONG,PYRAMIDINESTRAND]
                mutations_with_genomic_positions_df[SAMPLE] = mutations_with_genomic_positions_df[SAMPLE].astype(str)
                mutations_with_genomic_positions_df[CHROM] = mutations_with_genomic_positions_df[CHROM].astype(str)
                mutations_with_genomic_positions_df[START] = mutations_with_genomic_positions_df[START].astype(int)
                mutations_with_genomic_positions_df[MUTATIONLONG] = mutations_with_genomic_positions_df[MUTATIONLONG].astype(str)
                mutations_with_genomic_positions_df[PYRAMIDINESTRAND] = mutations_with_genomic_positions_df[PYRAMIDINESTRAND].astype(int)
                #Add new column
                # Add Context Column from T:TG[C>T]GC to  G[C>T]G
                # MatriXGenerator generates T:TG[C>T]GC
                # PCAWG_Matlab sbs probabilities has G[C>T]G
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

    chr_based_mutation_df = readChrBasedMutations(chr_based_mutation_filepath,mutation_type_context)

    if ((chr_based_mutation_df is not None) and (mutations_probabilities_df is not None)):

        ############################################################################
        #For PCAWG_Matlab
        #Convert CMDI-UK_SP116871_1 --> SP116871 # if(simNum>0): simNum=1 Simulation1
        #Convert CMDI-UK_SP116871 --> SP116871 # simNum=0 Original Data
        # chr_based_mutation_df[SAMPLE] = chr_based_mutation_df[SAMPLE].str.split('_', expand=True)[1]

        #For release
        if simNum>=1:
            chr_based_mutation_df[SAMPLE] = chr_based_mutation_df[SAMPLE].str.rsplit('_', 1, expand=True)[0]
        ############################################################################

        ############################################################################
        #For PCAWG Matlab Probabilities
        # merged_df = pd.merge(chr_based_mutation_df,mutations_probabilities_df, how='inner', left_on=[SAMPLE, MUTATION],right_on=[SAMPLE, MUTATION])

        # For release
        merged_df = pd.merge(chr_based_mutation_df,mutations_probabilities_df, how='inner', left_on=[SAMPLE, MUTATIONLONG],right_on=[SAMPLE, MUTATIONLONG])
        ############################################################################

        if ((merged_df is not None) and (chr_based_mutation_df.shape[0]!=merged_df.shape[0])):
            print('There is a situation. For simNum:%s chr:%s All %s mutations are not merged with signature probabilities'  %(simNum,chrShort,mutation_type_context))
            print('chr_based_mutation_df.shape for simNum:%s chr:%s %s' % (simNum, chrShort, mutation_type_context))
            print(chr_based_mutation_df.shape)
            print('merged_df.shape for simNum:%s chr:%s %s' % (simNum, chrShort, mutation_type_context))
            print(merged_df.shape)

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
def readMutationsWithGenomicPositions(snpsInputFile):
    mutationsWithGenomicPositionFilePath = os.path.join(snpsInputFile)

    # column_names_18 = ['Project', 'Sample', 'locID', 'Genome', 'mutType', 'Chrom', 'Start', 'End', 'Ref', 'Alt', 'Type', 'VarID', 'Strand', 'Gene', 'GeneID','ccdsID', 'TranscriptID',	'GeneType']
    # column_names_11 = ['Project', 'Sample', 'locID', 'Genome', 'mutType', 'Chrom', 'Start', 'End', 'Ref', 'Alt', 'Type']

    # Project Sample  locID   Genome  mutType Chrom Start   End     Ref     Alt     Type    VarID   Strand  Gene    GeneID  ccdsID  TranscriptID    GeneType
    # 21_BRCA_WGS     PD3851a .       GRCh37  SNV     1       809687  809688  G       C       SOMATIC 27610826
    #
    # Project Sample  locID   Genome  mutType Chrom   Start   End     Ref     Alt     Type
    #BRCA    PD10010a        CGP     GRCh37  SNP     X       4643309 4643309 G       A       SOMATIC

    # snps_with_genomic_positions_df = pd.read_table(mutationsWithGenomicPositionFilePath, sep="\t", header=0, names=column_names,dtype={'Chrom': str, 'Sample':str})
    snps_with_genomic_positions_df = pd.read_table(mutationsWithGenomicPositionFilePath, sep="\t", header=0, dtype={CHROM: str, SAMPLE: str})

    print('###########################################')
    print('SNPs information starts')
    print('snps_with_genomic_positions_df.head()')
    print(snps_with_genomic_positions_df.head())

    print('snps_with_genomic_positions_df.columns.values')
    print(snps_with_genomic_positions_df.columns.values)

    print('snps_with_genomic_positions_df.shape')
    print(snps_with_genomic_positions_df.shape)

    print('# of rows in snps_with_genomic_positions_df: %d' % len(snps_with_genomic_positions_df))

    # How many projects (cancer type)  are there?
    print('Unique projects (cancer type) in the snpsInputFile:')
    print(snps_with_genomic_positions_df[PROJECT].unique())
    print('# of projects in the snpsInputFile: %d' % len(snps_with_genomic_positions_df[PROJECT].unique()))
    print()

    #How many samples are there?
    print('Unique sample names in the snpsInputFile:')
    print(snps_with_genomic_positions_df[SAMPLE].unique())
    print('# of samples in the snpsInputFile: %d' %len(snps_with_genomic_positions_df[SAMPLE].unique()))
    print()

    print('Unique chroms in the snpsInputFile:')
    print(snps_with_genomic_positions_df[CHROM].unique())
    print('# of chroms in the snpsInputFile: %d' % len(snps_with_genomic_positions_df[CHROM].unique()))
    print()

    print('Unique genomes in the snpsInputFile:')
    print(snps_with_genomic_positions_df[GENOME].unique())
    print('# of genomes in the snpsInputFile: %d' % len(snps_with_genomic_positions_df[GENOME].unique()))
    print()

    print('Unique types in the snpsInputFile:')
    print(snps_with_genomic_positions_df[TYPE].unique())
    print('# of types in the snpsInputFile: %d' % len(snps_with_genomic_positions_df[TYPE].unique()))
    print()

    print('SNPs information ends')
    print('###########################################')

    return snps_with_genomic_positions_df
############################################################



############################################################
def convertStrandIntoNum(strand):
    if strand=='-':
        return -1
    else:
        return 1
############################################################

