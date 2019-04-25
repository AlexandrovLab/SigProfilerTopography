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
def readProbabilities(probabilitiesFile):

    probabilities_df = pd.read_table(probabilitiesFile,sep="\t", header=0)

    print('###########################################')
    print('Probabilities information starts')

    print('probabilities_df.shape')
    print(probabilities_df.shape)
    print('probabilities_df.head()')
    print(probabilities_df.head())

    if ('Sample' in probabilities_df.columns.values):
        print('Unique samples in probabilities_df')
        print(probabilities_df['Sample'].unique())
        print('# of unique samples in probabilities_df: %d' %(len(probabilities_df['Sample'].unique())))
        print()

    if ('Mutations' in probabilities_df.columns.values):
        print('Unique mutations in probabilities_df')
        print(probabilities_df['Mutations'].unique())
        print('# of unique mutations in probabilities_df: %d' %(len(probabilities_df['Mutations'].unique())))
        print()

    if ('MutationTypes' in probabilities_df.columns.values):
        print('Unique MutationTypes in probabilities_df')
        print(probabilities_df['MutationTypes'].unique())
        print('# of unique MutationTypes in probabilities_df: %d' %(len(probabilities_df['MutationTypes'].unique())))
        print()

    print('Probabilities information ends')
    print('###########################################')

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
def readChrBasedMutations(chr_based_mutation_filepath,mutationType):

    if (os.path.exists(chr_based_mutation_filepath)):
        try:
            mutations_with_genomic_positions_df = pd.read_table(chr_based_mutation_filepath, sep="\t", header=None)
        except pd.errors.EmptyDataError:
            mutations_with_genomic_positions_df = pd.DataFrame()

        if (not mutations_with_genomic_positions_df.empty):
            if (mutationType==DINUCS):
                mutations_with_genomic_positions_df.columns = [SAMPLE,CHROM,START,MUTATIONLONG,PYRAMIDINESTRAND]
                mutations_with_genomic_positions_df[SAMPLE] = mutations_with_genomic_positions_df[SAMPLE].astype(str)
                mutations_with_genomic_positions_df[CHROM] = mutations_with_genomic_positions_df[CHROM].astype(str)
                mutations_with_genomic_positions_df[START] = mutations_with_genomic_positions_df[START].astype(int)
                mutations_with_genomic_positions_df[MUTATIONLONG] = mutations_with_genomic_positions_df[MUTATIONLONG].astype(str)
                mutations_with_genomic_positions_df[PYRAMIDINESTRAND] = mutations_with_genomic_positions_df[PYRAMIDINESTRAND].astype(int)
                #Add new columns
                mutations_with_genomic_positions_df[TRANSCRIPTIONSTRAND] = mutations_with_genomic_positions_df[MUTATIONLONG].str[0]
                mutations_with_genomic_positions_df[MUTATION] = mutations_with_genomic_positions_df[MUTATIONLONG].str[4:9]
            elif (mutationType==INDELS):
                mutations_with_genomic_positions_df.columns = [SAMPLE,CHROM,START,MUTATIONLONG,REF,ALT,PYRAMIDINESTRAND]
                mutations_with_genomic_positions_df[SAMPLE] = mutations_with_genomic_positions_df[SAMPLE].astype(str)
                mutations_with_genomic_positions_df[CHROM] = mutations_with_genomic_positions_df[CHROM].astype(str)
                mutations_with_genomic_positions_df[START] = mutations_with_genomic_positions_df[START].astype(int)
                mutations_with_genomic_positions_df[MUTATIONLONG] = mutations_with_genomic_positions_df[MUTATIONLONG].astype(str)
                mutations_with_genomic_positions_df[PYRAMIDINESTRAND] = mutations_with_genomic_positions_df[PYRAMIDINESTRAND].astype(int)
                #Add new column
                mutations_with_genomic_positions_df[LENGTH] = mutations_with_genomic_positions_df[REF].apply(len)
                mutations_with_genomic_positions_df[TRANSCRIPTIONSTRAND] = mutations_with_genomic_positions_df[MUTATIONLONG].str[0]
                mutations_with_genomic_positions_df[MUTATION] = mutations_with_genomic_positions_df[MUTATIONLONG].str[2:]
                #order the columns make CONTEXT at the end.
                ordered_column_names = [SAMPLE,CHROM,START,MUTATIONLONG,REF,ALT,LENGTH,PYRAMIDINESTRAND,TRANSCRIPTIONSTRAND,MUTATION]
                mutations_with_genomic_positions_df = mutations_with_genomic_positions_df[ordered_column_names]
            elif(mutationType==SUBS):
                mutations_with_genomic_positions_df.columns = [SAMPLE,CHROM,START,MUTATIONLONG,PYRAMIDINESTRAND]
                mutations_with_genomic_positions_df[SAMPLE] = mutations_with_genomic_positions_df[SAMPLE].astype(str)
                mutations_with_genomic_positions_df[CHROM] = mutations_with_genomic_positions_df[CHROM].astype(str)
                mutations_with_genomic_positions_df[START] = mutations_with_genomic_positions_df[START].astype(int)
                mutations_with_genomic_positions_df[MUTATIONLONG] = mutations_with_genomic_positions_df[MUTATIONLONG].astype(str)
                mutations_with_genomic_positions_df[PYRAMIDINESTRAND] = mutations_with_genomic_positions_df[PYRAMIDINESTRAND].astype(int)
                #Add new column
                # Add Context Column from T:TG[C>T]GC to GCC
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
    mutationType = inputList[5]

    chr_based_mutation_df = readChrBasedMutations(chr_based_mutation_filepath,mutationType)

    #chr_based_dinuc_df columns ['Sample', 'Chrom', 'Start', 'Mutation','PyrimidineStrand']
    #dinucs_probabilities_df columns ['Sample', 'Mutation', 'DBS2', 'DBS4', 'DBS6', 'DBS7', 'DBS11']
    if ((chr_based_mutation_df is not None) and (mutations_probabilities_df is not None)):
        merged_df = pd.merge(chr_based_mutation_df,mutations_probabilities_df, how='inner', left_on=[SAMPLE, MUTATION],right_on=[SAMPLE, MUTATION])

        # print('chr_based_dinuc_df.shape')
        # print(chr_based_dinuc_df.shape[0])
        # print('dinucs_probabilities_df.shape')
        # print(dinucs_probabilities_df.shape)
        # print('merged_df.shape')
        # print(merged_df.shape[0])

        if (chr_based_mutation_df.shape[0]!=merged_df.shape[0]):
            print('There is a situation. For chr:%s All dinucs are not merged with dinucs signature probabilities' %(chrShort))

        if ((merged_df is not None) and (not merged_df.empty)):
            chrBasedMutationsFileName = 'chr%s_%s_for_topography.txt' %(chrShort,mutationType)
            chr_based_mutations_file_path = os.path.join(outputDir, jobname, DATA, CHRBASED, chrBasedMutationsFileName)

            # #After test uncomment
            # if ('MutationLong' in merged_df.columns.values):
            #     merged_df.drop(['MutationLong'], inplace=True, axis=1)

            #After merge
            if (mutationType==SUBS):
                merged_df[MUTATION] = merged_df[MUTATION].str[2:5]

            merged_df.to_csv(chr_based_mutations_file_path, sep='\t', header=True, index=False)
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



############################################################
def combineGenomicPositionsWithSignatureProbabilitiesParallel(chrBased_cancerTypeBased_Genomic_Positions_InputList):
    cancerType = chrBased_cancerTypeBased_Genomic_Positions_InputList[0]
    chrNumber = chrBased_cancerTypeBased_Genomic_Positions_InputList[1]
    chrBased_snps_df = chrBased_cancerTypeBased_Genomic_Positions_InputList[2]

    probabilitiles_df = chrBased_cancerTypeBased_Genomic_Positions_InputList[3]
    hg19_genome = chrBased_cancerTypeBased_Genomic_Positions_InputList[4]
    hg38_genome = chrBased_cancerTypeBased_Genomic_Positions_InputList[5]

    #####################################################
    # What apply method returns depend on the the return type in the apply method.
    # For example, if you return one value in apply, after your call, apply will return series.
    # If you return values in apply, after your call, apply will return dataframe.
    if STRAND in chrBased_snps_df.columns:
        chrBased_cancerTypeBased_df = chrBased_snps_df.apply(getMutationInformation,hg19_genome=hg19_genome,hg38_genome=hg38_genome,axis=1)
    else:
        chrBased_cancerTypeBased_df = chrBased_snps_df.apply(getMutationPyramidineStrandInformation,hg19_genome=hg19_genome,hg38_genome=hg38_genome,axis=1)


    print('For chr%s ' %chrNumber)
    result_df = pd.merge(chrBased_cancerTypeBased_df, probabilitiles_df, how='inner', left_on=[SAMPLE,MUTATIONS], right_on=[SAMPLE,MUTATIONS])

    #Drop the unnecessary columns after getMutationInformation
    result_df.drop([PROJECT, GENOME, REF, ALT, STRAND], inplace=True, axis=1,errors='ignore')

    columnNames = list(result_df.columns.values)
    indexofContext = columnNames.index(CONTEXT)
    signatureColumnNames = columnNames[indexofContext+1:]

    # Column order we want to have:
    ordered_column_names = [SAMPLE, CHROM, START, END, PYRAMIDINESTRAND, MUTATION, CONTEXT]
    ordered_column_names.extend(signatureColumnNames)

    result_df = result_df[ordered_column_names]
    #####################################################

    return result_df
############################################################

################################################################
#Read a sample file for a certain simulationNumber
#Returns snp_df and indel_df for each read sample file
def readSampleBasedSimulationBasedDF(sampleBasedSimulationsInputPath,sampleName,simulationNumber,mutationTypes):
    # sampleFileName =PD10010a_INDEL_96_10.txt

    if (len(mutationTypes)==2):
        sampleFileName = '%s_%s_%s_%d.txt' %(sampleName,mutationTypes[1],mutationTypes[0],simulationNumber)
    elif (len(mutationTypes)==1):
        sampleFileName = '%s_%s_%d.txt' %(sampleName,mutationTypes[0],simulationNumber)

    sampleBasedSimulationFilePath = os.path.join(sampleBasedSimulationsInputPath,sampleName,sampleFileName)

    # print('##############################################################')
    # print('sampleName:%s' % sampleName)
    #X       78540063        PD3851a T       A       .       Simulations     GRCh37  TTC     -1
    #19      33823827        PD3851a ACC     A       .       Simulations     GRCh37  2:Del:R:1       +1
    # X       112956612       PD3851a AT      A       .       Simulations     GRCh37  1:Del:T:0       +1

    #Please notice that later on I use column numbers to drop columns
    #column_names = ['chrNumber', 'start', 'sample', 'originalNucleotide', 'mutatedNucleotide','dot' ,'dataSource','genomeAssembly', 'context', 'strand']
    # dtypes = {'chrNumber':str,'start':int,'sample':str,'originalNucleotide':str,'mutatedNucleotide':str,'dot':str,'dataSource':str,'genomeAssembly':str,'context':str,'strand':int}
    dtypes = {0:str,1:int,2:str,3:str,4:str,5:str,6:str,7:str,8:str,9:int}

    # cancerBased_sampleBased_mutations_with_genomic_positions_df = pd.read_table(sampleBasedSimulationFilePath, sep="\t",header=None, names= column_names, dtype=dtypes)
    cancerBased_sampleBased_mutations_with_genomic_positions_df = pd.read_table(sampleBasedSimulationFilePath, sep="\t",header=None,dtype=dtypes)

    #If originalNucleotide and mutatedNucleotide have length of 1 then snp
    sampleBased_simulationBased_snp_df = cancerBased_sampleBased_mutations_with_genomic_positions_df[(cancerBased_sampleBased_mutations_with_genomic_positions_df.iloc[:,3].str.len()==1) & (cancerBased_sampleBased_mutations_with_genomic_positions_df.iloc[:,4].str.len()==1)]
    #If originalNucleotide or mutatedNucleotide have length not equal to 1 then indel
    sampleBased_simulationBased_indel_df = cancerBased_sampleBased_mutations_with_genomic_positions_df[(cancerBased_sampleBased_mutations_with_genomic_positions_df.iloc[:,3].str.len()!=1) | (cancerBased_sampleBased_mutations_with_genomic_positions_df.iloc[:,4].str.len()!=1)]

    #debug starts
    print('sampleBased_simulationBased_snp_df.shape')
    print(sampleBased_simulationBased_snp_df.shape)
    print('sampleBased_simulationBased_snp_df.columns')
    print(sampleBased_simulationBased_snp_df.columns)

    print('sampleBased_simulationBased_indel_df.shape')
    print(sampleBased_simulationBased_indel_df.shape)
    print('sampleBased_simulationBased_indel_df.columns')
    print(sampleBased_simulationBased_indel_df.columns)
    # #debug ends

    return sampleBased_simulationBased_snp_df,sampleBased_simulationBased_indel_df
################################################################


################################################################
def prepareSimulationBasedInputFilesForSigProfilerTopography(jobname,genome,mutationTypes,sigProfilerSimulatorSpecificDirName,numberofSimulations,probabilities_df):

    current_abs_path = os.path.dirname(os.path.realpath(__file__))

    # /oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/input/21BreastCancer/output/simulations/21BreastCancer_simulations_GRCh37_INDEL_96/
    sampleBasedSimulationsInputPath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, INPUT, jobname, OUTPUT, SIMULATIONS, sigProfilerSimulatorSpecificDirName)

    #/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/input/21BreastCancer/SIMULATIONS_FOR_TOPOGRAPHY
    output_path = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,INPUT,jobname,SIMULATIONS_FOR_TOPOGRAPHY)

    #read hg19_genome
    if (genome==HG19 or genome==GRCh37 ):
        human_genome = twobitreader.TwoBitFile(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, UCSCGENOME, HG19_2BIT))
    elif (genome==HG38 or genome==GRCh38):
        human_genome = twobitreader.TwoBitFile(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, UCSCGENOME, HG38_2BIT))

    #Let's not remove the existing prepared simulation data for SigProfilerTopography
    # if os.path.exists(output_path):
    #     os.system("rm -rf " + output_path)

    os.makedirs(output_path, exist_ok=True)

    # get the samples
    sampleDirNameList = os.listdir(sampleBasedSimulationsInputPath)

    #After simulations, we have sample based and simulation based files
    #We need to pool all the samples and convert these files into simulation based all samples together files
    for simNumber in range(1,numberofSimulations+1):
        print('###########################################################')
        print('Simulation Number: %s' %simNumber)

        simulationBased_allSamples_snp_filename = '%s_Sim%d_SNP.txt' %(jobname,simNumber)
        simulationBased_allSamples_indel_filename = '%s_Sim%d_INDEL.txt' %(jobname,simNumber)

        simulationBased_allSamples_snp_filepath = os.path.join(output_path,simulationBased_allSamples_snp_filename)
        simulationBased_allSamples_indel_filepath = os.path.join(output_path,simulationBased_allSamples_indel_filename)

        sampleBased_snp_df_list = []
        sampleBased_indel_df_list = []

        for sampleIndex, sampleName in enumerate(sampleDirNameList,1):
            # print('SampleIndex:%d' %(sampleIndex))
            sampleBased_simulationBased_snp_df, sampleBased_simulationBased_indel_df = readSampleBasedSimulationBasedDF(sampleBasedSimulationsInputPath,sampleName,simNumber,mutationTypes)
            sampleBased_snp_df_list.append(sampleBased_simulationBased_snp_df)
            sampleBased_indel_df_list.append(sampleBased_simulationBased_indel_df)

        print('Latest sampleIndex: %s' %(sampleIndex))

        # Create one df from list of dfs
        simulationBased_allSamples_snp_df = pd.concat(sampleBased_snp_df_list, axis=0)
        simulationBased_allSamples_indel_df = pd.concat(sampleBased_indel_df_list, axis=0)

        #Free up memory
        del sampleBased_snp_df_list
        del sampleBased_indel_df_list

        print('#############################################')
        print('simulationBased_allSamples_snp_df.shape')
        print(simulationBased_allSamples_snp_df.shape)

        print('simulationBased_allSamples_indel_df.shape')
        print(simulationBased_allSamples_indel_df.shape)
        print('#############################################')

        ##########################################################################################################
        ######################################## New way: Indels starts ##########################################
        ##########################################################################################################
        if ((simulationBased_allSamples_indel_df is not None) and (not simulationBased_allSamples_indel_df.empty)):

            print('########  before all ##########')
            print('simulationBased_allSamples_indel_df.shape')
            print(simulationBased_allSamples_indel_df.shape)

            print('simulationBased_allSamples_indel_df.head()')
            print(simulationBased_allSamples_indel_df.head())

            print('simulationBased_allSamples_indel_df first 2 rows')
            print(simulationBased_allSamples_indel_df.iloc[0, :])
            print(simulationBased_allSamples_indel_df.iloc[1, :])

            print('simulationBased_allSamples_indel_df.tail()')
            print(simulationBased_allSamples_indel_df.tail())

            # left here new columns
            # column_names = ['chrNumber', 'start', 'sample', 'originalNucleotide', 'mutatedNucleotide','dot' ,'dataSource','genomeAssembly', 'type', 'strand']
            # 1       198717710       PD10010a        ATTTG   A       .       Simulations     GRCh37  4:Del:M:2       +1

            simulationBased_allSamples_indel_df.columns = [CHROM, START, SAMPLE, REF, ALT, DOT, DATASOURCE, GENOME, MUTATION, STRAND]

            # Add Length column to the indels_df
            # simulationBased_allSamples_indel_df[LENGTH] = abs(simulationBased_allSamples_indel_df.iloc[:,3].str.len() - simulationBased_allSamples_indel_df.iloc[:,4].str.len()) +1
            simulationBased_allSamples_indel_df[LENGTH] = abs(simulationBased_allSamples_indel_df[REF].str.len() - simulationBased_allSamples_indel_df[ALT].str.len()) + 1

            # Add columns
            simulationBased_allSamples_indel_df[TYPE] = INDEL

            #End inclusive.
            # simulationBased_allSamples_indel_df[END] = simulationBased_allSamples_indel_df.iloc[:,1] + simulationBased_allSamples_indel_df[LENGTH] - 1
            simulationBased_allSamples_indel_df[END] = simulationBased_allSamples_indel_df[START] + simulationBased_allSamples_indel_df[LENGTH] - 1

            # Add a new column Category
            # simulationBased_allSamples_indel_df['Category'] = np.where(simulationBased_allSamples_indel_df['Length'] > 3, 'INDEL(>3bp)', 'INDEL(<=3bp)')
            simulationBased_allSamples_indel_df[CATEGORY] = np.where(simulationBased_allSamples_indel_df[LENGTH] > 3, INDEL_GT_3BP, INDEL_LTE_3BP)

            # drop the unnecessary columns
            # simulationBased_allSamples_indel_df.drop([5,6,7,8,9], inplace=True, axis=1)
            #TODO We may not drop MUTATION and STRAND
            simulationBased_allSamples_indel_df.drop([DOT,DATASOURCE,GENOME,MUTATION,STRAND], inplace=True, axis=1)

            #Provide column names
            simulationBased_allSamples_indel_df.columns = [CHROM, START, SAMPLE, REF, ALT, LENGTH, TYPE, END,  CATEGORY]

            #Reorder columns
            simulationBased_allSamples_indel_df = simulationBased_allSamples_indel_df[[SAMPLE, CHROM, START, END, REF, ALT, TYPE, LENGTH, CATEGORY]]

            print('############ after all ##########')
            print('simulationBased_allSamples_indel_df.shape')
            print(simulationBased_allSamples_indel_df.shape)

            print('simulationBased_allSamples_indel_df.head()')
            print(simulationBased_allSamples_indel_df.head())

            print('simulationBased_allSamples_indel_df first 2 rows')
            print(simulationBased_allSamples_indel_df.iloc[0,:])
            print(simulationBased_allSamples_indel_df.iloc[1,:])

            print('simulationBased_allSamples_indel_df.tail()')
            print(simulationBased_allSamples_indel_df.tail())

            # Save to a file
            simulationBased_allSamples_indel_df.to_csv(simulationBased_allSamples_indel_filepath, sep='\t', header=True, index=False)
        ##########################################################################################################
        ######################################## New way: Indels ends ############################################
        ##########################################################################################################


        ##########################################################################################################
        ######################################## New way: SNPs starts ############################################
        ##########################################################################################################
        if ((simulationBased_allSamples_snp_df is not None) and (not simulationBased_allSamples_snp_df.empty)):

            print('############# before all ############')
            print('simulationBased_allSamples_snp_df.shape')
            print(simulationBased_allSamples_snp_df.shape)

            print('simulationBased_allSamples_snp_df.head()')
            print(simulationBased_allSamples_snp_df.head())

            print('for debug: All pyrimidines? simulationBased_allSamples_snp_df[3].unique()')
            print(simulationBased_allSamples_snp_df[3].unique())

            print('simulationBased_allSamples_snp_df first 2 rows')
            print(simulationBased_allSamples_snp_df.iloc[0,:])
            print(simulationBased_allSamples_snp_df.iloc[1,:])

            print('simulationBased_allSamples_snp_df.tail()')
            print(simulationBased_allSamples_snp_df.tail())

            #Provide column names
            # simulationBased_allSamples_snp_df.columns = ['Chromosome', 'Start', 'Sample', 'Reference', 'Mutation','Dot','DataSource','Assembly', 'Context', 'Strand']
            # simulationBased_allSamples_snp_df = simulationBased_allSamples_snp_df.astype({'Chromosome': 'str', 'Start': 'int', 'Sample': 'str', 'Reference': 'str',
            #                                                                               'Mutation': 'str', 'Dot': 'str','DataSource': 'str', 'Assembly': 'str',
            #                                                                               'Context': 'str', 'Strand': 'int'})

            print('simulationBased_allSamples_snp_df.columns.values')
            print(simulationBased_allSamples_snp_df.columns.values)

            print('simulationBased_allSamples_snp_df.dtypes')
            print(simulationBased_allSamples_snp_df.dtypes)

            # left here new columns
            # column_names = ['chrNumber', 'start', 'sample', 'originalNucleotide', 'mutatedNucleotide','dot' ,'dataSource','genomeAssembly', 'type', 'strand']
            # X       25826390        PD10010a        T       A       .       Simulations     GRCh37  CTG     +1
            simulationBased_allSamples_snp_df.columns= [CHROM, START, SAMPLE, REF, ALT, DOT, DATASOURCE, GENOME, CONTEXT, PYRAMIDINESTRAND]


            # Add columns
            # simulationBased_allSamples_snp_df[END] = simulationBased_allSamples_snp_df.iloc[:, 1]
            simulationBased_allSamples_snp_df[END] = simulationBased_allSamples_snp_df[START]

            # Add Mutations column
            # For merging with probabilities
            # simulationBased_allSamples_snp_df[MUTATIONS] = simulationBased_allSamples_snp_df.iloc[:, 8].str[0] + '[' + \
            #                                                  simulationBased_allSamples_snp_df.iloc[:,3] + '>' + simulationBased_allSamples_snp_df.iloc[:,4] + ']' + \
            #                                                  simulationBased_allSamples_snp_df.iloc[:, 8].str[2]
            simulationBased_allSamples_snp_df[MUTATIONS] = simulationBased_allSamples_snp_df[CONTEXT].str[0] + '[' + \
                                                             simulationBased_allSamples_snp_df[REF] + '>' + simulationBased_allSamples_snp_df[ALT] + ']' + \
                                                             simulationBased_allSamples_snp_df[CONTEXT].str[2]


            # simulationBased_allSamples_snp_df[MUTATION] = simulationBased_allSamples_snp_df.iloc[:,3] + '>' + simulationBased_allSamples_snp_df.iloc[:,4]
            simulationBased_allSamples_snp_df[MUTATION] = simulationBased_allSamples_snp_df[REF] + '>' + simulationBased_allSamples_snp_df[ALT]

            #drop the unnecessary columns
            # simulationBased_allSamples_snp_df.drop([3,4,5,6,7], inplace=True, axis=1)
            simulationBased_allSamples_snp_df.drop([REF, ALT, DOT, DATASOURCE, GENOME], inplace=True, axis=1)

            #Give column names
            simulationBased_allSamples_snp_df.columns= [CHROM, START, SAMPLE, CONTEXT, PYRAMIDINESTRAND, END, MUTATIONS, MUTATION]

            #Reorder columns
            simulationBased_allSamples_snp_df = simulationBased_allSamples_snp_df[[SAMPLE, MUTATIONS, CHROM, START, END, PYRAMIDINESTRAND, MUTATION, CONTEXT]]

            print('################# after all ###############')
            print('simulationBased_allSamples_snp_df.shape')
            print(simulationBased_allSamples_snp_df.shape)

            print('simulationBased_allSamples_snp_df.head()')
            print(simulationBased_allSamples_snp_df.head())

            print('simulationBased_allSamples_snp_df first 2 rows')
            print(simulationBased_allSamples_snp_df.iloc[0,:])
            print(simulationBased_allSamples_snp_df.iloc[1,:])

            print('simulationBased_allSamples_snp_df.tail()')
            print(simulationBased_allSamples_snp_df.tail())

            print('probabilities_df.columns.values')
            print(probabilities_df.columns.values)

            print('simulationBased_allSamples_snp_df.columns.values')
            print(simulationBased_allSamples_snp_df.columns.values)

            print('for debug starts')
            print('Unique samples in probabilities_df')
            print(probabilities_df[SAMPLE].unique())
            print('Unique Mutations in probabilities_df')
            print(probabilities_df[MUTATIONS].unique())

            print('Unique samples in simulationBased_allSamples_snp_df')
            print(simulationBased_allSamples_snp_df[SAMPLE].unique())
            print('Unique Mutations in simulationBased_allSamples_snp_df')
            print(simulationBased_allSamples_snp_df[MUTATIONS].unique())
            print('for debug ends')

            # Merge with signature probabilities
            # In probabilities_df there is no 'Cancer Type' or 'Project' column, there is 'Sample' and 'Mutations' columns
            result_df = pd.merge(simulationBased_allSamples_snp_df, probabilities_df, how='inner',
                                 left_on=[SAMPLE, MUTATIONS],
                                 right_on=[SAMPLE, MUTATIONS])

            #Drop Mutations column after merge
            result_df.drop([MUTATIONS], inplace=True, axis=1)

            print('################# after merge ###############')
            print('result_df.shape')
            print(result_df.shape)

            print('result_df.head()')
            print(result_df.head())

            print('result_df first 2 rows')
            print(result_df.iloc[0, :])
            print(result_df.iloc[1, :])

            print('result_df.tail()')
            print(result_df.tail())

            # save to a file
            result_df.to_csv(simulationBased_allSamples_snp_filepath, sep='\t', header=True, index=False)
        ##########################################################################################################
        ######################################## New way: SNPs ends ##############################################
        ##########################################################################################################

################################################################
