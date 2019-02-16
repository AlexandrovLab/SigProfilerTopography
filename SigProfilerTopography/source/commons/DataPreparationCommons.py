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
    # Let's read the ';' at the dataSource column sometimes
    # mutation_df = pd.read_table(mutationsWithGenomicPositionFilePath, sep="\t", header=None, comment=';')

    print('###########################################')
    print('Probabilities information starts')

    print('probabilities_df.shape')
    print(probabilities_df.shape)
    print('probabilities_df.head()')
    print(probabilities_df.head())

    print('Unique samples in probabilities_df')
    print(probabilities_df['Sample'].unique())
    print('# of unique samples in probabilities_df: %d' %(len(probabilities_df['Sample'].unique())))
    print()

    print('Unique mutations in probabilities_df')
    print(probabilities_df['Mutations'].unique())
    print('# of unique mutations in probabilities_df: %d' %(len(probabilities_df['Mutations'].unique())))
    print()

    print('Probabilities information ends')
    print('###########################################')

    return probabilities_df
############################################################


############################################################
#Tested works correctly.
def get3Nucleotides(chromosomeShort,start,end,humanGenome):
    if chromosomeShort == 'MT':
        chromosomeShort = 'M'
    elif chromosomeShort == '23':
        chromosomeShort = 'X'
    elif chromosomeShort == '24':
        chromosomeShort = 'Y'

    chromosomeLong = 'chr' + chromosomeShort
    chrBased_humanGenome = humanGenome[chromosomeLong]
    seq = chrBased_humanGenome.get_slice(start-2, end + 1)
    seq = seq.upper()
    return seq
############################################################

############################################################
def complement(upperNucleotides):
    complemented = ''
    for upperNucleotide in upperNucleotides:
        complemented += complementDict[upperNucleotide]
    return complemented
############################################################


# ############################################################
# def getMutationInformationOld(row,hg19_genome):
#     originalNucleotide = row['originalNucleotide']
#     mutatedNucleotide = row['mutatedNucleotide']
#
#     pyramidineStrand = 1
#     # print('chrNumber: %s -- start: %d -- end: %d' %(chrNumber,start,end))
#     threeNucleotidesUpper = get3Nucleotides(row['chrNumber'], row['start'], row['end'], hg19_genome)
#
#     if (originalNucleotide=='G' or originalNucleotide=='A' or originalNucleotide=='g' or originalNucleotide=='a'):
#         originalNucleotide = complement(originalNucleotide.upper())
#         mutatedNucleotide = complement(mutatedNucleotide.upper())
#         threeNucleotidesUpper = complement(threeNucleotidesUpper)
#         pyramidineStrand = -1
#
#     mutationType = '%s>%s' %(originalNucleotide,mutatedNucleotide)
#     mutationSubtype = threeNucleotidesUpper
#
#     row['PyramidineStrand'] = pyramidineStrand
#     row['Mutation Type'] = mutationType
#     row['Mutation Subtype'] = mutationSubtype
#
#     # print('------- %s %s  %d' %(row['Mutation Type'],row['Mutation Subtype'],row['PyramidineStrand']))
#     return row
# ############################################################


#new code starts

############################################################
def getMutationInformation(row,hg19_genome, hg38_genome):
    originalNucleotide = row['Ref']
    mutatedNucleotide = row['Alt']
    pyramidineStrand = convertStrandIntoNum(row['Strand'])

    if row['Genome'] == 'GRCh37':
        threeNucleotidesUpper = get3Nucleotides(row['Chrom'], row['Start'], row['Start'], hg19_genome)
    elif row['Genome'] == 'GRCh38':
        threeNucleotidesUpper = get3Nucleotides(row['Chrom'], row['Start'], row['Start'], hg38_genome)

    if (originalNucleotide=='G' or originalNucleotide=='A' or originalNucleotide=='g' or originalNucleotide=='a'):
        originalNucleotide = complement(originalNucleotide.upper())
        mutatedNucleotide = complement(mutatedNucleotide.upper())
        threeNucleotidesUpper = complement(threeNucleotidesUpper)
        pyramidineStrand = -1*pyramidineStrand

    mutations = '%s[%s>%s]%s' %(threeNucleotidesUpper[0],originalNucleotide,mutatedNucleotide,threeNucleotidesUpper[2])

    row['Mutations'] = mutations
    row['PyramidineStrand'] = pyramidineStrand
    row['Mutation'] = '%s>%s' %(originalNucleotide,mutatedNucleotide)
    row['Context'] = threeNucleotidesUpper

    return row
############################################################


############################################################
def readIndelsandWriteWithGenomicPositions(jobname,indelsInputFile):
    # Project	Sample	locID	Genome	mutType	Chrom	Start	End	Ref	Alt	Type	VarID	Gene	GeneID	ccdsID	TranscriptID
    # 21_BRCA_WGS	PD3851a	.	GRCh37	INDEL	6	45960535	45960538	GTA	G	SOMATIC	67936580	CLIC5	ENSG00000112782	CCDS47438.1	ENST00000185206

    column_names = ['Project', 'Sample', 'locID', 'Genome', 'mutType', 'Chrom', 'Start', 'End', 'Ref', 'Alt', 'Type', 'VarID', 'Gene', 'GeneID', 'ccdsID', 'TranscriptID']
    indels_with_genomic_positions_df = pd.read_table(indelsInputFile, sep="\t", header=0, names=column_names,dtype={'Chrom': str, 'Sample':str})

    print('###########################################')
    print('Indels information starts')
    print('indels_with_genomic_positions_df.head()')
    print(indels_with_genomic_positions_df.head())

    print('indels_with_genomic_positions_df.columns.values')
    print(indels_with_genomic_positions_df.columns.values)

    print('indels_with_genomic_positions_df.shape')
    print(indels_with_genomic_positions_df.shape)

    print('# of rows in indels_with_genomic_positions_df: %d' % len(indels_with_genomic_positions_df))

    # How many projects (cancer type)  are there?
    print('Unique projects (cancer type) in the indelsInputFile:')
    print(indels_with_genomic_positions_df['Project'].unique())
    print('# of projects in the indelsInputFile: %d' % len(indels_with_genomic_positions_df['Project'].unique()))
    print()

    #How many samples are there?
    print('Unique sample names in the indelsInputFile:')
    print(indels_with_genomic_positions_df['Sample'].unique())
    print('# of samples in the indelsInputFile: %d' %len(indels_with_genomic_positions_df['Sample'].unique()))
    print()

    print('Unique chroms in the indelsInputFile:')
    print(indels_with_genomic_positions_df['Chrom'].unique())
    print('# of chroms in the indelsInputFile: %d' % len(indels_with_genomic_positions_df['Chrom'].unique()))
    print()

    print('Unique genomes in the indelsInputFile:')
    print(indels_with_genomic_positions_df['Genome'].unique())
    print('# of genomes in the indelsInputFile: %d' % len(indels_with_genomic_positions_df['Genome'].unique()))
    print()

    print('Unique types in the indelsInputFile:')
    print(indels_with_genomic_positions_df['Type'].unique())
    print('# of types in the indelsInputFile: %d' % len(indels_with_genomic_positions_df['Type'].unique()))
    print()

    print('Indels information ends')
    print('###########################################')

    # Drop the columns
    indels_with_genomic_positions_df.drop(['locID', 'Genome', 'Type', 'VarID', 'Gene', 'GeneID', 'ccdsID', 'TranscriptID'],inplace=True, axis=1)

    #Add the requested columns
    # Vectorized implementation on Pandas series
    indels_with_genomic_positions_df['Length'] = abs(indels_with_genomic_positions_df['Ref'].str.len() - indels_with_genomic_positions_df['Alt'].str.len())
    indels_with_genomic_positions_df['Category'] = np.where(indels_with_genomic_positions_df['Length'] > 3,'INDEL(>3bp)', 'INDEL(<=3bp)')

    # Rename
    indels_with_genomic_positions_df.rename(columns={'Start' : 'Position Start', 'End' : 'Position End','Ref':'Reference','Alt':'Mutation', 'mutType' : 'Type'}, inplace=True)


    #Order the columns
    # Column order we want to have:
    ordered_column_names = ['Project','Sample',	'Chrom',	'Position Start',	'Position End',	'Reference',	'Mutation',	'Type',	'Length','Category']
    indels_with_genomic_positions_df = indels_with_genomic_positions_df[ordered_column_names]


    #########################################################
    #Write indels_df by project
    indels_df_groupby_project = indels_with_genomic_positions_df.groupby('Project')

    numberofGroupsBasedonProject = len(indels_df_groupby_project)
    print('numberofGroupsBasedonProject:%d' %(numberofGroupsBasedonProject))


    for projectType, projectBased_indels_df  in indels_df_groupby_project:
        if (numberofGroupsBasedonProject==1):
            mutations_indels_filePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP,INPUT,jobname,jobname + '_indels_for_topography.txt')
        else:
            mutations_indels_filePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP,INPUT,jobname,jobname + '_' + projectType + '_indels_for_topography.txt')
        if (projectBased_indels_df is not None):
            projectBased_indels_df.to_csv(mutations_indels_filePath, sep='\t', header=True, index=False , columns = ['Sample',	'Chrom',	'Position Start',	'Position End',	'Reference',	'Mutation',	'Type',	'Length','Category'])
    #########################################################

############################################################

############################################################
def readMutationsWithGenomicPositions(snpsInputFile):
    mutationsWithGenomicPositionFilePath = os.path.join(snpsInputFile)

    column_names = ['Project', 'Sample', 'locID', 'Genome',	'mutType', 'Chrom',	'Start', 'End',	'Ref', 'Alt', 'Type', 'VarID', 'Strand', 'Gene', 'GeneID','ccdsID', 'TranscriptID',	'GeneType']
    snps_with_genomic_positions_df = pd.read_table(mutationsWithGenomicPositionFilePath, sep="\t", header=0, names=column_names,dtype={'Chrom': str, 'Sample':str})

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
    print(snps_with_genomic_positions_df['Project'].unique())
    print('# of projects in the snpsInputFile: %d' % len(snps_with_genomic_positions_df['Project'].unique()))
    print()

    #How many samples are there?
    print('Unique sample names in the snpsInputFile:')
    print(snps_with_genomic_positions_df['Sample'].unique())
    print('# of samples in the snpsInputFile: %d' %len(snps_with_genomic_positions_df['Sample'].unique()))
    print()

    print('Unique chroms in the snpsInputFile:')
    print(snps_with_genomic_positions_df['Chrom'].unique())
    print('# of chroms in the snpsInputFile: %d' % len(snps_with_genomic_positions_df['Chrom'].unique()))
    print()

    print('Unique genomes in the snpsInputFile:')
    print(snps_with_genomic_positions_df['Genome'].unique())
    print('# of genomes in the snpsInputFile: %d' % len(snps_with_genomic_positions_df['Genome'].unique()))
    print()

    print('Unique types in the snpsInputFile:')
    print(snps_with_genomic_positions_df['Type'].unique())
    print('# of types in the snpsInputFile: %d' % len(snps_with_genomic_positions_df['Type'].unique()))
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
    chrBased_cancerTypeBased_df = chrBased_snps_df.apply(getMutationInformation,hg19_genome=hg19_genome,hg38_genome=hg38_genome,axis=1)

    print('For chr%s ' %chrNumber)
    result_df = pd.merge(chrBased_cancerTypeBased_df, probabilitiles_df, how='inner', left_on=['Sample','Mutations'], right_on=['Sample','Mutations'])

    #Drop the unnecessary columns after getMutationInformation
    result_df.drop(['Project', 'Genome', 'Ref', 'Alt', 'Strand'], inplace=True, axis=1)

    # rename column names
    result_df.rename(columns={'Chrom': 'Chromosome'},inplace=True)

    columnNames = list(result_df.columns.values)
    indexofContext = columnNames.index('Context')
    signatureColumnNames = columnNames[indexofContext+1:]

    # Column order we want to have:
    ordered_column_names = ['Sample', 'Chromosome', 'Start', 'End', 'PyramidineStrand', 'Mutation', 'Context']
    ordered_column_names.extend(signatureColumnNames)

    result_df = result_df[ordered_column_names]
    #####################################################

    return result_df
############################################################


############################################################
#Assume that there can be different cancer types is snps_df
#Right now there is no cancer type column in probabilities_df
def mergeSNPsWithSignatureProbabilities(
        jobname,
        snps_df,
        probabilities_df,
        hg19_genome,
        hg38_genome):

    current_abs_path = os.path.abspath(os.path.dirname(__file__))
    os.makedirs(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, INPUT,jobname), exist_ok=True)

    ################################################################
    # Get signature column names after Mutations column
    signature_probabilities_columnNames = probabilities_df.columns.values
    mutationsIndex = list(signature_probabilities_columnNames).index('Mutations')
    signatureColumnNames = signature_probabilities_columnNames[mutationsIndex + 1:]
    print(signatureColumnNames)

    print('len(signatureColumnNames)')
    print(len(signatureColumnNames))
    ################################################################

    ################################################################
    print('#################################')
    snps_df_groupedby_project = snps_df.groupby(['Project'])
    print('#################################')
    ################################################################

    ############################################################
    ############### For each Cancer Type #######################
    ############################################################
    #sort each grouped_sub_df w.r.t. its size in ascending order
    snps_df_groupedbyproject_sorted = sorted(snps_df_groupedby_project,  # iterates pairs of (key, corresponding subDataFrame)
                            key=lambda x: len(x[1])  # sort by number of rows (len of subDataFrame)
                            )  # sort i.e. smallest first


    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)

    #We assume that we are doing analyses for a certain cancer type
    numberofGroupsBasedonProject = len(snps_df_groupedby_project)
    print('numberofGroupsBasedonProject:%d' %(numberofGroupsBasedonProject))

    #Do sequential for each cancerType
    #Start with the cancerType with sm allest number of samples
    for cancerType, cancerTypeBased_snps_df in snps_df_groupedbyproject_sorted:
        # if (cancerType == 'Skin-Melanoma'):
        print('################# For cancerType: %s starts ##########################' % cancerType)

        if (numberofGroupsBasedonProject==1):
            mutations_snps_snvs_filePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, INPUT,jobname,jobname + '_snps_for_topography.txt')
        else:
            mutations_snps_snvs_filePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP,ONE_DIRECTORY_UP, INPUT,jobname,jobname + '_' + cancerType + '_snps_for_topography.txt')


        # ################################################
        # #Sequential takes 10 minutes
        # chrBased_cancerTypeBased_df = cancerTypeBased_snps_df.apply(getMutationInformation, hg19_genome=hg19_genome,hg38_genome=hg38_genome, axis=1)
        # chrBased_cancerTypeBased_df.to_csv(mutations_snps_snvs_filePath, sep='\t', header=True, index=False)
        # ################################################

        ################################################
        # Parallel takes 3 minutes
        #group by chromosome
        chrBased_cancerTypeBased_grouped = cancerTypeBased_snps_df.groupby('Chrom')

        #Print number of groups
        print('len(chrBased_cancerTypeBased_grouped)')
        print(len(chrBased_cancerTypeBased_grouped))

        print('##############################################')
        print('cancerTypeBased_df.shape')
        print(cancerTypeBased_snps_df.shape)
        print('##############################################')

        print('##############################################')
        print('probabilitles_df.shape')
        print(probabilities_df.shape)
        print('##############################################')

        # Prepare the input starts
        chrBased_cancerTypeBased_Genomic_Positions_PoolInputList = []

        #Do in parallel for all chrosomomes of an cancer type
        for chrNumber, chrBased_cancerTypeBased_df in chrBased_cancerTypeBased_grouped:
            inputPoolList=[]
            inputPoolList.append(cancerType)
            inputPoolList.append(chrNumber)
            inputPoolList.append(chrBased_cancerTypeBased_df)

            #There is no chromosome information in signature probabilities dataframe
            inputPoolList.append(probabilities_df)
            inputPoolList.append(hg19_genome)
            inputPoolList.append(hg38_genome)
            chrBased_cancerTypeBased_Genomic_Positions_PoolInputList.append(inputPoolList)
        # Prepare the input ends

        #print number of elements in chrBased_cancerTypeBased_Genomic_Positions_PoolInputList
        print('len(chrBased_PoolInputList)')
        print(len(chrBased_cancerTypeBased_Genomic_Positions_PoolInputList))

        all_chrBased_combined_df_list = pool.map(combineGenomicPositionsWithSignatureProbabilitiesParallel,chrBased_cancerTypeBased_Genomic_Positions_PoolInputList)
        all_chroms_concatenated_df = pd.concat(all_chrBased_combined_df_list, axis=0)

        print('########################################')
        print('all_chroms_concatenated_df.shape')
        print(all_chroms_concatenated_df.shape)
        print('all_chroms_concatenated_df.columns.values')
        print(all_chroms_concatenated_df.columns.values)
        print('########################################')

        if (all_chroms_concatenated_df is not None):
            all_chroms_concatenated_df.to_csv(mutations_snps_snvs_filePath, sep='\t', header=True, index=False)
        ################################################

        print('################# For cancerType: %s ends ##########################' %cancerType)
        ########################################################

    ############################################################
    ############### For each Cancer Type #######################
    ############################################################

    pool.close()
    pool.join()
############################################################


################################################################
#Read a sample file for a certain simulationNumber
#Returns snp_df and indel_df for each read sample file
def readSampleBasedSimulationBasedDF(sampleBasedSimulationsInputPath,sampleName,simulationNumber,mutationTypes):
    # sampleFileName =PD10010a_INDEL_96_10.txt
    sampleFileName = '%s_%s_%s_%d.txt' %(sampleName,mutationTypes[1],mutationTypes[0],simulationNumber)
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
def prepareSimulationBasedInputFilesForSigProfilerTopography(jobname,genomeAssembly,mutationTypes,sigProfilerSimulatorSpecificDirName,numberofSimulations,probabilities_df):

    current_abs_path = os.path.dirname(os.path.realpath(__file__))

    #sigProfilerSimulatorSpecificDirName
    #21BreastCancer_simulations_GRCh37_INDEL_96

    # /oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/input/21BreastCancer/output/simulations/21BreastCancer_simulations_GRCh37_INDEL_96/
    sampleBasedSimulationsInputPath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, INPUT, jobname, OUTPUT, SIMULATIONS, sigProfilerSimulatorSpecificDirName)

    #/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/input/21BreastCancer/SIMULATIONS_FOR_TOPOGRAPHY
    output_path = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,INPUT,jobname,SIMULATIONS_FOR_TOPOGRAPHY)

    #read hg19_genome
    if (genomeAssembly==HG19):
        human_genome = twobitreader.TwoBitFile(os.path.join(current_abs_path, ONE_DIRECTORY_UP, LIB, UCSCGENOME, HG19_2BIT))
    elif (genomeAssembly==HG38):
        human_genome = twobitreader.TwoBitFile(os.path.join(current_abs_path, ONE_DIRECTORY_UP, LIB, UCSCGENOME, HG38_2BIT))

    if os.path.exists(output_path):
        os.system("rm -rf " + output_path)

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
        if (simulationBased_allSamples_indel_df is not None):

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
            # X       77873197        PD10010a        T       A       .       Simulations     GRCh37  TTT     +1



            # Add Length column to the indels_df
            # simulationBased_allSamples_indel_df['Length'] = simulationBased_allSamples_indel_df.apply(get_indel_length,axis=1)
            simulationBased_allSamples_indel_df['Length'] = abs(simulationBased_allSamples_indel_df.iloc[:,3].str.len() - simulationBased_allSamples_indel_df.iloc[:,4].str.len()) +1

            # Add columns
            simulationBased_allSamples_indel_df['Type'] = 'INDEL'

            # simulationBased_allSamples_indel_df['End'] = simulationBased_allSamples_indel_df['Start'] + simulationBased_allSamples_indel_df['Length']
            simulationBased_allSamples_indel_df['Position End'] = simulationBased_allSamples_indel_df.iloc[:,1] + simulationBased_allSamples_indel_df['Length'] - 1

            # Add a new column Category
            # simulationBased_allSamples_indel_df['Category'] = np.where(simulationBased_allSamples_indel_df['Length'] > 3, 'INDEL(>3bp)', 'INDEL(<=3bp)')
            simulationBased_allSamples_indel_df['Category'] = np.where(simulationBased_allSamples_indel_df['Length'] > 3, 'INDEL(>3bp)', 'INDEL(<=3bp)')

            # drop the unnecessary columns
            simulationBased_allSamples_indel_df.drop([5,6,7,8,9], inplace=True, axis=1)

            #Provide column names
            simulationBased_allSamples_indel_df.columns = ['Chrom', 'Position Start', 'Sample', 'Reference', 'Mutation', 'Length', 'Type', 'Position End',  'Category']

            #Reorder columns
            simulationBased_allSamples_indel_df = simulationBased_allSamples_indel_df[['Sample', 'Chrom', 'Position Start', 'Position End', 'Reference', 'Mutation', 'Type', 'Length', 'Category']]

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
        if (simulationBased_allSamples_snp_df is not None):

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
            #X       18050372        PD3851a C       A       .       Simulations     GRCh37  ACT     -1
            #X       77873197        PD10010a        T       A       .       Simulations     GRCh37  TTT     +1


            # Add columns
            simulationBased_allSamples_snp_df['End'] = simulationBased_allSamples_snp_df.iloc[:, 1]

            # Add Mutations column
            # simulationBased_allSamples_snp_df['Mutations'] = simulationBased_allSamples_snp_df.iloc[:,8] + '[' + simulationBased_allSamples_snp_df['Reference2Mutation'] + ']'+ simulationBased_allSamples_snp_df.iloc[:,8]
            simulationBased_allSamples_snp_df['Mutations'] = simulationBased_allSamples_snp_df.iloc[:, 8].str[0] + '[' + \
                                                             simulationBased_allSamples_snp_df.iloc[:,3] + '>' + simulationBased_allSamples_snp_df.iloc[:,4] + ']' + \
                                                             simulationBased_allSamples_snp_df.iloc[:, 8].str[2]

            simulationBased_allSamples_snp_df['Mutation'] = simulationBased_allSamples_snp_df.iloc[:,3] + '>' + simulationBased_allSamples_snp_df.iloc[:,4]


            #drop the unnecessary columns
            simulationBased_allSamples_snp_df.drop([3,4,5,6,7], inplace=True, axis=1)

            #Give column names
            simulationBased_allSamples_snp_df.columns= ['Chromosome', 'Start', 'Sample','Context', 'PyramidineStrand', 'End', 'Mutations', 'Mutation' ]

            #Reorder columns
            simulationBased_allSamples_snp_df = simulationBased_allSamples_snp_df[['Sample', 'Mutations', 'Chromosome', 'Start', 'End', 'PyramidineStrand', 'Mutation', 'Context']]

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
            print(probabilities_df['Sample'].unique())
            print('Unique Mutations in probabilities_df')
            print(probabilities_df['Mutations'].unique())

            print('Unique samples in simulationBased_allSamples_snp_df')
            print(simulationBased_allSamples_snp_df['Sample'].unique())
            print('Unique Mutations in simulationBased_allSamples_snp_df')
            print(simulationBased_allSamples_snp_df['Mutations'].unique())
            print('for debug ends')

            # Merge with signature probabilities
            # In probabilities_df there is no 'Cancer Type' or 'Project' column, there is 'Sample' and 'Mutations' columns
            result_df = pd.merge(simulationBased_allSamples_snp_df, probabilities_df, how='inner',
                                 left_on=['Sample', 'Mutations'],
                                 right_on=['Sample', 'Mutations'])

            #Drop Mutations column after merge
            result_df.drop(['Mutations'], inplace=True, axis=1)

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
