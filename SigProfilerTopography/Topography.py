# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018 Burcak Otlu

# #############################################################
# import sys
# import os
# current_abs_path = os.path.dirname(os.path.realpath(__file__))
# commonsPath = os.path.join(current_abs_path,'commons')
# sys.path.append(commonsPath)
# #############################################################

from SigProfilerTopography.source.commons.DataPreparationCommons import readMutationsWithGenomicPositions
from SigProfilerTopography.source.commons.DataPreparationCommons import readProbabilities
from SigProfilerTopography.source.commons.DataPreparationCommons import mergeSNPsWithSignatureProbabilities
from SigProfilerTopography.source.commons.DataPreparationCommons import readIndelsandWriteWithGenomicPositions
from SigProfilerTopography.source.commons.DataPreparationCommons import prepareSimulationBasedInputFilesForSigProfilerTopography
from SigProfilerTopography.source.commons.DataPreparationCommons import readChrBasedDinucsMergeWithProbabilitiesAndWrite

from SigProfilerTopography.source.commons.TopographyCommons import *

from SigProfilerTopography.source.commons.PartitionIndelsData import partitionIndelsData
from SigProfilerTopography.source.commons.PartitionSinglePointMutationsData import partitionSubsData

from SigProfilerTopography.source.commons.NucleosomeOccupancySignalCountArraysAndFigures import readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArrays
from SigProfilerTopography.source.nucleosomeoccupancy.NucleosomeOccupancyAnalysis import nucleosomeOccupancyAnalysis
from SigProfilerTopography.source.replicationtime.ReplicationTimeAnalysis import replicationTimeAnalysis
from SigProfilerTopography.source.replicationtime.ReplicationTimeAnalysis import readReplicationTimeDataAndWriteChrBasedReplicationTimeNPArrays


from SigProfilerTopography.source.replicationstrandbias.ReplicationStrandBiasAnalysis import replicationStrandBiasAnalysis
from SigProfilerTopography.source.transcriptionstrandbias.TranscriptionStrandBiasAnalysis import transcriptionStrandBiasAnalysis
from SigProfilerTopography.source.processivity.ProcessivityAnalysis import processivityAnalysis

from SigProfilerTopography.source.plotting.NucleosomeOccupancyAverageSignalFigures import nucleosomeOccupancyAverageSignalFigures
from SigProfilerTopography.source.plotting.ReplicationTimeNormalizedMutationDensityFigures import replicationTimeNormalizedMutationDensityFigures
from SigProfilerTopography.source.plotting.TranscriptionReplicationStrandBiasFigures import transcriptionReplicationStrandBiasFigures
from SigProfilerTopography.source.plotting.ProcessivityFigures import processivityFigures



############################################################
#CAn be move to DataPreparationCommons under /source/commons
#read chr based dinucs (provided by SigProfilerMatrixGenerator) and merge with probabilities (provided by SigProfilerTopography)
def prepareDinucsDataAfterMatrixGenerationAndExtractorForTopography(outputDir, jobname, matrix_generator_output_dir_path,dinucs_probabilities_file_path):
    os.makedirs(os.path.join(outputDir,jobname,DATA,CHRBASED),exist_ok=True)
    chrShortList = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
    if (os.path.exists(dinucs_probabilities_file_path)):
        dinucs_probabilities_df = readProbabilities(dinucs_probabilities_file_path)
        # dinucs_probabilities_df.columns.names [Sample Names    MutationTypes   DBS2    DBS4    DBS6    DBS7    DBS11]
        # dinucs_probabilities_df.columns.names [Sample    Mutation   DBS2    DBS4    DBS6    DBS7    DBS11]
        dinucs_probabilities_df.rename(columns={'Sample Names': 'Sample', 'MutationTypes': 'Mutation'}, inplace=True)

    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)

    poolInputList = []
    if (os.path.exists(matrix_generator_output_dir_path)):
        for chrShort in chrShortList:
            chr_based_dinuc_filename = '%s_seqinfo.txt' %(chrShort)
            chr_based_dinuc_filepath = os.path.join(matrix_generator_output_dir_path,chr_based_dinuc_filename)
            inputList=[]
            inputList.append(chrShort)
            inputList.append(outputDir)
            inputList.append(jobname)
            inputList.append(chr_based_dinuc_filepath)
            inputList.append(dinucs_probabilities_df)
            poolInputList.append(inputList)

        pool.map(readChrBasedDinucsMergeWithProbabilitiesAndWrite,poolInputList)

    pool.close()
    pool.join()

    # dinucs_probabilities_file_path
    # /oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/output/560_BRCA_WGS_DINUCS/DBS78/Suggested_Solution/Decomposed_Solution/Mutation_Probabilities.txt
    # Sample Names    MutationTypes   DBS2    DBS4    DBS6    DBS7    DBS11
    # PD10011a        AC>CA   0.0     0.04003180196828397     0.04286819373415273     0.887661695818471       0.029438308479092352
    # PD10011a        AC>CG   0.0     0.003750743595000353    0.01840082543246141     0.9452478818045245      0.03260054916801382

    #matrix_generator_output_directory
    #/ oasis / tscc / scratch / burcak / developer / python / SigProfilerTopography / SigProfilerTopography / input_for_matgen / BreastCancer560_subs_indels_dinucs / output / vcf_files / DINUC
    #chr based files
    # 10_seqinfo.txt  13_seqinfo.txt  16_seqinfo.txt  19_seqinfo.txt  21_seqinfo.txt  3_seqinfo.txt  6_seqinfo.txt  9_seqinfo.txt
    # 11_seqinfo.txt  14_seqinfo.txt  17_seqinfo.txt  1_seqinfo.txt   22_seqinfo.txt  4_seqinfo.txt  7_seqinfo.txt  X_seqinfo.txt
    # 12_seqinfo.txt  15_seqinfo.txt  18_seqinfo.txt  20_seqinfo.txt  2_seqinfo.txt   5_seqinfo.txt  8_seqinfo.txt  Y_seqinfo.txt
############################################################


############################################################
#Step1 read snp positions
#Step2 read probabilities
#Step3 merge snp positions with probabilities make ready for SigProfilerTopography
#Step4 read indels make ready for SigProfilerTopography
def prepareDataAfterExtractorForTopography(snpsInputFile,indelsInputFile,probabilitiesFile,jobname):
    current_abs_path = os.path.dirname(os.path.realpath(__file__))

    #############################################################################
    #Step1 read snps with genomic positions circa 28 million rows
    snps_df = readMutationsWithGenomicPositions(snpsInputFile)

    #Drop the unnecessary columns before the merge with probabilities
    snps_df.drop(['locID','mutType','Type','VarID','Gene','GeneID','ccdsID','TranscriptID','GeneType'], inplace=True, errors='ignore',axis=1)

    print('snps_df.columns.values')
    print(snps_df.columns.values)

    #If snps start and end are  not equal
    #We can make them equal here.
    if (not snps_df[START].equals(snps_df[END])):
        snps_df[END] = snps_df[START]
    #############################################################################

    #############################################################################
    # Step2 read sample based mutation based signature probabilities
    probabilities_df = readProbabilities(probabilitiesFile)
    #############################################################################

    ###################################################
    hg19_genome = twobitreader.TwoBitFile(os.path.join(current_abs_path, LIB, UCSCGENOME, 'hg19.2bit'))
    hg38_genome = twobitreader.TwoBitFile(os.path.join(current_abs_path, LIB, UCSCGENOME, 'hg38.2bit'))
    ###################################################

    #############################################################################
    # Step3
    mergeSNPsWithSignatureProbabilities(jobname,snps_df,probabilities_df,hg19_genome,hg38_genome)
    #############################################################################

    #############################################################################
    # Step4
    readIndelsandWriteWithGenomicPositions(jobname,indelsInputFile)
    #############################################################################

############################################################


################################################################
#Step1 read probabilities
#Step2 combine snps coming from all samples
#Step2 combine indels coming from all samples
#Step3 merge snp positions with probabilities make snps ready for SigProfilerTopography
#Step4 make indels ready for SigProfilerTopography
def prepareDataAfterSimulatorForTopography(jobname,genomeAssembly,mutationTypes,numberofSimulations,probabilitiesFile):

    # Where are simulations?
    #/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/input/21BreastCancer/output/simulations/21BreastCancer_simulations_GRCh37_INDEL_96/

    # jobname str default=NOTSET
    # numberofSimulations int default=0

    #Read probabilities
    probabilities_df = readProbabilities(probabilitiesFile)

    if (len(mutationTypes)==2):
        sigProfilerSimulatorSpecificDirName = '%s_simulations_%s_%s_%s' %(jobname,genomeAssembly,mutationTypes[1], mutationTypes[0])
    elif (len(mutationTypes)==1):
        sigProfilerSimulatorSpecificDirName = '%s_simulations_%s_%s' %(jobname,genomeAssembly,mutationTypes[0])

    prepareSimulationBasedInputFilesForSigProfilerTopography(jobname,genomeAssembly,mutationTypes,sigProfilerSimulatorSpecificDirName,numberofSimulations,probabilities_df)
################################################################


#######################################################
def download(genome):
    if (genome == GRCh37):
        os.makedirs(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, UCSCGENOME), exist_ok=True)
        filepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, UCSCGENOME, HG19_2BIT)
        downloadFromWeb(HG19_URL, filepath)
    elif (genome == GRCh38):
        os.makedirs(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, UCSCGENOME), exist_ok=True)
        filepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, UCSCGENOME, HG38_2BIT)
        downloadFromWeb(HG38_URL, filepath)

    # if ((genome == GRCh37) and (HG19_2BIT not in availableLibraryFilenamesList)):
    #     os.makedirs(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, UCSCGENOME), exist_ok=True)
    #     filepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, UCSCGENOME, HG19_2BIT)
    #     downloadFromWeb(HG19_URL, filepath)
    #     append2File(HG19_2BIT, AVAILABLE_LIBRARY_FILENAMES_PATH)
    # elif ((genome == GRCh38) and (HG38_2BIT not in availableLibraryFilenamesList)):
    #     os.makedirs(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, UCSCGENOME), exist_ok=True)
    #     filepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, UCSCGENOME, HG38_2BIT)
    #     downloadFromWeb(HG38_URL, filepath)
    #     append2File(HG38_2BIT, AVAILABLE_LIBRARY_FILENAMES_PATH)
#######################################################


#######################################################
def runNucleosomeOccupancyAnalyses(mutationTypes,genome,singlePointMutationsFilename,indelsFilename,outputDir,jobname,nucleosomeFilename,chromSizesDict,chromNamesList,availableLibraryFilenamesList):
    #############################################
    # NUCLEOSOMEOCCUPANCYANALYSIS
    # Delete the output/jobname/DATA/NUCLEOSOMEOCCUPANCY if exists
    jobnamePath = os.path.join(outputDir,jobname,DATA,NUCLEOSOMEOCCUPANCY)

    ################################################
    if (os.path.exists(jobnamePath)):
        try:
            shutil.rmtree(jobnamePath)
        except OSError as e:
            print('Error: %s - %s.' % (e.filename, e.strerror))
    ################################################

    nucleosomeFilename_woDir = os.path.basename(nucleosomeFilename)

    if (nucleosomeFilename_woDir not in availableLibraryFilenamesList):
        quantileValue = round(float(0.97), 2)
        # PartitionNucleosomeOccupancyData.partitionNucleosomeOccupancyData(jobname,nucleosomeFilename,quantileValue)
        readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArrays(genome,quantileValue,nucleosomeFilename)
        #append
        append2File(nucleosomeFilename_woDir,AVAILABLE_LIBRARY_FILENAMES_PATH)

    # computationType = COMPUTATION_ALL_CHROMOSOMES_PARALLEL
    computationType = COMPUTATION_CHROMOSOMES_SEQUENTIAL
    nucleosomeOccupancyAnalysis(mutationTypes,computationType,chromSizesDict,chromNamesList,outputDir,jobname,singlePointMutationsFilename,indelsFilename,nucleosomeFilename_woDir)
    ###############################################
#######################################################


#######################################################
def runReplicationTimeAnalysis(mutationTypes,genome,singlePointMutationsFilename,indelsFilename,outputDir,jobname,replicationTimeFilename,chromSizesDict,chromNamesList):
    #############################################
    # REPLICATIONTIME
    # Delete the output/jobname/DATA/REPLICATIONTIME if exists
    jobnamePath = os.path.join(outputDir, jobname, DATA, REPLICATIONTIME)

    ################################################
    if (os.path.exists(jobnamePath)):
        try:
            shutil.rmtree(jobnamePath)
        except OSError as e:
            print('Error: %s - %s.' % (e.filename, e.strerror))
    ################################################

    replicationTimeFilename_woDir = os.path.basename(replicationTimeFilename)

    #Please note that there are 3 options
    #Option1: Read data, prepare chrBased np arrays and load np arrays during runtime. Code below provides option1
    # if (replicationTimeFilename_woDir not in availableLibraryFilenamesList):
    #     readReplicationTimeDataAndWriteChrBasedReplicationTimeNPArrays(genome,chromNamesList,chromSizesDict,replicationTimeFilename)
    #     #append
    #     append2File(replicationTimeFilename_woDir,AVAILABLE_LIBRARY_FILENAMES_PATH)

    #Option2: Load offline prepared np arrays during runtime managby ed by replication_time_np_arrays_fill_runtime=False
    #Option2: Fill np array during runtime managed by replication_time_np_arrays_fill_runtime=True
    replication_time_np_arrays_fill_runtime = True
    computationType = COMPUTATION_ALL_CHROMOSOMES_PARALLEL
    # computationType = COMPUTATION_CHROMOSOMES_SEQUENTIAL_CHROMOSOME_SPLITS_PARALLEL
    # computationType = COMPUTATION_CHROMOSOMES_SEQUENTIAL # Took the longest time
    replicationTimeAnalysis(mutationTypes,computationType,replication_time_np_arrays_fill_runtime,genome,chromSizesDict,chromNamesList,outputDir,jobname,singlePointMutationsFilename,indelsFilename,replicationTimeFilename)
    ###############################################

#######################################################


#######################################################
def runReplicationStrandBiasAnalysis(mutationTypes,singlePointMutationsFilename,indelsFilename,outputDir,jobname,replicationTimeFilename,replicationTimeValleyFilename,replicationTimePeakFilename,chromSizesDict,chromNamesList):

    ###############################################
    # REPLICATIONSTRANDBIAS
    # Delete the output/jobname/DATA/REPLICATIONSTRANDBIAS if exists
    jobnamePath = os.path.join(outputDir,jobname,DATA,REPLICATIONSTRANDBIAS)

    ################################################
    if (os.path.exists(jobnamePath)):
        try:
            shutil.rmtree(jobnamePath)
        except OSError as e:
            print('Error: %s - %s.' % (e.filename, e.strerror))
    ################################################

    smoothedWaveletRepliseqDataFilename = replicationTimeFilename
    valleysBEDFilename = replicationTimeValleyFilename
    peaksBEDFilename = replicationTimePeakFilename

    if (singlePointMutationsFilename!=NOTSET):
        # computationType = COMPUTATION_ALL_CHROMOSOMES_PARALLEL
        computationType = COMPUTATION_CHROMOSOMES_SEQUENTIAL_CHROMOSOME_SPLITS_PARALLEL
        replicationStrandBiasAnalysis(mutationTypes,computationType,chromSizesDict,chromNamesList,outputDir,jobname,singlePointMutationsFilename,indelsFilename,smoothedWaveletRepliseqDataFilename,valleysBEDFilename,peaksBEDFilename)
    ###############################################

#######################################################

#######################################################
def runTranscriptionStradBiasAnalysis(mutationTypes,genome,singlePointMutationsFilename,indelsFilename,outputDir,jobname,chromSizesDict,chromNamesList):
    ###############################################
    # TRANSCRIPTIONSTRANDBIAS
    # Delete the output/jobname/DATA/TRANSCRIPTIONSTRANDBIAS if exists
    jobnamePath = os.path.join(outputDir,jobname,DATA,TRANSCRIPTIONSTRANDBIAS)

    ################################################
    if (os.path.exists(jobnamePath)):
        try:
            shutil.rmtree(jobnamePath)
        except OSError as e:
            print('Error: %s - %s.' % (e.filename, e.strerror))
    ################################################

    if (singlePointMutationsFilename!=NOTSET):
        computationType =COMPUTATION_ALL_CHROMOSOMES_PARALLEL
        # computationType = COMPUTATION_CHROMOSOMES_SEQUENTIAL_CHROMOSOME_SPLITS_PARALLEL
        transcriptionStrandBiasAnalysis(mutationTypes,computationType,genome,chromSizesDict,chromNamesList,outputDir,jobname, singlePointMutationsFilename,indelsFilename)
    ###############################################
#######################################################


#######################################################
def runProcessivityAnalysis(singlePointMutationsFilename,outputDir,jobname,chromNamesList):
    ###############################################
    # PROCESSIVITY
    # Delete the output/jobname/DATA/PROCESSIVITY if exists
    jobnamePath = os.path.join(outputDir,jobname,DATA,PROCESSIVITY)

    ###############################################
    if (os.path.exists(jobnamePath)):
        try:
            shutil.rmtree(jobnamePath)
        except OSError as e:
            print('Error: %s - %s.' % (e.filename, e.strerror))
    ###############################################

    #Internally Set
    considerProbabilityInProcessivityAnalysis = True

    # subprocess.call(['python', os.path.join(current_abs_path,SOURCE,PROCESSIVITY, 'ProcessivityAnalysis.py'),jobname,singlePointMutationsFilename,considerProbability])
    if (singlePointMutationsFilename != NOTSET):
        processivityAnalysis(chromNamesList,outputDir,jobname,singlePointMutationsFilename,considerProbabilityInProcessivityAnalysis)
    ###############################################

#######################################################


#######################################################
#Run SigProfilerTopography Analyses
def runAnalyses(genome,singlePointMutationsFilename,indelsFilename,outputDir,jobname,nucleosomeFilename,replicationTimeFilename,replicationTimeValleyFilename,replicationTimePeakFilename,mutationTypes=[SUBS, INDELS, DINUCS]):

    chromSizesDict = getChromSizesDict(genome)
    chromNamesList = list(chromSizesDict.keys())
    availableLibraryFilenamesList = getAvailableLibraryFilenamesList()

    print('Genome: %s' %(genome))
    print('Chromosome names: %s' %(chromNamesList))
    current_abs_path = os.path.dirname(os.path.realpath(__file__))
    print('current_abs_path: %s ' % current_abs_path)

    ##############################################
    #Partition the data (Single Point Mutations data and Indels data)
    # Delete the output/jobname/DATA/chrbased if exists
    jobnamePath = os.path.join(outputDir,jobname,DATA,CHRBASED)

    ######################################################
    # SigProfilerMatrixGenerator will provide chrbased subs, indels and dinucs mutation files.
    # Let's not delete the existing data
    if (os.path.exists(jobnamePath)):
        try:
            shutil.rmtree(jobnamePath)
        except OSError as e:
            print('Error: %s - %s.' % (e.filename, e.strerror))
    ######################################################

    if (indelsFilename!=NOTSET):
        partitionIndelsData(genome,outputDir,jobname,indelsFilename)
    if (singlePointMutationsFilename!=NOTSET):
        partitionSubsData(outputDir,jobname,singlePointMutationsFilename)
    #############################################


    ##############################################
    #Please note that singlePointMutationsFilename and indelsFilename are full paths.
    # They are all read and partitioned chr based
    # Get the fienames at the end
    singlePointMutationsFilename = os.path.basename(singlePointMutationsFilename)
    indelsFilename = os.path.basename(indelsFilename)
    ##############################################

    ##########################################################################################
    fill_dinucs_dictionaries_write(outputDir, jobname, chromNamesList)
    ##########################################################################################

    #############################################
    runNucleosomeOccupancyAnalyses(mutationTypes,genome,singlePointMutationsFilename,indelsFilename,outputDir,jobname,nucleosomeFilename,chromSizesDict,chromNamesList,availableLibraryFilenamesList)
    runReplicationTimeAnalysis(mutationTypes,genome,singlePointMutationsFilename,indelsFilename,outputDir,jobname,replicationTimeFilename,chromSizesDict,chromNamesList)
    runReplicationStrandBiasAnalysis(mutationTypes,singlePointMutationsFilename,indelsFilename,outputDir,jobname,replicationTimeFilename,replicationTimeValleyFilename,replicationTimePeakFilename,chromSizesDict,chromNamesList)
    runTranscriptionStradBiasAnalysis(mutationTypes,genome,singlePointMutationsFilename,indelsFilename,outputDir,jobname,chromSizesDict,chromNamesList)
    runProcessivityAnalysis(singlePointMutationsFilename,outputDir,jobname,chromNamesList)
    #############################################

#######################################################



##############################################################
# BONFERRONI_CORRECTION = 'BONFERRONI_CORRECTION'
# FDR_BH_CORRECTION = 'FDR_BH_CORRECTION'
#
# USING_POISSON_DISTRIBUTION = 'USING_POISSON_DISTRIBUTION'
# USING_NULL_DISTRIBUTION = 'USING_NULL_DISTRIBUTION'
# USING_GAUSSIAN_KDE = 'USING_GAUSSIAN_KDE'

#Plot Figures for the attainded data after SigProfilerTopography Analyses
def plotFigures(outputDir,jobname,numberofSimulations,multipleTesting,probabilityCalculation,mutationTypes=[SUBS, INDELS, DINUCS]):

    #Internally Set
    figureAugmentation = 'noaugmentation'

    current_abs_path = os.path.dirname(os.path.realpath(__file__))

    jobnamePath = os.path.join(current_abs_path,OUTPUT,jobname,FIGURE)
    print('Topography.py jobnamePath:%s ' %jobnamePath)
    if (os.path.exists(jobnamePath)):
        print('jobnamePath exists')

    ############################################################
    if (os.path.exists(jobnamePath)):
        try:
            shutil.rmtree(jobnamePath)
        except OSError as e:
            print('Error: %s - %s.' % (e.filename, e.strerror))
    ############################################################

    ############################################################
    nucleosomeOccupancyAverageSignalFigures(outputDir,jobname,figureAugmentation,numberofSimulations,mutationTypes)
    replicationTimeNormalizedMutationDensityFigures(outputDir,jobname,figureAugmentation,numberofSimulations,mutationTypes)
    transcriptionReplicationStrandBiasFigures(outputDir,jobname,figureAugmentation,numberofSimulations)
    processivityFigures(outputDir,jobname,numberofSimulations,multipleTesting,probabilityCalculation)
    ############################################################

##############################################################


# ##############################################################
# #To run local on laptop
# matrix_generator_output_dir_path = '/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/input_for_matgen/BreastCancer560_subs_indels_dinucs/output/vcf_files/DINUC'
# dinucs_probabilities_file_path = '/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/output/560_BRCA_WGS_DINUCS/DBS78/Suggested_Solution/Decomposed_Solution/Mutation_Probabilities.txt'
# prepareDinucsDataAfterMatrixGenerationAndExtractorForTopography(matrix_generator_output_dir_path,dinucs_probabilities_file_path)
# ##############################################################
