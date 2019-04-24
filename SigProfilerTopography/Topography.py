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

import SigProfilerMatrixGenerator as matgen_package
print(matgen_package.__path__[0])
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

from SigProfilerTopography.source.commons.DataPreparationCommons import readMutationsWithGenomicPositions
from SigProfilerTopography.source.commons.DataPreparationCommons import readProbabilities
from SigProfilerTopography.source.commons.DataPreparationCommons import mergeSNPsWithSignatureProbabilities
from SigProfilerTopography.source.commons.DataPreparationCommons import readIndelsandWriteWithGenomicPositions
from SigProfilerTopography.source.commons.DataPreparationCommons import prepareSimulationBasedInputFilesForSigProfilerTopography
from SigProfilerTopography.source.commons.DataPreparationCommons import readChrBasedMutationsMergeWithProbabilitiesAndWrite

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
def prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(outputDir, jobname,mutationType, matrix_generator_output_dir_path,mutations_probabilities_file_path):
    os.makedirs(os.path.join(outputDir,jobname,DATA,CHRBASED),exist_ok=True)
    chrShortList = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
    if (os.path.exists(mutations_probabilities_file_path)):
        mutations_probabilities_df = readProbabilities(mutations_probabilities_file_path)
        # dinucs_probabilities_df.columns.names [Sample Names    MutationTypes   DBS2    DBS4    DBS6    DBS7    DBS11]
        # dinucs_probabilities_df.columns.names [Sample    Mutation   DBS2    DBS4    DBS6    DBS7    DBS11]
        mutations_probabilities_df.rename(columns={'Sample Names': 'Sample', 'MutationTypes': 'Mutation'}, inplace=True)

    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)

    poolInputList = []
    if (os.path.exists(matrix_generator_output_dir_path)):
        for chrShort in chrShortList:
            chr_based_mutation_filename = '%s_seqinfo.txt' %(chrShort)
            chr_based_mutation_filepath = os.path.join(matrix_generator_output_dir_path,chr_based_mutation_filename)
            inputList=[]
            inputList.append(chrShort)
            inputList.append(outputDir)
            inputList.append(jobname)
            inputList.append(chr_based_mutation_filepath)
            inputList.append(mutations_probabilities_df)
            inputList.append(mutationType)
            poolInputList.append(inputList)

        pool.map(readChrBasedMutationsMergeWithProbabilitiesAndWrite,poolInputList)

    pool.close()
    pool.join()
############################################################


############################################################
#TODO To be depreceated. To be deleted.
#TODO All the auxiliary functions called from this function will be depreceated and then be deleted
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
def runNucleosomeOccupancyAnalyses(mutationTypes,genome,outputDir,jobname,nucleosomeFilename,chromSizesDict,chromNamesList,availableLibraryFilenamesList):
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

    nucleosomeFilename_wo_dir = os.path.basename(nucleosomeFilename)

    if (nucleosomeFilename_wo_dir not in availableLibraryFilenamesList):
        quantileValue = round(float(0.97), 2)
        # PartitionNucleosomeOccupancyData.partitionNucleosomeOccupancyData(jobname,nucleosomeFilename,quantileValue)
        readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArrays(genome,quantileValue,nucleosomeFilename)
        #append
        append2File(nucleosomeFilename_wo_dir,AVAILABLE_LIBRARY_FILENAMES_PATH)

    # computationType = COMPUTATION_ALL_CHROMOSOMES_PARALLEL
    computationType = COMPUTATION_CHROMOSOMES_SEQUENTIAL
    nucleosomeOccupancyAnalysis(mutationTypes,computationType,chromSizesDict,chromNamesList,outputDir,jobname,nucleosomeFilename)
    ###############################################
#######################################################


#######################################################
def runReplicationTimeAnalysis(mutationTypes,genome,outputDir,jobname,replicationTimeFilename,chromSizesDict,chromNamesList):
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
    replicationTimeAnalysis(mutationTypes,computationType,replication_time_np_arrays_fill_runtime,genome,chromSizesDict,chromNamesList,outputDir,jobname,replicationTimeFilename)
    ###############################################

#######################################################


#######################################################
def runReplicationStrandBiasAnalysis(mutationTypes,outputDir,jobname,replicationTimeFilename,replicationTimeValleyFilename,replicationTimePeakFilename,chromSizesDict,chromNamesList):

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

    # computationType = COMPUTATION_ALL_CHROMOSOMES_PARALLEL
    computationType = COMPUTATION_CHROMOSOMES_SEQUENTIAL_CHROMOSOME_SPLITS_PARALLEL
    replicationStrandBiasAnalysis(mutationTypes,computationType,chromSizesDict,chromNamesList,outputDir,jobname,smoothedWaveletRepliseqDataFilename,valleysBEDFilename,peaksBEDFilename)
    ###############################################

#######################################################

#######################################################
def runTranscriptionStradBiasAnalysis(mutationTypes,genome,outputDir,jobname,chromSizesDict,chromNamesList):
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

    computationType =COMPUTATION_ALL_CHROMOSOMES_PARALLEL
    # computationType = COMPUTATION_CHROMOSOMES_SEQUENTIAL_CHROMOSOME_SPLITS_PARALLEL
    useTranscriptionStrandColumn = False
    transcriptionStrandBiasAnalysis(mutationTypes,computationType,useTranscriptionStrandColumn,genome,chromSizesDict,chromNamesList,outputDir,jobname)
    ###############################################
#######################################################


#######################################################
def runProcessivityAnalysis(mutationTypes,outputDir,jobname,chromNamesList):
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

    processivityAnalysis(mutationTypes,chromNamesList,outputDir,jobname,considerProbabilityInProcessivityAnalysis)
    ###############################################

#######################################################


#######################################################
# inputDir ='/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/input_for_matgen/BreastCancer560_subs_indels_dinucs'
# outputDir = '/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/output_test/'
# jobname = 'BreastCancer560'

#Run SigProfilerTopography Analyses
#Former full path now only the filename with extension
# nucleosomeOccupancy = '/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/lib/nucleosome/wgEncodeSydhNsomeGm12878Sig.wig'
# replicationSignal = '/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.wig'
# replicationValley = '/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7ValleysRep1.bed'
# replicationPeak = '/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7PkRep1.bed'
# subs_probabilities_file_path = '/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/output/560_BRCA_WGS_DINUCS/SBS96/Suggested_Solution/Decomposed_Solution/Mutation_Probabilities.txt'
# indels_probabilities_file_path = '/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/output/560_BRCA_WGS_DINUCS/ID83/Suggested_Solution/Decomposed_Solution/Mutation_Probabilities.txt'
# dinucs_probabilities_file_path = '/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/output/560_BRCA_WGS_DINUCS/DBS78/Suggested_Solution/Decomposed_Solution/Mutation_Probabilities.txt'
def runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,subs_probabilities_file_path,indels_probabilities_file_path,dinucs_probabilities_file_path,nucleosomeFilename=DEFAULT_NUCLEOSOME_OCCUPANCY_FILE,replicationTimeFilename=DEFAULT_REPLICATION_TIME_SIGNAL_FILE,replicationTimeValleyFilename=DEFAULT_REPLICATION_TIME_VALLEY_FILE,replicationTimePeakFilename=DEFAULT_REPLICATION_TIME_PEAK_FILE,mutationTypes=[SUBS, INDELS, DINUCS]):
    chromSizesDict = getChromSizesDict(genome)
    chromNamesList = list(chromSizesDict.keys())
    availableLibraryFilenamesList = getAvailableLibraryFilenamesList()

    print('##################### SigProfilerTopography parameters ##########################')
    print('Genome: %s' %(genome))
    print('inputDir:%s' %inputDir)
    print('outputDir:%s' %outputDir)
    print('jobname:%s' %jobname)
    print('nucleosomeFilename:%s' %nucleosomeFilename)
    print('replicationTimeFilename:%s' % replicationTimeFilename)
    print('replicationTimeValleyFilename:%s' % replicationTimeValleyFilename)
    print('replicationTimePeakFilename:%s' % replicationTimePeakFilename)
    print('mutationTypes:%s' %mutationTypes)
    print('#################################################################################')

    #################################################################################
    print('Chromosome names: %s' %(chromNamesList))
    current_abs_path = os.path.dirname(os.path.realpath(__file__))
    print('current_abs_path: %s ' % current_abs_path)
    #################################################################################

    #################################################################################
    #######################  SigProfilerMatrixGenerator starts ######################
    #################################################################################
    #Prepare chrBased input files
    matGen.SigProfilerMatrixGeneratorFunc(jobname,genome,inputDir,plot=True, seqInfo=True)
    #################################################################################
    #######################  SigProfilerMatrixGenerator ends ########################
    #################################################################################


    #################################################################################
    ##################  Merge Files with Mutation Probabilities starts ##############
    #################################################################################
    #SUBS
    if (SUBS in mutationTypes):
        matrix_generator_output_dir_path = os.path.join(inputDir,'output','vcf_files','SNV')
        prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(outputDir,jobname,SUBS,matrix_generator_output_dir_path,subs_probabilities_file_path)

    #INDELS
    if (INDELS in mutationTypes):
        matrix_generator_output_dir_path =  os.path.join(inputDir,'output','vcf_files','ID')
        prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(outputDir,jobname,INDELS,matrix_generator_output_dir_path,indels_probabilities_file_path)

    #DINUCS
    if (DINUCS in mutationTypes):
        matrix_generator_output_dir_path = os.path.join(inputDir,'output','vcf_files','DBS')
        prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(outputDir,jobname,DINUCS,matrix_generator_output_dir_path,dinucs_probabilities_file_path)
    #################################################################################
    ##################  Merge Files with Mutation Probabilities ends ################
    #################################################################################

    #################################################################################
    fill_mutations_dictionaries_write(outputDir, jobname, chromNamesList,SUBS)
    fill_mutations_dictionaries_write(outputDir, jobname, chromNamesList,INDELS)
    fill_mutations_dictionaries_write(outputDir, jobname, chromNamesList,DINUCS)
    #################################################################################

    #################################################################################
    if (nucleosomeFilename == DEFAULT_NUCLEOSOME_OCCUPANCY_FILE):
        nucleosomeFilename = os.path.join(current_abs_path,LIB,NUCLEOSOMEOCCUPANCY,DEFAULT_NUCLEOSOME_OCCUPANCY_FILE)

    if (replicationTimeFilename == DEFAULT_REPLICATION_TIME_SIGNAL_FILE):
        replicationTimeFilename = os.path.join(current_abs_path,LIB,REPLICATION,DEFAULT_REPLICATION_TIME_SIGNAL_FILE)

    if (replicationTimeValleyFilename == DEFAULT_REPLICATION_TIME_VALLEY_FILE):
        replicationTimeValleyFilename = os.path.join(current_abs_path,LIB,REPLICATION,DEFAULT_REPLICATION_TIME_VALLEY_FILE)
    if (replicationTimePeakFilename == DEFAULT_REPLICATION_TIME_PEAK_FILE):
        replicationTimePeakFilename = os.path.join(current_abs_path,LIB,REPLICATION,DEFAULT_REPLICATION_TIME_PEAK_FILE)
    #################################################################################


    #################################################################################
    runNucleosomeOccupancyAnalyses(mutationTypes,genome,outputDir,jobname,nucleosomeFilename,chromSizesDict,chromNamesList,availableLibraryFilenamesList)
    runReplicationTimeAnalysis(mutationTypes,genome,outputDir,jobname,replicationTimeFilename,chromSizesDict,chromNamesList)
    runReplicationStrandBiasAnalysis(mutationTypes,outputDir,jobname,replicationTimeFilename,replicationTimeValleyFilename,replicationTimePeakFilename,chromSizesDict,chromNamesList)
    runTranscriptionStradBiasAnalysis(mutationTypes,genome,outputDir,jobname,chromSizesDict,chromNamesList)
    runProcessivityAnalysis(mutationTypes,outputDir,jobname,chromNamesList)
    #################################################################################

    #################################################################################
    plotFigures(outputDir, jobname, numofSimulations, 'BONFERRONI_CORRECTION', 'USING_POISSON_DISTRIBUTION')
    #################################################################################

#######################################################



##############################################################
# BONFERRONI_CORRECTION = 'BONFERRONI_CORRECTION'
# FDR_BH_CORRECTION = 'FDR_BH_CORRECTION'
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
