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

import time
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerSimulator import SigProfilerSimulator as simulator

from SigProfilerTopography.source.commons.DataPreparationCommons import readProbabilities
from SigProfilerTopography.source.commons.DataPreparationCommons import readChrBasedMutationsMergeWithProbabilitiesAndWrite

from SigProfilerTopography.source.commons.TopographyCommons import *

from SigProfilerTopography.source.nucleosomeoccupancy.NucleosomeOccupancyAnalysis import occupancyAnalysis
from SigProfilerTopography.source.replicationtime.ReplicationTimeAnalysis import replicationTimeAnalysis

from SigProfilerTopography.source.replicationstrandbias.ReplicationStrandBiasAnalysis import replicationStrandBiasAnalysis
from SigProfilerTopography.source.transcriptionstrandbias.TranscriptionStrandBiasAnalysis import transcriptionStrandBiasAnalysis
from SigProfilerTopography.source.processivity.ProcessivityAnalysis import processivityAnalysis

from SigProfilerTopography.source.plotting.NucleosomeOccupancyAverageSignalFigures import occupancyAverageSignalFigures
from SigProfilerTopography.source.plotting.ReplicationTimeNormalizedMutationDensityFigures import replicationTimeNormalizedMutationDensityFigures
from SigProfilerTopography.source.plotting.TranscriptionReplicationStrandBiasFigures import transcriptionReplicationStrandBiasFigures
from SigProfilerTopography.source.plotting.ProcessivityFigures import processivityFigures

import subprocess
import shutil
import logging

############################################################
#Can be move to DataPreparationCommons under /source/commons
#read chr based dinucs (provided by SigProfilerMatrixGenerator) and merge with probabilities (provided by SigProfilerTopography)
def prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList,inputDir,outputDir,jobname,mutation_type_context,mutations_probabilities_file_path,startSimNum, endSimNum,partialDirname,logger):

    #original matrix generator chrbased data will be under inputDir/output/vcf_files/SNV
    #original matrix generator chrbased data will be under inputDir/output/vcf_files/DBS
    #original matrix generator chrbased data will be under inputDir/output/vcf_files/ID

    #sim1 matrix generator chrbased data will be under inputDir/output/simulations/sim1/96/output/vcf_files/SNV
    #sim1 matrix generator chrbased data will be under inputDir/output/simulations/sim1/ID/output/vcf_files/ID
    #sim1 matrix generator chrbased data will be under inputDir/output/simulations/sim1/DBS/output/vcf_files/DBS

    os.makedirs(os.path.join(outputDir,jobname,DATA,CHRBASED),exist_ok=True)
    for simNum in range(1,endSimNum+1):
        simName = 'sim%d' % (simNum)
        os.makedirs(os.path.join(outputDir, jobname, DATA,CHRBASED,simName), exist_ok=True)

    if (os.path.exists(mutations_probabilities_file_path)):
        ##########################################################################################
        #Step1
        mutations_probabilities_df = readProbabilities(mutations_probabilities_file_path,logger)
        ##########################################################################################

        # print('mutations_probabilities_df.head()')
        # print(mutations_probabilities_df.head())

        # print('mutations_probabilities_df.columns.values')
        # print(mutations_probabilities_df.columns.values)

        ##########################################################################################
        #Step2
        #For release we will use SAMPLE as it is, no change in SAMPLE column is needed.

        #For PCAWG_Matlab
        # This statement below will be unnecessary
        # This is customized for  PCAWG_Matlab
        # To get rid of inconsistent cancer type names in sample column of chrbased mutation files and probabilities files
        # This behaviour can be parameterized
        # Breast-LobularCA_SP124191
        # mutations_probabilities_df[SAMPLE] = mutations_probabilities_df[SAMPLE].str.split('_',expand=True)[1]
        ##########################################################################################

        numofProcesses = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(numofProcesses)

        poolInputList = []

        for simNum in range(startSimNum,endSimNum+1):
            simName = 'sim%d' %(simNum)
            for chrShort in chromShortNamesList:
                chr_based_mutation_filename = '%s_seqinfo.txt' % (chrShort)
                if (simNum==0):
                    matrix_generator_output_dir_path = os.path.join(inputDir,'output','vcf_files',partialDirname)
                else:
                    matrix_generator_output_dir_path = os.path.join(inputDir,'output','simulations',simName,mutation_type_context,'output','vcf_files',partialDirname)

                if (os.path.exists(matrix_generator_output_dir_path)):
                    chr_based_mutation_filepath = os.path.join(matrix_generator_output_dir_path,chr_based_mutation_filename)
                    inputList = []
                    inputList.append(chrShort)
                    inputList.append(outputDir)
                    inputList.append(jobname)
                    inputList.append(chr_based_mutation_filepath)
                    inputList.append(mutations_probabilities_df)
                    inputList.append(mutation_type_context)
                    inputList.append(simNum)
                    poolInputList.append(inputList)

        pool.map(readChrBasedMutationsMergeWithProbabilitiesAndWrite, poolInputList)

        pool.close()
        pool.join()

    else:
        print('%s does not exist.' %(mutations_probabilities_file_path))
############################################################


#######################################################
def download_2bit_file(genome):
    if (genome == GRCh37):
        os.makedirs(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, UCSCGENOME), exist_ok=True)
        filepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, UCSCGENOME, HG19_2BIT)
        downloadFromWeb(HG19_URL, filepath)
    elif (genome == GRCh38):
        os.makedirs(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, UCSCGENOME), exist_ok=True)
        filepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, UCSCGENOME, HG38_2BIT)
        downloadFromWeb(HG38_URL, filepath)
    elif (genome == MM9):
        os.makedirs(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, UCSCGENOME), exist_ok=True)
        filepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, UCSCGENOME, MM9_2BIT)
        downloadFromWeb(MM9_URL, filepath)
    elif (genome == MM10):
        os.makedirs(os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, UCSCGENOME), exist_ok=True)
        filepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, UCSCGENOME, MM10_2BIT)
        downloadFromWeb(MM10_URL, filepath)
#######################################################

# #######################################################
# def download_bigwig2wig():
#     filepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME,BIGWIG2WIG)
#     downloadFromWeb(BIGWIG_TO_WIG_EXECUTABLE_LINUX_X86_64_URL,filepath)
#     os.chmod(filepath,0o744)
# #######################################################

########################################################
# bigWig2Wig executable is for linux/unix
# https://hgdownload.cse.ucsc.edu/admin/exe/
# At this address mac version is also provided but not for windows
def download_nucleosome_occupancy(cellLine):
    # bigWig2Wig_filepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME, BIGWIG2WIG)
    # os.chmod(bigWig2Wig_filepath,0o744)
    if (cellLine==GM12878):
        gm12878_bigWig_filepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME, ENCODE_NUCLEOSOME_GM12878_BIGWIG)
        downloadFromWeb(ENCODE_NUCLEOSOME_GM12878_BIGWIG_URL, gm12878_bigWig_filepath)
        # gm12878_wig_filepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME,ENCODE_NUCLEOSOME_GM12878_WIG)
        # subprocess.call([bigWig2Wig_filepath, gm12878_bigWig_filepath,gm12878_wig_filepath])
    elif (cellLine==K562):
        K562_bigWig_filepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME, ENCODE_NUCLEOSOME_K562_BIGWIG)
        downloadFromWeb(ENCODE_NUCLEOSOME_K562_BIGWIG_URL, K562_bigWig_filepath)
        # K562_wig_filepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME, ENCODE_NUCLEOSOME_K562_WIG)
        # subprocess.call([bigWig2Wig_filepath, K562_bigWig_filepath,K562_wig_filepath])
#######################################################


#######################################################
def runOccupancyAnalyses(outputDir,jobname,numofSimulations,sample_based,library_file_with_path,library_file_memo,chromSizesDict,chromNamesList,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,computation_type,occupancy_type,plusorMinus):

    # Check whether nucleosomeFilename_wo_dir is downloaded
    # if (os.path.exists(nucleosomeFilename_with_dir)):
    #     quantileValue = round(float(0.97), 2)
    #     package_lib_nucleosome_directory_path = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,NUCLEOSOME,CHRBASED)
    #
    #     if (not os.path.exists(package_lib_nucleosome_directory_path)):
    #         print('%s does not exists' %(package_lib_nucleosome_directory_path))
    #         readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArraysSequentially(genome,quantileValue,nucleosomeFilename_with_dir)
    #         # readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArrays(genome,quantileValue,nucleosomeFilename)
    #     elif (not os.listdir(package_lib_nucleosome_directory_path)):
    #         print('There is no file under: %s' %(package_lib_nucleosome_directory_path))
    #         readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArraysSequentially(genome,quantileValue,nucleosomeFilename_with_dir)
    #         # readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArrays(genome,quantileValue,nucleosomeFilename)

    if (not os.path.exists(library_file_with_path)):
        # nucleosomeFilename_wo_dir = os.path.basename(nucleosomeFilename_with_dir)
        # download_nucleosome_command = 'download_nucleosome_occupancy_convert_bigWig2wig(cellLine) command'
        print('There is no such file under %s' %(library_file_with_path))

    # computationType = COMPUTATION_CHROMOSOMES_SEQUENTIAL_ALL_SIMULATIONS_PARALLEL
    # computationType = COMPUTATION_CHROMOSOMES_SEQUENTIAL_SIMULATIONS_SEQUENTIAL
    # computationType = COMPUTATION_ALL_CHROMOSOMES_PARALLEL
    occupancyAnalysis(computation_type,occupancy_type,sample_based,plusorMinus,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,library_file_with_path,library_file_memo,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict)
    ###############################################
#######################################################


#######################################################
def runReplicationTimeAnalysis(genome,outputDir,jobname,numofSimulations,sample_based,replicationTimeFilename,chromSizesDict,chromNamesList,computation_type,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict):
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
    # computationType = COMPUTATION_ALL_CHROMOSOMES_PARALLEL
    # computationType = COMPUTATION_CHROMOSOMES_SEQUENTIAL_SIMULATIONS_SEQUENTIAL
    # computationType = COMPUTATION_CHROMOSOMES_SEQUENTIAL_ALL_SIMULATIONS_PARALLEL
    # computationType = COMPUTATION_CHROMOSOMES_SEQUENTIAL_CHROMOSOME_SPLITS_PARALLEL
    replicationTimeAnalysis(computation_type,sample_based,replication_time_np_arrays_fill_runtime,genome,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,replicationTimeFilename,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict)
    ###############################################

#######################################################


#######################################################
def runReplicationStrandBiasAnalysis(outputDir,jobname,numofSimulations,sample_based,replicationTimeFilename,replicationTimeValleyFilename,replicationTimePeakFilename,chromSizesDict,chromNamesList,computation_type,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict):

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
    # computationType = COMPUTATION_CHROMOSOMES_SEQUENTIAL_CHROMOSOME_SPLITS_PARALLEL
    # computationType = COMPUTATION_CHROMOSOMES_SEQUENTIAL_ALL_SIMULATIONS_PARALLEL
    # computationType = COMPUTATION_CHROMOSOMES_SEQUENTIAL_SIMULATIONS_SEQUENTIAL
    replicationStrandBiasAnalysis(computation_type,sample_based,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,smoothedWaveletRepliseqDataFilename,valleysBEDFilename,peaksBEDFilename,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict)
    ###############################################

#######################################################

#######################################################
def runTranscriptionStradBiasAnalysis(genome,outputDir,jobname,numofSimulations,sample_based,chromSizesDict,chromNamesList,computation_type,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict):
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

    # computation_type = COMPUTATION_ALL_CHROMOSOMES_PARALLEL
    # computation_type = COMPUTATION_CHROMOSOMES_SEQUENTIAL_CHROMOSOME_SPLITS_PARALLEL
    # computation_type = COMPUTATION_CHROMOSOMES_SEQUENTIAL_ALL_SIMULATIONS_PARALLEL
    # computation_type = COMPUTATION_CHROMOSOMES_SEQUENTIAL_SIMULATIONS_SEQUENTIAL
    # useTranscriptionStrandColumn = False
    useTranscriptionStrandColumn = True
    transcriptionStrandBiasAnalysis(computation_type,sample_based,useTranscriptionStrandColumn,genome,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict)
    ###############################################
#######################################################


#######################################################
def runProcessivityAnalysis(mutation_types_contexts,outputDir,jobname,numofSimulations,chromNamesList,signature2PropertiesListDict):
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
    # considerProbabilityInProcessivityAnalysis = False

    processivityAnalysis(mutation_types_contexts,chromNamesList,outputDir,jobname,numofSimulations,considerProbabilityInProcessivityAnalysis,signature2PropertiesListDict)
    ###############################################
#######################################################


#######################################################
def deleteOldData(outputDir,jobname,occupancy_type):
    #############################################
    # Delete the output/jobname/DATA/occupancy_type if exists
    jobnamePath = os.path.join(outputDir,jobname,DATA,occupancy_type)

    ################################################
    if (os.path.exists(jobnamePath)):
        try:
            shutil.rmtree(jobnamePath)
        except OSError as e:
            print('Error: %s - %s.' % (e.filename, e.strerror))
    ################################################
#######################################################

#######################################################
def deleteOldFigures(outputDir, jobname, occupancy_type):

    jobnamePath = os.path.join(outputDir, jobname, FIGURE, ALL, occupancy_type)
    print('Topography.py jobnamePath:%s ' %jobnamePath)

    ############################################################
    if (os.path.exists(jobnamePath)):
        try:
            shutil.rmtree(jobnamePath)
        except OSError as e:
            print('Error: %s - %s.' % (e.filename, e.strerror))
    ############################################################
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
def runAnalyses(genome,
                inputDir,
                outputDir,
                jobname,
                numofSimulations,
                sbs_probabilities_file_path=None,
                id_probabilities_file_path= None,
                dbs_probabilities_file_path=None,
                epigenomics_files=[DEFAULT_HISTONE_OCCUPANCY_FILE],
                epigenomics_files_memos=['Liver_H3K9me3'],
                nucleosome_file=DEFAULT_NUCLEOSOME_OCCUPANCY_FILE,
                replication_time_file=DEFAULT_REPLICATION_TIME_SIGNAL_FILE,
                replication_time_valley_file=DEFAULT_REPLICATION_TIME_VALLEY_FILE,
                replication_time_peak_file=DEFAULT_REPLICATION_TIME_PEAK_FILE,
                mutation_types_contexts=[SBS96,ID,DBS],
                computation_type=COMPUTATION_CHROMOSOMES_SEQUENTIAL_ALL_SIMULATIONS_PARALLEL,
                epigenomics=False,
                nucleosome=False,
                replication_time=False,
                strand_bias=False,
                processivity=False,
                sample_based=False,
                new_simulations_enforced=False,
                plot_figures=True):

    # ucsc hg19 chromosome names:
    # 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chrX', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr20', 'chrY', 'chr19', 'chr22', 'chr21', 'chrM'

    # ensembl GRCh37 chromosome names:
    # '1', '2', '3', '4', '5', '6', '7', 'X', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '20', 'Y', '19', '22', '21', 'MT'

    # default library files (nucleosome occupancy and replication time) are all in hg19
    # hg19 wgEncodeSydhNsomeGm12878Sig.bigWig from http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeSydhNsome
    # hg19 SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.wig
    # hg19 SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7PkRep1.bed
    # hg19 SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7ValleysRep1.bed

    chromSizesDict = getChromSizesDict(genome)
    chromNamesList = list(chromSizesDict.keys())
    chromShortNamesList=getShortNames(chromNamesList)

    ############################################################################################
    os.makedirs(os.path.join(outputDir,jobname),exist_ok=True)
    log_file=os.path.join(outputDir,jobname,'SigProfilerTopography.log')

    # set up logging to file - see previous section for more details
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename=log_file,
                        filemode='w')

    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)

    # set a format which is simpler for console use
    # formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    formatter = logging.Formatter('%(message)s')

    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)

    logger = logging.getLogger('SigProfilerTopography')
    ############################################################################################

    logger.info('#################################################################################')
    logger.info("--- SigProfilerTopography starts")
    logger.info('#################################################################################\n')

    logger.info('#################################################################################')
    logger.info('--- SigProfilerTopography parameters')
    logger.info('--- Genome: %s' %(genome))
    logger.info('--- inputDir:%s' %inputDir)
    logger.info('--- outputDir:%s' %outputDir)
    logger.info('--- jobname:%s' %jobname)
    logger.info('--- numofSimulations:%d' %numofSimulations)
    logger.info('--- epigenomics_files:%s' %epigenomics_files)
    logger.info('--- epigenomics_files_memos:%s' %epigenomics_files_memos)
    logger.info('--- nucleosome_file:%s' %nucleosome_file)
    logger.info('--- replication_time_file:%s' % replication_time_file)
    logger.info('--- replication_time_valley_file:%s' % replication_time_valley_file)
    logger.info('--- replication_time_peak_file:%s' % replication_time_peak_file)
    logger.info('--- mutation_types_contexts:%s' %mutation_types_contexts)
    logger.info('--- computation_type:%s' %computation_type)
    if sample_based:
        if epigenomics:
            logger.info('--- Epigenomics Sample Based Analysis.')
        if nucleosome:
            logger.info('--- Nucleosome Sample Based Analysis.')
        if replication_time:
            logger.info('--- Replication Time Sample Based Analysis.')
        if strand_bias:
            logger.info('--- Strand Bias Sample Based Analysis.')
        if processivity:
            logger.info('--- Processivity Analysis.')
    else:
        if epigenomics:
            logger.info('--- Epigenomics Analysis.')
        if nucleosome:
            logger.info('--- Nucleosome Analysis.')
        if replication_time:
            logger.info('--- Replication Time Analysis.')
        if strand_bias:
            print('--- Strand Bias Analysis.')
        if processivity:
            logger.info('--- Processivity Analysis.')
    logger.info('--- new_simulations_enforced:%s' %new_simulations_enforced)
    logger.info('--- plot_figures:%s' %plot_figures)
    logger.info('#################################################################################\n')

    #################################################################################
    logger.info('#################################################################################')
    logger.info('--- For Genome: %s' %(genome))
    logger.info('--- Chromosome names: %s' %(chromNamesList))
    logger.info('--- Chromosome short names: %s' % (chromShortNamesList))
    current_abs_path = os.path.dirname(os.path.realpath(__file__))
    logger.info('--- current_abs_path: %s ' % current_abs_path)
    logger.info('#################################################################################\n')
    #################################################################################


    ###################################################################################################
    #######################  SigProfilerMatrixGenerator for original data starts ######################
    ###################################################################################################
    #Run MatrixGenerator for original data: this call prepares chrBased input files for original data with mutation contexts
    logger.info('#################################################################################')
    logger.info('--- SigProfilerMatrixGenerator for original data')
    start_time = time.time()

    logger.info('For original data inputDir:%s' % (inputDir))
    matrices = matGen.SigProfilerMatrixGeneratorFunc(jobname,genome,inputDir,plot=False, seqInfo=True)
    # print('matrices')
    # print(matrices)

    # original matrix generator chrbased data will be under inputDir/output/vcf_files/SNV
    # original matrix generator chrbased data will be under inputDir/output/vcf_files/DBS
    # original matrix generator chrbased data will be under inputDir/output/vcf_files/ID

    logger.info("--- SigProfilerMatrixGenerator for original data: %s seconds ---" % (time.time() -  start_time))
    logger.info("--- SigProfilerMatrixGenerator for original data: %f minutess ---" %float((time.time() - start_time)/60))
    logger.info('#################################################################################\n')
    ###################################################################################################
    #######################  SigProfilerMatrixGenerator for original data ends ########################
    ###################################################################################################

    ####################################################################################################################
    ##################  Merge original chr based files with Mutation Probabilities starts ##############################
    ####################################################################################################################
    logger.info('#################################################################################')
    logger.info('--- Merge original chr based files with Mutation Probabilities starts')
    logger.info('#################################################################################')
    startSimNum = 0
    endSimNum = 0
    start_time = time.time()
    # SBS
    for mutation_type_context in mutation_types_contexts:
        if (mutation_type_context in SBS_CONTEXTS) and (sbs_probabilities_file_path is not None):
            logger.info('--- Merge %s mutations with probabilities for %s' %(mutation_type_context, sbs_probabilities_file_path))
            prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList, inputDir, outputDir,
                                                                               jobname, mutation_type_context,
                                                                               sbs_probabilities_file_path, startSimNum,
                                                                               endSimNum, 'SNV',logger)

    # ID
    if ((ID in mutation_types_contexts) and (id_probabilities_file_path is not None)):
        logger.info('--- Merge %s mutations with probabilities for %s' % (ID, id_probabilities_file_path))
        prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList, inputDir, outputDir,
                                                                           jobname, ID, id_probabilities_file_path,
                                                                           startSimNum, endSimNum, ID,logger)

    # DBS
    if ((DBS in mutation_types_contexts) and (dbs_probabilities_file_path is not None)):
        logger.info('--- Merge %s mutations with probabilities for %s' % (DBS, dbs_probabilities_file_path))
        prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList, inputDir, outputDir,
                                                                           jobname, DBS, dbs_probabilities_file_path,
                                                                           startSimNum, endSimNum, DBS,logger)

    logger.info("--- Merge original chr based files with Mutation Probabilities: %s seconds" % (time.time() - start_time))
    logger.info("--- Merge original chr based files with Mutation Probabilities: %f minutes" % (float((time.time() - start_time) / 60)))
    logger.info('--- Merge original chr based files with Mutation Probabilities ends')
    logger.info('#################################################################################\n')
    ####################################################################################################################
    ##################  Merge original chr based files with Mutation Probabilities ends ################################
    ####################################################################################################################

    ###################################################################################################################
    ###################################################################################################################
    ###################################################################################################################
    existsSimulations= doesSimulationsAlreadyExits(outputDir,jobname,numofSimulations)

    if existsSimulations:
        logger.info('#################################################################################')
        logger.info('--- %d simulations already exists' %(numofSimulations))
        if new_simulations_enforced:
            logger.info('--- new_simulations_enforced:%s' %(new_simulations_enforced))
            logger.info('--- New simulations will be generated')
        else:
            logger.info('--- new_simulations_enforced:%s' %(new_simulations_enforced))
            logger.info('--- Existing simulations will be used')
        logger.info('#################################################################################\n')

    if ((numofSimulations>0) and ((new_simulations_enforced) or (not existsSimulations))):

        ###################################################################################################
        ############################  SigProfilerSimulator for n simulations starts #######################
        ###################################################################################################
        logger.info('#################################################################################')
        logger.info('--- SigProfilerSimulator for %d simulations starts' %(numofSimulations))
        start_time = time.time()
        #Call SigProfilerSimulator separately for each mutation type context otherwise it counts DBS mutations also in SBS mutations
        # Topography uses same mutation types with Simulator
        # '96' or '384' for single base substitutions (Simulator 1536, or 3072)
        # 'DBS' for double base substitutions
        # 'ID' for indels
        for mutation_type_context in mutation_types_contexts:
            mutation_type_context_for_simulator = []
            mutation_type_context_for_simulator.append(mutation_type_context)
            # Please notice that Simulator reverse the given input mutationTypes_for_simulator
            logger.info('--- SigProfilerSimulator is running for %s' %(mutation_type_context))
            simulator.SigProfilerSimulator(jobname, inputDir, genome, mutation_type_context_for_simulator,simulations=numofSimulations)

        logger.info("--- SigProfilerSimulator for %d simulations: %s seconds" %(numofSimulations,(time.time() -  start_time)))
        logger.info("--- SigProfilerSimulator for %d simulations: %f minutes" %(numofSimulations,float((time.time()-start_time)/60)))
        logger.info('--- SigProfilerSimulator for %d simulations ends' %(numofSimulations))
        logger.info('#################################################################################\n')
        ###################################################################################################
        ############################  SigProfilerSimulator for n simulations ends #########################
        ###################################################################################################

        ###################################################################################################
        ########################### Create simN directories for MatrixGenerator starts ####################
        ###################################################################################################
        logger.info('#################################################################################')
        logger.info('--- Create directories for %d simulations under inputDir/output/simulations/' %(numofSimulations))
        start_time = time.time()
        #Create directories sim1 to SimN under inputDir/output/simulations/
        access_rights = 0o755
        for simNum in range(1,numofSimulations+1):
            try:
                simName = 'sim%d' %(simNum)
                simDir = os.path.join(inputDir,'output','simulations',simName)
                if (not os.path.exists(simDir)):
                    os.mkdir(simDir, access_rights)
                for mutation_type_context in mutation_types_contexts:
                    simDir = os.path.join(inputDir,'output','simulations',simName,mutation_type_context)
                    if (not os.path.exists(simDir)):
                        os.mkdir(simDir, access_rights)
            except OSError:
                logger.info("Creation of the directory %s failed" %simDir)
            else:
                logger.info("Successfully created the directory %s" %simDir)

        for mutation_type_context in mutation_types_contexts:
            dirName = '%s_simulations_%s_%s' %(jobname, genome,mutation_type_context)
            copyFromDir = os.path.join(inputDir,'output','simulations',dirName)
            copyToMainDir= os.path.join(inputDir,'output','simulations')
            copyMafFiles(copyFromDir,copyToMainDir,mutation_type_context,numofSimulations)
        logger.info("--- Create directories and copy files: %s seconds ---" %(time.time()-start_time))
        logger.info("--- Create directories and copy files: %f minutes ---" %(float((time.time()-start_time)/60)))
        logger.info('#################################################################################\n')
        ###################################################################################################
        ########################### Create simN directories for MatrixGenerator ends ######################
        ###################################################################################################

        ###################################################################################################
        ####################### Run MatrixGenerator for each simulation starts ############################
        ###################################################################################################
        logger.info('#################################################################################')
        logger.info('--- Run SigProfilerMatrixGenerator for each simulation starts')
        start_time = time.time()
        for simNum in range(1,numofSimulations+1):
            simName = 'sim%d' %(simNum)
            #For each simulation we are calling matrix generator separately for each mutation type context

            logger.info('--- SigProfilerMatrixGenerator is run for %s starts' %(simName))
            for mutation_type_context in mutation_types_contexts:
                simInputDir=  os.path.join(inputDir,'output','simulations',simName,mutation_type_context)
                logger.info('For %s: %s simInputDir:%s' %(mutation_type_context,simName,simInputDir))
                matrices = matGen.SigProfilerMatrixGeneratorFunc(jobname,genome,simInputDir,plot=False, seqInfo=True)
                # print('matrices')
                # print(matrices)
                logger.info('#####################################')
            logger.info('--- SigProfilerMatrixGenerator is run for %s ends\n' % (simName))
        #sim1 matrix generator chrbased data will be under inputDir/output/simulations/sim1/96/output/vcf_files/SNV
        #sim1 matrix generator chrbased data will be under inputDir/output/simulations/sim1/ID/output/vcf_files/ID
        #sim1 matrix generator chrbased data will be under inputDir/output/simulations/sim1/DBS/output/vcf_files/DBS

        #simN matrix generator chrbased data will be under inputDir/output/simulations/simN/96/output/vcf_files/SNV
        #simN matrix generator chrbased data will be under inputDir/output/simulations/simN/ID/output/vcf_files/ID
        #simN matrix generator chrbased data will be under inputDir/output/simulations/simN/DBS/output/vcf_files/DBS
        logger.info("--- Run MatrixGenerator for each simulation: %s seconds" %(time.time()-start_time))
        logger.info("--- Run MatrixGenerator for each simulation: %f minutes" %(float((time.time()-start_time)/60)))
        logger.info('--- Run SigProfilerMatrixGenerator for each simulation ends')
        logger.info('#################################################################################\n')
        ###################################################################################################
        ####################### Run MatrixGenerator for each simulation ends ##############################
        ###################################################################################################

        ####################################################################################################################
        ##################  Merge simulations chr based files with Mutation Probabilities starts ###########################
        ####################################################################################################################
        logger.info('#################################################################################')
        logger.info('--- Merge simulations chr based files with Mutation Probabilities starts')
        logger.info('#################################################################################')
        startSimNum=1
        endSimNum=numofSimulations
        start_time = time.time()
        #SBS
        for mutation_type_context in mutation_types_contexts:
            if (mutation_type_context in SBS_CONTEXTS) and (sbs_probabilities_file_path is not None):
                logger.info('--- Merge %s mutations with probabilities for %s' %(mutation_type_context,sbs_probabilities_file_path))
                prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList,inputDir,outputDir,jobname,mutation_type_context,sbs_probabilities_file_path,startSimNum,endSimNum,'SNV',logger)

        #ID
        if ((ID in mutation_types_contexts) and (id_probabilities_file_path is not None)):
            logger.info('--- Merge %s mutations with probabilities for %s' % (ID, id_probabilities_file_path))
            prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList,inputDir,outputDir,jobname,'ID',id_probabilities_file_path,startSimNum,endSimNum,'ID',logger)

        #DBS
        if ((DBS in mutation_types_contexts) and (dbs_probabilities_file_path is not None)):
            logger.info('--- Merge %s mutations with probabilities for %s' % (DBS,dbs_probabilities_file_path))
            prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList,inputDir,outputDir,jobname,'DBS',dbs_probabilities_file_path,startSimNum,endSimNum,'DBS',logger)

        logger.info("--- Merge simulations chr based files with Mutation Probabilities: %s seconds" %(time.time()-start_time))
        logger.info("--- Merge simulations chr based files with Mutation Probabilities: %f minutes" %(float((time.time()-start_time)/60)))
        logger.info('--- Merge simulations chr based files with Mutation Probabilities ends')
        logger.info('#################################################################################\n')
        ####################################################################################################################
        ##################  Merge simulations chr based files with Mutation Probabilities ends #############################
        ####################################################################################################################

    ###################################################################################################################
    ###################################################################################################################
    ###################################################################################################################

    #################################################################################
    #Please note that if you have chr based subs, indels, dinucs mutations combined with probabilities under data directory
    #You can comment the above code and run the rest of the code below.
    #################################################################################


    #################################################################################
    logger.info('#################################################################################')
    logger.info('--- Fill dictioaries using original data starts')
    # #################################################################################
    #For each signature we will find a cutoff value for mutations with average probability >=0.9
    #Our aim is to have at most 10% false positive rate in mutations
    #number of mutations >= 5K for subs signatures
    #number of mutations >= 1K for indels signatures
    #number of mutations >= 200 for dinuc signatures
    #If we can not satisfy this condition we will discard the signature

    cutoffs=[]
    for cufoff in np.arange(0.5,0.91,0.01):
        cutoffs.append("%.2f" %(cufoff))

    subsSignature2PropertiesListDict=fillCutoff2Signature2PropertiesListDictionary(outputDir,jobname,chromNamesList,SUBS,cutoffs)
    indelsSignature2PropertiesListDict=fillCutoff2Signature2PropertiesListDictionary(outputDir,jobname,chromNamesList,INDELS,cutoffs)
    dinucsSignature2PropertiesListDict=fillCutoff2Signature2PropertiesListDictionary(outputDir,jobname,chromNamesList,DINUCS,cutoffs)
    ##################################################################################

    #Create files
    createFiles(outputDir, jobname, MutationType2NumberofMutatiosDictFilename)

    #Initialize
    mutationType2NumberofMutationsDict = {}

    #Using original data
    fill_mutations_dictionaries_write(outputDir,jobname,chromNamesList,SUBS,mutationType2NumberofMutationsDict,subsSignature2PropertiesListDict)
    fill_mutations_dictionaries_write(outputDir,jobname,chromNamesList,INDELS,mutationType2NumberofMutationsDict,indelsSignature2PropertiesListDict)
    fill_mutations_dictionaries_write(outputDir,jobname,chromNamesList,DINUCS,mutationType2NumberofMutationsDict,dinucsSignature2PropertiesListDict)

    #We are writing number of mutations for each mutation type.
    # e.g.: {"SUBS": 3982196, "INDELS": 234731}
    appendDictionaryUnderDataDirectory(mutationType2NumberofMutationsDict,outputDir,jobname,MutationType2NumberofMutatiosDictFilename)
    logger.info('--- Fill dictioaries using original data ends')
    logger.info('#################################################################################\n')
    #################################################################################

    #################################################################################
    ################## Set full path library files starts ###########################
    #We need full path of the library files
    if ((len(epigenomics_files)==1) and (epigenomics_files[0]==DEFAULT_HISTONE_OCCUPANCY_FILE)):
        epigenomics_files[0] = os.path.join(current_abs_path,LIB,EPIGENOMICS,DEFAULT_HISTONE_OCCUPANCY_FILE)

    if (nucleosome_file== DEFAULT_NUCLEOSOME_OCCUPANCY_FILE):
        nucleosome_file = os.path.join(current_abs_path,LIB,NUCLEOSOME,DEFAULT_NUCLEOSOME_OCCUPANCY_FILE)

    if (replication_time_file == DEFAULT_REPLICATION_TIME_SIGNAL_FILE):
        replication_time_file = os.path.join(current_abs_path,LIB,REPLICATION,DEFAULT_REPLICATION_TIME_SIGNAL_FILE)

    if (replication_time_valley_file == DEFAULT_REPLICATION_TIME_VALLEY_FILE):
        replication_time_valley_file = os.path.join(current_abs_path,LIB,REPLICATION,DEFAULT_REPLICATION_TIME_VALLEY_FILE)

    if (replication_time_peak_file == DEFAULT_REPLICATION_TIME_PEAK_FILE):
        replication_time_peak_file = os.path.join(current_abs_path,LIB,REPLICATION,DEFAULT_REPLICATION_TIME_PEAK_FILE)
    ################## Set full path library files ends #############################
    #################################################################################

    ####################################################################################################################
    ################################### Run SigProfilerTopography Analysis starts ######################################
    ####################################################################################################################
    logger.info('#################################################################################')
    logger.info('--- Run SigProfilerTopography Analysis starts')
    plusorMinus_epigenomics = 2000
    plusorMinus_nucleosome = 1000

    if (epigenomics):
        #Epigenomics
        occupancy_type=EPIGENOMICSOCCUPANCY
        deleteOldData(outputDir,jobname,occupancy_type)

        for idx, epigenomics_file in enumerate(epigenomics_files):
            start_time = time.time()
            if idx<len(epigenomics_files_memos):
                epigenomics_file_memo= epigenomics_files_memos[idx]
            else:
                epigenomics_file_memo= None
            runOccupancyAnalyses(outputDir,jobname,numofSimulations,sample_based,epigenomics_file,epigenomics_file_memo,chromSizesDict,chromNamesList,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,computation_type,occupancy_type,plusorMinus_epigenomics)
            logger.info('#################################################################################')
            logger.info("--- Run Epigenomics Analyses: %s seconds --- %s" %((time.time()-start_time),computation_type))
            logger.info("--- Run Epigenomics Analyses: %f minutes --- %s" %(float((time.time()-start_time)/60),computation_type))
            logger.info('#################################################################################\n')

    if (nucleosome):
        #Nucleosome Occupancy
        occupancy_type = NUCLEOSOMEOCCUPANCY
        deleteOldData(outputDir,jobname,occupancy_type)

        start_time = time.time()
        runOccupancyAnalyses(outputDir,jobname,numofSimulations,sample_based,nucleosome_file,None,chromSizesDict,chromNamesList,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,computation_type,occupancy_type,plusorMinus_nucleosome)
        logger.info('#################################################################################')
        logger.info("--- Run Nucleosome Occupancy Analyses: %s seconds --- %s" %((time.time()-start_time),computation_type))
        logger.info("--- Run Nucleosome Occupancy Analyses: %f minutes --- %s" %(float((time.time()-start_time)/60),computation_type))
        logger.info('#################################################################################\n')

    if (replication_time):
        # Replication Time
        start_time = time.time()
        runReplicationTimeAnalysis(genome,outputDir,jobname,numofSimulations,sample_based,replication_time_file,chromSizesDict,chromNamesList,computation_type,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict)
        logger.info('#################################################################################')
        logger.info("--- Run Replication Time Analyses: %s seconds --- %s" %((time.time()-start_time),computation_type))
        logger.info("--- Run Replication Time Analyses: %f minutes --- %s" %(float((time.time()-start_time)/60),computation_type))
        logger.info('#################################################################################\n')

    if (strand_bias):
        # Replication Strand Bias
        start_time = time.time()
        runReplicationStrandBiasAnalysis(outputDir,jobname,numofSimulations,sample_based,replication_time_file,replication_time_valley_file,replication_time_peak_file,chromSizesDict,chromNamesList,computation_type,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict)
        logger.info('#################################################################################')
        logger.info("--- Run Replication Strand Bias Analyses: %s seconds --- %s" %((time.time()-start_time),computation_type))
        logger.info("--- Run Replication Strand Bias Analyses: %f minutes --- %s" %(float((time.time()-start_time)/60),computation_type))
        logger.info('#################################################################################\n')

        # Transcription Strand Bias
        start_time = time.time()
        runTranscriptionStradBiasAnalysis(genome,outputDir,jobname,numofSimulations,sample_based,chromSizesDict,chromNamesList,computation_type,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict)
        logger.info('#################################################################################')
        logger.info("--- Run Transcription Strand Bias Analyses: %s seconds --- %s" %((time.time()-start_time),computation_type))
        logger.info("--- Run Transcription Strand Bias Analyses: %f minutes --- %s" %(float((time.time()-start_time)/60),computation_type))
        logger.info('#################################################################################\n')

    if (processivity):
        # Processivity
        start_time = time.time()
        #TODO shall we condider only the signatures in subsSignature2PropertiesListDict?
        runProcessivityAnalysis(mutation_types_contexts,outputDir,jobname,numofSimulations,chromNamesList,subsSignature2PropertiesListDict)
        logger.info('#################################################################################')
        logger.info("--- Run Processivity Analyses: %s seconds ---" %(time.time()-start_time))
        logger.info("--- Run Processivity Analyses: %f minutes ---" %(float((time.time()-start_time)/60)))
        logger.info('#################################################################################\n')

    logger.info('#################################################################################\n')
    ####################################################################################################################
    ################################### Run SigProfilerTopography Analysis ends ########################################
    ####################################################################################################################

    ####################################################################################################################
    ############################################ Plot figures starts ###################################################
    ####################################################################################################################
    if (plot_figures):
        logger.info('#################################################################################')
        logger.info('--- Plot figures starts')
        start_time = time.time()
        # plotFigures(outputDir, jobname, numofSimulations, 'BONFERRONI_CORRECTION', 'USING_ONE_SAMPLE_TTEST')
        # plotFigures(outputDir, jobname, numofSimulations, 'BONFERRONI_CORRECTION', 'USING_NULL_DISTRIBUTION')
        plotFigures(outputDir, jobname, numofSimulations, sample_based,'USING_ZSCORE', 'USING_ZSCORE',mutation_types_contexts,epigenomics_files,epigenomics_files_memos,nucleosome_file,epigenomics,nucleosome,replication_time,strand_bias,processivity,plusorMinus_epigenomics,plusorMinus_nucleosome)
        logger.info('#################################################################################')
        logger.info("--- Plot Figures: %s seconds ---" %(time.time()-start_time))
        logger.info("--- Plot Figures: %f minutes ---" %(float((time.time()-start_time)/60)))
        logger.info('--- Plot figures ends')
        logger.info('#################################################################################\n')
    ####################################################################################################################
    ############################################ Plot figures ends #####################################################
    ####################################################################################################################

    logger.info('#################################################################################')
    logger.info("--- SigProfilerTopography ended successfully")
    logger.info("--- Thanks for using SigProfilerTopography")
    logger.info('#################################################################################\n')


#######################################################



##############################################################
# BONFERRONI_CORRECTION = 'BONFERRONI_CORRECTION'
# FDR_BH_CORRECTION = 'FDR_BH_CORRECTION'
# USING_POISSON_DISTRIBUTION = 'USING_POISSON_DISTRIBUTION'
# USING_NULL_DISTRIBUTION = 'USING_NULL_DISTRIBUTION'
# USING_GAUSSIAN_KDE = 'USING_GAUSSIAN_KDE'
#Plot Figures for the attainded data after SigProfilerTopography Analyses
def plotFigures(outputDir,jobname,numberofSimulations,sample_based,multipleTesting,probabilityCalculation,mutationTypes,epigenomics_files,epigenomics_files_memos,nucleosome_file,epigenomics,nucleosome,replication_time,strand_bias,processivity,plusOrMinus_epigenomics,plusOrMinus_nucleosome):

    #Internally Set
    figureAugmentation = 'noaugmentation'

    jobnameSamplesPath = os.path.join(outputDir,jobname,FIGURE,SAMPLES)
    print('Topography.py jobnameSamplesPath:%s ' %jobnameSamplesPath)

    ############################################################
    if (os.path.exists(jobnameSamplesPath)):
        try:
            shutil.rmtree(jobnameSamplesPath)
        except OSError as e:
            print('Error: %s - %s.' % (e.filename, e.strerror))
    ############################################################

    ############################################################
    if (epigenomics):
        occupancy_type=EPIGENOMICSOCCUPANCY
        deleteOldFigures(outputDir, jobname, occupancy_type)

        for idx, epigenomics_file in enumerate(epigenomics_files):
            epigenomics_file_basename = os.path.basename(epigenomics_file)

            if idx<len(epigenomics_files_memos):
                epigenomics_file_memo= epigenomics_files_memos[idx]
            else:
                epigenomics_file_memo= None
            occupancyAverageSignalFigures(outputDir, jobname, figureAugmentation, numberofSimulations,sample_based, mutationTypes,epigenomics_file_basename,epigenomics_file_memo,occupancy_type,plusOrMinus_epigenomics)
    if (nucleosome):
        occupancy_type=NUCLEOSOMEOCCUPANCY
        deleteOldFigures(outputDir, jobname, occupancy_type)
        nucleosome_file_basename = os.path.basename(nucleosome_file)
        occupancyAverageSignalFigures(outputDir,jobname,figureAugmentation,numberofSimulations,sample_based,mutationTypes,nucleosome_file_basename,None,occupancy_type,plusOrMinus_nucleosome)
    if (replication_time):
        replicationTimeNormalizedMutationDensityFigures(outputDir,jobname,figureAugmentation,numberofSimulations,sample_based,mutationTypes)
    if (strand_bias):
        transcriptionReplicationStrandBiasFigures(outputDir,jobname,figureAugmentation,numberofSimulations,sample_based)
    if (processivity):
        processivityFigures(outputDir,jobname,numberofSimulations,multipleTesting,probabilityCalculation)
    ############################################################

##############################################################



# # ##############################################################
# import os
#
# if __name__== "__main__":
#     genome= 'GRCh37'
#     inputDir ='/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/input_for_matgen/BreastCancer560_subs_indels_dinucs'
#     outputDir = os.path.join('C:\\','Users','burcak','Developer','Python','SigProfilerTopography','SigProfilerTopography','output_test')
#     jobname = 'BreastCancer560'
#     numberofSimulations = 2
#     subs_probabilities = os.path.join('C:\\','Users','burcak','Documents','DrLudmilAlexandrovLab','SigProfilerTopography','SigProfilerTopographyInput','Extractor','SBS_Mutation_Probabilities.txt')
#     indels_probabilities_file_path = os.path.join('C:\\','Users','burcak','Documents','DrLudmilAlexandrovLab','SigProfilerTopography','SigProfilerTopographyInput','Extractor','ID_Mutation_Probabilities.txt')
#     dinucs_probabilities_file_path = os.path.join('C:\\','Users','burcak','Documents','DrLudmilAlexandrovLab','SigProfilerTopography','SigProfilerTopographyInput','Extractor','DBS_Mutation_Probabilities.txt')
#     runAnalyses(genome,inputDir,outputDir,jobname,numberofSimulations,subs_probabilities_file_path=subs_probabilities,dinucs_probabilities_file_path=dinucs_probabilities_file_path,mutationTypes=['SUBS','DINUCS'])
# # ##############################################################