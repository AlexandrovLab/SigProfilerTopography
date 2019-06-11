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
import SigProfilerMatrixGenerator as matgen_package
# print(matgen_package.__path__[0])
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerSimulator import SigProfilerSimulator as simulator

from SigProfilerTopography.source.commons.DataPreparationCommons import readMutationsWithGenomicPositions
from SigProfilerTopography.source.commons.DataPreparationCommons import readProbabilities
from SigProfilerTopography.source.commons.DataPreparationCommons import readChrBasedMutationsMergeWithProbabilitiesAndWrite

from SigProfilerTopography.source.commons.TopographyCommons import *

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

import subprocess
import shutil


############################################################
#CAn be move to DataPreparationCommons under /source/commons
#read chr based dinucs (provided by SigProfilerMatrixGenerator) and merge with probabilities (provided by SigProfilerTopography)
def prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList,inputDir,outputDir,jobname,mutation_type_context,mutations_probabilities_file_path,numofSimulations,partialDirname):

    #original matrix generator chrbased data will be under inputDir/output/vcf_files/SNV
    #original matrix generator chrbased data will be under inputDir/output/vcf_files/DBS
    #original matrix generator chrbased data will be under inputDir/output/vcf_files/ID

    #new
    #sim1 matrix generator chrbased data will be under inputDir/output/simulations/sim1/96/output/vcf_files/SNV
    #sim1 matrix generator chrbased data will be under inputDir/output/simulations/sim1/ID/output/vcf_files/ID
    #sim1 matrix generator chrbased data will be under inputDir/output/simulations/sim1/DBS/output/vcf_files/DBS

    os.makedirs(os.path.join(outputDir,jobname,DATA,CHRBASED),exist_ok=True)
    for simNum in range(1,numofSimulations+1):
        simName = 'sim%d' % (simNum)
        os.makedirs(os.path.join(outputDir, jobname, DATA,CHRBASED,simName), exist_ok=True)

    if (os.path.exists(mutations_probabilities_file_path)):
        print('For mutation_type_context:%s mutations_probabilities_file_path:%s' %(mutation_type_context,mutations_probabilities_file_path))
        mutations_probabilities_df = readProbabilities(mutations_probabilities_file_path)
        # dinucs_probabilities_df.columns.names [Sample Names    MutationTypes   DBS2    DBS4    DBS6    DBS7    DBS11]
        # dinucs_probabilities_df.columns.names [Sample    Mutation   DBS2    DBS4    DBS6    DBS7    DBS11]
        mutations_probabilities_df.rename(columns={'Sample Names': 'Sample', 'MutationTypes': 'Mutation'}, inplace=True)

        numofProcesses = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(numofProcesses)

        poolInputList = []

        for simNum in range(0,numofSimulations+1):
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
def download_nucleosome_occupancy_convert_bigWig2wig(cellLine):
    bigWig2Wig_filepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME, BIGWIG2WIG)
    os.chmod(bigWig2Wig_filepath,0o744)
    if (cellLine==GM12878):
        gm12878_bigWig_filepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME, ENCODE_NUCLEOSOME_GM12878_BIGWIG)
        gm12878_wig_filepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME,ENCODE_NUCLEOSOME_GM12878_WIG)
        downloadFromWeb(ENCODE_NUCLEOSOME_GM12878_BIGWIG_URL, gm12878_bigWig_filepath)
        subprocess.call([bigWig2Wig_filepath, gm12878_bigWig_filepath,gm12878_wig_filepath])
    elif (cellLine==K562):
        K562_bigWig_filepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME, ENCODE_NUCLEOSOME_K562_BIGWIG)
        K562_wig_filepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME, ENCODE_NUCLEOSOME_K562_WIG)
        downloadFromWeb(ENCODE_NUCLEOSOME_K562_BIGWIG_URL, K562_bigWig_filepath)
        subprocess.call([bigWig2Wig_filepath, K562_bigWig_filepath,K562_wig_filepath])
#######################################################

#######################################################
def runNucleosomeOccupancyAnalyses(genome,outputDir,jobname,numofSimulations,nucleosomeFilename,chromSizesDict,chromNamesList,availableLibraryFilenamesList,computation_type):
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

        #Check whether nucleosomeFilename_wo_dir is downloaded
        if (os.path.exists(nucleosomeFilename)):
            quantileValue = round(float(0.97), 2)
            # PartitionNucleosomeOccupancyData.partitionNucleosomeOccupancyData(jobname,nucleosomeFilename,quantileValue)
            readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArrays(genome,quantileValue,nucleosomeFilename)
            #append
            append2File(nucleosomeFilename_wo_dir,AVAILABLE_LIBRARY_FILENAMES_PATH)
        else:
            download_nucleosome_command = 'download_nucleosome_occupancy_convert_bigWig2wig(cellLine) command'
            print('You need to download %s using %s' %(nucleosomeFilename_wo_dir,download_nucleosome_command))

    # computationType = COMPUTATION_CHROMOSOMES_SEQUENTIAL_ALL_SIMULATIONS_PARALLEL
    # computationType = COMPUTATION_CHROMOSOMES_SEQUENTIAL_SIMULATIONS_SEQUENTIAL
    # computationType = COMPUTATION_ALL_CHROMOSOMES_PARALLEL
    print('Nucleosome Occupancy Computation Type:%s' %(computation_type))
    nucleosomeOccupancyAnalysis(computation_type,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,nucleosomeFilename)
    ###############################################
#######################################################


#######################################################
def runReplicationTimeAnalysis(genome,outputDir,jobname,numofSimulations,replicationTimeFilename,chromSizesDict,chromNamesList,computation_type):
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
    print('Replication Time Analyis Computation Type:%s' %(computation_type))
    replicationTimeAnalysis(computation_type,replication_time_np_arrays_fill_runtime,genome,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,replicationTimeFilename)
    ###############################################

#######################################################


#######################################################
def runReplicationStrandBiasAnalysis(outputDir,jobname,numofSimulations,replicationTimeFilename,replicationTimeValleyFilename,replicationTimePeakFilename,chromSizesDict,chromNamesList,computation_type):

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
    replicationStrandBiasAnalysis(computation_type,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,smoothedWaveletRepliseqDataFilename,valleysBEDFilename,peaksBEDFilename)
    ###############################################

#######################################################

#######################################################
def runTranscriptionStradBiasAnalysis(genome,outputDir,jobname,numofSimulations,chromSizesDict,chromNamesList,computation_type):
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
    transcriptionStrandBiasAnalysis(computation_type,useTranscriptionStrandColumn,genome,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations)
    ###############################################
#######################################################


#######################################################
def runProcessivityAnalysis(mutation_types_contexts,outputDir,jobname,numofSimulations,chromNamesList):
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

    processivityAnalysis(mutation_types_contexts,chromNamesList,outputDir,jobname,numofSimulations,considerProbabilityInProcessivityAnalysis)
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
def runAnalyses(genome,inputDir, outputDir,jobname,numofSimulations,subs_probabilities_file_path=None,indels_probabilities_file_path= None,dinucs_probabilities_file_path=None,nucleosomeFilename=DEFAULT_NUCLEOSOME_OCCUPANCY_FILE,replicationTimeFilename=DEFAULT_REPLICATION_TIME_SIGNAL_FILE,replicationTimeValleyFilename=DEFAULT_REPLICATION_TIME_VALLEY_FILE,replicationTimePeakFilename=DEFAULT_REPLICATION_TIME_PEAK_FILE,mutation_types_contexts=[SBS96,ID,DBS],computation_type=COMPUTATION_CHROMOSOMES_SEQUENTIAL_ALL_SIMULATIONS_PARALLEL):

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
    print('mutation_types_contexts:%s' %mutation_types_contexts)
    print('computation_type:%s' %computation_type)
    print('#################################################################################')

    #################################################################################
    print('#################################################################################')
    print('Chromosome names: %s' %(chromNamesList))
    print('Chromosome short names: %s' % (chromShortNamesList))
    current_abs_path = os.path.dirname(os.path.realpath(__file__))
    print('current_abs_path: %s ' % current_abs_path)
    print('#################################################################################')
    #################################################################################

    ###################################################################################################
    #######################  SigProfilerMatrixGenerator for original data starts ######################
    ###################################################################################################
    #Run MatrixGenerator for original data: this call prepares chrBased input files for original data
    print('################## SigProfilerMatrixGenerator for original data #################')
    start_time = time.time()

    print('For original data inputDir:%s' % (inputDir))
    matrices = matGen.SigProfilerMatrixGeneratorFunc(jobname,genome,inputDir,plot=False, seqInfo=True)
    # print('matrices')
    # print(matrices)

    # original matrix generator chrbased data will be under inputDir/output/vcf_files/SNV
    # original matrix generator chrbased data will be under inputDir/output/vcf_files/DBS
    # original matrix generator chrbased data will be under inputDir/output/vcf_files/ID

    print("--- SigProfilerMatrixGenerator for original data: %s seconds ---" % (time.time() -  start_time))
    print("--- SigProfilerMatrixGenerator for original data: %f minutess ---" %float((time.time() - start_time)/60))
    print('#################################################################################')
    ###################################################################################################
    #######################  SigProfilerMatrixGenerator for original data ends ########################
    ###################################################################################################

    ###################################################################################################
    ############################  SigProfilerSimulator for n simulations starts #######################
    ###################################################################################################
    print('################## SigProfilerSimulator starts for %d simulations ###############' %(numofSimulations))
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
        print('SigProfilerSimulator is running for %s' %(mutation_type_context))
        simulator.SigProfilerSimulator(jobname, inputDir, genome, mutation_type_context_for_simulator,simulations=numofSimulations)
    print("--- SigProfilerSimulator for %d simulations: %s seconds ---" %(numofSimulations,(time.time() -  start_time)))
    print("--- SigProfilerSimulator for %d simulations: %f minutes ---" %(numofSimulations,float((time.time()-start_time)/60)))
    print('#################################################################################')
    ###################################################################################################
    ############################  SigProfilerSimulator for n simulations ends #########################
    ###################################################################################################

    ###################################################################################################
    ########################### Create simN directories for MatrixGenerator starts ####################
    ###################################################################################################
    print('################## Create directories for sim1 to sim%d under inputDir/output/simulations/ #################' %(numofSimulations))
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
            print("Creation of the directory %s failed" %simDir)
        else:
            print("Successfully created the directory %s" %simDir)
    print('#################################################################################')

    print('################### Copy sample files under corresponding simulation directory ######################')
    #Copy sample files under corresponding simulation directory
    # simulations sample based directories will be under inputDir/jobname/output/simulations/jobname_simulations_genome_mutationTypes_for_simulator[lastIndex]_..._mutationTypes_for_simulator[firstIndex]
    for mutation_type_context in mutation_types_contexts:
        dirName = '%s_simulations_%s' %(jobname, genome)
        dirName = dirName + '_' + mutation_type_context
        copyFromDir = os.path.join(inputDir,'output','simulations',dirName)
        copySampleFilesToCorrespondingSimulationDirectory(inputDir,copyFromDir,mutation_type_context)
    print("--- Create directories and copy files: %s seconds ---" %(time.time()-start_time))
    print("--- Create directories and copy files: %f minutes ---" %(float((time.time()-start_time)/60)))
    print('#################################################################################')
    ###################################################################################################
    ########################### Create simN directories for MatrixGenerator ends ######################
    ###################################################################################################

    ###################################################################################################
    ####################### Run MatrixGenerator for each simulation starts ############################
    ###################################################################################################
    print('################## Run SigProfilerMatrixGenerator for each simulation starts ####################')
    start_time = time.time()
    for simNum in range(1,numofSimulations+1):
        simName = 'sim%d' %(simNum)
        #For each simulation we are calling matrix generator separately for each mutation type context
        print('SigProfilerMatrixGenerator is run for %s starts' %(simName))
        for mutation_type_context in mutation_types_contexts:
            simInputDir=  os.path.join(inputDir,'output','simulations',simName,mutation_type_context)
            print('For %s: %s simInputDir:%s' %(mutation_type_context,simName,simInputDir))
            matrices = matGen.SigProfilerMatrixGeneratorFunc(jobname,genome,simInputDir,plot=False, seqInfo=True)
            # print('matrices')
            # print(matrices)
            print('#####################################')
        print('SigProfilerMatrixGenerator is run for %s ends' % (simName))
    #sim1 matrix generator chrbased data will be under inputDir/output/simulations/sim1/96/output/vcf_files/SNV
    #sim1 matrix generator chrbased data will be under inputDir/output/simulations/sim1/ID/output/vcf_files/ID
    #sim1 matrix generator chrbased data will be under inputDir/output/simulations/sim1/DBS/output/vcf_files/DBS

    #simN matrix generator chrbased data will be under inputDir/output/simulations/simN/96/output/vcf_files/SNV
    #simN matrix generator chrbased data will be under inputDir/output/simulations/simN/ID/output/vcf_files/ID
    #simN matrix generator chrbased data will be under inputDir/output/simulations/simN/DBS/output/vcf_files/DBS
    print("--- Run MatrixGenerator for each simulation: %s seconds ---" %(time.time()-start_time))
    print("--- Run MatrixGenerator for each simulation: %f minutes ---" %(float((time.time()-start_time)/60)))
    print('#######################################################################')
    ###################################################################################################
    ####################### Run MatrixGenerator for each simulation ends ##############################
    ###################################################################################################

    ####################################################################################################################
    ##################  Merge original and simulations chr based files with Mutation Probabilities starts ##############
    ####################################################################################################################
    print('##################  Merge original and simulations chr based files with Mutation Probabilities starts ####################')
    start_time = time.time()
    #SBS
    for mutation_type_context in mutation_types_contexts:
        if (mutation_type_context in SBS_CONTEXTS) and (subs_probabilities_file_path is not None):
            print('########## Merge mutations with probabilities for %s ##########' %(subs_probabilities_file_path))
            prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList,inputDir,outputDir,jobname,mutation_type_context,subs_probabilities_file_path,numofSimulations,'SNV')
            print('###############################################################')

    #ID
    if ((ID in mutation_types_contexts) and (indels_probabilities_file_path is not None)):
        print('Merge mutations with probabilities for %s' % (indels_probabilities_file_path))
        prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList,inputDir,outputDir,jobname,'ID',indels_probabilities_file_path,numofSimulations,'ID')
        print('###############################################################')

    #DBS
    if ((DBS in mutation_types_contexts) and (dinucs_probabilities_file_path is not None)):
        print('Merge mutations with probabilities for %s' % (dinucs_probabilities_file_path))
        prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList,inputDir,outputDir,jobname,'DBS',dinucs_probabilities_file_path,numofSimulations,'DBS')
        print('###############################################################')
    print("--- Merge original and simulations chr based files with Mutation Probabilities: %s seconds ---" %(time.time()-start_time))
    print("--- Merge original and simulations chr based files with Mutation Probabilities: %f minutes ---" %(float((time.time()-start_time)/60)))
    print('#######################################################################')
    ####################################################################################################################
    ##################  Merge original and simulations chr based files with Mutation Probabilities ends ################
    ####################################################################################################################

    #################################################################################
    print('####### Fill dictioaries using original data starts ###################')
    #Using original data
    fill_mutations_dictionaries_write(outputDir, jobname, chromNamesList,SUBS)
    fill_mutations_dictionaries_write(outputDir, jobname, chromNamesList,INDELS)
    fill_mutations_dictionaries_write(outputDir, jobname, chromNamesList,DINUCS)
    print('#######################################################################')
    #################################################################################

    #################################################################################
    print('############## Set full path library files starts #####################')
    #We need full path of the library files
    if (nucleosomeFilename == DEFAULT_NUCLEOSOME_OCCUPANCY_FILE):
        nucleosomeFilename = os.path.join(current_abs_path,LIB,NUCLEOSOME,DEFAULT_NUCLEOSOME_OCCUPANCY_FILE)

    if (replicationTimeFilename == DEFAULT_REPLICATION_TIME_SIGNAL_FILE):
        replicationTimeFilename = os.path.join(current_abs_path,LIB,REPLICATION,DEFAULT_REPLICATION_TIME_SIGNAL_FILE)

    if (replicationTimeValleyFilename == DEFAULT_REPLICATION_TIME_VALLEY_FILE):
        replicationTimeValleyFilename = os.path.join(current_abs_path,LIB,REPLICATION,DEFAULT_REPLICATION_TIME_VALLEY_FILE)

    if (replicationTimePeakFilename == DEFAULT_REPLICATION_TIME_PEAK_FILE):
        replicationTimePeakFilename = os.path.join(current_abs_path,LIB,REPLICATION,DEFAULT_REPLICATION_TIME_PEAK_FILE)
    print('#######################################################################')
    #################################################################################

    ####################################################################################################################
    ################################### Run SigProfilerTopography Analysis starts ######################################
    ####################################################################################################################
    print('################# Run SigProfilerTopography Analysis starts ###########')
    start_time = time.time()
    runNucleosomeOccupancyAnalyses(genome,outputDir,jobname,numofSimulations,nucleosomeFilename,chromSizesDict,chromNamesList,availableLibraryFilenamesList,computation_type)
    print("--- Run Nucleosome Occupancy Analyses: %s seconds --- %s" %((time.time()-start_time),computation_type))
    print("--- Run Nucleosome Occupancy Analyses: %f minutes --- %s" %(float((time.time()-start_time)/60),computation_type))

    start_time = time.time()
    runReplicationTimeAnalysis(genome,outputDir,jobname,numofSimulations,replicationTimeFilename,chromSizesDict,chromNamesList,computation_type)
    print("--- Run Replication Time Analyses: %s seconds --- %s" %((time.time()-start_time),computation_type))
    print("--- Run Replication Time Analyses: %f minutes --- %s" %(float((time.time()-start_time)/60),computation_type))

    start_time = time.time()
    runReplicationStrandBiasAnalysis(outputDir,jobname,numofSimulations,replicationTimeFilename,replicationTimeValleyFilename,replicationTimePeakFilename,chromSizesDict,chromNamesList,computation_type)
    print("--- Run Replication Strand Bias Analyses: %s seconds --- %s" %((time.time()-start_time),computation_type))
    print("--- Run Replication Strand Bias Analyses: %f minutes --- %s" %(float((time.time()-start_time)/60),computation_type))

    start_time = time.time()
    runTranscriptionStradBiasAnalysis(genome,outputDir,jobname,numofSimulations,chromSizesDict,chromNamesList,computation_type)
    print("--- Run Transcription Strand Bias Analyses: %s seconds --- %s" %((time.time()-start_time),computation_type))
    print("--- Run Transcription Strand Bias Analyses: %f minutes --- %s" %(float((time.time()-start_time)/60),computation_type))

    start_time = time.time()
    runProcessivityAnalysis(mutation_types_contexts,outputDir,jobname,numofSimulations,chromNamesList)
    print("--- Run Processivity Analyses: %s seconds ---" %(time.time()-start_time))
    print("--- Run Processivity Analyses: %f minutes ---" %(float((time.time()-start_time)/60)))
    print('#######################################################################')
    ####################################################################################################################
    ################################### Run SigProfilerTopography Analysis ends ########################################
    ####################################################################################################################

    ####################################################################################################################
    ############################################ Plot figures starts ###################################################
    ####################################################################################################################
    print('################ Plot figures starts ##################################')
    start_time = time.time()
    # plotFigures(outputDir, jobname, numofSimulations, 'BONFERRONI_CORRECTION', 'USING_ONE_SAMPLE_TTEST')
    # plotFigures(outputDir, jobname, numofSimulations, 'BONFERRONI_CORRECTION', 'USING_NULL_DISTRIBUTION')
    plotFigures(outputDir, jobname, numofSimulations, 'USING_ZSCORE', 'USING_ZSCORE')
    print("--- Plot Figures: %s seconds ---" %(time.time()-start_time))
    print("--- Plot Figures: %f minutes ---" %(float((time.time()-start_time)/60)))
    print('#######################################################################')
    ####################################################################################################################
    ############################################ Plot figures ends #####################################################
    ####################################################################################################################

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