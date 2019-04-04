# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018 Burcak Otlu

import shutil

# #############################################################
# import sys
# import os
#
# current_abs_path = os.path.dirname(os.path.realpath(__file__))
# print('Topography.py current_abs_path:%s' %(current_abs_path))
#
# commonsPath = os.path.join(current_abs_path, '..','commons')
# commonsPath = os.path.join(current_abs_path,'commons')
# sys.path.append(commonsPath)
# #############################################################

from SigProfilerTopography.source.commons.DataPreparationCommons import readMutationsWithGenomicPositions
from SigProfilerTopography.source.commons.DataPreparationCommons import readProbabilities
from SigProfilerTopography.source.commons.DataPreparationCommons import mergeSNPsWithSignatureProbabilities
from SigProfilerTopography.source.commons.DataPreparationCommons import readIndelsandWriteWithGenomicPositions
from SigProfilerTopography.source.commons.DataPreparationCommons import prepareSimulationBasedInputFilesForSigProfilerTopography

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


# #######################################################
# def downloadNuclesomeChrBasedSignalArrays(nucleosomeFilename):
#     left here
#     use wget
#     upload to a ftp site
#     download from that ftp site
#     https://stackoverflow.com/questions/25010369/wget-curl-large-file-from-google-drive/39225039#39225039
# #######################################################


#######################################################
#Run SigProfilerTopography Analyses
def runAnalyses(genome, singlePointMutationsFilename,indelsFilename,outputDir,jobname,nucleosomeFilename,replicationTimeFilename,replicationTimeValleyFilename,replicationTimePeakFilename):

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
    ##############################################


    ##############################################
    #Please note that singlePointMutationsFilename and indelsFilename are full paths.
    # They are all read and partitioned chr based
    # Get the fienames at the end
    singlePointMutationsFilename = os.path.basename(singlePointMutationsFilename)
    indelsFilename = os.path.basename(indelsFilename)
    ##############################################

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
    nucleosomeOccupancyAnalysis(computationType,chromSizesDict,chromNamesList,outputDir,jobname,singlePointMutationsFilename,indelsFilename,nucleosomeFilename_woDir)
    ###############################################


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
    replicationTimeAnalysis(computationType,replication_time_np_arrays_fill_runtime,genome,chromSizesDict,chromNamesList,outputDir,jobname,singlePointMutationsFilename,indelsFilename,replicationTimeFilename)
    ###############################################


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
        replicationStrandBiasAnalysis(computationType,chromSizesDict,chromNamesList,outputDir,jobname,singlePointMutationsFilename,indelsFilename,smoothedWaveletRepliseqDataFilename,valleysBEDFilename,peaksBEDFilename)
    ###############################################


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
        transcriptionStrandBiasAnalysis(computationType,genome,chromSizesDict,chromNamesList,outputDir,jobname, singlePointMutationsFilename,indelsFilename)
    ###############################################


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



##############################################################
# BONFERRONI_CORRECTION = 'BONFERRONI_CORRECTION'
# FDR_BH_CORRECTION = 'FDR_BH_CORRECTION'
#
# USING_POISSON_DISTRIBUTION = 'USING_POISSON_DISTRIBUTION'
# USING_NULL_DISTRIBUTION = 'USING_NULL_DISTRIBUTION'
# USING_GAUSSIAN_KDE = 'USING_GAUSSIAN_KDE'

#Plot Figures for the attainded data after SigProfilerTopography Analyses
def plotFigures(outputDir,jobname,numberofSimulations,multipleTesting,probabilityCalculation):

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
    nucleosomeOccupancyAverageSignalFigures(outputDir,jobname,figureAugmentation,numberofSimulations)
    replicationTimeNormalizedMutationDensityFigures(outputDir,jobname,figureAugmentation,numberofSimulations)
    transcriptionReplicationStrandBiasFigures(outputDir,jobname,figureAugmentation,numberofSimulations)
    processivityFigures(outputDir,jobname,numberofSimulations,multipleTesting,probabilityCalculation)
    ############################################################

##############################################################


# ##############################################################
# #This code for running SigProfilerTopography in local laptop from pycharm ide
# if __name__ == '__main__':
#     genome = 'GRCh37'
#     singlePointMutationsFilename = os.path.join('C:\\','Users','burcak', 'Developer','Python', 'SigProfilerTopography','SigProfilerTopography','input_test','BreastCancer560','560_BRCA_WGS_snps_for_topography.txt')
#     indelsFilename = os.path.join('C:\\','Users','burcak', 'Developer','Python', 'SigProfilerTopography','SigProfilerTopography','input_test','BreastCancer560','merged.txt')
#     outputDir = os.path.join('C:\\','Users','burcak', 'Developer','Python', 'SigProfilerTopography','SigProfilerTopography','output_test')
#     jobname = 'BreastCancer560'
#
#     nucleosomeFilename = os.path.join('C:\\','Users','burcak','Developer','Python','SigProfilerTopography','SigProfilerTopography','lib', 'nucleosome','wgEncodeSydhNsomeGm12878Sig.wig')
#     replicationTimeFilename = os.path.join('C:\\','Users','burcak','Developer','Python','SigProfilerTopography','SigProfilerTopography','lib', 'replication','GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.wig')
#     replicationTimeValleyFilename = os.path.join('C:\\','Users','burcak','Developer','Python','SigProfilerTopography','SigProfilerTopography','lib', 'replication','GSM923442_hg19_wgEncodeUwRepliSeqMcf7ValleysRep1.bed')
#     replicationTimePeakFilename = os.path.join('C:\\','Users','burcak','Developer','Python','SigProfilerTopography','SigProfilerTopography','lib', 'replication','GSM923442_hg19_wgEncodeUwRepliSeqMcf7PkRep1.bed')
#
#     runAnalyses(genome,singlePointMutationsFilename,indelsFilename,outputDir,jobname,nucleosomeFilename,replicationTimeFilename,replicationTimeValleyFilename,replicationTimePeakFilename)
#     plotFigures(outputDir, jobname, 0,'BONFERRONI_CORRECTION','USING_POISSON_DISTRIBUTION')
# ##############################################################

