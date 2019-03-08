# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018 Burcak Otlu

import sys
import os
import shutil
import twobitreader

#############################################################
current_abs_path = os.path.dirname(os.path.realpath(__file__))
commonsPath = os.path.join(current_abs_path,'source','commons')
sys.path.append(commonsPath)
#############################################################

from SigProfilerTopography.source.commons import TopographyCommons
from SigProfilerTopography.source.commons import DataPreparationCommons

from SigProfilerTopography.source.commons import PartitionIndelsData
from SigProfilerTopography.source.commons import PartitionSinglePointMutationsData
from SigProfilerTopography.source.commons import NucleosomeOccupancySignalCountArraysAndFigures
from SigProfilerTopography.source.nucleosomeoccupancy import NucleosomeOccupancyAnalysis
from SigProfilerTopography.source.replicationtime import ReplicationTimeAnalysis
from SigProfilerTopography.source.replicationstrandbias import ReplicationStrandBiasAnalysis
from SigProfilerTopography.source.transcriptionstrandbias import TranscriptionStrandBiasAnalysis
from SigProfilerTopography.source.processivity import ProcessivityAnalysis


from SigProfilerTopography.source.plotting.NucleosomeOccupancyAverageSignalFigures import *
from SigProfilerTopography.source.plotting.ReplicationTimeNormalizedMutationDensityFigures import *
from SigProfilerTopography.source.plotting.TranscriptionReplicationStrandBiasFigures import *
from SigProfilerTopography.source.plotting.ProcessivityFigures import *


############################################################
#Step1 read snp positions
#Step2 read probabilities
#Step3 merge snp positions with probabilities make ready for SigProfilerTopography
#Step4 read indels make ready for SigProfilerTopography
def prepareDataAfterExtractorForTopography(snpsInputFile,indelsInputFile,probabilitiesFile,jobname):

    current_abs_path = os.path.dirname(os.path.realpath(__file__))

    #############################################################################
    #Step1 read snps with genomic positions circa 28 million rows
    snps_df = DataPreparationCommons.readMutationsWithGenomicPositions(snpsInputFile)

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
    probabilities_df = DataPreparationCommons.readProbabilities(probabilitiesFile)
    #############################################################################

    ###################################################
    hg19_genome = twobitreader.TwoBitFile(os.path.join(current_abs_path, LIB, UCSCGENOME, 'hg19.2bit'))
    hg38_genome = twobitreader.TwoBitFile(os.path.join(current_abs_path, LIB, UCSCGENOME, 'hg38.2bit'))
    ###################################################

    #############################################################################
    # Step3
    DataPreparationCommons.mergeSNPsWithSignatureProbabilities(jobname,snps_df,probabilities_df,hg19_genome,hg38_genome)
    #############################################################################

    #############################################################################
    # Step4
    DataPreparationCommons.readIndelsandWriteWithGenomicPositions(jobname,indelsInputFile)
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
    probabilities_df = DataPreparationCommons.readProbabilities(probabilitiesFile)

    if (len(mutationTypes)==2):
        sigProfilerSimulatorSpecificDirName = '%s_simulations_%s_%s_%s' %(jobname,genomeAssembly,mutationTypes[1], mutationTypes[0])
    elif (len(mutationTypes)==1):
        sigProfilerSimulatorSpecificDirName = '%s_simulations_%s_%s' %(jobname,genomeAssembly,mutationTypes[0])

    DataPreparationCommons.prepareSimulationBasedInputFilesForSigProfilerTopography(jobname,genomeAssembly,mutationTypes,sigProfilerSimulatorSpecificDirName,numberofSimulations,probabilities_df)

################################################################



#######################################################
#Run SigProfilerTopography Analyses
def runAnalyses(genome, singlePointMutationsFilename,indelsFilename,jobname,nucleosomeFilename,replicationTimeFilename,replicationTimeValleyFilename,replicationTimePeakFilename):

    #Internally Set
    considerProbabilityInProcessivityAnalysis = True

    current_abs_path = os.path.dirname(os.path.realpath(__file__))
    print('current_abs_path: %s ' % current_abs_path)

    ##############################################
    #Partition the data (Single Point Mutations data and Indels data)
    # Delete the output/jobname/DATA/chrbased if exists
    jobnamePath = os.path.join(current_abs_path,OUTPUT,jobname,DATA,CHRBASED)
    print('sigProTopographyData.py jobnamePath:%s ' %jobnamePath)

    #######################################################
    #Let's not delete the existing data
    # if (os.path.exists(jobnamePath)):
    #     try:
    #         shutil.rmtree(jobnamePath)
    #     except OSError as e:
    #         print('Error: %s - %s.' % (e.filename, e.strerror))
    #######################################################


    if (indelsFilename!=TopographyCommons.NOTSET):
        PartitionIndelsData.partitionIndelsData(jobname,indelsFilename)
    if (singlePointMutationsFilename!=TopographyCommons.NOTSET):
        PartitionSinglePointMutationsData.partitionMutationsData(jobname,singlePointMutationsFilename)
    ##############################################

    #############################################
    print('current_abs_path: %s ' % current_abs_path)
    availableNucleosomeOccupancyFilesPath = os.path.join(current_abs_path,LIB,NUCLEOSOME,AVAILABLENUCLEOSOMEOCCUPANCYFILESNAME)

    if (os.path.exists(availableNucleosomeOccupancyFilesPath)):
        availableNucleosomeOccupancyFilesList = readAsAList(availableNucleosomeOccupancyFilesPath)

    if (nucleosomeFilename not in availableNucleosomeOccupancyFilesList):
        quantileValue = round(float(0.97), 2)
        # PartitionNucleosomeOccupancyData.partitionNucleosomeOccupancyData(jobname,nucleosomeFilename,quantileValue)
        NucleosomeOccupancySignalCountArraysAndFigures.readAllNucleosomeOccupancyDataAndWriteChrBasedSignalCountArrays(genome,quantileValue,nucleosomeFilename)
    #############################################

    ##############################################
    #Please note that singlePointMutationsFilename and indelsFilename are full paths.
    # Get the fienames at the end
    singlePointMutationsFilename = os.path.basename(singlePointMutationsFilename)
    indelsFilename = os.path.basename(indelsFilename)
    ##############################################

    ##############################################
    # NUCLEOSOMEOCCUPANCYANALYSIS
    # Delete the output/jobname/DATA/NUCLEOSOMEOCCUPANCY if exists
    jobnamePath = os.path.join(current_abs_path,OUTPUT, jobname,DATA, NUCLEOSOMEOCCUPANCY)

    ################################################
    if (os.path.exists(jobnamePath)):
        try:
            shutil.rmtree(jobnamePath)
        except OSError as e:
            print('Error: %s - %s.' % (e.filename, e.strerror))
    ################################################

    NucleosomeOccupancyAnalysis.nucleosomeOccupancyAnalysis(genome,jobname,singlePointMutationsFilename,indelsFilename,nucleosomeFilename)
    ###############################################

    # #############################################
    # # REPLICATIONTIME
    # # Delete the output/jobname/DATA/REPLICATIONTIME if exists
    # jobnamePath = os.path.join(current_abs_path, TopographyCommons.OUTPUT, jobname, TopographyCommons.DATA, TopographyCommons.REPLICATIONTIME)
    #
    # # ################################################
    # # if (os.path.exists(jobnamePath)):
    # #     try:
    # #         shutil.rmtree(jobnamePath)
    # #     except OSError as e:
    # #         print('Error: %s - %s.' % (e.filename, e.strerror))
    # # ################################################
    #
    # # ReplicationTime
    # print('current_abs_path: %s ' %current_abs_path)
    # # subprocess.call(['python', os.path.join(current_abs_path,SOURCE,REPLICATIONTIME, 'ReplicationTimeAnalysis.py'),jobname,singlePointMutationsFilename,indelsFilename,replicationTimeFilename])
    # ReplicationTimeAnalysis.replicationTimeAnalysis(jobname,singlePointMutationsFilename,indelsFilename,replicationTimeFilename)
    # ###############################################


    # ###############################################
    # # REPLICATIONSTRANDBIAS
    # # Delete the output/jobname/DATA/REPLICATIONSTRANDBIAS if exists
    # jobnamePath = os.path.join(current_abs_path,TopographyCommons.OUTPUT,jobname,TopographyCommons.DATA,TopographyCommons.REPLICATIONSTRANDBIAS)
    #
    # # ################################################
    # # if (os.path.exists(jobnamePath)):
    # #     try:
    # #         shutil.rmtree(jobnamePath)
    # #     except OSError as e:
    # #         print('Error: %s - %s.' % (e.filename, e.strerror))
    # # ################################################
    #
    # # ReplicationStrandBias
    # # subprocess.call(['python', os.path.join(current_abs_path,SOURCE,REPLICATIONSTRANDBIAS, 'ReplicationStrandBiasAnalysis.py'),jobname,singlePointMutationsFilename,replicationTimeFilename,replicationTimeValleyFilename, replicationTimePeakFilename,'0.5','0.5','0.01'])
    #
    # smoothedWaveletRepliseqDataFilename = replicationTimeFilename
    # valleysBEDFilename = replicationTimeValleyFilename
    # peaksBEDFilename = replicationTimePeakFilename
    #
    # startMutationProbability = round(float(0.5),2)
    # endMutationProbability = round(float(0.5),2)
    # step = round(float(0.01),2)
    #
    # if (singlePointMutationsFilename!=NOTSET):
    #     ReplicationStrandBiasAnalysis.replicationStrandBiasAnalysis(jobname,singlePointMutationsFilename, smoothedWaveletRepliseqDataFilename,valleysBEDFilename, peaksBEDFilename,startMutationProbability,endMutationProbability,step)
    # ###############################################


    # ###############################################
    # # TRANSCRIPTIONSTRANDBIAS
    # # Delete the output/jobname/DATA/TRANSCRIPTIONSTRANDBIAS if exists
    # jobnamePath = os.path.join(current_abs_path,OUTPUT,jobname,DATA,TRANSCRIPTIONSTRANDBIAS)
    #
    # # ################################################
    # # if (os.path.exists(jobnamePath)):
    # #     try:
    # #         shutil.rmtree(jobnamePath)
    # #     except OSError as e:
    # #         print('Error: %s - %s.' % (e.filename, e.strerror))
    # # ################################################
    #
    # # TranscriptionStrandBias
    # # subprocess.call(['python', os.path.join(current_abs_path,SOURCE,TRANSCRIPTIONSTRANDBIAS, 'TranscriptionStrandBiasAnalysis.py'),jobname,singlePointMutationsFilename,'0.5','0.5','0.01'])
    # startMutationProbability = round(float(0.5),2)
    # endMutationProbability = round(float(0.5),2)
    # step = round(float(0.01),2)
    #
    # if (singlePointMutationsFilename!=NOTSET):
    #     TranscriptionStrandBiasAnalysis.transcriptionStrandBiasAnalysis(jobname, singlePointMutationsFilename, startMutationProbability,endMutationProbability,step)
    # ###############################################


    # ###############################################
    # # PROCESSIVITY
    # # Delete the output/jobname/DATA/PROCESSIVITY if exists
    # jobnamePath = os.path.join(current_abs_path,TopographyCommons.OUTPUT,jobname,TopographyCommons.DATA,TopographyCommons.PROCESSIVITY)
    #
    # ################################################
    # # if (os.path.exists(jobnamePath)):
    # #     try:
    # #         shutil.rmtree(jobnamePath)
    # #     except OSError as e:
    # #         print('Error: %s - %s.' % (e.filename, e.strerror))
    # ################################################
    #
    # # TranscriptionStrandBias
    # # subprocess.call(['python', os.path.join(current_abs_path,SOURCE,PROCESSIVITY, 'ProcessivityAnalysis.py'),jobname,singlePointMutationsFilename,considerProbability])
    # if (singlePointMutationsFilename != NOTSET):
    #     ProcessivityAnalysis.processivityAnalysis(jobname,singlePointMutationsFilename,considerProbabilityInProcessivityAnalysis)
    # ###############################################

#######################################################



##############################################################
# BONFERRONI_CORRECTION = 'BONFERRONI_CORRECTION'
# FDR_BH_CORRECTION = 'FDR_BH_CORRECTION'
#
# USING_POISSON_DISTRIBUTION = 'USING_POISSON_DISTRIBUTION'
# USING_NULL_DISTRIBUTION = 'USING_NULL_DISTRIBUTION'
# USING_GAUSSIAN_KDE = 'USING_GAUSSIAN_KDE'

#Plot Figures for the attainded data after SigProfilerTopography Analyses
def plotFigures(jobname,numberofSimulations,multipleTesting,probabilityCalculation):

    #Internally Set
    figureAugmentation = 'noaugmentation'

    current_abs_path = os.path.dirname(os.path.realpath(__file__))
    print('current_abs_path: %s ' % current_abs_path)

    jobnamePath = os.path.join(current_abs_path,OUTPUT,jobname,FIGURE)
    print('SigProfilerTopography.py jobnamePath:%s ' %jobnamePath)
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


    # Notice that subprocess expects strings in the list
    # Plotting
    # subprocess.call(['python', os.path.join(current_abs_path, SOURCE, PLOTTING, 'NucleosomeOccupancyAverageSignalFigures.py'),jobname,figureAugmentation,str(numberofSimulations)])
    # subprocess.call(['python', os.path.join(current_abs_path, SOURCE, PLOTTING, 'ReplicationTimeNormalizedMutationDensityFigures.py'),jobname,figureAugmentation,str(numberofSimulations)])
    # subprocess.call(['python', os.path.join(current_abs_path, SOURCE, PLOTTING, 'TranscriptionReplicationStrandBiasFigures.py'),jobname,figureAugmentation,str(numberofSimulations)])
    # subprocess.call(['python', os.path.join(current_abs_path, SOURCE, PLOTTING, 'ProcessivityFigures.py'),jobname,str(numberofSimulations),multipleTesting,probabilityCalculation])

    nucleosomeOccupancyAverageSignalFigures(jobname,figureAugmentation,numberofSimulations)
    # replicationTimeNormalizedMutationDensityFigures(jobname,figureAugmentation,numberofSimulations)
    # transcriptionReplicationStrandBiasFigures(jobname,figureAugmentation,numberofSimulations)
    # processivityFigures(jobname,numberofSimulations,multipleTesting,probabilityCalculation)

##############################################################
