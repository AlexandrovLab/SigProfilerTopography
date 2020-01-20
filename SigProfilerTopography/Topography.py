# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018-2020 Burcak Otlu

# #############################################################
# import sys
# import os
# current_abs_path = os.path.dirname(os.path.realpath(__file__))
# commonsPath = os.path.join(current_abs_path,'commons')
# sys.path.append(commonsPath)
# #############################################################

import time
import numpy as np
import pandas as pd
import scipy
import statsmodels
import matplotlib as plt

import shutil
import platform
import multiprocessing

import SigProfilerMatrixGenerator as matrix_generator
MATRIX_GENERATOR_PATH=matrix_generator.__path__[0]

from SigProfilerTopography import version as topography_version
from SigProfilerMatrixGenerator import version as matrix_generator_version
from SigProfilerSimulator import version as simulator_version

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerSimulator import SigProfilerSimulator as simulator

from SigProfilerTopography.source.commons.TopographyCommons import readProbabilities
from SigProfilerTopography.source.commons.TopographyCommons import readChrBasedMutationsMergeWithProbabilitiesAndWrite

from SigProfilerTopography.source.commons.TopographyCommons import DATA
from SigProfilerTopography.source.commons.TopographyCommons import FIGURE
from SigProfilerTopography.source.commons.TopographyCommons import ALL
from SigProfilerTopography.source.commons.TopographyCommons import SAMPLES
from SigProfilerTopography.source.commons.TopographyCommons import SAMPLE

from SigProfilerTopography.source.commons.TopographyCommons import K562
from SigProfilerTopography.source.commons.TopographyCommons import MCF7

from SigProfilerTopography.source.commons.TopographyCommons import GM12878_NUCLEOSOME_OCCUPANCY_FILE
from SigProfilerTopography.source.commons.TopographyCommons import K562_NUCLEOSOME_OCCUPANCY_FILE

from SigProfilerTopography.source.commons.TopographyCommons import getNucleosomeFile
from SigProfilerTopography.source.commons.TopographyCommons import getReplicationTimeFiles

from SigProfilerTopography.source.commons.TopographyCommons import available_nucleosome_biosamples
from SigProfilerTopography.source.commons.TopographyCommons import available_replication_time_biosamples

from SigProfilerTopography.source.commons.TopographyCommons import NUCLEOSOME
from SigProfilerTopography.source.commons.TopographyCommons import REPLICATIONTIME
from SigProfilerTopography.source.commons.TopographyCommons import REPLICATIONSTRANDBIAS
from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIPTIONSTRANDBIAS
from SigProfilerTopography.source.commons.TopographyCommons import PROCESSIVITY
from SigProfilerTopography.source.commons.TopographyCommons import EPIGENOMICS
from SigProfilerTopography.source.commons.TopographyCommons import REPLICATION

from SigProfilerTopography.source.commons.TopographyCommons import EPIGENOMICSOCCUPANCY
from SigProfilerTopography.source.commons.TopographyCommons import NUCLEOSOMEOCCUPANCY

from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_HISTONE_OCCUPANCY_FILE1
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_HISTONE_OCCUPANCY_FILE2
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_HISTONE_OCCUPANCY_FILE3
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_HISTONE_OCCUPANCY_FILE4
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_HISTONE_OCCUPANCY_FILE5
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_HISTONE_OCCUPANCY_FILE6

from SigProfilerTopography.source.commons.TopographyCommons import BIOSAMPLE_UNDECLARED

from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC

from SigProfilerTopography.source.commons.TopographyCommons import SBS96
from SigProfilerTopography.source.commons.TopographyCommons import ID
from SigProfilerTopography.source.commons.TopographyCommons import DBS

from SigProfilerTopography.source.commons.TopographyCommons import SUBS
from SigProfilerTopography.source.commons.TopographyCommons import INDELS
from SigProfilerTopography.source.commons.TopographyCommons import DINUCS

from SigProfilerTopography.source.commons.TopographyCommons import SBS_CONTEXTS
from SigProfilerTopography.source.commons.TopographyCommons import SNV

from SigProfilerTopography.source.commons.TopographyCommons import CHRBASED

from SigProfilerTopography.source.commons.TopographyCommons import GRCh37
from SigProfilerTopography.source.commons.TopographyCommons import GRCh38
from SigProfilerTopography.source.commons.TopographyCommons import MM9
from SigProfilerTopography.source.commons.TopographyCommons import MM10

from SigProfilerTopography.source.commons.TopographyCommons import ONE_DIRECTORY_UP
from SigProfilerTopography.source.commons.TopographyCommons import LIB
from SigProfilerTopography.source.commons.TopographyCommons import UCSCGENOME

from SigProfilerTopography.source.commons.TopographyCommons import WIG
from SigProfilerTopography.source.commons.TopographyCommons import BED

from SigProfilerTopography.source.commons.TopographyCommons import current_abs_path

from SigProfilerTopography.source.commons.TopographyCommons import getChromSizesDict
from SigProfilerTopography.source.commons.TopographyCommons import getShortNames
from SigProfilerTopography.source.commons.TopographyCommons import doesSimulationsAlreadyExists
from SigProfilerTopography.source.commons.TopographyCommons import copyMafFiles
from SigProfilerTopography.source.commons.TopographyCommons import fillCutoff2Signature2PropertiesListDictionary
from SigProfilerTopography.source.commons.TopographyCommons import fill_mutations_dictionaries_write
from SigProfilerTopography.source.commons.TopographyCommons import appendDictionaryUnderDataDirectory

from SigProfilerTopography.source.commons.TopographyCommons import MutationType2NumberofMutatiosDictFilename
from SigProfilerTopography.source.commons.TopographyCommons import readDictionary

from SigProfilerTopography.source.nucleosomeoccupancy.ChrBasedSignalArrays import readBEDandWriteChromBasedSignalArrays

from SigProfilerTopography.source.commons.TopographyCommons import SubsSignature2PropertiesListDictFilename
from SigProfilerTopography.source.commons.TopographyCommons import IndelsSignature2PropertiesListDictFilename
from SigProfilerTopography.source.commons.TopographyCommons import DinucsSignature2PropertiesListDictFilename

from SigProfilerTopography.source.nucleosomeoccupancy.NucleosomeOccupancyAnalysis import occupancyAnalysis
from SigProfilerTopography.source.replicationtime.ReplicationTimeAnalysis import replicationTimeAnalysis
from SigProfilerTopography.source.replicationstrandbias.ReplicationStrandBiasAnalysis import replicationStrandBiasAnalysis
from SigProfilerTopography.source.transcriptionstrandbias.TranscriptionStrandBiasAnalysis import transcriptionStrandBiasAnalysis
from SigProfilerTopography.source.processivity.ProcessivityAnalysis import processivityAnalysis

from SigProfilerTopography.source.plotting.OccupancyAverageSignalFigures import occupancyAverageSignalFigures
from SigProfilerTopography.source.plotting.OccupancyAverageSignalFigures import plot_heatmaps
from SigProfilerTopography.source.plotting.ReplicationTimeNormalizedMutationDensityFigures import replicationTimeNormalizedMutationDensityFigures
from SigProfilerTopography.source.plotting.TranscriptionReplicationStrandBiasFigures import transcriptionReplicationStrandBiasFigures
from SigProfilerTopography.source.plotting.ProcessivityFigures import processivityFigures


############################################################
#Can be move to DataPreparationCommons under /source/commons
#read chr based dinucs (provided by SigProfilerMatrixGenerator) and merge with probabilities (provided by SigProfilerTopography)
def prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList,inputDir,outputDir,jobname,mutation_type_context,mutations_probabilities_file_path,startSimNum, endSimNum,partialDirname,PCAWG):

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
        mutations_probabilities_df = readProbabilities(mutations_probabilities_file_path)
        ##########################################################################################

        # print('mutations_probabilities_df.head()')
        # print(mutations_probabilities_df.head())

        # print('mutations_probabilities_df.columns.values')
        # print(mutations_probabilities_df.columns.values)

        ##########################################################################################
        #Step1 SigProfilerTopography Python Package
        #For Release we will use SAMPLE as it is, no change in SAMPLE column is needed.

        # For PCAWG_Matlab
        # This statement below is customized for  PCAWG_Matlab
        # To get rid of inconsistent cancer type names in sample column of chrbased mutation files and probabilities files
        # Breast-LobularCA_SP124191
        if PCAWG:
            mutations_probabilities_df[SAMPLE] = mutations_probabilities_df[SAMPLE].str.split('_',expand=True)[1]
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
                    inputList.append(PCAWG)
                    poolInputList.append(inputList)

        #TODO Right now it uses only one processor
        #TODO I guess this happens when sim data is big
        #TODO Use imap_unordered with big chunksize and monitor performance
        pool.map(readChrBasedMutationsMergeWithProbabilitiesAndWrite, poolInputList)

        pool.close()
        pool.join()

    else:
        #For Information
        print('%s does not exist.' %(mutations_probabilities_file_path))
############################################################

#######################################################
#JAN 9, 2020
def check_download_replication_time_files(replication_time_signal_file,replication_time_valley_file,replication_time_peak_file):
    current_abs_path = os.path.dirname(os.path.abspath(__file__))
    # print(current_abs_path)

    #These are currently full path, therefore convert them to filename
    replication_time_signal_file=os.path.basename(replication_time_signal_file)
    replication_time_valley_file=os.path.basename(replication_time_valley_file)
    replication_time_peak_file=os.path.basename(replication_time_peak_file)

    os.makedirs(os.path.join(current_abs_path,'lib','replication'),exist_ok=True)
    lib_replication_path = os.path.join(current_abs_path,'lib','replication')

    if os.path.isabs(lib_replication_path):
        # print('%s an absolute path.' %(lib_replication_path))
        os.chdir(lib_replication_path)

        replication_time_signal_file_path= os.path.join(lib_replication_path,replication_time_signal_file)
        replication_time_valley_file_path= os.path.join(lib_replication_path,replication_time_valley_file)
        replication_time_peak_file_path= os.path.join(lib_replication_path,replication_time_peak_file)

        if not os.path.exists(replication_time_signal_file_path):
            print('Does not exists: %s' %(replication_time_signal_file_path))
            try:
                # print('Downloading %s_signal_wgEncodeSydhNsome_%sSig.npy under %s' %(chrLong,cell_line,chrbased_npy_array_path))
                print('Downloading %s under %s' % (replication_time_signal_file, lib_replication_path))

                #wget -c Continue getting a partially-downloaded file
                #wget -nc  If a file is downloaded more than once in the same directory, the local file will be clobbered, or overwritten
                # cmd="bash -c '" + 'wget -r -l1 -c -nc --no-parent -nd -P ' + chrombased_npy_path + ' ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerTopography/lib/nucleosome/chrbased/' + filename + "'"
                cmd="bash -c '" + 'wget -r -l1 -c -nc --no-parent -nd ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerTopography/lib/replication/' + replication_time_signal_file + "'"
                # print(cmd)
                os.system(cmd)
            except:
                # print("The UCSD ftp site is not responding...pulling from sanger ftp now.")
                print("The ftp://alexandrovlab-ftp.ucsd.edu site is not responding...")

        if not os.path.exists(replication_time_valley_file_path):
            print('Does not exists: %s' %(replication_time_valley_file_path))
            try:
                # print('Downloading %s_signal_wgEncodeSydhNsome_%sSig.npy under %s' %(chrLong,cell_line,chrbased_npy_array_path))
                print('Downloading %s under %s' % (replication_time_valley_file, lib_replication_path))

                #wget -c Continue getting a partially-downloaded file
                #wget -nc  If a file is downloaded more than once in the same directory, the local file will be clobbered, or overwritten
                # cmd="bash -c '" + 'wget -r -l1 -c -nc --no-parent -nd -P ' + chrombased_npy_path + ' ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerTopography/lib/nucleosome/chrbased/' + filename + "'"
                cmd="bash -c '" + 'wget -r -l1 -c -nc --no-parent -nd ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerTopography/lib/replication/' + replication_time_valley_file + "'"
                # print(cmd)
                os.system(cmd)
            except:
                # print("The UCSD ftp site is not responding...pulling from sanger ftp now.")
                print("The ftp://alexandrovlab-ftp.ucsd.edu site is not responding...")

        if not os.path.exists(replication_time_peak_file_path):
            print('Does not exists: %s' %(replication_time_peak_file_path))
            try:
                # print('Downloading %s_signal_wgEncodeSydhNsome_%sSig.npy under %s' %(chrLong,cell_line,chrbased_npy_array_path))
                print('Downloading %s under %s' % (replication_time_peak_file, lib_replication_path))

                #wget -c Continue getting a partially-downloaded file
                #wget -nc  If a file is downloaded more than once in the same directory, the local file will be clobbered, or overwritten
                # cmd="bash -c '" + 'wget -r -l1 -c -nc --no-parent -nd -P ' + chrombased_npy_path + ' ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerTopography/lib/nucleosome/chrbased/' + filename + "'"
                cmd="bash -c '" + 'wget -r -l1 -c -nc --no-parent -nd ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerTopography/lib/replication/' + replication_time_peak_file + "'"
                # print(cmd)
                os.system(cmd)
            except:
                # print("The UCSD ftp site is not responding...pulling from sanger ftp now.")
                print("The ftp://alexandrovlab-ftp.ucsd.edu site is not responding...")

    else:
        #It has to be an absolute path
        print('%s is not an absolute path.' %(lib_replication_path))

    #go back
    os.chdir(current_abs_path)

#######################################################

#######################################################
#Nov25, 2019
# Download nucleosome occupancy chr based npy files from ftp alexandrovlab if they do not exists
# We are using this function if user is using our available nucleosome data for GM12878 adnd K562 cell lines
def check_download_chrbased_npy_nuclesome_files(nucleosome_file,chromNamesList):

    current_abs_path = os.path.dirname(os.path.abspath(__file__))
    # print(current_abs_path)

    os.makedirs(os.path.join(current_abs_path,'lib','nucleosome','chrbased'),exist_ok=True)
    chrombased_npy_path = os.path.join(current_abs_path,'lib','nucleosome','chrbased')
    # print(chrombased_npy_path)

    if os.path.isabs(chrombased_npy_path):
        # print('%s an absolute path.' %(chrombased_npy_path))
        os.chdir(chrombased_npy_path)

        nucleosome_filename_wo_extension = os.path.splitext(os.path.basename(nucleosome_file))[0]

        for chrLong in chromNamesList:
            # GM12878 and K562 comes from woman samples therefore there is no chrY
            if chrLong != 'chrY':
                # filename = '%s_signal_wgEncodeSydhNsome%sSig.npy' %(chrLong,cell_line)
                filename = '%s_signal_%s.npy' % (chrLong, nucleosome_filename_wo_extension)

                chrbased_npy_array_path = os.path.join(chrombased_npy_path, filename)
                if not os.path.exists(chrbased_npy_array_path):
                    print('Does not exists: %s' % (chrbased_npy_array_path))
                    try:
                        # print('Downloading %s_signal_wgEncodeSydhNsome_%sSig.npy under %s' %(chrLong,cell_line,chrbased_npy_array_path))
                        print('Downloading %s_signal_%s.npy under %s' % (
                        chrLong, nucleosome_filename_wo_extension, chrbased_npy_array_path))

                        # wget -c Continue getting a partially-downloaded file
                        # wget -nc  If a file is downloaded more than once in the same directory, the local file will be clobbered, or overwritten
                        # cmd="bash -c '" + 'wget -r -l1 -c -nc --no-parent -nd -P ' + chrombased_npy_path + ' ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerTopography/lib/nucleosome/chrbased/' + filename + "'"
                        cmd = "bash -c '" + 'wget -r -l1 -c -nc --no-parent -nd ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerTopography/lib/nucleosome/chrbased/' + filename + "'"
                        # print(cmd)
                        os.system(cmd)
                    except:
                        # print("The UCSD ftp site is not responding...pulling from sanger ftp now.")
                        print("The UCSD ftp site is not responding...")

    else:
        #It has to be an absolute path
        print('%s is not an absolute path.' %(chrombased_npy_path))


    #go back
    os.chdir(current_abs_path)
#######################################################


#######################################################

#######################################################
def runOccupancyAnalyses(genome,outputDir,jobname,numofSimulations,sample_based,library_file_with_path,library_file_memo,chromSizesDict,chromNamesList,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,computation_type,occupancy_type,plusorMinus,verbose):

    #We have to exclude for Topography provided nucleosome occupancy files
    exclude_from_check=False

    if (os.path.basename(library_file_with_path)==GM12878_NUCLEOSOME_OCCUPANCY_FILE) or (os.path.basename(library_file_with_path)==K562_NUCLEOSOME_OCCUPANCY_FILE):
        exclude_from_check=True

    if (not exclude_from_check) and (not os.path.exists(library_file_with_path)):
        print('There is no such file under %s' %(library_file_with_path))

    occupancyAnalysis(genome,computation_type,occupancy_type,sample_based,plusorMinus,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,library_file_with_path,library_file_memo,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,verbose)
#######################################################


#######################################################
def runReplicationTimeAnalysis(genome,outputDir,jobname,numofSimulations,sample_based,replicationTimeFilename,chromSizesDict,chromNamesList,computation_type,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,verbose,matrix_generator_path):

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

    #Option2: Load offline prepared np arrays during runtime managby ed by replication_time_np_arrays_fill_runtime=False
    #Option2: Fill np array during runtime managed by replication_time_np_arrays_fill_runtime=True
    replicationTimeAnalysis(computation_type,sample_based,genome,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,replicationTimeFilename,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,verbose,matrix_generator_path)
    ###############################################

#######################################################


#######################################################
def runReplicationStrandBiasAnalysis(outputDir,jobname,numofSimulations,sample_based,replicationTimeFilename,replicationTimeValleyFilename,replicationTimePeakFilename,chromSizesDict,chromNamesList,computation_type,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,verbose):

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

    replicationStrandBiasAnalysis(computation_type,sample_based,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,smoothedWaveletRepliseqDataFilename,valleysBEDFilename,peaksBEDFilename,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,verbose)
    ###############################################

#######################################################

#######################################################
def runTranscriptionStradBiasAnalysis(outputDir,jobname,numofSimulations,sample_based,chromSizesDict,chromNamesList,computation_type,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,verbose):
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

    transcriptionStrandBiasAnalysis(computation_type,sample_based,chromSizesDict,chromNamesList,outputDir,jobname,numofSimulations,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,verbose)
    ###############################################
#######################################################


#######################################################
def runProcessivityAnalysis(mutation_types_contexts,outputDir,jobname,numofSimulations,chromNamesList,computation_type,signature2PropertiesListDict,verbose):
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
    # computation_type=USING_APPLY_ASYNC
    processivityAnalysis(mutation_types_contexts,chromNamesList,computation_type,outputDir,jobname,numofSimulations,considerProbabilityInProcessivityAnalysis,signature2PropertiesListDict,verbose)
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
                sbs_probabilities=None,
                id_probabilities= None,
                dbs_probabilities=None,
                mutation_types_contexts=None,
                epigenomics_files=[DEFAULT_HISTONE_OCCUPANCY_FILE1,DEFAULT_HISTONE_OCCUPANCY_FILE2,DEFAULT_HISTONE_OCCUPANCY_FILE3,DEFAULT_HISTONE_OCCUPANCY_FILE4,DEFAULT_HISTONE_OCCUPANCY_FILE5,DEFAULT_HISTONE_OCCUPANCY_FILE6],
                epigenomics_files_memos=None,
                epigenomics_biosamples=None,
                nucleosome_biosample=K562,
                nucleosome_file=None,
                replication_time_biosample=MCF7,
                replication_time_signal_file=None,
                replication_time_valley_file=None,
                replication_time_peak_file=None,
                computation_type=USING_APPLY_ASYNC,
                epigenomics=False,
                nucleosome=False,
                replication_time=False,
                strand_bias=False,
                processivity=False,
                sample_based=False,
                new_simulations_enforced=True,
                plot_figures=True,
                full_mode=True,
                necessary_dictionaries_already_exists=False,
                average_probability=0.9,
                num_of_sbs_required=5000,
                num_of_id_required=1000,
                num_of_dbs_required=200,
                plusorMinus_epigenomics=2000,
                plusorMinus_nucleosome=1000,
                verbose=False,
                matrix_generator_path=MATRIX_GENERATOR_PATH,
                PCAWG=False,
                DATA_IS_READY_PLOT_FIGURES=False):

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

    ###################################################
    if mutation_types_contexts is None:
        mutation_types_contexts=[]
        if (sbs_probabilities is not None) and (os.path.exists(sbs_probabilities)):
            mutation_types_contexts.append(SBS96)
        if (id_probabilities is not None) and (os.path.exists(id_probabilities)):
            mutation_types_contexts.append(ID)
        if (dbs_probabilities is not None) and (os.path.exists(dbs_probabilities)):
            mutation_types_contexts.append(DBS)
    ###################################################

    print('#################################################################################')
    # print('--- %s' %platform.platform())
    # print('--- %s' %platform.system())
    #print("--- Operating System: %s" %(platform.uname()[0]))
    print("--- SigProfilerTopography starts")
    print('#################################################################################')

    print('#################################################################################')
    print("--- Operating System: %s" %(platform.platform()))
    print("--- Release: %s" %platform.uname()[2])
    print("--- Version: %s" %platform.uname()[3])
    print("--- Nodename: %s" %platform.uname()[1])
    print('#################################################################################')

    print('#################################################################################')
    print("--- Python and Package Versions")
    print("--- Python Version: %s" %(str(platform.sys.version_info.major) + "." + str(platform.sys.version_info.minor) + "." + str(platform.sys.version_info.micro)))
    print('--- SigProfilerTopography Version:%s' % topography_version.version)
    print("--- SigProfilerMatrixGenerator Version: %s" %matrix_generator_version.version)
    print("--- SigProfilerSimulator version: %s" %simulator_version.version)
    print("--- pandas version: %s" %pd.__version__)
    print("--- numpy version: %s" %np.__version__)
    print("--- statsmodels version: %s" %statsmodels.__version__)
    print("--- scipy version: %s" %scipy.__version__)
    print("--- matplotlib version: %s" %plt.__version__)
    print('#################################################################################\n')

    print('#################################################################################')
    print('--- SigProfilerTopography parameters')
    print('--- Genome: %s' %(genome))
    print('--- inputDir:%s' %inputDir)
    print('--- outputDir:%s' %outputDir)
    print('--- jobname:%s' %jobname)
    print('--- numofSimulations:%d' %numofSimulations)
    print('--- epigenomics_files:%s' %epigenomics_files)
    print('--- epigenomics_files_memos:%s' %epigenomics_files_memos)
    print('--- epigenomics_biosamples:%s' %epigenomics_biosamples)

    print('--- nucleosome_biosample:%s' %nucleosome_biosample)
    print('--- nucleosome_file:%s' % nucleosome_file)

    print('--- replication_time_biosample:%s' % replication_time_biosample)
    print('--- replication_time_signal_file:%s' % replication_time_signal_file)
    print('--- replication_time_valley_file:%s' % replication_time_valley_file)
    print('--- replication_time_peak_file:%s' % replication_time_peak_file)

    print('--- mutation_types_contexts:%s' %mutation_types_contexts)
    print('--- computation_type:%s' %computation_type)
    if sample_based:
        if epigenomics:
            print('--- Epigenomics Sample Based Analysis.')
        if nucleosome:
            print('--- Nucleosome Sample Based Analysis.')
        if replication_time:
            print('--- Replication Time Sample Based Analysis.')
        if strand_bias:
            print('--- Strand Bias Sample Based Analysis.')
        if processivity:
            print('--- Processivity Analysis.')
    else:
        if epigenomics:
            print('--- Epigenomics Analysis.')
        if nucleosome:
            print('--- Nucleosome Analysis.')
        if replication_time:
            print('--- Replication Time Analysis.')
        if strand_bias:
            print('--- Strand Bias Analysis.')
        if processivity:
            print('--- Processivity Analysis.')
    print('--- new_simulations_enforced:%s' %new_simulations_enforced)
    print('--- plot_figures:%s' %plot_figures)
    print('--- average mutation probability required %0.2f' %average_probability)
    print('--- minimum number of sbs mutations required: %d' %num_of_sbs_required)
    print('--- minimum number of id mutations required: %d' %num_of_id_required)
    print('--- minimum number of dbs mutations required: %d' %num_of_dbs_required)
    if epigenomics:
        print('--- number of bases considered before and after mutation start for epigenomics analysis: %d' %plusorMinus_epigenomics)
    if nucleosome:
        print('--- number of bases considered before and after mutation start for nucleosome occupancy analysis: %d' %plusorMinus_nucleosome)
    print('#################################################################################\n')

    print('#################################################################################')
    numofProcesses = multiprocessing.cpu_count()
    print('--- numofProcesses for multiprocessing: %d' %numofProcesses)
    print('#################################################################################\n')

    #################################################################################
    print('#################################################################################')
    print('--- For Genome: %s' %(genome))
    print('--- Chromosome names: %s' %(chromNamesList))
    print('--- Chromosome short names: %s' % (chromShortNamesList))
    current_abs_path = os.path.dirname(os.path.realpath(__file__))
    print('--- current_abs_path: %s ' % current_abs_path)
    print('#################################################################################\n')
    #################################################################################

    #######################################################################################################
    #####################################  Full mode starts ###############################################
    #######################################################################################################
    if (full_mode==True):

        ###################################################################################################
        #######################  SigProfilerMatrixGenerator for original data starts ######################
        ###################################################################################################
        # Run MatrixGenerator for original data: this call prepares chrBased input files for original data with mutation contexts
        print('#################################################################################')
        print('--- SigProfilerMatrixGenerator for original data')
        start_time = time.time()

        print('For original data inputDir:%s' % (inputDir))
        matrices = matGen.SigProfilerMatrixGeneratorFunc(jobname, genome, inputDir, plot=False, seqInfo=True)
        # print('matrices')
        # print(matrices)

        # original matrix generator chrbased data will be under inputDir/output/vcf_files/SNV
        # original matrix generator chrbased data will be under inputDir/output/vcf_files/DBS
        # original matrix generator chrbased data will be under inputDir/output/vcf_files/ID

        print("--- SigProfilerMatrixGenerator for original data: %s seconds ---" % (time.time() - start_time))
        print("--- SigProfilerMatrixGenerator for original data: %f minutess ---" % float((time.time() - start_time) / 60))
        print('#################################################################################\n')
        ###################################################################################################
        #######################  SigProfilerMatrixGenerator for original data ends ########################
        ###################################################################################################

        ####################################################################################################################
        ##################  Merge original chr based files with Mutation Probabilities starts ##############################
        ####################################################################################################################
        print('#################################################################################')
        print('--- Merge original chr based files with Mutation Probabilities starts')
        print('#################################################################################')
        startSimNum = 0
        endSimNum = 0
        start_time = time.time()
        # SBS
        for mutation_type_context in mutation_types_contexts:
            if (mutation_type_context in SBS_CONTEXTS) and (sbs_probabilities is not None):
                print('--- Merge %s mutations with probabilities for %s' % (
                mutation_type_context, sbs_probabilities))
                prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList, inputDir,
                                                                                   outputDir,
                                                                                   jobname, mutation_type_context,
                                                                                   sbs_probabilities,
                                                                                   startSimNum,
                                                                                   endSimNum, SNV,PCAWG)

        # ID
        if ((ID in mutation_types_contexts) and (id_probabilities is not None)):
            print('--- Merge %s mutations with probabilities for %s' % (ID, id_probabilities))
            prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList, inputDir, outputDir,
                                                                               jobname, ID, id_probabilities,
                                                                               startSimNum, endSimNum, ID,PCAWG)

        # DBS
        if ((DBS in mutation_types_contexts) and (dbs_probabilities is not None)):
            print('--- Merge %s mutations with probabilities for %s' % (DBS, dbs_probabilities))
            prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList, inputDir, outputDir,
                                                                               jobname, DBS,
                                                                               dbs_probabilities,
                                                                               startSimNum, endSimNum, DBS,PCAWG)

        print("--- Merge original chr based files with Mutation Probabilities: %s seconds" % (time.time() - start_time))
        print("--- Merge original chr based files with Mutation Probabilities: %f minutes" % (float((time.time() - start_time) / 60)))
        print('--- Merge original chr based files with Mutation Probabilities ends')
        print('#################################################################################\n')
        ####################################################################################################################
        ##################  Merge original chr based files with Mutation Probabilities ends ################################
        ####################################################################################################################

        ###################################################################################################################
        ###################################################################################################################
        ###################################################################################################################
        existsSimulations= doesSimulationsAlreadyExists(outputDir,jobname,numofSimulations)

        if existsSimulations:
            print('#################################################################################')
            print('--- %d simulations already exists' %(numofSimulations))
            if new_simulations_enforced:
                print('--- new_simulations_enforced:%s' %(new_simulations_enforced))
                print('--- New simulations will be generated')
            else:
                print('--- new_simulations_enforced:%s' %(new_simulations_enforced))
                print('--- Existing simulations will be used')
            print('#################################################################################\n')

        if ((numofSimulations>0) and ((new_simulations_enforced) or (not existsSimulations))):

            ###################################################################################################
            ############################  SigProfilerSimulator for n simulations starts #######################
            ###################################################################################################
            print('#################################################################################')
            print('--- SigProfilerSimulator for %d simulations starts' %(numofSimulations))
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
                print('--- SigProfilerSimulator is running for %s' %(mutation_type_context))
                simulator.SigProfilerSimulator(jobname, inputDir, genome, mutation_type_context_for_simulator,simulations=numofSimulations,chrom_based=True)

            print("--- SigProfilerSimulator for %d simulations: %s seconds" %(numofSimulations,(time.time() -  start_time)))
            print("--- SigProfilerSimulator for %d simulations: %f minutes" %(numofSimulations,float((time.time()-start_time)/60)))
            print('--- SigProfilerSimulator for %d simulations ends' %(numofSimulations))
            print('#################################################################################\n')
            ###################################################################################################
            ############################  SigProfilerSimulator for n simulations ends #########################
            ###################################################################################################

            ###################################################################################################
            ########################### Create simN directories for MatrixGenerator starts ####################
            ###################################################################################################
            print('#################################################################################')
            print('--- Create directories for %d simulations under inputDir/output/simulations/' %(numofSimulations))
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

            for mutation_type_context in mutation_types_contexts:
                # Simulator creates one maf file for each simulation for each mutation context
                # Simulator creates maf files under inputDir/output/simulations/jobname_simulations_GRCh37_96
                # Simulator creates maf files under inputDir/output/simulations/jobname_simulations_GRCh37_ID
                # Simulator creates maf files under inputDir/output/simulations/jobname_simulations_GRCh37_DBS
                dirName = '%s_simulations_%s_%s' %(jobname, genome,mutation_type_context)
                copyFromDir = os.path.join(inputDir,'output','simulations',dirName)
                copyToMainDir= os.path.join(inputDir,'output','simulations')

                # Topography copies these maf files  to inputDir/output/simulations/simX/mutation_type_context/X.maf
                copyMafFiles(copyFromDir,copyToMainDir,mutation_type_context,numofSimulations)
            print("--- Create directories and copy files: %s seconds ---" %(time.time()-start_time))
            print("--- Create directories and copy files: %f minutes ---" %(float((time.time()-start_time)/60)))
            print('#################################################################################\n')
            ###################################################################################################
            ########################### Create simN directories for MatrixGenerator ends ######################
            ###################################################################################################

            ###################################################################################################
            #Important note: Separate directory creation is necessary for Matrix Generator
            #inputDir/output/simulations/simX/96/X.maf
            #inputDir/output/simulations/simX/ID/X.maf
            #inputDir/output/simulations/simX/DBS/X.maf

            #enables MatrixGenerator to create chr based simulated data files under
            #sim1 matrix generator chrbased data will be under inputDir/output/simulations/simX/96/output/vcf_files/SNV
            #sim1 matrix generator chrbased data will be under inputDir/output/simulations/simX/ID/output/vcf_files/ID
            #sim1 matrix generator chrbased data will be under inputDir/output/simulations/simX/DBS/output/vcf_files/DBS

            #otherwise all simulations maf files will be under
            #inputDir/output/simulations/Skin-Melanoma_simulations_GRCh37_96
            #inputDir/output/simulations/Skin-Melanoma_simulations_GRCh37_DBS
            #inputDir/output/simulations/Skin-Melanoma_simulations_GRCh37_ID

            #Therefore running MatrixGenerator for each simulation will not be possible.
            ###################################################################################################

            ###################################################################################################
            ####################### Run MatrixGenerator for each simulation starts ############################
            ###################################################################################################
            print('#################################################################################')
            print('--- Run SigProfilerMatrixGenerator for each simulation starts')
            start_time = time.time()
            for simNum in range(1,numofSimulations+1):
                simName = 'sim%d' %(simNum)
                #For each simulation we are calling matrix generator separately for each mutation type context

                print('--- SigProfilerMatrixGenerator is run for %s starts' %(simName))
                for mutation_type_context in mutation_types_contexts:
                    simInputDir=  os.path.join(inputDir,'output','simulations',simName,mutation_type_context)
                    print('For %s: %s simInputDir:%s' %(mutation_type_context,simName,simInputDir))
                    matrices = matGen.SigProfilerMatrixGeneratorFunc(jobname,genome,simInputDir,plot=False, seqInfo=True)
                    # print('matrices')
                    # print(matrices)
                    print('#####################################')
                print('--- SigProfilerMatrixGenerator is run for %s ends\n' % (simName))
            #sim1 matrix generator chrbased data will be under inputDir/output/simulations/sim1/96/output/vcf_files/SNV
            #sim1 matrix generator chrbased data will be under inputDir/output/simulations/sim1/ID/output/vcf_files/ID
            #sim1 matrix generator chrbased data will be under inputDir/output/simulations/sim1/DBS/output/vcf_files/DBS

            #simN matrix generator chrbased data will be under inputDir/output/simulations/simN/96/output/vcf_files/SNV
            #simN matrix generator chrbased data will be under inputDir/output/simulations/simN/ID/output/vcf_files/ID
            #simN matrix generator chrbased data will be under inputDir/output/simulations/simN/DBS/output/vcf_files/DBS
            print("--- Run MatrixGenerator for each simulation: %s seconds" %(time.time()-start_time))
            print("--- Run MatrixGenerator for each simulation: %f minutes" %(float((time.time()-start_time)/60)))
            print('--- Run SigProfilerMatrixGenerator for each simulation ends')
            print('#################################################################################\n')
            ###################################################################################################
            ####################### Run MatrixGenerator for each simulation ends ##############################
            ###################################################################################################

            ####################################################################################################################
            ##################  Merge simulations chr based files with Mutation Probabilities starts ###########################
            ####################################################################################################################
            print('#################################################################################')
            print('--- Merge simulations chr based files with Mutation Probabilities starts')
            print('#################################################################################')
            startSimNum=1
            endSimNum=numofSimulations
            start_time = time.time()
            #SBS
            for mutation_type_context in mutation_types_contexts:
                if (mutation_type_context in SBS_CONTEXTS) and (sbs_probabilities is not None):
                    print('--- Merge %s mutations with probabilities for %s' %(mutation_type_context,sbs_probabilities))
                    prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList,inputDir,outputDir,jobname,mutation_type_context,sbs_probabilities,startSimNum,endSimNum,'SNV',PCAWG)

            #ID
            if ((ID in mutation_types_contexts) and (id_probabilities is not None)):
                print('--- Merge %s mutations with probabilities for %s' % (ID, id_probabilities))
                prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList,inputDir,outputDir,jobname,'ID',id_probabilities,startSimNum,endSimNum,'ID',PCAWG)

            #DBS
            if ((DBS in mutation_types_contexts) and (dbs_probabilities is not None)):
                print('--- Merge %s mutations with probabilities for %s' % (DBS,dbs_probabilities))
                prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList,inputDir,outputDir,jobname,'DBS',dbs_probabilities,startSimNum,endSimNum,'DBS',PCAWG)

            print("--- Merge simulations chr based files with Mutation Probabilities: %s seconds" %(time.time()-start_time))
            print("--- Merge simulations chr based files with Mutation Probabilities: %f minutes" %(float((time.time()-start_time)/60)))
            print('--- Merge simulations chr based files with Mutation Probabilities ends')
            print('#################################################################################\n')
            ####################################################################################################################
            ##################  Merge simulations chr based files with Mutation Probabilities ends #############################
            ####################################################################################################################

    #######################################################################################################
    #####################################  Full mode ends ###############################################
    #######################################################################################################


    #######################################################################################################
    #######################################################################################################
    #######################################################################################################
    if (necessary_dictionaries_already_exists):
        print('necessary_dictionaries_already_exists:%s' %(necessary_dictionaries_already_exists))
        subsSignature2PropertiesListDict=readDictionary(os.path.join(outputDir,jobname,DATA,SubsSignature2PropertiesListDictFilename))
        indelsSignature2PropertiesListDict= readDictionary(os.path.join(outputDir,jobname,DATA,IndelsSignature2PropertiesListDictFilename))
        dinucsSignature2PropertiesListDict=readDictionary(os.path.join(outputDir,jobname,DATA,DinucsSignature2PropertiesListDictFilename))
    else:
        #################################################################################
        print('#################################################################################')
        print('--- Fill dictionaries using original data starts')
        ##################################################################################
        # For each signature we will find a cutoff value for mutations with average probability >=0.9
        # Our aim is to have at most 10% false positive rate in mutations
        # number of mutations >= 5K for subs signatures
        # number of mutations >= 1K for indels signatures
        # number of mutations >= 200 for dinuc signatures
        # If we can not satisfy this condition we will discard the signature

        cutoffs = []
        for cufoff in np.arange(0.5, 0.91, 0.01):
            cutoffs.append("%.2f" % (cufoff))

        # We are reading original data to fill the signature2PropertiesListDict
        subsSignature2PropertiesListDict = fillCutoff2Signature2PropertiesListDictionary(outputDir, jobname,
                                                                                         chromNamesList, SUBS, cutoffs,
                                                                                         average_probability,
                                                                                         num_of_sbs_required,
                                                                                         num_of_id_required,
                                                                                         num_of_dbs_required)
        indelsSignature2PropertiesListDict = fillCutoff2Signature2PropertiesListDictionary(outputDir, jobname,
                                                                                           chromNamesList, INDELS,
                                                                                           cutoffs, average_probability,
                                                                                           num_of_sbs_required,
                                                                                           num_of_id_required,
                                                                                           num_of_dbs_required)
        dinucsSignature2PropertiesListDict = fillCutoff2Signature2PropertiesListDictionary(outputDir, jobname,
                                                                                           chromNamesList, DINUCS,
                                                                                           cutoffs, average_probability,
                                                                                           num_of_sbs_required,
                                                                                           num_of_id_required,
                                                                                           num_of_dbs_required)
        ##################################################################################

        ##################################################################################
        # We are reading original data again to fill the mutationType based, sample based and signature based dictionaries
        if sample_based:
            # Create files
            # createFiles(outputDir, jobname, MutationType2NumberofMutatiosDictFilename)

            # Initialize
            mutationType2NumberofMutationsDict = {}

            # Using original data
            fill_mutations_dictionaries_write(outputDir, jobname, chromNamesList, SUBS,
                                              mutationType2NumberofMutationsDict,
                                              subsSignature2PropertiesListDict, num_of_sbs_required, num_of_id_required,
                                              num_of_dbs_required)
            fill_mutations_dictionaries_write(outputDir, jobname, chromNamesList, INDELS,
                                              mutationType2NumberofMutationsDict,
                                              indelsSignature2PropertiesListDict, num_of_sbs_required,
                                              num_of_id_required,
                                              num_of_dbs_required)
            fill_mutations_dictionaries_write(outputDir, jobname, chromNamesList, DINUCS,
                                              mutationType2NumberofMutationsDict,
                                              dinucsSignature2PropertiesListDict, num_of_sbs_required,
                                              num_of_id_required,
                                              num_of_dbs_required)

            # We are writing number of mutations for each mutation type.
            # e.g.: {"SUBS": 3982196, "INDELS": 234731}
            appendDictionaryUnderDataDirectory(mutationType2NumberofMutationsDict, outputDir, jobname,MutationType2NumberofMutatiosDictFilename)
        ##################################################################################
        print('--- Fill dictionaries using original data ends')
        print('#################################################################################\n')
        #################################################################################

    #######################################################################################################
    #######################################################################################################
    #######################################################################################################


    #################################################################################
    ################## Set full path library files starts ###########################

    ###############################################
    # We need full path of the library files
    default_epigenomics_files = [DEFAULT_HISTONE_OCCUPANCY_FILE1, DEFAULT_HISTONE_OCCUPANCY_FILE2,
                                 DEFAULT_HISTONE_OCCUPANCY_FILE3, DEFAULT_HISTONE_OCCUPANCY_FILE4,
                                 DEFAULT_HISTONE_OCCUPANCY_FILE5, DEFAULT_HISTONE_OCCUPANCY_FILE6]
    if (set(epigenomics_files) == set(default_epigenomics_files)):
        # Set default
        # Make the order correctly
        epigenomics_files = default_epigenomics_files
        epigenomics_files_memos = ['H3K27me3_Breast_Epithelium', 'H3K36me3_Breast_Epithelium',
                                   'H3K9me3_Breast_Epithelium', 'H3K27ac_Breast_Epithelium',
                                   'H3K4me1_Breast_Epithelium', 'H3K4me3_Breast_Epithelium']
        epigenomics_biosamples = ['Breast_Epithelium']
        for file_index, filename in enumerate(epigenomics_files):
            epigenomics_files[file_index] = os.path.join(current_abs_path, LIB, EPIGENOMICS, filename)
    ###############################################

    ###############################################
    # We need set nucleosome_file
    # By default nucleosome_biosample=K562 and nucleosome_file is None
    # Here we set filename with extension not full path
    #There can be 2 cases:
    #Case 1 : nucleosome_biosample is an available nucleosome biosample and nucleosome_file is set as filename without fullpath
    #Case 2 : nucleosome_biosample is not an available nucleosome biosample and nucleosome_file is already set as filename with fullpath by the user

    #Case1: nucleosome_biosample is not None, nucleosome_file is a filename without fullpath
    if (nucleosome_biosample in available_nucleosome_biosamples):
        #Sets the filename without the full path
        nucleosome_file = getNucleosomeFile(nucleosome_biosample)
    #Case2: nucleosome_biosample is set to None, nucleosome_file is a filename with fullpath
    else:
        # We expect that user has provided nucleosome file with full path
        nucleosome_biosample = None
    ###############################################

    ###############################################
    #We need full path of the library files
    #By default replication_time_biosample=MCF7 and signal, valley, peak files are None
    if (replication_time_signal_file is None) and (replication_time_valley_file is None) and (replication_time_peak_file is None):
        replication_time_signal_file, replication_time_valley_file,replication_time_peak_file=getReplicationTimeFiles(replication_time_biosample)

    #User has provided replication time files
    #User must provide replication files with full paths
    elif (replication_time_signal_file is not None) and (replication_time_valley_file is not None) and (replication_time_peak_file is not None):
        replication_time_biosample=None
    ###############################################

    ################## Set full path library files ends #############################
    #################################################################################

    ####################################################################################################################
    ################################### Run SigProfilerTopography Analysis starts ######################################
    ####################################################################################################################
    print('#################################################################################')
    print('--- Run SigProfilerTopography Analysis starts')
    if (epigenomics):
        #Epigenomics
        occupancy_type=EPIGENOMICSOCCUPANCY
        deleteOldData(outputDir,jobname,occupancy_type)

        #data files are named using user provided epigenomics_files_memos or using epigenomics_file_memos_created
        epigenomics_file_memos_created=[]

        #Run for each epigenomics file
        for idx, epigenomics_file in enumerate(epigenomics_files):
            start_time = time.time()
            if (epigenomics_files_memos is not None) and (len(epigenomics_files_memos)==len(epigenomics_files)):
                epigenomics_file_memo= epigenomics_files_memos[idx]
            else:
                epigenomics_file_memo = os.path.splitext(os.path.basename(epigenomics_file))[0]
                epigenomics_file_memos_created.append(epigenomics_file_memo)

            runOccupancyAnalyses(genome,outputDir,jobname,numofSimulations,sample_based,epigenomics_file,epigenomics_file_memo,chromSizesDict,chromNamesList,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,computation_type,occupancy_type,plusorMinus_epigenomics,verbose)
            print('#################################################################################')
            print("--- Run Epigenomics Analyses: %s seconds --- %s" %((time.time()-start_time),epigenomics_file))
            print("--- Run Epigenomics Analyses: %f minutes --- %s" %(float((time.time()-start_time)/60),epigenomics_file))
            print('#################################################################################\n')

        #Used for plotting
        if (epigenomics_files_memos is None) or (len(epigenomics_files_memos) != len(epigenomics_files)):
            epigenomics_files_memos=epigenomics_file_memos_created

        if (epigenomics_biosamples is None) or (len(epigenomics_biosamples)==0):
            epigenomics_biosamples=[BIOSAMPLE_UNDECLARED]

    if (nucleosome):
        #Nucleosome Occupancy
        occupancy_type = NUCLEOSOMEOCCUPANCY
        deleteOldData(outputDir,jobname,occupancy_type)

        #For using SigProfilerTopography Nucleosme Files
        if nucleosome_biosample in available_nucleosome_biosamples:
            check_download_chrbased_npy_nuclesome_files(nucleosome_file,chromNamesList)

        start_time = time.time()
        runOccupancyAnalyses(genome,outputDir,jobname,numofSimulations,sample_based,nucleosome_file,None,chromSizesDict,chromNamesList,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,computation_type,occupancy_type,plusorMinus_nucleosome,verbose)
        print('#################################################################################')
        print("--- Run Nucleosome Occupancy Analyses: %s seconds --- %s" %((time.time()-start_time),nucleosome_file))
        print("--- Run Nucleosome Occupancy Analyses: %f minutes --- %s" %(float((time.time()-start_time)/60),nucleosome_file))
        print('#################################################################################\n')

    if (replication_time):
        # Replication Time
        # Required genome is already downloaded by matrix generator

        #For using SigProfilerTopography Provided Replication Time Files
        if replication_time_biosample in available_replication_time_biosamples:
            check_download_replication_time_files(replication_time_signal_file,replication_time_valley_file,replication_time_peak_file)

        start_time = time.time()
        runReplicationTimeAnalysis(genome,outputDir,jobname,numofSimulations,sample_based,replication_time_signal_file,chromSizesDict,chromNamesList,computation_type,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,verbose,matrix_generator_path)
        print('#################################################################################')
        print("--- Run Replication Time Analyses: %s seconds --- %s" %((time.time()-start_time),computation_type))
        print("--- Run Replication Time Analyses: %f minutes --- %s" %(float((time.time()-start_time)/60),computation_type))
        print('#################################################################################\n')

    if (strand_bias):

        # Replication Strand Bias
        start_time = time.time()
        runReplicationStrandBiasAnalysis(outputDir,jobname,numofSimulations,sample_based,replication_time_signal_file,replication_time_valley_file,replication_time_peak_file,chromSizesDict,chromNamesList,computation_type,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,verbose)
        print('#################################################################################')
        print("--- Run Replication Strand Bias Analyses: %s seconds --- %s" %((time.time()-start_time),computation_type))
        print("--- Run Replication Strand Bias Analyses: %f minutes --- %s" %(float((time.time()-start_time)/60),computation_type))
        print('#################################################################################\n')

        # Transcription Strand Bias
        start_time = time.time()
        runTranscriptionStradBiasAnalysis(outputDir,jobname,numofSimulations,sample_based,chromSizesDict,chromNamesList,computation_type,subsSignature2PropertiesListDict,indelsSignature2PropertiesListDict,dinucsSignature2PropertiesListDict,verbose)
        print('#################################################################################')
        print("--- Run Transcription Strand Bias Analyses: %s seconds --- %s" %((time.time()-start_time),computation_type))
        print("--- Run Transcription Strand Bias Analyses: %f minutes --- %s" %(float((time.time()-start_time)/60),computation_type))
        print('#################################################################################\n')

    if (processivity):
        # Processivity
        start_time = time.time()
        #TODO shall we consider only the signatures in subsSignature2PropertiesListDict?
        runProcessivityAnalysis(mutation_types_contexts,outputDir,jobname,numofSimulations,chromNamesList,computation_type,subsSignature2PropertiesListDict,verbose)
        print('#################################################################################')
        print("--- Run Processivity Analyses: %s seconds ---" %(time.time()-start_time))
        print("--- Run Processivity Analyses: %f minutes ---" %(float((time.time()-start_time)/60)))
        print('#################################################################################\n')

    print('#################################################################################\n')
    ####################################################################################################################
    ################################### Run SigProfilerTopography Analysis ends ########################################
    ####################################################################################################################

    ####################################################################################################################
    ############################################ Plot figures starts ###################################################
    ####################################################################################################################
    if (plot_figures):
        print('#################################################################################')
        print('--- Plot figures starts')
        start_time = time.time()
        # plotFigures(outputDir, jobname, numofSimulations, 'BONFERRONI_CORRECTION', 'USING_ONE_SAMPLE_TTEST')
        # plotFigures(outputDir, jobname, numofSimulations, 'BONFERRONI_CORRECTION', 'USING_NULL_DISTRIBUTION')
        plotFigures(outputDir, jobname, numofSimulations, sample_based,'USING_ZSCORE', 'USING_ZSCORE',mutation_types_contexts,epigenomics_files,epigenomics_files_memos,epigenomics_biosamples,nucleosome_file,epigenomics,nucleosome,replication_time,strand_bias,processivity,plusorMinus_epigenomics,plusorMinus_nucleosome,verbose,DATA_IS_READY_PLOT_FIGURES)
        print('#################################################################################')
        print("--- Plot Figures: %s seconds ---" %(time.time()-start_time))
        print("--- Plot Figures: %f minutes ---" %(float((time.time()-start_time)/60)))
        print('--- Plot figures ends')
        print('#################################################################################\n')
    ####################################################################################################################
    ############################################ Plot figures ends #####################################################
    ####################################################################################################################

    print('#################################################################################')
    print("--- SigProfilerTopography ended successfully")
    print("--- Thanks for using SigProfilerTopography")
    print('#################################################################################\n')


#######################################################



##############################################################
# BONFERRONI_CORRECTION = 'BONFERRONI_CORRECTION'
# FDR_BH_CORRECTION = 'FDR_BH_CORRECTION'
# USING_POISSON_DISTRIBUTION = 'USING_POISSON_DISTRIBUTION'
# USING_NULL_DISTRIBUTION = 'USING_NULL_DISTRIBUTION'
# USING_GAUSSIAN_KDE = 'USING_GAUSSIAN_KDE'
#Plot Figures for the attainded data after SigProfilerTopography Analyses
def plotFigures(outputDir,jobname,numberofSimulations,sample_based,multipleTesting,probabilityCalculation,mutationTypes,epigenomics_files,epigenomics_files_memos,epigenomics_biosamples,nucleosome_file,epigenomics,nucleosome,replication_time,strand_bias,processivity,plusOrMinus_epigenomics,plusOrMinus_nucleosome,verbose,DATA_IS_READY_PLOT_FIGURES):

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
    if (epigenomics or DATA_IS_READY_PLOT_FIGURES):
        occupancy_type=EPIGENOMICSOCCUPANCY
        deleteOldFigures(outputDir, jobname, occupancy_type)

        #Please note that epigenomics_file_memo is not None
        #If None it is created.
        for idx, epigenomics_file in enumerate(epigenomics_files):
            epigenomics_file_basename = os.path.basename(epigenomics_file)
            epigenomics_file_memo= epigenomics_files_memos[idx]
            occupancyAverageSignalFigures(outputDir, jobname, figureAugmentation, numberofSimulations,sample_based, mutationTypes,epigenomics_file_basename,epigenomics_file_memo,occupancy_type,plusOrMinus_epigenomics,verbose)

        plot_heatmaps(outputDir,jobname,numberofSimulations,epigenomics_files_memos,epigenomics_biosamples,occupancy_type)

    if (nucleosome or DATA_IS_READY_PLOT_FIGURES):
        occupancy_type=NUCLEOSOMEOCCUPANCY
        deleteOldFigures(outputDir, jobname, occupancy_type)
        nucleosome_file_basename = os.path.basename(nucleosome_file)
        occupancyAverageSignalFigures(outputDir,jobname,figureAugmentation,numberofSimulations,sample_based,mutationTypes,nucleosome_file_basename,None,occupancy_type,plusOrMinus_nucleosome,verbose)
    if (replication_time or DATA_IS_READY_PLOT_FIGURES):
        replicationTimeNormalizedMutationDensityFigures(outputDir,jobname,figureAugmentation,numberofSimulations,sample_based,mutationTypes)
    if (strand_bias or DATA_IS_READY_PLOT_FIGURES):
        transcriptionReplicationStrandBiasFigures(outputDir,jobname,figureAugmentation,numberofSimulations,sample_based)
    if (processivity or DATA_IS_READY_PLOT_FIGURES):
        processivityFigures(outputDir,jobname,numberofSimulations,multipleTesting,probabilityCalculation)
    ############################################################

##############################################################



##############################################################
#To run on laptob
import os

if __name__== "__main__":

    genome = 'GRCh37'
    jobname = 'Test-Skin-Melanoma'
    numberofSimulations = 2

    inputDir = '/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/input/PCAWG_Matlab_Clean/Skin-Melanoma/filtered/'
    outputDir = os.path.join('C:\\','Users','burcak','Developer','Python','SigProfilerTopography','SigProfilerTopography','output_test')

    sbs_probabilities_file_path = os.path.join('C:\\','Users','burcak','Developer','Python','SigProfilerTopography','SigProfilerTopography','output_for_extractor','PCAWG_Matlab','Skin-Melanoma_sbs96_mutation_probabilities.txt')
    id_probabilities_file_path = os.path.join('C:\\','Users','burcak','Developer','Python','SigProfilerTopography','SigProfilerTopography','output_for_extractor','PCAWG_Matlab','Skin-Melanoma_id83_mutation_probabilities.txt')
    dbs_probabilities_file_path = os.path.join('C:\\','Users','burcak','Developer','Python','SigProfilerTopography','SigProfilerTopography','output_for_extractor','PCAWG_Matlab','Skin-Melanoma_dbs_mutation_probabilities.txt')

    # user_provided_replication_time_file_path = os.path.join('C:\\','Users','burcak','Developer','Python','SigProfilerTopography','SigProfilerTopography','lib','replication','wgEncodeUwRepliSeqNhekWaveSignalRep1.wig')
    # user_provided_replication_time_valley_file_path = os.path.join('C:\\','Users','burcak','Developer','Python','SigProfilerTopography','SigProfilerTopography','lib','replication','wgEncodeUwRepliSeqNhekValleysRep1.bed')
    # user_provided_replication_time_peak_file_path = os.path.join('C:\\','Users','burcak','Developer','Python','SigProfilerTopography','SigProfilerTopography','lib','replication','wgEncodeUwRepliSeqNhekPkRep1.bed')

    # user_provided_nucleosome_file_path= os.path.join('C:\\','Users','burcak','Developer','Python','SigProfilerTopography','SigProfilerTopography','lib','nucleosome','wgEncodeSydhNsomeK562Sig.wig')
    user_provided_nucleosome_file_path = os.path.join('C:\\', 'Users', 'burcak', 'Developer', 'Python','SigProfilerTopography', 'SigProfilerTopography', 'lib','nucleosome', 'wgEncodeSydhNsomeGm12878Sig.wig')
    # user_provided_nucleosome_file_path= os.path.join('C:\\','Users','burcak','Developer','Python','SigProfilerTopography','SigProfilerTopography','lib','nucleosome','wgEncodeSydhNsomeGm12878Sig.bigWig')

    runAnalyses(genome, inputDir, outputDir, jobname, numberofSimulations,
                           sbs_probabilities=sbs_probabilities_file_path,
                           id_probabilities=id_probabilities_file_path,
                           dbs_probabilities=dbs_probabilities_file_path,
                            # nucleosome_biosample='K562',
                            # replication_time_biosample='NHEK',
                           # nucleosome_file=user_provided_nucleosome_file_path,
                           # replication_time_signal_file=user_provided_replication_time_file_path,
                           # replication_time_valley_file=user_provided_replication_time_valley_file_path,
                           # replication_time_peak_file=user_provided_replication_time_peak_file_path,
                           epigenomics=True, nucleosome=False, replication_time=False, strand_bias=False, processivity=False,
                           sample_based=False, new_simulations_enforced=False, full_mode=False, verbose=False,necessary_dictionaries_already_exists=True)
##############################################################