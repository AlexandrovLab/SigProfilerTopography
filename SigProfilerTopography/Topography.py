# !/usr/bin/env python3

# Author: burcakotlu

# Contact: burcakotlu@eng.ucsd.edu

# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018-2020 Burcak Otlu

import os
import sys
import math
import time
import numpy as np
import pandas as pd
import scipy
import statsmodels
import matplotlib as plt
import datetime

import shutil
import platform
import multiprocessing

import SigProfilerMatrixGenerator as matrix_generator

from SigProfilerMatrixGenerator import version as matrix_generator_version
from SigProfilerSimulator import version as simulator_version

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerSimulator import SigProfilerSimulator as simulator

from SigProfilerTopography import version as topography_version
from SigProfilerTopography.source.commons.TopographyCommons import readProbabilities
from SigProfilerTopography.source.commons.TopographyCommons import readChrBasedMutationsMergeWithProbabilitiesAndWrite

from SigProfilerTopography.source.commons.TopographyCommons import DATA
from SigProfilerTopography.source.commons.TopographyCommons import FIGURE
from SigProfilerTopography.source.commons.TopographyCommons import SAMPLE

from SigProfilerTopography.source.commons.TopographyCommons import K562
from SigProfilerTopography.source.commons.TopographyCommons import MCF7
from SigProfilerTopography.source.commons.TopographyCommons import IMR90

from SigProfilerTopography.source.commons.TopographyCommons import MEF
from SigProfilerTopography.source.commons.TopographyCommons import ESC

from SigProfilerTopography.source.commons.TopographyCommons import MM10
from SigProfilerTopography.source.commons.TopographyCommons import GRCh37
from SigProfilerTopography.source.commons.TopographyCommons import GRCh38

from SigProfilerTopography.source.commons.TopographyCommons import SIGPROFILERTOPOGRAPHY_DEFAULT_FILES

from SigProfilerTopography.source.commons.TopographyCommons import getNucleosomeFile
from SigProfilerTopography.source.commons.TopographyCommons import getReplicationTimeFiles

from SigProfilerTopography.source.commons.TopographyCommons import available_nucleosome_biosamples
from SigProfilerTopography.source.commons.TopographyCommons import available_replication_time_biosamples

from SigProfilerTopography.source.commons.TopographyCommons import EPIGENOMICSOCCUPANCY
from SigProfilerTopography.source.commons.TopographyCommons import NUCLEOSOMEOCCUPANCY
from SigProfilerTopography.source.commons.TopographyCommons import REPLICATIONTIME
from SigProfilerTopography.source.commons.TopographyCommons import REPLICATIONSTRANDBIAS
from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIPTIONSTRANDBIAS
from SigProfilerTopography.source.commons.TopographyCommons import PROCESSIVITY

from SigProfilerTopography.source.commons.TopographyCommons import EPIGENOMICS
from SigProfilerTopography.source.commons.TopographyCommons import STRANDBIAS

from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_H3K27ME3_OCCUPANCY_FILE
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_H3K36ME3_OCCUPANCY_FILE
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_H3K9ME3_OCCUPANCY_FILE
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_H3K27AC_OCCUPANCY_FILE
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_H3K4ME1_OCCUPANCY_FILE
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_H3K4ME3_OCCUPANCY_FILE
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_CTCF_OCCUPANCY_FILE
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_ATAC_SEQ_OCCUPANCY_FILE

from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_ATAC_SEQ_GRCh38_OCCUPANCY_FILE
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_H3K27ME3_GRCh38_OCCUPANCY_FILE
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_H3K36ME3_GRCh38_OCCUPANCY_FILE
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_H3K9ME3_GRCh38_OCCUPANCY_FILE
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_H3K27AC_GRCh38_OCCUPANCY_FILE
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_H3K4ME1_GRCh38_OCCUPANCY_FILE
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_H3K4ME3_GRCh38_OCCUPANCY_FILE
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_CTCF_GRCh38_OCCUPANCY_FILE

from SigProfilerTopography.source.commons.TopographyCommons import MM10_mmNuc0020101_GSM1004653_ESC_NUCLEOSOME_FILE
from SigProfilerTopography.source.commons.TopographyCommons import MM10_MEF_NUCLEOSOME_FILE

from SigProfilerTopography.source.commons.TopographyCommons import GM12878_NUCLEOSOME_OCCUPANCY_FILE

from SigProfilerTopography.source.commons.TopographyCommons import K562_NUCLEOSOME_OCCUPANCY_FILE
from SigProfilerTopography.source.commons.TopographyCommons import K562_GRCh38_NUCLEOSOME_OCCUPANCY_FILE

from SigProfilerTopography.source.commons.TopographyCommons import ENCFF575PMI_mm10_embryonic_facial_prominence_ATAC_seq
from SigProfilerTopography.source.commons.TopographyCommons import ENCFF993SRY_mm10_embryonic_fibroblast_H3K4me1
from SigProfilerTopography.source.commons.TopographyCommons import ENCFF912DNP_mm10_embryonic_fibroblast_H3K4me3
from SigProfilerTopography.source.commons.TopographyCommons import ENCFF611HDQ_mm10_embryonic_fibroblast_CTCF
from SigProfilerTopography.source.commons.TopographyCommons import ENCFF152DUV_mm10_embryonic_fibroblast_POLR2A
from SigProfilerTopography.source.commons.TopographyCommons import ENCFF114VLZ_mm10_embryonic_fibroblast_H3K27ac

from SigProfilerTopography.source.commons.TopographyCommons import SBS
from SigProfilerTopography.source.commons.TopographyCommons import DBS
from SigProfilerTopography.source.commons.TopographyCommons import ID

from SigProfilerTopography.source.commons.TopographyCommons import UNDECLARED

from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC
from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM
from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT
from SigProfilerTopography.source.commons.TopographyCommons import STRINGENT

from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_AVERAGE_PROBABILITY
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_CUTOFF
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_NUM_OF_SBS_REQUIRED
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_NUM_OF_DBS_REQUIRED
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_NUM_OF_ID_REQUIRED
from SigProfilerTopography.source.commons.TopographyCommons import DEFAULT_NUM_OF_REAL_DATA_OVERLAP_REQUIRED

from SigProfilerTopography.source.commons.TopographyCommons import CONSIDER_COUNT
from SigProfilerTopography.source.commons.TopographyCommons import CONSIDER_DISTANCE
from SigProfilerTopography.source.commons.TopographyCommons import CONSIDER_DISTANCE_ALL_SAMPLES_TOGETHER

from SigProfilerTopography.source.commons.TopographyCommons import MISSING_SIGNAL
from SigProfilerTopography.source.commons.TopographyCommons import NO_SIGNAL

from SigProfilerTopography.source.commons.TopographyCommons import SBS96
from SigProfilerTopography.source.commons.TopographyCommons import ID
from SigProfilerTopography.source.commons.TopographyCommons import DBS

from SigProfilerTopography.source.commons.TopographyCommons import SUBS
from SigProfilerTopography.source.commons.TopographyCommons import INDELS
from SigProfilerTopography.source.commons.TopographyCommons import DINUCS

from SigProfilerTopography.source.commons.TopographyCommons import SBS_CONTEXTS
from SigProfilerTopography.source.commons.TopographyCommons import SNV

from SigProfilerTopography.source.commons.TopographyCommons import CHRBASED
from SigProfilerTopography.source.commons.TopographyCommons import LIB

from SigProfilerTopography.source.commons.TopographyCommons import getChromSizesDict
from SigProfilerTopography.source.commons.TopographyCommons import getShortNames
from SigProfilerTopography.source.commons.TopographyCommons import copyMafFiles
from SigProfilerTopography.source.commons.TopographyCommons import fillCutoff2Signature2PropertiesListDictionary
from SigProfilerTopography.source.commons.TopographyCommons import fill_signature_number_of_mutations_df
from SigProfilerTopography.source.commons.TopographyCommons import fill_mutations_dictionaries_write
from SigProfilerTopography.source.commons.TopographyCommons import get_mutation_type_context_for_probabilities_file

from SigProfilerTopography.source.commons.TopographyCommons import Table_MutationType_NumberofMutations_NumberofSamples_SamplesList_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Table_ChrLong_NumberofMutations_Filename

from SigProfilerTopography.source.commons.TopographyCommons import Table_SBS_Signature_Cutoff_NumberofMutations_AverageProbability_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Table_DBS_Signature_Cutoff_NumberofMutations_AverageProbability_Filename
from SigProfilerTopography.source.commons.TopographyCommons import Table_ID_Signature_Cutoff_NumberofMutations_AverageProbability_Filename

from SigProfilerTopography.source.commons.TopographyCommons import NUMBER_OF_MUTATIONS_IN_EACH_SPLIT

from SigProfilerTopography.source.occupancy.OccupancyAnalysis import occupancyAnalysis
from SigProfilerTopography.source.replicationtime.ReplicationTimeAnalysis import replicationTimeAnalysis
from SigProfilerTopography.source.replicationstrandbias.ReplicationStrandBiasAnalysis import replicationStrandBiasAnalysis
from SigProfilerTopography.source.transcriptionstrandbias.TranscriptionStrandBiasAnalysis import transcriptionStrandBiasAnalysis
from SigProfilerTopography.source.processivity.ProcessivityAnalysis import processivityAnalysis

from SigProfilerTopography.source.annotation.Mutation_Annotation_Integration import mutation_annotation_replication_timing_integration
from SigProfilerTopography.source.annotation.Mutation_Annotation_Integration import mutation_annotation_replication_timing_integration_signature_specific

from SigProfilerTopography.source.plotting.OccupancyAverageSignalFigures import occupancyAverageSignalFigures
from SigProfilerTopography.source.plotting.OccupancyAverageSignalFigures import compute_fold_change_with_p_values_plot_heatmaps
from SigProfilerTopography.source.plotting.ReplicationTimeNormalizedMutationDensityFigures import replicationTimeNormalizedMutationDensityFigures
from SigProfilerTopography.source.plotting.TranscriptionReplicationStrandBiasFigures import transcriptionReplicationStrandBiasFiguresUsingDataframes
from SigProfilerTopography.source.plotting.ProcessivityFigures import processivityFigures

from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIBED_VERSUS_UNTRANSCRIBED
from SigProfilerTopography.source.commons.TopographyCommons import GENIC_VERSUS_INTERGENIC
from SigProfilerTopography.source.commons.TopographyCommons import LAGGING_VERSUS_LEADING
from SigProfilerTopography.source.commons.TopographyCommons import PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_TOOL

from SigProfilerTopography.source.commons.TopographyCommons import COMBINE_P_VALUES_METHOD_FISHER
from SigProfilerTopography.source.commons.TopographyCommons import WEIGHTED_AVERAGE_METHOD
from SigProfilerTopography.source.commons.TopographyCommons import COLORBAR_SEISMIC

from SigProfilerTopography.source.commons.TopographyCommons import natural_key

MATRIX_GENERATOR_PATH = matrix_generator.__path__[0]

# Can be moved to DataPreparationCommons under /source/commons
# Read chr based mutations (provided by SigProfilerMatrixGenerator) and merge with probabilities (provided by SigProfilerTopography)
def prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList,
                                                                       inputDir,
                                                                       outputDir,
                                                                       jobname,
                                                                       mutation_type_context,
                                                                       mutations_probabilities_file_path,
                                                                       mutation_type_context_for_probabilities,
                                                                       startSimNum,
                                                                       endSimNum,
                                                                       partialDirname,
                                                                       PCAWG,
                                                                       log_file,
                                                                       verbose):

    #original matrix generator chrbased data will be under inputDir/output/vcf_files/SNV
    #original matrix generator chrbased data will be under inputDir/output/vcf_files/DBS
    #original matrix generator chrbased data will be under inputDir/output/vcf_files/ID

    #sim1 matrix generator chrbased data will be under inputDir/output/simulations/sim1/96/output/vcf_files/SNV
    #sim1 matrix generator chrbased data will be under inputDir/output/simulations/sim1/ID/output/vcf_files/ID
    #sim1 matrix generator chrbased data will be under inputDir/output/simulations/sim1/DBS/output/vcf_files/DBS

    df_columns_contain_ordered_signatures = None

    os.makedirs(os.path.join(outputDir, jobname, DATA, CHRBASED), exist_ok=True)
    for simNum in range(1,endSimNum+1):
        simName = 'sim%d' % (simNum)
        os.makedirs(os.path.join(outputDir, jobname, DATA, CHRBASED, simName), exist_ok=True)

    if ((mutations_probabilities_file_path is not None) and (os.path.exists(mutations_probabilities_file_path))):
        mutations_probabilities_df = readProbabilities(mutations_probabilities_file_path, log_file, verbose)
        df_columns_contain_ordered_signatures = mutations_probabilities_df.columns.values

        if verbose:
            log_out = open(log_file, 'a')
            print('\tVerbose mutations_probabilities_df.head()', file=log_out)
            print('\tVerbose %s' %(mutations_probabilities_df.head()), file=log_out)
            print('\tVerbose mutations_probabilities_df.columns.values', file=log_out)
            print('\tVerbose %s' %(mutations_probabilities_df.columns.values), file=log_out)
            log_out.close()

        # Step1 SigProfilerTopography Python Package
        # For Release we will use SAMPLE as it is, no change in SAMPLE column is needed.

        # For PCAWG_Matlab
        # This statement below is customized for  PCAWG_Matlab
        # To get rid of inconsistent cancer type names in sample column of chrbased mutation files and probabilities files
        # Breast-LobularCA_SP124191
        if PCAWG:
            mutations_probabilities_df[SAMPLE] = mutations_probabilities_df[SAMPLE].str.split('_',expand=True)[1]

        ############################################################################################
        ##############################  pool.apply_async starts ####################################
        ############################################################################################

        numofProcesses = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=numofProcesses)

        jobs = []

        sim_nums = range(startSimNum,endSimNum+1)
        sim_num_chr_tuples = ((sim_num, chrShort) for sim_num in sim_nums for chrShort in chromShortNamesList)

        for simNum, chrShort in sim_num_chr_tuples:
            simName = 'sim%d' % (simNum)
            chr_based_mutation_filename = '%s_seqinfo.txt' % (chrShort)

            if (simNum == 0):
                matrix_generator_output_dir_path = os.path.join(inputDir, 'output', 'vcf_files', partialDirname)
                os.makedirs(matrix_generator_output_dir_path, exist_ok=True)
            else:
                matrix_generator_output_dir_path = os.path.join(inputDir, 'output', 'simulations', simName, mutation_type_context, 'output', 'vcf_files', partialDirname)
                os.makedirs(matrix_generator_output_dir_path, exist_ok=True)

            if (os.path.exists(matrix_generator_output_dir_path)):
                chr_based_mutation_filepath = os.path.join(matrix_generator_output_dir_path, chr_based_mutation_filename)
                inputList = []
                inputList.append(chrShort)
                inputList.append(outputDir)
                inputList.append(jobname)
                inputList.append(chr_based_mutation_filepath)
                inputList.append(mutations_probabilities_df)
                inputList.append(mutation_type_context_for_probabilities)
                inputList.append(mutation_type_context)
                inputList.append(simNum)
                inputList.append(PCAWG)
                inputList.append(log_file)
                jobs.append(pool.apply_async(readChrBasedMutationsMergeWithProbabilitiesAndWrite,
                                             args=(inputList,)))

        # wait for all jobs to finish
        for job in jobs:
            if verbose:
                log_out = open(log_file, 'a')
                print('\tVerbose merge chrom based mutations with probabilities worker pid %s job.get():%s ' % (str(os.getpid()), job.get()), file=log_out)
                log_out.close()

        pool.close()
        pool.join()

        ############################################################################################
        ##############################  pool.apply_async ends ######################################
        ############################################################################################

    elif ((mutations_probabilities_file_path is None) or (not (os.path.exists(mutations_probabilities_file_path)))):
        # For Information
        log_out = open(log_file, 'a')
        print('--- There is a situation/problem: mutations_probabilities_file_path:%s is None or does not exist.' %(mutations_probabilities_file_path), file=log_out)
        log_out.close()

        ############################################################################################
        ##############################  pool.apply_async starts ####################################
        ############################################################################################

        numofProcesses = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=numofProcesses)

        jobs = []

        sim_nums = range(startSimNum,endSimNum+1)
        sim_num_chr_tuples = ((sim_num, chrShort) for sim_num in sim_nums for chrShort in chromShortNamesList)

        for simNum, chrShort in sim_num_chr_tuples:
            simName = 'sim%d' % (simNum)
            chr_based_mutation_filename = '%s_seqinfo.txt' % (chrShort)
            if (simNum == 0):
                matrix_generator_output_dir_path = os.path.join(inputDir, 'output', 'vcf_files', partialDirname)
            else:
                matrix_generator_output_dir_path = os.path.join(inputDir, 'output', 'simulations', simName,mutation_type_context, 'output', 'vcf_files',partialDirname)
            if (os.path.exists(matrix_generator_output_dir_path)):
                chr_based_mutation_filepath = os.path.join(matrix_generator_output_dir_path,chr_based_mutation_filename)
                inputList = []
                inputList.append(chrShort)
                inputList.append(outputDir)
                inputList.append(jobname)
                inputList.append(chr_based_mutation_filepath)
                inputList.append(None)
                inputList.append(mutation_type_context_for_probabilities)
                inputList.append(mutation_type_context)
                inputList.append(simNum)
                inputList.append(PCAWG)
                inputList.append(log_file)
                jobs.append(pool.apply_async(readChrBasedMutationsMergeWithProbabilitiesAndWrite,
                                             args=(inputList,)))

        # wait for all jobs to finish
        for job in jobs:
            if verbose:
                log_out = open(log_file, 'a')
                print('\tVerbose  merge chrom based mutations with probabilities worker pid %s job.get():%s ' % (str(os.getpid()), job.get()), file=log_out)
                log_out.close()

        pool.close()
        pool.join()

        ############################################################################################
        ##############################  pool.apply_async ends ######################################
        ############################################################################################

    return df_columns_contain_ordered_signatures


def check_download_replication_time_files(replication_time_signal_file, replication_time_valley_file, replication_time_peak_file):
    current_abs_path = os.path.dirname(os.path.abspath(__file__))

    # These are currently full path, therefore convert them to filename
    replication_time_signal_file = os.path.basename(replication_time_signal_file)
    replication_time_valley_file = os.path.basename(replication_time_valley_file)
    replication_time_peak_file = os.path.basename(replication_time_peak_file)

    os.makedirs(os.path.join(current_abs_path,'lib','replication'),exist_ok=True)
    lib_replication_path = os.path.join(current_abs_path,'lib','replication')

    if os.path.isabs(lib_replication_path):
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
        # It has to be an absolute path
        print('%s is not an absolute path.' %(lib_replication_path))

    # go back
    os.chdir(current_abs_path)


def check_download_sample_probability_files():
    current_path = os.getcwd()
    os.makedirs(os.path.join(current_path, 'sample_probabilities'), exist_ok=True)
    sample_probability_files_path = os.path.join(current_path, 'sample_probabilities')

    probability_files = ['COSMIC_DBS78_Decomposed_Mutation_Probabilities.txt',
                        'COSMIC_SBS96_Decomposed_Mutation_Probabilities.txt']

    if os.path.isabs(sample_probability_files_path):
        os.chdir(sample_probability_files_path)

        for probability_filename in probability_files:
            probability_file_path = os.path.join(sample_probability_files_path, probability_filename)
            if not os.path.exists(probability_file_path):
                print('Does not exists: %s' % (probability_file_path))
            try:
                print('Downloading %s under %s' % (probability_filename, sample_probability_files_path))

                # wget -c Continue getting a partially-downloaded file
                # wget -nc  If a file is downloaded more than once in the same directory, the local file will be clobbered, or overwritten
                # cmd="bash -c '" + 'wget -r -l1 -c -nc --no-parent -nd -P ' + chrombased_npy_path + ' ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerTopography/lib/nucleosome/chrbased/' + filename + "'"

                # -r When included, the wget will recursively traverse subdirectories in order to obtain all content.
                # -l1 Limit recursion depth to a specific number of levels, by setting the <#> variable to the desired number.
                # -c option to resume a download
                # -nc, --no-clobber If a file is downloaded more than once in the same directory, Wget's behavior depends on a few options, including -nc.  In certain cases, the local file will be clobbered, or overwritten, upon repeated download.  In other cases it will be preserved.
                # -np, --no-parent Do not ever ascend to the parent directory when retrieving recursively.  This is a useful option, since it guarantees that only the files below a certain hierarchy will be downloaded.
                # -nd, --no-directories When included, directories will not be created. All files captured in the wget will be copied directly in to the active directory
                cmd = "bash -c '" + 'wget -r -l1 -c -nc --no-parent -nd ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerTopography/sample_probability_files/' + probability_filename + "'"
                print("cmd: %s" % cmd)
                os.system(cmd)
            except:
                print("The UCSD ftp site is not responding...")
    else:
        # It has to be an absolute path
        print('%s is not an absolute path.' % (sample_probability_files_path))

    # go back
    os.chdir(current_path)



def check_download_sample_vcf_files():
    current_path = os.getcwd()
    os.makedirs(os.path.join(current_path, 'sample_vcfs'), exist_ok=True)
    sample_vcf_files_path = os.path.join(current_path, 'sample_vcfs')

    vcf_files = ['PD4248a.vcf', 'PD4199a.vcf', 'PD4198a.vcf', 'PD4194a.vcf', 'PD4192a.vcf', 'PD4120a.vcf',
                'PD4116a.vcf', 'PD4115a.vcf', 'PD4109a.vcf', 'PD4107a.vcf', 'PD4103a.vcf', 'PD4088a.vcf',
                'PD4086a.vcf', 'PD4085a.vcf', 'PD4006a.vcf', 'PD4005a.vcf', 'PD3945a.vcf', 'PD3905a.vcf',
                'PD3904a.vcf', 'PD3890a.vcf', 'PD3851a.vcf']

    if os.path.isabs(sample_vcf_files_path):
        os.chdir(sample_vcf_files_path)

        for vcf_filename in vcf_files:
            vcf_file_path = os.path.join(sample_vcf_files_path, vcf_filename)
            if not os.path.exists(vcf_file_path):
                print('Does not exists: %s' % (vcf_file_path))
            try:
                print('Downloading %s under %s' % (vcf_filename, sample_vcf_files_path))

                # wget -c Continue getting a partially-downloaded file
                # wget -nc  If a file is downloaded more than once in the same directory, the local file will be clobbered, or overwritten
                # cmd="bash -c '" + 'wget -r -l1 -c -nc --no-parent -nd -P ' + chrombased_npy_path + ' ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerTopography/lib/nucleosome/chrbased/' + filename + "'"

                # -r When included, the wget will recursively traverse subdirectories in order to obtain all content.
                # -l1 Limit recursion depth to a specific number of levels, by setting the <#> variable to the desired number.
                # -c option to resume a download
                # -nc, --no-clobber If a file is downloaded more than once in the same directory, Wget's behavior depends on a few options, including -nc.  In certain cases, the local file will be clobbered, or overwritten, upon repeated download.  In other cases it will be preserved.
                # -np, --no-parent Do not ever ascend to the parent directory when retrieving recursively.  This is a useful option, since it guarantees that only the files below a certain hierarchy will be downloaded.
                # -nd, --no-directories When included, directories will not be created. All files captured in the wget will be copied directly in to the active directory
                cmd = "bash -c '" + 'wget -r -l1 -c -nc --no-parent -nd ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerTopography/sample_vcf_files/' + vcf_filename + "'"
                print("cmd: %s" % cmd)
                os.system(cmd)
            except:
                print("The UCSD ftp site is not responding...")
    else:
        # It has to be an absolute path
        print('%s is not an absolute path.' % (sample_vcf_files_path))

    # go back
    os.chdir(current_path)


def check_download_chrbased_npy_atac_seq_files(atac_seq_file, chromNamesList):
    current_abs_path = os.path.dirname(os.path.abspath(__file__))

    os.makedirs(os.path.join(current_abs_path, 'lib', 'epigenomics', 'chrbased'), exist_ok=True)
    chrombased_npy_path = os.path.join(current_abs_path, 'lib', 'epigenomics', 'chrbased')

    if os.path.isabs(chrombased_npy_path):
        os.chdir(chrombased_npy_path)

        atac_seq_filename_wo_extension = os.path.splitext(os.path.basename(atac_seq_file))[0]

        for chrLong in chromNamesList:
            filename = '%s_signal_%s.npy' % (chrLong, atac_seq_filename_wo_extension)

            chrbased_npy_array_path = os.path.join(chrombased_npy_path, filename)
            if not os.path.exists(chrbased_npy_array_path):
                print('Does not exists: %s' % (chrbased_npy_array_path))
                try:
                    print('Downloading %s under %s' % (filename, chrbased_npy_array_path))

                    # wget -c Continue getting a partially-downloaded file
                    # wget -nc  If a file is downloaded more than once in the same directory, the local file will be clobbered, or overwritten
                    # cmd="bash -c '" + 'wget -r -l1 -c -nc --no-parent -nd -P ' + chrombased_npy_path + ' ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerTopography/lib/nucleosome/chrbased/' + filename + "'"

                    # -r When included, the wget will recursively traverse subdirectories in order to obtain all content.
                    # -l1 Limit recursion depth to a specific number of levels, by setting the <#> variable to the desired number.
                    # -c option to resume a download
                    # -nc, --no-clobber If a file is downloaded more than once in the same directory, Wget's behavior depends on a few options, including -nc.  In certain cases, the local file will be clobbered, or overwritten, upon repeated download.  In other cases it will be preserved.
                    # -np, --no-parent Do not ever ascend to the parent directory when retrieving recursively.  This is a useful option, since it guarantees that only the files below a certain hierarchy will be downloaded.
                    # -nd, --no-directories When included, directories will not be created. All files captured in the wget will be copied directly in to the active directory
                    cmd = "bash -c '" + 'wget -r -l1 -c -nc --no-parent -nd ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerTopography/lib/epigenomics/chrbased/' + filename + "'"
                    print("cmd: %s" %cmd)
                    os.system(cmd)
                except:
                    # print("The UCSD ftp site is not responding...pulling from sanger ftp now.")
                    print("The UCSD ftp site is not responding...")

    else:
        #It has to be an absolute path
        print('%s is not an absolute path.' %(chrombased_npy_path))

    #go back
    os.chdir(current_abs_path)

# Download nucleosome occupancy chr based npy files from ftp alexandrovlab if they do not exists
# We are using this function if user is using our available nucleosome data for GM12878 adnd K562 cell lines
def check_download_chrbased_npy_nuclesome_files(nucleosome_file, chromNamesList):
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
        # It has to be an absolute path
        print('%s is not an absolute path.' %(chrombased_npy_path))

    # go back to current directory
    os.chdir(current_abs_path)


def install_default_nucleosome(genome):
    chromSizesDict = getChromSizesDict(genome)
    chromNamesList = list(chromSizesDict.keys())

    if genome == MM10:
        #Case1: File is not set, Biosample is not set
        nucleosome_file = MM10_mmNuc0020101_GSM1004653_ESC_NUCLEOSOME_FILE
        check_download_chrbased_npy_nuclesome_files(nucleosome_file, chromNamesList)

    elif genome == GRCh37:
        # Case1: File is not set, Biosample is not set
        nucleosome_file = K562_NUCLEOSOME_OCCUPANCY_FILE
        check_download_chrbased_npy_nuclesome_files(nucleosome_file, chromNamesList)

    elif genome == GRCh38:
        # Case1: File is not set, Biosample is not set
        nucleosome_file = K562_GRCh38_NUCLEOSOME_OCCUPANCY_FILE
        check_download_chrbased_npy_nuclesome_files(nucleosome_file, chromNamesList)


def install_default_atac_seq(genome):
    chromSizesDict = getChromSizesDict(genome)
    chromNamesList = list(chromSizesDict.keys())

    if genome == GRCh37:
        atac_seq_file = DEFAULT_ATAC_SEQ_OCCUPANCY_FILE
        check_download_chrbased_npy_atac_seq_files(atac_seq_file, chromNamesList)

    elif genome == GRCh38:
        atac_seq_file = DEFAULT_ATAC_SEQ_GRCh38_OCCUPANCY_FILE
        check_download_chrbased_npy_atac_seq_files(atac_seq_file, chromNamesList)

    elif genome == MM10:
        atac_seq_file = ENCFF575PMI_mm10_embryonic_facial_prominence_ATAC_seq
        check_download_chrbased_npy_atac_seq_files(atac_seq_file, chromNamesList)

def install_sample_vcf_files():
    # Download to where the SigProfilerTopography is run
    check_download_sample_vcf_files()


def install_sample_probability_files():
    # Download to where the SigProfilerTopography is run
    check_download_sample_probability_files()


# For Skin-Melanoma USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT is better
# For others USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM is better
def runOccupancyAnalyses(genome,
                         outputDir,
                         jobname,
                         numofSimulations,
                         job_tuples,
                         sample_based,
                         library_file_with_path,
                         library_file_memo,
                         chromSizesDict,
                         chromNamesList,
                         ordered_sbs_signatures_with_cutoffs_array,
                         ordered_dbs_signatures_with_cutoffs_array,
                         ordered_id_signatures_with_cutoffs_array,
                         ordered_sbs_signatures_cutoffs,
                         ordered_dbs_signatures_cutoffs,
                         ordered_id_signatures_cutoffs,
                         computation_type,
                         occupancy_type,
                         occupancy_calculation_type,
                         plusorMinus,
                         remove_outliers,
                         quantileValue,
                         discreet_mode,
                         default_cutoff,
                         parallel_mode,
                         log_file,
                         verbose):

    # common
    # computation_type = USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM
    # computation_type = USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT
    occupancyAnalysis(genome,
                      computation_type,
                      occupancy_type,
                      occupancy_calculation_type,
                      sample_based,
                      plusorMinus,
                      chromSizesDict,
                      chromNamesList,
                      outputDir,
                      jobname,
                      numofSimulations,
                      job_tuples,
                      library_file_with_path,
                      library_file_memo,
                      ordered_sbs_signatures_with_cutoffs_array,
                      ordered_dbs_signatures_with_cutoffs_array,
                      ordered_id_signatures_with_cutoffs_array,
                      ordered_sbs_signatures_cutoffs,
                      ordered_dbs_signatures_cutoffs,
                      ordered_id_signatures_cutoffs,
                      remove_outliers,
                      quantileValue,
                      discreet_mode,
                      default_cutoff,
                      parallel_mode,
                      log_file,
                      verbose)

def runReplicationTimeAnalysis(genome,
                               outputDir,
                               jobname,
                               numofSimulations,
                               job_tuples,
                               sample_based,
                               replicationTimeFilename,
                               chromSizesDict,
                               chromNamesList,
                               computation_type,
                               ordered_sbs_signatures_with_cutoffs,
                               ordered_dbs_signatures_with_cutoffs,
                               ordered_id_signatures_with_cutoffs,
                               ordered_sbs_signatures_cutoffs,
                               ordered_dbs_signatures_cutoffs,
                               ordered_id_signatures_cutoffs,
                               discreet_mode,
                               default_cutoff,
                               parallel_mode,
                               log_file,
                               verbose,
                               matrix_generator_path):

    # Fill np array during runtime managed by replication_time_np_arrays_fill_runtime=True
    # Supported computation types
    # computation_type = USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM
    # computation_type = USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT
    replicationTimeAnalysis(computation_type,
                            sample_based,
                            genome,
                            chromSizesDict,
                            chromNamesList,
                            outputDir,
                            jobname,
                            numofSimulations,
                            job_tuples,
                            replicationTimeFilename,
                            ordered_sbs_signatures_with_cutoffs,
                            ordered_dbs_signatures_with_cutoffs,
                            ordered_id_signatures_with_cutoffs,
                            ordered_sbs_signatures_cutoffs,
                            ordered_dbs_signatures_cutoffs,
                            ordered_id_signatures_cutoffs,
                            discreet_mode,
                            default_cutoff,
                            parallel_mode,
                            log_file,
                            verbose,
                            matrix_generator_path)


def runReplicationStrandBiasAnalysis(outputDir,
                                     jobname,
                                     numofSimulations,
                                     job_tuples,
                                     sample_based,
                                     all_samples_np_array,
                                     replicationTimeFilename,
                                     replicationTimeValleyFilename,
                                     replicationTimePeakFilename,
                                     chromSizesDict,
                                     chromNamesList,
                                     computation_type,
                                     ordered_sbs_signatures,
                                     ordered_dbs_signatures,
                                     ordered_id_signatures,
                                     ordered_sbs_signatures_cutoffs,
                                     ordered_dbs_signatures_cutoffs,
                                     ordered_id_signatures_cutoffs,
                                     discreet_mode,
                                     default_cutoff,
                                     parallel_mode,
                                     log_file,
                                     verbose):

    os.makedirs(os.path.join(outputDir,jobname,DATA,REPLICATIONSTRANDBIAS), exist_ok=True)

    smoothedWaveletRepliseqDataFilename = replicationTimeFilename
    valleysBEDFilename = replicationTimeValleyFilename
    peaksBEDFilename = replicationTimePeakFilename

    # Supported computation types
    # computation_type= USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM
    # computation_type =USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT
    replicationStrandBiasAnalysis(outputDir,
                                  jobname,
                                  numofSimulations,
                                  job_tuples,
                                  sample_based,
                                  all_samples_np_array,
                                  chromSizesDict,
                                  chromNamesList,
                                  computation_type,
                                  smoothedWaveletRepliseqDataFilename,
                                  valleysBEDFilename,
                                  peaksBEDFilename,
                                  ordered_sbs_signatures,
                                  ordered_dbs_signatures,
                                  ordered_id_signatures,
                                  ordered_sbs_signatures_cutoffs,
                                  ordered_dbs_signatures_cutoffs,
                                  ordered_id_signatures_cutoffs,
                                  discreet_mode,
                                  default_cutoff,
                                  parallel_mode,
                                  log_file,
                                  verbose)


def runTranscriptionStradBiasAnalysis(outputDir,
                                      jobname,
                                      numofSimulations,
                                      job_tuples,
                                      sample_based,
                                      all_samples_np_array,
                                      chromNamesList,
                                      computation_type,
                                      ordered_sbs_signatures,
                                      ordered_dbs_signatures,
                                      ordered_id_signatures,
                                      ordered_sbs_signatures_cutoffs,
                                      ordered_dbs_signatures_cutoffs,
                                      ordered_id_signatures_cutoffs,
                                      discreet_mode,
                                      default_cutoff,
                                      parallel_mode,
                                      log_file,
                                      verbose):

    os.makedirs(os.path.join(outputDir, jobname, DATA, TRANSCRIPTIONSTRANDBIAS), exist_ok=True)

    # Supported computation types
    # computation_type = USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM
    # computation_type = USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT
    transcriptionStrandBiasAnalysis(outputDir,
                                    jobname,
                                    numofSimulations,
                                    job_tuples,
                                    sample_based,
                                    all_samples_np_array,
                                    computation_type,
                                    chromNamesList,
                                    ordered_sbs_signatures,
                                    ordered_dbs_signatures,
                                    ordered_id_signatures,
                                    ordered_sbs_signatures_cutoffs,
                                    ordered_dbs_signatures_cutoffs,
                                    ordered_id_signatures_cutoffs,
                                    discreet_mode,
                                    default_cutoff,
                                    parallel_mode,
                                    log_file,
                                    verbose)


def runProcessivityAnalysis(mutation_types_contexts,
                            outputDir,
                            jobname,
                            numofSimulations,
                            chromNamesList,
                            processivity_calculation_type,
                            processivity_inter_mutational_distance,
                            subsSignature_cutoff_numberofmutations_averageprobability_df,
                            parallel_mode,
                            log_file,
                            verbose):

    os.makedirs(os.path.join(outputDir,jobname,DATA,PROCESSIVITY),exist_ok=True)

    # Internally Set
    considerProbabilityInProcessivityAnalysis = True

    processivityAnalysis(mutation_types_contexts,
                         chromNamesList,
                         processivity_calculation_type,
                         processivity_inter_mutational_distance,
                         outputDir,
                         jobname,
                         numofSimulations,
                         considerProbabilityInProcessivityAnalysis,
                         subsSignature_cutoff_numberofmutations_averageprobability_df,
                         parallel_mode,
                         log_file,
                         verbose)


def deleteOldData(outputDir,jobname, occupancy_type):
    # Delete the output/jobname/DATA/occupancy_type if exists
    jobnamePath = os.path.join(outputDir,jobname,DATA,occupancy_type)

    if (os.path.exists(jobnamePath)):
        try:
            shutil.rmtree(jobnamePath)
        except OSError as e:
            print('Error: %s - %s.' % (e.filename, e.strerror))


def deleteOldFigures(outputDir, jobname, occupancy_type):

    jobnamePath = os.path.join(outputDir, jobname, FIGURE, occupancy_type)
    print('Topography.py jobnamePath:%s ' %jobnamePath)

    if (os.path.exists(jobnamePath)):
        try:
            shutil.rmtree(jobnamePath)
        except OSError as e:
            print('Error: %s - %s.' % (e.filename, e.strerror))


# We assume that simulated data will have the same number_of_splits as the real data
def get_job_tuples(chrlong_numberofmutations_df,numofSimulations):

    job_tuples = []
    sim_nums = range(0, numofSimulations + 1)

    for chrLong in chrlong_numberofmutations_df['chrLong'].unique():
        number_of_mutations=int(chrlong_numberofmutations_df[chrlong_numberofmutations_df['chrLong']==chrLong]['number_of_mutations'].values[0])
        number_of_splits = math.ceil(number_of_mutations / NUMBER_OF_MUTATIONS_IN_EACH_SPLIT)
        split_indexes = range(0, number_of_splits)

        ###############################################################
        for sim_num in sim_nums:
            for split_index in split_indexes:
                job_tuples.append((chrLong, sim_num, split_index))
        ###############################################################

    return job_tuples


def get_all_signatures_array(ordered_all_sbs_signatures_wrt_probabilities_file_array, signature_starts_with):
    ordered_all_sbs_signatures = []

    if ordered_all_sbs_signatures_wrt_probabilities_file_array is not None:
        for i in ordered_all_sbs_signatures_wrt_probabilities_file_array:
            if i.startswith(signature_starts_with):
                ordered_all_sbs_signatures.append(i)

    return np.array(ordered_all_sbs_signatures)


#######################################################
# inputDir ='/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/input_for_matgen/BreastCancer560_subs_indels_dinucs'
# outputDir = '/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/output_test/'
# jobname = 'BreastCancer560'

# Run SigProfilerTopography Analyses
# Former full path now only the filename with extension
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
                gender = 'female', # 'male'
                sbs_probabilities = None,
                dbs_probabilities = None,
                id_probabilities = None,
                mutation_types_contexts = None,
                mutation_types_contexts_for_signature_probabilities = None,
                epigenomics_files = None,
                epigenomics_files_memos = None,
                epigenomics_biosamples = None,
                epigenomics_dna_elements = None,
                epigenomics_dir_name = None,
                nucleosome_biosample = None,
                nucleosome_file = None,
                replication_time_biosample = None,
                replication_time_signal_file = None,
                replication_time_valley_file = None,
                replication_time_peak_file = None,
                computation_type = USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM,
                epigenomics = False,
                nucleosome = False,
                replication_time = False,
                strand_bias = False,
                replication_strand_bias = False,
                transcription_strand_bias = False,
                processivity = False,
                mutation_annotation_integration=False,
                sample_based = False,
                plot_figures = True,
                step1_sim_data = True,
                step2_matgen_data = True,
                step3_prob_merged_data = True,
                step4_tables = True,
                discreet_mode = True, # discreet_mode = False means prob_mode
                show_all_signatures = False,
                average_probability = DEFAULT_AVERAGE_PROBABILITY, # applies for discreet_mode=True
                default_cutoff = DEFAULT_CUTOFF, # applies for discreet_mode=False
                parallel_mode = True,
                num_of_sbs_required = DEFAULT_NUM_OF_SBS_REQUIRED,
                num_of_dbs_required = DEFAULT_NUM_OF_DBS_REQUIRED,
                num_of_id_required = DEFAULT_NUM_OF_ID_REQUIRED,
                plusorMinus_epigenomics = 1000,
                plusorMinus_nucleosome = 1000,
                epigenomics_heatmap_significance_level = 0.05,
                processivity_significance_level = 0.05,
                verbose = False,
                matrix_generator_path = MATRIX_GENERATOR_PATH,
                PCAWG = False,
                plot_epigenomics = False,
                plot_nucleosome = False,
                plot_replication_time = False,
                plot_strand_bias = False,
                plot_replication_strand_bias = False,
                plot_transcription_strand_bias = False,
                plot_processivity = False,
                remove_outliers = False,
                quantileValue = 0.97,
                delete_old = False,
                plot_mode = PLOTTING_FOR_SIGPROFILERTOPOGRAPHY_TOOL,
                occupancy_calculation_type = MISSING_SIGNAL,
                processivity_calculation_type = CONSIDER_DISTANCE,
                processivity_inter_mutational_distance = 10000,
                combine_p_values_method = COMBINE_P_VALUES_METHOD_FISHER,
                fold_change_window_size = 100,
                num_of_real_data_avg_overlap = DEFAULT_NUM_OF_REAL_DATA_OVERLAP_REQUIRED):

    current_abs_path = os.path.dirname(os.path.realpath(__file__))

    chromSizesDict = getChromSizesDict(genome)
    chromNamesList = list(chromSizesDict.keys())
    chromShortNamesList=getShortNames(chromNamesList)

    # Filled in Step3
    # contains all the columns in order w.r.t. probabilities file
    ordered_all_sbs_signatures_wrt_probabilities_file_array = None
    ordered_all_dbs_signatures_wrt_probabilities_file_array = None
    ordered_all_id_signatures_wrt_probabilities_file_array = None

    if mutation_types_contexts is None:
        mutation_types_contexts=[]
        if (sbs_probabilities is not None):
            mutation_types_contexts.append(SBS96)
        if (id_probabilities is not None):
            mutation_types_contexts.append(ID)
        if (dbs_probabilities is not None):
            mutation_types_contexts.append(DBS)

    ############################## Log and Error Files #######################################
    time_stamp = datetime.date.today()

    current_hour = datetime.datetime.now().hour
    current_minute = datetime.datetime.now().minute
    current_second = datetime.datetime.now().second

    error_file = os.path.join(outputDir, jobname , 'logs', 'SigProfilerTopography_' + jobname + '_' + genome + '_' +
                              str(time_stamp) + '_' +
                              str(current_hour) + '-' + str(current_minute) + '-' + str(current_second) + '.err')

    log_file = os.path.join(outputDir, jobname, 'logs', 'SigProfilerTopography_' + jobname + '_' + genome + '_' +
                            str(time_stamp) + '_' +
                            str(current_hour) + '-' + str(current_minute) + '-' + str(current_second) + '.out')

    if not os.path.exists(os.path.join(outputDir, jobname, 'logs/')):
        os.makedirs(os.path.join(outputDir, jobname, 'logs/'))

    if os.path.exists(error_file):
        os.remove(error_file)
    if os.path.exists(log_file):
        os.remove(log_file)

    sys.stderr = open(error_file, 'w')
    log_out = open(log_file, 'w')
    log_out.write("THIS FILE CONTAINS THE METADATA ABOUT SYSTEM AND RUNTIME\n\n\n")
    log_out.write("-------System Info-------\n")

    log_out.write("Operating System Name: " + platform.platform() + "\n" +
                  "Operating System: " + platform.uname()[0] + "\n" +
                  "Nodename: " + platform.uname()[1] + "\n" +
                  "Release: " + platform.uname()[2] + "\n" +
                  "Version: " + platform.uname()[3] + "\n")

    log_out.write("\n-------Python and Package Versions------- \n")
    log_out.write("Python Version: " +
                  str(platform.sys.version_info.major) + "." +
                  str(platform.sys.version_info.minor) + "." +
                  str(platform.sys.version_info.micro) + "\n")
    log_out.write("\n-------Date and Time Data------- \n")
    tic = datetime.datetime.now()
    log_out.write("Date and Clock time when the execution started: " + str(tic) + "\n\n\n")

    ############################## Log and Error Files #######################################

    # If still None
    if mutation_types_contexts is None:
        print('--- There is a situation/problem: mutation_types_contexts is None.', file=log_out)
        print('--- mutation_types_contexts has to be set before SigProfilerTopography run.', file=log_out)

    if mutation_types_contexts_for_signature_probabilities is None:
        mutation_types_contexts_for_signature_probabilities=mutation_types_contexts

    if step1_sim_data:
        step2_matgen_data = True
        step3_prob_merged_data = True
        step4_tables = True
    elif step2_matgen_data:
        step3_prob_merged_data = True
        step4_tables = True
    elif step3_prob_merged_data:
        step4_tables = True

    if (average_probability!=DEFAULT_AVERAGE_PROBABILITY) or \
            (num_of_sbs_required!=DEFAULT_NUM_OF_SBS_REQUIRED) or \
            (num_of_dbs_required!=DEFAULT_NUM_OF_DBS_REQUIRED) or \
            (num_of_id_required!=DEFAULT_NUM_OF_ID_REQUIRED):
        step4_tables = True

    #################################################################################
    ################################## Setting starts ###############################
    ################## Set full path library files starts ###########################
    #################################################################################
    if genome is None:
        print('There is a situation/problem: Parameter genome:%s must be set for SigProfilerTopography Analysis.' %(genome), file=log_out)

    if strand_bias:
        replication_strand_bias = True
        transcription_strand_bias = True

    if plot_strand_bias:
        plot_replication_strand_bias = True
        plot_transcription_strand_bias = True

    # Epigenomics Occupancy
    # We need full path of the library files
    if (genome == GRCh37) and (epigenomics_files == None):
        epigenomics_files = [DEFAULT_ATAC_SEQ_OCCUPANCY_FILE,
                            DEFAULT_H3K27ME3_OCCUPANCY_FILE,
                            DEFAULT_H3K36ME3_OCCUPANCY_FILE,
                            DEFAULT_H3K9ME3_OCCUPANCY_FILE,
                            DEFAULT_H3K27AC_OCCUPANCY_FILE,
                            DEFAULT_H3K4ME1_OCCUPANCY_FILE,
                            DEFAULT_H3K4ME3_OCCUPANCY_FILE,
                            DEFAULT_CTCF_OCCUPANCY_FILE]

        epigenomics_files_memos = []
        for epigenomics_file in epigenomics_files:
            epigenomics_files_memos.append(os.path.splitext(os.path.basename(epigenomics_file))[0])

        # Defines columns in the heatmap
        # These strings must be within filenames (without file extension)
        # Order is not important
        epigenomics_dna_elements = ['H3K27me3', 'H3K36me3', 'H3K9me3', 'H3K27ac', 'H3K4me1', 'H3K4me3', 'CTCF', 'ATAC']

        # Defines rows in the detailed heatmap
        # These strings must be within filenames (without file extension)
        # Order is not important
        epigenomics_biosamples = ['breast_epithelium']

        for file_index, filename in enumerate(epigenomics_files):
            epigenomics_files[file_index] = os.path.join(current_abs_path, LIB, EPIGENOMICS, filename)
        # These must be under epigenomics under installed SigPofilerTopography

    elif (genome == GRCh38) and (epigenomics_files == None):
        epigenomics_files = [DEFAULT_ATAC_SEQ_GRCh38_OCCUPANCY_FILE,
                            DEFAULT_H3K27ME3_GRCh38_OCCUPANCY_FILE,
                            DEFAULT_H3K36ME3_GRCh38_OCCUPANCY_FILE,
                            DEFAULT_H3K9ME3_GRCh38_OCCUPANCY_FILE,
                            DEFAULT_H3K27AC_GRCh38_OCCUPANCY_FILE,
                            DEFAULT_H3K4ME1_GRCh38_OCCUPANCY_FILE,
                            DEFAULT_H3K4ME3_GRCh38_OCCUPANCY_FILE,
                            DEFAULT_CTCF_GRCh38_OCCUPANCY_FILE]

        epigenomics_files_memos = []
        for epigenomics_file in epigenomics_files:
            epigenomics_files_memos.append(os.path.splitext(os.path.basename(epigenomics_file))[0])

        # Defines columns in the heatmap
        # These strings must be within filenames (without file extension)
        # Order is not important
        epigenomics_dna_elements = ['H3K27me3', 'H3K36me3', 'H3K9me3', 'H3K27ac', 'H3K4me1', 'H3K4me3', 'CTCF', 'ATAC']

        # Defines rows in the detailed heatmap
        # These strings must be within filenames (without file extension)
        # Order is not important
        epigenomics_biosamples = ['lung']

        for file_index, filename in enumerate(epigenomics_files):
            epigenomics_files[file_index] = os.path.join(current_abs_path, LIB, EPIGENOMICS, filename)
        # These files must be under epigenomics under installed SigProfilerTopography

    elif (genome == MM10) and (epigenomics_files == None):
        epigenomics_files = [ENCFF575PMI_mm10_embryonic_facial_prominence_ATAC_seq,
                            ENCFF993SRY_mm10_embryonic_fibroblast_H3K4me1,
                            ENCFF912DNP_mm10_embryonic_fibroblast_H3K4me3,
                            ENCFF611HDQ_mm10_embryonic_fibroblast_CTCF,
                            ENCFF152DUV_mm10_embryonic_fibroblast_POLR2A,
                            ENCFF114VLZ_mm10_embryonic_fibroblast_H3K27ac]

        epigenomics_files_memos = []
        for epigenomics_file in epigenomics_files:
            epigenomics_files_memos.append(os.path.splitext(os.path.basename(epigenomics_file))[0])

        # Defines columns in the heatmap
        # These strings must be within filenames (without file extension)
        # Order is not important
        epigenomics_dna_elements = ['ATAC', 'H3K4me1', 'H3K4me3', 'CTCF', 'POLR2A', 'H3K27ac']

        # Defines rows in the detailed heatmap
        # These strings must be within filenames (without file extension)
        # Order is not important
        epigenomics_biosamples = ['embryonic_fibroblast']

        for file_index, filename in enumerate(epigenomics_files):
            epigenomics_files[file_index] = os.path.join(current_abs_path, LIB, EPIGENOMICS, filename)

    # Nucleosome Occupancy
    if genome == MM10:
        # Case1: File is not set, Biosample is not set
        if (nucleosome_file is None) and (nucleosome_biosample is None):
            nucleosome_biosample = ESC
            nucleosome_file = getNucleosomeFile(genome, nucleosome_biosample)
        # Case2: File is not set, Biosample is set
        elif (nucleosome_file is None) and (nucleosome_biosample is not None):
            if (nucleosome_biosample in available_nucleosome_biosamples):
                #Sets the filename without the full path
                nucleosome_file = getNucleosomeFile(genome, nucleosome_biosample)
        # Case3: nucleosome_file is a filename with fullpath (User provided) , biosample is not set
        elif ((nucleosome_file is not None) and (nucleosome_biosample is None)):
            # We expect that user has provided nucleosome file with full path
            nucleosome_biosample = UNDECLARED
        # Case4: nucleosome_file is a filename with fullpath (User provided), biosample is set
        # Do nothing use as it is

    elif genome == GRCh37:
        # Case1: File is not set, Biosample is not set
        if (nucleosome_file is None) and (nucleosome_biosample is None):
            nucleosome_biosample = K562
            nucleosome_file = getNucleosomeFile(genome, nucleosome_biosample)
        # Case2: File is not set, Biosample is set
        elif (nucleosome_file is None) and (nucleosome_biosample is not None):
            if (nucleosome_biosample in available_nucleosome_biosamples):
                #Sets the filename without the full path
                nucleosome_file = getNucleosomeFile(genome, nucleosome_biosample)
        # Case3: nucleosome_file is a filename with fullpath (User provided) , biosample is not set
        elif ((nucleosome_file is not None) and (nucleosome_biosample is None)):
            # We expect that user has provided nucleosome file with full path
            nucleosome_biosample = UNDECLARED
        # Case4: nucleosome_file is a filename with fullpath (User provided), biosample is set
        # Do nothing use as it is

    elif genome == GRCh38:
        # Case1: File is not set, Biosample is not set
        if (nucleosome_file is None) and (nucleosome_biosample is None):
            nucleosome_biosample = K562
            nucleosome_file = getNucleosomeFile(genome, nucleosome_biosample)
        # Case2: File is not set, Biosample is set
        elif (nucleosome_file is None) and (nucleosome_biosample is not None):
            if (nucleosome_biosample in available_nucleosome_biosamples):
                #Sets the filename without the full path
                nucleosome_file = getNucleosomeFile(genome, nucleosome_biosample)
        # Case3: nucleosome_file is a filename with fullpath (User provided) , biosample is not set
        elif ((nucleosome_file is not None) and (nucleosome_biosample is None)):
            # We expect that user has provided nucleosome file with full path
            nucleosome_biosample = UNDECLARED
        # Case4: nucleosome_file is a filename with fullpath (User provided), biosample is set
        # Do nothing use as it is

    # Replication Timing
    if genome == MM10:
        # Case1: Files are not set, Biosample is not set
        if (replication_time_signal_file is None) and (replication_time_valley_file is None) and (replication_time_peak_file is None) and (replication_time_biosample is None):
            replication_time_biosample = MEF
            # We only set replication_time_signal_file
            # replication_time_valley_file is None
            # replication_time_peak_file is None
            replication_time_signal_file, replication_time_valley_file,replication_time_peak_file = getReplicationTimeFiles(replication_time_biosample, genome)

    elif genome == GRCh37:
        # We need full path of the library files
        # By default replication_time_biosample=MCF7 and signal, valley, peak files are None
        # Case1: Files are not set, Biosample is not set
        if (replication_time_signal_file is None) and (replication_time_valley_file is None) and (replication_time_peak_file is None) and (replication_time_biosample is None):
            replication_time_biosample = MCF7
            replication_time_signal_file, replication_time_valley_file,replication_time_peak_file = getReplicationTimeFiles(replication_time_biosample, genome)
            if (replication_time or replication_strand_bias):
                # For using SigProfilerTopography Provided Replication Time Files
                check_download_replication_time_files(replication_time_signal_file, replication_time_valley_file,replication_time_peak_file)

        # Case2: Files are not set, Biosample is set
        elif (replication_time_signal_file is None) and (replication_time_valley_file is None) and (replication_time_peak_file is None) and (replication_time_biosample is not None):
            if (replication_time_biosample in available_replication_time_biosamples):
                replication_time_signal_file, replication_time_valley_file, replication_time_peak_file = getReplicationTimeFiles(replication_time_biosample, genome)
                if (replication_time or replication_strand_bias):
                    # For using SigProfilerTopography Provided Replication Time Files
                    check_download_replication_time_files(replication_time_signal_file, replication_time_valley_file,replication_time_peak_file)

        # Case3: replication_time_file is a filename with fullpath (User provided) , biosample is not set
        elif ((replication_time_signal_file is not None) or (replication_time_valley_file is not None) or (replication_time_peak_file is not None)) and (replication_time_biosample is None):
            replication_time_biosample = UNDECLARED
        # Case4: Files are set. Biosample is set. Use as it is. Do nothing.

    elif genome == GRCh38:
        # We need full path of the library files
        # By default replication_time_biosample = IMR90 and signal, valley, peak files are None
        # Case1: Files are not set, Biosample is not set
        if (replication_time_signal_file is None) and (replication_time_valley_file is None) and (replication_time_peak_file is None) and (replication_time_biosample is None):
            replication_time_biosample = IMR90
            replication_time_signal_file, replication_time_valley_file, replication_time_peak_file = getReplicationTimeFiles(replication_time_biosample, genome)
            if (replication_time or replication_strand_bias):
                # For using SigProfilerTopography Provided Replication Time Files
                check_download_replication_time_files(replication_time_signal_file, replication_time_valley_file, replication_time_peak_file)

        # Case4: Files are set. Biosample is set. Use as it is. Do nothing. Use them .
        pass

    # Data files are named using user provided epigenomics_files_memos or using epigenomics_file_memos_created
    epigenomics_file_memos_created = []

    # Run for each epigenomics file
    if (epigenomics_files_memos is None) or (len(epigenomics_files_memos) != len(epigenomics_files)):
        if epigenomics_files is not None:
            for idx, epigenomics_file in enumerate(epigenomics_files):
                epigenomics_file_memo = os.path.splitext(os.path.basename(epigenomics_file))[0]
                epigenomics_file_memos_created.append(epigenomics_file_memo)

    # Used for plotting
    if (epigenomics_files_memos is None) or (len(epigenomics_files_memos) != len(epigenomics_files)):
        epigenomics_files_memos = epigenomics_file_memos_created

    if (epigenomics_biosamples is None) or (len(epigenomics_biosamples) == 0):
        epigenomics_biosamples = [UNDECLARED]

    #################################################################################
    ################## Set full path library files ends #############################
    ################################## Setting ends #################################
    #################################################################################

    print('#################################################################################', file=log_out)
    print("--- SigProfilerTopography starts", file=log_out)
    print('#################################################################################', file=log_out)

    print('#################################################################################', file=log_out)
    print('--- SigProfilerTopography Version:%s' % topography_version.version, file=log_out)
    print("--- SigProfilerMatrixGenerator Version: %s" %matrix_generator_version.version, file=log_out)
    print("--- SigProfilerSimulator version: %s" %simulator_version.version, file=log_out)
    print("--- pandas version: %s" %pd.__version__, file=log_out)
    print("--- numpy version: %s" %np.__version__, file=log_out)
    print("--- statsmodels version: %s" %statsmodels.__version__, file=log_out)
    print("--- scipy version: %s" %scipy.__version__, file=log_out)
    print("--- matplotlib version: %s" %plt.__version__, file=log_out)
    print('#################################################################################\n', file=log_out)

    print('#################################################################################', file=log_out)
    print('--- SigProfilerTopography parameters', file=log_out)
    print('--- Genome: %s' %(genome), file=log_out)
    print('--- inputDir:%s' %inputDir, file=log_out)
    print('--- outputDir:%s' %outputDir, file=log_out)
    print('--- jobname:%s' %jobname, file=log_out)
    if (sbs_probabilities is not None):
        print('--- sbs_probabilities:%s' %sbs_probabilities, file=log_out)
    if (dbs_probabilities is not None):
        print('--- dbs_probabilities:%s' %dbs_probabilities, file=log_out)
    if (id_probabilities is not None):
        print('--- id_probabilities:%s' %id_probabilities, file=log_out)

    print('--- numofSimulations:%d' %numofSimulations, file=log_out)
    print('\n--- epigenomics_files:%s' %epigenomics_files, file=log_out)
    print('--- epigenomics_files_memos:%s' %epigenomics_files_memos, file=log_out)
    print('--- epigenomics_biosamples:%s' %epigenomics_biosamples, file=log_out)
    print('--- epigenomics_dna_elements:%s' %epigenomics_dna_elements, file=log_out)
    if epigenomics_files is not None:
        print('--- number of epigenomics_files:%d' %len(epigenomics_files), file=log_out)

    print('\n--- nucleosome_biosample:%s' %nucleosome_biosample, file=log_out)
    print('--- nucleosome_file:%s' % nucleosome_file, file=log_out)

    print('\n--- replication_time_biosample:%s' % replication_time_biosample, file=log_out)
    print('--- replication_time_signal_file:%s' % replication_time_signal_file, file=log_out)
    print('--- replication_time_valley_file:%s' % replication_time_valley_file, file=log_out)
    print('--- replication_time_peak_file:%s' % replication_time_peak_file, file=log_out)

    print('\n--- mutation_types_contexts:%s' %mutation_types_contexts, file=log_out)
    print('--- mutation_types_contexts_for_signature_probabilities:%s' %mutation_types_contexts_for_signature_probabilities, file=log_out)
    print('--- computation_type:%s' %computation_type, file=log_out)
    print('--- SigProfilerTopography run mode discreet:%s\n' %discreet_mode, file=log_out)
    if sample_based:
        print('--- Sample Based Analysis.', file=log_out)

    if epigenomics:
        print('--- Epigenomics Analysis.', file=log_out)
    if nucleosome:
        print('--- Nucleosome Analysis.', file=log_out)
    if replication_time:
        print('--- Replication Timing Analysis.', file=log_out)
    if (strand_bias or replication_strand_bias):
        print('--- Replication Strand Asymmetry Analysis.', file=log_out)
    if (strand_bias or transcription_strand_bias):
        print('--- Transcription Strand Asymmetry Analysis.', file=log_out)
    if processivity:
        print('--- Strand-coordinated Mutagenesis Analysis.', file=log_out)

    print('--- step1_sim_data:%s' %step1_sim_data, file=log_out)
    print('--- step2_matgen_data:%s' %step2_matgen_data, file=log_out)
    print('--- step3_prob_merged_data:%s' %step3_prob_merged_data, file=log_out)
    print('--- step4_tables:%s' %step4_tables, file=log_out)

    print('--- plot_figures:%s' %plot_figures, file=log_out)
    if discreet_mode:
        print('--- discreet_mode: %s' %discreet_mode, file=log_out)
        print('--- average mutation probability required: %0.2f' %average_probability, file=log_out)
    else:
        print('--- discreet_mode: %s' %discreet_mode, file=log_out)
        print('--- default_cutoff: %s' %default_cutoff, file=log_out)

    print('--- parallel_mode: %s' % parallel_mode, file=log_out)
    print('--- minimum number of sbs mutations required: %d' %num_of_sbs_required, file=log_out)
    print('--- minimum number of id mutations required: %d' %num_of_id_required, file=log_out)
    print('--- minimum number of dbs mutations required: %d' %num_of_dbs_required, file=log_out)
    if epigenomics:
        print('--- number of bases considered before and after mutation start for epigenomics analysis: %d' %plusorMinus_epigenomics, file=log_out)
    if nucleosome:
        print('--- number of bases considered before and after mutation start for nucleosome occupancy analysis: %d' %plusorMinus_nucleosome, file=log_out)
    print('#################################################################################\n', file=log_out)

    print('#################################################################################', file=log_out)
    numofProcesses = multiprocessing.cpu_count()
    print('--- numofProcesses for multiprocessing: %d' %numofProcesses, file=log_out)
    print('#################################################################################\n', file=log_out)


    print('#################################################################################', file=log_out)
    print('--- For Genome: %s' %(genome), file=log_out)
    print('--- Chromosome names: %s' %(chromNamesList), file=log_out)
    print('--- Chromosome short names: %s' % (chromShortNamesList), file=log_out)
    print('--- current_abs_path: %s ' % current_abs_path, file=log_out)
    print('#################################################################################\n', file=log_out)

    ###################################################################################################################
    ################################################# All Steps starts ################################################
    ###################################################################################################################

    ###################################################################################################
    ######################### SigProfilerMatrixGenerator for original data starts #####################
    ###################################################################################################
    if (step2_matgen_data):

        # Run MatrixGenerator for original data: this call prepares chrBased input files for original data with mutation contexts
        print('#################################################################################', file=log_out)
        print('--- SigProfilerMatrixGenerator for original data', file=log_out)
        start_time = time.time()

        print('For original data inputDir:%s' % (inputDir), file=log_out)
        matrices = matGen.SigProfilerMatrixGeneratorFunc(jobname, genome, inputDir, plot=False, seqInfo=True)

        # original matrix generator chrbased data will be under inputDir/output/vcf_files/SNV
        # original matrix generator chrbased data will be under inputDir/output/vcf_files/DBS
        # original matrix generator chrbased data will be under inputDir/output/vcf_files/ID

        print("--- SigProfilerMatrixGenerator for original data: %s seconds ---" % (time.time() - start_time), file=log_out)
        print("--- SigProfilerMatrixGenerator for original data: %f minutess ---" % float((time.time() - start_time) / 60), file=log_out)
        print('#################################################################################\n', file=log_out)
    ###################################################################################################
    ######################### SigProfilerMatrixGenerator for original data ends #######################
    ###################################################################################################

    ###################################################################################################################
    ################################## Step1 Simulations if any starts ################################################
    ###################################################################################################################
    if ((numofSimulations > 0) and (step1_sim_data)):

        ###################################################################################################
        ############################  SigProfilerSimulator for n simulations starts #######################
        ###################################################################################################
        print('#################################################################################', file=log_out)
        print('--- SigProfilerSimulator for %d simulations starts' %(numofSimulations), file=log_out)
        start_time = time.time()
        #Call SigProfilerSimulator separately for each mutation type context otherwise it counts DBS mutations also in SBS mutations
        # Topography uses same mutation types with Simulator
        # Acceptable contexts for Simulator include {'96', '384', '1536', '6144', 'DBS', 'ID', 'ID415'}.
        # '96' or '384' for single base substitutions (Simulator 1536, or 3072)
        # 'DBS' for double base substitutions
        # 'ID' for indels
        for mutation_type_context in mutation_types_contexts:
            mutation_type_context_for_simulator = []
            mutation_type_context_for_simulator.append(mutation_type_context)
            # Please notice that Simulator reverse the given input mutationTypes_for_simulator
            print('--- SigProfilerSimulator is running for %s' %(mutation_type_context), file=log_out)
            simulator.SigProfilerSimulator(jobname, inputDir, genome, mutation_type_context_for_simulator,
                                            gender=gender,
                                            simulations=numofSimulations,
                                            chrom_based=True)
        print("--- SigProfilerSimulator for %d simulations: %s seconds" %(numofSimulations,(time.time() -  start_time)), file=log_out)
        print("--- SigProfilerSimulator for %d simulations: %f minutes" %(numofSimulations,float((time.time()-start_time)/60)), file=log_out)
        print('--- SigProfilerSimulator for %d simulations ends' %(numofSimulations), file=log_out)
        print('#################################################################################\n', file=log_out)
        ###################################################################################################
        ############################  SigProfilerSimulator for n simulations ends #########################
        ###################################################################################################

    ###################################################################################################################
    ################################## Step1 Simulations if any ends ##################################################
    ###################################################################################################################

    ###################################################################################################################
    ################################## Step2 Matrix Generator for n simulations starts ################################
    ###################################################################################################################
    if (step2_matgen_data):

        if (numofSimulations > 0):
            ###################################################################################################
            ########################### Create simN directories for MatrixGenerator starts ####################
            ###################################################################################################

            print('#################################################################################', file=log_out)
            print('--- Create directories for %d simulations under %s/output/simulations/' %(numofSimulations,inputDir), file=log_out)
            start_time = time.time()
            # Create directories sim1 to SimN under inputDir/output/simulations/
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
                    print("Creation of the directory %s failed" %simDir, file=log_out)

            for mutation_type_context in mutation_types_contexts:
                # Simulator creates one maf file for each simulation for each mutation context
                # Simulator creates maf files under inputDir/output/simulations/jobname_simulations_GRCh37_96
                # Simulator creates maf files under inputDir/output/simulations/jobname_simulations_GRCh37_ID
                # Simulator creates maf files under inputDir/output/simulations/jobname_simulations_GRCh37_DBS
                dirName = '%s_simulations_%s_%s' %(jobname, genome,mutation_type_context)
                copyFromDir = os.path.join(inputDir,'output','simulations',dirName)
                copyToMainDir= os.path.join(inputDir,'output','simulations')

                # Topography copies these maf files to inputDir/output/simulations/simX/mutation_type_context/X.maf
                # So that, in the next step MatrixGenerator can create chrom based seqinfo text files for each X.maf file
                copyMafFiles(copyFromDir,copyToMainDir,mutation_type_context,numofSimulations)
            print("--- Create directories and copy files: %s seconds ---" %(time.time()-start_time), file=log_out)
            print("--- Create directories and copy files: %f minutes ---" %(float((time.time()-start_time)/60)), file=log_out)
            print('#################################################################################\n', file=log_out)

            ###################################################################################################
            ########################### Create simN directories for MatrixGenerator ends ######################
            ###################################################################################################

            ###################################################################################################
            # Important note: Separate directory creation is necessary for Matrix Generator
            # inputDir/output/simulations/simX/96/X.maf
            # inputDir/output/simulations/simX/ID/X.maf
            # inputDir/output/simulations/simX/DBS/X.maf

            # enables MatrixGenerator to create chr based simulated data files under
            # sim1 matrix generator chrbased data will be under inputDir/output/simulations/simX/96/output/vcf_files/SNV
            # sim1 matrix generator chrbased data will be under inputDir/output/simulations/simX/ID/output/vcf_files/ID
            # sim1 matrix generator chrbased data will be under inputDir/output/simulations/simX/DBS/output/vcf_files/DBS

            # otherwise all simulations maf files will be under
            # inputDir/output/simulations/Skin-Melanoma_simulations_GRCh37_96
            # inputDir/output/simulations/Skin-Melanoma_simulations_GRCh37_DBS
            # inputDir/output/simulations/Skin-Melanoma_simulations_GRCh37_ID
            # Then running MatrixGenerator for each simulation will not be possible.
            ###################################################################################################

            ###################################################################################################
            ####################### Run MatrixGenerator for each simulation starts ############################
            ###################################################################################################
            print('#################################################################################', file=log_out)
            print('--- Run SigProfilerMatrixGenerator for each simulation starts', file=log_out)
            start_time = time.time()
            for simNum in range(1,numofSimulations+1):
                simName = 'sim%d' %(simNum)
                #For each simulation we are calling matrix generator separately for each mutation type context

                print('--- SigProfilerMatrixGenerator is run for %s starts' %(simName), file=log_out)
                for mutation_type_context in mutation_types_contexts:
                    simInputDir=  os.path.join(inputDir,'output','simulations',simName,mutation_type_context)
                    print('--- For %s: %s simInputDir:%s' %(mutation_type_context,simName,simInputDir), file=log_out)
                    matrices = matGen.SigProfilerMatrixGeneratorFunc(jobname,genome,simInputDir,plot=False, seqInfo=True)
                print('--- SigProfilerMatrixGenerator is run for %s ends\n' % (simName), file=log_out)
            #sim1 matrix generator chrbased data will be under inputDir/output/simulations/sim1/96/output/vcf_files/SNV
            #sim1 matrix generator chrbased data will be under inputDir/output/simulations/sim1/ID/output/vcf_files/ID
            #sim1 matrix generator chrbased data will be under inputDir/output/simulations/sim1/DBS/output/vcf_files/DBS

            #simN matrix generator chrbased data will be under inputDir/output/simulations/simN/96/output/vcf_files/SNV
            #simN matrix generator chrbased data will be under inputDir/output/simulations/simN/ID/output/vcf_files/ID
            #simN matrix generator chrbased data will be under inputDir/output/simulations/simN/DBS/output/vcf_files/DBS
            print("--- Run MatrixGenerator for each simulation: %s seconds" %(time.time()-start_time), file=log_out)
            print("--- Run MatrixGenerator for each simulation: %f minutes" %(float((time.time()-start_time)/60)), file=log_out)
            print('--- Run SigProfilerMatrixGenerator for each simulation ends', file=log_out)
            print('#################################################################################\n', file=log_out)
            ###################################################################################################
            ####################### Run MatrixGenerator for each simulation ends ##############################
            ###################################################################################################

    ###################################################################################################################
    ################################## Step2 Matrix Generator for n simulations ends ##################################
    ###################################################################################################################

    log_out.close()

    # create DATA directory before step3 and step4
    os.makedirs(os.path.join(outputDir, jobname, DATA), exist_ok=True)

    ###################################################################################################################
    ########### Step3 Merge chrom based matrix generator generated files with probabilities starts ####################
    ###################################################################################################################
    if (step3_prob_merged_data):
        ####################################################################################################################
        ##################  Merge original chr based files with Mutation Probabilities starts ##############################
        ####################################################################################################################

        log_out = open(log_file, 'a')
        print('#################################################################################', file=log_out)
        print('--- Merge original chr based files with Mutation Probabilities starts', file=log_out)
        print('#################################################################################', file=log_out)
        log_out.close()

        startSimNum = 0
        endSimNum = 0
        start_time = time.time()
        # SBS
        for mutation_type_context in mutation_types_contexts:
            # if (mutation_type_context in SBS_CONTEXTS) and (sbs_probabilities is not None):
            if (mutation_type_context in SBS_CONTEXTS):
                mutation_type_context_for_probabilities = get_mutation_type_context_for_probabilities_file(mutation_types_contexts_for_signature_probabilities,SUBS)
                log_out = open(log_file, 'a')
                print('--- Merge %s context mutations with probabilities for %s' % (mutation_type_context, sbs_probabilities), file=log_out)
                log_out.close()
                ordered_all_sbs_signatures_wrt_probabilities_file_array = prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList,
                                                                                   inputDir,
                                                                                   outputDir,
                                                                                   jobname,
                                                                                   mutation_type_context,
                                                                                   sbs_probabilities,
                                                                                   mutation_type_context_for_probabilities,
                                                                                   startSimNum,
                                                                                   endSimNum,
                                                                                   SNV,
                                                                                   PCAWG,
                                                                                   log_file,
                                                                                   verbose)

        # ID
        # if ((ID in mutation_types_contexts) and (id_probabilities is not None)):
        if (ID in mutation_types_contexts):
            mutation_type_context_for_probabilities = get_mutation_type_context_for_probabilities_file(mutation_types_contexts_for_signature_probabilities, INDELS)

            log_out = open(log_file, 'a')
            print('--- Merge %s mutations with probabilities for %s' % (ID, id_probabilities), file=log_out)
            log_out.close()

            ordered_all_id_signatures_wrt_probabilities_file_array = prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList,
                                                                                inputDir,
                                                                                outputDir,
                                                                                jobname,
                                                                                ID,
                                                                                id_probabilities,
                                                                                mutation_type_context_for_probabilities,
                                                                                startSimNum,
                                                                                endSimNum,
                                                                                ID,
                                                                                PCAWG,
                                                                                log_file,
                                                                                verbose)

        # DBS
        # if ((DBS in mutation_types_contexts) and (dbs_probabilities is not None)):
        if (DBS in mutation_types_contexts):
            mutation_type_context_for_probabilities = get_mutation_type_context_for_probabilities_file(mutation_types_contexts_for_signature_probabilities, DINUCS)

            log_out = open(log_file, 'a')
            print('--- Merge %s mutations with probabilities for %s' % (DBS, dbs_probabilities), file=log_out)
            log_out.close()

            ordered_all_dbs_signatures_wrt_probabilities_file_array = prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList,
                                                                                inputDir,
                                                                                outputDir,
                                                                                jobname,
                                                                                DBS,
                                                                                dbs_probabilities,
                                                                                mutation_type_context_for_probabilities,
                                                                                startSimNum,
                                                                                endSimNum,
                                                                                DBS,
                                                                                PCAWG,
                                                                                log_file,
                                                                                verbose)


        log_out = open(log_file, 'a')
        print("--- Merge original chr based files with Mutation Probabilities: %s seconds" % (time.time() - start_time), file=log_out)
        print("--- Merge original chr based files with Mutation Probabilities: %f minutes" % (float((time.time() - start_time) / 60)), file=log_out)
        print('--- Merge original chr based files with Mutation Probabilities ends', file=log_out)
        print('#################################################################################\n', file=log_out)
        log_out.close()

        ####################################################################################################################
        ##################  Merge original chr based files with Mutation Probabilities ends ################################
        ####################################################################################################################

        ####################################################################################################################
        ##################  Merge simulations chr based files with Mutation Probabilities starts ###########################
        ####################################################################################################################
        if (numofSimulations > 0):

            log_out = open(log_file, 'a')
            print('#################################################################################', file=log_out)
            print('--- Merge simulations chr based files with Mutation Probabilities starts', file=log_out)
            print('#################################################################################', file=log_out)
            log_out.close()

            startSimNum = 1
            endSimNum = numofSimulations
            start_time = time.time()
            # SBS
            for mutation_type_context in mutation_types_contexts:
                # if (mutation_type_context in SBS_CONTEXTS) and (sbs_probabilities is not None):
                if (mutation_type_context in SBS_CONTEXTS):

                    log_out = open(log_file, 'a')
                    mutation_type_context_for_probabilities = get_mutation_type_context_for_probabilities_file(mutation_types_contexts_for_signature_probabilities, SUBS)
                    print('--- Merge %s mutations with probabilities for %s' %(mutation_type_context,sbs_probabilities), file=log_out)
                    log_out.close()

                    prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList,
                                                                                       inputDir,
                                                                                       outputDir,
                                                                                       jobname,
                                                                                       mutation_type_context,
                                                                                       sbs_probabilities,
                                                                                       mutation_type_context_for_probabilities,
                                                                                       startSimNum,
                                                                                       endSimNum,
                                                                                       SNV,
                                                                                       PCAWG,
                                                                                       log_file,
                                                                                       verbose)

            # ID
            # if ((ID in mutation_types_contexts) and (id_probabilities is not None)):
            if (ID in mutation_types_contexts):
                mutation_type_context_for_probabilities = get_mutation_type_context_for_probabilities_file(mutation_types_contexts_for_signature_probabilities, ID)

                log_out = open(log_file, 'a')
                print('--- Merge %s mutations with probabilities for %s' % (ID, id_probabilities), file=log_out)
                log_out.close()

                prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList,
                                                                                   inputDir,
                                                                                   outputDir,
                                                                                   jobname,
                                                                                   ID,
                                                                                   id_probabilities,
                                                                                   mutation_type_context_for_probabilities,
                                                                                   startSimNum,
                                                                                   endSimNum,
                                                                                   ID,
                                                                                   PCAWG,
                                                                                   log_file,
                                                                                   verbose)

            # DBS
            # if ((DBS in mutation_types_contexts) and (dbs_probabilities is not None)):
            if (DBS in mutation_types_contexts):
                mutation_type_context_for_probabilities = get_mutation_type_context_for_probabilities_file(mutation_types_contexts_for_signature_probabilities, DBS)

                log_out = open(log_file, 'a')
                print('--- Merge %s mutations with probabilities for %s' % (DBS,dbs_probabilities), file=log_out)
                log_out.close()

                prepareMutationsDataAfterMatrixGenerationAndExtractorForTopography(chromShortNamesList,
                                                                                   inputDir,
                                                                                   outputDir,
                                                                                   jobname,
                                                                                   DBS,
                                                                                   dbs_probabilities,
                                                                                   mutation_type_context_for_probabilities,
                                                                                   startSimNum,
                                                                                   endSimNum,
                                                                                   DBS,
                                                                                   PCAWG,
                                                                                   log_file,
                                                                                   verbose)

            log_out = open(log_file, 'a')
            print("--- Merge simulations chr based files with Mutation Probabilities: %s seconds" %(time.time()-start_time), file=log_out)
            print("--- Merge simulations chr based files with Mutation Probabilities: %f minutes" %(float((time.time()-start_time)/60)), file=log_out)
            print('--- Merge simulations chr based files with Mutation Probabilities ends', file=log_out)
            print('#################################################################################\n', file=log_out)
            log_out.close()
        ####################################################################################################################
        ##################  Merge simulations chr based files with Mutation Probabilities ends #############################
        ####################################################################################################################

    else:
        log_out = open(log_file, 'a')
        for mutation_type_context in mutation_types_contexts:
            if (mutation_type_context in SBS_CONTEXTS):
                if ((sbs_probabilities is not None) and (os.path.exists(sbs_probabilities))):
                    ordered_all_sbs_signatures_wrt_probabilities_file_array = pd.read_csv(sbs_probabilities, sep='\t', nrows=0).columns.values
                else:
                    filename = '%s_%s_for_topography.txt' % ('chr1', SUBS)
                    chrBasedMutationDFFilePath = os.path.join(outputDir, jobname, DATA, CHRBASED, filename)
                    if os.path.exists(chrBasedMutationDFFilePath):
                        ordered_all_sbs_signatures_wrt_probabilities_file_array = pd.read_csv(chrBasedMutationDFFilePath,sep='\t', nrows=0).columns.values
                        print('ordered_all_sbs_signatures_wrt_probabilities_file_array:%s' %(ordered_all_sbs_signatures_wrt_probabilities_file_array), file=log_out)
                    else:
                        print('There is a problem: ordered_all_sbs_signatures_wrt_probabilities_file_array is not filled.', file=log_out)

        if (DBS in mutation_types_contexts):
            if ((dbs_probabilities is not None) and (os.path.exists(dbs_probabilities))):
                ordered_all_dbs_signatures_wrt_probabilities_file_array = pd.read_csv(dbs_probabilities, sep='\t', nrows=0).columns.values
            else:
                filename = '%s_%s_for_topography.txt' % ('chr1', DINUCS)
                chrBasedMutationDFFilePath = os.path.join(outputDir, jobname, DATA, CHRBASED, filename)
                if os.path.exists(chrBasedMutationDFFilePath):
                    ordered_all_dbs_signatures_wrt_probabilities_file_array = pd.read_csv(chrBasedMutationDFFilePath, sep='\t', nrows=0).columns.values
                    print('ordered_all_dbs_signatures_wrt_probabilities_file_array:%s' %(ordered_all_dbs_signatures_wrt_probabilities_file_array), file=log_out)
                else:
                    print('There is a problem: ordered_all_dbs_signatures_wrt_probabilities_file_array is not filled.', file=log_out)

        if (ID in mutation_types_contexts):
            if ((id_probabilities is not None) and (os.path.exists(id_probabilities))):
                ordered_all_id_signatures_wrt_probabilities_file_array = pd.read_csv(id_probabilities,sep='\t', nrows=0).columns.values
            else:
                filename = '%s_%s_for_topography.txt' % ('chr1', INDELS)
                chrBasedMutationDFFilePath = os.path.join(outputDir, jobname, DATA, CHRBASED, filename)
                if os.path.exists(chrBasedMutationDFFilePath):
                    ordered_all_id_signatures_wrt_probabilities_file_array = pd.read_csv(chrBasedMutationDFFilePath, sep='\t', nrows=0).columns.values
                    print('ordered_all_id_signatures_wrt_probabilities_file_array:%s' %(ordered_all_id_signatures_wrt_probabilities_file_array), file=log_out)
                else:
                    print('There is a problem: ordered_all_id_signatures_wrt_probabilities_file_array is not filled.', file=log_out)
        log_out.close()
    ###################################################################################################################
    ########### Step# Merge chrom based matrix generator generated files with probabilities ends ######################
    ###################################################################################################################

    #######################################################################################################
    ################################### Step4 Fill Table Starts ###########################################
    #######################################################################################################
    # Step4 Initialize these dataframes as empty dataframe
    # Step4 We will fill these dataframes if there is the corresponding data
    subsSignature_cutoff_numberofmutations_averageprobability_df = pd.DataFrame()
    dinucsSignature_cutoff_numberofmutations_averageprobability_df = pd.DataFrame()
    indelsSignature_cutoff_numberofmutations_averageprobability_df = pd.DataFrame()

    # Fill these pandas dataframes
    # cancer_type signature  number_of_mutations average_probability samples_list len(samples_list) len(all_samples_list) percentage_of_samples
    sbs_signature_number_of_mutations_df = pd.DataFrame()
    dbs_signature_number_of_mutations_df = pd.DataFrame()
    id_signature_number_of_mutations_df = pd.DataFrame()

    mutationtype_numberofmutations_numberofsamples_sampleslist_df = pd.DataFrame()
    chrlong_numberofmutations_df = pd.DataFrame()

    if (step4_tables):

        log_out = open(log_file, 'a')
        print('#################################################################################', file=log_out)
        print('--- Fill tables/dictionaries using original data starts', file=log_out)
        start_time = time.time()
        log_out.close()

        # For each signature we will find a cutoff value for mutations with average probability >=0.9
        # Our aim is to have at most 10% false positive rate in mutations
        # number of mutations >= 5K for subs signatures
        # number of mutations >= 1K for indels signatures
        # number of mutations >= 200 for dinuc signatures
        # If we can not satisfy this condition we will discard the signature

        cutoffs = []
        for cufoff in np.arange(0.5, 0.91, 0.01):
            cutoffs.append("%.2f" % (cufoff))

        # Initialize
        # mutationType2PropertiesListDict: PropertiesList consists of [NumberofMutations NumberofSamples SamplesList]
        mutationType2PropertiesDict = {}
        chrLong2NumberofMutationsDict = {}

        for mutation_type_context in mutation_types_contexts:
            if (mutation_type_context in SBS_CONTEXTS):
                if discreet_mode:
                    # We are reading original data to fill the signature2PropertiesListDict
                    # We are writing all samples_mutations_cutoffs_tables and signature based decided samples_mutations_cutoffs_tables in table format.
                    subsSignature_cutoff_numberofmutations_averageprobability_df = fillCutoff2Signature2PropertiesListDictionary(
                        outputDir,
                        jobname,
                        chromNamesList,
                        SUBS,
                        cutoffs,
                        average_probability,
                        num_of_sbs_required,
                        num_of_id_required,
                        num_of_dbs_required,
                        mutationType2PropertiesDict,
                        chrLong2NumberofMutationsDict,
                        log_file)
                else:
                    # Probability Mode
                    subsSignature_cutoff_numberofmutations_averageprobability_df = fill_signature_number_of_mutations_df(outputDir,
                                                                                                jobname,
                                                                                                chromNamesList,
                                                                                                SUBS,
                                                                                                default_cutoff,
                                                                                                num_of_sbs_required,
                                                                                                mutationType2PropertiesDict,
                                                                                                chrLong2NumberofMutationsDict,
                                                                                                show_all_signatures)

                subsSignature_cutoff_numberofmutations_averageprobability_df.to_csv(os.path.join(outputDir, jobname, DATA,
                                     Table_SBS_Signature_Cutoff_NumberofMutations_AverageProbability_Filename), sep='\t', index=False)

                subsSignature_cutoff_numberofmutations_averageprobability_df.drop(['cancer_type', 'samples_list',
                                                                                       'len(samples_list)',
                                                                                       'len(all_samples_list)',
                                                                                       'percentage_of_samples'], inplace=True, axis=1)

        if (DBS in mutation_types_contexts):
            if discreet_mode:
                # We are reading original data to fill the signature2PropertiesListDict
                # We are writing all samples_mutations_cutoffs_tables and signature based decided samples_mutations_cutoffs_tables in table format.
                dinucsSignature_cutoff_numberofmutations_averageprobability_df = fillCutoff2Signature2PropertiesListDictionary(
                    outputDir,
                    jobname,
                    chromNamesList,
                    DINUCS,
                    cutoffs,
                    average_probability,
                    num_of_sbs_required,
                    num_of_id_required,
                    num_of_dbs_required,
                    mutationType2PropertiesDict,
                    chrLong2NumberofMutationsDict,
                    log_file)
            else:
                # Probability Mode
                dinucsSignature_cutoff_numberofmutations_averageprobability_df = fill_signature_number_of_mutations_df(outputDir,
                                                                                            jobname,
                                                                                            chromNamesList,
                                                                                            DINUCS,
                                                                                            default_cutoff,
                                                                                            num_of_dbs_required,
                                                                                            mutationType2PropertiesDict,
                                                                                            chrLong2NumberofMutationsDict,
                                                                                            show_all_signatures)

            dinucsSignature_cutoff_numberofmutations_averageprobability_df.to_csv(
                    os.path.join(outputDir, jobname, DATA,
                                 Table_DBS_Signature_Cutoff_NumberofMutations_AverageProbability_Filename),
                    sep='\t', index=False)

            dinucsSignature_cutoff_numberofmutations_averageprobability_df.drop(
                    ['cancer_type', 'samples_list', 'len(samples_list)', 'len(all_samples_list)',
                     'percentage_of_samples'], inplace=True, axis=1)


        if (ID in mutation_types_contexts):
            if discreet_mode:
                # We are reading original data to fill the signature2PropertiesListDict
                # We are writing all samples_mutations_cutoffs_tables and signature based decided samples_mutations_cutoffs_tables in table format.
                indelsSignature_cutoff_numberofmutations_averageprobability_df = fillCutoff2Signature2PropertiesListDictionary(
                    outputDir,
                    jobname,
                    chromNamesList,
                    INDELS,
                    cutoffs,
                    average_probability,
                    num_of_sbs_required,
                    num_of_id_required,
                    num_of_dbs_required,
                    mutationType2PropertiesDict,
                    chrLong2NumberofMutationsDict,
                    log_file)
            else:
                # Probability Mode
                indelsSignature_cutoff_numberofmutations_averageprobability_df = fill_signature_number_of_mutations_df(outputDir,
                                                                                            jobname,
                                                                                            chromNamesList,
                                                                                            INDELS,
                                                                                            default_cutoff,
                                                                                            num_of_id_required,
                                                                                            mutationType2PropertiesDict,
                                                                                            chrLong2NumberofMutationsDict,
                                                                                            show_all_signatures)

            indelsSignature_cutoff_numberofmutations_averageprobability_df.to_csv(
                    os.path.join(outputDir, jobname, DATA,
                                 Table_ID_Signature_Cutoff_NumberofMutations_AverageProbability_Filename),
                    sep='\t', index=False)

            indelsSignature_cutoff_numberofmutations_averageprobability_df.drop(
                    ['cancer_type', 'samples_list', 'len(samples_list)', 'len(all_samples_list)',
                     'percentage_of_samples'], inplace=True, axis=1)

        # Add the last row
        numberofMutations = 0
        all_samples = set()

        for mutation_type in mutationType2PropertiesDict:
            numberofMutations += mutationType2PropertiesDict[mutation_type]['number_of_mutations']
            samples_list = mutationType2PropertiesDict[mutation_type]['samples_list']
            all_samples = all_samples.union(samples_list)

        all_samples_list=list(all_samples)
        all_samples_list = sorted(all_samples_list, key=natural_key)

        log_out = open(log_file, 'a')
        print("--- Number of samples: %d" %len(all_samples_list), file=log_out)
        print("--- Samples: %s" %(all_samples_list), file=log_out)
        all_samples_np_array=np.array(all_samples_list)

        mutationType2PropertiesDict['All']={}
        mutationType2PropertiesDict['All']['number_of_mutations'] = numberofMutations
        mutationType2PropertiesDict['All']['number_of_samples'] = len(all_samples)
        mutationType2PropertiesDict['All']['samples_list'] = all_samples_list

        # Write mutationType2PropertiesListDict dictionary as a dataframe starts
        filePath = os.path.join(outputDir, jobname, DATA, Table_MutationType_NumberofMutations_NumberofSamples_SamplesList_Filename)

        L = sorted([(mutation_type, a['number_of_mutations'], a['number_of_samples'], a['samples_list'])
                    for mutation_type, a  in mutationType2PropertiesDict.items()])

        if L:
            mutationtype_numberofmutations_numberofsamples_sampleslist_df = pd.DataFrame(L, columns=['mutation_type', 'number_of_mutations', 'number_of_samples', 'samples_list'])
        # write this dataframe
        mutationtype_numberofmutations_numberofsamples_sampleslist_df.to_csv(filePath, sep='\t', header=True, index=False)
        # Write dictionary as a dataframe ends

        # Write chrLong2NumberofMutationsDict dictionary as a dataframe
        if chrLong2NumberofMutationsDict:
            filePath = os.path.join(outputDir, jobname, DATA, Table_ChrLong_NumberofMutations_Filename)

            L = sorted([(chrLong, number_of_mutations)
                        for chrLong, number_of_mutations in chrLong2NumberofMutationsDict.items()])

            if L:
                chrlong_numberofmutations_df = pd.DataFrame(L, columns=['chrLong', 'number_of_mutations'])
            # write this dataframe
            chrlong_numberofmutations_df.to_csv(filePath, sep='\t', header=True, index=False)

        # We are reading original data again to fill the mutationType based, sample based and signature based dictionaries
        # This part is deprecated
        if sample_based:
            # Using original data
            for mutation_type_context in mutation_types_contexts:
                if (mutation_type_context in SBS_CONTEXTS):
                    fill_mutations_dictionaries_write(outputDir, jobname, chromNamesList, SUBS,
                                                      subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                      num_of_sbs_required, num_of_id_required,
                                                      num_of_dbs_required)
            if (DBS in mutation_types_contexts):
                fill_mutations_dictionaries_write(outputDir, jobname, chromNamesList, DINUCS,
                                                  dinucsSignature_cutoff_numberofmutations_averageprobability_df,
                                                  num_of_sbs_required,
                                                  num_of_id_required,
                                                  num_of_dbs_required)

            if (ID in mutation_types_contexts):
                fill_mutations_dictionaries_write(outputDir, jobname, chromNamesList, INDELS,
                                                  indelsSignature_cutoff_numberofmutations_averageprobability_df,
                                                  num_of_sbs_required,
                                                  num_of_id_required,
                                                  num_of_dbs_required)

        print("--- Fill tables/dictionaries using original data: %s seconds" % (time.time() - start_time), file=log_out)
        print("--- Fill tables/dictionaries using original data: %f minutes" % (float((time.time() - start_time) / 60)), file=log_out)
        print('--- Fill tables/dictionaries using original data ends', file=log_out)
        print('#################################################################################\n', file=log_out)
        log_out.close()

    else:
        mutationtype_numberofmutations_numberofsamples_sampleslist_df = pd.read_csv(os.path.join(outputDir,jobname,DATA,Table_MutationType_NumberofMutations_NumberofSamples_SamplesList_Filename),sep='\t', header=0, dtype={'mutation_type':str, 'number_of_mutations':np.int32})

        all_samples_string = mutationtype_numberofmutations_numberofsamples_sampleslist_df[mutationtype_numberofmutations_numberofsamples_sampleslist_df['mutation_type']=='All']['samples_list'].values[0]
        all_samples_list = eval(all_samples_string)
        all_samples_list = sorted(all_samples_list, key=natural_key)
        all_samples_np_array = np.array(all_samples_list)

        log_out = open(log_file, 'a')
        print('sample_based:%s --- len(all_samples_list):%d --- all_samples_list:%s' %(sample_based,len(all_samples_list), all_samples_list), file=log_out)
        log_out.close()

        chrlong_numberofmutations_df = pd.read_csv(os.path.join(outputDir, jobname, DATA, Table_ChrLong_NumberofMutations_Filename), sep='\t', header=0, dtype={'chrLong': str, 'number_of_mutations': np.int32})
        for mutation_type_context in mutation_types_contexts:
            if (mutation_type_context in SBS_CONTEXTS):
                subsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(os.path.join(outputDir, jobname, DATA, Table_SBS_Signature_Cutoff_NumberofMutations_AverageProbability_Filename), sep='\t', header=0, dtype={'cutoff':np.float32,'signature':str, 'number_of_mutations':np.int32,'average_probability':np.float32})
        if (DBS in mutation_types_contexts):
            dinucsSignature_cutoff_numberofmutations_averageprobability_df = pd.read_csv(os.path.join(outputDir, jobname, DATA, Table_DBS_Signature_Cutoff_NumberofMutations_AverageProbability_Filename), sep='\t',header=0, dtype={'cutoff': np.float32, 'signature': str, 'number_of_mutations': np.int32,'average_probability': np.float32})
        if (ID in mutation_types_contexts):
            indelsSignature_cutoff_numberofmutations_averageprobability_df= pd.read_csv(os.path.join(outputDir,jobname,DATA, Table_ID_Signature_Cutoff_NumberofMutations_AverageProbability_Filename), sep='\t', header=0, dtype={'cutoff':np.float32,'signature':str, 'number_of_mutations':np.int32,'average_probability':np.float32})
    #######################################################################################################
    ################################### Step4 Fill Table ends #############################################
    #######################################################################################################

    ###################################################################################################################
    ################################################# All Steps ends ##################################################
    ###################################################################################################################

    # Fill numpy arrays with the signatures in cutoff files
    sbs_signatures_with_cutoffs = np.array([])
    dbs_signatures_with_cutoffs = np.array([])
    id_signatures_with_cutoffs = np.array([])

    # Fill ordered_signatures arrays w.r.t the order in probabilities file
    # cutoffs_df (e.g.: subsSignature_cutoff_numberofmutations_averageprobability_df) are filled in (Step4=True or False but full_mode=True) or full_mode=False
    # ordered_signatures_wrt_probabilities_file are filled in (Step3=True or False but full_mode=True) or full_mode=False
    # We are interested in the signatures in cutoffs_df
    # But user might have changed the order of lines in cutoffs_df
    # Therefore we are setting the order in signatures_array and signatures_cutoff_arrays w.r.t. probabilities file
    ordered_sbs_signatures_with_cutoffs = np.array([])
    ordered_dbs_signatures_with_cutoffs = np.array([])
    ordered_id_signatures_with_cutoffs = np.array([])

    # Fill the list with the cutoff values
    # Fill ordered_signatures_cutoffs
    ordered_sbs_signatures_cutoffs = []
    ordered_dbs_signatures_cutoffs = []
    ordered_id_signatures_cutoffs = []

    if not subsSignature_cutoff_numberofmutations_averageprobability_df.empty:
        sbs_signatures_with_cutoffs = subsSignature_cutoff_numberofmutations_averageprobability_df['signature'].values

    if not dinucsSignature_cutoff_numberofmutations_averageprobability_df.empty:
        dbs_signatures_with_cutoffs = dinucsSignature_cutoff_numberofmutations_averageprobability_df['signature'].values

    if not indelsSignature_cutoff_numberofmutations_averageprobability_df.empty:
        id_signatures_with_cutoffs = indelsSignature_cutoff_numberofmutations_averageprobability_df['signature'].values

    if ordered_all_sbs_signatures_wrt_probabilities_file_array is not None:
        df_columns_subs_signatures_mask_array = np.isin(ordered_all_sbs_signatures_wrt_probabilities_file_array, sbs_signatures_with_cutoffs)
        ordered_sbs_signatures_with_cutoffs = ordered_all_sbs_signatures_wrt_probabilities_file_array[df_columns_subs_signatures_mask_array]
        for signature in ordered_sbs_signatures_with_cutoffs:
            cutoff = subsSignature_cutoff_numberofmutations_averageprobability_df[subsSignature_cutoff_numberofmutations_averageprobability_df['signature'] == signature]['cutoff'].values[0]
            ordered_sbs_signatures_cutoffs.append(cutoff)

    if ordered_all_dbs_signatures_wrt_probabilities_file_array is not None:
        df_columns_dbs_signatures_mask_array = np.isin(ordered_all_dbs_signatures_wrt_probabilities_file_array, dbs_signatures_with_cutoffs)
        ordered_dbs_signatures_with_cutoffs = ordered_all_dbs_signatures_wrt_probabilities_file_array[df_columns_dbs_signatures_mask_array]
        for signature in ordered_dbs_signatures_with_cutoffs:
            cutoff = dinucsSignature_cutoff_numberofmutations_averageprobability_df[dinucsSignature_cutoff_numberofmutations_averageprobability_df['signature'] == signature]['cutoff'].values[0]
            ordered_dbs_signatures_cutoffs.append(cutoff)

    if ordered_all_id_signatures_wrt_probabilities_file_array is not None:
        df_columns_id_signatures_mask_array = np.isin(ordered_all_id_signatures_wrt_probabilities_file_array, id_signatures_with_cutoffs)
        ordered_id_signatures_with_cutoffs = ordered_all_id_signatures_wrt_probabilities_file_array[df_columns_id_signatures_mask_array]
        for signature in ordered_id_signatures_with_cutoffs:
            cutoff = indelsSignature_cutoff_numberofmutations_averageprobability_df[indelsSignature_cutoff_numberofmutations_averageprobability_df['signature'] == signature]['cutoff'].values[0]
            ordered_id_signatures_cutoffs.append(cutoff)

    ordered_sbs_signatures_cutoffs = np.array(ordered_sbs_signatures_cutoffs)
    ordered_dbs_signatures_cutoffs = np.array(ordered_dbs_signatures_cutoffs)
    ordered_id_signatures_cutoffs = np.array(ordered_id_signatures_cutoffs)

    # Get all signatures ordered array w.r.t. the probabilities file
    ordered_all_sbs_signatures_array = get_all_signatures_array(ordered_all_sbs_signatures_wrt_probabilities_file_array, SBS)
    ordered_all_dbs_signatures_array = get_all_signatures_array(ordered_all_dbs_signatures_wrt_probabilities_file_array, DBS)
    ordered_all_id_signatures_array = get_all_signatures_array(ordered_all_id_signatures_wrt_probabilities_file_array, ID)

    log_out = open(log_file, 'a')
    print('--- discreet_mode:', discreet_mode, file=log_out)
    print('--- ordered_all_sbs_signatures_array:', ordered_all_sbs_signatures_array, file=log_out)
    print('--- ordered_all_dbs_signatures_array:', ordered_all_dbs_signatures_array, file=log_out)
    print('--- ordered_all_id_signatures_array:', ordered_all_id_signatures_array, file=log_out)
    print('--- ordered_sbs_signatures_with_cutoffs:', ordered_sbs_signatures_with_cutoffs, file=log_out)
    print('--- ordered_dbs_signatures_with_cutoffs:', ordered_dbs_signatures_with_cutoffs, file=log_out)
    print('--- ordered_id_signatures_with_cutoffs:', ordered_id_signatures_with_cutoffs, file=log_out)

    if discreet_mode:
        print('--- ordered_sbs_signatures_cutoffs:', ordered_sbs_signatures_cutoffs, file=log_out)
        print('--- ordered_dbs_signatures_cutoffs:', ordered_dbs_signatures_cutoffs, file=log_out)
        print('--- ordered_id_signatures_cutoffs:', ordered_id_signatures_cutoffs, file=log_out)

    ####################################################################################################################
    ################################### Run SigProfilerTopography Analysis starts ######################################
    ####################################################################################################################
    print('#################################################################################', file=log_out)
    print('--- Run SigProfilerTopography Analysis starts', file=log_out)
    log_out.close()

    if (computation_type == USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT):
        job_tuples = get_job_tuples(chrlong_numberofmutations_df,numofSimulations)

    else:
        job_tuples = []

    if (nucleosome):
        # Nucleosome Occupancy
        occupancy_type = NUCLEOSOMEOCCUPANCY

        if delete_old:
            deleteOldData(outputDir,jobname,occupancy_type)

        start_time = time.time()

        runOccupancyAnalyses(genome,
                             outputDir,
                             jobname,
                             numofSimulations,
                             job_tuples,
                             sample_based,
                             nucleosome_file,
                             None,
                             chromSizesDict,
                             chromNamesList,
                             ordered_sbs_signatures_with_cutoffs,
                             ordered_dbs_signatures_with_cutoffs,
                             ordered_id_signatures_with_cutoffs,
                             ordered_sbs_signatures_cutoffs,
                             ordered_dbs_signatures_cutoffs,
                             ordered_id_signatures_cutoffs,
                             computation_type,
                             occupancy_type,
                             occupancy_calculation_type,
                             plusorMinus_nucleosome,
                             remove_outliers,
                             quantileValue,
                             discreet_mode,
                             default_cutoff,
                             parallel_mode,
                             log_file,
                             verbose)

        log_out = open(log_file, 'a')
        print('#################################################################################', file=log_out)
        print("--- Run Nucleosome Occupancy Analysis: %s seconds --- %s" %((time.time()-start_time),nucleosome_file), file=log_out)
        print("--- Run Nucleosome Occupancy Analysis: %f minutes --- %s" %(float((time.time()-start_time)/60),nucleosome_file), file=log_out)
        print('#################################################################################\n', file=log_out)
        log_out.close()

    if (replication_time):
        # Replication Time
        # Required genome is already downloaded by matrix generator
        if delete_old:
            deleteOldData(outputDir,jobname,REPLICATIONTIME)

        start_time = time.time()

        runReplicationTimeAnalysis(genome,
                                   outputDir,
                                   jobname,
                                   numofSimulations,
                                   job_tuples,
                                   sample_based,
                                   replication_time_signal_file,
                                   chromSizesDict,
                                   chromNamesList,
                                   computation_type,
                                   ordered_sbs_signatures_with_cutoffs,
                                   ordered_dbs_signatures_with_cutoffs,
                                   ordered_id_signatures_with_cutoffs,
                                   ordered_sbs_signatures_cutoffs,
                                   ordered_dbs_signatures_cutoffs,
                                   ordered_id_signatures_cutoffs,
                                   discreet_mode,
                                   default_cutoff,
                                   parallel_mode,
                                   log_file,
                                   verbose,
                                   matrix_generator_path)

        log_out = open(log_file, 'a')
        print('#################################################################################', file=log_out)
        print("--- Run Replication Timing Analysis: %s seconds --- %s" %((time.time()-start_time),computation_type), file=log_out)
        print("--- Run Replication Timing Analysis: %f minutes --- %s" %(float((time.time()-start_time)/60),computation_type), file=log_out)
        print('#################################################################################\n', file=log_out)
        log_out.close()


    if replication_strand_bias:
        # Replication Strand Bias
        if delete_old:
            deleteOldData(outputDir,jobname,REPLICATIONSTRANDBIAS)

        start_time = time.time()

        runReplicationStrandBiasAnalysis(outputDir,
                                         jobname,
                                         numofSimulations,
                                         job_tuples,
                                         sample_based,
                                         all_samples_np_array,
                                         replication_time_signal_file,
                                         replication_time_valley_file,
                                         replication_time_peak_file,
                                         chromSizesDict,
                                         chromNamesList,
                                         computation_type,
                                         ordered_sbs_signatures_with_cutoffs,
                                         ordered_dbs_signatures_with_cutoffs,
                                         ordered_id_signatures_with_cutoffs,
                                         ordered_sbs_signatures_cutoffs,
                                         ordered_dbs_signatures_cutoffs,
                                         ordered_id_signatures_cutoffs,
                                         discreet_mode,
                                         default_cutoff,
                                         parallel_mode,
                                         log_file,
                                         verbose)

        log_out = open(log_file, 'a')
        print('#################################################################################', file=log_out)
        print("--- Run Replication Strand Asymmetry Analysis: %s seconds --- %s" %((time.time()-start_time),computation_type), file=log_out)
        print("--- Run Replication Strand Asymmetry Analysis: %f minutes --- %s" %(float((time.time()-start_time)/60),computation_type), file=log_out)
        print('#################################################################################\n', file=log_out)
        log_out.close()

    if transcription_strand_bias:
        # Transcription Strand Bias
        if delete_old:
            deleteOldData(outputDir,jobname,TRANSCRIPTIONSTRANDBIAS)

        start_time = time.time()
        runTranscriptionStradBiasAnalysis(outputDir,
                                          jobname,
                                          numofSimulations,
                                          job_tuples,
                                          sample_based,
                                          all_samples_np_array,
                                          chromNamesList,
                                          computation_type,
                                          ordered_sbs_signatures_with_cutoffs,
                                          ordered_dbs_signatures_with_cutoffs,
                                          ordered_id_signatures_with_cutoffs,
                                          ordered_sbs_signatures_cutoffs,
                                          ordered_dbs_signatures_cutoffs,
                                          ordered_id_signatures_cutoffs,
                                          discreet_mode,
                                          default_cutoff,
                                          parallel_mode,
                                          log_file,
                                          verbose)

        log_out = open(log_file, 'a')
        print('#################################################################################', file=log_out)
        print("--- Run Transcription Strand Asymmetry Analysis: %s seconds --- %s" %((time.time()-start_time),computation_type), file=log_out)
        print("--- Run Transcription Strand Asymmetry Analysis: %f minutes --- %s" %(float((time.time()-start_time)/60),computation_type), file=log_out)
        print('#################################################################################\n', file=log_out)
        log_out.close()

    if (processivity):
        # Processivity
        if delete_old:
            deleteOldData(outputDir,jobname,PROCESSIVITY)

        start_time = time.time()

        # we can use subsSignature_cutoff_numberofmutations_averageprobability_df
        # either filled w.r.t. discreet mode or prob_mode_default_cutoff
        runProcessivityAnalysis(mutation_types_contexts,
                                outputDir,
                                jobname,
                                numofSimulations,
                                chromNamesList,
                                processivity_calculation_type,
                                processivity_inter_mutational_distance,
                                subsSignature_cutoff_numberofmutations_averageprobability_df,
                                parallel_mode,
                                log_file,
                                verbose)

        log_out = open(log_file, 'a')
        print('#################################################################################', file=log_out)
        print("--- Run Strand-coordinated Mutagenesis Analysis: %s seconds ---" %(time.time()-start_time), file=log_out)
        print("--- Run Strand-coordinated Mutagenesis Analysis: %f minutes ---" %(float((time.time()-start_time)/60)), file=log_out)
        print('#################################################################################\n', file=log_out)
        log_out.close()

    if (epigenomics):
        # Epigenomics
        # If there is  a user provided name use it as occupancy_type
        if (epigenomics_dir_name is not None):
            occupancy_type = epigenomics_dir_name
        else:
            occupancy_type = EPIGENOMICSOCCUPANCY

        if delete_old:
            deleteOldData(outputDir,jobname,occupancy_type)

        # Run for each epigenomics file
        for idx, epigenomics_file in enumerate(epigenomics_files):
            start_time = time.time()
            if (epigenomics_files_memos is not None) and (len(epigenomics_files_memos) == len(epigenomics_files)):
                epigenomics_file_memo = epigenomics_files_memos[idx]
            else:
                epigenomics_file_memo = os.path.splitext(os.path.basename(epigenomics_file))[0]

            runOccupancyAnalyses(genome,
                                 outputDir,
                                 jobname,
                                 numofSimulations,
                                 job_tuples,
                                 sample_based,
                                 epigenomics_file,
                                 epigenomics_file_memo,
                                 chromSizesDict,
                                 chromNamesList,
                                 ordered_sbs_signatures_with_cutoffs,
                                 ordered_dbs_signatures_with_cutoffs,
                                 ordered_id_signatures_with_cutoffs,
                                 ordered_sbs_signatures_cutoffs,
                                 ordered_dbs_signatures_cutoffs,
                                 ordered_id_signatures_cutoffs,
                                 computation_type,
                                 occupancy_type,
                                 occupancy_calculation_type,
                                 plusorMinus_epigenomics,
                                 remove_outliers,
                                 quantileValue,
                                 discreet_mode,
                                 default_cutoff,
                                 parallel_mode,
                                 log_file,
                                 verbose)

            log_out = open(log_file, 'a')
            print('#################################################################################', file=log_out)
            print("--- Run Epigenomics Analysis: %s seconds --- %s" %((time.time()-start_time),epigenomics_file), file=log_out)
            print("--- Run Epigenomics Analyses: %f minutes --- %s" %(float((time.time()-start_time)/60),epigenomics_file), file=log_out)
            print('#################################################################################\n', file=log_out)
            log_out.close()

    log_out = open(log_file, 'a')
    print('--- Run SigProfilerTopography Analysis ends', file=log_out)
    print('#################################################################################\n', file=log_out)
    log_out.close()
    ####################################################################################################################
    ################################### Run SigProfilerTopography Analysis ends ########################################
    ####################################################################################################################

    if mutation_annotation_integration:
        # mutation_annotation_replication_timing_integration
        ordered_all_consequences = mutation_annotation_replication_timing_integration(inputDir,
                                                           outputDir,
                                                           jobname,
                                                           cancer_type=jobname)

        mutation_annotation_replication_timing_integration_signature_specific(inputDir,
                                                        outputDir,
                                                        jobname,
                                                        ordered_all_consequences,
                                                        subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                        cancer_type=jobname)

    ####################################################################################################################
    ############################################ Plot figures starts ###################################################
    ####################################################################################################################
    if (plot_figures):

        log_out = open(log_file, 'a')
        print('#################################################################################', file=log_out)
        print('--- Plot figures starts', file=log_out)
        log_out.close()

        start_time = time.time()
        plotFigures(outputDir,
                    jobname,
                    numofSimulations,
                    sample_based,
                    mutation_types_contexts,
                    epigenomics_files,
                    epigenomics_files_memos,
                    epigenomics_biosamples,
                    epigenomics_dna_elements,
                    epigenomics_dir_name,
                    nucleosome_file,
                    nucleosome_biosample,
                    epigenomics,
                    nucleosome,
                    replication_time,
                    replication_strand_bias,
                    transcription_strand_bias,
                    processivity,
                    plusorMinus_epigenomics,
                    plusorMinus_nucleosome,
                    epigenomics_heatmap_significance_level,
                    processivity_significance_level,
                    log_file,
                    verbose,
                    plot_epigenomics,
                    plot_nucleosome,
                    plot_replication_time,
                    plot_replication_strand_bias,
                    plot_transcription_strand_bias,
                    plot_processivity,
                    delete_old,
                    plot_mode,
                    combine_p_values_method,
                    fold_change_window_size,
                    num_of_real_data_avg_overlap)

        log_out = open(log_file, 'a')
        print('#################################################################################', file=log_out)
        print("--- Plot Figures: %s seconds ---" %(time.time()-start_time), file=log_out)
        print("--- Plot Figures: %f minutes ---" %(float((time.time()-start_time)/60)), file=log_out)
        print('--- Plot figures ends', file=log_out)
        print('#################################################################################\n', file=log_out)
        log_out.close()
    ####################################################################################################################
    ############################################ Plot figures ends #####################################################
    ####################################################################################################################

    log_out = open(log_file, 'a')
    print('#################################################################################', file=log_out)
    print("--- SigProfilerTopography ended successfully", file=log_out)
    print("--- Thanks for using SigProfilerTopography", file=log_out)
    print('#################################################################################\n', file=log_out)
    log_out.close()
    sys.stderr.close()


# Plot figures for the attainded data after SigProfilerTopography Analyses
def plotFigures(outputDir,
                jobname,
                numberofSimulations,
                sample_based,
                mutation_types_contexts,
                epigenomics_files,
                epigenomics_files_memos,
                epigenomics_biosamples,
                epigenomics_dna_elements,
                epigenomics_dir_name,
                nucleosome_file,
                nucleosome_biosample,
                epigenomics,
                nucleosome,
                replication_time,
                replication_strand_bias,
                transcription_strand_bias,
                processivity,
                plusOrMinus_epigenomics,
                plusOrMinus_nucleosome,
                epigenomics_heatmap_significance_level,
                processivity_significance_level,
                log_file,
                verbose,
                plot_epigenomics,
                plot_nucleosome,
                plot_replication_time,
                plot_replication_strand_bias,
                plot_transcription_strand_bias,
                plot_processivity,
                delete_old,
                plot_mode,
                combine_p_values_method,
                fold_change_window_size,
                num_of_real_data_avg_overlap):

    if (nucleosome or plot_nucleosome):
        occupancy_type = NUCLEOSOMEOCCUPANCY
        if delete_old:
            deleteOldFigures(outputDir, jobname, occupancy_type)
        nucleosome_file_basename = os.path.basename(nucleosome_file)
        occupancyAverageSignalFigures(outputDir,
                                      jobname,
                                      numberofSimulations,
                                      sample_based,
                                      mutation_types_contexts,
                                      nucleosome_file_basename,
                                      None,
                                      occupancy_type,
                                      plusOrMinus_nucleosome,
                                      log_file,
                                      verbose,
                                      plot_mode)

        log_out = open(log_file, 'a')
        print("--- Plot nucleosome occupancy ends", file=log_out)
        log_out.close()

    if (replication_time or plot_replication_time):
        if delete_old:
            deleteOldFigures(outputDir, jobname, REPLICATIONTIME)

        replicationTimeNormalizedMutationDensityFigures(outputDir,
                                                        jobname,
                                                        numberofSimulations,
                                                        sample_based,
                                                        mutation_types_contexts,
                                                        plot_mode)
        log_out = open(log_file, 'a')
        print("--- Plot replication time ends", file=log_out)
        log_out.close()

    if ((replication_strand_bias and transcription_strand_bias) or (plot_replication_strand_bias and plot_transcription_strand_bias)):
        if delete_old:
            deleteOldFigures(outputDir, jobname, STRANDBIAS)
        strand_bias_list = [TRANSCRIBED_VERSUS_UNTRANSCRIBED,GENIC_VERSUS_INTERGENIC,LAGGING_VERSUS_LEADING]
        transcriptionReplicationStrandBiasFiguresUsingDataframes(outputDir, jobname, numberofSimulations,
                                                                 mutation_types_contexts, strand_bias_list, plot_mode,
                                                                 log_file)

        log_out = open(log_file, 'a')
        print("--- Plot strand bias ends", file=log_out)
        log_out.close()

    elif (replication_strand_bias or plot_replication_strand_bias):
        strand_bias_list = [LAGGING_VERSUS_LEADING]
        transcriptionReplicationStrandBiasFiguresUsingDataframes(outputDir, jobname, numberofSimulations,
                                                                 mutation_types_contexts, strand_bias_list, plot_mode,
                                                                 log_file)

        log_out = open(log_file, 'a')
        print("--- Plot strand bias ends", file=log_out)
        log_out.close()

    elif (transcription_strand_bias or plot_transcription_strand_bias):
        strand_bias_list = [TRANSCRIBED_VERSUS_UNTRANSCRIBED, GENIC_VERSUS_INTERGENIC]
        transcriptionReplicationStrandBiasFiguresUsingDataframes(outputDir, jobname, numberofSimulations,
                                                                 mutation_types_contexts, strand_bias_list, plot_mode,
                                                                 log_file)

        log_out = open(log_file, 'a')
        print("--- Plot strand bias ends", file=log_out)
        log_out.close()

    if (processivity or plot_processivity):
        if delete_old:
            deleteOldFigures(outputDir, jobname, PROCESSIVITY)

        processivityFigures(outputDir,
                            jobname,
                            numberofSimulations,
                            processivity_significance_level,
                            log_file,
                            verbose)

        log_out = open(log_file, 'a')
        print("--- Plot processivity ends", file=log_out)
        log_out.close()

    if (epigenomics or plot_epigenomics):
        if epigenomics_dir_name is not None:
            occupancy_type=epigenomics_dir_name
        else:
            occupancy_type=EPIGENOMICSOCCUPANCY

        if delete_old:
            deleteOldFigures(outputDir, jobname, occupancy_type)

        # Initiate the pool
        numofProcesses = multiprocessing.cpu_count()

        # For real runs uncomment
        pool = multiprocessing.Pool(numofProcesses)
        jobs=[]

        # Please note that epigenomics_file_memo is not None
        # If None then it is created from filename.
        for idx, epigenomics_file in enumerate(epigenomics_files):
            epigenomics_file_basename = os.path.basename(epigenomics_file)
            epigenomics_file_memo= epigenomics_files_memos[idx]
            jobs.append(pool.apply_async(occupancyAverageSignalFigures,
                                         args=(outputDir,
                                               jobname,
                                               numberofSimulations,
                                               sample_based,
                                               mutation_types_contexts,
                                               epigenomics_file_basename,
                                               epigenomics_file_memo,
                                               occupancy_type,
                                               plusOrMinus_epigenomics,
                                               log_file,
                                               verbose,
                                               plot_mode,)))


        # Wait for all jobs to finish
        for job in jobs:
            log_out = open(log_file, 'a')
            if verbose: print('\n\tVerbose %s Worker pid %s Plotting figures  job.get():%s ' %(occupancy_type,str(os.getpid()),job.get()), file=log_out)
            log_out.close()

        pool.close()
        pool.join()

        log_out = open(log_file, 'a')
        print("--- Plot epigenomics occupancy ends", file=log_out)
        log_out.close()

        # # sequential call for testing or debugging
        # occupancyAverageSignalFigures(outputDir,
        #                               jobname,
        #                               numberofSimulations,
        #                               sample_based,
        #                               mutation_types_contexts,
        #                               epigenomics_file_basename,
        #                               epigenomics_file_memo,
        #                               occupancy_type,
        #                               plusOrMinus_epigenomics,
        #                               verbose,
        #                               plot_mode)

        compute_fold_change_with_p_values_plot_heatmaps(combine_p_values_method,
                                              fold_change_window_size,
                                              num_of_real_data_avg_overlap,
                                              outputDir,
                                              jobname,
                                              numberofSimulations,
                                              mutation_types_contexts,
                                              nucleosome_file,
                                              nucleosome_biosample,
                                              epigenomics_files_memos,
                                              epigenomics_biosamples,
                                              epigenomics_dna_elements,
                                              plusOrMinus_epigenomics,
                                              plusOrMinus_nucleosome,
                                              epigenomics_heatmap_significance_level,
                                              log_file,
                                              verbose)

        log_out = open(log_file, 'a')
        print("--- Plot epigenomics heatmaps ends", file=log_out)
        log_out.close()



# To run on local machine
if __name__== "__main__":

    genome = 'GRCh37'
    jobname = 'Test-Skin-Melanoma'
    numberofSimulations = 2
    inputDir = '/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/input/PCAWG_Matlab_Clean/Skin-Melanoma/filtered/'
    outputDir = os.path.join('C:\\','Users','burcak','Developer','Python','SigProfilerTopography','SigProfilerTopography','output_test')
    sbs_probabilities_file_path = os.path.join('C:\\','Users','burcak','Developer','Python','SigProfilerTopography','SigProfilerTopography','output_for_extractor','PCAWG_Matlab','Skin-Melanoma_sbs96_mutation_probabilities.txt')
    id_probabilities_file_path = os.path.join('C:\\','Users','burcak','Developer','Python','SigProfilerTopography','SigProfilerTopography','output_for_extractor','PCAWG_Matlab','Skin-Melanoma_id83_mutation_probabilities.txt')
    dbs_probabilities_file_path = os.path.join('C:\\','Users','burcak','Developer','Python','SigProfilerTopography','SigProfilerTopography','output_for_extractor','PCAWG_Matlab','Skin-Melanoma_dbs_mutation_probabilities.txt')
    user_provided_nucleosome_file_path = os.path.join('C:\\', 'Users', 'burcak', 'Developer', 'Python','SigProfilerTopography', 'SigProfilerTopography', 'lib','nucleosome', 'wgEncodeSydhNsomeGm12878Sig.wig')

    runAnalyses(genome,
                inputDir,
                outputDir,
                jobname,
                numberofSimulations,
                sbs_probabilities=sbs_probabilities_file_path,
                id_probabilities=id_probabilities_file_path,
                dbs_probabilities=dbs_probabilities_file_path,
                epigenomics=True,
                nucleosome=False,
                replication_time=False,
                strand_bias=False,
                processivity=False,
                sample_based=False,
                verbose=False)

