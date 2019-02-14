# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time and strand bias.
# Copyright (C) 2018 Burcak Otlu
#
# SigProfilerTopography is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SigProfilerTopography is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SigProfilerTopography.  If not, see http://www.gnu.org/licenses/


#################################################################
#### This python code read the nucleosome occupancy filename ####
#### Remove the outliers or not #################################
#### Write the chromosome based nucleosome occupancy files ######
#################################################################


import multiprocessing
import os
import pandas as pd
import sys
import numpy as np

import matplotlib
BACKEND = 'Agg'
if matplotlib.get_backend().lower() != BACKEND.lower():
    # If backend is not set properly a call to describe will hang
    matplotlib.use(BACKEND)
from matplotlib import pyplot as plt

current_abs_path = os.path.abspath(os.path.dirname(__file__))
commonsPath = os.path.join(current_abs_path, '..','commons')
sys.path.append(commonsPath)

from TopographyCommons import OUTPUT
from TopographyCommons import NUCLEOSOMEOCCUPANCY
from TopographyCommons import FIGURE

from TopographyCommons import LIB
from TopographyCommons import NUCLEOSOME
from TopographyCommons import CHRBASED
from TopographyCommons import SIGNAL
from TopographyCommons import DISPLAY
from TopographyCommons import NODISPLAY


from TopographyCommons import ONE_DIRECTORY_UP
from TopographyCommons import ChrNamesInNucleosomesFilename

######################################################################
def  writeChrBasedNucleosome(chrNuclesomeList):
    chr = chrNuclesomeList[0]
    chr_based_nuclesome_df = chrNuclesomeList[1]
    name = chrNuclesomeList[2]

    current_abs_path = os.path.abspath(os.path.dirname(__file__))
    os.makedirs(os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,NUCLEOSOME,CHRBASED), exist_ok=True)

    print('writeChrBasedNucleosome:%s for %s starts' %(type(chr_based_nuclesome_df),chr))
    # write each chr based nucleosme to a text file

    #real file
    #wgEncodeSydhNsomeK562Sig.wig

    #test file
    #21 million rows
    #wgEncodeCrgMapabilityAlign100mer.wig

    #test file
    #411 rows
    #wgEncodeDacMapabilityConsensusExcludable.bed

    filename = '%s_%s' %(chr,name)

    partitionedNucleosmeFile = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,NUCLEOSOME,CHRBASED,filename)
    chr_based_nuclesome_df.to_csv(partitionedNucleosmeFile, header=None, index=None, sep='\t', mode='w')
    print('writeChrBasedNucleosome:%s for %s ends' % (type(chr_based_nuclesome_df), chr))
######################################################################



######################################################################
def readNucleosomeFile(filename,quantileValue,environment):
    #Read the file w.r.t. the current folder
    current_abs_path = os.path.abspath(os.path.dirname(__file__))

    #Real File: 595 million row file
    # 'wgEncodeSydhNsomeK562Sig.wig'
    #nucleosmeFile = os.path.join(current_abs_path, 'data',filename)
    # column_names = ['chrom', 'start', 'end', 'signal']

    #For Testing Purposes
    #21 million rows file
    # 'wgEncodeCrgMapabilityAlign100mer.wig'
    #nucleosmeFile = os.path.join(current_abs_path, 'data',filename)
    #column_names = ['chrom', 'start', 'end','signal']

    #For Testing Purposes
    #411 rows file
    #Small Bed File
    # 'wgEncodeDacMapabilityConsensusExcludable.bed'
    # column_names = ['chrom','start','end','info1','info2','info3']

    #e.g. sample lines for nucleosome occupancy
    # chr1    10180   10181   0.6
    # chr1    10181   10182   0.7
    # chr1    10182   10183   0.8
    # chr1    10183   10184   0.9
    # chr1    10184   10185   1.0
    # chr1    10185   10186   1.2
    # chr1    10186   10187   1.3

    nucleosmeFile = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,NUCLEOSOME,filename)
    column_names = ['chrom', 'start', 'end','signal']

    #original data signal is 0.8 float64 is enough for 0.
    # no need to make it 0.8000000000000000444 using float128
    #default float is float64 in python

    # nucleosome_df = pd.read_table(nucleosmeFile, sep="\t", header=None, comment='#', names=column_names, dtype={'chrom':str, 'start':np.int32, 'end':np.int32, 'signal':np.float32})
    nucleosome_df = pd.read_table(nucleosmeFile, sep="\t", header=None, comment='#', names=column_names,dtype={'chrom': str, 'start': np.int32, 'end': np.int32, 'signal': np.float64})
    # nucleosome_df = pd.read_table(nucleosmeFile, sep="\t", header=None, comment='#', names=column_names,dtype={'chrom': str, 'start': np.int32, 'end': np.int32, 'signal': np.float128})

    print('type(nucleosome_df):%s' % type(nucleosome_df))
    print(nucleosome_df.shape)
    print(nucleosome_df.shape[0])
    print(nucleosome_df.shape[1])
    print(nucleosome_df.head())

    current_abs_path = os.path.abspath(os.path.dirname(__file__))
    os.makedirs(os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,OUTPUT,FIGURE,NUCLEOSOMEOCCUPANCY), exist_ok=True)

    print('###### Describe Before #####')
    print(nucleosome_df[SIGNAL].describe())
    if (environment == DISPLAY):
        plt.hist(nucleosome_df[SIGNAL].astype(float), bins=100)
        #plt.show()
        figureName = 'OriginalNucleosomeOccupancyHistogram.png'
        figureFile = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,OUTPUT,FIGURE,NUCLEOSOMEOCCUPANCY,figureName)
        plt.savefig(figureFile)

    if (quantileValue<1.0):
        #remove the outliers
        q = nucleosome_df[SIGNAL].quantile(quantileValue)
        print('q:%f' % q)
        print('before %d' % (nucleosome_df.shape[0]))
        nucleosome_df = nucleosome_df[nucleosome_df[SIGNAL] < q]
        print('after %d' % (nucleosome_df.shape[0]))
        print('###### Describe After Outlier Removal #####')
        print(nucleosome_df[SIGNAL].describe())
        if (environment == DISPLAY):
            plt.hist(nucleosome_df[SIGNAL].astype(float), bins=100)
            # plt.show()
            figureName = 'AfterOutlierRemovalNucleosomeOccupancyHistogram.png'
            figureFile = os.path.join(current_abs_path,'..',OUTPUT,FIGURE,NUCLEOSOMEOCCUPANCY,figureName)
            plt.savefig(figureFile)

    ###################################################
    ######### for small bed file starts ###############
    ###################################################
    # # drop columns info1, info2, info3
    # nucleosome_df.drop(['info1', 'info2', 'info3'], axis=1, inplace=True)
    #
    # # for small bed file
    # # add signal column.
    # np.random.seed(0)
    # sLength =len(nucleosome_df['start'])
    # nucleosome_df['signal'] = np.random.randn(sLength)
    # print(nucleosome_df.head())
    #
    # print('type(nucleosome_df):%s' %type(nucleosome_df))
    # print(nucleosome_df.shape)
    # print(nucleosome_df.shape[0])
    # print(nucleosome_df.shape[1])
    # print(nucleosome_df.head())
    ###################################################
    ######### for small bed file ends #################
    ###################################################
    return nucleosome_df
######################################################################

######################################################################
# if __name__ == '__main__':
def partitionNucleosomeOccupancyData(nucleosomeFilename,quantileValue):

    # nucleosomeFilename = 'wgEncodeSydhNsomeGm12878Sig.wig'
    # nucleosomeFilename = 'wgEncodeCrgMapabilityAlign100mer.wig'
    # quantileValue = 0.97

    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)

    #############################################################################
    ########################## Read Nuclesome Data ##############################
    #############################################################################
    print('################# Nucleosome ##################')
    nucleosome_df = readNucleosomeFile(nucleosomeFilename,quantileValue,NODISPLAY)

    #########################################################
    ############### Write Unique Chrnames starts ############
    #########################################################
    uniqueChrNames = nucleosome_df['chrom'].unique()
    # The unique values are returned as a NumPy array

    # Write uniqueChrNames to a file
    ChrNamesFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME, ChrNamesInNucleosomesFilename)

    np.savetxt(ChrNamesFile,uniqueChrNames, delimiter='\t', fmt='%s')
    #########################################################
    ############### Write Unique Chrnames ends ##############
    #########################################################

    grouped_nucleosome_df = nucleosome_df.groupby('chrom')
    print('type(grouped_nucleosome_df):%s' %type(grouped_nucleosome_df))
    print(len(grouped_nucleosome_df))

    list = []

    for chr, chr_based_nuclesome_group in grouped_nucleosome_df:
        print(chr)
        chrNameDF = []
        chrNameDF.append(chr)
        chrNameDF.append(chr_based_nuclesome_group)
        chrNameDF.append(nucleosomeFilename)
        list.append(chrNameDF)

    pool.map(writeChrBasedNucleosome,list)
    pool.close()
    pool.join()
    #############################################################################
    #############################################################################
    #############################################################################

######################################################################


######################################################################
if __name__ == '__main__':
    partitionNucleosomeOccupancyData('wgEncodeCrgMapabilityAlign100mer.wig',0.97)
######################################################################
