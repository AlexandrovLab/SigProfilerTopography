# Depreceated
# Very slow becuase of unnesting
# Not used anymore

# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time and strand bias.
# Copyright (C) 2018 Burcak Otlu

###################################################################################
################ This python code read the nucleosome occupancy filename ##########
######################## Remove the outliers or not ###############################
###### Write the chromosome based average nucleosome occupancy signal arrays ######
###################################################################################


import os
import matplotlib
import sys
from sys import getsizeof

BACKEND = 'Agg'
if matplotlib.get_backend().lower() != BACKEND.lower():
    # If backend is not set properly a call to describe will hang
    matplotlib.use(BACKEND)

from matplotlib import pyplot as plt

#############################################################
current_abs_path = os.path.abspath(os.path.dirname(__file__))
#############################################################

commonsPath = os.path.join(current_abs_path, '..','commons')
sys.path.append(commonsPath)

from TopographyCommons import *

#########################################################################
def unnesting(df, explode):
    idx=df.index.repeat(df[explode[0]].str.len())
    df1=pd.concat([pd.DataFrame({x:np.concatenate(df[x].values)} )for x in explode],axis=1)
    df1.index=idx
    return df1.join(df.drop(explode,1),how='left')
#########################################################################



######################################################################
def computeSignalArrayAndCountArrayForEachSplit(chromSize,chrBased_nucleosome_split_df):
    print('for debug in computeSignalArrayAndCountArrayForEachSplit starts')

    print('before key column added')
    print('getsizeof(chrBased_nucleosome_split_df)')
    print(getsizeof(chrBased_nucleosome_split_df))

    chrBased_nucleosome_split_df['key'] = [list(range(x, y)) for x, y in zip(chrBased_nucleosome_split_df.start, chrBased_nucleosome_split_df.end)]

    print('after key column added')
    print('getsizeof(chrBased_nucleosome_split_df)')
    print(getsizeof(chrBased_nucleosome_split_df))

    chrBased_nucleosome_split_df.drop(['start','end'], axis=1, inplace=True)

    print('after drop column')
    print('getsizeof(chrBased_nucleosome_split_df)')
    print(getsizeof(chrBased_nucleosome_split_df))

    # Unnest and  generate signalArray and countArray
    chrBasedNucleosomeSplitDFSignalArray = unnesting(chrBased_nucleosome_split_df,['key']).groupby('key').signal.sum().reindex(list(range(chromSize))).fillna(0).values
    chrBasedNucleosomeSplitDFCountArray = unnesting(chrBased_nucleosome_split_df,['key']).groupby('key').signal.count().reindex(list(range(chromSize))).fillna(0).values
    print('getsizeof(chrBased_nucleosome_split_df): %d MB' %int(getsizeof(chrBased_nucleosome_split_df)/1000000))
    print('getsizeof(chrBasedNucleosomeSplitDFSignalArray): %d MB' %int(getsizeof(chrBasedNucleosomeSplitDFSignalArray)/1000000))
    print('getsizeof(chrBasedNucleosomeSplitDFCountArray): %d MB' %int(getsizeof(chrBasedNucleosomeSplitDFCountArray)/1000000))
    print('np.sum(chrBasedNucleosomeSplitDFSignalArray)')
    print(np.sum(chrBasedNucleosomeSplitDFSignalArray, dtype=np.float64))
    print('np.sum(chrBasedNucleosomeSplitDFCountArray)')
    print(np.sum(chrBasedNucleosomeSplitDFCountArray, dtype=np.float64))
    print('for debug in computeSignalArrayAndCountArrayForEachSplit ends')
    return chrBasedNucleosomeSplitDFSignalArray, chrBasedNucleosomeSplitDFCountArray
######################################################################


# ######################################################################
# def accumulateSignalArraysAndCountArraysAndComuputeAverageSignalArray(chromSize,listofSignalArrayAndCountArrayTuples):
#     accumulatedSignalArray = np.zeros(chromSize,dtype=np.float16)
#     accumulatedCountArray = np.zeros(chromSize, dtype=np.float16)
#
#     for tuple in listofSignalArrayAndCountArrayTuples:
#         accumulatedSignalArray += tuple[0]
#         accumulatedCountArray += tuple[1]
#
#     np.seterr(divide='ignore', invalid='ignore')
#
#     averageSignalArray = np.divide(accumulatedSignalArray, accumulatedCountArray, dtype=np.float16)
#     np.nan_to_num(averageSignalArray, copy=False)
#
#     return averageSignalArray
# ######################################################################


######################################################################
def  writeChrBasedNucleosomeOccupancyAverageSignalArray(chrLong, chromSize, chrBasedNuclesomeDF,nucleosomeFilename):

    numofProcesses = multiprocessing.cpu_count()

    print('for debug chrLong:%s chrSize:%d' %(chrLong,chromSize))
    os.makedirs(os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,NUCLEOSOME,CHRBASED), exist_ok=True)

    print('writeChrBasedNucleosome:%s for %s starts' %(nucleosomeFilename,chrLong))
    nucleosomeFilenameWoExtension = nucleosomeFilename[0:-4]
    filename = '%s_%s' %(chrLong,nucleosomeFilenameWoExtension)
    partitionedNucleosmeFile = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,NUCLEOSOME,CHRBASED,filename)

    chrBasedNucleosomeDFSplits = np.array_split(chrBasedNuclesomeDF,numofProcesses)

    accumulatedSignalArray = np.zeros(chromSize,dtype=np.float16)
    accumulatedCountArray = np.zeros(chromSize, dtype=np.float16)

    for idx, chrBasedNucleosomeDFSplit in enumerate(chrBasedNucleosomeDFSplits):
        print('for debug for split %d' %(idx))
        print('chromSize')
        print(chromSize)
        singalArray, countArray = computeSignalArrayAndCountArrayForEachSplit(chromSize,chrBasedNucleosomeDFSplit)
        accumulatedSignalArray += singalArray
        accumulatedCountArray += countArray

    np.seterr(divide='ignore', invalid='ignore')
    averageSignalArray = np.divide(accumulatedSignalArray, accumulatedCountArray, dtype=np.float16)
    np.nan_to_num(averageSignalArray, copy=False)

    np.save(partitionedNucleosmeFile, averageSignalArray)
    print('writeChrBasedNucleosome:%s for %s ends' % (nucleosomeFilename, chrLong))
######################################################################




######################################################################
def readChromBasedNucleosomeDF(chrLong,nucleosomeFilename):
    chrBasedNucleosmeFilename = '%s_%s' %(chrLong,nucleosomeFilename)
    chrBasedNucleosmeFilePath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME, CHRBASED, chrBasedNucleosmeFilename)

    column_names = ['chrom', 'start', 'end', 'signal']
    chrbased_nucleosome_df = pd.read_table(chrBasedNucleosmeFilePath, sep="\t", header=None, comment='#', names=column_names,dtype={'chrom': str, 'start': np.int32, 'end': np.int32, 'signal': np.float16})

    return chrbased_nucleosome_df
######################################################################



######################################################################
def readChrBasedNucleosomeOccupancyDataAndWriteAverageSignalArray(genomeassembly, nucleosomeFilename):

    # chrom2Size_dict = readChromSizes(HG19)
    if genomeassembly==HG19:
        chrom2Size_dict= {'chr1' : 249250621, 'chr2' : 243199373, 'chr3' : 198022430,'chr4': 191154276,'chr5': 180915260,'chr6':	171115067,'chr7':	159138663,'chrX':	155270560,'chr8':	146364022,'chr9':	141213431,
        'chr10':135534747, 'chr11':135006516, 'chr12':133851895, 'chr13':115169878, 'chr14':107349540, 'chr15':102531392, 'chr16':90354753, 'chr17':81195210, 'chr18':78077248, 'chr20':63025520, 'chrY':	59373566,
        'chr19':59128983, 'chr22':51304566, 'chr21':48129895, 'chrM' :16571}
    elif genomeassembly ==HG38:
        chrom2Size_dict = {'chr1':	248956422, 'chr2':	242193529, 'chr3':	198295559, 'chr4':	190214555, 'chr5':	181538259, 'chr6':	170805979, 'chr7':	159345973, 'chrX':	156040895, 'chr8':	145138636, 'chr9':	138394717,
                           'chr11':	135086622, 'chr10':	133797422, 'chr12':	133275309, 'chr13':	114364328, 'chr14':	107043718, 'chr15':	101991189, 'chr16':	90338345, 'chr17':	83257441, 'chr18':	80373285, 'chr20':	64444167,
                           'chr19':	58617616, 'chrY':	57227415, 'chr22':	50818468, 'chr21':	46709983, 'chrM':	16569}

    #read chromnames for this nucleosome data
    ChrNamesFilepath = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, NUCLEOSOME, ChrNamesInNucleosomesFilename)
    print('for debug ChrNamesFilepath: %s:' %(ChrNamesFilepath))

    if (os.path.exists(ChrNamesFilepath)):
        chromNames = np.loadtxt(ChrNamesFilepath,dtype=np.str)

    print('for debug type(chromNames): %s' %type(chromNames))

    for chrLong in chromNames:
        chromBasedNucleosomeDF = readChromBasedNucleosomeDF(chrLong,nucleosomeFilename)
        chromSize = chrom2Size_dict[chrLong]
        writeChrBasedNucleosomeOccupancyAverageSignalArray(chrLong, chromSize, chromBasedNucleosomeDF,nucleosomeFilename)
######################################################################

# ######################################################################
# if __name__== "__main__":
#     #Read chrom based nucleosome data and write the average signal array in npy format
#     # readChrBasedNucleosomeOccupancyDataAndWriteAverageSignalArray('hg19','wgEncodeSydhNsomeGm12878Sig.wig')
#     readChrBasedNucleosomeOccupancyDataAndWriteAverageSignalArray('hg19', 'wgEncodeCrgMapabilityAlign100mer.wig')
# #####################################################################
