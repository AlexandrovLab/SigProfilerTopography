#This python file enables epigenomics analysis of Mutational Signatures
#TODO1 we may try same HM using bed file and wig file and see whether it differs.
#TODO2 we may allow user to provide multiple epigenomics data such as more than one histone modifications and transcription factors
#However we have to plot each HM or TF separately in the figure since each one has a different meaning.

import os
from TopographyCommons import *


######################################################################
#BED files: end is not inclusive
#BED files: score - A score between 0 and 1000. See track lines, below, for ways to configure the display style of scored data.
def updateSignalArrays(data_row,signalArray):
    signalArray[data_row[start]:data_row[end]] += data_row[score]
######################################################################


######################################################################
# This is used right now.
def writeChrBasedEpigenomicsSignalArray(inputList):
    print('In write for each chromosome function')
    memory_usage()

    chrLong = inputList[0]
    chromSize = inputList[1]
    chrBasedEpigenomicsDF = inputList[2]
    epigenomicsFilename = inputList[3]

    epigenomicsFilename_wo_dir=os.path.basename(epigenomicsFilename)
    epigenomicsFilename_wo_dir_wo_extension = os.path.splitext(epigenomicsFilename_wo_dir)[0]

    print('For debug chrLong:%s chrSize:%d' %(chrLong,chromSize))
    os.makedirs(os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,EPIGENOMICS,CHRBASED), exist_ok=True)

    print('writeChrBasedEpigenomics:%s for %s starts' %(epigenomicsFilename,chrLong))

    #TODO We may only need 1 or 0  or very little number greater than 1
    signalArray = np.zeros(chromSize,dtype=np.float32)

    print('For %s signal array before sum' %(chrLong))
    print(np.sum(signalArray))

    print('Sum of score column in chrBasedEpigenomicsDF')
    print(chrBasedEpigenomicsDF[score].sum())

    chrBasedEpigenomicsDF.apply(updateSignalArrays, signalArray=signalArray, axis=1)

    print('For %s signal array after sum' % (chrLong))
    print(np.sum(signalArray))

    #############################  Save as npy starts ################################
    signalArrayFilename = '%s_signal_%s' %(chrLong,epigenomicsFilename_wo_dir_wo_extension)
    signalArrayFilenameText= '%s_signal_%s.txt' %(chrLong,epigenomicsFilename_wo_dir_wo_extension)
    chrBasedSignalEpigenomicsFile = os.path.join(current_abs_path,ONE_DIRECTORY_UP,ONE_DIRECTORY_UP,LIB,EPIGENOMICS,CHRBASED,signalArrayFilename)
    chrBasedSignalEpigenomicsFileText = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, LIB, EPIGENOMICS,CHRBASED,signalArrayFilenameText)

    print('Debug August 29, 2019 -- current_abs_path:%s' %current_abs_path)
    print('Debug August 29, 2019 -- chrBasedSignalEpigenomicsFile:%s' %chrBasedSignalEpigenomicsFile)
    np.save(chrBasedSignalEpigenomicsFile, signalArray)
    np.savetxt(chrBasedSignalEpigenomicsFileText, signalArray, delimiter='\t', fmt='%s')
    #############################  Save as npy ends ##################################

    print('writeChrBasedEpigenomics:%s for %s ends' % (epigenomicsFilename, chrLong))
######################################################################


######################################################################
#TODO Also try with wig file, do you get the same at the end.
#TODO Read with pybedtools, convert to pandas datafamre what is the memory usage?
#Make score type and signal array type same to get rid of  accumulation error
#Please notice that this is for bed file
def readEpigenomicsData(epigenomicsFilename):
    column_names = [chrom, start, end, name, score, strand]
    if os.path.exists(epigenomicsFilename):
        #TODO  update it
        epigenomics_df = pd.read_table(epigenomicsFilename, nrows=1)  # 2.25 GB
        ncols=epigenomics_df.shape[1]

        if (ncols<=3):
            print('There is no enough columns in this bed file')
        elif (ncols==4):
            print('SigProfilerTopogarphy assumes that score column is in the 4th column of this bed file and there is no header')
            epigenomics_df=pd.read_table('GSM1127065_UCSF-UBC.Breast_Fibroblast_Primary_Cells.H3K4me1.RM071.bed',header=None, usecols=[0, 1, 2, 3],names = [chrom,start,end,score],dtype={0: 'category', 1: np.int32, 2: np.int32, 3: np.float32})
        elif (ncols>=5):
            print('SigProfilerTopogarphy assumes that score column is in the 5th column of this bed file and there is no header')
            epigenomics_df=pd.read_table('GSM1127065_UCSF-UBC.Breast_Fibroblast_Primary_Cells.H3K4me1.RM071.bed',header=None, usecols=[0, 1, 2, 4], names = [chrom,start,end,score], dtype={0: 'category', 1: np.int32, 2: np.int32, 4: np.float32})

        epigenomics_df = pd.read_table(epigenomicsFilename, sep="\t", header=None, names=column_names, dtype={chrom: 'category', start: np.int32, end: np.int32, name:str, score:np.float32, strand: 'category'}) # 2.25 GB

        # print("epigenomics_df.dtypes")
        # print(epigenomics_df.dtypes)

        # print('epigenomics_df.columns.values')
        # print(epigenomics_df.columns.values)

        #drop name column
        # epigenomics_df.drop([name], inplace=True, axis=1) # 2.12 GB

    return epigenomics_df
######################################################################

######################################################################
def readAllEpigenomicsDataAndWriteChrBasedSignalArraysInParallel(genome, epigenomicsFilename):
    chromSizesDict = getChromSizesDict(genome)

    print('Before epigenomics data is loaded into memory')
    memory_usage()

    method="Using pandas read table"

    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)

    # column_names = [chrom, start, end, signal]
    if os.path.exists(epigenomicsFilename):
        epigenomics_df = readEpigenomicsData(epigenomicsFilename)

        print('epigenomics_df.head()')
        print(epigenomics_df.head())

        mem_usage = epigenomics_df.memory_usage(deep=True) / (1024 ** 2) # convert bytes to megabytes
        print('epigenomics_df.memory_usage(deep=True) / (1024 ** 2)')
        print(mem_usage)
        print('sys.getsizeof(epigenomics_df)/(1024*1024)')
        print(sys.getsizeof(epigenomics_df)/(1024*1024))

        epigenetics_df_grouped = epigenomics_df.groupby(chrom)

        poolInputList = []

        for chrLong, chromBasedEpigenomicsDF in epigenetics_df_grouped:
            chromSize = chromSizesDict[chrLong]
            inputList = []
            inputList.append(chrLong)
            inputList.append(chromSize)
            inputList.append(chromBasedEpigenomicsDF)
            inputList.append(epigenomicsFilename)
            poolInputList.append(inputList)

        pool.map(writeChrBasedEpigenomicsSignalArray, poolInputList)

        print('After all chr based files are written')
        memory_usage()
######################################################################


######################################################################
def readAllEpigenomicsDataAndWriteChrBasedSignalArraysSequentially(genome, epigenomicsFilename):
    chromSizesDict = getChromSizesDict(genome)

    print('Before epigenomics data is loaded into memory')
    memory_usage()

    method="Using pandas read table"

    # column_names = [chrom, start, end, signal]
    if os.path.exists(epigenomicsFilename):
        epigenomics_df = readEpigenomicsData(epigenomicsFilename)

        print('epigenomics_df.head()')
        print(epigenomics_df.head())

        mem_usage = epigenomics_df.memory_usage(deep=True) / 1024 ** 2 # convert bytes to megabytes
        print(mem_usage)
        print('in MBs')

        epigenetics_df_grouped = epigenomics_df.groupby(chrom)

        for chrLong, chromBasedEpigenomicsDF in epigenetics_df_grouped:
            chromSize = chromSizesDict[chrLong]
            inputList = []
            inputList.append(chrLong)
            inputList.append(chromSize)
            inputList.append(chromBasedEpigenomicsDF)
            inputList.append(epigenomicsFilename)
            writeChrBasedEpigenomicsSignalArray(inputList)

        print('After all chr based files are written')
        memory_usage()
######################################################################
