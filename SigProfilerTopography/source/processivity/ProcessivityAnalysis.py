# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018 Burcak Otlu

import os
import sys
from json import JSONEncoder


#############################################################
current_abs_path = os.path.abspath(os.path.dirname(__file__))
print('ProcessivityAnalysis.py current_abs_path:%s' %(current_abs_path))
#############################################################


commonsPath = os.path.join(current_abs_path, '..','commons')
sys.path.append(commonsPath)

from SigProfilerTopography.source.commons.TopographyCommons import *

#########################################################################
class ProcessiveGroupStructEncoder(JSONEncoder):
    def default(self, o):
        return o.__dict__
#########################################################################

#########################################################################
class ProcessiveGroupStruct(JSONEncoder):
    def __init__(self, processiveGroupLength, numberofProcessiveGroups, medianofNumberofProcessiveGroupsinMB):
        self.processiveGroupLength = processiveGroupLength
        self.numberofProcessiveGroups = numberofProcessiveGroups
        self.medianofNumberofProcessiveGroupsinMB = medianofNumberofProcessiveGroupsinMB
#########################################################################

#########################################################################
def myfunc(x,signatures):
    return (x.shape[0], x['Start'].iloc[-1] - x['Start'].iloc[0], *x[signatures].agg(max))
#########################################################################

#########################################################################
def myfuncUpdated(x):
    return (x.shape[0], x['Start'].iloc[-1] - x['Start'].iloc[0],  x['Signature'].iloc[0], x['Probability'].agg(max))
#########################################################################


#########################################################################
def fillDict(x,signatures,signature2ProcessiveGroupLength2DistanceListDict):
    for signature in signatures:
        if signature in signature2ProcessiveGroupLength2DistanceListDict:
            processiveGroupLength = x['ProcessiveGroupLength']
            distance = x['Distance']
            if processiveGroupLength in signature2ProcessiveGroupLength2DistanceListDict[signature]:
                if (x[signature]>=0.5):
                    signature2ProcessiveGroupLength2DistanceListDict[signature][processiveGroupLength].append(distance)
            else:
                if (x[signature]>=0.5):
                    signature2ProcessiveGroupLength2DistanceListDict[signature][processiveGroupLength] = [distance]
        else:
            if (x[signature] >= 0.5):
                processiveGroupLength = x['ProcessiveGroupLength']
                distance = x['Distance']
                signature2ProcessiveGroupLength2DistanceListDict[signature]= {}
                signature2ProcessiveGroupLength2DistanceListDict[signature][processiveGroupLength] = [distance]
#########################################################################


#########################################################################
def fillDictUpdated(x,signature2ProcessiveGroupLength2DistanceListDict, considerProbability):

    processiveGroupLength = x['ProcessiveGroupLength']
    distance = x['Distance']
    signature = x['Signature']
    probability = x['Probability']

    if (considerProbability and probability>=0.5):
        if signature in signature2ProcessiveGroupLength2DistanceListDict:
            if processiveGroupLength in signature2ProcessiveGroupLength2DistanceListDict[signature]:
                signature2ProcessiveGroupLength2DistanceListDict[signature][processiveGroupLength].append(distance)
            else:
                signature2ProcessiveGroupLength2DistanceListDict[signature][processiveGroupLength] = [distance]
        else:
            signature2ProcessiveGroupLength2DistanceListDict[signature]= {}
            signature2ProcessiveGroupLength2DistanceListDict[signature][processiveGroupLength] = [distance]
    elif (not considerProbability):
        if signature in signature2ProcessiveGroupLength2DistanceListDict:
            if processiveGroupLength in signature2ProcessiveGroupLength2DistanceListDict[signature]:
                signature2ProcessiveGroupLength2DistanceListDict[signature][processiveGroupLength].append(distance)
            else:
                signature2ProcessiveGroupLength2DistanceListDict[signature][processiveGroupLength] = [distance]
        else:
            signature2ProcessiveGroupLength2DistanceListDict[signature] = {}
            signature2ProcessiveGroupLength2DistanceListDict[signature][processiveGroupLength] = [distance]
#########################################################################



####################################################################################
def findProcessiveGroupsForInputList(inputList):

    # chrLong = inputList[0]
    # sample =  inputList[1]
    # signatures = inputList[3]
    sampleBased_chrBased_spms_df = inputList[0]
    considerProbability = inputList[1]

    #Sort it
    sampleBased_chrBased_spms_df.sort_values('Start', inplace=True)

    return findProcessiveGroupsUpdated(sampleBased_chrBased_spms_df,considerProbability)
####################################################################################

####################################################################################
#Old Way
def findProcessiveGroups(chrLong,sorted_spms_df,signatures):
    signature2ProcessiveGroupLength2DistanceListDict = {}

    if (sorted_spms_df is not None):
        sorted_spms_df['subgroup'] = (sorted_spms_df['Mutation'] != sorted_spms_df['Mutation'].shift(1)).cumsum()
        series_new = sorted_spms_df.groupby('subgroup', as_index=False).apply(myfunc, signatures=signatures)
        columns = ['ProcessiveGroupLength', 'Distance']
        columns.extend(signatures)
        df = series_new.apply(pd.Series)
        df.columns = columns
        df = df[df['ProcessiveGroupLength'].ne(1)]
        signature2ProcessiveGroupLength2DistanceListDict = {}
        df.apply(fillDict, signatures=signatures,signature2ProcessiveGroupLength2DistanceListDict=signature2ProcessiveGroupLength2DistanceListDict, axis=1)

    return signature2ProcessiveGroupLength2DistanceListDict
####################################################################################


####################################################################################
#New Way
def findProcessiveGroupsUpdated(sorted_sampleBased_chrBased_spms_df,considerProbability):
    signature2ProcessiveGroupLength2DistanceListDict = {}

    if (sorted_sampleBased_chrBased_spms_df is not None):
        sorted_sampleBased_chrBased_spms_df['subgroup'] = ((sorted_sampleBased_chrBased_spms_df['Mutation'] != sorted_sampleBased_chrBased_spms_df['Mutation'].shift(1)) |
                                                           (sorted_sampleBased_chrBased_spms_df['Signature'] != sorted_sampleBased_chrBased_spms_df['Signature'].shift(1))) .cumsum()


        series_new = sorted_sampleBased_chrBased_spms_df.groupby('subgroup', as_index=False).apply(myfuncUpdated)
        columns = ['ProcessiveGroupLength', 'Distance', 'Signature', 'Probability']
        df = series_new.apply(pd.Series)
        df.columns = columns
        df = df[df['ProcessiveGroupLength'].ne(1)]

        df.apply(fillDictUpdated,signature2ProcessiveGroupLength2DistanceListDict=signature2ProcessiveGroupLength2DistanceListDict, considerProbability=considerProbability , axis=1)

    return signature2ProcessiveGroupLength2DistanceListDict
####################################################################################

####################################################################################
def findMedians(signature2ProcessiveGroupLength2GroupSizeListDict):
    signature2ProcessiveGroupLength2PropertiesDict = {}

    for signature in sorted(signature2ProcessiveGroupLength2GroupSizeListDict.keys()):
        signature2ProcessiveGroupLength2PropertiesDict[signature] = {}
        for processiveGroupLength in sorted(signature2ProcessiveGroupLength2GroupSizeListDict[signature].keys()):

            groupSizeArray = np.array(signature2ProcessiveGroupLength2GroupSizeListDict[signature][processiveGroupLength])
            #if we remove the 0.0 size processive groups we will get rid of infinity
            groupSizeArray = groupSizeArray[groupSizeArray != 0]
            if (groupSizeArray.size>0):
                groupSizeArray = (1/groupSizeArray) * processiveGroupLength * 1000000
                median = np.median(groupSizeArray)
                numberofProcessiveGroups = len(groupSizeArray)

                processiveGroupStruct = ProcessiveGroupStruct(processiveGroupLength=processiveGroupLength,
                                                 numberofProcessiveGroups=numberofProcessiveGroups,
                                                 medianofNumberofProcessiveGroupsinMB = median)

                signature2ProcessiveGroupLength2PropertiesDict[signature][processiveGroupLength]= processiveGroupStruct

    return signature2ProcessiveGroupLength2PropertiesDict
####################################################################################


####################################################################################
def accumulateListofDicts(allSamples_chrBased_signature2ProcessiveGroupLength2DistanceListDict_list):
    allSamples_chrBased_signature2ProcessiveGroupLength2DistanceListDict = {}
    for sample_chrBased_signature2ProcessiveGroupLength2DistanceListDict in allSamples_chrBased_signature2ProcessiveGroupLength2DistanceListDict_list:
        accumulateDict(sample_chrBased_signature2ProcessiveGroupLength2DistanceListDict,allSamples_chrBased_signature2ProcessiveGroupLength2DistanceListDict)
    return  allSamples_chrBased_signature2ProcessiveGroupLength2DistanceListDict
####################################################################################


####################################################################################
def accumulateDict(small_signature2ProcessiveGroupLength2DistanceListDict,big_signature2ProcessiveGroupLength2DistanceListDict):
    for signature in small_signature2ProcessiveGroupLength2DistanceListDict:
        for processiveGroupLength in small_signature2ProcessiveGroupLength2DistanceListDict[signature]:
            if signature in big_signature2ProcessiveGroupLength2DistanceListDict:
                if processiveGroupLength in big_signature2ProcessiveGroupLength2DistanceListDict[signature]:
                    big_signature2ProcessiveGroupLength2DistanceListDict[signature][processiveGroupLength].extend(small_signature2ProcessiveGroupLength2DistanceListDict[signature][processiveGroupLength])
                else:
                    big_signature2ProcessiveGroupLength2DistanceListDict[signature][processiveGroupLength] = small_signature2ProcessiveGroupLength2DistanceListDict[signature][processiveGroupLength]
            else:
                big_signature2ProcessiveGroupLength2DistanceListDict[signature]={}
                big_signature2ProcessiveGroupLength2DistanceListDict[signature][processiveGroupLength] = small_signature2ProcessiveGroupLength2DistanceListDict[signature][processiveGroupLength]
####################################################################################

####################################################################################
def readSinglePointMutationsFindProcessivityGroupsWithMultiProcessing(jobname, singlePointMutationsFileName,considerProbability):

    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)

    if (singlePointMutationsFileName != NOTSET):
        # Load the chrnames in single point mutations data
        ChrNamesFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, DATA,ChrNamesInSPMsFilename)
        if (os.path.exists(ChrNamesFile)):
            chrNamesArray = np.loadtxt(ChrNamesFile, dtype=str, delimiter='\t')
            chrNamesInSPMs = chrNamesArray.tolist()

        signature2ProcessiveGroupLength2DistanceListDict = {}

        # read chrBased spms_df
        for chrShort in chrNamesInSPMs:
            chrLong = 'chr%s'%(chrShort)
            # print('%s %s' %(chrLong,chrShort))
            chrBased_spms_df = readChrBasedMutationDF(jobname, chrLong, singlePointMutationsFileName)

            if (chrBased_spms_df is not None):

                ###########################################
                columnNamesList = list(chrBased_spms_df.columns.values)
                contextIndex = columnNamesList.index('Context')
                # We assume that after the column named 'Context' there are the signature columns in tab separated way.
                signatures = columnNamesList[(contextIndex + 1):]
                ###########################################

                #delete unnecessary columns
                chrBased_spms_df.drop(['Chromosome','End', 'PyramidineStrand','Context'], inplace=True, axis=1)

                # df['Max'] = df.idxmax(axis=1).
                # left here
                # https://stackoverflow.com/questions/29919306/find-the-column-name-which-has-the-maximum-value-for-each-row

                #Each mutation will be assigned to the signature with the highest probability
                #Add new column
                chrBased_spms_df['Signature'] = (chrBased_spms_df.loc[:,signatures]).idxmax(axis=1)

                #Add new column. We take the that signatures's highest probability
                chrBased_spms_df['Probability'] = chrBased_spms_df[signatures].max(axis=1)

                #Drop the signatures columns
                chrBased_spms_df.drop(signatures,inplace=True,axis=1)

                # Prepare the poolInputList
                poolInputList = []
                sampleBased_chrBased_spms_df_grouped = chrBased_spms_df.groupby('Sample')

                for sample, sampleBased_chrBased_spms_df in sampleBased_chrBased_spms_df_grouped:
                    inputList = []
                    # inputList.append(chrLong)
                    # inputList.append(sample)
                    # inputList.append(signatures)
                    inputList.append(sampleBased_chrBased_spms_df)
                    inputList.append(considerProbability)
                    poolInputList.append(inputList)

                # Call the function findProcessiveGroups
                allSamples_chrBased_signature2ProcessiveGroupLength2DistanceListDict_list = pool.map(findProcessiveGroupsForInputList,poolInputList)

                #Accumuate all the list of dictionaries coming from all splits of each chromosome
                allSamples_chrBased_signature2ProcessiveGroupLength2DistanceListDict = accumulateListofDicts(allSamples_chrBased_signature2ProcessiveGroupLength2DistanceListDict_list)

                #Accumulate the dictionary coming from each chromosome
                accumulateDict(allSamples_chrBased_signature2ProcessiveGroupLength2DistanceListDict,signature2ProcessiveGroupLength2DistanceListDict)

        signature2ProcessiveGroupLength2PropertiesDict = findMedians(signature2ProcessiveGroupLength2DistanceListDict)
        writeDictionary(signature2ProcessiveGroupLength2PropertiesDict, jobname, 'Signature2ProcessiveGroupLength2PropertiesDict.txt', PROCESSIVITY, ProcessiveGroupStructEncoder)
####################################################################################


####################################################################################
def readSinglePointMutationsFindProcessivityGroupsWithoutMultiProcessing(jobname,singlePointMutationsFileName):
    if (singlePointMutationsFileName != NOTSET):

        signature2ProcessiveGroupLength2DistanceListDict = {}

        # Load the chrnames in single point mutations data
        ChrNamesFile = os.path.join(current_abs_path, ONE_DIRECTORY_UP, ONE_DIRECTORY_UP, OUTPUT, jobname, DATA,ChrNamesInSPMsFilename)
        if (os.path.exists(ChrNamesFile)):
            chrNamesArray = np.loadtxt(ChrNamesFile, dtype=str, delimiter='\t')
            chrNamesInSPMs = list(chrNamesArray)

        # read chrBased spms_df
        for chrShort in chrNamesInSPMs:
            chrLong = 'chr%s'%(chrShort)
            # print('%s %s' %(chrLong,chrShort))
            chrBased_spms_df = readChrBasedMutationDF(jobname, chrLong, singlePointMutationsFileName)
            # print(chrBased_spms_df.shape)

            ###########################################
            columnNamesList = list(chrBased_spms_df.columns.values)
            contextIndex = columnNamesList.index('Context')
            # We assume that after the column named 'Context' there are the signature columns in tab separated way.
            signatures = columnNamesList[(contextIndex + 1):]
            ###########################################

            #delete unnecessary columns
            chrBased_spms_df.drop(['Sample','Chromosome','End', 'PyramidineStrand','Context'], inplace=True, axis=1)

            #sort
            chrBased_spms_df.sort_values('Start', inplace=True)

            chrBased_signature2ProcessiveGroupLength2DistanceListDict = findProcessiveGroups(chrLong,chrBased_spms_df,signatures)

            accumulateDict(chrBased_signature2ProcessiveGroupLength2DistanceListDict, signature2ProcessiveGroupLength2DistanceListDict)

            writeDictionary(signature2ProcessiveGroupLength2DistanceListDict, jobname,'Signature2ProcessiveGroupLength2DistanceListDict.txt',PROCESSIVITY,None)

            findMedians(signature2ProcessiveGroupLength2DistanceListDict)
####################################################################################


####################################################################################
def convertStr2Bool(mystr):
    if (mystr=='True'):
        return True
    else:
        return False
####################################################################################



##################################################################################
def processivityAnalysis(jobname,singlePointMutationsFileName,considerProbability):
    print('########################## ProcessivityAnalysis starts ###########################')
    readSinglePointMutationsFindProcessivityGroupsWithMultiProcessing(jobname, singlePointMutationsFileName,considerProbability)
    print('########################## ProcessivityAnalysis ends #############################')
##################################################################################


