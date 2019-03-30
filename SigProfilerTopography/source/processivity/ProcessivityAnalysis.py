# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018 Burcak Otlu

from json import JSONEncoder
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
def myfuncUpdated(x):
    return (x.shape[0], x['Start'].iloc[-1] - x['Start'].iloc[0],  x['Signature'].iloc[0], x['Probability'].agg(max))
#########################################################################


#########################################################################
def fillDictUpdated(x,signature2ProcessiveGroupLength2DistanceListDict, considerProbability):

    processiveGroupLength = x['ProcessiveGroupLength']
    distance = x['Distance']
    signature = x['Signature']
    probability = x['Probability']

    if (considerProbability and (probability >= SUBSTITUTION_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD)):
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

    sampleBased_chrBased_spms_df = inputList[0]
    considerProbabilityInProcessivityAnalysis = inputList[1]

    #Sort it
    sampleBased_chrBased_spms_df.sort_values(START, inplace=True)

    return findProcessiveGroupsUpdated(sampleBased_chrBased_spms_df,considerProbabilityInProcessivityAnalysis)
####################################################################################



####################################################################################
#New Way
def findProcessiveGroupsUpdated(sorted_sampleBased_chrBased_spms_df,considerProbabilityInProcessivityAnalysis):
    signature2ProcessiveGroupLength2DistanceListDict = {}

    if (sorted_sampleBased_chrBased_spms_df is not None):
        sorted_sampleBased_chrBased_spms_df['subgroup'] = ((sorted_sampleBased_chrBased_spms_df['Mutation'] != sorted_sampleBased_chrBased_spms_df['Mutation'].shift(1)) |
                                                           (sorted_sampleBased_chrBased_spms_df['Signature'] != sorted_sampleBased_chrBased_spms_df['Signature'].shift(1))) .cumsum()


        series_new = sorted_sampleBased_chrBased_spms_df.groupby('subgroup', as_index=False).apply(myfuncUpdated)
        columns = ['ProcessiveGroupLength', 'Distance', 'Signature', 'Probability']
        df = series_new.apply(pd.Series)
        df.columns = columns
        df = df[df['ProcessiveGroupLength'].ne(1)]

        df.apply(fillDictUpdated,signature2ProcessiveGroupLength2DistanceListDict=signature2ProcessiveGroupLength2DistanceListDict, considerProbability=considerProbabilityInProcessivityAnalysis , axis=1)

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
def readSinglePointMutationsFindProcessivityGroupsWithMultiProcessing(chromNamesList,outputDir,jobname, singlePointMutationsFileName,considerProbabilityInProcessivityAnalysis):

    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)

    if (singlePointMutationsFileName != NOTSET):
        # Load the chrnames in single point mutations data


        signature2ProcessiveGroupLength2DistanceListDict = {}

        # read chrBased spms_df
        for chrLong in chromNamesList:
            # print('%s %s' %(chrLong,chrShort))
            chrBased_spms_df = readChrBasedSubsDF(outputDir,jobname, chrLong, singlePointMutationsFileName)

            if (chrBased_spms_df is not None):

                ###########################################
                columnNamesList = list(chrBased_spms_df.columns.values)
                contextIndex = columnNamesList.index('Context')
                # We assume that after the column named 'Context' there are the signature columns in tab separated way.
                signatures = columnNamesList[(contextIndex + 1):]
                ###########################################

                #delete unnecessary columns
                chrBased_spms_df.drop([CHROM, END, PYRAMIDINESTRAND,CONTEXT], inplace=True, errors='ignore',axis=1)

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
                sampleBased_chrBased_spms_df_grouped = chrBased_spms_df.groupby(SAMPLE)

                for sample, sampleBased_chrBased_spms_df in sampleBased_chrBased_spms_df_grouped:
                    inputList = []
                    # inputList.append(chrLong)
                    # inputList.append(sample)
                    # inputList.append(signatures)
                    inputList.append(sampleBased_chrBased_spms_df)
                    inputList.append(considerProbabilityInProcessivityAnalysis)
                    poolInputList.append(inputList)

                # Call the function findProcessiveGroups
                allSamples_chrBased_signature2ProcessiveGroupLength2DistanceListDict_list = pool.map(findProcessiveGroupsForInputList,poolInputList)

                #Accumuate all the list of dictionaries coming from all splits of each chromosome
                allSamples_chrBased_signature2ProcessiveGroupLength2DistanceListDict = accumulateListofDicts(allSamples_chrBased_signature2ProcessiveGroupLength2DistanceListDict_list)

                #Accumulate the dictionary coming from each chromosome
                accumulateDict(allSamples_chrBased_signature2ProcessiveGroupLength2DistanceListDict,signature2ProcessiveGroupLength2DistanceListDict)

        signature2ProcessiveGroupLength2PropertiesDict = findMedians(signature2ProcessiveGroupLength2DistanceListDict)
        writeDictionary(signature2ProcessiveGroupLength2PropertiesDict, outputDir,jobname, 'Signature2ProcessiveGroupLength2PropertiesDict.txt', PROCESSIVITY, ProcessiveGroupStructEncoder)
####################################################################################




####################################################################################
def convertStr2Bool(mystr):
    if (mystr=='True'):
        return True
    else:
        return False
####################################################################################



##################################################################################
def processivityAnalysis(chromNamesList,outputDir,jobname,singlePointMutationsFileName,considerProbabilityInProcessivityAnalysis):
    print('########################## ProcessivityAnalysis starts ###########################')
    readSinglePointMutationsFindProcessivityGroupsWithMultiProcessing(chromNamesList,outputDir,jobname, singlePointMutationsFileName,considerProbabilityInProcessivityAnalysis)
    print('########################## ProcessivityAnalysis ends #############################')
##################################################################################