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
#Return these columns  ['ProcessiveGroupLength', 'Distance', 'Signature', 'Probability']
def myfuncUpdated(x):
    return (x.shape[0], x['Start'].iloc[-1] - x['Start'].iloc[0],  x['Signature'].iloc[0], x['Probability'].agg(max))
#########################################################################


#########################################################################
def fillDict(x,signature2ProcessiveGroupLength2DistanceListDict, considerProbability,signature2PropertiesListDict):

    processiveGroupLength = x['ProcessiveGroupLength']
    distance = x['Distance']
    signature = x['Signature']
    probability = x['Probability']

    #Old way
    # if (considerProbability and (probability >= PROCESSIVITY_MUTATION_SIGNATURE_PROBABILITY_THRESHOLD)):

    #New way
    #We require signature to be in signature2PropertiesListDict
    if (signature in signature2PropertiesListDict):
        if (considerProbability and (probability>=float(signature2PropertiesListDict[signature][0]))):
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
    signature2PropertiesListDict=inputList[2]

    #Sort it
    sampleBased_chrBased_spms_df.sort_values(START, inplace=True)

    return findProcessiveGroups(sampleBased_chrBased_spms_df,considerProbabilityInProcessivityAnalysis,signature2PropertiesListDict)
####################################################################################



####################################################################################
def findProcessiveGroups(sorted_sampleBased_chrBased_spms_df,considerProbabilityInProcessivityAnalysis,signature2PropertiesListDict):
    signature2ProcessiveGroupLength2DistanceListDict = {}

    #They must be same type of mutation e.g.: T>A
    #They must be resulted from same signature
    #They must be on the same strand
    if (sorted_sampleBased_chrBased_spms_df is not None):
        sorted_sampleBased_chrBased_spms_df['subgroup'] = ((sorted_sampleBased_chrBased_spms_df[MUTATION] != sorted_sampleBased_chrBased_spms_df[MUTATION].shift(1)) |
                                                           (sorted_sampleBased_chrBased_spms_df['Signature'] != sorted_sampleBased_chrBased_spms_df['Signature'].shift(1)) |
                                                           (sorted_sampleBased_chrBased_spms_df[PYRAMIDINESTRAND] != sorted_sampleBased_chrBased_spms_df[PYRAMIDINESTRAND].shift(1))) .cumsum()


        series_new = sorted_sampleBased_chrBased_spms_df.groupby('subgroup', as_index=False).apply(myfuncUpdated)
        columns = ['ProcessiveGroupLength', 'Distance', 'Signature', 'Probability']
        df = series_new.apply(pd.Series)
        df.columns = columns
        df = df[df['ProcessiveGroupLength'].ne(1)]

        df.apply(fillDict,
                 signature2ProcessiveGroupLength2DistanceListDict=signature2ProcessiveGroupLength2DistanceListDict,
                 considerProbability=considerProbabilityInProcessivityAnalysis,
                 signature2PropertiesListDict=signature2PropertiesListDict,
                 axis=1)

    return signature2ProcessiveGroupLength2DistanceListDict
####################################################################################

####################################################################################
def findMedians(signature2ProcessiveGroupLength2GroupSizeListDict):
    signature2ProcessiveGroupLength2PropertiesDict = {}

    for signature in sorted(signature2ProcessiveGroupLength2GroupSizeListDict.keys()):
        signature2ProcessiveGroupLength2PropertiesDict[signature] = {}
        for processiveGroupLength in sorted(signature2ProcessiveGroupLength2GroupSizeListDict[signature].keys()):

            processiveGroupDistanceArray = np.array(signature2ProcessiveGroupLength2GroupSizeListDict[signature][processiveGroupLength])
            #if we remove the 0.0 size processive groups we will get rid of infinity

            processiveGroupDistanceArray = processiveGroupDistanceArray[processiveGroupDistanceArray != 0]
            if (processiveGroupDistanceArray.size>0):
                processiveGroupDistanceArray = (1/processiveGroupDistanceArray) * processiveGroupLength * 1000000
                median = np.median(processiveGroupDistanceArray)
                numberofProcessiveGroups = len(processiveGroupDistanceArray)

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
def readSinglePointMutationsFindProcessivityGroupsWithMultiProcessing(mutation_types_contexts,chromNamesList,outputDir,jobname,numofSimulations,considerProbabilityInProcessivityAnalysis,signature2PropertiesListDict):

    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numofProcesses)

    if ((SBS96 in mutation_types_contexts) or (SBS384 in mutation_types_contexts) or (SBS1536 in mutation_types_contexts) or (SBS3072 in mutation_types_contexts)):

        #####################################################################
        for simNum in range(0,numofSimulations+1):
            signature2ProcessiveGroupLength2DistanceListDict = {}

            # read chrBased spms_df
            for chrLong in chromNamesList:
                # print('%s %s' %(chrLong,chrShort))
                chrBased_spms_df = readChrBasedSubsDF(outputDir, jobname, chrLong, SUBS, simNum)

                if (chrBased_spms_df is not None):

                    ###########################################
                    chrBased_spms_df.drop([SIMULATION_NUMBER], inplace=True, errors='ignore', axis=1)
                    columnNamesList = list(chrBased_spms_df.columns.values)
                    # We assume that after the column named 'Context' there are the signature columns in tab separated way.
                    mutationIndex = columnNamesList.index(MUTATION)
                    signatures = columnNamesList[(mutationIndex + 1):]
                    ###########################################

                    # What are the columns in subs_df?
                    # Sample  Chrom   Start   MutationLong    PyramidineStrand        TranscriptionStrand     Mutation        SBS1    SBS2    SBS3    SBS13   SBS26   SBS40   SBS44

                    # delete unnecessary columns
                    # chrBased_spms_df.drop([CHROM, MUTATIONLONG, PYRAMIDINESTRAND, TRANSCRIPTIONSTRAND],inplace=True, errors='ignore', axis=1)
                    chrBased_spms_df.drop([CHROM, MUTATIONLONG, TRANSCRIPTIONSTRAND],inplace=True, errors='ignore', axis=1)

                    # df['Max'] = df.idxmax(axis=1).
                    # left here
                    # https://stackoverflow.com/questions/29919306/find-the-column-name-which-has-the-maximum-value-for-each-row

                    # Each mutation will be assigned to the signature with the highest probability
                    # Add new column
                    chrBased_spms_df['Signature'] = (chrBased_spms_df.loc[:, signatures]).idxmax(axis=1)

                    # Add new column. We take the that signatures's highest probability
                    chrBased_spms_df['Probability'] = chrBased_spms_df[signatures].max(axis=1)

                    # Drop the signatures columns
                    chrBased_spms_df.drop(signatures, inplace=True, axis=1)

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
                        inputList.append(signature2PropertiesListDict)
                        poolInputList.append(inputList)

                    # Call the function findProcessiveGroups
                    # Sequential for each simulation and chromosome
                    # Parallel for each sample
                    allSamples_chrBased_signature2ProcessiveGroupLength2DistanceListDict_list = pool.map(findProcessiveGroupsForInputList, poolInputList)

                    # Accumuate all the list of dictionaries coming from all samples
                    allSamples_chrBased_signature2ProcessiveGroupLength2DistanceListDict = accumulateListofDicts(allSamples_chrBased_signature2ProcessiveGroupLength2DistanceListDict_list)

                    # Accumulate the dictionary coming from each chromosome
                    accumulateDict(allSamples_chrBased_signature2ProcessiveGroupLength2DistanceListDict,signature2ProcessiveGroupLength2DistanceListDict)

            signature2ProcessiveGroupLength2PropertiesDict = findMedians(signature2ProcessiveGroupLength2DistanceListDict)
            filename = 'Sim%d_Signature2ProcessiveGroupLength2PropertiesDict.txt' %(simNum)

            # #For information
            # if (simNum==0):
            #     print('For original data signature2ProcessiveGroupLength2PropertiesDict')
            #     print(signature2ProcessiveGroupLength2PropertiesDict)

            writeDictionary(signature2ProcessiveGroupLength2PropertiesDict, outputDir, jobname,filename, PROCESSIVITY,ProcessiveGroupStructEncoder)

        #####################################################################

####################################################################################


##################################################################################
def processivityAnalysis(mutation_types_contexts,chromNamesList,outputDir,jobname,numofSimulations,considerProbabilityInProcessivityAnalysis,signature2PropertiesListDict):
    print('\n#################################################################################')
    print('--- ProcessivityAnalysis starts')
    readSinglePointMutationsFindProcessivityGroupsWithMultiProcessing(mutation_types_contexts,chromNamesList,outputDir,jobname,numofSimulations,considerProbabilityInProcessivityAnalysis,signature2PropertiesListDict)
    print('--- ProcessivityAnalysis ends')
    print('#################################################################################\n')
##################################################################################