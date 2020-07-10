# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018 Burcak Otlu

import multiprocessing
import numpy as np
import pandas as pd
import os

from json import JSONEncoder

from SigProfilerTopography.source.commons.TopographyCommons import PROCESSIVITY
from SigProfilerTopography.source.commons.TopographyCommons import DATA

from SigProfilerTopography.source.commons.TopographyCommons import START
from SigProfilerTopography.source.commons.TopographyCommons import MUTATION
from SigProfilerTopography.source.commons.TopographyCommons import PYRAMIDINESTRAND
from SigProfilerTopography.source.commons.TopographyCommons import SIMULATION_NUMBER
from SigProfilerTopography.source.commons.TopographyCommons import CHROM
from SigProfilerTopography.source.commons.TopographyCommons import MUTATIONLONG
from SigProfilerTopography.source.commons.TopographyCommons import TRANSCRIPTIONSTRAND
from SigProfilerTopography.source.commons.TopographyCommons import SAMPLE

from SigProfilerTopography.source.commons.TopographyCommons import SUBS

from SigProfilerTopography.source.commons.TopographyCommons import SBS96
from SigProfilerTopography.source.commons.TopographyCommons import SBS384
from SigProfilerTopography.source.commons.TopographyCommons import SBS1536
from SigProfilerTopography.source.commons.TopographyCommons import SBS3072

from SigProfilerTopography.source.commons.TopographyCommons import readChrBasedMutationsDF
from SigProfilerTopography.source.commons.TopographyCommons import writeDictionary

from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC
from SigProfilerTopography.source.commons.TopographyCommons import memory_usage

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
def fillDict(x,signature2ProcessiveGroupLength2DistanceListDict, considerProbability,signature_cutoff_numberofmutations_averageprobability_df):
    processiveGroupLength = x['ProcessiveGroupLength']
    distance = x['Distance']
    signature = x['Signature']
    probability = x['Probability']

    #We require signature to be in signature2PropertiesListDict
    if signature in signature_cutoff_numberofmutations_averageprobability_df['signature'].unique():
        cutoff=float(signature_cutoff_numberofmutations_averageprobability_df[signature_cutoff_numberofmutations_averageprobability_df['signature']==signature]['cutoff'].values[0])
        if (considerProbability and (probability>=cutoff)):
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
#DEC 22, 2019
def findProcessiveGroupsForApplySync(simNum,chrLong,sample,sampleBased_chrBased_spms_df,considerProbabilityInProcessivityAnalysis,subsSignature_cutoff_numberofmutations_averageprobability_df,verbose):
    #Sort it
    sampleBased_chrBased_spms_df.sort_values(START, inplace=True)

    return findProcessiveGroups(simNum,chrLong,sample,sampleBased_chrBased_spms_df,considerProbabilityInProcessivityAnalysis,subsSignature_cutoff_numberofmutations_averageprobability_df,verbose)
####################################################################################


####################################################################################
def findProcessiveGroupsForInputList(inputList):
    simNum=inputList[0]
    chrLong=inputList[1]
    sample=inputList[2]
    sampleBased_chrBased_spms_df = inputList[3]
    considerProbabilityInProcessivityAnalysis = inputList[4]
    subsSignature_cutoff_numberofmutations_averageprobability_df=inputList[5]
    verbose=inputList[6]

    #Sort it
    sampleBased_chrBased_spms_df.sort_values(START, inplace=True)

    return findProcessiveGroups(simNum,chrLong,sample,sampleBased_chrBased_spms_df,considerProbabilityInProcessivityAnalysis,subsSignature_cutoff_numberofmutations_averageprobability_df,verbose)
####################################################################################



####################################################################################
def findProcessiveGroups(simNum,chrLong,sample,sorted_sampleBased_chrBased_spms_df,considerProbabilityInProcessivityAnalysis,signature_cutoff_numberofmutations_averageprobability_df,verbose):
    signature2ProcessiveGroupLength2DistanceListDict = {}

    if verbose: print('\tVerbose Worker pid %s memory_usage %.2f MB simNum:%d chrLong:%s sample:%s findProcessiveGroups starts' %(str(os.getpid()), memory_usage(),simNum,chrLong,sample))
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
                 signature_cutoff_numberofmutations_averageprobability_df=signature_cutoff_numberofmutations_averageprobability_df,
                 axis=1)

    if verbose: print('\tVerbose Worker pid %s memory_usage %.2f MB simNum:%d chrLong:%s sample:%s findProcessiveGroups ends' %(str(os.getpid()), memory_usage(),simNum,chrLong,sample))

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

                processiveGroupDict={}
                processiveGroupDict['processiveGroupLength']=processiveGroupLength
                processiveGroupDict['numberofProcessiveGroups']=numberofProcessiveGroups
                processiveGroupDict['medianofNumberofProcessiveGroupsinMB'] = median

                # signature2ProcessiveGroupLength2PropertiesDict[signature][processiveGroupLength]= processiveGroupStruct
                signature2ProcessiveGroupLength2PropertiesDict[signature][processiveGroupLength] = processiveGroupDict

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
#This method can be customized for other dictionaries
def writeDictionaryAsADataframe(signature2ProcessiveGroupLength2PropertiesDict,filePath):
    L = sorted([(signature, processiveGroupLength, v1['numberofProcessiveGroups'],v1['medianofNumberofProcessiveGroupsinMB']) for signature, v in signature2ProcessiveGroupLength2PropertiesDict.items() for processiveGroupLength, v1 in v.items()])
    df = pd.DataFrame(L,columns=['Signature', 'Processsive_Group_Length', 'Number_of_Processive_Groups', 'Median_of_Number_of_Processive_Groups_in_MB'])

    #write this dataframe
    df.to_csv(filePath, sep='\t', header=True, index=False)
####################################################################################


####################################################################################
def readSinglePointMutationsFindProcessivityGroupsWithMultiProcessing(mutation_types_contexts,chromNamesList,computation_type,outputDir,jobname,numofSimulations,considerProbabilityInProcessivityAnalysis,subsSignature_cutoff_numberofmutations_averageprobability_df,verbose):

    if ((SBS96 in mutation_types_contexts) or (SBS384 in mutation_types_contexts) or (SBS1536 in mutation_types_contexts) or (SBS3072 in mutation_types_contexts)):

        #####################################################################
        if verbose: print('\tVerbose computation_type:%s' % (computation_type))

        if computation_type == USING_APPLY_ASYNC:
            for simNum in range(0,numofSimulations+1):

                ####################################################################################
                numofProcesses = multiprocessing.cpu_count()
                pool = multiprocessing.Pool(numofProcesses)
                ####################################################################################

                ####################################################################################
                jobs = []
                ####################################################################################

                signature2ProcessiveGroupLength2DistanceListDict = {}

                ####################################################################################
                def accumulate_apply_async_result(small_signature2ProcessiveGroupLength2DistanceListDict):
                    for signature in small_signature2ProcessiveGroupLength2DistanceListDict:
                        for processiveGroupLength in small_signature2ProcessiveGroupLength2DistanceListDict[signature]:
                            if verbose: print('\tVerbose Worker pid %s memory_usage %.2f MB accumulate_apply_async_result starts' % (str(os.getpid()), memory_usage()))

                            if signature in signature2ProcessiveGroupLength2DistanceListDict:
                                if processiveGroupLength in signature2ProcessiveGroupLength2DistanceListDict[signature]:
                                    signature2ProcessiveGroupLength2DistanceListDict[signature][processiveGroupLength].extend(small_signature2ProcessiveGroupLength2DistanceListDict[signature][processiveGroupLength])
                                else:
                                    signature2ProcessiveGroupLength2DistanceListDict[signature][processiveGroupLength] = small_signature2ProcessiveGroupLength2DistanceListDict[signature][processiveGroupLength]
                            else:
                                signature2ProcessiveGroupLength2DistanceListDict[signature] = {}
                                signature2ProcessiveGroupLength2DistanceListDict[signature][processiveGroupLength] = small_signature2ProcessiveGroupLength2DistanceListDict[signature][processiveGroupLength]
                ####################################################################################

                ####################################################################################
                for chrLong in chromNamesList:
                    chrBased_spms_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, SUBS, simNum)
                    if (chrBased_spms_df is not None):
                        ###########################################
                        chrBased_spms_df.drop([SIMULATION_NUMBER], inplace=True, errors='ignore', axis=1)
                        columnNamesList = list(chrBased_spms_df.columns.values)
                        # We assume that after the column named 'Context' there are the signature columns in tab separated way.
                        mutationIndex = columnNamesList.index(MUTATION)
                        signatures = columnNamesList[(mutationIndex + 1):]
                        ###########################################
                        chrBased_spms_df.drop([CHROM, MUTATIONLONG, TRANSCRIPTIONSTRAND],inplace=True, errors='ignore', axis=1)
                        chrBased_spms_df['Signature'] = (chrBased_spms_df.loc[:, signatures]).idxmax(axis=1)

                        # Add new column. We take the that signatures's highest probability
                        chrBased_spms_df['Probability'] = chrBased_spms_df[signatures].max(axis=1)

                        # Drop the signatures columns
                        chrBased_spms_df.drop(signatures, inplace=True, axis=1)

                        sampleBased_chrBased_spms_df_grouped = chrBased_spms_df.groupby(SAMPLE)

                        for sample, sampleBased_chrBased_spms_df in sampleBased_chrBased_spms_df_grouped:
                            jobs.append(pool.apply_async(findProcessiveGroupsForApplySync,(simNum,chrLong,sample,sampleBased_chrBased_spms_df,considerProbabilityInProcessivityAnalysis,subsSignature_cutoff_numberofmutations_averageprobability_df,verbose),callback=accumulate_apply_async_result))
                ####################################################################################

                ################################
                if verbose: print('\tVerbose Processivity len(jobs):%d\n' % (len(jobs)))

                # wait for all jobs to finish
                for job in jobs:
                    if verbose: print('\tVerbose Processivity Worker pid %s job.get():%s ' %(str(os.getpid()), job.get()))
                ################################

                ####################################################################################
                pool.close()
                pool.join()
                ####################################################################################

                signature2ProcessiveGroupLength2PropertiesDict = findMedians(signature2ProcessiveGroupLength2DistanceListDict)

                # ####################################################################################
                # filename = 'Sim%d_Signature2ProcessiveGroupLength2PropertiesDict.txt' % (simNum)
                # writeDictionary(signature2ProcessiveGroupLength2PropertiesDict, outputDir, jobname, filename,PROCESSIVITY, None)
                # ####################################################################################

                ####################################################################################
                filename = 'Sim%d_Processivity.txt' % (simNum)
                filePath = os.path.join(outputDir, jobname, DATA, PROCESSIVITY, filename)
                writeDictionaryAsADataframe(signature2ProcessiveGroupLength2PropertiesDict,filePath)
                ####################################################################################

        else:
            for simNum in range(0,numofSimulations+1):
                ####################################################################################
                numofProcesses = multiprocessing.cpu_count()
                pool = multiprocessing.Pool(numofProcesses)
                ####################################################################################

                signature2ProcessiveGroupLength2DistanceListDict = {}

                ####################################################################################
                # read chrBased spms_df
                for chrLong in chromNamesList:
                    # print('%s %s' %(chrLong,chrShort))
                    chrBased_spms_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, SUBS, simNum)

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
                            inputList.append(simNum)
                            inputList.append(chrLong)
                            inputList.append(sample)
                            inputList.append(sampleBased_chrBased_spms_df)
                            inputList.append(considerProbabilityInProcessivityAnalysis)
                            inputList.append(subsSignature_cutoff_numberofmutations_averageprobability_df)
                            inputList.append(verbose)
                            poolInputList.append(inputList)

                        # Call the function findProcessiveGroups
                        # Sequential for each simulation and chromosome
                        # Parallel for each sample
                        allSamples_chrBased_signature2ProcessiveGroupLength2DistanceListDict_list = pool.map(findProcessiveGroupsForInputList, poolInputList)

                        # Accumuate all the list of dictionaries coming from all samples
                        allSamples_chrBased_signature2ProcessiveGroupLength2DistanceListDict = accumulateListofDicts(allSamples_chrBased_signature2ProcessiveGroupLength2DistanceListDict_list)

                        # Accumulate the dictionary coming from each chromosome
                        accumulateDict(allSamples_chrBased_signature2ProcessiveGroupLength2DistanceListDict,signature2ProcessiveGroupLength2DistanceListDict)
                ####################################################################################

                ####################################################################################
                pool.close()
                pool.join()
                ####################################################################################

                signature2ProcessiveGroupLength2PropertiesDict = findMedians(signature2ProcessiveGroupLength2DistanceListDict)

                ####################################################################################
                filename = 'Sim%d_Processivity.txt' % (simNum)
                filePath = os.path.join(outputDir, jobname, DATA, PROCESSIVITY, filename)
                writeDictionaryAsADataframe(signature2ProcessiveGroupLength2PropertiesDict,filePath)
                ####################################################################################

    #####################################################################

####################################################################################


##################################################################################
def processivityAnalysis(mutation_types_contexts,chromNamesList,computation_type,outputDir,jobname,numofSimulations,considerProbabilityInProcessivityAnalysis,subsSignature_cutoff_numberofmutations_averageprobability_df,verbose):
    print('\n#################################################################################')
    print('--- ProcessivityAnalysis starts')
    readSinglePointMutationsFindProcessivityGroupsWithMultiProcessing(mutation_types_contexts,chromNamesList,computation_type,outputDir,jobname,numofSimulations,considerProbabilityInProcessivityAnalysis,subsSignature_cutoff_numberofmutations_averageprobability_df,verbose)
    print('--- ProcessivityAnalysis ends')
    print('#################################################################################\n')
##################################################################################