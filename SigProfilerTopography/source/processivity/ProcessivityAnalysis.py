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
import time


from json import JSONEncoder

from source.commons.TopographyCommons import PROCESSIVITY
from source.commons.TopographyCommons import DATA

from source.commons.TopographyCommons import START
from source.commons.TopographyCommons import MUTATION
from source.commons.TopographyCommons import PYRAMIDINESTRAND
from source.commons.TopographyCommons import SIMULATION_NUMBER
from source.commons.TopographyCommons import CHROM
from source.commons.TopographyCommons import MUTATIONLONG
from source.commons.TopographyCommons import TRANSCRIPTIONSTRAND
from source.commons.TopographyCommons import SAMPLE

from source.commons.TopographyCommons import SUBS

from source.commons.TopographyCommons import SBS96
from source.commons.TopographyCommons import SBS384
from source.commons.TopographyCommons import SBS1536
from source.commons.TopographyCommons import SBS3072

from source.commons.TopographyCommons import readChrBasedMutationsDF
from source.commons.TopographyCommons import writeDictionary

from source.commons.TopographyCommons import USING_APPLY_ASYNC
from source.commons.TopographyCommons import memory_usage

from source.commons.TopographyCommons import CONSIDER_COUNT
from source.commons.TopographyCommons import CONSIDER_DISTANCE
from source.commons.TopographyCommons import CONSIDER_DISTANCE_ALL_SAMPLES_TOGETHER


####################################################################################
#DEC 22, 2019
def findProcessiveGroupsForApplySyncWithDistance(simNum,
                                                 chrLong,
                                                 sample,
                                                 sampleBased_chrBased_spms_df,
                                                 considerProbabilityInProcessivityAnalysis,
                                                 subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                 verbose):
    #Sort it
    sampleBased_chrBased_spms_df.sort_values(START, inplace=True)

    return findProcessiveGroupsWithDistance(simNum,
                                            chrLong,
                                            sample,
                                            sampleBased_chrBased_spms_df,
                                            considerProbabilityInProcessivityAnalysis,
                                            subsSignature_cutoff_numberofmutations_averageprobability_df,
                                            verbose)
####################################################################################


####################################################################################
#SEP 17, 2020
#Did not work as expected
def findProcessiveGroupsForApplySyncWithDistanceAllSamples(simNum,chrLong,chrBased_spms_df,considerProbabilityInProcessivityAnalysis,subsSignature_cutoff_numberofmutations_averageprobability_df,verbose):
    #Sort it
    chrBased_spms_df.sort_values(START, inplace=True)

    return findProcessiveGroupsWithDistanceAllSamples(simNum,chrLong,chrBased_spms_df,considerProbabilityInProcessivityAnalysis,subsSignature_cutoff_numberofmutations_averageprobability_df,verbose)
####################################################################################

####################################################################################
#SEP 17, 2020
def findProcessiveGroupsForApplySyncWithCount(simNum,chrLong,sample,sampleBased_chrBased_spms_df,considerProbabilityInProcessivityAnalysis,subsSignature_cutoff_numberofmutations_averageprobability_df,verbose):
    #Sort it
    sampleBased_chrBased_spms_df.sort_values(START, inplace=True)

    return findProcessiveGroupsWithCount(simNum,chrLong,sample,sampleBased_chrBased_spms_df,considerProbabilityInProcessivityAnalysis,subsSignature_cutoff_numberofmutations_averageprobability_df,verbose)
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

    return findProcessiveGroupsForApplySyncWithDistance(simNum,chrLong,sample,sampleBased_chrBased_spms_df,considerProbabilityInProcessivityAnalysis,subsSignature_cutoff_numberofmutations_averageprobability_df,verbose)
####################################################################################


####################################################################################
def findProcessiveGroupsWithCount(simNum,chrLong,sample,sorted_sampleBased_chrBased_spms_df,considerProbabilityInProcessivityAnalysis,signature_cutoff_numberofmutations_averageprobability_df,verbose):

    my_dict={}

    if verbose: print('\tVerbose Worker pid %s memory_usage %.2f MB simNum:%d chrLong:%s sample:%s findProcessiveGroups starts' %(str(os.getpid()), memory_usage(),simNum,chrLong,sample))
    #They must be same type of mutation e.g.: T>A
    #They must be resulted from same signature
    #They must be on the same strand
    if (sorted_sampleBased_chrBased_spms_df is not None):

        ################################
        #As long as mutation, signature, and pyrimidine strands are the same continue to accumlate, if one of them is not the same calculate cumsum and start again
        sorted_sampleBased_chrBased_spms_df['subgroup'] = ((sorted_sampleBased_chrBased_spms_df[MUTATION].ne(sorted_sampleBased_chrBased_spms_df[MUTATION].shift())) |
                                                            (sorted_sampleBased_chrBased_spms_df['Signature'].ne(sorted_sampleBased_chrBased_spms_df['Signature'].shift())) |
                                                            (sorted_sampleBased_chrBased_spms_df[PYRAMIDINESTRAND].ne(sorted_sampleBased_chrBased_spms_df[PYRAMIDINESTRAND].shift()))) .cumsum()
        ################################


        ################################
        df = sorted_sampleBased_chrBased_spms_df.groupby("subgroup").agg(
            Signature=pd.NamedAgg(column='Signature', aggfunc="first"),
            ProcessiveGroupLength=pd.NamedAgg(column=MUTATION, aggfunc="count"),
            Probability=pd.NamedAgg(column='Probability', aggfunc="max"),)
        #agg().reset_index() can be used. not necessary
        ################################

        ################################
        # Remove rows with processive groups of length 1
        df = df[df['ProcessiveGroupLength'].ne(1)]
        ################################

        ################################
        if considerProbabilityInProcessivityAnalysis:
            # # Remove rows with processive groups of length 1
            df=df.loc[df['Probability'] >= df['Signature'].map(signature_cutoff_numberofmutations_averageprobability_df.set_index('signature')['cutoff'])]
        ################################

        ################################
        my_dict = df.groupby(['Signature', 'ProcessiveGroupLength']).size().to_dict()
        ################################

    if verbose: print('\tVerbose Worker pid %s memory_usage %.2f MB simNum:%d chrLong:%s sample:%s findProcessiveGroups ends' %(str(os.getpid()), memory_usage(),simNum,chrLong,sample))

    return my_dict
####################################################################################

####################################################################################
#Sep 17, 2020
#Did not work as expected
def findProcessiveGroupsWithDistanceAllSamples(simNum,chrLong,chrBased_spms_df,considerProbabilityInProcessivityAnalysis,signature_cutoff_numberofmutations_averageprobability_df,verbose):
    my_all_samples_dict={}

    if verbose: print('\tVerbose Worker pid %s memory_usage %.2f MB simNum:%d chrLong:%s findProcessiveGroups starts' %(str(os.getpid()), memory_usage(),simNum,chrLong))
    #They must be same type of mutation e.g.: T>A
    #They must be resulted from same signature
    #They must be on the same strand
    if (chrBased_spms_df is not None):

        for sample, sorted_sampleBased_chrBased_spms_df in chrBased_spms_df.groupby('Sample'):

            ################################
            #As long as mutation, signature, and pyrimidine strands are the same continue to accumlate, if one of them is not the same calculate cumsum and start again
            sorted_sampleBased_chrBased_spms_df['subgroup'] = ((sorted_sampleBased_chrBased_spms_df[MUTATION].ne(sorted_sampleBased_chrBased_spms_df[MUTATION].shift())) |
                                                                (sorted_sampleBased_chrBased_spms_df['Signature'].ne(sorted_sampleBased_chrBased_spms_df['Signature'].shift())) |
                                                                (sorted_sampleBased_chrBased_spms_df[PYRAMIDINESTRAND].ne(sorted_sampleBased_chrBased_spms_df[PYRAMIDINESTRAND].shift()))) .cumsum()
            ################################

            ################################
            df = sorted_sampleBased_chrBased_spms_df.groupby("subgroup").agg(
                Signature=pd.NamedAgg(column='Signature', aggfunc="first"),
                ProcessiveGroupLength=pd.NamedAgg(column=MUTATION, aggfunc="count"),
                Probability=pd.NamedAgg(column='Probability', aggfunc="max"),
                LastDistance=pd.NamedAgg(column='Start', aggfunc="last"),
                FirstDistance=pd.NamedAgg(column='Start', aggfunc="first"),
                ).assign(Distance=lambda x: x.pop('LastDistance') - x.pop('FirstDistance'))
            #agg().reset_index() can be used. not necessary
            ################################

            ################################
            # Remove rows with processive groups of length 1
            df = df[df['ProcessiveGroupLength'].ne(1)]
            ################################

            ################################
            if considerProbabilityInProcessivityAnalysis:
                # Remove rows with processive groups of length 1
                df=df.loc[df['Probability'] >= df['Signature'].map(signature_cutoff_numberofmutations_averageprobability_df.set_index('signature')['cutoff'])]
            ################################

            ################################
            my_dict = {k: np.hstack(v) for k, v in df.groupby(['Signature', 'ProcessiveGroupLength'])['Distance']}
            ################################

            ################################
            accumulateDict(my_dict,my_all_samples_dict)
            ################################

    if verbose: print('\tVerbose Worker pid %s memory_usage %.2f MB simNum:%d chrLong:%s findProcessiveGroups ends' %(str(os.getpid()), memory_usage(),simNum,chrLong))

    return my_all_samples_dict
####################################################################################

####################################################################################
def findProcessiveGroupsWithDistance(simNum,
                                     chrLong,
                                     sample,
                                     sorted_sampleBased_chrBased_spms_df,
                                     considerProbabilityInProcessivityAnalysis,
                                     signature_cutoff_numberofmutations_averageprobability_df,
                                     verbose):

    my_dict={}

    if verbose: print('\tVerbose Worker pid %s memory_usage %.2f MB simNum:%d chrLong:%s sample:%s findProcessiveGroups starts' %(str(os.getpid()), memory_usage(),simNum,chrLong,sample))
    #They must be same type of mutation e.g.: T>A
    #They must be resulted from same signature
    #They must be on the same strand
    if (sorted_sampleBased_chrBased_spms_df is not None):

        ################################
        #As long as mutation, signature, and pyrimidine strands are the same continue to accumlate, if one of them is not the same calculate cumsum and start again
        #If there are 4 consecutive rows with same mutatio, signature and pyrimidine strand, subgroup will be 4
        sorted_sampleBased_chrBased_spms_df['subgroup'] = ((sorted_sampleBased_chrBased_spms_df[MUTATION].ne(sorted_sampleBased_chrBased_spms_df[MUTATION].shift())) |
                                                            (sorted_sampleBased_chrBased_spms_df['Signature'].ne(sorted_sampleBased_chrBased_spms_df['Signature'].shift())) |
                                                            (sorted_sampleBased_chrBased_spms_df[PYRAMIDINESTRAND].ne(sorted_sampleBased_chrBased_spms_df[PYRAMIDINESTRAND].shift()))) .cumsum()
        ################################


        ################################
        df = sorted_sampleBased_chrBased_spms_df.groupby("subgroup").agg(
            Signature=pd.NamedAgg(column='Signature', aggfunc="first"),
            ProcessiveGroupLength=pd.NamedAgg(column=MUTATION, aggfunc="count"),
            Probability=pd.NamedAgg(column='Probability', aggfunc="max"),
            LastDistance=pd.NamedAgg(column='Start', aggfunc="last"),
            FirstDistance=pd.NamedAgg(column='Start', aggfunc="first"),
            ).assign(Distance=lambda x: x.pop('LastDistance') - x.pop('FirstDistance'))
        #agg().reset_index() can be used. not necessary
        ################################


        ################################
        # Remove rows with processive groups of length 1
        df = df[df['ProcessiveGroupLength'].ne(1)]
        ################################


        ################################
        if considerProbabilityInProcessivityAnalysis:
            df=df.loc[df['Probability'] >= df['Signature'].map(signature_cutoff_numberofmutations_averageprobability_df.set_index('signature')['cutoff'])]
        ################################

        ################################
        my_dict = {k: np.hstack(v) for k, v in df.groupby(['Signature', 'ProcessiveGroupLength'])['Distance']}
        ################################

    if verbose: print('\tVerbose Worker pid %s memory_usage %.2f MB simNum:%d chrLong:%s sample:%s findProcessiveGroups ends' %(str(os.getpid()), memory_usage(),simNum,chrLong,sample))

    return my_dict
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
                numberofProcessiveGroups = len(processiveGroupDistanceArray)
                median_distance_between_last_first_mutations = np.median(processiveGroupDistanceArray)
                median_distance_between_consecutive_mutations = np.median(processiveGroupDistanceArray/(processiveGroupLength-1))

                number_of_mutations_within_1MB_array = (1/processiveGroupDistanceArray) * processiveGroupLength * 1000000
                median_number_of_mutations_within_1MB = np.median(number_of_mutations_within_1MB_array)

                processiveGroupDict={}
                processiveGroupDict['processive_group_length']=processiveGroupLength
                processiveGroupDict['number_of_processive_groups']=numberofProcessiveGroups
                processiveGroupDict['median_distance_between_last_first_mutations'] = median_distance_between_last_first_mutations
                processiveGroupDict['median_distance_between_consecutive_mutations'] = median_distance_between_consecutive_mutations
                processiveGroupDict['median_number_of_mutations_within_1MB'] = median_number_of_mutations_within_1MB

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
    L = sorted([(signature,
                 processiveGroupLength,
                 v1['number_of_processive_groups'],
                 v1['median_distance_between_last_first_mutations'],
                 v1['median_distance_between_consecutive_mutations'],
                 v1['median_number_of_mutations_within_1MB']) for signature, v in signature2ProcessiveGroupLength2PropertiesDict.items() for processiveGroupLength, v1 in v.items()])
    df = pd.DataFrame(L,columns=['Signature',
                                 'Processsive_Group_Length',
                                 'Number_of_Processive_Groups',
                                 'Median_Distance_Between_Last_First_Mutations',
                                 'Median_Distance_Between_Consecutive_Mutations',
                                 'Median_Number_of_Mutations_Within_1MB'])

    #write this dataframe
    df.to_csv(filePath, sep='\t', header=True, index=False)
####################################################################################


####################################################################################
#This method can be customized for other dictionaries
def writeCountDictionaryAsADataframe(signature2ProcessiveGroupLength2CountDict,filePath):
    L = sorted([(signature, processiveGroupLength, count,0,0) for signature, v in signature2ProcessiveGroupLength2CountDict.items() for processiveGroupLength, count in v.items()])
    df = pd.DataFrame(L,columns=['Signature', 'Processsive_Group_Length', 'Number_of_Processive_Groups', 'Median_Distance', 'Median_of_Number_of_Processive_Groups_in_MB'])

    #write this dataframe
    df.to_csv(filePath, sep='\t', header=True, index=False)
####################################################################################


####################################################################################
def readSinglePointMutationsFindProcessivityGroupsWithMultiProcessing(mutation_types_contexts,
                                                                      chromNamesList,
                                                                      computation_type,
                                                                      processivity_calculation_type,
                                                                      outputDir,
                                                                      jobname,
                                                                      numofSimulations,
                                                                      considerProbabilityInProcessivityAnalysis,
                                                                      subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                      verbose):

    if ((SBS96 in mutation_types_contexts) or (SBS384 in mutation_types_contexts) or (SBS1536 in mutation_types_contexts) or (SBS3072 in mutation_types_contexts)):
        #####################################################################
        if verbose: print('\tVerbose computation_type:%s' % (computation_type))

        if computation_type == USING_APPLY_ASYNC:
            for simNum in range(0,numofSimulations+1):
                jobs = []
                numofProcesses = multiprocessing.cpu_count()
                pool = multiprocessing.Pool(processes=numofProcesses)

                signature2ProcessiveGroupLength2DistanceListDict = {}

                ####################################################################################
                def accumulate_apply_async_result_with_distance(my_dict):
                    if my_dict is not None:
                        for (signature,processive_group_length) in my_dict:
                            my_distance_array=my_dict[(signature,processive_group_length)]

                            if signature in signature2ProcessiveGroupLength2DistanceListDict:
                                if processive_group_length in signature2ProcessiveGroupLength2DistanceListDict[signature]:
                                    signature2ProcessiveGroupLength2DistanceListDict[signature][processive_group_length].extend(my_distance_array.tolist())
                                else:
                                    signature2ProcessiveGroupLength2DistanceListDict[signature][processive_group_length] = my_distance_array.tolist()
                            else:
                                signature2ProcessiveGroupLength2DistanceListDict[signature] = {}
                                signature2ProcessiveGroupLength2DistanceListDict[signature][processive_group_length] = my_distance_array.tolist()
                ####################################################################################


                ####################################################################################
                def accumulate_apply_async_result_with_count(my_dict):
                    if my_dict is not None:
                        for (signature,processive_group_length) in my_dict:
                            count=my_dict[(signature,processive_group_length)]

                            if signature in signature2ProcessiveGroupLength2DistanceListDict:
                                if processive_group_length in signature2ProcessiveGroupLength2DistanceListDict[signature]:
                                    signature2ProcessiveGroupLength2DistanceListDict[signature][processive_group_length]+=count
                                else:
                                    signature2ProcessiveGroupLength2DistanceListDict[signature][processive_group_length] = count
                            else:
                                signature2ProcessiveGroupLength2DistanceListDict[signature] = {}
                                signature2ProcessiveGroupLength2DistanceListDict[signature][processive_group_length] = count
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

                        # Add new column. We take the probability of signature with highest probability
                        chrBased_spms_df['Probability'] = chrBased_spms_df[signatures].max(axis=1)

                        # Drop the signatures columns
                        chrBased_spms_df.drop(signatures, inplace=True, axis=1)

                        sampleBased_chrBased_spms_df_grouped = chrBased_spms_df.groupby(SAMPLE)

                        if processivity_calculation_type==CONSIDER_COUNT:
                            # With Count
                            for sample, sampleBased_chrBased_spms_df in sampleBased_chrBased_spms_df_grouped:
                                jobs.append(pool.apply_async(findProcessiveGroupsForApplySyncWithCount,
                                                             args=(
                                                             simNum,
                                                             chrLong,
                                                             sample,
                                                             sampleBased_chrBased_spms_df,
                                                             considerProbabilityInProcessivityAnalysis,
                                                             subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                             verbose,),
                                                             callback=accumulate_apply_async_result_with_count))

                        elif processivity_calculation_type==CONSIDER_DISTANCE:
                            # With Distance
                            for sample, sampleBased_chrBased_spms_df in sampleBased_chrBased_spms_df_grouped:
                                jobs.append(pool.apply_async(findProcessiveGroupsForApplySyncWithDistance,
                                                             args=(simNum,
                                                                   chrLong,
                                                                   sample,
                                                                   sampleBased_chrBased_spms_df,
                                                                   considerProbabilityInProcessivityAnalysis,
                                                                   subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                   verbose,),
                                                             callback=accumulate_apply_async_result_with_distance))
                ####################################################################################

                if verbose: print('\tVerbose Processivity len(jobs):%d\n' % (len(jobs)))

                pool.close()
                pool.join()

                if (processivity_calculation_type == CONSIDER_COUNT):
                    ####################################################################################
                    filename = 'Sim%d_Processivity.txt' % (simNum)
                    filePath = os.path.join(outputDir, jobname, DATA, PROCESSIVITY, filename)
                    writeCountDictionaryAsADataframe(signature2ProcessiveGroupLength2DistanceListDict,filePath)
                    ####################################################################################
                elif (processivity_calculation_type == CONSIDER_DISTANCE) or (processivity_calculation_type == CONSIDER_DISTANCE_ALL_SAMPLES_TOGETHER):
                    ####################################################################################
                    filename = 'Sim%d_Processivity.txt' % (simNum)
                    filePath = os.path.join(outputDir, jobname, DATA, PROCESSIVITY, filename)
                    signature2ProcessiveGroupLength2PropertiesDict = findMedians(signature2ProcessiveGroupLength2DistanceListDict)
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
def processivityAnalysis(mutation_types_contexts,
                         chromNamesList,
                         computation_type,
                         processivity_calculation_type,
                         outputDir,
                         jobname,
                         numofSimulations,
                         considerProbabilityInProcessivityAnalysis,
                         subsSignature_cutoff_numberofmutations_averageprobability_df,
                         verbose):

    print('\n#################################################################################')
    print('--- ProcessivityAnalysis starts')
    print('processivity_calculation_type:%s' % processivity_calculation_type)

    readSinglePointMutationsFindProcessivityGroupsWithMultiProcessing(mutation_types_contexts,
                                                                      chromNamesList,
                                                                      computation_type,
                                                                      processivity_calculation_type,
                                                                      outputDir,
                                                                      jobname,
                                                                      numofSimulations,
                                                                      considerProbabilityInProcessivityAnalysis,
                                                                      subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                      verbose)
    print('--- ProcessivityAnalysis ends')
    print('#################################################################################\n')
##################################################################################
