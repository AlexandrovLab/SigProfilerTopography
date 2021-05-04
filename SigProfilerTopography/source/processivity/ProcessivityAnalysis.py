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

from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC
from SigProfilerTopography.source.commons.TopographyCommons import memory_usage

from SigProfilerTopography.source.commons.TopographyCommons import CONSIDER_COUNT
from SigProfilerTopography.source.commons.TopographyCommons import CONSIDER_DISTANCE
from SigProfilerTopography.source.commons.TopographyCommons import CONSIDER_DISTANCE_ALL_SAMPLES_TOGETHER

####################################################################################
#DEC 22, 2019
def findProcessiveGroupsForApplySyncWithDistance(inter_mutational_distance_for_processivity,
                                                 simNum,
                                                 chrLong,
                                                 sample,
                                                 sampleBased_chrBased_subs_df,
                                                 considerProbabilityInProcessivityAnalysis,
                                                 subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                 verbose):
    #Sort it
    sampleBased_chrBased_subs_df.sort_values(START, inplace=True)

    return findProcessiveGroupsWithDistance(inter_mutational_distance_for_processivity,
                                            simNum,
                                            chrLong,
                                            sample,
                                            sampleBased_chrBased_subs_df,
                                            considerProbabilityInProcessivityAnalysis,
                                            subsSignature_cutoff_numberofmutations_averageprobability_df,
                                            verbose)
####################################################################################


####################################################################################
def findProcessiveGroupsWithDistance(inter_mutational_distance_for_processivity,
                                     simNum,
                                     chrLong,
                                     sample,
                                     sorted_sampleBased_chrBased_subs_df,
                                     considerProbabilityInProcessivityAnalysis,
                                     signature_cutoff_numberofmutations_averageprobability_df,
                                     verbose):

    my_dict={}
    mutations_loci_df=None

    if verbose: print('\tVerbose Worker pid %s memory_usage %.2f MB simNum:%d chrLong:%s sample:%s findProcessiveGroups starts' %(str(os.getpid()), memory_usage(),simNum,chrLong,sample))
    # They must be coming from the same sample
    # They must be same type of mutation e.g.: T>A
    # They must be resulted from same signature
    # They must be on the same strand
    if (sorted_sampleBased_chrBased_subs_df is not None):

        #As long as mutation, signature, and pyrimidine strands are the same continue to accumulate, if one of them is not the same calculate cumsum and start again
        # e.g.: If there are 4 consecutive rows with same mutation, signature and pyrimidine strand, subgroup will be 4 for these 4 rows
        if inter_mutational_distance_for_processivity:
            sorted_sampleBased_chrBased_subs_df['subgroup'] = ((sorted_sampleBased_chrBased_subs_df[MUTATION].ne(sorted_sampleBased_chrBased_subs_df[MUTATION].shift())) |
                                                            (sorted_sampleBased_chrBased_subs_df['Signature'].ne(sorted_sampleBased_chrBased_subs_df['Signature'].shift())) |
                                                            (sorted_sampleBased_chrBased_subs_df[PYRAMIDINESTRAND].ne(sorted_sampleBased_chrBased_subs_df[PYRAMIDINESTRAND].shift()))|
                                                           (sorted_sampleBased_chrBased_subs_df[START] - sorted_sampleBased_chrBased_subs_df[START].shift() > inter_mutational_distance_for_processivity)) .cumsum()
        else:
            sorted_sampleBased_chrBased_subs_df['subgroup'] = ((sorted_sampleBased_chrBased_subs_df[MUTATION].ne(sorted_sampleBased_chrBased_subs_df[MUTATION].shift())) |
                                                            (sorted_sampleBased_chrBased_subs_df['Signature'].ne(sorted_sampleBased_chrBased_subs_df['Signature'].shift())) |
                                                            (sorted_sampleBased_chrBased_subs_df[PYRAMIDINESTRAND].ne(sorted_sampleBased_chrBased_subs_df[PYRAMIDINESTRAND].shift()))) .cumsum()


        ################################
        #Former way
        # df = sorted_sampleBased_chrBased_subs_df.groupby("subgroup").agg(
        #     Signature=pd.NamedAgg(column='Signature', aggfunc="first"),
        #     ProcessiveGroupLength=pd.NamedAgg(column=MUTATION, aggfunc="count"),
        #     Probability=pd.NamedAgg(column='Probability', aggfunc="max"),
        #     LastDistance=pd.NamedAgg(column='Start', aggfunc="last"),
        #     FirstDistance=pd.NamedAgg(column='Start', aggfunc="first"),
        #     ).assign(Distance=lambda x: x.pop('LastDistance') - x.pop('FirstDistance'))
        # #agg().reset_index() can be used. not necessary

        #Current way facilitates printing mutation start positions
        df = sorted_sampleBased_chrBased_subs_df.groupby('subgroup').agg(
            {'Sample':'first',
             'Signature': 'first',
             'Probability':'max',
             'Mutation': 'count',
             'PyramidineStrand': 'first',
             # 'Start': lambda x: x.unique().tolist(),
             # 'Start': lambda x: x.tolist(),
             'Start': list,
             }).rename(columns={'Start': 'Start_List', 'Mutation':'ProcessiveGroupLength'}).reset_index()

        df['FirstLoci'] = pd.DataFrame(df['Start_List'].values.tolist()).min(axis=1)
        df['LastLoci'] = pd.DataFrame(df['Start_List'].values.tolist()).max(axis=1)

        df['Distance'] = df['LastLoci'] - df['FirstLoci']
        df.drop(columns=['FirstLoci', 'LastLoci'], inplace=True)
        ################################

        ################################
        # Remove rows with processive groups of length 1
        df = df[df['ProcessiveGroupLength'].ne(1)]
        ################################

        ################################
        if considerProbabilityInProcessivityAnalysis:
            df=df.loc[df['Probability'] >= df['Signature'].map(signature_cutoff_numberofmutations_averageprobability_df.set_index('signature')['cutoff'])]
        ################################

        #'df columns: subgroup' 'Signature' 'Probability' 'ProcessiveGroupLength' 'PyramidineStrand' 'Start_List' 'Distance'
        mutations_loci_df = df.groupby(['Signature', 'ProcessiveGroupLength']).agg(
            {'Sample':'first',
             'Start_List': list,
             'Distance':list
             }).reset_index()

        mutations_loci_df['Chr']=chrLong
        mutations_loci_df = mutations_loci_df[['Sample', 'Chr', 'Signature', 'ProcessiveGroupLength', 'Start_List', 'Distance']]

        ################################
        my_dict = {k: np.hstack(v) for k, v in df.groupby(['Signature', 'ProcessiveGroupLength'])['Distance']}
        ################################

    if verbose: print('\tVerbose Worker pid %s memory_usage %.2f MB simNum:%d chrLong:%s sample:%s findProcessiveGroups ends' %(str(os.getpid()), memory_usage(),simNum,chrLong,sample))

    return (my_dict, mutations_loci_df)
    # return (my_dict,)
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
                 v1['median_number_of_mutations_within_1MB'])
                for signature, v in signature2ProcessiveGroupLength2PropertiesDict.items()
                for processiveGroupLength, v1 in v.items()])

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
def readSinglePointMutationsFindProcessivityGroupsWithMultiProcessing(mutation_types_contexts,
                                                                      chromNamesList,
                                                                      computation_type,
                                                                      inter_mutational_distance_for_processivity,
                                                                      outputDir,
                                                                      jobname,
                                                                      numofSimulations,
                                                                      considerProbabilityInProcessivityAnalysis,
                                                                      subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                      verbose):

    if ((SBS96 in mutation_types_contexts) or (SBS384 in mutation_types_contexts) or (SBS1536 in mutation_types_contexts) or (SBS3072 in mutation_types_contexts)):
        if verbose: print('\tVerbose computation_type:%s ' % (computation_type))

        for simNum in range(0,numofSimulations+1):
            jobs = []
            numofProcesses = multiprocessing.cpu_count()
            pool = multiprocessing.Pool(processes=numofProcesses)

            signature2ProcessiveGroupLength2DistanceListDict = {}
            mutations_loci_df_list=[]

            def accumulate_apply_async_result_with_distance(result_tuple):
                my_dict=result_tuple[0]
                mutations_loci_df=result_tuple[1]

                if (mutations_loci_df is not None) and (not mutations_loci_df.empty):
                    mutations_loci_df_list.append(mutations_loci_df)
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

            for chrLong in chromNamesList:
                chrBased_subs_df = readChrBasedMutationsDF(outputDir, jobname, chrLong, SUBS, simNum)
                if (chrBased_subs_df is not None):
                    ###########################################
                    chrBased_subs_df.drop([SIMULATION_NUMBER], inplace=True, errors='ignore', axis=1)
                    columnNamesList = list(chrBased_subs_df.columns.values)
                    # We assume that after the column named 'Context' there are the signature columns in tab separated way.
                    mutationIndex = columnNamesList.index(MUTATION)
                    signatures = columnNamesList[(mutationIndex + 1):]
                    ###########################################
                    chrBased_subs_df.drop([CHROM, MUTATIONLONG, TRANSCRIPTIONSTRAND],inplace=True, errors='ignore', axis=1)
                    chrBased_subs_df['Signature'] = (chrBased_subs_df.loc[:, signatures]).idxmax(axis=1)

                    # Add new column. We take the probability of signature with highest probability
                    chrBased_subs_df['Probability'] = chrBased_subs_df[signatures].max(axis=1)

                    # Drop the signatures columns
                    chrBased_subs_df.drop(signatures, inplace=True, axis=1)

                    sampleBased_chrBased_subs_df_grouped = chrBased_subs_df.groupby(SAMPLE)

                    # With Distance
                    for sample, sampleBased_chrBased_subs_df in sampleBased_chrBased_subs_df_grouped:
                        # Parallel version for real run
                        jobs.append(pool.apply_async(findProcessiveGroupsForApplySyncWithDistance,
                                                     args=(inter_mutational_distance_for_processivity,
                                                           simNum,
                                                           chrLong,
                                                           sample,
                                                           sampleBased_chrBased_subs_df,
                                                           considerProbabilityInProcessivityAnalysis,
                                                           subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                           verbose,),
                                                     callback=accumulate_apply_async_result_with_distance))

                        # Sequential version for testing/debugging
                        # result_tuple = findProcessiveGroupsForApplySyncWithDistance(inter_mutational_distance_for_processivity,
                        #                            simNum,
                        #                            chrLong,
                        #                            sample,
                        #                            sampleBased_chrBased_subs_df,
                        #                            considerProbabilityInProcessivityAnalysis,
                        #                            subsSignature_cutoff_numberofmutations_averageprobability_df,
                        #                            verbose)
                        # accumulate_apply_async_result_with_distance(result_tuple)

            if verbose: print('\tVerbose Processivity len(jobs):%d\n' % (len(jobs)))

            pool.close()
            pool.join()

            if len(mutations_loci_df_list)>0:
                all_mutations_loci_df=pd.concat(mutations_loci_df_list)
                all_mutations_loci_df.to_csv(os.path.join(outputDir, jobname, DATA, PROCESSIVITY, "Sim%d_Processive_Mutations_Loci.txt"%(simNum)), sep='\t', index=False, header=True)

            filename = 'Sim%d_Processivity.txt' % (simNum)
            filePath = os.path.join(outputDir, jobname, DATA, PROCESSIVITY, filename)
            signature2ProcessiveGroupLength2PropertiesDict = findMedians(signature2ProcessiveGroupLength2DistanceListDict)
            writeDictionaryAsADataframe(signature2ProcessiveGroupLength2PropertiesDict,filePath)
####################################################################################


##################################################################################
#Default value for processivity_calculation_type=CONSIDER_DISTANCE
def processivityAnalysis(mutation_types_contexts,
                         chromNamesList,
                         computation_type,
                         processivity_calculation_type,
                         inter_mutational_distance_for_processivity,
                         outputDir,
                         jobname,
                         numofSimulations,
                         considerProbabilityInProcessivityAnalysis,
                         subsSignature_cutoff_numberofmutations_averageprobability_df,
                         verbose):

    print('\n#################################################################################')
    print('--- ProcessivityAnalysis starts')
    print('processivity_calculation_type:%s' % processivity_calculation_type)
    print('inter_mutational_distance_for_processivity:%s' % inter_mutational_distance_for_processivity)

    readSinglePointMutationsFindProcessivityGroupsWithMultiProcessing(mutation_types_contexts,
                                                                      chromNamesList,
                                                                      computation_type,
                                                                      inter_mutational_distance_for_processivity,
                                                                      outputDir,
                                                                      jobname,
                                                                      numofSimulations,
                                                                      considerProbabilityInProcessivityAnalysis,
                                                                      subsSignature_cutoff_numberofmutations_averageprobability_df,
                                                                      verbose)
    print('--- ProcessivityAnalysis ends')
    print('#################################################################################\n')
##################################################################################