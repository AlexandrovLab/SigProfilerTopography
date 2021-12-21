# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018-2020 Burcak Otlu


# Version2
# This version use np.arrays
# Right now replication strand bias analysis works for single point mutations and signatures.
# This python code analyses the Replication Strand Bias

import multiprocessing
import numpy as np
import pandas as pd
import os

from SigProfilerTopography.source.commons.TopographyCommons import CHROM
from SigProfilerTopography.source.commons.TopographyCommons import START
from SigProfilerTopography.source.commons.TopographyCommons import END
from SigProfilerTopography.source.commons.TopographyCommons import SIGNAL

from SigProfilerTopography.source.commons.TopographyCommons import PYRAMIDINESTRAND
from SigProfilerTopography.source.commons.TopographyCommons import SAMPLE

from SigProfilerTopography.source.commons.TopographyCommons import TYPE

from SigProfilerTopography.source.commons.TopographyCommons import SUBS
from SigProfilerTopography.source.commons.TopographyCommons import INDELS
from SigProfilerTopography.source.commons.TopographyCommons import DINUCS

from SigProfilerTopography.source.commons.TopographyCommons import MUTATION
from SigProfilerTopography.source.commons.TopographyCommons import MUTATION_LONG
from SigProfilerTopography.source.commons.TopographyCommons import LENGTH

from SigProfilerTopography.source.commons.TopographyCommons import LEADING
from SigProfilerTopography.source.commons.TopographyCommons import LAGGING
from SigProfilerTopography.source.commons.TopographyCommons import REPLICATIONSTRANDBIAS

from SigProfilerTopography.source.commons.TopographyCommons import DATA
from SigProfilerTopography.source.commons.TopographyCommons import LIB
from SigProfilerTopography.source.commons.TopographyCommons import CHRBASED

from SigProfilerTopography.source.commons.TopographyCommons import readWig_with_fixedStep_variableStep

from SigProfilerTopography.source.occupancy.ChrBasedSignalArrays import readFileInBEDFormat
from SigProfilerTopography.source.commons.TopographyCommons import memory_usage

from SigProfilerTopography.source.commons.TopographyCommons import write_type_strand_bias_np_array_as_dataframe
from SigProfilerTopography.source.commons.TopographyCommons import write_signature_mutation_type_strand_bias_np_array_as_dataframe
from SigProfilerTopography.source.commons.TopographyCommons import write_sbs_signature_sbs96_mutation_type_replication_strand_bias
from SigProfilerTopography.source.commons.TopographyCommons import write_sample_based_strand1_strand2_as_dataframe

from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM
from SigProfilerTopography.source.commons.TopographyCommons import USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT

from SigProfilerTopography.source.commons.TopographyCommons import get_chrBased_simBased_combined_df_split
from SigProfilerTopography.source.commons.TopographyCommons import get_chrBased_simBased_combined_df
from SigProfilerTopography.source.commons.TopographyCommons import get_chrBased_simBased_dfs

from SigProfilerTopography.source.commons.TopographyCommons import decideFileType

from SigProfilerTopography.source.commons.TopographyCommons import C2A
from SigProfilerTopography.source.commons.TopographyCommons import C2G
from SigProfilerTopography.source.commons.TopographyCommons import C2T
from SigProfilerTopography.source.commons.TopographyCommons import T2A
from SigProfilerTopography.source.commons.TopographyCommons import T2C
from SigProfilerTopography.source.commons.TopographyCommons import T2G

#For Supp Fig2B
CHR10_THRESHOLD_START = 16400000
CHR10_THRESHOLD_END = 26400000

#For Supp Fig2A
CHR20_START = 36260000
CHR20_END = 36830000

# FOR FINDING TRANSITION ZONES (LEADING or LAGGING)
# THRESHOLD_CONSECUTIVE_LONG_STRETCH_LENGTH= 250000 #used in Supp Fig2B
# THRESHOLD_CONSECUTIVE_LONG_STRETCH_LENGTH= 150000
THRESHOLD_CONSECUTIVE_LONG_STRETCH_LENGTH= 10000

# THRESHOLD_DISCARD_LATEST_TRANSITION_ZONE = 100000 #used in Supp Fig2B
THRESHOLD_DISCARD_LATEST_TRANSITION_ZONE = 25000
# THRESHOLD_LATEST_TRANSITION_ZONE = 0


########################################################################
def checkForSameSignedSlopeBetweenConsecutivePeakandValley(chrLong,peakorValleyStart, peakorValleyEnd, chrBasedSmoothedWaveletReplicationTimeSignalDF):
    transitionZoneList =[]

    # print('################ checkForConsecutive starts ############ fromStart: %s toEnd: %s' %(peakorValleyStart,peakorValleyEnd))
    subset_df = chrBasedSmoothedWaveletReplicationTimeSignalDF[(chrBasedSmoothedWaveletReplicationTimeSignalDF[START]>=peakorValleyStart) & (chrBasedSmoothedWaveletReplicationTimeSignalDF[END]<=peakorValleyEnd)]

    consecutiveLength = 0
    formerRow= None
    formerSlopeDirection = None

    start = peakorValleyStart

    for index,row in subset_df.iterrows():
        if formerRow is None:
            #We read the row for the first time
            formerRow = row
            consecutiveLength += 1000
        else:
            slope = (row.get(SIGNAL) - formerRow.get(SIGNAL)) / 1000
            formerRow = row

            if (formerSlopeDirection is None):
                formerSlopeDirection = np.sign(slope)
                consecutiveLength += 1000
            elif (formerSlopeDirection==np.sign(slope)):
                consecutiveLength += 1000
            else:
                #They have different signs
                if (consecutiveLength>=THRESHOLD_CONSECUTIVE_LONG_STRETCH_LENGTH):
                    # print('Slope sign changed -- Found one: from %d to %d with %d bases with slope sign %s' %(start,((row.get('start') + row.get('end'))//2), consecutiveLength, formerSlopeDirection))
                    transitionZoneList.append((chrLong,start,(row.get(START) + row.get(END))//2,formerSlopeDirection,consecutiveLength))
                #initialize and start again
                consecutiveLength = 1000
                start = (row.get(START) + row.get(END))//2
                formerRow= row
                formerSlopeDirection= np.sign(slope)
                continue

            # print('slope: %f - np.sign(slope): %f -  consecutiveLength: %d ' %(slope,np.sign(slope),consecutiveLength))
            formerSlopeDirection = np.sign(slope)

    #This is for the last probable transition zone.
    if (consecutiveLength >= THRESHOLD_CONSECUTIVE_LONG_STRETCH_LENGTH):
        # print('After for loop ends, found one: from %d to %s with %d bases with slope sign %s' % (start, (row.get('start') + row.get('end'))//2, consecutiveLength, formerSlopeDirection))
        transitionZoneList.append((chrLong,start,(row.get(START) + row.get(END))//2,formerSlopeDirection,consecutiveLength))

    # print('################ checkForConsecutive ends ############ fromStart: %s toEnd: %s' % (peakorValleyStart,peakorValleyEnd))
    return transitionZoneList
########################################################################


########################################################################
# chr10_subset_wavelet_processed_df
#           chr     start       end   signal
# 265577  chr10  16400500  16401499  24.9438
# 265578  chr10  16401500  16402499  24.9585

# valleys_peaks_df
#         chr     start       end    type
# 415     chr10  16454500  16455500    Peak
# 415  chr10  16528500  16529500  Valley

def findLongStretchesofConsistentTransitionZones(chrLong,fromStart,toEnd,chrBasedSmoothedWaveletReplicationTimeSignalDF,valleys_peaks_df):
    transitionZonesList =[]
    for index,row in  valleys_peaks_df.iterrows():
        peakorValleyStart = row[START]
        peakorValleyEnd = row[END]
        peakorValleyMidpoint = (peakorValleyStart+peakorValleyEnd)//2

        type = row['type']
        if (type =='Peak'):
            if (peakorValleyMidpoint>fromStart):
                # print('from: %d - to: %d - difference: %d'  %(fromStart,peakorValleyMidpoint, (peakorValleyMidpoint-fromStart)))
                found = checkForSameSignedSlopeBetweenConsecutivePeakandValley(chrLong,fromStart, peakorValleyMidpoint, chrBasedSmoothedWaveletReplicationTimeSignalDF)
                transitionZonesList.extend(found)
                # print('found %s' %found)
            fromStart=peakorValleyMidpoint
        elif (type=='Valley'):
            valleyStart =row[START]
            valleyEnd = row[END]
            valleyMidpoint = (valleyStart+valleyEnd)//2
            # This is something special to valley
            newValleyStart1 = valleyMidpoint - THRESHOLD_DISCARD_LATEST_TRANSITION_ZONE
            newValleyStart2 = valleyMidpoint + THRESHOLD_DISCARD_LATEST_TRANSITION_ZONE
            if (newValleyStart1>fromStart):
                # print('from: %d - to: %d - difference: %d' % (fromStart, newValleyStart1, (newValleyStart1 - fromStart)))
                found = checkForSameSignedSlopeBetweenConsecutivePeakandValley(chrLong,fromStart, newValleyStart1,chrBasedSmoothedWaveletReplicationTimeSignalDF)
                transitionZonesList.extend(found)
                # print('found %s' % found)
            # bypass the genome region between newValleyStart1 and newValleyStart2
            fromStart = newValleyStart2
    #
    #For the last interval
    if (toEnd>fromStart):
        # print('last one from: %d - to: %d -difference: %d' %(fromStart,toEnd,(toEnd-fromStart)))
        found = checkForSameSignedSlopeBetweenConsecutivePeakandValley(chrLong,fromStart, toEnd, chrBasedSmoothedWaveletReplicationTimeSignalDF)
        transitionZonesList.extend(found)
        # print('found %s' %found)

    return transitionZonesList
########################################################################


########################################################################
#TODO Is (replicationStrand_row['end']+1) okey?
# We assume that there are no overlapping intervals with positive and negative slopes.
# To test it have one array for positive slope fill with 1s
#                one array for negative slope fill with -2a
#                add them if you habe any -1 that means that you contradict this assumption.
def fillReplicationStrandArray(replicationStrand_row,chrBased_replication_array):
    # e.g.: replicationStrand_row
    # chr chrX
    # start   154861998
    # end 155096999
    # slopeDirection  1 (1 means leading strand -1 means lagging strand on positive strand)
    # length  235000

    # labels = ['chr', 'start', 'end', 'slopeDirection', 'length']
    chrBased_replication_array[replicationStrand_row['start']:replicationStrand_row['end']+1] = replicationStrand_row['slopeDirection']
########################################################################


########################################################################
# July 28, 2020
# Using numpy arrays
#   if mutationPyramidineStrand and slope have the same sign increase LEADING STRAND count
#   else mutationPyramidineStrand and slope have the opposite sign increase LAGGING STRAND count
# sample_based for further usage
def searchAllMutationOnReplicationStrandArray_using_list_comprehension_using_numpy_array(
        mutation_row,
        my_type,
        chrBasedReplicationArray,
        SBS6_mutation_types_np_array,
        SBS96_mutation_types_np_array,
        ordered_signatures_cutoffs,
        df_columns_signatures_mask_array,
        SBS6_mutation_types_default_zeros_array,
        SBS96_mutation_types_default_zeros_array,
        subs_signatures_default_zeros_array,
        dinucs_signatures_default_zeros_array,
        indels_signatures_default_zeros_array,
        subs_signatures_SBS6_mutation_types_default_zeros_array,
        subs_signatures_SBS96_mutation_types_default_zeros_array,
        all_types_leading_np_array,
        all_types_lagging_np_array,
        subs_signature_SBS6_mutation_type_leading_np_array,
        subs_signature_SBS6_mutation_type_lagging_np_array,
        subs_signature_SBS96_mutation_type_leading_np_array,
        subs_signature_SBS96_mutation_type_lagging_np_array,
        all_samples_all_types_leading_np_array,
        all_samples_all_types_lagging_np_array,
        all_samples_subs_signature_mutation_type_leading_np_array,
        all_samples_subs_signature_mutation_type_lagging_np_array,
        sample_based,
        all_samples_np_array,
        is_discreet,
        df_columns):

    if sample_based:
        indexofSample = np.where(df_columns == SAMPLE)[0][0]
        mutation_sample = mutation_row[indexofSample]
        sample_index=np.where(all_samples_np_array == mutation_sample)[0][0]

    indexofStart = np.where(df_columns == START)[0][0]
    start = mutation_row[indexofStart]

    indexofPyrimidineStrand = np.where(df_columns == PYRAMIDINESTRAND)[0][0]
    pyramidineStrand = mutation_row[indexofPyrimidineStrand]

    subs_signature_SBS6_mutation_type_mask_array = subs_signatures_SBS6_mutation_types_default_zeros_array
    subs_signature_SBS96_mutation_type_mask_array = subs_signatures_SBS96_mutation_types_default_zeros_array

    probabilities = mutation_row[df_columns_signatures_mask_array]

    #############################################################################################################
    if(my_type==SUBS):
        end = start+1
        #e.g.: C>A

        indexofMutation = np.where(df_columns == MUTATION)[0][0]
        SBS6_mutation_type = mutation_row[indexofMutation]

        index_of_mutation_long = np.where(df_columns == MUTATION_LONG)[0][0]

        # e.g.: T:AA[C>A]AA
        mutation_type_long = mutation_row[index_of_mutation_long]
        SBS96_mutation_type = mutation_type_long[3:10]

        # six_mutation_types_mask_array.shape (6,)
        SBS6_mutation_types_mask_array = np.where(SBS6_mutation_types_np_array == SBS6_mutation_type, 1, 0)
        SBS96_mutation_types_mask_array = np.where(SBS96_mutation_types_np_array == SBS96_mutation_type, 1, 0)

        if is_discreet:
            # Convert True into 1, and False into 0
            # subs_signatures_mask_array.shape (num_of_subs_signatures,)
            threshold_mask_array = np.greater_equal(probabilities, ordered_signatures_cutoffs)
            subs_signatures_mask_array = threshold_mask_array.astype(int)
        else:
            subs_signatures_mask_array = np.array(probabilities).astype(float)

        #Concetanate
        all_types_mask_array = np.concatenate((SBS6_mutation_types_mask_array,
                                              SBS96_mutation_types_mask_array,
                                              subs_signatures_mask_array, # SUBS
                                              dinucs_signatures_default_zeros_array,
                                              indels_signatures_default_zeros_array), axis=None)

        # multiply subs_signatures_mask_array times six_mutation_types_mask_array
        # Add one more dimension to subs_signatures_mask_array and six_mutation_types_mask_array

        # subs_signatures_mask_array_2d.shape (1,num_of_subs_signatures)
        subs_signatures_mask_array_2d = np.expand_dims(subs_signatures_mask_array, axis=0)

        # six_mutation_types_mask_array_2d.shape (1,6)
        SBS6_mutation_types_mask_array_2d = np.expand_dims(SBS6_mutation_types_mask_array, axis=0)

        # SBS96_mutation_types_mask_array.shape (96,) --> (1,96)
        SBS96_mutation_types_mask_array_2d = np.expand_dims(SBS96_mutation_types_mask_array, axis=0)

        # to_be_accumulated_array.shape (num_of_subs_signatures,6)
        subs_signature_SBS6_mutation_type_mask_array = subs_signatures_mask_array_2d.T * SBS6_mutation_types_mask_array_2d

        # to_be_accumulated_array.shape (num_of_subs_signatures,96)
        subs_signature_SBS96_mutation_type_mask_array = subs_signatures_mask_array_2d.T * SBS96_mutation_types_mask_array_2d

    elif (my_type==DINUCS):
        end = start+2

        if is_discreet:
            # Convert True into 1, and False into 0
            threshold_mask_array = np.greater_equal(probabilities, ordered_signatures_cutoffs)
            dinucs_signatures_mask_array = threshold_mask_array.astype(int)
        else:
            dinucs_signatures_mask_array = np.array(probabilities).astype(float)

        #Concetanate
        all_types_mask_array= np.concatenate((SBS6_mutation_types_default_zeros_array,
                                              SBS96_mutation_types_default_zeros_array,
                                              subs_signatures_default_zeros_array,
                                              dinucs_signatures_mask_array, # DINUCS
                                              indels_signatures_default_zeros_array), axis=None)

    elif (my_type==INDELS):
        indexofLength = np.where(df_columns == LENGTH)[0][0]
        end = start+int(mutation_row[indexofLength])

        if is_discreet:
            # Convert True into 1, and False into 0
            threshold_mask_array = np.greater_equal(probabilities, ordered_signatures_cutoffs)
            indels_signatures_mask_array = threshold_mask_array.astype(int)
        else:
            indels_signatures_mask_array = np.array(probabilities).astype(float)

        #Concetanate
        all_types_mask_array= np.concatenate((SBS6_mutation_types_default_zeros_array,
                                              SBS96_mutation_types_default_zeros_array,
                                              subs_signatures_default_zeros_array,
                                              dinucs_signatures_default_zeros_array,
                                              indels_signatures_mask_array), axis=None)
    #############################################################################################################

    #############################################################################################################
    #if there is overlap with chrBasedReplicationArray
    slicedArray = chrBasedReplicationArray[int(start):int(end)]

    if (np.any(slicedArray)):
        #It must be full with at most -1 and +1
        uniqueValueArray = np.unique(slicedArray[np.nonzero(slicedArray)])

        # I expect the value of 1 (LEADING on the positive strand) or -1 (LAGGING on the positive strand) so size must be one.
        if (uniqueValueArray.size == 1):
            for uniqueValue in np.nditer(uniqueValueArray):
                # type(decileIndex) is numpy.ndarray
                slope = int(uniqueValue)
                #They have the same sign, multiplication (1,1) (-1,-1) must be 1
                if (slope*pyramidineStrand > 0):
                    all_types_leading_np_array += all_types_mask_array
                    subs_signature_SBS6_mutation_type_leading_np_array += subs_signature_SBS6_mutation_type_mask_array
                    subs_signature_SBS96_mutation_type_leading_np_array += subs_signature_SBS96_mutation_type_mask_array
                # They have the opposite sign, multiplication(1,-1) (-1,-)  must be -1
                elif (slope*pyramidineStrand < 0):
                    all_types_lagging_np_array += all_types_mask_array
                    subs_signature_SBS6_mutation_type_lagging_np_array += subs_signature_SBS6_mutation_type_mask_array
                    subs_signature_SBS96_mutation_type_lagging_np_array += subs_signature_SBS96_mutation_type_mask_array

        elif ((uniqueValueArray.size==2) and (pyramidineStrand!=0)):
            #Increment both LEADING and LAGGING
            all_types_leading_np_array += all_types_mask_array
            all_types_lagging_np_array += all_types_mask_array
            subs_signature_SBS6_mutation_type_leading_np_array += subs_signature_SBS6_mutation_type_mask_array
            subs_signature_SBS6_mutation_type_lagging_np_array += subs_signature_SBS6_mutation_type_mask_array
            subs_signature_SBS96_mutation_type_leading_np_array += subs_signature_SBS96_mutation_type_mask_array
            subs_signature_SBS96_mutation_type_lagging_np_array += subs_signature_SBS96_mutation_type_mask_array
        elif (uniqueValueArray.size>2):
            print('There is a situation!!!')
        else:
            print('There is a situation!!!')

    if sample_based:
        if (np.any(slicedArray)):
            # It must be full with at most -1 and +1
            uniqueValueArray = np.unique(slicedArray[np.nonzero(slicedArray)])

            # I expect the value of 1 (LEADING on the positive strand) or -1 (LAGGING on the positive strand) so size must be one.
            if (uniqueValueArray.size == 1):
                for uniqueValue in np.nditer(uniqueValueArray):
                    # type(decileIndex) is numpy.ndarray
                    slope = int(uniqueValue)
                    # They have the same sign, multiplication (1,1) (-1,-1) must be 1
                    if (slope * pyramidineStrand > 0):
                        all_samples_all_types_leading_np_array[sample_index] += all_types_mask_array
                        all_samples_subs_signature_mutation_type_leading_np_array[sample_index] += subs_signature_SBS6_mutation_type_mask_array
                    # They have the opposite sign, multiplication(1,-1) (-1,-)  must be -1
                    elif (slope * pyramidineStrand < 0):
                        all_samples_all_types_lagging_np_array[sample_index] += all_types_mask_array
                        all_samples_subs_signature_mutation_type_lagging_np_array[sample_index] += subs_signature_SBS6_mutation_type_mask_array
            elif ((uniqueValueArray.size == 2) and (pyramidineStrand != 0)):
                # Increment both LEADING and LAGGING
                all_samples_all_types_leading_np_array[sample_index] += all_types_mask_array
                all_samples_all_types_lagging_np_array[sample_index] += all_types_mask_array
                all_samples_subs_signature_mutation_type_leading_np_array[sample_index] += subs_signature_SBS6_mutation_type_mask_array
                all_samples_subs_signature_mutation_type_lagging_np_array[sample_index] += subs_signature_SBS6_mutation_type_mask_array
            elif (uniqueValueArray.size > 2):
                print('There is a situation!!!')
            else:
                print('There is a situation!!!')
    #############################################################################################################

########################################################################


########################################################################
# For df split
# July 28, 2020
# Using numpy arrays
#   if mutationPyramidineStrand and slope have the same sign increase LEADING STRAND count
#   else mutationPyramidineStrand and slope have the opposite sign increase LAGGING STRAND count
# sample_based for further usage
def searchAllMutationOnReplicationStrandArray_using_list_comprehension_using_numpy_array_for_df_split(
        mutation_row,
        chrBasedReplicationArray,
        six_mutation_types_np_array,
        ordered_sbs_signatures_cutoffs,
        ordered_dbs_signatures_cutoffs,
        ordered_id_signatures_cutoffs,
        df_columns_subs_signatures_mask_array,
        df_columns_dinucs_signatures_mask_array,
        df_columns_indels_signatures_mask_array,
        six_mutation_types_default_zeros_array,
        subs_signatures_default_zeros_array,
        dinucs_signatures_default_zeros_array,
        indels_signatures_default_zeros_array,
        subs_signatures_mutation_types_default_zeros_array,
        all_types_leading_np_array,
        all_types_lagging_np_array,
        subs_signature_mutation_type_leading_np_array,
        subs_signature_mutation_type_lagging_np_array,
        sample_based,
        is_discreet,
        df_columns):

    indexofStart = np.where(df_columns == START)[0][0]
    indexofPyrimidineStrand = np.where(df_columns == PYRAMIDINESTRAND)[0][0]
    indexofSample = np.where(df_columns == SAMPLE)[0][0]
    indexofType = np.where(df_columns == TYPE)[0][0]

    start = mutation_row[indexofStart]
    pyramidineStrand = mutation_row[indexofPyrimidineStrand]
    sample = mutation_row[indexofSample]
    my_type=mutation_row[indexofType]

    mutationType = None
    subs_signature_mutation_type_mask_array=subs_signatures_mutation_types_default_zeros_array

    #############################################################################################################
    if(my_type==SUBS):
        end = start+1
        #e.g.: C>A
        indexofMutation = np.where(df_columns == MUTATION)[0][0]
        mutationType = mutation_row[indexofMutation]

        #six_mutation_types_mask_array.shape (6,)
        six_mutation_types_mask_array= np.where(six_mutation_types_np_array == mutationType, 1, 0)

        probabilities = mutation_row[df_columns_subs_signatures_mask_array]

        if is_discreet:
            threshold_mask_array = np.greater_equal(probabilities, ordered_sbs_signatures_cutoffs)
            # Convert True into 1, and False into 0
            # subs_signatures_mask_array.shape (num_of_subs_signatures,)
            subs_signatures_mask_array = threshold_mask_array.astype(int)
        else:
            subs_signatures_mask_array = np.array(probabilities).astype(float)

        #Concetanate
        all_types_mask_array= np.concatenate((six_mutation_types_mask_array,
                                              subs_signatures_mask_array,
                                              dinucs_signatures_default_zeros_array,
                                              indels_signatures_default_zeros_array), axis=None)

        # Add one more dimension to subs_signatures_mask_array and six_mutation_types_mask_array
        # subs_signatures_mask_array_2d.shape (1,num_of_subs_signatures)
        subs_signatures_mask_array_2d = np.expand_dims(subs_signatures_mask_array, axis=0)

        # six_mutation_types_mask_array_2d.shape (1,6)
        six_mutation_types_mask_array_2d = np.expand_dims(six_mutation_types_mask_array, axis=0)

        # multiply subs_signatures_mask_array times six_mutation_types_mask_array
        subs_signature_mutation_type_mask_array = subs_signatures_mask_array_2d.T * six_mutation_types_mask_array_2d

    elif (my_type==DINUCS):
        end = start+2

        probabilities = mutation_row[df_columns_dinucs_signatures_mask_array]

        if is_discreet:
            threshold_mask_array = np.greater_equal(probabilities, ordered_dbs_signatures_cutoffs)
            # Convert True into 1, and False into 0
            dinucs_signatures_mask_array = threshold_mask_array.astype(int)
        else:
            dinucs_signatures_mask_array = np.array(probabilities).astype(float)

        #Concetanate
        all_types_mask_array= np.concatenate((six_mutation_types_default_zeros_array,
                                              subs_signatures_default_zeros_array,
                                              dinucs_signatures_mask_array,
                                              indels_signatures_default_zeros_array), axis=None)

    elif (my_type==INDELS):
        indexofLength = np.where(df_columns == LENGTH)[0][0]
        end = start+int(mutation_row[indexofLength])

        probabilities = mutation_row[df_columns_indels_signatures_mask_array]

        if is_discreet:
            threshold_mask_array = np.greater_equal(probabilities, ordered_id_signatures_cutoffs)
            # Convert True into 1, and False into 0
            indels_signatures_mask_array = threshold_mask_array.astype(int)
        else:
            indels_signatures_mask_array = np.array(probabilities).astype(float)

        #Concetanate
        all_types_mask_array= np.concatenate((six_mutation_types_default_zeros_array,
                                              subs_signatures_default_zeros_array,
                                              dinucs_signatures_default_zeros_array,
                                              indels_signatures_mask_array), axis=None)

    #############################################################################################################

    #############################################################################################################
    #if there is overlap with chrBasedReplicationArray
    slicedArray = chrBasedReplicationArray[int(start):int(end)]

    if (np.any(slicedArray)):
        #It must be full with at most -1 and +1
        uniqueValueArray = np.unique(slicedArray[np.nonzero(slicedArray)])

        # I expect the value of 1 (LEADING on the positive strand) or -1 (LAGGING on the positive strand) so size must be one.
        if (uniqueValueArray.size == 1):
            for uniqueValue in np.nditer(uniqueValueArray):
                # type(decileIndex) is numpy.ndarray
                slope = int(uniqueValue)
                #They have the same sign, multiplication (1,1) (-1,-1) must be 1
                if (slope*pyramidineStrand > 0):
                    all_types_leading_np_array += all_types_mask_array
                    subs_signature_mutation_type_leading_np_array += subs_signature_mutation_type_mask_array
                # They have the opposite sign, multiplication(1,-1) (-1,-)  must be -1
                elif (slope*pyramidineStrand < 0):
                    all_types_lagging_np_array += all_types_mask_array
                    subs_signature_mutation_type_lagging_np_array += subs_signature_mutation_type_mask_array
        elif ((uniqueValueArray.size==2) and (pyramidineStrand!=0)):
            #Increment both LEADING and LAGGING
            all_types_leading_np_array += all_types_mask_array
            all_types_lagging_np_array += all_types_mask_array
            subs_signature_mutation_type_leading_np_array += subs_signature_mutation_type_mask_array
            subs_signature_mutation_type_lagging_np_array += subs_signature_mutation_type_mask_array
        elif (uniqueValueArray.size>2):
            print('There is a situation!!!')
        else:
            print('There is a situation!!!')
    #############################################################################################################

########################################################################




########################################################################
#This code checks whether valleys and peaks are one after another, not two consecutive elements are both valley and peak.
def checkforValidness(chrBased_valleys_peaks_df):
    formerRowType = None

    for index, row in chrBased_valleys_peaks_df.iterrows():
        if formerRowType is None:
            formerRowType = row['type']
        elif (row['type']== formerRowType):
            return False
        else:
            formerRowType = row['type']

    return True
########################################################################


########################################################################
def get_chr_based_replication_strand_array_for_callback(chrLong,chromSize,repliseq_signal_df,valleys_df,peaks_df):
    chrBased_replication_array = get_chr_based_replication_strand_array(chrLong, chromSize, repliseq_signal_df, valleys_df, peaks_df)
    return (chrLong,chrBased_replication_array)
########################################################################


########################################################################
def get_chr_based_replication_strand_array(chrLong,chromSize,repliseq_signal_df,valleys_df,peaks_df):

    # Read chrBasedSmoothedWaveletReplicationTimeSignalDF
    chrBased_SmoothedWaveletReplicationTimeSignal_df = repliseq_signal_df[repliseq_signal_df[CHROM] == chrLong]

    chrBasedValleysDF = valleys_df[valleys_df[CHROM] == chrLong].copy()
    chrBasedValleysDF['type'] = 'Valley'
    chrBasedValleysDF.astype(dtype={START: int, END: int})

    chrBasedPeaksDF = peaks_df[peaks_df[CHROM] == chrLong].copy()
    chrBasedPeaksDF['type'] = 'Peak'
    chrBasedPeaksDF.astype(dtype={START: int, END: int})

    # Concat Peaks and Valleys
    chrBased_valleys_peaks_df = pd.concat([chrBasedValleysDF, chrBasedPeaksDF], axis=0)

    # Sort Valleys and peaks
    chrBased_valleys_peaks_df.sort_values(START, inplace=True)

    if ((chrBased_SmoothedWaveletReplicationTimeSignal_df is not None) and (not chrBased_SmoothedWaveletReplicationTimeSignal_df.empty) and (checkforValidness(chrBased_valleys_peaks_df))):
        chrBased_replication_array = fill_chr_based_replication_strand_array(chrLong,
                                                                             chromSize,
                                                                             chrBased_SmoothedWaveletReplicationTimeSignal_df,
                                                                             chrBased_valleys_peaks_df)
        return chrBased_replication_array
    else:
        return None
########################################################################


########################################################################
def fill_chr_based_replication_strand_array(chrLong,
                                            chromSize,
                                            chrBasedSmoothedWaveletReplicationTimeSignalDF,
                                            chrBased_valleys_peaks_df):
    # +1 means leading strand, -1 means lagging strand
    # we will fill this array using smoothedSignal, peaks and valleys for each chromosome
    chrBased_replication_array = np.zeros(chromSize, dtype=np.int8)

    firstIndex = chrBasedSmoothedWaveletReplicationTimeSignalDF.index[0]
    lastIndex = chrBasedSmoothedWaveletReplicationTimeSignalDF.index[-1]

    startColumnIndex = chrBasedSmoothedWaveletReplicationTimeSignalDF.columns.get_loc(START)
    endColumnIndex = chrBasedSmoothedWaveletReplicationTimeSignalDF.columns.get_loc(END)

    start = chrBasedSmoothedWaveletReplicationTimeSignalDF.iloc[0, startColumnIndex]  # get the first row start
    end = chrBasedSmoothedWaveletReplicationTimeSignalDF.iloc[-1, endColumnIndex]  # get the last row end

    # Step1 Find the transition zones
    chrBasedTransitionZonesList = findLongStretchesofConsistentTransitionZones(chrLong,
                                                                               start,
                                                                               end,
                                                                               chrBasedSmoothedWaveletReplicationTimeSignalDF,
                                                                               chrBased_valleys_peaks_df)

    labels = ['chr', 'start', 'end', 'slopeDirection', 'length']
    chrBasedTransitionZonesDF = pd.DataFrame.from_records(chrBasedTransitionZonesList, columns=labels)

    # Step2 Fill the replication array using transition zones
    chrBasedTransitionZonesDF.apply(fillReplicationStrandArray, chrBased_replication_array=chrBased_replication_array,axis=1)

    return chrBased_replication_array
########################################################################


########################################################################
def read_repliseq_dataframes(smoothedWaveletRepliseqDataFilename,valleysBEDFilename,peaksBEDFilename):

    ################### Read the Smoothed Wavelet Replication Time Signal starts ###########################
    #new way, JAN 7, 2020
    file_extension = os.path.splitext(os.path.basename(smoothedWaveletRepliseqDataFilename))[1]
    if (file_extension.lower() == '.wig'):
        isFileTypeBEDGRAPH = decideFileType(smoothedWaveletRepliseqDataFilename)
        if isFileTypeBEDGRAPH:
            repliseq_wavelet_signal_df=pd.read_csv(smoothedWaveletRepliseqDataFilename, sep='\t', comment='#', header=None, names=[CHROM,START,END,SIGNAL])
        else:
            repliseq_wavelet_signal_df=readWig_with_fixedStep_variableStep(smoothedWaveletRepliseqDataFilename)
    elif (file_extension.lower()=='.bedgraph'):
        repliseq_wavelet_signal_df = pd.read_csv(smoothedWaveletRepliseqDataFilename, sep='\t', comment='#',header=None, names=[CHROM, START, END, SIGNAL])

    print('Chromosome names in replication time signal data: %s' % (repliseq_wavelet_signal_df[CHROM].unique()))
    # print('repliseq_wavelet_signal_df[chr].unique')
    # print(repliseq_wavelet_signal_df['chr'].unique())
    ################### Read the Smoothed Wavelet Replication Time Signal ends #############################


    ############## Read the Valleys and Peaks starts #######################################
    discard_signal=True
    valleys_df= readFileInBEDFormat(valleysBEDFilename,discard_signal)
    valleys_df[END] = valleys_df[END] - 1
    print('Chromosome names in replication time valleys data: %s' % (valleys_df[CHROM].unique()))

    peaks_df = readFileInBEDFormat(peaksBEDFilename,discard_signal)
    peaks_df[END] = peaks_df[END] - 1
    print('Chromosome names in replication time peaks data: %s' % (peaks_df[CHROM].unique()))
    ############## Read the Valleys and Peaks ends ########################################

    return repliseq_wavelet_signal_df, valleys_df, peaks_df
########################################################################


########################################################################
# For df_split
# Using Numpy Arrays
# Search for all mutatitions
def searchAllMutationsOnReplicationStrandArray_for_df_split(chrBased_simBased_combined_df_split,
                                            chrBased_replication_array,
                                            sim_num,
                                            six_mutation_types_np_array,
                                            ordered_sbs_signatures,
                                            ordered_dbs_signatures,
                                            ordered_id_signatures,
                                            ordered_sbs_signatures_cutoffs,
                                            ordered_dbs_signatures_cutoffs,
                                            ordered_id_signatures_cutoffs,
                                            all_types_leading_np_array,
                                            all_types_lagging_np_array,
                                            subs_signature_mutation_type_leading_np_array,
                                            subs_signature_mutation_type_lagging_np_array,
                                            sample_based,
                                            is_discreet,
                                            verbose):

    ################################################################################
    if ((chrBased_simBased_combined_df_split is not None) and (not chrBased_simBased_combined_df_split.empty)):
        if verbose: print('\tVerbose Worker pid %s SBS searchMutationd_comOnReplicationStrandArray_simulations_integrated starts %s MB' % (str(os.getpid()), memory_usage()))

        #df_columns numpy array
        df_columns = chrBased_simBased_combined_df_split.columns.values

        ###############################################################################
        ################################ Initialization ###############################
        ###############################################################################
        df_columns_subs_signatures_mask_array = np.isin(df_columns,ordered_sbs_signatures)
        df_columns_dinucs_signatures_mask_array = np.isin(df_columns,ordered_dbs_signatures)
        df_columns_indels_signatures_mask_array = np.isin(df_columns,ordered_id_signatures)

        six_mutation_types_default_zeros_array= np.zeros(six_mutation_types_np_array.size) # dtype=int
        subs_signatures_default_zeros_array = np.zeros(ordered_sbs_signatures.size) # dtype=int
        dinucs_signatures_default_zeros_array = np.zeros(ordered_dbs_signatures.size) # dtype=int
        indels_signatures_default_zeros_array = np.zeros(ordered_id_signatures.size) # dtype=int
        subs_signatures_mutation_types_default_zeros_array = np.zeros((ordered_sbs_signatures.size, six_mutation_types_np_array.size)) # dtype=int
        ###############################################################################
        ################################ Initialization ###############################
        ###############################################################################

        ##############################################################################################
        #In list comprehesion, mutation_row becomes <class 'numpy.ndarray'>
        #In apply, mutation_row becomes <class 'pandas.core.series.Series'>
        #Therefore, list comprehesion is adopted.
        [searchAllMutationOnReplicationStrandArray_using_list_comprehension_using_numpy_array_for_df_split(mutation_row,
                                                                            chrBased_replication_array,
                                                                            six_mutation_types_np_array,
                                                                            ordered_sbs_signatures_cutoffs,
                                                                            ordered_dbs_signatures_cutoffs,
                                                                            ordered_id_signatures_cutoffs,
                                                                            df_columns_subs_signatures_mask_array,
                                                                            df_columns_dinucs_signatures_mask_array,
                                                                            df_columns_indels_signatures_mask_array,
                                                                            six_mutation_types_default_zeros_array,
                                                                            subs_signatures_default_zeros_array,
                                                                            dinucs_signatures_default_zeros_array,
                                                                            indels_signatures_default_zeros_array,
                                                                            subs_signatures_mutation_types_default_zeros_array,
                                                                            all_types_leading_np_array,
                                                                            all_types_lagging_np_array,
                                                                            subs_signature_mutation_type_leading_np_array,
                                                                            subs_signature_mutation_type_lagging_np_array,
                                                                            sample_based,
                                                                            is_discreet,
                                                                            df_columns) for mutation_row in chrBased_simBased_combined_df_split.values]
        ##############################################################################################

        if verbose: print('\tVerbose Worker pid %s SBS searchMutationOnReplicationStrandArray_simulations_integrated ends %s MB' % (str(os.getpid()), memory_usage()))
    ################################################################################

    return(sim_num,
           all_types_leading_np_array,
           all_types_lagging_np_array,
           subs_signature_mutation_type_leading_np_array,
           subs_signature_mutation_type_lagging_np_array)
########################################################################




########################################################################
# Using Numpy Arrays
# Search for all mutations
def searchAllMutationsOnReplicationStrandArray(chrBased_simBased_subs_df,
                                            chrBased_simBased_dinucs_df,
                                            chrBased_simBased_indels_df,
                                            chrBased_replication_array,
                                            sim_num,
                                            SBS6_mutation_types_np_array,
                                            SBS96_mutation_types_np_array,
                                            ordered_all_sbs_signatures_array,
                                            ordered_all_dbs_signatures_array,
                                            ordered_all_id_signatures_array,
                                            ordered_sbs_signatures,
                                            ordered_dbs_signatures,
                                            ordered_id_signatures,
                                            ordered_sbs_signatures_cutoffs,
                                            ordered_dbs_signatures_cutoffs,
                                            ordered_id_signatures_cutoffs,
                                            all_types_leading_np_array,
                                            all_types_lagging_np_array,
                                            subs_signature_SBS6_mutation_type_leading_np_array,
                                            subs_signature_SBS6_mutation_type_lagging_np_array,
                                            subs_signature_SBS96_mutation_type_leading_np_array,
                                            subs_signature_SBS96_mutation_type_lagging_np_array,
                                            all_samples_all_types_leading_np_array,
                                            all_samples_all_types_lagging_np_array,
                                            all_samples_subs_signature_mutation_type_leading_np_array,
                                            all_samples_subs_signature_mutation_type_lagging_np_array,
                                            sample_based,
                                            all_samples_np_array,
                                            is_discreet,
                                            verbose):

    if is_discreet:
        number_of_sbs_signatures = ordered_sbs_signatures.size
        number_of_dbs_signatures = ordered_dbs_signatures.size
        number_of_id_signatures = ordered_id_signatures.size
    else:
        number_of_sbs_signatures = ordered_all_sbs_signatures_array.size
        number_of_dbs_signatures = ordered_all_dbs_signatures_array.size
        number_of_id_signatures = ordered_all_id_signatures_array.size

    ###############################################################################
    ################################ Initialization ###############################
    ###############################################################################
    SBS6_mutation_types_default_zeros_array = np.zeros(SBS6_mutation_types_np_array.size) # dtype=int
    SBS96_mutation_types_default_zeros_array = np.zeros(SBS96_mutation_types_np_array.size) # dtype=int
    subs_signatures_default_zeros_array = np.zeros(number_of_sbs_signatures) # dtype=int
    dinucs_signatures_default_zeros_array = np.zeros(number_of_dbs_signatures) # dtype=int
    indels_signatures_default_zeros_array = np.zeros(number_of_id_signatures) # dtype=int
    subs_signatures_SBS6_mutation_types_default_zeros_array = np.zeros((number_of_sbs_signatures, SBS6_mutation_types_np_array.size)) # dtype=int
    subs_signatures_SBS96_mutation_types_default_zeros_array = np.zeros((number_of_sbs_signatures, SBS96_mutation_types_np_array.size)) # dtype=int
    ###############################################################################
    ################################ Initialization ###############################
    ###############################################################################

    ################################################################################
    #SUBS
    if ((chrBased_simBased_subs_df is not None) and (not chrBased_simBased_subs_df.empty)):
        if verbose: print('\tVerbose Worker pid %s SBS searchMutationd_comOnReplicationStrandArray_simulations_integrated starts %s MB' % (str(os.getpid()), memory_usage()))

        # df_columns numpy array
        df_columns = chrBased_simBased_subs_df.columns.values

        if is_discreet:
            df_columns_subs_signatures_mask_array = np.isin(df_columns, ordered_sbs_signatures)
        else:
            df_columns_subs_signatures_mask_array = np.isin(df_columns, ordered_all_sbs_signatures_array)

        ##############################################################################################
        # In list comprehesion, mutation_row becomes <class 'numpy.ndarray'>
        # In apply, mutation_row becomes <class 'pandas.core.series.Series'>
        # Therefore, list comprehesion is adopted.
        [searchAllMutationOnReplicationStrandArray_using_list_comprehension_using_numpy_array(mutation_row,
                                                                            SUBS,
                                                                            chrBased_replication_array,
                                                                            SBS6_mutation_types_np_array,
                                                                            SBS96_mutation_types_np_array,
                                                                            ordered_sbs_signatures_cutoffs,
                                                                            df_columns_subs_signatures_mask_array,
                                                                            SBS6_mutation_types_default_zeros_array,
                                                                            SBS96_mutation_types_default_zeros_array,
                                                                            subs_signatures_default_zeros_array,
                                                                            dinucs_signatures_default_zeros_array,
                                                                            indels_signatures_default_zeros_array,
                                                                            subs_signatures_SBS6_mutation_types_default_zeros_array,
                                                                            subs_signatures_SBS96_mutation_types_default_zeros_array,
                                                                            all_types_leading_np_array,
                                                                            all_types_lagging_np_array,
                                                                            subs_signature_SBS6_mutation_type_leading_np_array,
                                                                            subs_signature_SBS6_mutation_type_lagging_np_array,
                                                                            subs_signature_SBS96_mutation_type_leading_np_array,
                                                                            subs_signature_SBS96_mutation_type_lagging_np_array,
                                                                            all_samples_all_types_leading_np_array,
                                                                            all_samples_all_types_lagging_np_array,
                                                                            all_samples_subs_signature_mutation_type_leading_np_array,
                                                                            all_samples_subs_signature_mutation_type_lagging_np_array,
                                                                            sample_based,
                                                                            all_samples_np_array,
                                                                            is_discreet,
                                                                            df_columns) for mutation_row in chrBased_simBased_subs_df.values]
        ##############################################################################################

    ################################################################################
    # DINUCS
    if ((chrBased_simBased_dinucs_df is not None) and (not chrBased_simBased_dinucs_df.empty)):
        if verbose: print('\tVerbose Worker pid %s SBS searchMutationd_comOnReplicationStrandArray_simulations_integrated starts %s MB' % (str(os.getpid()), memory_usage()))

        # df_columns numpy array
        df_columns = chrBased_simBased_dinucs_df.columns.values

        if is_discreet:
            df_columns_dinucs_signatures_mask_array = np.isin(df_columns, ordered_dbs_signatures)
        else:
            df_columns_dinucs_signatures_mask_array = np.isin(df_columns, ordered_all_dbs_signatures_array)

        ##############################################################################################
        # In list comprehesion, mutation_row becomes <class 'numpy.ndarray'>
        # In apply, mutation_row becomes <class 'pandas.core.series.Series'>
        # Therefore, list comprehesion is adopted.
        [searchAllMutationOnReplicationStrandArray_using_list_comprehension_using_numpy_array(mutation_row,
                                                                                              DINUCS,
                                                                                              chrBased_replication_array,
                                                                                              SBS6_mutation_types_np_array,
                                                                                              SBS96_mutation_types_np_array,
                                                                                              ordered_dbs_signatures_cutoffs,
                                                                                              df_columns_dinucs_signatures_mask_array,
                                                                                              SBS6_mutation_types_default_zeros_array,
                                                                                              SBS96_mutation_types_default_zeros_array,
                                                                                              subs_signatures_default_zeros_array,
                                                                                              dinucs_signatures_default_zeros_array,
                                                                                              indels_signatures_default_zeros_array,
                                                                                              subs_signatures_SBS6_mutation_types_default_zeros_array,
                                                                                              subs_signatures_SBS96_mutation_types_default_zeros_array,
                                                                                              all_types_leading_np_array,
                                                                                              all_types_lagging_np_array,
                                                                                              subs_signature_SBS6_mutation_type_leading_np_array,
                                                                                              subs_signature_SBS6_mutation_type_lagging_np_array,
                                                                                              subs_signature_SBS96_mutation_type_leading_np_array,
                                                                                              subs_signature_SBS96_mutation_type_lagging_np_array,
                                                                                              all_samples_all_types_leading_np_array,
                                                                                              all_samples_all_types_lagging_np_array,
                                                                                              all_samples_subs_signature_mutation_type_leading_np_array,
                                                                                              all_samples_subs_signature_mutation_type_lagging_np_array,
                                                                                              sample_based,
                                                                                              all_samples_np_array,
                                                                                              is_discreet,
                                                                                              df_columns) for mutation_row in chrBased_simBased_dinucs_df.values]
        ##############################################################################################

    ################################################################################
    # INDELS
    if ((chrBased_simBased_indels_df is not None) and (not chrBased_simBased_indels_df.empty)):
        if verbose: print('\tVerbose Worker pid %s SBS searchMutationd_comOnReplicationStrandArray_simulations_integrated starts %s MB' % (str(os.getpid()), memory_usage()))

        # df_columns numpy array
        df_columns = chrBased_simBased_indels_df.columns.values

        if is_discreet:
            df_columns_indels_signatures_mask_array = np.isin(df_columns, ordered_id_signatures)
        else:
            df_columns_indels_signatures_mask_array = np.isin(df_columns, ordered_all_id_signatures_array)

        ##############################################################################################
        # In list comprehesion, mutation_row becomes <class 'numpy.ndarray'>
        # In apply, mutation_row becomes <class 'pandas.core.series.Series'>
        # Therefore, list comprehesion is adopted.
        [searchAllMutationOnReplicationStrandArray_using_list_comprehension_using_numpy_array(mutation_row,
                                                                                              INDELS,
                                                                                              chrBased_replication_array,
                                                                                              SBS6_mutation_types_np_array,
                                                                                              SBS96_mutation_types_np_array,
                                                                                              ordered_id_signatures_cutoffs,
                                                                                              df_columns_indels_signatures_mask_array,
                                                                                              SBS6_mutation_types_default_zeros_array,
                                                                                              SBS96_mutation_types_default_zeros_array,
                                                                                              subs_signatures_default_zeros_array,
                                                                                              dinucs_signatures_default_zeros_array,
                                                                                              indels_signatures_default_zeros_array,
                                                                                              subs_signatures_SBS6_mutation_types_default_zeros_array,
                                                                                              subs_signatures_SBS96_mutation_types_default_zeros_array,
                                                                                              all_types_leading_np_array,
                                                                                              all_types_lagging_np_array,
                                                                                              subs_signature_SBS6_mutation_type_leading_np_array,
                                                                                              subs_signature_SBS6_mutation_type_lagging_np_array,
                                                                                              subs_signature_SBS96_mutation_type_leading_np_array,
                                                                                              subs_signature_SBS96_mutation_type_lagging_np_array,
                                                                                              all_samples_all_types_leading_np_array,
                                                                                              all_samples_all_types_lagging_np_array,
                                                                                              all_samples_subs_signature_mutation_type_leading_np_array,
                                                                                              all_samples_subs_signature_mutation_type_lagging_np_array,
                                                                                              sample_based,
                                                                                              all_samples_np_array,
                                                                                              is_discreet,
                                                                                              df_columns) for mutation_row in chrBased_simBased_indels_df.values]
        ##############################################################################################

        if verbose: print('\tVerbose Worker pid %s SBS searchMutationOnReplicationStrandArray_simulations_integrated ends %s MB' % (str(os.getpid()), memory_usage()))
    ################################################################################

    return(sim_num,
           all_types_leading_np_array,
           all_types_lagging_np_array,
           subs_signature_SBS6_mutation_type_leading_np_array,
           subs_signature_SBS6_mutation_type_lagging_np_array,
           subs_signature_SBS96_mutation_type_leading_np_array,
           subs_signature_SBS96_mutation_type_lagging_np_array,
           all_samples_all_types_leading_np_array,
           all_samples_all_types_lagging_np_array,
           all_samples_subs_signature_mutation_type_leading_np_array,
           all_samples_subs_signature_mutation_type_lagging_np_array)
########################################################################


########################################################################
def searchAllMutationsOnReplicationStrandArray_simbased_chrombased_splitbased(outputDir,
                                                                              jobname,
                                                                              chrLong,
                                                                              simNum,
                                                                              splitIndex,
                                                                              six_mutation_types_np_array,
                                                                              ordered_sbs_signatures,
                                                                              ordered_dbs_signatures,
                                                                              ordered_id_signatures,
                                                                              ordered_sbs_signatures_cutoffs,
                                                                              ordered_dbs_signatures_cutoffs,
                                                                              ordered_id_signatures_cutoffs,
                                                                              all_types_np_array_size,
                                                                              sample_based,
                                                                              is_discreet,
                                                                              verbose):

    chr_based_replication_time_file_name = '%s_replication_time.npy' % (chrLong)
    chr_based_replication_time_file_path = os.path.join(outputDir, jobname, DATA, REPLICATIONSTRANDBIAS, LIB, CHRBASED,chr_based_replication_time_file_name)

    # Initialization
    all_types_leading_np_array = np.zeros((all_types_np_array_size)) # dtype=int
    all_types_lagging_np_array = np.zeros((all_types_np_array_size)) # dtype=int
    subs_signature_mutation_type_leading_np_array = np.zeros((ordered_sbs_signatures.size, six_mutation_types_np_array.size)) # dtype=int
    subs_signature_mutation_type_lagging_np_array = np.zeros((ordered_sbs_signatures.size, six_mutation_types_np_array.size)) # dtype=int

    if (os.path.exists(chr_based_replication_time_file_path)):
        chrBased_replication_array = np.load(chr_based_replication_time_file_path)
    else:
        chrBased_replication_array=None

    if chrBased_replication_array is not None:
        chrBased_simBased_combined_df_split = get_chrBased_simBased_combined_df_split(outputDir, jobname, chrLong,simNum, splitIndex)

        return searchAllMutationsOnReplicationStrandArray_for_df_split(chrBased_simBased_combined_df_split,
                                                          chrBased_replication_array,
                                                          simNum,
                                                          six_mutation_types_np_array,
                                                          ordered_sbs_signatures,
                                                          ordered_dbs_signatures,
                                                          ordered_id_signatures,
                                                          ordered_sbs_signatures_cutoffs,
                                                          ordered_dbs_signatures_cutoffs,
                                                          ordered_id_signatures_cutoffs,
                                                          all_types_leading_np_array,
                                                          all_types_lagging_np_array,
                                                          subs_signature_mutation_type_leading_np_array,
                                                          subs_signature_mutation_type_lagging_np_array,
                                                          sample_based,
                                                          is_discreet,
                                                          verbose)
    else:
        return (simNum,
                all_types_leading_np_array,
                all_types_lagging_np_array,
                subs_signature_mutation_type_leading_np_array,
                subs_signature_mutation_type_lagging_np_array)
########################################################################


########################################################################
# April 30, 2020
# Read chromBased and simBased combined (SBS, DBS and ID) dataframe in the process
def searchAllMutationsOnReplicationStrandArray_simbased_chrombased(outputDir,
                                                                   jobname,
                                                                   chrLong,
                                                                   simNum,
                                                                   SBS6_mutation_types_np_array,
                                                                   SBS96_mutation_types_np_array,
                                                                   ordered_all_sbs_signatures_array,
                                                                   ordered_all_dbs_signatures_array,
                                                                   ordered_all_id_signatures_array,
                                                                   ordered_sbs_signatures,
                                                                   ordered_dbs_signatures,
                                                                   ordered_id_signatures,
                                                                   ordered_sbs_signatures_cutoffs,
                                                                   ordered_dbs_signatures_cutoffs,
                                                                   ordered_id_signatures_cutoffs,
                                                                   all_types_np_array_size,
                                                                   sample_based,
                                                                   all_samples_np_array,
                                                                   is_discreet,
                                                                   verbose):

    chr_based_replication_time_file_name = '%s_replication_time.npy' % (chrLong)
    chr_based_replication_time_file_path = os.path.join(outputDir, jobname, DATA, REPLICATIONSTRANDBIAS, LIB, CHRBASED,chr_based_replication_time_file_name)

    if is_discreet:
        number_of_sbs_signatures = ordered_sbs_signatures.size
    else:
        number_of_sbs_signatures = ordered_all_sbs_signatures_array.size

    # Initialization
    all_types_leading_np_array = np.zeros((all_types_np_array_size)) # dtype=int
    all_types_lagging_np_array = np.zeros((all_types_np_array_size)) # dtype=int
    subs_signature_SBS6_mutation_type_leading_np_array = np.zeros((number_of_sbs_signatures, SBS6_mutation_types_np_array.size)) # dtype=int
    subs_signature_SBS6_mutation_type_lagging_np_array = np.zeros((number_of_sbs_signatures, SBS6_mutation_types_np_array.size)) # dtype=int
    subs_signature_SBS96_mutation_type_leading_np_array = np.zeros((number_of_sbs_signatures, SBS96_mutation_types_np_array.size)) # dtype=int
    subs_signature_SBS96_mutation_type_lagging_np_array = np.zeros((number_of_sbs_signatures, SBS96_mutation_types_np_array.size)) # dtype=int

    all_samples_np_array_size = all_samples_np_array.size
    all_samples_all_types_leading_np_array = np.zeros((all_samples_np_array_size, all_types_np_array_size)) # dtype=int
    all_samples_all_types_lagging_np_array = np.zeros((all_samples_np_array_size, all_types_np_array_size)) # dtype=int
    all_samples_subs_signature_mutation_type_leading_np_array = np.zeros((all_samples_np_array_size, number_of_sbs_signatures, SBS6_mutation_types_np_array.size)) # dtype=int
    all_samples_subs_signature_mutation_type_lagging_np_array = np.zeros((all_samples_np_array_size, number_of_sbs_signatures, SBS6_mutation_types_np_array.size)) # dtype=int

    if (os.path.exists(chr_based_replication_time_file_path)):
        chrBased_replication_array = np.load(chr_based_replication_time_file_path)
    else:
        chrBased_replication_array=None

    if chrBased_replication_array is not None:
        chrBased_simBased_subs_df, chrBased_simBased_dinucs_df, chrBased_simBased_indels_df = get_chrBased_simBased_dfs(outputDir, jobname, chrLong, simNum)

        return searchAllMutationsOnReplicationStrandArray(chrBased_simBased_subs_df,
                                                          chrBased_simBased_dinucs_df,
                                                          chrBased_simBased_indels_df,
                                                          chrBased_replication_array,
                                                          simNum,
                                                          SBS6_mutation_types_np_array,
                                                          SBS96_mutation_types_np_array,
                                                          ordered_all_sbs_signatures_array,
                                                          ordered_all_dbs_signatures_array,
                                                          ordered_all_id_signatures_array,
                                                          ordered_sbs_signatures,
                                                          ordered_dbs_signatures,
                                                          ordered_id_signatures,
                                                          ordered_sbs_signatures_cutoffs,
                                                          ordered_dbs_signatures_cutoffs,
                                                          ordered_id_signatures_cutoffs,
                                                          all_types_leading_np_array,
                                                          all_types_lagging_np_array,
                                                          subs_signature_SBS6_mutation_type_leading_np_array,
                                                          subs_signature_SBS6_mutation_type_lagging_np_array,
                                                          subs_signature_SBS96_mutation_type_leading_np_array,
                                                          subs_signature_SBS96_mutation_type_lagging_np_array,
                                                          all_samples_all_types_leading_np_array,
                                                          all_samples_all_types_lagging_np_array,
                                                          all_samples_subs_signature_mutation_type_leading_np_array,
                                                          all_samples_subs_signature_mutation_type_lagging_np_array,
                                                          sample_based,
                                                          all_samples_np_array,
                                                          is_discreet,
                                                          verbose)
    else:
        return (simNum,
                all_types_leading_np_array,
                all_types_lagging_np_array,
                subs_signature_SBS6_mutation_type_leading_np_array,
                subs_signature_SBS6_mutation_type_lagging_np_array,
                subs_signature_SBS96_mutation_type_leading_np_array,
                subs_signature_SBS96_mutation_type_lagging_np_array,
                all_samples_all_types_leading_np_array,
                all_samples_all_types_lagging_np_array,
                all_samples_subs_signature_mutation_type_leading_np_array,
                all_samples_subs_signature_mutation_type_lagging_np_array)
########################################################################


########################################################################
def read_create_write_replication_time_array_in_parallel(outputDir,jobname,chromNamesList,chromSizesDict,smoothedWaveletRepliseqDataFilename,valleysBEDFilename,peaksBEDFilename,verbose):

    repliseq_signal_df, valleys_df, peaks_df = read_repliseq_dataframes(smoothedWaveletRepliseqDataFilename,valleysBEDFilename,peaksBEDFilename)

    ################################
    def write_chrom_based_replication_array(result_tuple):
        chrLong=result_tuple[0]
        chrBased_replication_array=result_tuple[1]
        if (chrBased_replication_array is not None):
            os.makedirs(os.path.join(outputDir, jobname, DATA, REPLICATIONSTRANDBIAS, LIB, CHRBASED), exist_ok=True)
            #File name without extension
            chr_based_replication_time_file_name='%s_replication_time' %(chrLong)
            chr_based_replication_time_file_path = os.path.join(outputDir,jobname,DATA,REPLICATIONSTRANDBIAS,LIB,CHRBASED,chr_based_replication_time_file_name)
            np.save(chr_based_replication_time_file_path, chrBased_replication_array)
    ################################

    ################################
    numofProcesses = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=numofProcesses)
    ################################

    ################################
    jobs = []
    ################################

    for chrLong in chromNamesList:
        chromSize = chromSizesDict[chrLong]
        jobs.append(pool.apply_async(get_chr_based_replication_strand_array_for_callback,
                                 args=(chrLong, chromSize, repliseq_signal_df,valleys_df, peaks_df,),
                                 callback=write_chrom_based_replication_array))

    ##############################################################################
    # wait for all jobs to finish
    for job in jobs:
        if verbose: print('\tVerbose Write Chrom Based Replication Time Array for Replicatio Strand Bias Analysis Worker pid %s job.get():%s ' % (str(os.getpid()), job.get()))
    ##############################################################################

    ################################
    pool.close()
    pool.join()
    ################################

########################################################################


########################################################################
# pool.imap_unordered : Fills poolInputList (jobs), sends poolInputList and accumulates results one by one.
# Slowest one but completes without memory problem.
# Can be updated and tested to read chrom based sim based mutations data and chrom based replication array in the worker process.
#
# pool.apply_async:  USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT
#
# pool.apply_async:  USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT_USING_POOL_INPUT_LIST
# You fill your pool input list therefore you decide how many jobs  to send at once.
# Faster than imap_unordered. Low memory usage.
# Can be updated to read chrom based sim based mutations data and chrom based replication array in the worker process.
# If USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM does not end because of memory error, this can be used.
# All 28/28 processes are running. When the jobs in the pool input list are finishing some processes waits others to finish.
#
# pool.apply_async:  USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM
# For each possible (chrLong,simNum) couple read the data and array on the worker process
# Fastest, consumes more memory than others. 22/28 processes are running. For Combined_PACWG_nonPCAWG Skin_Melanoma after 1 hour all 28/28 running.
def replicationStrandBiasAnalysis(outputDir,
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
                                  ordered_all_sbs_signatures_array,
                                  ordered_all_dbs_signatures_array,
                                  ordered_all_id_signatures_array,
                                  ordered_sbs_signatures,
                                  ordered_dbs_signatures,
                                  ordered_id_signatures,
                                  ordered_sbs_signatures_cutoffs,
                                  ordered_dbs_signatures_cutoffs,
                                  ordered_id_signatures_cutoffs,
                                  is_discreet,
                                  verbose):

    print('\n#################################################################################')
    print('--- ReplicationStrandBias Analysis starts')

    ###############################################
    #April 30, 2020
    read_create_write_replication_time_array_in_parallel(outputDir,jobname,chromNamesList,chromSizesDict,smoothedWaveletRepliseqDataFilename,valleysBEDFilename,peaksBEDFilename,verbose)
    ###############################################

    #########################################################################################
    SBS6_mutation_types_np_array = np.array([C2A, C2G, C2T, T2A, T2C, T2G])

    nucleotides = ['A', 'C', 'G', 'T']
    mutations = ['[C>A]', '[C>G]', '[C>T]', '[T>A]', '[T>C]', '[T>G]']
    SBS96_mutation_types_np_array = np.array([nucleotide_left + middle + nucleotide_right
                                                for nucleotide_left in nucleotides
                                                    for nucleotide_right in nucleotides
                                                        for middle in mutations])

    if is_discreet:
        all_types_np_array = np.concatenate((SBS6_mutation_types_np_array,
                                         SBS96_mutation_types_np_array,
                                         ordered_sbs_signatures,
                                         ordered_dbs_signatures,
                                         ordered_id_signatures), axis=None)

        number_of_sbs_signatures = ordered_sbs_signatures.size
        sbs_signatures = ordered_sbs_signatures

    else:
        all_types_np_array = np.concatenate((SBS6_mutation_types_np_array,
                                         SBS96_mutation_types_np_array,
                                         ordered_all_sbs_signatures_array,
                                         ordered_all_dbs_signatures_array,
                                         ordered_all_id_signatures_array), axis=None)
        sbs_signatures = ordered_all_sbs_signatures_array


        number_of_sbs_signatures = ordered_all_sbs_signatures_array.size


    all_types_np_array_size = all_types_np_array.size
    #########################################################################################

    #########################################################################################
    # Initialization
    all_sims_all_types_leading_np_array = np.zeros((numofSimulations+1, all_types_np_array_size)) # dtype=int
    all_sims_all_types_lagging_np_array = np.zeros((numofSimulations+1, all_types_np_array_size)) # dtype=int
    all_sims_subs_signature_SBS6_mutation_type_leading_np_array = np.zeros((numofSimulations+1, number_of_sbs_signatures, SBS6_mutation_types_np_array.size)) # dtype=int
    all_sims_subs_signature_SBS6_mutation_type_lagging_np_array= np.zeros((numofSimulations+1, number_of_sbs_signatures, SBS6_mutation_types_np_array.size)) # dtype=int
    all_sims_subs_signature_SBS96_mutation_type_leading_np_array = np.zeros((numofSimulations+1, number_of_sbs_signatures, SBS96_mutation_types_np_array.size)) # dtype=int
    all_sims_subs_signature_SBS96_mutation_type_lagging_np_array= np.zeros((numofSimulations+1, number_of_sbs_signatures, SBS96_mutation_types_np_array.size)) # dtype=int
    #########################################################################################

    #########################################################################################
    if sample_based:
        # Initialization for accumulated arrays
        all_samples_np_array_size = all_samples_np_array.size
        all_sims_all_samples_all_types_leading_np_array = np.zeros((numofSimulations + 1, all_samples_np_array_size, all_types_np_array_size)) # dtype=int
        all_sims_all_samples_all_types_lagging_np_array = np.zeros((numofSimulations + 1, all_samples_np_array_size, all_types_np_array_size)) # dtype=int
        all_sims_all_samples_subs_signature_mutation_type_leading_np_array = np.zeros((numofSimulations + 1, all_samples_np_array_size, number_of_sbs_signatures, SBS6_mutation_types_np_array.size)) # dtype=int
        all_sims_all_samples_subs_signature_mutation_type_lagging_np_array = np.zeros((numofSimulations + 1, all_samples_np_array_size, number_of_sbs_signatures, SBS6_mutation_types_np_array.size)) # dtype=int
    #########################################################################################

    #########################################################################################
    # Accumulate Numpy Arrays
    # July 27, 2020
    def accumulate_np_arrays(result_tuple):
        sim_num = result_tuple[0]
        all_types_leading_np_array = result_tuple[1]
        all_types_lagging_np_array = result_tuple[2]
        subs_signature_SBS6_mutation_type_leading_np_array = result_tuple[3]
        subs_signature_SBS6_mutation_type_lagging_np_array = result_tuple[4]
        subs_signature_SBS96_mutation_type_leading_np_array = result_tuple[5]
        subs_signature_SBS96_mutation_type_lagging_np_array = result_tuple[6]

        all_sims_all_types_leading_np_array[sim_num] += all_types_leading_np_array
        all_sims_all_types_lagging_np_array[sim_num] += all_types_lagging_np_array
        all_sims_subs_signature_SBS6_mutation_type_leading_np_array[sim_num] += subs_signature_SBS6_mutation_type_leading_np_array
        all_sims_subs_signature_SBS6_mutation_type_lagging_np_array[sim_num] += subs_signature_SBS6_mutation_type_lagging_np_array
        all_sims_subs_signature_SBS96_mutation_type_leading_np_array[sim_num] += subs_signature_SBS96_mutation_type_leading_np_array
        all_sims_subs_signature_SBS96_mutation_type_lagging_np_array[sim_num] += subs_signature_SBS96_mutation_type_lagging_np_array
        # print('MONITOR ACCUMULATE', flush=True)

        if sample_based:
            #Jan 21, 2021
            all_samples_all_types_leading_np_array = result_tuple[5]
            all_samples_all_types_lagging_np_array = result_tuple[6]
            all_samples_subs_signature_mutation_type_leading_np_array = result_tuple[7]
            all_samples_subs_signature_mutation_type_lagging_np_array = result_tuple[8]

            all_sims_all_samples_all_types_leading_np_array[sim_num] += all_samples_all_types_leading_np_array
            all_sims_all_samples_all_types_lagging_np_array[sim_num] += all_samples_all_types_lagging_np_array
            all_sims_all_samples_subs_signature_mutation_type_leading_np_array[sim_num] += all_samples_subs_signature_mutation_type_leading_np_array
            all_sims_all_samples_subs_signature_mutation_type_lagging_np_array[sim_num] += all_samples_subs_signature_mutation_type_lagging_np_array
    #########################################################################################

    jobs = []

    ###############################################################################
    #April 30, 2020
    #read chrom based sim based mutations data and chrom based replication time data in each worker process
    if (computation_type==USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM):
        sim_nums = range(0, numofSimulations + 1)
        sim_num_chr_tuples = ((sim_num, chrLong) for sim_num in sim_nums for chrLong in chromNamesList)

        numofProcesses = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=numofProcesses)

        for simNum, chrLong in sim_num_chr_tuples:
            jobs.append(pool.apply_async(searchAllMutationsOnReplicationStrandArray_simbased_chrombased,
                                    args=(outputDir,
                                          jobname,
                                          chrLong,
                                          simNum,
                                          SBS6_mutation_types_np_array,
                                          SBS96_mutation_types_np_array,
                                          ordered_all_sbs_signatures_array,
                                          ordered_all_dbs_signatures_array,
                                          ordered_all_id_signatures_array,
                                          ordered_sbs_signatures,
                                          ordered_dbs_signatures,
                                          ordered_id_signatures,
                                          ordered_sbs_signatures_cutoffs,
                                          ordered_dbs_signatures_cutoffs,
                                          ordered_id_signatures_cutoffs,
                                          all_types_np_array_size,
                                          sample_based,
                                          all_samples_np_array,
                                          is_discreet,
                                          verbose,),
                                    callback=accumulate_np_arrays))
        pool.close()
        pool.join()
    ###############################################################################

    ###############################################################################
    #July 23, 2020 starts
    elif (computation_type==USING_APPLY_ASYNC_FOR_EACH_CHROM_AND_SIM_SPLIT):
        numofProcesses = multiprocessing.cpu_count()
        pool = multiprocessing.Pool(processes=numofProcesses)

        for chrLong, simNum, splitIndex in job_tuples:
            jobs.append(pool.apply_async(searchAllMutationsOnReplicationStrandArray_simbased_chrombased_splitbased,
                                    args=(outputDir,
                                          jobname,
                                          chrLong,
                                          simNum,
                                          splitIndex,
                                          SBS6_mutation_types_np_array,
                                          ordered_sbs_signatures,
                                          ordered_dbs_signatures,
                                          ordered_id_signatures,
                                          ordered_sbs_signatures_cutoffs,
                                          ordered_dbs_signatures_cutoffs,
                                          ordered_id_signatures_cutoffs,
                                          all_types_np_array_size,
                                          sample_based,
                                          is_discreet,
                                          verbose,),
                                    callback=accumulate_np_arrays))

        pool.close()
        pool.join()
    #July 23, 2020 ends
    ###############################################################################


    ############################################################################################################
    #####################################       Output starts      #############################################
    ############################################################################################################

    strand_bias = REPLICATIONSTRANDBIAS
    replication_strands = [LAGGING,LEADING]

    SBS6_np_arrays_list =[  all_sims_subs_signature_SBS6_mutation_type_lagging_np_array,
                            all_sims_subs_signature_SBS6_mutation_type_leading_np_array]

    # Write files for real data
    # For each signature all_sims_subs_signature_SBS96_mutation_type_lagging_np_array
    # For each signature all_sims_subs_signature_SBS96_mutation_type_leading_np_array
    write_sbs_signature_sbs96_mutation_type_replication_strand_bias(all_sims_subs_signature_SBS96_mutation_type_lagging_np_array[0],
                                                                    all_sims_subs_signature_SBS96_mutation_type_leading_np_array[1],
                                                                    SBS96_mutation_types_np_array,
                                                                    sbs_signatures,
                                                                    strand_bias,
                                                                    outputDir,
                                                                    jobname)

    write_signature_mutation_type_strand_bias_np_array_as_dataframe(SBS6_np_arrays_list,
                                                                    SBS6_mutation_types_np_array,
                                                                    sbs_signatures,
                                                                    strand_bias,
                                                                    replication_strands,
                                                                    outputDir,
                                                                    jobname)

    all_sims_all_types_strand_np_arrays_list =[all_sims_all_types_lagging_np_array,
                                               all_sims_all_types_leading_np_array]

    write_type_strand_bias_np_array_as_dataframe(all_sims_all_types_strand_np_arrays_list,
                                                all_types_np_array,
                                                strand_bias,
                                                replication_strands,
                                                outputDir,
                                                jobname)


    if sample_based:
        write_sample_based_strand1_strand2_as_dataframe(outputDir,
                                                        jobname,
                                                        numofSimulations,
                                                        strand_bias,
                                                        all_samples_np_array,
                                                        all_types_np_array,
                                                        all_sims_all_samples_all_types_lagging_np_array,
                                                        all_sims_all_samples_all_types_leading_np_array)

    ############################################################################################################
    #####################################       Output ends      ###############################################
    ############################################################################################################

    print('--- ReplicationStrandBias Analysis ends')
    print('#################################################################################\n')
########################################################################



