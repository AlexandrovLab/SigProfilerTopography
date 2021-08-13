# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018 Burcak Otlu

#Read the processivityDict for original data
#Read the processivityDict for simulations
#Plot the signatures versus processivity group lengths figure where the circle radius shows the number of processive groups
# and  the color represents the signficance of number of processive groups in original data w.r.t. simulations data

import os
import sys
import shutil
import statsmodels.stats.multitest
import math
import  numpy as np
import pandas as pd

import matplotlib
BACKEND = 'Agg'
if matplotlib.get_backend().lower() != BACKEND.lower():
    # If backend is not set properly a call to describe will hang
    matplotlib.use(BACKEND)

from matplotlib import pyplot as plt

import matplotlib as mpl
import matplotlib.cm as cm
from statsmodels.stats.weightstats import ztest
from matplotlib.colors import Normalize

from SigProfilerTopography.source.commons.TopographyCommons import DATA
from SigProfilerTopography.source.commons.TopographyCommons import FIGURE
from SigProfilerTopography.source.commons.TopographyCommons import PROCESSIVITY
from SigProfilerTopography.source.commons.TopographyCommons import TABLES
from SigProfilerTopography.source.commons.TopographyCommons import USING_ONE_SAMPLE_TTEST
from SigProfilerTopography.source.commons.TopographyCommons import USING_NULL_DISTRIBUTION
from SigProfilerTopography.source.commons.TopographyCommons import USING_GAUSSIAN_KDE
from SigProfilerTopography.source.commons.TopographyCommons import USING_ZSCORE
from SigProfilerTopography.source.commons.TopographyCommons import FDR_BH_CORRECTION
from SigProfilerTopography.source.commons.TopographyCommons import BONFERRONI_CORRECTION

from SigProfilerTopography.source.commons.TopographyCommons import readDictionary
from SigProfilerTopography.source.commons.TopographyCommons import natural_key

plt.rcParams.update({'figure.max_open_warning': 0})

# SigProfilerTopography PROCESSIVITY CONSTRAINTS
PROCESSIVITY_SIGNIFICANCE_LEVEL = 0.01
MINIMUM_REQUIRED_PROCESSIVE_GROUP_LENGTH = 2
MINIMUM_REQUIRED_NUMBER_OF_PROCESSIVE_GROUPS = 2

###################################################################
def readSimulationBasedDictionaries(outputDir,jobname,numberofSimulations):
    simulation2Signature2ProcessiveGroupLength2PropertiesDict = {}

    for simNum in range(1,numberofSimulations+1):
        filename = 'Sim%d_Signature2ProcessiveGroupLength2PropertiesDict.txt' % (simNum)
        simulationFilePath = os.path.join(outputDir,jobname,DATA,PROCESSIVITY,filename)

        signature2ProcessiveGroupLength2PropertiesDict =  readDictionary(simulationFilePath)

        if (signature2ProcessiveGroupLength2PropertiesDict is not None):
            simulation2Signature2ProcessiveGroupLength2PropertiesDict[simNum] = signature2ProcessiveGroupLength2PropertiesDict

    return simulation2Signature2ProcessiveGroupLength2PropertiesDict
###################################################################

###################################################################
def plot_color_bar(outputDir,jobname,norm):

    # Make a figure and axes with dimensions as desired.
    fig = plt.figure(figsize=(4, 8))
    ax = fig.add_axes([0.05, 0.05, 0.15, 0.9])

    cmap = cm.get_cmap('YlOrRd')  # Looks better good
    # cb=fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)

    cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                     norm=norm,
                                     spacing='proportional',
                                     orientation='vertical')

    # cb = plt.colorbar(cmap=cmap,ax=ax,orientation='vertical')  # this works because of the scatter
    cb.ax.set_xticklabels(cb.ax.get_xticklabels(), fontsize=20)
    cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=20)

    font = mpl.font_manager.FontProperties(size=30)
    cbax = cb.ax
    text_x = cbax.xaxis.label
    text_y = cbax.yaxis.label
    text_x.set_font_properties(font)
    text_y.set_font_properties(font)

    cb.set_label("-log10\n  (q-value)", horizontalalignment='right', rotation=0, labelpad=150)

    ##################################################################################
    #create the directory if it does not exists
    os.makedirs(os.path.join(outputDir, jobname, FIGURE, PROCESSIVITY), exist_ok=True)
    filename = '%s_Processivity_ColorBar.png' %(jobname)

    figFile = os.path.join(outputDir, jobname, FIGURE, PROCESSIVITY, filename)
    fig.savefig(figFile)
    plt.cla()
    plt.close(fig)
    ##################################################################################

###################################################################

###################################################################
def set_radius(df):
    df['radius'] = df['log10_number_of_processive_groups'] / df['log10_number_of_processive_groups'].max() * 0.48
    return df
###################################################################

###################################################################
def plotRelationshipBetweenSignaturesandProcessiveGroupLengthsUsingDataframes(outputDir,jobname,processivity_df,numberofSimulations,verbose):

    # processivity_df: columns below
    # Simulation_Number
    # Signature
    # Processsive_Group_Length
    # Number_of_Processive_Groups
    # Median_Distance_Between_Last_First_Mutations
    # Median_Distance_Between_Consecutive_Mutations
    # Median_Number_of_Mutations_Within_1MB

    ####################################################
    # Get the list of signatures using original
    signatures_list = processivity_df[processivity_df['Simulation_Number']==0]['Signature'].unique()
    sorted_signature_list = sorted(signatures_list,reverse=True,key=natural_key)
    ####################################################

    ####################################################
    # Get the list of processive group lengths using original
    processsive_group_length_list = processivity_df[processivity_df['Simulation_Number']==0]['Processsive_Group_Length'].unique()
    sorted_processsive_group_length_list = sorted(processsive_group_length_list, key=int)
    ####################################################

    signature_processive_group_length_properties_df = pd.DataFrame(columns=["signature",
                                                                            "processive_group_length",
                                                                            "number_of_processive_groups",
                                                                            "log10_number_of_processive_groups",
                                                                            "radius",
                                                                            "avg_sims",
                                                                            "min_sims",
                                                                            "max_sims",
                                                                            "mean_sims",
                                                                            "std_sims",
                                                                            "pvalue",
                                                                            "qvalue",
                                                                            "minus_log10_qvalue", #will be used in coloring
                                                                            "zscore", # for information only
                                                                            "expected_number_of_processive_groups"])
    for signature in sorted_signature_list:
        for processive_group_length in sorted_processsive_group_length_list:

            if processivity_df[
                (processivity_df['Simulation_Number'] == 0) &
                (processivity_df['Signature'] == signature) &
                (processivity_df['Processsive_Group_Length'] == processive_group_length)].values.any():

                expected_number_of_processive_groups = []

                number_of_processive_groups = processivity_df[(processivity_df['Simulation_Number'] == 0) &
                                                              (processivity_df['Signature'] == signature) &
                                                              (processivity_df['Processsive_Group_Length'] == processive_group_length)]['Number_of_Processive_Groups'].values[0]


                if processivity_df[
                    (processivity_df['Simulation_Number'] != 0) &
                    (processivity_df['Signature'] == signature) &
                    (processivity_df['Processsive_Group_Length'] == processive_group_length)].values.any():

                    expected_number_of_processive_groups = processivity_df[(processivity_df['Simulation_Number'] != 0) &
                                                                  (processivity_df['Signature'] == signature) &
                                                                  (processivity_df['Processsive_Group_Length'] == processive_group_length)]['Number_of_Processive_Groups'].values.tolist()

                signature_processive_group_length_properties_df = signature_processive_group_length_properties_df.append(
                    {"signature": signature,
                     "processive_group_length": processive_group_length,
                     "number_of_processive_groups": number_of_processive_groups,
                     "log10_number_of_processive_groups": np.nan,
                     "radius": np.nan,
                     "avg_sims": np.nan,
                     "min_sims": np.nan,
                     "max_sims": np.nan,
                     "mean_sims": np.nan,
                     "std_sims": np.nan,
                     "pvalue": np.nan,
                     "qvalue": np.nan,
                     "minus_log10_qvalue": np.nan,
                     "zscore": np.nan,
                     "expected_number_of_processive_groups": expected_number_of_processive_groups}, ignore_index=True)


    ##########################################################################################
    ############################# Calculate p-values starts ##################################
    ##########################################################################################
    #p values
    all_p_values = []
    names_list = []

    for signature in sorted_signature_list:
        for processive_group_length in sorted_processsive_group_length_list:
            if signature_processive_group_length_properties_df[(signature_processive_group_length_properties_df['signature']==signature) &
                                                               (signature_processive_group_length_properties_df['processive_group_length']==processive_group_length)].values.any():

                observed_value = signature_processive_group_length_properties_df[
                    (signature_processive_group_length_properties_df['signature'] == signature) &
                    (signature_processive_group_length_properties_df['processive_group_length'] == processive_group_length)]['number_of_processive_groups'].values[0]

                expected_values = signature_processive_group_length_properties_df[
                    (signature_processive_group_length_properties_df['signature'] == signature) &
                    (signature_processive_group_length_properties_df['processive_group_length'] == processive_group_length)]['expected_number_of_processive_groups'].values[0]

                zscore = None
                if (not np.isnan(observed_value)) and (len(expected_values)>0 and np.count_nonzero(expected_values)>0):
                    # zstat, pvalue = ztest(expectedValues, value=observedValue) results in very small p-values therefore we are not calling in this way.
                    try:
                        zstat, pvalue = ztest(expected_values, [observed_value], alternative='smaller')
                    except FloatingPointError:
                        print(signature,' observed_value: ', observed_value, ' expected_values: ', expected_values, ' FloatingPointError: divide by zero encountered in double_scalars')

                    # Please note
                    # If pvalue is np.nan e.g.: due to a few expected values like only one [1]
                    # Then there must be cases when you may want to manually set minus_log10_qvalue to np.inf
                    if np.isnan(pvalue):
                        signature_processive_group_length_properties_df.loc[
                            ((signature_processive_group_length_properties_df['signature'] == signature) &
                             (signature_processive_group_length_properties_df['processive_group_length'] == processive_group_length)), 'minus_log10_qvalue'] = np.inf

                    if (pvalue is not None) and (not np.isnan(pvalue)):
                        all_p_values.append(pvalue)
                        names_list.append((signature,processive_group_length))

                    avg_sims=sum(expected_values)/len(expected_values)
                    min_sims=min(expected_values)
                    max_sims=max(expected_values)
                    mean_sims = np.mean(expected_values)
                    std_sims = np.std(expected_values)

                    if (std_sims > 0):
                        zscore = (observed_value - mean_sims)/std_sims

                    signature_processive_group_length_properties_df.loc[
                        ((signature_processive_group_length_properties_df['signature'] == signature) &
                        (signature_processive_group_length_properties_df['processive_group_length'] == processive_group_length)), 'avg_sims'] = avg_sims

                    signature_processive_group_length_properties_df.loc[
                        ((signature_processive_group_length_properties_df['signature'] == signature) &
                        (signature_processive_group_length_properties_df['processive_group_length'] == processive_group_length)), 'min_sims'] = min_sims

                    signature_processive_group_length_properties_df.loc[
                        ((signature_processive_group_length_properties_df['signature'] == signature) &
                        (signature_processive_group_length_properties_df['processive_group_length'] == processive_group_length)), 'max_sims'] = max_sims

                    signature_processive_group_length_properties_df.loc[
                        ((signature_processive_group_length_properties_df['signature'] == signature) &
                        (signature_processive_group_length_properties_df['processive_group_length'] == processive_group_length)), 'mean_sims'] = mean_sims

                    signature_processive_group_length_properties_df.loc[
                        ((signature_processive_group_length_properties_df['signature'] == signature) &
                        (signature_processive_group_length_properties_df['processive_group_length'] == processive_group_length)), 'std_sims'] = std_sims

                    signature_processive_group_length_properties_df.loc[
                        ((signature_processive_group_length_properties_df['signature'] == signature) &
                        (signature_processive_group_length_properties_df['processive_group_length'] == processive_group_length)), 'pvalue'] = pvalue

                    signature_processive_group_length_properties_df.loc[
                        ((signature_processive_group_length_properties_df['signature'] == signature) &
                        (signature_processive_group_length_properties_df['processive_group_length'] == processive_group_length)), 'zscore'] = zscore

                # elif (not np.isnan(observed_value)) and (observed_value >= MINIMUM_REQUIRED_NUMBER_OF_PROCESSIVE_GROUPS) and (len(expected_values) == 0):
                elif (not np.isnan(observed_value)) and ((len(expected_values)>0 and np.count_nonzero(expected_values)==0)):

                    # manually set
                    signature_processive_group_length_properties_df.loc[
                        ((signature_processive_group_length_properties_df['signature'] == signature) &
                         (signature_processive_group_length_properties_df['processive_group_length'] == processive_group_length)), 'minus_log10_qvalue'] = np.inf

                # elif (not np.isnan(observed_value)) and (observed_value >= MINIMUM_REQUIRED_NUMBER_OF_PROCESSIVE_GROUPS) and (np.count_nonzero(expected_values) == 0):
                elif (not np.isnan(observed_value)) and (len(expected_values) == 0):
                    # manually set
                    signature_processive_group_length_properties_df.loc[
                        ((signature_processive_group_length_properties_df['signature'] == signature) &
                         (signature_processive_group_length_properties_df['processive_group_length'] == processive_group_length)), 'minus_log10_qvalue'] = np.inf

    ##########################################################################################
    ############################# Calculate p-values ends ####################################
    ##########################################################################################

    all_p_values_array = np.asarray(all_p_values)
    all_FDR_BH_adjusted_p_values=None

    #FDR BH Multiple Testing Correction
    try:
        rejected, all_FDR_BH_adjusted_p_values, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(all_p_values_array, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    except ZeroDivisionError:
        print('ZeroDivisionError during statsmodels.stats.multitest.multipletests')
        print('all_p_values_array: %s' %(all_p_values_array))

    if all_FDR_BH_adjusted_p_values is not None:
        minus_log10_all_FDR_BH_adjusted_p_values = [-math.log10(q_value) if (q_value > 0) else np.inf for q_value in all_FDR_BH_adjusted_p_values]
    else:
        minus_log10_all_FDR_BH_adjusted_p_values = []

    # Get the corrected p values in an order
    for index, (signature,processive_group_length) in  enumerate(names_list,0):
        qvalue = all_FDR_BH_adjusted_p_values[index]
        minus_log10_qvalue = minus_log10_all_FDR_BH_adjusted_p_values[index]

        if signature_processive_group_length_properties_df[
            ((signature_processive_group_length_properties_df['signature'] == signature) &
                 (signature_processive_group_length_properties_df['processive_group_length'] == processive_group_length))].values.any():

            signature_processive_group_length_properties_df.loc[
                ((signature_processive_group_length_properties_df['signature'] == signature) &
                 (signature_processive_group_length_properties_df['processive_group_length'] == processive_group_length)), 'qvalue'] = qvalue

            signature_processive_group_length_properties_df.loc[
                ((signature_processive_group_length_properties_df['signature'] == signature) &
                 (signature_processive_group_length_properties_df['processive_group_length'] == processive_group_length)), 'minus_log10_qvalue'] = minus_log10_qvalue

    # modify starts
    # Filter the rows where processive_group_length >= MINIMUM_REQUIRED_PROCESSIVE_GROUP_LENGTH
    signature_processive_group_length_properties_df = signature_processive_group_length_properties_df[
        signature_processive_group_length_properties_df['processive_group_length'] >= MINIMUM_REQUIRED_PROCESSIVE_GROUP_LENGTH]

    signature_processive_group_length_properties_df['log10_number_of_processive_groups']=\
        np.log(signature_processive_group_length_properties_df['number_of_processive_groups'].replace(0,np.nan))

    # To show avg_number_of_processive_groups=1
    if MINIMUM_REQUIRED_NUMBER_OF_PROCESSIVE_GROUPS == 1:
        signature_processive_group_length_properties_df.loc[(signature_processive_group_length_properties_df['number_of_processive_groups'] == 1), 'log10_number_of_processive_groups'] = np.log(2) / 2

    # Here we set radius
    signature_processive_group_length_properties_df = \
        signature_processive_group_length_properties_df.groupby('signature').apply(lambda df: set_radius(df))
    # modify ends

    # Get the highest processive group length with a nonzero radius
    if (len(signature_processive_group_length_properties_df.index)>0):
        max_processive_group_length = signature_processive_group_length_properties_df[
            (round(signature_processive_group_length_properties_df['radius'],2)>0) &
            (signature_processive_group_length_properties_df['number_of_processive_groups'] >= MINIMUM_REQUIRED_NUMBER_OF_PROCESSIVE_GROUPS)]['processive_group_length'].max()

        # Update sorted_processsive_group_length_list
        processsive_group_length_list = signature_processive_group_length_properties_df[
            (round(signature_processive_group_length_properties_df['radius'], 2) > 0) &
            (signature_processive_group_length_properties_df['number_of_processive_groups'] >= MINIMUM_REQUIRED_NUMBER_OF_PROCESSIVE_GROUPS)]['processive_group_length'].unique()
        sorted_processsive_group_length_list = sorted(processsive_group_length_list, key=int)

        # Update sorted_signature_list
        signatures_list = signature_processive_group_length_properties_df[
            (round(signature_processive_group_length_properties_df['radius'], 2) > 0) &
            (signature_processive_group_length_properties_df['number_of_processive_groups'] >= MINIMUM_REQUIRED_NUMBER_OF_PROCESSIVE_GROUPS)]['signature'].unique()
        sorted_signature_list = sorted(signatures_list, reverse=True, key=natural_key)

        if verbose: print('\tVerbose #############################################')
        if verbose: print('\tVerbose len(all_p_values):%d\n all_p_values: %s' % (len(all_p_values), all_p_values))

        if verbose: print('\tVerbose #############################################')
        if verbose:
            if (all_FDR_BH_adjusted_p_values is not None): print('\tVerbose len(all_FDR_BH_adjusted_p_values):%d\n all_FDR_BH_adjusted_p_values: %s' % (len(all_FDR_BH_adjusted_p_values), all_FDR_BH_adjusted_p_values))
        if verbose: print('\tVerbose len(minus_log10_all_FDR_BH_adjusted_p_values):%d\n minus_log10_all_FDR_BH_adjusted_p_values:%s' % (len(minus_log10_all_FDR_BH_adjusted_p_values), minus_log10_all_FDR_BH_adjusted_p_values))

        # Plotting starts
        # create the directory if it does not exists
        os.makedirs(os.path.join(outputDir, jobname, FIGURE, PROCESSIVITY), exist_ok=True)
        os.makedirs(os.path.join(outputDir, jobname, FIGURE, PROCESSIVITY, TABLES), exist_ok=True)

        # Plot processivity figure
        plot_processivity_figure(outputDir,
                                 jobname,
                                 numberofSimulations,
                                 sorted_signature_list,
                                 sorted_processsive_group_length_list,
                                 max_processive_group_length,
                                 signature_processive_group_length_properties_df,
                                 verbose)

    # Write dataframe
    filePath = os.path.join(outputDir, jobname, FIGURE, PROCESSIVITY, TABLES, '%s_Signatures_Processivity.txt' % (jobname))
    signature_processive_group_length_properties_df.to_csv(filePath, sep='\t', header=True, index=False)

    # Append  all_mutations_loci_df' to 'signature_processive_group_length_properties_df'
    all_mutations_loci_df = pd.read_csv(os.path.join(outputDir, jobname, DATA, PROCESSIVITY, "Sim0_Processive_Mutations_Loci.txt"), sep='\t', header=0)

    f = open(filePath,"a")
    f.write("\n")
    all_mutations_loci_df.to_csv(f, sep='\t', header=True, index=False)
    f.close()
###################################################################



###################################################################
def plot_processivity_figure(outputDir,
                             jobname,
                             numberofSimulations,
                             sorted_signature_list,
                             sorted_processsive_group_length_list,
                             max_processive_group_length,
                             signature_processive_group_length_properties_df,
                             verbose):

    index = None
    if ((len(sorted_processsive_group_length_list) > 0) and (max_processive_group_length > 0)):
        # Find index of maxProcessiveGroupLength in sortedProcessiveGroupLengthList
        index = sorted_processsive_group_length_list.index(max_processive_group_length)
        if verbose: print('\tVerbose sortedProcessiveGroupLengthList[index]:%s' % (sorted_processsive_group_length_list[index]))
        if verbose: print('\tVerbose ##########################################')

    if verbose: print('\tVerbose maxProcessiveGroupLength:%d len(sortedSignatureList):%d ' % (max_processive_group_length, len(sorted_signature_list)))

    if (len(sorted_signature_list)>0):
        plot1, panel1 = plt.subplots(figsize=(20+1.5*len(sorted_processsive_group_length_list), 10+1.5*len(sorted_signature_list)))
        plt.rc('axes', edgecolor='lightgray')

        #make aspect ratio square
        panel1.set_aspect(1.0)

        #set title
        panel1.text(0.1, 1.2, jobname,horizontalalignment='center', verticalalignment='top', fontsize=60, fontweight='bold', fontname='Arial',transform=panel1.transAxes)

        #To get rid of  UserWarning: Attempting to set identical left==right results in singular transformations; automatically expanding.
        if (len(sorted_processsive_group_length_list)>1):
            panel1.set_xlim([1,index+1])
            panel1.set_xticks(np.arange(0,index+2,1))
        else:
            panel1.set_xlim([0,len(sorted_processsive_group_length_list)])
            panel1.set_xticks(np.arange(0,len(sorted_processsive_group_length_list)+1,1))

        if (len(sorted_signature_list)>1):
            panel1.set_ylim([1, len(sorted_signature_list)])
        else:
            panel1.set_ylim([0, len(sorted_signature_list)])

        panel1.set_yticks(np.arange(0, len(sorted_signature_list) + 1, 1))

        cmap = cm.get_cmap('YlOrRd')  # Looks better good
        v_min = 2
        v_max = 20
        #Very important: You have to normalize
        norm = mpl.colors.Normalize(vmin=v_min, vmax=v_max)

        if not signature_processive_group_length_properties_df.empty:
            # Plot the circles with color
            for signature_index, signature in enumerate(sorted_signature_list):
                for processive_group_length_index, processive_group_length in enumerate(sorted_processsive_group_length_list):
                    number_of_processive_groups=np.nan
                    radius = np.nan
                    color = np.nan

                    if (signature_processive_group_length_properties_df[
                        (signature_processive_group_length_properties_df['signature'] == signature) &
                        (signature_processive_group_length_properties_df['processive_group_length'] == processive_group_length)]['number_of_processive_groups'].values.any()):

                        number_of_processive_groups= signature_processive_group_length_properties_df[
                            (signature_processive_group_length_properties_df['signature'] == signature) &
                            (signature_processive_group_length_properties_df['processive_group_length'] == processive_group_length)]['number_of_processive_groups'].values[0]

                    if (signature_processive_group_length_properties_df[
                        (signature_processive_group_length_properties_df['signature'] == signature) &
                        (signature_processive_group_length_properties_df['processive_group_length'] == processive_group_length)]['radius'].values.any()):

                        radius = signature_processive_group_length_properties_df[
                            (signature_processive_group_length_properties_df['signature'] == signature) &
                            (signature_processive_group_length_properties_df['processive_group_length'] == processive_group_length)]['radius'].values[0]

                    if (signature_processive_group_length_properties_df[
                        (signature_processive_group_length_properties_df['signature'] == signature) &
                        (signature_processive_group_length_properties_df['processive_group_length'] == processive_group_length)]['minus_log10_qvalue'].values.any()):

                        color = signature_processive_group_length_properties_df[
                            (signature_processive_group_length_properties_df['signature'] == signature) &
                            (signature_processive_group_length_properties_df['processive_group_length'] == processive_group_length)]['minus_log10_qvalue'].values[0]

                    if ((not np.isnan(number_of_processive_groups)) and (number_of_processive_groups >= 5)
                            and (not np.isnan(radius)) and (radius > 0) and (not np.isnan(color))):
                        # Very important: You have to norm
                        circle = plt.Circle((processive_group_length_index + 0.5, signature_index + 0.5), radius,color=cmap(norm(color)), fill=True)
                        panel1.add_patch(circle)
                    elif ((not np.isnan(number_of_processive_groups)) and (number_of_processive_groups >= 5)
                          and (not np.isnan(radius)) and (radius > 0) and np.isnan(color)):
                        circle = plt.Circle((processive_group_length_index + 0.5, signature_index + 0.5), radius,color="g", fill=True)
                        panel1.add_patch(circle)

        panel1.set_facecolor('white')
        #When there are subplots, this is needed.
        panel1.grid(color='black')

        for edge, spine in panel1.spines.items():
            spine.set_visible(True)
            spine.set_color('black')

        xlabels=None
        if (index is not None):
            xlabels = sorted_processsive_group_length_list[0:index+1]
        ylabels = sorted_signature_list

        # Put the color bar if there are simulations
        if (numberofSimulations>0):
            #Vertical/Horizontal Colorbar
            # Vertical volorbar to the right
            cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap))  # this works because of the scatter
            # cb.ax.set_ylabel("-log10 (q-value)", va="bottom", rotation=-90, labelpad=25)
            cb.ax.set_ylabel("-log10 (q-value)", fontsize=50, labelpad=25)

            # Horizontal colorbar to the bottom
            # cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), orientation='horizontal')  # this works because of the scatter
            # cb.ax.set_xlabel("colorbar label", fontsize=50, labelpad=25)

            #common for horizontal colorbar and vertical colorbar
            cbax = cb.ax
            cbax.tick_params(labelsize=40)
            text_x = cbax.xaxis.label
            text_y = cbax.yaxis.label
            font = mpl.font_manager.FontProperties(size=40)
            text_x.set_font_properties(font)
            text_y.set_font_properties(font)

        # CODE GOES HERE TO CENTER X-AXIS LABELS...
        panel1.set_xticklabels([])
        mticks = panel1.get_xticks()
        panel1.set_xticks((mticks[:-1] + mticks[1:]) / 2, minor=True)
        panel1.tick_params(axis='x', which='minor', length=0,labelsize=50)

        if xlabels is not None:
            panel1.set_xticklabels(xlabels, minor=True)

        panel1.xaxis.set_ticks_position('top')

        plt.tick_params(
            axis='x',  # changes apply to the x-axis
            which='major',  # both major and minor ticks are affected
            bottom=False,  # ticks along the bottom edge are off
            top=False)  # labels along the bottom edge are off

        # CODE GOES HERE TO CENTER Y-AXIS LABELS...
        panel1.set_yticklabels([])
        mticks = panel1.get_yticks()
        panel1.set_yticks((mticks[:-1] + mticks[1:]) / 2, minor=True)
        panel1.tick_params(axis='y', which='minor', length=0,labelsize=50)
        panel1.set_yticklabels(ylabels, minor=True) # fontsize

        plt.tick_params(
            axis='y',  # changes apply to the x-axis
            which='major',  # both major and minor ticks are affected
            left=False)  # labels along the bottom edge are off

        filename = '%s_Processivity.png' %(jobname)
        figFile = os.path.join(outputDir, jobname, FIGURE, PROCESSIVITY, filename)
        plot1.savefig(figFile,dpi=100, bbox_inches="tight")

        plt.cla()
        plt.close(plot1)
###################################################################

###################################################################
def processivityFigures(outputDir,jobname,numberofSimulations,verbose):

    jobnamePath = os.path.join(outputDir,jobname,FIGURE,PROCESSIVITY)
    if verbose: print('\tVerbose Topography.py jobnamePath:%s ' %jobnamePath)

    processivity_table_file_list=[]
    processivity_df_list=[]
    #Fill signature_processive_group_length_number_of_processive_groups_median_number_of_processive_groups_in_MB_df

    for simNum in range(0,numberofSimulations+1):
        filename="Sim%d_Processivity.txt" %(simNum)
        filepath=os.path.join(outputDir,jobname,DATA,PROCESSIVITY,filename)
        if os.path.exists(filepath):
            processivity_table_file_list.append((simNum,filepath))

    for (simNum,processivity_table_file) in processivity_table_file_list:
        processivity_df = pd.read_csv(processivity_table_file, header=0, sep='\t')
        processivity_df['Simulation_Number']=simNum
        processivity_df=processivity_df[['Simulation_Number',
                                         'Signature',
                                         'Processsive_Group_Length',
                                         'Number_of_Processive_Groups',
                                         'Median_Distance_Between_Last_First_Mutations',
                                         'Median_Distance_Between_Consecutive_Mutations',
                                         'Median_Number_of_Mutations_Within_1MB']]
        processivity_df_list.append(processivity_df)

    # Vertically combine dfs
    processivity_df = pd.concat(processivity_df_list)

    # Using dataframes
    plotRelationshipBetweenSignaturesandProcessiveGroupLengthsUsingDataframes(outputDir,jobname,processivity_df,numberofSimulations,verbose)
###################################################################