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

SIGNIFICANCE_LEVEL=0.01

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
def plotRelationshipBetweenSignaturesandProcessiveGroupLengthsUsingDataframes(outputDir,jobname,processivity_df,numberofSimulations,verbose):

    #processivity_df
    #'Simulation_Number', 'Signature', 'Processsive_Group_Length', 'Number_of_Processive_Groups', 'Median_of_Number_of_Processive_Groups_in_MB'

    #This dictionary is filled and written as a dataframe at the end.
    signature2ProcessiveGroupLength2ProcessiveGroupPropertiesDict={}

    ####################################################
    #Get the list of signatures using original
    signatures_list=processivity_df[processivity_df['Simulation_Number']==0]['Signature'].unique()
    sorted_signature_list = sorted(signatures_list,reverse=True,key=natural_key)
    ####################################################

    ####################################################
    #Get the list of processive group lengths using original
    processsive_group_length_list = processivity_df[processivity_df['Simulation_Number']==0]['Processsive_Group_Length'].unique()
    sorted_processsive_group_length_list = sorted(processsive_group_length_list, key=int)
    ####################################################

    ###################################################################
    ############### Fill this radius dataframe starts #################
    ###################################################################
    list_of_lists = []

    for signature in sorted_signature_list:
        # Normalize number of groups within a signature
        numberofProcessiveGroupsList = []
        numberofProcessiveGroupsList_in_log10 = []

        for processiveGroupLength in sorted_processsive_group_length_list:
            if (processivity_df[(processivity_df['Simulation_Number']==0) & (processivity_df['Signature']==signature) & (processivity_df['Processsive_Group_Length']==processiveGroupLength)]['Number_of_Processive_Groups'].values.size>0):
                number_of_processive_groups= processivity_df[(processivity_df['Simulation_Number']==0) & (processivity_df['Signature']==signature) & (processivity_df['Processsive_Group_Length']==processiveGroupLength)]['Number_of_Processive_Groups'].values[0]
                if (number_of_processive_groups >= 5):
                    numberofProcessiveGroupsList.append(number_of_processive_groups)
                    numberofProcessiveGroupsList_in_log10.append(math.log10(number_of_processive_groups))
                else:
                    numberofProcessiveGroupsList.append(0)
                    numberofProcessiveGroupsList_in_log10.append(0)

        # If numberofProcessiveGroupsList is not empty list
        if (numberofProcessiveGroupsList_in_log10 and (max(numberofProcessiveGroupsList_in_log10) > 0)):
            normalizedNumberofProcessiveGroupsList = [i / max(numberofProcessiveGroupsList_in_log10) for i in numberofProcessiveGroupsList_in_log10]
            radiusNumberofProcessiveGroupsList = [i * 0.48 for i in normalizedNumberofProcessiveGroupsList]

            # Fill dictionary using radiusNormalizedNumberofProcessiveGroupsList
            radiusIndex = 0
            for processiveGroupLength in sorted_processsive_group_length_list:
                if (processivity_df[(processivity_df['Simulation_Number']==0) & (processivity_df['Signature']==signature) & (processivity_df['Processsive_Group_Length'] == processiveGroupLength)]['Number_of_Processive_Groups'].values.size > 0):
                    number_of_processive_groups=processivity_df[(processivity_df['Simulation_Number'] == 0) & (processivity_df['Signature'] == signature) & (processivity_df['Processsive_Group_Length'] == processiveGroupLength)]['Number_of_Processive_Groups'].values[0]
                    radius = radiusNumberofProcessiveGroupsList[radiusIndex]
                    list_of_lists.append([signature, processiveGroupLength, number_of_processive_groups, radius])
                    radiusIndex += 1

    signature_radius_df=pd.DataFrame(list_of_lists, columns=['Signature', 'Processsive_Group_Length', 'Number_of_Processive_Groups', 'Radius'])
    ###################################################################
    ############### Fill this radius dataframe ends ###################
    ###################################################################


    ##########################################################################################
    ############################# Calculate p-values starts ##################################
    ##########################################################################################
    # p values
    all_p_values = []
    all_p_values_element_names = []

    for signature in sorted_signature_list:
        for processiveGroupLength in sorted_processsive_group_length_list:

            if (signature_radius_df[(signature_radius_df['Signature']==signature) & (signature_radius_df['Processsive_Group_Length']==processiveGroupLength)]['Radius'].values.size>0):
                radius=signature_radius_df[(signature_radius_df['Signature']==signature) & (signature_radius_df['Processsive_Group_Length']==processiveGroupLength)]['Radius'].values[0]

            if (processivity_df[(processivity_df['Simulation_Number']==0) & (processivity_df['Signature']==signature) & (processivity_df['Processsive_Group_Length']==processiveGroupLength)]['Number_of_Processive_Groups'].values.size>0):
                observedValue=processivity_df[(processivity_df['Simulation_Number'] == 0) & (processivity_df['Signature'] == signature) & (processivity_df['Processsive_Group_Length'] == processiveGroupLength)]['Number_of_Processive_Groups'].values[0]

                expectedValues = []
                for simNum in range(1,numberofSimulations+1):
                    if (processivity_df[(processivity_df['Simulation_Number'] == simNum) & (processivity_df['Signature'] == signature) & (processivity_df['Processsive_Group_Length'] == processiveGroupLength)]['Number_of_Processive_Groups'].values.size > 0):
                        expectedValues.append(processivity_df[(processivity_df['Simulation_Number'] == simNum) & (processivity_df['Signature'] == signature) & (processivity_df['Processsive_Group_Length'] == processiveGroupLength)]['Number_of_Processive_Groups'].values[0])
                    else:
                        expectedValues.append(0)

                ############################################################################
                mean_sims = None
                min_sims = None
                max_sims = None
                std_sims = None
                pvalue = None
                zscore = None

                if expectedValues and (len(expectedValues) > 0) and (np.any(expectedValues)):
                    # zstat, pvalue = ztest(expectedValues, value=observedValue) results in very small p-values therefore we are not calling in this way.
                    try:
                        # zstat, pvalue = ztest(expectedValues, [observedValue])
                        #alternative hypothesis: difference between mean of expected values and mean of observed values is less than zero.
                        zstat, pvalue = ztest(expectedValues, [observedValue],alternative = 'smaller')
                    except FloatingPointError:
                        print('FloatingPointError zstat, pvalue = ztest(%s, [%d])' % (expectedValues, observedValue))
                    # if pvalue is None or np.nan qvalue can not be calculated.
                    if (pvalue is not None) and (not np.isnan(pvalue)):
                        all_p_values.append(pvalue)
                        all_p_values_element_names.append((signature, processiveGroupLength))

                if (expectedValues and len(expectedValues) > 0):
                    mean_sims = np.mean(expectedValues)
                    min_sims = np.min(expectedValues)
                    max_sims = np.max(expectedValues)
                    std_sims = np.std(expectedValues)

                    if (std_sims > 0):
                        zscore = (observedValue - mean_sims) / std_sims

                processiveGroupProperties = {}
                processiveGroupProperties['processive_group_length'] = processiveGroupLength
                processiveGroupProperties['number_of_processive_groups'] = observedValue
                processiveGroupProperties['radius'] = radius
                processiveGroupProperties['mean_sims'] = mean_sims
                processiveGroupProperties['min_sims'] = min_sims
                processiveGroupProperties['max_sims'] = max_sims
                processiveGroupProperties['std_sims'] = std_sims
                processiveGroupProperties['pvalue'] = pvalue
                processiveGroupProperties['zscore'] = zscore
                processiveGroupProperties['expectedValues'] = expectedValues

                if signature in signature2ProcessiveGroupLength2ProcessiveGroupPropertiesDict:
                    signature2ProcessiveGroupLength2ProcessiveGroupPropertiesDict[signature][processiveGroupLength] = processiveGroupProperties
                else:
                    signature2ProcessiveGroupLength2ProcessiveGroupPropertiesDict[signature] = {}
                    signature2ProcessiveGroupLength2ProcessiveGroupPropertiesDict[signature][processiveGroupLength] = processiveGroupProperties
                ############################################################################

    ##########################################################################################
    ############################# Calculate p-values ends ####################################
    ##########################################################################################


    ####################################################################################
    ############################# Calculate q-values starts ############################
    ####################################################################################

    ##########################################################################################
    all_p_values_array = np.asarray(all_p_values)
    all_FDR_BH_adjusted_p_values = None

    # FDR BH Multiple Testing Correction
    try:
        rejected, all_FDR_BH_adjusted_p_values, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(all_p_values_array, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
    except ZeroDivisionError:
        print('ZeroDivisionError during statsmodels.stats.multitest.multipletests')
        print('all_p_values_array: %s' %(all_p_values_array))
    ##########################################################################################


    if all_FDR_BH_adjusted_p_values is not None:
        minus_log10_all_FDR_BH_adjusted_p_values = [-math.log10(q_value) if (q_value > 0 and q_value < SIGNIFICANCE_LEVEL) else np.nan for q_value in all_FDR_BH_adjusted_p_values]
    else:
        minus_log10_all_FDR_BH_adjusted_p_values = []

    if verbose: print('\tVerbose #############################################')
    if verbose: print('\tVerbose len(all_p_values):%d\n all_p_values: %s' % (len(all_p_values), all_p_values))

    if verbose: print('\tVerbose #############################################')
    if verbose:
        if (all_FDR_BH_adjusted_p_values is not None):
            print('\tVerbose len(all_FDR_BH_adjusted_p_values):%d\n all_FDR_BH_adjusted_p_values: %s' % (len(all_FDR_BH_adjusted_p_values), all_FDR_BH_adjusted_p_values))

    if verbose: print('\tVerbose #############################################')
    if verbose: print('\tVerbose len(minus_log10_all_FDR_BH_adjusted_p_values):%d\n minus_log10_all_FDR_BH_adjusted_p_values:%s' % (len(minus_log10_all_FDR_BH_adjusted_p_values), minus_log10_all_FDR_BH_adjusted_p_values))
    if verbose: print('\tVerbose #############################################')

    #######################################################################################
    #######################  Get the corrected q values in an order starts ################
    #######################################################################################
    for element_index, all_p_values_element_name in enumerate(all_p_values_element_names):
        q_value = all_FDR_BH_adjusted_p_values[element_index]
        (signature, processiveGroupLength) = all_p_values_element_name

        signature2ProcessiveGroupLength2ProcessiveGroupPropertiesDict[signature][processiveGroupLength]['qvalue'] = q_value
        if q_value is not None and q_value == 0:
            signature2ProcessiveGroupLength2ProcessiveGroupPropertiesDict[signature][processiveGroupLength]['minus_log10_qvalue'] = np.inf
        elif q_value is not None:
            signature2ProcessiveGroupLength2ProcessiveGroupPropertiesDict[signature][processiveGroupLength]['minus_log10_qvalue'] = -math.log10(q_value)
        else:
            signature2ProcessiveGroupLength2ProcessiveGroupPropertiesDict[signature][processiveGroupLength]['minus_log10_qvalue'] = None
    #######################################################################################
    #######################  Get the corrected q values in an order ends ##################
    #######################################################################################

    # For the (signature,processiveGroupLength) tuples where pvalue can not be calculated using ztest e.g.: expected values are all zero
    # We need to set qvalue and minus_log10_qvalue as None
    # So that writeDictionaryAsADataframe does not give error.
    for signature in signature2ProcessiveGroupLength2ProcessiveGroupPropertiesDict:
        for processiveGroupLength in signature2ProcessiveGroupLength2ProcessiveGroupPropertiesDict[signature]:
            if (signature, processiveGroupLength) not in all_p_values_element_names:
                signature2ProcessiveGroupLength2ProcessiveGroupPropertiesDict[signature][processiveGroupLength]['qvalue'] = None
                signature2ProcessiveGroupLength2ProcessiveGroupPropertiesDict[signature][processiveGroupLength]['minus_log10_qvalue'] = None

    if verbose: print('\tVerbose ################################################################')
    ####################################################################################
    ############################# Calculate q-values ends ############################
    ####################################################################################


    ###################################################################
    ############### For information starts ############################
    ###################################################################
    # Get the highest processive group length with a nonzero radius
    #'Signature', 'Processsive_Group_Length', 'Number_of_Processive_Groups', 'Radius'

    processive_group_length_nparray = signature_radius_df[(signature_radius_df['Radius'] > 0)]['Processsive_Group_Length'].values

    if (processive_group_length_nparray.size>0):
        maxProcessiveGroupLength= max(processive_group_length_nparray)
        if verbose: print('\tVerbose Processivity plot will be for this maxProcessiveGroupLength:%d' % (maxProcessiveGroupLength))
        if verbose: print('\tVerbose len(sortedProcessiveGroupLengthList):%d' % (len(sorted_processsive_group_length_list)))

    ###################################################################
    ############### For information ends ##############################
    ###################################################################

    ###################################################################
    ############### Plotting starts ###################################
    ###################################################################
    if (processive_group_length_nparray.size>0):
        # create the directory if it does not exists
        os.makedirs(os.path.join(outputDir, jobname, FIGURE, PROCESSIVITY), exist_ok=True)
        os.makedirs(os.path.join(outputDir, jobname, FIGURE, PROCESSIVITY, TABLES), exist_ok=True)

        #plot processivity figure
        plot_processivity_figure(outputDir,
                                 jobname,
                                 numberofSimulations,
                                 sorted_signature_list,
                                 sorted_processsive_group_length_list,
                                 maxProcessiveGroupLength,
                                 signature_radius_df,
                                 signature2ProcessiveGroupLength2ProcessiveGroupPropertiesDict,
                                 verbose)

        #write accompanying processivity file
        filePath = os.path.join(outputDir, jobname, FIGURE, PROCESSIVITY, TABLES, '%s_Signatures_Processivity.txt' % (jobname))
        writeDictionaryAsADataframe(jobname, signature2ProcessiveGroupLength2ProcessiveGroupPropertiesDict, filePath)
    ###################################################################
    ############### Plotting ends #####################################
    ###################################################################

###################################################################



###################################################################
def plot_processivity_figure(outputDir,
                             jobname,
                             numberofSimulations,
                             sorted_signature_list,
                             sorted_processsive_group_length_list,
                             maxProcessiveGroupLength,
                             signature_radius_df,
                             signature2ProcessiveGroupLength2ProcessiveGroupPropertiesDict,
                             verbose):

    index = None
    if ((len(sorted_processsive_group_length_list) > 0) and (maxProcessiveGroupLength > 0)):
        # Find index of maxProcessiveGroupLength in sortedProcessiveGroupLengthList
        index = sorted_processsive_group_length_list.index(maxProcessiveGroupLength)
        if verbose: print('\tVerbose sortedProcessiveGroupLengthList[index]:%s' % (sorted_processsive_group_length_list[index]))
        if verbose: print('\tVerbose ##########################################')

    if verbose: print('\tVerbose maxProcessiveGroupLength:%d len(sortedSignatureList):%d ' % (maxProcessiveGroupLength, len(sorted_signature_list)))

    if (len(sorted_signature_list)>0):
        plot1, panel1 = plt.subplots(figsize=(20+1.5*len(sorted_processsive_group_length_list), 10+1.5*len(sorted_signature_list)))
        plt.rc('axes', edgecolor='lightgray')

        #make aspect ratio square
        panel1.set_aspect(1.0)

        #set title
        # panel1.text(len(sorted_processsive_group_length_list)*3, len(sorted_signature_list)+2.5, jobname,  horizontalalignment='center', fontsize=60, fontweight='bold', fontname='Arial')
        panel1.text(0.1, 1.2, jobname,horizontalalignment='center', verticalalignment='top', fontsize=60, fontweight='bold', fontname='Arial',transform=panel1.transAxes)

        #To get rid of  UserWarning: Attempting to set identical left==right results in singular transformations; automatically expanding.
        if (len(sorted_processsive_group_length_list)>1):
            panel1.set_xlim([1,index+1])
            panel1.set_xticks(np.arange(0,index+2,1))
        else:
            panel1.set_xlim([0,len(sorted_processsive_group_length_list)])
            panel1.set_xticks(np.arange(0,len(sorted_processsive_group_length_list),1))

        if (len(sorted_signature_list)>1):
            panel1.set_ylim([1, len(sorted_signature_list)])
        else:
            panel1.set_ylim([0, len(sorted_signature_list)])

        panel1.set_yticks(np.arange(0, len(sorted_signature_list) + 1, 1))
        #######################################################################

        #######################################################################
        cmap = cm.get_cmap('YlOrRd')  # Looks better good
        v_min = 2
        v_max = 20
        #Very important: You have to normalize
        norm = mpl.colors.Normalize(vmin=v_min, vmax=v_max)
        #######################################################################

        ##########################################################################################
        if (numberofSimulations>0):
            #Plot the circles with color
            for sigIndex, signature in enumerate(sorted_signature_list):
                for lengthIndex, processiveGroupLength in enumerate(sorted_processsive_group_length_list):
                    if (signature_radius_df[(signature_radius_df['Signature'] == signature) & (signature_radius_df['Processsive_Group_Length'] == processiveGroupLength)]['Radius'].values > 0):
                        radius = signature_radius_df[(signature_radius_df['Signature'] == signature) & (signature_radius_df['Processsive_Group_Length'] == processiveGroupLength)]['Radius'].values[0]
                        color=None
                        if 'minus_log10_qvalue' in signature2ProcessiveGroupLength2ProcessiveGroupPropertiesDict[signature][processiveGroupLength]:
                            color=signature2ProcessiveGroupLength2ProcessiveGroupPropertiesDict[signature][processiveGroupLength]['minus_log10_qvalue']

                        if ((radius is not None) and (radius>0) and color):
                            #Very important: You have to norm
                            circle = plt.Circle((lengthIndex + 0.5, sigIndex + 0.5), radius, color=cmap(norm(color)), fill=True)
                            panel1.add_patch(circle)

        else:
            #There is  no simulation data therefore no p values
            #Plot the circles without color
            for sigIndex, signature in enumerate(sorted_signature_list):
                for lengthIndex, processiveGroupLength in enumerate(sorted_processsive_group_length_list):
                    if (signature_radius_df[(signature_radius_df['Signature'] == signature) & (signature_radius_df['Processsive_Group_Length'] == processiveGroupLength)]['Radius'].values > 0):
                        radius = signature_radius_df[(signature_radius_df['Signature'] == signature) & (signature_radius_df['Processsive_Group_Length'] == processiveGroupLength)]['Radius'].values[0]
                        circle = plt.Circle((lengthIndex + 0.5, sigIndex + 0.5),radius,color="g", fill=True)
                        panel1.add_patch(circle)
        ##########################################################################################

        panel1.set_facecolor('white')
        #When there are subplots, this is needed.
        panel1.grid(color='black')

        panel1.grid(color='black')

        for edge, spine in panel1.spines.items():
            spine.set_visible(True)
            spine.set_color('black')

        xlabels=None
        if (index is not None):
            xlabels = sorted_processsive_group_length_list[0:index+1]
        ylabels = sorted_signature_list
        ##########################################################################################

        ################### Put the color bar if there are simulations starts ###################
        if (numberofSimulations>0):

            #########################################################################################
            #Vertical Colorbar
            #Used for scatter plot
            x = []
            y = []
            c = []

            for lengthIndex, processiveGroupLength in enumerate(sorted_processsive_group_length_list):
                x.append(lengthIndex)
                y.append(lengthIndex)
                c.append(0.5)

            #This code defines the ticks on the color bar
            # plot the scatter plot
            sc = plt.scatter(x, y, s=0, c=c, cmap=cmap, vmin=v_min, vmax=v_max, edgecolors='black')
            cb = plt.colorbar(sc)  # this works because of the scatter

            # cb.ax.set_ylabel("-log10 (q-value)", va="bottom", rotation=-90, labelpad=25)
            cb.ax.set_ylabel("-log10 (q-value)", fontsize=50, labelpad=25)
            #########################################################################################

            ##########################################################################################
            #common for horizontal colorbar and vertical colorbar
            cbax = cb.ax
            cbax.tick_params(labelsize=40)
            text_x = cbax.xaxis.label
            text_y = cbax.yaxis.label
            font = mpl.font_manager.FontProperties(size=40)
            text_x.set_font_properties(font)
            text_y.set_font_properties(font)
            ##########################################################################################

        ################### Put the color bar if there are simulations ends #####################

        ##################################################################################
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
        ##################################################################################

        ##################################################################################
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
        ##################################################################################

        ##################################################################################
        filename = '%s_Processivity.png' %(jobname)
        figFile = os.path.join(outputDir, jobname, FIGURE, PROCESSIVITY, filename)
        plot1.savefig(figFile,dpi=100, bbox_inches="tight")
        # plot1.tight_layout()

        plt.cla()
        plt.close(plot1)
        ##################################################################################

###################################################################


####################################################################################
# processiveGroupProperties['processive_group_length'] = processiveGroupLength
# processiveGroupProperties['number_of_processive_groups'] = observedValue
# processiveGroupProperties['mean_sims'] = mean_sims
# processiveGroupProperties['min_sims'] = min_sims
# processiveGroupProperties['max_sims'] = max_sims
# processiveGroupProperties['std_sims'] = std_sims
# processiveGroupProperties['pvalue'] = pvalue
# processiveGroupProperties['minus_log10_pvalue'] = minus_log10_pvalue
# processiveGroupProperties['zscore'] = zscore
# processiveGroupProperties['expectedValues'] = expectedValues
def writeDictionaryAsADataframe(jobname,signature2ProcessiveGroupLength2ProcessiveGroupPropertiesDict,filePath):
    L = sorted([(jobname,signature, processiveGroupLength,v1['number_of_processive_groups'],v1['radius'],v1['mean_sims'],v1['min_sims'],v1['max_sims'],v1['std_sims'],v1['pvalue'],v1['qvalue'],v1['minus_log10_qvalue'],v1['zscore'],v1['expectedValues']) for signature, v in signature2ProcessiveGroupLength2ProcessiveGroupPropertiesDict.items() for processiveGroupLength, v1 in v.items()])
    df = pd.DataFrame(L,columns=['tissue','signature', 'processsive_group_length','number_of_processive_groups','radius','mean_sims','min_sims','max_sims','std_sims','pvalue','qvalue','minus_log10_qvalue','zscore','expectedValues'])

    #write this dataframe
    df.to_csv(filePath, sep='\t', header=True, index=False)
####################################################################################


###################################################################
def processivityFigures(outputDir,jobname,numberofSimulations,verbose):

    jobnamePath = os.path.join(outputDir,jobname,FIGURE,PROCESSIVITY)
    if verbose: print('\tVerbose Topography.py jobnamePath:%s ' %jobnamePath)

    ############################################################
    processivity_table_file_list=[]
    processivity_df_list=[]
    #Fill signature_processive_group_length_number_of_processive_groups_median_number_of_processive_groups_in_MB_df

    for simNum in range(0,numberofSimulations+1):
        filename="Sim%d_Processivity.txt" %(simNum)
        filepath=os.path.join(outputDir,jobname,DATA,PROCESSIVITY,filename)
        if os.path.exists(filepath):
            processivity_table_file_list.append((simNum,filepath))

    for (simNum,processivity_table_file) in processivity_table_file_list:
        #Signature       Processsive_Group_Length        Number_of_Processive_Groups     Median_of_Number_of_Processive_Groups_in_MB
        processivity_df=pd.read_csv(processivity_table_file, header=0, sep='\t')
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
    #signature_processive_group_length_number_of_processive_groups_median_number_of_processive_groups_in_MB_df
    processivity_df = pd.concat(processivity_df_list)

    #New way using dataframes
    plotRelationshipBetweenSignaturesandProcessiveGroupLengthsUsingDataframes(outputDir,jobname,processivity_df,numberofSimulations,verbose)
    ############################################################

###################################################################
