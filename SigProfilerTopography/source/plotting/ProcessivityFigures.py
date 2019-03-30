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

import scipy.stats
import statsmodels.stats.multitest

import matplotlib
BACKEND = 'Agg'
if matplotlib.get_backend().lower() != BACKEND.lower():
    # If backend is not set properly a call to describe will hang
    matplotlib.use(BACKEND)

from matplotlib import pyplot as plt

import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from scipy.stats import poisson
from SigProfilerTopography.source.commons.TopographyCommons import *

###################################################################
def readSimulationBasedDictionaries(outputDir,jobname,numberofSimulations):
    simulation2Signature2ProcessiveGroupLength2PropertiesDict = {}

    for simNum in range(1,numberofSimulations+1):
        simJobName = '%s_Sim%d' % (jobname, simNum)
        simulationFilePath = os.path.join(outputDir,simJobName,DATA,PROCESSIVITY,'Signature2ProcessiveGroupLength2PropertiesDict.txt')

        signature2ProcessiveGroupLength2PropertiesDict =  readDictionary(simulationFilePath)

        if (signature2ProcessiveGroupLength2PropertiesDict is not None):
            simulation2Signature2ProcessiveGroupLength2PropertiesDict[simNum] = signature2ProcessiveGroupLength2PropertiesDict

    return simulation2Signature2ProcessiveGroupLength2PropertiesDict
###################################################################



#consider p values seems ro change nothing? can it be? Check this out.

###################################################################
def plotRelationshipBetweenSignaturesandProcessiveGroupLengths(outputDir,jobname,originalSignature2ProcessiveGroupLength2PropertiesDict,
                                                                simulation2Signature2ProcessiveGroupLength2PropertiesDict,
                                                                multipleTestingCorrection,
                                                                pValueCalculation):
    print('########################################')
    print('pValueCalculation:%s\tmultipleTestingCorrection:%s' %(pValueCalculation,multipleTestingCorrection))
    print('########################################')

    ####################################################
    numberofCorrectedPValuesGreaterThanZero = 0
    numberofCorrectedPValuesEqualsToZero = 0
    numberofCorrectedPValuesLessThanZero = 0
    ####################################################

    ####################################################
    #Get the list of signatures using original
    signature_dict_keys = originalSignature2ProcessiveGroupLength2PropertiesDict.keys()
    sortedSignatureList = sorted(signature_dict_keys, reverse=True)
    ####################################################

    ####################################################
    #Get the list of processive group lengths using original
    processiveGroupLengthSet = set()
    for signature in signature_dict_keys:
        processiveGroupLengthSet = processiveGroupLengthSet.union(set(originalSignature2ProcessiveGroupLength2PropertiesDict[signature].keys()))

    #Convert set to list
    processiveGroupLengthList = list(processiveGroupLengthSet)
    sortedProcessiveGroupLengthList = sorted(processiveGroupLengthList, key=int)
    ####################################################


    ####################################################
    #Fill the number of cells in the table
    numberofMultipleTests = 0
    ####################################################

    ####################################################
    #Fill this radius dictionary
    signature2ProcessiveGroupLength2RadiusDict = {}

    for signature in sortedSignatureList:
        #Normalize number of groups within in signature
        numberofProcessiveGroupsList = []
        for processiveGroupLength in sortedProcessiveGroupLengthList:
            if signature in originalSignature2ProcessiveGroupLength2PropertiesDict:
                if processiveGroupLength in originalSignature2ProcessiveGroupLength2PropertiesDict[signature]:
                    numberofMultipleTests +=1
                    #Take log10
                    numberofProcessiveGroupsList.append(math.log10(originalSignature2ProcessiveGroupLength2PropertiesDict[signature][processiveGroupLength]['numberofProcessiveGroups']))

        #If numberofProcessiveGroupsList is not empty list
        if (numberofProcessiveGroupsList and (max(numberofProcessiveGroupsList)>0)):
            normalizedNumberofProcessiveGroupsList = [float(i) / max(numberofProcessiveGroupsList) for i in numberofProcessiveGroupsList]
            atMost0Point9NumberofProcessiveGroupsList = [i*0.9 for i in normalizedNumberofProcessiveGroupsList]
            radiusNumberofProcessiveGroupsList = [i/2 for i in atMost0Point9NumberofProcessiveGroupsList]

            #Fill dictionary using radiusNormalizedNumberofProcessiveGroupsList
            radiusIndex = 0
            for processiveGroupLength in sortedProcessiveGroupLengthList:
                if signature in originalSignature2ProcessiveGroupLength2PropertiesDict:
                    if processiveGroupLength in originalSignature2ProcessiveGroupLength2PropertiesDict[signature]:
                        radius = radiusNumberofProcessiveGroupsList[radiusIndex]
                        if signature in signature2ProcessiveGroupLength2RadiusDict:
                            signature2ProcessiveGroupLength2RadiusDict[signature][processiveGroupLength] = radius
                        else:
                            signature2ProcessiveGroupLength2RadiusDict[signature] = {}
                            signature2ProcessiveGroupLength2RadiusDict[signature][processiveGroupLength] = radius
                        radiusIndex +=1

    ####################################################


    ####################################################
    #Are these dummy variables? Yes.
    x = []
    y = []
    c = []
    ####################################################

    # create a new figure
    fig = plt.figure(figsize=(45, 15), dpi=300)
    ax = plt.gca()

    #######################################################################
    #To get rid of  UserWarning: Attempting to set identical left==right results in singular transformations; automatically expanding.
    if (len(sortedProcessiveGroupLengthList)>1):
        plt.xlim([1, len(sortedProcessiveGroupLengthList)])
    else:
        plt.xlim([0, len(sortedProcessiveGroupLengthList)])

    if (len(sortedSignatureList)>1):
        plt.ylim([1, len(sortedSignatureList)])
    else:
        plt.ylim([0, len(sortedSignatureList)])
    #######################################################################

    ax.set_xticks(np.arange(0, len(sortedProcessiveGroupLengthList) + 1, 1))
    ax.set_yticks(np.arange(0, len(sortedSignatureList) + 1, 1))

    cmap = cm.get_cmap('Greens')

    # Set the Colour Scale
    vmin = 0
    vmax = 7
    norm = Normalize(vmin,vmax,clip=True)

    #p values
    all_p_values = []

    # Fill all_p_values
    ##########################################################################################
    ##########################################################################################
    ##########################################################################################
    for signature in sortedSignatureList:
        for processiveGroupLength in sortedProcessiveGroupLengthList:

            if signature in signature2ProcessiveGroupLength2RadiusDict:
                if processiveGroupLength in signature2ProcessiveGroupLength2RadiusDict[signature]:

                    observedValue = originalSignature2ProcessiveGroupLength2PropertiesDict[signature][processiveGroupLength]['numberofProcessiveGroups']

                    ##########################################################################
                    if (simulation2Signature2ProcessiveGroupLength2PropertiesDict is not None):

                        ########################## Fill expected values ############################
                        expectedValues = []

                        for simulation in simulation2Signature2ProcessiveGroupLength2PropertiesDict.keys():
                            # print('for debug simulation2Signature2ProcessiveGroupLength2PropertiesDict[simulation]')
                            # print(simulation2Signature2ProcessiveGroupLength2PropertiesDict[simulation])

                            if signature in simulation2Signature2ProcessiveGroupLength2PropertiesDict[simulation]:
                                if processiveGroupLength in simulation2Signature2ProcessiveGroupLength2PropertiesDict[simulation][signature]:
                                    expectedValues.append(simulation2Signature2ProcessiveGroupLength2PropertiesDict[simulation][signature][processiveGroupLength]['numberofProcessiveGroups'])
                                else:
                                    expectedValues.append(0)
                            else:
                                expectedValues.append(0)
                        ########################## Fill expected values ############################


                        if (pValueCalculation == USING_POISSON_DISTRIBUTION):
                            #################################################
                            #Using Poisson Distribution
                            mu = sum(expectedValues)/len(expectedValues)
                            prob = poisson.pmf(observedValue, mu)
                            #################################################

                        elif (pValueCalculation == USING_NULL_DISTRIBUTION):
                            #################################################
                            #Using Null Distribution
                            numberofExpectedValuesGreaterThanObservedValue = sum(expectedValue >= observedValue for expectedValue in expectedValues)
                            prob = numberofExpectedValuesGreaterThanObservedValue / len(expectedValues)
                            #################################################

                        elif (pValueCalculation == USING_GAUSSIAN_KDE):
                            #################################################
                            #Using Gaussian KDE
                            try:
                                kde = scipy.stats.gaussian_kde(expectedValues)
                                prob = kde.pdf(observedValue)
                            except np.linalg.linalg.LinAlgError:
                                print('Check this. np.linalg.linalg.LinAlgError')
                                prob = 1
                            except ValueError:
                                print('Check this. ValueError')
                                prob = 1
                            #################################################

                        all_p_values.append(prob)
                    ##########################################################################
    ##########################################################################################
    ##########################################################################################
    ##########################################################################################

    ####################################################################################
    ######  simulation2Signature2ProcessiveGroupLength2PropertiesDict is not None ######
    ####################################################################################
    if (simulation2Signature2ProcessiveGroupLength2PropertiesDict is not None):

        ##########################################################################################
        all_p_values_array = np.asarray(all_p_values)

        if (multipleTestingCorrection==FDR_BH_CORRECTION):
            #FDR BH Multiple Testing Correction
            try:
                rejected, all_FDR_BH_adjusted_p_values, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(all_p_values_array, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
            except ZeroDivisionError:
                print('ZeroDivisionError during statsmodels.stats.multitest.multipletests')
                print('for debug ZeroDivisionError, all_p_values_array:')
                print(all_p_values_array)
        elif (multipleTestingCorrection==BONFERRONI_CORRECTION):
            #Bonferroni Corrected P Values
            all_Bonferroni_corrected_p_values = all_p_values_array * numberofMultipleTests
        ##########################################################################################


        ##########################################################################################
        #Plot the circles with color
        correctedPValueIndex = 0
        for sigIndex, signature in enumerate(sortedSignatureList):
            for lengthIndex, processiveGroupLength in enumerate(sortedProcessiveGroupLengthList):

                if signature in signature2ProcessiveGroupLength2RadiusDict:
                    if processiveGroupLength in signature2ProcessiveGroupLength2RadiusDict[signature]:

                        radius = signature2ProcessiveGroupLength2RadiusDict[signature][processiveGroupLength]

                        if (multipleTestingCorrection==FDR_BH_CORRECTION):
                            correctedPValue = all_FDR_BH_adjusted_p_values[correctedPValueIndex]
                        elif (multipleTestingCorrection==BONFERRONI_CORRECTION):
                            correctedPValue = all_Bonferroni_corrected_p_values[correctedPValueIndex]

                        correctedPValueIndex +=1

                        if (correctedPValue >= 1):
                            correctedPValue = 1

                        # print('correctedPValue: %f' % (correctedPValue))

                        if (correctedPValue > 0):
                            color = -1 * math.log10(correctedPValue)
                            numberofCorrectedPValuesGreaterThanZero += 1
                        elif (correctedPValue < 0):
                            color = 0
                            numberofCorrectedPValuesLessThanZero += 1
                        else:
                            color = 7
                            numberofCorrectedPValuesEqualsToZero += 1

                        color = norm(color) * vmax
                        circle = plt.Circle((lengthIndex + 0.5, sigIndex + 0.5), radius, color=cmap(color), fill=True)
                        ax.add_artist(circle)
        ##########################################################################################

    ####################################################################################
    ######  simulation2Signature2ProcessiveGroupLength2PropertiesDict is not None ######
    ####################################################################################


    else:

    #There is  no simulation data therefore no p values
    #######################################################################################
    ######  simulation2Signature2ProcessiveGroupLength2PropertiesDict is None starts ######
    #######################################################################################

        ##########################################################################################
        #Plot the circles without color
        for sigIndex, signature in enumerate(sortedSignatureList):
            for lengthIndex, processiveGroupLength in enumerate(sortedProcessiveGroupLengthList):

                if signature in signature2ProcessiveGroupLength2RadiusDict:
                    if processiveGroupLength in signature2ProcessiveGroupLength2RadiusDict[signature]:

                        radius = signature2ProcessiveGroupLength2RadiusDict[signature][processiveGroupLength]
                        circle = plt.Circle((lengthIndex + 0.5, sigIndex + 0.5), radius,color="g", fill=True)
                        ax.add_artist(circle)
        ##########################################################################################

    #######################################################################################
    ######  simulation2Signature2ProcessiveGroupLength2PropertiesDict is None ends ########
    #######################################################################################


    ##########################################################################################
    for lengthIndex, processiveGroupLength in enumerate(processiveGroupLengthList):
        x.append(lengthIndex)
        y.append(lengthIndex)
        c.append(0.5)

    ax.set_aspect(1.0)  # make aspect ratio square

    # plot the scatter plot
    sc = plt.scatter(x, y, s=0, c=c, cmap=cmap, vmin=vmin, vmax=vmax, edgecolors='black')
    ax.set_facecolor('white')
    plt.grid(color='black')

    for edge, spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_color('black')

    xlabels = sortedProcessiveGroupLengthList
    ylabels = sortedSignatureList
    ##########################################################################################

    ################### Put the color bar if there are simulations starts ###################
    if (simulation2Signature2ProcessiveGroupLength2PropertiesDict is not None):
        cb = plt.colorbar(sc)  # this works because of the scatter
        cb.ax.set_xticklabels(cb.ax.get_xticklabels(), fontsize=20)
        cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=20)

        font = mpl.font_manager.FontProperties(size=20)
        cbax = cb.ax
        text_x = cbax.xaxis.label
        text_y = cbax.yaxis.label
        text_x.set_font_properties(font)
        text_y.set_font_properties(font)

        cb.set_label("-log10\n  (P value)", horizontalalignment='right', rotation=0, labelpad=80)
    ################### Put the color bar if there are simulations ends #####################


    ##################################################################################
    #Put a title
    plt.title('%s with %s' %(pValueCalculation,multipleTestingCorrection), fontsize=30, y=1.1)
    ##################################################################################

    ##################################################################################
    # CODE GOES HERE TO CENTER X-AXIS LABELS...
    ax.set_xticklabels([])
    mticks = ax.get_xticks()
    ax.set_xticks((mticks[:-1] + mticks[1:]) / 2, minor=True)
    ax.tick_params(axis='x', which='minor', length=0,labelsize=20)
    ax.set_xticklabels(xlabels, minor=True)

    ax.xaxis.set_ticks_position('top')

    plt.tick_params(
        axis='x',  # changes apply to the x-axis
        which='major',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False)  # labels along the bottom edge are off
    ##################################################################################

    ##################################################################################
    # CODE GOES HERE TO CENTER Y-AXIS LABELS...
    ax.set_yticklabels([])
    mticks = ax.get_yticks()
    ax.set_yticks((mticks[:-1] + mticks[1:]) / 2, minor=True)
    ax.tick_params(axis='y', which='minor', length=0,labelsize=20)
    ax.set_yticklabels(ylabels, minor=True) # fontsize

    plt.tick_params(
        axis='y',  # changes apply to the x-axis
        which='major',  # both major and minor ticks are affected
        left=False)  # labels along the bottom edge are off
    ##################################################################################


    ##################################################################################
    #create the directory if it does not exists
    os.makedirs(os.path.join(outputDir, jobname, FIGURE, ALL, PROCESSIVITY), exist_ok=True)
    filename = 'Processivity_Using_%s_%s.png' %(pValueCalculation,multipleTestingCorrection)

    figFile = os.path.join(outputDir, jobname, FIGURE, ALL, PROCESSIVITY, filename)
    fig.savefig(figFile)
    plt.close(fig)
    ##################################################################################

###################################################################



###################################################################
def processivityFigures(outputDir,jobname,numberofSimulations,multipleTestingCorrection,probabilityCalculation):

    simulation2Signature2ProcessiveGroupLength2PropertiesDict = None
    simulation2Signature2ProcessiveGroupLength2PropertiesDict = None
    originalSignature2ProcessiveGroupLength2PropertiesDict = None

    ############################################################
    if (numberofSimulations > 0):
        simulation2Signature2ProcessiveGroupLength2PropertiesDict = readSimulationBasedDictionaries(outputDir,jobname,numberofSimulations)
    ############################################################

    ############################################################
    originalSignature2ProcessiveGroupLength2PropertiesDictFilePath = os.path.join(outputDir,jobname,DATA,PROCESSIVITY,
                                                                                  'Signature2ProcessiveGroupLength2PropertiesDict.txt')
    originalSignature2ProcessiveGroupLength2PropertiesDict = readDictionary(originalSignature2ProcessiveGroupLength2PropertiesDictFilePath)
    ############################################################


    ############################################################
    if (numberofSimulations > 0):
        print('Which simulations are available?: simulation2Signature2ProcessiveGroupLength2PropertiesDict.keys()')
        print(simulation2Signature2ProcessiveGroupLength2PropertiesDict.keys())
        print('Number of simulations: len(simulation2Signature2ProcessiveGroupLength2PropertiesDict.keys())')
        print(len(simulation2Signature2ProcessiveGroupLength2PropertiesDict.keys()))
    ############################################################


    ############################################################
    # possible multipleTestingCorrection
    # BONFERRONI_CORRECTION
    # FDR_BH_CORRECTION

    # possible probabilityCalculation
    # USING_POISSON_DISTRIBUTION
    # USING_GAUSSIAN_KDE
    # USING_NULL_DISTRIBUTION

    plotRelationshipBetweenSignaturesandProcessiveGroupLengths(outputDir,jobname,originalSignature2ProcessiveGroupLength2PropertiesDict,simulation2Signature2ProcessiveGroupLength2PropertiesDict, multipleTestingCorrection,probabilityCalculation)
    ############################################################

###################################################################
