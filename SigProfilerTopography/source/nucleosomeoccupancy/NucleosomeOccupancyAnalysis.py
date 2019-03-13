# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time, strand bias and processivity.
# Copyright (C) 2018 Burcak Otlu

import sys
import os

#############################################################
current_abs_path = os.path.abspath(os.path.dirname(__file__))
print('NucleosomeOccupancyAnalysis.py current_abs_path:%s' %(current_abs_path))
#############################################################

commonsPath = os.path.join(current_abs_path, '..','commons')
sys.path.append(commonsPath)

from SigProfilerTopography.source.commons import TopographyCommons
from SigProfilerTopography.source.nucleosomeoccupancy import NucleosomeOccupancyAnalysis_SPMs_SignatureBased
from SigProfilerTopography.source.nucleosomeoccupancy import NucleosomeOccupancyAnalysis_Subs_Indels
from SigProfilerTopography.source.nucleosomeoccupancy import NucleosomeOccupancyAnalysis_Indels


##############################################################################################################
#main function
def nucleosomeOccupancyAnalysis(genome,outputDir,jobname,singlePointMutationsFilename, indelsFilename, nucleosomeFilename):
    print('########################## NucleosomeOccupancyAnalysis starts ##########################')

    print('#################### NucleosomeOccupancyAnalysis system arguments: #####################')
    print('jobname:%s singlePointMutationsFilename:%s indelsFilename:%s nucleosomeFilename:%s' %(jobname, singlePointMutationsFilename, indelsFilename, nucleosomeFilename))

    # Case 1: SPMs
    if (singlePointMutationsFilename!=TopographyCommons.NOTSET and indelsFilename==TopographyCommons.NOTSET):
        # subprocess.call(['python',os.path.join(current_abs_path,'NucleosomeOccupancyAnalysis_SPMs_SignatureBased.py'),jobname,singlePointMutationFilename,nucleosomeFilename])
        NucleosomeOccupancyAnalysis_SPMs_SignatureBased.nucleosomeOccupancyAnalysis_SPMs_SignatureBased(jobname,singlePointMutationsFilename,nucleosomeFilename)

    # Case 2: SPMs and Indels
    elif (singlePointMutationsFilename!=TopographyCommons.NOTSET and indelsFilename!=TopographyCommons.NOTSET):
        NucleosomeOccupancyAnalysis_Subs_Indels.nucleosome_occupancy_analysis_subs_indels_all_chroms_parallel(genome, outputDir, jobname, singlePointMutationsFilename, indelsFilename, nucleosomeFilename)

    # Case 3: Indels
    elif (singlePointMutationsFilename==TopographyCommons.NOTSET and indelsFilename!=TopographyCommons.NOTSET):
        # subprocess.call(['python',os.path.join(current_abs_path,'NucleosomeOccupancyAnalysis_Indels.py'),jobname,indelsFilename,nucleosomeFilename])
        NucleosomeOccupancyAnalysis_Indels.nucleosomeOccupancyAnalysis_Indels(jobname,indelsFilename,nucleosomeFilename)

    print('########################## NucleosomeOccupancyAnalysis ends ############################')
##############################################################################################################