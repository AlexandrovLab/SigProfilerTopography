# This source code file is a part of SigProfilerTopography
# SigProfilerTopography is a tool included as part of the SigProfiler
# computational framework for comprehensive analysis of mutational
# signatures from next-generation sequencing of cancer genomes.
# SigProfilerTopography provides the downstream data analysis of
# mutations and extracted mutational signatures w.r.t.
# nucleosome occupancy, replication time and strand bias.
# Copyright (C) 2018 Burcak Otlu


from SigProfilerTopography.source.commons.TopographyCommons import readIndelsAndWriteChrBasedParallel

###############################################################
def partitionIndelsData(genome,outputDir,jobname,indelsFilename):
    # Parallel
    readIndelsAndWriteChrBasedParallel(genome,outputDir,jobname,indelsFilename)
###############################################################
