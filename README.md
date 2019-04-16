# SigProfilerTopography
SigProfilerTopography provides topography analyses for substitutions, dinucleotides and indels for each sample and all samples pooled.


**QUICK START GUIDE**

This section will guide you through the minimum steps required to run SigProfilerTopography:
1. Install the python package using pip:
```
$ pip install SigProfilerTopography
```

2. If you have installed SigProfilerTopography before, upgrade using pip:
```
$ pip install SigProfilerTopography --upgrade
```

3. Import SigProfilerTopography as follows:
```
$ python
>> from SigProfilerTopography import Topography as topography
```
3. Within a python session, you can run the topography analyses as follows:
```
>> genome= 'GRCh37'
>> jobname = 'BreastCancer560'
>> inputDir = '/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/input_test/%s' %(jobname)
>> outputDir = '/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/output_test/'
>> snpsForTopography= '%s/560_BRCA_WGS_snps_for_topography.txt' %(inputDir)
>> indelsForTopography= '%s/560_BRCA_WGS_indels_for_topography.txt' %(inputDir)

>> nucleosomeOccupancy = '/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/lib/nucleosome/wgEncodeSydhNsomeGm12878Sig.wig'
>> replicationSignal = '/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.wig'
>> replicationValley = '/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7ValleysRep1.bed'
>> replicationPeak = '/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7PkRep1.bed'

>> topography.runAnalyses(genome,snpsForTopography,indelsForTopography,outputDir,jobname,nucleosomeOccupancy,replicationSignal,replicationValley,replicationPeak)
```

4. Within a python session, you can plot the topography figures as follows:
```
>> numofSimulations = 0
>> jobname = 'BreastCancer560'
>> outputDir = '/oasis/tscc/scratch/burcak/developer/python/SigProfilerTopography/SigProfilerTopography/output_test/'

>> topography.plotFigures(outputDir,jobname,numofSimulations,'BONFERRONI_CORRECTION','USING_POISSON_DISTRIBUTION')
```

**INPUT FILE FORMAT**

This tool currently supports simple text file format. The user must provide input files with their paths.

**SAMPLE INPUT FILES**
Download sample snps and indels input data from
https://drive.google.com/open?id=1CZh_oLPmje5hcpb1x0w-64Nklf9d51ZX


**LIBRARY**

This tool uses ENCODE provided files for topography analysis such as nucleosome occupancy and replcation time.
You can also provide your local nucleosome occupancy (.wig format) and replication time (WaveSignal in .wig and Pk and Valleys in bed format) files with their paths.

**LIBRARY NUCLEOSOME OCCUPANCY**

Step1: Download nucleosome occupancy data wgEncodeSydhNsomeGm12878Sig.bigWig from
http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeSydhNsome

Step2: Convert .bigWig file into .wig

Step3: Provide wgEncodeSydhNsomeGm12878Sig.wig file under SigProfilerTopography/lib/nucleosome/wgEncodeSydhNsomeGm12878Sig.wig

**LIBRARY REPLICATION TIME**

When you install SigProfilerTopography python package, SigProfilerTopography downloads replication time data for Mcf7 cell lines under
```
SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.wig
SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7PkRep1.bed
SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7ValleysRep1.bed
```

[comment]: <Step1: Download GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.bigWig,>
[comment]: <GSM923442_hg19_wgEncodeUwRepliSeqMcf7ValleysRep1.bed.gz,>
[comment]: <GSM923442_hg19_wgEncodeUwRepliSeqMcf7PkRep1.bed.gz from>
[comment]: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM923442>
[comment]: <Step2: Convert .bed.gz into .bed>
[comment]: <Step3: Convert .bigWig file into .wig>
[comment]: <Step4: Provide these files under>
[comment]: <SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.wig>
[comment]: <SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7PkRep1.bed>
[comment]: <SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7ValleysRep1.bed>


**LIBRARY TRANSCRIPTS**

When you install SigProfilerTopography python package, SigProfilerTopography downloads transcripts for GRCh37 under
```
SigProfilerTopography/lib/transcripts/GRCh37_transcripts.txt
```
[comment]: <Transcripts>
[comment]: <Step1: Download GRCh37_transcripts.txt from>
[comment]: <https://drive.google.com/open?id=1TSyV_wA5pbPYg2g7M63m4QEp0bd7nYLB>
[comment]: <Step2: Provide GRCh37_transcripts.txt under>
[comment]: <SigProfilerTopography/lib/transcripts/GRCh37_transcripts.txt>

**LIBRARY HG19 and HG38 2bit files**

Within a python session, you can download the human genome data as follows:
```
$ python
>> from SigProfilerTopography import Topography as topography
```

For GRCh37,  you can download hg19.2bit as follows:
```
>> topography.download('GRCh37')
```

For GRCh38,  you can download hg38.2bit as follows:
```
>> topography.download('GRCh38')
```
Corresponding .2bit file will be downloaded under:
```
SigProfilerTopography/lib/ucscgenome/hg19.2bit
or
SigProfilerTopography/lib/ucscgenome/hg38.2bit
```


[comment]: <Step1: Download hg19.2bit and hg38.2bit from>
[comment]: <http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit and>
[comment]: <http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit>
[comment]: <Step2: Provide hg19.2bit and hg38.2bit under>
[comment]: <SigProfilerTopography/lib/ucscgenome/hg19.2bit>
[comment]: <SigProfilerTopography/lib/ucscgenome/hg38.2bit>


**COPYRIGHT**

This software and its documentation are copyright 2018 as a part of the SigProfiler project.
The SigProfilerTopography framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

**CONTACT INFORMATION**

Please address any queries or bug reports to Burcak Otlu at burcakotlu@eng.ucsd.edu


[comment]: <https://dillinger.io/>