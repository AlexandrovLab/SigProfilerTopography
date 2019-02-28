# SigProfilerTopography
SigProfilerTopography provides topography analyses for substitutions, dinucleotides and indels for each sample and all samples pooled.


**QUICK START GUIDE**

This section will guide you through the minimum steps required to run SigProfilerTopography:
1. Install the python package using pip:
```
$ pip install SigProfilerTopography
```
2. Import SigProfilerTopography as follows:
```
$ python
>> from SigProfilerTopography import SigProfilerTopography as topography
>> import os
```


3. From within a python session, you can run the topography analyses as follows:
```
$ python3
>>current_abs_path = os.path.abspath(os.path.dirname(__file__))
>>jobname = '21BreastCancer'
>>dataDir = '%s/SigProfilerTopography/input/%s' %(current_abs_path,jobname)
>>snpsForTopography= '%s/%s_snps_for_topography.txt' %(dataDir,jobname)
>>indelsForTopography= '%s/%s_indels_for_topography.txt' %(dataDir,jobname)
>>topography.runAnalyses(snpsForTopography,indelsForTopography,jobname,'wgEncodeSydhNsomeGm12878Sig.wig','GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.wig','GSM923442_hg19_wgEncodeUwRepliSeqMcf7ValleysRep1.bed','GSM923442_hg19_wgEncodeUwRepliSeqMcf7PkRep1.bed')

```


**INPUT FILE FORMAT**

This tool currently supports simple text file format. The user must provide input files with their paths.

**SAMPLE INPUT FILES**
Download sample snps and indels input data from
https://drive.google.com/open?id=1CZh_oLPmje5hcpb1x0w-64Nklf9d51ZX


**LIBRARY**

This tool uses ENCODE provided files for topography analysis such as nucleosome occupancy and replcation time.
These files have to provided under

SigProfilerTopography/lib/

By the way
SigProfilerTopography/lib/ and SigProfilerTopography/source/ must be at the same level.

**LIBRARY NUCLEOSOME OCCUPANCY**

Step1: Download nucleosome occupancy data wgEncodeSydhNsomeGm12878Sig.bigWig from
http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeSydhNsome

Step2: Convert .bigWig file into .wig

Step3: Provide wgEncodeSydhNsomeGm12878Sig.wig file under SigProfilerTopography/lib/nucleosome/wgEncodeSydhNsomeGm12878Sig.wig

**LIBRARY REPLICATION TIME**

Step1: Download GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.bigWig,
GSM923442_hg19_wgEncodeUwRepliSeqMcf7ValleysRep1.bed.gz,
GSM923442_hg19_wgEncodeUwRepliSeqMcf7PkRep1.bed.gz from
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM923442

Step2: Convert .bed.gz into .bed
Step3: Convert .bigWig file into .wig

Step4: Provide these files under
SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.wig
SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7PkRep1.bed
SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7ValleysRep1.bed


**LIBRARY TRANSCRIPTS**

Transcripts
Step1: Download GRCh37_transcripts.txt from
https://drive.google.com/open?id=1TSyV_wA5pbPYg2g7M63m4QEp0bd7nYLB

Step2: Provide GRCh37_transcripts.txt under
SigProfilerTopography/lib/transcripts/GRCh37_transcripts.txt

**LIBRARY HG19 and HG38 twobit files**

Step1: Download hg19.2bit and hg38.2bit from
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit and
http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit

Step2: Provide hg19.2bit and hg38.2bit under
SigProfilerTopography/lib/ucscgenome/hg19.2bit
SigProfilerTopography/lib/ucscgenome/hg38.2bit


**COPYRIGHT**

This software and its documentation are copyright 2018 as a part of the sigProfiler project. The SigProfilerTopography framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

**CONTACT INFORMATION**

Please address any queries or bug reports to Burcak Otlu at burcakotlu@eng.ucsd.edu