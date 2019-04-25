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

3. SigProfilerTopography requires SigProf ilerMatrixGenerator. Install branch of this SigProfilerMatrixGenerator.
```
$ git clone --single-branch --branch Development https://github.com/AlexandrovLab/SigProfilerMatrixGenerator.git
$ cd SigProfilerMatrixGenerator
$ pip install .
```
4. Install your desired reference genome from the command line/terminal as follows (available reference genomes are: GRCh37, GRCh38, mm9, and mm10):
```
$ python
>> from SigProfilerMatrixGenerator import install as genInstall
>> genInstall.install('GRCh37')
```

5. Import SigProfilerTopography as follows:
```
$ python
>> from SigProfilerTopography import Topography as topography
```

6. Within a python session, you can run the topography analyses as follows: 
This call also plots topography output figures.
```
>> genome= 'GRCh37'
>> inputDir = '.../from/googledrive/you/can/download/sample/input/under/matrixgenerator/'
>> outputDir = '.../as/you/wish/output/'
>> jobname = 'BreastCancer560'
>> numofSimulations = 0
>> subs_probabilities_file_path = '.../from/googledrive/you/can/download/sample/input/under/extractor/SBS96_Mutation_Probabilities.txt'
>> indels_probabilities_file_path = '.../from/googledrive/you/can/download/sample/input/under/extractor/ID83_Mutation_Probabilities.txt'
>> dinucs_probabilities_file_path = '.../from/googledrive/you/can/download/sample/input/under/extractor/DBS78_Mutation_Probabilities.txt'
>> topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,subs_probabilities_file_path,indels_probabilities_file_path,dinucs_probabilities_file_path)
```

**INPUT FILE FORMAT**

This tool currently supports formats (maf, vcf, simple text file, and ICGC) that are supported by SigProfilerMatrixGenerator. The user must provide input files with their paths.

**SAMPLE INPUT FILES**
Download sample input files from
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
Unless you set different files, these installed files are used for topography replication time and replication strand bias analyses.
```
**LIBRARY TRANSCRIPTS**

When you install SigProfilerTopography python package, SigProfilerTopography downloads transcripts for GRCh37 under
```
SigProfilerTopography/lib/transcripts/GRCh37_transcripts.txt
For GRCh37,  you can download hg19.2bit as follows:
```
**LIBRARY HG19 and HG38 2bit files**

Within a python session, you can download the human genome data as follows:
```

For GRCh37,  you can download hg19.2bit as follows:
```
>> $ python
>> from SigProfilerTopography import Topography as topography
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

**COPYRIGHT**

This software and its documentation are copyright 2018 as a part of the SigProfiler project.
The SigProfilerTopography framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

**CONTACT INFORMATION**

Please address any queries or bug reports to Burcak Otlu at burcakotlu@eng.ucsd.edu
