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

3. SigProfilerTopography requires SigProfilerMatrixGenerator and SigProfilerSimulator. Please install them as follows:
```
$ pip install SigProfilerMatrixGenerator
$ pip install SigProfilerSimulator
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

6. Within a python session, you need to download corresponding 2bit file for your genome.
Command below will download hg19.2bit for 'GRCh37' and hg38.2bit for 'GRCh38'.
```
>> genome= 'GRCh37'
>> topography.download_2bit_file(genome)
```

7. Within a python session, you can run the topography analyses as follows:
You must provide the mutation types that you want to carry out topography analyses for in mutationTypes list.
You must provide the corresponding probabilities files in subs_probabilities_file_path, indels_probabilities_file_path and dinucs_probabilities_file_path accordingly.
For example, if you want to carry out topography analyses only for substitution (one base) and dinucleotide (two base) mutations then you must supply subs_probabilities_file_path and dinucs_probabilities_file_path with mutation_types_contexts=['96', 'DBS'].
'96' for substitutions and 'DBS' for dinucleotides and 'ID' for indels.
This call also plots topography output figures.
```
>> genome= 'GRCh37'
>> inputDir = '.../from/googledrive/you/can/download/sample/input/under/matrixgenerator/'
>> outputDir = '.../as/you/wish/output/'
>> jobname = 'BreastCancer560'
>> numberofSimulations = 2
>> subs_probabilities = '.../from/googledrive/you/can/download/sample/input/under/extractor/SBS96_Mutation_Probabilities.txt'
>> dinucs_probabilities = '.../from/googledrive/you/can/download/sample/input/under/extractor/DBS78_Mutation_Probabilities.txt'
topography.runAnalyses(genome,inputDir,outputDir,jobname,numberofSimulations,subs_probabilities_file_path=subs_probabilities,dinucs_probabilities_file_path=dinucs_probabilities,mutation_types_contexts=['96','DBS'])
```

**INPUT FILE FORMAT**

This tool currently supports formats (maf, vcf, simple text file, and ICGC) that are supported by SigProfilerMatrixGenerator. The user must provide input files with their paths.

**SAMPLE INPUT FILES**
Download sample input files from
https://drive.google.com/open?id=1CZh_oLPmje5hcpb1x0w-64Nklf9d51ZX

**LIBRARY**

This tool uses ENCODE provided files for topography analysis such as nucleosome occupancy and replcation time.
You can also provide your local nucleosome occupancy (.bigWig or .bigBed format) and replication time (WaveSignal in .wig and Pk and Valleys in bed formats) files with their paths.

**NUCLEOSOME OCCUPANCY**

You can provide your nucleosome occupancy data file as follows.
```
>> user_provided_nucleosome_data_file_path = '.../user_provided_nucleosome.bigWig'
>> user_provided_nucleosome_data_file_path = '.../user_provided_nucleosome.bigBed'

>> topography.runAnalyses(genome,inputDir,outputDir,jobname,numberofSimulations,subs_probabilities_file_path=subs_probabilities,dinucs_probabilities_file_path=dinucs_probabilities,nucleosomeFilename=user_provided_nucleosome_data_file_path,mutation_types_contexts=['96','DBS'])
```

**LIBRARY REPLICATION TIME**

When you install SigProfilerTopography python package, SigProfilerTopography downloads replication time data for Mcf7 cell lines under
```
SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.wig
SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7PkRep1.bed
SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7ValleysRep1.bed
```
Unless you set different files, these installed files are used for topography replication time and replication strand bias analyses.
You can provide your replication data file as follows.
```
>> user_provided_replication_time_file_path = '.../user_provided_replication_time.wig'
>> user_provided_replication_time_valley_file_path = '.../user_provided_replication_time_valley.bed'
>> user_provided_replication_time_peak_file_path = '.../user_provided_replication_time_peak.bed'

>> topography.runAnalyses(genome,inputDir,outputDir,jobname,numberofSimulations,subs_probabilities_file_path=subs_probabilities,dinucs_probabilities_file_path=dinucs_probabilities,replicationTimeFilename=user_provided_replication_time_file_path, replicationTimeValleyFilename=user_provided_replication_time_valley_file_path, replicationTimePeakFilename=user_provided_replication_time_peak_file_path,mutation_types_contexts=['96','DBS'])
```

**LIBRARY TRANSCRIPTS**

When you install SigProfilerTopography python package, SigProfilerTopography downloads transcripts for GRCh37 under
```
SigProfilerTopography/lib/transcripts/GRCh37_transcripts.txt
```


**COPYRIGHT**

This software and its documentation are copyright 2018 as a part of the SigProfiler project.
The SigProfilerTopography framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

**CONTACT INFORMATION**

Please address any queries or bug reports to Burcak Otlu at burcakotlu@eng.ucsd.edu
