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
>> genInstall.install('GRCh37', rsync=False, bash=True)
```
This will install the human 37 assembly as a reference genome. You may install as many genomes as you wish. If you have a firewall on your server, you may need to install rsync and use the rsync=True parameter. Similarly, if you do not have bash, 
use bash=False.
```
5. Import SigProfilerTopography as follows:
```
$ python
>> from SigProfilerTopography import Topography as topography
```
6. Download default nucleosome occupancy bigWig file as follows:
```
>> cell_line='GM12878'
>> topography.download_nucleosome_occupancy(cell_line)
```
7. Within a python session, you need to download corresponding 2bit file for your genome.
Command below will download hg19.2bit for 'GRCh37' and hg38.2bit for 'GRCh38'.
```
>> genome= 'GRCh37'
>> topography.download_2bit_file(genome)
```
8. Within a python session, you can run the topography analyses as follows:
You must provide the mutation types that you want to carry out topography analyses for in mutation_types_contexts list.
You must provide the corresponding probabilities files in sbs_probabilities_file_path, id_probabilities_file_path and dbs_probabilities_file_path accordingly.
For example, if you want to carry out topography analyses only for substitution (one base) and dinucleotide (two base) mutations then you must supply subs_probabilities_file_path and dinucs_probabilities_file_path with mutation_types_contexts=['96', 'DBS'].
'96' for substitutions and 'DBS' for dinucleotides and 'ID' for indels.
This call also plots topography output figures.
```
>> genome= 'GRCh37'
>> inputDir = '.../your_input_dir/'
>> outputDir = '.../your_output_dir/'
>> jobname = 'any_job_name'
>> numofSimulations = 2
>> sbs_probabilities = '.../SBS96_Mutation_Probabilities_that_comes_from_SigProfilerExtractor.txt'
>> id_probabilities = '.../ID83_Mutation_Probabilities_that_comes_from_SigProfilerExtractor.txt'
>> dbs_probabilities = '.../DBS78_Mutation_Probabilities_that_comes_from_SigProfilerExtractor.txt'
topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities_file_path=sbs_probabilities,id_probabilities_file_path=id_probabilities,dbs_probabilities_file_path=dbs_probabilities,mutation_types_contexts=['96','ID','DBS'],epigenomics=True,nucleosome=True,replication_time=True,strand_bias=True,processivity=True,sample_based=False)
```

**INPUT FILE FORMAT**

This tool currently supports formats (maf, vcf, simple text file, and ICGC) that are supported by SigProfilerMatrixGenerator. The user must provide input files under inputDir.

**LIBRARY**

This tool uses ENCODE provided files for topography analysis such as histone modifications narrow peaks, nucleosome occupancy and replcation time.
You can provide your local histone modifications and nucleosome occupancy files in .bigWig or .bigBed formats and replication time files: WaveSignal in .wig and Peaks and Valleys in bed formats with their paths.

**NUCLEOSOME OCCUPANCY**

You can provide your nucleosome occupancy data file as follows.
```

>> user_provided_nucleosome_data_file_path = '.../user_provided_nucleosome.bigWig'
or
>> user_provided_nucleosome_data_file_path = '.../user_provided_nucleosome.bigBed'

>> topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities_file_path=sbs_probabilities,dbs_probabilities_file_path=dbs_probabilities,nucleosome_file=user_provided_nucleosome_data_file_path,mutation_types_contexts=['96','DBS'],nucleosome=True)
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

>> topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities_file_path=sbs_probabilities,dbs_probabilities_file_path=dbs_probabilities,replication_time_file=user_provided_replication_time_file_path, replication_time_valley_file=user_provided_replication_time_valley_file_path, replication_time_peak_file=user_provided_replication_time_peak_file_path,mutation_types_contexts=['96','DBS'],replication_time=True,strand_bias=True)
```


**COPYRIGHT**

This software and its documentation are copyright 2018 as a part of the SigProfiler project.
The SigProfilerTopography framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

**CONTACT INFORMATION**

Please address any queries or bug reports to Burcak Otlu at burcakotlu@eng.ucsd.edu
