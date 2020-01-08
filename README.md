# SigProfilerTopography
SigProfilerTopography provides topography analyses for mutations such as

- SBS (Single Base Substitutions)
- ID (indels)
- DBS (Double Base Substitutions)


SigProfilerTopography carries out following analyses:

- Histone Occupancy
- Nucleosome Occupancy
- Replication Time
- Replication Strand Bias
- Transcription Strand Bias
- Processivity

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
4. Install your desired reference genome from the command line/terminal.
Currently SigProfilerTopography only supports GRCh37 (hg19), therefore install GRCh37 reference genome as follows: 
```
$ python
>> from SigProfilerMatrixGenerator import install as genInstall
>> genInstall.install('GRCh37', rsync=False, bash=True)
```
This will install the human 37 assembly as a reference genome. 
If you have a firewall on your server, you may need to install rsync and use the rsync=True parameter. Similarly, if you do not have bash,  use bash=False.

5. Import SigProfilerTopography as follows:
```
$ python
>>> from SigProfilerTopography import Topography as topography
```

6. Within a python session, you can run the topography analyses as follows:

	You must provide the mutation types in`mutation_types_contexts`  list for the ones you want to carry out topography analyses.
You must provide the corresponding probabilities files in `sbs_probabilities_file_path`, `id_probabilities_file_path` and `dbs_probabilities_file_path` accordingly.
These probabilities files must be output probabilities files of SigProfilerExtractor.

	For example, if you want to carry out topography analyses only for single base substitutions and dinucleotides then you must supply `subs_probabilities_file_path` and `dinucs_probabilities_file_path` with `mutation_types_contexts=['96', 'DBS']`.

	By the way, `96` stands for substitutions and `DBS` stand for dinucleotides and `ID` stand for indels. This call also plots topography output figures.
```
>>> genome= 'GRCh37'
>>> inputDir = '.../your_input_dir/'
>>> outputDir = '.../your_output_dir/'
>>> jobname = 'your_job_name'
>>> numofSimulations = 2
>>> sbs_probabilities = '.../SBS96_Mutation_Probabilities_that_comes_from_SigProfilerExtractor.txt'
>>> id_probabilities = '.../ID83_Mutation_Probabilities_that_comes_from_SigProfilerExtractor.txt'
>>> dbs_probabilities = '.../DBS78_Mutation_Probabilities_that_comes_from_SigProfilerExtractor.txt'
>>>topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities_file_path=sbs_probabilities,id_probabilities_file_path=id_probabilities,dbs_probabilities_file_path=dbs_probabilities,mutation_types_contexts=['96','ID','DBS'],epigenomics=True,nucleosome=True,replication_time=True,strand_bias=True,processivity=True)
```

**INPUT FILE FORMAT**

SigProfilerTopography uses formats (maf, vcf, simple text file, and ICGC) that are supported by SigProfilerMatrixGenerator. The user must provide input files under inputDir.

**LIBRARY**

SigProfilerTopography uses ENCODE provided files for topography analyses such as histone occupancy, nucleosome occupancy, replication time and replication strand bias.

You can provide your local *Histone DNA binding files* in  **bed** format and replication time files: *signal* in **wig** and *peaks* and *valleys* in **bed** format with their paths.

**HISTONE OCCUPANCY**

By default, SigProfilerTopograpy makes use of these 6 Histone DNA binding files for histone occupancy analyses.
                
1. ENCFF291WFP_H3K27me3_breast_epithelium.bed
2. ENCFF906MJM_H3K36me3_breast_epithelium.bed
3. ENCFF065FJK_H3K9me3_breast_epithelium.bed
4. ENCFF154XFN_H3K27ac_breast_epithelium.bed
5. ENCFF336DDM_H3K4me1_breast_epithelium.bed
6. ENCFF065TIH_H3K4me3_breast_epithelium.bed
                


If you prefer to use other Histone DNA binding files, you have to provide them with their full path and then include them in list `epigenomics_files_list`.

You can have as many Histone DNA binding files as you want.

`histone_dna_binding_file1=.../full_path_to/file1.bed`
`histone_dna_binding_file2=.../full_path_to/file2.bed`
`epigenomics_files_list=[histone_dna_binding_file1,histone_dna_binding_file2]`

Then you need to provide `epigenomics_files_list` in the `runAnalyses` call as follows:

`>>>topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities_file_path=sbs_probabilities,id_probabilities_file_path=id_probabilities,dbs_probabilities_file_path=dbs_probabilities,mutation_types_contexts=['96','ID','DBS'],epigenomics=True,nucleosome=True,replication_time=True,strand_bias=True,processivity=True,epigenomics_files=epigenomics_files_list)`

**NUCLEOSOME OCCUPANCY**

We obtained micrococcal nuclease sequencing (MNase-seq) data for GM12878 and K562 cell lines from ENCODE.

By default, SigProfilerTopography makes use of MNase-seq of GM12878 cell line for nucleosome occupancy analysis.

If you want to run SigProfilerTopography using  MNase-seq of K562 cell line, you have to include
`nucleosome_file='wgEncodeSydhNsomeK562Sig.bigWig'`  in the `runAnalyses` call as follows:

`>>>topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities_file_path=sbs_probabilities,id_probabilities_file_path=id_probabilities,dbs_probabilities_file_path=dbs_probabilities,mutation_types_contexts=['96','ID','DBS'],epigenomics=True,nucleosome=True,replication_time=True,strand_bias=True,processivity=True,nucleosome_file='wgEncodeSydhNsomeK562Sig.bigWig')`


**REPLICATION TIME** and **REPLICATION STRAND BIAS**

When you install SigProfilerTopography, SigProfilerTopography downloads replication time data for Mcf7 cell lines under
```
SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.wig
SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7PkRep1.bed
SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7ValleysRep1.bed
```
Unless you set different files, these installed files are used for replication time and replication strand bias analyses.

You can provide your replication data file as follows:
```
>>> user_provided_replication_time_file_path = '.../user_provided_replication_time.wig'
>>> user_provided_replication_time_valley_file_path = '.../user_provided_replication_time_valley.bed'
>>> user_provided_replication_time_peak_file_path = '.../user_provided_replication_time_peak.bed'
>>> topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities_file_path=sbs_probabilities,dbs_probabilities_file_path=dbs_probabilities,replication_time_file=user_provided_replication_time_file_path, replication_time_valley_file=user_provided_replication_time_valley_file_path, replication_time_peak_file=user_provided_replication_time_peak_file_path,mutation_types_contexts=['96','DBS'],replication_time=True,strand_bias=True)
```
**SIGPROFILERTOPOGRAPHY PARAMETERS**

`num_of_sbs_required` `num_of_id_required` `num_of_dbs_required` 
`average probabilty`

                

+ By default, we require at least 5000, 1000, 200 mutations with average probability of 0.9 for a SBS, ID, DBS signature, respectively, to be considered and analysed.

    * `average_probability=0.9`
    * `num_of_sbs_required=5000`
    * `num_of_id_required=1000`
	* `num_of_dbs_required=200`

+ However, as an example these parameters can be relaxed in the `runAnalyses` call as follows:

`>>>topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities_file_path=sbs_probabilities,id_probabilities_file_path=id_probabilities,dbs_probabilities_file_path=dbs_probabilities,mutation_types_contexts=['96','ID','DBS'],epigenomics=True,nucleosome=True,replication_time=True,strand_bias=True,processivity=True,average_probability=0.75,num_of_sbs_required=2000,num_of_id_required=1000,num_of_dbs_required=200)`


**COPYRIGHT**

This software and its documentation are copyright 2020 as a part of the SigProfiler project.
The SigProfilerTopography framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

**CONTACT INFORMATION**

Please address any queries or bug reports to Burcak Otlu at burcakotlu@eng.ucsd.edu
