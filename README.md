[![License](https://img.shields.io/badge/License-BSD\%202--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause)

# SigProfilerTopography
SigProfilerTopography provides topography analyses for mutations such as

- Single Base Substitutions (SBS)
- Double Base Substitutions (DBS)
- Insertions and deletions, indels (ID)


and carries out following analyses:

- Epigenomics Occupancy (e.g.: Histone Modifications, Transcription Factors, ATAC-seq for open chromatin regions)
- Nucleosome Occupancy
- Replication Time
- Replication Strand Bias
- Transcription Strand Bias
- Processivity

**PREREQUISITES**

SigProfilerTopography is written in PYTHON, and it also requires the following software with the given version (or newer):

- WGET version 1.9 or RSYNC if you have a firewall

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

3. Install default nucleosome occupancy (Mnase-seq) and open chromatin regions (ATAC-seq) files:
```
$ python
>> from SigProfilerTopography import Topography as topography
>> topography.install_default_nucleosome('GRCh37')
>> topography.install_default_atac_seq('GRCh37')
```

4. SigProfilerTopography requires SigProfilerMatrixGenerator and SigProfilerSimulator. Please install them as follows:
```
$ pip install SigProfilerMatrixGenerator
$ pip install SigProfilerSimulator
```
5. Install your desired reference genome from the command line/terminal.
SigProfilerTopography supports GRCh37 (hg19), GRCh38 (hg38), mm9 and mm10.
For GRCh37 (hg19) install GRCh37 reference genome as follows:
```
$ python
>> from SigProfilerMatrixGenerator import install as genInstall
>> genInstall.install('GRCh37', rsync=False, bash=True)
```
This will install the human 37 assembly as a reference genome. 
If you have a firewall on your server, you may need to install rsync and use the rsync=True parameter. Similarly, if you do not have bash,  use bash=False.

6. You can download 21 sample vcf files (in GRCh37) and corresponding probability files using SigProfilerTopography.
```
$ python
>> from SigProfilerTopography import Topography as topography
>> topography.install_sample_vcf_files()
>> topography.install_sample_probability_files()

# This will download 21 sample vcf files under current_working_directory/sample_vcfs/.
# This will download sample probabilities files under current_working_directory/sample_probabilities/.
```

7. Within a python session, you can run the topography analyses for these 21 sample vcf files as follows:

	+ You must provide the corresponding probabilities files in `sbs_probabilities`, `id_probabilities` and `dbs_probabilities` accordingly.
		* For example, if you want to carry out topography analyses only for single base substitutions and dinucleotides then you must supply only `sbs_probabilities` and `dbs_probabilities`

	+ These probabilities files must be output probabilities files coming from SigProfilerExtractor.

	+ The run below will also plot the resulting figures.

```
from SigProfilerTopography import Topography as topography

def main_function():
 genome= 'GRCh37'
 inputDir = '/path/to/sample_vcfs/'
 outputDir = '/path/to/output_dir/'
 jobname = 'sample_21_vcfs'
 numofSimulations = 10
 sbs_probabilities_file_path = '/path/to/sample_probabilities/COSMIC_SBS96_Decomposed_Mutation_Probabilities.txt'
 dbs_probabilities_file_path = '/path/to/sample_probabilities/COSMIC_DBS78_Decomposed_Mutation_Probabilities.txt'
 topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities=sbs_probabilities_file_path,dbs_probabilities=dbs_probabilities_file_path,epigenomics=True,nucleosome=True,replication_time=True,strand_bias=True,processivity=True)

if __name__ == "__main__":
 main_function()
```

**INPUT FILE FORMAT**

SigProfilerTopography uses formats (maf, vcf, simple text file, and ICGC) that are supported by SigProfilerMatrixGenerator. The user must provide input files under inputDir.

**LIBRARY**

SigProfilerTopography uses ENCODE provided files for topography analyses such as histone occupancy, nucleosome occupancy, replication time and replication strand bias.


+ **HISTONE OCCUPANCY**

	By default, SigProfilerTopograpy makes use of these 6 Histone DNA binding files for histone occupancy analyses.

		ENCFF291WFP_H3K27me3_breast_epithelium.bed
		ENCFF906MJM_H3K36me3_breast_epithelium.bed
		ENCFF065FJK_H3K9me3_breast_epithelium.bed
		ENCFF154XFN_H3K27ac_breast_epithelium.bed
		ENCFF336DDM_H3K4me1_breast_epithelium.bed
		ENCFF065TIH_H3K4me3_breast_epithelium.bed


+ **NUCLEOSOME OCCUPANCY**

	+ We obtained micrococcal nuclease sequencing (MNase-seq) data for GM12878 and K562 cell lines from ENCODE.

	+ By default, SigProfilerTopography makes use of MNase-seq of K562 cell line for nucleosome occupancy analysis.

	+ If you want to run SigProfilerTopography using  MNase-seq of GM12878 cell line, you have to include `nucleosome_biosample='GM12878'`  in the `runAnalyses` call as follows:
		
		`>>>topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities=sbs_probabilities_file_path,dbs_probabilities=dbs_probabilities_file_path,id_probabilities=id_probabilities_file_path,epigenomics=True,nucleosome=True,replication_time=True,strand_bias=True,processivity=True,nucleosome_biosample='GM12878')`

	+ SigProfilerTopography downloads offline prepared chrom based signal arrays from **ftp://alexandrovlab-ftp.ucsd.edu/**  under *.../SigProfilerTopography/lib/nucleosome/chrbased/*  for the `nucleosome_biosample` you have set  which requires ~6 GB of storage.

+ **REPLICATION TIME** and **REPLICATION STRAND BIAS**

	+ By default, SigProfilerTopography carries out replication time and replication strand bias analyses using Repli-seq of MCF7 cell line.

	+ SigProfilerTopography provides replication time files for the biosamples listed in the table below:
	
	                    
		| Biosample | Organism  | Tissue | Cell Type | Disease |
		| ------------- | ------------- | ------------- | ------------- | ------------- |
		| MCF7 | human  | breast  | mammary | cancer |
		| HEPG2 | human  | liver  |  | cancer |
		| HELAS3 | human  | cervix  |  | cancer |
		| SKNSH | human  | brain  |  | cancer |
		| K562 | human  | bone marrow  |  | cancer |
		| IMR90 | human  | lung  | fibroblast | normal |
		| NHEK | human  | skin  | keratinocyte | normal |
		| BJ | human  | skin  | fibroblast | normal |
		| HUVEC | human  | umbilical vein  | endothelial | normal |
		| BG02ES | human  |   | embryonic stem cell | None reported  |
		| GM12878 | human  | blood  | B-Lymphocyte | normal |
		| GM06990 | human  | blood  | B-Lymphocyte | Unknown |
		| GM12801 | human  | blood | B-Lymphocyte | Unknown  |
		| GM12812 | human  | blood | B-Lymphocyte | Unknown |
		| GM12813 | human  | blood | B-Lymphocyte | Unknown |
                    

	
	+ SigProfilerTopography downloads replication time files from **ftp://alexandrovlab-ftp.ucsd.edu/**  under *.../SigProfilerTopography/lib/replication/*  for the `replication_time_biosample` you have set which requires ~25 MB of storage.

	+ If you want to run SigProfilerTopography using  Repli-seq of any available biosamples e.g.: NHEK,  then you have to include `replication_time_biosample='NHEK'`  in the `runAnalyses` call as follows:

		`>>>topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities=sbs_probabilities_file_path,dbs_probabilities=dbs_probabilities_file_path,id_probabilities=id_probabilities_file_path,epigenomics=True,nucleosome=True,replication_time=True,strand_bias=True,processivity=True,replication_time_biosample='NHEK')`
                    

**USER PROVIDED LIBRARY FILES**
+ **HISTONE OCCUPANCY**

	+ If you prefer to use other Histone DNA binding files, you have to provide them with their full path and then include them in `epigenomics_files`

	+ You can have as many Histone DNA binding files as you want.

		`histone_dna_binding_file1=.../path_to/file1.bed`
	
		`histone_dna_binding_file2=.../path_to/file2.bed`
	
		`epigenomics_files_list=[histone_dna_binding_file1,histone_dna_binding_file2]`

	+ Then you need to provide `epigenomics_files` in the `runAnalyses` call as follows:

		`>>>topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities=sbs_probabilities_file_path,dbs_probabilities=dbs_probabilities_file_path,id_probabilities=id_probabilities_file_path,epigenomics=True,nucleosome=True,replication_time=True,strand_bias=True,processivity=True,epigenomics_files=epigenomics_files_list)`

	
+ **NUCLEOSOME OCCUPANCY**
	+ SigProfilerTopography enables user to provide `nucleosome_file` to be used in nucleosome occupancy analysis.

	+ You can provide your nucleosome data file as follows:
	
		` user_nucleosome_file = '.../path_to/nucleosome.wig'`

	+ Then you need to set `nucleosome_file` in the `runAnalyses` call as follows:
	
		`>>>topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities=sbs_probabilities_file_path,dbs_probabilities=dbs_probabilities_file_path, nucleosome_file=user_nucleosome_file,nucleosome=True)`

+ **REPLICATION TIME** and **REPLICATION STRAND BIAS**
	+ You can provide your replication data files as follows:
	
		` user_replication_time_signal_file = '.../path_to/replication_time_signal.wig'`
	
		`user_replication_time_valley_file = '.../path_to/replication_time_valley.bed'`
	
		`user_replication_time_peak_file = '.../path_to/replication_time_peak.bed'`

	+ Then you need to set `replication_time_signal_file`, `replication_time_valley_file`, and `replication_time_peak_file`, in the `runAnalyses` call as follows:
	
		`>>>topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities=sbs_probabilities_file_path,dbs_probabilities=dbs_probabilities_file_path,replication_time_signal_file=user_replication_time_signal_file, replication_time_valley_file=user_replication_time_valley_file, replication_time_peak_file=user_replication_time_peak_file,replication_time=True,strand_bias=True)`

**SIGPROFILERTOPOGRAPHY PARAMETERS**

                    
| Parameter | Default Value  |
| ------------- | ------------- |
| average probability | 0.75  |
| num_of_sbs_required | 2000  |
| num_of_dbs_required | 200  |
| num_of_id_required | 1000  | |
                    


+ By default, we require at least 2000, 200, 1000 mutations with average probability of 0.75 for a SBS, DBS, ID signature, respectively, to be considered and analysed.

+ However, as an example these parameters can be relaxed in the `runAnalyses` call as follows:

	`>>>topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities=sbs_probabilities_file_path,dbs_probabilities=dbs_probabilities_file_path,id_probabilities=id_probabilities_file_path,epigenomics=True,nucleosome=True,replication_time=True,strand_bias=True,processivity=True,average_probability=0.7,num_of_sbs_required=1500,num_of_dbs_required=200,num_of_id_required=1000)`
	
+ By default, SigProfilerTopography appends `'96'`, `'DBS'`, and `'ID'` to  `mutation_types_contexts` if `sbs_probabilities`, `dbs_probabilities`,  and  `id_probabilities` are set respectively for topography analyses.

+ If you want to analyse SBS signatures using `'288'` context instead of `'96'` context then you have to set mutation_types_contexts=`['288', 'DBS', 'ID']` in the `runAnalyses` call as follows:
	
	
`>>>topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities=sbs_probabilities_file_path,dbs_probabilities=dbs_probabilities_file_path,id_probabilities=id_probabilities_file_path,epigenomics=True,nucleosome=True,replication_time=True,strand_bias=True,processivity=True,mutation_types_contexts=['288', 'DBS', 'ID'])`

+ If you want to analyse only SBS signatures using `'288'` context then you have to set `mutation_types_contexts=['288']` in the `runAnalyses` call as follows: 
	
	
`>>>topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities=sbs_probabilities_file_path,epigenomics=True,nucleosome=True,replication_time=True,strand_bias=True,processivity=True,mutation_types_contexts=['288'])`

+ If you don't have `'sbs_probabilities'`, `'dbs_probabilities'` and `'id_probabilities'`, you can still achieve aggregated analyses by setting mutation_types_contexts=`['96', 'DBS', 'ID']` in the `runAnalyses` call as follows: 
	
	
`>>>topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,epigenomics=True,nucleosome=True,replication_time=True,strand_bias=True,processivity=True,mutation_types_contexts=['96','DBS','ID'])`



**COPYRIGHT**

This software and its documentation are copyright 2018-2020 as a part of the SigProfiler project.
The SigProfilerTopography framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

**CONTACT INFORMATION**

Please address any queries or bug reports to Burcak Otlu at burcakotlu@eng.ucsd.edu
