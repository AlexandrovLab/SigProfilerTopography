[![License](https://img.shields.io/badge/License-BSD\%202--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause)

![schematic](SigProfilerTopography.png)

# SigProfilerTopography
SigProfilerTopography provides topography analyses for mutations such as

- Single Base Substitutions (SBS)
- Doublet Base Substitutions (DBS)
- Small insertions and deletions, indels (ID)

and carries out following analyses:

- Epigenomics Occupancy (e.g.: Histone Modifications, Transcription Factors, Open Chromatin Regions)
- Nucleosome Occupancy
- Replication Timing
- Replication Strand Asymmetry
- Transcription Strand Asymmetry
- Genic versus Intergenic Regions
- Strand-coordinated Mutagenesis

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

3. Install default nucleosome occupancy (Mnase-seq), open chromatin regions (ATAC-seq) and replication timing (Repli-seq) files for the genome of interest as follows:

	For GRCh37, Mnase-seq of K562, ATAC-seq of breast epithelium tissue and Repli-seq of MCF7 breast cancer tissue are used as default files:
	```
	$ python
	>> from SigProfilerTopography import Topography as topography
	>> topography.install_nucleosome('GRCh37')
	>> topography.install_atac_seq('GRCh37')
	>> topography.install_repli_seq('GRCh37')
	```

	For GRCh38, Mnase-seq of K562, ATAC-seq of left lung tissue and Repli-seq of IMR90 lung tissue are used as default files:
	```
	$ python
	>> from SigProfilerTopography import Topography as topography
	>> topography.install_nucleosome('GRCh38')
	>> topography.install_atac_seq('GRCh38')
	>> topography.install_repli_seq('GRCh38')
	```

	For mm10, Mnase-seq of mouse embryonic stem cells (GSM1004653), ATAC-seq of embryonic_facial_prominence, Repli-Seq of endoterm are used as default files:
	```
	$ python
	>> from SigProfilerTopography import Topography as topography
	>> topography.install_nucleosome('mm10')
	>> topography.install_atac_seq('mm10')
	>> topography.install_repli_seq('mm10')
	```

4. For GRCh37 and GRCh38, SigProfilerTopography provides Repli-seq files of the biosamples listed in the table below:


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
	| HEK293  | human  | embryonic kidney | epithelial | Unknown |

	To carry out replication timing analysis for e.g.: GRCh38 and biosample of HEK293, first install Repli-seq files as follows:

	```
	$ python
	>> from SigProfilerTopography import Topography as topography
	>> topography.install_repli_seq('GRCh38', 'HEK293')
	```

5. SigProfilerTopography requires SigProfilerMatrixGenerator and SigProfilerSimulator. Please install them as follows:
	```
	$ pip install SigProfilerMatrixGenerator
	$ pip install SigProfilerSimulator
	```
6. Install your desired reference genome from the command line/terminal. SigProfilerTopography supports GRCh37 (hg19), GRCh38 (hg38), mm9 and mm10. 

	For GRCh37 (hg19) install GRCh37 reference genome as follows:

	````
	$ python
	>> from SigProfilerMatrixGenerator import install as genInstall
	>> genInstall.install('GRCh37', rsync=False, bash=True)
	````
	This will install the human GRCh37 assembly as a reference genome. 
	
	If you have a firewall on your server, you may need to install rsync and use the rsync=True parameter. Similarly, if you do not have bash, use bash=False.


7. You can download 21 sample vcf files (in GRCh37) and corresponding probability files using SigProfilerTopography.

	```
	$ python
	>> from SigProfilerTopography import Topography as topography
	>> topography.install_sample_vcf_files()
	>> topography.install_sample_probability_files()

	# This will download 21 sample vcf files under current_working_directory/sample_vcfs/.
	# This will download sample probabilities files under current_working_directory/sample_probabilities/.
	```

8. Within a python session, you can run the topography analyses for these 21 sample vcf files as follows:

	+ You must provide the corresponding probabilities files in `sbs_probabilities`, `dbs_probabilities` and `id_probabilities` accordingly.
		* For example, if you want to carry out topography analyses only for single base substitutions and doublet base substitutions then you must only set `sbs_probabilities` and `dbs_probabilities` parameters.

	+ These probabilities files must be output probabilities files coming from SigProfilerExtractor.

	+ The run below will also plot the resulting figures.

		````
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
		````

**INPUT FILE FORMAT**

SigProfilerTopography uses formats (maf, vcf, simple text file, and ICGC) that are supported by SigProfilerMatrixGenerator. The user must provide input files under inputDir.

**LIBRARY**

SigProfilerTopography uses ENCODE provided files for topography analyses such as histone occupancy, nucleosome occupancy, replication timing and replication strand asymmetry.


+ **EPIGENOMICS OCCUPANCY**

	For GRCh37, by default SigProfilerTopograpy makes use of 7 ChIP-seq files of histone post-translational modification sites and CTCF binding sites.

		ENCFF291WFP_H3K27me3_breast_epithelium.bed
		ENCFF906MJM_H3K36me3_breast_epithelium.bed
		ENCFF065FJK_H3K9me3_breast_epithelium.bed
		ENCFF154XFN_H3K27ac_breast_epithelium.bed
		ENCFF336DDM_H3K4me1_breast_epithelium.bed
		ENCFF065TIH_H3K4me3_breast_epithelium.bed
        ENCFF782GCQ_CTCF_breast_epithelium.bed


	For GRCh38, by default SigProfilerTopograpy makes use of 7 ChIP-seq files of histone post-translational modification sites and CTCF binding sites.

		ENCFF367LMT_H3K27me3_lower-lobe-of-left-lung.bed
		ENCFF078AMJ_H3K36me3_upper-lobe-of-left-lung.bed
		ENCFF185EKW_H3K9me3_left-lung.bed
		ENCFF846LAZ_H3K27ac_left-lung.bed
		ENCFF908EUN__H3K4me1_lung.bed
		ENCFF069GQK_H3K4me3_lung.bed
        ENCFF061UVF_CTCF_upper-lobe-of-left-lung.bed


	For mm10, by default SigProfilerTopograpy makes use of 5 ChIP-seq files of histone post-translational modification sites, CTCF and POLR2A binding sites.

		ENCFF114VLZ_H3K27ac_embryonic_fibroblast.bed
		ENCFF993SRY_H3K4me1_embryonic_fibroblast.bed
		ENCFF912DNP_H3K4me3_embryonic_fibroblast.bed
        ENCFF611HDQ_CTCF_embryonic_fibroblast.bed
        ENCFF152DUV_POLR2A_embryonic_fibroblast.bed

+ **NUCLEOSOME OCCUPANCY**

	+ We obtained micrococcal nuclease sequencing (MNase-seq) data for GM12878 and K562 cell lines from ENCODE for GRCh37 and GRCh38 .

	+ By default, SigProfilerTopography makes use of MNase-seq of K562 cell line for nucleosome occupancy analysis.

	+ If you want to run SigProfilerTopography using MNase-seq of GM12878 cell line, you need to first install nucleosome occupancy data for the genome of interest e.g.: GRCh38 as follows:

		````
		$ python
		>> from SigProfilerTopography import Topography as topography
		>> topography.install_nucleosome('GRCh38', 'GM12878')
		````

	+ SigProfilerTopography downloads chrom based signal arrays from **ftp://alexandrovlab-ftp.ucsd.edu/**  under *.../SigProfilerTopography/lib/nucleosome/chrbased/*  for the `nucleosome_biosample` of interest which requires ~6 GB of storage.

	+ Then you have to include `nucleosome_biosample='GM12878'` in the `runAnalyses` call as follows:

		````
		>>>topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities=sbs_probabilities_file_path,dbs_probabilities=dbs_probabilities_file_path,id_probabilities=id_probabilities_file_path,epigenomics=True,nucleosome=True,replication_time=True,strand_bias=True,processivity=True,nucleosome_biosample='GM12878')	
		````


+ **REPLICATION TIMING** and **REPLICATION STRAND ASYMMETRY**

	+ By default, SigProfilerTopography carries out replication timing and replication strand asymmetry analyses using Repli-seq of MCF7 and IMR90 cell line for GRCh37 and GRCh38, respectively.

	+ If you want to run SigProfilerTopography with Repli-seq of HELAS3 cell line, you need to first install replication timing data for the genome of interest e.g.: GRCh38 as follows:
		````
		$ python
		>> from SigProfilerTopography import Topography as topography
		>> topography.install_repli_seq('GRCh38', 'HELAS3')
		````
	
	+ Then you have to include `replication_time_biosample='HELAS3'` in the `runAnalyses` call as follows:
		            
		````
		>>>topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities=sbs_probabilities_file_path,dbs_probabilities=dbs_probabilities_file_path,id_probabilities=id_probabilities_file_path,epigenomics=True,nucleosome=True,replication_time=True,strand_bias=True,processivity=True,replication_time_biosample='HELAS3')
		````
		                    
	+ SigProfilerTopography downloads replication timing files from **ftp://alexandrovlab-ftp.ucsd.edu/**  under *.../SigProfilerTopography/lib/replication/*  for the `replication_time_biosample` of interest which requires ~20-100 MB of storage.

                    

**USER PROVIDED LIBRARY FILES**
+ **EPIGENOMICS OCCUPANCY**

	+ If you prefer to use your own library files, then you have to provide them with their full path and then include them in `epigenomics_files` parameter.
	
	+ You can have as many files as you want.
		````
		histone_dna_binding_file1=.../path_to/file1.bed
		histone_dna_binding_file2=.../path_to/file2.bed
		ctcf_dna_binding_file1=.../path_to/file3.bed
		ctcf_dna_binding_file2=.../path_to/file4.bed
		atac_seq_open_chromatin_region_file1=.../path_to/file5.bed
		atac_seq_open_chromatin_region_file2=.../path_to/file6.bed

		epigenomics_files_list=[histone_dna_binding_file1, histone_dna_binding_file2, ctcf_dna_binding_file1, ctcf_dna_binding_file2, atac_seq_open_chromatin_region_file1, atac_seq_open_chromatin_region_file2]
		````

	+ Then you need to provide `epigenomics_files` in the `runAnalyses` call as follows:
			`>>>topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities=sbs_probabilities_file_path,dbs_probabilities=dbs_probabilities_file_path,id_probabilities=id_probabilities_file_path,epigenomics=True,nucleosome=True,replication_time=True,strand_bias=True,processivity=True,epigenomics_files=epigenomics_files_list)`

	
+ **NUCLEOSOME OCCUPANCY**
	+ SigProfilerTopography enables user to provide `nucleosome_file` to be used in nucleosome occupancy analysis.

	+ You can provide your nucleosome data file as follows:
		````
		user_nucleosome_file = '.../path_to/nucleosome.wig'
		````

	+ Then you need to set `nucleosome_file` in the `runAnalyses` call as follows:
	
		````
		>>>topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities=sbs_probabilities_file_path,dbs_probabilities=dbs_probabilities_file_path,nucleosome_file=user_nucleosome_file,nucleosome=True)
		````

+ **REPLICATION TIMING** and **REPLICATION STRAND ASYMMETRY**
	+ You can provide your replication timing files as follows:
		````
		user_replication_time_signal_file = .../path_to/replication_time_signal.wig'
		user_replication_time_valley_file = '.../path_to/replication_time_valley.bed'
		user_replication_time_peak_file = '.../path_to/replication_time_peak.bed'
		````

	+ Then you need to set `replication_time_signal_file`, `replication_time_valley_file`, and `replication_time_peak_file`parameters in the `runAnalyses` call as follows:
	
		```
		>>>topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities=sbs_probabilities_file_path,dbs_probabilities=dbs_probabilities_file_path,replication_time_signal_file=user_replication_time_signal_file,replication_time_valley_file=user_replication_time_valley_file,replication_time_peak_file=user_replication_time_peak_file,replication_time=True,strand_bias=True)
		```

**SIGPROFILERTOPOGRAPHY PARAMETERS**

                  
| Parameter | Default Value  |
| ------------- | ------------- |
| average probability | 0.75  |
| num_of_sbs_required | 2000  |
| num_of_dbs_required | 200  |
| num_of_id_required | 1000  |
                    


+ By default, we require at least 2000, 200, 1000 mutations with average probability of 0.75 for a SBS, DBS, ID signature, respectively, to be considered and analysed.

+ However, as an example these parameters can be relaxed in the `runAnalyses` call as follows:

	````
	>>>topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities=sbs_probabilities_file_path,dbs_probabilities=dbs_probabilities_file_path,id_probabilities=id_probabilities_file_path,epigenomics=True,nucleosome=True,replication_time=True,strand_bias=True,processivity=True,average_probability=0.7,num_of_sbs_required=1500,num_of_dbs_required=200,num_of_id_required=1000)
	````

+ You can run topography analyses with `discreet_mode = False` and allow each mutation with signature specific probability >= 0.5 (`default_cutoff` parameter) considered in the analyses in the `runAnalyses` call as follows:

	````
	>>>topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities=sbs_probabilities_file_path,dbs_probabilities=dbs_probabilities_file_path,id_probabilities=id_probabilities_file_path,epigenomics=True,nucleosome=True,replication_time=True,strand_bias=True,processivity=True,disreet_mode=False)	
	````
	
+ By default, SigProfilerTopography expects `sbs_probabilities` in any of these contexts: SBS_6, SBS_24, SBS_96, SBS_192, SBS_288, SBS_384, SBS_1536 or SBS_6144; `dbs_probabilities`,  in DBS_78 and  `id_probabilities` in ID_83 contexts, respectively for topography analyses.

+ SigProfilerTopography automatically detects the mutation type context in the 2nd column of `sbs_probabilities` . However, you can also set `sigprofiler_extractor_sbs_mutation_context` e.g.: as SBS_288 in the `runAnalyses` call as follows:
		
	````
	>>>topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,sbs_probabilities=sbs_probabilities_file_path,dbs_probabilities=dbs_probabilities_file_path,id_probabilities=id_probabilities_file_path,epigenomics=True,nucleosome=True,replication_time=True,strand_bias=True,processivity=True,sigprofiler_extractor_sbs_mutation_context='288')
	````

+ If you don't have `'sbs_probabilities'`, `'dbs_probabilities'` and `'id_probabilities'`, you can still achieve topography analyses in aggregated mutations mode as follows: 
	
	
	```
	>>>topography.runAnalyses(genome,inputDir,outputDir,jobname,numofSimulations,epigenomics=True,nucleosome=True,replication_time=True,strand_bias=True,processivity=True)
	```



**COPYRIGHT**

This software and its documentation are copyright 2018-2022 as a part of the SigProfiler project.
The SigProfilerTopography framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

**CONTACT INFORMATION**

Please address any queries or bug reports to Burcak Otlu at burcakotlu@eng.ucsd.edu

