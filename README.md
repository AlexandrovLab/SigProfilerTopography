[![License](https://img.shields.io/badge/License-BSD\%202--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause)
[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://osf.io/5unby/wiki/home/)
 [![Build Status](https://app.travis-ci.com/AlexandrovLab/SigProfilerTopography.svg?branch=master)](https://app.travis-ci.com/AlexandrovLab/SigProfilerTopography)
[![Uptime Robot status](https://img.shields.io/uptimerobot/status/m795312784-02766a79f207f67626cef289)](https://stats.uptimerobot.com/jjqW4Ulymx)


![schematic](SigProfilerTopography.png)

# SigProfilerTopography

SigProfilerTopography allows evaluating the effect of chromatin organization, histone modifications, transcription factor binding, DNA replication, and DNA transcription on the activities of different mutational processes. SigProfilerTopography elucidates the unique topographical characteristics of mutational signatures. The tool seamlessly integrates with other SigProfiler tools including SigProfilerMatrixGenerator, SigProfilerSimulator, and SigProfilerAssignment. Detailed documentation can be found at: https://osf.io/5unby/wiki/home/

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

The framework is written in PYTHON, however, it also requires the following software with the given versions (or newer):

  * PYTHON          version 3.8 or newer
  * WGET               version 1.9  or RSYNC if you have a firewall

**QUICK START GUIDE**

This section will guide you through the minimum steps required to run SigProfilerTopography:

1. For most recent stable PyPI version of this tool, install the python package using pip:
	```
	$ pip install SigProfilerTopography
	```
	If you have installed SigProfilerTopography before, upgrade using pip:
	```
	$ pip install SigProfilerTopography --upgrade
	```

<!---
```To install the current version of this Github repo, git clone this repo or download the zip file. Unzip the contents of SigProfilerTopography-master.zip or the zip file of a corresponding branch.

	In the command line, please run the following:
	```bash
	$ cd SigProfilerTopography-master
	$ pip install .
	```
	```
-->


2. Imports the example data that is provided by SigProfilerTopography. This data can be used to run the example program and ensure that the environment is set up.

	```python 
	>>> from SigProfilerTopography import Topography as topography
	>>> topography.install_example_data()
	```

	Imports `21BRCA.zip` under the current working directory. Once `21BRCA.zip` has been downloaded, unzip the file. The unzipped `21BRCA` folder contains two folders: `21BRCA_vcfs` and `21BRCA_probabilities`. The folder `21BRCA_vcfs` contains 21 VCF files (one per each breast cancer sample) in GRCh37 and 21BRCA_probabilities` contains probability matrix files for single base substitutions and doublet base substitutions.

3. Install your desired reference genome from the command line/terminal as follows (available reference genomes are: GRCh37, GRCh38, mm9, and mm10):

	```python
	$ python
	>>> from SigProfilerMatrixGenerator import install as genInstall
	>>> genInstall.install('GRCh37')
	```
	
	This will install the human 37 assembly as a reference genome. 

4. Imports the nucleosome library file that is necessary for nucleosome occupancy analyses. Next, choose the genome that you would like to import:  

	```python
	>>> from SigProfilerTopography import Topography as topography
	>>> topography.install_nucleosome("GRCh37")
	```
	By default, `install_nucleosome` imports nucleosome data of `K562` cell line for GRCh37 and GRCh38 genome assemblies.

5. Imports the open chromatin library file that is necessary for epigenomics analyses. Next, choose the genome that you would like to import:  
	
	```python
	>>> from SigProfilerTopography import Topography as topography
	>>> topography.install_atac_seq("GRCh37")
	```
	By default,  `install_atac_seq` imports open chromatin data of `breast epithelium` tissue for GRCh37 and `left lung` tissue for GRCh38.


6. Imports the replication timing library file that is necessary for replication timing analyses. Next, choose the genome that you would like to import:  

	```python
	>>> from SigProfilerTopography import Topography as topography
	>>> topography.install_repli_seq("GRCh37")
	```
	By default, `install_repli_seq` imports replication time data of `MCF7` and `IMR90` for GRCh37 and GRCh38, respectively.
	
7. Conducts topography analyses for your samples. Here is an example of a call to  `runAnalyses`  that generates all of the different analyses.

	```python 
	>>> from SigProfilerTopography import Topography as topography

	>>> genome = "GRCh37"
	>>> inputDir = "path/to/21BRCA_vcfs"
	>>> outputDir = "path/to/results"
	>>> jobname = "21BRCA_SPT"
	>>> numofSimulations = 5

	>>> topography.runAnalyses(genome, 
	                   inputDir, 
	                   outputDir, 
	                   jobname, 
	                   numofSimulations, 
	                   epigenomics=True,
	                   nucleosome=True, 
	                   replication_time=True, 
	                   strand_bias=True, 
	                   processivity=True)
	```

	If probability files are not provided, SigProfilerTopography utilizes SigProfilerAssignment by default to attribute the activities of known reference mutational signatures from the Catalogue Of Somatic Mutations In Cancer (COSMIC) database to each examined sample.

8. Here is an example of a call to  `runAnalyses`  with probability files using the 21 VCF files located in the subfolder `21BRCA_vcfs` as input and providing the probability files in the subfolder `21BRCA_probabilities`.

	```python 
	>>> from SigProfilerTopography import Topography as topography

	>>> genome = "GRCh37"
	>>> inputDir = "path/to/21BRCA_vcfs"
	>>> outputDir = "path/to/results"
	>>> jobname = "21BRCA_SPT_with_probability_matrices"
	>>> numofSimulations = 5
	>>> sbs_probability_file = "path/to/21BRCA_probabilities/COSMIC_SBS96_Decomposed_Mutation_Probabilities.txt"
	>>> dbs_probability_file = "path/to/21BRCA_probabilities/COSMIC_DBS78_Decomposed_Mutation_Probabilities.txt"

	>>> topography.runAnalyses(genome, 
	                   inputDir, 
	                   outputDir, 
	                   jobname, 
	                   numofSimulations, 
	                   sbs_probabilities = sbs_probability_file,
	                   dbs_probabilities = dbs_probability_file,
	                   epigenomics=True,
	                   nucleosome=True, 
	                   replication_time=True, 
	                   strand_bias=True, 
	                   processivity=True)
	```

	SigProfilerTopography utilizes probability matrix files containing the probability of each signature to cause a specific mutation type in a cancer sample. 

View the table below for the full list of `runAnalyses` parameters.

**PARAMETERS**
| Category | Parameter | Variable Type | Parameter Description |
| ------ | ----------- | ----------- | ----------- |
| Required |  |  |  |
|  | **genome** | String | The reference genome used for the topography analyses. Accepted values include: {"GRCh37", "GRCh38", "mm10"}.  |
| | **inputDir** | String | The path to the directory containing the input files. SigProfilerTopography accepts all input files that SigProfilerMatriXGenerator can process. |
|  | **outputDir** | String | The path of the directory where the output will be saved. If this directory doesn't exist, a new one will be created. |
| |  **jobname** | String | The name of the directory containing all of the outputs under `outputDir/jobname`. If this directory doesn't exist, a new one will be created.   |
| | **numofSimulations** | Integer | The number of simulations to be created. |
| Optional |  |  |  |
| | **epigenomics** | Boolean | Generate epigenomics analysis when True. By default, this is set to False. |
| | **nucleosome** | Boolean | Generate nucleosome occupancy analysis when True. By default, this is set to False. |
|  | **replication_time** | Boolean | Generate replication timing analysis when True. By default, this is set to False. |
| | **strand_bias** | Boolean | Generate replication and transcription strand asymmetry analysis when True. By default, this is set to False. |
| | **replication_strand_bias** | Boolean | Generate replication strand asymmetry analysis when True. By default, this is set to False. |
| | **transcription_strand_bias** | Boolean | Generate transcription strand asymmetry analysis (including genic versus intergenic regions) when True. By default, this is set to False. |
| | **processivity** | Boolean | Generate strand-coordinated mutagenesis when True. By default, this is set to False. |
| | **epigenomics_files** | List of Strings | Python list of paths for each epigenomics library file utilized in the epigenomics analysis. By default, epigenomics files of open chromatin, CTCF and histone modifications attained from "breast_epithelium" and "lung" tissue are utilized for GRCh37 and GRCh38, respectively. |
| | **epigenomics_dna_elements** | List of Strings | Python list of unique DNA element names for the epigenomics files utilized in the epigenomics analysis. Each DNA element name must be contained in at least one epigenomics library filename. E.g., DNA element is 'CTCF' for the epigenomics file of 'ENCFF782GCQ_breast_epithelium_Normal_CTCF-human.bed'. By default, DNA elements of ['H3K27me3', 'H3K36me3', 'H3K9me3', 'H3K27ac', 'H3K4me1', 'H3K4me3', 'CTCF', 'ATAC'] are utilized for GRCh37 and GRCh38. If user provided `epigenomics_files` is provided, then `epigenomics_dna_elements` is mandatory. |
| | **epigenomics_biosamples** | List of Strings | Python list of unique biosample names for the epigenomics files utilized in the epigenomics analyses. Each biosample name must be contained in at least one epigenomics library filename. E.g., biosample is 'breast_epithelium' for the epigenomics file of 'ENCFF782GCQ_breast_epithelium_Normal_CTCF-human.bed'. By default, "breast_epithelium" and "lung" biosamples are utilized for GRCh37 and GRCh38, respectively. Biosamples are shown in the epigenomics heatmaps if `plot_detailed_epigemomics_heatmaps` is set to True. |
| | **nucleosome_biosample** | String | Biosample that will be used for nucleosome occupancy analysis. Analysis can be done by using either K562 or GM12878 cell line from ENCODE. By default, the K562 cell line is used for GRCh37 and GRCh38. |
| | **nucleosome_file** | String | The path to the nucleosome occupancy library file that will be used for the analysis. By default, nucleosome occupancy file (MNase-seq) of K562 cell line is used for GRCh37 and GRCh38. |
| | **replication_time_biosample** | String | Biosample that will be used to carry out replication timing and replication strand asymmetry analyses. By default, MCF7 and IMR90 cell lines are utilized for GRCh37 and GRCh38, respectively. For the complete list of available replication time biosamples, refer to the Replication Time Biosamples table below. |
| | **replication_time_signal_file** | String | The path to the replication time signal file. By default, replication time signal file (wig file) of MCF7 and IMR90 cell lines are utilized for GRCh37 and GRCh38, respectively. |
| | **replication_time_valley_file** | String | The path to the replication time valley file. By default, replication time valley file (bed file) of MCF7 and IMR90 cell lines are utilized for GRCh37 and GRCh38, respectively. |
| | **replication_time_peak_file** | String | The path to the replication time peak file. By default, replication time peak file (bed file) of MCF7 and IMR90 cell lines are utilized for GRCh37 and GRCh38, respectively. |
| | **samples_of_interest** | List of Strings | Conduct topography analyses for these samples of interest only. By default, it is set to None and topography analyses are carried out for all samples. |
| | **discreet_mode** | Boolean | Each mutation contributes to the topography analyses either with 1 or 0 when True; otherwise, each mutation contributes with its probability when False. By default, this is set to True.   |
| | **average_probability** | Float | The average probability of the mutations assigned to a SBS, DBS, and ID signature. By default, it is set to 0.90. The `average_probability` applies when `discreet_mode` is True. We set signature specific cutoffs, such that for the mutations satisfying mutation_signature_probability >= cutoff, average probability of these mutations must be at least 0.90. |
| | **num_of_sbs_required** | Integer | The minimum required number of mutations for a SBS signature. The `num_of_sbs_required` applies when `discreet_mode` is True or when `discreet_mode` is False and `show_all_signatures` is False. By default, it is set to 2000. |
| | **num_of_dbs_required** | Integer | The minimum required number of mutations for a DBS signature. The `num_of_dbs_required` applies when `discreet_mode` is True or when `discreet_mode` is False and `show_all_signatures` is False. By default, it is set to 200. |
| | **num_of_id_required** | Integer | The minimum required number of mutations for a ID signature. The `num_of_id_required` applies when `discreet_mode` is True or when `discreet_mode` is False and `show_all_signatures` is False. By default, it is set to 1000. |
| | **exceptional_signatures** | Dictionary | The dictionary of exceptional signatures. The `exceptional_signatures` applies when `discreet_mode` is True. E.g., `exceptional_signatures` = {"SBS32" : 0.63} is a Python dictionary where key is a *mutational signature* and value is an *average probability*.  Exceptional signatures are included in the topography analyses if they satisfy `num_of_sbs_required`, `num_of_dbs_required`, and `num_of_id_required` constraints with `average_probability` >= given *average probability*. |
| | **default_cutoff** | Float | The `default_cutoff` applies for all signatures when `discreet_mode` is False. Mutations satisfying mutation_signature_probability >= `default_cutoff` are considered in the topography analyses with their probability. By default, it is set to 0.5. |
| | **show_all_signatures** | Boolean | The `show_all_signatures` applies when `discreet_mode` is False. All signatures are considered in the topography analyses when True, otherwise signatures satisfying `num_of_sbs_required`, `num_of_dbs_required`, and `num_of_id_required` are considered in the topography analyses when False. By default, it is set to True. |
| | **plot_figures** | Boolean | Generate plots displaying the results of all topography analyses when True. By default, this is set to True. |
| | **plot_epigenomics** | Boolean | Generate epigenomics heatmaps and occupancy plots when True. By default, this is set to False. |
| | **plot_nucleosome** | Boolean | Generate nucleosome occupancy plots when True. By default, this is set to False. |
| | **plot_replication_time** | Boolean | Generate replication timing plots when True. By default, this is set to False. |
| | **plot_strand_bias** | Boolean | Generate replication strand asymmetry, transcription strand asymmetry, genic versus intergenic regions plots when True. By default, this is set to False. |
| | **plot_replication_strand_bias** | Boolean | Generate replication strand asymmetry plots when True. By default, this is set to False. |
| | **plot_transcription_strand_bias** | Boolean | Generate transcription strand asymmetry and genic versus intergenic regions plots when True. By default, this is set to False. |
| | **plot_processivity** | Boolean | Generate strand-coordinated mutagenesis plots when True. By default, this is set to False. |
| | **step1_matgen_real_data** | Boolean | Run SigProfilerMatrixGenerator to generate matrices for the real mutations when True. By default, this is set to True. |
| | **step2_gen_sim_data** | Boolean | Run SigProfilerSimulator to generate simulated mutations when True. By default, this is set to True. |
| | **step3_matgen_sim_data** | Boolean | Run SigProfilerMatrixGenerator to generate matrices for the simulated mutations when True. By default, this is set to True. |
| | **step4_merge_prob_data** | Boolean | Merge real and simulated mutations with the probabilities files when True. By default, this is set to True. |
| | **step5_gen_tables** | Boolean | Generate tables for providing information on mutational signatures, cutoffs, number of mutations and average probability when True. By default, this is set to True. |
| | **sbs_probabilities** | String | The path to the probabilities matrix file. The probabilities matrix includes the probabilities of each mutation type in each sample. The first column lists all the samples, the second column lists all the mutation types, and the following columns list the calculated probability value for the respective SBS signatures where the sum of each row is 1. The probabilities file can be in *SBS_6*,  *SBS_24* *SBS_96*, *SBS_192*, *SBS_288*, *SBS_384*, *SBS_1536*, or *SBS_6144* context produced by mutational signature extractor. |
| | **dbs_probabilities** | String | The path to the probabilities matrix file. The probabilities matrix includes the probabilities of each mutation type in each sample. The first column lists all the samples, the second column lists all the mutation types, and the following columns list the calculated probability value for the respective DBS signatures where the sum of each row is 1. The probabilities file in DBS-78 context produced by mutational signature extractor. |
| | **id_probabilities** | String |  The path to the probabilities matrix file. The probabilities matrix includes the probabilities of each mutation type in each sample. The first column lists all the samples, the second column lists all the mutation types, and the following columns list the calculated probability value for the respective ID signatures where the sum of each row is 1. The probabilities file in ID-83 context produced by mutational signature extractor. |
| | **sbs_signatures** | String |  The path to the signatures matrix file. The signatures matrix contains the distribution of mutation types in the SBS mutational signatures. The first column lists all of the mutation types. e.g., There are 96 possible mutations that are considered for the SBS-96 context. The following columns are the SBS signatures. The sum of each column is 1, and each value in a column indicates the proportion of a mutational context in the signature. |
| | **dbs_signatatures** | String | The path to the signatures matrix file. The signatures matrix contains the distribution of mutation types in the DBS mutational signatures. The first column lists all of the mutation types. e.g., There are 78 possible mutations that are considered for the DBS-78 context. The following columns are the DBS signatures. The sum of each column is 1, and each value in a column indicates the proportion of a mutational context in the signature. |
| | **id_signatures** | String | The path to the signatures matrix file. The signatures matrix contains the distribution of mutation types in the ID mutational signatures. The first column lists all of the mutation types. e.g., There are 83 possible mutations that are considered for the ID-83 context. The following columns are the ID signatures. The sum of each column is 1, and each value in a column indicates the proportion of a mutational context in the signature. |
| | **sbs_activities** | String | The path to the activities matrix file. The activity matrix for the selected SBS signatures. The first column lists all of the samples and the second and the following columns list the calculated activity value (number of mutations) for the respective SBS signatures. |
| | **dbs_activities** | String | The path to the activities matrix file. The activity matrix for the selected DBS signatures. The first column lists all of the samples and the second and the following columns list the calculated activity value (number of mutations) for the respective DBS signatures. |
| | **id_activities** | String | The path to the activities matrix file. The activity matrix for the selected ID signatures. The first column lists all of the samples and the second and the following columns list the calculated activity value (number of mutations) for the respective ID signatures. |
| | **verbose** | Boolean | Set to True for detailed debugging messages. By default, this is set to False. |
| | **parallel_mode** | Boolean | Set to True for running SigProfilerTopography using multiprocessing. By default, this is set to True. |
| | **plusorMinus_epigenomics** | Integer | The number of bases considered before and after mutation start for epigenomics occupancy analysis. |
| | **plusorMinus_nucleosome** | Integer | The number of bases considered before and after  mutation start for nucleosome occupancy analysis. |
| | **epigenomics_heatmap_<br>significance_level** | Float | Corrected p-values <= `epigenomics_heatmap_significance_level` are considered statistically significant. By default, this is set to 0.05. |
| | **fold_change_window_size** | Integer | In epigenomics analysis, fold change of real versus simulated mutations is calculated for the window size centered at the mutation start. E.g., for window size of 100 bases, ± 50 bases are considered before and after mutation start. By default, this is set to 100. |
| | **num_of_avg_overlap** | Integer | The minimum required average number of overlaps between the mutations and the regions outlined in the epigenomics files. By default, set to 100.|
| | **plot_detailed_epigemomics_<br>heatmaps** | Boolean | Plot detailed epigenomics heatmaps when True. By default, set to False. |
| | **remove_dna_elements_with_all_<br>nans_in_epigemomics_heatmaps** | Boolean | Remove the DNA elements from the epigenomics heatmap if no result exists. By default, set to True. |
| | **odds_ratio_cutoff** | Float | Strand asymmetries with odd ratio >= `odds_ratio_cutoff` are shown in the strand asymmetry circle plots. By default, set to 1.1. |
| | **percentage_of_real_<br>mutations_cutoff** | Float | Strand asymmetries of the SBS signatures with percentage of the mutations >= `percentage_of_real_mutations_cutoff` are shown in the plots. By default, set to 5. |
| | **ylim_multiplier** | Float | Multiply the y-axis view limits with `ylim_multiplier` in strand asymmetry bar plots. By default, set to 1.25. |
| | **processivity_inter_<br>mutational_distance** | Integer | Consecutive mutations with distance <= `processivity_inter_mutational_distance` are considered for the strand-coordinated mutagenesis. By default, set to 10000. |
| | **processivity_significance_level** | Float | Corrected p-values <= `processivity_significance_level` are considered statistically significant for strand coordinated mutagenesis. By default, this is set to 0.05. |
| | **exome** | Boolean | SigProfilerSimulator simulates on the exome of the reference genome. By default, set to None. |
| | **updating** | Boolean | SigProfilerSimulator updates the chromosome with each mutation. By default, set to False. |
| | **bed_file** | String | SigProfilerSimulator simulates on custom regions of the genome. Requires the full path to the BED file. By default, set to None. |
| | **overlap** | Boolean | SigProfilerSimulator allows overlapping of mutations along the chromosome. By default, set to False. |
| | **gender** | String | SigProfilerSimulator simulates male or female genomes. By default, set to 'female'. |
| | **seed_file** | String | SigProfilerSimulator uses this path to user defined seeds. One seed is required per processor. Uses a built in file by default. By default, this is set to None. |
| | **noisePoisson** | Boolean | SigProfilerSimulator adds poisson noise to the simulations. By default, set to False. |
| | **noiseUniform** | Integer | SigProfilerSimulator adds a noise dependent on a +/- allowance of noise (e.g., noiseUniform=5 allows +/-2.5% of mutations for each mutation type). By default, this is set to  0. |
| | **cushion** | Integer | SigProfilerSimulator allows cushion when simulating on the exome or targetted panel. By default, this is set to 100 base pairs.  |
| | **region** | String | For SigProfilerSimulator. Path to targetted region panel for simulated on a user-defined region. Default is whole-genome simulations. |
| | **vcf** | Boolean | SigProfilerSimulator outputs simulated samples as vcf files with one file per iteration per sample when True. SigProfilerSimulator outputs all samples from an iteration into a single maf file when False. By default, this is set to False. |
| | **mask** | String | For SigProfilerSimulator. Path to probability mask file. A mask file format is tab-separated with the following required columns: Chromosome, Start, End, Probability. Note: Mask parameter does not support exome data where bed_file flag is set to true, and the following header fields are required: Chromosome, Start, End, Probability. By default, this is set to None. |
|  |  |  | |


**SigProfilerTopography Output**
To learn about the output, please visit https://osf.io/5unby/wiki/home/
  
<!---

**Replication Time Biosamples**
 For GRCh37 and GRCh38, SigProfilerTopography provides Repli-seq files of the biosamples listed in the table below, which are the valid parameter values for `replication_time_biosample`.


| Biosample | Organism | Tissue | Cell Type | Diseases |
| --------- | --------- | -------- | -------- | -------- |
| **MCF7** | human | breast | mammary | Cancer |
| **HEPG2** | human | liver | liver cells (hepatocytes) | Cancer |
| **HELAS3** | human | cervix | epithelial-like cervical cells | Cancer |
| **SKNSH** | human | brain | neuronal-like cells  | Cancer |
| **K562** | human | bone marrow | lymphoblast cells | Cancer |
| **IMR90** | human | lung | fibroblast | Normal |
| **NHEK** | human | skin | keratinocyte | Normal |
| **BJ** | human | skin | fibroblast | Normal |
| **HUVEC** | human | skin | fibroblast | Normal |
| **BG02ES** | human | early developmental stage of an embryo, not from a differentiated tissue | embyronic stem cell | None reported |
| **GM12878** | human | blood | B-Lymphocyte | Normal |
| **GM06990** | human | blood | B-Lymphocyte | Unknown |
| **GM12801** | human | blood |  B-Lymphocyte | Unknown |
| **GM12812** | human | blood | B-Lymphocyte | Unknown |
| **GM12813** | human | blood |  B-Lymphocyte | Unknown |
| | | | | |

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
	
		
		````
		>>> topography.runAnalyses(genome,     inputDir,outputDir,jobname,numofSimulations,sbs_probabilities=sbs_probabilities_file_path,dbs_probabilities=dbs_probabilities_file_path,replication_time_signal_file=user_replication_time_signal_file,replication_time_valley_file=user_replication_time_valley_file,replication_time_peak_file=user_replication_time_peak_file,replication_time=True,strand_bias=True)
		````
-->


## <a name="citation"></a> Citation
Otlu B, Alexandrov LB: Evaluating topography of mutational signatures with SigProfilerTopography. __BioRxiv 2024__, [https://doi.org/10.1101/2024.01.08.574683](https://doi.org/10.1101/2024.01.08.574683).

Otlu B, Diaz-Gay M, Vermes I, Bergstrom EN, Zhivagui M, Barnes M, Alexandrov LB: Topography of mutational signatures in human cancer. __Cell Rep 2023__,  [https://doi.org/10.1016/j.celrep.2023.112930](https://doi.org/10.1016/j.celrep.2023.112930).


## <a name="copyright"></a> Copyright
This software and its documentation are copyright 2018 as a part of the SigProfiler project. The SigProfilerTopography framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.


## <a name="contact"></a> Contact Information
Please address any queries or bug reports to Burcak Otlu at burcakotlu@eng.ucsd.edu

