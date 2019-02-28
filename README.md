# SigProfilerTopography
SigProfilerTopography provides topography analyses for substitutions, dinucleotides and indels for each sample and all samples pooled.


**QUICK START GUIDE**

This section will guide you through the minimum steps required to run SigProfilerTopography:
1. Install the python package using pip:
```
                          pip install SigProfilerTopography
```
2. Import SigProfilerTopography as follows:
```
$ python
>> from SigProfilerTopography import SigProfilerTopography as topography
>> import os
```


3. Place your vcf files in your desired output folder. It is recommended that you name this folder based on your project's name
4. From within a python session, you can now generate the matrices as follows:
```
$ python3
>>current_abs_path = os.path.abspath(os.path.dirname(__file__))
>>jobname = '21BreastCancer'
>>dataDir = '%s/SigProfilerTopography/input/%s' %(current_abs_path,jobname)
>>snpsForTopography= '%s/%s_snps_for_topography.txt' %(dataDir,jobname)
>>indelsForTopography= '%s/%s_indels_for_topography.txt' %(dataDir,jobname)
>>topography.runAnalyses(snpsForTopography,indelsForTopography,jobname,'wgEncodeSydhNsomeGm12878Sig.wig','GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.wig','GSM923442_hg19_wgEncodeUwRepliSeqMcf7ValleysRep1.bed','GSM923442_hg19_wgEncodeUwRepliSeqMcf7PkRep1.bed')

```
  The layout of the parameters are as follows:

      SigProfilerMatrixGeneratorFunc(project, reference_genome, path_to_input_files, plot)

  where project, reference_genome, and path_to_input_files must be strings (surrounded by quotation marks, ex: "test") and plot is a boolean argument (True or False)


**INPUT FILE FORMAT**

This tool currently supports simple text file format. The user must provide variant data with its path as in the case of snpsForTopography and indelsForTopography in the sample run.

**LIBRARY**

This tool uses ENCODE provided files for topography analysis such as nucleosome occupancy, replcation time files.
These files have to provided under

SigProfilerTopography/lib/

By the way

SigProfilerTopography/lib/ and SigProfilerTopography/source/ are at the same level.

**LIBRARY NUCLEOSOME OCCUPANCY**

SigProfilerTopography/lib/nucleosome/wgEncodeSydhNsomeGm12878Sig.wig

**LIBRARY REPLICATION TIME**

SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7WaveSignalRep1.wig
SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7PkRep1.bed
SigProfilerTopography/lib/replication/GSM923442_hg19_wgEncodeUwRepliSeqMcf7ValleysRep1.bed

**LIBRARY TRANSCRIPTS**

SigProfilerTopography/lib/transcripts/GRCh37_transcripts.txt

**LIBRARY HG19 and HG38 twobit files **

SigProfilerTopography/lib/ucscgenome/hg19.2bit
SigProfilerTopography/lib/ucscgenome/hg38.2bit


**COPYRIGHT**

This software and its documentation are copyright 2018 as a part of the sigProfiler project. The SigProfilerTopography framework is free software and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

**CONTACT INFORMATION**

Please address any queries or bug reports to Burcak Otlu at burcakotlu@eng.ucsd.edu