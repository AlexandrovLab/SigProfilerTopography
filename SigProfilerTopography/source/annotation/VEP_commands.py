# !/usr/bin/env python3

# Author: burcakotlu

# Contact: burcakotlu@eng.ucsd.edu

"""
This python code prepares bash script consisting of  VEP command lines to run VEP
for all the samples under the given directory.

"""

import os

# LICA-FR_SP97771.snv_mnv_filtered.vcf
# LICA-FR_SP97771.indel_filtered.vcf

def get_vcf_samples(vcf_dir, substring=None):
    samples = []

    if os.path.exists(vcf_dir):
        files_list = os.listdir(vcf_dir)
        if substring:
            samples = [file[:-4] for file in files_list if (file.endswith('.vcf') and (file.find(substring) != -1))]
        else:
            samples = [file[:-4] for file in files_list if file.endswith('.vcf') ]

    return samples


def write_bash_script(cancer_type, mutation_type, samples, vcf_dir, vep_dir):
    bash_file_name = '%s_%s_VEP_run_commands.sh' %(cancer_type, mutation_type)
    f = open(bash_file_name, "w")

    f.write('#!/bin/bash\n')
    f.write('# Bash script\n')
    f.write('\n')

    f.write('cd $PBS_O_WORKDIR\n')
    f.write('export MODULEPATH=/projects/builder-group/jpg/modulefiles/applications:$MODULEPATH\n')
    f.write('module load vep\n')
    f.write('\n')

    for sample in samples:
        text = 'vep -i %s/%s.vcf -o %s/%s.txt -offline --dir_cache /projects/builder-group/data/vep-96' %(vcf_dir, sample, vep_dir, sample)
        f.write(text + '\n\n')
    f.close()


# cancer_types = ['Liver-HCC', 'Lung-AdenoCA']
cancer_types = ['Lung-AdenoCA']

for cancer_type in cancer_types:
    pcawg_vcf_dir = os.path.join('/restricted/alexandrov-group/burcak/data/PCAWG/' + cancer_type + '/filtered')
    pcawg_vep_dir = os.path.join('/restricted/alexandrov-group/burcak/data/PCAWG/' + cancer_type + '/vep_files')

    os.makedirs(pcawg_vep_dir, exist_ok=True)

    snv_mnv_samples = get_vcf_samples(pcawg_vcf_dir, 'snv_mnv_filtered')
    indel_samples = get_vcf_samples(pcawg_vcf_dir, 'indel_filtered')

    write_bash_script(cancer_type, 'snv_mnv', snv_mnv_samples, pcawg_vcf_dir, pcawg_vep_dir)
    write_bash_script(cancer_type, 'indel', indel_samples, pcawg_vcf_dir, pcawg_vep_dir)

