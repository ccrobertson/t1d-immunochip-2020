#!/bin/bash

group=$1

keeplist=${fam_assoc}/family_analysis_keep_list_${group}.txt
dir=${fam_assoc}/${group}

mkdir -p ${dir}
cd ${dir}


awk -v file=${keeplist} 'BEGIN {while (getline<file) {a[$2]=1}} $3 in a {print $3,$3,$2,$3}' ${phenofile} > update_fid_${group}.txt
awk -v file=${keeplist} 'BEGIN {while (getline<file) {a[$2]=1}} $3 in a {print $2,$3,$4,$5}' ${phenofile} > update_parents_${group}.txt
awk -v file=${keeplist} 'BEGIN {while (getline<file) {a[$2]=1}} $3 in a {print $2,$3,$7}' ${phenofile} > update_case_${group}.txt
awk -v file=${keeplist} 'BEGIN {while (getline<file) {a[$2]=1}} $3 in a {print $2,$3,$6}' ${phenofile} > update_sex_${group}.txt
awk -v file=${keeplist} 'BEGIN {a[1]=2; a[2]=1; a[0]=0} {print $1, $2, a[$3]}'  update_case_${group}.txt >  swap_cc_status_${group}.txt


sbatch --array=1-22 ${scripts}/assoc/convert_to_plink.slurm ${group}
