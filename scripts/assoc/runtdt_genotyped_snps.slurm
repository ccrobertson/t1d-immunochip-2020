#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account rich_immunochip_impute
#SBATCH --partition=largemem
#SBATCH --mem=20000


group=$1

keeplist=${fam_assoc}/family_analysis_keep_list_${group}.txt

cd ${fam_assoc}/${group}

### Run TDT in genotyped data
#extract family data set
plink --bfile ${genodat} --keep ${keeplist} --make-bed --out genotyped
#swap case-control status
awk 'BEGIN {a[1]=2; a[2]=1; a[0]=0; a[-9]=0} {print $1, $2, a[$6]}'  genotyped.fam >  genotyped_swap_cc_status.txt
plink --bfile genotyped --pheno genotyped_swap_cc_status.txt --make-bed --out genotyped_WRONG
#run basic tdt analysis
plink --bfile genotyped --tdt --ci 0.95 --out genotyped
plink --bfile genotyped_WRONG --ci 0.95 --tdt --out genotyped_WRONG
