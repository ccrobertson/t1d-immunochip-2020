#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account rich_immunochip_impute
#SBATCH --partition=standard
#SBATCH --mem=20000



region=$1

awk -v reg=${region} '$2==reg' ${qtl}/chromQTL/chromQTL_unstim_cis.txt > ${qtl}/chromQTL/chromQTL_unstim_cis_${region}.txt
