#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account rich_immunochip_impute
#SBATCH --partition=largemem
#SBATCH --mem=20000

group=$1
dir=${fam_assoc}/${group}

Rscript ${scripts}/assoc/family_based_gwas_plots.R ${group} ${dir}
