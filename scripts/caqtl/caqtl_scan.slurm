#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --account rich_immunochip_impute
#SBATCH --partition=largemem
#SBATCH --mem=500000


#module load gcc/7.1.0  openmpi/3.1.4
#module load R/3.6.1

module load gcc/7.1.0
module load R/3.5.1

#cd /nv/vol185/T1DGC/T1DGC_DP3/analysis
#R -e "rmarkdown::render('analysis_first_pass.Rmd')"

Rscript ${scripts}/caqtl/caqtl_scan.R
