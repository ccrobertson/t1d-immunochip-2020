#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --account rich_immunochip_impute
#SBATCH --partition=standard
#SBATCH --mem=50000


module load gcc/7.1.0  openmpi/3.1.4
module load R/3.6.1

Rscript ${scripts}/caqtl/conditional_caqtl.R
