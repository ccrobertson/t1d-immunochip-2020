#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account rich_immunochip_impute
#SBATCH --partition=standard
#SBATCH --mem=10000


module load gcc/7.1.0  openmpi/3.1.4
module load R/3.6.1

### Format and subset data for matrix EQTL
Rscript ${scripts}/caqtl/prepare_data_for_MatrixEQTL.R
