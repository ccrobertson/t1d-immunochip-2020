#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account rich_immunochip_impute
#SBATCH --partition=standard
#SBATCH --mem=200000


module load gcc/7.1.0
module load R/3.5.1

group=$1
refpanel=$2
chr=${SLURM_ARRAY_TASK_ID}

echo "CALCULATING ACCURACY WITHIN ICHIP REGIONS ONLY"
Rscript ${scripts}/imputation_accuracy_coverage/calculate_imputation_accuracy.R ${group} ${refpanel} ichip ${chr}
date

echo "CALCULATING ACCURACY IN ALL REGIONS"
Rscript ${scripts}/imputation_accuracy_coverage/calculate_imputation_accuracy.R ${group} ${refpanel} all ${chr}
date

echo "DONE"
