#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --account rich_immunochip_impute
#SBATCH --partition=standard
#SBATCH --mem=10000


dataset=$1

module load peer
cd ${freeze}
python -u ${scripts}/caqtl/peer_correction.py ${dataset}
