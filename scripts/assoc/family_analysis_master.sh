#!/bin/bash

group=$1

mkdir -p family_analysis


### Convert data to plink and update fam files
echo "CONVERT TO PLINK"
bash ${scripts}/assoc/convert_to_plink.sh ${group}


### Run TDT in imputed data
echo "RUN TDT"
bash ${scripts}/assoc/runtdt.sh ${group}


### Calculate Mendelian Error rates
echo "RUN KING"
sbatch ${scripts}/assoc/snp_qc_mi.slurm ${group}
