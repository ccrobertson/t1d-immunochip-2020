#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account rich_immunochip_impute
#SBATCH --partition=standard
#SBATCH --mem=20000


gene=$1


echo -e "STARTING ${gene}"
date

zcat ${PUBLIC_DATA}/eQTLGen/cis-eQTLs_full_20180905.txt.gz | head | awk 'NR==1' > ${PUBLIC_DATA}/eQTLGen/cis-eQTLs_full_20180905__${gene}.txt
zcat ${PUBLIC_DATA}/eQTLGen/cis-eQTLs_full_20180905.txt.gz | awk -v gene=${gene} '$9==gene' >> ${PUBLIC_DATA}/eQTLGen/cis-eQTLs_full_20180905__${gene}.txt

echo -e "FINISHED ${gene}"
date
echo -e " "
echo -e " "
