#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account rich_immunochip_impute
#SBATCH --partition=standard
#SBATCH --mem=50000


date

getCombinedFile () {
	local group=$1
	local refpanel=$2
	local subset=$3
	local
	outfile=${impbench}/imputation_vs_wgs_${refpanel}_${group}_${subset}_regions.txt
	awk 'NR==1' ${impbench}/imputation_vs_wgs_${refpanel}_${group}_chr1_${subset}_regions.txt > ${outfile}
	for i in {1..22}; do
    echo "chr${i}"
		awk 'NR>1' ${impbench}/imputation_vs_wgs_${refpanel}_${group}_chr${i}_${subset}_regions.txt  >> ${outfile}
	done
}

getCombinedFile AFR IMPUTED_TOPMED ichip
getCombinedFile AFR IMPUTED_1KG ichip

getCombinedFile AMR IMPUTED_TOPMED ichip
getCombinedFile AMR IMPUTED_1KG ichip

getCombinedFile EUR IMPUTED_TOPMED ichip
getCombinedFile EUR IMPUTED_1KG ichip

#getCombinedFile AFR IMPUTED_TOPMED all
#getCombinedFile AFR IMPUTED_1KG all

#getCombinedFile AMR IMPUTED_TOPMED all
#getCombinedFile AMR IMPUTED_1KG all

echo "DONE"
date
