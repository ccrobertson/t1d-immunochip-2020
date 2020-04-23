# Imputation coverage and accuracy

Get WashU WGS data
	```bash

	#find overlapping subjects and create keeplists
	bash ${scripts}/imputation_accuracy_coverage/get_washu_mega_overlap_and_keeplists.sh
	#extract mega subjects and filter for MAF>0.005
  sbatch --array=1-22 --output=${impbench}/filter_washu_files_AFR_%a.log ${scripts}/imputation_accuracy_coverage/filter_washu_files.slurm AFR
	sbatch --array=1-22 --output=${impbench}/filter_washu_files_AMR_%a.log ${scripts}/imputation_accuracy_coverage/filter_washu_files.slurm AMR
	sbatch --array=1-22 --output=${impbench}/filter_washu_files_EUR_%a.log ${scripts}/imputation_accuracy_coverage/filter_washu_files.slurm EUR

	#some of the filtering timed out so had to rerun
	sbatch --output=${impbench}/filter_washu_files_AFR_chr2_restart.log ${impbench}/filter_washu_files_AFR_chr2.slurm
	sbatch --output=${impbench}/filter_washu_files_AMR_chr1_restart.log ${impbench}/filter_washu_files_AMR_chr1.slurm
  ```

Prep data sets
	```bash
	prepDataForAccuracy () {
		local group=$1
		sbatch --array=1-22 --output=${impbench}/prep_washu_data_${group}_%a.log ${scripts}/imputation_accuracy_coverage/prep_washu_data.slurm ${group}
		#sbatch --array=1-22 --output=${impbench}/prep_mega_data_IMPUTED_TOPMED_${group}_%a.log ${scripts}/imputation_accuracy_coverage/prep_mega_data.slurm ${group} IMPUTED_TOPMED gr38
		#sbatch --array=1-22 --output=${impbench}/prep_mega_data_IMPUTED_1KG_${group}_%a.log ${scripts}/imputation_accuracy_coverage/prep_mega_data.slurm ${group} IMPUTED_1KG gr37
	}
	prepDataForAccuracy AFR
	prepDataForAccuracy	AMR
	prepDataForAccuracy	EUR
	prepDataForAccuracy	FIN

	NOTE FOR AFR - ERROR in CHR2 and CHR14 due to duplicate variant IDs in 1000G imputed file (possibly created during liftover?)
	#Manual fix:
	mv mega_IMPUTED_1KG_AFR_chr2_ichip_regions.vcf.gz mega_IMPUTED_1KG_AFR_chr2_ichip_regions_DUPS.vcf.gz
	mv mega_IMPUTED_1KG_AFR_chr14_ichip_regions.vcf.gz mega_IMPUTED_1KG_AFR_chr14_ichip_regions_DUPS.vcf.gz
	zcat mega_IMPUTED_1KG_AFR_chr2_ichip_regions_DUPS.vcf.gz | grep -v "2:61390913:C:<CN0>" | bgzip > mega_IMPUTED_1KG_AFR_chr2_ichip_regions.vcf.gz
	zcat mega_IMPUTED_1KG_AFR_chr14_ichip_regions_DUPS.vcf.gz  | grep -v "14:88504832:\|14:98480754:" | bgzip > mega_IMPUTED_1KG_AFR_chr14_ichip_regions.vcf.gz
	sbatch --array=2,14 --output=${impbench}/calculate_imputation_accuracy_IMPUTED_1KG_AFR_%a.log ${scripts}/imputation_accuracy_coverage/calculate_imputation_accuracy.slurm AFR IMPUTED_1KG
	```

Compute accuracy
	```bash
	calcAccuracy () {
	local group=$1
	sbatch --array=1-22 --output=${impbench}/calculate_imputation_accuracy_IMPUTED_TOPMED_${group}_%a.log ${scripts}/imputation_accuracy_coverage/calculate_imputation_accuracy.slurm ${group} IMPUTED_TOPMED
	sbatch --array=1-22 --output=${impbench}/calculate_imputation_accuracy_IMPUTED_1KG_${group}_%a.log ${scripts}/imputation_accuracy_coverage/calculate_imputation_accuracy.slurm ${group} IMPUTED_1KG
	}
	calcAccuracy AFR
	calcAccuracy AMR
	calcAccuracy EUR

	#some jobs ran out of memory or timed out, re-running these here on largemem
	sbatch --array=1-12 --output=${impbench}/calculate_imputation_accuracy_IMPUTED_TOPMED_AFR_%a_restart.log ${impbench}/calculate_imputation_accuracy_restart.slurm AFR IMPUTED_TOPMED
	sbatch --array=2 --output=${impbench}/calculate_imputation_accuracy_IMPUTED_TOPMED_AMR_%a_restart.log ${impbench}/calculate_imputation_accuracy_restart.slurm AMR IMPUTED_TOPMED
	```

Count genotyped variants in ichip regions
	```bash
	plink1.9 --bfile ${genodat}_filtered_for_hwe --freq --out ${genodat}_filtered_for_hwe
	awk 'NR>1{print $2, $3, $4, $1}' ${PROJECT_MEGA}/define_ichip_regions/ichip_regions_gr37.txt > ${PROJECT_MEGA}/define_ichip_regions/ichip_regions_gr37_plinkformat.txt
	plink1.9 --bfile ${genodat}_filtered_for_hwe --extract ${PROJECT_MEGA}/define_ichip_regions/ichip_regions_gr37_plinkformat.txt --range --make-bed --out ${genodat}_filtered_for_hwe_ichip_regions
	plink1.9 --bfile ${genodat}_filtered_for_hwe_ichip_regions --freq --out ${genodat}_filtered_for_hwe_ichip_regions
	awk '$5>0.005' ${genodat}_filtered_for_hwe_ichip_regions.frq | wc -l
	```

Count imputed variants in ichip regions
	```bash
	bash ${scripts}/imputation_accuracy_coverage/extract_ichip_regions.sh EUR
	bash ${scripts}/imputation_accuracy_coverage/extract_ichip_regions.sh AFR
	bash ${scripts}/imputation_accuracy_coverage/extract_ichip_regions.sh FIN
	bash ${scripts}/imputation_accuracy_coverage/extract_ichip_regions.sh AMR
	```

	```R
	for (group in c("AFR","EUR","FIN","EAS","AMR")) {
	  d = read.table(paste0("count_imputed_variants_in_ichip_regions_",group,".txt"))
	  print(paste(group, sum(d$V2), sum(d$V2)/84923))
	}
	```

Summarize accuracy results
  ```bash
	#merge results across chromosomes
	sbatch --output=${impbench}/merge_accuracy_results.log ${scripts}/imputation_accuracy_coverage/merge_accuracy_results.slurm
	#compare and visualize results for between refpanels
  Rscript ${scripts}/imputation_accuracy_coverage/summarize_accuracy.R ichip IMPUTED_TOPMED
  Rscript ${scripts}/imputation_accuracy_coverage/compare_accuracy_topmed_vs_1kg.R
  ```

