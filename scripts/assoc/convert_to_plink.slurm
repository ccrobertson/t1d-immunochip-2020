#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account rich_immunochip_impute
#SBATCH --partition=standard
#SBATCH --mem=40000


group=$1

chr=${SLURM_ARRAY_TASK_ID}

if [ ! -f convert_finished_chr${chr}.flag ]; then

	echo -e "STARTING CHR${chr}"
	vcffile=${filtered}/${group}/relateds/chr${chr}.filter_maf_gt_0.005_rsq_gt_0.8.relateds.dose.vcf.gz
	prefix=chr${chr}.filter_maf_gt_0.005_rsq_gt_0.8.relateds

	#convert to plink
	echo "CONVERTING TO PLINK"
	plink2 --memory 15000 --vcf ${vcffile} --double-id --make-bed --out ${prefix}
	echo " "

	#update FID
	echo "UPDATING FID"
	plink --memory 15000 --bfile ${prefix} --update-ids update_fid_${group}.txt --make-bed --out ${prefix}2
	echo " "

	#extract complete families
	echo "EXTRACTING COMPLETE TRIOS"
	plink --memory 15000 --bfile ${prefix}2 --keep update_parents_${group}.txt --make-bed --out ${prefix}2b
	echo " "

	#update parent IDs
	echo "UPDATING PARENT IDS"
	plink --memory 15000 --bfile ${prefix}2b --update-parents update_parents_${group}.txt --make-bed --out ${prefix}3
	echo " "

	#update sex
	echo "UPDATING SEX"
	plink --memory 15000 --bfile ${prefix}3 --update-sex update_sex_${group}.txt --make-bed --out ${prefix}4
	echo " "

	#update t1d status
	echo "UPDATING T1D STATUS"
	plink --memory 15000 --bfile ${prefix}4 --pheno update_case_${group}.txt --make-bed --out ${prefix}_with_pheno
	echo " "

	if [ -f ${prefix}_with_pheno.bim ]; then

		echo -e "FINISHED CHR${chr}"

		rm -f ${prefix}2.*
		rm -f ${prefix}2b.*
		rm -f ${prefix}3.*
		rm -f ${prefix}4.*
		if [ -f convert_error_chr${chr}.flag ]; then
			mv convert_error_chr${chr}.flag convert_finished_chr${chr}.flag
		else
			> convert_finished_chr${chr}.flag
		fi

	else
		> convert_error_chr${chr}.flag
	fi

fi
