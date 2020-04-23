#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account rich_immunochip_impute
#SBATCH --partition=standard
#SBATCH --mem=50000

chr=${SLURM_ARRAY_TASK_ID}

group=$1

while [ ! -f convert_finished_chr${chr}.flag ]; do
	sleep 120
	echo "SLEEPING"
done

if [ ! -f runtdt_finished_chr${chr}.flag ]; then

	echo -e "STARTING CHR${chr}"
	prefix=chr${chr}.filter_maf_gt_0.005_rsq_gt_0.8.relateds

	echo -e "SWAPPING CASE-CONTROL STATUS"
	awk 'BEGIN {a[1]=2; a[2]=1; a[0]=0; a[-9]=0} {print $1, $2, a[$6]}'  ${prefix}_with_pheno.fam >  swap_cc_status_${chr}.txt
	plink --memory 15000 --bfile ${prefix}_with_pheno --pheno swap_cc_status_${chr}.txt --make-bed --out ${prefix}_WRONG


	echo "RUNNING TDT ON TRIOS"
	plink --memory 15000 --bfile ${prefix}_with_pheno --tdt --ci 0.95 --out ${prefix}_with_pheno
	echo " "

	echo "RUNNING TDT ON NEG CONTROL"
	plink --memory 15000 --bfile ${prefix}_WRONG --tdt --ci 0.95 --out ${prefix}_WRONG
	echo " "

	wait
	echo -e "FINISHED CHR${chr}"

	if [ -f ${prefix}_with_pheno.tdt ]; then

		if [ -f runtdt_error_chr${chr}.flag ]; then
			mv runtdt_error_chr${chr}.flag runtdt_finished_chr${chr}.flag
		else
			> runtdt_finished_chr${chr}.flag
		fi

	else
		> runtdt_error_chr${chr}.flag
	fi

fi
