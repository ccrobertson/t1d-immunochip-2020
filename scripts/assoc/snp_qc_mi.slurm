#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account rich_immunochip_impute
#SBATCH --partition=standard
#SBATCH --mem=50000


group=$1
cd ${fam_assoc}/${group}

while [ ! -f convert_finished_chr2.flag ]; do
	echo "SLEEPING 1"
	sleep 120
done

while [ $(ls -1 convert_finished_chr*.flag | wc -l) -lt 22 ]; do
	echo "SLEEPING 2"
	sleep 120
done



if [ ! -f king_snpqcbySNP.txt ]; then
	echo -e "START ${group}"

	#Merge chromosomes
	echo "MERGING PLINK FILES"
	> plinkfiles.txt
	for i in {1..22}; do
		echo "chr${i}.filter_maf_gt_0.005_rsq_gt_0.8.relateds_with_pheno" >> plinkfiles.txt
	done
	plink --memory 50000 --merge-list plinkfiles.txt --make-bed --out chrALL.filter_maf_gt_0.005_rsq_gt_0.8.relateds_with_pheno
	echo " "

	#SNP and SAMPLE QC in KING
	echo "RUNNING KING SNP QC"
	king -b chrALL.filter_maf_gt_0.005_rsq_gt_0.8.relateds_with_pheno.bed --cluster --bySNP --prefix king_snpqc
	echo " "

	#Filter out Err_InHomPO > 0.01 OR Err_InHetTrio > 0.01
	awk 'NR==1{for (i=1; i<=NF; i++) if ($i~/Err_InPO|Err_InTrio|Err_InHomPO|Err_InHetTrio/) {a[$i]=i}} {print $1, $a["Err_InPO"], $a["Err_InHomPO"], $a["Err_InTrio"], $a["Err_InHetTrio"]}' king_snpqcbySNP.txt | awk '$3>0.01 || $5>0.01 {print $1}' > snps_with_mendelInconsistencies_gt1pct.txt
	awk 'NR==1{for (i=1; i<=NF; i++) if ($i~/Err_InPO|Err_InTrio|Err_InHomPO|Err_InHetTrio/) {a[$i]=i}} {print $1, $a["Err_InPO"], $a["Err_InHomPO"], $a["Err_InTrio"], $a["Err_InHetTrio"]}' king_snpqcbySNP.txt | awk '$3>0.05 || $5>0.05 {print $1}' > snps_with_mendelInconsistencies_gt5pct.txt

	#Break into chromosome
	for ((i=1; i<=22; i++)); do
	awk -v chr=${i} 'NR==1 || $2==chr' king_snpqcbySNP.txt > king_snpqcbySNP_chr${i}.txt
	done

	echo -e "FINISHED ${group}"
fi
