#!/bin/bash
set -euxo pipefail

group=$1
cd /m/CPHG/MEGA/IMPUTED_TOPMED/results/${group}/filtered


echo -e "########################################################################################################################################"
echo -e
echo -e "# Usage: bash merge_chromosomes.sh {group}"
echo -e "#"
echo -e "#"
echo -e "#"
echo -e "#"
echo -e "Performed on:"
date
echo -e
echo -e "########################################################################################################################################"


softpath=/h4/t1/users/ccr5ju/software
plink=${softpath}/bin/plink1.9

merged=${group}_filtered_chrAll

#echo -n > merge_list.txt
#for chr in {1..22}; do
	#	echo "chr${chr}_filtered" >> merge_list.txt
#done

ls -l | grep .bim | awk '{print $9}' | sed 's/.bim//g' > merge_list.txt

${plink} --merge-list merge_list.txt --make-bed --out ${merged}
