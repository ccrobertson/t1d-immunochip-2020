set -euxo pipefail

awk -v file1=${data_derived}/combined_pheno_updated3.txt 'BEGIN { while( getline <file1 ) a[$3]=$1} {print $1, $2, a[$2]$1, $2}'  ${data_derived}/prelimqcC.fam > ${data_derived}/update-ids_cohort.txt 
plink1.9 --bfile ${data_derived}/prelimqcC --update-ids ${data_derived}/update-ids_cohort.txt --make-bed  --out ${data_derived}/prelimqc_FINAL

