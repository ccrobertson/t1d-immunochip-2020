set -euxo pipefail

#remove subjects with sex inconsistencies
awk '{print $1, $2}' ${qclists}/SEXQC2/SEXCHECK2/non_missing_sex_error.txt > ${qclists}/SEXQC2/sex_inconsistencies_final.txt
plink1.9 --bfile ${relqc}/relatedqc_FINAL --remove ${qclists}/SEXQC2/sex_inconsistencies_final.txt --make-bed --out ${data_derived}/prelimqcA 

#infer sex from genotypes in when sex is missing
awk '{print $1, $2, $4}' ${qclists}/SEXQC2/SEXCHECK2/missing_sex_imputed.txt >  ${qclists}/SEXQC2/update_sex_final.txt
plink1.9 --bfile ${data_derived}/prelimqcA --update-sex ${qclists}/SEXQC2/update_sex_final.txt --make-bed --out ${data_derived}/prelimqcB

#remove subjects not in phenofile
awk -v file1=${relqc}/combined_pheno_updated.txt 'BEGIN { while( getline <file1 ) a[$3]=1} !($2 in a) {print $1, $2}'  ${data_derived}/prelimqcB.fam > ${data_derived}/missing_pheno.txt
plink1.9 --bfile ${data_derived}/prelimqcB --remove ${data_derived}/missing_pheno.txt --make-bed --out ${data_derived}/prelimqcC 

#update phenotype file: extract phenofile for only genotyped subjects passing QC AND update sex
perl ${scripts}/update_phenofile_sex.pl ${relqc}/combined_pheno_updated.txt ${qclists}/SEXQC2/update_sex_final.txt > ${data_derived}/combined_pheno_updated2.txt
awk -v file1=${data_derived}/prelimqcC.fam 'BEGIN { while( getline <file1 ) a[$2]=1} $3 in a' ${data_derived}/combined_pheno_updated2.txt > ${data_derived}/combined_pheno_updated3.txt



