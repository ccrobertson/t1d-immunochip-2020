set -euxo pipefail

##############################################################################################
#
#    Add parents and phenotypes
#
##############################################################################################	
comment "Add parents and phenotypes (sex and T1D status)" 

	
#update parents
awk -v file=${phenofile}  'BEGIN { while(getline<file) {a[$2"\t"$3]=$4"\t"$5;} } {if ($1"\t"$2 in a) {print $1,$2,a[$1"\t"$2]}}' ${data_derived}/${nickname}_raw_nophenoI.fam > ${data_derived}/update_parents.txt

#update sex (1=male; 2=female; other=unknown)
awk -v file=${phenofile}  'BEGIN { while(getline<file) {a[$2"\t"$3]=$6} } {if ($1"\t"$2 in a) {print $1,$2,a[$1"\t"$2]}}' ${data_derived}/${nickname}_raw_nophenoI.fam > ${data_derived}/update_sex.txt

#update t1d status
awk -v file=${phenofile}  'BEGIN { while(getline<file) {a[$2"\t"$3]=$7} } {if ($1"\t"$2 in a) {print $1,$2,a[$1"\t"$2]}}' ${data_derived}/${nickname}_raw_nophenoI.fam > ${data_derived}/update_t1d.txt

plink1.9 --bfile ${data_derived}/${nickname}_raw_nophenoI --update-parents ${data_derived}/update_parents.txt --make-bed --out ${data_derived}/${nickname}_raw_nophenoJ 
plink1.9 --bfile ${data_derived}/${nickname}_raw_nophenoJ --update-sex ${data_derived}/update_sex.txt --make-bed --out ${data_derived}/${nickname}_raw_nophenoK 
plink1.9 --bfile ${data_derived}/${nickname}_raw_nophenoK  --pheno ${data_derived}/update_t1d.txt --make-bed --out ${data_derived}/${nickname}_rawA 



