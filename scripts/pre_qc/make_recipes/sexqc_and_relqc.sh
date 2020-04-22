set -euxo pipefail

##############################################################################################
#
#    Sex and Relationship QC
#
##############################################################################################
comment "Sex and Relationship QC" 


cd ${prelim}

#preliminary sex QC (SEXQC1)
echo -e "See ${sexqc}/sexqc1.log"  
bash ${scripts}/sexqc2.sh ${data_derived}/${nickname}_rawE SEXQC1 > ${sexqc}/sexqc1.log 2>&1

#relationship QC
echo -e "See ${relqc}/relatqc2.log" 
bash ${scripts}/relatedqc2.sh ${data_derived}/${nickname}_rawE > ${relqc}/relatqc2.log 2>&1


#copy relationship QC'd data
cp ${relqc}/relatedqc_FINAL.bed ${data_derived}/${nickname}_rawF.bed
cp ${relqc}/relatedqc_FINAL.bim ${data_derived}/${nickname}_rawF.bim
cp ${relqc}/relatedqc_FINAL.fam ${data_derived}/${nickname}_rawF.fam
cp ${relqc}/combined_pheno_updated.txt ${data_derived}/combined_pheno_updated.txt

#second sex QC (SEXQC2)
echo -e "See ${sexqc}/sexqc2.log" 
bash ${scripts}/sexqc2.sh ${data_derived}/${nickname}_rawF SEXQC2 > ${sexqc}/sexqc2.log 2>&1
wc -l ${sexqc}/SEXQC*/SEXCHECK2/non_missing_sex_error.txt 

#remove subjects with sex inconsistencies
awk '{print $1, $2}' ${qclists}/SEXQC2/SEXCHECK2/non_missing_sex_error.txt > ${qclists}/SEXQC2/sex_inconsistencies_final.txt
plink1.9 --bfile ${data_derived}/${nickname}_rawF --remove ${qclists}/SEXQC2/sex_inconsistencies_final.txt --make-bed --out ${data_derived}/${nickname}_rawG 

#infer sex from genotypes in when sex is missing
awk '{print $1, $2, $4}' ${qclists}/SEXQC2/SEXCHECK2/missing_sex_imputed.txt >  ${qclists}/SEXQC2/update_sex_final.txt
plink1.9 --bfile ${data_derived}/${nickname}_rawG --update-sex ${qclists}/SEXQC2/update_sex_final.txt --make-bed --out ${data_derived}/${nickname}_rawH 

#remove subjects not in phenofile
awk -v file1=${data_derived}/combined_pheno_updated.txt 'BEGIN { while( getline <file1 ) a[$3]=1} !($2 in a) {print $1, $2}'  ${data_derived}/${nickname}_rawH.fam > ${data_derived}/missing_pheno.txt
plink1.9 --bfile ${data_derived}/${nickname}_rawH --remove ${data_derived}/missing_pheno.txt --make-bed --out ${data_derived}/${nickname}_rawI  

#update phenotype file: extract phenofile for only genotyped subjects passing QC AND update sex
perl ${scripts}/update_phenofile_sex.pl ${data_derived}/combined_pheno_updated.txt ${qclists}/SEXQC2/update_sex_final.txt > ${data_derived}/combined_pheno_updated2.txt
awk -v file1=${data_derived}/${nickname}_rawI.fam 'BEGIN { while( getline <file1 ) a[$2]=1} $3 in a' ${data_derived}/combined_pheno_updated2.txt > ${data_derived}/combined_pheno_updated3.txt

#add cohort info to family ids
awk -v file1=${data_derived}/combined_pheno_updated3.txt 'BEGIN { while( getline <file1 ) a[$3]=$1} {print $1, $2, a[$2]$1, $2}'  ${data_derived}/${nickname}_rawI.fam > ${data_derived}/update-ids_cohort.txt 
plink1.9 --bfile ${data_derived}/${nickname}_rawI --update-ids ${data_derived}/update-ids_cohort.txt --make-bed  --out ${data_derived}/${nickname}_rawJ 



