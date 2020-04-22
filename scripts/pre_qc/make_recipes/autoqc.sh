set -euxo pipefail

##############################################################################################
#
#    Final SNP QC
#
##############################################################################################
comment "SNP QC" 

cd ${data_derived}

#autoQC		
king213 -b ${data_derived}/prelimqc_FINAL.bed --autoQC --callrateM 0.98 --callrateN 0.98 --prefix ${snpqc}/${nickname}  
plink1.9 --bfile ${data_derived}/prelimqc_FINAL --remove ${snpqc}/${nickname}_autoQC_sampletoberemoved.txt --exclude ${snpqc}/${nickname}_autoQC_snptoberemoved.txt --make-bed --out ${data_derived}/prelimqc_FINAL2

#get rid of fids
awk '{print $1, $2, $2, $2}' ${data_derived}/prelimqc_FINAL2.fam > ${snpqc}/remove_fids.txt
plink1.9 --bfile ${data_derived}/prelimqc_FINAL2 --update-ids ${snpqc}/remove_fids.txt --make-bed --out ${snpqc}/snpqc1 

#remove snps with low concordance between duplicates
#Err_InMZ > 0.01
#Err_InPO > 0.01
#Err_InTrio > 0.01
#Err_InHomPO > 0.1
#Err_InHetTrio > 0.1
king213 -b ${snpqc}/snpqc1.bed --cluster --bySNP --prefix ${snpqc}/snpqc1 
awk '{print $1, $15, $20, $25, $21, $26}' ${snpqc}/snpqc1bySNP.txt | awk '$2>0.01 || $3>0.01 || $4>0.01 || $5>0.1 || $6>0.1' > ${snpqc}/MI_error_snps_to_remove.txt
plink1.9 --bfile ${data_derived}/prelimqc_FINAL2 --exclude ${snpqc}/MI_error_snps_to_remove.txt --make-bed --out ${clean}/${nickname}_b37_clean 

#update pheno file
cp ${data_derived}/combined_pheno_updated3.txt ${clean}/${nickname}_b37_clean_pheno.txt
awk -v file=${clean}/${nickname}_b37_clean.fam 'BEGIN { while ( getline<file ) {a[$1"\t"$2]=1} } ($1$2"\t"$3 in a) {print $0}' ${clean}/${nickname}_b37_clean_pheno.txt | sort | uniq > ${clean}/${nickname}_b37_clean_pheno_genotypedOnly.txt



