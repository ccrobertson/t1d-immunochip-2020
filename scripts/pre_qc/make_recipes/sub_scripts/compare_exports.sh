#!/bin/bash

set -o errexit
set -o pipefail

cd ${data_derived}

#if [ ]; then
### Compare allele frequencies across exports
getAFbyExport() {
	local export1=$1
	local export2=$2
	label=${export2}_vs_${export1}
	pheno1=pheno_${label}.txt
	awk -v e1=${export1} -v e2=${export2} 'BEGIN {print "FID", "IID", "Export"} NR>1 {if ($3==e1) {print $1, $2, "1"} else if ($3==e2) {print $1, $2, "2"} else {print $1, $2, "-9"}}' genome_studio_export_map.txt > ${pheno1} 
	plink1.9 --bfile mega_raw_nophenoA --pheno ${pheno1} --make-bed --out testing_af_diff_${label}
	plink1.9 --bfile testing_af_diff_${label} --assoc --out testing_af_diff_${label} --allow-no-sex
}
getAFbyExport ${nickname_1} ${nickname_2}
getAFbyExport ${nickname_1} ${nickname_3}
getAFbyExport ${nickname_1} ${nickname_4}


### Compare SNP Mendelian error rates across exports (based on inferred relationships)
#remove family ids from merged data
genodat=mega_rawD
awk '{print $1, $2, $2, $2}' ${genodat}.fam > remove_fids_2.txt
plink1.9 --bfile ${genodat} --update-ids remove_fids_2.txt --make-bed --out ${genodat}_nofamids
getMERatebyExport() {
	local export=$1
	#extract export subjects
	awk -v export=${export} '$3==export' ${data_derived}/genome_studio_export_map.txt > ${export}_keeplist 
	plink1.9 --bfile ${genodat}_nofamids --keep ${export}_keeplist --make-bed --out ${genodat}_nofamids_${export}
	#run snp qc
	king213 -b ${genodat}_nofamids_${export}.bed --cluster --bySNP --prefix testing_mendelError_${export}
}
getMERatebyExport ${nickname_1}
getMERatebyExport ${nickname_2}
getMERatebyExport ${nickname_3}
getMERatebyExport ${nickname_4}


### Compare relationship error rates within exports (how many pedigree errors?)
king213 -b ${genodat}.bed --bysample --prefix testing_MendelErrors


### Relationship errors in families from multiple exports
#find families in multiple exports
perl ${scripts}/check_exports.pl ${phenofile} ${data_derived}/genome_studio_export_map.txt pheno_with_export.txt > families_in_multiple_exports.txt
awk -v file=families_in_multiple_exports.txt 'BEGIN { while(getline<file) a[$1]=1} $2 in a {print $2, $3, $8}' pheno_with_export.txt > families_in_multiple_exports_iids.txt
plink1.9 --bfile ${genodat} --keep families_in_multiple_exports_iids.txt --make-bed --out ${genodat}_multiple_exports
#run relationship qc
king213 -b ${genodat}_multiple_exports.bed --related --prefix testing_RelErrors

#fi 

# GENERATE RMARKDOWN REPORT
cp ${scripts}/export_bias.Rmd .  
R -e "rmarkdown::render('export_bias.Rmd')" 



