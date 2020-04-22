set -euxo pipefail

##############################################################################################
#
#    Prelim sample QC
#
##############################################################################################
comment "Preliminary sample QC" 



cd ${prelim}

echo "Relaxed sample filter -- remove samples with call rate <80%" 
plink1.9 --bfile ${data_derived}/${nickname}_rawD --missing --out ${data_derived}/${nickname}_rawD 
awk 'NR>1 && $6>0.2 {print $1, $2}' ${data_derived}/${nickname}_rawD.imiss > ${qclists}/relaxed_sample_filter.txt
plink1.9 --bfile ${data_derived}/${nickname}_rawD --remove ${qclists}/relaxed_sample_filter.txt --make-bed --out ${data_derived}/${nickname}_rawE 




