set -euxo pipefail

##############################################################################################
#
#    PRE-QC: Remove samples
#
##############################################################################################
#
#	Here we are removing samples that (a priori) we do not want to include in 
#	analysis -- either because we do not have permission to use them or because the 
#	study dropped them from their cohorts (could be for a number of reasons)
#

comment "Remove samples based on pre-qc drop lists" 


### remove TEDDY samples
grep "TEDDY" ${phenofile} | awk '{print $2, $3, $1}' > ${data_derived}/sampletoberemoved_TEDDY.txt
cat  ${data_derived}/sampletoberemoved_TEDDY.txt > ${data_derived}/sampletoberemoved_preqc.txt

### remove T1DGC one-way switches
awk '{print $1, $2, "T1DGC_one_way_switches"}' ${data_derived}/t1dgc_switches_to_drop.txt  >> ${data_derived}/sampletoberemoved_preqc.txt

### remove nPOD samples
grep "nPOD" ${phenofile} | awk '{print $2, $3, $1}' > ${data_derived}/sampletoberemoved_nPOD.txt
cat  ${data_derived}/sampletoberemoved_nPOD.txt >> ${data_derived}/sampletoberemoved_preqc.txt

### remove Existing_Source_to_T1DGC:She samples
grep "Existing_Source_to_T1DGC:She" ${phenofile} | awk '{print $2, $3, $1}' > ${data_derived}/sampletoberemoved_OTHER.txt
cat  ${data_derived}/sampletoberemoved_OTHER.txt >> ${data_derived}/sampletoberemoved_preqc.txt

### remove DK samples flagged during q/c
awk 'NR>1' ${data_raw}/sample_qc_files/sampletoberemoved_DAN.txt >> ${data_derived}/sampletoberemoved_preqc.txt

### remove TrialNet samples flagged during q/c
cat ${data_raw}/sample_qc_files/sampletoberemoved_TrialNet.txt >> ${data_derived}/sampletoberemoved_preqc.txt

### remove samples flagged during first round of T1DGC ImmunoChip_Project
cat ${data_raw}/sample_qc_files/T1DGC_SAMPLEQC_DROP_LIST_201105.tsv >> ${data_derived}/sampletoberemoved_preqc.txt

### remove Denver samples that are actually T2D subjects
cat ${data_raw}/sample_qc_files/Denver_drop_list.txt >> ${data_derived}/sampletoberemoved_preqc.txt

### update famids in preqc droplist to match pedigree file
#awk -v file=${data_derived}/${nickname}_rawA.fam 'BEGIN {while(getline<file) {a[$2]=$1}} ($2 in a) && !($1==a[$2]) {print $1, a[$2], $2, $3}' ${data_derived}/sampletoberemoved_preqc.txt
awk -v file=${data_derived}/${nickname}_rawA.fam 'BEGIN {while(getline<file) {a[$2]=$1}} $2 in a {print a[$2], $2, $3}' ${data_derived}/sampletoberemoved_preqc.txt | sort | uniq  > ${data_derived}/sampletoberemoved_preqc_updated_fids.txt

### remove all preqc samples using sampletoberemoved_preqc.txt
plink1.9 --bfile ${data_derived}/${nickname}_rawA --remove ${data_derived}/sampletoberemoved_preqc_updated_fids.txt --make-bed  --out ${data_derived}/${nickname}_rawD 




