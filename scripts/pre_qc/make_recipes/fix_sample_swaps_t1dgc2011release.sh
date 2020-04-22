set -euxo pipefail

##############################################################################################
#
#    Fix sample switches based on May 2011 T1DGC release switch lists
#
##############################################################################################
comment "Fix sample switches based on May 2011 T1DGC release switch lists" 


### remove carriage return
sed 's/\r//g'  ${data_raw}/sample_qc_files/T1DGC_SAMPLEQC_SWITCH_LIST_201105.txt > ${data_derived}/T1DGC_SAMPLEQC_SWITCH_LIST_201105_A.txt

### note there are some sample switches that only map one direction (so not really "switches" but rather mapping one id to another)
### for now will drop them
awk 'NR>1{print $1"\n"$2}' ${data_derived}/T1DGC_SAMPLEQC_SWITCH_LIST_201105_A.txt | LANG=en_EN sort | uniq -u | LANG=en_EN join -1 1 -2 3 - <(LANG=en_EN sort -k3 ${phenofile}) | awk '{print $3, $1}' > ${data_derived}/t1dgc_switches_to_drop.txt 
awk -v file=${data_derived}/t1dgc_switches_to_drop.txt 'BEGIN {while (getline<file) {a[$2]=0;} } NR>1 && !($1 in a || $2 in a)' ${data_derived}/T1DGC_SAMPLEQC_SWITCH_LIST_201105_A.txt > ${data_derived}/T1DGC_SAMPLEQC_SWITCH_LIST_201105_B.txt 

### get family ids for sample switches
LANG=en_EN join -1 3 -2 1 <( LANG=en_EN sort -k3 ${phenofile} ) <( LANG=en_EN sort -k1 ${data_derived}/T1DGC_SAMPLEQC_SWITCH_LIST_201105_B.txt ) | awk '{print $3, $1, $8, $8"X"}' > ${data_derived}/T1DGC_SAMPLEQC_SWITCH_LIST_201105_C.txt
LANG=en_EN join -1 3 -2 2 <( LANG=en_EN sort -k3 ${phenofile} ) <( LANG=en_EN sort -k2 ${data_derived}/T1DGC_SAMPLEQC_SWITCH_LIST_201105_B.txt ) | awk '{print $1, $1"X", $3, $1}' > ${data_derived}/T1DGC_SAMPLEQC_SWITCH_LIST_201105_D.txt

plink1.9 --bfile ${data_derived}/${nickname}_raw_nophenoF --update-ids ${data_derived}/T1DGC_SAMPLEQC_SWITCH_LIST_201105_C.txt --make-bed  --out ${data_derived}/${nickname}_raw_nophenoG 
plink1.9 --bfile ${data_derived}/${nickname}_raw_nophenoG --update-ids ${data_derived}/T1DGC_SAMPLEQC_SWITCH_LIST_201105_D.txt --make-bed  --out ${data_derived}/${nickname}_raw_nophenoI 



