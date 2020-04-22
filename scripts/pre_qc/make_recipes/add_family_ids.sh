#!/bin/bash

set -euxo pipefail

##############################################################################################
#
#    Add family ids
#
##############################################################################################	
comment "Add family ids" 

	

#update famids 
awk -v file=${phenofile}  'BEGIN { while(getline<file) {fid[$3]=$2;} } {if ($2 in fid) {print $1,$2,fid[$2],$2}}' ${data_derived}/${nickname}_raw_nophenoE.fam > ${data_derived}/update_famids.txt
plink1.9 --bfile ${data_derived}/${nickname}_raw_nophenoE --update-ids ${data_derived}/update_famids.txt --make-bed --out ${data_derived}/${nickname}_raw_nophenoF 


