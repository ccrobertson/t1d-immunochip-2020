#!/bin/bash

#run from ~/prelim_QC/RELATEDQC directory

famid=$1


echo -e "PEDIGREE"
grep -w ${famid} relatedqc3.fam 

echo -e " "
echo -e "KINSHIP"
grep -w ${famid} relatedqc3.kin
grep -w ${famid} relatedqc3.kin0 

echo -e " "
echo -e "MASTER KINSHIP"
grep -w ${famid} relatedqc3.MKIN

echo -e " "
echo -e "SEX ERRORS"
grep -w ${famid} ../QC_lists/SEXQC1/SEXCHECK2/non_missing_sex_error.txt


#poplus	FIN-1790	FIN-1790.2	FIN-1790	FIN-1790.4
#poplus	FIN-1790	FIN-1790.4	FIN-1790	FIN-1790.2
