#!/bin/bash

set -o errexit
set -o pipefail

date

###############################################################################
# Usage: sh  relatedqc.sh {inputfile} {phenofile} {scripts_dir} {relqc_dir} {qclists_dir}
##############################################################################

plinkraw=$1
prefix=relatedqc


###############################################################################
echo -e
echo -e "\n------------------"
echo "Retrieving data: copying files to ${relqc}"
echo -e "------------------\n"
echo -e


cp ${plinkraw}.bim ${relqc}/${prefix}1.bim
cp ${plinkraw}.bed ${relqc}/${prefix}1.bed
cp ${plinkraw}.fam ${relqc}/${prefix}1.fam
cp ${plinkraw}.log ${relqc}/${prefix}1.log


cd ${relqc}


###############################################################################
echo -e
echo -e "\n------------------"
echo "Setting up: defining functions"
echo -e "------------------\n"
echo -e

#This file contains
#FID1 IID1 --> family and individual ids for first person
#FID2 IID2 --> family and individual ids for second person
#Phi --> kinship coefficient predicted by pedigree file
#PEDREL --> relationship predicated by pedigree file
#SNPREL --> relationship inferred from SNP data (by king --related)
#PEDCHILD --> which subject is the child when the pedigree predicts parent-offspring relationship (otherwise missing)
createMasterKinship () {
	local prefix=$1
	local pheno=$2
	
	echo -e "Creating master kinship from ${prefix}.fam, ${prefix}.kin, and ${prefix}.kin0"
	echo -e " " 
	n_snps=$(grep "pass filters and QC" ${prefix}.log | awk '{print $1}')
	echo -e "${n_snps} SNPs used in kinship analysis"
	n_snp_thresh=$(echo "${n_snps} * 0.75" | bc)
	echo -e "${n_snp_thresh} SNPs required to trust SNP-based relationship inference between expected unrelateds"

	echo -e "FID1\tIID1\tFID2\tIID2\tPhi\tSNPREL" > ${prefix}.MKIN_tmp
	awk -v n_snp_thresh=${n_snp_thresh} 'NR>1 {print $1, $2, $1, $3, $6, $15}' ${prefix}.kin >> ${prefix}.MKIN_tmp
	awk -v n_snp_thresh=${n_snp_thresh} '$5>n_snp_thresh && NR>1 {print $1, $2, $3, $4, "UN", $14}' ${prefix}.kin0 >> ${prefix}.MKIN_tmp
	
	#note, we must use the updated fam files from each iteration to get parent info 
	perl ${scripts}/extract_ped_info.pl ${prefix}.MKIN_tmp ${prefix}.fam ${pheno} > ${prefix}.MKIN
	
	echo -e "$(awk '$6!=$7' ${prefix}.MKIN | wc -l) errors detected in ${prefix}.MKIN"
	echo -e " " 

}


#####################################################################################################################################
echo -e
echo -e "\n------------------------------------------------------------------------------------------------------------"
echo "STEP 1: prelim-SNP/Sample QC -- remove missing non-Y SNPs and all monomorphic SNPs + remove sex_anomalies.txt"
echo -e "------------------------------------------------------------------------------------------------------------\n"
echo -e
#####################################################################################################################################
		
#run king snp qc
king213 -b ${prefix}1.bed --bySNP --prefix ${prefix}1

awk '$11 < 0.95 && $2 != "Y"' ${prefix}1bySNP.txt | awk '{print $1, "CallRateLessThan95"}' > snptoberemoved.txt
awk '$8+$9==0 || $9+$10==0' ${prefix}1bySNP.txt  | awk 'NR>1 {print $1, "Monomorphic"}'  >> snptoberemoved.txt
awk '{print $1, "PoorYsnp"}' ${qclists}/SEXQC1/poorYSNP_remove.txt >> snptoberemoved.txt
cat ${qclists}/SEXQC1/SEXCHECK1/sex_anomalies.txt ${qclists}/SEXQC1/SEXCHECK2/sex_anomalies.txt > sex_anomalies_to_remove.txt
plink1.9 --bfile ${prefix}1 --exclude snptoberemoved.txt --remove sex_anomalies_to_remove.txt --make-bed --out ${prefix}2



#####################################################################################################################################
echo -e
echo -e "\n------------------------------------------------------------------------------------------------------------"
echo "STEP 2: Within-cohort relationship QC"
echo -e "------------------------------------------------------------------------------------------------------------\n"
echo -e
#####################################################################################################################################


##########################################
echo -e
echo -e "\n------------------"
echo "Round 1 (obvious duplicates)"
echo -e "------------------\n"
echo -e

#run king relationship qc						 
king213 -b ${prefix}2.bed --related --prefix ${prefix}2 --degree 2 --rplot
createMasterKinship ${prefix}2 ${phenofile}

#extract within cohort info
awk '$10==$11' ${prefix}2.MKIN > ${prefix}2_withinCohort.MKIN
echo -e "Within cohort errors:"
awk 'NR>1 && $6!=$7 {print $6":"$7}' ${prefix}2_withinCohort.MKIN | sort | uniq -c

#obvious duplicates
awk '$6=="UN" && $7=="DUP/MZ"' ${prefix}2_withinCohort.MKIN | awk '{gsub(/\_2$/,"",$2); gsub(/\_2$/,"",$4); print $1, $2, $3, $4}' | awk '$2==$4 {print $1, $2, $3, $4"_2"}' > duplicates.txt
awk '$6=="UN" && $7=="DUP/MZ"' ${prefix}2_withinCohort.MKIN | awk '{gsub(/\-2$/,"",$2); gsub(/\-2$/,"",$4); print $1, $2, $3, $4}' | awk '$2==$4 {print $1, $2, $3, $4"-2"}' >> duplicates.txt

#obvious duplicate drop list
awk '{print $3, $4}' duplicates.txt > duplicates_drop.txt
#remove obvious duplicates
plink1.9 --bfile ${prefix}2 --remove  duplicates_drop.txt --make-bed --out ${prefix}3
																																																																																							

##########################################
echo -e
echo -e "\n------------------"
echo "Round 2 (within family swaps)"
echo -e "------------------\n"
echo -e



	
#run king relationship qc
king213 -b ${prefix}3.bed --related --prefix ${prefix}3 --degree 2 --rplot
createMasterKinship ${prefix}3 ${phenofile}

#extract within family info
awk '$1==$3' ${prefix}3.MKIN > ${prefix}3_withinFamily.MKIN

#find within family swaps 
#awk 'NR==1 || $3>0' ${qclists}/SEXQC1/SEXCHECK2/sexcheck2_sexupdate.txt > ${qclists}/SEXQC1/SEXCHECK2/sexcheck2_sexupdate_nonmissing.txt
#perl ${scripts}/within_family_swaps3.pl ${prefix}3_withinFamily.MKIN ${qclists}/SEXQC1/SEXCHECK2/sexcheck2_sexupdate_nonmissing.txt ${prefix}3.fam
perl ${scripts}/within_family_swaps3.pl ${prefix}3_withinFamily.MKIN ${qclists}/SEXQC1/SEXCHECK2/sexcheck2_sexupdate.txt ${prefix}3.fam

echo -e "Breakdown in family categories"
cat family_categories.txt | awk '{print $1}'| sort | uniq -c
echo -e " "

#make sure each update file is unique 
sort within_family_swaps_sexupdate.txt | uniq > within_family_swaps_sexupdate2.txt
sort within_family_swaps_parentupdate.txt | uniq > within_family_swaps_parentupdate2.txt
	
#create update file for parent swaps
awk '{print $2, $3, $4"XXX", $5"XXX"}' within_family_parent_swaps.txt > within_family_parent_swaps_update1
awk '{print $4"XXX", $5"XXX", $4, $5}' within_family_parent_swaps.txt > within_family_parent_swaps_update2

#create update file for sibling swaps
awk '{print $2, $3, $4"XXX", $5"XXX"}' within_family_sibling_swaps.txt > within_family_sibling_swaps_update1
awk '{print $4"XXX", $5"XXX", $4, $5}' within_family_sibling_swaps.txt > within_family_sibling_swaps_update2

#create update file for parent-child swaps
awk '{print $2, $3, $4"XXX", $5"XXX"}' within_family_parentchild_swaps.txt > within_family_parentchild_swaps_update1
awk '{print $4"XXX", $5"XXX", $4, $5}' within_family_parentchild_swaps.txt > within_family_parentchild_swaps_update2

#update plink files to fix iid swaps
plink1.9 --bfile ${prefix}3 --update-ids within_family_parent_swaps_update1 --make-bed --out ${prefix}4a
plink1.9 --bfile ${prefix}4a --update-ids within_family_parent_swaps_update2 --make-bed --out ${prefix}4b

plink1.9 --bfile ${prefix}4b --update-ids within_family_sibling_swaps_update1 --make-bed --out ${prefix}4c
plink1.9 --bfile ${prefix}4c --update-ids within_family_sibling_swaps_update2 --make-bed --out ${prefix}4d

plink1.9 --bfile ${prefix}4d --update-ids within_family_parentchild_swaps_update1 --make-bed --out ${prefix}4e
plink1.9 --bfile ${prefix}4e --update-ids within_family_parentchild_swaps_update2 --make-bed --out ${prefix}4f

plink1.9 --bfile ${prefix}4f --update-parents within_family_swaps_parentupdate2.txt --make-bed --out ${prefix}4g
plink1.9 --bfile ${prefix}4g --update-sex within_family_swaps_sexupdate2.txt --make-bed --out ${prefix}4



##########################################
echo -e
echo -e "\n------------------"
echo "Round 3 (between family swaps)"
echo -e "------------------\n"
echo -e


#run king relationship qc	
king213 -b ${prefix}4.bed --related --prefix ${prefix}4 --degree 2 --rplot
createMasterKinship ${prefix}4 ${phenofile}
 
#extract within cohort info
awk '$10==$11' ${prefix}4.MKIN > ${prefix}4_withinCohort.MKIN

#find families with unexpected relatedness
awk '$1!=$3 && ($7=="PO" || $7=="FS" || $7=="DUP/MZ") {print $1, $3}' ${prefix}4_withinCohort.MKIN > family_swaps 
awk -v file=family_swaps 'BEGIN {while (getline<file) {a[$1]=1; a[$2]=1} } ($1 in a && $3 in a)' ${prefix}4_withinCohort.MKIN > ${prefix}4_familySwaps.MKIN

#find between family swaps
#awk 'NR==1 || $3>0' ${qclists}/SEXQC1/SEXCHECK2/sexcheck2_sexupdate.txt > ${qclists}/SEXQC1/SEXCHECK2/sexcheck2_sexupdate_nonmissing.txt
#perl ${scripts}/between_family_swaps3.pl ${prefix}4_familySwaps.MKIN family_swaps ${qclists}/SEXQC1/SEXCHECK2/sexcheck2_sexupdate_nonmissing.txt ${prefix}4.fam
perl ${scripts}/between_family_swaps3.pl ${prefix}4_familySwaps.MKIN family_swaps ${qclists}/SEXQC1/SEXCHECK2/sexcheck2_sexupdate.txt ${prefix}4.fam

#create update file for child swaps
awk '{print $2, $3, $5"XXX", $6"XXX"}' between_family_child_swaps.txt > between_family_child_swaps_update1
awk '{print $5"XXX", $6"XXX", $5, $6}' between_family_child_swaps.txt > between_family_child_swaps_update2

#create update file for parent swaps
awk '{print $2, $3, $5"XXX", $6"XXX"}' between_family_parent_swaps.txt > between_family_parent_swaps_update1
awk '{print $5"XXX", $6"XXX", $5, $6}' between_family_parent_swaps.txt > between_family_parent_swaps_update2

#create update file for parent-child swaps
awk '{print $2, $3, $5"XXX", $6"XXX"}' between_family_parent_child_swaps.txt > between_family_parent_child_swaps_update1 
awk '{print $5"XXX", $6"XXX", $5, $6}' between_family_parent_child_swaps.txt > between_family_parent_child_swaps_update2

#update plink files
plink1.9 --bfile ${prefix}4 --update-ids between_family_child_swaps_update1 --make-bed --out ${prefix}5a
plink1.9 --bfile ${prefix}5a --update-ids between_family_child_swaps_update2 --make-bed --out ${prefix}5b

plink1.9 --bfile ${prefix}5b --update-ids between_family_parent_swaps_update1 --make-bed --out ${prefix}5c
plink1.9 --bfile ${prefix}5c --update-ids between_family_parent_swaps_update2 --make-bed --out ${prefix}5d

plink1.9 --bfile ${prefix}5d --update-ids between_family_parent_child_swaps_update1 --make-bed --out ${prefix}5e
plink1.9 --bfile ${prefix}5e --update-ids between_family_parent_child_swaps_update2 --make-bed --out ${prefix}5f

plink1.9 --bfile ${prefix}5f --update-parents between_family_swaps_parentupdate.txt --make-bed --out ${prefix}5g
plink1.9 --bfile ${prefix}5g --update-sex between_family_swaps_sexupdate.txt --make-bed --out ${prefix}5







##########################################
echo -e
echo -e "\n------------------"
echo "Round 4 (remaining duplicates)"
echo -e "------------------\n"
echo -e



#run king relationship qc
king213 -b ${prefix}5.bed --related --prefix ${prefix}5 --degree 2 --rplot
createMasterKinship ${prefix}5 ${phenofile}
awk '$10==$11' ${prefix}5.MKIN > ${prefix}5_withinCohort.MKIN

#true twins
awk '$1==$3 && $6=="FS" && $7=="DUP/MZ"' ${prefix}5_withinCohort.MKIN | grep -v "_2" > true_twins.txt

#probable duplicates within families
awk 'NR>1 && $1==$3 && $6=="UN" && $7=="DUP/MZ" {print $1"\t"$2"\n"$3"\t"$4}' ${prefix}5_withinCohort.MKIN > probable_duplicates.txt
awk '$1==$3 && $6=="FS" && $7=="DUP/MZ" {print $1"\t"$2"\n"$3"\t"$4}' ${prefix}5_withinCohort.MKIN | grep "_2" >> probable_duplicates.txt
awk '$1==$3 && $6=="PO" && $7=="DUP/MZ" {print $1"\t"$2"\n"$3"\t"$4}' ${prefix}5_withinCohort.MKIN  >> probable_duplicates.txt

#remove probable duplicates
plink1.9 --bfile ${prefix}5 --remove probable_duplicates.txt --make-bed --out ${prefix}6a

#run king relationship qc
king213 -b ${prefix}6a.bed --related --prefix ${prefix}6a --degree 2 --rplot
createMasterKinship ${prefix}6a ${phenofile}
awk '$10==$11' ${prefix}6a.MKIN > ${prefix}6a_withinCohort.MKIN

#probable duplicates between families
awk 'NR>1 && $1!=$3 && $6=="UN" && $7=="DUP/MZ" {print $1, $3}' ${prefix}6a_withinCohort.MKIN > duplicate_families.txt
perl ${scripts}/decide_dup_family_to_drop.pl ${prefix}6a.fam duplicate_families.txt > duplicate_families_drop_list.txt

#remove one family from each pair of duplicate families
plink1.9 --bfile ${prefix}6a --remove duplicate_families_drop_list.txt --make-bed --out ${prefix}6




##########################################
echo -e
echo -e "\n------------------"
echo "Round 5 (wrong parents)" 
echo -e "------------------\n"
echo -e

	
#run king relationship qc
king213 -b ${prefix}6.bed --related --prefix ${prefix}6 --degree 2 --rplot
createMasterKinship ${prefix}6 ${phenofile}
awk '$10==$11' ${prefix}6.MKIN > ${prefix}6_withinCohort.MKIN

#get PO errors (expected parent-offspring pair but were unrelated, or vice versa) 
awk '($6=="UN" && $7=="PO") || ($6=="PO" && $7=="UN")' ${prefix}6_withinCohort.MKIN > po_errors

#flag individuals who are expected to be parents (according to the pedigree) but are neither the parent or child of anyone in the cohort (according to the genotype data)
perl ${scripts}/flag_wrong_parents.pl po_errors > parent_errors

#create parent update files -- this was surprisingly challenging !!!
awk -v file=${prefix}5.fam 'BEGIN {while (getline<file) {fat[$1"\t"$2]=$3; mot[$1"\t"$2]=$4}} {if($3==fat[$1"\t"$2]) {print $1,$2, "0", mot[$1"\t"$2]}}' parent_errors > update_fat
awk -v file=${prefix}5.fam  'BEGIN {while (getline<file) {fat[$1"\t"$2]=$3; mot[$1"\t"$2]=$4}} {if($3==mot[$1"\t"$2]) {print $1,$2, fat[$1"\t"$2], "0"}}' parent_errors > update_mot
awk '{print $1, $2}' parent_errors | sort | uniq > update_merged
perl ${scripts}/merge_parent_updates.pl update_fat update_mot update_merged > parent_errors_update.txt

#update parent status
plink1.9 --bfile ${prefix}6 --update-parents parent_errors_update.txt --make-bed --out ${prefix}7





##########################################
echo -e
echo -e "\n------------------"
echo "Manually check remaining family errors"
echo -e "------------------\n"
echo -e



#run king relationship qc
king213 -b ${prefix}7.bed --related --prefix ${prefix}7 --degree 2 --rplot
createMasterKinship ${prefix}7 ${phenofile}
awk '$10==$11' ${prefix}7.MKIN > ${prefix}7_withinCohort.MKIN


#flag remaining errors for manual inspection and removal
#awk 'NR>1 && $6!=$7 && $6!="UN" {print $0, $6":"$7}' ${prefix}7_withinCohort.MKIN | awk '$12!="FS:DUP/MZ" && $12!="HS:2nd" && $12!="HS:FS"' 
#awk 'NR>1 && $6!=$7 && $6!="UN" {print $0, $6":"$7}' ${prefix}7_withinCohort.MKIN | awk '$12!="FS:DUP/MZ" && $12!="HS:2nd" && $12!="HS:FS" {print $1, $2"\n",$3, $4}' > manual_drop_list.txt
awk 'NR>1 && $6!=$7 && $6!="UN" {print $0, $6":"$7}' ${prefix}7_withinCohort.MKIN | awk '$12!="FS:DUP/MZ" && $12!="HS:2nd"' 
awk 'NR>1 && $6!=$7 && $6!="UN" {print $0, $6":"$7}' ${prefix}7_withinCohort.MKIN | awk '$12!="FS:DUP/MZ" && $12!="HS:2nd" {print $1, $2"\n",$3, $4}' > manual_drop_list.txt


#update plink file
plink1.9 --bfile ${prefix}7 --remove manual_drop_list.txt --make-bed --out ${prefix}8



##########################################
echo -e
echo -e "\n------------------"
echo "Update phenofile to match QC'd data"
echo -e "------------------\n"
echo -e

phenoupdated=combined_pheno_updated.txt

perl ${scripts}/update_phenofile_parents.pl ${phenofile} parent_errors_update.txt > ${phenoupdated}

  

##########################################
echo -e
echo -e "\n------------------"
echo "Select one family from each related pair WITHIN a cohort"
echo -e "------------------\n"
echo -e



#run king relationship qc
king213 -b ${prefix}8.bed --related --prefix ${prefix}8 --degree 2 --rplot
createMasterKinship ${prefix}8 ${phenoupdated}
awk '$10==$11' ${prefix}8.MKIN > ${prefix}8_withinCohort.MKIN

#count errors remaining by type
awk 'NR>1 && $6!=$7 {print $6":"$7}' ${prefix}8_withinCohort.MKIN | sort | uniq -c

#choose one family from each pair to remove
#awk 'NR>1 && $6=="UN" && ($7=="PO" || $7=="FS" || $7=="DUP/MZ") {if (!($1"\t"$3 in a) && !($3"\t"$1 in a)) {a[$1"\t"$3]=1; print $1, $3}}' ${prefix}8_withinCohort.MKIN > unexpected_related_families.txt
awk 'NR>1 && $6=="UN" && ($7=="PO" || $7=="FS" || $7=="DUP/MZ" || $7=="2nd") {if (!($1"\t"$3 in a) && !($3"\t"$1 in a)) {a[$1"\t"$3]=1; print $1, $3}}' ${prefix}8_withinCohort.MKIN > unexpected_related_families.txt
perl ${scripts}/decide_dup_family_to_drop.pl ${prefix}8.fam unexpected_related_families.txt > unexpected_related_families_drop_list.txt 

#remove families from plink files
plink1.9 --bfile ${prefix}8 --remove unexpected_related_families_drop_list.txt --make-bed --out ${prefix}9


#####################################################################################################################################
echo -e
echo -e "\n------------------------------------------------------------------------------------------------------------"
echo "STEP 3: Between-cohort relationship QC"
echo -e "------------------------------------------------------------------------------------------------------------\n"
echo -e
#####################################################################################################################################


##########################################
echo -e
echo -e "\n------------------"
echo "Select one family from each related pair BETWEEN cohorts"
echo -e "------------------\n"
echo -e


#run king relationship qc
king213 -b ${prefix}9.bed --related --prefix ${prefix}9 --degree 2 --rplot
createMasterKinship ${prefix}9 ${phenoupdated}
awk '$10==$11' ${prefix}9.MKIN > ${prefix}9_withinCohort.MKIN
awk '$10!=$11' ${prefix}9.MKIN > ${prefix}9_betweenCohort.MKIN

#count errors remaining by type
echo "Remaining within cohort errors:"
awk 'NR>1 && $6!=$7 {print $6":"$7}' ${prefix}9_withinCohort.MKIN | sort | uniq -c

echo "Between cohort relationships (all are unexpected!):"
awk 'NR>1 && $7!="UN" {print $7}' ${prefix}9_betweenCohort.MKIN | sort | uniq -c

#decide which families to drop
#awk 'NR>1 && $6=="UN" && ($7=="PO" || $7=="FS" || $7=="DUP/MZ" || $7=="2nd") {if (!($1"\t"$3 in a) && !($3"\t"$1 in a)) {a[$1"\t"$3]=1; print $1, $3, $10":"$11}}' ${prefix}9.MKIN > unexpected_related_families2.txt
awk 'NR>1 && ($7=="PO" || $7=="FS" || $7=="DUP/MZ" || $7=="2nd") {if (!($10$1"\t"$11$3 in a) && !($11$3"\t"$10$1 in a)) {a[$10$1"\t"$11$3]=1; print $10$1, $11$3, $10":"$11}}' ${prefix}9_betweenCohort.MKIN > unexpected_related_families2.txt
perl ${scripts}/decide_dup_family_to_drop_between_cohort.pl ${prefix}9.fam unexpected_related_families2.txt ${phenoupdated} > unexpected_related_families_drop_list2.txt 

#remove subjects from dropped families 
plink1.9 --bfile ${prefix}9 --remove unexpected_related_families_drop_list2.txt --make-bed --out ${prefix}_FINAL




##########################################
echo -e
echo -e "\n------------------"
echo "FINAL RELATIONSHIP CHECK"
echo -e "------------------\n"
echo -e



#run king relationship qc
king213 -b ${prefix}_FINAL.bed --related --prefix ${prefix}_FINAL --degree 2 --rplot
createMasterKinship ${prefix}_FINAL ${phenoupdated}
awk '$10==$11' ${prefix}_FINAL.MKIN > ${prefix}_FINAL_withinCohort.MKIN
awk '$10!=$11' ${prefix}_FINAL.MKIN > ${prefix}_FINAL_betweenCohort.MKIN

#count errors remaining by type
echo "Remaining within cohort errors:"
awk 'NR>1 && $6!=$7 {print $6":"$7}' ${prefix}_FINAL_withinCohort.MKIN | sort | uniq -c

echo "Between cohort relationships (all are unexpected!):"
awk 'NR>1 && $7!="UN" {print $6":"$7}' ${prefix}_FINAL_betweenCohort.MKIN | sort | uniq -c


##########################################
echo -e
echo -e "\n------------------"
echo "Copying relationship QC drop and swap lists to ${qclists}"
echo -e "------------------\n"
echo -e

awk '{print $0, "duplicate"}' duplicates_drop.txt > samples_removed.txt
awk '{print $0, "duplicate"}' probable_duplicates.txt  >> samples_removed.txt
awk '{print $0, "duplicate_family"}' duplicate_families_drop_list.txt  >> samples_removed.txt
awk '{print $1, $2, "manual"}' manual_drop_list.txt >> samples_removed.txt
awk '{print $1, $2, "overlapping_family_within_cohort"}' unexpected_related_families_drop_list.txt >> samples_removed.txt
awk '{print $1, $2, "overlapping_family_between_cohort"}' unexpected_related_families_drop_list2.txt >> samples_removed.txt

awk '{print $0, "within_family:parent_swap"}' within_family_parent_swaps.txt > samples_swapped.txt
awk '{print $0, "within_family:sibling_swap"}' within_family_sibling_swaps.txt >> samples_swapped.txt
awk '{print $0, "within_family:parent-child_swap"}' within_family_parentchild_swaps.txt >> samples_swapped.txt
awk '{print $0, "between_family:parent_swap"}' between_family_parent_swaps.txt >> samples_swapped.txt 
awk '{print $0, "between_family:child_swap"}' between_family_child_swaps.txt >> samples_swapped.txt
awk '{print $0, "between_family:parentchild_swap"}' between_family_parent_child_swaps.txt >> samples_swapped.txt

awk '{print $0, "wrong_parent"}' parent_errors_update.txt > pedigree_updates.txt

mkdir -p ${qclists}/RELQC

cp  samples_removed.txt ${qclists}/RELQC
cp  samples_swapped.txt ${qclists}/RELQC
cp  pedigree_updates.txt ${qclists}/RELQC

echo -e "Finished Relationship QC"
