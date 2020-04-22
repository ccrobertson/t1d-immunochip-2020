set -euxo pipefail

##############################################################################################
#
#    Remove duplicate SNPs
#
##############################################################################################
#	
#	Many snps on ichip are duplicated. In other words, the exact same chr:pos:ref:alt is found
#	multiple times with two different ids (typically a standard rs id and then an ichip-specific id
#	e.g. 1:67303323:A:G is found as both 1kg_1_67541594 and rs12131065
#
#	Some of these were included intentionally for QC purposes (i think)
#	I am taking the one with the lowest missiningness and removing the other
#
#
comment "Remove duplicate SNPs" 

cd ${data_derived}

### Calculate missingness
echo -e " " 
echo -e "Now calculating snp and sample-level missingness in merged file"  
plink1.9 --bfile ${data_derived}/${nickname}_raw_nophenoA --missing --out ${data_derived}/${nickname}_raw_nophenoA 		
	
### Identify duplicate variant ids
echo -e " " 
echo -e "Now finding duplicate snp with higher missingness"  
awk 'BEGIN {a[$1":"$4":"$5":"$6]=0} {if (a[$1":"$4":"$5":"$6]==0) {a[$1":"$4":"$5":"$6]=$2} else {print $1":"$4":"$5":"$6, a[$1":"$4":"$5":"$6], $2}}' ${data_derived}/${nickname}_raw_nophenoA.bim > ${data_derived}/duplicate_snps_chrpos.txt

#for each duplicate take variant with lower missingness
awk -v file=${data_derived}/${nickname}_raw_nophenoA.lmiss 'BEGIN {while ( getline<file ) {a[$2]=$5} } {print $1, $2, a[$2], $3, a[$3]}' ${data_derived}/duplicate_snps_chrpos.txt > ${data_derived}/duplicate_snp_missingness.txt
awk -v file=${data_derived}/${nickname}_raw_nophenoA.lmiss 'BEGIN {while ( getline<file ) {a[$2]=$5} } a[$2]>=a[$3] {print $2} a[$3]>a[$2] {print $3}' ${data_derived}/duplicate_snps_chrpos.txt > ${data_derived}/duplicate_snps_with_higher_missingness.txt
	
#for each pair of duplicate variants, calculate genotype concordance
echo -e " " 
echo -e "Now finding duplicate snp pairs with low concordance"  
awk '{print $2"\n"$3}' ${data_derived}/duplicate_snps_chrpos.txt > ${data_derived}/duplicate_snps_all.txt
plink1.9 --bfile ${data_derived}/${nickname}_raw_nophenoA --extract ${data_derived}/duplicate_snps_all.txt --make-bed --out ${data_derived}/${nickname}_raw_nophenoA_duplicates_only  
plink1.9 --bfile ${data_derived}/${nickname}_raw_nophenoA_duplicates_only --recode A --out ${data_derived}/${nickname}_raw_nophenoA_duplicates_only  
Rscript ${scripts}/duplicate_snp_concordance.R ${data_derived}/${nickname}_raw_nophenoA_duplicates_only.raw ${data_derived}/duplicate_snps_low_concordance.txt 0.95

#duplicate snps with poor genotype concordance (<95%) accross samples will be removed
cat <(sed 's/^X1kg/1kg/g; s/_[ATGC]$//g' ${data_derived}/duplicate_snps_low_concordance.txt) ${data_derived}/duplicate_snps_with_higher_missingness.txt | sort | uniq >  ${data_derived}/duplicate_snps_to_remove.txt

#remove duplicates from plink file
echo -e " " 
echo -e "Now removing duplicate with higher missingness or duplicate pairs with low concordance."  
plink1.9 --bfile ${data_derived}/${nickname}_raw_nophenoA --exclude ${data_derived}/duplicate_snps_to_remove.txt  --out ${data_derived}/${nickname}_raw_nophenoC --make-bed 




