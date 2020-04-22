set -euxo pipefail

##############################################################################################
#
#    REPORT QC STATISTICS
#
##############################################################################################

comment "SAMPLE QC SUMMARY STATISTICS" 

STARTSS=$(wc -l <(awk 'NR>1' ${data_derived}/genome_studio_export_map.txt) | awk '{print $1}') 
echo -e "Starting sample size: ${STARTSS}" 
awk 'NR>1 {print $3}' ${data_derived}/genome_studio_export_map.txt | sort | uniq -c
echo -e " "

PREDROPS=$(wc -l <(awk '{print $2}' ${data_derived}/sampletoberemoved_preqc_updated_fids.txt | sort | uniq) | awk '{print $1}') 
echo -e "1 - Dropped based on pre-defined droplists: ${PREDROPS}" 

LOWCALL=$(wc -l ${qclists}/relaxed_sample_filter.txt | awk '{print $1}') 
echo -e "2 - Dropped based on call rate <80%: ${LOWCALL}" 

SEXANOM=$(wc -l ${relqc}/sex_anomalies_to_remove.txt | awk '{print $1}') 
echo -e "3 - Dropped due to sex anomalies: ${SEXANOM}" 

RELDROPS=$(wc -l <(awk '{print $2}' ${qclists}/RELQC/samples_removed.txt | sort | uniq) | awk '{print $1}') 
echo -e "4 - Dropped during relationship QC: ${RELDROPS}" 

SEXINC=$(wc -l ${qclists}/SEXQC2/sex_inconsistencies_final.txt | awk '{print $1}')
echo -e "5 - Dropped due to sex inconsistencies (reported sex inconsistent with observed genotypes): ${SEXINC}" 

NOPHENO=$(wc -l ${data_derived}/missing_pheno.txt | awk '{print $1}')
echo -e "6 - Dropped because subject is missing from phenotype file: ${NOPHENO}" 

AUTODROP=$(wc -l <(awk 'NR>1' ${snpqc}/mega_autoQC_sampletoberemoved.txt) | awk '{print $1}')
echo -e "7 - Dropped in king213 autoQC: ${AUTODROP}"  
echo -e " "

TOTALSAMPLEDROP=$(expr ${PREDROPS} + ${LOWCALL} + ${SEXANOM} + ${RELDROPS} + ${SEXINC} + ${NOPHENO} + ${AUTODROP})
FINALSS=$(wc -l ${clean}/${nickname}_b37_clean.fam | awk '{print $1}') 
echo -e "Total samples dropped:" ${TOTALSAMPLEDROP} 
echo -e "Final sample size: ${FINALSS}" 
echo -e " "

> ${qclists}/sampletoberemoved_all.txt
awk '{print $2, "preqc:"$3}' ${data_derived}/sampletoberemoved_preqc_updated_fids.txt | sort | uniq >> ${qclists}/sampletoberemoved_all.txt
awk '{print $2, "call_rate_lt80"}' ${qclists}/relaxed_sample_filter.txt >> ${qclists}/sampletoberemoved_all.txt
awk '{print $2, "sex_anomalies:"$3}' ${relqc}/sex_anomalies_to_remove.txt >> ${qclists}/sampletoberemoved_all.txt
awk '{print $2, "relqc:"$3}' ${relqc}/samples_removed.txt | sort | uniq >> ${qclists}/sampletoberemoved_all.txt
awk '{print $2, "sex_error"}' ${qclists}/SEXQC2/sex_inconsistencies_final.txt >> ${qclists}/sampletoberemoved_all.txt
awk '{print $2, "missing_pheno"}' ${data_derived}/missing_pheno.txt >> ${qclists}/sampletoberemoved_all.txt
awk 'NR>1{print $2, "autoqc:"$3}' ${snpqc}/mega_autoQC_sampletoberemoved.txt >> ${qclists}/sampletoberemoved_all.txt

> ${qclists}/sampleswapped_all.txt
awk '{print $1, $2, $3, $4, "updateID_6"}' ${data_derived}/updateID_0.txt >> ${qclists}/sampleswapped_all.txt
awk '{print "XXX", $1, "XXX", $2, "T1DGC_SAMPLEQC_SWITCH_LIST_201105"}' ${data_derived}/T1DGC_SAMPLEQC_SWITCH_LIST_201105_B.txt >> ${qclists}/sampleswapped_all.txt
awk '$0~/within_family/ {print $2, $3, $4, $5, "RELQC:"$6} $0~/between_family/ {print $2, $3, $5, $6,"RELQC:"$7}' ${qclists}/RELQC/samples_swapped.txt >> ${qclists}/sampleswapped_all.txt

if [ $(expr ${STARTSS} - ${FINALSS}) -eq ${TOTALSAMPLEDROP} ]; then
	echo " "
else 
	echo "WARNING: Something in Sample stats doesn't add up"
	echo " "
fi

echo -e "For list of all samples removed, see: ${qclists}/sampletoberemoved_all.txt"
echo -e "For list of all sample swaps, see: ${qclists}/sampleswapped_all.txt"


comment "Relationship errors"
echo " "
echo "Pedigree Codes:"
echo "FS means both parents are in pedigree and shared by offspring"
echo "HS means only one parent in pedigree is shared between offspring (so could be FS but second parent id is missing)"
echo "UN means no 1st degree relationship is defined in the pedigree"
echo " "
echo "Inferred Codes are based on IBD sharing, as defined by KING 2.1.3"
echo " "
echo "REMAINING INCONSISTENCIES AFTER QC: (Format is Pedigree:Inferred)"
awk 'NR>1 && $6!=$7 {print $6":"$7}' ${relqc}/relatedqc_FINAL.MKIN | sort | uniq -c
echo " "
echo " "
awk 'NR==1 ||  $6":"$7=="HS:FS"' ${relqc}/relatedqc_FINAL.MKIN > ${relqc}/halfSib_vs_fullSib.txt
echo -e "See ${relqc}/halfSib_vs_fullSib.txt for HS:FS errors"
echo " "



comment "SNP QC SUMMARY STATISTICS" 

STARTSNP=$(wc -l ${data_derived}/${nickname}_main_raw/${nickname}_main0.bim | awk '{print $1}') 
echo -e "Starting number of SNPs: ${STARTSNP}" 
echo -e " "

awk -v file=${data_derived}/${nickname}_main_raw/${nickname}_main1.bim 'BEGIN {while(getline<file) {a[$2]=1}} ($1 in a)' ${data_derived}/mydrop > ${data_derived}/mydrop_nodups
UNMAPPED=$(wc -l ${data_derived}/mydrop_nodups | awk '{print $1}') 
echo -e "1 - Dropped due to non-unique probe mapping: ${UNMAPPED}" 

DUPSNP=$(wc -l ${data_derived}/duplicate_snps_to_remove.txt | awk '{print $1}') 
echo -e "2 - Duplicate SNPs removed: ${DUPSNP}" 

PRELIMSNP=$(wc -l <(awk '{print $1}' ${relqc}/snptoberemoved.txt | sort | uniq) | awk '{print $1}')
echo -e "3 -  Monomorphic and low quality SNPs: ${PRELIMSNP}" 

AUTODROPSNP=$(wc -l <(awk 'NR>1' ${snpqc}/mega_autoQC_snptoberemoved.txt) | awk '{print $1}')
echo -e "4 - SNPs removed in king213 autoQC: ${AUTODROPSNP}" 

MISNP=$(wc -l <(awk 'NR>1' ${snpqc}/MI_error_snps_to_remove.txt) | awk '{print $1}')
echo -e "5 - SNPs removed due to Mendelian Inconsistencies: ${MISNP}" 
echo -e " "

TOTALSNPDROP=$(expr ${UNMAPPED} + ${DUPSNP} + ${PRELIMSNP} + ${AUTODROPSNP} + ${MISNP})
FINALSNP=$(wc -l ${clean}/${nickname}_b37_clean.bim | awk '{print $1}') 
echo -e "Total SNPs dropped: ${TOTALSNPDROP}"
echo -e "Final number of SNPs: ${FINALSNP}" 
echo -e " "

> ${qclists}/snpstoberemoved_all.txt
awk '{print $1, "unmappable"}' ${data_derived}/mydrop_nodups >> ${qclists}/snpstoberemoved_all.txt
awk '{print $1, "duplicate"}' ${data_derived}/duplicate_snps_to_remove.txt >> ${qclists}/snpstoberemoved_all.txt
awk '{print $1, "relqc:"$2}' ${relqc}/snptoberemoved.txt >> ${qclists}/snpstoberemoved_all.txt
awk 'NR>1 {print $1, "autoqc:"$2}' ${snpqc}/mega_autoQC_snptoberemoved.txt >> ${qclists}/snpstoberemoved_all.txt
awk 'NR==1 {for(i=1; i<=NF; i++){H[i]=$i}} NR>1 {out[1]="MendelError(see_third_column)";  for(i=2; i<=NF; i++) {out[i]=H[i]":"$i}; print $1, out[1], "{"out[2]","out[3]","out[4]","out[5]","out[6]"}"}' ${snpqc}/MI_error_snps_to_remove.txt >> ${qclists}/snpstoberemoved_all.txt

if [ $(expr ${STARTSNP} - ${FINALSNP}) -eq ${TOTALSNPDROP} ]; then
	echo " "
else 
	echo "WARNING: Something in SNP stats doesn't add up"
	echo " "
fi

echo -e "For list of all SNPs removed, see: ${qclists}/snpstoberemoved_all.txt"








comment "Track number of samples and snps through pipeline"

numSampleAndSNP () {
	local file=$1
	numSNP=$(wc -l $file.bim | awk '{print $1}')
	numSample=$(wc -l $file.fam | awk '{print $1}')
	echo -e "SNPs: ${numSNP}"
	echo -e "Samples: ${numSample}" 
	echo " "
}


echo -e "Start"
numSampleAndSNP ${data_derived}/${nickname}_main_raw/${nickname}_main0

echo -e "Update dup sample ids"
numSampleAndSNP ${data_derived}/${nickname}_main_raw/${nickname}_main1

echo -e "Drop unmappable probes"
numSampleAndSNP ${data_derived}/${nickname}_main_raw/${nickname}_main2

echo -e "Merge plink files"
numSampleAndSNP ${data_derived}/${nickname}_raw_nophenoA

echo -e "Remove duplicate snps"
numSampleAndSNP ${data_derived}/${nickname}_raw_nophenoC

echo -e "Add parents and phenotypes"
numSampleAndSNP ${data_derived}/${nickname}_rawA

echo -e "Remove samples based on predefined lists"
numSampleAndSNP ${data_derived}/${nickname}_rawD

echo -e "Remove samples with call rate <80%"
numSampleAndSNP ${data_derived}/${nickname}_rawE

echo -e "Relqc"
numSampleAndSNP ${relqc}/relatedqc_FINAL

echo -e "Sexqc + other updates"
numSampleAndSNP ${data_derived}/prelimqc_FINAL

echo -e "Autoqc"
numSampleAndSNP ${data_derived}/prelimqc_FINAL2

echo -e "Mendel errors (snp filter)"
numSampleAndSNP ${clean}/${nickname}_b37_clean
