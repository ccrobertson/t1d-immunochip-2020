outer_project_folder=/data1/ccr5ju/MEGA/pre_qc/Mega_06Aug2018
data_raw=${outer_project_folder}/data_raw
data_derived=${outer_project_folder}/data_derived
clean=${outer_project_folder}/data_clean
prelim=${outer_project_folder}/prelim_QC
qclists=${prelim}/QC_lists
sexqc=${prelim}/SEXQC
relqc=${prelim}/RELATEDQC
snpqc=${prelim}/SNPQC

cd /m/CPHG/MEGA/mega_genotyped/release3

#compile droplist
out=dropped_samples.txt
> ${out}
awk -v file=${data_derived}/sampletoberemoved_preqc_updated_fids.txt '{print $1, $2, "provided_droplists"}' ${data_derived}/sampletoberemoved_preqc.txt >> ${out}
awk '{print $1, $2, "call_rate_lt_80_pct"}' ${qclists}/relaxed_sample_filter.txt >> ${out}   
awk '{print $1, $2, $3}' ${relqc}/sex_anomalies_to_remove.txt >> ${out}
awk '{print $1, $2, $3}' ${qclists}/RELQC/samples_removed.txt  >> ${out}
awk '{print $1, $2, "sex_inconsistency2"}' ${qclists}/SEXQC2/sex_inconsistencies_final.txt >> ${out}
awk '{print $1, $2, "missing_in_pheno_files"}' ${data_derived}/missing_pheno.txt >> ${out}
awk 'NR> 1{print $1, $2, "king_autoqc:"$3}' ${snpqc}/mega_autoQC_sampletoberemoved.txt >> ${out}
awk '{print $3}' dropped_samples.txt | sort | uniq -c

awk 'BEGIN {count[$2]=0} {if (count[$2]>0) {print $1, $2, reason[$2]";"$3}; count[$2]++; reason[$2]=$3}' dropped_samples.txt 


#compile swaplist
awk '{print $2, $3, $5, $6, $7}' ${prelim}/QC_lists/RELQC/samples_swapped.txt > swapped_samples.txt


#other lists
${data_derived}/updateID_2.txt
${data_raw}/T1DGC_SAMPLEQC_SWITCH_LIST_201105.txt

awk -v file=${data_derived}/sampletoberemoved_preqc.txt  'BEGIN {while(getline<file) {a[$2]=$0}} $2 in a {print $1, $2, a[$2]}' /m/CPHG/MEGA/mega_genotyped/release3/mega_b37_clean.fam





cd /data1/ccr5ju/MEGA

cat /m/jdrfdn_scratch/users/projects/IMCHIP/20190315/07132017-ImmunochipPlates760-768_Sample_sheet.csv | awk 'FS="," {print $1}' | awk 'NR>9' > additions_samples.txt

grep -f additions_samples.txt --word-regexp processed_for_assoc/family_analysis_keep_list.txt | wc -l
#12

grep -f additions_samples.txt --word-regexp /m/CPHG/MEGA/cc_analysis/unrelateds_all.txt | wc -l
#370

grep -f additions_samples.txt --word-regexp /m/CPHG/MEGA/mega_genotyped/release3/dropped_samples.txt | awk '{print $3}' | sort | uniq -c
#438


### break down of problem samples INCLUDED in analyses
getCounts () {
	group=$1
	grep --word-regexp -f additions_samples.txt processed_for_imputation_hwe_filter/mega_b37_${group}_cluster.txt > included_${group}  
	fb=$(grep --word-regexp -f  <(awk '{print $2}' included_${group}) processed_for_assoc/family_analysis_keep_list_${group}.txt | wc -l)
	cc=$(grep --word-regexp -f <(awk '{print $2}' included_${group}) /m/CPHG/MEGA/cc_analysis/unrelateds_all.txt | wc -l)
	echo "$group $fb $cc"
}
getCounts "EUR"
getCounts "AMR"
getCounts "AFR"
getCounts "FIN"
getCounts "EAS"

### breakdown of problem samples EXCLUDED from analyses
#awk '{print $2}' raw/release3/mega_b37_clean.fam > passed_qc
#cat /m/CPHG/MEGA/cc_analysis/unrelateds_all.txt processed_for_assoc/family_analysis_keep_list_???.txt | awk 'NR>1 {print $2}' > included_in_assoc
#grep -v --word-regexp -f included_in_assoc  additions_samples.txt > excluded_assoc
#grep --word-regexp -f <(awk '{print $2}' /m/CPHG/MEGA/mega_genotyped/release3/dropped_samples.txt)  additions_samples.txt > excluded_qc
awk -v file1=/m/CPHG/MEGA/mega_genotyped/release3/dropped_samples.txt -v file2=${data_derived}/combined_pheno.txt  'BEGIN {while(getline<file1) {a[$2]=$0}; while(getline<file2) {cohort[$3]=$1}} $1 in a {print cohort[$1],a[$1]}' additions_samples.txt > additions_samples_excluded.txt
#why were they excluded?
awk '{print $4}'  additions_samples_excluded.txt | sort | uniq -c
#what cohorts do sex inconsistencies come from?
awk '$4=="sex_inconsistency2" {print $1}'  additions_samples_excluded.txt | sort | uniq -c

awk -v file=additions_samples.txt 'BEGIN {while(getline<file) {a[$1]=$0}} $2 in a' raw/release3/mega_b37_clean.fam > additions_samples_included.txt
awk -v file=additions_samples_excluded.txt 'BEGIN {while(getline<file) {a[$3]=$0}} $2 in a' raw/release3/mega_b37_clean.fam

#are there any samples in drop list that made it into clean set
awk -v file=/m/CPHG/MEGA/mega_genotyped/release3/dropped_samples.txt 'BEGIN {while(getline<file) {a[$2]=$0}} $2 in a' raw/release3/mega_b37_clean.fam
