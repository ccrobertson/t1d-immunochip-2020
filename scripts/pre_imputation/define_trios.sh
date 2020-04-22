#!/bin/bash

set -euxo pipefail

nickname=${project}_${release_build}

comment "Determining trios for TDT analysis"
date

### Define trios for tdt (note, we need to do this on a batch basis because a few families have been broken into different batches)
perl ${scripts}/pre_imputation/define_trios_by_batch.pl FIN ${phenofile} ${proc_for_imp}/${nickname}_FIN_cluster.txt > ${fam_assoc}/family_analysis_keep_list_FIN.txt
perl ${scripts}/pre_imputation/define_trios_by_batch.pl EUR ${phenofile} ${proc_for_imp}/${nickname}_EUR_cluster.txt > ${fam_assoc}/family_analysis_keep_list_EUR.txt
perl ${scripts}/pre_imputation/define_trios_by_batch.pl AFR ${phenofile} ${proc_for_imp}/${nickname}_AFR_cluster.txt > ${fam_assoc}/family_analysis_keep_list_AFR.txt
perl ${scripts}/pre_imputation/define_trios_by_batch.pl EAS ${phenofile} ${proc_for_imp}/${nickname}_EAS_cluster.txt > ${fam_assoc}/family_analysis_keep_list_EAS.txt
perl ${scripts}/pre_imputation/define_trios_by_batch.pl AMR ${phenofile} ${proc_for_imp}/${nickname}_AMR_cluster.txt > ${fam_assoc}/family_analysis_keep_list_AMR.txt

cat ${fam_assoc}/family_analysis_keep_list_???.txt > ${fam_assoc}/family_analysis_keep_list.txt
awk '{print $2}' ${fam_assoc}/family_analysis_keep_list.txt > ${fam_assoc}/family_analysis_keep_list_iids.txt

#get affected offspring ids
awk 'BEGIN {print "group","famid","iid","cohort"}' > ${fam_assoc}/affected_offspring.txt
perl ${scripts}/pre_imputation/define_trios_count_affecteds_by_batch.pl FIN ${phenofile} ${proc_for_imp}/${nickname}_FIN_cluster.txt | awk '{print "FIN",$0}' >> ${fam_assoc}/affected_offspring.txt
perl ${scripts}/pre_imputation/define_trios_count_affecteds_by_batch.pl EUR ${phenofile} ${proc_for_imp}/${nickname}_EUR_cluster.txt | awk '{print "EUR",$0}' >> ${fam_assoc}/affected_offspring.txt
perl ${scripts}/pre_imputation/define_trios_count_affecteds_by_batch.pl AFR ${phenofile} ${proc_for_imp}/${nickname}_AFR_cluster.txt | awk '{print "AFR",$0}' >> ${fam_assoc}/affected_offspring.txt
perl ${scripts}/pre_imputation/define_trios_count_affecteds_by_batch.pl EAS ${phenofile} ${proc_for_imp}/${nickname}_EAS_cluster.txt | awk '{print "EAS",$0}' >> ${fam_assoc}/affected_offspring.txt
perl ${scripts}/pre_imputation/define_trios_count_affecteds_by_batch.pl AMR ${phenofile} ${proc_for_imp}/${nickname}_AMR_cluster.txt | awk '{print "AMR",$0}' >> ${fam_assoc}/affected_offspring.txt

#count affected trios by ancestry batch
awk 'BEGIN {print "group","affected_trios"}' > ${fam_assoc}/trio_count.txt
wc -l ${fam_assoc}/affected_offspring.txt | awk '{print $1, "all"}' >> ${fam_assoc}/trio_count.txt
awk 'NR>1 {print $1}' ${fam_assoc}/affected_offspring.txt | sort | uniq -c >> ${fam_assoc}/trio_count.txt

cat ${fam_assoc}/trio_count.txt
