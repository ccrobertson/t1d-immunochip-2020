#!/bin/bash
set -euxo pipefail


password1=$1
password2=$2
password3=$3
mafthresh=$4
rsqthresh=$5

WORK=${results}

mkdir -p ${WORK}/EUR ${WORK}/EUR1 ${WORK}/EUR2 ${WORK}/EUR3


echo -e "1. DOWNLOADING"
cd ${WORK}/EUR1
#bash download_EUR1.sh
cd ${WORK}/EUR2
#bash download_EUR2.sh
cd ${WORK}/EUR3
#bash download_EUR3.sh


cd ${WORK}
echo -e "2. UNZIPPING"
#bash ${scripts}/post_imputation/unzipfiles2.sh EUR1 "${password1}" > EUR1/unzip_EUR1.log 2>&1
#bash ${scripts}/post_imputation/unzipfiles2.sh EUR2 "${password2}" > EUR2/unzip_EUR2.log 2>&1
#bash ${scripts}/post_imputation/unzipfiles2.sh EUR3 "${password3}" > EUR3/unzip_EUR3.log 2>&1



echo -e "3. MERGING"
#bash ${scripts}/post_imputation/merge_EUR_info_files.sh > merge_EUR_info_files.log 2>&1
bash ${scripts}/post_imputation/merge_EUR_batches.sh > merge_EUR_batches_2.log 2>&1


echo -e "4. FILTERING FOR MAF>${mafthresh} AND RSQ>${rsqthresh}"
bash ${scripts}/post_imputation/filter_for_maf_rsq_EUR.sh EUR ${mafthresh} ${rsqthresh} dose

cd ../results_filtered
echo -e "5. CREATING ANALYSIS DATA SETS"
bash ${scripts}/post_imputation/create_analysis_vcfs.sh EUR filter_maf_gt_${mafthresh}_rsq_gt_${rsqthresh}


#Combine info files across chromosomes
cd ${WORK}/EUR
head -n 1 chr1.merged_info > chrALL.merged_info
for i in {1..22}; do
  echo "chr${i}"
  awk 'NR>1' chr${i}.merged_info >> chrALL.merged_info
done
