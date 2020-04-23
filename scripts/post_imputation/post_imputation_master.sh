#!/bin/bash
set -euxo pipefail

group=$1
password=$2
mafthresh=$3
rsqthresh=$4

WORK=${results}

mkdir -p ${WORK}/${group}

cd ${WORK}/${group}

if [ ! -f download_${group}.complete ]; then
	echo -e "1. DOWNLOADING"
	bash download_${group}.sh
	> download_${group}.complete
fi

cd ${WORK}
echo -e "2. UNZIPPING"
if [ ! -f ${group}/unzip_${group}.complete ]; then
	bash ${scripts}/post_imputation/unzipfiles2.sh ${group} "${password}" > ${group}/unzip_${group}.log 2>&1
	success=$(grep "Everything is Ok" ${group}/unzip_${group}.log | wc -l)
	if [ ${success} -eq "22" ]; then
		>${group}/unzip_${group}.complete
	else
		>${group}/unzip_${group}.error
	fi
fi

echo -e "3. FILTERING FOR MAF>${mafthresh} AND RSQ>${rsqthresh}"
bash ${scripts}/post_imputation/filter_for_maf_rsq.sh ${group} ${mafthresh} ${rsqthresh} dose


cd ${filtered}
echo -e "4. CREATING ANALYSIS DATA SETS"
bash ${scripts}/post_imputation/create_analysis_vcfs.sh ${group} filter_maf_gt_${mafthresh}_rsq_gt_${rsqthresh}
