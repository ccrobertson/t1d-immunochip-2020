#!/bin/bash


getCombinedInfo () {
	local group=$1
  echo ${group}
	zcat ${results}/${group}/chr1.info.gz | awk 'NR==1' > ${results}/${group}/chrALL.info
	for i in {1..22}; do
    echo "chr${i}"
		zcat ${results}/${group}/chr${i}.info.gz | awk 'NR>1' >> ${results}/${group}/chrALL.info
	done
  awk '$8=="Genotyped"' ${results}/${group}/chrALL.info > ${results}/${group}/typed_variants.txt
  echo " "
}

getCombinedInfo EUR1
getCombinedInfo EUR2
getCombinedInfo EUR3
getCombinedInfo EAS
getCombinedInfo AFR
getCombinedInfo AMR
getCombinedInfo FIN
