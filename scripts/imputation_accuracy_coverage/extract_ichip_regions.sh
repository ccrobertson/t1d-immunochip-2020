#!/bin/bash

module load htslib

group=$1
tallyfile=${impbench}/count_imputed_variants_in_ichip_regions_${group}.txt

> ${tallyfile}
for i in {1..22}; do
  echo "chr${i}"
  tabix -h -R ${ichip_regions}/ichip_regions_gr38_chr${i}.txt ${filtered}/${group}/all/chr${i}.filter_maf_gt_0.005_rsq_gt_0.8.dose.vcf.gz | bgzip > ${filtered}/${group}/all/chr${i}.filter_maf_gt_0.005_rsq_gt_0.8_ichip_regions.dose.vcf.gz
  no_variants=$(zcat ${filtered}/${group}/all/chr${i}.filter_maf_gt_0.005_rsq_gt_0.8_ichip_regions.dose.vcf.gz | grep -v '^#' | wc -l)
  echo "chr${i} ${no_variants}" >> ${tallyfile}
done
