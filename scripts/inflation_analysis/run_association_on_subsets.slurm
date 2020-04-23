#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account rich_immunochip_impute
#SBATCH --partition=standard
#SBATCH --mem=2000

group=$1
subset=${SLURM_ARRAY_TASK_ID}

set -euxo pipefail
module load bcftools

echo "RUNNING SNPTEST ASSOCIATION"
getAssocResults () {
  local group=$1
  local subset=$2
  sample_list_ordered=${inflation}/subdata/${subset}/${group}_sub${subset}samps_list_ordered
  sample_file_subset=${inflation}/subdata/${subset}/${group}_sub${subset}samps_snptest
  bcftools query --list-samples ${inflation}/subdata/${subset}/chr1.${group}_sub${subset}.unrelateds.dose.vcf.gz > ${sample_list_ordered}
  awk 'NR<=2' ${inflation}/${group}samps_snptest > ${sample_file_subset}
  awk -v file=${inflation}/${group}samps_snptest 'BEGIN {while(getline<file) {a[$2]=$0}} {print a[$1]}' ${sample_list_ordered} >> ${sample_file_subset}
  for ((chr=1; chr<=22; chr+=1)); do
    echo -e "chr${chr}"
    snptest_v2.5.4-beta3 -data ${inflation}/subdata/${subset}/chr${chr}.${group}_sub${subset}.unrelateds.dose.vcf.gz ${sample_file_subset} -filetype vcf -include_samples ${sample_list_ordered} -genotype_field GP -pheno t1d -cov_all_continuous -frequentist 1 -method newml -o ${inflation}/subanalyses/${subset}/chr${chr}.${group}_sub${subset}.tbl -log ${inflation}/subanalyses/${subset}/chr${chr}.${group}_sub${subset}.log
    Rscript ${scripts}/inflation_analysis/convert_to_metal.R ${inflation}/subanalyses/${subset}/chr${chr}.${group}_sub${subset}.tbl ${inflation}/subanalyses/${subset}/metal_chr${chr}.${group}_sub${subset}.tbl
  done
}

getAssocResults ${group} ${subset}

echo "COMBINING RESULTS ACROSS CHROMOSOMES"
awk 'NR==1' ${inflation}/subanalyses/${subset}/metal_chr1.${group}_sub${subset}.tbl > ${inflation}/subanalyses/${subset}/metal_${group}_sub${subset}.tbl
for ((chr=1; chr<=22; chr+=1)); do
  echo -e "chr${chr}"
  awk 'NR>1' ${inflation}/subanalyses/${subset}/metal_chr${chr}.${group}_sub${subset}.tbl >> ${inflation}/subanalyses/${subset}/metal_${group}_sub${subset}.tbl
done
