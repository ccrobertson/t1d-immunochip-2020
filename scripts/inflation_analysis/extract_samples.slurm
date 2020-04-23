#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account rich_immunochip_impute
#SBATCH --partition=standard
#SBATCH --mem=2000


group=$1
subset=$2
chr=$3

module load htslib
module load bcftools/1.9

date
vcfin=${filtered}/${group}/unrelateds/chr${chr}.filter_maf_gt_0.005_rsq_gt_0.8.unrelateds.dose.vcf.gz
bcftools view --samples-file ${inflation}/${group}_sub${subset}.txt --force-samples --output-type z --output-file ${inflation}/subdata/${subset}/chr${chr}.${group}_sub${subset}.unrelateds.dose.vcf.gz ${vcfin}
tabix -p vcf ${inflation}/subdata/${subset}/chr${chr}.${group}_sub${subset}.unrelateds.dose.vcf.gz

echo "DONE"
date
