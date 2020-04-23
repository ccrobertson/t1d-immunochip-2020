#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account rich_immunochip_impute
#SBATCH --partition=standard
#SBATCH --mem=100000

chr=${SLURM_ARRAY_TASK_ID}

module load bcftools/1.9

set -euxo pipefail


#extract subjects
extractSubjects() {
  local group=$1
  echo "extracting chr${chr} for ${group}"
  bcftools view --samples-file ${freeze}/atac_samples_iids_${group}.txt -O z -o ${freeze}/mega_genotypes/mega_${group}_chr${chr}.vcf.gz ${filtered}/${group}/all/chr${chr}.filter_maf_gt_0.005_rsq_gt_0.8.dose.vcf.gz
  tabix -f -p vcf ${freeze}/mega_genotypes/mega_${group}_chr${chr}.vcf.gz
}
extractSubjects AFR
extractSubjects EUR

wait

echo "merging chr${chr}"
bcftools merge ${freeze}/mega_genotypes/mega_AFR_chr${chr}.vcf.gz ${freeze}/mega_genotypes/mega_EUR_chr${chr}.vcf.gz -O z -o ${freeze}/mega_genotypes/mega_chr${chr}.vcf.gz
tabix -f -p vcf ${freeze}/mega_genotypes/mega_chr${chr}.vcf.gz
