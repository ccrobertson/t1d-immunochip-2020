#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account rich_immunochip_impute
#SBATCH --partition=standard
#SBATCH --mem=20000

module load htslib
module load bcftools

group=$1
chr=${SLURM_ARRAY_TASK_ID}

invcf=${impbench}/washu_${group}_chr${chr}_filtered.gz
outprefix=${impbench}/washu_${group}_chr${chr}_filtered

echo "Extracting subjects"
bcftools view --samples-file ${impbench}/washu_batch1_washu_keeplist_${group}.txt --force-samples --output-type z --output-file ${outprefix}.vcf.gz ${invcf}

echo "Indexing ${outprefix}.vcf.gz"
tabix -f -p vcf ${outprefix}.vcf.gz

echo "Extracting ichip regions"
tabix -h -R ${PROJECT_MEGA}/define_ichip_regions/ichip_regions_gr38_chr${chr}.txt ${outprefix}.vcf.gz | bgzip > ${outprefix}_ichip_regions.vcf.gz

echo "Indexing ${outprefix}_ichip_regions.vcf.gz"
tabix -f -p vcf ${outprefix}_ichip_regions.vcf.gz

echo "DONE"
