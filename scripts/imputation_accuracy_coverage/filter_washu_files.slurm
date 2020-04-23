#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account rich_immunochip_impute
#SBATCH --partition=largemem
#SBATCH --mem=200000

chr=${SLURM_ARRAY_TASK_ID}
group=$1

echo -e "START CHR${chr}"

#Usage: bash filter_washu_files.sh {input vcf} {output vcf}
invcf=/nv/vol185/T1DGC/MGI_WashU/vcf-batch1/concatenated.cchr${chr}.ab_filtered.vcf.gz
outvcf=${impbench}/washu_${group}_chr${chr}_filtered.gz

module load bcftools

echo "Extracting MEGA subjects"
bcftools view --samples-file <(awk -v g=${group} '$5==g {print $1}' ${impbench}/washu_batch1_in_mega.txt) --force-samples --output-type z --output-file ${outvcf}_tmp1.gz ${invcf}

echo "Filtering for MAF>0.005"
bcftools view --min-af 0.005 --output-type z --output-file ${outvcf} ${outvcf}_tmp1.gz

echo "Removing tmp file"
rm -f ${outvcf}_tmp1.gz

echo -e "FINISHED CHR${chr}"
