#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account rich_immunochip_impute
#SBATCH --partition=standard
#SBATCH --mem=10000


sample=$1
SQB=$2
group=$3

bam=${PROJECT_DP3}/processed/${SQB}/output/results_pipeline/${sample}/aligned_hg38/${sample}_sort_dedup.bam
vcf=${filtered}/${group}/all/chr1.filter_maf_gt_0.005_rsq_gt_0.8.dose.vcf.gz

QTLtools mbv --bam ${bam} --vcf ${vcf} --out ${freeze}/sample_match_full_cohort/${sample}.match_${group}.txt
