<!-- TOC depthFrom:1 depthTo:3 withLinks:1 updateOnSave:1 orderedList:1 -->

1. [Pre-processing](#pre-processing)
2. [Pre-imputation analysis](#pre-imputation-analysis)
3. [Define ichip regions](#define-ichip-regions)
4. [Post imputation processing](#post-imputation-processing)
5. [Imputation coverage and accuracy](#imputation-coverage-and-accuracy)
6. [Family-based association analysis](#family-based-association-analysis)
7. [Case-control association analysis](#case-control-association-analysis)
8. [Meta-analysis](#meta-analysis)
9. [Inflation comparison](#inflation-comparison)
10. [Fine-mapping](#fine-mapping)
11. [Variant annotation](#variant-annotation)
12. [QTL colocalisation](#qtl-colocalisation)

<!-- /TOC -->

# Summary

Analysis pipeline for "Fine-mapping, trans-ancestral and genomic analyses identify causal variants, cells, genes and drug targets for type 1 diabetes"

Pipeline authors: Jamie Inshaw & Cassie Robertson (ccr5ju@virginia.edu)


## iii. Pre-processing

Map ImmunoChip probe sequences to hg19
  ```bash
  ${scripts}/probe_alignment/README_map_probes
  ```

Run pre-processing pipeline
  ```bash
  cd ${PROJECT_MEGA}/scripts/pre_qc
  make --quiet > make.log
  ```

Move to final directory
  ```bash
  cp ${PROJECT_MEGA}/pre_qc/Mega_15Jun2018/data_clean/* ${PROJECT_MEGA}/mega_genotyped/release2
  cp ${PROJECT_MEGA}/pre_qc/Mega_06Aug2018/data_clean/* ${PROJECT_MEGA}/mega_genotyped/release3
  cp ${PROJECT_MEGA}/pre_qc/Mega_10April2019/data_clean/* ${PROJECT_MEGA}/mega_genotyped/release4
  ```

Post pre-processing releases:

- release1: 2017-05-26 (b37, 56693 subjects)
- release2: 2018-06-15 (b38, 58187 subjects)
- release3 (added cbr): 2018-08-09 (b37, 61104 subjects)
- release4 (fixed additions export): 2019-04-06 (b37, 61427 subjects)



## iv. Pre-imputation analysis




Cohort summary tables for manuscript
  ```bash
  Rscript ${scripts}/manuscript/summary_tables.R
  ```


## v. Define ichip regions
(insert script from Jamie)
  ```bash
  bash ${scripts}/define_ichip_regions.sh
  ```

Add cytobands to the ichip region file
  ```bash
  Rscript ${scripts}/define_ichip_regions/add_cytobands.R
  ```

Add updated ichip genotyped variant counts (based only on those passing QC)

## vi. Post imputation processing
  ```bash
  bash ${scripts}/post_imputation/post_imputation_master.sh {group} {password} {mafthresh} {rsqthresh}
  bash ${scripts}/post_imputation/post_imputation_master_EUR.sh {password1} {password2} {password3} {mafthresh} {rsqthresh}
  ```

	Get combined info files
	```bash
	bash ${scripts}/post_imputation/combine_info_files.sh
	```

## vii. Imputation coverage and accuracy

Get WashU WGS data
```bash

	#find overlapping subjects and create keeplists
	bash ${scripts}/imputation_accuracy_coverage/get_washu_mega_overlap_and_keeplists.sh
	#extract mega subjects and filter for MAF>0.005
  sbatch --array=1-22 --output=${impbench}/filter_washu_files_AFR_%a.log ${scripts}/imputation_accuracy_coverage/filter_washu_files.slurm AFR
	sbatch --array=1-22 --output=${impbench}/filter_washu_files_AMR_%a.log ${scripts}/imputation_accuracy_coverage/filter_washu_files.slurm AMR
	sbatch --array=1-22 --output=${impbench}/filter_washu_files_EUR_%a.log ${scripts}/imputation_accuracy_coverage/filter_washu_files.slurm EUR

	#some of the filtering timed out so had to rerun
	sbatch --output=${impbench}/filter_washu_files_AFR_chr2_restart.log ${impbench}/filter_washu_files_AFR_chr2.slurm
	sbatch --output=${impbench}/filter_washu_files_AMR_chr1_restart.log ${impbench}/filter_washu_files_AMR_chr1.slurm
  ```

Prep data sets
```bash
prepDataForAccuracy () {
	local group=$1
	sbatch --array=1-22 --output=${impbench}/prep_washu_data_${group}_%a.log ${scripts}/imputation_accuracy_coverage/prep_washu_data.slurm ${group}
	#sbatch --array=1-22 --output=${impbench}/prep_mega_data_IMPUTED_TOPMED_${group}_%a.log ${scripts}/imputation_accuracy_coverage/prep_mega_data.slurm ${group} IMPUTED_TOPMED gr38
	#sbatch --array=1-22 --output=${impbench}/prep_mega_data_IMPUTED_1KG_${group}_%a.log ${scripts}/imputation_accuracy_coverage/prep_mega_data.slurm ${group} IMPUTED_1KG gr37
}
prepDataForAccuracy AFR
prepDataForAccuracy	AMR
prepDataForAccuracy	EUR
prepDataForAccuracy	FIN

#NOTE FOR AFR - ERROR in CHR2 and CHR14 due to duplicate variant IDs in 1000G imputed file (possibly created during liftover?)
#Manual fix:
mv mega_IMPUTED_1KG_AFR_chr2_ichip_regions.vcf.gz mega_IMPUTED_1KG_AFR_chr2_ichip_regions_DUPS.vcf.gz
mv mega_IMPUTED_1KG_AFR_chr14_ichip_regions.vcf.gz mega_IMPUTED_1KG_AFR_chr14_ichip_regions_DUPS.vcf.gz
zcat mega_IMPUTED_1KG_AFR_chr2_ichip_regions_DUPS.vcf.gz | grep -v "2:61390913:C:<CN0>" | bgzip > mega_IMPUTED_1KG_AFR_chr2_ichip_regions.vcf.gz
zcat mega_IMPUTED_1KG_AFR_chr14_ichip_regions_DUPS.vcf.gz  | grep -v "14:88504832:\|14:98480754:" | bgzip > mega_IMPUTED_1KG_AFR_chr14_ichip_regions.vcf.gz
sbatch --array=2,14 --output=${impbench}/calculate_imputation_accuracy_IMPUTED_1KG_AFR_%a.log ${scripts}/imputation_accuracy_coverage/calculate_imputation_accuracy.slurm AFR IMPUTED_1KG
```

Compute accuracy
```bash
calcAccuracy () {
local group=$1
sbatch --array=1-22 --output=${impbench}/calculate_imputation_accuracy_IMPUTED_TOPMED_${group}_%a.log ${scripts}/imputation_accuracy_coverage/calculate_imputation_accuracy.slurm ${group} IMPUTED_TOPMED
sbatch --array=1-22 --output=${impbench}/calculate_imputation_accuracy_IMPUTED_1KG_${group}_%a.log ${scripts}/imputation_accuracy_coverage/calculate_imputation_accuracy.slurm ${group} IMPUTED_1KG
}
calcAccuracy AFR
calcAccuracy AMR
calcAccuracy EUR

#some jobs ran out of memory or timed out, re-running these here on largemem
sbatch --array=1-12 --output=${impbench}/calculate_imputation_accuracy_IMPUTED_TOPMED_AFR_%a_restart.log ${impbench}/calculate_imputation_accuracy_restart.slurm AFR IMPUTED_TOPMED
sbatch --array=2 --output=${impbench}/calculate_imputation_accuracy_IMPUTED_TOPMED_AMR_%a_restart.log ${impbench}/calculate_imputation_accuracy_restart.slurm AMR IMPUTED_TOPMED
```

Count genotyped variants in ichip regions
```bash
plink1.9 --bfile ${genodat}_filtered_for_hwe --freq --out ${genodat}_filtered_for_hwe
awk 'NR>1{print $2, $3, $4, $1}' ${PROJECT_MEGA}/define_ichip_regions/ichip_regions_gr37.txt > ${PROJECT_MEGA}/define_ichip_regions/ichip_regions_gr37_plinkformat.txt
plink1.9 --bfile ${genodat}_filtered_for_hwe --extract ${PROJECT_MEGA}/define_ichip_regions/ichip_regions_gr37_plinkformat.txt --range --make-bed --out ${genodat}_filtered_for_hwe_ichip_regions
plink1.9 --bfile ${genodat}_filtered_for_hwe_ichip_regions --freq --out ${genodat}_filtered_for_hwe_ichip_regions
awk '$5>0.005' ${genodat}_filtered_for_hwe_ichip_regions.frq | wc -l
```

Count imputed variants in ichip regions
```bash
bash ${scripts}/imputation_accuracy_coverage/extract_ichip_regions.sh EUR
bash ${scripts}/imputation_accuracy_coverage/extract_ichip_regions.sh AFR
bash ${scripts}/imputation_accuracy_coverage/extract_ichip_regions.sh FIN
bash ${scripts}/imputation_accuracy_coverage/extract_ichip_regions.sh AMR
```

```R
for (group in c("AFR","EUR","FIN","EAS","AMR")) {
  d = read.table(paste0("count_imputed_variants_in_ichip_regions_",group,".txt"))
  print(paste(group, sum(d$V2), sum(d$V2)/84923))
}
```

Summarize accuracy results
  ```bash
	#merge results across chromosomes
	sbatch --output=${impbench}/merge_accuracy_results.log ${scripts}/imputation_accuracy_coverage/merge_accuracy_results.slurm
	#compare and visualize results for between refpanels
  Rscript ${scripts}/imputation_accuracy_coverage/summarize_accuracy.R ichip IMPUTED_TOPMED
  Rscript ${scripts}/imputation_accuracy_coverage/compare_accuracy_topmed_vs_1kg.R
  ```



## viii. Family-based association analysis

Run TDT pipeline
  ```bash
  bash ${scripts}/assoc/family_analysis_master.sh {group}
  ```

Plot TDT results across data sets
  ```bash
  sbatch ${scripts}/assoc/family_analysis_plots.slurm {group}
  ```

## ix. Case-control association analysis
See README for scripts and running order in:
```
${scripts}/assoc_case_control/README.md
```

## x. Meta-analysis
See README for scripts and running order in:
```
${scripts}/assoc_case_control/README.md
```


Replication - get plot comparing association statistics in case-control and family analysis
  ```bash
  Rscript ${scripts}/replication.R
  ```

## xi. Inflation comparison
Generate family-only meta-analysis summary statistics
  ```bash
  cd ${meta}
  metal ${scripts}/inflation_analysis/metal_script_for_fam_only_metaanaylsis
  mv METAANALYSIS1.TBL METAANALYSIS_fam_only_vcf_newqc_1.TBL
  mv METAANALYSIS1.TBL.info METAANALYSIS_fam_only_vcf_newqc_1.info
  ```

Create five sub-sampled case-control data sets
  ```bash
  for i in {1..5}; do
    mkdir -p ${inflation}/subdata/${i}
    mkdir -p ${inflation}/subanalyses/${i}
  done
  Rscript ${scripts}/inflation_analysis/subset_cohorts.R
  bash ${scripts}/inflation_analysis/subset_cohorts.sh > subset_cohorts.log
  ```

Run SNPTEST on sub-sampled data
  ```bash
  sbatch --array=1-5 ${scripts}/inflation_analysis/run_association_on_subsets.slurm AFR
  sbatch --array=1-5 ${scripts}/inflation_analysis/run_association_on_subsets.slurm AMR
  sbatch --array=1-5 ${scripts}/inflation_analysis/run_association_on_subsets.slurm EUR
  sbatch --array=1-5 ${scripts}/inflation_analysis/run_association_on_subsets.slurm FIN
  ```

Compare sub-sampled association statistic distributions to family-based association stats with equivalent power
  ```bash
  Rscript ${scripts}/inflation_analysis/visualize_inflation.R
  ```

Calculate number of informative trios across MAF spectrum
  ```R
  getNum = function(p, n) { 2*p*(1-p)*n }
  nums = NULL
  mafs = seq(0,0.5, by=0.001)
  for (i in 1:length(mafs)) { nums[i]=getNum(mafs[i],4578)}
  d = data.frame(mafs, nums)
  plot(d$mafs, d$nums, xlab="MAF", ylab="Number of informative trios")
  ```

## xii. Fine-mapping
See README for scripts and running order, which includes GUESSFM (using EUR only) and PAINTOR (using all ancestry groups) in:
```
${scripts}/fine_mapping/README.md
```


## xiii. Variant annotation
Convert to ANNOVAR format
```bash
Rscript ${scripts}/fine_mapping/prep_for_annovar.R
```

Annotate with ANNOVAR
```bash
ANNOVAR=${resources}/annovar

perl ${ANNOVAR}/table_annovar.pl ${annot}/annovar_input_file.txt ${ANNOVAR}/humandb/ -buildver hg38 -out ${annot}/annovar_output -protocol refGene,ensGene,cytoBand,1000g2015aug_eur,1000g2015aug_afr,avsnp150,dbnsfp30a -operation g,g,r,f,f,f,f -nastring . -csvout

perl ${ANNOVAR}/table_annovar.pl ${annot}/annovar_input_file_credible_sets_with_proxies.txt ${ANNOVAR}/humandb/ -buildver hg38 -out ${annot}/annovar_output_credible_sets_with_proxies -protocol refGene,ensGene,cytoBand,1000g2015aug_afr,1000g2015aug_eur,avsnp150,dbnsfp30a -operation g,g,r,f,f,f,f -nastring . -csvout
```

Annotations for all the variants included in our analysis are in ${annot}/annovar_output.hg38_multianno.csv

Annotations for credible sets are in ${annot}/annovar_output_credible_sets_with_proxies.hg38_multianno.csv

Create annotated fine-mapping table (for manuscript)
```bash
Rscript ${scripts}/add_annotation_to_credible_set_table.R
```

## cis-caQTL scan

Link to DP3 data
```bash
linkToProcessedData () {
	local freezename=$1
	mkdir -p ${qtl}/${freezename}
	ln -s -f ${PROJECT_DP3}/analysis/${freezename}/processed_merged_peaks/processed_sampleDF.rds ${qtl}/${freezename}/processed_sampleDF.rds
	ln -s -f ${PROJECT_DP3}/analysis/${freezename}/processed_merged_peaks/processed_sampleDF.txt ${qtl}/${freezename}/processed_sampleDF.txt
	ln -s -f ${PROJECT_DP3}/analysis/${freezename}/processed_merged_peaks/processed_countMat_clean_voomnormalized.rds ${qtl}/${freezename}/processed_countMat_clean_voomnormalized.rds
	ln -s -f ${PROJECT_DP3}/analysis/${freezename}/processed_merged_peaks/processed_countMat_autosome.rds ${qtl}/${freezename}/processed_countMat_autosome.rds
}

linkToProcessedData "freeze1_SQB1to8"
linkToProcessedData "freeze2_SQB1to16"

```

1. extract genotype data for subjects with ATAC-seq
  ```bash
	#figure out which samples to extract
	awk 'NR>1 {print $1}' ${freeze}/processed_sampleDF.txt | sort | uniq > ${freeze}/atac_samples_iids.txt
	getSampleIIDs() {
	  local group=$1
	  bcftools query --list-samples ${filtered}/${group}/all/chr22.filter_maf_gt_0.005_rsq_gt_0.8.dose.vcf.gz | awk -v file=${freeze}/atac_samples_iids.txt 'BEGIN {while(getline<file){a[$1]=1}} {i=1; while(i<=NF) { if ($i in a) {print $i}; i++}}' > ${freeze}/atac_samples_iids_${group}.txt
	}
	getSampleIIDs AFR
	getSampleIIDs EUR

	#extract samples for each chromosome
	mkdir -p ${freeze}/mega_genotypes
  sbatch --output=${freeze}/extract_mega_data_%a.log --array=1-22 ${scripts}/caqtl/extract_mega_data.slurm

  #concatenate into one vcf
  >${freeze}/mega_genotypes/mega_filelist.txt
  for i in {1..22}; do
      echo "${freeze}/mega_genotypes/mega_chr${i}.vcf.gz" >> ${freeze}/mega_genotypes/mega_filelist.txt
  done
  bcftools concat -f ${freeze}/mega_genotypes/mega_filelist.txt -O z -o ${freeze}/mega_genotypes/mega.vcf.gz
  tabix -f -p vcf ${freeze}/mega_genotypes/mega.vcf.gz
  ```

2. Make sure samples match genotypes
This script goes through the processed_sampleDF.txt file, and for each row/sample, it runs QTLtools MBV and flags the sample if the top match is not the same subject or if top match is <99%.
```bash
mkdir -p ${freeze}/sample_match

Rscript ${scripts}/caqtl/sample_matching_run_QTLtools.R

#note, this script needs to be run interactively still
Rscript ${scripts}/caqtl/sample_matching_master.R
```

3. Prepare data for analysis
```bash
sbatch --output=${freeze}/prep_genotypes.log ${scripts}/caqtl/prep_genotypes.slurm
sbatch --output=${freeze}/prep_data.log ${scripts}/caqtl/prep_data.slurm
```

4. Run peer algorithm on data sets
```bash
sbatch --output=${freeze}/peer_correction_DATASET_test.log ${scripts}/caqtl/peer_correction.slurm DATASET_test
sbatch --output=${freeze}/peer_correction_DATASET_un.log ${scripts}/caqtl/peer_correction.slurm DATASET_un
sbatch --output=${freeze}/peer_correction_DATASET_stim.log ${scripts}/caqtl/peer_correction.slurm DATASET_stim
sbatch --output=${freeze}/peer_correction_DATASET_un_EUR.log ${scripts}/caqtl/peer_correction.slurm DATASET_un_EUR
sbatch --output=${freeze}/peer_correction_DATASET_un_AFR.log ${scripts}/caqtl/peer_correction.slurm DATASET_un_AFR
```

5. run chromatin QTL scan
```bash
sbatch --output ${freeze}/caqtl_scan.log ${scripts}/caqtl/caqtl_scan.slurm
```

6. conditional and joint analysis for significant peaks
```bash
Rscript ${scripts}/caqtl/meta_qtl.R
sbatch --output ${freeze}/conditional_caqtl.log ${scripts}/caqtl/conditional_caqtl.slurm
```

## cis-QTL colocalisation

Run coloc
```bash
mkdir -p ${freeze}/locuscompare
Rscript ${scripts}/fine_mapping/colocalisation.R


### Troubleshoot differences between freezes in IKZF3 region
#Look for peaks
getRegion = function(x){
  s = strsplit(x,split="_")[[1]]
  n = length(s)
  c(s[1], s[(n-1)],s[n])
}
REGmat = do.call("rbind",lapply(row.names(dat_un[["counts"]]),getRegion))
REG = data.frame(geneid=row.names(dat_un[["counts"]]), chr=as.numeric(gsub("chr","",REGmat[,1])), left=as.numeric(REGmat[,2]),  right=as.numeric(REGmat[,3]))
chr17_regions = REG[REG$chr==17 & (REG$left>(39866973-5000) & REG$left<(39868416+5000) | REG$left>(39939542-5000) & REG$left<(39940558+5000)),]
# geneid chr     left    right
# 19357 chr17_39863823_39864542  17 39863823 39864542
# 19358 chr17_39866973_39868400  17 39866973 39868400  <--- this is the stim-responsive peak
# 19361 chr17_39939542_39940558  17 39939542 39940558

res_EUR = read.table("caqtl_scan_unstim_cis_EUR.txt", header=TRUE)
res_EUR[res_EUR$gene=="chr17_39866973_39868400",]

res_AFR = read.table("caqtl_scan_unstim_cis_AFR.txt", header=TRUE)
res_AFR[res_AFR$gene=="chr17_39866973_39868400",]

res_AFR[res_AFR$gene=="chr17_39866973_39868400",]$SNP %in% res_EUR[res_EUR$gene=="chr17_39866973_39868400",]$SNP

res_meta = read.table("caqtl_scan_unstim_cis_meta.txt", header=TRUE)
res_meta[res_meta$gene=="chr17_39866973_39868400",]

res_EUR[res_EUR$SNP=="chr17_39883308_TAACA_T",]
res_AFR[res_AFR$SNP=="chr17_39883308_TAACA_T",]
```

Calculate allelic imbalance at colocalized caQTLs
```bash
bash ${scripts}/fine_mapping/allelic_imbalance.R

## Trick pipeline into running on rs705705 (in perfect LD with rs705704)
cp rs705704/rs705704_HET.txt rs705705/rs705705_HET.txt
cp rs705704/rs705704_HOMREF.txt rs705705/rs705705_HOMREF.txt
cp rs705704/rs705704_HOMALT.txt rs705705/rs705705_HOMALT.txt
```

Create colocalisation table
```bash
Rscript ${scripts}/fine_mapping/merge_coloc_ecaviar_results.R
```


## Generate fine-mapping region plots

Get annotation tracks from BLUEPRINT and ENCODE
```bash
cd ${PUBLIC_DATA}/EncodeReadmap
wget http://dcc.blueprint-epigenome.eu/#/md/secondary_analysis/Segmentation_of_ChIP-Seq_data_20140811

cd ${PUBLIC_DATA}/Blueprint
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/
```

Visualize credible set overlap with tracks in each region (for Supplemental Data)
```bash
Rscript ${scripts}/fine_mapping/generate_region_tracks.R
```

Create caQTL tracks

ANKRD55: rs10213692_chr5_56146422_T_C
BACH2: rs72928038_chr6_90267049_G_A
CENPW: rs9388486_chr6_126340008_T_C
RP11-6101.1: rs3902659_chr14_98021110_G_A

For each SNP, extract region by genotype and create mean track.
```bash
bash ${scripts}/fine_mapping/chromqtl_tracks_get_region_by_genotype.sh
```
