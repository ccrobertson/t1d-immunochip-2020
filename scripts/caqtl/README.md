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

