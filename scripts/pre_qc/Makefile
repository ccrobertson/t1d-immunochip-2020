
qctitle=Mega_10April2019
qclocation=${PROJECT_MEGA}/pre_qc

#GNU Make tutorial 
#https://bsmith89.github.io/make-bml/#makefile-basics-45-minutes

### Define global environment variables (inherited by all subprocesses)
export resources=${PROJECT_MEGA}/resources/

export preqc_folder=${qclocation}/${qctitle}
export scripts=${preqc_folder}/scripts
export data_raw=${preqc_folder}/data_raw
export data_derived=${preqc_folder}/data_derived
export clean=${preqc_folder}/data_clean
export prelim=${preqc_folder}/prelim_QC
export qclists=${prelim}/QC_lists
export sexqc=${prelim}/SEXQC
export relqc=${prelim}/RELATEDQC
export snpqc=${prelim}/SNPQC
export phenofile=${data_derived}/combined_pheno.txt
export logs=${preqc_folder}/logs

export nickname=mega
export nickname_1=mega_main
export nickname_2=mega_cbr
export nickname_3=mega_jhu
export nickname_4=mega_additions
export nicknames=${nickname_1} ${nickname_2} ${nickname_3} ${nickname_4}
export PATH=/apps/software/standard/compiler/gcc/7.1.0/R/3.5.1/bin:/apps/software/standard/core/htslib/1.9/bin:/apps/software/standard/core/anaconda/5.2.0-py3.6:/apps/software/standard/core/anaconda/5.2.0-py3.6/bin:/apps/software/standard/compiler/gcc/7.1.0/bedops/2.4.1/bin:/apps/software/standard/core/gcc/7.1.0/bin:/apps/software/standard/core/java/1.8.0:/apps/software/standard/core/java/1.8.0/bin:/usr/lib64/qt-3.3/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/opt/slurm/current/bin:/opt/slurm/current/sbin:/opt/singularity/current/bin:/opt/rci/bin:/opt/rci/sbin:/opt/nhc/current/sbin:/share/rci_apps/common/bin:/share/resources/HPCtools:/nm/t1/users/ccr5ju/software/bin:/nm/t1/users/ccr5ju/software/pigz-2.4:/home/ccr5ju/.local/bin:/home/ccr5ju/bin:/apps/software/standard/core/rstudio/1.0.143/bin/pandoc


### Define local envorinment variables (only accessible within Makefile)
RECIPES=${PROJECT_MEGA}/scripts/pre_qc/make_recipes
SUB_SCRIPTS=${RECIPES}/sub_scripts


.PHONY: all clean updated_scripts

all: updated_scripts ${data_raw}/phenotype_files/combined_pheno_additional.txt ${phenofile} ${data_derived}/mydrop ${data_derived}/export_bias.html ${clean}/${nickname}_b37_clean.qcreport
	@echo "FINISHED"

### Setup for qc pipeline
${data_raw}/phenotype_files/combined_pheno_additional.txt ${data_raw}/sample_qc_files/T1DGC_SAMPLEQC_SWITCH_LIST_201105.txt ${data_raw}/snpqc_files/mymapInfo: 
	@echo "retrieving data and scripts"
	mkdir -p ${logs}
	bash ${RECIPES}/setup.sh > ${logs}/setup.log 2>&1

### Retrieve scripts -- run this every time
updated_scripts: ${data_raw}/phenotype_files/combined_pheno_additional.txt
	@echo "retrieving sub_scripts"
	cp -u ${SUB_SCRIPTS}/* ${scripts}
	chmod u=rwx ${scripts}/*

### Create combined phenofile
${phenofile}: ${data_raw}/phenotype_files/combined_pheno_additional.txt
	@echo "creating phenofile"
	bash ${RECIPES}/create_combined_pheno.sh > ${logs}/create_combined_pheno.log 2>&1

### Generate Plink Files [ only if explicitly requested ]
${data_derived}/mega_main_raw/mega_main0.bed ${data_derived}/genome_studio_export_map.txt: 
	@echo "generating plink files"
	bash ${RECIPES}/generate_plink_files.sh > ${logs}/generate_plink_files.log 2>&1

### Updating duplicate sample ids 
${data_derived}/mega_main_raw/mega_main1.bed: ${data_derived}/mega_main_raw/mega_main0.bed
	@echo "updating duplicate sample ids"
	bash ${RECIPES}/update_dup_sample_ids.sh > ${logs}/update_dup_sample_ids.log 2>&1
	
### Update SNP map
${data_derived}/mydrop: ${data_derived}/mega_main_raw/mega_main0.bed ${data_raw}/snpqc_files/mymapInfo
	@echo "generating map update files"
	bash ${RECIPES}/generate_snp_position_update_files.sh > ${logs}/generate_snp_position_update_files.log 2>&1

### Update SNP positions
${data_derived}/mega_main_raw/mega_main4.bed: ${data_derived}/mega_main_raw/mega_main1.bed ${data_derived}/mydrop
	@echo "updating map"
	bash ${RECIPES}/update_snp_positions.sh > ${logs}/update_snp_positions.log 2>&1

### Merge plink files
${data_derived}/${nickname}_raw_nophenoA.bed: ${data_derived}/mega_main_raw/mega_main4.bed
	@echo "merging plink files"
	bash ${RECIPES}/merge_plink_files.sh > ${logs}/merge_plink_files.log 2>&1

### Compare exports	
${data_derived}/export_bias.html: ${data_derived}/mega_rawD.bed ${data_derived}/mega_raw_nophenoA.bed ${data_derived}/genome_studio_export_map.txt 
	@echo "comparing qc metrics across exports"
	bash ${RECIPES}/compare_export_files.sh	> ${logs}/compare_export_files.log 2>&1
	
### Remove duplicate snps
${data_derived}/${nickname}_raw_nophenoC.bed: ${data_derived}/${nickname}_raw_nophenoA.bed
	@echo "removing duplicate snps"
	bash ${RECIPES}/remove_dup_snps.sh > ${logs}/remove_dup_snps.log 2>&1
	
### Update rsids to b150 
${data_derived}/${nickname}_raw_nophenoD.bed: ${data_derived}/${nickname}_raw_nophenoC.bed ${data_raw}/snpqc_files/mymapInfo
	@echo "updating rsids to b150"
	bash ${RECIPES}/update_rsid_to_b150.sh > ${logs}/update_rsid_to_b150.log 2>&1

### Update sample ids using provided maps
${data_derived}/${nickname}_raw_nophenoE.bed: ${data_derived}/${nickname}_raw_nophenoD.bed ${data_raw}/sample_qc_files/T1DGC_SAMPLEQC_SWITCH_LIST_201105.txt
	@echo "updating sample ids using provided maps"
	bash ${RECIPES}/update_sample_ids_with_provided_maps.sh > ${logs}/update_sample_ids_with_provided_maps.log 2>&1

### Add family ids
${data_derived}/${nickname}_raw_nophenoF.bed: ${data_derived}/${nickname}_raw_nophenoE.bed ${phenofile}
	@echo "adding family ids"
	bash ${RECIPES}/add_family_ids.sh > ${logs}/add_family_ids.log 2>&1

### Fix sample switches based on May 2011 T1DGC release
${data_derived}/${nickname}_raw_nophenoI.bed: ${data_derived}/${nickname}_raw_nophenoF.bed ${data_raw}/sample_qc_files/T1DGC_SAMPLEQC_SWITCH_LIST_201105.txt
	@echo "fixing sample switches based on May 2011 release"
	bash ${RECIPES}/fix_sample_swaps_t1dgc2011release.sh > ${logs}/fix_sample_swaps_t1dgc2011release.log 2>&1

### Add parents and phenotypes
${data_derived}/${nickname}_rawA.bed: ${data_derived}/${nickname}_raw_nophenoI.bed ${phenofile}
	@echo "adding parents and phenotypes"
	bash ${RECIPES}/add_parents_and_phenotypes.sh > ${logs}/add_parents_and_phenotypes.log 2>&1
	
###	Remove samples based on predefined lists
${data_derived}/${nickname}_rawD.bed: ${data_derived}/${nickname}_rawA.bed ${data_raw}/sample_qc_files/T1DGC_SAMPLEQC_SWITCH_LIST_201105.txt
	@echo "removing samples based on predefined lists"
	bash ${RECIPES}/remove_samples_in_provided_droplists.sh > ${logs}/remove_samples_in_provided_droplists.log 2>&1
	
###	Preliminary sample QC
${data_derived}/${nickname}_rawE.bed: ${data_derived}/${nickname}_rawD.bed
	@echo "preliminary sample qc, filtering for call rate<80%"
	bash ${RECIPES}/preliminary_sample_qc.sh > ${logs}/preliminary_sample_qc.log 2>&1

### Sex QC1
${sexqc}/SEXQC1/SEXCHECK1/sex_anomalies.txt: ${data_derived}/${nickname}_rawE.bed 
	@echo "running sex qc pipeline (round 1)"
	bash ${RECIPES}/sexqc2.sh ${data_derived}/${nickname}_rawE SEXQC1 > ${logs}/sexqc1.log 2>&1

### Relationship QC
${relqc}/relatedqc_FINAL.bed: ${data_derived}/${nickname}_rawE.bed ${sexqc}/SEXQC1/SEXCHECK1/sex_anomalies.txt
	@echo "running relationship qc pipeline"
	bash ${RECIPES}/relatedqc2.sh ${data_derived}/${nickname}_rawE > ${logs}/relatqc2.log 2>&1

### Sex QC2
${sexqc}/SEXQC2/SEXCHECK2/sex_anomalies.txt: ${relqc}/relatedqc_FINAL.bed 
	@echo "running sex qc pipeline (round 2)"
	bash ${RECIPES}/sexqc2.sh ${relqc}/relatedqc_FINAL SEXQC2 > ${sexqc}/sexqc2.log 2>&1

### Update data set with SEXQC2 results
${data_derived}/prelimqcC.bed: ${sexqc}/SEXQC2/SEXCHECK2/sex_anomalies.txt
	@echo "update based on sex qc round 2 results"
	bash ${RECIPES}/update_with_relqc_and_sexqc_lists.sh > ${logs}/update_with_relqc_and_sexqc_lists.log 2>&1

### Add cohort ids to fids
${data_derived}/prelimqc_FINAL.bed: ${data_derived}/prelimqcC.bed
	@echo "adding cohort to fids"
	bash ${RECIPES}/add_cohort_to_fids.sh > ${logs}/add_cohort_to_fids.log 2>&1
		
### Auto QC and Mendel Error check
${clean}/${nickname}_b37_clean.bed: ${data_derived}/prelimqc_FINAL.bed
	@echo "running auto qc"
	bash ${RECIPES}/autoqc.sh > ${logs}/autoqc.log 2>&1

### Visualize remaining pedigree errors
${relqc}/KING_MIerrorplot.pdf: ${clean}/${nickname}_b37_clean.bed
	@echo "visualizing remaining pedigree errors"
	king22 -b ${clean}/${nickname}_b37_clean.bed --prefix ${relqc}/KING --related --degree 2 --rplot 

### Sex QC3
${sexqc}/SEXQC3/SEXCHECK2/sex_anomalies.txt: ${clean}/${nickname}_b37_clean.bed
	@echo "running sex qc pipeline (round 3 - final)" 
	bash ${RECIPES}/sexqc2.sh ${clean}/${nickname}_b37_clean SEXQC3 > ${sexqc}/sexqc3.log 2>&1 
	
### Summarize QC statistics
${clean}/${nickname}_b37_clean.qcreport: ${clean}/${nickname}_b37_clean.bed ${sexqc}/SEXQC3/SEXCHECK2/sex_anomalies.txt
	@echo "summarizing qc steps"
	bash  ${RECIPES}/report_qc_statistics.sh > ${clean}/${nickname}_b37_clean.qcreport 
	
### Clean up [ only if explicitly requested ]
clean:	
	@echo "cleaning up"
	bash ${RECIPES}/cleanup.sh ${preqc_folder}/tmp_trash
	