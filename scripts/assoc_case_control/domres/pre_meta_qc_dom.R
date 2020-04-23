#pre_meta_qc_dom.R

#this script generates shell scripts to perform association tests for each different cohort, directly from the VCF file 
#and taking into account the uncertainty in the genotype imputation.
#Testing the dominance model for each SNP

extension="final_collection"
cohorts=c("AMR","EUR","AFR","FIN")
vcfpath="/nv/vol185/MEGA/release4/IMPUTED_TOPMED/results_filtered/"
dobbypath=paste0("/scratch/ji3x/",extension,"_gwas/vcf/")
outdir<-paste0("/scratch/ji3x/",extension,"_gwas/META/")

#Now lets generate shell scripts that will test for association for each cohort:
ext="final_collection"
dothescripts<-function(cohort, chr){
sink(file=paste0("~/programs/",ext,"/META/snptest_scripts/dom/",cohort,"_chr",chr))
cat(paste0("~/software/snptest_v2.5.4-beta3_linux_x86_64_static/snptest_v2.5.4-beta3 -data ",
vcfpath,cohort,"/unrelateds/chr",chr,".filter_maf_gt_0.005_rsq_gt_0.8.unrelateds.dose.vcf.gz ",dobbypath,
cohort,"/samps_snptest -filetype vcf -include_samples ",dobbypath,cohort,
"/samples.txt -genotype_field GP -pheno t1d -cov_all_continuous -frequentist 2 -method newml -o ",outdir,
cohort,"/dom/chr_",chr))
sink()
sink(file=paste0("~/programs/",ext,"/META/snptest_scripts/dom/",cohort,"_chr",chr,".sh"))
cat(paste0("#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH -A rich_immunochip_impute
#SBATCH -p standard
#SBATCH -n 1

bash ~/programs/",ext,"/META/snptest_scripts/dom/",cohort,"_chr",chr,"
"))
sink()

system(paste0("sbatch ~/programs/",ext,"/META/snptest_scripts/dom/",cohort,"_chr",chr,".sh"))
}


lapply(c(1:22),dothescripts,cohort="EUR")
lapply(c(1:22),dothescripts,cohort="FIN")
lapply(c(1:22),dothescripts,cohort="AMR")
lapply(c(1:22),dothescripts,cohort="AFR")



