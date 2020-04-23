#cond_setup_5pc_diff_dens.R

#set up the conditional analyses by producing subsets of VCFs to perform conditional regressions with.
library(parallel)
extension="final_collection"
vcfpath="/nv/vol185/MEGA/release4/IMPUTED_TOPMED/results_filtered/"
cohorts=c("AMR","EUR","AFR","FIN")
outputdir="~/output/final_collection_gwas/META/"
path=paste0("/scratch/ji3x/",extension,"_gwas/vcf/conditional/")

signals<-read.table(file=paste0(outputdir,"/pmeta_5pcs_vcf_5pc_diff_dens.txt"),as.is=T, sep=",",header=T)
#remove those there due to LD (think just one at MAPT):
signals$chromosome<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\1",signals$Marker)
signals$chromosome<-gsub("chr","",signals$chromosome)
signals$chromosome<-as.numeric(signals$chromosome)
signals$position<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\2",signals$Marker)
signals$position<-as.numeric(signals$position)


#use BCF tools to keep only the SNPs of interest
ext="final_collection"
getreg<-function(snp1, chr, min, max, cohort){
sink(file=paste0("~/programs/",ext,"/META/conditional/",cohort,"/",snp1,"_get.sh"))
cat(paste0("#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH -A rich_immunochip_impute
#SBATCH -p standard
#SBATCH -n 1
#",snp1,"_get.sh\n\n"))
cat(paste0("/home/ji3x/software/bcftools/bcftools view -r chr",chr,
":",min,"-",max," ",vcfpath, cohort,"/unrelateds/chr",chr,".filter_maf_gt_0.005_rsq_gt_0.8.unrelateds.dose.vcf.gz > ",
path,cohort,"/",snp1,"out.vcf\n"))
cat(paste0("bgzip -f ",path,cohort,"/",snp1,"out.vcf\n"))
cat(paste0("tabix -f ",path,cohort,"/",snp1,"out.vcf.gz\n"))
sink()
message(paste0("getting region for ",snp1," - ",cohort))
system(paste0("sbatch ~/programs/",ext,"/META/conditional/",cohort,"/",snp1,"_get.sh"))
}



f<-function(snp){
s<-signals[signals$Marker==snp,]
chr<-s$chromosome
min<-max(c(0,(s$position-750000)))
max<-s$position+750000
pos<-s$position
snp1<-gsub(":","_",snp)
lapply(cohorts, getreg, snp1=snp1, chr=chr, min=min, max=max)
}
lapply(signals$Marker, FUN=f)


