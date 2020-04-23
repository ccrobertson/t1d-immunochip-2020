#guessfm_readin_script_gen.R

#generating scripts to run guessfm_readin at each locus:

extension="final_collection"
orig<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/pmeta_fam_ind_risk_regions_vcf_5pc_diff.txt"),header=T,sep=",",as.is=T)
orig<-orig[orig$ichip=="yes",]
orig$mem<-ifelse(orig$Marker %in% c("chr10:6052734:C:T","chr4:122322439:A:ATC","chr11:2160994:A:T"),100000,
ifelse(orig$Marker %in% c("chr2:203874196:G:A"),200000,
ifelse(orig$Marker %in% c("chr16:11097543:G:A","chr17:39910119:T:C"),48000,24000)))
orig$time<-ifelse(orig$Marker %in% c("chr10:6052734:C:T","chr4:122322439:A:ATC"),"24:00:00",
ifelse(orig$Marker %in% c("chr2:203874196:G:A"),"7-00:00:00",
ifelse(orig$Marker %in% c("chr16:11097543:G:A"),"16:00:00",
ifelse(orig$Marker %in% c("chr11:2163618:G:T"),"15:00:00","8:00:00"))))


genscript<-function(snp, time, mem){
sink(file=paste0("~/programs/",extension,"/guessscripts/readin/",gsub(":",".",snp),"_5pc_diff.sh"))
n<-which(orig$Marker==snp)
cat(paste0("#!/bin/bash
#SBATCH --time=",time,"
#SBATCH -A rich_immunochip_impute
#SBATCH -p standard
#SBATCH -n 1
#SBATCH --mem=",format(mem,scientific=F),"

module load gcc R/3.5.1
"))

cat(paste0("Rscript ~/programs/",extension,"/guessfm_readin_5pc_diff.R ",n,"\n"))
sink()
system(paste0("sbatch ~/programs/",extension,"/guessscripts/readin/",gsub(":",".",snp),"_5pc_diff.sh"))
}

mapply(genscript, snp=orig$Marker,time=orig$time, mem=orig$mem)
