#post_imputation_qc_5pc_diff_dens.R

#This script performs post-imputation QC on the set of analysed SNPs:
#then generates METAL scripts to run the fixed effects meta analysis 
library(GenomicRanges)

#format for each cohort, 
extension="final_collection"
cohorts=c("AFR","EUR","AMR", "FIN")
outdir<-paste0("/scratch/ji3x/",extension,"_gwas/META/")
metdir<-paste0("/scratch/ji3x/",extension,"_gwas/META/META/")

#define ichip regions:
ichip<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/rawdat/dens_and_ichip_regs.txt"),header=T,as.is=T,sep="\t")
ichip<-GRanges(seqnames=ichip$seqnames,
IRanges(ichip$start, end=ichip$end))

#and add the SNPs that were on the iChip to this set:
ichip1<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/rawdat/topmed_imputation_input_variants_b38_sorted_excludingMismatches.txt"),header=F,sep=" ",as.is=T)
ichip1<-GRanges(seqnames=paste0("chr",ichip1$V1),
IRanges(ichip1$V2,end=ichip1$V2))


k<-setdiff(ichip1,ichip)
allichip<-c(ichip, k)

combineallsnps<-function(cohort){
l<-lapply(c(1:22),combinethesnps,cohort=cohort)
l<-do.call("rbind",l)
return(l)
}

combinethesnps<-function(cohort, chrom){
r<-read.table(file=paste0(outdir,"/",cohort,"/chr_",chrom), header=T, as.is=T)
cols<-colnames(r)
r$allele1<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\3",r$rsid)
r$allele2<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\4",r$rsid)
r$chromosome<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\1",r$rsid)
r$chromosome<-gsub("chr","",r$chromosome)
r$chromosome<-as.numeric(r$chromosome)

#ensure the effect directions are all alligned to TOPMED reference panel:
w<-which(r$alleleA==r$allele2 & r$alleleB==r$allele1)
r[w,"aAdummy"]<-r[w,"alleleA"]
r[w,"aBdummy"]<-r[w,"alleleB"]
r[w,"alleleA"]<-r[w,"aBdummy"]
r[w,"alleleB"]<-r[w,"aAdummy"]
r[w,"frequentist_add_beta_1.add.t1d.1"]<-r[w,"frequentist_add_beta_1.add.t1d.1"]*-1
r[w,"abdummy"]<-r[w,"controls_AB"]
r[w,"aadummy"]<-r[w,"controls_AA"]
r[w,"bbdummy"]<-r[w,"controls_BB"]
r[w,"controls_BB"]<-r[w,"aadummy"]
r[w,"controls_AA"]<-r[w,"bbdummy"]
r$AF<-(r$controls_AB+(2*r$controls_BB))/((2*r$controls_AA) + 2*r$controls_AB + (2*r$controls_BB))

#remove the results of those with MAF<0.005:
message(paste0("Removing SNPs with MAF<0.005"))
r[r$all_maf<0.005,"frequentist_add_beta_1.add.t1d.1"]<-NA
r[r$all_maf<0.005,"frequentist_add_wald_pvalue_1"]<-NA
r[r$all_maf<0.005,"frequentist_add_se_1"]<-NA
r[r$all_maf<0.005,"all_info"]<-NA

#remove the SNPs with imputation information score<0.8:
message(paste0("Removing SNPs with info score<0.8 (overall or in cases or controls"))
w1<-which(r$all_info<0.8 | r$cases_info<0.8 | r$controls_info<0.8)
r[w1,"frequentist_add_beta_1.add.t1d.1"]<-NA
r[w1,"frequentist_add_wald_pvalue_1"]<-NA
r[w1,"frequentist_add_se_1"]<-NA
r[w1,"all_info"]<-NA

#remove SNPs with imputation score difference >5% between cases and controls:
message(paste0("Removing SNPs with info score >5% different between cases and controls"))
r$diff<-abs(r$cases_info-r$controls_info)
w2<-which(r$diff>0.05)
r[w2,"frequentist_add_beta_1.add.t1d.1"]<-NA
r[w2,"frequentist_add_wald_pvalue_1"]<-NA
r[w2,"frequentist_add_se_1"]<-NA
r[w2,"all_info"]<-NA
#remove imputed  SNPs outside ichip regions:
r1<-GRanges(seqnames=paste0("chr",r$chromosome),
IRanges(r$position,end=r$position),
MarkerName=r$rsid)

ot<-mergeByOverlaps(r1,allichip)
ot<-as.data.frame(ot)
w<-which(!r$rsid %in% ot$MarkerName)
r[w,"frequentist_add_beta_1.add.t1d.1"]<-NA
r[w,"frequentist_add_wald_pvalue_1"]<-NA
r[w,"frequentist_add_se_1"]<-NA
r[w,"all_info"]<-NA
r<-r[,cols]
return(r)
}
results<-lapply(cohorts, combineallsnps)
names(results)<-cohorts



#Get them into metal format and write the files out:
metalthem<-function(cohort){
co<-results[[cohort]]
co$RefAllele<-co$alleleA
co$NonRefAllele<-co$alleleB
co$BETA<-co$frequentist_add_beta_1.add.t1d.1
co$SE<-co$frequentist_add_se_1
co$SNP<-co$rsid
df<-co[,c("SNP","BETA","SE","NonRefAllele", "RefAllele")]
write.table(df,paste0(metdir,"metal_results_",cohort,"_5pc_diff_dens.tbl"), sep=" ", col.names=T, row.names=F, quote=F)
return(df)
}
r<-lapply(cohorts,metalthem)



#now generate a file to use in METAL:
generatefile<-function(extension){
sink(file=paste0("~/programs/",extension,"/META/metal_script_5pc_diff_dens"))
cat(paste0("#metal_script_5pc_diff_dens

#THIS SCRIPT EXECUTES AN ANALYSIS OF THE SEVEN INDEPENDENT GROUPS

#LOAD THE SEVEN INPUT FILES

# === DESCRIBE AND PROCESS THE FIRST INPUT FILE ===
MARKER SNP
ALLELE RefAllele NonRefAllele
EFFECT BETA
STDERR SE 
SCHEME STDERR 
"))

cat(paste0("
PROCESS ",metdir,"/metal_results_EUR_5pc_diff_dens.tbl 

#All files are of the same format so can process immediately:
PROCESS ",metdir,"/metal_results_AMR_5pc_diff_dens.tbl
PROCESS ",metdir,"/metal_results_FIN_5pc_diff_dens.tbl
PROCESS ",metdir,"/metal_results_AFR_5pc_diff_dens.tbl

OUTFILE \"",metdir,"/META/METAANALYSIS_initial_5pc_diff_dens\"
ANALYZE

QUIT"))
sink()
}
generatefile("final_collection")


#finally the file to run the metal script:

generatescript<-function(extension){
sink(file=paste0("~/programs/",extension,"/META/run_metal_5pc_diff_dens.sh"))
cat(paste0("#!/bin/bash
#run_metal_5pc_diff_dens.sh
"))

cat(paste0("/home/ji3x/software/generic-metal/metal metal_script_5pc_diff_dens

#Currently can't get the results to output to where i'd like them to,
#so I am moving here:
mv METAANALYSIS1.TBL METAANALYSIS_5pc_diff_dens.TBL
mv METAANALYSIS1.TBL.info METAANALYSIS_5pc_diff_dens.info

mv METAANALYSIS_5pc_diff_dens.TBL ",metdir,"
mv METAANALYSIS_5pc_diff_dens.info ",metdir,"
"))
sink()
system(paste0("chmod a=rwx ~/programs/",extension,"/META/run_metal_5pc_diff_dens.sh"))
}
generatescript("final_collection")


