#guessfm_setup_5pc_diff.R
#generate region VCFs
#only for hits in iChip regions, as agreed with Cassie and runs the GUESS stochastic search in the region

library(vcfR)
library(reshape2)
library(GUESSFM)
library(parallel)
library(filesstrings)
library(snpStats)
library(GenomicRanges)
library(rtracklayer)


extension="final_collection"
vcfpath<-"/nv/vol185/MEGA/release4/IMPUTED_TOPMED/results_filtered/EUR/unrelateds/"
orig<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/pmeta_fam_ind_risk_regions_vcf_5pc_diff.txt"),header=T,sep=",",as.is=T)
orig<-orig[orig$ichip=="yes",]

#load the indepedent case-control results files to see which SNPs were included/excluded in that analysis (should be consistent with fine mapping).
load(file=paste0("/scratch/ji3x/",extension,"_gwas/resultsmeta_5pcs_vcf_5pc_diff.RData"))
#looking at europeans only, so just keeoping SNPs that passed european QC:
res<-res[substr(res$Direction,1,1)!="?",]


#get ichip reginos defined as only examining variants in these densely genotyped regions (otherwise could be misleading as 
#tag variants for a haplotype can have very high posterior probability of causailty due to absense of other SNPs on htat haplotype in the
#set of SNPs under consideration).
ichip<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/rawdat/ichip_regions_gr38.txt"),header=T,as.is=T)

chipreg<-GRanges(seqnames=ichip$chromosome,
IRanges(ichip$start, end=ichip$end),
region=ichip$region)

#take a 1.5 megabase region around the index SNPs:
getreg<-function(index){
o<-orig[orig$Marker==index,]
chr<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\1",o$Marker)
chr<-as.numeric(gsub("chr","",chr))
position<-as.numeric(sub("(^.*)[:](.*)[:](.*)[:](.*)","\\2",o$Marker))
min<-position-750000
max<-position+750000
r<-res[res$chromosome==chr & res$position>min & res$position<max,]
#remove those not in iChip regions:
rreg<-GRanges(seqnames=paste0("chr",r$chromosome),
IRanges(r$position, end=r$position), MarkerName=r$MarkerName)
rr<-mergeByOverlaps(chipreg,rreg)
r1<-rr$MarkerName
write.table(r1,file=paste0("/scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/snps_",gsub(":",".",index),".txt"), col.names=F, row.names=F, quote=F)


#create file that SnpStats can read into R for use in GUESSFM from VCF file:
system(paste0("/home/ji3x/software/bcftools/bcftools view -r chr",chr,":",min,"-",max," -o /scratch/ji3x/",extension,
"_gwas/guessfm_5pc_diff/input/",gsub(":",".",index),".vcf -O v ",vcfpath,"/chr",chr,
".filter_maf_gt_0.005_rsq_gt_0.8.unrelateds.dose.vcf.gz"))

columns<-system(paste0("tail -n 1 /scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/",gsub(":",".",index),".vcf | awk -F \'\\t\' \'{print NF;exit}\'"),intern=T)

system(noquote(paste0("awk \'{for (i=10;i<=",columns,";i++)sub(/.+:/,\"\",$i);print}\' /scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/",
gsub(":",".",index),".vcf > /scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/try",gsub(":",".",index))))
system(noquote(paste0("awk -F \' \' \'{gsub(/,/,\" \");print}\' /scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/try",
gsub(":",".",index)," >/scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/try",gsub(":",".",index),"_1")))
system(noquote(paste0("cat /scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/try",
gsub(":",".",index),"_1 | sed \'/^#/ d\' > /scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/try",gsub(":",".",index),"_2")))
system(noquote(paste0("awk -F \' \' \'{t=$3; s=$2; $1=t; $2=t; $3=s;print}\' /scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/try",
gsub(":",".",index),"_2 > /scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/try",gsub(":",".",index),"_3")))
system(noquote(paste0("awk -F \' \' \'{$6=$7=$8=$9=\"\";print}\'  /scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/try",
gsub(":",".",index),"_3 > /scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/see_",gsub(":",".",index),".impute")))
system(noquote(paste0("rm /scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/try",gsub(":",".",index),"*")))


system(paste0("bgzip -f /scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/",gsub(":",".",index),".vcf"))
system(paste0("tabix -f /scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/",gsub(":",".",index),".vcf.gz"))
t1<-system(paste0("tabix -H /scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/",gsub(":",".",index),".vcf.gz"),intern=T)
t1<-t1[[length(t1)]]
t1<-strsplit(t1,split="\t")
t1<-t1[[1]]
t1<-t1[10:length(t1)]
DATA<-read.impute(file=paste0("/scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/see_",gsub(":",".",index),".impute"), rownames=t1)

#keep only high quality variants:
DATA<-DATA[,colnames(DATA) %in% r1]

#Get phenotype data:
init<-read.table(file="/nv/vol185/MEGA/mega_genotyped/release4/mega_b37_clean_pheno_genotypedOnly.txt", as.is=T)
pcs<-read.table(file="/nv/vol185/MEGA/release4/pca/mega_pca_b37_pruned_unrelated_EUR_pca_proj_onto_controls.sscore",as.is=T)
rownames(pcs)<-pcs$V2
rownames(init)<-init$V3
pcs<-pcs[t1,]
init<-init[t1,]
pcs<-pcs[,c(6,7,8,9,10)]
colnames(pcs)<-c("PC1","PC2","PC3","PC4","PC5")
init<-init[,c(2,7)]
colnames(init)<-c("ID","t1d")
init$ID<-rownames(init)
all<-cbind(init,pcs)
all<-all[all$ID %in% t1,]
all<-all[!duplicated(all$ID),]
all$t1d<-ifelse(all$t1d==2,1,ifelse(all$t1d==1,0,NA))
all<-all[!is.na(all$t1d),]
DATA<-DATA[rownames(DATA) %in% all$ID,]
DATA<-DATA[all$ID,]
Y<-data.frame(outcome=all$t1d)
covariates<-data.frame(pc1=all$PC1, pc2=all$PC2,
pc3=all$PC3, pc4=all$PC4, pc5=all$PC5)
rownames(Y)<-rownames(all)
rownames(covariates)<-rownames(all)

#remove missing
all<-all[!is.na(all$t1d),]
DATA<-DATA[all$ID,]

#filter poorly imputed SNPs:
message("Removing poorly imputed SNPs.")
message("Input matrix has ",ncol(DATA)," SNPs.")
cs <- col.summary(DATA)
wh <- which(cs[,"MAF"]<0.005 | cs[,"Call.rate"]<0.99 | cs[,"Certain.calls"]<0.75 | abs(cs[,"z.HWE"])>20 | is.na(cs[,"z.HWE"]))
if(length(wh)) {
  message("Dropping ",length(wh)," SNPs with |z.HWE|>20, MAF < 0.005 or call rate <0.99")
  DATA <- DATA[,-wh]
}



#Run bayesian variable selection via GUESS
dir.create(file.path(paste0("/scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/",gsub(":",".",index),"/")), showWarnings = FALSE)
mydir <-paste0("/scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/",gsub(":",".",index),"/")
save(Y, DATA, covariates, file=paste0(mydir, "data.RData"))
load(file=paste0(mydir, "data.RData"))


com<-run.bvs(X=DATA,Y=Y[,"outcome"],gdir=paste0(mydir),
        guess.command="/home/ji3x/software/GUESS_v1.1/Main/GUESS",
        nexp=3,                # expected number of causal variants, an overestimate
        nsave=10000,            # number of best models to save
        tag.r2=0.99,            #R2 tag value
        family="binomial",
        covars=covariates, run=FALSE)

sink(file=paste0("~/programs/",extension,"/guessscripts/",gsub(":",".",index),"_5pc_diff.sh"))
cat(paste0("#!/bin/bash
"))
if(index=="chr10:6052734:C:T"){
cat(paste0("#SBATCH --time=12:00:00
"))
}
if(index!="chr10:6052734:C:T"){
cat(paste0("#SBATCH --time=8:00:00
"))
}
cat(paste0("#SBATCH -A rich_immunochip_impute
#SBATCH -p standard
#SBATCH -n 1
#SBATCH --mem=18000

module load gcc 

#",gsub(":",".",index),".sh
",com,"\n"))
sink()
system(paste0("sbatch ~/programs/",extension,"/guessscripts/",gsub(":",".",index),"_5pc_diff.sh"))
}

mclapply(orig$Marker[1:nrow(orig)], getreg, mc.cores=4)



