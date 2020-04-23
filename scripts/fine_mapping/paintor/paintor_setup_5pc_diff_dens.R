#paintor_setup_5pc_diff_dens.R
#create the input files for paintor and the bash scripts to run the regions associatied on the iChip.
library(parallel)
library(GUESSFM)

extension="final_collection"
orig<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/pmeta_fam_ind_risk_regions_vcf_5pc_diff.txt"),header=T,sep=",",as.is=T)
orig<-orig[orig$ichip=="yes",]
#remove duplicate due to long LD block:
orig<-orig[!(orig$Marker %in% c("chr12:112153882:G:A","chr12:112730563:A:G","chr7:50406053:G:A",
"chr1:113302202:A:C","chr17:46751565:G:A","chr14:68286876:C:T","chr12:111066229:T:C",
"chr17:40623690:G:A")),]
orig$snpnum=c(1:nrow(orig))

#load the indepedent case-control results files to see which SNPs were included/excluded in that analysis (should be consistent with fine mapping.
r1<-get(load(file=paste0("/scratch/ji3x/",extension,"_gwas/resultsmeta_5pcs_vcf_5pc_diff.RData")))
res<-r1[substr(r1$Direction,1,1)!="?",]


#the below function wont work unless guessfm_setup has been run, relies on files created in that program:
getreg<-function(index,snpnum){

#examining those regions with PP>0.5 of 1 causal variant in the region in GUESSFM
mydir <-paste0("/scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/",gsub(":",".",index),"/")
dd <- read.snpmod(mydir)
pp <- pp.nsnp(dd,plot=FALSE,expected=3, overdispersion = 1.00000001)
stop1<-(pp$trait[2]<0.5)
if(stop1==T){
message(paste0(index," - not proceeding as pp 1 causal variant<0.5"))
out<-NULL
return(out)
}
if(stop1==F){
message(paste0(index," - proceeding as pp 1 causal variant>=0.5"))
#first job is to create z score files for each cohort that has a SNP with an association of at least < 5 * 10^-4 at the index SNP (so not prioritising noise)
getzs<-function(cohort,cohortnum,index){
l<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/snps_",gsub(":",".",index),".txt"),as.is=T)
l<-data.frame(V1=l[!l$V1 %in% c("chr22:29825212:G:A", "chr22:29874081:C:G","chr22:30226880:A:T","chr22:30223185:T:G"),])
reses<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/META/META/METAANALYSIS_",cohort,"_dens.TBL"),header=T)
reses<-reses[reses$MarkerName %in% l$V1,]
reses$Allele1<-toupper(reses$Allele1)
reses$Allele2<-toupper(reses$Allele2)
reses$A1<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\3",reses$MarkerName)
reses$A2<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\4",reses$MarkerName)
w<-which(reses$Allele1==reses$A2 & reses$Allele2==reses$A1)
reses[w,"Effect"]<-reses[w,"Effect"]*-1
reses[,paste0("ZSCORE.P",cohortnum)]<-reses$Effect/reses$StdErr
reses[,paste0("p_",cohortnum)]<-reses$P.value
mi<-min(reses[reses$MarkerName==index,paste0("p_",cohortnum)],na.rm=T)
reses<-reses[,c("MarkerName",paste0("ZSCORE.P",cohortnum))]
#take those included in the guessfm analysis.
reses<-reses[reses$MarkerName %in% res$MarkerName,]
if(mi>5*10^-4 | is.na(mi)){
reses<-NULL
}
return(reses)
}
k<-mapply(getzs,cohort=c("EUR","AFR","FIN"),cohortnum=c(1,2,3),index=c(index,index,index),SIMPLIFY=F)
#take those included in the guessfm analysis.


if(is.null(k[[2]]) & is.null(k[[3]])){
ks<-k[[1]]
}
if(is.null(k[[1]]) & is.null(k[[3]])){
ks<-k[[2]]
}
if(is.null(k[[1]]) & is.null(k[[2]])){
ks<-k[[3]]
}
stop<-(is.null(k[[2]]) & is.null(k[[3]]) |
is.null(k[[1]]) & is.null(k[[3]]) |
is.null(k[[1]]) & is.null(k[[2]]))
if(stop==T){
message(paste0(index," - not proceeding as only association in 1 ancestry group"))
out<-NULL
return(out)
}
if(stop==F){
message(paste0(index," - proceeding as association in >1 ancestry group"))
if (!is.null(k[[1]]) & !is.null(k[[2]])){
ks<-merge(k[[1]],k[[2]],by="MarkerName")
}
if(is.null(k[[2]]) & !is.null(k[[1]]) & !is.null(k[[3]])){
ks<-merge(k[[1]],k[[3]],by="MarkerName")
}
if(is.null(k[[1]]) & !is.null(k[[2]]) & !is.null(k[[3]])){
ks<-merge(k[[2]],k[[3]],by="MarkerName")
}
if(!is.null(k[[1]]) & !is.null(k[[2]]) & !is.null(k[[3]])){
ks<-merge(ks,k[[3]],by="MarkerName")
}
ks$chr<-gsub(":.*","",ks$MarkerName)
ks$position<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\2",ks$MarkerName)
zs<-colnames(ks)[substr(colnames(ks),1,2)=="ZS"]
ks$position<-as.numeric(ks$position)
ks<-ks[order(ks$position),]
ks<-as.data.frame(ks[,c("chr","position","MarkerName",zs)])
#KEEP ONLY THOSE INCLUDED IN GUESSFM ANALSYSIS:
dat<-get(load(file=paste0("/scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/",gsub(":",".",index),"/data.RData"))[[2]])
ks<-ks[ks$MarkerName %in% colnames(dat),]
if(index=="chr6:126343187:C:T"){
ks<-ks[c(1:700),]
}
if(index=="chr17:39910119:T:C"){
ks<-ks[ks$position<40500000,]
}
if(index=="chr17:39910119:T:C"){
ks<-ks[ks$position<56200000,]
}
write.table(ks,file=paste0("/scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/Locus",snpnum), sep=" ", 
col.names=T, row.names=F, quote=F)

#and to keep these SNPs only for the LD calculation:
kn<-ks[,c("chr","position")]
write.table(kn,file=paste0("/scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/",index,".keep"),col.names=F,row.names=F,quote=F,sep="\t")

#and in case of multi allelic variants in some populations:
ks$SNPID<-ks$MarkerName
ks$rsid<-ks$MarkerName
ks$a1<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\3",ks$rsid)
ks$a2<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\4",ks$rsid)
kn1<-ks[,c("SNPID","rsid","chr","position","a1","a2")]
colnames(kn1)<-c("SNPID","rsid","chromosome","position","alleleA","alleleB")
write.table(kn1,file=paste0("/scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/",gsub(":",".",index),".keep1"),col.names=T,row.names=F,quote=F,sep=" ")

#now for the LD estimation - going to use LDSTORE to do this quickly and efficiently
#convert the VCF file to bgen:
#which cohorts to run this for:
w<-c(!is.null(k[[1]]),!is.null(k[[2]]),!is.null(k[[3]]))
df<-data.frame(cohort=c("EUR","AFR","FIN"),cohortnum=c(1,2,3),do=w)
cohorts<-df[df$do==TRUE,"cohort"]
cohorts<-as.character(cohorts)
cohortnums<-df[df$do==TRUE,"cohortnum"]
chr<-gsub("chr","",ks$chr[1])
#produce dosage files to convet to bgen for LD calculation:
genvcf<-function(cohort){
system(paste0("/home/ji3x/software/bcftools/bcftools view -R /scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/",index,
".keep -o /scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/",cohort,
"/",gsub(":",".",index),".vcf -O v /nv/vol185/MEGA/release4/IMPUTED_TOPMED/results_filtered/",cohort,"/unrelateds/chr",chr,
".filter_maf_gt_0.005_rsq_gt_0.8.unrelateds.dose.vcf.gz"))
system(paste0("/home/ji3x/software/bcftools/bcftools view /scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/",cohort,
"/",gsub(":",".",index),".vcf | awk \'/^#/{print}; !/^#/{if (!uniq[$3]++) print}\' > /scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/",cohort,
"/",gsub(":",".",index),"_nodup.vcf"))
}
lapply(cohorts,genvcf)


#now for the LD estimation - going to use LDSTORE to do this quickly and efficiently
#convert the VCF file to bgen:

genbgen<-function(cohort){
system(paste0("/home/ji3x/software/gavinband-qctool-ba5eaa44a62f/bin/qctool_v2.0.1 -g /scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/",cohort,
"/",gsub(":",".",index),"_nodup.vcf -filetype vcf -og /scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/",cohort,
"/",gsub(":",".",index),".bgen -ofiletype bgen"))
system(paste0("/home/ji3x/software/gavinband-qctool-ba5eaa44a62f/bin/qctool_v2.0.1 -g /scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/",cohort,
"/",gsub(":",".",index),".bgen -incl-variants /scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/",gsub(":",".",index),".keep1 -og /scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/",cohort,
"/",gsub(":",".",index),"_filt.bgen"))
}
lapply(cohorts,genbgen)

genbcor<-function(cohort,cohortnum,snpnum){
system(paste0("~/software/ldstore_v1.1_x86_64/ldstore --bgen /scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/",cohort,
"/",gsub(":",".",index),"_filt.bgen --bcor /scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/",cohort,
"/",gsub(":",".",index),".bcor --ld-thold 0 --n-threads 1"))
system(paste0("~/software/ldstore_v1.1_x86_64/ldstore --bcor /scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/",cohort,
"/",gsub(":",".",index),".bcor_1 --matrix /scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/Locus",snpnum,".LD",cohortnum))
}
mcmapply(genbcor,cohort=cohorts,cohortnum=cohortnums, snpnum=c(rep(snpnum,nrow(df[df$do==TRUE,]))),SIMPLIFY=F, mc.cores=8)

#the annotation file (all 0s):
annot<-data.frame(dummy=c(rep(0,nrow(ks))))
write.table(annot, file=paste0("/scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/Locus",snpnum,".annotations"), quote=F, col.names=T, row.names=F)

#output file to help generate the locus file PAINTOR needs to run:
out<-data.frame(snp=index, locusnum=snpnum)
out$inc1<-(w[[1]]==TRUE)
out$inc2<-(w[[2]]==TRUE)
out$inc3<-(w[[3]]==TRUE)
message(paste0("done ",index))
}
}
return(out)
}

forscript<-mcmapply(getreg, index=orig$Marker, snpnum=orig$snpnum,SIMPLIFY=F,mc.cores=4)
forscript<-do.call("rbind",forscript)
save(forscript, file=paste0("/scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/forscript.RData"))

load(file=paste0("/scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/forscript.RData"))

#finally the shell script to run PAINTOR (through the cluster):
combsrun<-function(inc1,inc2,inc3,name,filename,zs,lds){
f<-forscript[forscript$inc1==inc1 & forscript$inc2==inc2 & forscript$inc3==inc3,]
f$l<-paste0("Locus",f$locusnum)
write.table(f$l,file=paste0("/scratch/ji3x/",extension,"_gwas/input_paintor_5pc_diff_",filename),quote=F,col.names=F,row.names=F)


#having to do a number of scripts depending on how many of the cohorts have a SNP with some association signal:
sink(file=paste0("~/programs/",extension,"/paintorscripts/",name,"_5pc_diff.sh"))
cat(paste0("#!/bin/bash
#SBATCH -A rich_immunochip_impute
#SBATCH --time=12:00:00
#SBATCH -p standard
#SBATCH -n 1
#SBATCH --mem=18000


module load gcc gsl/2.4

"))

cat(paste0("~/software/paintor/PAINTOR_V3.0/PAINTOR -input /scratch/ji3x/",extension,"_gwas/input_paintor_5pc_diff_",filename," ", 
"-in /scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/ -out /scratch/ji3x/",extension,"_gwas/paintor_results_5pc_diff/ ",
"-Zhead ",zs," -LDname ",lds," -enumerate 2\n"))
sink()
}
if(nrow(forscript[forscript$inc1==T & forscript$inc2==T & forscript$inc3==T,])>0){
combsrun(inc1=T,inc2=T,inc3=T, name="three",filename="3",zs="ZSCORE.P1,ZSCORE.P2,ZSCORE.P3",lds="LD1,LD2,LD3")
}
if(nrow(forscript[forscript$inc1==T & forscript$inc2==T	& forscript$inc3==F,])>0){
combsrun(inc1=T,inc2=T,inc3=F, name="oneandtwo",filename="1_and_2",zs="ZSCORE.P1,ZSCORE.P2",lds="LD1,LD2")
}
if(nrow(forscript[forscript$inc1==T & forscript$inc2==F & forscript$inc3==T,])>0){
combsrun(inc1=T,inc2=F,inc3=T, name="oneandthree",filename="1_and_3",zs="ZSCORE.P1,ZSCORE.P3",lds="LD1,LD3")
}


