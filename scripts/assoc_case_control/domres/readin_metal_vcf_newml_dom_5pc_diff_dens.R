#readin_metal_vcf_newml_dom_5pc_diff_dens.R

#reads in results from initial meta run:

library(snpStats)
library(annotSnpStats)
library(qqman)
#change the name "extension" below in future analyses to generate same output in different location.

extension="final_collection"
cohorts=c("AFR","EUR","AMR", "FIN")
outdir<-paste0("/scratch/ji3x/",extension,"_gwas/META/")
metdir<-paste0("/scratch/ji3x/",extension,"_gwas/META/META/")
dobbydir<-paste0("/scratch/ji3x/",extension,"_gwas/vcf/")
outputdir<-"~/output/final_collection_gwas/META/"


res<-read.table(file=paste0(metdir,"/dom/METAANALYSIS_5pcs_vcf_5pc_diff_dens.TBL"), header=T, as.is=T)

#allign to topmed:
res$a0<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\3",res$MarkerName)
res$a1<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\4",res$MarkerName)
res$Allele1<-toupper(res$Allele1)
res$Allele2<-toupper(res$Allele2)
w<-which(res$Allele1==res$a1)
res[w,"Effect"]<-res[w,"Effect"]*-1
res$Allele1<-res$a0
res$Allele2<-res$a1


res$snp.name<-res$MarkerName
res$chisq <- qchisq(1-res$P.value,1)
res$position<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\2",res$snp.name)
res$chromosome<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\1",res$snp.name)
res$chromosome<-gsub("chr","",res$chromosome)
res$position<-as.numeric(res$position)
res$chromosome<-as.numeric(res$chromosome)
res$snp.name=res$MarkerName
res$z<-res$Effect/res$StdErr

#remove MHC and INS regions then calculate lambda1000:
nonsig<-res[!(res$chromosome==6 & res$position>25000000 & res$position<35000000) & !(res$chromosome==11 & res$position>2000000 & res$position<3000000),]
lam<-median(nonsig[,"chisq"], na.rm=T)/qchisq(0.5,1)

readitin<-function(name){
samps<-read.table(file=paste0(dobbydir,"/",name,"/samples.txt"), header=F, as.is=T)
all<-read.table(file="/nv/vol185/MEGA/mega_genotyped/release4/mega_b37_clean_pheno_genotypedOnly.txt", header=F, as.is=T)
r<-all[all$V2 %in% samps$V1,]
t<-table(r[,6])
r<-data.frame(n0=t[1],n1=t[2])
return(r)
}
reads<-lapply(cohorts, readitin)
reads<-do.call("rbind",reads)
n0<-sum(reads[,1])
n1<-sum(reads[,2])


lam1000<-1+(lam-1)*((1/n0) + (1/n1))/(2/1000)

save(res,file=paste0(metdir,"/dom/resultsmeta_5pcs_vcf_5pc_diff_dens.RData"))

load(file=paste0(metdir,"/dom/resultsmeta_5pcs_vcf_5pc_diff_dens.RData"))
res<-res[!is.na(res$z),]
chr<-paste0("chr", res$chromosome)
start<-res$position
snp.name<-res$MarkerName

#remove MHC:
res<-res[!(res$chromosome==6 & res$position>25000000 & res$position<35000000),]
res$z<-res$Effect/res$StdErr
res$abz<-abs(res$z)
getmin<-function(frame){
min<-frame[frame$abz==max(frame$abz),]
#in situations with identical effect sizes, choose either (for now)...
min<-min[1,]
return(min)
}

dropnext<-function(min, results){
droplow<-min$position-500000
drophigh<-min$position+500000
chrom<-min$chromosome
m<-results[!(results$chromosome==chrom & results$position>droplow & results$position<drophigh),]
return(m)
}

n<-getmin(res)
dropped<-dropnext(n,res)
both<-n

while(n$P.value<5*10^-8){
n<-getmin(dropped)
dropped<-dropnext(n,dropped)
both<-rbind(both, n)
}

both<-both[order(both$chromosome, both$position),]
both$or<-exp(both$Effect)


#get the summary stats for each cohort:
getco<-function(cohort){
co<-read.table(file=paste0(metdir,"/dom/metal_results_",cohort,"_5pc_diff_dens.tbl"),
header=T, as.is=T)
rownames(co)<-co$SNP
co<-co[both$MarkerName,]
co[,paste0("a0_",cohort)]<-co$RefAllele
co[,paste0("a1_",cohort)]<-co$NonRefAllele
w<-which(co$RefAllele!=both$Allele1)
co[w,"BETA"]<-co[w,"BETA"]*-1
co[,paste0("beta_",cohort)]<-co$BETA
co<-co[,c(1,8)]
return(co)
}
indiv<-lapply(cohorts,getco)

out<-cbind(both, indiv[[1]])
out<-cbind(out, indiv[[2]])
out<-cbind(out, indiv[[3]])
out<-cbind(out, indiv[[4]])

#create regions file for BCF tools to extract SNPs needed:
o<-out[,c("chromosome","position")]
write.table(o, file=paste0("/scratch/ji3x/full_cohort_gwas/vars/signewmlr2.txt"), sep="\t",
col.names=F, row.names=F, quote=F)

system(paste0("/home/ji3x/software/bcftools/bcftools view -R /scratch/ji3x/full_cohort_gwas",
"/vars/signewmlr2.txt /scratch/ji3x/full_cohort_gwas",
"/vars/Kaviar-160204-Public-hg38-trim.vcf.gz > /scratch/ji3x/full_cohort_gwas",
"/vars//signewmlr2_ids.txt"))

vars<-read.table(file=paste0("/scratch/ji3x/full_cohort_gwas",
"/vars/signewmlr2_ids.txt"),skip=40, header=T, comment.char="",as.is=T)
vars$MarkerName<-paste0("chr",vars$X.CHROM,":",vars$POS,":",vars$REF,":",vars$ALT)

vars<-vars[vars$ID!=".",]
see<-strsplit(vars$ALT,split=",")
l<-lengths(see)
alts<-str_split_fixed(vars$ALT, ",", max(l)+1)
colnames(alts)<-paste0("alt",c(1:ncol(alts)))
vars<-cbind(vars,alts)
for(i in c(1:(max(l)+1))){
vars[,paste0("MarkerName",i)]<-paste0("chr",vars$X.CHROM,":",vars$POS,":",vars$REF,":",vars[,paste0("alt",i)])
}

vars$MarkerName<-ifelse(vars$MarkerName1 %in% out$MarkerName,vars$MarkerName1,
ifelse(vars$MarkerName2 %in% out$MakerName,vars$MarkerName2,
ifelse(vars$MarkerName3 %in% out$MarkerName,vars$MarkerName3,NA)))

vars<-vars[!is.na(vars$MarkerName),]
vars<-vars[!duplicated(vars$MarkerName),]
vars<-vars[,c("MarkerName","ID")]
out<-merge(out,vars,by="MarkerName",all.x=T)
out$MarkerName<-as.character(out$MarkerName)
out$ID<-ifelse(is.na(out$ID), out$MarkerName,out$ID)


#include additive p-value in Table:
load(file=paste0(metdir,"/resultsmeta_5pcs_vcf_5pc_diff_dens.RData"))
r<-res[res$MarkerName %in% out$MarkerName,]
r$p_add<-r$P.value
r<-r[,c("MarkerName","p_add")]
out<-merge(out, r, by="MarkerName",all.x=TRUE)
out<-out[order(out$chromosome, out$position),]


#remove those there due to LD:

out<-out[!out$Marker %in% c("chr12:112246417:A:T","chr12:112753113:T:C",
"chr17:46751565:G:A"),]

#exporting the regions to a word table:
sink(file=paste0(outputdir,"/dom/pmeta_5pcs_vcf_5pc_diff_dens.txt"))
cat("Marker,ID,Beta,Pmeta,BetaAFR,BetaEUR,BetaAMR,BetaFIN,padditive \n")
for(i in 1:nrow(out)){
cat(paste0(out[i,"MarkerName"],",",out[i,"ID"],",",
round(out[i,"Effect"],digits=3),",",
format(out[i,"P.value"],digits=3, scientific=T),",",
round(out[i,"beta_AFR"],digits=3),",",
round(out[i,"beta_EUR"],digits=3),",",
round(out[i,"beta_AMR"],digits=3),",",
round(out[i,"beta_FIN"],digits=3),",",
format(out[i,"p_add"], digits=3, scientific=T),"\n"))
}
sink()


