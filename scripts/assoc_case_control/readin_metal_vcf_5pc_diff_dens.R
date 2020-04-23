#readin_metal_vcf_5pc_diff_dens.R

#reads in results from initial meta run:

library(snpStats)
#library(annotSnpStats)
library(qqman)
library(GenomicRanges)
#change the name "extension" below in future analyses to generate same output in different location.
extension="final_collection"
metdir<-paste0("/scratch/ji3x/",extension,"_gwas/META/META/")
dobbypath="/scratch/ji3x/final_collection_gwas/vcf/"
outputdir<-"~/output/final_collection_gwas/META/"

extension="final_collection"
cohorts=c("AFR","EUR","AMR", "FIN")
res<-read.table(file=paste0(metdir,"/METAANALYSIS_5pc_diff_dens.TBL"), header=T, as.is=T)

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
samps<-read.table(file=paste0(dobbypath,"/",name,"/samples.txt"), header=F, as.is=T)
all<-read.table(file="/nv/vol185/MEGA/mega_genotyped/release4/mega_b37_clean_pheno_genotypedOnly.txt", header=F, as.is=T)
r<-all[all$V2 %in% samps$V1,]
r<-r[r$V7!=0,]
t<-table(r[,7])
r<-data.frame(n0=t[1],n1=t[2])
return(r)
}
reads<-lapply(cohorts, readitin)
reads<-do.call("rbind",reads)
n0<-sum(reads[,1])
n1<-sum(reads[,2])


lam1000<-1+(lam-1)*((1/n0) + (1/n1))/(2/1000)


save(res,file=paste0(metdir,"/resultsmeta_5pcs_vcf_5pc_diff_dens.RData"))

load(file=paste0(metdir,"/resultsmeta_5pcs_vcf_5pc_diff_dens.RData"))

dataframenona<-res[!is.na(res$P.value) & res$P.value!=0,]
dataframenona$logp<-log10(dataframenona$P.value)*-1
yval<-max(dataframenona$logp)
yval<-yval*0.8
pdf(file=paste0(outputdir,"/manhattan_qq_initial_5pcs_5pc_diff_dens.pdf"),
    width=7, height=7)
par(mfrow=c(2,1))
manhattan(dataframenona, p="P.value", bp="position", chr="chromosome",cex=0.5)
qq(dataframenona$P.value,cex=0.5)
text(x=1,y=yval, labels=expression(paste(lambda,"=")))
text(x=1.3, y=yval, labels=paste0(round(lam,2)))
dev.off()

tiff(file=paste0(outputdir,"/manhattan_qq_initial_5pcs_5pc_diff_dens.tiff"),
    width=7, height=7, units="in", res=400)
par(mfrow=c(2,1))
manhattan(dataframenona, p="P.value", bp="position", chr="chromosome",cex=0.5)
qq(dataframenona$P.value,cex=0.5)
text(x=1,y=yval, labels=expression(paste(lambda,"=")))
text(x=1.3, y=yval, labels=paste0(round(lam,2)))
dev.off()

png(file=paste0(outputdir,"/manhattan_qq_initial_5pcs_5pc_diff.png"),
    width=7, height=7, units="in", res=400)
par(mfrow=c(2,1))
manhattan(dataframenona, p="P.value", bp="position", chr="chromosome",cex=0.5)
qq(dataframenona$P.value,cex=0.5)
text(x=1,y=yval, labels=expression(paste(lambda,"=")))
text(x=1.3, y=yval, labels=paste0(round(lam,2)))
dev.off()

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
co<-read.table(file=paste0(metdir,"metal_results_",cohort,"_5pc_diff_dens.tbl"),
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
"/vars/signewmlr2_ids.txt"),skip=40, header=T, comment.char="")
vars$MarkerName<-paste0("chr",vars$X.CHROM,":",vars$POS,":",vars$REF,":",vars$ALT)
vars<-vars[,c("MarkerName","ID")]

out<-merge(out,vars,by="MarkerName", all.x=T)
out<-out[order(out$chromosome,out$position),]

out$ID<-as.character(out$ID)
out$ID<-ifelse(is.na(out$ID),".",out$ID)


out$ID<-as.character(out$ID)
out$ID<-ifelse(is.na(out$ID),".",out$ID)

#which of these are in ichip regions?:
ichip<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/rawdat/dens_and_ichip_regs.txt"),header=T,as.is=T,sep="\t")
ichip<-GRanges(seqnames=ichip$seqnames,
IRanges(ichip$start, end=ichip$end))

outr<-GRanges(seqnames=paste0("chr",out$chromosome),
IRanges(out$position,end=out$position),
MarkerName=out$MarkerName)

outr<-mergeByOverlaps(outr,ichip)
outr<-as.data.frame(outr)

out$ichipreg<-ifelse(out$MarkerName %in% outr$MarkerName,"yes","no")
#removing those variants there due to LD alone:
out<-out[!out$MarkerName %in% c("chr1:113298204:C:T",
"chr3:45821775:T:G","chr12:112153882:G:A","chr12:112731156:C:A",
"chr17:46751565:G:A"),]

out<-out[out$P.value<5*10^-8,]
out$ID<-ifelse(out$MarkerName=="chr22:30026118:C:G","rs41171",out$ID)
#exporting the regions to a word table:
sink(file=paste0(outputdir,"/pmeta_5pcs_vcf_5pc_diff_dens.txt"))
cat("Marker,ID,ichipreg,Beta,SE,Pmeta,BetaAFR,BetaEUR,BetaAMR,BetaFIN \n")
for(i in 1:nrow(out)){
cat(paste0(out[i,"MarkerName"],",",out[i,"ID"],",",
out[i,"ichipreg"],",",
round(out[i,"Effect"],digits=3),",",
round(out[i,"StdErr"],digits=3),",",
format(out[i,"P.value"],digits=3, scientific=T),",",
round(out[i,"beta_AFR"],digits=3),",",
round(out[i,"beta_EUR"],digits=3),",",
round(out[i,"beta_AMR"],digits=3),",",
round(out[i,"beta_FIN"],digits=3),"\n"))
}
sink()


