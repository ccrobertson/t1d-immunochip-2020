#fdr_check_dens.R
#seeing how many associations at FDR<0.1 (BY):
library(GenomicRanges)
library(stringr)
extension="final_collection"
metdir<-paste0("/scratch/ji3x/",extension,"_gwas/META/META/")
outdir<-paste0("/scratch/ji3x/",extension,"_gwas/META/")
outputdir<-paste0("~/output/",extension,"_gwas/META/")

load(file=paste0(metdir,"/resultsmeta_ind_fam_tdt_5pc_diff_dens.RData"))

chr<-paste0("chr", res$chromosome)
start<-res$position
snp.name<-res$MarkerName

#remove MHC:
res<-res[!(res$chromosome==6 & res$position>25000000 & res$position<35000000),]
res$z<-res$Effect/res$StdErr
res$abz<-abs(res$z)

res$pfdr<-p.adjust(res$P.value,method="BY")


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

while(n$pfdr<0.01){
n<-getmin(dropped)
dropped<-dropnext(n,dropped)
both<-rbind(both, n)
}


both<-both[order(both$chromosome, both$position),]
both$or<-exp(both$Effect)


#any that weren't significant in the all individuals analysis?
resoriginal<-get(load(file=paste0(metdir,"/resultsmeta_5pcs_vcf_5pc_diff_dens.RData")))
resoriginal$porig<-resoriginal$P.value
resoriginal$efforig<-resoriginal$Effect
resoriginal<-resoriginal[,c("MarkerName","porig","efforig")]

sigregions<-merge(both, resoriginal, by="MarkerName", all.x=T)
sigregions<-sigregions[order(sigregions$chromosome, sigregions$position),]

#create regions file for BCF tools to extract SNPs needed:
o<-sigregions[,c("chromosome","position")]
write.table(o, file=paste0("/scratch/ji3x/full_cohort_gwas/vars/signewfr2.txt"), sep="\t",
col.names=F, row.names=F, quote=F)

system(paste0("~/software/bcftools/bcftools view -R /scratch/ji3x/full_cohort_gwas/vars/signewfr2.txt",
" /scratch/ji3x/full_cohort_gwas/vars/Kaviar-160204-Public-hg38-trim.vcf.gz > /scratch/ji3x/full_cohort_gwas/vars/signewfr2_ids.txt"))

vars<-read.table(file=paste0("/scratch/ji3x/full_cohort_gwas/vars/signewfr2_ids.txt"),skip=40, header=T, comment.char="",as.is=T)
vars$refcoms<-str_count(vars$REF,",")
vars$altcoms<-str_count(vars$ALT,",")
vars<-vars[vars$ID!=".",]
see<-strsplit(vars$ALT,split=",")
l<-lengths(see)
alts<-str_split_fixed(vars$ALT, ",", max(l)+1)
colnames(alts)<-paste0("alt",c(1:ncol(alts)))
vars<-cbind(vars,alts)
for(i in c(1:max(l))){
vars[,paste0("MarkerName",i)]<-paste0("chr",vars$X.CHROM,":",vars$POS,":",vars$REF,":",vars[,paste0("alt",i)])
}

vars$MarkerName<-ifelse(vars$MarkerName1 %in% sigregions$MarkerName,vars$MarkerName1,
ifelse(vars$MarkerName2 %in% sigregions$MarkerName,vars$MarkerName2,
ifelse(vars$MarkerName3 %in% sigregions$MarkerName,vars$MarkerName3,NA)))


vars<-vars[!is.na(vars$MarkerName),]
vars<-vars[!duplicated(vars$MarkerName),]
vars<-vars[,c("MarkerName","ID")]
sigregions<-merge(sigregions,vars,by="MarkerName",all.x=T)
sigregions$MakerName<-as.character(sigregions$MarkerName)
sigregions$ID<-ifelse(is.na(sigregions$ID), sigregions$MarkerName,sigregions$ID)


sigregions<-sigregions[sigregions$pfdr<0.01,]


#are these SNPs in iChip regions?:
#the ones in iChip regions will be brought forward for conditional regression analyses + fine mapping.
ichip<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/rawdat/dens_and_ichip_regs.txt"),header=T,as.is=T,sep="\t")
ichip<-GRanges(seqnames=ichip$seqnames,
IRanges(ichip$start, end=ichip$end))

outr<-GRanges(seqnames=paste0("chr",sigregions$chromosome),
IRanges(sigregions$position,end=sigregions$position),
MarkerName=sigregions$MarkerName)

outr<-mergeByOverlaps(outr,ichip)
outr<-as.data.frame(outr)

sigregions$ichipreg<-ifelse(sigregions$MarkerName %in% outr$MarkerName,"yes","no")


#and get european AF of alternative allele:
getaf<-function(snp){
s<-sigregions[sigregions$MarkerName==snp,]
s<-s[order(s$ID),]
s<-s[nrow(s),]
chrom<-s$chromosome
r<-read.table(file=paste0(outdir,"/EUR/chr_",chrom), header=T, as.is=T)
cols<-colnames(r)
r<-r[r$rsid %in% s$MarkerName,]
if(nrow(r)==0){
s$AF=NA
}
if(nrow(r)>0){
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
s$AF<-r$AF
}
return(s)
}
sigregions<-lapply(sigregions$MarkerName,getaf)
sigregions<-do.call("rbind",sigregions)


#remove those likely due to LD:
sigregions<-sigregions[!sigregions$MarkerName %in% c("chr1:113298204:C:T",
"chr3:45821775:T:G","chr12:110181544:G:C",
"chr12:111066229:T:C","chr12:112153882:G:A",
"chr12:112731156:C:A","chr17:46751565:G:A",
"chr22:29511412:C:T"),]

#exporting the regions to a word table:
sink(file=paste0(outputdir,"/pmeta_fam_ind_risk_regions_fdr_by_0.01_dens.txt"))
cat("Marker,ID,ichip,EURAF,Pcasecontrol,Beta,Pmeta\n")
for(i in 1:nrow(sigregions)){
cat(paste0(sigregions[i,"MarkerName"],",",sigregions[i,"ID"],",",
sigregions[i,"ichipreg"],",",round(sigregions[i,"AF"],digits=3),",",
format(sigregions[i,"porig"],digits=3, scientific=T),",",
sigregions[i,"Effect"],",",
format(sigregions[i,"P.value"],digits=3,scientific=T),"\n"))
}
sink()




