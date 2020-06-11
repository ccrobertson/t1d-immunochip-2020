#supplementary_tables_generate_dens_new.R
#generating a number of supplementary Tables for the manuscript
#this uis an updated set of supplementary tables, comining some together as we have too many supplementary tables for nature genetics.
#1) All Univariable results, including Phase I, Phase II, and overall results
#2) Univariable hits + conditionally independent results to <5*10^-8 for Phase I (this includes the Phase II results for the univariable associations (missing for conditional)
#3) Univariable hits (<5*10^-8) and FDR (<0.01) hits for Phase I and II meta-analysis
#4) Dominent and recessive hits, including Phase I and AIC of additive and dominent models from european collections only

library(GenomicRanges)
library(stringr)
##############################################################################
#1) All	Univariable results, including Phase I,	Phase II, and overall results#
##############################################################################
extension="final_collection"
cohorts=c("EUR","FIN","AFR","AMR")
cohorts1=c("EUR","FIN","AFR","AMR","EAS")
cassiedir<-"/nv/vol185/MEGA/release4/IMPUTED_TOPMED/family_analysis/"
metdir<-paste0("/scratch/ji3x/",extension,"_gwas/META/META/")
outdir<-paste0("/scratch/ji3x/",extension,"_gwas/META/")
outputdir<-paste0("~/output/",extension,"_gwas/META/")

#define ichip regions:
ichip<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/rawdat/dens_and_ichip_regs.txt"),header=T,as.is=T,sep="\t")
ichip<-GRanges(seqnames=ichip$seqnames,
IRanges(ichip$start, end=ichip$end))

#and add the SNPs that were on the iChip to this set:
ichip1<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/rawdat/topmed_imputation_input_variants_b38_sorted_excludingMismatches.txt"),header=F, as.is=T)
ichip1<-GRanges(seqnames=paste0("chr",ichip1$V1),
IRanges(ichip1$V2, end=ichip1$V2))

k<-setdiff(ichip1,ichip)
allichip<-c(ichip, k)


load(file=paste0(metdir,"/resultsmeta_5pcs_vcf_5pc_diff_dens.RData"))
rownames(res)<-res$MarkerName
res$alleleA<-res$Allele1
res$alleleB<-res$Allele2
res$EffectphaseI<-res$Effect
res$SEphaseI<-res$StdErr
res$pphaseI<-res$P.value
res$rsid<-res$MarkerName 

res<-res[,c("rsid","chromosome","position","alleleA","alleleB","EffectphaseI","SEphaseI","pphaseI")]

#for each variant, get the required info:effect, se and af for each collection,
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
r$AIC<-14-(2*r$frequentist_add_ll)
r<-r[,c(cols,"AF","AIC")]
return(r)
}
combineallsnps<-function(cohort){
l<-lapply(c(1:22),combinethesnps,cohort=cohort)
l<-do.call("rbind",l)
return(l)
}


getinfo<-function(cohort){
init<-combineallsnps(cohort)
init<-init[init$rsid %in% res$rsid,]
init[,paste0("Effect_",cohort)]<-init$frequentist_add_beta_1.add.t1d.1
init[,paste0("SE_",cohort)]<-init$frequentist_add_se_1
init[,paste0("AF_",cohort)]<-init$AF
if(cohort=="EUR"){
init[,paste0("AIC_",cohort)]<-init$AIC
}
rownames(init)<-init$rsid
if(cohort!="EUR"){
keeps<-c(paste0("Effect_",cohort),paste0("SE_",cohort),paste0("AF_",cohort))
}
if(cohort=="EUR"){
keeps<-c(paste0("Effect_",cohort),paste0("SE_",cohort),paste0("AF_",cohort),paste0("AIC_",cohort))
}
init<-init[,c("rsid","chromosome","position","alleleA","alleleB",keeps)]
return(init)
}
indep<-lapply(cohorts,getinfo)

ind<-merge(res,indep[[1]],all.x=T, by=c("rsid","chromosome","position","alleleA","alleleB"))
ind<-merge(ind,indep[[2]],all.x=T, by=c("rsid","chromosome","position","alleleA","alleleB"))
ind<-merge(ind,indep[[3]],all.x=T, by=c("rsid","chromosome","position","alleleA","alleleB"))
ind<-merge(ind,indep[[4]],all.x=T, by=c("rsid","chromosome","position","alleleA","alleleB"))

#Now get the family results:
readinfam<-function(cohort, chromosome){
r<-read.table(paste0(cassiedir,"/",cohort,"/chr",chromosome,".filter_maf_gt_0.005_rsq_gt_0.8.relateds_with_pheno.tdt"), header=T)
return(r)
}
getallfam<-function(cohort){
g<-lapply(c(1:22),cohort=cohort, readinfam)
g<-do.call("rbind",g)
return(g)
}

getbetas<-function(cohort){
fam<-getallfam(cohort)
fam$A1<-as.character(fam$A1)
fam$A2<-as.character(fam$A2)
fam$allele1<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\3",fam$SNP)
fam$allele2<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\4",fam$SNP)
fam$chromosome<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\1",fam$SNP)
fam$chromosome<-gsub("chr","",fam$chromosome)
fam$chromosome<-as.numeric(fam$chromosome)
fam$position<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\2",fam$SNP)
fam$position<-as.numeric(fam$position)

#ensure the effect directions are all alligned to TOPMED reference panel:
w<-which(fam$A2==fam$allele2 & fam$A1==fam$allele1)
fam[w,"A1dummy"]<-fam[w,"A2"]
fam[w,"A2dummy"]<-fam[w,"A1"]
fam[w,"A1"]<-fam[w,"A1dummy"]
fam[w,"A2"]<-fam[w,"A2dummy"]
fam[,paste0("Effect_fam",cohort)]<-log(fam$OR)
fam[w,paste0("Effect_fam",cohort)]<-fam[w,paste0("Effect_fam",cohort)]*-1
fam[,paste0("SE_fam",cohort)]<-sqrt((1/fam$T)+(1/fam$U))
fam[,paste0("SE_fam",cohort)]<-ifelse(is.infinite(fam[,paste0("SE_fam",cohort)]),NA,fam[,paste0("SE_fam",cohort)])
men<-read.table(file=paste0(cassiedir,"/",cohort,"/snps_with_mendelInconsistencies_gt1pct.txt"), header=T, as.is=T)
fam<-fam[!fam$SNP %in% men$SNP,]
fam$alleleA<-fam$A2
fam$alleleB<-fam$A1
keeps<-c(paste0("Effect_fam",cohort),paste0("SE_fam",cohort))
fam$rsid<-fam$SNP
fam<-fam[,c("rsid","chromosome","position","alleleA","alleleB",keeps)]
return(fam)
}
fams<-lapply(cohorts1, getbetas)
names(fams)<-cohorts1


all<-merge(ind, fams[[1]],all.x=T, by=c("rsid","chromosome","position","alleleA","alleleB"))
all<-merge(all,	fams[[2]],all.x=T, by=c("rsid","chromosome","position","alleleA","alleleB"))
all<-merge(all,	fams[[3]],all.x=T, by=c("rsid","chromosome","position","alleleA","alleleB"))
all<-merge(all,	fams[[4]],all.x=T, by=c("rsid","chromosome","position","alleleA","alleleB"))
all<-merge(all,	fams[[5]],all.x=T, by=c("rsid","chromosome","position","alleleA","alleleB"))

#phase II results:
load(file=paste0(metdir,"/resultsmeta_fam_only_dens.RData"))
res$effect_phaseII<-res$Effect
res$SE_phaseII<-res$StdErr
res$rsid<-res$MarkerName
res$alleleA<-res$Allele1
res$alleleB<-res$Allele2
res$p_phaseII<-res$P.value
res<-res[c("rsid","chromosome","position","alleleA","alleleB","effect_phaseII","SE_phaseII","p_phaseII")]
all<-merge(all, res,all.x=T, by=c("rsid","chromosome","position","alleleA","alleleB"))


#finally the PhaseIandII results:
load(file=paste0(metdir,"/resultsmeta_ind_fam_tdt_5pc_diff_dens.RData"))
res$rsid<-res$MarkerName
res$EffectphaseIandII<-res$Effect
res$SEphaseIandII<-res$StdErr
res$pphaseIandII<-res$P.value
res$alleleA<-res$Allele1
res$alleleB<-res$Allele2
res<-res[,c("rsid","chromosome","position","alleleA","alleleB","EffectphaseIandII","SEphaseIandII","pphaseIandII")]
final<-merge(all,res,by=c("rsid","chromosome","position","alleleA","alleleB"),all.x=T)


#get RSIDs where they are present:
df<-final
o<-df[,c("chromosome","position")]
write.table(o, file=paste0("/scratch/ji3x/full_cohort_gwas/vars/signewfr2.txt"), sep="\t",
col.names=F, row.names=F, quote=F)

system(paste0("/home/ji3x/software/bcftools/bcftools view -T /scratch/ji3x/full_cohort_gwas",
"/vars/signewfr2.txt /scratch/ji3x/full_cohort_gwas/",
"/vars/Kaviar-160204-Public-hg38-trim.vcf.gz > /scratch/ji3x/full_cohort_gwas",
"/vars/signewfr2_ids.txt"))

vars<-read.table(file=paste0("/scratch/ji3x/full_cohort_gwas",
"/vars/signewfr2_ids.txt"),skip=40, header=T, comment.char="",as.is=T)
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

vars$MarkerName<-ifelse(vars$MarkerName1 %in% df$rsid,vars$MarkerName1,
ifelse(vars$MarkerName2 %in% df$rsid,vars$MarkerName2,
ifelse(vars$MarkerName3 %in% df$rsid,vars$MarkerName3,
ifelse(vars$MarkerName4 %in% df$rsid,vars$MarkerName4,
ifelse(vars$MarkerName5 %in% df$rsid,vars$MarkerName5,
ifelse(vars$MarkerName6 %in% df$rsid,vars$MarkerName6,
ifelse(vars$MarkerName7 %in% df$rsid,vars$MarkerName7,
ifelse(vars$MarkerName8 %in% df$rsid,vars$MarkerName8,
ifelse(vars$MarkerName9 %in% df$rsid,vars$MarkerName9,
ifelse(vars$MarkerName10 %in% df$rsid,vars$MarkerName10,
ifelse(vars$MarkerName11 %in% df$rsid,vars$MarkerName11,
ifelse(vars$MarkerName12 %in% df$rsid,vars$MarkerName12,
ifelse(vars$MarkerName13 %in% df$rsid,vars$MarkerName13,
ifelse(vars$MarkerName14 %in% df$rsid,vars$MarkerName14,
ifelse(vars$MarkerName15 %in% df$rsid,vars$MarkerName15,
ifelse(vars$MarkerName16 %in% df$rsid,vars$MarkerName16,
ifelse(vars$MarkerName17 %in% df$rsid,vars$MarkerName17,
ifelse(vars$MarkerName18 %in% df$rsid,vars$MarkerName18,
ifelse(vars$MarkerName19 %in% df$rsid,vars$MarkerName19,
ifelse(vars$MarkerName20 %in% df$rsid,vars$MarkerName20,
ifelse(vars$MarkerName21 %in% df$rsid,vars$MarkerName21,NA)))))))))))))))))))))

vars<-vars[!is.na(vars$MarkerName),]
vars<-vars[!duplicated(vars$MarkerName),]
vars$rsid<-vars$MarkerName
vars<-vars[,c("rsid","ID")]
df<-merge(df,vars,by="rsid",all.x=T)
df$rsid<-as.character(df$rsid)
df$ID<-ifelse(is.na(df$ID), df$rsid,df$ID)

final<-df
final$MarkerName<-final$rsid

ichip<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/rawdat/dens_and_ichip_regs.txt"),header=T,as.is=T,sep="\t")
ichip<-GRanges(seqnames=ichip$seqnames,
IRanges(ichip$start, end=ichip$end))

outr<-GRanges(seqnames=paste0("chr",final$chromosome),
IRanges(final$position,end=final$position),
MarkerName=final$MarkerName)

outr<-mergeByOverlaps(outr,ichip)
outr<-as.data.frame(outr)

final$ichipreg<-ifelse(final$MarkerName %in% outr$MarkerName,"yes","no")


fi<-final[,c("MarkerName","ichipreg","chromosome","position","ID","alleleA","alleleB",
"Effect_EUR","SE_EUR","AF_EUR",
"Effect_FIN","SE_FIN","AF_FIN",
"Effect_AFR","SE_AFR","AF_AFR",
"Effect_AMR","SE_AMR","AF_AMR",
"EffectphaseI","SEphaseI","pphaseI",
"Effect_famEUR","SE_famEUR",
"Effect_famFIN","SE_famFIN",
"Effect_famAFR","SE_famAFR",
"Effect_famAMR","SE_famAMR",
"Effect_famEAS","SE_famEAS",
"effect_phaseII","SE_phaseII",
"p_phaseII",
"EffectphaseIandII",
"SEphaseIandII",
"pphaseIandII")]
fi<-fi[order(fi$chromosome, fi$position),]
fi$AF_EUR<-ifelse(fi$MarkerName=="chr1:93145882:G:A",0.0001,fi$AF_EUR)
write.table(fi, file=paste0(outputdir,"/summaries/all_univariable_dens.txt"),
sep=",",col.names=T,row.names=F, quote=F)
save(final, fi, file=paste0(outputdir,"/summaries/robs.RData"))

#####################
#2) Univariable hits#
#####################
load(file=paste0(outputdir,"/summaries/robs.RData"))
hits<-read.table(file=paste0(outputdir,"pmeta_fam_ind_risk_regions_vcf_5pc_diff_dens.txt"),header=T, as.is=T,sep=",")
hits$EURAF<-ifelse(hits$Marker=="chr1:93145882:G:A",0.0001,hits$EURAF)

f<-fi[fi$MarkerName %in% hits$Marker,]
f<-f[order(f$chromosome, f$position),]
write.table(f, file=paste0(outputdir,"/summaries/genomewide_hits_univariable_dens.txt"),
sep=",",col.names=T,row.names=F, quote=F)




####################
#Hits from Phase I#
###################
hitsp1<-read.table(file=paste0(outputdir,"/pmeta_5pcs_vcf_5pc_diff_dens.txt"),header=T, as.is=T,sep=",")

f1<-fi[fi$MarkerName %in% hitsp1$Marker,]
f1<-f1[order(f1$chromosome, f1$position),]
f1$condrsid=""

#as well as the conditional hits:
cond<-read.table(file=paste0("~/output/",extension,"_gwas/META/conditional/init_detail_diff_dens_genomewide.txt"),header=T,as.is=T,sep="\t")
cond$MarkerName<-cond$rsid
f1<-f1[,colnames(f1) %in% colnames(cond)]
cond<-cond[,colnames(cond) %in% colnames(f1)]
cond<-cond[!cond$ID %in% f1$ID,]
f1<-rbind(f1,cond)
ichip<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/rawdat/dens_and_ichip_regs.txt"),header=T,as.is=T,sep="\t")
ichip<-GRanges(seqnames=ichip$seqnames,
IRanges(ichip$start, end=ichip$end))

outr<-GRanges(seqnames=paste0("chr",f1$chromosome),
IRanges(f1$position,end=f1$position),
MarkerName=f1$MarkerName)

outr<-mergeByOverlaps(outr,ichip)
outr<-as.data.frame(outr)

f1$ichipreg<-ifelse(f1$MarkerName %in% outr$MarkerName,"yes","no")
f1<-f1[!f1$ID %in% "rs78636954",]
f1$condno<-str_count(f1$condrsid, "\\+")+1
f1<-f1[order(f1$chromosome,f1$position,-f1$condno),]
f1$gene<-ifelse(f1$ID=="rs2269240","PGM1",ifelse(f1$ID=="rs2476601","PTPN22",
ifelse(f1$ID=="rs2229238","IL6R",ifelse(f1$ID=="rs78037977","FASLG",
ifelse(f1$ID=="rs6690988","CAMSAP2",ifelse(f1$ID=="rs3024493","IL10",
ifelse(f1$ID=="rs55730781","IL10",ifelse(f1$ID=="rs11120029","TATDN3",
ifelse(f1$ID=="rs10169963","TRIB2",ifelse(f1$ID=="rs2111485","IFIH1",
ifelse(f1$ID=="rs72871627","IFIH1",ifelse(f1$ID=="rs35667974","IFIH1",
ifelse(f1$ID=="rs6434435","STAT4",ifelse(f1$ID=="rs3087243","CTLA4",
ifelse(f1$ID=="rs112165453","CTLA4",ifelse(f1$ID=="rs376119912","CTLA4",
ifelse(f1$ID=="rs113341849","CCR5",ifelse(f1$ID=="rs113881148","TMEM175",
ifelse(f1$ID=="rs6848239","RBPJ",ifelse(f1$ID=="rs337637","KLF3",
ifelse(f1$ID=="rs70950892","IL2/IL21",ifelse(f1$ID=="rs2611215","CPE",
ifelse(f1$ID=="rs2287900","IL7R",ifelse(f1$ID=="rs1876142","PTGER4",NA))))))))))))))))))))))))
f1$gene<-ifelse(f1$ID=="rs10213692","ANKRD55",ifelse(f1$ID=="rs72928038","BACH2",
ifelse(f1$ID=="rs17754780","CENPW",ifelse(f1$ID=="rs12665429","TNFAIP3",
ifelse(f1$ID=="rs212408","TAGAP",ifelse(f1$ID=="rs73069533","SKAP2",
ifelse(f1$ID=="rs1018942","COBL",ifelse(f1$ID=="rs6944602","IKZF1",
ifelse(f1$ID=="rs1624088","CTSB",ifelse(f1$ID=="rs1574285","GLIS3",
ifelse(f1$ID=="rs10815223","PLGRKT",ifelse(f1$ID=="rs61839660","IL2RA",
ifelse(f1$ID=="rs35285258","IL2RA",ifelse(f1$ID=="rs56179589","IL2RA",
ifelse(f1$ID=="rs6602437","IL2RA",ifelse(f1$ID=="rs947474","IL2RA",
ifelse(f1$ID=="rs722988","NRP1",ifelse(f1$ID=="rs2018705","RNLS",
ifelse(f1$ID=="rs689","INS",ifelse(f1$ID=="rs7119275","INS",
ifelse(f1$ID=="rs10836367","SLC1A2",ifelse(f1$ID=="rs79538630","CD5/CD6",
ifelse(f1$ID=="rs607703","FLI1",ifelse(f1$ID=="rs917911","CD69",
ifelse(f1$ID=="rs7313065","ITGB7",ifelse(f1$ID=="rs34415530","IKZF4",
ifelse(f1$ID=="rs77508451","IKZF4",ifelse(f1$ID=="rs7310615","SH2B3",
ifelse(f1$ID=="rs9557217","GPR183",ifelse(f1$ID=="rs174213","ZFP36L1",
ifelse(f1$ID=="rs911263","ZFP36L1",ifelse(f1$ID=="rs4900384","AP000997.1",
ifelse(f1$ID=="rs56994090","MEG3",ifelse(f1$ID=="rs2289702","CTSH",
ifelse(f1$ID=="rs12927355","DEXI",ifelse(f1$ID=="rs66718203","DEXI",
ifelse(f1$ID=="rs4238595","UMOD",ifelse(f1$ID=="rs151174","IL27",f1$gene))))))))))))))))))))))))))))))))))))))
f1$gene<-ifelse(f1$ID=="rs72802365","BCAR1",ifelse(f1$ID=="rs56380902","GSDMB",
ifelse(f1$ID=="rs8070723","MAPT2",ifelse(f1$ID=="rs80262450","PTPN2",
ifelse(f1$ID=="rs62097857","PTPN2",ifelse(f1$ID=="rs12969657","CD226",
ifelse(f1$ID=="rs34536443","TYK2",ifelse(f1$ID=="rs12720356","TYK2",
ifelse(f1$ID=="rs314675","PRKD2",ifelse(f1$ID=="rs601338","FUT2",
ifelse(f1$ID=="rs6043409","SIRPG",ifelse(f1$ID=="rs6072275","LPIN3",
ifelse(f1$ID=="rs11203203","UBASH3A",ifelse(f1$ID=="rs58911644","ICOSLG",
ifelse(f1$ID=="rs41171","ASCC2/LIF",ifelse(f1$ID=="rs229531","C1QTNF6",f1$gene))))))))))))))))

f1<-f1[,c("MarkerName","gene","ichipreg","chromosome","position","ID",
"alleleA","alleleB","condrsid","Effect_EUR","SE_EUR","AF_EUR",
"Effect_FIN","SE_FIN","AF_FIN","Effect_AFR","SE_AFR","AF_AFR",
"Effect_AMR","SE_AMR","AF_AMR","EffectphaseI","SEphaseI","pphaseI")]
for (k in c("Effect_EUR","SE_EUR","AF_EUR","Effect_FIN","SE_FIN","AF_FIN",
"Effect_AFR","SE_AFR","AF_AFR","Effect_AMR","SE_AMR","AF_AMR","EffectphaseI","SEphaseI")){
f1[,k]<-round(f1[,k],digits=3)
}
f1$pphaseI<-format(f1$pphaseI,digits=3,SCIENTIFIC=T)
write.table(f1, file=paste0(outputdir,"/summaries/genomewide_hits_univariable_p1_dens_plus_cond.txt"),
sep=",",col.names=T,row.names=F, quote=F)



library(ggplot2)
g<-ggplot(f1, aes(EffectphaseI, effect_phaseII)) + geom_point() +
geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0)) + 
geom_abline(aes(intercept=0, slope=1)) + 
scale_x_continuous(name="Phase I log-odds ratio") + 
scale_y_continuous(name="Phase II log-odds ratio")
png(file=paste0(outputdir,"/phase1_vs_phase2_logodds_dens.png"),
height=20, width=20, units="cm", res=300)
g
dev.off()

#now removing any hits where the EUR associaiton statistics aren't present in both:
feur<-f1[!is.na(f1$Effect_EUR) & !is.na(f1$Effect_famEUR),]
g<-ggplot(feur, aes(EffectphaseI, effect_phaseII)) + geom_point() +
geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=0)) +
geom_abline(aes(intercept=0, slope=1)) +
scale_x_continuous(name="Phase I log-odds ratio") +
scale_y_continuous(name="Phase II log-odds ratio")
png(file=paste0(outputdir,"/phase1_vs_phase2_logodds_eur_dens.png"),
height=20, width=20, units="cm", res=300)
g
dev.off()



######################
#FDR associated 1% BY#
######################
hitsfdr<-read.table(file=paste0(outputdir,"/pmeta_fam_ind_risk_regions_fdr_by_0.01_dens.txt"),header=T,as.is=T,sep=",")

f2<-fi[fi$MarkerName %in% hitsfdr$Marker,]
f2<-f2[order(f2$chromosome, f2$position),]
write.table(f2, file=paste0(outputdir,"/summaries/genomewide_hits_univariable_fdr_dens.txt"),
sep=",",col.names=T,row.names=F, quote=F)

f2$abseffect<-abs(f2$EffectphaseIandII)
f2$or<-exp(f2$abseffect)
f3<-f2[order(f2$pphaseIandII),]
f3$ord<-c(1:nrow(f3))
f3$variant<-ifelse(f3$ID==".",f3$MarkerName,f3$ID)
f3$MarkerName<-as.factor(f3$MarkerName)
f3$MarkerName<-reorder(f3$MarkerName, f3$ord)
f3$maf_eur<-ifelse(f3$AF_EUR>0.5,1-f3$AF_EUR,f3$AF_EUR)
bonferroni<-0.05/nrow(fi)
w<-which(f3$pphaseIandII>bonferroni)
w<-w[1]
bon<-w-0.5
png(file=paste0(outputdir,"/effect_sizes_genomewide_fdr_dens.png"),height=20, width=30, units='cm', res=500)
ggplot(data=f3, aes(as.factor(rsid),or, colour=maf_eur))	+ geom_point() + theme(axis.text.x=element_text(angle=90)) +
scale_x_discrete(name="Variant") + 
scale_y_continuous(name="Absolute odds ratio") +
geom_vline(aes(xintercept=74.5), colour="red",linetype="dashed") +
geom_vline(aes(xintercept=bon), colour="red",linetype="dotted") +
dev.off()

f3$variant<-as.factor(f3$variant)
f3$variant<-reorder(f3$variant, f3$ord)
pdf(file=paste0(outputdir,"/effect_sizes_genomewide_fdr_dens.pdf"),width=12,height=7)
ggplot(data=f3, aes(as.factor(variant),or, colour=maf_eur))        + geom_point() + theme(axis.text.x=element_text(angle=90)) +
scale_x_discrete(name="Variant") +
scale_y_continuous(name="Absolute odds ratio") +
scale_colour_continuous(name="EUR MAF") +
geom_vline(aes(xintercept=78.5), colour="red",linetype="dashed") +
geom_vline(aes(xintercept=bon), colour="red",linetype="dotted") +
theme(axis.text.x=element_text(size=7))
dev.off()

png(file=paste0(outputdir,"/effect_sizes_genomewide_fdr_dens.png"),width=30,height=20, units='cm', res=500)
ggplot(data=f3, aes(as.factor(variant),or, colour=maf_eur))        + geom_point() + theme(axis.text.x=element_text(angle=90)) +
scale_x_discrete(name="Variant") +
scale_y_continuous(name="Absolute odds ratio") +
scale_colour_continuous(name="EUR MAF") +
geom_vline(aes(xintercept=78.5), colour="red",linetype="dashed") +
geom_vline(aes(xintercept=bon), colour="red",linetype="dotted") +
theme(axis.text.x=element_text(size=7))
dev.off()



#compare effect size and AF of variants associated to genomewide and FDR1%:
f2$genomewide<-(f2$MarkerName %in% f$MarkerName)
library(plyr)
library(dplyr)
f2$abseffect<-abs(f2$EffectphaseIandII)
f2$eurmaf<-ifelse(f2$AF_EUR>0.5,1-f2$AF_EUR,f2$AF_EUR)
sum<- f2 %>%
group_by(genomewide) %>%
summarise(meff=median(abseffect), lowereff=quantile(abseffect,prob=0.25),
uppereff=quantile(abseffect,prob=0.75), mafeur=median(eurmaf),
lowemaf=quantile(eurmaf, prob=0.25),uppermaf=quantile(eurmaf, prob=0.75))


###########################################################################################################################
#3) Dominent univariable results, including Phase I and AIC of additive and dominent models from european collections only#
###########################################################################################################################

combinethesnpsresdom<-function(cohort, chrom, resdomname, resdom){
r<-read.table(file=paste0(outdir,"/",cohort,"/",resdomname,"/chr_",chrom), header=T, as.is=T)
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
r[w,paste0("frequentist_",resdom,"_beta_1.",resdom,".t1d.1")]<-r[w,paste0("frequentist_",resdom,"_beta_1.",resdom,".t1d.1")]*-1
r[w,"abdummy"]<-r[w,"controls_AB"]
r[w,"aadummy"]<-r[w,"controls_AA"]
r[w,"bbdummy"]<-r[w,"controls_BB"]
r[w,"controls_BB"]<-r[w,"aadummy"]
r[w,"controls_AA"]<-r[w,"bbdummy"]
r$AF<-(r$controls_AB+(2*r$controls_BB))/((2*r$controls_AA) + 2*r$controls_AB + (2*r$controls_BB))

#remove the results of those with MAF<0.005:
message(paste0("Removing SNPs with MAF<0.005"))
r[r$all_maf<0.005,paste0("frequentist_",resdom,"_beta_1.",resdom,".t1d.1")]<-NA
r[r$all_maf<0.005,paste0("frequentist_",resdom,"_wald_pvalue_1")]<-NA
r[r$all_maf<0.005,paste0("frequentist_",resdom,"_se_1")]<-NA
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
r[w2,paste0("frequentist_",resdom,"_beta_1.",resdom,".t1d.1")]<-NA
r[w2,paste0("frequentist_",resdom,"_wald_pvalue_1")]<-NA
r[w2,paste0("frequentist_",resdom,"_se_1")]<-NA
r[w2,"all_info"]<-NA
#remove imputed  SNPs outside ichip regions:
r1<-GRanges(seqnames=paste0("chr",r$chromosome),
IRanges(r$position,end=r$position),
MarkerName=r$rsid)

ot<-mergeByOverlaps(r1,allichip)
ot<-as.data.frame(ot)
w<-which(!r$rsid %in% ot$MarkerName)
r[w,paste0("frequentist_",resdom,"_beta_1.",resdom,".t1d.1")]<-NA
r[w,paste0("frequentist_",resdom,"_wald_pvalue_1")]<-NA
r[w,paste0("frequentist_",resdom,"_se_1")]<-NA
r[w,"all_info"]<-NA
r$AIC<-14-(2*r[,paste0("frequentist_",resdom,"_ll")])
r<-r[,c(cols,"AF","AIC")]
return(r)
}

combineallsnpsresdom<-function(cohort,resdomname, resdom){
l<-lapply(c(1:22),combinethesnpsresdom,cohort=cohort,resdomname=resdomname, resdom=resdom)
l<-do.call("rbind",l)
return(l)
}

getitresdom<-function(cohort, resdomname, resdom){
init<-combineallsnpsresdom(cohort,resdomname=resdomname, resdom=resdom)
init<-init[init$rsid %in% res$rsid,]
init[,paste0("Effect_",cohort)]<-init[,paste0("frequentist_",resdom,"_beta_1.",resdom,".t1d.1")]
init[,paste0("SE_",cohort)]<-init[,paste0("frequentist_",resdom,"_se_1")]
init[,paste0("AF_",cohort)]<-init[,paste0("AF")]
if(cohort=="EUR"){
init[,paste0("AIC_",cohort)]<-init[,paste0("AIC")]
}
rownames(init)<-init$rsid
if(cohort=="EUR"){
keeps<-c(paste0("Effect_",cohort),paste0("SE_",cohort),paste0("AIC_",cohort),paste0("AF_",cohort))
}
if(cohort!="EUR"){
keeps<-c(paste0("Effect_",cohort),paste0("SE_",cohort),paste0("AF_",cohort))
}
init<-init[,c("rsid","chromosome","position","alleleA","alleleB",keeps)]
return(init)
}

load(file=paste0(metdir,"/resultsmeta_5pcs_vcf_5pc_diff_dens.RData"))
res$Effectmeta<-res$Effect
res$SEmeta<-res$StdErr
res$pmeta<-res$P.value
res$alleleA<-res$Allele1
res$alleleB<-res$Allele2
res$rsid<-res$MarkerName
res<-res[,c("rsid","chromosome","position",
"alleleA","alleleB","Effectmeta","SEmeta","pmeta")]

dom<-lapply(cohorts,getitresdom,resdom="dom", resdomname="dom")

allcomb<-merge(res,dom[[1]], all.x=T, by=c("rsid","chromosome","position","alleleA","alleleB"))
allcomb<-merge(allcomb,dom[[2]], all.x=T, by=c("rsid","chromosome","position","alleleA","alleleB"))
allcomb<-merge(allcomb,dom[[3]], all.x=T, by=c("rsid","chromosome","position","alleleA","alleleB"))
allcomb<-merge(allcomb,dom[[4]], all.x=T, by=c("rsid","chromosome","position","alleleA","alleleB"))
colnames(allcomb)[colnames(allcomb) %in% "AIC_EUR"]<-"AIC_EUR_DOM"

fin<-final[,c("MarkerName","position","chromosome","ID","alleleA","alleleB","AIC_EUR")]
fin$rsid<-fin$MarkerName
allcomb<-merge(allcomb, fin, by=c("rsid","position","chromosome","alleleA","alleleB"), all.x=T)
colnames(allcomb)[colnames(allcomb) %in% "AIC_EUR"]<-"AIC_EUR_ADD"

allcomb<-allcomb[,c("rsid","chromosome","position","ID","alleleA","alleleB",
"Effect_EUR","SE_EUR","AF_EUR","Effect_FIN","SE_FIN","AF_FIN","Effect_AFR","SE_AFR","AF_AFR",
"Effect_AMR","SE_AMR","AF_AMR","Effectmeta","SEmeta","pmeta","AIC_EUR_DOM", "AIC_EUR_ADD")]

allcomb<-allcomb[order(allcomb$chromosome, allcomb$position),]
write.table(allcomb, file=paste0(outputdir,"/summaries/all_univariable_dom_dens.txt"),
sep=",",col.names=T,row.names=F, quote=F)

##################
#4) Dominent hits#
##################
domhit<-read.table(file=paste0(outputdir,"/dom/pmeta_5pcs_vcf_5pc_diff_dens_fdr.txt"),header=T,
as.is=T, sep=",")
domhit<-domhit[!duplicated(domhit$Marker),]
allhitsdom<-allcomb[allcomb$rsid %in% domhit$Marker,]
allhitsdom<-allhitsdom[order(allhitsdom$chromosome, allhitsdom$position),]
allhitsdom<-allhitsdom[!is.na(allhitsdom$AIC_EUR_DOM),]
allhitsdom<-allhitsdom[allhitsdom$AIC_EUR_DOM<allhitsdom$AIC_EUR_ADD,]
w<-which(allhitsdom$AF_EUR>0.5)
allhitsdom[w,"AF_EUR"]<-1-allhitsdom[w,"AF_EUR"]
allhitsdom[w,"AF_FIN"]<-1-allhitsdom[w,"AF_FIN"]
allhitsdom[w,"AF_AFR"]<-1-allhitsdom[w,"AF_AFR"]
allhitsdom[w,"AF_AMR"]<-1-allhitsdom[w,"AF_AMR"]
allhitsdom[w,"Effect_EUR"]<-allhitsdom[w,"Effect_EUR"]*-1
allhitsdom[w,"Effect_FIN"]<-allhitsdom[w,"Effect_FIN"]*-1
allhitsdom[w,"Effect_AFR"]<-allhitsdom[w,"Effect_AFR"]*-1
allhitsdom[w,"Effect_AMR"]<-allhitsdom[w,"Effect_AMR"]*-1
allhitsdom[w,"Effectmeta"]<-allhitsdom[w,"Effectmeta"]*-1
allhitsdom[w,"adum"]<-allhitsdom[w,"alleleB"]
allhitsdom[w,"bdum"]<-allhitsdom[w,"alleleA"]
allhitsdom[w,"alleleA"]<-allhitsdom[w,"adum"]
allhitsdom[w,"alleleB"]<-allhitsdom[w,"bdum"]
allhitsdom[,!(colnames(allhitsdom) %in% c("adum","bdum"))]
allhitsdom$type<-"Dominant"
write.table(allhitsdom, file=paste0(outputdir,"/summaries/all_hits_dom_dens_fdr.txt"),
sep=",",col.names=T,row.names=F, quote=F)


#############################################################################################################################
#5) Recessive univariable results, including Phase I and AIC of additive and recessive models from european collections only#
#############################################################################################################################

load(file=paste0(metdir,"/res/resultsmeta_5pcs_vcf_5pc_diff_dens.RData"))
res$Effectmeta<-res$Effect
res$SEmeta<-res$StdErr
res$pmeta<-res$P.value
res$alleleA<-res$Allele1
res$alleleB<-res$Allele2
res$rsid<-res$MarkerName
res<-res[,c("rsid","chromosome","position",
"alleleA","alleleB","Effectmeta","SEmeta","pmeta")]

rec<-lapply(cohorts,getitresdom, resdom="rec",resdomname="res")

allcomb1<-merge(res,rec[[1]], all.x=T, by=c("rsid","chromosome","position","alleleA","alleleB"))
allcomb1<-merge(allcomb1,rec[[2]], all.x=T, by=c("rsid","chromosome","position","alleleA","alleleB"))
allcomb1<-merge(allcomb1,rec[[3]], all.x=T, by=c("rsid","chromosome","position","alleleA","alleleB"))
allcomb1<-merge(allcomb1,rec[[4]], all.x=T, by=c("rsid","chromosome","position","alleleA","alleleB"))
colnames(allcomb1)[colnames(allcomb1) %in% "AIC_EUR"]<-"AIC_EUR_REC"


allcomb1<-merge(allcomb1, fin, by=c("rsid","position","chromosome","alleleA","alleleB"), all.x=T)
colnames(allcomb1)[colnames(allcomb1) %in% "AIC_EUR"]<-"AIC_EUR_ADD"

allcomb1<-allcomb1[,c("rsid","chromosome","position","ID","alleleA","alleleB",
"Effect_EUR","SE_EUR","AF_EUR","Effect_FIN","SE_FIN","AF_FIN","Effect_AFR","SE_AFR","AF_AFR",
"Effect_AMR","SE_AMR","AF_AMR","Effectmeta","SEmeta","pmeta","AIC_EUR_REC", "AIC_EUR_ADD")]
allcomb1<-allcomb1[order(allcomb1$chromosome, allcomb1$position),]
write.table(allcomb1, file=paste0(outputdir,"/summaries/all_univariable_rec_dens.txt"),
sep=",",col.names=T,row.names=F, quote=F)

###################
#6) Recessive hits#
###################
rechit<-read.table(file=paste0(outputdir,"/res/pmeta_5pcs_vcf_5pc_diff_dens_fdr.txt"),header=T,
as.is=T, sep=",")
rechit<-rechit[!duplicated(rechit$Marker),]

allhitsrec<-allcomb1[allcomb1$rsid %in% rechit$Marker,]
allhitsrec<-allhitsrec[order(allhitsrec$chromosome, allhitsrec$position),]
allhitsrec<-allhitsrec[allhitsrec$AIC_EUR_REC<allhitsrec$AIC_EUR_ADD,]
w<-which(allhitsrec$AF_EUR>0.5)
allhitsrec[w,"AF_EUR"]<-1-allhitsrec[w,"AF_EUR"]
allhitsrec[w,"AF_FIN"]<-1-allhitsrec[w,"AF_FIN"]
allhitsrec[w,"AF_AFR"]<-1-allhitsrec[w,"AF_AFR"]
allhitsrec[w,"AF_AMR"]<-1-allhitsrec[w,"AF_AMR"]
allhitsrec[w,"Effect_EUR"]<-allhitsrec[w,"Effect_EUR"]*-1
allhitsrec[w,"Effect_FIN"]<-allhitsrec[w,"Effect_FIN"]*-1
allhitsrec[w,"Effect_AFR"]<-allhitsrec[w,"Effect_AFR"]*-1
allhitsrec[w,"Effect_AMR"]<-allhitsrec[w,"Effect_AMR"]*-1
allhitsrec[w,"Effectmeta"]<-allhitsrec[w,"Effectmeta"]*-1
allhitsrec[w,"adum"]<-allhitsrec[w,"alleleB"]
allhitsrec[w,"bdum"]<-allhitsrec[w,"alleleA"]
allhitsrec[w,"alleleA"]<-allhitsrec[w,"adum"]
allhitsrec[w,"alleleB"]<-allhitsrec[w,"bdum"]
allhitsrec[,!(colnames(allhitsrec) %in% c("adum","bdum"))]
allhitsrec$type<-"Recessive"
write.table(allhitsrec, file=paste0(outputdir,"/summaries/all_hits_rec_dens_fdr.txt"),
sep=",",col.names=T,row.names=F, quote=F)

colnames(allhitsrec)[22]<-"AIC_EUR_DOM"

allhitsdom<-rbind(allhitsdom, allhitsrec)
allhitsdom<-allhitsdom[order(allhitsdom$chromosome, allhitsdom$position),]
write.table(allhitsdom, file=paste0(outputdir,"/summaries/all_hits_dom_dens_fdr.txt"),
sep=",",col.names=T,row.names=F, quote=F)

