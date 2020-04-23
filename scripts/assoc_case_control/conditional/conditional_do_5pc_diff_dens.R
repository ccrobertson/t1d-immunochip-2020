#conditional_do_5pc_diff_dens.R

#generate the SNPTEST scripts which will carry out the stepwise conditional regressions
library(parallel)
extension="final_collection"
vcfpath="/nv/vol185/MEGA/release4/IMPUTED_TOPMED/results_filtered/"
cohorts=c("AMR","EUR","AFR","FIN")
path=paste0("/scratch/ji3x/",extension,"_gwas/vcf/conditional/")
outputdir<-paste0("~/output/",extension,"_gwas/META/")
metdir<-paste0("/scratch/ji3x/",extension,"_gwas/META/META/")
dobbydir<-paste0("/scratch/ji3x/",extension,"_gwas/vcf/")
resdir<-paste0("/scratch/ji3x/",extension,"_gwas/META/")

signals<-read.table(file=paste0(outputdir,"/pmeta_5pcs_vcf_5pc_diff_dens.txt"),as.is=T, sep=",",header=T)

signals$chromosome<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\1",signals$Marker)
signals$chromosome<-gsub("chr","",signals$chromosome)
signals$chromosome<-as.numeric(signals$chromosome)
signals$position<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\2",signals$Marker)
signals$position<-as.numeric(signals$position)

load(file=paste0(metdir,"/resultsmeta_5pcs_vcf_5pc_diff_dens.RData"))
#generate the snptest scripts:
ext="final_collection"
snptestscripts<-function(snp,cohort, adj=NULL){
s<-signals[signals$Marker==snp,]
snp1<-gsub(":","_",snp)
snp2<-signals[signals$Marker==snp,"Marker"]
if (is.null(adj)){
sink(file=paste0("~/programs/",ext,"/META/conditional/snptest_scripts/",cohort,"/",snp1,"_do_5pc_diff_dens.sh"))
cat(paste0("#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH -A rich_immunochip_impute
#SBATCH -p standard
#SBATCH -n 1

/home/ji3x/software/snptest_v2.5.4-beta3_linux_x86_64_static/snptest_v2.5.4-beta3 -data ",
path,cohort,"/",snp1,"out.vcf.gz",
" ",dobbydir,"/",cohort,"/samps_snptest -filetype vcf -genotype_field GP -pheno t1d ",
"-frequentist 1 -method newml -cov_names PC1 PC2 PC3 PC4 PC5 -condition_on ",snp2,
" add -o ",path,"/results/",cohort,"/",snp1,"_5pc_diff_dens\n"))
sink()
system(paste0("sbatch ~/programs/",ext,"/META/conditional/snptest_scripts/",cohort,"/",snp1,"_do_5pc_diff_dens.sh"))
}
if (!is.null(adj)){
j<-strsplit(adj, split="+", fixed=T)
n<-length(j[[1]])
adjusted<-paste0(j[[1]]," add")
adjusted<-paste(adjusted, collapse=" ")
sink(file=paste0("~/programs/",ext,"/META/conditional/snptest_scripts/",cohort,"/",snp1,"_do_5pc_diff_dens_",n,".sh"))
cat(paste0("#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH -A rich_immunochip_impute
#SBATCH -p standard
#SBATCH -n 1

/home/ji3x/software/snptest_v2.5.4-beta3_linux_x86_64_static/snptest_v2.5.4-beta3 -data ",
path,cohort,"/",snp1,"out.vcf.gz",
" ",dobbydir,"/",cohort,"/samps_snptest -filetype vcf -genotype_field GP -pheno t1d ",
"-frequentist 1 -method newml -cov_names PC1 PC2 PC3 PC4 PC5 -condition_on ",adjusted,
" -o ",path,"/results/",cohort,"/",snp1,"_5pc_diff_dens_",n,"\n"))
sink()
system(paste0("sbatch /home/ji3x/programs/",ext,"/META/conditional/snptest_scripts/",cohort,"/",snp1,"_do_5pc_diff_dens_",n,".sh"))
}
}

#second signals?
lapply(signals$Marker,snptestscripts, cohort="EUR")
lapply(signals$Marker,snptestscripts, cohort="AMR")
lapply(signals$Marker,snptestscripts, cohort="AFR")
lapply(signals$Marker,snptestscripts, cohort="FIN")

stayInLoop=TRUE
while(stayInLoop) {
Sys.sleep(30)
a<-system("squeue -u ji3x | wc -l", intern=T)
stayInLoop <- a != "1"
message(paste0("Still ",a,"jobs..."))
}

readmetain<-function(snp,cohort, adj=NULL){
snp1<-gsub(":","_",snp)
if(is.null(adj)){
g<-tryCatch({r<-read.table(file=paste0(path,"results/",cohort,"/",snp1,"_5pc_diff_dens"), header=T, as.is=T)
r$conditional=snp
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
#filter out those that failed QC:
r<-r[r$rsid %in% res$MarkerName,]
#and those already adjusted for:
r<-r[!r$rsid %in% snp,]
#and re do the QC quickly:
r$diff<-abs(r$cases_info-r$controls_info)
w<-which(r$diff>0.05 | r$all_info<0.8 | r$controls_info<0.8 | r$cases_info<0.8)
r[w,"frequentist_add_beta_1.add.t1d.1"]<-NA
r[w,"frequentist_add_se_1"]<-NA
r[w,"frequentist_add_wald_pvalue_1"]<-NA
g<-r[,colnames(r) %in% c("rsid", "chromosome", "position", "alleleA", "alleleB",
"all_info", "all_maf", "cases_maf", "controls_maf", "frequentist_add_wald_pvalue_1",
"frequentist_add_beta_1.add.t1d.1","frequentist_add_se_1","AF")]
g$weight<-1/(g$frequentist_add_se_1^2)
g$betawt<-g$weight*g$frequentist_add_beta_1.add.t1d.1
g$or<-exp(g$frequentist_add_beta_1.add.t1d.1)
colnames(g)[!colnames(g) %in% c("rsid","position","chromosome","alleleA","alleleB")]<-
paste0(colnames(g)[!colnames(g) %in% c("rsid","position","chromosome","alleleA","alleleB")],"_",cohort)
return(g)
}, warning = function(w) {
    g<-data.frame(rsid=NA, chromosome=NA,  position=NA, alleleA=NA, alleleB=NA)
g[,paste0("all_maf_",cohort)]=NA
g[,paste0("all_info_",cohort)]=NA
g[,paste0("cases_maf_",cohort)]=NA
g[,paste0("controls_maf_",cohort)]=NA
g[,paste0("frequentist_add_beta_1.add.t1d.1_",cohort)]=NA
g[,paste0("frequentist_add_se_1_",cohort)]=NA
g[,paste0("frequentist_add_wald_pvalue_1_",cohort)]=NA
g[,paste0("AF_",cohort)]=NA
g[,paste0("weight_",cohort)]=NA
g[,paste0("betawt_",cohort)]=NA
g[,paste0("or_",cohort)]=NA
return(g)
}, error = function(e) {
    g<-data.frame(rsid=NA, chromosome=NA,  position=NA, alleleA=NA, alleleB=NA)
g[,paste0("all_maf_",cohort)]=NA
g[,paste0("all_info_",cohort)]=NA
g[,paste0("cases_maf_",cohort)]=NA
g[,paste0("controls_maf_",cohort)]=NA
g[,paste0("frequentist_add_beta_1.add.t1d.1_",cohort)]=NA
g[,paste0("frequentist_add_se_1_",cohort)]=NA
g[,paste0("frequentist_add_wald_pvalue_1_",cohort)]=NA
g[,paste0("AF_",cohort)]=NA
g[,paste0("weight_",cohort)]=NA
g[,paste0("betawt_",cohort)]=NA
g[,paste0("or_",cohort)]=NA
return(g)
})
}
if(!is.null(adj)){
j<-strsplit(adj, split="+", fixed=T)
n<-length(j[[1]])

g<-tryCatch({r<-read.table(file=paste0(path,"results/",cohort,"/",snp1,"_5pc_diff_dens_",n),header=T, as.is=T)
r$conditional=adj
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
#filter	out those that failed QC:
r<-r[r$rsid %in% res$MarkerName,]
#and re do the QC quickly:
r$diff<-abs(r$cases_info-r$controls_info)
w<-which(r$diff>0.05 | r$all_info<0.8 | r$cases_info<0.8 | r$controls_info<0.8)
r[w,"frequentist_add_beta_1.add.t1d.1"]<-NA
r[w,"frequentist_add_se_1"]<-NA
r[w,"frequentist_add_wald_pvalue_1"]<-NA

#and remove those already adjusted for:
a<-strsplit(adj, split="+",fixed=TRUE)
r<-r[!r$rsid %in% c(snp,a[[1]]),]
g<-r[,colnames(r) %in% c("rsid", "chromosome", "position", "alleleA", "alleleB",
"all_info", "all_maf", "cases_maf", "controls_maf", "frequentist_add_wald_pvalue_1",
"frequentist_add_beta_1.add.t1d.1","frequentist_add_se_1","AF")]
g$weight<-1/(g$frequentist_add_se_1^2)
g$betawt<-g$weight*g$frequentist_add_beta_1.add.t1d.1
g$or<-exp(g$frequentist_add_beta_1.add.t1d.1)
colnames(g)[!colnames(g) %in% c("rsid","position","chromosome","alleleA","alleleB")]<-
paste0(colnames(g)[!colnames(g) %in% c("rsid","position","chromosome","alleleA","alleleB")],"_",cohort)
return(g)
}, warning = function(w) {
g<-data.frame(rsid=NA, chromosome=NA,  position=NA, alleleA=NA, alleleB=NA)
g[,paste0("all_maf_",cohort)]=NA
g[,paste0("all_info_",cohort)]=NA
g[,paste0("cases_maf_",cohort)]=NA
g[,paste0("controls_maf_",cohort)]=NA
g[,paste0("frequentist_add_beta_1.add.t1d.1_",cohort)]=NA
g[,paste0("frequentist_add_se_1_",cohort)]=NA
g[,paste0("frequentist_add_wald_pvalue_1_",cohort)]=NA
g[,paste0("AF_",cohort)]=NA
g[,paste0("weight_",cohort)]=NA
g[,paste0("betawt_",cohort)]=NA
g[,paste0("or_",cohort)]=NA
return(g)
}, error = function(e) {
g<-data.frame(rsid=NA, chromosome=NA,  position=NA, alleleA=NA, alleleB=NA)
g[,paste0("all_maf_",cohort)]=NA
g[,paste0("all_info_",cohort)]=NA 
g[,paste0("cases_maf_",cohort)]=NA 
g[,paste0("controls_maf_",cohort)]=NA 
g[,paste0("frequentist_add_beta_1.add.t1d.1_",cohort)]=NA 
g[,paste0("frequentist_add_se_1_",cohort)]=NA 
g[,paste0("frequentist_add_wald_pvalue_1_",cohort)]=NA 
g[,paste0("AF_",cohort)]=NA 
g[,paste0("weight_",cohort)]=NA 
g[,paste0("betawt_",cohort)]=NA 
g[,paste0("or_",cohort)]=NA 
return(g)
})
}
return(g)
}

metathem<-function(snp, adj=NULL){
if(is.null(adj)){
eur<-readmetain("EUR", snp=snp)
amr<-readmetain("AMR",snp=snp)
afr<-readmetain("AFR",snp=snp)
fin<-readmetain("FIN",snp=snp)
}
if(!is.null(adj)){
j<-strsplit(adj, split="+", fixed=T)
n<-length(j[[1]])
eur<-readmetain("EUR",snp=snp,adj=adj)
amr<-readmetain("AMR",snp=snp,adj=adj)
afr<-readmetain("AFR",snp=snp,adj=adj)
fin<-readmetain("FIN",snp=snp,adj=adj)
}
meta<-merge(eur,amr,by=c("rsid","chromosome","position","alleleA","alleleB"), all.x=T)
meta<-merge(meta,afr,by=c("rsid","chromosome","position","alleleA","alleleB"), all=T)
meta<-merge(meta,fin,by=c("rsid","chromosome","position","alleleA","alleleB"), all=T)

meta$weightsum<-rowSums(meta[,c("weight_EUR","weight_AMR", "weight_AFR","weight_FIN")], na.rm=TRUE)
meta$seall<-sqrt(1/meta$weightsum)
meta$sumbeta<-rowSums(meta[,c("betawt_EUR","betawt_AMR", "betawt_AFR","betawt_FIN")], na.rm=TRUE)
meta$beta<-meta$sumbeta/meta$weightsum
meta$z<-meta$beta/meta$seall
meta$pmeta<-2*pnorm(-abs(meta$z))
m<-meta[meta$pmeta==min(meta$pmeta,na.rm=T),]
m<-m[!is.na(m$rsid),]
m$or<-exp(m$beta)
m$p<-m$pmeta
if (!is.null(adj)){
m$conditional<-paste0(adj)
}
if (is.null(adj)){
m$conditional<-gsub("_",":",snp)
}
m$index<-snp
for(co in c("EUR","FIN","AFR","AMR")){
m[,paste0("Effect_",co)]<-m[,paste0("frequentist_add_beta_1.add.t1d.1_",co)]
m[,paste0("SE_",co)]<-m[,paste0("frequentist_add_se_1_",co)]
m[,paste0("AF_",co)]<-m[,paste0("AF_",co)]
}
m$EffectphaseI<-m$beta
m$SEphaseI<-m$seall
m$pphaseI<-m$p
m<-m[,c("rsid","chromosome","position","alleleA","alleleB","or_EUR",
"or_AMR", "or_AFR","or_FIN","or","p","conditional","index", "all_maf_EUR",
"all_maf_AMR","all_maf_AFR","all_maf_FIN","Effect_EUR","SE_EUR","AF_EUR",
"Effect_FIN","SE_FIN","AF_FIN","Effect_AFR","SE_AFR","AF_AFR",
"Effect_AMR","SE_AMR","AF_AMR","EffectphaseI","SEphaseI","pphaseI")]
message(paste0("DONE ",snp))
return(m)
}

s<-mclapply(gsub(":","_",signals$Marker), metathem, mc.cores=4)
s<-do.call("rbind", s)
save(s, file=paste0(path,"/results//META/signals_2_5pc_diff_dens.RData"))
load(file=paste0(path,"/results/META/signals_2_5pc_diff_dens.RData"))


while(min(s$p,na.rm=T)<5*10^-6){
#have a look at third signals:
stillsig<-s[s$p<5*10^-6,"index"]
adjfor<-paste0(s[s$p<5*10^-6,"conditional"],"+",s[s$p<5*10^-6,"rsid"])
j1<-strsplit(adjfor[1], split="+", fixed=T)
n1<-length(j1[[1]])+1
mapply(snptestscripts, cohort=rep("EUR",length(adjfor)), adj=adjfor, snp=stillsig)
mapply(snptestscripts, cohort=rep("AMR",length(adjfor)), adj=adjfor, snp=stillsig)
mapply(snptestscripts, cohort=rep("AFR",length(adjfor)), adj=adjfor, snp=stillsig)
mapply(snptestscripts, cohort=rep("FIN",length(adjfor)), adj=adjfor, snp=stillsig)

stayInLoop=TRUE
while(stayInLoop) {
Sys.sleep(30)
a<-system("squeue -u ji3x | wc -l", intern=T)
stayInLoop <- a != "1"
message(paste0("Still ",a,"jobs..."))
}


s<-mcmapply(metathem,snp=gsub(":","_",signals[gsub(":","_",signals$Marker) %in% s[s$p<5*10^-6,"index"],"Marker"]), adj=adjfor, SIMPLIFY=FALSE, mc.cores=4)
s<-do.call("rbind", s)
save(s, file=paste0(path,"/results/META/signals_",n1,"_5pc_diff_dens.RData"))
message(paste0("Done ",n1))
}


#get info for the unadjusted analysis:
signals<-read.table(file=paste0(outputdir,"/pmeta_5pcs_vcf_5pc_diff_dens.txt"),as.is=T, sep=",",header=T)
signals$chromosome<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\1",signals$Marker)
signals$chromosome<-gsub("chr","",signals$chromosome)
signals$chromosome<-as.numeric(signals$chromosome)
signals$position<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\2",signals$Marker)
signals$position<-as.numeric(signals$position)

getres<-function(cohort){
getall<-function(chrom){
r1<-read.table(file=paste0(resdir,"/",cohort,"/chr_",chrom), header=T, as.is=T)
return(r1)
}
r<-lapply(c(1:22),getall)
r<-do.call("rbind",r)
r<-r[r$rsid %in% signals$Marker,]
r[,paste0("or_",cohort)]<-exp(r$frequentist_add_beta_1.add.t1d.1)
r[,paste0("all_maf_",cohort)]<-r$all_maf
r[,paste0("Effect_",cohort)]<-r$frequentist_add_beta_1.add.t1d.1
r[,paste0("SE_",cohort)]<-r$frequentist_add_se_1
r[,paste0("AF_",cohort)]<-(r$controls_AB+(2*r$controls_BB))/((2*r$controls_AA) + 2*r$controls_AB + (2*r$controls_BB))
r<-r[,c("rsid",paste0("or_",cohort),paste0("Effect_",cohort),paste0("SE_",cohort),paste0("AF_",cohort),paste0("all_maf_",cohort))]
return(r)
}
a<-lapply(cohorts,getres)
keep<-merge(a[[1]],a[[2]],by="rsid",all=T)
keep<-merge(keep,a[[3]],by="rsid",all=T) 
keep<-merge(keep,a[[4]],by="rsid",all=T) 

res$Marker<-res$MarkerName
res$rsid<-res$MarkerName
keep<-merge(keep,res,by="rsid",all.x=T)
keep$conditional<-""
signals$p<-signals$Pmeta
signals$or<-exp(signals$Beta)
signals$pphaseI<-signals$p
signals$EffectphaseI<-signals$Beta
signals$SEphaseI<-signals$SE
leads<-signals[,c("Marker","ID","p","or","EffectphaseI","SEphaseI","pphaseI")]
keep$Marker<-keep$MarkerName
keep<-merge(keep,leads,by="Marker")
keep$index<-keep$rsid
keep$alleleA<-keep$Allele1
keep$alleleB<-keep$Allele2

keep<-keep[,c("rsid", "chromosome",
"position", "alleleA", "alleleB", "or_EUR",
"or_AMR", "or_AFR", "or_FIN", "or", "p", "conditional", 
"index","all_maf_EUR","all_maf_AMR",
"all_maf_AFR","all_maf_FIN","Effect_EUR","SE_EUR","AF_EUR",
"Effect_FIN","SE_FIN","AF_FIN","Effect_AFR","SE_AFR","AF_AFR",
"Effect_AMR","SE_AMR","AF_AMR","EffectphaseI","SEphaseI","pphaseI")]


#load all conditional results:
getc<-function(num){
g<-get(load(file=paste0(path,"/results/META/signals_",num,"_5pc_diff_dens.RData")))
return(g)
}
nums<-system(paste0("ls ",path,"/results/META/ | grep dens.RData | wc -l"),intern=T)
nums<-as.numeric(nums)
all<-lapply(c(2:nums),getc)
all<-do.call("rbind",all)

all<-rbind(keep,all)


all$chromosome<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\1",all$rsid)
all$chromosome<-gsub("chr","",all$chromosome)
all$chromosome<-as.numeric(all$chromosome)
all<-all[order(all$chromosome, all$position),]
#remove one there due to LD and boundaries:
allg<-all[all$p<5*10^-8,]
allg<-allg[!(allg$rsid=="chr7:50406053:G:A" & allg$conditional==""),]
allg<-allg[!(allg$rsid=="chr7:50970505:G:A" & allg$conditional=="chr7:50406053:G:A"),]
allg<-allg[!(allg$rsid=="chr14:68286876:C:T" & allg$conditional==""),]
allg<-allg[!(allg$rsid=="chr14:68804010:T:C" & allg$conditional=="chr14:68286876:C:T"),]

write.table(allg, file=paste0(outputdir,"/conditional/init_5pc_diff_dens_genomwide.txt"), sep="\t",
col.names=T, row.names=F, quote=F)

alls<-all[all$p<5*10^-6,]
alls<-alls[!(alls$rsid=="chr7:50406053:G:A" & alls$conditional==""),]
alls<-alls[!(alls$rsid=="chr7:50970505:G:A" & alls$conditional=="chr7:50406053:G:A"),]
alls<-alls[!(alls$rsid=="chr7:51022213:G:A" & alls$conditional=="chr7:50406053:G:A+chr7:50970505:G:A"),]
alls<-alls[!(alls$rsid=="chr14:68286876:C:T" & alls$conditional==""),]
alls<-alls[!(alls$rsid=="chr14:68804010:T:C" & alls$conditional=="chr14:68286876:C:T"),]

write.table(alls, file=paste0(outputdir,"/conditional/init_5pc_diff_dens_minus6.txt"), sep="\t",
col.names=T, row.names=F, quote=F)

#and now a nicer table:
#get rsids:
o<-all[,c("chromosome","position")]
write.table(o, file=paste0("/scratch/ji3x/full_cohort_gwas/vars/signewmlr2.txt"), sep="\t",
col.names=F, row.names=F, quote=F)

system(paste0("/home/ji3x/software/bcftools/bcftools view -R /scratch/ji3x/full_cohort_gwas",
"/vars/signewmlr2.txt /scratch/ji3x/full_cohort_gwas",
"/vars/Kaviar-160204-Public-hg38-trim.vcf.gz > /scratch/ji3x/full_cohort_gwas",
"/vars//signewmlr2_ids.txt"))

vars<-read.table(file=paste0("/scratch/ji3x/full_cohort_gwas",
"/vars/signewmlr2_ids.txt"),skip=40, header=T, comment.char="")
vars$rsid<-paste0("chr",vars$X.CHROM,":",vars$POS,":",vars$REF,":",vars$ALT)
vars<-vars[,c("rsid","ID")]


#do this for the 'conditional' column too:
getcond<-function(i){
f<-all[i,]
cond<-f$conditional
if(cond==""){
out=""
}
if(cond!=""){
c1<-strsplit(cond, split="+",fixed=TRUE)
c1<-c1[[1]]
o<-data.frame(chromosome=sub("(^.*)[:](.*)[:](.*)[:](.*)","\\1",c1),
position=sub("(^.*)[:](.*)[:](.*)[:](.*)","\\2",c1))
o$chromosome<-gsub("chr","",o$chromosome)
write.table(o, file=paste0("/scratch/ji3x/full_cohort_gwas",
"/vars/sigcond.txt"), sep="\t",
col.names=F, row.names=F, quote=F)
system(paste0("/home/ji3x/software/bcftools/bcftools view -R /scratch/ji3x/full_cohort_gwas",
"/vars/sigcond.txt /scratch/ji3x/full_cohort_gwas",
"/vars/Kaviar-160204-Public-hg38-trim.vcf.gz > /scratch/ji3x/full_cohort_gwas",
"/vars/sigcond_ids.txt"))

vars<-read.table(file=paste0("/scratch/ji3x/full_cohort_gwas",
"/vars/sigcond_ids.txt"),as.is=T,skip=40, header=T, comment.char="")
vars$MarkerName<-paste0("chr",vars$X.CHROM,":",vars$POS,":",vars$REF,":",vars$ALT)
vars<-vars[,c("MarkerName","ID")]
vars$ID<-ifelse(vars$ID==".",vars$MarkerName,vars$ID)
out=paste0(vars$ID,collapse="+")
}
return(out)
}
conds<-lapply(c(1:nrow(all)),getcond)
conds<-do.call("rbind", conds)
all$condrsid<-conds
all<-all[all$p<5*10^-6,]
all$condrsid<-ifelse(all$condrsid=="rs61839660+rs56179589+chr10:6062411:TTTTTTTTTTT:TTTTTTTTT,TTTTTTTTTTTT+rs35285258","rs61839660+rs56179589+rs35285258",
ifelse(all$condrsid=="rs61839660+rs56179589+chr10:6062411:TTTTTTTTTTT:TTTTTTTTT,TTTTTTTTTTTT+rs35285258+rs6602437","rs61839660+rs56179589+rs35285258+rs6602437",
ifelse(all$condrsid=="rs61839660+rs56179589+chr10:6062411:TTTTTTTTTTT:TTTTTTTTT,TTTTTTTTTTTT+rs35285258+rs6602437+chr10:6348487:CG:CA,TG+rs947474",
"rs61839660+rs56179589+rs35285258+rs6602437+rs947474",
ifelse(all$condrsid==paste0("rs71459332+rs77508451+chr12:56050824:TCTTCAATAATAAATTAACCTTAGCTTACTGTAACTTTT:TCTTCAATAATAAATTAACCTTAGTTTACTGTAACTTTT+",
"chr12:56050847:GC:GT+rs34415530+chr12:56050848:CTTA:TTTA+chr12:56050848:CTTACTGTAACTTTTTTTTTCTTTC:TTTACTGTAACTTTTTTTTTCTTTC"),
"rs71459332+rs77508451+rs34415530",
ifelse(all$condrsid==paste0("rs77508451+chr12:56050824:TCTTCAATAATAAATTAACCTTAGCTTACTGTAACTTTT:TCTTCAATAATAAATTAACCTTAGTTTACTGTAACTTTT+",
"chr12:56050847:GC:GT+rs34415530+chr12:56050848:CTTA:TTTA+chr12:56050848:CTTACTGTAACTTTTTTTTTCTTTC:TTTACTGTAACTTTTTTTTTCTTTC"),
"rs77508451+rs34415530",
ifelse(all$condrsid==paste0("chr12:56050824:TCTTCAATAATAAATTAACCTTAGCTTACTGTAACTTTT:TCTTCAATAATAAATTAACCTTAGTTTACTGTAACTTTT+",
"chr12:56050847:GC:GT+rs34415530+chr12:56050848:CTTA:TTTA+chr12:56050848:CTTACTGTAACTTTTTTTTTCTTTC:TTTACTGTAACTTTTTTTTTCTTTC"),
"rs34415530",ifelse(all$condrsid=="rs56380902+chr17:39910119:TG:CA","rs56380902",
ifelse(all$condrsid=="rs112165453+rs3087243+rs376119912+chr2:203888383:CT:C+rs796311366","rs112165453+rs3087243+rs376119912+rs796311366",
ifelse(all$condrsid=="chr20:1635559:CA:CG+rs6043409","rs6043409",all$condrsid)))))))))

all<-all[!(all$rsid=="chr7:50406053:G:A" & all$conditional==""),]
all<-all[!(all$rsid=="chr7:50970505:G:A" & all$conditional=="chr7:50406053:G:A"),]
all<-all[!(all$rsid=="chr7:51022213:G:A" & all$conditional=="chr7:50960703:T:C+chr7:50406053:G:A"),]
all<-all[!(all$rsid=="chr14:68286876:C:T" & all$conditional==""),]
all<-all[!(all$rsid=="chr14:68804010:T:C" & all$conditional=="chr14:68286876:C:T"),]
all<-merge(all,vars,by="rsid",all.x=TRUE)
all$index<-gsub("_",":",all$index)
all$length<-nchar(all$conditional)
all$posindex<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\2",all$index)
all$posindex<-as.numeric(all$posindex)
all<-all[order(all$chromosome, all$posindex, all$length),]
all$ID<-as.character(all$ID)
all$ID<-ifelse(is.na(all$ID),".",all$ID)
all$ID<-ifelse(all$rsid=="chr2:203888383:CT:C","rs376119912",all$ID)
all$ID<-ifelse(all$rsid=="chr16:28494339:G:C","rs151234",all$ID)
all$ID<-ifelse(all$rsid=="chr22:30026118:C:G","rs41171",all$ID)


allg<-all[all$p<5*10^-8,]
alls<-all[all$p<5*10^-6,]

addline<-function(i,df){
cat(paste0(df[i,"rsid"],";",
df[i,"ID"],";",
df[i,"condrsid"],";",
round(as(df[i,"or_EUR"],"numeric"),digits=3),";",
round(as(df[i,"or_AFR"],"numeric"),digits=3),";",
round(as(df[i,"or_AMR"],"numeric"),digits=3),";",
round(as(df[i,"or_FIN"],"numeric"),digits=3),";",
round(as(df[i,"or"],"numeric"),digits=3),";",	
format(df[i,"p"],scientific=TRUE, digits=3),"\n"))
}


sink(file=paste0(outputdir,"/conditional/init_tidy_5pc_diff_dens_genomewide.txt"))
cat("MarkerName;rsid;conditional on; OREUR;ORAFR;ORAMR;ORFIN;ORMETA;pMETA\n")
invisible(lapply(c(1:nrow(allg)),addline,df=allg))
sink()

sink(file=paste0(outputdir,"/conditional/init_tidy_5pc_diff_dens_minus6.txt"))
cat("MarkerName;rsid;conditional on; OREUR;ORAFR;ORAMR;ORFIN;ORMETA;pMETA\n")
invisible(lapply(c(1:nrow(alls)),addline,df=alls))
sink()



#another export to use in the supplementary table generation:
out<-allg[,c("rsid","ID","chromosome","position","alleleA", "alleleB",
"condrsid","Effect_EUR", "SE_EUR","AF_EUR","Effect_FIN","SE_FIN",
"AF_FIN","Effect_AFR","SE_AFR","AF_AFR","Effect_AMR","SE_AMR","AF_AMR",
"EffectphaseI","SEphaseI","pphaseI")]

write.table(out, file=paste0(outputdir,"/conditional/init_detail_diff_dens_genomewide.txt"),col.names=T, row.names=F,quote=F, sep="\t")
