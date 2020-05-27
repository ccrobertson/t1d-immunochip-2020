#guessfm_sum_5pc_diff.R
#summarising the guessfm results in the europeans (newqc_1):
library(GUESSFM)
library(stringr)

extension="final_collection"
orig<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/pmeta_fam_ind_risk_regions_vcf_5pc_diff.txt"),header=T,sep=",",as.is=T)
orig<-orig[orig$ichip=="yes",]
orig<-orig[!(orig$Marker %in% c("chr12:112153882:G:A","chr12:112730563:A:G","chr7:50406053:G:A",
"chr1:113302202:A:C","chr17:46751565:G:A","chr14:68286876:C:T","chr12:111066229:T:C",
"chr17:40623690:G:A")),]

#first export all creidble SNPs for David and Tony
getall<-function(snp){
if (file.exists(paste0("/home/ji3x/output/",extension,"_gwas/guessfm_5pc_diff/",gsub(":",".",snp),"/credible_snps.txt"))){
eur<-read.table(file=paste0("~/output/",extension,"_gwas/guessfm_5pc_diff/",gsub(":",".",snp),"/credible_snps.txt"),
header=T, as.is=T, sep="\t")
eur$region<-paste0("chr",eur$chromosome[1],":",min(eur$position),"-",max(eur$position))
eur<-eur[,c("region","MarkerName","ID","tag","Allele1","Allele2","Effect","StdErr","P.value","Direction","chromosome","position",
"pp","ppsum")]
}
if (!file.exists(paste0("/home/ji3x/output/",extension,"_gwas/guessfm_5pc_diff/",gsub(":",".",snp),"/credible_snps.txt"))){
eur<-NULL
}
return(eur)
}
d<-lapply(orig$Marker,getall)
d<-do.call("rbind",d)

write.table(d, file=paste0("~/output/",extension,"_gwas/guessfm_5pc_diff/all_creds.txt"),sep="\t", col.names=T, row.names=F, quote=F)


df<-d[,-3]
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


vars$MarkerName<-ifelse(vars$MarkerName1 %in% df$MarkerName,vars$MarkerName1,
ifelse(vars$MarkerName2 %in% df$MarkerName,vars$MarkerName2,
ifelse(vars$MarkerName3 %in% df$MarkerName,vars$MarkerName3,
ifelse(vars$MarkerName4 %in% df$MarkerName,vars$MarkerName4,
ifelse(vars$MarkerName5 %in% df$MarkerName,vars$MarkerName5,NA)))))

vars<-vars[!is.na(vars$MarkerName),]
vars<-vars[!duplicated(vars$MarkerName),]
vars$rsid<-vars$MarkerName
vars<-vars[,c("MarkerName","ID")]
df<-merge(df,vars,by="MarkerName",all.x=T)
df$MarkerName<-as.character(df$MarkerName)
df$ID<-ifelse(is.na(df$ID), df$MarkerName,df$ID)

d<-df
d<-d[order(d$chromosome, d$position),]

write.table(d, file=paste0("~/output/",extension,"_gwas/guessfm_5pc_diff/all_creds_idcorr.txt"),sep="\t", col.names=T, row.names=F, quote=F)


d<-read.table(file=paste0("~/output/",extension,"_gwas/guessfm_5pc_diff/all_creds_idcorr.txt"),sep="\t",header=T,as.is=T)
#and gr38 for tony (he needs both 38 and 37:
#only those with ppsum>0.8
frame1<-d[d$ppsum>0.8 & !is.na(d$ppsum),]
frame1$chromosome<-paste0("chr",frame1$chromosome)
frame1<-frame1[,c("chromosome","position","ID")]
frame1$positionend<-frame1$position
frame1<-frame1[,c("chromosome","positionend","position","ID")]
write.table(frame1, file=paste0("~/output/",extension,"_gwas/guessfm_5pc_diff/creds_38_pp_0.8.txt"),
col.names=F, row.names=F, quote=F, sep="\t")

#now all of them:
frame<-d[,c("chromosome","position", "tag","ID","ppsum")]
frame$chromosome<-paste0("chr",frame$chromosome)
frame$positionend<-frame$position
frame<-frame[,c("chromosome","position", "positionend", "tag","ID","ppsum")]
write.table(frame, file=paste0("~/output/",extension,"_gwas/guessfm_5pc_diff/all_creds_38_tag.txt"),
col.names=F, row.names=F, quote=F, sep="\t")

frame$tag1<-sub("(^.*)[.](.*)[.](.*)[.](.*)","\\2",frame$tag)
frame<-frame[,c("chromosome","position", "positionend", "tag1")]
write.table(frame, file=paste0("~/output/",extension,"_gwas/guessfm_5pc_diff/all_creds_38_tag_bed.txt"),
col.names=F, row.names=F, quote=F, sep="\t")

#and by tag group:
tagsum<-function(tag){
f<-frame[frame$tag==tag,]
write.table(f, file=paste0("~/output/",extension,"_gwas/guessfm_5pc_diff/creds/38/creds_",tag,"_38.txt"),
col.names=F, row.names=F, quote=F, sep="\t")
}
lapply(unique(frame$tag),tagsum)





#export in gr37 for Tony:
library(GenomicRanges)
library(rtracklayer)
chr <- paste0("chr", d$chromosome)
start <- d$position
snp <- d$MarkerName
d$snp.name <- d$MarkerName
d$pos <- c(1:nrow(d))
chain.file.path <- paste0("/scratch/ji3x/liftover_allign/hg38ToHg19.over.chain")
input <- GRanges(seqname = Rle(chr), ranges = IRanges(start = start, 
     end = start), snp.name = snp)
c <- import.chain(chain.file.path)
alligned <- unlist(liftOver(input, c))
names(mcols(alligned)) <- "snp.name"
updated <- data.frame(snp.name = alligned@elementMetadata@listData$snp.name)
updated[, paste0("position37")] <- alligned@ranges@start
frame <- merge(d, updated, by = "snp.name")
frame <- frame[order(frame$pos), ]
frame$pos <- NULL
frame$chromosome<-paste0("chr",frame$chromosome)
#only those with ppsum>0.8
frame1<-frame[frame$ppsum>0.8,]
frame1<-frame1[,c("chromosome","position37","ID")]
frame1$positionend<-frame1$position37
frame1<-frame1[,c("chromosome","positionend","position37","ID")]
write.table(frame1, file=paste0("~/output/",extension,"_gwas/guessfm_5pc_diff/creds_37_pp_0.8.txt"), 
col.names=F, row.names=F, quote=F, sep="\t")

#now all of them:
frame<-frame[,c("chromosome","position37", "tag","ID","ppsum")]
frame$positionend<-frame$position37
frame<-frame[,c("chromosome","position37", "positionend", "tag","ID","ppsum")]
write.table(frame, file=paste0("~/output/",extension,"_gwas/guessfm_5pc_diff/all_creds_37_tag.txt"),
col.names=F, row.names=F, quote=F, sep="\t")

frame$tag1<-sub("(^.*)[.](.*)[.](.*)[.](.*)","\\2",frame$tag)
frame<-frame[,c("chromosome","position37", "positionend", "tag1")]
write.table(frame, file=paste0("~/output/",extension,"_gwas/guessfm_5pc_diff/all_creds_37_tag_bed.txt"),
col.names=F, row.names=F, quote=F, sep="\t")


#and by tag group:
tagsum<-function(tag){
f<-frame[frame$tag==tag,]
write.table(f, file=paste0("~/output/",extension,"_gwas/guessfm_5pc_diff/creds/37/creds_",tag,".txt"),
col.names=F, row.names=F, quote=F, sep="\t")
}
lapply(unique(frame$tag),tagsum)


getres<-function(snp){
if (file.exists(paste0("/home/ji3x/output/",extension,"_gwas/guessfm_5pc_diff/",gsub(":",".",snp),"/credible_snps.txt"))){
eurdd<-read.snpmod(paste0("/scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/",gsub(":",".",snp)))
eur<-read.table(file=paste0("~/output/",extension,"_gwas/guessfm_5pc_diff/",gsub(":",".",snp),"/credible_snps.txt"),
header=T, as.is=T, sep="\t")

postnumeur<-pp.nsnp(eurdd,plot=FALSE,expected=3, overdispersion = 1.00000001)
nprobe<-as.numeric(names(postnumeur$trait[postnumeur$trait==max(postnumeur$trait)]))

out<-data.frame(snp=snp,nsig=nprobe, ncreds=nrow(eur))
}
if (!file.exists(paste0("/home/ji3x/output/",extension,"_gwas/guessfm_5pc_diff/",gsub(":",".",snp),"/credible_snps.txt"))){
out<-data.frame(snp=snp,nsig=NA, ncreds=NA)
}
return(out)
}

l<-lapply(orig$Marker,getres)
l<-do.call("rbind", l)

write.table(l, file=paste0("~/output/",extension,"_gwas/guessfm_5pc_diff/summary_nsigs.txt"),sep="\t", col.names=T, row.names=F, quote=F)


#now per group:
groupswithless<-function(snp){
if (file.exists(paste0("/home/ji3x/output/",extension,"_gwas/guessfm_5pc_diff/",gsub(":",".",snp),"/credible_snps.txt"))){
eur<-read.table(file=paste0("~/output/",extension,"_gwas/guessfm_5pc_diff/",gsub(":",".",snp),"/credible_snps.txt"),
header=T, as.is=T, sep="\t")
}

j<-length(unique(eur$tag))
tags<-unique(eur$tag)
getsize<-function(tag,marker){
g<-eur[eur$tag==tag,]
n<-nrow(g)
out<-data.frame(marker=marker,tag=tag, n=n, ppsum=g$ppsum[1])
return(out)
}
l<-lapply(tags, getsize, marker=snp)
l<-do.call("rbind",l)
}
if (!file.exists(paste0("/home/ji3x/output/",extension,"_gwas/guessfm_5pc_diff/",gsub(":",".",snp),"/credible_snps.txt"))){
l<-data.frame(marker=snp, tag=NA, n=NA, ppsum=NA)
}
return(l)
}

groups<-lapply(orig$Marker,groupswithless)
groups<-do.call("rbind", groups)


promisinggroups<-groups[groups$n<20 & !is.na(groups$n) & groups$ppsum>0.8 & !is.na(groups$ppsum),]
promisinggroups$ppsum<-round(promisinggroups$ppsum,digits=3)
write.table(promisinggroups, file=paste0("~/output/",extension,"_gwas/guessfm_5pc_diff/promising_ngroups.txt"),sep="\t", col.names=T, row.names=F, quote=F)
getcreds<-function(snp, tag){
j<-read.table(file=paste0("~/output/",extension,"_gwas/guessfm_5pc_diff/",gsub(":",".",snp),"/credible_snps.txt"),
header=T, as.is=T, sep="\t")
j<-j[j$tag==tag,]
return(j)
}
promisingcreds<-mapply(getcreds, snp=promisinggroups$marker, tag=promisinggroups$tag, SIMPLIFY=FALSE)
promisingcreds<-do.call("rbind", promisingcreds)

write.table(promisingcreds, file=paste0("~/output/",extension,"_gwas/guessfm_5pc_diff/promising_summary.txt"),sep="\t", col.names=T, row.names=F, quote=F)




#highlight suspcious groups as those without any actually genotyped SNPs on.
genod<-read.table(file="/nv/vol185/MEGA/release4/processed_for_imputation/topmed_alignment/mega_b37-updated.bim",header=F,as.is=T)
colnames(genod)<-c("chromosome","snp.name","cM","position","A1","A2")
#lift to build 38:
library(GenomicRanges)
library(rtracklayer)
chr <- paste0("chr", genod$chromosome)
start <- genod$position
snp <- genod$snp.name
genod$pos <- c(1:nrow(genod))
chain.file.path <- paste0("/scratch/ji3x/liftover_allign/hg19ToHg38.over.chain")
input <- GRanges(seqname = Rle(chr), ranges = IRanges(start = start,
     end = start), snp.name = snp)
c <- import.chain(chain.file.path)
alligned <- unlist(liftOver(input, c))
names(mcols(alligned)) <- "snp.name"
updated <- data.frame(snp.name = alligned@elementMetadata@listData$snp.name)
updated[, paste0("position38")] <- alligned@ranges@start
frame <- merge(genod, updated, by = "snp.name")
frame <- frame[order(frame$pos), ]
frame$pos <- NULL

d<-d[d$ppsum>0.5,]
checkifany<-function(group){
g<-d[d$tag==group,]
chr<-g$chromosome[1]
f<-frame[frame$chromosome==chr,]
n<-nrow(g[g$position %in% f$position38,])
if(n==0){
out<-data.frame(group=group, anygenod="No")
}
if(n>0){
out<-data.frame(group=group, anygenod="Yes")
}
return(out)
}
see<-lapply(unique(d$tag),checkifany)
see<-do.call("rbind",see)



#bar graph summarising the number of variuants per group per region - prior to inclusion of excluded variants due to imputation issues (i.e. not figure 1 in the manuscript)
d<-read.table(file=paste0("~/output/",extension,"_gwas/guessfm_5pc_diff/all_creds_idcorr.txt"),sep="\t", header=T)
regs<-d[d$ppsum>0.5,]
regs$ID<-as.character(regs$ID)
library(plyr)
library(dplyr)

regs$tagit<-ifelse(gsub(":",".",regs$MarkerName)==regs$tag,1,0)
regs$tag1<-ifelse(regs$tagit==1,regs$ID,NA)

doall<-function(tag){
p<-regs[regs$tag==tag,]
p<-p[order(p$tag1),]
p$tag1<-p$tag1[1]
return(p)
}
regs<-lapply(unique(regs$tag),doall)
regs<-do.call("rbind",regs)
regs$genename<-ifelse(regs$region=="chr1:113761186-113834946","PTPN22",
ifelse(regs$region=="chr1:172699548-173286670","FASLG",
ifelse(regs$region=="chr1:192513038-192575969","RGS1",
ifelse(regs$region=="chr1:206629095-206812410","IL10",
ifelse(regs$region=="chr2:100024524-100254400","AFF3",
ifelse(regs$region=="chr2:162157475-162469787","IFIH1",
ifelse(regs$region=="chr2:191051012-191137812","STAT4",
ifelse(regs$region=="chr2:203708579-203922837","CTLA4",
ifelse(regs$region=="chr3:45888730-46440558","CCR5",
ifelse(regs$region=="chr4:26083858-26127088","RBPJ",
ifelse(regs$region=="chr4:122110339-122620186","IL2_IL21",
ifelse(regs$region=="chr5:35820791-35924646","IL7R",
ifelse(regs$region=="chr5:40345563-40623244","PTGER4",
ifelse(regs$region=="chr5:56141024-56148856","ANKRD55",
ifelse(regs$region=="chr6:90227175-90296024","BACH2",
ifelse(regs$region=="chr6:126337897-126517459","CENPW",
ifelse(regs$region=="chr6:137639876-137923679","TNFAIP3",
ifelse(regs$region=="chr6:158904748-159103161","TAGAP",
ifelse(regs$region=="chr7:26655307-27157700","SKAP2",
ifelse(regs$region=="chr7:28098574-28216621","JAZF1",
ifelse(regs$region=="chr7:50221511-51009338","COBL_IKZF1",
ifelse(regs$region=="chr9:4282536-4296430","GLIS3",
ifelse(regs$region=="chr10:6028576-6489093","IL2RA",
ifelse(regs$region=="chr10:88263276-88291560","RNLS",
ifelse(regs$region=="chr11:2157974-2177435","INS",
ifelse(regs$region=="chr11:60960006-61039249","CD5_CD6",
ifelse(regs$region=="chr12:9753255-9757853","CD69",
ifelse(regs$region=="chr12:55976127-56360038","IKZF4",
ifelse(regs$region=="chr12:111427245-111569952","SH2B3",
ifelse(regs$region=="chr13:42268918-42364338","AKAP11",
ifelse(regs$region=="chr13:99389908-99442617","GPR183",
ifelse(regs$region=="chr14:68792124-68847342","ZFP36L1",
ifelse(regs$region=="chr14:97924282-98032614","14q32.2",
ifelse(regs$region=="chr14:100831117-100855421","MEG3",
ifelse(regs$region=="chr15:78942128-78944951","CTSH",
ifelse(regs$region=="chr16:11070710-11371022","DEXI",
ifelse(regs$region=="chr16:28326722-28986790","IL27",
ifelse(regs$region=="chr16:75206985-75252586","BCAR1",
ifelse(regs$region=="chr17:39747478-40657572","GSDMB_ORMDL3_IKZF3",
ifelse(regs$region=="chr18:12773554-12884344","PTPN2",
ifelse(regs$region=="chr18:69846569-69876452","CD226",
ifelse(regs$region=="chr19:10352442-10473649","TYK2",
ifelse(regs$region=="chr19:46651227-46753285","PRKD2",
ifelse(regs$region=="chr19:48682911-48714803","FUT2",
ifelse(regs$region=="chr20:1629306-1699463","SIRPG",
ifelse(regs$region=="chr21:42400164-42444409","UBASH3A",
ifelse(regs$region=="chr21:44201347-44280294","ICOSLG",
ifelse(regs$region=="chr22:29775042-30195180","ASCC2_LIF",
ifelse(regs$region=="chr22:37148205-37256722","C1QTNF6",NA)))))))))))))))))))))))))))))))))))))))))))))))))

regs<-regs[order(regs$chromosome, regs$position),]
regs$ord<-c(1:nrow(regs))
regs$genename<-as.factor(regs$genename)
regs$genename<-reorder(regs$genename, regs$ord)
j<- regs %>%
group_by(genename, tag1) %>%
summarise(n=n(),pp=max(ppsum))
j$pp<-ifelse(j$pp>1,1,j$pp)

library(RColorBrewer)
library(ggplot2)
#highlight the UBASH3A groups for Figure 1
cols<-c(rep("black",nrow(j)))
names(cols)<-j$tag1
cols<-ifelse(names(cols)=="rs7276555","green",
ifelse(names(cols)=="rs13048049","red",
ifelse(names(cols)=="rs9984852","blue",cols)))
names(cols)<-j$tag1
j$tag2<-as.factor(j$tag1)

names(cols)<-j$tag2

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, 1),name="Group\nPosterior\nProbability")
sc1 <- scale_fill_gradientn(colours = myPalette(100), limits=c(0, 1),name="Group\nPosterior\nProbability")
names(cols)<-NULL
png(file="~/output/final_collection_gwas/guessfm_5pc_diff/cred_group_sizes_sm.png", width=10, height=5, units="cm",res=800)
ggplot(data=j, aes(tag1,n, colour=pp,fill=pp)) + geom_bar(stat="identity") +
facet_grid(. ~ genename, space="free_x", scales="free_x", switch="x") +
#theme_classic() +
  theme(#strip.placement = "outside",
        panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",size=0.2),
        strip.background = element_rect(fill=NA, colour="black",size=0.2),
        panel.spacing.x=unit(0,"cm"), axis.text.x=element_text(angle=90,size=3,colour=cols),
        strip.text.x = element_text(angle = 90,size=3),
        legend.text=element_text(size=3), legend.title=element_text(size=3), 
        axis.title.x=element_text(size=5), axis.title.y=element_text(size=5),
        axis.text.y=element_text(size=3),
	axis.ticks = element_line(size = 0.2),legend.key.size = unit(0.15, "cm"),
	legend.position=c(0.5,0.8)) +
scale_x_discrete(name="Candidate gene(s) and group index variant") +
scale_y_continuous(name="Number of variants in group") + sc + sc1
dev.off()

png(file="~/output/final_collection_gwas/guessfm_5pc_diff/cred_group_sizes.png", width=25, height=20, units="cm",res=400)
ggplot(data=j, aes(tag1,n, colour=pp,fill=pp)) + geom_bar(stat="identity") + 
facet_grid(. ~ genename, space="free_x", scales="free_x", switch="x") +
theme_classic() +
  theme(#strip.placement = "outside",
        strip.background = element_rect(fill=NA, colour="grey50"),
        panel.spacing.x=unit(0,"cm"), axis.text.x=element_text(angle=90),
	strip.text.x = element_text(angle = 90,size=8)) +
scale_x_discrete(name="Candidate gene(s) and group index variant") +
scale_y_continuous(name="Number of variants in group") + sc + sc1 
dev.off()


#and what is the predicted number of variants per region?
library(GUESSFM)
f<-function(index){
mydir <-paste0("/scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/",gsub(":",".",index),"/")
dd <- read.snpmod(mydir)
pp <- pp.nsnp(dd,plot=FALSE,expected=3, overdispersion = 1.00000001)
pps<-data.frame(nums=names(pp$trait),prob=pp$trait)
alls<-data.frame(nums=c(0:10),ps=0)
out<-merge(alls,pps, by="nums",all.x=T)
out$probs<-ifelse(is.na(out$prob),out$ps,out$prob)
o<-out[out$probs==max(out$probs),]
o<-data.frame(index, o$nums)
return(o)
}
l<-lapply(orig$Marker,f)
l<-do.call("rbind",l)

