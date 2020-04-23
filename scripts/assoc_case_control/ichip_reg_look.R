#ichip_look_reg.R

#some exploratory analyses examining where the genotyped variants are on the iChip, but the main purpose is to
#save the ImmunoChip or densely genotyped regions on the ImmunoChip, only variants within these regions are examined in the 
#genetic discovery analysis

library(GenomicRanges)
extension="final_collection"
cohorts=c("EUR")

#define ichip regions:
ichip<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/rawdat/ichip_regions_gr38.txt"),header=T,as.is=T,sep="\t")

#and add the SNPs that were on the iChip to this set:
ichip1<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/rawdat/topmed_imputation_input_variants_b38_sorted_excludingMismatches.txt"),header=F, as.is=T)


library(ggplot2)
library(ggbio)
#plot for each chromosome, genotyped variants vs. ichip region definition:
plotchrom<-function(chrom){
chr<-ichip[ichip$chromosome==paste0("chr",chrom),]

genochr<-ichip1[ichip1$V1==chrom,]

one<-ggplot(data=chr, aes(x=chr$start, xend=chr$end, y=0, yend=0)) + geom_segment() +
theme(axis.text.y=element_blank(),
axis.title.y=element_blank(),
axis.ticks.y=element_blank())

two<-ggplot(data=genochr, aes(x=genochr$V2, xend=genochr$V2, y=1, yend=-1)) + geom_segment(size=0.01) +
theme(axis.text.y=element_blank(),
axis.title.y=element_blank(),axis.ticks.y=element_blank())


t<-tracks(one,two)
ggsave(t, file=paste0("~/output/final_collection_gwas/ichip_regs/Chr",chrom,".png"),height=10, width=20,units="cm", dpi=200)
}
g<-lapply(c(1:22), plotchrom)



#for each region, zoom in to 500kb either side and see if looks approximately correct definition.


plotreg<-function(reg){
regs<-ichip[ichip$region==reg,]
min<-regs$start-500000
max<-regs$end+500000
chrom<-as.numeric(gsub("chr","",regs$chromosome))
genochr<-ichip1[ichip1$V1==chrom,]

one<-ggplot(data=regs, aes(x=start, xend=end, y=0, yend=0)) + geom_segment() +
theme(axis.text.y=element_blank(),
axis.title.y=element_blank(), 
axis.ticks.y=element_blank()) + coord_cartesian(xlim=c(min,max))

two<-ggplot(data=genochr, aes(x=V2, xend=V2, y=1, yend=-1)) + geom_segment(size=0.01) +
theme(axis.text.y=element_blank(),
axis.title.y=element_blank(),axis.ticks.y=element_blank()) + coord_cartesian(xlim=c(min,max))


t<-tracks(one,two)
ggsave(t, file=paste0("~/output/final_collection_gwas/ichip_regs/region",reg,".png"),height=10, width=20,units="cm", dpi=200)
}
g<-lapply(c(1:nrow(ichip)), plotreg)

ichipext<-ichip
ichipext$start<-ichipext$start-50000
ichipext$end<-ichipext$end+50000

plotregplus50<-function(reg){
regs<-ichipext[ichipext$region==reg,]
min<-regs$start-500000
max<-regs$end+500000
chrom<-as.numeric(gsub("chr","",regs$chromosome))
genochr<-ichip1[ichip1$V1==chrom,]

one<-ggplot(data=regs, aes(x=start, xend=end, y=0, yend=0)) + geom_segment() +
theme(axis.text.y=element_blank(),
axis.title.y=element_blank(),
axis.ticks.y=element_blank()) + coord_cartesian(xlim=c(min,max))

two<-ggplot(data=genochr, aes(x=V2, xend=V2, y=1, yend=-1)) + geom_segment(size=0.01) +
theme(axis.text.y=element_blank(),
axis.title.y=element_blank(),axis.ticks.y=element_blank()) + coord_cartesian(xlim=c(min,max))


t<-tracks(one,two)
ggsave(t, file=paste0("~/output/final_collection_gwas/ichip_regs/region",reg,"_ext.png"),height=10, width=20,units="cm", dpi=200)
}
g<-lapply(c(1:nrow(ichip)), plotregplus50)





ichip<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/rawdat/ichip_regions_gr38.txt"),header=T,as.is=T,sep="\t")
ichipr<-GRanges(seqnames=ichip$chromosome,
IRanges(ichip$start, end=ichip$end),
region=ichip$region)

#and add the SNPs that were on the iChip to this set:
ichip1<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/rawdat/topmed_imputation_input_variants_b38_sorted_excludingMismatches.txt"),header=F, as.is=T)

ichips<-GRanges(seqnames=paste0("chr",ichip1$V1),
IRanges(ichip1$V2,end=ichip1$V2))

getsegs<-function(chrom){
chr<-ichip1[ichip1$V1==chrom,]
max<-max(chr$V2)
div<-max/500000
l<-div-round(div,digits=0)
if(l>0){
max<-max+500000
}
j<-seq(0,max,by=500000)
ob<-j[2:length(j)]
j<-j[-length(j)]
out<-data.frame(seqnames=rep(paste0("chr",chrom),length(ob)),
start=j, end=ob)
return(out)
}

segsall<-lapply(c(1:22),getsegs)
segsall<-do.call("rbind",segsall)
segsall$seg<-c(1:nrow(segsall))
obs<-GRanges(seqnames=segsall$seqnames,
IRanges(segsall$start, end=segsall$end),
seg=segsall$seg)

#get SNP density in 500 kb windows?
getdens<-function(seg){
segs<-obs[obs$seg==seg,]
chipreg<-sum(countOverlaps(ichips, segs))
dens<-chipreg/500000
df<-data.frame(seg=seg, dens=dens, sum=chipreg)
return(df)
}
densities<-lapply(c(1:nrow(segsall)),getdens)
densities<-do.call("rbind",densities)

png(file="~/output/final_collection_gwas/ichip_regs/densities_all.png",
height=20, width=20,units="cm",res=200)
ggplot(densities, aes(sum)) + geom_histogram()
dev.off()


#pick all regions with >50 variants and see how they line up with the ichip regions:
hund<-densities[densities$sum>=50,]
segsh<-segsall[segsall$seg %in% hund$seg,]

densreg<-GRanges(seqnames=segsh$seqnames,
IRanges(segsh$start,end=segsh$end),
seq=segsh$seg)

inichip<-subsetByOverlaps(densreg,ichipr)

dregtab<-as.data.frame(densreg)
dregtab$nvar=hund$sum
dregtab$chromosome<-dregtab$seqnames
dregtab$ichip<-ifelse(dregtab$seq %in% inichip$seq,T,F)
dregtab<-dregtab[,c("chromosome","start","end","nvar","ichip")]
write.table(dregtab, sep="\t", quote=F, col.names=T, row.names=F, file="~/output/final_collection_gwas/META/dens_regs.txt")

#plot for each chromosome, genotyped variants vs. ichip region definition:
chromichipdensreg<-function(chrom){
chr<-ichip[ichip$chromosome==paste0("chr",chrom),]

dreg<-as.data.frame(densreg)
genochr<-dreg[dreg$seqnames==paste0("chr",chrom),]

one<-ggplot(data=chr, aes(x=chr$start, xend=chr$end, y=0, yend=0)) + geom_segment() +
theme(axis.text.y=element_blank(),
axis.title.y=element_blank(),
axis.ticks.y=element_blank())

two<-ggplot(data=genochr, aes(x=genochr$start, xend=genochr$end, y=0, yend=0)) + geom_segment() +
theme(axis.text.y=element_blank(),
axis.title.y=element_blank(),
axis.ticks.y=element_blank())

genochr1<-ichip1[ichip1$V1==chrom,]

three<-ggplot(data=genochr1, aes(x=V2, xend=V2, y=1, yend=-1)) + geom_segment(size=0.01) +
theme(axis.text.y=element_blank(),
axis.title.y=element_blank(),axis.ticks.y=element_blank())


t<-tracks(one,two,three)
ggsave(t, file=paste0("~/output/final_collection_gwas/ichip_regs/Chr",chrom,"_dens_ichip.png"),height=10, width=20,units="cm", dpi=200)
}
g<-lapply(c(1:22), chromichipdensreg)


#dens + ichip region look
plotregboth<-function(reg){
regs<-ichipext[ichipext$region==reg,]
min<-regs$start-500000
max<-regs$end+500000
chrom<-as.numeric(gsub("chr","",regs$chromosome))
genochr<-ichip1[ichip1$V1==chrom,]

dreg<-as.data.frame(densreg)
genochr1<-dreg[dreg$seqnames==paste0("chr",chrom),]

one<-ggplot(data=regs, aes(x=start, xend=end, y=0, yend=0)) + geom_segment() +
theme(axis.text.y=element_blank(),
axis.title.y=element_blank(),
axis.ticks.y=element_blank()) + coord_cartesian(xlim=c(min,max))

onea<-ggplot(data=genochr1, aes(x=start, xend=end, y=0, yend=0)) + geom_segment() +
theme(axis.text.y=element_blank(),
axis.title.y=element_blank(),
axis.ticks.y=element_blank()) + coord_cartesian(xlim=c(min,max))

two<-ggplot(data=genochr, aes(x=V2, xend=V2, y=1, yend=-1)) + geom_segment(size=0.01) +
theme(axis.text.y=element_blank(),
axis.title.y=element_blank(),axis.ticks.y=element_blank()) + coord_cartesian(xlim=c(min,max))


t<-tracks(one,onea,two)
ggsave(t, file=paste0("~/output/final_collection_gwas/ichip_regs/both_region",reg,".png"),height=10, width=20,units="cm", dpi=200)
}

g<-lapply(c(1:nrow(ichip)), plotregboth)

dreg<-as.data.frame(densreg)
colnames(ichipext)[2]<-"seqnames"
newregs<-rbind(as.data.frame(dreg)[,c("seqnames","start","end")],ichipext[,c("seqnames","start","end")])
write.table(newregs, file="/scratch/ji3x/final_collection_gwas/rawdat/dens_and_ichip_regs.txt",sep="\t", quote=F, col.names=T, row.names=F)


chromallreg<-function(chrom){
chr<-newregs[newregs$seqnames==paste0("chr",chrom),]


one<-ggplot(data=chr, aes(x=chr$start, xend=chr$end, y=0, yend=0)) + geom_segment() +
theme(axis.text.y=element_blank(),
axis.title.y=element_blank(),
axis.ticks.y=element_blank())

genochr1<-ichip1[ichip1$V1==chrom,]

two<-ggplot(data=genochr1, aes(x=V2, xend=V2, y=1, yend=-1)) + geom_segment(size=0.01) +
theme(axis.text.y=element_blank(),
axis.title.y=element_blank(),axis.ticks.y=element_blank())


t<-tracks(one,two)
ggsave(t, file=paste0("~/output/final_collection_gwas/ichip_regs/allregs",chrom,"_ichip.png"),height=10, width=20,units="cm", dpi=200)
}
g<-lapply(c(1:22), chromallreg)

