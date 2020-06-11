#ld_overall_newb_genedens_isletup_incld.R
library(jimisc)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(plyr)
library(dplyr)

#finding regions in 1000 genomes with haplotypes of similar size to the MEGA credible SNPs.
#in addition I want them to have roughly the same gene density as iChip regions.
tmpdir<-"/well/todd/users/jinshaw/mega/"
outdir<-"/well/todd/users/jinshaw/output/uva/atac/"

creds<-read.table(file=paste0(tmpdir,"/all_creds_plus_ld.txt"),header=T, as.is=T)
creds<- creds %>%
group_by(tag) %>%
mutate(ppsums=max(ppsum,na.rm=T))

creds<-as.data.frame(creds)
creds$ppsum<-creds$ppsums
creds$ID<-ifelse(creds$ID==".",creds$MarkerName,creds$ID)

liftthem<-function (frame, chain, snp.name = "snp.name", chromosome = "chromosome",
    position = "position", updateto)
{
    if (is.character(frame[, chromosome])) {
        if (substr(1, 3, frame[1, chromosome]) == "chr")
            chr <- paste0(frame[, chromosome])
    }
    chr <- paste0("chr", frame[, chromosome])
    start <- frame[, position]
    snp <- frame[, snp.name]
    frame[, "snp.name"] <- frame[, snp.name]
    frame[, "pos"] <- c(1:nrow(frame))
    chain.file.path <- paste0("/scratch/ji3x/liftover_allign/",
        chain)
    input <- GRanges(seqname = Rle(chr), ranges = IRanges(start = start,
        end = start), snp.name = snp)
    c <- import.chain(chain.file.path)
    alligned <- unlist(liftOver(input, c))
    names(mcols(alligned)) <- "snp.name"
    updated <- data.frame(snp.name = alligned@elementMetadata@listData$snp.name)
    updated[, paste0("position", updateto)] <- alligned@ranges@start
    frame <- merge(frame, updated, by = "snp.name")
    frame <- frame[order(frame$pos), ]
    frame$pos <- NULL
    return(frame)
}

#lift over to build 37:
creds<-liftthem(creds,chain="hg38ToHg19.over.chain", snp.name = "ID", chromosome = "chromosome",
    position = "position", updateto="37")
creds<-creds[order(creds$chromosome, creds$position),]

creds<-creds[creds$ppsum>0.8,]
cred<- creds %>%
group_by(tag) %>%
summarise(min=min(position37),max=max(position37), chrom=min(chromosome), n=n())
cred$cat<-ifelse(cred$n<10,10,
ifelse(cred$n>=10 & cred$n<20,20,
ifelse(cred$n>=20 & cred$n<50,50,
ifelse(cred$n>=50 & cred$n<75,75,
ifelse(cred$n>=75 & cred$n<100,100,
ifelse(cred$n>=100 & cred$n<150,150,
ifelse(cred$n>=150 & cred$n<260,260,NA)))))))

credg<-GRanges(seqnames=paste0("chr",cred$chrom),
IRanges(cred$min, end=cred$max), tag=cred$tag,
n=cred$n, cat=cred$cat)

#get gene density of the credible snp sets:
library(humarray)
library(Homo.sapiens)
human.genes = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

credg$totgenes=countOverlaps(credg, human.genes)


#now get group size:
getgroup<-function(tag){
cr<-creds[creds$tag==tag,]
j<-data.frame(tag=tag, region=cr$region[1],
n=nrow(cr))
return(j)
}
all<-lapply(unique(creds$tag),getgroup)
all<-do.call("rbind",all)

all$cat<-ifelse(all$n<10,10,
ifelse(all$n>=10 & all$n<20,20,
ifelse(all$n>=20 & all$n<50,50,
ifelse(all$n>=50 & all$n<75,75,
ifelse(all$n>=75 & all$n<100,100,
ifelse(all$n>=100 & all$n<150,150,
ifelse(all$n>=150 & all$n<260,260,NA)))))))


#want to identify other regions across the genome with haplotypes this size.
load(file=paste0(tmpdir,"/all_lds.RData"))
g<-allhaps %>%
group_by(RSID1) %>%
summarise(N=n(),
start=min(c(position1,position2)),
end=max(c(position1,position2)),
chrom=min(chromosome))


#compare to the same chunk size:
t2<-table(all$cat)

g$cat<-ifelse(g$N<10,10,
ifelse(g$N>=10 & g$N<20,20,
ifelse(g$N>=20 & g$N<50,50,
ifelse(g$N>=50 & g$N<75,75,
ifelse(g$N>=75 & g$N<100,100,
ifelse(g$N>=100 & g$N<150,150,
ifelse(g$N>=150 & g$N<260,260,
ifelse(g$N>=260,500,NA))))))))

g1<-GRanges(seqnames=paste0("chr",g$chrom),
IRanges(g$start, g$end),
RSID1=g$RSID1,cat=g$cat, n=g$N)
g1$totgene<-countOverlaps(g1, human.genes)


cat10<-subset(g1, cat==10)
cat20<-subset(g1, cat==20)
cat50<-subset(g1, cat==50)
cat75<-subset(g1, cat==75)
cat100<-subset(g1, cat==100)
cat150<-subset(g1, cat==150)
cat260<-subset(g1, cat==260)


dosamp<-function(tag){
one<-credg[credg$tag==tag,]
if(one$cat==10){
df=cat10
}
if(one$cat==20){
df=cat20
}
if(one$cat==50){
df=cat50
}
if(one$cat==75){
df=cat75
}
if(one$cat==100){
df=cat100
}
if(one$cat==150){
df=cat150
}
if(one$cat==260){
df=cat260
}
df$diff<-abs(df$totgene-one$totgenes)
df$diffsize<-abs(df$n-one$n)
catb<-sort(df,by=~diff + diffsize)
catb<-catb[c(1:1000),]
h<-runif(1)
select=ceiling(1000*h)
out<-catb[select,]
return(out)
}
dosamp1<-function(num){
j1<-lapply(unique(creds$tag),dosamp)
j1<-Reduce(union, c(j1[[1]],j1[[2]],j1[[3]],j1[[4]],
j1[[5]],j1[[6]],j1[[7]],j1[[8]],
j1[[9]],j1[[10]],j1[[11]],j1[[12]],
j1[[13]],j1[[14]],j1[[15]],j1[[16]],
j1[[17]],j1[[18]],j1[[19]],j1[[20]],
j1[[21]],j1[[22]],j1[[23]],j1[[24]],
j1[[25]],j1[[26]],j1[[27]],j1[[28]],
j1[[29]],j1[[30]],j1[[31]],j1[[32]],
j1[[33]],j1[[34]],j1[[35]],j1[[36]],
j1[[37]],j1[[38]],j1[[39]],j1[[40]],
j1[[41]],j1[[42]],j1[[43]],j1[[44]],
j1[[45]],j1[[46]],j1[[47]],j1[[48]],
j1[[49]],j1[[50]],j1[[51]],j1[[52]],
j1[[53]],j1[[54]],j1[[55]],j1[[56]],
j1[[57]],j1[[58]],j1[[59]],j1[[60]],
j1[[61]],j1[[62]],j1[[63]],j1[[64]],
j1[[65]],j1[[66]],j1[[67]],j1[[68]],
j1[[69]],j1[[70]],j1[[71]],j1[[72]],
j1[[73]]))
j1$num=num
message(paste0("Done ",num))
return(j1)
}

js<-mclapply(c(1:100),dosamp1,mc.cores=4)


#now generate the dataset of variants in roughly equivalent LD. Save this rather than ALL the SNPs in LD so that:
#a) the same set of SNPs are sampled for each cell type and 
#b) won't use up loads of memory

savethesamps<-function(num){
allhaps1<-allhaps[allhaps$RSID1 %in% js[[num]]$RSID1,]
allhaps1$chrpos1<-paste0(allhaps1$chromosome,":",allhaps1$position1)
allhaps1$chrpos2<-paste0(allhaps1$chromosome,":",allhaps1$position2)

allhaps2<-data.frame(RSID=c(unique(allhaps1$RSID1),allhaps1$RSID2),
chrpos=c(unique(allhaps1$chrpos1),allhaps1$chrpos2))
allhaps2$chromosome<-sub("(^.*)[:](.*)","\\1",allhaps2$chrpos)
allhaps2$position<-sub("(^.*)[:](.*)","\\2",allhaps2$chrpos)

allhaps2$position<-as.numeric(allhaps2$position)
allhaps2$chromosome<-as.numeric(allhaps2$chromosome)
allhaps2<-liftthem(allhaps2,chain="hg19ToHg38.over.chain", snp.name = "RSID", chromosome = "chromosome",
    position = "position", updateto="38")
allhaps2$position37<-allhaps2$position
allhaps2$position<-allhaps2$position38
message(paste0("Done ",num))
return(allhaps2)
}
sampsthem<-mclapply(c(1:100),savethesamps,mc.cores=6)
save(js, sampsthem,file=paste0(tmpdir,"/ld_samples_genedens_incld.RData"))

#load the random samples of 100 in:
load(file=paste0(tmpdir,"/ld_samples_genedens_incld.RData"))

#create the GRanges object for the credible SNPs
cr1<-GRanges(seqnames=paste0("chr",creds$chromosome),
IRanges(creds$position,end=creds$position))

lookatthem<-function(loc,folder,stub){
peaks<-read.table(file=paste0(loc,folder,"/",stub,".bed.gz"),as.is=T)
#GRanges object for the ATAC-seq peaks
p<-GRanges(seqnames=peaks$V1,
IRanges(peaks$V2,end=peaks$V3))

m<-subsetByOverlaps(cr1,p)
reg<-nrow(as.data.frame(m))
n<-nrow(creds)

#now perform the enrichment analyses:
getenrich<-function(num,m,n,reg){
#GRanges object for the randomly sampled:
crnull<-GRanges(seqnames=paste0("chr",sampsthem[[num]]$chromosome),
IRanges(sampsthem[[num]]$position,end=sampsthem[[num]]$position))

mnull<-subsetByOverlaps(crnull, p)
regnull<-nrow(as.data.frame(mnull))
nnull<-nrow(sampsthem[[num]])

mat<-matrix(c(reg,regnull,n,nnull),ncol=2)
f<-fisher.test(mat)
z<-ifelse(f$estimate>1,
qnorm(f$p.value/2,lower.tail=F),
ifelse(f$estimate<1,
qnorm(f$p.value/2,lower.tail=F)*-1,0))
out<-data.frame(stub=stub,num=num,overlap=reg, denom=n,
overlapnull=regnull, denomnull=nnull, z=z, p=f$p.value)
message(paste0("Done ",num," - ",stub)) 
return(out)
}
outs<-lapply(c(1:100),getenrich,m=m,n=n,reg=reg)
outs<-do.call("rbind",outs)
outs$stub<-paste0(folder,"_",stub)
outs1<-data.frame(stub=paste0(folder,"-",stub),
overlap=outs$overlap[1], denom=outs$denom[1],
overlapmean_null=mean(outs$overlapnull), 
denommean_null=mean(outs$denomnull),
zmean=mean(outs$z))
outs1$p=2*pnorm(-abs(outs1$zmean))
message(paste0("Done ",folder," - ",stub))
return(outs1)
}


types=c("Bulk_B-S","Effector_CD4pos_T-U", "Immature_NK-U","Memory_Teffs-U","Naive_B-U", "Plasmablasts-U", "Th2_precursors-S",
"Bulk_B-U","Effector_memory_CD8pos_T-S","Mature_NK-S", "Memory_Tregs-S", "Naive_CD8_T-S", "Regulatory_T-S", "Th2_precursors-U",
"CD8pos_T-S","Effector_memory_CD8pos_T-U","Mature_NK-U", "Memory_Tregs-U", "Naive_CD8_T-U", "Regulatory_T-U",
"CD8pos_T-U","Follicular_T_Helper-S", "Mem_B-S","Monocytes-S", "Naive_Teffs-S", "Th17_precursors-S", "pDCs-U",
"Central_memory_CD8pos_T-S", "Follicular_T_Helper-U", "Mem_B-U", "Monocytes-U", "Naive_Teffs-U",  "Th17_precursors-U",
"Central_memory_CD8pos_T-U", "Gamma_delta_T-S", "Memory_NK-U", "Myeloid_DCs-U", "Naive_Tregs-S",  "Th1_precursors-S",
"Effector_CD4pos_T-S", "Gamma_delta_T-U", "Memory_Teffs-S",  "Naive_B-S", "Naive_Tregs-U",  "Th1_precursors-U",
"ACT_IL4IL21_new","NACT_IL4IL21_new","ACT","NACT","HumanCardiacFibroblast_Adult","HumanCardiacFibroblast_Fetal","EndoC-BetaH1_control",
"EndoC-BetaH1_cytokine-treated","Human-islet_control","Human-islet_cytokine-treated")

folders=c(rep("Stanford",times=45),rep("BCell",times=2),rep("CD4Cambridge",times=2), rep("HFC",times=2),
rep("Islets_new",times=4))

outs<-mcmapply(lookatthem,folder=folders, stub=types,loc=c(rep(paste0(tmpdir,"/atac"),55)),
 SIMPLIFY=FALSE,mc.cores=6)
outs<-do.call("rbind",outs)
outs$stub<-as.character(outs$stub)
outs$stim<-ifelse(substr(outs$stub,nchar(outs$stub),nchar(outs$stub))=="S","Stim",
ifelse(grepl("-ACT",outs$stub),"Stim",
ifelse(grepl("_cytokine-treated",outs$stub),"Stim","Unstim")))
outs$stubs<-gsub("-S","",outs$stub)
outs$stubs<-gsub("-U","",outs$stubs)
outs$stubs<-gsub("-ACT_new","",outs$stubs)
outs$stubs<-gsub("-NACT_new","",outs$stubs)
outs$stubs<-gsub("-ACT","",outs$stubs)
outs$stubs<-gsub("-NACT","",outs$stubs)
outs$stubs<-gsub("_IL4IL21_new","",outs$stubs)
outs$stubs<-gsub("_IL4IL21_new","",outs$stubs)
outs$stubs<-gsub("_cytokine-treated","",outs$stubs)
outs$stubs<-gsub("_control","",outs$stubs)

outs$stubs<-ifelse(outs$stubs=="BCell","Oxford-BCell",
ifelse(outs$stubs=="CD4Cambridge","Oxford-CD4",
ifelse(outs$stubs=="HFC-HumanCardiacFibroblast_Adult","Jonsson-Cardiac_Fibroblasts-a",
ifelse(outs$stubs=="HFC-HumanCardiacFibroblast_Fetal","Jonsson-Cardiac_Fibroblasts-f",
ifelse(outs$stubs=="Islets_new-EndoC-BetaH1", "Ramos-Rodriguez-EndoC_BetaH1",
ifelse(outs$stubs=="Islets_new-Human-islet", "Ramos-Rodriguez-Islets", outs$stubs))))))

outs$stubs<-gsub("Stanford","Calderon",outs$stubs)




g<-ggplot(data=outs, aes(y=zmean, x=as.factor(stubs),fill=stim)) + geom_bar(stat="identity",position=position_dodge(width=1)) + 
coord_flip() + geom_hline(aes(yintercept=1.959964)) + scale_y_continuous(name="Mean enrichment z-score from 100 iterations") + 
scale_x_discrete(name="Cell type") +
theme(panel.background = element_blank(),
      axis.line=element_line())

ggsave(g, file=paste0(outdir,"enrich_all_newb_genedens_isletup_incld.png"),
width=20, height=20, units="cm", dpi=300)
ggsave(g, file=paste0(outdir,"/enrich_all_newb_genedens_isletup_incld.pdf"),
width=7, height=8)
outs$logp<-log10(outs$p)*-1
bon<-0.05/nrow(outs)
bonline<-log10(bon)*-1
outs<-outs[order(-outs$logp),]
outs$ords<-c(1:nrow(outs))
outs <- outs %>%
group_by(stubs) %>%
mutate(m=min(ords))
outs$stub<-as.factor(outs$stubs)
outs$stub<-reorder(outs$stub,-outs$m)
g1<-ggplot(outs, aes(y=logp, x=as.factor(stub),fill=stim)) + geom_bar(stat="identity",position=position_dodge(width=1)) + 
coord_flip() + geom_hline(aes(yintercept=bonline),colour="red",linetype="dashed") + 
scale_y_continuous(name=expression(paste("-",log[10]," enrichment p-value from 100 iterations")),breaks=c(0,2,4,6,8)) +
scale_x_discrete(name="Cell type") +
scale_fill_discrete(name="Condition", labels=c("Stimulated", "Unstimulated")) +
theme(panel.background = element_blank(),
      axis.line=element_line())

ggsave(g1, file=paste0(outdir,"/enrich_all_newb_p_genedens_isletup_incld.png"),
width=20, height=20, units="cm", dpi=300)
ggsave(g1, file=paste0(outdir,"/enrich_all_newb_p_genedens_isletup_incld.pdf"),
width=7, height=8)
save(outs, file=paste0(outdir,"/enrichments_all_newb_genedens_isletup_incld.RData"))




#get the list of variants for each cell type:
creds$variant<-ifelse(creds$ID==".",creds$snp.name, creds$ID)
cr1<-GRanges(seqnames=paste0("chr",creds$chromosome),
IRanges(creds$position,end=creds$position),variant=creds$variant)

getsnps<-function(loc,folder, stub){
peaks<-read.table(file=paste0(loc,folder,"/",stub,".bed.gz"),as.is=T)
#GRanges object for the ATAC-seq peaks
p<-GRanges(seqnames=peaks$V1,
IRanges(peaks$V2,end=peaks$V3))

m<-subsetByOverlaps(cr1,p)
mout<-as.data.frame(m)
write.table(mout, file=paste0(outdir,"/snps/",folder,"-",stub,"_incld.txt"),
sep="\t", col.names=T,row.names=F, quote=F)
}



types=c("Bulk_B-S","Effector_CD4pos_T-U", "Immature_NK-U","Memory_Teffs-U","Naive_B-U", "Plasmablasts-U", "Th2_precursors-S",
"Bulk_B-U","Effector_memory_CD8pos_T-S","Mature_NK-S", "Memory_Tregs-S", "Naive_CD8_T-S", "Regulatory_T-S", "Th2_precursors-U",
"CD8pos_T-S","Effector_memory_CD8pos_T-U","Mature_NK-U", "Memory_Tregs-U", "Naive_CD8_T-U", "Regulatory_T-U",
"CD8pos_T-U","Follicular_T_Helper-S", "Mem_B-S","Monocytes-S", "Naive_Teffs-S", "Th17_precursors-S", "pDCs-U",
"Central_memory_CD8pos_T-S", "Follicular_T_Helper-U", "Mem_B-U", "Monocytes-U", "Naive_Teffs-U",  "Th17_precursors-U",
"Central_memory_CD8pos_T-U", "Gamma_delta_T-S", "Memory_NK-U", "Myeloid_DCs-U", "Naive_Tregs-S",  "Th1_precursors-S",
"Effector_CD4pos_T-S", "Gamma_delta_T-U", "Memory_Teffs-S",  "Naive_B-S", "Naive_Tregs-U",  "Th1_precursors-U",
"ACT_IL4IL21_new","NACT_IL4IL21_new","ACT","NACT","HumanCardiacFibroblast_Adult","HumanCardiacFibroblast_Fetal",
"EndoC-BetaH1_control","EndoC-BetaH1_cytokine-treated","Human-islet_control","Human-islet_cytokine-treated")

folders=c(rep("Stanford",times=45),rep("BCell",times=2),rep("CD4Cambridge",times=2), rep("HFC",times=2),
rep("Islets_new",times=4))

mcmapply(getsnps,folder=folders, stub=types, loc=c(rep(paste0(tmpdir,"/atac/"),55)),mc.cores=6)





#now generate a figure showing whether any/how many variants are in open chromatin in each group in each cell type
cred<- creds %>%
group_by(tag) %>%
summarise(min=min(position),max=max(position), chrom=min(chromosome), n=n())

getregsums<-function(loc,folder,stub){
peaks<-read.table(file=paste0(loc,folder,"/",stub,".bed.gz"),as.is=T)
#GRanges object for the ATAC-seq peaks
p<-GRanges(seqnames=peaks$V1,
IRanges(peaks$V2,end=peaks$V3))

doforall<-function(tag,p, stub, folder){
cr1<-creds[creds$tag==tag,]
cr1<-GRanges(seqnames=paste0("chr",cr1$chrom),
IRanges(cr1$position,end=cr1$position))

m<-subsetByOverlaps(cr1,p)
reg<-nrow(as.data.frame(m))
out<-data.frame(tag=tag,stub=paste0(folder,"_",stub),number=reg)
return(out)
}
outs1<-lapply(cred$tag,doforall, p=p, stub=stub, folder=folder)
outs1<-do.call("rbind",outs1)

return(outs1)
}

types=c("Bulk_B-S","Effector_CD4pos_T-U", "Immature_NK-U","Memory_Teffs-U","Naive_B-U", "Plasmablasts-U", "Th2_precursors-S",
"Bulk_B-U","Effector_memory_CD8pos_T-S","Mature_NK-S", "Memory_Tregs-S", "Naive_CD8_T-S", "Regulatory_T-S", "Th2_precursors-U",
"CD8pos_T-S","Effector_memory_CD8pos_T-U","Mature_NK-U", "Memory_Tregs-U", "Naive_CD8_T-U", "Regulatory_T-U",
"CD8pos_T-U","Follicular_T_Helper-S", "Mem_B-S","Monocytes-S", "Naive_Teffs-S", "Th17_precursors-S", "pDCs-U",
"Central_memory_CD8pos_T-S", "Follicular_T_Helper-U", "Mem_B-U", "Monocytes-U", "Naive_Teffs-U",  "Th17_precursors-U",
"Central_memory_CD8pos_T-U", "Gamma_delta_T-S", "Memory_NK-U", "Myeloid_DCs-U", "Naive_Tregs-S",  "Th1_precursors-S",
"Effector_CD4pos_T-S", "Gamma_delta_T-U", "Memory_Teffs-S",  "Naive_B-S", "Naive_Tregs-U",  "Th1_precursors-U",
"ACT_IL4IL21_new","NACT_IL4IL21_new","ACT","NACT","HumanCardiacFibroblast_Adult","HumanCardiacFibroblast_Fetal","EndoC-BetaH1_control",
"EndoC-BetaH1_cytokine-treated","Human-islet_control","Human-islet_cytokine-treated")

folders=c(rep("Stanford",times=45),rep("BCell",times=2),rep("CD4Cambridge",times=2), rep("HFC",times=2),
rep("Islets_new",times=4))

outs<-mcmapply(getregsums,folder=folders, stub=types,loc=c(rep(paste0(tmpdir,"/atac/",55)), SIMPLIFY=FALSE,mc.cores=6)
outs<-do.call("rbind",outs)


outs$stub<-as.character(outs$stub)
outs$stim<-ifelse(substr(outs$stub,nchar(outs$stub),nchar(outs$stub))=="S","Stim",
ifelse(grepl("-ACT",outs$stub),"Stim",
ifelse(grepl("_cytokine-treated",outs$stub),"Stim","Unstim")))
outs$stubs<-gsub("-S","",outs$stub)
outs$stubs<-gsub("-U","",outs$stubs)
outs$stubs<-gsub("-ACT_new","",outs$stubs)
outs$stubs<-gsub("-NACT_new","",outs$stubs)
outs$stubs<-gsub("-ACT","",outs$stubs)
outs$stubs<-gsub("-NACT","",outs$stubs)
outs$stubs<-gsub("_IL4IL21_new","",outs$stubs)
outs$stubs<-gsub("_IL4IL21_new","",outs$stubs)
outs$stubs<-gsub("_cytokine-treated","",outs$stubs)
outs$stubs<-gsub("_control","",outs$stubs)

outs$stubs<-ifelse(outs$stubs=="BCell","Oxford-BCell",
ifelse(outs$stubs=="CD4Cambridge","Oxford-CD4",
ifelse(outs$stubs=="HFC-HumanCardiacFibroblast_Adult","Jonsson-Cardiac_Fibroblasts-a",
ifelse(outs$stubs=="HFC-HumanCardiacFibroblast_Fetal","Jonsson-Cardiac_Fibroblasts-f",
ifelse(outs$stubs=="Islets_new-EndoC-BetaH1", "Ramos-Rodriguez-EndoC_BetaH1",
ifelse(outs$stubs=="Islets_new-Human-islet", "Ramos-Rodriguez-Islets", outs$stubs))))))

outs$stubs<-paste0(outs$stubs,"-",outs$stim)

outs$yesno<-ifelse(outs$number>0,1,0)
p <- ggplot(outs, aes(tag, stubs)) + geom_tile(aes(fill = as.factor(number))) + 
theme(axis.text.x=element_text(angle=90))
ggsave(p, file=paste0(outdir,"/groups_open_up.png"),height=20, width=30,units="cm",dpi=400)

p <- ggplot(outs, aes(tag, stubs)) + geom_tile(aes(fill = as.factor(yesno))) +
theme(axis.text.x=element_text(angle=90))
ggsave(p, file=paste0(outdir,"/groups_open_up_yesno.png"),height=20, width=30,units="cm",dpi=400)
