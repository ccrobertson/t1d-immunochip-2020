#expr_stims_newb_genedens_isletup_incld.R

#creates a list of FDR<0.01 differentially open peaks then tests enrichment of cell types within stimulation or non-stimulation specific peaks.

library(jimisc)
library(GenomicRanges)
library(rtracklayer)
library(DESeq2)
library(ggplot2)
library(plyr)
library(dplyr)

tmpdir<-"/well/todd/users/jinshaw/mega/"
outdir<-"/well/todd/users/jinshaw/output/uva/atac/"

#read in the MEGA credible SNPs:
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


creds<-liftthem(creds,chain="hg38ToHg19.over.chain", snp.name = "ID", chromosome = "chromosome",
    position = "position", updateto="37")
creds<-creds[order(creds$chromosome, creds$position),]
creds<-creds[creds$ppsum>0.8,]


counts<-read.table(file=paste0(tmpdir,"/ATACseq/counts_all.csv.gz"),header=T,as.is=T,sep=",")
counts<-counts[substr(counts$Position,1,2)!="GL" & substr(counts$Position,1,2)!="KI",]
counts$chromosome<-sub("(.*)[:](.*)[-](.*)","\\1",counts$Position)
counts$posstart<-sub("(.*)[:](.*)[-](.*)","\\2",counts$Position)
counts$posstop<-sub("(.*)[:](.*)[-](.*)","\\3",counts$Position)

counts$posstart<-as.numeric(counts$posstart)
counts$posstop<-as.numeric(counts$posstop)

regions<-GRanges(seqnames=counts$chromosome,
IRanges(counts$posstart,end=counts$posstop),
name=counts$Position)


credi<-GRanges(seqnames=paste0("chr",creds$chromosome),
IRanges(creds$position,end=creds$position),
region=creds$region,snp.name=creds$snp.name,
ID=creds$ID)


#just look at the Stanford + in-house generated stuff to start with:
rownames(counts)<-counts$Position
stan<-counts[,grepl("Stanford|Bcell|CD4cam|ISislets",colnames(counts))]
stan<-stan[,!grepl("old",colnames(stan))]
#remove B cells stimulated with different stimulus:
stan<-stan[,!colnames(stan) %in% c("Bcell.S15500030_ACT_new","Bcell.S15500030_NACT_new","Bcell.S15500037_ACT_new","Bcell.S15500037_NACT_new")]

#for each cell type and each peak, examinine differential expression between stim and unstim conditions
#also naive (unstim) vs memory (unstim) cell types:

coldata<-data.frame(o=colnames(stan))
library(stringr)
coldata$o<-as.character(coldata$o)
coldata$stim<-grepl("\\.ACT|_ACT|\\.S$|treated",coldata$o)
coldata$cell<-sub("(^.*)[\\.](.*)[\\.](.*)[\\.](.*)","\\3",coldata$o)
coldata$cell<-ifelse(substr(coldata$cell,1,6)=="CD4cam","In-house_CD4",
ifelse(substr(coldata$cell,1,5)=="Bcell","In-house_B",
ifelse(substr(coldata$cell,1,6)=="BetaH1"| substr(coldata$cell,16,21)=="BetaH1","EndoC_BetaH1",
ifelse(substr(coldata$cell,1,5)=="islet"| substr(coldata$cell,16,20)=="islet","Islets",coldata$cell))))


coldata$ID<-sub("(^.*)[\\.](.*)[\\.](.*)[\\.](.*)","\\2",coldata$o)
coldata$ID<-ifelse(grepl("Donor1",coldata$ID),"Donor1",
ifelse(grepl("Donor2",coldata$ID),"Donor2",coldata$ID))
coldata$ID<-ifelse(coldata$cell=="In-house_B",sub("(^.*)[\\.](.*)[_](.*)[_](.*)","\\2",coldata$ID),coldata$ID)
coldata$ID<-ifelse(coldata$cell=="EndoC_BetaH1",sub("(^.*)[_](.*)[_](.*)","\\3",coldata$o),coldata$ID)
coldata$ID<-ifelse(coldata$cell=="Islets",sub("(^.*)[_](.*)[_](.*)","\\3",coldata$o),coldata$ID)

rownames(coldata)<-coldata$o


#for each cell type and region, test differential openess:
#(can only do this when two conditions present):
keepboth<-function(cell){
p<-coldata[coldata$cell==cell,]
t<-table(p$stim)
if(length(t)==2){
p$test=TRUE
}
if(length(t)!=2){
p$test=FALSE
}
return(p)
}

coldata<-lapply(unique(coldata$cell),keepboth)
coldata<-do.call("rbind",coldata)
coldata<-coldata[coldata$test==TRUE,]


testdeseq<-function(cell){
names<-rownames(stan)
coldat1<-coldata[grepl(paste0("^",cell),coldata$cell),]
stan1<-stan[,rownames(coldat1)]
dds <- DESeqDataSetFromMatrix(countData = stan1,
                              colData = coldat1,
                              design= ~ stim)

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds)
res<-as.data.frame(res)
res$celltype=cell
res$region=rownames(res)
res<-res[!is.na(res$log2FoldChange),]
res<-res[!is.na(res$pvalue),]
res$pfdr<-p.adjust(res$pvalue,method="BH")
res<-res[res$pfdr<0.01,]
res$start<-as.numeric(sub("(^.*)[:](.*)[-](.*)","\\2",res$region))
res$stop<-as.numeric(sub("(^.*)[:](.*)[-](.*)","\\3",res$region))
res$enrich<-ifelse(res$log2FoldChange>0,TRUE,ifelse(res$log2FoldChange<0,FALSE,NA))
ranges<-GRanges(seqnames=sub("(^.*)[:](.*)[-](.*)","\\1",res$region),
IRanges(res$start, end=res$stop),enrich=res$enrich)

message(paste0("Done ",cell))
return(ranges)
}
results<-lapply(unique(coldata$cell),testdeseq)
names(results)<-unique(coldata$cell)

#now using these bed files, want to examine enrichemnt looking specifically at the differentially-open regions between states.
stimget<-function(cell,enrich){
stim<-results[[cell]][results[[cell]]$enrich==enrich,]
return(stim)
}
stims<-lapply(unique(coldata$cell),stimget, enrich=TRUE)
names(stims)<-unique(coldata$cell)
unstims<-lapply(unique(coldata$cell),stimget, enrich=FALSE)
names(unstims)<-unique(coldata$cell)
save(stims, unstims, file=paste0(tmpdir,"/enrichments/stimunstim_newb_isletup.RData"))

load(file=paste0(tmpdir,"/enrichments/stimunstim_newb_isletup.RData"))
#number of differentially-associated peaks by cell type:
getnum<-function(cell){
s<-stims[[cell]]
u<-unstims[[cell]]
out<-data.frame(cell=cell, stim=nrow(as.data.frame(s)),
unstim=nrow(as.data.frame(u)))
out$percstim=(out$stim/nrow(counts))*100
out$percunstim=(out$unstim/nrow(counts))*100
return(out)
}
numbs<-lapply(unique(coldata$cell),getnum)
numbs<-do.call("rbind",numbs)
numbs$cell<-as.character(numbs$cell)
numbs$stub<-ifelse(numbs$cell=="In-house_CD4","Oxford-CD4",
ifelse(numbs$cell=="In-house_B","Oxford-B",
ifelse(numbs$cell=="EndoC_BetaH1","Ramos-Rodriguez-EndoC_BetaH1",
ifelse(numbs$cell=="Islets","Ramos-Rodriguez-Islets",paste0("Calderon-",numbs$cell)))))
numbs<-numbs[!(numbs$stim==0 & numbs$unstim==0),]
numbs<-numbs[order(-numbs$stim),]

sink(file=paste0(outdir,"/stim_unstim_diff_peaks_newb_isletup.txt"))
cat(paste0("Cell type;Activation specific peaks (%);Non activation specific peaks (%)\n"))
for (i in c(1:nrow(numbs))){
cat(paste0(numbs$stub[i],";",numbs$stim[i]," (",round(numbs$percstim[i],digits=2),"%);",
numbs$unstim[i]," (",round(numbs$percunstim[i],digits=2),"%)\n"))
}
sink()

load(file=paste0(tmpdir,"/ld_samples_genedens_incld.RData"))

#create the GRanges object for the credible SNPs
cr1<-GRanges(seqnames=paste0("chr",creds$chromosome),
IRanges(creds$position,end=creds$position))

lookatthem<-function(cell, stimunstim){

#GRanges object for the ATAC-seq peaks
p<-stimunstim[[cell]]

m<-subsetByOverlaps(cr1,p)
reg<-nrow(as.data.frame(m))
n<-nrow(creds)

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
out<-data.frame(cell=cell,num=num,overlap=reg, denom=n,
overlapnull=regnull, denomnull=nnull, z=z, p=f$p.value)
out$peaks<-nrow(as.data.frame(p))
message(paste0("Done ",num," - ",cell))
return(out)
}
outs<-lapply(c(1:100),getenrich,m=m,n=n,reg=reg)
outs<-do.call("rbind",outs)
outs$cell<-paste0(cell)
outs1<-data.frame(stub=paste0(cell),peaks=outs$peaks[1],
overlap=outs$overlap[1], denom=outs$denom[1],
overlapmean_null=mean(outs$overlapnull),
denommean_null=mean(outs$denomnull),
zmean=mean(outs$z))
outs1$p=2*pnorm(-abs(outs1$zmean))
message(paste0("Done - ",cell))
return(outs1)
}


outs<-mclapply(unique(coldata$cell[coldata$cell %in% numbs$cell]),stimunstim=stims,lookatthem,mc.cores=6)
outs<-do.call("rbind",outs)
outsun<-mclapply(unique(coldata$cell[coldata$cell %in% numbs$cell]),stimunstim=unstims,lookatthem,mc.cores=6)
outsun<-do.call("rbind",outsun)

bon<-0.05/nrow(outs)
bonline<-log10(bon)*-1
library(gridExtra)
outs$stub<-ifelse(outs$stub=="In-house_CD4","Oxford-CD4",
ifelse(outs$stub=="In-house_B","Oxford-BCell",
ifelse(outs$stub=="EndoC_BetaH1","Ramos-Rodriguez-EndoC_BetaH1",
ifelse(outs$stub=="Islets","Ramos-Rodriguez-Islets",paste0("Calderon-",outs$stub)))))
outs$stim="Stimulation-specific"

outsun$stub<-ifelse(outsun$stub=="In-house_CD4","Oxford-CD4",
ifelse(outsun$stub=="In-house_B","Oxford-BCell",
ifelse(outsun$stub=="EndoC_BetaH1","Ramos-Rodriguez-EndoC_BetaH1",
ifelse(outsun$stub=="Islets","Ramos-Rodriguez-Islets",paste0("Calderon-",outsun$stub)))))
outsun$stim="Unstimulated-specific"

both<-rbind(outs,outsun)
bon<-0.05/nrow(both)
bonline<-log10(bon)*-1
save(outs,outsun,both, file=paste0(tmpdir,"/atac/enrichments_all_stim_unstim_peaks_newb_genedens_isletup_incld.RData"))


both$logp<-log10(both$p)*-1
both<-both[order(-both$logp),]
both$ords<-c(1:nrow(both))
both <- both %>%
group_by(stub) %>%
mutate(m=min(ords))
both$stub<-as.factor(both$stub)
both$stub<-reorder(both$stub,-both$m)
g<-ggplot(data=both, aes(y=logp, x=as.factor(stub),fill=as.factor(stim))) + geom_bar(stat="identity",position=position_dodge(width=1)) +
coord_flip(ylim=c(0,9)) + geom_hline(aes(yintercept=bonline),colour="red",linetype="dashed") + scale_y_continuous(name=expression(paste("-",log[10]," enrichment p-value from 100 iterations")),
breaks=c(0,2,4,6,8)) +
scale_x_discrete(name="Cell type") +
scale_fill_discrete(name="Specificity") +
theme(axis.line=element_line(),
panel.background = element_blank())

ggsave(g, file=paste0(outdir,"/enrich_all_stim_unstim_peaks_logp_newb_genedens_isletup_incld.png"),
width=20, height=20, units="cm", dpi=300)
ggsave(g, file=paste0(outdir,"/enrich_all_stim_unstim_peaks_logp_newb_genedens_isletup_incld.pdf"),
width=7, height=8)

#get stimulation specific overlaps:
load(file=paste0(tmpdir,"/enrichments/stimunstim_newb_isletup.RData"))
creds$variant<-ifelse(creds$ID==".",creds$snp.name, creds$ID)
cr1<-GRanges(seqnames=paste0("chr",creds$chromosome),
IRanges(creds$position,end=creds$position),variant=creds$variant)


getsnps<-function(cell, stimunstim){
#GRanges object for the ATAC-seq peaks
p<-stimunstim[[cell]]

m<-subsetByOverlaps(cr1,p)
mout<-as.data.frame(m)
if(as.data.frame(p[1,])[6]==TRUE){
write.table(mout, file=paste0(outdir,"/snps/stim/",cell,"_newb_isletup_incld.txt"),
sep="\t", col.names=T,row.names=F, quote=F)
}
if(as.data.frame(p[1,])[6]==FALSE){
write.table(mout, file=paste0(outdir,"/snps/unstim/",cell,"_newb_isletup_incld.txt"),
sep="\t", col.names=T,row.names=F, quote=F)
}
message(paste0("Done ",cell))
}

getsnps("In-house_B",stimunstim=stims)
getsnps("In-house_B",stimunstim=unstims)





#examine how many stimlation-specific peak in each group:
load(file=paste0(tmpdir,"/enrichments/stimunstim_newb_isletup.RData"))
creds<-read.table(file=paste0(tmpdir,"/all_creds_plus_ld.txt"),header=T, as.is=T)
creds<- creds %>%
group_by(tag) %>%
mutate(ppsums=max(ppsum,na.rm=T))

creds<-as.data.frame(creds)
creds$ppsum<-creds$ppsums
creds$ID<-ifelse(creds$ID==".",creds$MarkerName,creds$ID)

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


cred<- creds %>%
group_by(tag) %>%
summarise(min=min(position),max=max(position), chrom=min(chromosome), n=n())



getregsums<-function(set,stub){
if(set=="stim"){
peaks<-stims[[stub]]
}
if(set=="unstim"){
peaks<-unstims[[stub]]
}


doforall<-function(tag,peaks, stub){
cr1<-creds[creds$tag==tag,]
cr1<-GRanges(seqnames=paste0("chr",cr1$chrom),
IRanges(cr1$position,end=cr1$position))

m<-subsetByOverlaps(cr1,peaks)
reg<-nrow(as.data.frame(m))
out<-data.frame(tag=tag,stub=paste0(stub),number=reg)
return(out)
}
outs1<-lapply(cred$tag,doforall, peaks=peaks, stub=stub)
outs1<-do.call("rbind",outs1)

return(outs1)
}

types<-names(stims)

outs<-mcmapply(getregsums, stub=types,set=c(rep("stim",length(types))), SIMPLIFY=FALSE,mc.cores=6)
outs<-do.call("rbind",outs)



outs$stubs<-ifelse(outs$stub=="In-house_CD4","Oxford-CD4",
ifelse(outs$stub=="In-house_B","Oxford-BCell",
ifelse(outs$stub %in% c("Naive_Teffs","Memory_Tregs","Bulk_B","Central_memory_CD8pos_T",
"CD8pos_T","Th1_precursors","Effector_CD4pos_T","Naive_B","Monocytes","Gamma_delta_T",
"Naive_CD8_T","Follicular_T_Helper","Regulatory_T","Effector_memory_CD8pos_T",
"Th2_precursors","Memory_Teffs","Mature_NK","Th17_precursors","Mem_B",
"Naive_Tregs"),paste0("Calderon-",outs$stub),
ifelse(outs$stub %in% c("EndoC_BetaH1","Islets"),paste0("Ramos-Rodriguez-",outs$stub),NA))))

outs$yesno<-ifelse(outs$number>0,1,0)

p <- ggplot(outs, aes(tag, stubs)) + geom_tile(aes(fill = as.factor(number))) +
theme(axis.text.x=element_text(angle=90))
ggsave(p, file=paste0(outdir,"/groups_open_stim_up.png"),height=20, width=30,units="cm",dpi=400)

p <- ggplot(outs, aes(tag, stubs)) + geom_tile(aes(fill = as.factor(yesno))) +
theme(axis.text.x=element_text(angle=90))
ggsave(p, file=paste0(outdir,"/groups_open_stim_up_yesno.png"),height=20, width=30,units="cm",dpi=400)


outs<-mcmapply(getregsums, stub=types,set=c(rep("unstim",length(types))), SIMPLIFY=FALSE,mc.cores=6)
outs<-do.call("rbind",outs)



outs$stubs<-ifelse(outs$stub=="In-house_CD4","Oxford-CD4",
ifelse(outs$stub=="In-house_B","Oxford-BCell",
ifelse(outs$stub %in% c("Naive_Teffs","Memory_Tregs","Bulk_B","Central_memory_CD8pos_T",
"CD8pos_T","Th1_precursors","Effector_CD4pos_T","Naive_B","Monocytes","Gamma_delta_T",
"Naive_CD8_T","Follicular_T_Helper","Regulatory_T","Effector_memory_CD8pos_T",
"Th2_precursors","Memory_Teffs","Mature_NK","Th17_precursors","Mem_B",
"Naive_Tregs"),paste0("Calderon-",outs$stub),
ifelse(outs$stub %in% c("EndoC_BetaH1","Islets"),paste0("Ramos-Rodriguez-",outs$stub),NA))))


outs$yesno<-ifelse(outs$number>0,1,0)
p <- ggplot(outs, aes(tag, stubs)) + geom_tile(aes(fill = as.factor(number))) +
theme(axis.text.x=element_text(angle=90))
ggsave(p, file=paste0(outdir,"groups_open_unstim_up.png"),height=20, width=30,units="cm",dpi=400)


p <- ggplot(outs, aes(tag, stubs)) + geom_tile(aes(fill = as.factor(yesno))) +
theme(axis.text.x=element_text(angle=90))
ggsave(p, file=paste0(outdir,"/groups_open_unstim_up_yesno.png"),height=20, width=30,units="cm",dpi=400)
