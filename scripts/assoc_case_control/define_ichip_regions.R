#define_ichip_regions.R

#read in ichip regions hg19 and liftover all SNPs... define the regions in hg38:
library(GenomicRanges)
library(rtracklayer)
library(humarray)

#load the b36 regions
data(iChipRegionsB36)
#now load the SNPs we used as input for the imputation:
extension="final_collection"

snps<-read.table(file="/nv/vol185/MEGA/pre_qc/Mega_9April2019/data_derived/mega_rawD.bim", header=F, as.is=T)
colnames(snps)<-c("chromosome","snp.name","distance","position37","a1","a2")
#lift over to b36:

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

#get snps in all build positions:
snps<-liftthem(snps,"hg19ToHg18.over.chain",position="position37",updateto="36")
snps<-liftthem(snps,"hg19ToHg38.over.chain",position="position37", updateto="38")

#make it a GRanges object:
snps1<-GRanges(paste0("chr",seqnames=snps$chromosome),
ranges=IRanges(snps$position36, end=snps$position36),
snpname=snps$snp.name,
position37=snps$position37,
position38=snps$position38)

ichip<-as(iChipRegionsB36,"GRanges")
elementMetadata(ichip)[[ "region" ]] <- c(1:188)

out<-mergeByOverlaps(ichip,snps1)
out<-as.data.frame(out)

regiondef<-function(reg){
o<-out[out$region==reg,]
m<-min(o$position38)
ma<-max(o$position38)
ch<-o$ichip.seqnames[1]
t<-data.frame(region=reg, chromosome=ch, start=m, end=ma,nsnp=nrow(o))
return(t)
}
regions<-lapply(c(1:188),regiondef)
regions<-do.call("rbind", regions)

write.table(regions, file=paste0("/scratch/ji3x/",extension,"_gwas/rawdat/ichip_regions_gr38.txt"), col.names=T, row.names=F, quote=F, sep="\t")
