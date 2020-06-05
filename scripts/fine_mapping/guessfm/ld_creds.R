#ld_creds.R

library(GenomicRanges)
library(rtracklayer)

extension="final_collection"
creddir<-paste0("~/output/",extension,"_gwas/guessfm_5pc_diff/")
tmpdir<-"/well/todd/users/jinshaw/mega/"

#read in the MEGA credible SNPs:
creds<-read.table(file=paste0(creddir,"/all_creds_idcorr.txt"),header=T, as.is=T)


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


creds<-liftthem(creds,chain="hg38ToHg19.over.chain", snp.name = "MarkerName", chromosome = "chromosome",
    position = "position", updateto="37")
creds<-creds[order(creds$chromosome, creds$position),]



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
ifelse(all$n>=150 & all$n<250,250,NA)))))))

#load all LDs with Pearson correlation>0.9
load(file=paste0(tmpdir,"/all_lds_9.RData"))
allhaps9$snp1<-paste0(allhaps9$chromosome,":",allhaps9$position1)
allhaps9$snp2<-paste0(allhaps9$chromosome,":",allhaps9$position2)
creds$snp<-paste0(creds$chromosome,":",creds$position37)
creds$index<-gsub("\\.",":",creds$tag)
creds$ta<-ifelse(creds$snp.name==creds$index,creds$snp,NA)

dota<-function(snp){
t<-creds[creds$tag==snp,]
t<-t[order(t$ta),]
t$ta<-ifelse(is.na(t$ta),t$ta[1],t$ta)
return(t)
}
cred<-lapply(unique(creds$tag),dota)
cred<-do.call("rbind",cred)

getlds<-function(snp){
t<-cred[cred$tag==snp,]
t1<-allhaps9[allhaps9$snp1 %in% t$ta,]
t2<-allhaps9[allhaps9$snp2 %in%	t$ta,]
if(nrow(t1)>0 | nrow(t2)>0){
others<-t1[!duplicated(t1$RSID1),]
others2<-t1[!duplicated(t1$RSID2),]
others3<-t2[!duplicated(t2$RSID1),]
others4<-t2[!duplicated(t2$RSID2),]
all<-data.frame(ID=c(others$RSID1,
others2$RSID2,others3$RSID1,others4$RSID2),
chromosome=c(others$chromosome,
others2$chromosome,others3$chromosome,others4$chromosome),
position=c(others$position1,
others2$position2,others3$position1,others4$position2))
all<-all[!duplicated(all$ID),]
all<-all[!is.na(all$ID),]
all<-all[!all$ID %in% t$ID,]
if(nrow(all)>0){
all<-liftthem(all,chain="hg19ToHg38.over.chain", snp.name = "ID", chromosome = "chromosome",
    position = "position", updateto="38")
all$position=all$position38
all$tag=snp
}
}
if(nrow(t1)==0 & nrow(t2)==0){
all<-NULL
}
message(paste0("Done ",snp))
return(all)
}

j<-mclapply(unique(cred$tag),getlds,mc.cores=6)
j<-do.call("rbind",j)
j$excluded<-TRUE
#remove those variants included in the stochastic search and not prioritised:
load(file="/well/todd/users/jinshaw/mega/resultsmeta_5pcs_vcf_5pc_diff.RData")
res$chrpos<-paste0("chr",res$chromosome,":",res$position)
#keep variants not excluded from EUR collection:
res<-res[substr(res$Direction,1,1)!="?",]
j$chrpos<-paste0("chr",j$chromosome,":",j$position38)
j<-j[!j$chrpos %in% res$chrpos,]

j<-j[,c("ID","chromosome","position","tag","excluded")]
creds$excluded=FALSE

creds<-merge(creds,j,by=c("ID","chromosome","position","tag","excluded"),all=T)
creds<-creds[order(creds$tag,creds$chromosome,creds$position),]

docred<-function(snp){
t<-creds[creds$tag==snp,]
t<-t[order(t$region),]
t$region<-ifelse(is.na(t$region),t$region[1],t$region)
return(t)
}
cred1<-lapply(unique(creds$tag),docred)
cred1<-do.call("rbind",cred1)

cred1<-cred1[,c("region","MarkerName","ID","tag","Allele1","Allele2","Effect","StdErr","P.value","Direction","chromosome","position","pp","ppsum","excluded")]

cred1<-cred1[order(cred1$chromosome, cred1$position),]
library(dplyr)

number<-length(unique(cred1$tag))
assignnum<-function(tag,num){
g<-cred1[cred1$tag==tag,]
g$id=num
return(g)
}
cred2<-mapply(assignnum, tag=unique(cred1$tag),num=c(1:number),SIMPLIFY=FALSE)
cred2<-do.call("rbind",cred2)

cred2<-cred2[order(cred2$id, cred2$chromosome, cred2$position),]

cred2<-cred2[,c("region","MarkerName","ID","tag","Allele1","Allele2","Effect","StdErr","P.value","Direction","chromosome","position","pp","ppsum","excluded")]
write.table(cred2,file=paste0(tmpdir,"/all_creds_plus_ld.txt"),col.names=T,row.names=F,
sep="\t",quote=F)

