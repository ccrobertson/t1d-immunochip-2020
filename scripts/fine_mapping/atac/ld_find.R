#ld_find.R

#finding regions in 1000 genomes with haplotypes of similar size to the MEGA credible SNPs.

library(GenomicRanges)
library(rtracklayer)

tgdir<-"/well/1000G/WTCHG/1000GP_Phase3/"
tmpdir<-"/well/todd/users/jinshaw/mega/"
progdir<-"~/programs/uva/"
ldstore<-"~/software/ldstore_v1.1_x86_64/ldstore"
qctool<-"/apps/well/qctool/2.0.1/qctool"

folder=as.character(commandArgs(trailingOnly = TRUE)[1])
stub=as.character(commandArgs(trailingOnly = TRUE)[2])


message(paste0("Folder=",folder,", stub=",stub))
#want to identify other regions across the genome with haplotypes this size (in europeans).
eur<-read.table(file=paste0(tgdir,"/1000GP_Phase3.sample"),as.is=T,header=T)
eur1<-eur[eur$GROUP=="EUR",]
ids<-eur1$ID
write.table(ids,file=paste0(tmpdir,"/eurs.txt"),
col.names=F,row.names=F,quote=F)

eur$sex=ifelse(eur$SEX=="male","M",ifelse(eur$SEX=="female","F",NA))
eur<-eur[,c("ID","sex")]
cols<-data.frame(ID="0",sex="D")

eurs<-rbind(cols,eur)
write.table(eurs, col.names=T, row.names=F, quote=F, sep=" ", file=paste0(tmpdir,"/eurs_sample.txt"))

#generate bgen for each chromosome:
dochrom<-function(chrom){
sink(file=paste0(progdir,"/ldget/chr",chrom,".sh"))
cat(paste0("#!/bin/bash

#$ -cwd -V
#$ -N chrom_",chrom, " -j y
#$ -P todd.prjc -q short.qc
#$ -pe shmem 3

OMP_NUM_THREADS=5
"))
cat(paste0(qctool," -g ",tgdir,"/ALL.chr",chrom,
".phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz -vcf-genotype-field GT -og ",tmpdir,"/chr",chrom,".bgen\n"))
#filter only europeans:
cat(paste0("echo filtering Europeans...\n"))
cat(paste0(qctool," -g ",tmpdir,"/chr",chrom,
".bgen -s ",tmpdir,"/eurs_sample.txt -incl-samples ",tmpdir,"/eurs.txt -og ",tmpdir,"/chr",chrom,"_eur.bgen\n\n"))
cat(paste0("echo exclude rare variants:\n"))
cat(paste0(qctool," -g ",tmpdir,"/chr",chrom,
"_eur.bgen -threshold 0.9 -snp-stats -osnp ",tmpdir,"/snp-stats_",chrom,".txt\n\n"))

cat(paste0("echo keep only those with MAF<0.01 so we can eventually filter them out:\n"))
cat(paste0("Rscript ",progdir,"filterthem.R ",chrom,"\n\n"))

cat(paste0(qctool," -g ",tmpdir,"/chr",chrom,
"_eur.bgen -og ",tmpdir,"/chr_",chrom,"_subsetted.bgen -excl-positions ",tmpdir,"/exclude_",chrom,".txt\n"))

cat(paste0("echo finding correlations up r 0.8...\n"))
cat(paste0(ldstore," --bgen ",tmpdir,"/chr_",chrom,
"_subsetted.bgen --bcor ",tmpdir,"/chr",chrom,".bcor --ld-thold 0.8 --n-threads 5\n"))
cat(paste0(ldstore," --bcor ",tmpdir,"/chr",chrom,".bcor --merge 5\n"))
cat(paste0(ldstore," --bcor ",tmpdir,"/chr",chrom,
".bcor --table ",tmpdir,"/chr",chrom,".matrix\n"))


cat(paste0("echo finding correlations up r 0.9...\n"))
cat(paste0(ldstore," --bgen ",tmpdir,"/chr_",chrom,
"_subsetted.bgen --bcor ",tmpdir,"/chr",chrom,"_9.bcor --ld-thold 0.9 --n-threads 5\n"))
cat(paste0(ldstore," --bcor ",tmpdir,"/chr",chrom,"_9.bcor --merge 5\n"))
cat(paste0(ldstore," --bcor ",tmpdir,"/chr",chrom,
"_9.bcor --table ",tmpdir,"/chr",chrom,"_9.matrix\n"))
sink()
system(paste0("qsub ",progdir,"/ldget/chr",chrom,".sh"))
}
lapply(c(1:22),dochrom)



readinall<-function(chrom){
get<-read.table(file=paste0(tmpdir,"/chr",chrom,".matrix"),header=T,as.is=T)
return(get)
}
allhaps<-lapply(c(1:22),readinall)
allhaps<-do.call("rbind", allhaps)

allhaps<-allhaps[allhaps$RSID1!="." & allhaps$RSID2!=".",]

library(plyr)
library(dplyr)
g<-allhaps %>%
group_by(RSID1) %>%
summarise(N=n())


save(allhaps, g, file=paste0(tmpdir,"/all_lds.RData"))


#same for the 0.9 data:
readinall9<-function(chrom){
get<-read.table(file=paste0(tmpdir,"/chr",chrom,"_9.matrix"),header=T,as.is=T)
return(get)
}
allhaps9<-lapply(c(1:22),readinall9)
allhaps9<-do.call("rbind", allhaps9)

g9<-allhaps9 %>%
group_by(RSID1) %>%
summarise(N=n())
save(allhaps9, g9, file=paste0(tmpdir,"/all_lds_9.RData"))


