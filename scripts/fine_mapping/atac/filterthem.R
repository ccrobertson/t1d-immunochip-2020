#filterthem.R
#get low frequency variants to exclude for each chromosome:


args=commandArgs(trailingOnly=TRUE)
tmpdir<-"/well/todd/users/jinshaw/mega/"

r<-read.table(file=paste0(tmpdir,"/snp-stats_",args,".txt"),header=T,as.is=T)
r<-r[r$minor_allele_frequency<0.01 | r$rsid=="." | is.na(r$rsid) | is.na(r$position),]
r<-r[!is.na(r$position),]
pos<-paste0(args,":",r$position)
write.table(pos,file=paste0(tmpdir,"/exclude_",args,".txt"),col.names=F, row.names=F,
quote=F,sep="\t")
