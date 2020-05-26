#ldget_for_cassie.R

#cassie needs the allele 1 and allele 2 info for
#non credible vaiarnts in LD with the credible variants
library(GenomicRanges)
library(stringr)

tmpdir<-"/well/todd/users/jinshaw/mega/"
bcftools<-"/home/ji3x/software/bcftools/bcftools"
kavdir<-"/scratch/ji3x/full_cohort_gwas/vars/"

#read in all variants:
r<-read.table(file=paste0(tmpdir,"/all_creds_plus_ld.txt"),header=T,as.is=T)
issues<-r[is.na(r$Allele1),]



#get RSIDs where they are present:
df<-issues
o<-df[,c("chromosome","position")]
write.table(o, file=paste0(tmpdir,"/signewfr2.txt"), sep="\t",
col.names=F, row.names=F, quote=F)

system(paste0(bcftools," view -T ",tmpdir,"/signewfr2.txt ",kavdir,"/Kaviar-160204-Public-hg38-trim.vcf.gz > ",
tmpdir,"/signewfr2_ids.txt"))

vars<-read.table(file=paste0(tmpdir,"/signewfr2_ids.txt"),skip=40, header=T, comment.char="",as.is=T)
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

vars$AF<-gsub(";.*","",vars$INFO)
vars$AF1<-gsub("AF=","",vars$AF)
see<-strsplit(vars$AF1,split=",")
l<-lengths(see)
altaf<-str_split_fixed(vars$AF1, ",", max(l)+1)
colnames(altaf)<-paste0("af",c(1:ncol(altaf)))
vars<-cbind(vars,altaf)
for(i in c(1:max(l))){
vars[,paste0("af",i)]<-as.character(vars[,paste0("af",i)])
vars[,paste0("af",i)]<-ifelse(vars[,paste0("af",i)]=="",NA,vars[,paste0("af",i)])
vars[,paste0("af",i)]<-as.numeric(vars[,paste0("af",i)])
}
vars$freq<-NA
for(j in c(1:nrow(vars))){
vars[j,"freq"]<-which(c(vars[j,"af1"],
vars[j,"af2"],vars[j,"af3"],
vars[j,"af4"],vars[j,"af5"])==max(c(vars[j,"af1"],
vars[j,"af2"],vars[j,"af3"],
vars[j,"af4"],vars[j,"af5"]),na.rm=T))
}
vars$MarkerName<-ifelse(vars$freq==1,vars$MarkerName1,
ifelse(vars$freq==2,vars$MarkerName2,
ifelse(vars$freq==3,vars$MarkerName3,
ifelse(vars$freq==4,vars$MarkerName4,
ifelse(vars$freq==5,vars$MarkerName5,NA)))))

vars<-vars[vars$ID %in%	issues$ID,] 


issues1<-issues[issues$ID %in% vars$ID,]
rownames(vars)<-vars$ID
vars<-vars[issues1$ID,]
vars$allele1<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\3",vars$MarkerName)
vars$allele2<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\4",vars$MarkerName)

issues1$Allele1<-vars$allele1
issues1$Allele2<-vars$allele2


ord<-paste0(r$chromosome,":",r$position)

r1<-r[!r$ID %in% issues1$ID,]


r<-rbind(r1,issues1)
rownames(r)<-paste0(r$chromosome,":",r$position)
r<-r[ord,]
write.table(r,file=paste0(tmpdir,"/all_creds_plus_ld_a12.txt"),
col.names=T,quote=F,row.names=F,sep="\t")
#this was shared with Cassie to run annovar with
