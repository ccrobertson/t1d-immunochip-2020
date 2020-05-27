#ldget_for_cassie.R

#cassie needs the allele 1 and allele 2 info for
#non credible vaiarnts in LD with the credible variants
library(GenomicRanges)
library(stringr)

tmpdir<-"/well/todd/users/jinshaw/mega/"
bcftools<-"/home/ji3x/software/bcftools/bcftools"
kavdir<-"/scratch/ji3x/full_cohort_gwas/vars/"
outdir<-"~/output/final_collection_gwas/guessfm_5pc_diff/"

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



#bar graph summarising the number of variuants per group per region - currently used in Figure 1 in the manuscript
d<-read.table(file=paste0(tmpdir,"/all_creds_plus_ld_a12.txt"),sep="\t", header=T)
getppsum<-function(tag){
t<-d[d$tag==tag,]
m<-max(t$ppsum,na.rm=T)
t$ppsum<-ifelse(t$excluded==TRUE,m,t$ppsum)
return(t)
}
d<-lapply(unique(d$tag),getppsum)
d<-do.call("rbind",d)

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
png(file=paste0(outdir,"/cred_group_sizes_sm_incld.png"), width=10, height=5, units="cm",res=800)
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

png(file=paste0(outdir,"/cred_group_sizes_incld.png"), width=25, height=20, units="cm",res=400)
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

