#paintor_readin_5pc_diff_up.R
#Reads in the paintor results and produces plots
library(ggplot2)
library(ggbio)

extension="final_collection"

orig<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/pmeta_fam_ind_risk_regions_vcf_5pc_diff.txt"),header=T,sep=",",as.is=T)
orig<-orig[orig$ichip=="yes",]
#remove duplicate due to long LD block:
orig<-orig[!(orig$Marker %in% c("chr12:112153882:G:A","chr12:112730563:A:G","chr7:50406053:G:A",
"chr1:113302202:A:C","chr17:46751565:G:A","chr14:68286876:C:T","chr12:111066229:T:C",
"chr17:40623690:G:A")),]

orig$snpnum=c(1:nrow(orig))

orig$gene<-ifelse(orig$Marker=="chr1:92358141:G:C","RPAP2",
ifelse(orig$Marker=="chr1:113834946:A:G","PTPN22",
ifelse(orig$Marker=="chr1:172746562:A:G","FASLG",
ifelse(orig$Marker=="chr1:192570207:G:A","RGS1",
ifelse(orig$Marker=="chr1:206769068:C:T","IL10",
ifelse(orig$Marker=="chr2:100147438:G:T","AFF3",
ifelse(orig$Marker=="chr2:162254026:A:G","IFIH1",
ifelse(orig$Marker=="chr2:191105394:C:G","STAT4",
ifelse(orig$Marker=="chr2:203874196:G:A","CTLA4",
ifelse(orig$Marker=="chr3:46342713:G:A","CCR5",
ifelse(orig$Marker=="chr4:26083889:G:A","RBPJ",       
ifelse(orig$Marker=="chr4:122322439:A:ATC","IL2/IL21",
ifelse(orig$Marker=="chr5:35888101:G:A","IL7R",
ifelse(orig$Marker=="chr5:40521603:G:T","TTC33",
ifelse(orig$Marker=="chr5:56146422:T:C","ANKRD55",
ifelse(orig$Marker=="chr6:424915:C:A","IRF4",
ifelse(orig$Marker=="chr6:90267049:G:A","BACH2",
ifelse(orig$Marker=="chr6:126343187:C:T","CENPW",
ifelse(orig$Marker=="chr6:137682468:T:C","TNFAIP3",
ifelse(orig$Marker=="chr6:159049210:G:T","TAGAP",
ifelse(orig$Marker=="chr7:26852706:G:A","SKAP2",
ifelse(orig$Marker=="chr7:28102567:G:T","JAZF1",
ifelse(orig$Marker=="chr7:50970505:G:A","IKZF1/COBL",
ifelse(orig$Marker=="chr9:4283137:G:T","GLIS3",
ifelse(orig$Marker=="chr10:6052734:C:T","IL2RA",
ifelse(orig$Marker=="chr10:88287593:T:G","RNLS",
ifelse(orig$Marker=="chr11:2163618:G:T","INS",
ifelse(orig$Marker=="chr11:60961822:G:T","CD5/CD6",
ifelse(orig$Marker=="chr12:9757568:G:A","CD69",
ifelse(orig$Marker=="chr12:56050848:C:T","IKZF4",
ifelse(orig$Marker=="chr12:111569952:C:T","SH2B3",
ifelse(orig$Marker=="chr13:42343795:C:T","AKAP11",
ifelse(orig$Marker=="chr13:99411911:G:A","GPR183",
ifelse(orig$Marker=="chr14:68792124:C:G","ZFP36L1",
ifelse(orig$Marker=="chr14:98032614:A:G","-",
ifelse(orig$Marker=="chr14:100840110:T:C","MEG3",NA))))))))))))))))))))))))))))))))))))
orig$gene<-ifelse(orig$Marker=="chr15:38670555:G:T","RASGRP1",
ifelse(orig$Marker=="chr15:78944951:C:T","CTSH",
ifelse(orig$Marker=="chr16:11097543:G:A","DEXI",
ifelse(orig$Marker=="chr16:28496748:C:T","IL27",
ifelse(orig$Marker=="chr16:75212137:G:C","BCAR1",
ifelse(orig$Marker=="chr17:39910119:T:C","GSDMB/CCR7",
ifelse(orig$Marker=="chr18:12818923:G:A","PTPN2",
ifelse(orig$Marker=="chr18:69869260:C:T","CD226",
ifelse(orig$Marker=="chr19:10317045:T:A","TYK2",
ifelse(orig$Marker=="chr19:46705224:T:C","PRKD2",
ifelse(orig$Marker=="chr19:48703417:G:A","FUT2",
ifelse(orig$Marker=="chr20:1634898:T:C","SIRPG",
ifelse(orig$Marker=="chr21:42416077:G:A","UBASH3A",
ifelse(orig$Marker=="chr21:44204668:C:T","ICOSLG",
ifelse(orig$Marker=="chr22:30055497:C:T","ASCC2/LIF",
ifelse(orig$Marker=="chr22:37189696:T:C","C1QTNF6",orig$gene))))))))))))))))

load(file=paste0("/scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/forscript.RData"))


getmanhats<-function(number){
if(file.exists(paste0("/scratch/ji3x/final_collection_gwas/input_paintor_5pc_diff_3"))){
three<-read.table(file="/scratch/ji3x/final_collection_gwas/input_paintor_5pc_diff_3",header=F,as.is=T)
}
if(file.exists(paste0("/scratch/ji3x/final_collection_gwas/input_paintor_5pc_diff_1_and_2"))){
oneandtwo<-read.table(file="/scratch/ji3x/final_collection_gwas/input_paintor_5pc_diff_1_and_2",header=F,as.is=T)
}
if(file.exists(paste0("/scratch/ji3x/final_collection_gwas/input_paintor_5pc_diff_1_and_3"))){
oneandthree<-read.table(file="/scratch/ji3x/final_collection_gwas/input_paintor_5pc_diff_1_and_3",header=F,as.is=T)
}
init<-read.table(file=paste0("/scratch/ji3x/final_collection_gwas/paintor_5pc_diff/Locus",number),header=T, as.is=T)
paintor<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/paintor_results_5pc_diff/Locus",number,".results"),header=T,as.is=T)
top<-paintor[paintor$Posterior_Prob==max(paintor$Posterior_Prob),]
top<-top[abs(top$ZSCORE.P1)==max(abs(top$ZSCORE.P1)),]
top<-top[1,]

snp<-orig[number,"Marker"]
system(paste0("~/software/ldstore_v1.1_x86_64/ldstore --bcor /scratch/ji3x/",
extension,"_gwas/paintor_5pc_diff/EUR/",gsub(":",".",snp),".bcor_1 --matrix /scratch/ji3x/",
extension,"_gwas/paintor_5pc_diff/EUR/",gsub(":",".",snp),".matrix"))
system(paste0("~/software/ldstore_v1.1_x86_64/ldstore --bcor /scratch/ji3x/",
extension,"_gwas/paintor_5pc_diff/EUR/",gsub(":",".",snp),".bcor_1 --meta /scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/EUR/",gsub(":",".",snp),".meta"))
snps<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/EUR/",gsub(":",".",snp),".meta"),header=T,as.is=T)
mat<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/EUR/",gsub(":",".",snp),".matrix"))
colnames(mat)<-snps$RSID
rownames(mat)<-snps$RSID

ld<-mat[,colnames(mat)==top$MarkerName]
names(ld)<-snps$RSID
ld<-ld[init$MarkerName]

one<-ggplot(data=init, aes(position,abs(ZSCORE.P1), colour=abs(ld))) + geom_point() +
scale_color_gradient(low = "black", high = "red", name="EUR LD to lead variant") +
scale_y_continuous(name="Absolute z-score (EUR)")

ch<-grep(paste0("Locus",number,"$"),oneandtwo$V1)
ch1<-grep(paste0("Locus",number,"$"),oneandthree$V1)
ch2<-grep(paste0("Locus",number,"$"),three$V1)
if(length(ch)==0 & length(ch1)==0 & length(ch2)==0){
}
if(length(ch)==1 | length(ch2)==1){
system(paste0("~/software/ldstore_v1.1_x86_64/ldstore --bcor /scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/AFR/",gsub(":",".",snp),".bcor_1 --matrix ",
"/scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/AFR/",gsub(":",".",snp),".matrix"))
system(paste0("~/software/ldstore_v1.1_x86_64/ldstore --bcor /scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/AFR/",gsub(":",".",snp),".bcor_1 --meta ",
"/scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/AFR/",gsub(":",".",snp),".meta"))
snpsa<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/AFR/",gsub(":",".",snp),".meta"),header=T,as.is=T)
mata<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/AFR/",gsub(":",".",snp),".matrix"))
colnames(mata)<-snpsa$RSID
rownames(mata)<-snpsa$RSID

ld<-mata[,colnames(mata)==top$MarkerName]
init$lda<-ld
two<-ggplot(data=init, aes(position,abs(ZSCORE.P2),colour=abs(lda))) + geom_point() +
scale_color_gradient(low = "black", high = "red", name="AFR LD to lead variant") +
scale_y_continuous(name="Absolute z-score (AFR)")
}
if(length(ch1)==1 | length(ch2)==1){
system(paste0("~/software/ldstore_v1.1_x86_64/ldstore --bcor /scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/FIN/",gsub(":",".",snp),".bcor_1 --matrix ",
"/scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/FIN/",gsub(":",".",snp),".matrix"))
system(paste0("~/software/ldstore_v1.1_x86_64/ldstore --bcor /scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/FIN/",gsub(":",".",snp),".bcor_1 --meta ",
"/scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/FIN/",gsub(":",".",snp),".meta"))
snpsf<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/FIN/",gsub(":",".",snp),".meta"),header=T,as.is=T)
matf<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/paintor_5pc_diff/FIN/",gsub(":",".",snp),".matrix"))
colnames(matf)<-snpsf$RSID
rownames(matf)<-snpsf$RSID

ld<-matf[,colnames(matf)==top$MarkerName]
init$ldf<-ld
three<-ggplot(data=init, aes(position,abs(ZSCORE.P3),colour=abs(ldf))) + geom_point() + 
scale_color_gradient(low = "black", high = "red", name="FIN LD to lead variant") +
scale_y_continuous(name="Absolute z-score (FIN)")
}
ld<-mat[,colnames(mat)==top$MarkerName]
names(ld)<-snps$RSID
ld<-ld[paintor$MarkerName]
paintor$ld1<-ld

paint<-ggplot(data=paintor, aes(position, Posterior_Prob,colour=abs(ld1))) + geom_point() +
scale_color_gradient(low = "black", high = "red", name="EUR LD to lead variant") +
scale_y_continuous(name="PAINTOR posterior probability")

if(length(ch2)==0 & length(ch1)==0 & length(ch)==0){
t<-tracks(one,paint)
}
if(length(ch2)==1){
t<-tracks(one,two,three,paint)
}
if(length(ch)==1){
t<-tracks(one,two,paint)
}
if(length(ch1)==1){
t<-tracks(one,three,paint)
}
ggsave(t, file=paste0("~/output/final_collection_gwas/paintor_5pc_diff/",snp,"_manhat_paint_dens.png"),
height=25, width=20, units="cm", dpi=200)
message(paste0("Done ",snp))
}
lapply(forscript$locusnum,getmanhats)

#generate supplemtary tables for manuscript:
gensup<-function(number){
snp=orig[orig$snpnum==number,"Marker"]
guessfm<-read.table(file=paste0("/home/ji3x/output/final_collection_gwas/guessfm_5pc_diff/",gsub(":",".",snp),
"/credible_snps.txt"),header=T, as.is=T)
paintor<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/paintor_results_5pc_diff/Locus",number,".results"),header=T,as.is=T)

gu<-merge(guessfm,paintor,by="MarkerName",all.x=T)
gu<-gu[order(-gu$Posterior_Prob),]
if("ZSCORE.P2" %in% colnames(gu) & !"ZSCORE.P3" %in% colnames(gu)){
gu<-gu[,c("MarkerName","ID","position.x",
"pp","ppsum","ZSCORE.P1", "ZSCORE.P2","Posterior_Prob")]

colnames(gu)<-c("Marker","rsid","position","ppguessfm","ppsumguessfm","zeur","zafr","pppaintor")
gu$ppguessfm<-round(gu$ppguessfm,digits=3)
gu$ppsumguessfm<-round(gu$ppsumguessfm,digits=3)
gu$pppaintor<-round(gu$pppaintor,digits=3)
gu$zafr<-round(gu$zafr,digits=3)
gu$zeur<-round(gu$zeur,digits=3)
write.table(gu,file=paste0("~/output/",extension,"_gwas/paintor_5pc_diff/comp_",snp,"_dens"),
col.names=T,row.names=F,quote=F,sep="\t")
}
if(!"ZSCORE.P2" %in% colnames(gu) & "ZSCORE.P3"	%in% colnames(gu)){
gu<-gu[,c("MarkerName","ID","position.x",
"pp","ppsum","ZSCORE.P1", "ZSCORE.P3","Posterior_Prob")]

colnames(gu)<-c("Marker","rsid","position","ppguessfm","ppsumguessfm","zeur","zfin","pppaintor")
gu$ppguessfm<-round(gu$ppguessfm,digits=3)
gu$ppsumguessfm<-round(gu$ppsumguessfm,digits=3)
gu$pppaintor<-round(gu$pppaintor,digits=3)
gu$zfin<-round(gu$zfin,digits=3)
gu$zeur<-round(gu$zeur,digits=3)
write.table(gu,file=paste0("~/output/",extension,"_gwas/paintor_5pc_diff/comp_",snp,"_dens"),
col.names=T,row.names=F,quote=F,sep="\t")
}
if("ZSCORE.P2" %in% colnames(gu) & "ZSCORE.P3" %in% colnames(gu)){
gu<-gu[,c("MarkerName","ID","position.x",
"pp","ppsum","ZSCORE.P1","ZSCORE.P2", "ZSCORE.P3","Posterior_Prob")]

colnames(gu)<-c("Marker","rsid","position","ppguessfm","ppsumguessfm","zeur","zafr","zfin","pppaintor")
gu$ppguessfm<-round(gu$ppguessfm,digits=3)
gu$ppsumguessfm<-round(gu$ppsumguessfm,digits=3)
gu$pppaintor<-round(gu$pppaintor,digits=3)
gu$zafr<-round(gu$zafr,digits=3)
gu$zfin<-round(gu$zfin,digits=3)
gu$zeur<-round(gu$zeur,digits=3)
write.table(gu,file=paste0("~/output/",extension,"_gwas/paintor_5pc_diff/comp_",snp,"_dens"),
col.names=T,row.names=F,quote=F,sep="\t")
}
}
lapply(forscript$locusnum,gensup)
