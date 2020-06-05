#paintor_readin_5pc_diff_up.R
#Reads in the paintor results and produces plots
library(ggplot2)
library(ggbio)
library(GenomicRanges)
library(rtracklayer)
library(RColorBrewer)
library(biomaRt)
library(Gviz)




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

#script to add gene track:
ensembl=useMart("ENSEMBL_MART_ENSEMBL")
datasets=listDatasets(ensembl)
human_ensembl = useDataset("hsapiens_gene_ensembl", mart=ensembl)

geneTrack = function(region, gene_list=NULL, colour=NULL) {

  #define colour
  if (is.null(colour)) {
    geneColor = "#999999"
  } else {
    geneColor = colour
  }

  #get exons annotations
  exonData = data.frame(geneTrack <- Gviz::BiomartGeneRegionTrack(biomart = human_ensembl, genome = "hg38", chromosome = as.numeric(gsub("chr","",region$chr)), start = region$start, end = region$end, name = "ENSEMBL", transcriptAnnotation="symbol")@range)
  exonData$exon_uniq = paste(exonData$exon, exonData$feature, sep=":")
  exonData$features_in_transcript = sapply(exonData$transcript, function(x) {sum(exonData$transcript==x)})
  exonData$longest_transcript = sapply(exonData$symbol, function(x) {d = exonData[exonData$symbol==x,c("transcript","features_in_transcript")]; d[order(d$features_in_transcript, decreasing=TRUE),]$transcript[1]})

  #subset exons for plotting
  exons = exonData[exonData$transcript==exonData$longest_transcript & exonData$feature %in% c("lncRNA","miRNA","protein_coding","utr3","utr5"),]
  exons$plotStart = as.numeric(sapply(exons$start, function(x) max(x,region[["start"]])))
  exons$plotEnd = sapply(exons$end, function(x) min(x,region[["end"]]))
  exons$width = exons$plotEnd - exons$plotStart
  exons$ymin = -0.25
  exons$ymax = 0.25
  exons$ymin[exons$feature=="protein_coding"] <- -0.5
  exons$ymax[exons$feature=="protein_coding"] <- 0.5

  #create genes data frame
  genes = data.frame(
    symbol=unique(exons$symbol),
    t(sapply(unique(exons$symbol), function(x){d=exons[exons$symbol==x,]; return(c(min(d$start, na.rm=TRUE), max(d$end, na.rm=TRUE)))})),
    strand=sapply(unique(exons$symbol), function(x){d=exons[exons$symbol==x,]; return(d$strand[1])}))
  names(genes) <- c("symbol","start","end", "strand")
  genes$plotStart = as.numeric(sapply(genes$start, function(x) max(x,region[["start"]])))
  genes$plotEnd = sapply(genes$end, function(x) min(x,region[["end"]]))
  genes$width = genes$plotEnd - genes$plotStart

  #filter based on boundaries or provided gene list
  if (is.null(gene_list)) {
    exons$include = exons$width>0
    genes$include = genes$width>0
  } else {
    exons$include = exons$symbol%in%gene_list & exons$width>0
    genes$include = genes$symbol%in%gene_list & genes$width>0
  }
  exons = exons[exons$include,]
  genes = genes[genes$include,]

  #create arrowtick data frame
  arrowtickpos = seq(region[['start']],region[['end']], length=20)
  arrows = NULL
  for (i in 1:dim(genes)[1]) {
    if (sum(arrowtickpos>genes$plotStart[i] & arrowtickpos<genes$plotEnd[i])>0) {
      d = data.frame(symbol=genes$symbol[i], strand=genes$strand[i],start=arrowtickpos[arrowtickpos>genes$plotStart[i] & arrowtickpos<genes$plotEnd[i]])
      if (nrow(d)>3) {
        d_trimmed = d[2:(nrow(d)-1),]
      } else {
        d_trimmed = d
      }
      arrows = rbind(arrows, d_trimmed)
    }
  }
  arrows$end[arrows$strand=="-"] = arrows$start[arrows$strand=="-"]-1
  arrows$end[arrows$strand=="+"] = arrows$start[arrows$strand=="+"]+1

  if(!is.null(arrows)) {
  ggplot() +
    #plot exons
    geom_rect(data=exons, aes(xmin=plotStart, xmax=plotEnd, ymin=ymin, ymax=ymax), colour=geneColor, fill=geneColor, alpha=1) +
    #connect exons
    geom_rect(data=genes, aes(xmin=plotStart, xmax=plotEnd, ymin=-0.005, ymax=0.005), colour=geneColor, fill=geneColor, alpha=1) +
    #add arrow ticks
    geom_segment(data=arrows, aes(x=start, xend=end, y=0, yend=0), lineend="butt", linejoin="mitre", arrow = arrow(angle=20, length = unit(0.5, "lines"), ends="last", type="closed"), colour=geneColor) +
    #tweak aesthetics
    scale_x_continuous(limits = c(region[["start"]],region[["end"]]), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1,1), expand = c(0, 0)) +
    facet_grid(symbol~., margins=FALSE, switch="y") +
    ylab(" ") + xlab(" ") +
    theme(
      strip.text.y = element_blank(),
      strip.background = element_blank(),
      panel.spacing.y=unit(0,"lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid=element_line(size=0.5),
      panel.background = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      legend.position = "none")
  } else {
    ggplot() +
      #plot exons
      geom_rect(data=exons, aes(xmin=plotStart, xmax=plotEnd, ymin=ymin, ymax=ymax), colour=geneColor, fill=geneColor, alpha=1) +
      #connect exons
      geom_rect(data=genes, aes(xmin=plotStart, xmax=plotEnd, ymin=-0.005, ymax=0.005), colour=geneColor, fill=geneColor, alpha=1) +
      #tweak aesthetics
      scale_x_continuous(limits = c(region[["start"]],region[["end"]]), expand = c(0, 0)) +
      scale_y_continuous(limits = c(-1,1), expand = c(0, 0)) +
      facet_grid(symbol~., margins=FALSE, switch="y") +
      ylab(" ") + xlab(" ") +
      theme(
        strip.text.y = element_blank(),
        strip.background = element_blank(),
        panel.spacing.y=unit(0,"lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid=element_line(size=0.5),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")
  }
}






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

reg<-list(chr=init$chr, start=min(init$position), end=max(init$position))
genes<-geneTrack(reg)

one<-ggplot(data=init, aes(position,abs(ZSCORE.P1), colour=abs(ld))) + geom_point(size=0.2) +
scale_color_gradient(low = "black", high = "red", name="EUR LD to \nlead variant") +
scale_y_continuous(name="Absolute z-score (EUR)") + 
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid=element_line(size=0.5),
      panel.background = element_blank(),
      axis.line.x = element_line(size=0.5),
      axis.line.y = element_line(size=0.5),
      axis.text.y=element_text(size=4),
      axis.title.y=element_text(size=6),
      legend.title=element_text(size=6),
      legend.text=element_text(size=4),
      legend.key.size=unit(0.6,"lines")) 

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
two<-ggplot(data=init, aes(position,abs(ZSCORE.P2),colour=abs(lda))) + geom_point(size=0.2) +
scale_color_gradient(low = "black", high = "red", name="AFR LD to \nlead variant") +
scale_y_continuous(name="Absolute z-score (AFR)") +
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid=element_line(size=0.5),
      panel.background = element_blank(),
      axis.line.x = element_line(size=0.5),
      axis.line.y = element_line(size=0.5),
      axis.text.y=element_text(size=4),
      axis.title.y=element_text(size=6),
      legend.title=element_text(size=6), 
      legend.text=element_text(size=4),
      legend.key.size=unit(0.6,"lines"))
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
three<-ggplot(data=init, aes(position,abs(ZSCORE.P3),colour=abs(ldf))) + geom_point(size=0.2) + 
scale_color_gradient(low = "black", high = "red", name="FIN LD to \nlead variant") +
scale_y_continuous(name="Absolute z-score (FIN)") + 
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid=element_line(size=0.5),
      panel.background = element_blank(),
      axis.line.x = element_line(size=0.5),
      axis.line.y = element_line(size=0.5),
      axis.text.y=element_text(size=4),
      axis.title.y=element_text(size=6),
      legend.title=element_text(size=6), 
      legend.text=element_text(size=4),
      legend.key.size=unit(0.6,"lines"))
}
ld<-mat[,colnames(mat)==top$MarkerName]
names(ld)<-snps$RSID
ld<-ld[paintor$MarkerName]
paintor$ld1<-ld

chrpos_ticks_breaks = seq(from = min(init$position), to = max(init$position), by = round((max(init$position)-min(init$position))/4, digits=0))
chrpos_ticks_labels = paste0(round(chrpos_ticks_breaks/1e6, digits=3),"Mb")
paint<-ggplot(data=paintor, aes(position, Posterior_Prob,colour=abs(ld1))) + geom_point(size=0.2) +
scale_color_gradient(low = "black", high = "red", name="EUR LD to \nlead variant") +
scale_y_continuous(name="PAINTOR posterior probability") + 
scale_x_continuous(breaks=chrpos_ticks_breaks, labels=chrpos_ticks_labels) +
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid=element_line(size=0.5),
      panel.background = element_blank(),
      axis.line.x = element_line(size=0.5),
      axis.line.y = element_line(size=0.5),
      axis.text.x=element_text(size=4),
      axis.text.y=element_text(size=4),
      axis.title.y=element_text(size=6),
      legend.title=element_text(size=6), 
      legend.text=element_text(size=4),
      legend.key.size=unit(0.6,"lines"))

if(length(ch2)==0 & length(ch1)==0 & length(ch)==0){
t<-tracks(genes,one,paint, heights=c(2,4,4)) + theme(axis.ticks.x=element_line(size=0.1))
}
if(length(ch2)==1){
t<-tracks(genes,one,two,three,paint, heights=c(1,2,2,2,3)) + theme(axis.ticks.x=element_line(size=0.1))
}
if(length(ch)==1){
t<-tracks(genes,one,two,paint,heights=c(1,3,3,3)) + theme(axis.ticks.x=element_line(size=0.1))
}
if(length(ch1)==1){
t<-tracks(genes,one,three,paint,heights=c(1,3,3,3)) + theme(axis.ticks.x=element_line(size=0.1))
}
ggsave(t, file=paste0("~/output/final_collection_gwas/paintor_5pc_diff/",snp,"_manhat_paint_dens.png"),
height=10, width=7.5, units="cm", dpi=800)
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









#overlay on publically available ATAC seq data
#doing this in depth for the RBPJ region (figure 2 in the manuscript).
getAtacTrack = function(sample, dir, region) {
  regionGR <- GRanges(c(region[["chr"]]), IRanges(start=region[["start"]], end=region[["end"]]))
  bw = file.path(dir, paste0(sample,".bw"))
  bw_reg = data.frame(import.bw(bw, which = regionGR))
  names(bw_reg) <- c("chr", "start","end","width","strand","score")
  bw_reg$type = sample
  return(bw_reg)
}

extractAtacRegion = function(samples, dir, region) {
  tracks = lapply(samples, getAtacTrack, dir, region)
  dat_reg = do.call("rbind",tracks)
  dat_reg$type = factor(dat_reg$type, levels=samples)
  dat_reg$type_Num = as.numeric(dat_reg$type)
  return(dat_reg)
}




pubdir<-"/nv/vol185/T1DGC/PublicData/"
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

finemap<-read.table(file=paste0("~/output/",extension,"_gwas/guessfm_5pc_diff/all_creds_idcorr.txt"),sep="\t", header=T)

finemap_dat = data.frame(chr=paste0("chr",finemap$chromosome), start=finemap$position, end=finemap$position, ppsum=finemap$ppsum, tag=finemap$tag, marker=finemap$ID)
finemap_reg = finemap_dat[finemap_dat$tag=="chr4.26108491.G.A",]
finemap_reg$chr<-as.character(finemap_reg$chr)

stanford_dir = paste0(pubdir,"/Calderon/atac_data/processed_by_Flores/bigwigs")
stanford_samples =  c(
  "B_cells_S", "B_cells_U",  "Naive_B-S",  "Naive_B-U",  "Mem_B-S",  "Mem_B-U",
  "CD4_Teff_S",  "CD4_Teff_U",  "Naive_Teffs-S",  "Naive_Teffs-U",  "Memory_Teffs-S",  "Memory_Teffs-U",  "Memory_Tregs-S",  "Memory_Tregs-U",
  "Th1_precursors_S",  "Th1_precursors_U",  "Th2_precursors_S",  "Th2_precursors_U",  "Th17_precursors_S",  "Th17_precursors_U",  "Tfh-S",  "Tfh-U",
  "CD8pos_S",  "CD8pos_U",  "Naive_CD8_T-S",  "Naive_CD8_T-U",  "CD8cm_S",  "CD8cm_U",  "CD8em_S",  "CD8em_U",
  "Regulatory_T-S",  "Regulatory_T-U",
  "Gamma_delta_T-S", "Gamma_delta_T-U",
  "Immature_NK-U", "Mature_NK-S",  "Mature_NK-U",  "Memory_NK-U",
  "Monocytes-S",  "Monocytes-U",  "Myeloid_DCs-U",  "pDCs")

region_broad = list(chr=finemap_reg$chr, start=min(finemap_reg[,"start"])-6000, end=max(finemap_reg[,"start"])-17500)

stanford_reg = extractAtacRegion(stanford_samples, stanford_dir, region_broad)
stanford_reg_sub = stanford_reg[stanford_reg$type %in% c("B_cells_S","B_cells_U","CD4_Teff_S","CD4_Teff_U","CD8pos_S","CD8pos_U"),]






snpTrack<-ggplot(data=finemap_reg, aes(x=start,xend=start, y=0,yend=1)) + geom_segment(size=0.1) +
facet_grid(tag~., margins=FALSE, switch="y") +
scale_x_continuous(limits = c(region_broad[["start"]],region_broad[["end"]]), expand = c(0, 0)) +
theme(strip.text.y = element_blank(),
      strip.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid=element_line(size=0.5),
      panel.background = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      legend.position = "none",
      axis.title.y=element_blank())


genet<-geneTrack(region_broad,gene_list=c("LINC02357","RBPJ"))


getAtacPlot = function(dat_reg, region, title=NULL) {
  B_types = c("B_cells_S", "B_cells_U",  "Naive_B-S",  "Naive_B-U",  "Mem_B-S",  "Mem_B-U")
  CD8_types = c("CD8pos_S",  "CD8pos_U",  "Naive_CD8_T-S",  "Naive_CD8_T-U",  "CD8cm_S",  "CD8cm_U",  "CD8em_S",  "CD8em_U")
  CD4_types = c("CD4_Teff_S",  "CD4_Teff_U",  "Naive_Teffs-S",  "Naive_Teffs-U",  "Memory_Teffs-S",  "Memory_Teffs-U","Th1_precursors_S",  "Th1_precursors_U",  "Th2_precursors_S",  "Th2_precursors_U",  "Th17_precursors_S",  "Th17_precursors_U",  "Tfh-S",  "Tfh-U")
  Regulatory_T_types = c("Memory_Tregs-S", "Memory_Tregs-U","Regulatory_T-S",  "Regulatory_T-U")
  GD_T_types = c("Gamma_delta_T-S","Gamma_delta_T-U")
  NK_types = c("Immature_NK-U", "Mature_NK-S",  "Mature_NK-U",  "Memory_NK-U")
  innate_types = c("Monocytes-S","Monocytes-U","Myeloid_DCs-U", "pDCs")
  atac_cell_types_ordered = c(B_types, CD8_types, CD4_types, Regulatory_T_types, GD_T_types, NK_types, innate_types)
  dat_reg$plotType = factor(dat_reg$type, levels=atac_cell_types_ordered)
  dat_reg$plotCat[dat_reg$type %in% B_types] <- "B"
  dat_reg$plotCat[dat_reg$type %in% CD8_types] <- "CD8"
  dat_reg$plotCat[dat_reg$type %in% CD4_types] <- "CD4"
  dat_reg$plotCat[dat_reg$type %in% Regulatory_T_types] <- "Regulatory"
  dat_reg$plotCat[dat_reg$type %in% GD_T_types] <- "Gamma-delta"
  dat_reg$plotCat[dat_reg$type %in% NK_types] <- "NK"
  dat_reg$plotCat[dat_reg$type %in% innate_types] <- "Innate"
  atac_cell_categories_ordered = c("B", "CD8",  "CD4",  "Regulatory", "Gamma-delta",  "NK", "Innate")
  atac_cell_categories_colors = c(brewer.pal(9, "Reds")[5],  #B
                                  brewer.pal(9, "Blues")[6],#CD8
                                  brewer.pal(9, "Blues")[8], #CD4
                                  brewer.pal(9, "Greys")[5], #Treg
                                  brewer.pal(9, "Greys")[5], #Gamma-delta
                                  brewer.pal(9, "Purples")[6],  #NK
                                  brewer.pal(9,"YlOrRd")[5]) #Innate
  dat_reg$plotCat = factor(dat_reg$plotCat, levels=atac_cell_categories_ordered)
  dat_reg$score[dat_reg$score>100] <- 100
  chrpos_ticks_breaks = seq(from = region[["start"]], to = region[["end"]], by = round((region[["end"]]-region[["start"]])/4, digits=0))
  chrpos_ticks_labels = paste0(round(chrpos_ticks_breaks/1e6, digits=3),"Mb")
  ggplot(data=dat_reg,aes(xmin=start, xmax=end, ymin=0, ymax=score, fill=plotCat)) + geom_rect() +
    scale_x_continuous(limits = c(region[["start"]],region[["end"]]), expand = c(0, 0), breaks=chrpos_ticks_breaks, labels=chrpos_ticks_labels) +
    #scale_y_continuous(limits = c(0, 1.25*max(dat_reg$score)), expand = c(0,0)) +
    scale_fill_manual(values=atac_cell_categories_colors) +
    facet_grid(plotType~., margins=FALSE, switch="y") +
    ggtitle(title) +
    theme(
      strip.text.y = element_blank(),
      strip.background = element_blank(),
      panel.spacing.y = unit(0.1,"lines"),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_text(size=4), 
      axis.text.y=element_blank(),axis.ticks.y=element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust=0.5, face="bold"))
}


plot1<-getAtacPlot(stanford_reg_sub, region_broad) 
library(ggbio)
plots<-tracks(snpTrack,genet, plot1,heights=c(0.4,1,8.6)) + geom_vline(xintercept=26083858, size=0.08, colour="lightblue") + geom_vline(xintercept=26083889, size=0.08, colour="lightblue") +
geom_vline(xintercept=26084947, size=0.08, colour="lightblue") + geom_vline(xintercept=26086506, size=0.08, colour="lightblue") + geom_vline(xintercept=26093692, size=0.08, colour="lightblue") + 
theme(axis.ticks.x=element_line(size=0.1))
ggsave(plots, file=paste0("~/output/",extension,"_gwas/paintor_5pc_diff/atac_bigwigs_1_up.png"),dpi=800, height=10, width=7.5,units="cm")


