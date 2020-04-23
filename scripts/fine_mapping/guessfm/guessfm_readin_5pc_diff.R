#guessfm_readin_5pc_diff.R

#This reads in the GUESS stochastic search and does the post-processing analysis to generate	'credible variants'.

library(ggplot2)
library(ggbio)
library(GUESSFM)
library(R2GUESS)
library(reshape2)
library(ggplot2)
library(snpStats)
library(speedglm)
library(parallel)
library(GenomicRanges)


args = commandArgs(trailingOnly=TRUE)

#define regions:
extension="final_collection"
orig<-read.table(file=paste0("/scratch/ji3x/",extension,"_gwas/pmeta_fam_ind_risk_regions_vcf_5pc_diff.txt"),header=T,sep=",",as.is=T)
orig<-orig[orig$ichip=="yes",]

#set paths
mydir <-paste0("/scratch/ji3x/",extension,"_gwas/guessfm_5pc_diff/input/")
outdir<-paste0("~/output/",extension,"_gwas/guessfm_5pc_diff/")

#load previous results which will contain meta data on each SNP:
load(file=paste0("/scratch/ji3x/",extension,"_gwas/resultsmeta_ind_fam_tdt_5pc_diff.RData"))
genes<-read.table(file="/scratch/ji3x/final_collection_gwas/rawdat/refFlat.txt.gz",as.is=T)
genes<-genes[!duplicated(genes$V1),]
genes<-GRanges(seqnames=genes$V3,
ranges=IRanges(genes$V5,end=genes$V6),
gene=genes$V1)


readitallin<-function(snp){
load(file=paste0(mydir, gsub(":",".",snp),"/data.RData"))

mydir<-paste0(mydir,gsub(":",".",snp),"/")
outdir<-paste0(outdir,gsub(":",".",snp),"/")
list.files(mydir)

#create outdir if it doesn't exist:
if(!dir.exists(paste0(outdir)))
dir.create(paste0(outdir),recursive=T)


## read output with R2GUESS and run a basic qc plot
if (!file.exists(paste0(outdir,"/model_diagnostics.png"))){
ess <- ess.read(f=paste0(mydir,"out_55000_55000"))
png(file=paste0(outdir,"/model_diagnostics.png"),
res=200, height=20, width=20, units="cm")
k<-R2GUESS:::plot.ESS(ess)
dev.off()

png(file=paste0(outdir,"/model_convergence.png"),
res=200, height=20, width=20, units="cm")
par(mfrow=c(1,2))
k1<-R2GUESS::check.convergence(ess)
dev.off()
}

## now read output with GUESSFM
if(file.exists(paste0(mydir,"out_55000_55000_sweeps_output_best_visited_models.txt"))){
system(paste0("mv ", mydir, "out_55000_55000_sweeps_output_best_visited_models.txt ",mydir,"out_55000_55000_output_best_visited_models.txt"))
system(paste0("mv ",mydir,"out_55000_55000_sweeps_output_marg_prob_incl.txt ", mydir,"out_55000_55000_output_marg_prob_incl.txt"))
}
dd <- read.snpmod(mydir)

#Examine posterior of number of SNPs in model:

png(file=paste0(outdir,"/model_size.png"),
res=200, height=20, width=20, units="cm")
pp <- pp.nsnp(dd,plot=TRUE,expected=3, overdispersion = 1.00000001)
dev.off()

## examine the best models and SNPs with greatest marginal support within the tagged data.
best.models(dd)
best.snps(dd)


#Now expand these tag SNPs:
load(paste0(mydir,"/tags.RData"))
dx <- expand.tags(dd, tags)

best<-best.models(dx, cpp.thr = 0.99)

save(best, dx, file=paste0(mydir,"/best_models.RData"))
load(file=paste0(mydir,"/best_models.RData"))



colnames(DATA)<-gsub(":",".", colnames(DATA))
Y<-as.matrix(Y)

#taking residuals from logistic regression including covriates as outcome... it's an approximation
abf <- abf.calc(y=Y,x=DATA,q=covariates,models=best$str,family="binomial", verbose=TRUE, approx.lm=T)
abf$nsnp=ncol(DATA)

sm.all <- abf2snpmod(abf,expected=3, overdispersion = 1.00000001)
sp.all <-snp.picker(d=sm.all, data=DATA)

save(sm.all,sp.all,file=paste0(mydir,"/smsp_residuals.RData"))

#Grouping:
load(file=paste0(mydir,"/smsp_residuals.RData"))

groups <- as(sp.all,"groups")
l<-res[gsub(":",".",res$MarkerName) %in% colnames(DATA),]
rownames(l)<-gsub(":",".",l$MarkerName)
if(length(groups)>0){
summx <- guess.summ(sm.all,groups=groups,snps=l,position="position")
summx <- scalepos(summx,position="position")

save(sm.all,sp.all,groups,l,summx, file=paste0(mydir,"summx.RData"))
load(file=paste0(mydir,"summx.RData"))

#get the variant names where available:
#create regions file for BCF tools to extract SNPs needed:
o<-summx[,c("chromosome","position")]
write.table(o, file=paste0("/scratch/ji3x/full_cohort_gwas/vars/sig_",gsub(":",".",snp),".txt"), sep="\t",
col.names=F, row.names=F, quote=F)

system(paste0("/home/ji3x/software/bcftools/bcftools view -R /scratch/ji3x/full_cohort",
"_gwas/vars/sig_",gsub(":",".",snp),".txt /scratch/ji3x/full_cohort",
"_gwas/vars/Kaviar-160204-Public-hg38-trim.vcf.gz > /scratch/ji3x/full_cohort",
"_gwas/vars/sig_",gsub(":",".",snp),"_ids.txt"))

vars<-read.table(file=paste0("/scratch/ji3x/full_cohort",
"_gwas/vars/sig_",gsub(":",".",snp),"_ids.txt"),skip=40, header=T, comment.char="", as.is=T,colClasses = c("character"))

vars$MarkerName<-paste0("chr",vars$X.CHROM,":",vars$POS,":",vars$REF,":",vars$ALT)
vars<-vars[,c("MarkerName","ID")]
vars<-vars[vars$ID!=".",]
summx<-merge(summx,vars,by="MarkerName", all.x=T)
summx<-summx[order(summx$chromosome, summx$position),]
rownames(summx)<-summx$snp
summx$ID<-ifelse(summx$snp=="chr2.162279995.C.G","rs35337543", 
ifelse(summx$snp=="chr12.56042145.C.G","rs1131017",
ifelse(summx$snp=="chr14.68796357.G.A","rs6573857",
ifelse(summx$snp=="chr14.98011018.A.T","rs79835824",
ifelse(summx$snp=="chr16.11190252.A.G","rs248832",
ifelse(summx$snp=="chr16.28494339.G.C","rs151234",
ifelse(summx$snp=="chr21.44280294.T.A","rs9671675",
ifelse(summx$snp=="chr21.42418534.G.A","rs13048049",summx$ID))))))))
summx$ID<-ifelse(is.na(summx$ID),".",summx$ID)
write.table(summx, file=paste0(outdir,"credible_snps.txt"), quote=F, row.names=F, col.names=T, sep="\t")


mod2group <- function(str,groups) {  
  x.snps <- strsplit(str,"%")
  x.groups <- lapply(x.snps,snpin,groups)
  G <- sapply(x.groups,function(g) {
    if(is.null(g))
      return("N")
    match <- apply(g,1,any)
    ret <- numeric(nrow(g))
    if(any(match))
      ret[match] <- apply(g[match,,drop=FALSE],1,which)
    return(paste(sort(ret),collapse="-"))
  })
  G[str=="1"] <- "-" # distinguish null model from ungrouped snps
  return(G)
}
numorneg <- function(x) {
  suppressWarnings(n <- as(x,"numeric"))
  if(any(is.na(n)))
    n[is.na(n)] <- -1
  return(n)
}
gsumm <- function(x,groups) {
  n <- max(numorneg(unlist(strsplit(x$group,"-"))))
  counts <- tapply(x$PP,x$group,sum)
  df <- data.frame(pattern=names(counts), PP=counts)
  df <- df[order(df$PP,decreasing=TRUE), ]
  df$ymax <- cumsum(df$PP)
  df$ymin <- c(0,df$ymax[-nrow(df)])
  cn <- lapply(strsplit(rownames(df),"-"),numorneg)  
  df2 <- df[rep(1:nrow(df), times=sapply(cn,length)),]
  df2$xmin <- unlist(cn)
  df2$xmax <- df2$xmin+1
  return(df2)
}


pattern.plot1<-function (SM, groups) 
{
    if (!is.list(SM)) 
        SM <- list(trait = SM)
    BM <- lapply(SM, function(x) best.models(x, cpp.thr = 0.99)[, 
        c("str", "PP")])
    for (i in seq_along(BM)) BM[[i]]$group <- mod2group(BM[[i]]$str, 
        groups)
    G <- lapply(BM, gsumm)
    if (is.null(names(G))) 
        names(G) <- paste0("trait", seq_along(G))
    for (i in names(G)) G[[i]]$trait <- i
    G <- do.call("rbind", G)
    G$xmin <- G$xmin - 1
    G$xmax <- G$xmax - 1
    labs<-summx[summx$snp %in% groups@tags,]
    labs<-labs[groups@tags,]
    la<-labs$ID
    do<-grepl("^-",G$pattern)
    cols=c("grey","black","red","blue","green","purple","yellow","brown","pink","turquoise","khaki1","slateblue1","deepskyblue1","linen","firebrick1")
    cols<-cols[c(1:(length(la)+2))]
    names(cols)<-c("None","Other",la)
    labs<-c("None","Other",la)
    if(!TRUE %in% do){
    cols<-cols[2:length(cols)]
    labs<-labs[2:length(labs)]
    }
    p <- ggplot(G, aes(xmin = xmin, xmax = xmax, ymin = ymin, 
        ymax = ymax)) + geom_rect() + geom_vline(aes(xintercept = xmin), 
        col = "gray60", size = 0.2) + geom_hline(aes(yintercept = ymax), 
        col = "gray60", size = 0.2) + scale_x_continuous(breaks = sort(unique(G$xmin)) + 
        0.5, labels = labs, limits = c(min(G$xmin), 
        max(G$xmin) + 1), expand = c(0, 0)) + scale_y_continuous(breaks = c(0, 
        1), expand = c(0, 0), limits = c(0, 1)) + theme_bw() + 
        theme(panel.grid = element_blank(), legend.position = "none", 
            panel.spacing = unit(1, "lines"),
            axis.text.x = element_text(colour=cols[labs],angle = 45,
            vjust = 1, hjust = 1)) + xlab("SNP group index") + ylab("Cumulative Model Posterior Probility") + 
        scale_fill_manual(values = c(`FALSE` = "grey", `TRUE` = "grey20"))
        return(p)
}

pat<-pattern.plot1(sm.all,groups)
ggsave(pat,file=paste0(outdir,"/cum_prob_plot_upd.png"),
dpi=200, height=20, width=20, units="cm")
ggsave(pat,file=paste0(outdir,"/cum_prob_plot_upd.pdf"))

signal.plot1<-function (summ, groups, w = 0.2, highlight = NULL) 
{
    if (!("x.scale" %in% colnames(summ)) || any(is.na(summ$x.scale))) {
        stop("Missing x co-ordinates.  Do you need to run scalepos or fix some missing values?")
    }
    summ$xmin <- summ$xmin.scale - w
    summ$xmax <- summ$xmax.scale + w
    summ$x <- summ$x.scale
    summ$tag <- as.factor(summ$tag)
    ord<-data.frame(tag=groups@tags,ord=c(1:length(groups@tags)))
    summ<-merge(summ,ord,by="tag")
    summ$tag<-reorder(summ$tag,summ$ord)
    cols=c("red","blue","green","purple","yellow","brown","pink","turquoise","khaki1","slateblue1","deepskyblue1","linen","firebrick1")
    cols<-cols[c(1:nrow(ord))]
    names(cols)<-ord$tag
    summ$ID<-ifelse(summ$ID==".",summ$snp,summ$ID)
    p <- ggplot(summ, aes(x = x, col = tag)) + geom_vline(aes(xintercept = x, 
        col = tag), alpha = 0.5, size = 0.2) + geom_hline(yintercept = 1, 
        col = "grey", size = 0.2, linetype = "dashed") + geom_point(aes(y = pp), 
        size = 1) + ylab("PP") + theme(legend.position = "none", 
        strip.text.y = element_text(size = 9, angle = 0), panel.grid = element_blank()) + 
        geom_segment(data = subset(summ, !is.na(xmin)), mapping = aes(x = xmin, 
            xend = xmax, y = ppsum, yend = ppsum, col = tag), 
            size = 0.5) + geom_rect(data = subset(summ, !is.na(x.min)), 
        mapping = aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = ppsum, 
            fill = tag), alpha = 0.1, size = 0.2) +
	scale_fill_manual(values=cols) +
	scale_colour_manual(values=cols) +
	annotate("text",x=summ$x, y=rep(c(0.2,0.8),times=1000)[c(1:nrow(summx))] , label=summ$ID, colour="black", angle=90,size=2)
        return(p)
}

signals<-signal.plot1(summx,groups)

ggchr1<-function (summ, groups) 
{
    ord<-data.frame(tag=groups@tags,ord=c(1:length(groups@tags)))
    summ<-merge(summ,ord,by="tag")
    summ$tag<-reorder(summ$tag,summ$ord)
    cols=c("red","blue","green","purple","yellow","brown","pink","turquoise","khaki1","slateblue1","deepskyblue1","linen","firebrick1")
    cols<-cols[c(1:nrow(ord))]
    names(cols)<-ord$tag
    summ$ID<-ifelse(summ$ID==".",summ$snp,summ$ID)

    ggplot(summ, aes(x = x.scale, xend = position.plot, y = 0, 
        yend = 1, colour = tag)) + ggplot2::geom_segment() + 
	scale_colour_manual(values=cols) +
        theme(legend.position = "none", strip.text.y = element_text(size = 9, 
            angle = 0), axis.title.y = element_blank(), axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(), panel.grid = element_blank())
}


chr<-ggchr1(summx,groups)
if(nrow(summx)>1)
lds<-ggld(DATA, summx)


if (nrow(summx)>1){
out<-tracks(chr,signals,lds,heights=c(1,2,1))
}
if(nrow(summx)==1){
out<-tracks(chr,signals,heights=c(1,2))
}
ggsave(out,file=paste0(outdir,"/init_out_upd.png"),
dpi=200, height=25, width=20, units="cm")
chromosome<-sub("(^.*)[:](.*)[:](.*)[:](.*)","\\1",snp)
g<-GRanges(seqnames=c(paste0(chromosome)),
ranges=IRanges(min(summx$position),end=max(summx$position)))
o1<-subsetByOverlaps(genes,g)
o1<-as.data.frame(o1)
if(nrow(o1)>0){
m<-(0.2*nrow(o1))
o1$y=c(seq(0.2,m,0.2))
}
if(nrow(o1)==0){
g<-GRanges(seqnames=c(paste0("chr",chromosome)),
ranges=IRanges(min(summx$position)-200000,end=max(summx$position)+200000))
o1<-subsetByOverlaps(genes,g)
o1<-as.data.frame(o1)
if(nrow(o1)>0){
m<-(0.2*nrow(o1))
o1$y=c(seq(0.2,m,0.2))
}
}
if(nrow(o1)>0){
huesg = c("green","purple","turquoise1","orange","slategrey","tan","pink","lightblue")
huesg =huesg[1:nrow(o1)]
o1$gene<-as.factor(o1$gene)
names(huesg)<-levels(o1$gene)


t1<-ggplot(data=o1, aes(x=start, xend=end, y=y, yend=y,colour=as.factor(gene))) + geom_segment(arrow=arrow(length = unit(0.03, "npc"))) +
theme(axis.text.y=element_blank(),
axis.ticks.y=element_blank(),
axis.title.y=element_blank(),
axis.text.x=element_blank(),
axis.title.x=element_blank(),
axis.ticks.x=element_blank()) +
scale_colour_manual(name="Gene name",values=c(huesg)) + scale_y_continuous(limits=c(0, (m+0.2)))
if (nrow(summx)>1){
out1<-tracks(t1,chr,signals,lds,heights=c(2,1,2,1))
}
if(nrow(o1)==0){
out1<-tracks(chr,signals,lds,heights=c(2,1,2))
}
ggsave(out1,file=paste0(outdir,"/init_out_genes_upd.png"),
dpi=200, height=25, width=20, units="cm")
ggsave(out1, file=paste0(outdir,"/init_out_genes_upd.pdf"))
}

}
if(length(groups)==0){
}
message(paste0("Done ",snp))
}

readitallin(orig$Marker[as.numeric(args[1])])



