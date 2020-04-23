options(stringsAsFactors=FALSE)
library(snpStats)
library(qqman)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
group = args[1]
dir=args[2]

cat("START ",group,"\n")

setwd(dir)

### Get genotyped results
cat("Getting genotyped tdt results...\n")
geno1=read.table(paste("genotyped.tdt",sep=""), header=TRUE)
geno2=read.table(paste("genotyped_WRONG.tdt",sep=""),header=TRUE)

### Get imputed results
cat("Getting TOPMed imputed tdt results...\n")
topmedList1 = list()
topmedList2 = list()
for (i in 1:22) {
	filename1=paste("chr",i,".filter_maf_gt_0.005_rsq_gt_0.8.relateds_with_pheno.tdt",sep="")
	filename2=paste("chr",i,".filter_maf_gt_0.005_rsq_gt_0.8.relateds_WRONG.tdt",sep="")
	topmedList1[[i]] = read.table(filename1, header=TRUE)
	topmedList2[[i]] = read.table(filename2, header=TRUE)
}
topmed1 = do.call("rbind", topmedList1)
topmed2 = do.call("rbind", topmedList2)

#write.table(topmed1, file="chrALL.filter_maf_gt_0.005_rsq_gt_0.8.relateds_with_pheno.tdt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

#Filter TOPMed and 1000G for Mendelian Inconsistencies
droplistMI = scan("snps_with_mendelInconsistencies_gt5pct.txt", what="character")
#droplistMI = scan("snps_with_mendelInconsistencies_gt1pct.txt", what="character")
topmed1_filtered= topmed1[!topmed1$SNP%in%droplistMI,]
topmed2_filtered= topmed2[!topmed2$SNP%in%droplistMI,]

cat("Creating plots...\n")

### QQ plots
jpeg(paste("tdt_results_",group,"_1.jpg",sep=""))
par(mfrow=c(3,2))
qq(geno1$P, main="Genotyped (Affected trios)")
qq(geno2$P, main="Genotyped (Unaffected trios)")
qq(topmed1$P, main="TOPMed Imputed (Affected trios)")
qq(topmed2$P, main="TOPMed Imputed (Unaffected trios)")
qq(topmed1_filtered$P, main="TOPMed Imputed Filtered (Affected trios)")
qq(topmed2_filtered$P, main="TOPMed Imputed Filtered (Unaffected trios)")
dev.off()


### Manhattan plots
jpeg(paste("tdt_results_",group,"_2.jpg",sep=""))
par(mfrow=c(3,1))
manhattan(geno1, ylim=c(0,20), main="Genotyped (Affected trios)")
manhattan(thousg1, ylim=c(0,20), main="1000G Imputed (Affected trios)")
manhattan(topmed1, ylim=c(0,20), main="TOPMed Imputed (Affected trios)")
dev.off()

jpeg(paste("tdt_results_",group,"_3.jpg",sep=""))
par(mfrow=c(3,1))
manhattan(geno2, ylim=c(0,20), main="Genotyped (Unaffected trios)")
manhattan(thousg2, ylim=c(0,20), main="1000G Imputed (Unaffected trios)")
manhattan(topmed2, ylim=c(0,20), main="TOPMed Imputed (Unaffected trios)")
dev.off()


### Direction of effect in Affected vs Unaffected trios
jpeg(paste("tdt_results_",group,"_4.jpg",sep=""))
par(mfrow=c(2,2))
plot(log10(geno1$OR),log10(geno2$OR), xlab="Affected trios", ylab="Unaffected trios", main="Genotyped")
plot(log10(thousg1$OR),log10(thousg2$OR), xlab="Affected trios", ylab="Unaffected trios", main="1000G Imputed")
plot(log10(topmed1$OR),log10(topmed2$OR), xlab="Affected trios", ylab="Unaffected trios", main="TOPMed Imputed")
dev.off()

### Manhattan plots after filtering
jpeg(paste("tdt_results_",group,"_5.jpg",sep=""))
par(mfrow=c(2,1))
manhattan(topmed1_filtered, ylim=c(0,20), main="TOPMed Imputed (Affected trios)\nFiltered for MI")
manhattan(topmed2_filtered, ylim=c(0,20), main="TOPMed Imputed (Unaffected trios)\nFiltered for MI")
dev.off()

jpeg(paste("tdt_results_",group,"_6.jpg",sep=""))
par(mfrow=c(2,1))
manhattan(thousg1_filtered, ylim=c(0,20), main="1000G Imputed (Affected trios)\nFiltered for MI")
manhattan(thousg2_filtered, ylim=c(0,20), main="1000G Imputed (Unaffected trios)\nFiltered for MI")
dev.off()


### QQ plots after filtering
jpeg(paste("tdt_results_",group,"_1B.jpg",sep=""))
par(mfrow=c(2,2))
qq(thousg1_filtered$P, main="1000G Imputed (Affected trios)\nFiltered for MI")
qq(thousg2_filtered$P, main="1000G Imputed (Unaffected trios)\nFiltered for MI")
qq(topmed1_filtered$P, main="TOPMed Imputed (Affected trios)\nFiltered for MI")
qq(topmed2_filtered$P, main="TOPMed Imputed (Unaffected trios)\nFiltered for MI")
dev.off()


### PLOT SNP QC INFO
mi=read.table("king_snpqcbySNP.txt", header=TRUE)
info=read.table("../../results_filtered/infofiles/FIN_chrAll.maf_gt_0.005.info.gz", header=TRUE)
minfo = merge(mi, info, by="SNP")

mi_tg = read.table("king_snpqc_thousgbySNP.txt", header=TRUE)
info_tg = readRDS("../../../IMPUTED_1KG/results/vcffiles/mymega_chrAll.filtered.info.rds")
info_tg$SNP2 = paste(info_tg$SNP,info_tg$REF.0., info_tg$ALT.1., sep=":")
minfo_tg = merge(mi_tg, info_tg, by.x="SNP", by.y="SNP2")

jpeg(paste("tdt_results_",group,"_7.jpg",sep=""))
par(mfrow=c(4,2))
plot(minfo$Rsq, minfo$Err_InPO, xlab="R-squared", ylab="Error Rate in PO pairs", main="TOPMed Imputed", ylim=c(0,0.08))
plot(minfo_tg$Rsq, minfo_tg$Err_InPO, xlab="R-squared", ylab="Error Rate in PO pairs", main="1000G Imputed", ylim=c(0,0.08))
plot(minfo$MAF, minfo$Err_InPO, xlab="MAF", ylab="Error Rate in PO pairs", main="TOPMed Imputed", ylim=c(0,0.08))
plot(minfo_tg$MAF, minfo_tg$Err_InPO, xlab="MAF", ylab="Error Rate in PO pairs", main="1000G Imputed", ylim=c(0,0.08))
hist(minfo$Rsq[minfo$Err_InPO==0],xlab="R-squared", main="TOPMed Imputed (Error Rate == 0)")
hist(minfo_tg$Rsq[minfo$Err_InPO==0],xlab="R-squared", main="1000G Imputed (Error Rate == 0)")
hist(minfo$Rsq[minfo$Err_InPO>0],xlab="R-squared", main="TOPMed Imputed (Error Rate > 0)")
hist(minfo_tg$Rsq[minfo$Err_InPO>0],xlab="R-squared", main="1000G Imputed (Error Rate > 0)")
dev.off()

cat("FINISHED ",group,"\n")



##2D density plots in ggplot --> https://www.r-graph-gallery.com/2d-density-plot-with-ggplot2/
#ggplot(minfo, aes(x=Rsq, y=Err_InPO) ) +
#		geom_bin2d(bins=70) +
#		theme_bw()
#
#ggplot(minfo[minfo$Err_InPO>0,], aes(x=Rsq, y=Err_InPO) ) +
#		geom_hex(bins=70) +
#		theme_bw()


#### TDT
#args = commandArgs(trailingOnly=TRUE)
#filename = args[1]
#d = read.table(filename, header=TRUE)
##d = read.table("relateds_snpqc.tdt", header=TRUE)
##d = read.table("mega_b37_TRIOS_snpqc.tdt", header=TRUE)
#
###Get lead snp on chromosome 8
##d[d$CHR==8 & d$P==min(d$P[d$CHR==8]),]
##d[d$CHR==8 & d$P<5e-5,]
#
#d$P[d$P<1e-50] <- 1e-50
#
#getQCplots = function(d, title) {
#	manhattan(d)
#	CHISQ <- qchisq(1-d$P,1)
#	gif = median(CHISQ)/qchisq(0.5,1)
#	qq(d$P, main = paste(title,"\n GIF=", round(gif, digits=2)))
#
#	#gif_1k = 1 + (gif-1)*((1/ncase + 1/ncontrol)/(1/1000 + 1/1000))
#	#qq(d$P, main = paste(title,"\n GIF=", round(gif, digits=2), "\n GIF_1000=",round(gif_1k, digits=2)))
#}
#
#
#pdf(file=paste("plots/",filename,".pdf", sep=""))
#	getQCplots(d, title="Basic TDT Test P-values")
#	d_noMHC = d[!((d$CHR==6 & d$BP>20000000 & d$BP<40000000)|(d$CHR==11 & d$BP>1000000 & d$BP<3000000)|(d$CHR==1 & d$BP>112000000 & d$BP<115000000)),]
#	getQCplots(d_noMHC, title="Basic TDT Test P-values (removed MHC, INS, and PTPN22 regions)")
#dev.off()




### GDT
#gd = read.table("gdt.tbl", header=TRUE)
#gd$PVALUE[gd$PVALUE<1e-20] <- 1e-20
#manhattan(gd, chr="CHR", p="PVALUE", bp="POS")
#
#CHISQ <- qchisq(1-gd$PVALUE,1)
#gif = median(CHISQ)/qchisq(0.5,1)


#group="FIN"
#plinkresults1=paste("/nv/vol185/MEGA/IMPUTED_TOPMED/family_analysis/",group,"/mega_b37_TRIOS_",group,".tdt", sep="")
#plinkresults2=paste("/nv/vol185/MEGA/IMPUTED_TOPMED/family_analysis/",group,"/mega_b37_TRIOS_",group,"_WRONG.tdt", sep="")
#
#dplink1=read.table(plinkresults1, header=TRUE)
#dplink2=read.table(plinkresults2, header=TRUE)
#sum(!dplink1$SNP==dplink2$SNP)
#
#head(dplink1); head(dplink2)
#
#plot(dplink1$OR,dplink2$OR)
#plot(log(dplink1$OR),log(dplink2$OR),xlab="true trios", ylab="unaffected trios")
#plot(log(dplink1$OR[dplink1$P<1e-5]),log(dplink2$OR[dplink1$P<1e-5]))
#
#qq(dplink1$P)
#qq(dplink2$P)
