options(stringsAsFactors=FALSE)
setwd(Sys.getenv("impbench"))
library(vcfR)
library(ggplot2)
library(gridExtra)
library(MASS)
library(dplyr)
library(GenomicRanges)

args = commandArgs(trailingOnly=TRUE)
subset = args[1]
refpanel = args[2]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		Read in accuracy tables
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
eur.df = read.table(paste0("imputation_vs_wgs_",refpanel,"_EUR_",subset,"_regions.txt"), header=TRUE)
amr.df = read.table(paste0("imputation_vs_wgs_",refpanel,"_AMR_",subset,"_regions.txt"), header=TRUE)
afr.df = read.table(paste0("imputation_vs_wgs_",refpanel,"_AFR_",subset,"_regions.txt"), header=TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		Plot relationship between true R2 and estimated R2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plotTrueVsEst = function(df, maf_min, title) {
	df_keep = df[df$wgs==1 & df$mega==1 & df$wgs_maf>maf_min & df$mega_maf>0,]
	true_thresh = 0.5
	smoothScatter(df_keep[,c("imp_Rsq2", "genomic_Rsq")], xlab="Imputation R-squared", ylab="True R-squared", main=title)
	abline(v=0.8)
	segments(x0=0.8, x1=1.2, y0=true_thresh, y1=true_thresh, lty=2)
	legend(0.85,0.97, paste(round(mean(df_keep$genomic_Rsq[df_keep$imp_Rsq2>0.8]>true_thresh),4)*100, "%", sep=""), bty="n")
}

pdf(paste0("true_vs_estimated_rsq_",refpanel,".pdf"))
par(mfrow=c(2,2))
plotTrueVsEst(afr.df, 0.005, "African-American (N=1,411)")
plotTrueVsEst(amr.df, 0.005, "Admixed (N=641)")
plotTrueVsEst(eur.df, 0.005, "European (N=95)")
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		Plot distribution of trueRsq among variants with impRsq>0.8
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plotDistributions = function(df, maf_min, Ersq_thresh, title) {
	df_keep = df[df$mega==1 & df$mega_maf>maf_min & df$imp_Rsq2>Ersq_thresh,]
	ggplot(df_keep) + geom_histogram(aes(x=genomic_Rsq)) + xlab("True R-squared") + ggtitle(title)
}
pdf(paste0("genomicRsq_dist_in_snps_to_keep_",refpanel,".pdf"))
plotDistributions(eur.df, 0.005, 0.8, "EUR")
plotDistributions(afr.df, 0.005, 0.8, "AFR")
plotDistributions(amr.df, 0.005, 0.8, "AMR")
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		Calculate coverage by R2 threshold
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getCoverage = function(df, maf_min, maf_max=1) {
	df_common = df[df$wgs==1 & df$wgs_maf>maf_min & df$wgs_maf<maf_max,]
	thresholds = seq(0.5,1, by=0.005)
	coverage_Trsq = NULL
	coverage_Ersq = NULL
	for (i in 1:length(thresholds)) {
		t = thresholds[i]
		coverage_Trsq[i] = sum(df_common$genomic_Rsq>t, na.rm=TRUE)/sum(df_common$wgs==1)
		coverage_Ersq[i] = sum(df_common$imp_Rsq2>t, na.rm=TRUE)/sum(df_common$wgs==1)
	}
	summDF = data.frame(thresholds, coverage_Trsq, coverage_Ersq)
	summDF
}

getCoverageTable = function(maf_min, maf_max) {
	summDF_all = rbind(data.frame(group="AFR",getCoverage(afr.df, maf_min, maf_max)), data.frame(group="AMR",getCoverage(amr.df, maf_min, maf_max)))
	write.table(summDF_all, file=paste0("coverage_by_imprsq_threshold.maf_",maf_min,"-",maf_max,"_",refpanel,".txt"), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
}

getCoverageTable(maf_min=0, maf_max=1)
getCoverageTable(maf_min=0.005, maf_max=1)
getCoverageTable(maf_min=0, maf_max=0.005)
getCoverageTable(maf_min=0.005, maf_max=0.01)
getCoverageTable(maf_min=0.01, maf_max=0.05)
getCoverageTable(maf_min=0.05, maf_max=1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate coverage by region
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getCoveragebyRegion = function(df, Rsq_thresh, maf_min) {
	df$keep = NA
	df$keep[df$mega==1 & df$imp_Rsq2>Rsq_thresh & df$mega_maf>maf_min] <- 1
	df$keep[df$mega==1 & (df$imp_Rsq2<=Rsq_thresh | df$mega_maf<=maf_min)] <- 0
	df$keep[df$mega==0] <- 0
	df$TYPED[df$mega==0] <- 0
	summary(df$keep)
	df.wgs = df[df$wgs==1,]
	df.wgs.common = df.wgs[df.wgs$wgs_maf>maf_min,]
	summary(df.wgs.common$keep)
	coverage_by_region = df.wgs.common %>% group_by(region) %>% summarize(prop_kept=mean(keep), prop_typed=mean(TYPED), total_kept=sum(keep, na.rm=TRUE), total_typed=sum(TYPED), total_wgs=length(keep))
	as.data.frame(coverage_by_region)
}
getCoveragebyRegion(afr.df, Rsq_thresh=0.8, maf_min=0.005)
getCoveragebyRegion(amr.df, Rsq_thresh=0.8, maf_min=0.005)
getCoveragebyRegion(eur.df, Rsq_thresh=0.8, maf_min=0.005)

regionSumm = rbind(
		data.frame(group="AFR", getCoveragebyRegion(afr.df, Rsq_thresh=0.8, maf_min=0.005)),
		data.frame(group="AMR", getCoveragebyRegion(amr.df, Rsq_thresh=0.8, maf_min=0.005)),
		data.frame(group="EUR", getCoveragebyRegion(eur.df, Rsq_thresh=0.8, maf_min=0.005)))

write.table(regionSumm, file=paste0("coverage_by_region_rsq",refpanel,".txt"), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
