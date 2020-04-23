setwd("/nv/vol185/MEGA/IMPUTED_TOPMED_HWE/coverage")
library(ggplot2)
library(gridExtra)
library(ggExtra)
library(dplyr)
library(ggpubr)
library(rtf)

d = read.table("mendel_error_vs_imp_rsq.txt", header=TRUE)
d$cc_exclude = d$Err_InHetTrio>0.05 | d$Err_InHomPO>0.05
d$fam_exclude =d$Err_InHetTrio>0.01 | d$Err_InHomPO>0.01

getPercent = function(x) {
	round(mean(x)*100, digits=2)
}

sumTableTotalExcluded <- d %>% group_by(MAF_bin) %>%
		summarise(
				PropExcludeCC = getPercent(cc_exclude),
				PropExcludeFam = getPercent(fam_exclude),
				NumExcludeCC = sum(cc_exclude),
				NumExcludeFam = sum(fam_exclude)
				)
sumTableTotalExcluded$CC_sum = paste0(sumTableTotalExcluded$NumExcludeCC, " (", sumTableTotalExcluded$PropExcludeCC,"%)")
sumTableTotalExcluded$Fam_sum = paste0(sumTableTotalExcluded$NumExcludeFam, " (", sumTableTotalExcluded$PropExcludeFam,"%)")
		
sumTableByCategory <- d %>% group_by(MAF_bin) %>%
		summarise(
				PropHetTrio_gt_5 = getPercent(Err_InHetTrio>0.05),
				PropHetTrio_gt_1 = getPercent(Err_InHetTrio>0.01),
				PropHomPO_gt_5 = getPercent(Err_InHomPO>0.05),
				PropHomPO_gt_1 = getPercent(Err_InHomPO>0.01),
				NumHetTrio_gt_5 = sum(Err_InHetTrio>0.05),
				NumHetTrio_gt_1 = sum(Err_InHetTrio>0.01),
				NumHomPO_gt_5 = sum(Err_InHomPO>0.05),
				NumHomPO_gt_1 = sum(Err_InHomPO>0.01)
		)
sumTableByCategory$HTgt5 = paste0(sumTableByCategory$NumHetTrio_gt_5, " (", sumTableByCategory$PropHetTrio_gt_5,"%)")
sumTableByCategory$HTgt1 = paste0(sumTableByCategory$NumHetTrio_gt_1, " (", sumTableByCategory$PropHetTrio_gt_1,"%)")
sumTableByCategory$HPgt5 = paste0(sumTableByCategory$NumHomPO_gt_5, " (", sumTableByCategory$PropHomPO_gt_5,"%)")
sumTableByCategory$HPgt1 = paste0(sumTableByCategory$NumHomPO_gt_1, " (", sumTableByCategory$PropHomPO_gt_1,"%)")


rtffile <- RTF("mendel_error_vs_maf_bin.doc")
addParagraph(rtffile, "Supplementary Table ## - Variants with imputation R-squared>0.8 that were subsequently excluded for Mendel Errors")
addTable(rtffile, sumTableTotalExcluded[c("MAF_bin","CC_sum","Fam_sum")])
addParagraph(rtffile, "")

addParagraph(rtffile, "Supplementary Table ## - Variants with imputation R-squared>0.8 that were subsequently excluded for Mendel Errors")
addTable(rtffile, sumTableByCategory[,c("MAF_bin", "HTgt5","HTgt1", "HPgt5","HPgt1")])
done(rtffile)
		
		
		
#PLOTS SHOWING MENDEL ERROR RATES VS IMPUTATION RSQ 
png("mendel_error_vs_imp_rsq_CCTHRESH.png")
p1 = ggplot(d) + geom_point(aes(x=rsq, y=Err_InHetTrio*100), alpha=0.5, size=0.5) +
		ylab("Heterozygous Trio Error Rate (%)") +
		xlab("Imputation R-squared") +
		facet_grid(.~MAF_bin) + theme_classic() +
		geom_hline(yintercept=5, colour="red")  +
		ylim(0,45)
p2 = ggplot(d) + geom_point(aes(x=rsq, y=Err_InHomPO*100), alpha=0.5, size=0.5) +
		ylab("Homozygous PO Error Rate (%)") +
		xlab("Imputation R-squared") +
		facet_grid(.~MAF_bin) + theme_classic() +
		geom_hline(yintercept=5, colour="red") +
		ylim(0,45)
grid.arrange(p1, p2, nrow=2)
dev.off()


png("mendel_error_vs_imp_rsq_FAMTHRESH.png")
p1 = ggplot(d) + geom_point(aes(x=rsq, y=Err_InHetTrio*100), alpha=0.5, size=0.5) +
		ylab("Heterozygous Trio Error Rate (%)") +
		xlab("Imputation R-squared") +
		facet_grid(.~MAF_bin) + theme_classic() +
		geom_hline(yintercept=1, colour="red") +
		ylim(0,45)
p2 = ggplot(d) + geom_point(aes(x=rsq, y=Err_InHomPO*100), alpha=0.5, size=0.5) +
		ylab("Homozygous PO Error Rate (%)") +
		xlab("Imputation R-squared") +
		facet_grid(.~MAF_bin) + theme_classic() +
		geom_hline(yintercept=1, colour="red") +
		ylim(0,45)
grid.arrange(p1, p2, nrow=2)
dev.off()



### PLOTS WITH MARGINAL HISTOGRAMS
getPlot1 = function(d, thresh) {
	p = ggplot(d) + geom_point(aes(x=rsq, y=Err_InHetTrio*100), alpha=0.5, size=0.5) +
			ylab("Heterozygous Trio Error Rate (%)") +
			xlab("Imputation R-squared") +
			theme_classic() +
			geom_hline(yintercept=thresh, colour="red") +
			ylim(0,50)
	return(p)
}
getPlot2 = function(d, thresh) {
	p = ggplot(d) + geom_point(aes(x=rsq, y=Err_InHomPO*100), alpha=0.5, size=0.5) +
			ylab("Homozygous PO Error Rate (%)") +
			xlab("Imputation R-squared") +
			theme_classic() +
			geom_hline(yintercept=thresh, colour="red") +
			ylim(0,50)
	return(p)
}
getPlotGrid = function(thresh) {
	p1 = ggMarginal(getPlot1(d[d$MAF_bin=="MAF 0.5-1%",], thresh), type="histogram")
	p2 = ggMarginal(getPlot1(d[d$MAF_bin=="MAF 1-5%",], thresh), type="histogram")
	p3 = ggMarginal(getPlot1(d[d$MAF_bin=="MAF >5%",], thresh), type="histogram")
	p4 = ggMarginal(getPlot2(d[d$MAF_bin=="MAF 0.5-1%",], thresh), type="histogram")
	p5 = ggMarginal(getPlot2(d[d$MAF_bin=="MAF 1-5%",], thresh), type="histogram")
	p6 = ggMarginal(getPlot2(d[d$MAF_bin=="MAF >5%",],thresh), type="histogram")
	grid.arrange(p1, p2, p3, p4, p5, p6, nrow=2)
}

png("mendel_error_vs_imp_rsq_FAMTHRESH_histograms.png", width=720, height=480)
getPlotGrid(1)
dev.off()

png("mendel_error_vs_imp_rsq_CCTHRESH_histograms.png", width=720, height=480)
getPlotGrid(5)
dev.off()

