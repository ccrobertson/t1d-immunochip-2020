library(ggplot2)
library(dplyr)
setwd(Sys.getenv("impbench"))

# #Compare coverage by threshold
# topmed.coverage = read.table("coverage_by_imprsq_threshold_IMPUTED_TOPMED.txt", header=TRUE)
# kg.coverage = read.table("coverage_by_imprsq_threshold_IMPUTED_1KG.txt", header=TRUE)
# coverage = rbind(data.frame(panel="topmed", topmed.coverage), data.frame(panel="kg",kg.coverage))
# coverage$panel = factor(coverage$panel, levels=c("topmed","kg"), labels=c("TOPMed", "1000G Phase 3"))
# pdf("compare_coverage_by_threshold.pdf", width=12, height=8)
# ggplot(coverage) + geom_point(aes(x=thresholds, y=coverage_Ersq, color=panel)) +
# 		facet_grid(.~group) +
# 		xlab("Imputation R-squared filtering threshold") + ylab("Percent of variants in WGS") +
# 		theme_bw()
# dev.off()


#Compare coverage by threshold across different MAF spectrums
topmed.coverage.rare = read.table("coverage_by_imprsq_threshold.maf_0-0.005_IMPUTED_TOPMED.txt", header=TRUE)
topmed.coverage.lowfreq1 = read.table("coverage_by_imprsq_threshold.maf_0.005-0.01_IMPUTED_TOPMED.txt", header=TRUE)
topmed.coverage.lowfreq2 = read.table("coverage_by_imprsq_threshold.maf_0.01-0.05_IMPUTED_TOPMED.txt", header=TRUE)
topmed.coverage.common = read.table("coverage_by_imprsq_threshold.maf_0.05-1_IMPUTED_TOPMED.txt", header=TRUE)

kg.coverage.rare = read.table("coverage_by_imprsq_threshold.maf_0-0.005_IMPUTED_1KG.txt", header=TRUE)
kg.coverage.lowfreq1 = read.table("coverage_by_imprsq_threshold.maf_0.005-0.01_IMPUTED_1KG.txt", header=TRUE)
kg.coverage.lowfreq2 = read.table("coverage_by_imprsq_threshold.maf_0.01-0.05_IMPUTED_1KG.txt", header=TRUE)
kg.coverage.common = read.table("coverage_by_imprsq_threshold.maf_0.05-1_IMPUTED_1KG.txt", header=TRUE)

coverage_by_maf_range = rbind(
		data.frame(panel="topmed", maf_range="0-0.005", topmed.coverage.rare),
		data.frame(panel="topmed", maf_range="0.005-0.01", topmed.coverage.lowfreq1),
		data.frame(panel="topmed", maf_range="0.01-0.05", topmed.coverage.lowfreq2),
		data.frame(panel="topmed", maf_range="0.05-1", topmed.coverage.common),
		data.frame(panel="kg",maf_range="0-0.005", kg.coverage.rare),
		data.frame(panel="kg",maf_range="0.005-0.01", kg.coverage.lowfreq1),
		data.frame(panel="kg",maf_range="0.01-0.05", kg.coverage.lowfreq2),
		data.frame(panel="kg",maf_range="0.05-1", kg.coverage.common))
coverage_by_maf_range$panel = factor(coverage_by_maf_range$panel, levels=c("topmed","kg"), labels=c("TOPMed", "1000G Phase 3"))
coverage_by_maf_range$maf_range = factor(coverage_by_maf_range$maf_range, levels=c("0-0.005", "0.005-0.01", "0.01-0.05","0.05-1"), labels=c("MAF < 0.5%", "MAF 0.5 - 1%","MAF 1 - 5%","MAF > 5%"))

pdf("compare_coverage_by_threshold_and_maf_range.pdf", width=12, height=8)
coverage_by_maf_range$group_factor = factor(coverage_by_maf_range$group, levels=c("AFR","AMR","EUR"), labels=c("African-American (N=1,411)", "Admixed (N=641)", "European"))
ggplot(coverage_by_maf_range[coverage_by_maf_range$group%in%c("AFR","AMR"),]) + geom_point(aes(x=thresholds, y=coverage_Ersq, color=panel)) +
		facet_grid(maf_range~group_factor) +
		xlab("Imputation R-squared filtering threshold") + ylab("Percent of variants in WGS") +
		theme_bw() +
		theme(legend.title=element_blank())
# ggplot(coverage_by_maf_range[!coverage_by_maf_range$maf_range=="MAF < 0.5%" & coverage_by_maf_range$group%in%c("AFR","AMR"),]) + geom_point(aes(x=thresholds, y=coverage_Ersq, color=panel)) +
# 		facet_grid(maf_range~group) +
# 		xlab("Imputation R-squared filtering threshold") + ylab("Percent of variants in WGS") +
# 		theme_bw() +
# 		theme(legend.title=element_blank())
dev.off()

#Compare coverage by region
topmed.regions = read.table("coverage_by_region_rsqIMPUTED_TOPMED.txt", header=TRUE)
kg.regions = read.table("coverage_by_region_rsqIMPUTED_1KG.txt", header=TRUE)
coverage_by_region = rbind(data.frame(panel="topmed", topmed.regions), data.frame(panel="kg", kg.regions))
coverage_by_region$panel = factor(coverage_by_region$panel, levels=c("topmed","kg"), labels=c("TOPMed", "1000G Phase 3"))
pdf("compare_coverage_by_region.pdf", width=12, height=8)
ggplot(coverage_by_region) + geom_boxplot(aes(x=group, y=prop_kept, fill=panel)) + theme_bw() +
		xlab("") + ylab("Percent of variants in WGS") +
		ggtitle("Coverage of variants with MAF>0.005 after filtering for Rsq>0.8")
ggplot(coverage_by_region) + geom_point(aes(x=prop_typed, y=prop_kept, color=panel)) +
		facet_grid(panel~group) + ylim(0,1) + xlim(0,1) +
		geom_abline(intercept=0, slope=1) + theme_bw() +
		ggtitle("Proportion of variants in WGS by region") +
		xlab("Immunochip") +
		ylab("Imputation")
ggplot(coverage_by_region) + geom_point(aes(x=total_wgs, y=total_kept, color=panel)) +
		facet_grid(panel~group) +
		geom_abline(intercept=0, slope=1) + theme_bw() +
		ggtitle("Total number of variants in region") + xlab("WGS") + ylab("Imputation")
dev.off()
