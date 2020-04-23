library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)

setwd(Sys.getenv("inflation"))

#get family results
load(paste0(Sys.getenv("meta"),"/combined_res_newqc_5pc_diff.RData"))
csq_fam = (f1$beta_EUR_fams/f1$se_EUR_fams)^2
csq_fam = csq_fam[!is.na(csq_fam)]
median(csq_fam)/qchisq(0.5,1)

#remove MHC (chr6:20-40MB) and INS (chr11:1-3MB) region variants
remove_MHC_INS_PTP = function(vars) {
	varDF = as.data.frame(str_split_fixed(vars, fixed(":"), 4))
	names(varDF) <- c("CHR","BP","REF","ALT")
	mhc_vars = vars[varDF$CHR=="chr6" & varDF$BP>20000000 & varDF$BP<40000000]
	ins_vars = vars[varDF$CHR=="chr11" & varDF$BP>1000000 & varDF$BP<3000000]
	ptp_vars = vars[varDF$CHR=="chr1" & varDF$BP>112000000 & varDF$BP<115000000]
	keep = vars[!(vars %in% c(mhc_vars, ins_vars, ptp_vars))]
	return(keep)
}

#get list of typed variants
typed_vars = list()
for (group in c("EUR1", "EUR2","EUR3","AMR","AFR","FIN")) {
	info = read.table(paste0(Sys.getenv("results"),"/",group,"/typed_variants.txt"), header=TRUE)
	typed_vars[[group]] = info$SNP[info$Genotyped=="Genotyped"]
}
typed_vars[["EUR"]] = unique(c(typed_vars[["EUR1"]], typed_vars[["EUR2"]], typed_vars[["EUR3"]]))


#compute lambda and plot qq from metal results
getLambda = function(d, snp_col="SNP", beta_col="BETA", se_col="SE", title, sample_size) {

	chisq_all = (d[,beta_col]/d[,se_col])^2
	chisq_obs = sort(chisq_all[!is.na(chisq_all)])
	n_vars=length(chisq_obs)
	n_cases=sample_size/2
	n_controls=sample_size/2
	chisq_null = sort(rchisq(n=n_vars, df=1, ncp=0))
	lambda = round(median(chisq_obs)/qchisq(0.5,1), digits=2)
	lambda_1k = round(1 + (lambda-1)*((1/n_cases + 1/n_controls)/(1/1000 + 1/1000)), digits=2)
	dat = data.frame(chisq_null, chisq_obs)
	plottitle=paste0(title, "\n lambda =", lambda, "\nlambda_1000 =", lambda_1k)
	p = ggplot(dat) + geom_point(aes(chisq_null, chisq_obs), size=0.5) +
			theme_classic() +
			geom_abline(intercept=0, slope=1) +
			ggtitle(plottitle)
	return(list(lambda, p))
}

#compare fams to five subsets
getComparison = function(group, nickname, remove_mhc=FALSE, typed_only=FALSE, intersectOrUnion) {

	png(paste0("inflation_analysis_",nickname,".png"), width=720, height=480)
	sample_size <- as.numeric(unlist(strsplit(system(paste0("wc -l ",group,"_sub1.txt"), intern=TRUE),split=" "))[1])

	#get intersection between fam and cc variants
	datfam_tmp=read.table(paste0(Sys.getenv("meta"),"/meta_file_fam_",group,"_newqc_1.tbl"), header=TRUE)
	datcc_tmp = read.table(paste0("subanalyses/",sub,"/metal_",group,"_sub",sub,".tbl"), header=TRUE)
	varintersect = intersect(datfam_tmp$SNP[!is.na(datfam_tmp$BETA)], datcc_tmp$SNP[!is.na(datcc_tmp$BETA)])
	varunion = union(datfam_tmp$SNP[!is.na(datfam_tmp$BETA)], datcc_tmp$SNP[!is.na(datcc_tmp$BETA)])

	#process fam data
	datfam_1=read.table(paste0(Sys.getenv("meta"),"/meta_file_fam_",group,"_newqc_1.tbl"), header=TRUE)
	if (typed_only) {
		typed_list=typed_vars[[group]]
		datfam_2 = datfam_1[datfam_1[,"SNP"]%in%typed_list,]
	} else {
		datfam_2 = datfam_1
	}
	if (remove_mhc) {
		keep = remove_MHC_INS_PTP(datfam_2[,"SNP"])
		datfam_3 = datfam_2[datfam_2[,"SNP"]%in%keep,]
	} else {
		datfam_3 = datfam_2
	}
	if (intersectOrUnion=="intersect") {
		datfam = datfam_3[datfam_3$SNP %in% varintersect,]
	} else if (intersectOrUnion=="union") {
		datfam = datfam_3[datfam_3$SNP %in% varunion,]
	}
	pFam=getLambda(datfam,title=paste0(group," fams"), sample_size=sample_size)

	#process cc data
	pCC=list()
	for (sub in 1:5) {
		datcc_1 = read.table(paste0("subanalyses/",sub,"/metal_",group,"_sub",sub,".tbl"), header=TRUE)
		if (typed_only) {
			typed_list=typed_vars[[group]]
			datcc_2 = datcc_1[datcc_1[,"SNP"]%in%typed_list,]
		} else {
			datcc_2 = datcc_1
		}
		if (remove_mhc) {
			keep = remove_MHC_INS_PTP(datcc_2[,"SNP"])
			datcc_3 = datcc_2[datcc_2[,"SNP"]%in%keep,]
		} else {
			datcc_3 = datcc_2
		}
		if (intersectOrUnion=="intersect") {
			datcc = datcc_3[datcc_3$SNP %in% varintersect,]
		} else if (intersectOrUnion=="union") {
			datcc = datcc_3[datcc_3$SNP %in% varunion,]
		}
		pCC[[sub]] = getLambda(datcc,title=paste0(group, " CC - sub",sub), sample_size=sample_size)

	}

	#print plots
	plottitle=paste0(sample_size/2," cases, ", sample_size/2, " controls \n versus \n", sample_size/2, " trios")
	grid.arrange(pFam[[2]], pCC[[1]][[2]],pCC[[2]][[2]],pCC[[3]][[2]],pCC[[4]][[2]],pCC[[5]][[2]], nrow=2, top=plottitle)

	dev.off()
	return(list(sample_size/2, pFam, pCC))
}

summarizeResults = function(remove_mhc, typed_only, intersectOrUnion) {

	if (remove_mhc) {
		p1 = "noMHC"
	} else {
		p1 = "withMHC"
	}
	if (typed_only) {
		p2 = "genotyped"
	} else {
		p2 = "imputed"
	}
	p3 = intersectOrUnion

	prefix = paste(p1, p2, p3, sep="_")
	print(prefix)

	inflationResults = list()
	for (group in c("EUR","AFR","FIN","AMR")) {
		print(group)
		inflationResults[[group]] = getComparison(group, paste(group, prefix, sep="_"), remove_mhc=remove_mhc, typed_only=typed_only, intersectOrUnion=intersectOrUnion)
	}

	for (group in c("EUR","AFR","FIN","AMR")) {
		vec = unlist(lapply(inflationResults[[group]][[3]], function(x) {x[[1]]}))
		inflationResults[[group]][["lambda_cc"]] = paste0(round(mean(vec), digits=2), " (", round(sd(vec), digits=2), ")")
	}

	data.frame( sample_size=unlist(lapply(inflationResults, function(x) {x[[1]]})),
		lambda_fam=unlist(lapply(inflationResults, function(x) {x[[2]][[1]]})),
		lambda_cc=unlist(lapply(inflationResults, function(x) {x[["lambda_cc"]]}))
	)

}


summarizeResults(remove_mhc=TRUE, typed_only=TRUE, intersectOrUnion="intersect")
summarizeResults(remove_mhc=TRUE, typed_only=FALSE, intersectOrUnion="intersect")




# getComparison("FIN", "FIN_noMHC_genotyped", remove_mhc=TRUE, typed_only=TRUE)
# getComparison("AFR", "AFR_noMHC_genotyped", remove_mhc=TRUE, typed_only=TRUE)
# getComparison("AMR", "AMR_noMHC_genotyped", remove_mhc=TRUE, typed_only=TRUE)
#
# getComparison("FIN", "FIN_withMHC_genotyped", remove_mhc=FALSE, typed_only=TRUE)
# getComparison("AFR", "AFR_withMHC_genotyped", remove_mhc=FALSE, typed_only=TRUE)
# getComparison("AMR", "AMR_withMHC_genotyped", remove_mhc=FALSE, typed_only=TRUE)
#
#
# getComparison("EUR", "EUR_noMHC_imputed", remove_mhc=TRUE, typed_only=FALSE)
# getComparison("EUR", "EUR_withMHC_imputed", remove_mhc=FALSE, typed_only=FALSE)
