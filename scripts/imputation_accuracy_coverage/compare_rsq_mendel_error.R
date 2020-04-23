
setwd("/nv/vol185/MEGA/IMPUTED_TOPMED_HWE/coverage")

#get data sets
dataSets=list()
for (i in 1:22) {
	
	cat("chr",i,"\n")
	#ALL: imputed variants with rsq>0.3
	info = read.table(paste0("../results/EUR/chr",i,".merged_info"), header=TRUE)
	row.names(info) <- info$SNP
	
	#ICHIP REGIONS: union of (imputed variants with rsq>0.3) + (wgs variants) 
	rsq = read.table(paste0("imputation_vs_wgs_EUR_chr",i,"_ichip.txt"), header=TRUE)
	row.names(rsq) <- rsq$variant
	
	#ALL: imputed variants with maf>0.005 and rsq>0.8
	mend =  read.table(paste0("../family_analysis/EUR/king_snpqcbySNP_chr",i,".txt"), header=TRUE)
	row.names(mend) <- mend$SNP
	
	#WANT INTERSECTION: imputed variants in ichip regions with rsq>0.8 & maf>0.005
	vlist1 = row.names(info)[info$MAF>0.005 & info$rsq>0.8] #imputed variants with maf>0.005 & rsq>0.8
	vlist2 = row.names(rsq)[rsq$ichip==1] #imputed variants in ichip regions
	keep = intersect(vlist1, vlist2)
	
	#compare maf and rsq estimates based on subset used for wgs vs imp comparison to those from full EUR cohort
	#d = data.frame(SNP=keep, mend[keep,c("Err_InHetTrio","Err_InHomPO")], info[keep,c("MAF","rsq")], rsq[keep,c("ichip_maf","imp_Rsq1","imp_Rsq2")])
	#plot(d$rsq, d$imp_Rsq1)
	#plot(d$MAF, d$ichip_maf)
	
	#Visualize relationship between imputation Rsq and Mendel Error Rates
	d = data.frame(SNP=keep, mend[keep,c("Err_InHetTrio","Err_InHomPO")], info[keep,c("MAF","rsq")])
	dataSets[[i]] = d

}
d = do.call("rbind",dataSets)
d$MAF_bin = NA
d$MAF_bin[d$MAF<=0.01] <- "MAF 0.5-1%"
d$MAF_bin[d$MAF>0.01 & d$MAF<=0.05] <- "MAF 1-5%"
d$MAF_bin[d$MAF>0.05] <- "MAF >5%"
d$MAF_bin = factor(d$MAF_bin, levels=c("MAF 0.5-1%", "MAF 1-5%", "MAF >5%"))
write.table(d, "mendel_error_vs_imp_rsq.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")


