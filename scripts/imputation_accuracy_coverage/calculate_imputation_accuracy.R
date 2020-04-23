options(stringsAsFactors=FALSE)

library(vcfR)
library(ggplot2)

setwd(Sys.getenv("impbench"))

args = commandArgs(trailingOnly=TRUE)

group = args[1]
refpanel = args[2]
subset = args[3]
chr = paste0("chr",args[4])


### Read data
if (subset == "all") {

	cat("STARTING: calculating WGS concordance for ", subset," sites\n")
	wgs = read.vcfR(paste("washu",group,chr,"filtered.vcf.gz", sep="_"))
	mega = read.vcfR(paste("mega",refpanel,group,chr,"gr38.vcf.gz",sep="_"))

} else if (subset == "ichip") {

	cat("STARTING: calculating WGS concordance for ", subset," sites\n")
	wgs = read.vcfR(paste("washu",group,chr,"filtered_ichip_regions.vcf.gz",sep="_"))
	mega = read.vcfR(paste("mega",refpanel,group,chr,"ichip_regions.vcf.gz",sep="_"))
}

mi_error = read.table(paste0(Sys.getenv("fam_assoc"),"/",group,"/king_snpqcbySNP_",chr,".txt"), header=TRUE)
row.names(mi_error) = mi_error$SNP

### Extract VCF components
#get genotypes
wgs_gt <- extract.gt(wgs, element="GT")
mega_gt <- extract.gt(mega, element="GT")
mega_ds <- extract.gt(mega, element="DS", as.numeric=TRUE)

#get variant info
mega_chrpos = as.data.frame(getFIX(mega))
row.names(mega_chrpos) = as.character(mega_chrpos$ID)
mega_chrpos$REF = gsub(":","_", mega_chrpos$REF)
mega_chrpos$ALT = gsub(":","_", mega_chrpos$ALT)
mega_chrpos$b38_variant_id = paste(mega_chrpos$CHROM, mega_chrpos$POS, mega_chrpos$REF, mega_chrpos$ALT, sep=":")

### Convert to genotype matrix to counts (ASSUMES biallelic sites)
wgs_gt[wgs_gt=="./."] <- NA
wgs_gt[wgs_gt=="0/."] <- NA
wgs_gt[wgs_gt=="./0"] <- NA
wgs_gt[wgs_gt=="1/."] <- NA
wgs_gt[wgs_gt=="./1"] <- NA
wgs_gt[wgs_gt=="0/0"] <- 0
wgs_gt[wgs_gt=="0/1"] <- 1
wgs_gt[wgs_gt=="1/0"] <- 1
wgs_gt[wgs_gt=="1/1"] <- 2

mega_gt[mega_gt=="0|0"] <- 0
mega_gt[mega_gt=="0|1"] <- 1
mega_gt[mega_gt=="1|0"] <- 1
mega_gt[mega_gt=="1|1"] <- 2


### Make genotype matrices compatable
#update variant names
row.names(wgs_gt) <- gsub("-",":",row.names(wgs_gt))
row.names(mega_gt) <- mega_chrpos[row.names(mega_gt), "b38_variant_id"]
row.names(mega_ds) <- mega_chrpos[row.names(mega_ds), "b38_variant_id"]
common_variants = row.names(mega_gt)[row.names(mega_gt) %in% row.names(wgs_gt)]

#update washu subject names (NOTE -- should include a KEY here to maximize washu sample size!!!!)
colnames(wgs_gt) <- gsub("H_VD-","", colnames(wgs_gt))
common_subjects = colnames(mega_gt) [colnames(mega_gt) %in% colnames(wgs_gt)]



### Get info data on each data set
getMissing = function(x) {
	sum(is.na(x))/length(x)
}
getAF = function(x) {
	x = as.numeric(x)
	sum(x, na.rm=TRUE)/(2*sum(!is.na(x)))
}
getMAF = function(x) {
	x = as.numeric(x)
	p = getAF(x)
	min(p, 1-p)
}
getImpRsq = function(x) {
	x = as.numeric(x)
	p = getAF(x)
	if (var(x,na.rm=TRUE)>0) {
		rsq = var(x, na.rm=TRUE)/(2*p*(1-p))
	} else {
		rsq = NA
	}
	min(rsq,1)
}

wgs_maf = apply(wgs_gt, MARGIN=1, FUN=getMAF)
wgs_af = apply(wgs_gt, MARGIN=1, FUN=getAF)
wgs_miss = apply(wgs_gt, MARGIN=1, FUN=getMissing)
wgs_info = data.frame(variant=row.names(wgs_gt), wgs_maf, wgs_af, wgs_miss, wgs=1)

mega_maf = apply(mega_ds, MARGIN=1, FUN=getMAF)
mega_af = apply(mega_ds, MARGIN=1, FUN=getAF)
mega_imprsq1 = apply(mega_ds, MARGIN=1, FUN=getImpRsq)
mega_imprsq2 = INFO2df(mega)$R2
mega_typed=as.numeric(grepl("TYPED",getINFO(mega)))
mega_mi_hettrio=mi_error[row.names(mega_gt),"Err_InHetTrio"]
mega_mi_hompo=mi_error[row.names(mega_gt),"Err_InHomPO"]
mega_info = data.frame(variant=row.names(mega_gt), mega_maf, mega_af, imp_Rsq1=mega_imprsq1, imp_Rsq2=mega_imprsq2, mega=1, TYPED=mega_typed, mega_mi_hettrio, mega_mi_hompo)

merged_info = merge(wgs_info, mega_info, by="variant", all=TRUE)
chrpos = do.call("rbind",sapply(merged_info$variant, strsplit, split=":"))
merged_info$chr = as.character(chrpos[,1])
merged_info$position = as.numeric(chrpos[,2])
merged_info$ref_allele = as.character(chrpos[,3])
merged_info$alt_allele = as.character(chrpos[,4])
merged_info$mega[is.na(merged_info$mega)] <- 0
merged_info$wgs[is.na(merged_info$wgs)] <- 0


### Calculate empirical genomic Rsq
getCorrelation = function(var, trueGT, impDS) {

	if (var %in% row.names(trueGT) & var %in% row.names(impDS)) {
		x1 = as.numeric(trueGT[var,])
		x2 = as.numeric(impDS[var,])
		if (sd(x1, na.rm=TRUE)>0 & sd(x2, na.rm=TRUE)>0) {
			r = cor(x1, x2, use="complete.obs")
		} else {
			r = NA
		}

	} else {
		r = NA
	}
	r

}
merged_info$genomic_R = sapply(merged_info$variant, FUN=getCorrelation, trueGT = wgs_gt[common_variants, common_subjects], impDS = mega_ds[common_variants, common_subjects])
merged_info$genomic_Rsq = merged_info$genomic_R^2

#assign each variant to a region
regions_all = read.table(paste0(Sys.getenv("PROJECT_MEGA"),"/define_ichip_regions/ichip_regions_gr38.txt"), header=TRUE)
regions = regions_all[regions_all$chromosome==chr,]
names(regions) <- c("region","chr","start","end")
for (i in 1:nrow(regions)) {
	cat(i,"\n")
	REGION_LABEL = regions[i,"region"]
	START = regions[i,"start"]
	END = regions[i,"end"]
	merged_info$region[merged_info$position>=START & merged_info$position<=END] <- REGION_LABEL
}
table(merged_info$region)
sum(is.na(merged_info$region))


### Save data
save(merged_info, file=paste("imputation_vs_wgs_",refpanel,"_",group,"_",chr,"_",subset,"_regions.RData",sep=""))
write.table(merged_info, file=paste("imputation_vs_wgs_",refpanel,"_",group,"_",chr,"_",subset,"_regions.txt",sep=""), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)


### Are some of the non-overlapping variants due to strand swaps? NO, doesn't look like it, REF allele matches for all sites where chro:pos matches betwen imputed and wgs data
#wgs_only = merged_info[merged_info$ichip==0 & merged_info$wgs==1,]
#ichip_only = merged_info[merged_info$ichip==1 & merged_info$wgs==0,]
#sum(wgs_only$position %in% ichip_only$position)
#
#wgs_FIX = as.data.frame(getFIX(wgs))
#wgs_FIX$chrpos = paste(wgs_FIX$CHROM, wgs_FIX$POS, sep=":")
#wgs_FIX$chrposref = paste(wgs_FIX$CHROM, wgs_FIX$POS, wgs_FIX$REF, sep=":")
#wgs_FIX$chrposrefalt = paste(wgs_FIX$CHROM, wgs_FIX$POS, wgs_FIX$REF, wgs_FIX$ALT, sep=":")
#
#ichip_FIX = as.data.frame(getFIX(ichip))
#ichip_FIX$chrpos = paste(ichip_FIX$CHROM, ichip_FIX$POS, sep=":")
#ichip_FIX$chrposref = paste(ichip_FIX$CHROM, ichip_FIX$POS, ichip_FIX$REF, sep=":")
#ichip_FIX$chrposrefalt = paste(ichip_FIX$CHROM, ichip_FIX$POS, ichip_FIX$REF, ichip_FIX$ALT, sep=":")
#
#sum(ichip_FIX$chrpos %in% wgs_FIX$chrpos)
#sum(ichip_FIX$chrposref %in% wgs_FIX$chrposref)
#sum(ichip_FIX$chrposrefalt %in% wgs_FIX$chrposrefalt)

### What were some good SNPs excluded in initial filter? Probably due to differences in MAF in full cohort vs subset used for WGS
#good = ichip_info[ichip_info$imp_Rsq2>0.8 & ichip_info$ichip_maf>0.01,]
#good[is.na(good$ichip_mi_hettrio),]
#ichip_info["chr22:39376024:T:C",]
#merged_info[merged_info$variant=="chr22:39376024:T:C",]
