library(GenomicRanges)
setwd(Sys.getenv("freeze"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get peaks with significant meta caQTLs that are T1D credible variants
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
res_meta = read.table("caqtl_scan_unstim_cis_meta.txt", header=TRUE)
finemap = readRDS(file=paste0(Sys.getenv("finemap"),"/annotated_credible_sets.rds"))
credible_variants = gsub(":","_",finemap$MarkerName[finemap$ppsum>0.8])
significant_peaks = unique(res_meta[res_meta$SNP%in%credible_variants & res_meta$pval.fixed<5e-5,"gene"])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get datasets
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat_un = readRDS("DATASET_un.rds")
#dat_stim = readRDS("DATASET_stim.rds")
dat_un_EUR = readRDS("DATASET_un_EUR.rds")
dat_un_AFR = readRDS("DATASET_un_AFR.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create variant lookup object
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (file.exists("mega_variant_info_granges_object.rds")) {

  variantInfoGR = readRDS(file="mega_variant_info_granges_object.rds")

} else {

  variants = colnames(dat_un[["dosages"]])
  variantInfo = do.call("rbind",lapply(variants, function(x) {y = unlist(strsplit(x, split="_")); if (length(y)==4) {return(c(x, unlist(y)))}}))
  annot = read.csv(paste0(Sys.getenv("annot"), "/annovar_output.hg38_multianno.csv"))
  annot$AltID = paste0("chr",annot$Chr,"_",annot$Start,"_",annot$Ref,"_",annot$Alt)
  annot = annot[!duplicated(annot$AltID),]
  row.names(annot) = annot$AltID
  variantInfoGR = GRanges(seqnames = variantInfo[,2], ranges = IRanges(start=as.numeric(variantInfo[,3]), end=as.numeric(variantInfo[,3])), id=variantInfo[,1], rsid = annot[variantInfo[,1],"avsnp150"])
  saveRDS(variantInfoGR, file="mega_variant_info_granges_object.rds")

}





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Conditional analysis functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getPeakObject = function(peak) {
  vals = strsplit(peak, split="_")[[1]]
  list(label=peak, chr=vals[[1]], start=as.numeric(vals[[2]]), end=as.numeric(vals[[3]]))
}

assoc = function(model) {
  coeffMat = summary(lm(model))$coefficients
  if ("x" %in% row.names(coeffMat)) {
    out = coeffMat["x",]
  } else {
    out = rep(NA, dim(coeffMat)[2])
  }
}

getProxies = function(index, X) {
  ld = apply(X, 2, function(x) {cor(x,X[,index], use="complete.obs")^2})
  return(names(ld[ld>0.9 & !is.na(ld)]))
}


conditionalAnalysisOneRound = function(peak, dat_sub, flank = 5e5, condition_snps = NULL, covariates) {
  #peak is a peak object as defined by function "getPeakObject"
  #rows of sampleDF_sub and genoDF_sub must be in the same order as columns of countMat_sub

  #sampleDF_sub = dat_sub[["samples"]]
  sampleDF_sub = data.frame(dat_sub[["samples"]], dat_sub[["pheno"]], dat_sub[["peer"]])
  countMat_sub = dat_sub[["counts"]]
  genoDF_sub = dat_sub[["dosages"]]

  #extract variants in the region surrounding peak
  region = GRanges(seqnames = peak[["chr"]], ranges = IRanges(start=peak[["start"]]-flank, end=peak[["end"]]+flank))
  variant_overlap = subsetByOverlaps(variantInfoGR, region, ignore.strand=TRUE)$id

  #deal with colinearity with condition snps
  if (length(condition_snps)>0) {
    #ldmat = cor(X, use="complete.obs")
    #X_cond = data.frame(genoDF_sub[,condition_snps])
    X = genoDF_sub[,variant_overlap]
    proxies = lapply(condition_snps, getProxies, X)
    variant_keeplist = variant_overlap[!variant_overlap %in% c(condition_snps, unlist(proxies))]
  } else {
    variant_keeplist = variant_overlap[!variant_overlap %in% condition_snps]
  }

  #setup data objects
  X_cond = data.frame(genoDF_sub[,condition_snps]); names(X_cond) <- condition_snps
  X = data.frame(genoDF_sub[,variant_keeplist]); names(X) <- variant_keeplist
  Y = countMat_sub[peak[["label"]],]
  cov = sampleDF_sub[,covariates]

  #define regression model
  if (length(condition_snps)>0) {
    lm_formula = paste0("Y ~ x + ",paste(paste0("X_cond$",names(X_cond)),collapse=" + ")," + ",paste(paste0("cov$",names(cov)),collapse=" + "))
  } else {
    lm_formula = paste0("Y ~ x + ", paste(paste0("cov$",names(cov)),collapse=" + "))
  }

  #run model on each snp
  OUT = data.frame(t(apply(X, 2, function(x) {assoc(as.formula(lm_formula))})))
  OUT$formula = lm_formula
  names(OUT) <- c("Beta","StdErr","tvalue","pvalue", "formula")
  return(OUT[!is.na(OUT$Beta),])
}

conditionalAnalysisStepwise = function(peak, dat_sub, covariates) {

  #initialize
  cat(peak$label,"\n")
  output = list()

  #run first round
  round=1; cat("round",round,"\n")
  condition_snps = NULL
  res = conditionalAnalysisOneRound(peak=peak, dat_sub=dat_sub, covariates=covariates)
  index = row.names(res)[res$pvalue==min(res$pvalue)][1]

  #run additional rounds while significant
  while (res[index,"pvalue"]<1e-3) {
    output[[index]] = res[index,]
    round=round+1; cat("round",round,"\n")
    condition_snps = unique(c(index, condition_snps))
    res = conditionalAnalysisOneRound(peak=peak, dat_sub=dat_sub, condition_snps=condition_snps, covariates=covariates)
    index = row.names(res)[res$pvalue==min(res$pvalue)][1]
  }

  #return results
  if (length(output)>0) {
    outDF = data.frame(peak=peak$label, do.call("rbind",output))
    outDF$SNP = row.names(outDF)
    OUT = outDF[,c("peak","SNP","Beta","StdErr","tvalue","pvalue","formula")]
  } else {
    OUT = NULL
  }
  return(list(peak=peak, stepwise=OUT))
}

conditionalAnalysisJointModel = function(result, dat_sub, covariates) {

  peak = result[["peak"]]
  output = list()
  if (!is.null(result[["stepwise"]])) {
    index_variants = result[["stepwise"]]$SNP
    if (length(index_variants)>1) {
      for (i in 1:length(index_variants)) {
        index = index_variants[i]
        condition_snps = index_variants[!index_variants==index]
        res = conditionalAnalysisOneRound(peak=peak, dat_sub=dat_sub, condition_snps=condition_snps, covariates=covariates)
        res$peak = peak$label
        res$SNP = row.names(res)
        output[[index]] = res[,c("peak","SNP","Beta","StdErr","tvalue","pvalue","formula")]
      }
    } else {
      index = index_variants[1]
      res = conditionalAnalysisOneRound(peak=peak, dat_sub=dat_sub, covariates=covariates)
      res$peak = peak$label
      res$SNP = row.names(res)
      output[[index]] = res[,c("peak","SNP","Beta","StdErr","tvalue","pvalue","formula")]
    }
  }
  return(output)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run stepwise conditional analysis for significant peaks
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
peaks = lapply(significant_peaks, getPeakObject)
stepwise_results_un = lapply(peaks, conditionalAnalysisStepwise, dat_sub=dat_un, covariates=c("PC1","PC2", "age_sample", "TSS_Score"))
stepwise_results_un_EUR = lapply(peaks, conditionalAnalysisStepwise, dat_sub=dat_un_EUR, covariates=c("PC1_ancestry_specific","PC2_ancestry_specific", "age_sample", "TSS_Score"))
stepwise_results_un_AFR = lapply(peaks, conditionalAnalysisStepwise, dat_sub=dat_un_AFR, covariates=c("PC1_ancestry_specific","PC2_ancestry_specific", "age_sample", "TSS_Score"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fit joint models and save summary stats by index variant
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
joint_results_un = unlist(lapply(stepwise_results_un, conditionalAnalysisJointModel, dat_sub=dat_un, covariates=c("PC1","PC2", "age_sample", "TSS_Score")), recursive=FALSE)
saveRDS(joint_results_un, file="joint_caqtl_models_un.rds")


# TRY JOINT MODEL IN EUR AND AFR SAMPLES USING INDEX VARIANTS FROM META-ANALYSIS
joint_results_un_EUR = unlist(lapply(stepwise_results_un, conditionalAnalysisJointModel, dat_sub=dat_un_EUR, covariates=c("PC1_ancestry_specific","PC2_ancestry_specific", "age_sample", "TSS_Score")), recursive=FALSE)
saveRDS(joint_results_un_EUR, file="joint_caqtl_models_EUR.rds")

joint_results_un_AFR = unlist(lapply(stepwise_results_un, conditionalAnalysisJointModel, dat_sub=dat_un_AFR, covariates=c("PC1_ancestry_specific","PC2_ancestry_specific", "age_sample", "TSS_Score")), recursive=FALSE)
saveRDS(joint_results_un_AFR, file="joint_caqtl_models_AFR.rds")
