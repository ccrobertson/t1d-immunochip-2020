setwd(Sys.getenv("freeze"))

#Install Jamie's branch of coloc
#library(devtools)
#install_github("jinshaw16/coloc")

library(coloc)
library(locuscomparer)


# Calculate r2 between each pair of caqtl and t1d index variants
# For each peak-t1d_index pair with r2>0.5
#   - extract jointmodel caqtl summary stats for region
#   - extract jointmodel t1d summary stats for region (from conditional results files, or unconditional results files, depending on number of causal variants in region)
#   - define overlapping variant set between caqtl and t1d
#   - test for colocalisation


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get T1D index variants and assoc data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
finemap = readRDS(file=paste0(Sys.getenv("finemap"),"/annotated_credible_sets.rds"))
t1d_index_variants = gsub(":","_",unique(finemap$tag[finemap$ppsum>0.8]))

annot =  read.csv(paste0(Sys.getenv("annot"),"/annovar_output.hg38_multianno.csv"))
annot$MarkerName = paste(paste0("chr",annot$Chr),annot$Start, annot$Ref,annot$Alt, sep=":")
annot$MarkerName_swap = paste(paste0("chr",annot$Chr),annot$Start, annot$Alt,annot$Ref, sep=":")

#get gwas results from each independent analysis (broken down by ancestry and case-control vs family)
load(paste0(Sys.getenv("meta"),"/combined_res_newqc_5pc_diff.RData"))
f1 = merge(f1, annot[,c("Chr","Start","avsnp150")], by.x=c("chromosome","position"), by.y=c("Chr","Start"))
cond_dir = paste0(Sys.getenv("meta"),"/conditional_coloc_1")
cond_files = list.files(cond_dir)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get caqtl index sets
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
results_qtl_un = readRDS(file="joint_caqtl_models_un.rds")
results_qtl_AFR = readRDS(file="joint_caqtl_models_AFR.rds")
results_qtl_EUR = readRDS(file="joint_caqtl_models_EUR.rds")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get datasets
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat_un = readRDS("DATASET_un.rds")
dat_un_EUR = readRDS("DATASET_un_EUR.rds")
dat_un_AFR = readRDS("DATASET_un_AFR.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Colocalisation functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getLD = function(snpList1, snpList2, dat_sub) {
  genoDF = dat_sub[["genotypes"]]
  pairwiseLD = NULL
  for (snp1 in snpList1) {
    for(snp2 in snpList2) {
      if (snp1 %in% names(genoDF) & snp2 %in% names(genoDF)) {
        rsq = (cor(genoDF[,snp1], genoDF[,snp2], use="complete.obs"))^2
      } else {
        rsq = NA
      }
      pairwiseLD = rbind(pairwiseLD, data.frame(snp1=snp1, snp2=snp2, rsq=rsq))
    }
  }
  return(pairwiseLD)
}

getRsid = function (x) {
  x = gsub("_",":", x)
  if (x %in% finemap$MarkerName) {
    if (!finemap$ID[!is.na(finemap$MarkerName) & finemap$MarkerName==x]==".") {
      rsid = finemap$ID[!is.na(finemap$MarkerName) & finemap$MarkerName==x][1]
    } else {
      rsid = finemap$MarkerID[!is.na(finemap$MarkerName) & finemap$MarkerName==x][1]
    }
  } else if (x %in% annot$MarkerName) {
    rsid = annot$avsnp150[annot$MarkerName==x][1]
  } else if (x %in% annot$MarkerName_swap) {
    rsid = annot$avsnp150[annot$MarkerName_swap==x][1]
  } else {
    rsid = x
  }
  return(rsid)
}

getIndexPairs = function(gwas_index_variants, qtl_index_variants, ld_dat, threshold=0.5) {
  ldDF = getLD(gwas_index_variants, qtl_index_variants, ld_dat)
  index_pairs = ldDF[ldDF$rsq>threshold & !is.na(ldDF$rsq),]
  index_pairs$rsid1 = sapply(index_pairs$snp1, getRsid)
  index_pairs$rsid2 = sapply(index_pairs$snp2, getRsid)
  return(index_pairs)
}

getEUR_AF = function(x) {
  x = gsub("_",":", x)
  if (x %in% annot$MarkerName) {
    if (!annot$X1000g2015aug_eur[annot$MarkerName==x][1]==".") {
      AF = as.numeric(annot$X1000g2015aug_eur[annot$MarkerName==x][1])
    } else {
      AF = NA
    }
  } else if (x %in% annot$MarkerName_swap) {
    if (!annot$X1000g2015aug_eur[annot$MarkerName_swap==x][1]==".") {
      AF = 1-as.numeric(annot$X1000g2015aug_eur[annot$MarkerName_swap==x][1])
    } else {
      AF = NA
    }

  } else {
    AF = NA
  }
  return(AF)
}

getColocDat = function(index_gwas, index_qtl, dat_qtl, results_qtl) {

  #get gwas suummary stats
  if (getRsid(index_gwas) %in% cond_files) {
    res_gwas = read.table(paste0(cond_dir,"/",getRsid(index_gwas)), comment.char="#", header=TRUE)
    res_gwas$MarkerName = res_gwas$rsid
    res_gwas$chromosome = sapply(res_gwas$MarkerName, function(x) {gsub("chr","",strsplit(x, split=":")[[1]][1])})
    res_gwas$beta = res_gwas$frequentist_add_beta_1.add.t1d.1
    res_gwas$SE = res_gwas$frequentist_add_se_1
    res_gwas$source = "COND"
  } else {
    res_gwas = f1
    res_gwas$MarkerName = res_gwas$rsid
    res_gwas$beta = res_gwas$beta_EUR
    res_gwas$SE = res_gwas$se_EUR
    res_gwas$source = "META"
  }
  res_gwas = res_gwas[!is.na(res_gwas$beta) & !is.na(res_gwas$SE),c("MarkerName", "chromosome", "position","alleleA","alleleB","beta","SE", "source")]
  res_gwas$varbeta = (res_gwas$SE)^2
  res_gwas$pval = pchisq((res_gwas$beta/res_gwas$SE)^2, df=1,lower.tail=FALSE)

  #get qtl summary stats
  res_qtl = results_qtl[[index_qtl]]
  res_qtl$MarkerName = gsub("_",":",res_qtl$SNP)
  res_qtl$chromosome = sapply(res_qtl$SNP, function(x) {gsub("chr","",strsplit(x, split="_")[[1]][1])})
  res_qtl$position = sapply(res_qtl$SNP, function(x) {strsplit(x, split="_")[[1]][2]})
  res_qtl$alleleA = sapply(res_qtl$SNP, function(x) {strsplit(x, split="_")[[1]][3]})
  res_qtl$alleleB = sapply(res_qtl$SNP, function(x) {strsplit(x, split="_")[[1]][4]})
  res_qtl$beta = res_qtl$Beta
  res_qtl$SE = res_qtl$StdErr
  res_qtl = res_qtl[!is.na(res_qtl$beta) & !is.na(res_qtl$SE),c("MarkerName", "peak","chromosome", "position","alleleA","alleleB","beta","SE")]
  res_qtl$varbeta = (res_qtl$SE)^2
  res_qtl$pval = pchisq((res_qtl$beta/res_qtl$SE)^2, df=1,lower.tail=FALSE)

  #combine
  combdat = merge(res_gwas, res_qtl, by=c("MarkerName","chromosome","position","alleleA","alleleB"), suffixes=c(".pheno",".qtl"))
  combdat$match = paste(paste0("chr",combdat$chromosome), combdat$position, combdat$alleleA, combdat$alleleB,sep=":")==combdat$MarkerName
  combdat$rsid = sapply(combdat$MarkerName, getRsid)
  combdat$rsid[is.na(combdat$rsid) | combdat$rsid=="."] <- combdat$MarkerName[is.na(combdat$rsid) | combdat$rsid=="."]

  #run coloc
  coloc.res <- coloc.abf(p12=5e-6, dataset1=list(N=33601, beta=combdat$beta.pheno, varbeta=combdat$varbeta.pheno, type="cc", s=0.5, snp=combdat$rsid),
                      dataset2=list(N=nrow(dat_qtl[["samples"]]), beta=combdat$beta.qtl, varbeta=combdat$varbeta.qtl, type="quant", snp=combdat$rsid, sdY=sd(dat_qtl[["counts"]][combdat$peak[1],])))
  sorted_results = coloc.res$results[order(coloc.res$results$SNP.PP.H4, decreasing=TRUE),]
  return(list(index_gwas=index_gwas, index_qtl=index_qtl, coloc_summary = coloc.res$summary, coloc_results = sorted_results, coloc_input = combdat, feature=res_qtl$peak[1]))

}

getEQTL_ColocDat = function(index_gwas, gene) {

  #get gwas summary stats
  if (getRsid(index_gwas) %in% cond_files) {
    res_gwas = read.table(paste0(cond_dir,"/",getRsid(index_gwas)), comment.char="#", header=TRUE, sep="")
    res_gwas$MarkerName = res_gwas$rsid
    res_gwas$rsid = sapply(res_gwas$MarkerName, getRsid)
    res_gwas$chromosome = sapply(res_gwas$MarkerName, function(x) {gsub("chr","",strsplit(x, split=":")[[1]][1])})
    res_gwas$beta = res_gwas$frequentist_add_beta_1.add.t1d.1
    res_gwas$SE = res_gwas$frequentist_add_se_1
    res_gwas$AF = res_gwas$all_maf
    res_gwas$source = "COND"
  } else {
    res_gwas = f1
    res_gwas$MarkerName = res_gwas$rsid
    res_gwas$rsid = res_gwas$avsnp150
    res_gwas$beta = res_gwas$beta_EUR
    res_gwas$SE = res_gwas$se_EUR
    res_gwas$AF = res_gwas$AF_EUR
    res_gwas$source = "META"
  }
  res_gwas = res_gwas[!is.na(res_gwas$beta) & !is.na(res_gwas$SE),]
  res_gwas$zscore = res_gwas$beta/res_gwas$SE
  res_gwas$pval = pchisq((res_gwas$beta/res_gwas$SE)^2, df=1,lower.tail=FALSE)
  res_gwas = res_gwas[,c("rsid","MarkerName", "chromosome", "position","alleleA","alleleB","AF","zscore","pval", "source")]

  #get qtl summary stats
  res_qtl = read.table(paste0(Sys.getenv("PUBLIC_DATA"),"/eQTLGen/cis-eQTLs_full_20180905__",gene,".txt"), header=TRUE, sep="")
  res_qtl$rsid = res_qtl$SNP
  res_qtl$chromosome = res_qtl$SNPChr
  res_qtl$position_hg19 = res_qtl$SNPPos
  granges_obj_hg19 = GRanges(seqnames = paste0("chr",res_qtl$chromosome), ranges = IRanges(start=res_qtl$position_hg19-1, end=res_qtl$position_hg19), id=res_qtl$SNP)
  lifted = data.frame(liftOver(x=granges_obj_hg19, chain=chain_hg19_to_hg38))
  res_qtl$position_hg38 = sapply(res_qtl$rsid, function(x){lifted[lifted$id==x,"end"]})
  res_qtl$position = res_qtl$position_hg38
  res_qtl$pval = res_qtl$Pvalue

  #fix alleles
  res_qtl$alleleA = res_qtl$OtherAllele
  res_qtl$alleleB = res_qtl$AssessedAllele
  res_qtl$MarkerName = paste(paste0("chr",res_qtl$chromosome),res_qtl$position_hg38,res_qtl$alleleA,res_qtl$alleleB,sep=":")
  res_qtl$MarkerName_swap = paste(paste0("chr",res_qtl$chromosome),res_qtl$position_hg38,res_qtl$alleleB,res_qtl$alleleA,sep=":")

  res_qtl$swap = res_qtl$MarkerName_swap %in% res_gwas$MarkerName

  res_qtl$MarkerName_new[!res_qtl$swap] <- res_qtl$MarkerName[!res_qtl$swap]
  res_qtl$MarkerName_new[res_qtl$swap] <- res_qtl$MarkerName_swap[res_qtl$swap]
  res_qtl$MarkerName_old = res_qtl$MarkerName

  res_qtl$alleleA_new[!res_qtl$swap] <- res_qtl$alleleA[!res_qtl$swap]
  res_qtl$alleleA_new[res_qtl$swap] <- res_qtl$alleleB[res_qtl$swap]

  res_qtl$alleleB_new[!res_qtl$swap] <- res_qtl$alleleB[!res_qtl$swap]
  res_qtl$alleleB_new[res_qtl$swap] <- res_qtl$alleleA[res_qtl$swap]

  res_qtl$zscore[!res_qtl$swap] <- res_qtl$Zscore[!res_qtl$swap]
  res_qtl$zscore[res_qtl$swap] <- (-1)*res_qtl$Zscore[res_qtl$swap]

  res_qtl = res_qtl[,c("MarkerName_new", "rsid","zscore","pval","MarkerName_old", "swap", "Zscore","NrSamples")]


  #combine
  combdat = merge(res_gwas, res_qtl, by.x="MarkerName", by.y="MarkerName_new", suffixes=c(".pheno",".qtl"))
  combdat$rsid = combdat$rsid.pheno
  #combdat = combdat[!is.na(combdat$rsid.pheno),]

  #run coloc
  # coloc.res <- coloc.abf(dataset1=list(N=33601, pvalues=combdat$pval.pheno, type="cc", s=0.5, snp=combdat$rsid.pheno, MAF=combdat$AF),
  #                     dataset2=list(N=combdat$NrSamples, pvalues=combdat$pval.qtl, type="quant", snp=combdat$rsid.pheno, MAF=combdat$AF))
  coloc.res <- coloc.abf(p12=5e-6, dataset1=list(N=33601, zs=combdat$zscore.pheno, type="cc", s=0.5, snp=combdat$rsid.pheno, MAF=combdat$AF),
                  dataset2=list(N=combdat$NrSamples, zs=combdat$zscore.qtl, type="quant", snp=combdat$rsid.pheno, MAF=combdat$AF))
  sorted_results = coloc.res$results[order(coloc.res$results$SNP.PP.H4, decreasing=TRUE),]

  return(list(index_gwas=index_gwas, coloc_summary = coloc.res$summary, coloc_results = sorted_results, coloc_input = combdat, feature=gene))
  #return(list(index_gwas=index_gwas, res_gwas, res_qtl, combdat2))
  #return(list(index_gwas=index_gwas, res_gwas, res_qtl, combdat))

}
#getEQTL_ColocDat(index_gwas="chr6:90267049:G:A", gene="BACH2")

getLCPlot = function(coloc_out_obj, group, ylab=NULL) {

  index_gwas = coloc_out_obj[["index_gwas"]]
  #index_qtl = coloc_out_obj[["index_qtl"]]
  feature = coloc_out_obj[["feature"]]
  combdat = coloc_out_obj[["coloc_input"]]
  label = gsub(":","_",paste0(index_gwas,"__",feature,"__",group))
  cat(label,"\n")

  #create lcformat files
  if (!dir.exists("locuscompare")) {dir.create("locuscompare")}
  gwas_fn = paste0("locuscompare/lcformat_",label,"_gwas.txt")
  gwas_d = combdat[,c("rsid","pval.pheno")]
  names(gwas_d) <- c("rsid","pval")
  #cat("GWAS:",class(gwas_d),"\n")
  write.table(gwas_d, file=gwas_fn, col.names=TRUE, row.names=FALSE, quote=FALSE)

  qtl_fn = paste0("locuscompare/lcformat_",label,"_qtl.txt")
  qtl_d = combdat[,c("rsid","pval.qtl")]
  names(qtl_d) <- c("rsid","pval")
  #cat("QTL:",class(qtl_d),"\n")
  write.table(qtl_d, file=qtl_fn, col.names=TRUE, row.names=FALSE, quote=FALSE)

  #run locuscompare
  png(paste0("locuscompare/locuscompare_chromqtls_",label,".png"), width = 860, height = 480)
  p = locuscompare(in_fn1 = gwas_fn, in_fn2 = qtl_fn, title = 'GWAS', title2 = ylab)
  print(p)
  dev.off()

}

getPeakObject = function(peak) {
  vals = strsplit(peak, split="_")[[1]]
  out = list(label=peak, chr=vals[[1]], start=as.numeric(vals[[2]]), end=as.numeric(vals[[3]]))
  out[["prettyLabel"]] = paste0(out$chr, ":",out$start,"-",out$end)
  return(out)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run colocalisation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Due to lack of power in EUR only analysis we did the following:
# 1. defined index variants based on conditional analysis in EUR and AFR analyzed jointly (adjusting for global PCs projected onto 1000G)
#       - only 2/10 regions in this analysis had secondary signals at p<1e-3 (chr16_28914657_T_C and chr2_203687587_T_C)
# 2. for each of the 12 index variants found in (1), we fit EUR-specific models adjusting for EUR-specific PCs and conditioning on other index SNPs in the region (only for two regions with secondary signals)
# 3. for each index variant in LD with a GUESSFM index variant (rsq>0.5 in caQTL EUR data set), we formally tested for colocalisation using EUR-specifc caQTL stats or pooled AFR and EUR caQTL stats

index_pairs_un = getIndexPairs(gwas_index_variants=t1d_index_variants, qtl_index_variants=names(results_qtl_un), ld_dat=dat_un_EUR, threshold=0.25)
#index_pairs_un_all = getIndexPairs(gwas_index_variants=t1d_index_variants, qtl_index_variants=names(results_qtl_un), ld_dat=dat_un, threshold=0.25)

coloc_out_un = apply(index_pairs_un, 1, function(x) {getColocDat(index_gwas = x["snp1"], index_qtl = x["snp2"], dat_qtl = dat_un, results_qtl = results_qtl_un)} )
lapply(coloc_out_un, getLCPlot, group="un", ylab="caQTL")

coloc_out_EUR = apply(index_pairs_un, 1, function(x) {getColocDat(index_gwas = x["snp1"], index_qtl = x["snp2"], dat_qtl = dat_un_EUR, results_qtl = results_qtl_EUR)} )
lapply(coloc_out_EUR, getLCPlot, group="EUR", ylab="caQTL")

index_pairs_AFR = getIndexPairs(gwas_index_variants=t1d_index_variants, qtl_index_variants=names(results_qtl_AFR), ld_dat=dat_un_AFR)
# coloc_out_AFR = apply(index_pairs_AFR, 1, function(x) {getColocDat(index_gwas = x["snp1"], index_qtl = x["snp2"], dat_qtl = dat_un_AFR, results_qtl = results_qtl_AFR)} )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prep whole Blood eQTL data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library("liftOver")
chain_hg19_to_hg38 = import.chain(con=paste0(Sys.getenv("resources"),"/hg19ToHg38.over.chain"))

#read in genome-wide significant eqtls
wb_eqtl = readRDS(file=file.path(Sys.getenv("qtl"),"wholeblood_eQTL/cis-eQTLs_significant_FDR_lt_0.001.rds"))
credible_vars = unlist(strsplit(finemap[finemap$tag_rsid %in% index_pairs_un$rsid1, "rsid_merged"], split=";"))

#get significant eqtls overlapping T1D credible regions
wb_eqtl_overlap = wb_eqtl[wb_eqtl$SNP %in% credible_vars,]
saveRDS(wb_eqtl_overlap, file="cis-eQTLs_significant_FDR_lt_0.001_coloc_index_vars.rds")

#extract complete eqtl summary stats for each region
extractEqtlData = function(gene) {
  cat("extracting eQTL statistics for",gene,"\n")
  command = paste("sbatch",paste0("--output=extract_eqtl_results_by_gene_",gene,".log"),paste0(Sys.getenv("scripts"),"/fine_mapping/extract_eqtl_results_by_gene.slurm"), gene)
  cat(command,"\n")
  system(command)
}
genes = unique(wb_eqtl_overlap$GeneSymbol)
for ( i in 1:length(genes)) {
  genefile = paste0(Sys.getenv("PUBLIC_DATA"),"/eQTLGen/cis-eQTLs_full_20180905__",genes[i],".txt")
  if (!file.exists(genefile)) {
    cat("creating", genefile,"\n")
    extractEqtlData(genes[i])
  }
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Whole Blood eQTL colocalisation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
matchIndexToGenes = function(index_gwas) {
  cred_vars = unlist(strsplit(finemap[finemap$tag==index_gwas,"rsid_merged"], split=";"))
  unique(wb_eqtl_overlap[wb_eqtl_overlap$SNP%in%cred_vars,"GeneSymbol"])
}
matchIndexToGenes("chr12:56047884:TA:T")

runEQTL = function(index) {
  genes = matchIndexToGenes(index)
  cat(index,"-->",paste(genes,collapse=" "),"\n")
  results = list()
  for (gene in genes) {
    cat("Running colocalisation for",gene,"\n")
    results[[gene]] = getEQTL_ColocDat(index_gwas=index, gene=gene)
    cat(" \n")
  }
  return(results)
}

# run coloc
coloc_eqtl_out = lapply(gsub("_",":",index_pairs_un$snp1), runEQTL)
names(coloc_eqtl_out)<- index_pairs_un$snp1
saveRDS(coloc_eqtl_out, file="eqtl_colocalisation_stats.rds")

# generate lc plots
for ( i in 1:length(coloc_eqtl_out)) {
  lapply(coloc_eqtl_out[[i]], getLCPlot, group="EUR", ylab="eQTL")
}


summarizeEQTL = function(index) {
  cat(index,"\n")
  out = coloc_eqtl_out[[index]]
  genes = names(out)
  if (length(genes)>0) {
    index_reformat = gsub("_",":",index)
    if (index_reformat %in% out[[1]]$coloc_input$MarkerName) {
      zscores = sapply(genes, function(x){coloc_input=out[[x]]$coloc_input; coloc_input[coloc_input$MarkerName==index_reformat,"zscore.qtl"]})
      pvals = sapply(genes, function(x){coloc_input=out[[x]]$coloc_input; coloc_input[coloc_input$MarkerName==index_reformat,"pval.qtl"]})
    } else {
      zscores = NA
      pvals = NA
    }
    PP.H4.abfs = sapply(genes, function(x) {out[[x]]$coloc_summary["PP.H4.abf"]})
    df = data.frame(index=index, gene=genes, zscore=zscores, pvalue=pvals, PP.H4.abf=as.numeric(PP.H4.abfs))
    df = df[df$PP.H4.abf>0.8,]
    if (nrow(df)>0) {
      df$OUT = paste0(df$gene, " (z=", format(df$zscore, digit=2, trim=TRUE), "; P=", format(df$pvalue, digits=2, trim=TRUE), "; PP.H4.abf=", format(round(df$PP.H4.abf, digit=2), scientific=FALSE), ")")
      OUTSTRING = paste(df$OUT, collapse="\n")
    } else {
      OUTSTRING = NA
    }
  } else {
    OUTSTRING = NA
  }
  return(OUTSTRING)
}

# summarizeEQTL(index_pairs_un$snp1[1])
# summarizeEQTL("chr12_56047884_TA_T")
# summarizeEQTL("chr14_98032614_A_G")
# lapply(index_pairs_un$snp1, summarizeEQTL)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Collect colocalisation data into a table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## For each index pair, we want
# - gwas index variant
# - index variant alleleB freq
# - gwas assoc stats
# - peak
# - qtl index variant
# - ld with gwas index
# - combined caQTL stats: beta, pval, PP.H4.abf
# - EUR caQTL stats: beta, pval, PP.H4.abf
# - GUESSFM-prioritized variants in peak
# - eQTL stats
#     - gene; eqtl index variant; eqtl beta; eqtl pvalue

collectData = function(coloc_res, coloc_res_EUR) {
  coloc_input = coloc_res[["coloc_input"]]
  coloc_input_EUR = coloc_res_EUR[["coloc_input"]]

  index_gwas = gsub("_",":",coloc_res[["index_gwas"]])
  gwas_alleles = coloc_input[coloc_input$MarkerName==index_gwas,c("alleleA","alleleB")]
  gwas_AF = getEUR_AF(index_gwas)
  gwas_stats = coloc_input[coloc_input$MarkerName==index_gwas,c("beta.pheno","pval.pheno")]

  index_qtl = gsub("_",":",coloc_res[["index_qtl"]])
  index_ld = getLD(coloc_res[["index_gwas"]], coloc_res[["index_qtl"]], dat_sub=dat_un_EUR)["rsq"]
  qtl_stats =  c(coloc_input[coloc_input$MarkerName==index_qtl,c("beta.qtl","pval.qtl")], coloc_res[["coloc_summary"]]["PP.H4.abf"])
  qtl_stats_EUR =  c(coloc_input_EUR[coloc_input_EUR$MarkerName==index_qtl,c("beta.qtl","pval.qtl")], coloc_res_EUR[["coloc_summary"]]["PP.H4.abf"])
  peak = getPeakObject(coloc_res[["coloc_input"]]$peak[1])
  dat_snps_in_group = data.frame(finemap[finemap$tag==index_gwas,])
  snps_in_peak = dat_snps_in_group$ID[paste0("chr",dat_snps_in_group$chromosome)==peak$chr & dat_snps_in_group$position>=peak$start & dat_snps_in_group$position<=peak$end]
  eqtl_stats = summarizeEQTL(gsub(":","_",index_gwas))

  #deal with missingness
  #gwas_stats = lapply(gwas_stats, function(x) if (length(x)==0){NA}else{x}); names(gwas_stats) <- c("beta.pheno","pval.pheno")
  #eqtl_stats = lapply(eqtl_stats, function(x) if (length(x)==0){NA}else{x})
  qtl_stats = lapply(qtl_stats, function(x) if (length(x)==0){NA}else{x})
  qtl_stats_EUR = lapply(qtl_stats_EUR, function(x) if (length(x)==0){NA}else{x})

  OUT = data.frame(
    gwas_lead=paste(getRsid(index_gwas),index_gwas, sep=";"),
    gwas_alleles,
    AF_alleleB = gwas_AF,
    gwas_stats,
    eqtl_stats,
    peak = peak$prettyLabel,
    qtl_lead=paste(getRsid(index_qtl), index_qtl, sep=";"),
    qtl_ld = index_ld,
    qtl_stats,
    qtl_stats_EUR,
    snps_in_peak = paste(snps_in_peak, collapse=";")
  )
  return(OUT)
}


ivals = 1:nrow(index_pairs_un); #ivals = ivals[!index_pairs_un$snp1%in%c("chr12_56034460_G_A","chr16_75218429_G_A")]
#ivals = 1:(nrow(index_pairs_un)-1); #ivals = ivals[!(index_pairs_un$snp1 == "chr12_56034460_G_A" | index_pairs_un$snp2=="chr16_75213493_A_G")]
OUTTABLE = do.call("rbind",lapply(ivals, function(i) {collectData(coloc_out_un[[i]],coloc_out_EUR[[i]])}))
saveRDS(OUTTABLE, file="caqtl_colocalisation_stats.rds")

coloc_out_un_SUMMARY = data.frame(do.call("rbind",lapply(coloc_out_un, function(x) c(x[["index_gwas"]],x[["index_qtl"]],x[["coloc_summary"]]))))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save to table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source(paste0(Sys.getenv("scripts"),"/manuscript/format_tables.R"))
OUTTABLE_TOP = OUTTABLE[OUTTABLE$PP.H4.abf>0.75 | OUTTABLE$PP.H4.abf.1>0.75,]


library(rtf)
rtffile <- RTF("colocalisation_caqtl_tophits.doc")
addParagraph(rtffile, "Colocalisation between T1D association signals and caQTLs in CD4+ T cells")
addTable(rtffile, prettyTable(OUTTABLE_TOP))
addParagraph(rtffile, "Note - rs796916887 was merged into rs61555617 on October 12, 2018")
done(rtffile)

write.table(OUTTABLE_TOP,file="colocalisation_caqtl_tophits.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Sensitivity analysis for BACH2 association (removing homozygous alt individuals)
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# samples_rs729 = dat_un[["samples"]]
# geno_rs729 = dat_un[["genotypes"]][,"chr6_90267049_G_A"]
# counts_rs729 = dat_un[["counts"]]["chr6_90266736_90267780",]
# dat_rs729 = data.frame(samples_rs729, gt=geno_rs729,count=counts_rs729)
# summary(lm(data=dat_rs729, formula=count~gt+PC1+PC2))
# summary(lm(data=dat_rs729[dat_rs729$gt<2,], formula=count~gt+PC1+PC2))
#
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Look up for ANKRD55
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# f1_ank = f1[f1$avsnp150%in%c("rs71624119","rs7731626"),]
# f1_ank$Z_AFR = f1_ank$beta_AFR/f1_ank$se_AFR
# f1_ank$Z_EUR = f1_ank$beta_EUR/f1_ank$se_EUR
# f1_ank$Z_AMR = f1_ank$beta_AMR/f1_ank$se_AMR
# f1_ank$Z_FIN = f1_ank$beta_FIN/f1_ank$se_FIN
#
# #grep 'chr5:56144903:G:A\|chr5:56148856:G:A'
#
# #Variants ordered by genomic position
# chr5:56141024:T:A --> rs6873385
# chr5:56142753:C:A --> rs6859219
# chr5:56143024:C:T --> rs10065637
# chr5:56144903:G:A --> rs71624119
# chr5:56146422:T:C --> rs10213692
# chr5:56148856:G:A --> rs7731626
#
# finemap[!is.na(finemap$MarkerName) & finemap$MarkerName=="chr5:56142753:C:A","rsid_merged"]
# finemap[!is.na(finemap$MarkerName) & finemap$MarkerName=="chr5:56143024:C:T","rsid_merged"]
# finemap[!is.na(finemap$MarkerName) & finemap$MarkerName=="chr5:56148856:G:A","rsid_merged"]
