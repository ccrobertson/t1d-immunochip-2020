#Pi_dens.R

#trying Pi on the FDR associated variants:
#Hai Fang, Nature Genetics, 2019.

library(Pi)
library(GenomicRanges)
library(ggforce)

tmpdir<-"/well/todd/users/jinshaw/mega/"
#Read in the T1D summary stats:
hits<-read.table(file=paste0(tmpdir,"sup_t5.txt.gz"), header=T, as.is=T, sep=",")
hits<-hits[!is.na(hits$Effect_EUR),]

#remove MHC:
hits<-hits[!(hits$chromosome==6 & hits$position>25000000 & hits$position<35000000),]

#for now, keeping those in ImmunoChip regions, since these are likely the peak in the signal - no guarentee of this in regions outside iChip:
ichip<-read.table(file=paste0(tmpdir,"/dens_and_ichip_regs.txt"),header=T,as.is=T)

hitsr<-GRanges(seqnames=paste0("chr",as.character(hits$chromosome)),
IRanges(hits$position,end=hits$position))
ichips<-GRanges(seqnames=as.character(ichip$seqnames),
IRanges(ichip$start,end=ichip$end))

l<-subsetByOverlaps(hitsr, ichips)

hitsichip<-hits[hits$position %in% as.data.frame(l)$start,]
hits<-hitsichip
hits1<-hits[,c("ID","pphaseIandII","alleleB","alleleA","EffectphaseIandII","SEphaseIandII","MarkerName","AF_EUR")]
colnames(hits1)<-c("snp","p","effect","other","b","se","rsid","AF_EUR")

hits1<-hits1[hits1$snp!=".",]

#if duplicate variants exists, keep the most associated:
library(plyr)
library(dplyr)
hits1<-hits1 %>%
group_by(snp) %>%
filter(p==min(p))

hits1<-as.data.frame(hits1)

doitall<-function(name,eqtls,score.cap){

include.LD <- 'EUR'
LD.r2 <- 0.8
LD.customised <- NULL
significance.threshold <- 5e-05
score.cap=score.cap
distance.max <- 20000
decay.kernel <- "constant"
decay.exponent <- 2
eQTL.customised <- NULL
include.HiC <- c("Monocytes", "Fetal_thymus", "Naive_CD4_T_cells", 
        "Total_CD4_T_cells", "Activated_total_CD4_T_cells", "Nonactivated_total_CD4_T_cells", 
        "Naive_CD8_T_cells", "Total_CD8_T_cells", "Naive_B_cells", 
        "Total_B_cells")

include.eQTL<-eqtls

include.TAD="GM12878"
GR.SNP <- "dbSNP_GWAS"
GR.Gene <- "UCSC_knownGene"
cdf.function <- "empirical"
STRING_only<-NA
scoring.scheme <- 'max'
network <- "STRING_high"
network.customised <- NULL
weighted <- FALSE
normalise <- "laplacian"
normalise.affinity.matrix <- "none"
restart <- 0.7
parallel <- TRUE
multicores <- NULL
verbose <- TRUE
RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
seeds.inclusive <- TRUE
STRING.only <- NA

data<-hits1
data$suggestive=T
fdrhits<-read.table(file=paste0(tmp,"/genomewide_hits_univariable_fdr_dens.txt"),header=T,as.is=T,sep=",")
fdrhits$rsid<-fdrhits$MarkerName
data$suggestive<-ifelse(data$rsid %in% fdrhits$rsid,F,data$suggestive)

#however, LD function doesnt seem to work with indels or variants without RSID's. I am therefore going to 
#manually change the variants to define the regions to be SNPs and have RSIDs. 
#changing chr1:116738074:C:T to rs112360546
#changing rs376644292 to rs7566146
#changing rs78636954 to rs6534347
#changing rs140418154 to rs138650988
#changing rs375341506 to rs72793372

data$suggestive<-ifelse(data$snp=="chr1:116738074:C:T",T,
ifelse(data$snp=="rs376644292",T,
ifelse(data$snp=="rs78636954",T,
ifelse(data$snp=="rs140418154",T,
ifelse(data$snp=="rs375341506",T,
ifelse(data$snp=="rs112360546",F,
ifelse(data$snp=="rs7566146",F,
ifelse(data$snp=="rs6534347",F,
ifelse(data$snp=="rs138650988",F,
ifelse(data$snp=="rs72793372",F,data$suggestive))))))))))

data<-data[,-7]
data<-data[!is.na(data$AF_EUR),]

#edit functions to be able to include eqlgen data:
xPierSNPsAdvABF1<-function (data, include.LD = NA, LD.customised = NULL, LD.r2 = 0.8, 
    significance.threshold = 1e-05, score.cap = 10, distance.max = 2000, 
    decay.kernel = c("slow", "constant", "linear", "rapid"), 
    decay.exponent = 2, GR.SNP = c("dbSNP_GWAS", "dbSNP_Common", 
        "dbSNP_Single"), GR.Gene = c("UCSC_knownGene", "UCSC_knownCanonical"), 
    include.TAD = c("none", "GM12878", "IMR90", "MSC", "TRO", 
        "H1", "MES", "NPC"), include.eQTL = c("CD14", "LPS2", 
        "LPS24", "IFN", "Bcell", "NK", "Neutrophil", "CD4", "CD8", 
        "Blood", "Monocyte", "shared_CD14", "shared_LPS2", "shared_LPS24", 
        "shared_IFN"), include.HiC = NA, cdf.function = c("empirical", 
        "exponential"), scoring.scheme = c("max", "sum", "sequential"), 
    network = c("STRING_highest", "STRING_high", "STRING_medium", 
        "STRING_low", "PCommonsUN_high", "PCommonsUN_medium", 
        "PCommonsDN_high", "PCommonsDN_medium", "PCommonsDN_Reactome", 
        "PCommonsDN_KEGG", "PCommonsDN_HumanCyc", "PCommonsDN_PID", 
        "PCommonsDN_PANTHER", "PCommonsDN_ReconX", "PCommonsDN_TRANSFAC", 
        "PCommonsDN_PhosphoSite", "PCommonsDN_CTD", "KEGG", "KEGG_metabolism", 
        "KEGG_genetic", "KEGG_environmental", "KEGG_cellular", 
        "KEGG_organismal", "KEGG_disease", "REACTOME"), STRING.only = c(NA, 
        "neighborhood_score", "fusion_score", "cooccurence_score", 
        "coexpression_score", "experimental_score", "database_score", 
        "textmining_score")[1], weighted = FALSE, network.customised = NULL, 
    seeds.inclusive = TRUE, normalise = c("laplacian", "row", 
        "column", "none"), restart = 0.7, normalise.affinity.matrix = c("none", 
        "quantile"), parallel = TRUE, multicores = NULL, verbose = TRUE, 
    verbose.details = FALSE, RData.location = "http://galahad.well.ox.ac.uk/bigdata", 
    ...) 
{
    startT <- Sys.time()
    if (verbose) {
        message(paste(c("Start at ", as.character(startT)), collapse = ""), 
            appendLF = TRUE)
        message("", appendLF = TRUE)
    }
    decay.kernel <- match.arg(decay.kernel)
    cdf.function <- match.arg(cdf.function)
    scoring.scheme <- match.arg(scoring.scheme)
    network <- match.arg(network)
    normalise <- match.arg(normalise)
    normalise.affinity.matrix <- match.arg(normalise.affinity.matrix)
    suggestive <- NULL
    data_significant <- subset(data, !suggestive)
    if (verbose == FALSE) {
        verbose.details <- FALSE
    }
    if (verbose) {
        now <- Sys.time()
        message(sprintf("Preparing the distance predictor (%s) ...", 
            as.character(now)), appendLF = TRUE)
    }
    relative.importance <- c(1, 0, 0)
    pNode_distance <- xPierSNPs(data = data_significant, include.LD = include.LD, 
        LD.customised = LD.customised, LD.r2 = LD.r2, significance.threshold = significance.threshold, 
        score.cap = score.cap, distance.max = distance.max, decay.kernel = decay.kernel, 
        decay.exponent = decay.exponent, GR.SNP = GR.SNP, GR.Gene = GR.Gene, 
        include.TAD = include.TAD, include.eQTL = NA, eQTL.customised = NULL, 
        include.HiC = NA, cdf.function = cdf.function, relative.importance = relative.importance, 
        scoring.scheme = scoring.scheme, network = network, weighted = weighted, 
        network.customised = network.customised, seeds.inclusive = seeds.inclusive, 
        normalise = normalise, restart = restart, normalise.affinity.matrix = normalise.affinity.matrix, 
        parallel = parallel, multicores = multicores, verbose = verbose.details, 
        RData.location = RData.location)
    ls_pNode_distance <- list(pNode_distance)
    names(ls_pNode_distance) <- paste("nGene_", distance.max, 
        "_", decay.kernel, sep = "")
    ls_pNode_eQTL <- NULL
    include.eQTLs <- include.eQTL[!is.na(include.eQTL)]
    if (length(include.eQTLs) > 0) {
        names(include.eQTLs) <- include.eQTLs
        ls_pNode_eQTL <- lapply(include.eQTLs, function(x) {
            if (verbose) {
                startT <- Sys.time()
                message(sprintf("Preparing the eQTL predictor '%s' through ABF (%s) ...", 
                  x, as.character(startT)), appendLF = TRUE)
            }
            suppressMessages(pNode <- xPierABF1(data = data, eqtl = x, 
                network = network, STRING.only = STRING.only, 
                weighted = weighted, network.customised = network.customised, 
                seeds.inclusive = seeds.inclusive, normalise = normalise, 
                restart = restart, normalise.affinity.matrix = normalise.affinity.matrix, 
                
parallel = parallel, multicores = multicores, 
                verbose = verbose, RData.location = RData.location,cutoff.pgwas=1e-04, 
                ...))
            if (verbose & is.null(pNode)) {
                message(sprintf("\tNote: this predictor '%s' is NULL", 
                  x), appendLF = TRUE)
            }
            if (verbose) {
                endT <- Sys.time()
                runTime <- as.numeric(difftime(strptime(endT, 
                  "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), 
                  units = "secs"))
                message(paste(c("Runtime in total is: ", runTime, 
                  " secs\n"), collapse = ""), appendLF = TRUE)
            }
            return(pNode)
        })
        names(ls_pNode_eQTL) <- paste("eGene_", names(ls_pNode_eQTL), 
            sep = "")
    }
    include.HiCs <- include.HiC[!is.na(include.HiC)]
    if (length(include.HiCs) > 0) {
        names(include.HiCs) <- include.HiCs
        ls_pNode_HiC <- lapply(include.HiCs, function(x) {
            if (verbose) {
                message(sprintf("Preparing the HiC predictor '%s' (%s) ...", 
                  x, as.character(Sys.time())), appendLF = TRUE)
            }
            relative.importance <- c(0, 0, 1)
            pNode <- xPierSNPs(data = data_significant, include.LD = include.LD, 
                LD.customised = LD.customised, LD.r2 = LD.r2, 
                significance.threshold = significance.threshold, 
                score.cap = score.cap, distance.max = distance.max, 
                decay.kernel = decay.kernel, decay.exponent = decay.exponent, 
                GR.SNP = GR.SNP, GR.Gene = GR.Gene, include.eQTL = NA, 
                eQTL.customised = NULL, include.HiC = x, cdf.function = cdf.function, 
                relative.importance = relative.importance, scoring.scheme = scoring.scheme, 
                network = network, weighted = weighted, network.customised = network.customised, 
                seeds.inclusive = seeds.inclusive, normalise = normalise, 
                restart = restart, normalise.affinity.matrix = normalise.affinity.matrix, 
                parallel = parallel, multicores = multicores, 
                verbose = verbose.details, RData.location = RData.location)
            if (verbose & is.null(pNode)) {
                message(sprintf("\tNote: this predictor '%s' has NULL", 
                  x), appendLF = TRUE)
            }
            return(pNode)
        })
        names(ls_pNode_HiC) <- paste("cGene_", names(ls_pNode_HiC), 
            sep = "")
    }
    else {
        ls_pNode_HiC <- NULL
    }
    ls_pNode <- c(ls_pNode_distance, ls_pNode_eQTL, ls_pNode_HiC)
    endT <- Sys.time()
    if (verbose) {
        message(paste(c("\nFinish at ", as.character(endT)), 
            collapse = ""), appendLF = TRUE)
    }
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), 
        strptime(startT, "%Y-%m-%d %H:%M:%S"), units = "secs"))
    message(paste(c("Runtime in total is: ", runTime, " secs\n"), 
        collapse = ""), appendLF = TRUE)
    invisible(ls_pNode)
}



approx.bf.z <- function(z,f,type, N, s, suffix=NULL) {
  if(type=="quant") {
    sd.prior <- 0.15
    V <- Var.data(f, N)
  } else {
    sd.prior <- 0.2
    V <- Var.data.cc(f, N, s)
  }
  ## Shrinkage factor: ratio of the prior variance to the total variance
  r <- sd.prior^2 / (sd.prior^2 + V)
  ## Approximate BF  # I want ln scale to compare in log natural scale with LR diff
  lABF = 0.5 * (log(1-r) + (r * z^2))
  ret <- data.frame(V,z,r,lABF)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep=".")
  return(ret)
}

logdiff <- function(x,y) {
  my.max <- max(x,y)                              ##take out the maximum value in log form
  my.res <- my.max + log(exp(x - my.max ) - exp(y-my.max))
  return(my.res)
}

logsum <- function(x) {
  my.max <- max(x)                              ##take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}

Var.data.cc <- function(f, N, s) {
  1 / (2 * N * f * (1 - f) * s * (1 - s))
}

Var.data <- function(f, N) {
  1 / (2 * N * f * (1 - f))
}

sdY.est <- function(vbeta, maf, n) {
    warning("estimating sdY from maf and varbeta, please directly supply sdY if known")
    oneover <- 1/vbeta
    nvx <- 2 * n * maf * (1-maf)
    m <- lm(nvx ~ oneover - 1)
    cf <- coef(m)[['oneover']]
    if(cf < 0)
	stop("estimated sdY is negative - this can happen with small datasets, or those with errors.  A reasonable estimate of sdY is required to continue.")
    return(sqrt(cf))
}

combine.abf <- function(l1, l2, p1, p2, p12) {
  lsum <- l1 + l2
  lH0.abf <- 0
  lH1.abf <- log(p1) + logsum(l1)
  lH2.abf <- log(p2) + logsum(l2)
  lH3.abf <- log(p1) + log(p2) + logdiff(logsum(l1) + logsum(l2), logsum(lsum))
  lH4.abf <- log(p12) + logsum(lsum)

  all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3.abf, lH4.abf)
  my.denom.log.abf <- logsum(all.abf)
  pp.abf <- exp(all.abf - my.denom.log.abf)
  names(pp.abf) <- paste("PP.H", (1:length(pp.abf)) - 1, ".abf", sep = "")
  print(signif(pp.abf,3))
  print(paste("PP abf for shared variant: ", signif(pp.abf["PP.H4.abf"],3)*100 , '%', sep=''))
  return(pp.abf)
}

approx.bf.estimates <- function (z, V, type, suffix=NULL, sdY=1) {
  sd.prior <- if (type == "quant") { 0.15*sdY } else { 0.2 }
  r <- sd.prior^2/(sd.prior^2 + V)
  lABF = 0.5 * (log(1 - r) + (r * z^2))
  ret <- data.frame(V, z, r, lABF)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep = ".")
  return(ret)
}

process.dataset <- function(d, suffix) {
  #message('Processing dataset')

  nd <- names(d)
  if (! 'type' %in% nd)
    stop("dataset ",suffix,": ",'The variable type must be set, otherwise the Bayes factors cannot be computed')

  if(!(d$type[1] %in% c("quant","cc")))
      stop("dataset ",suffix,": ","type must be quant or cc")

  if(d$type[1]=="cc") {
      if(! "s" %in% nd)
          stop("dataset ",suffix,": ","please give s, proportion of samples who are cases")
      if("pvalues" %in% nd && !( "MAF" %in% nd))
          stop("dataset ",suffix,": ","please give MAF if using p values")
      if(d$s<=0 || d$s>=1)
          stop("dataset ",suffix,": ","s must be between 0 and 1")
  }

  if(d$type[1]=="quant") {
      if(!("sdY" %in% nd || ("MAF" %in% nd && "N" %in% nd )))
          stop("dataset ",suffix,": ","must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated")
  }

  if("beta" %in% nd && "varbeta" %in% nd) {  ## use beta/varbeta.  sdY should be estimated by now for quant
    if(length(d$beta) != length(d$varbeta))
      stop("dataset ",suffix,": ","Length of the beta vectors and variance vectors must match")
    if(!("snp" %in% nd))
      d$snp <- sprintf("SNP.%s",1:length(d$beta))
    if(length(d$snp) != length(d$beta))
      stop("dataset ",suffix,": ","Length of snp names and beta vectors must match")

    if(d$type[1]=="quant" && !('sdY' %in% nd))
          d$sdY <- sdY.est(d$varbeta, d$MAF, d$N)
   df <- approx.bf.estimates(z=d$beta/sqrt(d$varbeta),
                              V=d$varbeta, type=d$type, suffix=suffix, sdY=d$sdY)
    df$snp <- as.character(d$snp)
    return(df)
}
  if("z" %in% nd & "MAF" %in% nd & "N" %in% nd) { ## no beta/varbeta: use p value / MAF approximation
    if (length(d$z) != length(d$MAF))
      stop('Length of the z-value vectors and MAF vector must match')
    if(!("snp" %in% nd))
      d$snp <- sprintf("SNP.%s",1:length(d$z))
    df <- data.frame(z = d$z,
                     MAF = d$MAF,
                     snp=as.character(d$snp))
    colnames(df)[-3] <- paste(colnames(df)[-3], suffix, sep=".")
#   df <- subset(df, df$MAF>0 & df$z>0) # all p values and MAF > 0
    abf <- approx.bf.z(z=df$z, f=df$MAF, type=d$type, N=d$N, s=d$s, suffix=suffix)
    df <- cbind(df, abf)
    return(df)
  }
  else{
  stop("Must give, as a minimum, one of:\n(beta, varbeta, type, sdY)\n(beta, varbeta, type, MAF)\n(z, MAF, N, type)")
}
}


coloc.abf <- function(dataset1, dataset2, MAF=NULL,
                      p1=1e-4, p2=1e-4, p12=1e-5) {

  if(!is.list(dataset1) || !is.list(dataset2))
    stop("dataset1 and dataset2 must be lists.")
  if(!("MAF" %in% names(dataset1)) & !is.null(MAF))
    dataset1$MAF <- MAF
  if(!("MAF" %in% names(dataset2)) & !is.null(MAF))
    dataset2$MAF <- MAF

  df1 <- process.dataset(d=dataset1, suffix="df1")
  df2 <- process.dataset(d=dataset2, suffix="df2")
  merged.df <- merge(df1,df2)

   if(!nrow(merged.df))
    stop("dataset1 and dataset2 should contain the same snps in the same order, or should contain snp names through which the common snps can be identified")

  merged.df$internal.sum.lABF <- with(merged.df, lABF.df1 + lABF.df2)
  ## add SNP.PP.H4 - post prob that each SNP is THE causal variant for a shared signal
  my.denom.log.abf <- logsum(merged.df$internal.sum.lABF)
  merged.df$SNP.PP.H4 <- exp(merged.df$internal.sum.lABF - my.denom.log.abf)


##############################

  pp.abf <- combine.abf(merged.df$lABF.df1, merged.df$lABF.df2, p1, p2, p12)
  common.snps <- nrow(merged.df)
  results <- c(nsnps=common.snps, pp.abf)

  output<-list(summary=results, results=merged.df)
  return(output)
}








xPierABF1<-function (data, eqtl = c("CD14", "LPS2", "LPS24", "IFN", "Bcell", 
    "NK", "Neutrophil", "CD4", "CD8", "Blood", "Monocyte", "shared_CD14", 
    "shared_LPS2", "shared_LPS24", "shared_IFN","eQTLGen"), prior.eqtl = 1e-04, 
    prior.gwas = 1e-04, prior.both = 1e-05, cutoff.H4 = 0.8, 
    cutoff.pgwas = 1e-05, network = c("STRING_highest", "STRING_high", 
        "STRING_medium", "STRING_low", "PCommonsUN_high", "PCommonsUN_medium", 
        "PCommonsDN_high", "PCommonsDN_medium", "PCommonsDN_Reactome", 
        "PCommonsDN_KEGG", "PCommonsDN_HumanCyc", "PCommonsDN_PID", 
        "PCommonsDN_PANTHER", "PCommonsDN_ReconX", "PCommonsDN_TRANSFAC", 
        "PCommonsDN_PhosphoSite", "PCommonsDN_CTD", "KEGG", "KEGG_metabolism", 
        "KEGG_genetic", "KEGG_environmental", "KEGG_cellular", 
        "KEGG_organismal", "KEGG_disease"), STRING.only = c(NA, 
        "neighborhood_score", "fusion_score", "cooccurence_score", 
        "coexpression_score", "experimental_score", "database_score", 
        "textmining_score")[1], weighted = FALSE, network.customised = NULL, 
    seeds.inclusive = TRUE, normalise = c("laplacian", "row", 
        "column", "none"), restart = 0.7, normalise.affinity.matrix = c("none", 
        "quantile"), parallel = TRUE, multicores = NULL, verbose = TRUE, 
    RData.location = "http://galahad.well.ox.ac.uk/bigdata") 
{
    startT <- Sys.time()
    if (verbose) {
        message(paste(c("Start at ", as.character(startT)), collapse = ""), 
            appendLF = TRUE)
        message("", appendLF = TRUE)
    }
    eqtl <- match.arg(eqtl)
    network <- match.arg(network)
    normalise <- match.arg(normalise)
    normalise.affinity.matrix <- match.arg(normalise.affinity.matrix)
    pNode <- NULL
    se <- b <- effect <- other <- p <- NULL
    pp_ABF <- p_GWAS <- p_eQTL <- H4 <- NULL
    if (class(data) == "data.frame") {
        if (all(c("snp", "effect", "other", "b", "p", "se") %in% 
            colnames(data))) {
            summary_gwas <- data %>% dplyr::filter(!is.na(se), 
                se != 0, b != 0, effect != "", other != "") %>% 
                dplyr::arrange(p)
            summary_gwas <- summary_gwas[!duplicated(summary_gwas$snp), 
                ]
        }
        else {
            warnings("The input data.frame does not contain required columns: c('snp','effect','other','b','p','se').\n")
            return(NULL)
        }
    }
    vec_N_eqtl <- c(414, 261, 322, 367, 286, 245, 101, 293, 283, 
        5311, 287, 228, 228, 228, 228)
    names(vec_N_eqtl) <- c("CD14", "LPS2", "LPS24", "IFN", "Bcell", 
        "NK", "Neutrophil", "CD4", "CD8", "Blood", "Monocyte", 
        "shared_CD14", "shared_LPS2", "shared_LPS24", "shared_IFN")
    N_eqtl <- vec_N_eqtl[eqtl]
    if (0) {
        JK_cohort_xMEdb <- xRDataLoader("JK_cohort_xMEdb", verbose = F, 
            RData.location = RData.location)
        ind <- match(JK_cohort_xMEdb$context, eqtl)
        summary_eqtl <- JK_cohort_xMEdb[!is.na(ind), ]
    }
   if (eqtl=="eQTLGen"){
        #summary_eqtl<-read.table(file="/well/todd/users/jinshaw/eqtl/eqtlgen/cis-eQTLs_full_20180905.txt.gz",header=T,as.is=T)
        #summary_eqtl<-summary_eqtl[summary_eqtl$Pvalue<0.0005,]
	#save(summary_eqtl,file="/well/todd/users/jinshaw/eqtl/eqtlgen/cis-eQTLs_assoc.RData")
        load(file="/well/todd/users/jinshaw/eqtl/eqtlgen/cis-eQTLs_assoc.RData")
	summary_eqtl$snps<-summary_eqtl$SNP
	summary_eqtl$context<-"Whole blood"
	summary_eqtl$gene<-summary_eqtl$Gene
	summary_eqtl$Symbol<-summary_eqtl$GeneSymbol
	summary_eqtl$effect_allele<-summary_eqtl$AssessedAllele
	summary_eqtl$other_allele<-summary_eqtl$OtherAllele
	summary_eqtl$N<-summary_eqtl$NrSamples        
	summary_eqtl$mode<-"cis"
	summary_eqtl$ProbeID<-summary_eqtl$Gene
        summary_eqtl$snp_cse<-paste0("chr",summary_eqtl$SNPChr,":",summary_eqtl$SNPPos,"-",summary_eqtl$SNPPos)
	summary_eqtl$gene_cse<-paste0("chr",summary_eqtl$GeneChr,":",summary_eqtl$GenePos,"-",summary_eqtl$GenePos,":-")
	summary_eqtl$pvalue<-summary_eqtl$Pvalue
}
    else {
        summary_eqtl <- xRDataLoader(paste0("JK_cohort_xMEdb_", 
            eqtl), verbose = F, RData.location = RData.location)
    }
    ind <- match(summary_gwas$snp, summary_eqtl$snps)
    summary_gwas <- summary_gwas[!is.na(ind), ]
    ls_gene <- split(x = summary_eqtl, f = summary_eqtl$gene)
    ls_df_output <- lapply(1:length(ls_gene), function(i) {
        if (verbose) {
            if (i%%1000 == 0) {
                message(sprintf("Analysing %d (%d) (%s)", i, 
                  length(ls_gene), as.character(Sys.time())), 
                  appendLF = T)
            }
            message(sprintf("Analysing %d (%d) (%s)", i, length(ls_gene), 
                as.character(Sys.time())), appendLF = T)
        }
        df_output <- NULL
        df_eqtl <- ls_gene[[i]]
        ind <- match(summary_gwas$snp, df_eqtl$snps)
        if (sum(!is.na(ind)) >= 3) {
            df_gwas <- summary_gwas[!is.na(ind), ]
            df_eqtl <- df_eqtl[ind[!is.na(ind)], ]
            df_gwas$b_corrected <- df_gwas$b
            ind <- which(df_gwas$effect == df_eqtl$other_allele & 
                df_gwas$other == df_eqtl$effect_allele)
            df_gwas$b_corrected[ind] <- -1 * df_gwas$b[ind]
            ind <- df_gwas$effect != df_eqtl$effect_allele & 
                df_gwas$effect != df_eqtl$other_allele
            if (sum(ind) > 0) {
                df_gwas <- df_gwas[!ind, ]
                df_eqtl <- df_eqtl[!ind, ]
            }
            if (nrow(df_gwas) >= 3) {
		if(eqtl!="eQTLGen"){
                eqtl.summary <- list(beta = df_eqtl$beta, varbeta = df_eqtl$se^2, 
                  N = N_eqtl, MAF = df_eqtl$effect_maf, type = "quant", 
                  snp = df_eqtl$snps)
                gwas.summary <- list(beta = df_gwas$b_corrected, 
                  varbeta = df_gwas$se^2, type = "cc", snp = df_gwas$snp)
                res <- XGR::xMEabf(eqtl.summary, gwas.summary, 
                  prior.eqtl = prior.eqtl, prior.gwas = prior.gwas, 
                  prior.both = prior.both)
}
		if(eqtl=="eQTLGen"){
			df_eqtl$type="quant"
			df_eqtl$MAF<-ifelse(df_gwas$AF_EUR>0.5,1-df_gwas$AF_EUR,df_gwas$AF_EUR)
			df_eqtl$z<-df_eqtl$Zscore
			df_eqtl$snp<-df_eqtl$snps
			t<-process.dataset(df_eqtl,suffix="quant")
			df_gwas$s=0.4
			df_gwas$type="cc"
			df_gwas$z<-df_gwas$b_corrected/df_gwas$se
			df_gwas$MAF<-ifelse(df_gwas$AF_EUR>0.5,1-df_gwas$AF_EUR,df_gwas$AF_EUR)
			df_gwas$N<-59974
			t1<-process.dataset(df_gwas,suffix="cc")
		  
 		res <- coloc.abf(df_eqtl, df_gwas)
}
                df_res <- res$results
                ind <- match(df_res$snp, df_eqtl$snps)
                df_res <- cbind(df_res, df_eqtl[ind, ], df_gwas[ind, 
                  ])
if(eqtl!="eQTLGen"){
                df_output <- df_res[, c("context", "mode", "ProbeID", 
                  "Symbol", "gene_cse", "snps", "snp_cse", "effect_allele", 
                  "other_allele", "b_corrected", "beta", "p", 
                  "pvalue", "SNP.PP.H4")]
                colnames(df_output) <- c("context", "mode", "ProbeID", 
                  "Symbol", "gene_cse", "snps", "snp_cse", "A1", 
                  "A2", "b_GWAS", "b_eQTL", "p_GWAS", "p_eQTL", 
                  "pp_ABF")
#                df_output$b_ABF <- df_output$b_GWAS/df_output$b_eQTL
                df_output <- df_output %>% dplyr::arrange(-pp_ABF, 
                  p_GWAS, p_eQTL)
}
if(eqtl=="eQTLGen"){
		df_output <- df_res[, c("context", "mode", "ProbeID",
                  "Symbol", "gene_cse", "snps", "snp_cse", "effect_allele",
                  "other_allele", "b_corrected", "z", "p",
                  "pvalue", "SNP.PP.H4")]
                colnames(df_output) <- c("context", "mode", "ProbeID",
                  "Symbol", "gene_cse", "snps", "snp_cse", "A1",
                  "A2", "b_GWAS", "z_eQTL", "p_GWAS", "p_eQTL",
                  "pp_ABF")
#                df_output$b_ABF <- df_output$b_GWAS/df_output$b_eQTL
                df_output <- df_output %>% dplyr::arrange(-pp_ABF,
                  p_GWAS, p_eQTL)
}

                df_output$nsnps <- res$summary[1]
                df_output$H0 <- res$summary[2]
                df_output$H1 <- res$summary[3]
                df_output$H2 <- res$summary[4]
                df_output$H3 <- res$summary[5]
                df_output$H4 <- res$summary[6]
            }
        }
        df_output
    })
    df_output <- do.call(rbind, ls_df_output)
    ind <- which(!duplicated(df_output$ProbeID))
    df <- df_output[ind, ]
    df <- subset(df, H4 > cutoff.H4 & p_GWAS < cutoff.pgwas)
    df_evidence <- df %>% dplyr::arrange(-H4, p_GWAS, p_eQTL)
    df <- df %>% dplyr::arrange(-H4, p_GWAS, p_eQTL)
    ind <- which(!duplicated(df$Symbol))
    df <- df[ind, ]
    data_subset <- df[, c("Symbol", "pp_ABF")]
    if (nrow(data_subset) != 0) {
        pNode <- suppressMessages(xPierGenes(data = data_subset, 
            network = network, STRING.only = STRING.only, weighted = weighted, 
            network.customised = network.customised, seeds.inclusive = seeds.inclusive, 
            normalise = normalise, restart = restart, normalise.affinity.matrix = normalise.affinity.matrix, 
            parallel = parallel, multicores = multicores, verbose = verbose, 
            RData.location = RData.location))
        pNode$evidence <- df_evidence
    }
    endT <- Sys.time()
    if (verbose) {
        message(paste(c("\nFinish at ", as.character(endT)), 
            collapse = ""), appendLF = TRUE)
    }
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), 
        strptime(startT, "%Y-%m-%d %H:%M:%S"), units = "secs"))
    message(paste(c("Runtime in total is: ", runTime, " secs\n"), 
        collapse = ""), appendLF = TRUE)
    invisible(pNode)
}




ls_pNode_genomic <- xPierSNPsAdvABF1(data, include.LD=include.LD, LD.customised=LD.customised, LD.r2=LD.r2, 
significance.threshold=significance.threshold, score.cap=score.cap, distance.max=distance.max, decay.kernel=decay.kernel, 
decay.exponent=decay.exponent, GR.SNP=GR.SNP, GR.Gene=GR.Gene, include.TAD=include.TAD, include.eQTL=include.eQTL, 
include.HiC=include.HiC, cdf.function=cdf.function, scoring.scheme=scoring.scheme, network=network, STRING.only=STRING.only, 
weighted=weighted, network.customised=network.customised, seeds.inclusive=seeds.inclusive, normalise=normalise, restart=restart, 
normalise.affinity.matrix=normalise.affinity.matrix, parallel=parallel, multicores=multicores, verbose=verbose, RData.location=RData.location)




data.file <- file.path(RData.location, "iAnno.txt")
iA <- read.delim(data.file, header=TRUE, stringsAsFactors=FALSE)[,c("Symbol","OMIM","Phenotype","Function")]
colnames(iA) <- c("Symbol","dGene","pGene","fGene")
ls_pNode_anno <- lapply(2:4, function(j){
    data_anno <- subset(data.frame(seed=iA$Symbol,weight=iA[,j],stringsAsFactors=F), weight>0)
    pNode <- xPierAnno(data_anno, list_pNode=ls_pNode_genomic, network=network, STRING.only=STRING.only, 
weighted=weighted, network.customised=network.customised, seeds.inclusive=seeds.inclusive, 
normalise=normalise, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix, 
parallel=parallel, multicores=multicores, verbose=verbose, RData.location=RData.location)
})
names(ls_pNode_anno) <- colnames(iA)[2:4]
## bring together both predictors
ls_pNode <- c(ls_pNode_anno, ls_pNode_genomic)

dTarget <- xPierMatrix(ls_pNode, displayBy="pvalue", aggregateBy="fishers", RData.location=RData.location)

save(dTarget,file=paste0(tmpdir,"/Pi/",name,".RData"))
}

#doitall(name="noeqtlgen_score_10",eqtls=c("CD14", "LPS2", "LPS24", "IFN", "Bcell", 
#        "NK", "Neutrophil", "CD4", "CD8", "Blood", "Monocyte"),score.cap=10)
#doitall(name="noeqtlgen_score_20_dens",eqtls=c("CD14", "LPS2", "LPS24", "IFN", "Bcell",
#        "NK", "Neutrophil", "CD4", "CD8", "Blood", "Monocyte"),score.cap=20)
#doitall(name="noeqtlgen_score_30",eqtls=c("CD14", "LPS2", "LPS24", "IFN", "Bcell",
#        "NK", "Neutrophil", "CD4", "CD8", "Blood", "Monocyte"),score.cap=30)
#doitall(name="noeqtlgen_score_40",eqtls=c("CD14", "LPS2", "LPS24", "IFN", "Bcell",
#        "NK", "Neutrophil", "CD4", "CD8", "Blood", "Monocyte"),score.cap=40)
#doitall(name="noeqtlgen_score_50",eqtls=c("CD14", "LPS2", "LPS24", "IFN", "Bcell",
#        "NK", "Neutrophil", "CD4", "CD8", "Blood", "Monocyte"),score.cap=50)


#doitall(name="eqtlgen_score_10",eqtls=c("CD14", "LPS2", "LPS24", "IFN", "Bcell",
#        "NK", "Neutrophil", "CD4", "CD8", "Blood", "Monocyte","eQTLGen"),score.cap=10)
doitall(name="eqtlgen_score_20_dens",eqtls=c("CD14", "LPS2", "LPS24", "IFN", "Bcell",
        "NK", "Neutrophil", "CD4", "CD8", "Blood", "Monocyte","eQTLGen"),score.cap=20)
#doitall(name="eqtlgen_score_30",eqtls=c("CD14", "LPS2", "LPS24", "IFN", "Bcell",
#        "NK", "Neutrophil", "CD4", "CD8", "Blood", "Monocyte","eQTLGen"),score.cap=30)
#doitall(name="eqtlgen_score_40",eqtls=c("CD14", "LPS2", "LPS24", "IFN", "Bcell",
#        "NK", "Neutrophil", "CD4", "CD8", "Blood", "Monocyte","eQTLGen"),score.cap=40)
#doitall(name="eqtlgen_score_50",eqtls=c("CD14", "LPS2", "LPS24", "IFN", "Bcell",
#        "NK", "Neutrophil", "CD4", "CD8", "Blood", "Monocyte","eQTLGen"),score.cap=50)


#doitall(name="eqtlgen_cd4_cd8_score_10",eqtls=c("CD4", "CD8","eQTLGen"),score.cap=10)
#doitall(name="eqtlgen_cd4_cd8_score_20",eqtls=c("CD4", "CD8","eQTLGen"),score.cap=20)
#doitall(name="eqtlgen_cd4_cd8_score_30",eqtls=c("CD4", "CD8","eQTLGen"),score.cap=30)
#doitall(name="eqtlgen_cd4_cd8_score_40",eqtls=c("CD4", "CD8","eQTLGen"),score.cap=40)
#doitall(name="eqtlgen_cd4_cd8_score_50",eqtls=c("CD4", "CD8","eQTLGen"),score.cap=50)
