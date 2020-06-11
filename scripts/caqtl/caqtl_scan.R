### Setup

#NOTE: must use github version of MatrixEQTL.
#CRAN version throws an error when using noFDRsaveMemory = TRUE
#see https://github.com/andreyshabalin/MatrixEQTL/issues/8 for details
#install_github("andreyshabalin/MatrixEQTL")
library(MatrixEQTL)

setwd(Sys.getenv("freeze"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		Read in data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat_un = readRDS("DATASET_un.rds")
dat_stim = readRDS("DATASET_stim.rds")
dat_un_EUR = readRDS("DATASET_un_EUR.rds")
dat_un_AFR = readRDS("DATASET_un_AFR.rds")
dat_test = readRDS("DATASET_test.rds")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		Run MatrixEQTL
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
runQTLScan = function(dat_sub, outfile_cis, outfile_trans, pvOutputThreshold_cis, pvOutputThreshold_trans, covariates) {

  genodat = dat_sub[["dosages"]]
  countMat = dat_sub[["counts"]]
  sampleDF = data.frame(dat_sub[["samples"]], dat_sub[["pheno"]], dat_sub[["peer"]])
  sampleDF$selection = sampleDF$SQB %in% c("ATAC_UF_SQB101", "ATAC_UF_SQB102", "ATAC_UF_SQB103", "ATAC_UF_SQB104", "ATAC_UF_SQB105", "ATAC_UF_SQB106")


	#set up input matrices
	X = t(genodat)
	COV=t(sampleDF[,covariates])
	#SNPlength = sapply(row.names(X), function(x){ length(unlist(strsplit(x, split="_")))})
	SNPmat = t(sapply(row.names(X), function(x){ unlist(strsplit(x, split="_"))[c(1,2)]}))
	SNP = data.frame(snpid=row.names(SNPmat), chr=as.numeric(gsub("chr","",SNPmat[,1])), pos=as.numeric(SNPmat[,2]))
	getRegion = function(x){
		s = strsplit(x,split="_")[[1]]
		n = length(s)
		c(s[1], s[(n-1)],s[n])
	}
	REGmat = do.call("rbind",lapply(row.names(countMat),getRegion))
	REG = data.frame(geneid=row.names(countMat), chr=as.numeric(gsub("chr","",REGmat[,1])), left=as.numeric(REGmat[,2]),  right=as.numeric(REGmat[,3]))

	## Load genotype data
	snps = SlicedData$new();
	snps$CreateFromMatrix(X);

	## Load gene expression data
	gene = SlicedData$new();
	gene$CreateFromMatrix(countMat);

	## Load covariates
	cvrt = SlicedData$new();
	cvrt$CreateFromMatrix(COV);

	## Run the analysis
	#run all: pvOutputThreshold > 0, pvOutputThreshold.cis = 0
	#local only: pvOutputThreshold = 0, pvOutputThreshold.cis > 0
	me = Matrix_eQTL_main(
	snps = snps,
	gene = gene,
	cvrt = cvrt,
	output_file_name = outfile_trans,
	pvOutputThreshold = pvOutputThreshold_trans,
	useModel = modelLINEAR, # modelANOVA, modelLINEAR, or modelLINEAR_CROSS,
	errorCovariance = numeric(), # Error covariance matrix - Set to numeric() for identity.
	verbose = TRUE,
	output_file_name.cis = outfile_cis,
	pvOutputThreshold.cis = pvOutputThreshold_cis,
	snpspos = SNP,
	genepos = REG,
	cisDist = 1e6, # Distance for local gene-SNP pairs
	pvalue.hist = "qqplot",
	min.pv.by.genesnp = FALSE,
	noFDRsaveMemory = TRUE)
	cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');

	## Plot the Q-Q plot of local and distant p-values
	plot(me)
}



#cis-eQTLs
pdf("caqtl_scan_cis.pdf")
#runQTLScan(dat_test, outfile_cis="caqtl_scan_unstim_cis.txt", outfile_trans=NULL, pvOutputThreshold_cis=0.05, pvOutputThreshold_trans=0, covariates=c("PC1","PC2",paste0("PF",seq(2,10))))
#runQTLScan(dat_un, outfile_cis="caqtl_scan_unstim_cis.txt", outfile_trans=NULL, pvOutputThreshold_cis=0.05, pvOutputThreshold_trans=0, covariates=c("PC1","PC2",paste0("PF",seq(2,10))))
#runQTLScan(dat_stim, outfile_cis="caqtl_scan_stim_cis.txt", outfile_trans=NULL, pvOutputThreshold_cis=0.05, pvOutputThreshold_trans=0, covariates=c("PC1","PC2", paste0("PF",seq(2,10))))
#runQTLScan(dat_un_EUR, outfile_cis="caqtl_scan_unstim_cis_EUR.txt", outfile_trans=NULL, pvOutputThreshold_cis=0.05, pvOutputThreshold_trans=0, covariates=c("PC1_ancestry_specific","PC2_ancestry_specific", paste0("PF",seq(2,10))))
#runQTLScan(dat_un_AFR, outfile_cis="caqtl_scan_unstim_cis_AFR.txt", outfile_trans=NULL, pvOutputThreshold_cis=0.05, pvOutputThreshold_trans=0, covariates=c("PC1_ancestry_specific","PC2_ancestry_specific", paste0("PF",seq(2,10))))

runQTLScan(dat_un_EUR, outfile_cis="caqtl_scan_unstim_cis_EUR.txt", outfile_trans=NULL, pvOutputThreshold_cis=0.05, pvOutputThreshold_trans=0, covariates=c("PC1_ancestry_specific","PC2_ancestry_specific", "age_sample", "TSS_Score", "selection"))
runQTLScan(dat_un_AFR, outfile_cis="caqtl_scan_unstim_cis_AFR.txt", outfile_trans=NULL, pvOutputThreshold_cis=0.05, pvOutputThreshold_trans=0, covariates=c("PC1_ancestry_specific","PC2_ancestry_specific", "age_sample", "TSS_Score", "selection"))
dev.off()

#trans-eQTLs
# pdf("caqtl_scan_trans.pdf")
# runQTLScan(dat_test, outfile_cis=NULL, outfile_trans="caqtl_scan_trans_test.txt", pvOutputThreshold_cis=0, pvOutputThreshold_trans=1e-6,  covariates=c("PC1","PC2","age_sample",paste0("PF",seq(10))))
# dev.off()
