###### ###### ###### ###### ###### ###### 
###### PCA analysis of array data  ###### 
###### ###### ###### ###### ###### ###### 

options(stringsAsFactor=FALSE)
library("snpStats")
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
bed = args[1]
bim = args[2]
fam = args[3]
nickname = args[4]
#setwd("/data1/ccr5ju/MEGA/pre_qc/Mega_1Mar2018/raw")

pdf(paste("pca/pca_",nickname,".pdf", sep=""))
sample <- read.plink(bed=bed, bim=bim, fam=fam)

snps = sample$genotypes
fam = sample$fam
map = sample$map

snps
head(fam)
head(map)

#note, in vignette
#snps.10 = sample$genotypes
#snp.support = sample$map
#sample.support = sample$fam

#thinning snps (take every 10th snp)
use <- seq(1, ncol(snps), 10)
snps.pca <- snps[,use]

#calculate eigenvectors
xxmat <- xxt(snps.pca, correct.for.missing=FALSE)
evv <- eigen(xxmat, symmetric=TRUE)

#take first 5 pcs
pcs <- evv$vectors[,1:10]
evals <- evv$values[1:10]
evals

saveRDS(pcs, paste("pca/pcs_",nickname,".rds", sep=""))

#plot pcs
plot(pcs[,1], pcs[,2], xlab="PC1", ylab="PC2", main=nickname)

dev.off()

#####################

#
#d1 = read.table("mega_main_raw/mega_main0.frq", header=TRUE)
#d2 = read.table("mega_additions_raw/mega_additions0.frq", header=TRUE)
#d3 = read.table("mega_cbr_raw/mega_cbr0.frq", header=TRUE)
#
#dim(d1)
#dim(d2)
#dim(d3)
#
##Check intersection of SNP ids between files
#d1[!d1$SNP%in%d2$SNP,]
#d2[!d2$SNP%in%d1$SNP,]
#
#d1[!d1$SNP%in%d3$SNP,]
#d3[!d3$SNP%in%d1$SNP,]
#
#
##Merge files by SNP id
#row.names(d1) = d1$SNP
#row.names(d2) = d2$SNP
#row.names(d3) = d3$SNP
#m = data.frame(SNP=d1$SNP, MAF_main=d1[d1$SNP,"MAF"], MAF_add=d2[d1$SNP,"MAF"],MAF_cbr=d3[d1$SNP,"MAF"])
#dim(m); head(m)
#
##Plot AF 
#ggplot(m, aes(x=MAF_main, y=MAF_add)) + geom_point()