options(stringsAsFactors=FALSE)

args = commandArgs(trailingOnly=TRUE)
prefix = args[1]
xHetThresh = args[2]
yThreshLow = args[3]
yThreshHigh = args[4]

fileBySample = paste(prefix,"bySample.txt", sep="")
d = read.table(fileBySample, header=TRUE)
d$sex_color = NA
d$sex_color[d$SEX==1] <- "red"
d$sex_color[d$SEX==2] <- "blue"

pdf(paste(prefix,"pdf",sep="."))
plot(d$xHeterozygosity, d$N_ySNP, xlab="X Heterozygosity", ylab="non-missing Y SNPs", col=d$sex_color)
abline(v=xHetThresh)
abline(h=yThreshHigh)
abline(h=yThreshLow)
dev.off()

