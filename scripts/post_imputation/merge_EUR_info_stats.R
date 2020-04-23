#args = c(2, "/m/CPHG/MEGA/IMPUTED_1KG_rel3/results/EUR1/chr2.info.gz","/m/CPHG/MEGA/IMPUTED_1KG_rel3/results/EUR2/chr2.info.gz", "/m/CPHG/MEGA/IMPUTED_1KG_rel3/results/EUR3/chr2.info.gz", "EUR")
#args = c(1, "/m/CPHG/MEGA/IMPUTED_1KG_rel3/results/EUR1/chr1.info.gz","/m/CPHG/MEGA/IMPUTED_1KG_rel3/results/EUR2/chr1.info.gz", "/m/CPHG/MEGA/IMPUTED_1KG_rel3/results/EUR3/chr1.info.gz", "EUR")
args = commandArgs(trailingOnly=TRUE)
chr=args[1]
file1=args[2]
file2=args[3]
file3=args[4]
outdir=args[5]

getData = function(file) {
	d = read.table(file, comment.char="~", header=TRUE)
	d$Rsq = as.numeric(d$Rsq)
	d$ALT_frq = as.numeric(d$ALT_Frq)
	if (sum(duplicated(d$SNP))>0) {
		d2 = d[!duplicated(d$SNP),]
		row.names(d2) <- d2$SNP
		return(d2)
	} else {
		row.names(d) <- d$SNP
		return(d)
	}
}

d1 = getData(file1)
d2 = getData(file2)
d3 = getData(file3)

#d1 = read.table(file1, comment.char="~", header=TRUE)
#d1 = d1[unique(d1$SNP),]
#row.names(d1) <- d1$SNP

#d2 = read.table(file2, comment.char="~", header=TRUE)
#d2 = d2[unique(d2$SNP),]
#row.names(d2) <- d2$SNP

#d3 = read.table(file3, comment.char="~", header=TRUE)
#d3 = d3[unique(d3$SNP),]
#row.names(d3) <- d3$SNP


#extract intersecting variants (variants with MAF>0 in all batches)
common_vars = Reduce(intersect, list(d1$SNP[d1$MAF>0], d2$SNP[d2$MAF>0], d3$SNP[d3$MAF>0] ))
df = data.frame(d1[common_vars,c("ALT_Frq","Rsq")], d2[common_vars,c("ALT_Frq","Rsq")], d3[common_vars,c("ALT_Frq","Rsq")])
names(df) <- c("p1", "rsq1","p2","rsq2","p3","rsq3")
df$ALT_Frq = df$p
df$SNP = row.names(df)

#calculate pooled alt allele freq
n1 = 15000
n2 = 15000
n3 = 17049
N = sum(n1+n2+n3)
df$p = with(df, (n1*p1 + n2*p2 + n3*p3)/N)
df$rsq = with(df, (n1*(rsq1*p1*(1-p1)+(p1-p)^2) + n2*(rsq2*p2*(1-p2)+(p2-p)^2) + n3*(rsq3*p3*(1-p3)+(p3-p)^2))/(N*p*(1-p)) )
getMAF = function(p) {
	min(p, 1-p)
}
df$ALT_Frq = df$p
df$MAF = unlist(lapply(df$ALT_Frq, getMAF))

write.table(df[,c("SNP","ALT_Frq","MAF","rsq")], file=paste(outdir,"/chr",chr,".merged_info",sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")





#par(mfrow=c(2,2))
#with(df, plot(rsq1, rsq))
#with(df, plot(rsq2, rsq))
#with(df, plot(rsq3, rsq))


###Compare frequencies across other ethnicities batches
#dAfr_1 = read.table("FIN/chr22.info.gz", comment.char="~", header=TRUE)
#dAfr = dAfr_1[dAfr_1$Rsq>0.3,]
#row.names(dAfr) <- dAfr$SNP
#
#common_vars_afr = Reduce(intersect, list(d1$SNP[d1$MAF>0], d2$SNP[d2$MAF>0], d3$SNP[d3$MAF>0], dAfr$SNP[dAfr$MAF>0]))
#df2 = data.frame(d1[common_vars_afr,c("ALT_Frq","Rsq")], d2[common_vars_afr,c("ALT_Frq","Rsq")], d3[common_vars_afr,c("ALT_Frq","Rsq")], dAfr[common_vars_afr, c("ALT_Frq", "Rsq")])
#names(df2) <- c("p1", "rsq1","p2","rsq2","p3","rsq3", "p4","rsq4")
#
#
#par(mfrow=c(2,2))
#with(df2, plot(p1,p2))
#with(df2, plot(p3,p2))
#with(df2, plot(p1,p3))
#with(df2, plot(p1,p4))
