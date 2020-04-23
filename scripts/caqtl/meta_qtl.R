library(meta)
setwd(Sys.getenv("freeze"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Meta-analysis of caQTL EUR and AFR results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
res_EUR = read.table("caqtl_scan_unstim_cis_EUR.txt", header=TRUE)
res_EUR$SE = res_EUR$beta/res_EUR$t.stat
res_AFR = read.table("caqtl_scan_unstim_cis_AFR.txt", header=TRUE)
res_AFR$SE = res_AFR$beta/res_AFR$t.stat

res_merged = merge(res_EUR, res_AFR, by=c("SNP","gene"), suffixes=c("_EUR","_AFR"))

metaAnalysis = function(x, groups) {
  fields_beta = paste("beta",groups, sep="_")
  fields_SE = paste("SE",groups, sep="_")
  out = metagen(TE=as.numeric(x[fields_beta]), seTE=as.numeric(x[fields_SE]), comb.fixed=TRUE)
  unlist(out[c("TE.fixed", "seTE.fixed", "zval.fixed", "pval.fixed", "pval.Q")])
}
res_meta = data.frame(res_merged, t(apply(res_merged, 1, metaAnalysis, groups=c("EUR","AFR"))))
write.table(res_meta, file="caqtl_scan_unstim_cis_meta.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
