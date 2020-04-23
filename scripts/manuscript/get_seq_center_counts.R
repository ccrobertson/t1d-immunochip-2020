
tracking = read.csv("raw_data/10022015-Updated_Complete_Immuno_Sample_Sheet.csv", skip=8)

pheno = read.table("mega_genotyped/release4/mega_b37_clean_pheno_genotypedOnly.txt", header=FALSE)
names(pheno) <- c("Cohort","FID","IID","FAT","MOT","Sex", "T1D")
update = read.table("raw_data/phenotype_files/updateID_6.txt",header=FALSE)
pheno$IID_old = sapply(pheno$IID, function(x) {
  if (x %in% update$V4) {
    old=update$V2[update$V4==x][1]
  } else {
    old=NA
  };
  return(old)
})


sum(!(pheno$IID %in% tracking$Sample_ID | pheno$IID_old %in% tracking$Sample_ID))

write.csv(pheno[!(pheno$IID %in% tracking$Sample_ID | pheno$IID_old %in% tracking$Sample_ID),c("Cohort","FID","IID")], file="mega_not_on_sample_sheet.csv", row.names=FALSE)

pheno$IID_merged =
pheno$GenoCenter = apply(pheno, MARGIN=1, FUN=function(x) {
  if (x["IID"] %in% tracking$Sample_ID) {
    return(tracking$GenotypingCenter[tracking$Sample_ID==x["IID"]][1])
  } else if (x["IID_old"] %in% tracking$Sample_ID) {
    return(tracking$GenotypingCenter[tracking$Sample_ID==x["IID_old"]][1])
  } else {
    return(NA)
  }
})
table(pheno$GenoCenter)

Cambridge 2941
Feinstein Institute 1811
Sanger 4347
UVA 52219
