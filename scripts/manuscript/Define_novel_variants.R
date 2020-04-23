setwd("/Users/Cassie/Box\ Sync/Projects/MEGA")

### Get data
gwas_cat = read.table("gwas-association-downloaded_2020-03-06-EFO_0001359-withChildTraits.txt", header=TRUE, sep="\t", comment.char="~")
mega_res = read.table("mega_meta_analysis_suptab7.txt", header=TRUE, sep="\t", comment.char="~")
mega_res_phase1 = read.table("mega_meta_analysis_suptab6.txt", header=TRUE, sep="\t", comment.char="~")


### Filter for true T1D studies
t1d_studies = c("GCST008377", "GCST005536","GCST001394", "GCST001255", "GCST000539", "GCST000392", "GCST000258","GCST000244","GCST000141","GCST000054", "GCST000054","GCST000043","GCST000038")
gwas_cat_t1d = gwas_cat[gwas_cat$STUDY.ACCESSION%in%t1d_studies & gwas_cat$P.VALUE<=5e-8,]

### Define regions aroung GWAS Catalog lead snps
gwas_cat_t1d$CHR_ID[gwas_cat_t1d$CHR_ID=="X"] <- 23
gwas_cat_t1d$region_chr = as.numeric(gwas_cat_t1d$CHR_ID)
gwas_cat_t1d$region_start = as.numeric(gwas_cat_t1d$CHR_POS) - 250000
gwas_cat_t1d$region_end = as.numeric(gwas_cat_t1d$CHR_POS) + 250000


### Get novel regions for combined analysis
novel_phase1 = NULL
for (i in 1:nrow(mega_res_phase1)) {
  chromosome=mega_res_phase1$chromosome[i]
  position=mega_res_phase1$position[i]
  in_region = NULL
  for (k in 1:nrow(gwas_cat_t1d)) {
    in_region = c(in_region, (gwas_cat_t1d[k,"region_chr"]==chromosome & gwas_cat_t1d[k,"region_start"]<position & gwas_cat_t1d[k,"region_end"]>position))
  }
  cat(chromosome, position, sum(in_region),"\n")
  if (sum(in_region)==0) {
    novel_phase1 = c(novel_phase1, 1)
  } else {
    novel_phase1 = c(novel_phase1, 0)
  }
}
sum(novel_phase1)
mega_res_phase1$NOVEL[novel_phase1==1]<-1
mega_res_phase1$NOVEL[novel_phase1==0]<-0
write.table(mega_res_phase1, file="mega_meta_analysis_suptab6_updated.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")


### Get novel regions for combined analysis
novel = NULL
for (i in 1:nrow(mega_res)) {
  chromosome=mega_res$chromosome[i]
  position=mega_res$position[i]
  in_region = NULL
  for (k in 1:nrow(gwas_cat_t1d)) {
    in_region = c(in_region, (gwas_cat_t1d[k,"region_chr"]==chromosome & gwas_cat_t1d[k,"region_start"]<position & gwas_cat_t1d[k,"region_end"]>position))
  }
  cat(chromosome, position, sum(in_region),"\n")
  if (sum(in_region)==0) {
    novel = c(novel, 1)
  } else {
    novel = c(novel, 0)
  }
}

sum(novel[mega_res$Genomewide.FDR=="Genomewide"])
mega_res$NOVEL[novel==1]<-1
mega_res$NOVEL[novel==0]<-0
mega_res$NOVEL[mega_res$Genomewide.FDR=="FDR"]<-NA
write.table(mega_res, file="mega_meta_analysis_suptab7_updated.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")


##Look at just the novel ones
write.table(mega_res[mega_res$Genomewide.FDR=="Genomewide" & mega_res$NOVEL==1,], file="mega_meta_analysis_suptab7_novel_only.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
