
### Get GWAS catalog data
resources = "/Users/Cassie/Box Sync/Projects/MEGA/"
gwas_cat = read.table(file.path(resources,"gwas-association-downloaded_2020-03-06-EFO_0001359-withChildTraits.tsv"), header=TRUE, sep="\t", comment.char="~")

### Get current tables
current_version="/Users/Cassie/UVA_Dropbox/Dropbox/mega_ng/v1.6"
mega_res_phase1 = read.csv(file.path(current_version,"sup_tables_all_ST6.csv"), header=TRUE, comment.char="~")
mega_res_meta = read.table(file.path(current_version,"sup_tables_all_ST8.txt"), header=TRUE, sep="\t", comment.char="~")


### Filter for true T1D studies
t1d_studies = c("GCST008377", "GCST005536","GCST001394", "GCST001255", "GCST000539", "GCST000392", "GCST000258","GCST000244","GCST000141","GCST000054", "GCST000054","GCST000043","GCST000038")
gwas_cat_t1d = gwas_cat[gwas_cat$STUDY.ACCESSION%in%t1d_studies & gwas_cat$P.VALUE<=5e-8,]

### Define regions aroung GWAS Catalog lead snps
gwas_cat_t1d$CHR_ID[gwas_cat_t1d$CHR_ID=="X"] <- 23
gwas_cat_t1d$region_chr = as.numeric(gwas_cat_t1d$CHR_ID)
gwas_cat_t1d$region_start = as.numeric(gwas_cat_t1d$CHR_POS) - 250000
gwas_cat_t1d$region_end = as.numeric(gwas_cat_t1d$CHR_POS) + 250000


### Define function
getNovel = function(mega_results) {
  novel = NULL
  for (i in 1:nrow(mega_results)) {
    chromosome=mega_results$chromosome[i]
    position=mega_results$position[i]
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
  return(novel)
}

### Get novel region indicator for each table
mega_res_phase1$Novel = getNovel(mega_res_phase1)
mega_res_meta$Novel = getNovel(mega_res_meta)
mega_res_meta$Novel[mega_res_meta$Genomewide.FDR=="FDR"] <- NA


### Make sure novel column is consistent across tables
merged1 = merge(mega_res_phase1, mega_res_domrec, by="MarkerName")
merged2 = merge(mega_res_phase1, mega_res_meta, by="MarkerName")
merged2$NOVEL.x == merged2$Novel.y
sum(!merged1$Novel.x == merged1$Novel.y)
sum(!merged2$Novel.x == merged2$Novel.y)


### Write to file
write.table(mega_res_phase1, file=file.path(current_version,"sup_tables_all_ST6_updated.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
write.table(mega_res_meta, file=file.path(current_version,"sup_tables_all_ST8_updated.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
