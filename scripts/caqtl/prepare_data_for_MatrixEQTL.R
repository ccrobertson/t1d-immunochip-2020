
setwd(Sys.getenv("freeze"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		Read in processed atac and sample meta data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sampleDF = readRDS("processed_sampleDF.rds")
countMat = readRDS("processed_countMat_clean_voomnormalized.rds")
doseDF = readRDS("mega_genotypes/mega_dosage_matrix_filtered.rds")
genoDF = readRDS("mega_genotypes/mega_genotype_matrix_filtered.rds")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		Remove (or fix?) sample swaps
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
swaps = read.table("sample_swaps.txt", header=TRUE)
sampleDF = sampleDF[!sampleDF$sample %in% swaps$sample,]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#		Get T1DGC data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
t1dgc_dir = "/nv/vol145/users/projects/T1DGC_SC_2011_release/T1DGC.2011.03_Phenotypic_Current/Analyses\ Files"
t1dgc_pheno = read.csv(paste0(t1dgc_dir,"/T1DGC.2011.03_Phenotypic_Current_Repository.csv"))
t1dgc_labs = read.csv(paste0(t1dgc_dir,"/T1DGC.2011.03_Phenotypic_Current_Labs.csv"))
t1dgc_hla = read.csv(paste0(t1dgc_dir,"/T1DGC.2011.03_Phenotypic_Current_HLA.csv"), colClasses="character")
t1dgc_forms = read.csv(paste0(t1dgc_dir,"/T1DGC.2011.03_Phenotypic_Current_Forms.csv"))

row.names(t1dgc_pheno) = t1dgc_pheno$analytic_id
row.names(t1dgc_labs) = t1dgc_labs$analytic_id
row.names(t1dgc_hla) = t1dgc_hla$analytic_id
row.names(t1dgc_forms) = t1dgc_forms$analytic_id

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare data sets -- separated by condition and ancestry
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getDatSubObject = function(samples_to_keep, prefix) {

  #get core data sets
  sampleDF_sub = sampleDF[sampleDF$sample %in% samples_to_keep,]
  countMat_sub = countMat[,sampleDF_sub$sample]
  genoDF_sub = genoDF[sampleDF_sub$iid,]
  doseDF_sub = doseDF[sampleDF_sub$iid,]
  dat_sub = list(samples=sampleDF_sub, counts=countMat_sub, genotypes=genoDF_sub, dosages=doseDF_sub)

  #add in t1dgc info
  dat_sub[["pheno"]] = t1dgc_pheno[dat_sub$sample$iid,]
  dat_sub[["labs"]] = t1dgc_labs[dat_sub$sample$iid,]
  dat_sub[["hla"]] = t1dgc_hla[dat_sub$sample$iid,]
  dat_sub[["forms"]] = t1dgc_forms[dat_sub$sample$iid,]

  #save files to use for peers
  write.table(t(dat_sub$counts), paste0(prefix,".counts.csv"), row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)
  write.table(data.frame(dat_sub$samples$T1D, dat_sub$pheno$age_sample), paste0(prefix,".samples.csv"), row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)

  #save object for R analysis
  saveRDS(dat_sub, file=paste0(prefix,".rds"))

  #return
  return(dat_sub)
}

samples_unstim = sampleDF$sample[sampleDF$iid %in% row.names(genoDF) & sampleDF$condition=="u"]
dat_un = getDatSubObject(samples_unstim, "DATASET_un")

samples_stim = sampleDF$sample[sampleDF$iid %in% row.names(genoDF) & sampleDF$condition=="s"]
dat_stim = getDatSubObject(samples_stim, "DATASET_stim")

samples_unstim_EUR = sampleDF$sample[sampleDF$iid %in% row.names(genoDF) & sampleDF$condition=="u" & sampleDF$cluster_label=="EUR"]
dat_un_EUR = getDatSubObject(samples_unstim_EUR, "DATASET_un_EUR")

samples_unstim_AFR = sampleDF$sample[sampleDF$iid %in% row.names(genoDF) & sampleDF$condition=="u" & sampleDF$cluster_label=="AFR"]
dat_un_AFR = getDatSubObject(samples_unstim_AFR, "DATASET_un_AFR")


### Create data set for testing
dat_test = dat_un
dat_test$counts = dat_test$counts[sample(1:nrow(dat_test$counts), 100),]
saveRDS(dat_test, file="DATASET_test.rds")
write.table(t(dat_test$counts), "DATASET_test.counts.csv", row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)
write.table(data.frame(dat_test$samples$T1D, dat_test$pheno$age_sample), "DATASET_test.samples.csv", row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)
