library(dplyr)

### Read sample QC and info files
sampleDF = readRDS(paste0(Sys.getenv("freeze"), "/processed_sampleDF.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Check for matches
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Now check whether sample is best match
getMatch = function(sample) {
  sample = sampleDF[sampleDF$sample==sample, "sample"]
  SQB = sampleDF[sampleDF$sample==sample, "SQB"]
  iid = sampleDF[sampleDF$sample==sample, "iid"]
  d = read.table(paste0(Sys.getenv("freeze"),"/sample_match/",sample,".match.txt"), header=TRUE)

  #get stats for individual with matching subject id
  if (iid %in% d$SampleID) {
    iid_n = d[d$SampleID==iid,"n_het_covered"]
    iid_percent = d[d$SampleID==iid,"perc_het_consistent"]
  } else {
    iid_n = NA
    iid_percent = NA
  }

  d_ordered = d[order(d$perc_het_consistent, decreasing=TRUE),]
  topmatch_iid = d_ordered[1,"SampleID"]
  topmatch_n = d_ordered[1,"n_het_covered"]
  topmatch_percent = d_ordered[1,"perc_het_consistent"]

  return(data.frame(sample, SQB, iid, iid_n, iid_percent, topmatch_iid, topmatch_n, topmatch_percent))
}

### Extract match stats
matches = do.call("rbind", lapply(sampleDF$sample, getMatch))
write.table(matches, file=paste0(Sys.getenv("freeze"),"/sample_matches.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

### Find correct samples
pass = matches[matches$iid==matches$topmatch_iid & matches$topmatch_percent>0.9,]
summary(pass$topmatch_percent)

### Find samples where top match is not consisent with IID
fail = matches[!(matches$iid==matches$topmatch_iid & matches$topmatch_percent>0.9),]
write.table(fail, file=paste0(Sys.getenv("freeze"),"/samples_no_match.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Investigate Swaps
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### For potential swaps, run agains the full mega EUR and AFR vcfs
dir.create(paste0(Sys.getenv("freeze"),"/sample_match_full_cohort"))
for (i in 1:nrow(fail)) {
  sample = fail[i, "sample"]
  SQB = fail[i, "SQB"]

  command_template1 = "sbatch --output=${freeze}/sample_match_full_cohort/${sample}.match_EUR.log ${scripts}/caqtl/sample_matching_full_cohort.slurm ${sample} ${SQB} EUR"
  command_template2 = "sbatch --output=${freeze}/sample_match_full_cohort/${sample}.match_AFR.log ${scripts}/caqtl/sample_matching_full_cohort.slurm ${sample} ${SQB} AFR"

  command1 = command_template1 %>%
    gsub("${freeze}", Sys.getenv("freeze"), ., fixed=TRUE) %>%
    gsub("${scripts}", Sys.getenv("scripts"), ., fixed=TRUE) %>%
    gsub("${sample}", sample, ., fixed=TRUE) %>%
    gsub("${SQB}", SQB, ., fixed=TRUE)
  command2 = command_template2 %>%
    gsub("${freeze}", Sys.getenv("freeze"), ., fixed=TRUE) %>%
    gsub("${scripts}", Sys.getenv("scripts"), ., fixed=TRUE) %>%
    gsub("${sample}", sample, ., fixed=TRUE) %>%
    gsub("${SQB}", SQB, ., fixed=TRUE)

  cat(command1,"\n")
  cat(command2,"\n")
  system(command1)
  system(command2)
}


getMatch_FullCohort = function(sample) {
  sample = sampleDF[sampleDF$sample==sample, "sample"]
  SQB = sampleDF[sampleDF$sample==sample, "SQB"]
  iid = sampleDF[sampleDF$sample==sample, "iid"]
  d1 = read.table(paste0(Sys.getenv("freeze"),"/sample_match_full_cohort/",sample,".match_EUR.txt"), header=TRUE)
  d2 = read.table(paste0(Sys.getenv("freeze"),"/sample_match_full_cohort/",sample,".match_AFR.txt"), header=TRUE)
  d = rbind(d1, d2)

  #get stats for individual with matching subject id
  if (iid %in% d$SampleID) {
    iid_n = d[d$SampleID==iid,"n_het_covered"]
    iid_percent = d[d$SampleID==iid,"perc_het_consistent"]
  } else {
    iid_n = NA
    iid_percent = NA
  }

  d_ordered = d[order(d$perc_het_consistent, decreasing=TRUE),]
  topmatch_iid = d_ordered[1,"SampleID"]
  topmatch_n = d_ordered[1,"n_het_covered"]
  topmatch_percent = d_ordered[1,"perc_het_consistent"]

  return(data.frame(sample, SQB, iid, iid_n, iid_percent, topmatch_iid, topmatch_n, topmatch_percent))
}

swap_matches = do.call("rbind", lapply(fail$sample, getMatch_FullCohort))
swaps = swap_matches[swap_matches$topmatch_percent>0.9,]
write.table(swaps, file=paste0(Sys.getenv("freeze"),"/sample_swaps.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Investigate missing subjects
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Look into subjects missing from VCF
missing = fail[!fail$iid %in% swaps$iid,]

#removed in mega QC
sample_droplist = read.table(paste0(Sys.getenv("PROJECT_MEGA"),"/pre_qc/Mega_10April2019/prelim_QC/QC_lists/sampletoberemoved_all.txt"), header=FALSE)
names(sample_droplist) <- c("iid","reason")
missing_dropped = sample_droplist[sample_droplist$iid %in% missing$iid,]

paste(missing_dropped[missing_dropped$reason=="autoqc:MissingMoreThan2","iid"], collapse=",")
paste(missing_dropped[missing_dropped$reason=="relqc:overlapping_family_within_cohort","iid"], collapse=",")

#non EUR or AFR ancestry
missing_not_dropped = missing[!missing$iid %in% sample_droplist$iid,]
pcs = read.table(paste0(Sys.getenv("pca"),"/mega_cluster_assignments_and_pcs.txt"), header=TRUE, sep="\t", comment.char="~")
pcs[pcs$IID %in% missing_not_dropped$iid,]
