library(dplyr)

### Read sample QC and info files
sampleDF = readRDS(paste0(Sys.getenv("freeze"), "/processed_sampleDF.rds"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Run QTLtools on QTL VCF
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir.create(paste0(Sys.getenv("freeze"),"/sample_match"))

### For each sample run QTLtools
for (i in 1:nrow(sampleDF)) {
  sample = sampleDF[i, "sample"]
  SQB = sampleDF[i, "SQB"]

  command_template = "sbatch --output=${freeze}/sample_match/${sample}.match.log ${scripts}/caqtl/sample_matching.slurm ${sample} ${SQB}"
  command = command_template %>%
    gsub("${freeze}", Sys.getenv("freeze"), ., fixed=TRUE) %>%
    gsub("${scripts}", Sys.getenv("scripts"), ., fixed=TRUE) %>%
    gsub("${sample}", sample, ., fixed=TRUE) %>%
    gsub("${SQB}", SQB, ., fixed=TRUE)

  cat(command,"\n")
  system(command)
}
