library(vcfR)
setwd(paste0(Sys.getenv("freeze"),"/mega_genotypes"))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Convert to matrix
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#read in vcf
vcf = read.vcfR("mega.vcf.gz")
gt_t = extract.gt(vcf, element = 'GT')
ds_t = extract.gt(vcf, element = 'DS')

#exclude variants that aren't biallelic
gt_t[!gt_t %in% c("0|0", "0|1","1|0", "1|1")] <- NA
ds_t[!gt_t %in% c("0|0", "0|1","1|0", "1|1")] <- NA

#convert gt to allele counts
gt_t[gt_t=="0|0"] <- 0
gt_t[gt_t=="0|1"] <- 1
gt_t[gt_t=="1|0"] <- 1
gt_t[gt_t=="1|1"] <- 2


#transpose
gt = as.data.frame(apply(t(gt_t), MARGIN=2, FUN=as.numeric))
rownames(gt) <- colnames(gt_t)
colnames(gt) <- gsub(":","_",rownames(gt_t))

ds = as.data.frame(apply(t(ds_t), MARGIN=2, FUN=as.numeric))
rownames(ds) <- colnames(ds_t)
colnames(ds) <- gsub(":","_",rownames(ds_t))

#save
saveRDS(gt, file="mega_genotype_matrix.rds")
saveRDS(ds, file="mega_dosage_matrix.rds")

gt = readRDS("mega_genotype_matrix.rds")
ds = readRDS("mega_dosage_matrix.rds")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Filter variants
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#removing variants with missingness
missingness = apply(gt,2, function(x) {sum(is.na(x))})
variants_nonmissing = names(missingness[missingness==0])

#exclude variants dropped in association QC --> where the first character in the “direction" column==“?” - I.e. Failed QC in EUR
load(paste0(Sys.getenv("meta"),"/resultsmeta_5pcs_vcf_5pc_diff.RData"))
res$MarkerName2 = gsub(":","_",res$MarkerName)
res$Direction_1 = sapply(res$Direction, function(x) {substring(x, 1,1)})
res_keep = res[!res$Direction_1=="?",]
gt_keep = gt[,names(gt)%in%intersect(res_keep$MarkerName2,variants_nonmissing)]
ds_keep = ds[,names(ds)%in%intersect(res_keep$MarkerName2,variants_nonmissing)]

#save
saveRDS(gt_keep, file="mega_genotype_matrix_filtered.rds")
saveRDS(ds_keep, file="mega_dosage_matrix_filtered.rds")
