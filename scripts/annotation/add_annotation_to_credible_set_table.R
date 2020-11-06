setwd(Sys.getenv("finemap"))


### Get data
dat = read.table("credible_sets_with_proxies_from_Jamie.txt", header=TRUE)
#dat = read.table("all_creds_plus_ld_a12.txt", header=TRUE)
#dat_old = read.table("credible_sets_with_proxies_from_Jamie_20191001.txt", header=TRUE)
dat$ChrPos = paste0("chr",dat$chromosome, ":",dat$position)
annot = read.csv(paste0(Sys.getenv("annot"), "/annovar_output_credible_sets_with_proxies.hg38_multianno.csv"))
annot$ChrPos = paste0("chr",annot$Chr, ":",annot$End)

### Merge
#make sure chr:pos is a reasonable unique identifier for merging data sets
sum(!annot$ChrPos %in% dat$ChrPos)
sum(!dat$ChrPos %in% annot$ChrPos)
sum(duplicated(annot$ChrPos))
datnew = merge(dat, annot, by="ChrPos")

### Check for swapped alleles and flip beta when necessary
datnew$Allele1_new <-  datnew$Ref
datnew$Allele2_new <-  datnew$Alt
datnew$swap = (!is.na(datnew$Allele1) & !is.na(datnew$Allele2) & (datnew$Allele1==datnew$Allele2_new) & (datnew$Allele2==datnew$Allele1_new))
datnew$Effect_new[datnew$swap==TRUE] <-  (-1)*datnew$Effect[datnew$swap==TRUE]
datnew$Effect_new[datnew$swap==FALSE] <-  datnew$Effect[datnew$swap==FALSE]

### Reformat variables
datnew$OddsRatio = exp(datnew$Effect_new)
datnew$AF_1000G_AFR = datnew$X1000g2015aug_afr
datnew$AF_1000G_EUR = datnew$X1000g2015aug_eur
datnew$rsid_merged = sapply(datnew$ChrPos, function(x) {
    ID1 = datnew[datnew$ChrPos==x & !is.na(datnew$ChrPos),"ID"];
    ID2 = datnew[datnew$ChrPos==x & !is.na(datnew$ChrPos),"avsnp150"];
    if (ID1==ID2) {rsid=ID2}
    else if (ID1==".") {rsid=ID2}
    else if (ID2==".") {rsid=ID1}
    else {rsid=ID2}
  }
)
datnew$tag=gsub(".",":",datnew$tag, fixed=TRUE)
#datnew$tag_rsid = sapply(datnew$tag, function(x) { rsid = datnew$avsnp150[datnew$MarkerName==x & !is.na(datnew$MarkerName)]; if (rsid==".") {tag=x} else {tag=rsid}; return(tag)} )
datnew$tag_rsid = sapply(datnew$tag, function(x) { rsid = datnew$rsid_merged[datnew$MarkerName==x & !is.na(datnew$MarkerName)]; if (rsid==".") {tag=x} else {tag=rsid}; return(tag)} )
datnew$ppsum = sapply(datnew$tag, function(x) {datnew$ppsum[datnew$MarkerName==x & !is.na(datnew$MarkerName)]})
datnew$credset_size = sapply(datnew$tag, function(x) {d = datnew[datnew$tag==x,]; return(dim(d)[1])})

### Order data.frame by chromosome and region
datnew$tag_pos = sapply(datnew$tag, function(x) {datnew$position[datnew$MarkerName==x & !is.na(datnew$MarkerName)]})
datnew$region_pos = sapply(datnew$tag, function(x) {reg = datnew$region[datnew$MarkerName==x & !is.na(datnew$MarkerName)];  as.numeric(strsplit(reg, split="[-:]", perl=TRUE)[[1]][2])})
datnew_ordered = datnew[order(datnew$chromosome, datnew$region_pos, datnew$tag_pos, datnew$position),]

### Generate tables for saving
fields = c("region","region_pos","MarkerName","ID","avsnp150","tag","tag_rsid","tag_pos","credset_size","chromosome","position","Allele1_new","Allele2_new","AF_1000G_AFR","AF_1000G_EUR","OddsRatio","Effect_new","StdErr","P.value","ppsum","pp","excluded","cytoBand","Func.ensGene", "Func.refGene", "Gene.ensGene","Gene.refGene","ExonicFunc.ensGene","ExonicFunc.refGene", "AAChange.ensGene", "AAChange.refGene")
fields_manuscript = c("region","MarkerName","rsid_merged","tag","tag_rsid","credset_size","chromosome","position","Allele1_new","Allele2_new","AF_1000G_AFR","AF_1000G_EUR","OddsRatio","Effect_new","StdErr","P.value","ppsum","pp","excluded","cytoBand","Func.ensGene", "Func.refGene", "Gene.ensGene","Gene.refGene","ExonicFunc.ensGene","ExonicFunc.refGene", "AAChange.ensGene", "AAChange.refGene")
OUTTABLE_all = datnew_ordered[,fields_manuscript]
OUTTABLE_protein_altering = OUTTABLE_all[OUTTABLE_all$ppsum>0.5 & (OUTTABLE_all$ExonicFunc.refGene %in% c("stopgain","nonsynonymous SNV","frameshift substitution") | OUTTABLE_all$Func.refGene=="splicing" | OUTTABLE_all$ExonicFunc.ensGene %in% c("stopgain","nonsynonymous SNV","frameshift substitution") | OUTTABLE_all$Func.ensGene=="splicing"), ]
OUTTABLE_top_candidates = OUTTABLE_protein_altering[!is.na(OUTTABLE_protein_altering$pp) & OUTTABLE_protein_altering$pp>0.1,]

### Save to files (Supplementary Tables 11 & 13)
write.table(OUTTABLE_all, file="credible_sets_with_proxies_from_Jamie_ANNOTATED.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(OUTTABLE_protein_altering, file="credible_sets_with_proxies_from_Jamie_ANNOTATED_protein_altering.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(OUTTABLE_top_candidates, file="credible_sets_with_proxies_from_Jamie_ANNOTATED_protein_altering_top_candidates.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

#rsync -v $rivanna:/nv/vol185/MEGA/release4/IMPUTED_TOPMED/fine_mapping/credible_sets_with_proxies_from_Jamie_ANNOTATED.txt .
#rsync -v $rivanna:/nv/vol185/MEGA/release4/IMPUTED_TOPMED/fine_mapping/credible_sets_with_proxies_from_Jamie_ANNOTATED_protein_altering.txt .
#rsync -v $rivanna:/nv/vol185/MEGA/release4/IMPUTED_TOPMED/fine_mapping/credible_sets_with_proxies_from_Jamie_ANNOTATED_protein_altering_top_candidates.txt .


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## save credible set rds object for later use
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
finemap = read.table("credible_sets_with_proxies_from_Jamie_ANNOTATED.txt", header=TRUE, comment.char="~", sep="\t")
finemap$chr = paste0("chr", finemap$chromosome)
finemap$MarkerID = finemap$ID
finemap$MarkerID[finemap$MarkerID=="."] <- finemap$MarkerName[finemap$MarkerID=="."]
finemap$include_line=(finemap$ppsum>0.8 & finemap$credset_size<=5)
saveRDS(finemap, file="annotated_credible_sets.rds")
