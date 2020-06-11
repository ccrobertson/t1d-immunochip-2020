#install.packages("bedr")
library(bedr)
#module load htslib

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Prep annovar file for all meta-analysis variants
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load(paste0(Sys.getenv('meta'),'/resultsmeta_5pcs_vcf_5pc_diff.RData'))
annovar_input = res[order(res$chromosome, res$position),c('chromosome', 'position','position','Allele1','Allele2')]

annovar_input$swap = apply(annovar_input, 1, function(x) {nchar(x["Allele1"])>nchar(x["Allele2"])})
annovar_input$Allele1_swapped = apply(annovar_input, 1, function(x) {if (x["swap"]==TRUE) {x["Allele2"]} else {x["Allele1"]} })
annovar_input$Allele2_swapped = apply(annovar_input, 1, function(x) {if (x["swap"]==TRUE) {x["Allele1"]} else {x["Allele2"]} })

write.table(annovar_input[,c('chromosome','position','position.1','Allele1_swapped','Allele2_swapped')],file=paste0(Sys.getenv('annot'),'/annovar_input_file.txt'),row.names=FALSE, col.names=FALSE, quote=FALSE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Prep annovar file for credible variants
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Read data
d = read.table(paste0(Sys.getenv('finemap'),'/credible_sets_with_proxies_from_Jamie.txt'), header=TRUE)

### Add ref and alt alleles to variants missing them
#get snps with ref/alt missing
snps_missing_alleles = d[d$exclude & is.na(d$Allele1),c("ID","chromosome","position")]
snps_missing_alleles$query = paste0(snps_missing_alleles$chromosome,":",snps_missing_alleles$position-1,"-",snps_missing_alleles$position)
snps_missing_alleles$file_1000g = paste0(Sys.getenv("resources"),"/1000g_b38/ALL.chr",snps_missing_alleles$chromosome,"_GRCh38_sites.20170504.vcf.gz")

#extract ref/alt for each snpfrom 1000g b38 files
info_1000g = do.call("rbind",apply(snps_missing_alleles, 1, function(x) { tabix(region=x["query"], file.name=x["file_1000g"], check.chr=FALSE)}))

#first merge on rsid
m1 = merge(snps_missing_alleles, info_1000g, by=c("ID"))
m1$ID.old = m1$ID
m1$ID.new = m1$ID

#merge everything else on chrpos
m2 = merge(snps_missing_alleles[!snps_missing_alleles$ID %in% m1$ID,], info_1000g, by.x=c("chromosome","position"), by.y=c("CHROM","POS"),suffixes=c(".old",".new"))
m2$ID.new[m2$ID.new=="rs9582295;rs55960850"] <- "rs55960850"

#combine
update = rbind(m1[,c("ID.old","ID.new","chromosome","position","REF","ALT")], m2[,c("ID.old","ID.new","chromosome","position","REF","ALT")])
row.names(update) = update$ID.old

#add updated REF and ALT to fine-mapping dataframe
d$ID_new = apply(d, 1 , function(x) { snp=x["ID"]; if(snp  %in% row.names(update) ) { update[snp, "ID.new"]} else {snp }})
d$Allele1_new = apply(d, 1 , function(x) { snp=x["ID"]; if(snp  %in% row.names(update) ) { update[snp, "REF"]} else {x["Allele1"] }})
d$Allele2_new = apply(d, 1 , function(x) { snp=x["ID"]; if(snp  %in% row.names(update) ) { update[snp, "ALT"]} else {x["Allele2"]}})

### Swap allele for indels (want smaller allele to be Allele1)
d$swap = apply(d, 1, function(x) {nchar(x["Allele1_new"])>nchar(x["Allele2_new"])})
d$Allele1_new_swapped = apply(d, 1, function(x) {if (x["swap"]==TRUE) {x["Allele2_new"]} else {x["Allele1_new"]} })
d$Allele2_new_swapped = apply(d, 1, function(x) {if (x["swap"]==TRUE) {x["Allele1_new"]} else {x["Allele2_new"]} })


### find alleles based on dbsnp 151 for proxy snps (this didn't work -- the variants that we are missing are also not in this dbsnp151 database)
# snps <- SNPlocs.Hsapiens.dbSNP151.GRCh38
# map = as.data.frame(snpsById(snps, ids=d$ID[!d$ID=="." & d$excluded], ifnotfound="drop"))
# map$Alleles = IUPAC_CODE_MAP[map$alleles_as_ambig]
# map$Allele1 = sapply(map$Alleles, function(x) {if(nchar(x)==2) {A=substr(x,1,1)} else {A=NA}; return(A)})
# map$Allele2 = sapply(map$Alleles, function(x) {if(nchar(x)==2) {A=substr(x,2,2)} else {A=NA}; return(A)})
# m = merge(d, map, all=TRUE, by.x="ID", by.y="RefSNP_id")
# m[is.na(m$Allele1.x),c("ID","MarkerName","Allele1.y","Allele2.y")]

### Write to file in annovar format
annovar_input = d[,c('chromosome','position','position','Allele1_new_swapped','Allele2_new_swapped')]
write.table(annovar_input, file=paste0(Sys.getenv('annot'),'/annovar_input_file_credible_sets_with_proxies.txt'),row.names=FALSE, col.names=FALSE, quote=FALSE)
