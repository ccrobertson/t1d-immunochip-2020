library(vcfR)
library(ggplot2)
library(reshape2)
library(lme4)
library(cowplot)


setwd(Sys.getenv("freeze"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get data sets
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat_un = readRDS("DATASET_un.rds")
finemap = readRDS(file=paste0(Sys.getenv("finemap"),"/annotated_credible_sets.rds"))
annot =  read.csv(paste0(Sys.getenv("annot"),"/annovar_output.hg38_multianno.csv"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get caqtl summary stats
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
results_qtl_un = readRDS(file="joint_caqtl_models_un.rds")
#results_qtl_AFR = readRDS(file="joint_caqtl_models_AFR.rds")
#results_qtl_EUR = readRDS(file="joint_caqtl_models_EUR.rds")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Collect variant info
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getVariantInfo = function(rsid) {
  if(rsid %in% finemap$MarkerID) {
    altid = finemap$MarkerName[finemap$MarkerID==rsid]
    params = unlist(strsplit(altid, split=":"))
    chr = params[1]
    start = params[2]
    end = params[2]
    ref = params[3]
    alt = params[4]
  } else if (rsid %in% annot$avsnp150) {
    params = as.list(annot[annot$avsnp150==rsid,c("Chr","Start","End","Ref","Alt")])
    chr = paste0("chr",params[[1]])
    start = params[[2]]
    end = params[[3]]
    ref = params[[4]]
    alt = params[[5]]
    altid = paste(chr, start, ref, alt, sep=":")
  } else {
    chr=NA
    start=NA
    end=NA
    ref=NA
    alt=NA
    altid=NA
  }
  return(list(rsid=rsid, altid=altid, chr=chr, start=start, end=end, ref=ref, alt=alt))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Collect sample info
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getSampleInfo = function(variant, dat_sub) {
  sampleDF = dat_sub[["samples"]]
  geno_mat = dat_sub[["genotypes"]]

  #create directory to store sample info
  vardir = paste0(Sys.getenv("freeze"),"/allelic_bias/",variant[["rsid"]])
  if (!dir.exists(vardir)) {
    dir.create(vardir, showWarnings=TRUE, recursive=TRUE)
  }

  #extract sample iids for each genotype
  varColName = gsub(":","_", variant[["altid"]])
  iids0 = row.names(geno_mat)[geno_mat[,varColName]==0]
  iids1 = row.names(geno_mat)[geno_mat[,varColName]==1]
  iids2 = row.names(geno_mat)[geno_mat[,varColName]==2]

  #create file listing sample info for each genotype
  maptorun = sampleDF[,c("SQB", "iid", "sample")]
  row.names(maptorun) = maptorun$iid
  write.table(maptorun[iids0,], file=paste0(vardir, "/", variant[["rsid"]],"_HOMREF.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
  write.table(maptorun[iids1,], file=paste0(vardir, "/", variant[["rsid"]],"_HET.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
  write.table(maptorun[iids2,], file=paste0(vardir, "/", variant[["rsid"]],"_HOMALT.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Count alleles
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
countAlleles = function(variant) {
  command = paste("bash",paste0(Sys.getenv("scripts"),"/fine_mapping/allelic_bias.sh"), variant[["rsid"]], variant[["chr"]], variant[["start"]], variant[["ref"]], variant[["alt"]])
  cat(command,"\n")
  system(command)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot allelic imbalance
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getAlleleImbalance = function(variant) {
  vardir = paste0(Sys.getenv("freeze"),"/allelic_bias/",variant[["rsid"]])
	d = read.table(paste(vardir, "allelic_bias_HET_summary.txt", sep="/"), header=TRUE)
  d$Sample = row.names(d)
  d$Total = d[,variant[["ref"]]]+d[,variant[["alt"]]]
  d$Prop = d[,variant[["alt"]]]/d$Total
  dlong = melt(d, id=c("Sample","Prop","Total"), value.name="Count", variable.name="Allele")

  #test for overall statistical significance
  #model1 <- glmer(Count ~ (1|Sample) + Allele, family = "poisson", data=dlong)
  #model2 <- glmer(Prop ~ (1|Sample) + Allele, weights=Total, family = "binomial", data=dlong)
  if (sum(dlong$Total>5)>5) {
    dlong_filtered = dlong[dlong$Total>5,]
    model1 <- glmer(Count ~ (1|Sample) + Allele, family = "poisson", data=dlong_filtered)
    model2 <- glmer(Count/Total ~ (1|Sample) + Allele, weights=Total, family = "binomial", data=dlong_filtered)
  } else {
    model1 = NA
    model2 = NA
  }
  return(list(variant=variant, d=d, dlong=dlong, modelPoisson = model1, modelBinomial = model2))
}

plotAlleleImbalance = function(aiobject) {
  #sample_order = dlong[order(d$Total), "Sample"]
  variant = aiobject[["variant"]]
  dlong = aiobject[["dlong"]]

  #plot sample-level data
  vardir = paste0(Sys.getenv("freeze"),"/allelic_bias/",variant[["rsid"]])
  pdf(paste0(vardir,"/allelic_imbalance.pdf"))
  dlong$CountPlot = apply(dlong, 1, function (x) { allele=x["Allele"]; count=as.numeric(x["Count"]); if (allele==variant[["ref"]]) {return(-count)} else {return(count)}})
  ymax = round(max(dlong$Count)/5)*5
  ybreaks = seq(-ymax, ymax, by=5)
  p = ggplot(dlong, aes(x = Sample, y = CountPlot, group = Allele, fill = Allele)) +
    geom_bar(stat = "identity", width = 0.75) +
    coord_flip() +
    scale_x_discrete(limits = sort(unique(dlong$Sample))) +
    scale_y_continuous(breaks = ybreaks,
                       labels = abs(ybreaks),
                     limits = c(min(ybreaks)-5, max(ybreaks)+5)) +
    labs(x = "Sample", y = "Read count", title = variant[["rsid"]]) +
    scale_fill_manual(values=c("red", "blue"),
                      name="",
                      breaks=c(variant[["ref"]], variant[["alt"]]),
                      labels=c(variant[["ref"]], variant[["alt"]])) +
    theme_cowplot() + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "top")
  print(p)
  dev.off()
  return(p)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Combine into pipeline
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prepAndCount = function(rsid, dat_sub) {
  cat("STARTING ", rsid, "\n")
  cat("getting variant info ...\n")
  variant = getVariantInfo(rsid)
  cat("getting sample info ...\n")
  getSampleInfo(variant, dat_sub)
  cat("counting alleles ...\n")
  countAlleles(variant)
  cat("FINISHED ", rsid, "\n")
}


summarizeAndPlot = function(rsid) {
  cat("STARTING ", rsid, "\n")
  cat("getting variant info ...\n")
  variant = getVariantInfo(rsid)
  cat("summarzing allele counts ...\n")
  aiobject = getAlleleImbalance(variant)
  cat("plotting allele counts ...\n")
  plotAlleleImbalance(aiobject)
  cat("FINISHED ", rsid, "\n")
  return(aiobject)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run pipeline
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Get snps in peaks to run pipeline on
coloc_table = readRDS(file=paste0(Sys.getenv("freeze"),"/caqtl_colocalisation_stats.rds"))
coloc_table = coloc_table[coloc_table$PP.H4.abf>0.8,]
snps_in_peak = unlist(strsplit(paste(coloc_table$snps_in_peak[!coloc_table$snps_in_peak==""], collapse=";"), split=";"))
#add rs7731626 as it is most likely causal variant in ANKRD55 locus
#snps_in_peak = c(snps_in_peak,"rs7731626")

### Run pipeline for each snp in peak
lapply(snps_in_peak[!snps_in_peak=="rs705705"], prepAndCount, dat_un)
alleleImbalance_out = lapply(snps_in_peak[!snps_in_peak=="rs705705"], summarizeAndPlot)
alleleImbalance_plots = lapply(alleleImbalance_out, plotAlleleImbalance)

#prepAndCount("rs11160429", dat_un)
#bash ${scripts}/fine_mapping/allelic_bias.sh rs11160429 chr14 98019103 G C

### Check to see if models converge
lapply(alleleImbalance_out, function(x) {summary(x$modelPoisson)})
lapply(alleleImbalance_out, function(x) {summary(x$modelBinomial)})
lapply(alleleImbalance_out, function(x) {if (is.na(x$modelPoisson)) {cat(x$variant$rsid,":","no model\n")} else {cat(x$variant$rsid,":",x$modelPoisson@optinfo$conv$lme4$messages,"\n")}})
lapply(alleleImbalance_out, function(x) {if (is.na(x$modelPoisson)) {cat(x$variant$rsid,":","no model\n")} else {cat(x$variant$rsid,":",x$modelBinomial@optinfo$conv$lme4$messages,"\n")}})

### Run pipeline for special cases
#run for rs705705 (this variant was not imputed in MEGA but is in near perfect LD with rs705704)
rs705705_varobj = list(rsid="rs705705", altid="chr12:56041720:G:C", chr="chr12", start=56041720, end=56041720, ref="G", alt="C")
countAlleles(rs705705_varobj)
rs705705_out = getAlleleImbalance(rs705705_varobj)
rs705705_plot = plotAlleleImbalance(rs705705_out)

#run for rs6908626
prepAndCount("rs72928038", dat_un)
summarizeAndPlot("rs72928038")

prepAndCount("rs6908626", dat_un)
summarizeAndPlot("rs6908626")

prepAndCount("rs7731626", dat_un)
summarizeAndPlot("rs7731626")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Collect data on allelic imbalance for table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## For each snp in OUTTABLE$snps_in_peak, we want
#   - tag_snp
#   - snp: rsid, markername, alleleA, alleleB, AF_alleleB
#   - qtl_stats: beta, pval, PP.H4.abf
#   - number_of_heterozygotes
#   - total_read_count: alleleA
#   - total_read_count: alleleB
#   - alleic_imabalance_stats: prop (alleleB/(alleleA+alleleB)), pvalue (poisson model)

snpToColocDict = list()
for (i in 1:nrow(coloc_table)) {
  snps = unlist(strsplit(coloc_table[i, "snps_in_peak"], split=";"))
  if (length(snps>0)) {
    for (j in 1:length(snps)) {
      for (k in 1:ncol(coloc_table)) {
        snpToColocDict[[snps[j]]][[names(coloc_table)[k]]] = coloc_table[i,names(coloc_table)[k]]
      }
    }
  }
}
snpToColoc = function(rsid, colName) {
  snpToColocDict[[rsid]][[colName]]
}

collectDataForTable = function(aiobject, results) {

  variant = aiobject[["variant"]]
  snp_rsid = variant[["rsid"]]
  snp_altid = variant[["altid"]]

  GUESSFM_tag = paste(finemap[!is.na(finemap$MarkerName) & finemap$MarkerName==snp_altid,c("tag_rsid","tag")], collapse=";")
  caqtl_index =
  peak = snpToColoc(snp_rsid, "peak")
  snp = paste(snp_rsid, snp_altid,sep=";")
  alleleA = variant[["ref"]]
  alleleB = variant[["alt"]]
  qtl_lead = gsub(":","_",unlist(strsplit(snpToColoc(snp_rsid, "qtl_lead"), split=";"))[2])
  qtl_dat = results_qtl_un[[qtl_lead]]
  qtl_stats = qtl_dat[gsub(":","_", snp_altid),c("Beta","pvalue")]
  names(qtl_stats) <- paste0("caQTL_", names(qtl_stats))
  qtl_coloc = snpToColoc(snp_rsid, "PP.H4.abf")
  no_hets = sum(aiobject$d$Total>0)
  total_ref = sum(aiobject$d[,alleleA])
  total_alt = sum(aiobject$d[,alleleB])
  prop_alt = total_alt/(total_ref+total_alt)
  if (class(aiobject$modelPoisson)=="glmerMod") {
    pval_poisson = summary(aiobject$modelPoisson)$coefficients[paste0("Allele",alleleB),"Pr(>|z|)"]
    pval_binomial = summary(aiobject$modelBinomial)$coefficients[paste0("Allele",alleleB),"Pr(>|z|)"]
  } else {
    pval_poisson = NA
    pval_binomial = NA
  }
  OUT = data.frame(snp_in_peak=snp,
    GUESSFM_tag=GUESSFM_tag, peak=peak, ref=alleleA, alt=alleleB,
    qtl_stats, PP.H4.abf = qtl_coloc,
    number_heterozygous=no_hets, REF=total_ref, ALT=total_alt, Prop=prop_alt, pvalue_AllelicBias=pval_poisson, pval_binomial)
  return(OUT)
}

OUTTABLE = do.call("rbind", lapply(alleleImbalance_out, collectDataForTable))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save to table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source(paste0(Sys.getenv("scripts"),"/manuscript/format_tables.R"))
library(rtf)
rtffile <- RTF("allelic_imbalance_tophits.doc")
addParagraph(rtffile, "Allelic bias in heterozygous individuals for GUESSFM-prioritised variants overlapping caQTL peaks")
addTable(rtffile, prettyTable(OUTTABLE[OUTTABLE$PP.H4.abf>0.8,]))
addParagraph(rtffile, "In all cases where the caQTL effect is positive (Beta>0, meaning the alternative allele correlates with greater accessibility),
we see more than 50 percent of reads from heterozygous individuals containing the alternative allele.
Likewise, we see that less than 50 percent of reads contain the alternative allele for caQTL variants with a negative effect on accessibility (Beta<0).")
done(rtffile)




# ### BACH2
# runAll("rs72928038")
# runAll("rs6908626")
#
# ### CENPW
# "rs9388486"
#
# ### ANKRD55
# #credible set
# "rs10213692"
# "rs6859219"
# "rs10065637"
# #negative controls (snps in peak)
# "rs10214316"
# "rs78213348"
# "rs12515884"
# "rs160924"
#
# ### RP11-6101.1
# "rs11628807"
# "rs4383076"
# "rs11628876"
# "rs11160429"
#
# ### IKZF4/RPS26
# "rs705704"
# "rs34415530"
