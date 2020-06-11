library(ggplot2)
library(cowplot)
library(gridExtra)
library(RColorBrewer)
library(rtracklayer) #use this for importing bigwigs
library(locuscomparer)
library(reshape2)
library(lme4)

pubdir = Sys.getenv('PUBLIC_DATA')
setwd(Sys.getenv('finemap'))


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Input parameters
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
finemap = readRDS(file="annotated_credible_sets.rds")

getVariantInfo = function(rsid) {
  altid = finemap$MarkerName[finemap$MarkerID==rsid]
  params = unlist(strsplit(altid, split=":"))
  return(list(rsid=rsid, altid=altid, chr=params[1], start=params[2], end=params[2], ref=params[3], alt=params[4]))
}

variant1 = getVariantInfo("rs72928038")
variant2 = getVariantInfo("rs6908626")

#find peaks around variants
dp3_dir = Sys.getenv('freeze')
dp3_un = readRDS(file.path(dp3_dir,"DATASET_un.rds"))
dp3_counts = dp3_un[["counts"]]
dp3_geno = dp3_un[["genotypes"]]
dp3_sampleinfo = dp3_un[["samples"]]

peaks_all = data.frame(do.call("rbind", sapply(row.names(dp3_counts), function(x) {y = unlist(strsplit(x, split="_")); if (length(y)==3) { return(y) } } )))
names(peaks_all) <- c("chromosome","start","end")
peaks_all$start = as.numeric(peaks_all$start)
peaks_all$end = as.numeric(peaks_all$end)
peaks_all[peaks_all$chromosome=="chr6" & peaks_all$start<variant1[["start"]] & peaks_all$end>variant1[["start"]],]
peaks_all[peaks_all$chromosome=="chr6" & peaks_all$start<variant2[["start"]] & peaks_all$end>variant2[["start"]],]

peak1 = list(label="chr6_90266766_90267747", chr="chr6", start=90266766, end=90267747)
peak2 = list(label="chr6_90294665_90297341", chr="chr6", start=90294665, end=90297341)



region_broad = list(chr=variant1[["chr"]], start=min(finemap[finemap$tag_rsid==variant1[["rsid"]],"position"])-5000, end=max(finemap[finemap$tag_rsid==variant1[["rsid"]],"position"])+5000)
region_narrow1 = list(chr=peak1[["chr"]], start=peak1[["start"]]-2000, end=peak1[["end"]]+2000)
region_narrow2 = list(chr=peak2[["chr"]], start=peak2[["start"]]-2000, end=peak2[["end"]]+2000)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Source plotting functions
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source(paste0(Sys.getenv("scripts"),"/fine_mapping/plotting_functions.R"))
source(paste0(Sys.getenv("scripts"),"/fine_mapping/gene_track.R"))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get credible set data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
finemap = readRDS(file="annotated_credible_sets.rds")
finemap_dat = data.frame(chr=finemap$chr, start=finemap$position, end=finemap$position, ppsum=finemap$ppsum, tag=finemap$tag, tag_rsid=finemap$tag_rsid, marker=finemap$MarkerID, include_line=finemap$include_line)
finemap_reg = finemap_dat[finemap_dat$tag_rsid==variant1[["rsid"]],]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get Blueprint data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
blu_dat = readRDS(file="Blueprint_chromHMM_tracks_dataframe.rds")
blukey_tissue = readRDS(file="Blueprint_chromHMM_tracks_tissue_key.rds")
blukey_state = readRDS(file="Blueprint_chromHMM_tracks_state_key.rds")
blu_reg = extractBedRegion(blu_dat, region_broad)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get Calderon ATAC-seq
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stanford_dir = paste0(pubdir,"/Calderon/atac_data/processed_by_Flores/bigwigs")
stanford_samples =  c(
  "B_cells_S", "B_cells_U",  "Naive_B-S",  "Naive_B-U",  "Mem_B-S",  "Mem_B-U",
  "CD4_Teff_S",  "CD4_Teff_U",  "Naive_Teffs-S",  "Naive_Teffs-U",  "Memory_Teffs-S",  "Memory_Teffs-U",  "Memory_Tregs-S",  "Memory_Tregs-U",
  "Th1_precursors_S",  "Th1_precursors_U",  "Th2_precursors_S",  "Th2_precursors_U",  "Th17_precursors_S",  "Th17_precursors_U",  "Tfh-S",  "Tfh-U",
  "CD8pos_S",  "CD8pos_U",  "Naive_CD8_T-S",  "Naive_CD8_T-U",  "CD8cm_S",  "CD8cm_U",  "CD8em_S",  "CD8em_U",
  "Regulatory_T-S",  "Regulatory_T-U",
  "Gamma_delta_T-S", "Gamma_delta_T-U",
  "Immature_NK-U", "Mature_NK-S",  "Mature_NK-U",  "Memory_NK-U",
  "Monocytes-S",  "Monocytes-U",  "Myeloid_DCs-U",  "pDCs")
stanford_reg = extractAtacRegion(stanford_samples, stanford_dir, region_broad)
stanford_reg_sub = stanford_reg[stanford_reg$type %in% c("Naive_Teffs-U","Naive_CD8_T-U","Naive_Teffs-S","Naive_CD8_T-S"),]

pdf("bach2_stanford_atac.pdf")
getAtacPlot(stanford_reg[stanford_reg$type %in% c("CD4_Teff_U","CD4_Teff_S","Naive_Teffs-U","Naive_Teffs-S","CD8pos_U","CD8pos_S"),], region_narrow1, title=variant1[["rsid"]])
dev.off()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## PCHi-C data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Read in PCHi-C data from Javierre et al
contactDat = read.table(paste0(pubdir, "/Javierre_pchic/PCHiC_peak_matrix_cutoff5.tsv"), header=TRUE)

### Extract just lines where bait or target is in BACH2
contactDat_BACH2 = contactDat[(!is.na(contactDat$baitName) & contactDat$baitName=="BACH2") | (!is.na(contactDat$oeName) & contactDat$oeName=="BACH2"),]
contactDat_BACH2$ID = paste0("N",seq(1:nrow(contactDat_BACH2)))

### LiftOver to hg38
library("liftOver")
liftoverPCHicDat = function(dat) {
  chain_hg19_to_hg38 = import.chain(con=paste0(Sys.getenv("resources"),"/hg19ToHg38.over.chain"))
  granges_obj_hg19_baits = GRanges(seqnames = paste0("chr",dat$baitChr), ranges = IRanges(start=dat$baitStart, end=dat$baitEnd), ID=dat$ID)
  granges_obj_hg19_targets = GRanges(seqnames = paste0("chr",dat$oeChr), ranges = IRanges(start=dat$oeStart, end=dat$oeEnd), ID=dat$ID)
  lifted_baits = data.frame(liftOver(x=granges_obj_hg19_baits, chain=chain_hg19_to_hg38))
  lifted_targets = data.frame(liftOver(x=granges_obj_hg19_targets, chain=chain_hg19_to_hg38))
  lifted_coords = merge(lifted_baits, lifted_targets, by="ID",suffixes=c("Bait","Target"))
  lifted_coords_with_scores = merge(lifted_coords, dat[,names(dat)[!names(dat) %in% c("baitChr","baitStart","baitEnd","oeChr","oeStart","oeEnd")]], by="ID")
  return(lifted_coords_with_scores)
}
contactDat_BACH2_hg38 = liftoverPCHicDat(contactDat_BACH2)

### Get interactions overlapping variant
getPCHicOverlap = function(dat, variant) {
  var_pos = as.numeric(variant[["start"]])
  dat_overlap = dat[(dat$startBait<var_pos & dat$endBait>var_pos) | (dat$startTarget<var_pos & dat$endTarget>var_pos),]
  return(dat_overlap)
}
interactions_to_plot = getPCHicOverlap(contactDat_BACH2_hg38,variant1)
interactions_to_plot$midBait = mean(c(interactions_to_plot$startBait, interactions_to_plot$endBait))
interactions_to_plot$midTarget = mean(c(interactions_to_plot$startTarget, interactions_to_plot$endTarget))

### Plot bait, target,and interaction
getInteractionPlot = function(interactions_to_plot, region) {
  ggplot(interactions_to_plot) +
    geom_rect(aes(xmin=startBait,xmax=endBait, ymin=0, ymax=1)) +
    geom_rect(aes(xmin=startTarget,xmax=endTarget, ymin=0, ymax=1)) +
    scale_x_continuous(limits = c(region[["start"]],region[["end"]]), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0,5), expand = c(0, 0)) +
    geom_curve(
      aes(x=midTarget, xend=startBait, y=1, yend=1.2),
      curvature=-0.2, #negative = left-hand curves; positive = right-hand curves
      lineend="butt",
      arrow = arrow(angle=20, length = unit(0.5, "lines"), ends="last", type="closed")) +
      #arrow = arrow(length = unit(0.03, "npc")), size=10) +
    #theme_cowplot() +
    theme(
      plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       panel.border = element_blank(),
       panel.background = element_blank(),
      axis.title.x = element_blank(),
       axis.title.y = element_blank(),
       axis.text.x = element_blank(),
       axis.text.y = element_blank(),
       axis.ticks = element_blank(),
       axis.line = element_blank())
}
pdf("bach2_pchic.pdf")
getInteractionPlot(interactions_to_plot,region_broad)
dev.off()


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## caQTL plots
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
caQTLTrack = function(dp3_reg, region, variant, ylimit) {
   ggplot(dp3_reg) +
      geom_rect(data=dp3_reg, aes(xmin=start-1, xmax=end+1, ymin=0, ymax=score)) +
      scale_x_continuous(limits = c(region[["start"]],region[["end"]]), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, ylimit), expand = c(0,0)) +
      facet_grid(genotypeFactor~., margins=FALSE, switch="y") +
      labs(x = " ", y = " ", title = variant[["rsid"]]) +
      theme_cowplot(12) +
      theme(
        strip.text.y = element_text(size = 15, colour = "black", angle = 180),
        strip.background = element_blank(),
        panel.spacing.y=unit(0.1,"lines"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype="solid", colour="grey", size=0.3),
        panel.grid.minor = element_blank(),
        panel.background=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="none")
}

caQTLBoxPlot = function(dp3_countdat, title=NULL) {
  ggplot(dp3_countdat[!is.na(dp3_countdat$genotype),]) +
    geom_boxplot(aes(x=genotypeFactor, y=peak_count), fill="#999999", outlier.shape = NA) +
    #geom_dotplot(aes(x=genotypeFactor, y=peak_count, fill=genotypeFactor), binaxis='y', stackdir='center') +
    geom_jitter(aes(x=genotypeFactor, y=peak_count), height=0, width=0.1) +
    theme_cowplot(12) +
    ggtitle(title) +
    theme(axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.title.x=element_blank(),
      axis.text.x=element_text(size=15),
      axis.title.y=element_text(size=15),
      legend.position = "none",
      plot.title = element_text(hjust=0, face="bold")) +
    ylab("Peak accessibility")
}

tracks_dir = paste0(Sys.getenv("qtl"),"/tracks")
createDP3trackDat = function(variant, peak, region, ylimit) {
  dp3_samples = paste(variant[["rsid"]], c("homref","het","homalt"),"mean", sep="_")
  dp3_reg = extractAtacRegion(dp3_samples, paste(tracks_dir,variant[["rsid"]],sep="/"), region)
  variant[["genotypes"]] = c(paste0(variant[["ref"]],variant[["ref"]]), paste0(variant[["ref"]],variant[["alt"]]), paste0(variant[["alt"]],variant[["alt"]]))
  dp3_reg$genotypeFactor = factor(dp3_reg$type, levels=dp3_samples, labels=variant[["genotypes"]])
  dp3_reg$score[dp3_reg$score>ylimit] <- ylimit

  samples_to_keep = dp3_sampleinfo[dp3_sampleinfo$condition=="u",]
  dp3_countdat = data.frame(genotype=dp3_geno[samples_to_keep$iid,gsub(":","_",variant[["altid"]])],peak_count=dp3_counts[peak[["label"]],samples_to_keep$sample], condition=samples_to_keep$condition)
  dp3_countdat$genotypeFactor = factor(dp3_countdat$genotype, levels=c(0,1,2), labels=variant[["genotypes"]])

  p_caqtl_track = caQTLTrack(dp3_reg, region, variant, ylimit) + snp_lines
  p_caqtl_boxplot = caQTLBoxPlot(dp3_countdat, title=variant[["rsid"]])
  return(list(trackDat=dp3_reg, countDat=dp3_countdat, trackPlot=p_caqtl_track, countPlot=p_caqtl_boxplot))
}

# variant=variant1
# peak=peak1
# region=region_narrow1
# ylimit=0.4
# rm(variant, peak, region, ylimit)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Allelic imbalance
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### NOTE this function should be the same as in allelic_imbalance.R
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

### NOTE this plotting function has been modified to remove titles
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
    theme_cowplot() + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "top", plot.title=element_blank())
  print(p)
  dev.off()
  return(p)
}



aiPipeline = function(variant) {
  aiobject = getAlleleImbalance(variant)
  plotAlleleImbalance(aiobject)
}

pdf("bach2_allele_imbalance.pdf")
aiPipeline(variant1)
aiPipeline(variant2)
dev.off()

pdf("rs705704_allelic_imbalance.pdf")
getAlleleImbalance(getVariantInfo("rs705704"))
dev.off()


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Colocalisation
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# getLCPlots = function(peak, gene) {
#   gwas_fn = paste0(Sys.getenv("qtl"),"/locuscompare/lcformat_",peak,"_gwas.txt")
#   caqtl_fn = paste0(Sys.getenv("qtl"),"/locuscompare/lcformat_",peak,"_qtl.txt")
#   eqtl_fn = paste0(Sys.getenv("qtl"),"/locuscompare/lcformat_",gene,"_qtl.txt")
#   p1 = locuscompare(in_fn1 = gwas_fn, in_fn2 = caqtl_fn, title = 'GWAS -log10(P)', title2 = 'caQTL -log10(P)', combine=FALSE)
#   p2 = locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'GWAS -log10(P)', title2 = 'BACH2 eQTL -log10(P)', combine=FALSE)
#   return(list(p1, p2))
# }

getLCPlot = function(label, ylab) {
  gwas_fn = paste0(Sys.getenv("freeze"),"/locuscompare/lcformat_",label,"_gwas.txt")
  qtl_fn = paste0(Sys.getenv("freeze"),"/locuscompare/lcformat_",label,"_qtl.txt")
  locuscompare(in_fn1 = gwas_fn, in_fn2 = qtl_fn, title = 'GWAS', title2 = ylab, combine=FALSE)
}
pdf("bach2_coloc.pdf")
getLCPlot(label="chr6_90267049_G_A__chr6_90266766_90267747__un", ylab="chr6:90266766-90267747 caQTL")
getLCPlot(label="chr6_90267049_G_A__BACH2__EUR", ylab="BACH2 eQTL")
dev.off()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create plots
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
snp_lines = geom_vline(xintercept=finemap_reg[finemap_reg$include_line,"start"],linetype = "longdash", size=0.2, alpha=1)
p_credset = getCredSetPlot2(finemap_reg, region_broad)
p_gene = geneTrack(region_broad, gene_list=c("BACH2")) + snp_lines
p_hic = getInteractionPlot(interactions_to_plot,region_broad) + snp_lines
p_qtl_1 = createDP3trackDat(variant1, peak1, region_narrow1, ylimit=0.4)
p_qtl_2 = createDP3trackDat(variant2, peak2, region_narrow2, ylimit=1)
p_ai_1 = aiPipeline(variant1)
p_ai_2 = aiPipeline(variant2)
p_atac = getAtacPlot(stanford_reg[stanford_reg$type %in% c("CD4_Teff_U","CD4_Teff_S","Naive_Teffs-U","Naive_Teffs-S","CD8pos_U","CD8pos_S"),], region_narrow1, title=variant1[["rsid"]]) + snp_lines



pdf("bach2_plot.pdf", width=11, height=11)
p_main1 = plot_grid(p_credset,p_gene, getChromHMMPlot2(blu_reg, region_broad, blukey_tissue, blukey_state, snp_lines),p_hic,blankPlot, nrow=5, ncol=1, align="v", axis="lr", rel_heights=c(1,1,4,1,0.5))
p_main2 = plot_grid(p_qtl_1[["countPlot"]], p_qtl_2[["countPlot"]],p_ai_1, p_ai_2, nrow=2, ncol=2, align="h", axis="b")
p_main3 = plot_grid(p_atac, blankPlot, nrow=2, rel_heights=c(10,1))
p_main4 = plot_grid(
  getLCPlot(label="chr6_90267049_G_A__chr6_90266766_90267747__un", ylab="caQTL")[[1]],
  getLCPlot(label="chr6_90267049_G_A__BACH2__EUR", ylab="BACH2 eQTL")[[1]],
  nrow=1, ncol=2, align="h", axis="b")
grid.arrange(p_main1, p_main2, p_main3, p_main4, ncol=2, nrow=4, layout_matrix = cbind(c(1,1,2,2), c(3,3,3,4)))
dev.off()


pdf("bach2_chromhmm.pdf")
plot_grid(getChromHMMPlot2(blu_reg, region_broad, blukey_tissue, blukey_state, snp_lines))
dev.off()
# pdf("bach2_plot.pdf", width=11, height=11)
# p_main1 = plot_grid(p_credset,p_gene,getChromHMMPlot2(blu_reg, region_broad, blukey_tissue, blukey_state, snp_lines),p_atac, nrow=4, ncol=1, align="v", axis="lr", rel_heights=c(1,1,3,6))
# p_main2 = plot_grid(p_qtl_1[["trackPlot"]],p_qtl_2[["trackPlot"]], p_qtl_1[["countPlot"]], p_qtl_2[["countPlot"]], nrow=2, ncol=2, align="h", axis="b")
# p_main3 = plot_grid(p_ai_1, p_ai_2, nrow=1, ncol=2, align="h")
# p_main4 = plot_grid(getLCPlot(peak1[["label"]], "caQTL")[[1]], getLCPlot("BACH2", "BACH2 eQTL")[[1]], nrow=1, ncol=2, align="h", axis="b")
# #plot_grid(p_main1, p_main2, p_main3, p_main4, nrow=1, ncol=4)
# grid.arrange(p_main1, p_main2, p_main3, p_main4, ncol=2, nrow=3, layout_matrix = cbind(c(1,1,1,1), c(2,2,3,4)))
# dev.off()
#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Look up for ANKRD55 region
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
variant1 = list(rsid="rs7731626", altid="chr5:56148856:G:A", chr="chr5", start=56148856, end=56148856, ref="G", alt="A")
peak1 = list(label="chr5_56147972_56149111", chr="chr5", start=56147972, end=56149111)
region_narrow1 = list(chr=peak1[["chr"]], start=peak1[["start"]]-8000, end=peak1[["end"]]+1000)
region_broad = region_narrow1
stanford_reg = extractAtacRegion(stanford_samples, stanford_dir, region_broad)
stanford_reg_sub = stanford_reg[stanford_reg$type %in% c("Naive_Teffs-U","Naive_CD8_T-U","Naive_Teffs-S","Naive_CD8_T-S"),]
snp_positions = c(finemap_dat[finemap_dat$tag_rsid=="rs71624119","start"], variant1[["start"]])
snp_lines = geom_vline(xintercept=snp_positions,linetype = "longdash", size=0.2, alpha=1)
p_atac = getAtacPlot(stanford_reg[stanford_reg$type %in% c("CD4_Teff_U","CD4_Teff_S","Naive_Teffs-U","Naive_Teffs-S","CD8pos_U","CD8pos_S"),], region_narrow1, title=" ") +
  snp_lines +
  annotate("rect", xmin=peak1[["start"]], xmax=peak1[["end"]], ymin=0, ymax=Inf, alpha=0.2, fill="red")

pdf("ANKRD55_rs7731626_vs_rs71624119_stanford_atac.pdf")
print(p_atac)
dev.off()
