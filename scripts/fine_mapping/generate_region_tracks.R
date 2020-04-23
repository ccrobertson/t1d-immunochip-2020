library(ggplot2)
library(cowplot)
library(rtracklayer)
library(RColorBrewer)

pubdir = Sys.getenv('PUBLIC_DATA')
setwd(Sys.getenv('finemap'))
cbbPalette <- c("#000000", "#E69F00", "#999999", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Data pre-processing:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#source(paste0(Sys.getenv("scripts"),"/fine_mapping/plotting_data_processing.R"))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get association stats
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load(paste0(Sys.getenv('meta'),"/resultsmeta_ind_fam_tdt_5pc_diff.RData"))
res$chr = paste0("chr",res$chromosome)
mega_dat = data.frame(chr=res$chr, start=res$position, end=res$position, P=res$P.value)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get credible sets
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
finemap = readRDS(file="annotated_credible_sets.rds")
paintor_vars = c(
  "rs10517086","rs34185821","rs35944082",
  "rs12969657","rs1865761","rs17207042",
  "rs67126256","rs2045258","rs9388486")
finemap$include_paintor_line = finemap$ID %in% paintor_vars

finemap_dat = data.frame(chr=finemap$chr, start=finemap$position, end=finemap$position, ppsum=finemap$ppsum, tag=finemap$tag, tag_rsid=finemap$tag_rsid, include_line=finemap$include_line, include_paintor_line=finemap$include_paintor_line)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get ATAC-seq file paths
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stanford_dir = paste0(pubdir,"/Calderon/atac_data/processed_by_Flores/bigwigs")
stanford_samples =  c(
  "B_cells_S", "B_cells_U",  "Naive_B-S",  "Naive_B-U",  "Mem_B-S",  "Mem_B-U",
  "CD4_Teff_S",  "CD4_Teff_U",  "Naive_Teffs-S",  "Naive_Teffs-U",  "Memory_Teffs-S",  "Memory_Teffs-U",
  "Th1_precursors_S",  "Th1_precursors_U",  "Th2_precursors_S",  "Th2_precursors_U",  "Th17_precursors_S",  "Th17_precursors_U",  "Tfh-S",  "Tfh-U",
  "CD8pos_S",  "CD8pos_U",  "Naive_CD8_T-S",  "Naive_CD8_T-U",  "CD8cm_S",  "CD8cm_U",  "CD8em_S",  "CD8em_U",
  "Regulatory_T-S",  "Regulatory_T-U", "Memory_Tregs-S",  "Memory_Tregs-U",
  "Gamma_delta_T-S", "Gamma_delta_T-U",
  "Immature_NK-U", "Mature_NK-S",  "Mature_NK-U",  "Memory_NK-U",
  "Monocytes-S",  "Monocytes-U",  "Myeloid_DCs-U",  "pDCs")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get EncodeRoadmap data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
enc_dat = readRDS(file="EncodeRoadmap_chromHMM_tracks_dataframe.rds")
enckey_tissue  = readRDS(file="EncodeRoadmap_chromHMM_tracks_tissue_key.rds")
enckey_state = readRDS(file="EncodeRoadmap_chromHMM_tracks_state_key.rds")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get Blueprint data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
blu_dat = readRDS(file="Blueprint_chromHMM_tracks_dataframe.rds")
blukey_tissue = readRDS(file="Blueprint_chromHMM_tracks_tissue_key.rds")
blukey_state = readRDS(file="Blueprint_chromHMM_tracks_state_key.rds")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Source plotting functions
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source(paste0(Sys.getenv("scripts"),"/fine_mapping/plotting_functions.R"))
source(paste0(Sys.getenv("scripts"),"/fine_mapping/gene_track.R"))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Main plotting function
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plotRegion_main = function(region, paintor=FALSE) {

  cat(region[["chr"]],":",region[["start"]],"-",region[["end"]],"\n")

  cat("extracting data for region...\n")
  enc_reg = extractBedRegion(enc_dat, region)
  blu_reg = extractBedRegion(blu_dat, region)
  mega_reg = extractBedRegion(mega_dat, region)
  finemap_reg = extractBedRegion(finemap_dat, region)
  stanford_reg = extractAtacRegion(stanford_samples, stanford_dir, region)

  if (paintor) {
    snp_lines = geom_vline(xintercept=finemap_reg[finemap_reg$include_paintor_line,"start"],linetype = "longdash", size=0.2, alpha=1, colour="red")
  } else {
    snp_lines = geom_vline(xintercept=finemap_reg[finemap_reg$include_line,"start"],linetype = "longdash", size=0.1, alpha=1)
  }

  cat("generating plots... ")
  cat(" p1")
  p1 = getCredSetPlot(finemap_reg, region) + snp_lines
  cat(" p2")
  p2 = geneTrack(region, colour=cbbPalette[2]) + snp_lines
  cat(" p3")
  #p3 = getChromHMMPlot(enc_reg, region, enckey_tissue, enckey_state, snp_lines)
  cat(" p4")
  #p4 = getChromHMMPlot(blu_reg, region, blukey_tissue, blukey_state, snp_lines)
  cat(" p5\n")
  p5 = getAtacPlot(stanford_reg, region) + snp_lines

  cat("combining and printing plots...\n")
  if (paintor) {
    pdffile = paste0("region_plots/",region[["chr"]],"_",region[["start"]],"_",region[["end"]],"_paintor.pdf")
  } else {
    pdffile = paste0("region_plots/",region[["chr"]],"_",region[["start"]],"_",region[["end"]],".pdf")
  }

  pdf(pdffile, width=8.5, height=11)
  p_chmm = plot_grid(p1,p2,getChromHMMPlot(enc_reg, region, enckey_tissue, enckey_state, snp_lines),getChromHMMPlot(blu_reg, region, blukey_tissue, blukey_state, snp_lines), nrow=4, ncol=1, align="v", axis="lr", rel_heights=c(1.2+0.8*length(unique(finemap_reg$tag)),2, 5,7))
  print(p_chmm)
  p_atac = plot_grid(p1,p2,p5, nrow=3, ncol=1, align="v", axis="lr", rel_heights=c(1.2+0.8*length(unique(finemap_reg$tag)),2,12))
  print(p_atac)
  dev.off()

}

plotRegion_main(region = list(chr="chr10", start=5988475,end=6127208))
plotRegion_main(region = list(chr="chr5", start=56140058, end=56148447))
plotRegion_main(region = list(chr="chr6",start=90250000,end=90310000))
plotRegion_main(region = list(chr="chr6",start=126327897,end=126527459))



### Generate plots with PAINTOR lines
plotRegion_main(region = list(chr="chr4", start=26073858, end=26137088), paintor=TRUE)
plotRegion_main(region = list(chr="chr6", start=126327897, end=126527459), paintor=TRUE)
plotRegion_main(region = list(chr="chr18", start=69836569, end=69886452), paintor=TRUE)

### Generate plots with caQTL peaks


### Generate custom plots in regions of interest
plotRegion_main(region = list(chr="chr17",start=39747478-2000,end=39924612+2000))
plotRegion_main(region = list(chr="chr17",start=40598769-1000,end=40657572+1000))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Plot all regions
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
regions = lapply(unique(finemap$region), function(x) {d=finemap[finemap$region==x,]; return(list(chr=paste0("chr",d$chromosome[1]), start=min(d$position)-10000, end=max(d$position)+10000))} )
lapply(regions, plotRegion_main)
