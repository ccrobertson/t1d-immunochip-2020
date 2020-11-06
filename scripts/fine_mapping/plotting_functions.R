##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Functions for creating plots in generate_region_plots.R and bach2_plots.R
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Blank
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
blankPlot <- ggplot()+geom_blank(aes(1,1))+
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


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## General function for extracting regions from data frames where first 3 cols are bed format
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extractBedRegion = function(dat, region) {
  dat_reg = dat[dat$chr==region[["chr"]] & dat$end>=region[["start"]] & dat$start<=region[["end"]],]
  dat_reg$plotStart = sapply(dat_reg$start, function(x) max(x,region[["start"]]))
  dat_reg$plotEnd = sapply(dat_reg$end, function(x) min(x,region[["end"]]))
  return(dat_reg)
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ATAC-seq plotting function
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

getAtacTrack = function(sample, dir, region) {
  regionGR <- GRanges(c(region[["chr"]]), IRanges(start=region[["start"]], end=region[["end"]]))
  bw = file.path(dir, paste0(sample,".bw"))
  bw_reg = data.frame(import.bw(bw, which = regionGR))
  names(bw_reg) <- c("chr", "start","end","width","strand","score")
  bw_reg$type = sample
  return(bw_reg)
}

extractAtacRegion = function(samples, dir, region) {
  tracks = lapply(samples, getAtacTrack, dir, region)
  dat_reg = do.call("rbind",tracks)
  dat_reg$type = factor(dat_reg$type, levels=samples)
  dat_reg$type_Num = as.numeric(dat_reg$type)
  return(dat_reg)
}


getAtacPlot = function(dat_reg, region, title=NULL) {
  B_types = c("B_cells_S", "B_cells_U",  "Naive_B-S",  "Naive_B-U",  "Mem_B-S",  "Mem_B-U")
  CD8_types = c("CD8pos_S",  "CD8pos_U",  "Naive_CD8_T-S",  "Naive_CD8_T-U",  "CD8cm_S",  "CD8cm_U",  "CD8em_S",  "CD8em_U")
  CD4_types = c("CD4_Teff_S",  "CD4_Teff_U",  "Naive_Teffs-S",  "Naive_Teffs-U",  "Memory_Teffs-S",  "Memory_Teffs-U","Th1_precursors_S",  "Th1_precursors_U",  "Th2_precursors_S",  "Th2_precursors_U",  "Th17_precursors_S",  "Th17_precursors_U",  "Tfh-S",  "Tfh-U")
  Regulatory_T_types = c("Memory_Tregs-S", "Memory_Tregs-U","Regulatory_T-S",  "Regulatory_T-U")
  GD_T_types = c("Gamma_delta_T-S","Gamma_delta_T-U")
  NK_types = c("Immature_NK-U", "Mature_NK-S",  "Mature_NK-U",  "Memory_NK-U")
  innate_types = c("Monocytes-S","Monocytes-U","Myeloid_DCs-U", "pDCs")
  atac_cell_types_ordered = c(B_types, CD8_types, CD4_types, Regulatory_T_types, GD_T_types, NK_types, innate_types)
  dat_reg$plotType = factor(dat_reg$type, levels=atac_cell_types_ordered)
  dat_reg$plotCat[dat_reg$type %in% B_types] <- "B"
  dat_reg$plotCat[dat_reg$type %in% CD8_types] <- "CD8"
  dat_reg$plotCat[dat_reg$type %in% CD4_types] <- "CD4"
  dat_reg$plotCat[dat_reg$type %in% Regulatory_T_types] <- "Regulatory"
  dat_reg$plotCat[dat_reg$type %in% GD_T_types] <- "Gamma-delta"
  dat_reg$plotCat[dat_reg$type %in% NK_types] <- "NK"
  dat_reg$plotCat[dat_reg$type %in% innate_types] <- "Innate"
  atac_cell_categories_ordered = c("B", "CD8",  "CD4",  "Regulatory", "Gamma-delta",  "NK", "Innate")
  atac_cell_categories_colors = c(brewer.pal(9, "Reds")[5],  #B
                                  brewer.pal(9, "Blues")[6],#CD8
                                  brewer.pal(9, "Blues")[8], #CD4
                                  brewer.pal(9, "Greys")[5], #Treg
                                  brewer.pal(9, "Greys")[5], #Gamma-delta
                                  brewer.pal(9, "Purples")[6],  #NK
                                  brewer.pal(9,"YlOrRd")[5]) #Innate
  dat_reg$plotCat = factor(dat_reg$plotCat, levels=atac_cell_categories_ordered)
  dat_reg$score[dat_reg$score>100] <- 100
  chrpos_ticks_breaks = seq(from = region[["start"]], to = region[["end"]], by = round((region[["end"]]-region[["start"]])/4, digits=0))
  chrpos_ticks_labels = paste0(round(chrpos_ticks_breaks/1e3, digits=2),"Kb")
  ggplot(dat_reg) + geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=score, fill=plotCat)) +
    scale_x_continuous(limits = c(region[["start"]],region[["end"]]), expand = c(0, 0), breaks=chrpos_ticks_breaks, labels=chrpos_ticks_labels) +
    #scale_y_continuous(limits = c(0, 1.25*max(dat_reg$score)), expand = c(0,0)) +
    scale_fill_manual(values=atac_cell_categories_colors[atac_cell_categories_ordered %in% unique(dat_reg$plotCat)]) +
    facet_grid(plotType~., margins=FALSE, switch="y") +
    ggtitle(title) +
    theme(
      strip.text.y = element_text(size = 8, colour = "black", angle = 180),
      strip.background = element_blank(),
      panel.spacing.y = unit(0.1,"lines"),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_text(size=8),
      #axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),axis.ticks.y=element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust=0.5, face="bold"))
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## chromHMM plotting function
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getChromHMMPlot = function(dat_reg, region, tissue_key, state_key, snp_lines) {

  row.names(state_key) = state_key$code
  dat_reg$state_label = state_key[dat_reg$state,"label"]
  state_key_nodups = state_key[!duplicated(state_key$label) & state_key$label %in% unique(dat_reg$state_label),]
  dat_reg$plotState = factor(dat_reg$state_label, levels=state_key_nodups$label, labels=state_key_nodups$label)
  dat_reg$plotTissue = factor(dat_reg$tissue, levels=levels(tissue_key$category))

  p = ggplot() +
    geom_rect(data=dat_reg, aes(xmin=plotStart,xmax=plotEnd, ymin=(EID_Num-1), ymax=EID_Num, fill=plotState)) +
    scale_x_continuous(limits = c(region[["start"]],region[["end"]]), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(breaks=state_key_nodups$label, values=state_key_nodups$color_hex) +
    facet_grid(plotTissue~., scales="free_y", space="free_y", margins=FALSE, switch="y") +
    #facet_grid(plotTissue~., scales="free_y", margins=FALSE, switch="y") +
    theme(
      #strip.background = element_rect(),
      strip.text.y = element_text(size = 8, colour = "black", angle = 180),
      panel.spacing.y = unit(0,"inches"),
      panel.background = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      legend.position = "none") +
    snp_lines

  g <- ggplot_gtable(ggplot_build(p))
  strip_both <- which(grepl('strip-', g$layout$name))
  tissue_colors = levels(tissue_key$COLOR)
  k <- 1
  for (i in strip_both) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- tissue_colors[k]
    k <- k+1
  }
  g
}
# pdf("chrom_hmm_track.pdf")
# grid.draw(getChromHMMPlot(enc_reg, region, enckey_tissue, enckey_state))
# grid.draw(getChromHMMPlot(blu_reg, region, blukey_tissue, blukey_state))
# ggplot() + geom_rect(data=blukey_tissue, aes(xmin=0, xmax=1,ymin=0,ymax=1, fill=category)) + facet_wrap(~category) +
# scale_fill_manual(values=levels(blukey_tissue$COLOR))
# dev.off()



getChromHMMPlot2 = function(dat_reg, region, tissue_key, state_key, snp_lines) {

  row.names(state_key) = state_key$code
  dat_reg$state_label = state_key[dat_reg$state,"label"]
  state_key_nodups = state_key[!duplicated(state_key$label) & state_key$label %in% unique(dat_reg$state_label),]
  dat_reg$plotState = factor(dat_reg$state_label, levels=state_key_nodups$label, labels=state_key_nodups$label)
  dat_reg$plotTissue = factor(dat_reg$tissue, levels=levels(tissue_key$category))
  chrpos_ticks_breaks = seq(from = region[["start"]], to = region[["end"]], by = round((region[["end"]]-region[["start"]])/4, digits=0))
  chrpos_ticks_labels = paste0(round(chrpos_ticks_breaks/1e6, digits=3),"Mb")
  p = ggplot() +
    geom_rect(data=dat_reg, aes(xmin=plotStart,xmax=plotEnd, ymin=(EID_Num-1), ymax=EID_Num, fill=plotState)) +
    #scale_x_continuous(limits = c(region[["start"]],region[["end"]]), expand = c(0, 0)) +
    scale_x_continuous(limits = c(region[["start"]],region[["end"]]), expand = c(0, 0), breaks=chrpos_ticks_breaks, labels=chrpos_ticks_labels) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(breaks=state_key_nodups$label, values=state_key_nodups$color_hex) +
    facet_grid(plotTissue~., scales="free_y", margins=FALSE, switch="y") +
    #facet_grid(plotTissue~., scales="free_y", margins=FALSE, switch="y") +
    theme(
      #strip.background = element_rect(),
      strip.text.y = element_text(size = 8, colour = "black", angle = 180),
      panel.spacing.y = unit(0,"inches"),
      panel.background = element_blank(),
      axis.text.x = element_text(size=8),
      #axis.ticks.x = element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      legend.position = "none") +
    snp_lines

  g <- ggplot_gtable(ggplot_build(p))
  strip_both <- which(grepl('strip-', g$layout$name))
  tissue_colors = levels(tissue_key$COLOR)
  k <- 1
  for (i in strip_both) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- tissue_colors[k]
    k <- k+1
  }
  g
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Credible set plotting function
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getCredSetPlot = function(dat_reg, region) {
  dat_reg$ppsum = round(dat_reg$ppsum, digits=3)
  tags = unique(dat_reg[order(dat_reg$ppsum),"tag"])
  ppsums = unique(dat_reg[order(dat_reg$ppsum),"ppsum"])
  bandsize= (region[["end"]]-region[["start"]])/500
  dat_reg$tag_num = as.numeric(factor(dat_reg$tag, levels=tags))
  ggplot(dat_reg) + geom_rect(aes(xmin=start-bandsize, xmax=end+bandsize, ymin=0.25, ymax=0.75, fill=ppsum)) +
    facet_grid(tag_rsid~., scales="free_y", space="free_y", margins=FALSE, switch="y") +
    scale_x_continuous(limits = c(region[["start"]],region[["end"]]), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
    scale_fill_gradient(limits=c(0,1), name="Credible Set Posterior Probability", low="white", high="black",
      guide=guide_colourbar(barheight=0.5, barwidth=5, title.position="left", direction = "horizontal", ticks=FALSE, title.theme=element_text(size=7), label.theme=element_text(size=5))) +
    theme(
      strip.text.y = element_text(size = 8, colour = "black", angle = 180),
      strip.background = element_blank(),
      panel.border = element_blank(),
      panel.spacing.y = unit(0,"inches"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background=element_blank(),
      axis.text.x = element_blank(), axis.ticks.x=element_blank(),
      axis.text.y = element_blank(), axis.ticks.y=element_blank(),
      legend.position = "top")
}

getCredSetPlot2 = function(dat_reg, region) {
  dat_reg$ppsum = round(dat_reg$ppsum, digits=3)
  tags = unique(dat_reg[order(dat_reg$ppsum),"tag"])
  ppsums = unique(dat_reg[order(dat_reg$ppsum),"ppsum"])
  bandsize= (region[["end"]]-region[["start"]])/500
  dat_reg$tag_num = as.numeric(factor(dat_reg$tag, levels=tags))
  ggplot(dat_reg) + geom_rect(aes(xmin=start-bandsize, xmax=end+bandsize, ymin=0.25, ymax=0.75, fill=ppsum)) +
  geom_text(aes(x=end, y=1, label=marker), hjust="middle", vjust="bottom", size=3) +
    facet_grid(tag_rsid~., scales="free_y", space="free_y", margins=FALSE, switch="y") +
    scale_x_continuous(limits = c(region[["start"]],region[["end"]]), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0,1.5), expand = c(0, 0)) +
    scale_fill_gradient(limits=c(0,1), name="Credible Set Posterior Probability", low="white", high="black",
      guide=guide_colourbar(barheight=0.5, barwidth=5, title.position="left", direction = "horizontal", ticks=FALSE, title.theme=element_text(size=7), label.theme=element_text(size=5))) +
    theme(
      strip.text.y = element_blank(),
      strip.background = element_blank(),
      panel.border = element_blank(),
      panel.spacing.y = unit(0,"inches"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background=element_blank(),
      axis.text.x = element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(),
      axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
      legend.position = "none")
}
#
# pdf("bach2_plot.pdf", width=11, height=11)
# getCredSetPlot2(finemap_reg, region_broad)
# dev.off()
