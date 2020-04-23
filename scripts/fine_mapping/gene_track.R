cat("loading biomaRt library\n")
library(biomaRt)

cat("loading Gviz library\n")
library(Gviz)

cat("loading ggplot2 library\n")
library(ggplot2)

cat("extracting Ensembl Homo sapiens data\n")
ensembl=useMart("ENSEMBL_MART_ENSEMBL")
datasets=listDatasets(ensembl)
human_ensembl = useDataset("hsapiens_gene_ensembl", mart=ensembl)

cat("defining geneTrack plotting function\n")
geneTrack = function(region, gene_list=NULL, colour=NULL) {

  #define colour
  if (is.null(colour)) {
    geneColor = "#999999"
  } else {
    geneColor = colour
  }

  #get exons annotations
  exonData = data.frame(geneTrack <- Gviz::BiomartGeneRegionTrack(biomart = human_ensembl, genome = "hg38", chromosome = as.numeric(gsub("chr","",region$chr)), start = region$start, end = region$end, name = "ENSEMBL", transcriptAnnotation="symbol")@range)
  exonData$exon_uniq = paste(exonData$exon, exonData$feature, sep=":")
  exonData$features_in_transcript = sapply(exonData$transcript, function(x) {sum(exonData$transcript==x)})
  exonData$longest_transcript = sapply(exonData$symbol, function(x) {d = exonData[exonData$symbol==x,c("transcript","features_in_transcript")]; d[order(d$features_in_transcript, decreasing=TRUE),]$transcript[1]})

  #subset exons for plotting
  exons = exonData[exonData$transcript==exonData$longest_transcript & exonData$feature %in% c("lncRNA","miRNA","protein_coding","utr3","utr5"),]
  exons$plotStart = as.numeric(sapply(exons$start, function(x) max(x,region[["start"]])))
  exons$plotEnd = sapply(exons$end, function(x) min(x,region[["end"]]))
  exons$width = exons$plotEnd - exons$plotStart
  exons$ymin = -0.25
  exons$ymax = 0.25
  exons$ymin[exons$feature=="protein_coding"] <- -0.5
  exons$ymax[exons$feature=="protein_coding"] <- 0.5

  #create genes data frame
  genes = data.frame(
    symbol=unique(exons$symbol),
    t(sapply(unique(exons$symbol), function(x){d=exons[exons$symbol==x,]; return(c(min(d$start, na.rm=TRUE), max(d$end, na.rm=TRUE)))})),
    strand=sapply(unique(exons$symbol), function(x){d=exons[exons$symbol==x,]; return(d$strand[1])}))
  names(genes) <- c("symbol","start","end", "strand")
  genes$plotStart = as.numeric(sapply(genes$start, function(x) max(x,region[["start"]])))
  genes$plotEnd = sapply(genes$end, function(x) min(x,region[["end"]]))
  genes$width = genes$plotEnd - genes$plotStart

  #filter based on boundaries or provided gene list
  if (is.null(gene_list)) {
    exons$include = exons$width>0
    genes$include = genes$width>0
  } else {
    exons$include = exons$symbol%in%gene_list & exons$width>0
    genes$include = genes$symbol%in%gene_list & genes$width>0
  }
  exons = exons[exons$include,]
  genes = genes[genes$include,]

  #create arrowtick data frame
  arrowtickpos = seq(region[['start']],region[['end']], length=20)
  arrows = NULL
  for (i in 1:dim(genes)[1]) {
    if (sum(arrowtickpos>genes$plotStart[i] & arrowtickpos<genes$plotEnd[i])>0) {
      d = data.frame(symbol=genes$symbol[i], strand=genes$strand[i],start=arrowtickpos[arrowtickpos>genes$plotStart[i] & arrowtickpos<genes$plotEnd[i]])
      if (nrow(d)>3) {
        d_trimmed = d[2:(nrow(d)-1),]
      } else {
        d_trimmed = d
      }
      arrows = rbind(arrows, d_trimmed)
    }
  }
  arrows$end[arrows$strand=="-"] = arrows$start[arrows$strand=="-"]-1
  arrows$end[arrows$strand=="+"] = arrows$start[arrows$strand=="+"]+1

  if(!is.null(arrows)) {
  ggplot() +
    #plot exons
    geom_rect(data=exons, aes(xmin=plotStart, xmax=plotEnd, ymin=ymin, ymax=ymax), colour=geneColor, fill=geneColor, alpha=1) +
    #connect exons
    geom_rect(data=genes, aes(xmin=plotStart, xmax=plotEnd, ymin=-0.005, ymax=0.005), colour=geneColor, fill=geneColor, alpha=1) +
    #add arrow ticks
    geom_segment(data=arrows, aes(x=start, xend=end, y=0, yend=0), lineend="butt", linejoin="mitre", arrow = arrow(angle=20, length = unit(0.5, "lines"), ends="last", type="closed"), colour=geneColor) +
    #tweak aesthetics
    scale_x_continuous(limits = c(region[["start"]],region[["end"]]), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1,1), expand = c(0, 0)) +
    facet_grid(symbol~., margins=FALSE, switch="y") +
    ylab(" ") + xlab(" ") +
    theme(
      strip.text.y = element_text(size = 8, colour = "black", angle = 180),
      strip.background = element_blank(),
      panel.spacing.y=unit(0,"lines"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid=element_line(size=0.5),
      panel.background = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      legend.position = "none")
  } else {
    ggplot() +
      #plot exons
      geom_rect(data=exons, aes(xmin=plotStart, xmax=plotEnd, ymin=ymin, ymax=ymax), colour=geneColor, fill=geneColor, alpha=1) +
      #connect exons
      geom_rect(data=genes, aes(xmin=plotStart, xmax=plotEnd, ymin=-0.005, ymax=0.005), colour=geneColor, fill=geneColor, alpha=1) +
      #tweak aesthetics
      scale_x_continuous(limits = c(region[["start"]],region[["end"]]), expand = c(0, 0)) +
      scale_y_continuous(limits = c(-1,1), expand = c(0, 0)) +
      facet_grid(symbol~., margins=FALSE, switch="y") +
      ylab(" ") + xlab(" ") +
      theme(
        strip.text.y = element_text(size = 8, colour = "black", angle = 180),
        strip.background = element_blank(),
        panel.spacing.y=unit(0,"lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid=element_line(size=0.5),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none")
  }
}



# pdf("gene_track.pdf", width=10, height=2)
# geneTrack(region=list(chr="chr10", start=5988475,end=6127208))
# geneTrack(region=list(chr="chr10", start=5988475,end=6127208), gene_list=c("IL2RA", "RBM17"), colour=cbbPalette[2])
# geneTrack(region=list(chr="chr2", start=162147475,end=162479787), colour=cbbPalette[2])
# geneTrack(region=list(chr="chr2", start=162147475,end=162479787), gene_list=c("IFIH1"))
# dev.off()
