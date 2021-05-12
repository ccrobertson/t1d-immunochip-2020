library(ggplot2)
library(cowplot)

## Get Pi data object
load("/nv/vol185/MEGA/Jamies_files/pi/eqtlgen_score_20_dens.RData")


## Gather eQTL evidence across cell types
eGeneDataList = list()
k=0
for (i in 5:length(dTarget$list_pNode)) {
  k = k+1
  eGeneDataList[[k]] = dTarget$list_pNode[[i]]$evidence
}
#fix header name for eQTLGen data set
names(eGeneDataList[[12]]) <- names(eGeneDataList[[11]])
#combine into one data.frame
eGeneData = do.call("rbind", eGeneDataList)

## Plot eGene -- direction of effects
targetsDF = read.table("supp_table_17.txt", header=TRUE)
eGeneDataTargets = eGeneData[eGeneData$Symbol%in%targetsDF$Gene,]
eGeneDataTargets$sign_GWAS = sign(eGeneDataTargets$b_GWAS)
eGeneDataTargets$sign_eQTL = sign(eGeneDataTargets$b_eQTL)
eGeneDataTargets$relationship = eGeneDataTargets$sign_GWAS * eGeneDataTargets$sign_eQTL
eGeneDataTargets$celltype = factor(eGeneDataTargets$context,
  levels=c("Whole blood", "Blood","Neutrophil","CD14","Monocyte","LPS2","LPS24","IFN","NK","Bcell","CD4","CD8"),
  labels=c("Whole blood (Vosa)", "Blood (Westra)","Neutrophil","Macrophage (CD14)","Monocyte","Monocyte + LPS 2h","Monocyte + LPS 24h","Monocyte + IFN 24h","Natural Killer","B cell","CD4 T cell","CD8 T cell"))

## Plot eGene direction of effects by cell type
pdf("pi_directionality_plots.pdf")
ggplot(eGeneDataTargets, aes(x=celltype, y=Symbol)) +
  geom_tile(aes(fill=relationship*pp_ABF)) +
  scale_fill_gradient2(midpoint=0,low="blue",high="red", mid="white") +
  theme_cowplot() +
  theme(axis.ticks=element_blank(), axis.text.x = element_text(angle = 90), legend.position="bottom", legend.title=element_blank()) +
  xlab(" ") + ylab(" ")
dev.off()
