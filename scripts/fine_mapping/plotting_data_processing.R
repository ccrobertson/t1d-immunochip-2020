setwd(Sys.getenv('finemap'))
pubdir = Sys.getenv('PUBLIC_DATA')
#source()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Credible set data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
finemap = read.table("credible_sets_with_proxies_from_Jamie_ANNOTATED.txt", header=TRUE, comment.char="~", sep="\t")
finemap$chr = paste0("chr", finemap$chromosome)
finemap$MarkerID = finemap$ID
finemap$MarkerID[finemap$MarkerID=="."] <- finemap$MarkerName[finemap$MarkerID=="."]
finemap$include_line=(finemap$ppsum>0.8 & finemap$credset_size<=5)
saveRDS(finemap, file="annotated_credible_sets.rds")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## EncodeRoadmap data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
encdir = paste0(pubdir,"/EncodeRoadMap")
enckey_tissue = read.table(paste0(encdir,"/Roadmap.metadata.qc.jul2013_EncodeRoadmap_for_MEGA_figures.txt"), header=TRUE, sep="\t", comment.char="~")
enckey_tissue$category = factor(enckey_tissue$Category_MEGA, levels=unique(enckey_tissue$Category_MEGA))
getPalette = colorRampPalette(brewer.pal(12, "Set3"))
enc_tissue_colors = getPalette(length(levels(enckey_tissue$category)))
enckey_tissue$COLOR = factor(enckey_tissue$Category_MEGA, levels=levels(enckey_tissue$category), labels=enc_tissue_colors)

enckey_state_labels = read.table(paste0(encdir,"/labelmap_15_coreMarks.tab"), header=FALSE)
names(enckey_state_labels) <- c("code","label")

enckey_state_colors = read.table(paste0(encdir,"/colormap_15_coreMarks.tab"), header=FALSE)
names(enckey_state_colors) <- c("code","color_rgb")
enckey_state_colors$code = paste0("E",enckey_state_colors$code)

enckey_state = merge(enckey_state_labels, enckey_state_colors, by="code")
row.names(enckey_state) = enckey_state$label

rgb2hex <- function(x) {
  vec = as.numeric(unlist(strsplit(x, split=",")))
  rgb(vec[1],vec[2],vec[3], maxColorValue = 255)
}
enckey_state$color_hex = sapply(enckey_state$color_rgb, rgb2hex)

getEncodeTrack = function(eid) {
  bed = read.table(paste0(encdir,"/",eid,"_15_coreMarks_hg38lift_segments.bed.gz"), header=FALSE)
  names(bed) <- c("chr","start","end","state")
  bed$EID = eid
  bed$tissue = enckey_tissue[enckey_tissue$EID==eid,"category"]
  return(bed)
}

eids = enckey_tissue$EID
enctracks = lapply(eids, getEncodeTrack)
enc_dat_tmp = do.call("rbind", enctracks)
enc_dat = enc_dat_tmp[!enc_dat_tmp$tissue=="Other",]
enc_dat$EID_Num = as.numeric(factor(enc_dat$EID, levels=enckey_tissue$EID))
saveRDS(enc_dat, file="EncodeRoadmap_chromHMM_tracks_dataframe.rds")
saveRDS(enckey_tissue, file="EncodeRoadmap_chromHMM_tracks_tissue_key.rds")
saveRDS(enckey_state, file="EncodeRoadmap_chromHMM_tracks_state_key.rds")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Blueprint data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bludir = paste0(pubdir,"/Blueprint/hg38")
blukey_tissue = read.table(paste0(bludir,"/sample_key.txt"), sep="\t", header=TRUE)
blukey_tissue = blukey_tissue[order(blukey_tissue$major_catNum, blukey_tissue$category, decreasing=TRUE),]

blukey_state = read.table(paste0(bludir,"/state_key_reformatted.txt"))
names(blukey_state) <- c("code","label_12state")
blukey_state$label = NA
blukey_state$label[blukey_state$label_12state=="Repressed"] <- "ReprPC"
blukey_state$label[blukey_state$label_12state=="Repressed_Polycomb_Low_signal_H3K27me3"] <- "ReprPCWk"
blukey_state$label[blukey_state$label_12state=="Low_signal"] <- "Quies"
blukey_state$label[blukey_state$label_12state=="Heterochromatin_High_Signal_H3K9me3"] <- "Het"
blukey_state$label[blukey_state$label_12state=="Transcription_High_signal_H3K36me3"] <- "Tx"
blukey_state$label[blukey_state$label_12state=="Transcription_Low_signal_H3K36me3"] <- "TxWk"
blukey_state$label[blukey_state$label_12state=="Genic_Enhancer_High_Signal_H3K4me1_&_H3K36me3"] <- "EnhG"
blukey_state$label[blukey_state$label_12state=="Enhancer_High_Signal_H3K4me1"] <- "Enh"
blukey_state$label[blukey_state$label_12state=="Active_Enhancer_High_Signal_H3K4me1_&_H3K27Ac"] <- "Enh"
blukey_state$label[blukey_state$label_12state=="Distal_Active_Promoter_(2Kb)_High_Signal_H3K4me3_&_H3K27Ac_&_H3K4me1"] <- "TssAFlnk"
blukey_state$label[blukey_state$label_12state=="Active_TSS_High_Signal_H3K4me3_&_H3K4me1"] <- "TssA"
blukey_state$label[blukey_state$label_12state=="Active_TSS_High_Signal_H3K4me3_&_H3K27Ac"] <- "TssA"

blukey_state$code_15state = enckey_state[blukey_state$label,"code"]
blukey_state$color_hex = enckey_state[blukey_state$label,"color_hex"]

getBluTrack = function(eid) {
  bed = read.table(paste0(bludir,"/",eid,"_12_12_Blueprint_release_201608_segments.bed"), header=FALSE)
  names(bed) <- c("chr","start","end","state")
  bed$EID = eid
  bed$tissue = blukey_tissue[blukey_tissue$EID==eid,"category"]
  return(bed)
}
blu_eids = blukey_tissue$EID
blutracks = lapply(blu_eids, getBluTrack)
blu_dat_tmp = do.call("rbind", blutracks)
blu_dat_tmp$EID_Num = as.numeric(factor(blu_dat_tmp$EID, levels=blukey_tissue$EID))

#drop tissues we won't use
tissues_to_drop = c("Endothelial cell","Endothelial progenitor","Erythroblast","Megakaryocyte","Mesenchymal stem cell","Osteoclast")
samples_to_drop = blukey_tissue[blukey_tissue$label %in% c("band_form_neutrophil","neutrophilic_metamyelocyte","neutrophilic_myelocyte","segmented_neutrophil_of_bone_marrow"),"EID"]
blu_dat = blu_dat_tmp[!blu_dat_tmp$tissue %in% tissues_to_drop & !blu_dat_tmp$EID %in% samples_to_drop,]
blukey_tissue = blukey_tissue[!blukey_tissue$category %in% tissues_to_drop,]

#create color palettes and factor variables for plotting
B_palette = brewer.pal(9, "Reds")[3:6]
T_palette = brewer.pal(9, "Blues")[3:5]
NK_palette = brewer.pal(9, "Purples")[4] #brewer.pal requires at least n=3, so taking middle color
innate_palette = brewer.pal(9, "YlOrRd")[1:5]
blu_cell_types_ordered = c("Naive B cell", "Memory B cell", "Germinal B cell", "Plasma cell",
  "CD8 positive T cell", "CD4 positive T cell", "Memory T cell",
  "Natural Killer cell",
  "Dendritic cell", "Monocyte","Macrophage", "Eosinophil", "Neutrophil")
blukey_tissue$COLOR= factor(blukey_tissue$category,
  levels=blu_cell_types_ordered,
  labels=c(B_palette, T_palette, NK_palette, innate_palette))
blukey_tissue$category = factor(blukey_tissue$category, levels=blu_cell_types_ordered)

saveRDS(blu_dat, file="Blueprint_chromHMM_tracks_dataframe.rds")
saveRDS(blukey_tissue, file="Blueprint_chromHMM_tracks_tissue_key.rds")
saveRDS(blukey_state, file="Blueprint_chromHMM_tracks_state_key.rds")
