# Generate BACH2 fine-mapping plots

Get annotation tracks from BLUEPRINT and ENCODE
```bash
cd ${PUBLIC_DATA}/EncodeReadmap
wget http://dcc.blueprint-epigenome.eu/#/md/secondary_analysis/Segmentation_of_ChIP-Seq_data_20140811

cd ${PUBLIC_DATA}/Blueprint
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/
```

Plot BACH2 region (Figure 3)
```bash
Rscript ${scripts}/fine_mapping/bach2_plots.R
```


Visualize credible set overlap with tracks in each region (for Supplemental Data)
```bash
Rscript ${scripts}/fine_mapping/generate_region_tracks.R
```


View chromHmm colors
```R
key = readRDS("Blueprint_chromHMM_tracks_state_key.rds")
key = key[!duplicated(key$label),]
key$plotState = factor(key$label, levels=key$label, labels=key$label)
key$x = 1
key$y = 1
pdf("chromhmm_blueprint_colors.pdf")
ggplot(data=key) + geom_rect(aes(xmin=x, xmax=x+1, ymin=y, ymax=y+1, fill=plotState)) + scale_fill_manual(breaks=key$label, values=key$color_hex) + facet_wrap(vars(label))
ggplot(data=key) + geom_rect(aes(xmin=x, xmax=x+1, ymin=y, ymax=y+1, fill=plotState)) + scale_fill_manual(breaks=key$label, values=key$color_hex) + facet_wrap(vars(label_12state))
dev.off()
```
