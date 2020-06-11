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
