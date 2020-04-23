# Generate BACH2 fine-mapping plots

Get annotation tracks from BLUEPRINT and ENCODE
	```bash
	cd ${PUBLIC_DATA}/EncodeReadmap
	wget http://dcc.blueprint-epigenome.eu/#/md/secondary_analysis/Segmentation_of_ChIP-Seq_data_20140811

	cd ${PUBLIC_DATA}/Blueprint
	wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/
	```

Visualize credible set overlap with tracks in each region (for Supplemental Data)
	```bash
	Rscript ${scripts}/fine_mapping/generate_region_tracks.R
	```

For each SNP, extract region by genotype and create mean track.
	```bash
	bash ${scripts}/fine_mapping/chromqtl_tracks_get_region_by_genotype.sh
	```

