# cis-QTL colocalisation

Run coloc
	```bash
	mkdir -p ${freeze}/locuscompare
	Rscript ${scripts}/fine_mapping/colocalisation.R
	```

Calculate allelic imbalance at colocalized caQTLs
	```bash
	bash ${scripts}/fine_mapping/allelic_imbalance.R

	## Trick pipeline into running on rs705705 (in perfect LD with rs705704)
	cp rs705704/rs705704_HET.txt rs705705/rs705705_HET.txt
	cp rs705704/rs705704_HOMREF.txt rs705705/rs705705_HOMREF.txt
	cp rs705704/rs705704_HOMALT.txt rs705705/rs705705_HOMALT.txt
	```

How many credible variants in peaks were genotyped on ImmunoChip?
```bash
awk '$2=="rs7731626"' ${data}/mega_release4_filtered_for_hwe.bim
awk '$2=="rs72928038"' ${data}/mega_release4_filtered_for_hwe.bim
awk '$2=="rs9388486"' ${data}/mega_release4_filtered_for_hwe.bim
awk '$2=="rs705704" || $2=="rs705705"' ${data}/mega_release4_filtered_for_hwe.bim
awk '$2=="rs11628807" || $2=="rs4383076" || $2=="rs11628876" || $2=="rs11160429"' ${data}/mega_release4_filtered_for_hwe.bim
```
