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

Create colocalisation table
	```bash
	Rscript ${scripts}/fine_mapping/merge_coloc_ecaviar_results.R
	```

