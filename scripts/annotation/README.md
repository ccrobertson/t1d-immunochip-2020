	# Variant annotation
	Convert to ANNOVAR format
		```bash
		Rscript ${scripts}/fine_mapping/prep_for_annovar.R
		```

	Annotate with ANNOVAR
		```bash
		ANNOVAR=${resources}/annovar

		perl ${ANNOVAR}/table_annovar.pl ${annot}/annovar_input_file.txt ${ANNOVAR}/humandb/ -buildver hg38 -out ${annot}/annovar_output -protocol refGene,ensGene,cytoBand,1000g2015aug_eur,1000g2015aug_afr,avsnp150,dbnsfp30a -operation g,g,r,f,f,f,f -nastring . -csvout

		perl ${ANNOVAR}/table_annovar.pl ${annot}/annovar_input_file_credible_sets_with_proxies.txt ${ANNOVAR}/humandb/ -buildver hg38 -out ${annot}/annovar_output_credible_sets_with_proxies -protocol refGene,ensGene,cytoBand,1000g2015aug_afr,1000g2015aug_eur,avsnp150,dbnsfp30a -operation g,g,r,f,f,f,f -nastring . -csvout
		```

	Annotations for all the variants included in our analysis are in ${annot}/annovar_output.hg38_multianno.csv

	Annotations for credible sets are in ${annot}/annovar_output_credible_sets_with_proxies.hg38_multianno.csv

	Create annotated fine-mapping table (for manuscript)
		```bash
		Rscript ${scripts}/add_annotation_to_credible_set_table.R
		```

