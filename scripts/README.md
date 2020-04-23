<!-- TOC depthFrom:1 depthTo:1 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Background and contact](#background-and-contact)
- [Quality control and ancestry analyses](#quality-control-and-ancestry-analyses)
- [Imputation](#imputation)
- [Association analysis](#association-analysis)
- [Fine-mapping and enrichment analysis](#fine-mapping-and-enrichment-analysis)
- [QTL analyses](#qtl-analyses)
- [Other analyses](#other-analyses)

<!-- /TOC -->

# Background and contact
This is the analysis pipeline for "Fine-mapping, trans-ancestral and genomic analyses identify causal variants, cells, genes and drug targets for type 1 diabetes."

For questions, please contact the pipeline authors: Jamie Inshaw & Cassie Robertson (ccr5ju@virginia.edu)


# Quality control and ancestry analyses
* Map ImmunoChip probe sequences to hg19 (see probe_alignment/README_map_probes)
* Run the pre-processing QC pipeline (see pre_qc/Makefile)
* Genotype-based ancestry analyses (see pca/README.md)


# Imputation
* Prepare data for imputaton (see pre_imputation/README.md)
* Download, unzip, and filter imputation results (see post_imputation/README.md)
* Assess accuracy of imputation  (see imputation_accuracy_coverage/README.md)


# Association analysis
* Family-based association analysis (see assoc/README.md)
* Case-control association analysis (see assoc_case_control/README.md)
* Meta-analysis (assoc_case_control/README.md)
* Inflation comparison between case-control and family-based analysis	(see inflation_analysis/README.md)


# Fine-mapping and enrichment analysis
* Fine-mapping of T1D regions with GUESSFM (using EUR only) and PAINTOR (using all ancestry groups) (see fine_mapping/README.md)
* Haplotype analysis (Jamie please add info)
* Enrichment for credible variants in ATAC-seq peaks (Jamie please add info)


# QTL analyses
* cis-caQTL scan (see caqtl/README.md)
* cis-QTL colocalisation with T1D (see fine-mapping/README_coloc.md)
* BACH2 fine-mapping plots (see fine-mapping/README_bach2.md)


# Other analyses
* Define ichip and dense regions (Jamie please add info)
* Variant annotation (see annotation/README.md)
* Cohort summary tables (see manuscript/summary_tables.R)







scripts/fine_mapping/README.md
scripts/imputation_accuracy_coverage/README.md
scripts/inflation_analysis/README.md
scripts/pca/README.md
scripts/post_imputation/README.md
scripts/pre_imputation/README.md
