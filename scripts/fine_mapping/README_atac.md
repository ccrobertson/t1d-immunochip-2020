# ATAC-seq enrichment analyses using credible sets from GUESSFM analysis


```
1) ./atac/ld_find.R
```
Takes all the 1000 Genomes European population data (publically available) and identifies all 'haplotypes' (variants in r2>0.9 with each other) across the whole genome.
These are used in the enrichment analysis later.
NOTE: this is also used in the ./guessfm/ analysis to find variants dropped due to failing QC.

```
2) ./atac/ld_overall_newb_genedens_isletup_incld.R
```
Performs enrichment analysis: are more T1D credible variants in open chromatin in each cell type relative to what would be expected by chance.


```
3) ./atac/expr_stims_newb_genedens_isletup_incld.R

```
Performs enrichment analysis: are more T1D credible variants in ATAC-seq peaks in differentially-open peaks after stimulation or differentially-open peaks under resting conditions, than expected by chance.
It first identifies differentially-open peaks under resting and stimulated conditions, then carries out the enrichment analysis.
