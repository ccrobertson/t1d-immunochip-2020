Results from fine-mapping T1D-associated ImmunoChip regions using the software [GUESSFM](https://github.com/chr1swallace/GUESSFM) are in these directory subdirectories.

For easy access to the image files, we recommend cloning the repository to your local computer, which will allow you to open any image (.png) files with your preferred desktop image viewer.
You can do so by navigaating to the desired directory and running the following command in your terminal:
```
git clone https://github.com/ccrobertson/t1d-immunochip-2020
```

Each directory shows a fine mapped region, named either by a proximal candidate causal gene or a cytoband if no known candidate causal genes are in the region.
As standard, the directory will contain the following figures:
1)	A haplotype analysis. The name of this file is the chromosome, position, reference and alternate allele from the index variant in the univariable analysis. The first panel contains estimates of each haplotype on T1D risk, containing index variants from GUESSFM groups with group posterior probability >0.8. The second panel contains the frequency estimates of these haplotypes in the MEGAT1D population. The third panel shows the haplotype, with black squares representing the major allele and white squares representing the minor allele at that variant. 
2)	A cumulative probability plot, named “cum_prob_plot_upd.png”. This shows the proportion of times a variant from the group indexed by the column name were selected in the stochastic search. Black regions are model visits where a variant from this group were included in the model, white regions reflect models where no variant from that group were included in the model. If two or more columns are black simultaneously, this reflects model visits where more than one variant were included in the model.
3)	GUESSFM output, named “init_out_genes_upd.png “. The fourth panel shows the linkage disequilibrium between GUESSM prioritised variants. The third panel shows the GUESSFM prioritised variants, coloured by group (these colours match those in the haplotype analysis), the dots are the variant posterior probability and the shaded region shows the group posterior probability. The second panel maps the location of the variants in the third panel to their position in the genome. The first panel shows the location of any protein coding genes in the region (genome build 38).

In four directories, there are five figures, where there is evidence of a ‘tag’ variant, which has a strong univariable association statistic due to being in linkage disequilibrium with two or more disease associated variants, but does not appear to be disease associated when present in isolation without the disease-associated variants. These figures are called the same as the haplotype analysis name, but with”_inc_tag” at the end. The following directories contain five figure and therefore have evidence of a tag variant:

- CTLA4
- IL2/IL21
- MEG3
- UBASH3A

In three directories, there are three figures. These are regions where there are no groups with posterior probability >0.8 and therefore no haplotype analysis:

- RASGRP1
- IRF4
- RPAP2


