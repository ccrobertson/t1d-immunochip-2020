# haplotype analyses using credible sets from GUESSFM analysis


```
1) ./haplotype/haplotype_init.R
```
Performs haplotype analysis using all index variants from GUESSFM groups with PP>0.5 (except IL2RA, which is particularly complex, so threshold of 0.8 chosen).


```
2) ./haplotype/haplo_init_tags.R
```
Re-performs the haplotype analysis, but this time if the index variant from the univaribale analysis has fine-mapping PP<0.5, it is included in the 
haplotype analysis to show that when that index variant appears on it's own without the other fine-mapping signals, there is no effect of the 
univariable analysis index variant on T1D risk.


