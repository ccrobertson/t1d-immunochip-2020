# Fine mapping scripts and running order


## GUESSFM using European cases and controls only
```
1) ./guessfm/guessfm_setup_5pc_diff.R 
```
Generate region VCFs only for hits in iChip regions, as agreed with Cassie and runs the GUESS stochastic search in the region


```
2) ./guessfm/guessfm_readin_script_gen.R 
```
Generates and submits batch scripts to run guessfm_readin_5pc_diff.R at each locus. This reads in the GUESS stochastic search and does the post-processing analysis to generate 'credible variants'.


```
3) ./guessfm/ld_creds.R 

```
Identify all variants in LD (r2>0.9) with credible variants that may have been removed due to failing imputation QC or other reasons.
NOTE: you will have to run ./atac/ld_find.R first for this script to work as it relies on haplotypes identified in that script.


```
4) ./guessfm/ldget_for_cassie.R
```
This script creates a file for Cassie to run through annovar to identify protein coding variants.


## Haplotype analysis accompanying the GUESSFM analyses using European cases and controls only
```
5) ./haplotypes/haplo_init.R
```
Haplotype analysis of the GUESSFM regions.


```
6) ./haplotypes/haplo_init_tags.R
```
Haplotype analysis of the GUESSFM regions including the tag variants in the haplotype plots whenever the index variant has PP<0.5.


## PAINTOR using all ancestry groups

```
7) ./paintor/meta_by_ancestry_5pc_diff_dens.R
```
Meta-analysis done for each ancestry group to generate ancestry sprecific summary stats for PAINTOR trans ethnic analysis.


```
8) ./paintor/paintor_setup_5pc_diff_dens.R
```
Creates the input files for paintor and the bash scripts to run the regions associatied on the iChip.
Requires (LDSTORE www.christianbenner.com)


```
9) ./paintor/paintor_readin_5pc_diff_up.R
```
Reads in the paintor results and produces plots.


