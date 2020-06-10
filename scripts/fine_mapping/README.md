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
Also requires GUESS to be installed from http://www.bgx.org.uk/software/guess.html

```
3) ./guessfm/guessfm_sum_5pc_diff.R
```
Produces various summaries or credible variants prior to including variants in high LD with the credible variants that failed imputation QC.


```
4) ./guessfm/ld_creds.R 

```
Identify all variants in LD (r2>0.9) with credible variants that may have been removed due to failing imputation QC or other reasons.
NOTE: you will have to run ./atac/ld_find.R first for this script to work as it relies on haplotypes identified in that script.


```
5) ./guessfm/ldget_for_cassie.R
```
This script creates a file for Cassie to run through annovar to identify protein coding variants. 
Also produces the summary of variants in all GUESSFM groups used in Figure 1 of manuscript.



## Haplotype analysis accompanying the GUESSFM analyses using European cases and controls only
```
6) ./haplotypes/haplo_init.R
```
Haplotype analysis of the GUESSFM regions.


```
7) ./haplotypes/haplo_init_tags.R
```
Haplotype analysis of the GUESSFM regions including the tag variants in the haplotype plots whenever the index variant has PP<0.5. The UBASH3A figure is also used in Figure 1 of the manuscript.

```
8) ./guessfms/ubash3a_look_up.R
```
Further examintion of the UBASH3A region - produces plots used in Figure 1 of the manuscript.



## PAINTOR using all ancestry groups

```
9) ./paintor/meta_by_ancestry_5pc_diff_dens.R
```
Meta-analysis done for each ancestry group to generate ancestry sprecific summary stats for PAINTOR trans ethnic analysis.


```
10) ./paintor/paintor_setup_5pc_diff_dens.R
```
Creates the input files for paintor and the bash scripts to run the regions associatied on the iChip.
Requires (LDSTORE www.christianbenner.com)


```
11) ./paintor/paintor_readin_5pc_diff_up.R
```
Reads in the paintor results and produces plots. Produces supplementary tables and plots for manuscript, as well as Figure 2 in manuscript.


