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


## Haplotype analysis accompanying the GUESSFM analyses using European cases and controls only
```
3) ./haplotypes/haplo_init.R
```
Haplotype analysis of the GUESSFM regions.


```
4) ./haplotypes/haplo_init_tags.R
```
Haplotype analysis of the GUESSFM regions including the tag variants in the haplotype plots whenever the index variant has PP<0.5.


## PAINTOR using all ancestry groups

```
5) ./paintor/meta_by_ancestry_5pc_diff_dens.R
```
Meta-analysis done for each ancestry group to generate ancestry sprecific summary stats for PAINTOR trans ethnic analysis.


```
6) ./paintor/paintor_setup_5pc_diff_dens.R
```
Creates the input files for paintor and the bash scripts to run the regions associatied on the iChip.
Requires (LDSTORE www.christianbenner.com)


```
7) ./paintor/paintor_readin_5pc_diff_up.R
```
Reads in the paintor results and produces plots.

## Other - Cassie, please complete?
