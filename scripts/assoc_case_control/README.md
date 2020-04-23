# Case-control meta-analysis scripts and running order


## Phase I analysis - case-control analyses only (no families)
```
1) define_ichip_regions.R
```
read in ichip regions hg19 and liftover all SNPs... define the regions in hg38:


```
2) ichip_reg_look.R
```
Some exploratory analyses examining where the genotyped variants are on the iChip, but	the main purpose is to save the ImmunoChip or densely genotyped regions on the ImmunoChip, only variants within these regions are examined in the genetic discovery analysis.


```
3) pre_meta_qc.R
```
This script generates shell scripts to perform association tests for each different cohort, directly from the VCF file, and taking into account the uncertainty in the genotype imputation.


```
4) post_imputation_qc_5pc_diff_dens.R
```
This script performs post-imputation QC on the set of analysed SNPs then generates METAL scripts to run the fixed effects meta analysis.


```
5) readin_metal_vcf_5pc_diff_dens.R 
```
Reads in results from initial meta run.


```
6) ./conditional/cond_setup_5pc_diff_dens.R 
```
Set up the conditional analyses by producing subsets of VCFs to perform conditional regressions with.


```
7) ./conditional/conditional_do_5pc_diff_dens.R
```
Generates the SNPTEST scripts which will carry out the stepwise conditional regressions.


## Phase I and II meta-analyis - case-control and families
```
8) meta_fams_inds_vcf_5pc_diff_dens.R  
```
Inverse variance-weighted meta analysis combining the family TDT results with the indepedent individual results.


## Dominant and recessive mode of inheritance analyses
```
9) ./domres/pre_meta_qc_dom.R 
```
This script generates shell scripts to perform association tests for each different cohort, directly from the VCF file and taking into account the uncertainty in the genotype imputation, assuming a dominant mode of inheritance.


```
10) ./domres/post_imputation_qc_newml_dom_5pc_diff_dens.R
```
This script performs post-imputation QC on the set of analysed SNPs assuming a dominant mode of inheritance.


```
11) ./domres/readin_metal_vcf_newml_dom_5pc_diff_dens.R
```
Reads in results from initial meta-analysis run with dominant mode of inheritance.


```
12) ./domres/pre_meta_qc_res.R
```
This script generates shell scripts to perform association tests for each different cohort, directly from the VCF file and taking into account the uncertainty in the genotype imputation, assuming a recessive mode of inheritance.


```
13) ./domres/post_imputation_qc_newml_res_5pc_diff_dens.R
```
This script performs post-imputation QC on the set of analysed SNPs assuming a recessive mode of inheritance.



```
14) ./domres/readin_metal_vcf_newml_res_5pc_diff_dens.R
```
Reads in results from initial meta-analysis run with recessive mode of inheritance.


# Generating some supplementary Tables
```
15) supplementary_tables_generate_dens_new.R
```
Generating a number of supplementary Tables for the manuscript.
