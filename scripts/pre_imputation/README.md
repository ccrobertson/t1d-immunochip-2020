# Prep data for imputation

Get the set of all unrelated cases and controls
```bash
awk 'BEGIN {print "fid", "iid"}' > ${cc_assoc}/unrelateds_all.txt
cat ${pca}/${nicknamepca}_pruned_notrios_???unrelated.txt >> ${cc_assoc}/unrelateds_all.txt
```

Run imputation preparation pipeline
```bash
bash ${scripts}/pre_imputation/prep_data_for_imputation.sh > ${logs}/prep_data_for_imputation.log 2>&1
```

### TOPMED IMPUTATION PREP STATISTICS:
* Start with 162,289 variants
* Remove 1,280 on X, Y or MT chromosomes
* Remove 9,437 because they have no position match (there are no variants in HRC with same position)
* Remove 8,289 because the allele frequency different between our EUR unrelated controls and HRC is >0.2
* Remove 2,169 because they are palindromic (A/T or C/G) with MAF>0.4
* Remove 407 because the alleles annotated on the chip do not match the alleles annotated in the HRC (e.g. ichip says its an A/T snp but HRC says its A/G -- flipping strands or ref/alt alleles cannot reconcile these as being the same variant)
* Remove 1 variant because it is a duplicate (same chr:pos:ref:alt variant appears twice)
* Finish with 140706 variants

### 1000G IMPUTATION PREP STATISTICS:
* Start with 162,289 variants
* Remove 1,280 on X, Y or MT chromosomes
* Remove 6,620 because they have no position match (there are no variants in 1000G with same position)
* Remove 8,326 because the allele frequency is different between our EUR unrelated controls and 1000G EUR population is >0.2
* Remove 2,191  because they are palindromic (A/T or C/G) with MAF>0.4
* Remove 1,607 because the alleles annotated on the chip do not match the alleles annotated in the 1000G (e.g. ichip says its an A/T snp but 1000G says its A/G -- flipping strands or ref/alt alleles cannot reconcile these as being the same variant)
* Remove 1 variant because it is a duplicate (same chr:pos:ref:alt variant appears twice)
* Finish with 142,264 variants
