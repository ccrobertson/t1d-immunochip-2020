Ancestry analysis
  ```bash
  bash ${scripts}/pca/pca_on_1000g.sh > ${logs}/pca_analysis.log
  Rscript ${scripts}/pca/define_ancestry_clusters.R ${pca} ${phenofile} ${resources}/1000genomes_populations_superpopulations.txt
  Rscript ${scripts}/pca/plot_1000g_projected_pcas.R ${pca} ${phenofile} ${resources}/1000genomes_populations_superpopulations.txt
	```

Define trios and case-control analysis groups ([kmeans approach](https://bmcgenet.biomedcentral.com/articles/10.1186/1471-2156-11-108))
  ```bash
  bash ${scripts}/pre_imputation/define_trios.sh >${logs}/define_analysis_groups.log
	```

Case-control ancestry plots
```bash
sbatch --output=${pca}/pca_case_control_AFR.log ${scripts}/pca/pca_case_control.slurm AFR
sbatch --output=${pca}/pca_case_control_AMR.log ${scripts}/pca/pca_case_control.slurm AMR
sbatch --output=${pca}/pca_case_control_FIN.log ${scripts}/pca/pca_case_control.slurm FIN
sbatch --output=${pca}/pca_case_control_EUR.log ${scripts}/pca/pca_case_control.slurm EUR

#plot PC1 vs PC2 colored by case/control status (in unrelateds)
Rscript ${scripts}/pca/plot_case_control_projected_pcs.R ${pca} AFR ${phenofile} mega_pca_b37
Rscript ${scripts}/pca/plot_case_control_projected_pcs.R ${pca} AMR ${phenofile} mega_pca_b37
Rscript ${scripts}/pca/plot_case_control_projected_pcs.R ${pca} EUR ${phenofile} mega_pca_b37
Rscript ${scripts}/pca/plot_case_control_projected_pcs.R ${pca} FIN ${phenofile} mega_pca_b37

Rscript ${scripts}/pca/plot_case_control_combine_plots.R

#get case controls
awk 'BEGIN {print "fid", "iid"}' > ${cc_assoc}/unrelateds_all.txt
cat ${pca}/${nicknamepca}_pruned_notrios_???unrelated.txt >> ${cc_assoc}/unrelateds_all.txt
```

Prep data for imputation
  ```bash
  bash ${scripts}/pre_imputation/prep_data_for_imputation.sh > ${logs}/prep_data_for_imputation.log 2>&1
  ```

TOPMED IMPUTATION PREP STATISTICS:
* Start with 162,289 variants
* Remove 1,280 on X, Y or MT chromosomes
* Remove 9,437 because they have no position match (there are no variants in HRC with same position)
* Remove 8,289 because the allele frequency different between our EUR unrelated controls and HRC is >0.2
* Remove 2,169 because they are palindromic (A/T or C/G) with MAF>0.4
* Remove 407 because the alleles annotated on the chip do not match the alleles annotated in the HRC (e.g. ichip says its an A/T snp but HRC says its A/G -- flipping strands or ref/alt alleles cannot reconcile these as being the same variant)
* Remove 1 variant because it is a duplicate (same chr:pos:ref:alt variant appears twice)
* Finish with 140706 variants

1000G IMPUTATION PREP STATISTICS:
* Start with 162,289 variants
* Remove 1,280 on X, Y or MT chromosomes
* Remove 6,620 because they have no position match (there are no variants in 1000G with same position)
* Remove 8,326 because the allele frequency is different between our EUR unrelated controls and 1000G EUR population is >0.2
* Remove 2,191  because they are palindromic (A/T or C/G) with MAF>0.4
* Remove 1,607 because the alleles annotated on the chip do not match the alleles annotated in the 1000G (e.g. ichip says its an A/T snp but 1000G says its A/G -- flipping strands or ref/alt alleles cannot reconcile these as being the same variant)
* Remove 1 variant because it is a duplicate (same chr:pos:ref:alt variant appears twice)
* Finish with 142,264 variants

