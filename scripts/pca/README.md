# Ancestry analysis

### Define ancestry groups for stratification
Calculate principal components on full cohort projected on to 1000G Phase 3.
```bash
bash ${scripts}/pca/pca_on_1000g.sh > ${logs}/pca_analysis.log
```

Now use k-means clustering of the top 10 projected PCs to determine gentically similar subgroups for association analysis. This type of kmeans approach has been suggested [previously](https://bmcgenet.biomedcentral.com/articles/10.1186/1471-2156-11-108))
```bash
Rscript ${scripts}/pca/define_ancestry_clusters.R ${pca} ${phenofile} ${resources}/1000genomes_populations_superpopulations.txt
Rscript ${scripts}/pca/plot_1000g_projected_pcs.R ${pca} ${phenofile} ${resources}/1000genomes_populations_superpopulations.txt
```

### Define trios and case-control analysis groups
```bash
bash ${scripts}/pre_imputation/define_trios.sh >${logs}/define_analysis_groups.log
```

### Calculate subgroup-specific PCs
Calculate subgroup-specific principal components to include as covariates in logistic regression models for case-control analysis (which will be stratified by ancestry group)
```bash
sbatch --output=${pca}/pca_case_control_AFR.log ${scripts}/pca/pca_case_control.slurm AFR
sbatch --output=${pca}/pca_case_control_AMR.log ${scripts}/pca/pca_case_control.slurm AMR
sbatch --output=${pca}/pca_case_control_FIN.log ${scripts}/pca/pca_case_control.slurm FIN
sbatch --output=${pca}/pca_case_control_EUR.log ${scripts}/pca/pca_case_control.slurm EUR
```

Plot PC1 vs PC2 colored by case/control status (in unrelateds)
```bash
Rscript ${scripts}/pca/plot_case_control_projected_pcs.R ${pca} AFR ${phenofile} mega_pca_b37
Rscript ${scripts}/pca/plot_case_control_projected_pcs.R ${pca} AMR ${phenofile} mega_pca_b37
Rscript ${scripts}/pca/plot_case_control_projected_pcs.R ${pca} EUR ${phenofile} mega_pca_b37
Rscript ${scripts}/pca/plot_case_control_projected_pcs.R ${pca} FIN ${phenofile} mega_pca_b37

Rscript ${scripts}/pca/plot_case_control_combine_plots.R
```
