# Inflation comparison
Generate family-only meta-analysis summary statistics
  ```bash
  cd ${meta}
  metal ${scripts}/inflation_analysis/metal_script_for_fam_only_metaanaylsis
  mv METAANALYSIS1.TBL METAANALYSIS_fam_only_vcf_newqc_1.TBL
  mv METAANALYSIS1.TBL.info METAANALYSIS_fam_only_vcf_newqc_1.info
  ```

Create five sub-sampled case-control data sets
  ```bash
  for i in {1..5}; do
    mkdir -p ${inflation}/subdata/${i}
    mkdir -p ${inflation}/subanalyses/${i}
  done
  Rscript ${scripts}/inflation_analysis/subset_cohorts.R
  bash ${scripts}/inflation_analysis/subset_cohorts.sh > subset_cohorts.log
  ```

Run SNPTEST on sub-sampled data
  ```bash
  sbatch --array=1-5 ${scripts}/inflation_analysis/run_association_on_subsets.slurm AFR
  sbatch --array=1-5 ${scripts}/inflation_analysis/run_association_on_subsets.slurm AMR
  sbatch --array=1-5 ${scripts}/inflation_analysis/run_association_on_subsets.slurm EUR
  sbatch --array=1-5 ${scripts}/inflation_analysis/run_association_on_subsets.slurm FIN
  ```

Compare sub-sampled association statistic distributions to family-based association stats with equivalent power
  ```bash
  Rscript ${scripts}/inflation_analysis/visualize_inflation.R
  ```

Calculate number of informative trios across MAF spectrum
  ```R
  getNum = function(p, n) { 2*p*(1-p)*n }
  nums = NULL
  mafs = seq(0,0.5, by=0.001)
  for (i in 1:length(mafs)) { nums[i]=getNum(mafs[i],4578)}
  d = data.frame(mafs, nums)
  plot(d$mafs, d$nums, xlab="MAF", ylab="Number of informative trios")
  ```

