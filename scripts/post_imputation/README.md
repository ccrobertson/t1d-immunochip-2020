# Post imputation processing
  ```bash
  bash ${scripts}/post_imputation/post_imputation_master.sh {group} {password} {mafthresh} {rsqthresh}
  bash ${scripts}/post_imputation/post_imputation_master_EUR.sh {password1} {password2} {password3} {mafthresh} {rsqthresh}
  ```

Get combined info files
```bash
bash ${scripts}/post_imputation/combine_info_files.sh
```




# Fix header on 1000G imputation files
Note: the following files are missing a vcf header. This was the case when the were downloaded
from the imputation server and unzipped -- perhaps due to a bug on the server?
I checked to make sure the number of variants in these files matches the corresponding info file
and is approximately the same as the number of variants in the same chromosome for the other EUR groups
e.g. EUR2/chr4.dose.gz has 130271 variants and EUR1/chr4.dose.gz has 128735 variants
so it doesn't seem like the files were truncated during download
```bash
/nv/vol185/MEGA/release4/IMPUTED_1KG/results/EUR2/chr4.dose.vcf.gz
/nv/vol185/MEGA/release4/IMPUTED_1KG/results/EUR2/chr5.dose.vcf.gz
/nv/vol185/MEGA/release4/IMPUTED_1KG/results/EUR2/chr12.dose.vcf.gz
/nv/vol185/MEGA/release4/IMPUTED_1KG/results/EUR3/chr9.dose.vcf.gz
```

Create header for each batch
```bash
cd /nv/vol185/MEGA/release4/IMPUTED_1KG/results
zcat EUR1/chr1.dose.vcf.gz | head -n 100 | awk '$1~/^#/' | bgzip > EUR1/EUR1_header.gz
zcat EUR2/chr1.dose.vcf.gz | head -n 100 | awk '$1~/^#/' | bgzip > EUR2/EUR2_header.gz
zcat EUR3/chr1.dose.vcf.gz | head -n 100 | awk '$1~/^#/' | bgzip > EUR3/EUR3_header.gz
```

Add header to files
```bash
cd /nv/vol185/MEGA/release4/IMPUTED_1KG/results/EUR2
mv chr4.dose.vcf.gz chr4.dose_no_header.vcf.gz
mv chr5.dose.vcf.gz chr5.dose_no_header.vcf.gz
mv chr12.dose.vcf.gz chr12.dose_no_header.vcf.gz
cat EUR2_header.gz chr4.dose_no_header.vcf.gz > chr4.dose.vcf.gz
cat EUR2_header.gz chr5.dose_no_header.vcf.gz > chr5.dose.vcf.gz
cat EUR2_header.gz chr12.dose_no_header.vcf.gz > chr12.dose.vcf.gz

cd /nv/vol185/MEGA/release4/IMPUTED_1KG/results/EUR3
mv chr9.dose.vcf.gz chr9.dose_no_header.vcf.gz
cat EUR3_header.gz chr9.dose_no_header.vcf.gz > chr9.dose.vcf.gz
```

Rerun merging for chr 4, 5, 9, and 12
```bash
module load htslib
cd /nv/vol185/MEGA/release4/IMPUTED_1KG/results
bash ${scripts}/post_imputation/merge_EUR_batches.sh > merge_EUR_batches_3.log 2>&1

sbatch --output merge_chr4.log ${scripts}/post_imputation/merge_chr4.slurm
```
