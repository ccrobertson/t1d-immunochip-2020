#!/bin/bash

################################################################################
#
#	 Usage: bash pca_analysis.sh {bfile} {phenofile} {build}
#
################################################################################

set -euxo pipefail

#input variables
raw=${genodat}
build=${release_build}
nickname=${project}_pca_${release_build}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	PREP DATA SETS FOR PCA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
date


#### Remove high LD regions (b37 coordinates obtained from https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD))
comment "Removing high LD regions"
if [ "${build}" != 'b37' ]; then
	#liftover regions to correct build
	bash ${scripts}/pca/liftover_high_ld_regions.sh
fi
#add label for each region to ld region file (required format for plink "--exclude range" option)
awk '{print $1, $2, $3, "R"$NR}' ${resources}/high_ld_regions_${build}.txt > ${resources}/high_ld_regions_${build}.range.txt
plink1.9 --bfile ${raw} --exclude range ${resources}/high_ld_regions_${build}.range.txt --make-bed --out ${pca}/${nickname}_dropped_highLD



#### Prune variants --> pairs of variants within 50kb window must have r2<0.2
comment "Pruning variants for R2<0.2 within 50kb window"
plink1.9 --bfile ${pca}/${nickname}_dropped_highLD --geno 0.05 --maf 0.01 --indep-pairwise 50 5 0.2 --autosome --biallelic-only strict --make-bed --out ${pca}/${nickname}_dropped_highLD
plink1.9 --bfile ${pca}/${nickname}_dropped_highLD --extract ${pca}/${nickname}_dropped_highLD.prune.in --make-bed --out ${pca}/${nickname}_pruned



#### Prep 1000G Data
comment "Filtering 1000G data for pca analysis"
#extract pruned SNPs from 1000G
awk -F ' '  '{print $1, $4, $4}' ${pca}/${nickname}_pruned.bim | sed 's/ /\t/g' > ${pca}/${nickname}_pruned.tab
for i in {1..22}; do
	echo "Starting filtering procedure for Chromsome ${i}"

	if [ "${build}" == 'b37' ]; then
		raw_1000g=${resources}/1000g_b37/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
	elif [ "${build}" == 'b38' ]; then
		raw_1000g=${resources}/1000g_b38/ALL.chr${i}_GRCh38.genotypes.20170504.vcf.gz
	else
		echo -e "Please download appropriate reference panel for build ${build}"
		exit
	fi
	bcftools view --regions-file ${pca}/${nickname}_pruned.tab ${raw_1000g} -m2 -M2 -v snps  > ${pca}/1000_genomes_${build}_pruned_overlap_chr${i}.vcf
	echo "Completed Chromosome ${i}"
done

#convert 1000G to plink and merge chromosomes
comment "Converting 1000G vcfs to plink"
echo -n > ${pca}/1000_genomes_${build}_merge_list.txt
for i in {1..22}; do
	#bcftools norm -Ou -m -any ${pca}/1000_genomes_${build}_pruned_overlap_chr${i}.vcf | bcftools norm -Ou -f ${resources}/human_g1k_v37.fasta | bcftools annotate -Oz -x ID -I +'%CHROM:%POS:%REF:%ALT' > ${pca}/1000_genomes_${build}_pruned_overlap_chr${i}_stripped.vcf.gz
	bcftools norm -Ou -m -any ${pca}/1000_genomes_${build}_pruned_overlap_chr${i}.vcf | bcftools annotate -Oz -x ID -I +'%CHROM:%POS:%REF:%ALT' > ${pca}/1000_genomes_${build}_pruned_overlap_chr${i}_stripped.vcf.gz
	plink2 --vcf ${pca}/1000_genomes_${build}_pruned_overlap_chr${i}_stripped.vcf.gz --const-fid --make-bed --out ${pca}/1000_genomes_${build}_pruned_overlap_chr${i}
	echo -e ${pca}/1000_genomes_${build}_pruned_overlap_chr${i} >> ${pca}/1000_genomes_${build}_merge_list.txt
done
plink1.9 --merge-list ${pca}/1000_genomes_${build}_merge_list.txt --out ${pca}/1000_genomes_${build}_pruned_overlap_ALL_tmp

#filter for maf
plink1.9 --bfile ${pca}/1000_genomes_${build}_pruned_overlap_ALL_tmp --maf 0.05 --make-bed --out ${pca}/1000_genomes_${build}_pruned_overlap_ALL




#### Prep MEGA Data
#update variant ids in mega to CHR:POS:REF:ALT
awk '{print $2, $1":"$4":"$6":"$5}' ${pca}/${nickname}_pruned.bim > ${pca}/update_variants_ids.txt
plink1.9 --bfile ${pca}/${nickname}_pruned --update-name ${pca}/update_variants_ids.txt --make-bed --out ${pca}/${nickname}_pruned_updatedIDs

#extract intersection of variants from mega
awk '{print $2}' ${pca}/1000_genomes_${build}_pruned_overlap_ALL.bim > ${pca}/variant_keep_list1
plink1.9 --bfile ${pca}/${nickname}_pruned_updatedIDs --extract ${pca}/variant_keep_list1 --make-bed --out ${pca}/${nickname}_pruned_FINAL

#extract intersection of variants from 1000g
awk '{print $2}' ${pca}/${nickname}_pruned_FINAL.bim > ${pca}/variant_keep_list2
plink1.9 --bfile ${pca}/1000_genomes_${build}_pruned_overlap_ALL --extract ${pca}/variant_keep_list2 --make-bed --out ${pca}/1000_genomes_${build}_pruned_overlap_ALL_FINAL





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Create function for projecting pcas onto reference cohort
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#Usage: calculateProjectedPCAs {1000g file} {1000g nickname} {mega file} {mega nickname}
calculateProjectedPCAs () {

	#calculate principal components in reference cohort
	plink2 --bfile ${pca}/$1 --freq --pca var-wts --out ${pca}/pca_$2 --threads 32

	#project mega onto reference pcs
	plink2 --bfile ${pca}/$3 --read-freq ${pca}/pca_$2.afreq --score ${pca}/pca_$2.eigenvec.var 2 3 \
	header-read no-mean-imputation variance-standardize --score-col-nums 5-14 --out ${pca}/pca_proj_${4}_onto_$2 --threads 32

	#get 1000 genomes pcs on the same scale as the projected mega pcs
	plink2 --bfile ${pca}/$1 --read-freq ${pca}/pca_$2.afreq --score ${pca}/pca_$2.eigenvec.var 2 3 \
	header-read no-mean-imputation variance-standardize --score-col-nums 5-14 --out ${pca}/pca_proj_$2 --threads 32
}




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	FULL COHORT PCA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comment "PCA for entire cohort"
calculateProjectedPCAs 1000_genomes_${build}_pruned_overlap_ALL_FINAL 1000g ${nickname}_pruned_FINAL ${project}_all



date
