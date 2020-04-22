#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --account rich_immunochip_impute
#SBATCH --partition=standard
#SBATCH --mem=50000

nicknamepca=${project}_pca_${release_build}
nickname=${project}_${release_build}

group=$1

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	RUN PCA FOR EACH CLUSTER
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comment "Run PCA for each cluster"
date

#extract trios from pruned data set
plink1.9 --bfile ${pca}/${nicknamepca}_pruned --keep ${fam_assoc}/family_analysis_keep_list.txt --make-bed --out ${pca}/${nicknamepca}_pruned_TRIOS
plink1.9 --bfile ${pca}/${nicknamepca}_pruned --remove ${fam_assoc}/family_analysis_keep_list.txt --make-bed --out ${pca}/${nicknamepca}_pruned_notrios


runPCA_with_projection () {
	local group=$1

	comment "Starting ${group}"

	#extract subjects in cluster
	plink1.9 --bfile ${pca}/${nicknamepca}_pruned_notrios --keep ${proc_for_imp}/${nickname}_${group}_cluster.txt --make-bed --out ${pca}/${nicknamepca}_pruned_notrios_${group}
	plink1.9 --bfile ${pca}/${nicknamepca}_pruned --keep ${proc_for_imp}/${nickname}_${group}_cluster.txt --make-bed --out ${pca}/${nicknamepca}_pruned_all_${group}

	#extract unrelated set
	king -b ${pca}/${nicknamepca}_pruned_notrios_${group}.bed --unrelated --prefix ${pca}/${nicknamepca}_pruned_notrios_${group}
	plink1.9 --bfile ${pca}/${nicknamepca}_pruned_notrios_${group} --keep ${pca}/${nicknamepca}_pruned_notrios_${group}unrelated.txt --make-bed --out ${pca}/${nicknamepca}_pruned_unrelated_${group}

	#get controls
	awk '$6==1 {print $1, $2}' ${pca}/${nicknamepca}_pruned_unrelated_${group}.fam > ${pca}/${nicknamepca}_pruned_unrelated_${group}_controls.txt
	plink1.9 --bfile ${pca}/${nicknamepca}_pruned_unrelated_${group} --keep ${pca}/${nicknamepca}_pruned_unrelated_${group}_controls.txt --make-bed --out ${pca}/${nicknamepca}_pruned_unrelated_${group}_controls

	#calculate principal components in controls
	plink2 --bfile ${pca}/${nicknamepca}_pruned_unrelated_${group}_controls --maf 0.05 --freq --pca var-wts --out ${pca}/${nicknamepca}_pruned_unrelated_${group}_controls --threads 32

	#project unrelated subjects onto control PCA space
	plink2 --bfile ${pca}/${nicknamepca}_pruned_unrelated_${group} --read-freq ${pca}/${nicknamepca}_pruned_unrelated_${group}_controls.afreq --score ${pca}/${nicknamepca}_pruned_unrelated_${group}_controls.eigenvec.var 2 3 \
	header-read no-mean-imputation variance-standardize --score-col-nums 5-14 --out ${pca}/${nicknamepca}_pruned_unrelated_${group}_pca_proj_onto_controls --threads 32

	#project all subjects onto control PCA space
	plink2 --bfile ${pca}/${nicknamepca}_pruned_all_${group} --read-freq ${pca}/${nicknamepca}_pruned_unrelated_${group}_controls.afreq --score ${pca}/${nicknamepca}_pruned_unrelated_${group}_controls.eigenvec.var 2 3 \
	header-read no-mean-imputation variance-standardize --score-col-nums 5-14 --out ${pca}/${nicknamepca}_pruned_all_${group}_pca_proj_onto_controls --threads 32

	echo -e "Finished ${group}"
	date
	echo -e " "
	echo -e " "
}

runPCA_with_projection ${group}
