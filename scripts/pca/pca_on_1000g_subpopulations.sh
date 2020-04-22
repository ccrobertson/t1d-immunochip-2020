#!/bin/bash

################################################################################
#
#	 Usage: bash pca_analysis.sh {bfile} {phenofile} {build} {prefix}
#
################################################################################

set -euxo pipefail
source ${PROJECT_MEGA}/scripts/configure_env.sh


#input variables
raw=${genodat}
build="b37"
prefix=mega_pca
nickname=${prefix}_${build}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Create function for projecting pcas onto reference cohort
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#Usage: calculateProjectedPCAs {1000g file} {1000g nickname} {mega file} {mega nickname} 
calculateProjectedPCAs () {

	#calculate principal components in reference cohort
	${plink2} --bfile ${pca}/$1 --freq --pca var-wts --out ${pca}/pca_$2 --threads 32

	#project mega onto reference pcs
	${plink2} --bfile ${pca}/$3 --read-freq ${pca}/pca_$2.afreq --score ${pca}/pca_$2.eigenvec.var 2 3 \
	header-read no-mean-imputation variance-standardize --score-col-nums 5-14 --out ${pca}/pca_proj_${4}_onto_$2 --threads 32
		
	#get 1000 genomes pcs on the same scale as the projected mega pcs
	${plink2} --bfile ${pca}/$1 --read-freq ${pca}/pca_$2.afreq --score ${pca}/pca_$2.eigenvec.var 2 3 \
	header-read no-mean-imputation variance-standardize --score-col-nums 5-14 --out ${pca}/pca_proj_$2 --threads 32				
}


		
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	HISPANIC PCA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comment "Admixed-American/Hispanic specific PCA"


#find column number for self-reported hispanic 
awk 'BEGIN {FS=","} NR=1 {for (i = 1; i <= NF; i++) { if ($i == "Hispanic") {print i}} }' /m/jdrfdn_scratch/users/projects/IMCHIP/T1DGC.2011.03_Golden_Pedigree/T1DGC.2011.03_Golden_Pedigree/Resources/T1DGC.2011.03_Resources.csv 

#extract hispanic subjects 
awk 'BEGIN {FS=","} $19=="1" {print $6}' /m/jdrfdn_scratch/users/projects/IMCHIP/T1DGC.2011.03_Golden_Pedigree/T1DGC.2011.03_Golden_Pedigree/Resources/T1DGC.2011.03_Resources.csv > ${pca}/hispanic_subjects.txt 
awk -v file=${pca}/hispanic_subjects.txt 'BEGIN { while( getline<file) {a[$1]=1}} $3 in a {print $1$2,$3}' ${phenofile} > ${pca}/hispanics_keeplist.txt
${plink} --bfile ${pca}/${nickname}_pruned_FINAL --keep ${pca}/hispanics_keeplist.txt --make-bed --out ${pca}/${nickname}_pruned_HISPANICS

#get AMR 1000G subjects
awk '$3=="AMR" {print "0",$1}' ${resources}/1000genomes_populations_superpopulations.txt > ${pca}/1000g_AMR_ids.txt
${plink} --bfile ${pca}/1000_genomes_${build}_pruned_overlap_ALL_FINAL --keep ${pca}/1000g_AMR_ids.txt --maf 0.05 --make-bed --out ${pca}/1000_genomes_${build}_pruned_overlap_ALL_AMR


calculateProjectedPCAs 1000_genomes_${build}_pruned_overlap_ALL_AMR 1000g_amr ${nickname}_pruned_HISPANICS mega_hisp
calculateProjectedPCAs 1000_genomes_${build}_pruned_overlap_ALL_AMR 1000g_amr ${nickname}_pruned_ALL_FINAL mega_all




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	EUROPEAN PCA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comment "European-specific PCA"

	
#get EUR 1000G subjects
awk '$3=="EUR" {print "0",$1}' ${resources}/1000genomes_populations_superpopulations.txt > ${pca}/1000g_EUR_ids.txt
${plink} --bfile ${pca}/1000_genomes_${build}_pruned_overlap_ALL_FINAL --keep ${pca}/1000g_EUR_ids.txt --maf 0.05 --make-bed --out ${pca}/1000_genomes_${build}_pruned_overlap_ALL_EUR
 
#get EUR mega subjects
awk 'NR>1 && $4=="EUR" {print $1, $2}' ${pca}/mega_cluster_assignments.txt > ${pca}/mega_cluster_EUR.txt
${plink} --bfile ${pca}/${nickname}_pruned_FINAL --keep ${pca}/mega_cluster_EUR.txt --make-bed --out ${pca}/${nickname}_pruned_EUR


calculateProjectedPCAs 1000_genomes_${build}_pruned_overlap_ALL_EUR 1000g_eur ${nickname}_pruned_EUR mega_eur
 


