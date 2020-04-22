#!/bin/bash
##### Using Will Rayner's script (http://www.well.ox.ac.uk/~wrayner/tools/)

set -euxo pipefail

nickname=${project}_${release_build}

cd ${proc_for_imp}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comment	"Define ancestry groups"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
makeAncestryCluster () {
	cluster=$1
	awk -v cluster=${cluster} '$4==cluster {print $1, $2}' ${pca}/${project}_cluster_assignments.txt > ${nickname}_${cluster}_cluster.txt
}

makeAncestryCluster "EUR"
makeAncestryCluster "AFR"
makeAncestryCluster "FIN"
makeAncestryCluster "EAS"
makeAncestryCluster "AMR"

#Divide EUR into 20k batches
awk 'NR<=15000' ${nickname}_EUR_cluster.txt > ${nickname}_EUR1_cluster.txt
awk 'NR>15000 && NR<=30000' ${nickname}_EUR_cluster.txt > ${nickname}_EUR2_cluster.txt
awk 'NR>30000' ${nickname}_EUR_cluster.txt > ${nickname}_EUR3_cluster.txt


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comment	"Define European controls for alignment"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Define EUR subjects for alignment --> intersection of EUR cluster and UK control cohorts: BC58,cbr_controls,NIMH,UKBS
awk -v file=${nickname}_EUR_cluster.txt 'BEGIN { while (getline<file) {a[$1"\t"$2]=1} } ($1~/BC58|cbr_controls|NIMH|UKBS/) && ($1$2"\t"$3 in a) {print $1$2,$3}' ${phenofile} > ${nickname}_EUR_for_alignment.txt
plink --bfile ${genodat} --keep ${nickname}_EUR_for_alignment.txt --make-bed --out ${nickname}_EUR_for_alignment_tmp
#confirm that subjects are unrelated
#${king} -b ${nickname}_EUR_for_alignment.bed --unrelated --degree 2 --prefix ${nickname}_EUR_for_alignment



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comment	"Filter for deviation from HWE"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plink --bfile ${nickname}_EUR_for_alignment_tmp --hardy --out ${nickname}_EUR_for_alignment_tmp
awk '$3=="UNAFF" && $1<23' ${nickname}_EUR_for_alignment_tmp.hwe| awk '$9<0.00005 {print $2}' > snps_failing_hwe.txt
plink --bfile ${nickname}_EUR_for_alignment_tmp --exclude snps_failing_hwe.txt --make-bed --out ${nickname}_EUR_for_alignment
plink --bfile ${genodat} --exclude snps_failing_hwe.txt --make-bed --out ${genodat}_filtered_for_hwe

#generate frequency file
plink --bfile ${nickname}_EUR_for_alignment --freq --out ${nickname}_EUR_for_alignment_frq


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	Define functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

makeClusteredDatasets () {
	cluster=$1
	mkdir -p ${cluster}
	for chr in {1..23}; do
	plink --bfile ${nickname}-updated-chr${chr} --keep-allele-order --keep ../${nickname}_${cluster}_cluster.txt --make-bed --out ${cluster}/${nickname}-updated-chr${chr}_${cluster}_cluster
	done
}

convertToVCF () {
	cluster=$1
	for chr in {1..23}; do
		plink --bfile ${cluster}/${nickname}-updated-chr${chr}_${cluster}_cluster --keep-allele-order --recode vcf-iid --out ${cluster}/${nickname}-updated-chr${chr}_${cluster}_cluster
		bgzip ${cluster}/${nickname}-updated-chr${chr}_${cluster}_cluster.vcf

	done
}

prepForImputation () {

	refpanel=$1

	#
	comment "Running Rayner ${refpanel} alignment script"
	if [ ${refpanel} == "HRC" ]; then
		perl ${resources}/HRC-1000G-check-bim-NoReadKey.pl -b ${proc_for_imp}/${nickname}_EUR_for_alignment.bim -f ${proc_for_imp}/${nickname}_EUR_for_alignment_frq.frq -r ${resources}/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h
	elif [ ${refpanel} == "1000G" ]; then
		perl ${resources}/HRC-1000G-check-bim-NoReadKey.pl -b ${proc_for_imp}/${nickname}_EUR_for_alignment.bim -f ${proc_for_imp}/${nickname}_EUR_for_alignment_frq.frq -r ${resources}/1000GP_Phase3_combined.legend.gz -g -p EUR
	fi

	#
	comment "Run Rayner commands"
	sed "1s|${nickname}_EUR_for_alignment|${genodat}_filtered_for_hwe|" Run-plink.sh > Run-plink_B.sh
	sed "s|plink|plink|g;s|--bfile ${nickname}_EUR_for_alignment-updated|--bfile ${nickname}-updated|g;s|--out ${nickname}_EUR_for_alignment-updated|--out ${nickname}-updated|g" Run-plink_B.sh > Run-plink_C.sh

	bash Run-plink_C.sh

	#
	comment	"Break into clusters"
	makeClusteredDatasets "AFR"
	makeClusteredDatasets "FIN"
	makeClusteredDatasets "EAS"
	makeClusteredDatasets "AMR"
	makeClusteredDatasets "EUR1"
	makeClusteredDatasets "EUR2"
	makeClusteredDatasets "EUR3"

	#
	comment	"Convert to VCF"
	convertToVCF "AFR"
	convertToVCF "FIN"
	convertToVCF "EAS"
	convertToVCF "AMR"
	convertToVCF "EUR1"
	convertToVCF "EUR2"
	convertToVCF "EUR3"

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comment	"Align data to HRC"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### HRC alignment (for imputation to topmed)
mkdir -p ${proc_for_imp}/topmed_alignment
cd ${proc_for_imp}/topmed_alignment
prepForImputation "HRC"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
comment	"Align data to 1000G"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### 1000G alignment (for imputation to 1000G)
mkdir -p ${proc_for_imp}/kg_alignment
cd ${proc_for_imp}/kg_alignment
prepForImputation "1000G"
