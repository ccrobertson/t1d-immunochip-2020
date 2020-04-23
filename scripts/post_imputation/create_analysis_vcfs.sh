#!/bin/bash
set -euxo pipefail

group=$1
suffix=$2
dir=${filtered}/${group}
dir_all=${dir}/all
dir_rel=${dir}/relateds
dir_unrel=${dir}/unrelateds
dir_log=${dir}/log

mkdir -p ${dir_rel} ${dir_unrel} ${dir_log}

getDataSets () {

	local chr=$1
	echo -e "STARTING CHR${chr}"
	date

	#extract families
	echo -e "Creating trio data set..."
	families=${fam_assoc}/family_analysis_keep_list_iids.txt
	bcftools view --samples-file ${families} --force-samples --output-type z --output-file ${dir_rel}/chr${chr}.${suffix}.relateds.dose.vcf.gz ${dir_all}/chr${chr}.${suffix}.dose.vcf.gz

	#extract unrelateds
	echo -e "Creating cc data set..."
	unrelateds=${pca}/${project}_pca_${release_build}_pruned_notrios_${group}unrelated.txt
	bcftools view --samples-file <(awk '{print $2}' ${unrelateds}) --force-samples --output-type z --output-file ${dir_unrel}/chr${chr}.${suffix}.unrelateds.dose.vcf.gz ${dir_all}/chr${chr}.${suffix}.dose.vcf.gz

	#index files
	echo -e "Indexing data sets..."
	tabix -p vcf ${dir_rel}/chr${chr}.${suffix}.relateds.dose.vcf.gz
	tabix -p vcf ${dir_unrel}/chr${chr}.${suffix}.unrelateds.dose.vcf.gz

	echo -e "FINISHED CHR${chr}"
	date
}

#getDataSets "22" > ${dir_log}/getDataSets_chr22.log 2>&1 &

for ((i=1; i<=22; i+=1)); do
	echo -e "chr${i}"
	getDataSets "${i}"  > ${dir_log}/getDataSets.${suffix}.chr${i}.log 2>&1 &
done
wait

echo "DONE"
