#!/bin/bash
set -euxo pipefail

dir=${results}
outdir=${dir}/EUR
mkdir -p ${outdir}

findAndMergeStats () {

	local chr=$1
	echo -e "STARTING CHR${chr}"
	date

	#Get file for each batch
	infofiles=()
	batches=("EUR1" "EUR2" "EUR3")
	for batch in ${batches[@]}; do
		if [ -f ${dir}/${batch}/chr${chr}.dose.vcf.gz ]; then
			filedir=${dir}/${batch}

		elif [ -f ${dir}/${batch}B/chr${chr}.dose.vcf.gz ]; then
			filedir=${dir}/${batch}B

		else
			echo -e "Missing ${batch} chr${chr}.dose.vcf.gz"
			exit
		fi
		infofiles+=("${filedir}/chr${chr}.info.gz")
	done
	echo "INFO FILES:"
	echo "${infofiles[@]}"
	echo " "

	#Get merged info file
	echo -e "Calculating merged rsq stats"
	Rscript ${scripts}/post_imputation/merge_EUR_info_stats.R ${chr} ${infofiles[@]} ${outdir}
	awk 'NR>1 {print $1}' ${outdir}/chr${chr}.merged_info > ${outdir}/chr${chr}_variants_to_keep.txt
	date

	echo -e "FINISHED CHR${chr}"
	date
	echo -e " "
}



for i in {1..22}; do
	findAndMergeStats "${i}"
done
wait

echo "DONE"
