#!/bin/bash
set -euxo pipefail

#Usage: merge_EUR_batches.sh

dir=${results}
outdir=${dir}/EUR

findAndMerge () {

	local chr=$1
	echo -e "STARTING CHR${chr}"
	date

	outfiles=()
	batches=("EUR1" "EUR2" "EUR3")
	for batch in ${batches[@]}; do

		#Find file for batch
		if [ -f ${dir}/${batch}/chr${chr}.dose.vcf.gz ]; then
			filedir=${dir}/${batch}
		elif [ -f ${dir}/${batch}B/chr${chr}.dose.vcf.gz ]; then
			filedir=${dir}/${batch}B
		else
			echo -e "Missing ${batch} chr${chr}.dose.vcf.gz"
			exit
		fi
		infile=${filedir}/chr${chr}.dose.vcf.gz
		outfile=${filedir}/chr${chr}_shared_only.dose.vcf.gz

		#Filter files
		echo -e "Filtering ${infile} for shared variants..."
		zcat ${infile} | awk -v file=${outdir}/chr${chr}.merged_info 'BEGIN { while( getline<file ) {a[$1]=1} } $1~/^#/ || $3 in a' | bgzip > ${outfile}
		outfiles+=("${outfile}")

		echo "Indexing ${outfile} ..."
		tabix -p vcf ${outfile}


	done

	wait
	echo -e "Merging files"
	bcftools merge ${outfiles[@]} --output-type z --output ${outdir}/chr${chr}.dose.vcf.gz

	echo -e "FINISHED CHR${chr}"
	date
	echo -e " "
}


#for ((i=1; i<=22; i+=1)); do
#	echo -e "chr${i}"
#	findAndMerge "${i}"  > ${outdir}/findAndMerge_chr${i}.log 2>&1
#done
#wait

findAndMerge "4"  > ${outdir}/findAndMerge_chr4_3.log 2>&1
findAndMerge "5"  > ${outdir}/findAndMerge_chr5_3.log 2>&1
findAndMerge "9"  > ${outdir}/findAndMerge_chr9_3.log 2>&1
findAndMerge "12"  > ${outdir}/findAndMerge_chr12_3.log 2>&1
wait

echo "DONE"
