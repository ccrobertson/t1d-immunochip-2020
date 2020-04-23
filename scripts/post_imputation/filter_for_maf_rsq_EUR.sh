#!/bin/bash
set -euxo pipefail

#Usage: filter_for_maf_rsq_EUR.sh {group} {maf_thresh} {rsq_thresh} {infile_suffix}
date

group=$1
maf_thresh=$2
rsq_thresh=$3
infile_suffix=$4
outfile_suffix=filter_maf_gt_${maf_thresh}_rsq_gt_${rsq_thresh}


INPUT_DIR=${results}/${group}
OUTPUT_DIR=${filtered}/${group}/all
mkdir -p ${filtered}/${group}
mkdir -p ${filtered}/${group}/all
cd ${OUTPUT_DIR}



echo -e "Input files: ${INPUT_DIR}/chr*.${infile_suffix}.vcf.gz"
echo -e "MAF threshold: ${maf_thresh}"
echo -e "Rsq threshold: ${rsq_thresh}"
echo -e "Output files: ${OUTPUT_DIR}/chr*.${outfile_suffix}.vcf.gz"
echo -e " "

processData () {

	local chr=$1
	echo -e "STARTING CHR${chr}"

	infile=${INPUT_DIR}/chr${chr}.${infile_suffix}.vcf.gz
	info=${INPUT_DIR}/chr${chr}.merged_info
	keepfile=${OUTPUT_DIR}/chr${chr}.${outfile_suffix}.merged_info
	outfile=${OUTPUT_DIR}/chr${chr}.${outfile_suffix}.dose.vcf.gz

	if [ -f ${infile} ]; then

		echo -e "Filtering ${infile} for MAF>${maf_thresh} and Rsq>${rsq_thresh}..."
		awk -v maf=${maf_thresh} -v rsq=${rsq_thresh} '$3>maf && $4>rsq' ${info} > ${keepfile}
		zcat ${infile} | awk -v file=${keepfile} 'BEGIN { while( getline<file ) {a[$1]=1} } $1~/^#/ || $3 in a' | bgzip > ${outfile}

		echo "Indexing ${outfile} ..."
		tabix -p vcf ${outfile}
		echo "Indexing complete"
		date

	fi

	echo -e "FINISHED CHR${chr}"
	date
	echo " "

}


for ((i=1; i<=22; i+=1)); do
	echo -e "chr${i}"
	processData "${i}" > ${OUTPUT_DIR}/chr${i}.${outfile_suffix}.log 2>&1 &
done
wait

#processData "22"

echo "DONE"
date
