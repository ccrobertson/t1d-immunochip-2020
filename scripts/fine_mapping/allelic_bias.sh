#!/bin/bash


snp=$1
chr=$2
pos=$3
a1=$4
a2=$5

vardir=${freeze}/allelic_bias/${snp}

start=$( expr ${pos} - 1)
end=${pos}
echo -e "${chr}\t${start}\t${end}\tAF=0.500;AD::${a1}=10,${a2}=10" > ${vardir}/count_hetalleles_sites_${snp}.bed

calculateAndSummarize () {
	local group=$1
	echo -e "Sample\t${a1}\t${a2}" > ${vardir}/allelic_bias_${group}_summary.txt
	> ${vardir}/allelic_bias_${group}_errors.txt
	while read line; do
		ary=($line)
		SQB=${ary[0]}
		sample=${ary[2]}
		echo "Sample: ${sample}"
		bampath=${PROJECT_DP3}/processed/${SQB}/output/results_pipeline/${sample}/aligned_hg38/${sample}_sort_dedup.bam
		python ${scripts}/fine_mapping/count_hetalleles.py --bam ${bampath} --hetsites_bed ${vardir}/count_hetalleles_sites_${snp}.bed > ${vardir}/allelic_bias_${group}_${sample}_out.bed 2> ${vardir}/allelic_bias_${group}_${sample}_err.bed
		cat ${vardir}/allelic_bias_${group}_${sample}_out.bed | sed 's/;/\t/g' | awk '{print $6}' | gawk -v sample=${sample} -v A1=${a1} -v A2=${a2} -v pattern="${a1}=([[:digit:]]+),${a2}=([[:digit:]]+)" 'match($0, pattern, a) {print sample, a[1], a[2]}' >> ${vardir}/allelic_bias_${group}_summary.txt
		cat ${vardir}/allelic_bias_${group}_${sample}_err.bed | sed 's/;/\t/g' | awk '{print $6}' | gawk -v sample=${sample} -v A1=${a1} -v A2=${a2} -v pattern="${a1}=([[:digit:]]+),${a2}=([[:digit:]]+)" 'match($0, pattern, a) {print sample, a[1], a[2]}' >> ${vardir}/allelic_bias_${group}_errors.txt
		rm -f ${vardir}/allelic_bias_${group}_${sample}_out.bed
		rm -f ${vardir}/allelic_bias_${group}_${sample}_err.bed
	done < ${vardir}/${snp}_${group}.txt
}

echo "Getting HET counts"
calculateAndSummarize HET
echo "Getting HOMALT counts"
calculateAndSummarize HOMALT
echo "Getting HOMREF counts"
calculateAndSummarize HOMREF
