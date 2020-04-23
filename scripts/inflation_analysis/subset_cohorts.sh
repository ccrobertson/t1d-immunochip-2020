#!/bin/bash


module load htslib


getDataSubSets () {

	group=$1
	subset=$2
	echo ${subset}
	for ((i=1; i<=22; i+=1)); do
		echo -e "chr${i}"
		sbatch --output	${inflation}/subdata/${subset}/extract_${group}_sub${subset}.chr${i}.log ${scripts}/inflation_analysis/extract_samples.slurm ${group} ${subset} ${i}
	done
	echo " "

}


getDataSets () {

	group=$1
	getDataSubSets ${group} 1
	getDataSubSets ${group} 2
	getDataSubSets ${group} 3
	getDataSubSets ${group} 4
	getDataSubSets ${group} 5

}

getDataSets "AFR"
getDataSets "EUR"
getDataSets "AMR"
getDataSets "FIN"
