set -euxo pipefail

##############################################################################################
#
#    Generate Plink Files
#
##############################################################################################
comment "Generate plink files from Genome Studio exports" 
echo -e "See ${data_derived}/${nickname}_makeraw.log" 


DATAFILE_1=/data1/so4g/ic_mega/10142015-Complete_Immunochip_Full_Data_Table.txt
export nickname_1=mega_main

DATAFILE_2=/data1/so4g/cbr/07122018-Cambridge_Immuno_Additions_Full_Data_Table.txt
export nickname_2=mega_cbr

DATAFILE_3=/m/jdrfdn_scratch/users/projects/IMCHIP/20190320_JHU/03192019-GSLcreated_JohnsHopkins_Project_Full_Data_Table.txt
export nickname_3=mega_jhu

DATAFILE_4=/m/jdrfdn_scratch/users/projects/IMCHIP/20190315/03152019-Immunoplates760-768_Full_Data_Table.txt
export nickname_4=mega_additions


nicknames=(${nickname_1} ${nickname_2} ${nickname_3} ${nickname_4})
DATAFILES=(${DATAFILE_1} ${DATAFILE_2} ${DATAFILE_3} ${DATAFILE_4})



### Create plink files from the full data exports (only if they don't already exist)
cd ${data_derived}
echo -e $(date) > ${data_derived}/${nickname}_makeraw.log

len=${#nicknames[@]}
for (( i=0; i<$len; i++ )); do
	nickname_X=${nicknames[${i}]}
	DATAFILE_X=${DATAFILES[${i}]}
	if [ ! -f ${data_derived}/${nickname_X}_raw/${nickname_X}0.fam ]; then
		echo -e "Making ${nickname_X}_raw" 
		echo ${DATAFILE_X} 
		echo ""  
		bash ${scripts}/convert_genome_studio_raw.sh ${DATAFILE_X} ${nickname_X} >> ${data_derived}/${nickname}_makeraw.log
	fi
done


#check for samples with missing phenotypes
for nickname_X in ${nicknames[@]}; do
	awk -v file1=${phenofile} 'BEGIN { while( getline <file1 ) a[$3]=1} !($2 in a) {print $2}'  ${data_derived}/${nickname_X}_raw/${nickname_X}0.fam > ${data_raw}/missing_pheno_${nickname_X}.txt
done



