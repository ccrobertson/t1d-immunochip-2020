set -euxo pipefail

##############################################################################################
#
#    Update duplicate sample ids 
#
##############################################################################################
#
#	This must be done before merging files or performing any other commands in plink
#	because plink will throw an error if it finds more than one sample with same ID
#
#
comment "Update duplicate sample ids" 


cd ${data_derived}
export nicknames=(${nickname_1} ${nickname_2} ${nickname_3} ${nickname_4})

############# CHECK FOR DUPLICATES (and triplicates) within and between files before merging
echo -e "Now checking for duplicate and triplicate sample IDs." 

echo -n > ${nickname}0_sample_ids.txt
for nickname_X in ${nicknames[@]}; do
	awk -v nickname=${nickname_X} '{print $0,nickname}' ${nickname_X}_raw/${nickname_X}0.fam >> ${nickname}0_sample_ids.txt
done
		
#get duplicate update file
awk '{OFS="\t"; print $1,$2, $7}' ${nickname}0_sample_ids.txt | awk 'BEGIN {a[$2]=0} {a[$2]++; if(a[$2]==2) {print $1,$2,$1"_2",$2"_2", $3}}' > duplicate_sample_update.txt
for nickname_X in ${nicknames[@]}; do
	awk -v nickname=${nickname_X} '$5==nickname {print $1,$2,$3,$4}' duplicate_sample_update.txt > ${nickname_X}_raw/${nickname_X}0_duplicate_sample_updates.txt
done

#get triplicate update file
awk '{OFS="\t"; print $1,$2, $7}' ${nickname}0_sample_ids.txt | awk '{a[$2]++; if(a[$2]==3) {print $1,$2,$1"_3",$2"_3", $3}}' > triplicate_sample_update.txt
for nickname_X in ${nicknames[@]}; do	
	awk -v nickname=${nickname_X} '$5==nickname {print $1,$2,$3,$4}' triplicate_sample_update.txt > ${nickname_X}_raw/${nickname_X}0_triplicate_sample_updates.txt
done


############### UPDATE PLINK FILES
for nickname_X in ${nicknames[@]}; do	
	awk -v file=${nickname_X}_raw/${nickname_X}0_duplicate_sample_updates.txt 'BEGIN { while( getline <file ) a[$1"\t"$2]=$3"\t"$4 } {b[$1"\t"$2]++; if( ($1"\t"$2 in a) && (b[$1"\t"$2]==1)) {print a[$1"\t"$2]"\t"$3"\t"$4"\t"$5"\t"$6} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' ${nickname_X}_raw/${nickname_X}0.fam > ${nickname_X}_raw/${nickname_X}0_update1.fam
	awk -v file=${nickname_X}_raw/${nickname_X}0_triplicate_sample_updates.txt 'BEGIN { while( getline <file ) a[$1"\t"$2]=$3"\t"$4 } {b[$1"\t"$2]++; if( ($1"\t"$2 in a) && (b[$1"\t"$2]==1)) {print a[$1"\t"$2]"\t"$3"\t"$4"\t"$5"\t"$6} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' ${nickname_X}_raw/${nickname_X}0_update1.fam > ${nickname_X}_raw/${nickname_X}0_update2.fam
	plink1.9 --bed ${nickname_X}_raw/${nickname_X}0.bed --bim ${nickname_X}_raw/${nickname_X}0.bim --fam ${nickname_X}_raw/${nickname_X}0_update2.fam --make-bed --out ${nickname_X}_raw/${nickname_X}1 
done


############### Create export indicator file --> Which subject came from which export?
echo -e "FID\tIID\tExport" > ${data_derived}/genome_studio_export_map.txt
for nickname_X in ${nicknames[@]}; do	
	
	awk -v export=${nickname_X} '{print $1, $2, export}' ${data_derived}/${nickname_X}_raw/${nickname_X}1.fam >> ${data_derived}/genome_studio_export_map.txt

done



