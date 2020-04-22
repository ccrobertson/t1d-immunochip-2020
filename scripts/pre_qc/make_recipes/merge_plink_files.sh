set -euxo pipefail

##############################################################################################
#
#    Merge plink files
#
##############################################################################################
comment "Merge plink files" 


cd ${data_derived}
echo -e "Now merging plink files"  
echo -n > plink_merge_list.txt
for nickname_X in ${nicknames[@]}; do	
	echo ${nickname_X}_raw/${nickname_X}4 >> plink_merge_list.txt
done
plink1.9 --merge-list plink_merge_list.txt --out ${nickname}_raw_nophenoA 




