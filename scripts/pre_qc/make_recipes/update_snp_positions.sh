set -euxo pipefail

##############################################################################################
#
#    Update SNP positions (by aligning probe sequences to hg19)
#
##############################################################################################
#
#	This must be done before merging files because not all exports were generated with
#	the same manifest, which means different files have different snps positions
#
#	Code and files used to generate mymapInfo and mydrop:
# 
#	${PROJECT_MEGA}/scripts/probe_alignment/README_map_probes
#	

comment "Update SNP positions to hg19" 

cd ${data_derived}

for nickname_X in ${nicknames[@]}; do
	
	echo "Updating ${nickname_X}:"  
	
	#### USING MY MAPPING
	echo -e "1. Dropping indels and otherwise unmapped variants"  
	plink1.9 --bfile ${nickname_X}_raw/${nickname_X}1 --exclude ${data_derived}/mydrop --out ${nickname_X}_raw/${nickname_X}2 --make-bed  
	
	echo -e "2. Updating chromosome"  
	plink1.9 --bfile ${nickname_X}_raw/${nickname_X}2 --update-chr ${data_derived}/mymap_update_chr --out ${nickname_X}_raw/${nickname_X}3 --make-bed  
	
	echo -e "3. Updating position"  
	plink1.9 --bfile ${nickname_X}_raw/${nickname_X}3 --update-map ${data_derived}/mymap_update_pos --out ${nickname_X}_raw/${nickname_X}4 --make-bed  

	
	#### USING QUINLAN MAPPING
#	echo -e "1. Updating chromosome"    
#	plink1.9 --bfile ${nickname_X}_raw/${nickname_X}1 --update-chr ${data_derived}/quinlan_update_chr --out ${nickname_X}_raw/${nickname_X}3 --make-bed  
#
#	echo -e "2. Updating position"    
#	plink1.9 --bfile ${nickname_X}_raw/${nickname_X}3 --update-map ${data_derived}/quinlan_update_pos --out ${nickname_X}_raw/${nickname_X}4 --make-bed  
	
	
	#### USING QUINLAN MAPPING -- NO CHR UPDATE
	#echo -e "1. Dropping indels and otherwise unmapped variants"  
	#plink1.9 --bfile ${nickname_X}_raw/${nickname_X}1 --exclude ${data_derived}/mydrop --out ${nickname_X}_raw/${nickname_X}3 --make-bed    

	#echo -e "2. Updating position"    
	#plink1.9 --bfile ${nickname_X}_raw/${nickname_X}3 --update-map ${data_derived}/quinlan_update_pos --out ${nickname_X}_raw/${nickname_X}4 --make-bed  
	
	
	echo -e " " 
	
done



