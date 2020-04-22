set -euxo pipefail


##############################################################################################
#
#    Create Combined Phenotype file
#
##############################################################################################
#	See README_file_dict for more details on files origins
#
#	MEGA PHENOTYPE FILES
#--------------------------------------------------------------------------
#	Cohort		File (unless specified, path is /h4/t1/projects/IMCHIP/IMCHIP_phenotype/data_pheno/*)
#--------------------------------------------------------------------------										
#	nPOD		immunochip_nPOD_sex_race_20130906_pheno.txt 
#	T1DGC		T1DGC_2011_03_Resources_20151029_v1.txt
#	Trialnet	Trialnet_query1_pheno.txt
#	TEDDY		teddy_family_mapping_rev_jul15_pheno.txt
#	other		samplewith_2_pheno_20151030.txt <--- these are duplicates, there are several versions of this file. this seems to be the most recent
#	???			combined_pheno_additional.txt  <--- not sure where this file comes from (see README_file_dict)
# 	CBR			/data1/so4g/cbr/cbrmicro_ic_inclusion_list_20171003.txt
#	Additions	/data1/so4g/ic_mega/ic_additional/IC-mega-additional_samples_20180305.csv
#
#

#create formatted phenofile for CBR controls
echo -e "Cohort\tFamilyID\tSampleID\tFather\tMother\tSex\tT1D" > ${data_derived}/cbr_phenofile_July2018.txt 
awk 'BEGIN {OFS="\t"}; {print "cbr_controls", $1, $1, "0", "0", "0", "1"}' ${data_raw}/phenotype_files/cbr-controls-2018-06-28_formatted.txt >> ${data_derived}/cbr_phenofile_July2018.txt 

#generate combined file from the multiple phenotype files
echo -e "Cohort\tFamilyID\tSampleID\tFather\tMother\tSex\tT1D" > ${phenofile}_tmp

#s/\r//g <-- fixes carriage return
#gsub(/-9/,"0",$6) <--- replaces -9 with 0 in $6 column
#awk 'BEGIN {FS="\t"; OFS="\t"};NR>1{print "nPOD",$1, $2, $3, $4, $5, $6}' ${data_raw}/phenotype_files/nPOD_Ped_info_file.csv | sed 's/ /_/g' >> ${phenofile}_tmp
awk 'BEGIN {FS="\t"; OFS="\t"};NR>1{gsub(/-9/,"0",$6); gsub(/-9/,"0",$7); print $1, $2, $3, $4, $5, $6,$7}' ${data_raw}/phenotype_files/immunochip_nPOD_sex_race_20130906_pheno.txt | sed 's/ /_/g;s/\r//g' >> ${phenofile}_tmp
awk 'BEGIN {FS="\t"; OFS="\t"};NR>1{gsub(/-9/,"0",$6); gsub(/-9/,"0",$7); print $1, $2, $3, $4, $5, $6, $7}' ${data_raw}/phenotype_files/T1DGC_2011_03_Resources_20151029_v1.txt | sed 's/ /_/g;s/\:_/:/g;s/\r//g' >> ${phenofile}_tmp
awk 'BEGIN {FS="\t"; OFS="\t"};NR>1{gsub(/-9/,"0",$6); gsub(/-9/,"0",$7); print $1, $2, $3, $4, $5, $6, $7}' ${data_raw}/phenotype_files/Trialnet_query1_pheno.txt | sed 's/ /_/g;s/\r//g' >> ${phenofile}_tmp
awk 'BEGIN {FS="\t"; OFS="\t"};NR>1{gsub(/-9/,"0",$8); gsub(/-9/,"0",$7); print $1, $2, $3, $6, $5, $8, $7}' ${data_raw}/phenotype_files/teddy_family_mapping_rev_jul15_pheno.txt | sed 's/ /_/g;s/\r//g' >> ${phenofile}_tmp
awk 'BEGIN {FS="\t"; OFS="\t"};NR>1{gsub(/-9/,"0",$6); gsub(/-9/,"0",$7); print $1, $2, $3, $4, $5, $6, $7}' ${data_raw}/phenotype_files/samplewith_2_pheno_20151030.txt | sed 's/ /_/g;s/\r//g' >> ${phenofile}_tmp
awk 'BEGIN {FS="\t"; OFS="\t"};{gsub(/-9/,"0",$6); gsub(/-9/,"0",$7); print $1, $2, $3, $4, $5, $6, $7}' ${data_raw}/phenotype_files/combined_pheno_additional.txt | sed 's/ /_/g;s/\r//g' >> ${phenofile}_tmp
awk 'BEGIN {FS="\t"; OFS="\t"};NR>1{gsub(/-9/,"0",$6); gsub(/-9/,"0",$7); print $1, $2, $3, $4, $5, $6, $7}' ${data_derived}/cbr_phenofile_July2018.txt | sed 's/ /_/g;s/\r//g' >> ${phenofile}_tmp
awk 'BEGIN {FS=","; OFS="\t"};NR>1{gsub(/-9/,"0",$6); gsub(/-9/,"0",$7); print $1, $2, $3, $4, $5, $6, $7}' ${data_raw}/phenotype_files/IC-mega-additional_samples_20180305.csv | sed 's/ /_/g;s/\r//g' >> ${phenofile}_tmp

#fix cohort ids --> define T1DGC groups
awk 'BEGIN {FS="\t"; OFS="\t"};NR==1 {print $1, $2, $3, $4, $5, $6,$7}; NR>1 {\
gsub(/^AP$/,"T1DGC-AP",$1); gsub(/^T1DGC:Asia-Pacific_DNA_Repository$/,"T1DGC-AP",$1); \
gsub(/^EUR$/,"T1DGC-EUR",$1); gsub(/^Existing_Source_to_T1DGC:_SWEDISH$/,"T1DGC-EUR",$1); gsub(/^Existing_Source_to_T1DGC:SWEDISH$/,"T1DGC-EUR",$1); gsub(/^T1DGC:European_DNA_Repository$/,"T1DGC-EUR",$1); \
gsub(/^NA$/,"T1DGC-NA",$1); gsub(/^T1DGC:North_American_DNA_Repository$/,"T1DGC-NA",$1); \
gsub(/^UK$/,"T1DGC-UK",$1); gsub(/^T1DGC:United_Kingdom_DNA_Repository$/,"T1DGC-UK",$1); gsub(/^T1DGC-WARREN$/,"T1DGC-UK",$1); \
gsub(/^58BC$/,"BC58",$1); gsub(/^B58C$/,"BC58",$1); gsub(/^Sanger_BC58$/,"BC58",$1); \
print $1, $2, $3, $4, $5, $6,$7}' ${phenofile}_tmp > ${phenofile}_tmp2

#assign families with missing subcohort to correct T1DGC subcohort 
#awk '$1=="T1DGC" {print $2}' ${phenofile}_tmp2 > ${data_derived}/T1DGC_unassigned.txt
#awk -v file=${data_derived}/T1DGC_unassigned.txt 'BEGIN { while(getline<file){a[$1]=1} } $2 in a {print $2, $1}' ${phenofile}_tmp2 
awk 'BEGIN {FS="\t"; OFS="\t"} { if ($2=="278002") $1="T1DGC-EUR"; if ($2=="489652") $1="T1DGC-NA"; if ($2=="531896") $1="T1DGC-UK"; ; print $0}' ${phenofile}_tmp2 > ${phenofile}

