set -euxo pipefail

##############################################################################################
#
#    Update sample ids using provided sample id maps (from cohort studies)
#
##############################################################################################	
comment "Update ids using provided sample id maps" 


	
	
echo -e "Now updating IDs to match phenotype file, and removing necessary samples." 

#fix carriage return, update ids when multiple old ids map to same new id (i.e. update duplicate ids by appending "_2")
sed 's/\r//g' ${data_raw}/phenotype_files/updateID_6.txt | awk '$1!=$3 || $2!=$4' | sort -u > ${data_derived}/updateID_0.txt
awk -v file=${data_derived}/updateID_0.txt 'BEGIN { while( getline<file ) { a[$1"\t"$2]=$3"\t"$4; } }   { if ($1"\t"$2 in a) {print $1"\t"$2"\t"a[$1"\t"$2]} else {print $1"\t"$2"\t"$1"\t"$2} } ' ${data_derived}/${nickname}_raw_nophenoD.fam > ${data_derived}/updateID_1.txt
awk 'BEGIN {a[$4]=0} { a[$4]++; if(a[$4]==2) {print $1,$2,$3"_2",$4"_2"} else {print $1,$2,$3,$4}}'  ${data_derived}/updateID_1.txt >  ${data_derived}/updateID_2.txt
plink1.9 --bfile ${data_derived}/${nickname}_raw_nophenoD --update-ids ${data_derived}/updateID_2.txt --make-bed --out ${data_derived}/${nickname}_raw_nophenoE 

#how many ids are missing in phenofile now?
awk -v file1=${phenofile} 'BEGIN { while( getline <file1 ) a[$3]=1} !($2 in a) {print $2}'  ${data_derived}/${nickname}_raw_nophenoE.fam > ${data_raw}/phenotype_files/missing_pheno_${nickname}_raw_nophenoE.txt



