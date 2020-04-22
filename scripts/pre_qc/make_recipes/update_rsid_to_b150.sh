set -euxo pipefail

##############################################################################################
#
#    Update rs ids to b150
#
##############################################################################################
#
#
#	https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/
#	https://www.biostars.org/p/349284/
#
#	#### get dbSNP 150 file ####
#	cd ${resources}
#	wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/All_20170710.vcf.gz
#	tabix -p vcf All_20170710.vcf.gz

comment "Update rsids to dbsnp150"

cd ${data_derived}

#use tabix to extract snps at updated position in mega data set
input=${data_derived}/${nickname}_raw_nophenoC.bim
vcf=${resources}/All_20170710.vcf.gz
awk  '{print $1, $4"-"$4}' ${input} | awk 'BEGIN {mapToXY[23]="X"; mapToXY[24]="Y"; mapToXY[25]="PAR"; mapToXY[26]="M"} {if ($1>22) {print mapToXY[$1], $2} else {print $0}}' > ichip_snp_positions.txt
tabix ${vcf} -R ichip_snp_positions.txt > ichip_snp_positions_to_dbsnp.txt

#create update file to recode variant ids in mega data
python ${scripts}/update_rsids.py ${input} > rsid_map
awk 'NR>1 {if ($4=="NA") {print $2, $1} else {print $2, $4}}' rsid_map > rsid_update.txt


#update rsids
plink1.9 --bfile ${data_derived}/${nickname}_raw_nophenoC --update-name rsid_update.txt --make-bed --out ${data_derived}/${nickname}_raw_nophenoD
