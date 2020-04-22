#!/bin/bash

set -o errexit
set -o pipefail

###############################################################################
# Usage: sh  sexqc.sh {inputfile}
##############################################################################
#(1=male; 2=female; other=unknown)
 
plinkraw=$1
sexqc_folder=${sexqc}/$2
qclists_folder=${qclists}/$2
prefix=sexqc

echo -e "########################################################################################################################################"
echo -e
echo -e "Sex QC: checking for sex chromosome anomalies and inconsistencies between reported sex and genotype data"
echo -e
echo -e "Input files: ${plinkraw}"
echo -e "Performed on:"
date
echo -e
echo -e "########################################################################################################################################"


##############################################################################################
echo -e
echo -e "\n------------------"
echo "Retrieving data: copying ${plinkraw} to ${sexqc}"
echo -e "------------------\n"
echo -e

mkdir -p ${sexqc_folder}
mkdir -p ${qclists_folder}

cp ${plinkraw}.bim ${sexqc_folder}/${prefix}1.bim
cp ${plinkraw}.bed ${sexqc_folder}/${prefix}1.bed
cp ${plinkraw}.fam ${sexqc_folder}/${prefix}1.fam
cp ${plinkraw}.log ${sexqc_folder}/${prefix}1.log

cd ${sexqc_folder}



##############################################################################################
echo -e
echo -e "\n------------------"
echo "STEP 1: prelim-SNP QC -- remove missing non-Y SNPs and all monomorphic SNPs"
echo -e "------------------\n"
echo -e
	

king213 -b ${prefix}1.bed --bySNP --prefix ${prefix}1

#find non-Y SNPs with call rate <0.8
awk '$11 < 0.8 && $2 != "Y"' ${prefix}1bySNP.txt | awk '{print $1}' > missingSNP.txt

#find monomorphic snps
awk '$8+$9==0 || $9+$10==0' ${prefix}1bySNP.txt | awk '{print $1}' > monomorphicSNP.txt

awk '{print $1, "CallRateLessThan80"}' missingSNP.txt >> snptoberemoved.txt
awk 'NR>1{print $1, "Monomorphic"}' monomorphicSNP.txt >> snptoberemoved.txt

plink1.9 --bfile ${prefix}1 --exclude snptoberemoved.txt --make-bed --out ${prefix}2


#######################################################################################################
echo -e
echo -e "\n------------------"
echo "Define sexCheck function"
echo -e "------------------\n"
echo -e


sexCheck () {
	local infile=$1
	local dir=$2
	local prefix=$3
	
	king213 -b ${infile}.bed --bysample --prefix ${dir}/${prefix} 
	
	#Calculate y SNP counts
	ysnps=$(awk '$1==24' ${infile}.bim | wc -l)
	half_ysnps=$(echo ${ysnps} / 2 | bc )
	third_ysnps=$(echo ${ysnps} / 3 | bc )
	twothird_ysnps=$(echo "2 * ${third_ysnps}" | bc)
	sixth_ysnps=$(echo ${ysnps} / 6 | bc)
	fivesixth_ysnps=$(echo ${ysnps} - ${sixth_ysnps} | bc)
	
	#Define thresholds for detecting anomalies
	xHetThresh=0.08
	yThreshLow=${third_ysnps}
	yThreshHigh=${twothird_ysnps} 
	
	cd ${dir}
	
	#Visualize 
	Rscript ${scripts}/visualize_sexqc.R ${prefix} ${xHetThresh} ${yThreshLow} ${yThreshHigh}
	
	#Kleinfelter (XXY)
	#X-heterozygosity is greater than 8% (so probably more than 1 X chromosome) but more than 2/3 y snps are non-missing (so there is probably a y chromosome)
	#awk 'NR==1 {for(i = 1; i <= NF; i++) {a[$i]=i}} {print $a["xHeterozygosity"], $a["xHeterozygosity"]}' SEXCHECK1/sexcheck1bySample.txt |head
	awk -v yThreshHigh=${yThreshHigh} -v xHetThresh=${xHetThresh} 'NR>1 && $10 > xHetThresh && $11 > yThreshHigh' ${prefix}bySample.txt | awk '{print $1,$2}' > XXYerror.nb

	#Midsex (possibly mosaic, possibly Keinfelters) ==> ambiguous frequency of non-missing y snps
	awk -v yThreshLow=${yThreshLow} -v yThreshHigh=${yThreshHigh} 'NR>1 && $11 > yThreshLow && $11 < yThreshHigh' ${prefix}bySample.txt | awk '{print $1,$2}' > midsexerror.nb

	#Turner Syndrome (X0) or other error ==> low X-heterozygosity (so probably only one X chromosome) but less than 1/3 of y-snps are non-missing (so there probably isn't a y chromosome) 
	awk -v yThreshLow=${yThreshLow} -v xHetThresh=${xHetThresh} 'NR>1 && $10 < xHetThresh && $11 < yThreshLow' ${prefix}bySample.txt | awk '{print $1,$2}' > X0error.nb
	
	#Flagged sex chromosome anomalies
	awk '{print $1, $2, "MidSexError"}' midsexerror.nb > sex_anomalies.txt
	awk '{print $1, $2, "X0error"}' X0error.nb >> sex_anomalies.txt
	awk '{print $1, $2, "XXYerror"}' XXYerror.nb >> sex_anomalies.txt
	
	echo -e "FID\tIID\tPEDSEX\tSNPSEX\tF\tN_ySNP" > ${prefix}_sexupdate.txt
	awk -v xHetThresh=${xHetThresh} -v yThreshLow=${yThreshLow} -v yThreshHigh=${yThreshHigh} 'NR>1 {if ($10>xHetThresh && $11<yThreshLow) {sex="2"} else if ($10<xHetThresh && $11>yThreshHigh) {sex="1"} else {sex="0"}; print $1, $2, $5, sex, $10, $11}'  ${prefix}bySample.txt >> ${prefix}_sexupdate.txt		
	
	#echo -e "FID\tIID\tPEDSEX\tSNPSEX\tF\tN_ySNP" > sexupdate.txt
	#awk -v xHetThresh=0.08 -v yThreshLow=300 -v yThreshHigh=800 'NR>1 {if ($10>xHetThresh && $11<yThreshLow) {sex="2"} else if ($10<xHetThresh && $11>yThreshHigh) {sex="1"} else {sex="0"}; print $1, $2, $5, sex, $10, $11}'  SEXQC/sexqc5bySample.txt >> sexupdate.txt

	awk '$3==0 && $4==0 {print $1, $2, $3, $4}' ${prefix}_sexupdate.txt > missing_sex_ambiguous.txt
	awk '$3==0 && $4!=0 {print $1, $2, $3, $4}' ${prefix}_sexupdate.txt > missing_sex_imputed.txt
	awk '$3!=0 && $4==0 {print $1, $2, $3, $4}' ${prefix}_sexupdate.txt > non_missing_sex_ambiguous.txt
	awk '($3==1 && $4==2) || ($3==2 && $4==1) {print $1, $2, $3, $4}' ${prefix}_sexupdate.txt > non_missing_sex_error.txt

	cd ..	
}


#######################################################################################################
echo -e
echo -e "\n------------------"
echo "STEP 2: flag sex anomalies"
echo -e "------------------\n"
echo -e


mkdir -p SEXCHECK1
sexCheck ${prefix}2 SEXCHECK1 sexcheck1

#remove anomalies
plink1.9 --bfile ${prefix}2 --remove SEXCHECK1/sex_anomalies.txt --make-bed --out ${prefix}3


##############################################################################################
echo -e
echo -e "\n------------------"
echo "STEP 3: Y-SNP filtering"
echo -e "------------------\n"
echo -e


#calculate snp-level stats
king213 -b ${prefix}3.bed --bySNP --prefix ${prefix}3

#Determine the relative ratio of males to males+females
males=$(awk '$5==1' ${prefix}3.fam | wc -l)
females=$(awk '$5==2' ${prefix}3.fam | wc -l)
mratio=$(echo ${males}/$((${males}+${females}))|bc -l)
ratio=$(echo ${mratio}*1.1|bc -l)
echo "The ratio of males:(males+females) is " ${mratio} ", which means we will cut Y SNPs with call rates greater than " ${ratio}

#Find Y snps with call rate greater than expected based on proportion male
awk -v r=${ratio} '$2=="Y" && $11 > r' ${prefix}3bySNP.txt | awk '{print $1, "HighCallRate"}' > poorYSNP.txt

#Find Y snps with high call rates in females
plink1.9 --bfile ${prefix}3 --filter-females --out ${prefix}3female --make-bed
king213 -b ${prefix}3female.bed --bySNP --prefix ${prefix}3female
awk '$2=="Y" && $11 > 0.04' ${prefix}3femalebySNP.txt | awk '{print $1 , "HighCallRateInFemales"}' >> poorYSNP.txt

#Find Y snps with low call rates in males
plink1.9 --bfile ${prefix}3 --filter-males --out ${prefix}3male --make-bed
king213 -b ${prefix}3male.bed --bySNP --prefix ${prefix}3male
awk '$2=="Y" && $11 < 0.80' ${prefix}3malebySNP.txt | awk '{print $1 , "LowCallRateInMales"}' >> poorYSNP.txt


#Find Y snps with high heterozygosity (should be homozygous everywhere!)
#any Y SNPs with a number of heterozygous samples greater than one tenth the number of males, should be removed
males=$(awk '$5==1' ${prefix}3.fam | wc -l)
males10=$(echo "$males *.1" | bc -l)
echo "Any Y SNPs with a number of heterozygous samples greater than one tenth the number of males, should be removed. We have" ${males} "males"
awk -v males10=${males10} '$2=="Y" && $9>males10' ${prefix}3bySNP.txt | awk '{print $1, "HighHeterozygosity"}' >> poorYSNP.txt
#awk 'BEGIN {a[$1]=0} {if (a[$1]==0) {a[$1]=$2} else {a[$1]=a[$1]","$2}} END {for (i in a) {print i, a[i]}}' poorYSNP.txt

#remove poor Y snps
awk '{print $1}' poorYSNP.txt | sort | uniq > poorYSNP_remove.txt 
plink1.9 --bfile ${prefix}3 --exclude poorYSNP_remove.txt  --make-bed --out ${prefix}4 


#######################################################################################################
echo -e
echo -e "\n------------------"
echo "STEP 4: flag sex inconsistencies and impute sex when missing"
echo -e "------------------\n"
echo -e

mkdir -p SEXCHECK2
sexCheck ${prefix}4 SEXCHECK2 sexcheck2
	 		 		 		 		 		 		 		 		 		 		 		 		 	 		 		 		 		 		 		 		 		 		 		 		 	 		 		 		 		 		 		 		 		 		 		 		 	 		 		 		 		 		 		 		 		 		 		 		 	 		 		 		 		 		 	

##############################################################################################
echo -e
echo -e "\n------------------"
echo "STEP 5: copy sex qc lists to ${qclists_folder}"
echo -e "------------------\n"
echo -e


cp poorYSNP_remove.txt ${qclists_folder}/
cp -r SEXCHECK1 ${qclists_folder}/
cp -f -r SEXCHECK2 ${qclists_folder}/


