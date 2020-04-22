#!/bin/bash
#!/usr/bin/bash

#Cassie Robertson
#January 2018
#This script is adapted from Wei-Min's script: /t121/jdrfdn/projects/JDRFDN/workfiles/Wei-Min/script/makedata.sh 
#CHANGES:
#plink --> plink1.9
#removed --no-web flag from plink commands


### This script transforms output from Genome Studio (Full Data Tables in Top Allele format) into PLINK format.


##############################################################################
# Usage: sh  make_data.sh {DATAFILE} {nickname}
##############################################################################

##################################################################################################################
#Arguments
# DATAFILE : the path to the full data table exported from Genome Studio.
#		It is expected that the table will contain more than one column per sample. Including Top Allele and GType
# nickname : string, the alias for the project, used in all filenames
#
########################################################################################################



# Definining/Reading user specified script parameters/variables/values
DATAFILE=$1
nickname=$2



echo -e "########################################################################################################################################"
echo -e
echo -e "# Title: MAKE DATA - Standard GWAS QC Pipeline"
echo -e "#"
 echo -e "#"
echo -e "# Disc: Transforms output from Genome Studio into plink files."
echo -e "#"
echo -e "# Note: PLINK output files will NOT contain sex, status, or family information, those are added during pheno_inc.sh"
echo -e "# Note: DATAFILE expected to be in Top Allele format."
echo -e "# Note: DATAFILE expected to contain a GType column."
echo -e "# Note: DATAFILE expected to have more than one column per sample."
echo -e "# Note: The first two '.Top Allele' columns are expected to occur within the first 20 columns of DATAFILE."
echo -e "#"
echo -e "# Usage: sh  ~/SHcode/make_data.sh {DATAFILE} {nickname}"
echo -e "#"
echo -e "# See Script for option/parameter details"
echo -e "#"
echo -e "# by:  jbonnie"
echo -e "# date: 10.25.13"
echo -e
echo -e "Full Data Table: ${DATAFILE}"
echo -e "Study Alias: ${nickname}"
echo -e "Performed on:"
date
echo -e
echo -e "########################################################################################################################################"




project_folder=$(pwd)


#It's easier to make variables to indicate these subfolders and files

raw_folder=${project_folder}/${nickname}_raw

echo -e
echo -e "\n------------------"
echo "Creating Raw Data Folder"
echo -e "------------------\n"
echo -e
mkdir ${raw_folder}
cd ${raw_folder}

echo -e
echo -e "\n------------------"
echo "Storing column numbers from the Full Data Table &"
echo "Determining Uniform Number of Columns for Each Sample"
echo -e "------------------\n"
echo -e


#Lets read some columns numbers!


#Locate the first two occurances of the pattern ".Top Allele", and then find the difference.
#This tells us the uniform number of columns between each ".Top Allele"!

topcols=$(cut -f1-20 $DATAFILE| sed "s/\r//g" | head -n1| sed 's/\t/\n/g'| awk '$0 ~/.Top Allele/ {print NR}')
top1=$( echo ${topcols} | awk '{print $1}')
top2=$( echo ${topcols} | awk '{print $2}')
dif=$((${top2} - ${top1}))

#Locate the first occurance of the pattern ".Top Allele", that is the first data column of interest.
#top1=$(cut -f1-8 $DATAFILE| sed "s/\r//g" | head -n1| sed 's/\t/\n/g'| awk '$0 ~/.Top Allele/ {print NR}')

#The second occurance of the pattern ".Top Allele" cannot be more than 10 columns in, right? Lets Find it!
#Once we know where it is we will know the uniform number of columns between each ".Top Allele"!
#dif=$(cut -f$((${top1}+1))-10 $DATAFILE| sed "s/\r//g" | head -n1| sed 's/\t/\n/g'| awk '$0 ~/.Top Allele/ {print NR}')
#top2=$((${top1} + ${dif}))

#We assume that the preliminary columns also occur within the first 20 columns

chr=$(cut -f1-20 $DATAFILE| sed "s/\r//g" | head -n1| sed 's/\t/\n/g'| awk '$0 ~/^Chr$/ {print NR}')
pos=$(cut -f1-20 $DATAFILE| sed "s/\r//g" | head -n1| sed 's/\t/\n/g'| awk '$0 ~/^Position$/ {print NR}')
snp=$(cut -f1-20 $DATAFILE| sed "s/\r//g" | head -n1| sed 's/\t/\n/g'| awk '$0 ~/^Name$/ {print NR}')


echo -e
echo -e "Column Numbers"
echo -e "Chromosome : ${chr}"
echo -e "SNP : ${snp}"
echo -e "SNP POSITION : ${pos}"
echo -e "First Top Allele : ${top1}"
echo -e "Second Top Allele : ${top2}"
echo -e "Difference (Number of Columns per Sample) : ${dif}"
echo -e




echo -e
echo -e "\n------------------"
echo "Creating PLINK TPED: Transferring Allele Data"
echo -e "------------------\n"
echo -e


## For each .Top Allele column, print the allele information, replacing "--" with "0 0" and adding spaces between the alleles
## Also, for some reason, flip the alleles, this was how it was done in inherited skeleton, so it has been left that way.

awk 'NR>1' $DATAFILE | sed "s/\r//g" |\
  awk -v top1=${top1} -v top2=${top2} -v dif=${dif} '{printf("%s", $top1); for(i=top2;i<=NF;i+=dif) printf(" %s", $i);printf("\n");}' |\
  sed 's/--/0 0/g' |\
 sed 's/\([A-Z]\)\([A-Z]\)/\2 \1/g' > topAllele_tped.tmp


echo -e
echo -e "\n------------------"
echo "Creating PLINK TPED: Transferring SNP Information"
echo -e "------------------\n"
echo -e


# The first four columns of the tped hold the SNP info from the table with "0" place holders in the "Genetic Distance" column (column 3)

awk 'NR>1' $DATAFILE | sed "s/\r//g" | awk -v chr=${chr} -v pos=${pos} -v snp=${snp} '{printf("%s %s 0 %d\n", $chr,$snp,$pos);}' > snpInfo_tped.tmp

echo -e
echo -e "\n------------------"
echo "Creating PLINK TPED: Pasting TMP files"
echo -e "------------------\n"
echo -e

# create tped by pasting snpInfo_tped.tmp and topAllele_tped.tmp
paste -d" " snpInfo_tped.tmp topAllele_tped.tmp > data.tped



echo -e
echo -e "\n------------------"
echo "Creating PLINK TFAM"
echo -e "------------------\n"
echo -e

# -9 means unknown/missing affection status
# take header, replace tabs with newlines, and strip windows newlines; take only the lines containing 'GType' and then strip '.GType' to leave only sampleIDs;
# NOT SURE WHAT THIS IF STATEMENT WOULD CATCH.... maybe a colname that is "FID IID.GType"... 
# regardless, outside of the if statement the IID is doubled to also be FID and parents and affection status are defaulted to their missing values

head -1 $DATAFILE | sed "s/\r//g" | sed "s/\t/\n/g" |\
 grep "GType" | sed "s/.GType//g" |\
 awk '{if(NF==2){print $1,$2,0,0,0,-9}else{print $1,$1,0,0,0,-9}}' > data.tfam




echo -e
echo -e "\n------------------"
echo "Creating PLINK BINARY FILES"
echo -e "------------------\n"
echo -e


plink1.9 --tfile data --make-bed --out ${nickname}0


echo -e
echo -e "\n------------------"
echo "Tar/Zip"
echo -e "------------------\n"
echo -e

#rmlist=$(ls data.* *.tmp *.hh *.nosex)
#rm -f $rmlist

tar --create -f raw.tar *

gzip raw.tar


echo 'convert_genome_studio_raw.sh completed: '$(date)
