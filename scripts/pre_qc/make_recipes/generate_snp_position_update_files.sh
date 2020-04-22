#!/bin/bash

set -o errexit
set -o pipefail

cd ${data_derived}

#create mappable and unmappable probe lists (note, "unmappable" includes indels)
awk 'BEGIN {FS="\t"} NR>2 && ($3=="[NA, NA, NA, NA]" || $2=="[D/I]" || $2=="[I/D]")' ${data_raw}/snpqc_files/mymapInfo > my_unmappable 
awk 'BEGIN {FS="\t"} NR>2 && !($3=="[NA, NA, NA, NA]" || $2=="[D/I]" || $2=="[I/D]")' ${data_raw}/snpqc_files/mymapInfo > my_mappable 

#create variant exclusion file
vartodrop=mydrop
> ${vartodrop}

#add variants that don't match beteen mega exports (there is only one -- 'seq-VH-220Z' in main is 'seq-VH-22Z' in jhu)
awk -v file=mega_jhu_raw/mega_jhu0.bim  'BEGIN {while(getline<file){a[$2]=1}} !($2 in a) {print $2}' mega_main_raw/mega_main0.bim >> ${vartodrop}
awk -v file=mega_main_raw/mega_main0.bim  'BEGIN {while(getline<file){a[$2]=1}} !($2 in a) {print $2}' mega_jhu_raw/mega_jhu0.bim >> ${vartodrop}

#add variants without a unique mapping
awk '{print $1}' my_unmappable >> ${vartodrop} 

#create MASTER map file
cat my_mappable | sed 's/\[//g; s/\]//g' | awk 'BEGIN {FS="\t"} {split($3, a, ","); print $1, a[2], a[4]}' | sed -e 's/[ \t]+/\t/g' | awk -v file=${vartodrop} 'BEGIN {while(getline<file){a[$1]=1}} !($1 in a)' > mymap  

#extract variants where map is different from manifest/bim (i.e. only positions that need to be updated)
cat mymap | awk -v file=mega_main_raw/mega_main0.bim 'BEGIN {while(getline<file){chrpos[$2]=$1":"$4}} {gsub("chr","", $2); if (chrpos[$1]!=$2":"$3)  print $0, chrpos[$1]}' > mymap_update

#create chr and pos map update files
awk '{print $1, $2}' mymap_update > mymap_update_chr 
awk '{print $1, $3}' mymap_update > mymap_update_pos
awk 'BEGIN {FS="\t"} NR>1 {print $7, $1}' ${data_raw}/snpqc_files/quinlin_mapInfo > quinlan_update_chr
awk 'BEGIN {FS="\t"} NR>1 {print $7, $3}' ${data_raw}/snpqc_files/quinlin_mapInfo > quinlan_update_pos








#LOOK AT SOME OF THE DISCREPANCIES BETWEEN MYMAP AND QUINLIN MAP
#cat my_mappable | sed 's/\[//g; s/\]//g' | awk 'BEGIN {FS="\t"} {split($4, a, ","); print $1, a[2], a[4], "chr"$6, $8}' | awk '$2":"$3!=$4":"$5 {print $1, $2":"$3, $4":"$5}' > mymap_vs_quinlin_discrepancies.txt 
#findProbeSeq () {
#	local var=$1
#	line=$(awk -v var=${var} '$1==">"var {print NR+1}' manifest_A_alleleA.fa)
#	head -n $line manifest_A_alleleA.fa | tail -n 2
#}
#findProbeSeq "imm_1_2475164"
#findProbeSeq "rs45540544"
#
#grep -f <(cat /nv/vol185/MEGA/probe_alignment/mymap_vs_quinlin_discrepancies.txt | grep chrX| awk '{print $1}') manifest_A_alleleA.psl | awk '{print $14, $16, $17}' > problem_chrX.txt
#awk '$1=="chrX" && $2>60001 && $2<2699520' problem_chrX.txt > par1_X.txt
#awk '$1=="chrX" && $2>154931044 && $2<155260560' problem_chrX.txt > par2_X.txt
#
#grep -f <(cat /nv/vol185/MEGA/probe_alignment/mymap_vs_quinlin_discrepancies.txt | grep chrY| awk '{print $1}') manifest_A_alleleA.psl | awk '{print $14, $16, $17}' > problem_chrY.txt
#awk '$1=="chrY" && $2>10001 && $2<2649520' problem_chrY.txt > par1_Y.txt
#awk '$1=="chrY" && $2>59034050 && $2<59363566' problem_chrY.txt > par2_Y.txt
