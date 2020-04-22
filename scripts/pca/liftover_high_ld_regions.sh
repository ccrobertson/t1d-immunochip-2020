#!/bin/bash



resources=/data1/ccr5ju/MEGA/resources
software=/h4/t1/users/ccr5ju/software/

liftOver=${software}/liftOver/liftOver
chain=${software}/liftOver/hg19ToHg38.over.chain.gz
tmp=${resources}/tmp
mkdir ${tmp}

#liftover start positions
awk '{print "chr"$1, $2-1, $2, "R"NR}' ${resources}/high_ld_regions_b37.txt > ${tmp}/high_ld_regions_b37_start.bed 
${liftOver} ${tmp}/high_ld_regions_b37_start.bed ${chain} ${tmp}/high_ld_regions_b38_start.bed ${tmp}/high_ld_regions_b38_start.err -multiple -minMatch=0.99

#liftover end positions
awk '{print "chr"$1, $3-1, $3, "R"NR}' ${resources}/high_ld_regions_b37.txt > ${tmp}/high_ld_regions_b37_end.bed  
${liftOver} ${tmp}/high_ld_regions_b37_end.bed ${chain} ${tmp}/high_ld_regions_b38_end.bed ${tmp}/high_ld_regions_b38_end.err -multiple -minMatch=0.99


#merge start and end positions
join -14 -24 ${tmp}/high_ld_regions_b38_start.bed ${tmp}/high_ld_regions_b38_end.bed | awk '{gsub("chr","",$2); print $2, $3, $7, $1}' > ${resources}/high_ld_regions_b38.txt 
