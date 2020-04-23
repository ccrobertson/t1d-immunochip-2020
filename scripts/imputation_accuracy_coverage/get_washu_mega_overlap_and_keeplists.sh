
#get corresponding mega ids, cohorts, and ancestry cluster
awk -v file1=${washu}/washu_batch1_sample_sheet.txt  -v file2=${pca}/mega_cluster_assignments.txt 'BEGIN { while( getline<file1 ) {a[$2]=$1}; while( getline<file2 ) {b[$2]=$4} } $3 in a {print "H_VD-"a[$3], $2, $3, $1, b[$3]}' ${phenofile} > ${impbench}/washu_batch1_in_mega.txt
awk -v file1=${washu}/washu_batch2_sample_sheet.txt  -v file2=${pca}/mega_cluster_assignments.txt 'BEGIN { while( getline<file1 ) {a[$2]=$1}; while( getline<file2 ) {b[$2]=$4} } $3 in a {print "H_VD-"a[$3], $2, $3, $1, b[$3]}' ${phenofile} > ${impbench}/washu_batch2_in_mega.txt
awk '{print $5}' ${impbench}/washu_batch1_in_mega.txt| sort | uniq -c
awk '{print $5}' ${impbench}/washu_batch2_in_mega.txt| sort | uniq -c

#get keep lists (for MEGA) by ancestry
awk '$5=="EUR" {print $3}' ${impbench}/washu_batch1_in_mega.txt > ${impbench}/washu_batch1_mega_keeplist_EUR.txt
awk '$5=="AMR" {print $3}' ${impbench}/washu_batch1_in_mega.txt > ${impbench}/washu_batch1_mega_keeplist_AMR.txt
awk '$5=="AFR" {print $3}' ${impbench}/washu_batch1_in_mega.txt > ${impbench}/washu_batch1_mega_keeplist_AFR.txt

#get keep lists (for WashU) by ancestry
awk '$5=="EUR" {print $1}' ${impbench}/washu_batch1_in_mega.txt > ${impbench}/washu_batch1_washu_keeplist_EUR.txt
awk '$5=="AMR" {print $1}' ${impbench}/washu_batch1_in_mega.txt > ${impbench}/washu_batch1_washu_keeplist_AMR.txt
awk '$5=="AFR" {print $1}' ${impbench}/washu_batch1_in_mega.txt > ${impbench}/washu_batch1_washu_keeplist_AFR.txt
