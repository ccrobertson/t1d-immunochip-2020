#haplo_init.R
#haplotype analysis of the GUESSFM regions
library(snpStatsWriter)
library(mice)
library(plyr)
library(dplyr)
library(reshape)
library(ggplot2)
library(ggbio)
library(GUESSFM)

## NOTE must run: module load bcftools
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Set data input and output paths
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

jamie_dir = "/scratch/ccr5ju/MEGA/Jamies_files"
out_dir = "/scratch/ccr5ju/MEGA/halotypes_lookup"
cred_dir = "/scratch/ccr5ju/MEGA/release4/IMPUTED_TOPMED/fine_mapping/REDO"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Define function for running haplotype analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("/project/CPHG/T1D/MEGA/scripts/fine_mapping/haplotypes/haplotype_functions.R")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Run on CTLA4 region only
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
haplorun(index = "chr2:203874196:G:A", base=NULL, thr=0.004)

#haplotype analysis 1 =     rs231779, rs11681201, rs78847176, rs76562827,rs62184029, rs7369876 and rs117701653
taghap1 = c("rs231779", "rs11681201", "rs78847176", "rs76562827","rs62184029", "rs7369876", "rs117701653")
ss1 = haplorun_custom(index = "chr2:203874196:G:A", base="AATGCCC", thr=0.004, taghap=taghap1, hapnickname="taghap1")

#haplotype analysis 2 =  rs231779, rs231775, rs78847176, rs76562827, rs62184029, rs7369876 and rs117701653
taghap2 = c("rs231779", "rs231775", "rs78847176", "rs76562827", "rs62184029", "rs7369876", "rs117701653")
ss2 = haplorun_custom(index = "chr2:203874196:G:A", base=NULL, thr=0.001, taghap=taghap2, hapnickname="taghap2")

taghap2 = c("rs231779", "rs231775", "rs78847176", "rs76562827", "rs62184029", "rs7369876", "rs117701653")
ss2_collaped = haplorun_custom(index = "chr2:203874196:G:A", base=NULL, thr=0.004, taghap=taghap2, hapnickname="taghap2_collapsed")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Run haplotype analysis across T1D regions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#get significant regions
orig<-read.table(file=paste0(jamie_dir, "/pmeta_fam_ind_risk_regions_vcf_5pc_diff.txt"),header=T,sep=",",as.is=T)
orig<-orig[orig$ichip=="yes",]

#remove duplicate due to long LD block:
orig<-orig[!(orig$Marker %in% c("chr12:112153882:G:A","chr12:112730563:A:G","chr7:50406053:G:A",
"chr1:113302202:A:C","chr17:46751565:G:A","chr14:68286876:C:T","chr12:111066229:T:C",
"chr17:40623690:G:A")),]

#run for each index snp
#lapply(orig$Marker,haplorun, base=NULL, thr=0.004)
