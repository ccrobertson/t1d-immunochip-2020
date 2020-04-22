##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
options(stringsAsFactors=FALSE)
library(ggplot2)
library(gridExtra)


args = commandArgs(trailingOnly=TRUE)

dir=args[1]
phenofile=args[2]
superpopulations=args[3]
project=args[4]

setwd(dir)

#Get MEGA phenofile
mega_pheno = read.table(phenofile)
names(mega_pheno) <- c("cohort","FID","IID", "FAT","MOT","SEX","T1D")
row.names(mega_pheno) = mega_pheno$IID


#Get 1000G population codes
df.loc = read.table(superpopulations, header=FALSE)
names(df.loc) <- c("IID", "region", "continent")
#define AMR subregion variable
df.loc$amr_region = df.loc$region
#df.loc$amr_region[!df.loc$continent == "AMR"] <- df.loc$continent[!df.loc$continent == "AMR"]
df.loc$amr_region[!df.loc$continent == "AMR"] <- "Non-AMR"
row.names(df.loc) = df.loc$IID


#Define MEGA analysis "Groups"
mega_pheno$group = "unassigned"
mega_pheno$group[mega_pheno$cohort%in% c("BC58","BRI","BDA","GRID","Sanger_BC58","T1DGC-WARREN","T1DGC-UK","UKBS","cbr_controls")] <- "GB"
mega_pheno$group[mega_pheno$cohort%in% c("DAN-SDC","IC_Cases_Steno","IC_Cases_HSG")] <- "DAN"
mega_pheno$group[mega_pheno$cohort%in% c("NIMH", "Trialnet", "UCSF_AA", "BRI","EDIC","GoKinD","NYCP","SEARCH","T1DGC-NA","UAB","CLEAR","UCHSC2")] <- "USA"
mega_pheno$group[mega_pheno$cohort%in% c("T1DGC-EUR")] <- "EUR"
mega_pheno$group[mega_pheno$cohort%in% c("MILWAUKEE", "IDDMGEN","T1DGEN")] <- "FIN"
mega_pheno$group[mega_pheno$cohort%in% c("T1DGC-NI", "T1DGC-YH")] <- "NI"
mega_pheno$group[mega_pheno$cohort%in% c("T1DGC-AP")] <- "AP"
mega_pheno$group[mega_pheno$cohort%in% c("T1DGC")] <- "Other T1DGC"
mega_pheno$group[mega_pheno$cohort%in% c("HapMap")] <- "HapMap"
sum(is.na(mega_pheno$group))

### ALL SUBJECTS

#get MEGA PCs
mega_all =  read.table(paste0("pca_proj_",project,"_all_onto_1000g.sscore"))
names(mega_all) = c("FID", "IID","PHENO", "CNT1", "CNT2", "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10" )
mega_all$cohort = mega_pheno[mega_all$IID, "cohort"]
mega_all$group = mega_pheno[mega_all$IID, "group"]

#get 1000G PCs
thous_g = read.table("pca_proj_1000g.sscore")
names(thous_g) = c("FID", "IID","CNT1", "CNT2","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10" )
thous_g$continent = df.loc[thous_g$IID,"continent"]



#Define groups for imputation
set.seed(43)
clusters5 = kmeans(mega_all[,c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")], 5)
mega_all$clusters5 = as.factor(clusters5$cluster)
tbl_c5 = table(mega_all$clusters5, mega_all$group)
#tbl_c5 = table(mega_all$clusters5, mega_all$cohort)
mega_all$cluster_label = NA
mega_all$cluster_label[mega_all$clusters5==order(tbl_c5[,"USA"], decreasing=TRUE)[3]] <- "AMR"
mega_all$cluster_label[mega_all$clusters5==order(tbl_c5[,"FIN"], decreasing=TRUE)[1]] <- "FIN"
mega_all$cluster_label[mega_all$clusters5==order(tbl_c5[,"GB"], decreasing=TRUE)[1]] <- "EUR"
mega_all$cluster_label[mega_all$clusters5==order(tbl_c5[,"USA"], decreasing=TRUE)[2]] <- "AFR"
mega_all$cluster_label[mega_all$clusters5==order(tbl_c5[,"AP"], decreasing=TRUE)[3]] <- "EAS"

write.table(mega_all, file=paste0(project,"_cluster_assignments_and_pcs.txt"),quote=FALSE,row.names=FALSE, col.names=TRUE, sep="\t")
write.table(mega_all[,c("FID","IID","clusters5","cluster_label")], file=paste0(project,"_cluster_assignments.txt"),quote=FALSE,row.names=FALSE, col.names=TRUE)



#Plots of kmeans clustering
getClusterPlot = function(x) {
	set.seed(43)
	clusters = kmeans(mega_all[,c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")], x)
	mega_all$clusters = as.factor(clusters$cluster)
	p = ggplot() +
			geom_point(data=mega_all, aes(x=PC1, y=PC2, color=clusters)) +
			ggtitle("K-means clustering of MEGA subjects based on top 10 PCs") +
			annotate("text", x=-0.1, y=0.15, label = paste("K =", x), size=10)
	print(p)
	table(mega_all$clusters, mega_all$group)
}

for (x in 3:8) {
	png(paste0("pca_clustering_k",x,".png"))
	getClusterPlot(x)
	dev.off()
}
