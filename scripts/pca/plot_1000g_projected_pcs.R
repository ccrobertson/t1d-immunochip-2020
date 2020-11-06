options(stringsAsFactors=FALSE)
library(ggplot2)
library(gridExtra)
library(svglite)
library(ggpubr)

args = commandArgs(trailingOnly=TRUE)
cdir=args[1]
phenofile=args[2]
superpopulations=args[3]

setwd(cdir)

#get MEGA PCs
mega_all =  read.table(paste0("mega_cluster_assignments_and_pcs.txt"), header=TRUE, sep="\t")

#Get 1000G population codes
df.loc = read.table(superpopulations, header=FALSE)
names(df.loc) <- c("IID", "region", "continent")
#define AMR subregion variable
df.loc$amr_region = df.loc$region
#df.loc$amr_region[!df.loc$continent == "AMR"] <- df.loc$continent[!df.loc$continent == "AMR"]
df.loc$amr_region[!df.loc$continent == "AMR"] <- "Non-AMR"
row.names(df.loc) = df.loc$IID

#get 1000G PCs
thous_g = read.table("pca_proj_1000g.sscore")
names(thous_g) = c("FID", "IID","CNT1", "CNT2","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10" )
thous_g$continent = df.loc[thous_g$IID,"continent"]


#plot clusters (right = 1000g on top; left = mega on top)
getSideBySidePlot = function(cluster, title=NULL) {
	pointsize=0.3
	p1 = ggplot() +
			geom_point(data=mega_all[mega_all$cluster_label==cluster,], aes(x=PC1, y=PC2), size=pointsize) +
			geom_point(data=thous_g, aes(x=PC1, y=PC2, color=continent), size=pointsize) +
			theme_classic() + theme(legend.position="none")

	p2 = ggplot() +
			geom_point(data=thous_g, aes(x=PC1, y=PC2, color=continent), size=pointsize) +
			geom_point(data=mega_all[mega_all$cluster_label==cluster,], aes(x=PC1, y=PC2), size=pointsize) +
			theme_classic() + theme(legend.position="none")
	#pdf(paste0("supp_figure_",cluster,".pdf"), height=6, width=12)
	ppair = arrangeGrob(p1,p2, nrow=1, top=title)
	#, common.legend = TRUE, legend = "bottom"
	return(ppair)
	#dev.off()
}


thous_g$continent2 = factor(thous_g$continent, levels=c("AFR","AMR","EAS","EUR","SAS"), labels=c("1000G_AFR","1000G_AMR","1000G_EAS","1000G_EUR","1000G_SAS"))
p_legend = ggplot() +
		geom_point(data=thous_g, aes(x=PC1, y=PC2, color=continent2)) +
		theme(
				legend.text=element_text(size=12),
				legend.background = element_rect(colour = NA),
				legend.key=element_blank()) +
		labs(color="1000 Genomes\nContinental Populations")

png(paste0("ancestry_all.png"), height=480, width=640)
p_final = grid.arrange(getSideBySidePlot("AFR","African Admixed (AFR)"),
		getSideBySidePlot("AMR","Other Admixed (AMR)"),
		getSideBySidePlot("EAS","East Asian (EAS)"),
		getSideBySidePlot("FIN","Finnish (FIN)"),
		getSideBySidePlot("EUR","European (EUR)"),
		get_legend(p_legend),
		nrow=3)
p_final
dev.off()


#ggsave(file="ancestry_all.svg", plot=p_final, width=12, height=9)
#ggsave(file="ancestry_legend.svg", plot=p_legend, width=10, height=10)



##1000 Genomes only (no legend)
#p0a = ggplot() +
#		#geom_point(data=mega_all, aes(x=PC1, y=PC2)) +
#		geom_point(data=thous_g, aes(x=PC1, y=PC2, color=continent)) +
#		theme_classic() + theme(legend.position="none")
#pdf("supp_figure_1000g.pdf", height=6, width=6)
#print(p0a)
#dev.off()
#
##1000 Genomes with legend -- for use in inkscape
#p0b = ggplot() +
#		#geom_point(data=mega_all, aes(x=PC1, y=PC2)) +
#		geom_point(data=thous_g, aes(x=PC1, y=PC2, color=continent)) +
#		theme_classic()
#pdf("supp_figure_1000g_legend.pdf", height=6, width=6)
#print(p0b)
#dev.off()
#

#clusters=c("AMR","AFR","EAS","EUR","FIN")
#ppairs=list()
#for (i in 1:length(clusters)) {
#	cat(clusters[i],"\n")
#	ppairs[[i]] = getSideBySidePlot(clusters[i])
#}

#getSideBySidePlot = function(cluster) {
#	pointsize=0.5
#	p1 = ggplot() +
#			geom_point(data=mega_all[mega_all$cluster_label==cluster,], aes(x=PC1, y=PC2), size=pointsize) +
#			geom_point(data=thous_g, aes(x=PC1, y=PC2, color=continent), size=pointsize) +
#			theme_classic() + theme(legend.position="none")
#
#	p2 = ggplot() +
#			geom_point(data=thous_g, aes(x=PC1, y=PC2, color=continent), size=pointsize) +
#			geom_point(data=mega_all[mega_all$cluster_label==cluster,], aes(x=PC1, y=PC2), size=pointsize) +
#			theme_classic() + theme(legend.position="none")
#	list(p1, p2)
#}
#plotsByAncestry = list()
#for (i in 1:length(clusters)) {
#	cat(clusters[i],"\n")
#	plotsByAncestry[[i]] = getSideBySidePlot(clusters[i])
#}
#
#png("supp_figure_all.png", height=10, width=4)
#grid.arrange(plotsByAncestry[[1]][[1]], plotsByAncestry[[1]][[2]],
#		plotsByAncestry[[2]][[1]], plotsByAncestry[[2]][[2]],
#		plotsByAncestry[[3]][[1]], plotsByAncestry[[3]][[2]],
#		plotsByAncestry[[4]][[1]], plotsByAncestry[[4]][[2]],
#		plotsByAncestry[[5]][[1]], plotsByAncestry[[5]][[2]], nrow=5)
#
#dev.off()


#p2 = ggplot() +
#		geom_point(data=mega_all[mega_all$cluster_label=="AFR",], aes(x=PC1, y=PC2)) +
#		geom_point(data=thous_g, aes(x=PC1, y=PC2, color=continent)) +
#		theme_bw() + theme(legend.position="none")
#
#p3 = ggplot() +
#		geom_point(data=mega_all[mega_all$cluster_label=="EAS",], aes(x=PC1, y=PC2)) +
#		geom_point(data=thous_g, aes(x=PC1, y=PC2, color=continent)) +
#		theme_bw() + theme(legend.position="none")
#
#p4 = ggplot() +
#		geom_point(data=mega_all[mega_all$cluster_label=="EUR",], aes(x=PC1, y=PC2)) +
#		geom_point(data=thous_g, aes(x=PC1, y=PC2, color=continent)) +
#		theme_bw() + theme(legend.position="none")
#
#p5 = ggplot() +
#		geom_point(data=mega_all[mega_all$cluster_label=="FIN",], aes(x=PC1, y=PC2)) +
#		geom_point(data=thous_g, aes(x=PC1, y=PC2, color=continent)) +
#		theme_bw() + theme(legend.position="none")
#
#
#pdf("supp_figure.pdf", height=8, width=12)
#grid.arrange(p0,p1,p2,p3,p4,p5, nrow=2)
#dev.off()
#
#
##Plot of each cluster
#p0 = ggplot() +
#		geom_point(data=mega_all, aes(x=PC1, y=PC2, color=cluster_label)) +
#		ggtitle("All participants ") +
#		theme(legend.position="none",
#				axis.line = element_line(colour = "black"),
#				panel.grid.major = element_blank(),
#				panel.grid.minor = element_blank(),
#				panel.border = element_blank(),
#				panel.background = element_blank(),
#				axis.text = element_text(size=10),
#				axis.title = element_text(size=15),
#				title = element_text(size=20))
#
#p1 = ggplot() +
#		geom_point(data=mega_all, aes(x=PC1, y=PC2, color=cluster_label)) +
#		geom_point(data=mega_all[mega_all$cluster_label=="AMR",], aes(x=PC1, y=PC2)) +
#		ggtitle("Admixed American/\nSouth Asian") +
#		theme(legend.position="none",
#				axis.line = element_line(colour = "black"),
#				panel.grid.major = element_blank(),
#				panel.grid.minor = element_blank(),
#				panel.border = element_blank(),
#				panel.background = element_blank(),
#				axis.text = element_text(size=10),
#				axis.title = element_text(size=15),
#				title = element_text(size=20))
#
#p2 = ggplot() +
#		geom_point(data=mega_all, aes(x=PC1, y=PC2, color=cluster_label)) +
#		geom_point(data=mega_all[mega_all$cluster_label=="AFR",], aes(x=PC1, y=PC2)) +
#		ggtitle("African American") +
#		theme(legend.position="none",
#				axis.line = element_line(colour = "black"),
#				panel.grid.major = element_blank(),
#				panel.grid.minor = element_blank(),
#				panel.border = element_blank(),
#				panel.background = element_blank(),
#				axis.text = element_text(size=10),
#				axis.title = element_text(size=15),
#				title = element_text(size=20))
#
#p3 = ggplot() +
#		geom_point(data=mega_all, aes(x=PC1, y=PC2, color=cluster_label)) +
#		geom_point(data=mega_all[mega_all$cluster_label=="EAS",], aes(x=PC1, y=PC2)) +
#		ggtitle("East Asian") +
#		theme(legend.position="none",
#				axis.line = element_line(colour = "black"),
#				panel.grid.major = element_blank(),
#				panel.grid.minor = element_blank(),
#				panel.border = element_blank(),
#				panel.background = element_blank(),
#				axis.text = element_text(size=10),
#				axis.title = element_text(size=15),
#				title = element_text(size=20))
#
#p4 = ggplot() +
#		geom_point(data=mega_all, aes(x=PC1, y=PC2, color=cluster_label)) +
#		geom_point(data=mega_all[mega_all$cluster_label=="FIN",], aes(x=PC1, y=PC2)) +
#		ggtitle("Finnish") +
#		theme(legend.position="none",
#				axis.line = element_line(colour = "black"),
#				panel.grid.major = element_blank(),
#				panel.grid.minor = element_blank(),
#				panel.border = element_blank(),
#				panel.background = element_blank(),
#				axis.text = element_text(size=10),
#				axis.title = element_text(size=15),
#				title = element_text(size=20))
#
#p5 = ggplot() +
#		geom_point(data=mega_all, aes(x=PC1, y=PC2, color=cluster_label)) +
#		geom_point(data=mega_all[mega_all$cluster_label=="EUR",], aes(x=PC1, y=PC2)) +
#		ggtitle("European") +
#		theme(legend.position="none",
#				axis.line = element_line(colour = "black"),
#				panel.grid.major = element_blank(),
#				panel.grid.minor = element_blank(),
#				panel.border = element_blank(),
#				panel.background = element_blank(),
#				axis.text = element_text(size=10),
#				axis.title = element_text(size=15),
#				title = element_text(size=20))
#
#pdf('pca_by_cluster.pdf', width=12, height=7)
#grid.arrange(p0, p1,p3,p2,p4,p5, nrow=2)
#dev.off()
#
#
#
#
#
#
#
#
##full projection
#pdf("pca_all_onto_1000g.pdf")
#ggplot() +
#		geom_point(data=mega_all, aes(x=PC1, y=PC2)) +
#		geom_point(data=thous_g, aes(x=PC1, y=PC2, color=continent)) +
#		theme_bw()
#dev.off()
#
#
##Plots for each cohort
#pdf("pca_by_cohort.pdf")
#cohort_names = unique(mega_all$cohort)
#for (i in 1:length(cohort_names)) {
#	cohort=cohort_names[i]
#	p = ggplot() +
#			geom_point(data=thous_g, aes(x=PC1, y=PC2, color=continent)) +
#			geom_point(data=mega_all[mega_all$cohort==cohort,], aes(x=PC1, y=PC2)) +
#			ggtitle(cohort)
#	print(p)
#}
#dev.off()
#
##Plots for each analysis group
#pdf("pca_by_group.pdf")
#group_names = unique(mega_all$group)
#for (i in 1:length(group_names)) {
#	group=group_names[i]
#	p = ggplot() +
#			geom_point(data=thous_g, aes(x=PC1, y=PC2, color=continent)) +
#			geom_point(data=mega_all[mega_all$group==group,], aes(x=PC1, y=PC2)) +
#			ggtitle(group)
#	print(p)
#}
#dev.off()
#
#
#
#####EUROPEANS
#mega_eur =  read.table("pca_proj_mega_eur_onto_1000g_eur.sscore")
#names(mega_eur) = c("FID", "IID","PHENO", "CNT1", "CNT2", "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10" )
#mega_eur$cohort = mega_pheno[mega_eur$IID, "cohort"]
#mega_eur$group = mega_pheno[mega_eur$IID, "group"]
#
##get 1000G PCs
#thous_g_eur = read.table("pca_proj_1000g_eur.sscore")
#names(thous_g_eur) = c("FID", "IID","CNT1", "CNT2","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10" )
#thous_g_eur$region = df.loc[thous_g_eur$IID,"region"]
#
##full projection
#getClusterPlot2 = function(x, df) {
#	set.seed(15)
#	clusters = kmeans(df[,c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")], x)
#	df$clusters = as.factor(clusters$cluster)
#	p = ggplot() +
#			geom_point(data=thous_g_eur, aes(x=PC1, y=PC2)) +
#			geom_point(data=df, aes(x=PC1, y=PC2, color=clusters)) +
#			ggtitle("K-means clustering of MEGA subjects based on top 10 PCs") +
#			annotate("text", x=-0.1, y=0.15, label = paste("K =", x), size=10)
#	print(p)
#	table(df$clusters, df$group)
#}
#pdf("pca_clusters_EUR_only.pdf")
#ggplot() +
#		geom_point(data=mega_eur, aes(x=PC1, y=PC2)) +
#		geom_point(data=thous_g_eur, aes(x=PC1, y=PC2, color=region))
#
#ggplot() + geom_point(data=thous_g_eur, aes(x=PC1, y=PC2, color=region))
#getClusterPlot2(2, mega_eur)
#getClusterPlot2(3, mega_eur)
#getClusterPlot2(4, mega_eur)
#getClusterPlot2(5, mega_eur)
#dev.off()
#
#set.seed(15)
#eur_clusters2 = kmeans(mega_eur[,c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")], 2)
#table(eur_clusters2$cluster)
#output_eur_clusters = data.frame(mega_eur[,c("FID","IID")], cluster=eur_clusters2$cluster)
#output_eur_clusters$group = NA
#output_eur_clusters$group[output_eur_clusters$cluster==1] <- "SOUTHERN"
#output_eur_clusters$group[output_eur_clusters$cluster==2] <- "NORTHERN"
#write.table(output_eur_clusters, file="mega_eur_cluster_assignments.txt", quote=FALSE, row.names=FALSE, col.names=TRUE)
#
#pdf("pca_south_vs_north_eur.pdf")
#southern_europeans = output_eur_clusters[output_eur_clusters$group=="SOUTHERN", "IID"]
#ggplot() +
#		geom_point(data=thous_g, aes(x=PC1, y=PC2, color=continent)) +
#		geom_point(data=mega_all[mega_all$IID%in%southern_europeans,], aes(x=PC1, y=PC2))
#
#northern_europeans = output_eur_clusters[output_eur_clusters$group=="NORTHERN", "IID"]
#ggplot() +
#		geom_point(data=thous_g, aes(x=PC1, y=PC2, color=continent)) +
#		geom_point(data=mega_all[mega_all$IID%in%northern_europeans,], aes(x=PC1, y=PC2))
#dev.off()
#
#
#
#
#
#
#### HISPANICS
##pca_proj_hispanics.sscore
##pca_proj_1000g_onto_hispanics.sscore
##get MEGA self-reported hispanic status
#hisp = scan("hispanic_subjects.txt")
#
#
#
##cd /Users/Cassie/Box\ Sync/Rich\ Lab/MEGA/ancestry/
##scp ccr5ju@dobby.cphg.virginia.edu:/data1/ccr5ju/MEGA/pca/pca_*.pdf .
#
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
#
###readinpcas_allinds.R
##library(ggplot2)
##
##
###examining the proportion of individuals in our netire dataset that look like they may be of African ancestry:
##labs = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "P10")
##
###get location data:
##loc<-read.table(file="/data1/ji3x/1000genomes/1000genomes_locations.csv", sep=",", as.is=T)
##
##
###have a look at the PCs from the ic65k cohort:
##vecs1<-read.table(file="/data1/ji3x/qc/SNPQC/imputed_gwas/pca_proj_ic65k_allinds.sscore", header=F, as.is=T)
##colnames(vecs1)<-c("FID","IID","PHENO1","NMISS_ALLELE_CT","NAMED_ALLELE_DOSAGE_SUM",labs)
##
##vecs1$country="Projected"
##vecs1$pop="Projected"
##
###merge with info file to establish membership of which group:
##samp<-read.table(file="/data1/ji3x/processed/megaic_with_cohort.fam", header=T, as.is=T)
##vecs1$uniqueID<-paste0(vecs1$FID,".",vecs1$IID)
##
###### VISUALIZE ANALYSIS GROUPS
##ggplot(data=vecsboth, aes(PC2, PC1, color=group)) + geom_point() +
##		scale_x_continuous(name="PC2") + scale_y_continuous(name="PC1") +
##		geom_hline(yintercept=0) +
##		geom_vline(xintercept=-0.04) +
##		geom_vline(xintercept=0.05)
##
##ggplot(data=vecsboth, aes(PC2, PC1, color=group)) + geom_point() +
##		scale_x_continuous(name="PC2") + scale_y_continuous(name="PC1") +
##		facet_grid(.~group)
##
##
###### VISUALIZE COHORTS
##vecsboth<-merge(vecs1, samp, by="uniqueID",all.x=T)
##vecsboth$group<-ifelse(vecsboth$Cohort %in% c("58BC","B58C","BRI","BDA","GRID","Sanger_BC58","T1DGC-WARREN","UK","UKBS"), "GB",
##		ifelse(vecsboth$Cohort %in% c("DAN-SDC"), "DAN",
##				ifelse(vecsboth$Cohort %in% c("NIMH", "Trialnet", "UCSF_AA", "BRI","EDIC","GoKinD","NYCP","SEARCH") | is.na(vecsboth$Cohort), "USA",
##						ifelse(vecsboth$Cohort %in% c("EUR", "Existing_Source_to_T1DGC:SWEDISH"), "EUR",
##								ifelse(vecsboth$Cohort %in% c("MILWAUKEE", "IDDMGEN"), "FIN",
##										ifelse(vecsboth$Cohort %in% c("T1DGC-NI", "T1DGC-YH"),"NI",
##												ifelse(vecsboth$Cohort %in% c("AP"), "AP", NA)))))))
##
##ggplot(data=vecsboth, aes(PC2, PC1, color=Cohort)) + geom_point() +
##		scale_x_continuous(name="PC2") + scale_y_continuous(name="PC1") +
##		facet_grid(.~group)
##
##
###### VISUALIZE RACE
###create self-reported hispanic variable
##t1dgc_phenos = read.csv("/m/jdrfdn_scratch/users/projects/IMCHIP/T1DGC.2011.03_Golden_Pedigree/T1DGC.2011.03_Golden_Pedigree/Resources/T1DGC.2011.03_Resources.csv")
##hispanic_subjects = t1dgc_phenos$Analytic.ID[t1dgc_phenos$Hispanic==1]
##vecsboth$hispanic = NA
#vecsboth$hispanic[vecsboth$IID%in%hispanic_subjects] <- 1
#vecsboth$hispanic[!vecsboth$IID%in%hispanic_subjects] <- 0
#
##plot hispanics
#ggplot(data=vecsboth) +
#		geom_point(aes(PC2, PC1, color=as.factor(hispanic))) +
#		geom_point(data = vecsboth[vecsboth$hispanic==1,], aes(PC2, PC1, color=as.factor(hispanic))) +
#		scale_x_continuous(name="PC2") + scale_y_continuous(name="PC1")
#
##hclustOut = hclust(dist(vecsboth[,labs]))
#
#
#
##from this projection, it looks like anyone with PC1<0 is more similar to African than any other ancestry:
#nottg<-vecs1[vecs1$country=="Projected",]
#africans<-nottg[nottg$PC1<0,]
#message(paste0(nrow(africans), " Africans"))
#
##south asians going to take as PC1>0 & -0.04<PC2<0.05
#southasians<-nottg[nottg$PC1>=0 & nottg$PC2>-0.04 & nottg$PC2<=0.05,]
#message(paste0(nrow(southasians), " South Asians"))
#
##east asians going to take as PC1>0 & PC2>0.05
#eastasians<-nottg[nottg$PC1>=0 & nottg$PC2>0.05,]
#message(paste0(nrow(eastasians), " East Asians"))
#
#
##europeans going to take as PC1>0 & PC2<=-0.04
#europeans<-nottg[nottg$PC1>0 & nottg$PC2<=-0.04,]
#message(paste0(nrow(europeans), " Europeans"))
#
#
#


#
