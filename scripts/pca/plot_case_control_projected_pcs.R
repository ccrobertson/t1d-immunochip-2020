library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
dir = args[1]
group = args[2]
phenofile = args[3]
nickname = args[4]

setwd(dir)

eigenscorefile = paste0(nickname,"_pruned_unrelated_",group,"_pca_proj_onto_controls.sscore")
eigenvecvarfile = paste0(nickname,"_pruned_unrelated_",group,"_controls.eigenvec.var")
eigenvalfile = paste0(nickname,"_pruned_unrelated_",group,"_controls.eigenval")

mega_pheno = read.table(phenofile)
names(mega_pheno) <- c("cohort","FID","IID", "FAT","MOT","SEX","T1D")
row.names(mega_pheno) = mega_pheno$IID

pdf(paste("case_control_pca_plots_", group,".pdf", sep=""))

d = read.table(eigenscorefile, comment.char="~", header=TRUE)
d$T1D = factor(mega_pheno[d$IID, "T1D"], levels=c(1,2), labels=c("Unaffected", "Affected"))
d$PC1 = d$PC1_AVG
d$PC2 = d$PC2_AVG
d$cohort = mega_pheno[d$IID, "cohort"]

#plot by case control status
#ggplot() + geom_point(data=d, aes(x=PC1, y=PC2, color=T1D)) + ggtitle(group)

p1 = ggplot() +
		geom_point(data=d[d$T1D=="Unaffected",], aes(x=PC1, y=PC2, color=T1D), size=0.5) + ggtitle(group) +
		geom_point(data=d[d$T1D=="Affected",], aes(x=PC1, y=PC2, color=T1D), size=0.5) + ggtitle(group) +
		theme_bw() + theme(legend.position = "none")
print(p1)

p2 = ggplot() +
		geom_point(data=d[d$T1D=="Affected",], aes(x=PC1, y=PC2, color=T1D), size=0.5) + ggtitle(" ") +
		geom_point(data=d[d$T1D=="Unaffected",], aes(x=PC1, y=PC2, color=T1D),size=0.5) + ggtitle(" ") +
		theme_bw() + theme(legend.position = "none")
print(p2)

#plot by cohort
ggplot() + geom_point(data=d, aes(x=PC1, y=PC2, color=cohort)) + ggtitle(group)


#plot variance explained (scree plot)
varweights = read.table(eigenvecvarfile, comment.char="~", header=TRUE)
varsum = dim(varweights)[1]
eigenvals =  scan(eigenvalfile)
df = data.frame(pc=seq(1:10), var_explained=eigenvals/varsum)
p3 = ggplot() +
		geom_line(data=df, aes(x=pc, y=var_explained)) +
		geom_point(data=df, aes(x=pc, y=var_explained)) +
		scale_x_continuous(breaks=seq(1:10)) +
		xlab("Principal Component") +
		ylab("Variance explained") +
		theme_bw() + ggtitle(" ")
print(p3)

dev.off()

plots=list(p1,p2,p3)
save(plots, file=paste("case_control_pca_plots_", group,".RData", sep=""))
