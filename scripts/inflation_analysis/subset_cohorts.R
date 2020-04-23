### Generate keep lists for subsetted case-control analyses

setwd(Sys.getenv("inflation"))

#get cluster assignments
clusters = read.table("../../pca/mega_cluster_assignments_and_pcs.txt", header=TRUE, sep="\t")
#get set of ids included in case control analyses
unrelateds = read.table("../cc_analysis/unrelateds_all.txt", header=TRUE)

keep = clusters[clusters$IID %in% unrelateds$iid & !is.na(clusters$PHENO), c("IID","cluster_label", "PHENO")]

getRandomSet = function(ancestry, n, seed) {

	set.seed(seed)
	source_cases = keep$IID[keep$cluster_label==ancestry & keep$PHENO==2]
	source_controls = keep$IID[keep$cluster_label==ancestry & keep$PHENO==1]
	return(c(sample(source_cases, size=n, replace=FALSE), sample(source_controls, size=n, replace=FALSE)))

}

getKeepLists = function(ancestry, n_cases, n_subs) {
	seeds = seq(1,n_subs)
	subs = list()
	for (i in 1:n_subs) {
		subs[[i]] = getRandomSet(ancestry, n=n_cases, seed=seeds[i])
		write.table(subs[[i]], file=paste0(ancestry,"_sub",i,".txt"), col.names=FALSE, row.names=FALSE, quote=FALSE)
	}
	subs
}


#get trio counts
trio_counts = read.table(paste0(Sys.getenv("fam_assoc"),"/trio_count.txt"), header=FALSE, skip=1)
names(trio_counts) = c("affected_trios","group")
trio_counts_dict = list()
for (i in 1:nrow(trio_counts)) {
	trio_counts_dict[[trio_counts$group[i]]] = trio_counts$affected_trios[i]
}

## EUR --> 4758 cases, 4758 controls
subs_EUR = getKeepLists("EUR",n_cases=trio_counts_dict[["EUR"]],n_subs=5)

## FIN --> 812 cases, 812 controls
subs_FIN = getKeepLists("FIN",n_cases=trio_counts_dict[["FIN"]],n_subs=5)

## AFR --> 122 cases, 122 controls
subs_AFR = getKeepLists("AFR",n_cases=trio_counts_dict[["AFR"]],n_subs=5)

## AMR --> 294 cases, 294 controls
subs_AMR = getKeepLists("AMR",n_cases=trio_counts_dict[["AMR"]],n_subs=5)
