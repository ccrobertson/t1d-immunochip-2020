##get atac data
#sample_ids = scan("../Calderon_et_al_ATAC/cell_types.txt", what="character")
#getBedGraph = function(sample) {
#	file = paste("../Calderon_et_al_ATAC/dds_ds", sample, "v1.bedGraph", sep=".")
#	cat(file,"\n")
#	type=unlist(strsplit(sample, split=".", fixed=TRUE))[1]
#	stim=unlist(strsplit(sample, split=".", fixed=TRUE))[2]
#	data=read.table(file, header=FALSE, skip=1)
#	names(data) <- c("chr","start","end",sample)
#	list(type=type, stim=stim, data=data)
#}
#bedGraphs = lapply(sample_ids, getBedGraph)
#names(bedGraphs) <- sample_ids
#bedGraphsCombinedList = lapply(bedGraphs, `[[`, 3)
#bedGraphsCombinedDF = Reduce(function(x, y) merge(x, y, by = c("chr","start","end"), all = TRUE), bedGraphsCombinedList)
#cell_types = unlist(lapply(bedGraphs, `[[`, 1))
#stim_groups = unlist(lapply(bedGraphs, `[[`, 2))
#
#write.table(bedGraphsCombinedDF, file="../Calderon_et_al_ATAC/combined_bedGraph.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
#save(sample_ids, cell_types, stim_groups, bedGraphs, bedGraphsCombinedDF, file="../Calderon_et_al_ATAC/combined_bedGraph.RData")
#load(file="../Calderon_et_al_ATAC/combined_bedGraph.RData")


#load newest data
setwd("/nv/vol185/MEGA/atac_data/20190127/bedGraphs")
config = read.table("sample_config.txt", header=TRUE)
bedGraphs = list()
sample_ids = NULL
for (i in 1:nrow(config)) {
	group = config[i, "group"]
	sample_id = config[i, "sample_id"]
	file_path = config[i, "file_path"]
	cat(file_path,"\n")
	data=read.table(file_path, header=FALSE, skip=1)
	names(data) <- c("chr","start","end",sample_id)
	sample_ids[i] = sample_id
	bedGraphs[[i]] = data
}
names(bedGraphs) <- sample_ids
bedGraphsCombinedDF = Reduce(function(x, y) merge(x, y, by = c("chr","start","end"), all = TRUE), bedGraphs)
save(sample_ids, bedGraphs, bedGraphsCombinedDF, file="combined_bedGraph.RData")