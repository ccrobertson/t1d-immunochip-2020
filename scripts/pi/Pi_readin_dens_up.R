#Pi_readin_dens_up.R

library(Pi)
tmpdir<-"/well/todd//users/jinshaw/mega/"
outdir<-"/well/todd/users/jinshaw/output/uva/Pi/"

load(file=paste0(tmpdir,"/Pi/eqtlgen_score_20_dens.RData"))

#show the networks for genes prioritised in the top 50 here that are not listed in top 150 in the Fang analysis:
pri<-dTarget$priority
tp<-pri[1:50,]


#compare to Fang results previously:
r<-read.table(file=paste0(tmpdir,"/Pi/fang_ranks.txt"),header=T,as.is=T)
r$name<-r$Gene
r<-r[,c(2,3)]

tp<-merge(tp,r, by="name",all.x=T)
tp$Rank<-ifelse(is.na(tp$Rank),151,tp$Rank)
tp$Rank<-as.character(tp$Rank)
tp$Rank<-ifelse(tp$Rank=="151",">150",tp$Rank)
tp<-tp[order(tp$rank),]

sink(file=paste0(outdir,"/top50_dens.txt"))
cat("Gene&rating&seed&nGene&cGene&eGene&dGene&pGene&fGene&fang_rank\\\\\n")
for(i in c(1:50)){
cat(paste0(tp[i,"name"],"&",tp[i,"rating"],"&",tp[i,"seed"],"&",tp[i,"nGene"],"&",
tp[i,"cGene"],"&",tp[i,"eGene"],"&",tp[i,"dGene"],"&",tp[i,"pGene"],"&",tp[i,"fGene"],"&",tp[i,"Rank"],"\\\\\n"))
}
sink()




xVisEvidence(dTarget,nodes="UBA52")
savePlot(filename=paste0(outdir,"/EBA52_network_dens"),type="png")
dev.off()
dev.off()

xVisEvidence(dTarget,nodes="IL2RA")
savePlot(filename=paste0(outdir,"/IL2RA_network_dens"),type="png")
dev.off()
dev.off()


