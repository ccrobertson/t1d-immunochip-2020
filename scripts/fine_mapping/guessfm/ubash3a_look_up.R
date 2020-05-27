#ubash3a_look_up.R
#updated the AIC plot as suggested by Linda, just include single SNP and 3 SNP model
#examining the UBASH3A region and it's complexity. Used to produce plots in Figure 1 in the manuscript.

library(snpStats)
library(annotSnpStats)
library(ggplot2)
extension="final_collection"
jamiedir<-paste0("/scratch/ji3x/",extension,"_gwas/META/EUR/")

load(file="/scratch/ji3x/final_collection_gwas/guessfm_5pc_diff/input/chr21.42416077.G.A/data.RData")


r<-read.table(file=paste0(jamiedir,"chr_21"),header=T,as.is=T)
r<-r[r$rsid %in% colnames(DATA),]
r$logp<-log10(r$frequentist_add_wald_pvalue_1)*-1

r$highlight<-ifelse(r$rsid=="chr21:42416077:G:A", "rs11203203",
ifelse(r$rsid=="chr21:42419803:T:C","rs7276555",
ifelse(r$rsid=="chr21:42418534:G:A","rs13048049",
ifelse(r$rsid=="chr21:42408836:T:C","rs9984852","other"))))
cols=c("grey","green","red","blue","black","yellow","brown","pink","turquoise","khaki1","slateblue1","deepskyblue1","linen","firebrick1")
cols<-cols[1:5]
names(cols)<-c("rs11203203","rs7276555","rs13048049","rs9984852","other")


g<-ggplot(data=r, aes(position, logp, colour=as.factor(highlight))) + geom_point(size=0.8) + 
theme(legend.position="none",
panel.background = element_blank(),
panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
axis.line = element_line(colour = "black"),
axis.text.y=element_text(size=4),
axis.title.y=element_text(size=5),
axis.text.x=element_text(size=4),
axis.title.x=element_text(size=5)) +
scale_colour_manual(values=cols) + scale_y_continuous(name=expression("-log"[10]*" p-value")) +
scale_x_continuous(name="Position on chromosome 21")
ggsave(g, file="/home/ji3x/output/final_collection_gwas/META/UBASH3A/init_manhattan.png",dpi=800, height=5,width=5,units="cm")


geno<-DATA[,c("chr21:42416077:G:A","chr21:42419803:T:C","chr21:42418534:G:A","chr21:42408836:T:C")]
cs<-col.summary(geno)
w<-which(cs$RAF>0.5)
if(length(w)>0)
geno<-switch.alleles(geno,snps=w)
geno<-as(geno,"numeric")

samps<-cbind(Y,covariates)
samps<-cbind(samps,geno)
colnames(samps)[7:10]<-c("rs11203203","rs7276555","rs13048049","rs9984852")


getaic<-function(mod){

aic=AIC(glm(data=samps, formula=as.formula(paste0("outcome ~ pc1 + pc2     + pc3 + pc4 + pc5",mod)),family="binomial"))
out<-data.frame(model=mod, aic=aic)
return(out)
}
mods<-c("","+rs7276555","+rs13048049","+rs9984852","+rs11203203","+rs9984852+rs13048049+rs7276555")

outs<-lapply(mods, getaic)
outs<-do.call("rbind",outs)
outs$model<-substr(outs$model,2,1000)
outs$ord<-c(1:nrow(outs))

g1<-ggplot(data=outs, aes(ord, aic)) + geom_point(size=0.8) +
#annotate("text",y=45195,x=1,label="no variants",colour="black",angle=90,size=1) +
#annotate("text",y=45195,x=2,label="rs7276555",colour="green",angle=90,size=1) +
#annotate("text",y=45195,x=3,label="rs13048049",colour="red",angle=90,size=1) +
#annotate("text",y=45195,x=4,label="rs9984852",colour="blue",angle=90,size=1) +
#annotate("text",y=45195,x=5,label="rs11203203",colour="grey",angle=90,size=1) +
#annotate("text",y=45171,x=6,label="rs11203203",colour="grey",angle=90,size=1) +
#annotate("text",y=45183,x=6,label="+",colour="black",angle=90,size=1) +
#annotate("text",y=45195,x=6,label="rs7276555",colour="green",angle=90,size=1) +
#annotate("text",y=45171,x=7,label="rs11203203",colour="grey",angle=90,size=1) +
#annotate("text",y=45183,x=7,label="+",colour="black",angle=90,size=1) +
#annotate("text",y=45195,x=7,label="rs13048049",colour="red",angle=90,size=1) +
#annotate("text",y=45171,x=8,label="rs9984852",colour="blue",angle=90,size=1) +
#annotate("text",y=45183,x=8,label="+",colour="black",angle=90,size=1) +
#annotate("text",y=45195,x=8,label="rs7276555",colour="green",angle=90,size=1) +
#annotate("text",y=45171,x=9,label="rs9984852",colour="blue",angle=90,size=1) +
#annotate("text",y=45183,x=9,label="+",colour="black",angle=90,size=1) +
#annotate("text",y=45195,x=9,label="rs13048049",colour="red",angle=90,size=1) +
#annotate("text",y=45147,x=10,label="rs11203203",colour="grey",angle=90,size=1) +
#annotate("text",y=45159,x=10,label="+",colour="black",angle=90,size=1) +
#annotate("text",y=45171,x=10,label="rs7276555",colour="green",angle=90,size=1) +
#annotate("text",y=45183,x=10,label="+",colour="black",angle=90,size=1) +
#annotate("text",y=45195,x=10,label="rs13048049",colour="red",angle=90,size=1) +
#annotate("text",y=45147,x=11,label="rs9984852",colour="blue",angle=90,size=1) +
#annotate("text",y=45159,x=11,label="+",colour="black",angle=90,size=1) +
#annotate("text",y=45171,x=11,label="rs7276555",colour="green",angle=90,size=1) +
#annotate("text",y=45183,x=11,label="+",colour="black",angle=90,size=1) +
#annotate("text",y=45195,x=11,label="rs13048049",colour="red",angle=90,size=1) +
theme(panel.background = element_blank(),
panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
axis.line = element_line(colour = "black"),axis.title.x=element_blank(),
axis.text.x=element_blank(),axis.ticks.x=element_blank(),
axis.text.y=element_text(size=3),
axis.title.y=element_text(size=5)) +
scale_y_continuous(name="Model AIC")

ggsave(g1, file="/home/ji3x/output/final_collection_gwas/META/UBASH3A/aic_comp_up.png",dpi=800, height=5,width=5,units="cm")
ggsave(g1, file="/home/ji3x/output/final_collection_gwas/META/UBASH3A/aic_comp_dpi100_up.png",dpi=100, height=5,width=5,units="cm")






g1<-ggplot(data=outs, aes(ord, aic)) + geom_point(size=0.7) +
theme(panel.background = element_blank(),
panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
axis.line = element_line(colour = "black"),axis.title.x=element_blank(),
axis.text.x=element_blank(),axis.ticks.x=element_blank(),
axis.text.y=element_text(size=4),
axis.title.y=element_text(size=5)) +
scale_y_continuous(name="Model AIC")

ggsave(g1, file="/home/ji3x/output/final_collection_gwas/META/UBASH3A/aic_comp_nolabels_up.png",dpi=800, height=5,width=5,units="cm")
