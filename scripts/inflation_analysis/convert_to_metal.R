

args = commandArgs(trailingOnly=TRUE)
filein = args[1]
fileout = args[2]

metalthem<-function(filein, fileout){
  d = read.table(filein, header=TRUE)
  d$RefAllele<-d$alleleA
  d$NonRefAllele<-d$alleleB
  d$BETA<-d$frequentist_add_beta_1.add.t1d.1
  d$SE<-d$frequentist_add_se_1
  d$SNP<-d$rsid
  dout <-d[,c("SNP","BETA","SE","NonRefAllele", "RefAllele")]
  write.table(dout,file=fileout, sep=" ", col.names=T, row.names=F, quote=F)
}

metalthem(filein,fileout)
