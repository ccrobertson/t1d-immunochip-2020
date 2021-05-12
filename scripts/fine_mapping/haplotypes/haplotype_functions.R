#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Define function for reading SNPHAP output
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hap.read <- function(f) {
  d <- read.table(f,sep="\t",as.is=TRUE,header=TRUE)
  d <- d[,-c(2,3)]
  d$hap <- apply(d[,-c(1,ncol(d)),drop=FALSE],1,paste,collapse="")
  d
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Define function for multiple imputation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

model.mi <- function(dir,df,thr=0.004,maxi=NULL,haps.pattern, covars=c("pc1","pc2","pc3","pc4","pc5"), phenotype="outcome",base=NULL, family="binomial") {

  # dir=out_subdir
  # df=samples
  # thr=0.004
  # maxi=NULL
  # haps.pattern="init.out2"
  # covars=c("pc1","pc2","pc3","pc4","pc5")
  # phenotype="outcome"
  # base=NULL
  # family="binomial"

  ## fix haplotype base and set of haps to use
  haps <- hap.read(file.path(dir,haps.pattern))
  tt <- tapply(haps$Prob,haps$hap,sum)
  tt <- 100*tt/sum(tt)
  use <- names(tt)[tt/sum(tt)>thr]
  if(!is.null(base)) {
    if(!(base %in% use)) {
      message("haplotypes found:")
      cat(use,sep="\n")
      stop("requested base haplotype not found: ", base)
    }
  } else {
    base <- names(tt)[which.max(tt)]
  }
  lev <- names(sort(tt[use],decreasing=TRUE))
  message(length(lev)," haplotypes found meeting frequency threshold ",thr)
  lev <- c(lev,"OTHER")
  wh <- which(lev==base)
  if(wh!=1)
    lev <- c(lev[wh],lev[-wh])

  ## model formula (adjusting for the population structure variables)
  f <- "y ~ hap + pc1 + pc2 + pc3 + pc4 + pc5"

  ## read MI and fit models
  files <- list.files(dir,pattern=paste0(haps.pattern,"."),full=TRUE)
  if(!is.null(maxi) && length(files)>maxi)
    files <- files[1:maxi]
  message("reading ",length(files)," multiply imputed datasets from ",dir)
  fits <- lapply(seq_along(files), function(i) {
    cat(".")
    haps <- hap.read(files[[i]])
    haps <- haps[ haps$id %in% rownames(df), ]
    m <- match(haps$id, rownames(df))
    if(length(covars))
      haps <- cbind(haps, df[m,covars,drop=FALSE])

    haps$y <- df[ m, phenotype ]
    ## if(length(kk<-setdiff(use,haps$hap)))
    ##   stop(length(kk)," haps not found")
    haps$hap[ !(haps$hap %in% lev) ] <- "OTHER"
    ##    haps <- subset(haps, hap %in% use)
    haps$hap <- factor(haps$hap, levels=lev)
    haps<-haps[!is.na(haps$y),]
    glm(data=haps, as.formula(f),family=family)
  })
  ## make mira
  ## print(table(sapply(lapply(fits,coefficients),length)))
  message("averaging results over the MI fits")
  mi <- as.mira(fits)
  ss <- as.data.frame(summary(pool(mi)))
  null.value=0

  rownames(ss) <- ss$term
  rownames(ss)[1] <- paste0("hap",base)

  ss$Fq <- tt[sub("hap","",rownames(ss))]
  ss$estimate[1] <- NA
  ss$std.error[1] <- NA
  ss$statistic[1] <- NA
  ss$p.value[1]<-NA
  ss$p <- ss$p.value
  ss$beta <- ss$estimate
  ss$se<-ss$std.error
  ss <- ss[,c("beta","se","p","Fq")]
  return(ss)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Define function for plotting results -- averaging across imputations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plotHaplotypes = function(ss, thr, l, cols) {

  hsnps<-rownames(l)
  hsnps1<-l$ID
  hsnps1<-ifelse(hsnps1==".",hsnps,hsnps1)
  l$dbsnp<-l$cols
  l$allele.A<-l$a0
  l$allele.B<-l$a1
  l$dbsnp=l$MarkerName

  results<-ss
  results <- subset(results,Fq>(thr*100))
  results$ref<-ifelse(rownames(results)==rownames(ss)[1],1,0)
  results<-results[order(-results$ref,-results$Fq),]
  results$ord<-c(1:nrow(results))
  results$out="t1d"
  results$hap<-reorder(results$hap, results$ord)
  #results$hap<-results$ord

  hlab <- results[!duplicated(results[,"hap"]),]
  hlab <- within(hlab, {
    x <- as.numeric(hlab$hap) + 0.25
    hap <- sub("hap|dip","",hap)
  }) %>%
    subset(., hap!="OTHER" & Fq>thr*100, select=c("hap","x"))
  hlab <- cbind(hlab, strsplit(hlab$hap,"") %>% do.call("rbind",.))

  h <- melt(hlab,c("hap","x"))
  h$snp <- hsnps[as.numeric(h$variable)]
  h$dbsnp<-h$snp
  l<-l[,!colnames(l) %in% "snp"]
  h<-join(h,l, by="dbsnp")
  h$major<-(as.character(h$Allele1)==as.character(h$value))
  h$y <- as.numeric(h$variable)

  s<-data.frame(dbsnp=hsnps)
  g<-data.frame(group=h$variable, dbsnp=h$dbsnp)
  s<-merge(s, g, by="dbsnp", all.x=TRUE)
  s<-unique(s[,1:2])
  rownames(s)<-s$dbsnp
  s<-s[l$dbsnp,]

  #plot tiles showing haplotypes
  hplot <- ggplot(h,aes(x=x,y=y,label=value,colour=major,xmin=x-0.5,xmax=x+0.5,ymin=y-0.5,ymax=y+0.5,fill=major)) +
    geom_rect() +
    geom_text(size=5) +
    theme_minimal() +
    scale_y_continuous(breaks=seq_along(hsnps1),labels=gsub("\\..*","",hsnps1),limits=c(0.5,length(hsnps1) + 0.5)) +
    scale_fill_manual("",values=c("white","grey20"), labels=c("Minor allele", "Major allele")) +
    scale_colour_manual("",values=c("grey20","white"), labels=c("Minor allele", "Major allele")) + ylab("SNP") +
    scale_x_continuous("Haplotype",breaks=as.numeric(hlab$hap)+0.25,
                       labels=sub("hap","",hlab$hap)) +
    theme(legend.position="none", legend.key.size=unit(0.5,"cm"), title=element_text(size=10)) +
    theme(axis.text=element_text(size=12),
          axis.title.y=element_text(angle=90, vjust=0.5),
    axis.text.y=element_text(colour=cols[hsnps1]))

  #plot haplotype effect size and confidence interval
  #calculate CIs
  m<-results[!duplicated(results[,"hap"]),]
  m$ymin <- m$beta - 1.96*m$se
  m$ymax <- m$beta + 1.96*m$se
  vlines <- unique(as.numeric(m$hap)) - 0.25
  oplot <- ggplot(m, aes(x=(as.numeric(hap)+0.25),y=beta,ymin=ymin,ymax=ymax)) + geom_pointrange(position = position_dodge(width=0.75)) +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=vlines,colour="grey") +
    labs(colour="Outcome") +
    scale_x_continuous("Haplotype",breaks=as.numeric(hlab$hap)+0.25,
                       labels=sub("hap","",hlab$hap)) +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=90,size=10,family="monospace"),
          legend.position="top", axis.title=element_text(size=10),
          axis.title.y=element_text(angle=90, vjust=0.5)) +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major = element_blank()) +
    ylab("log-odds ratio (95% CI)")

  #plot haplotype frequencies
  labs<-c(round(m$Fq/100,2))
  getxs<-function(out){
    h<-m[m$out==out,]
    h$xs<-c(1:nrow(h))
    if (out=="t1d"){
    h$xs<-h$xs+0.5
    }
    if (out=="controls"){
    h$xs<-h$xs+0.25
    }
    if (out=="cases"){
    h$xs<-h$xs
    }
    return(h)
  }
  m<-lapply(unique(m$out),getxs)
  m<-do.call("rbind",m)
  labs<-c(round(m$Fq/100,3))
  xs<-c(1:nrow(m))
  xs<-xs+0.25
  fplot<- ggplot(m, aes(x=as.numeric(hap), y=0)) +
    annotate("text", label=labs, x=xs, y=0, size=3, angle=90) +
    scale_y_continuous(name="Frequency", breaks=c()) + theme_minimal() +
    theme(axis.title.y=element_text(angle=90, size=10, vjust=0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_vline(xintercept=vlines,colour="grey")

  #combine plots
  t<-tracks(oplot, fplot,  hplot, heights=c(12,6,30))
  return(t)
}


plotHaplotypes_nocolor = function(ss, thr, l) {

  hsnps<-rownames(l)
  hsnps1<-l$ID
  hsnps1<-ifelse(hsnps1==".",hsnps,hsnps1)
  l$dbsnp<-l$cols
  l$allele.A<-l$a0
  l$allele.B<-l$a1
  l$dbsnp=l$MarkerName

  results<-ss
  results <- subset(results,Fq>(thr*100))
  results$ref<-ifelse(rownames(results)==rownames(ss)[1],1,0)
  results<-results[order(-results$ref,-results$Fq),]
  results$ord<-c(1:nrow(results))
  results$out="t1d"
  results$hap<-reorder(results$hap, results$ord)
  #results$hap<-results$ord

  hlab <- results[!duplicated(results[,"hap"]),]
  hlab <- within(hlab, {
    x <- as.numeric(hlab$hap) + 0.25
    hap <- sub("hap|dip","",hap)
  }) %>%
    subset(., hap!="OTHER" & Fq>thr*100, select=c("hap","x"))
  hlab <- cbind(hlab, strsplit(hlab$hap,"") %>% do.call("rbind",.))

  h <- melt(hlab,c("hap","x"))
  h$snp <- hsnps[as.numeric(h$variable)]
  h$dbsnp<-h$snp
  l<-l[,!colnames(l) %in% "snp"]
  h<-join(h,l, by="dbsnp")
  h$major<-(as.character(h$Allele1)==as.character(h$value))
  h$y <- as.numeric(h$variable)

  s<-data.frame(dbsnp=hsnps)
  g<-data.frame(group=h$variable, dbsnp=h$dbsnp)
  s<-merge(s, g, by="dbsnp", all.x=TRUE)
  s<-unique(s[,1:2])
  rownames(s)<-s$dbsnp
  s<-s[l$dbsnp,]

  #plot tiles showing haplotypes
  hplot <- ggplot(h,aes(x=x,y=y,label=value,colour=major,xmin=x-0.5,xmax=x+0.5,ymin=y-0.5,ymax=y+0.5,fill=major)) +
    geom_rect() +
    geom_text(size=5) +
    theme_minimal() +
    scale_y_continuous(breaks=seq_along(hsnps1),labels=gsub("\\..*","",hsnps1),limits=c(0.5,length(hsnps1) + 0.5)) +
    scale_fill_manual("",values=c("white","grey20"), labels=c("Minor allele", "Major allele")) +
    scale_colour_manual("",values=c("grey20","white"), labels=c("Minor allele", "Major allele")) + ylab("SNP") +
    scale_x_continuous("Haplotype",breaks=as.numeric(hlab$hap)+0.25,
                       labels=sub("hap","",hlab$hap)) +
    theme(legend.position="none", legend.key.size=unit(0.5,"cm"), title=element_text(size=10)) +
    theme(axis.text=element_text(size=12),
          axis.title.y=element_text(angle=90, vjust=0.5))
    #axis.text.y=element_text(colour=cols[hsnps1]))

  #plot haplotype effect size and confidence interval
  #calculate CIs
  m<-results[!duplicated(results[,"hap"]),]
  m$ymin <- m$beta - 1.96*m$se
  m$ymax <- m$beta + 1.96*m$se
  vlines <- unique(as.numeric(m$hap)) - 0.25
  oplot <- ggplot(m, aes(x=(as.numeric(hap)+0.25),y=beta,ymin=ymin,ymax=ymax)) + geom_pointrange(position = position_dodge(width=0.75)) +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=vlines,colour="grey") +
    labs(colour="Outcome") +
    scale_x_continuous("Haplotype",breaks=as.numeric(hlab$hap)+0.25,
                       labels=sub("hap","",hlab$hap)) +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=90,size=10,family="monospace"),
          legend.position="top", axis.title=element_text(size=10),
          axis.title.y=element_text(angle=90, vjust=0.5)) +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major = element_blank()) +
    ylab("log-odds ratio (95% CI)")

  #plot haplotype frequencies
  labs<-c(round(m$Fq/100,2))
  getxs<-function(out){
    h<-m[m$out==out,]
    h$xs<-c(1:nrow(h))
    if (out=="t1d"){
    h$xs<-h$xs+0.5
    }
    if (out=="controls"){
    h$xs<-h$xs+0.25
    }
    if (out=="cases"){
    h$xs<-h$xs
    }
    return(h)
  }
  m<-lapply(unique(m$out),getxs)
  m<-do.call("rbind",m)
  labs<-c(round(m$Fq/100,3))
  xs<-c(1:nrow(m))
  xs<-xs+0.25
  fplot<- ggplot(m, aes(x=as.numeric(hap), y=0)) +
    annotate("text", label=labs, x=xs, y=0, size=3, angle=90) +
    scale_y_continuous(name="Frequency", breaks=c()) + theme_minimal() +
    theme(axis.title.y=element_text(angle=90, size=10, vjust=0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_vline(xintercept=vlines,colour="grey")

  #combine plots
  t<-tracks(oplot, fplot,  hplot, heights=c(12,6,30))
  return(t)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Master function to run it all
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

haplorun<-function(index, base, thr){

  #index = "chr2:203874196:G:A"
  #base = NULL
  #thr=0.004
  out_subdir = file.path(out_dir,gsub(":",".",index),"guessfmtags")
  dir.create(out_subdir)

  #load genotype matrix and phenotype data (objects: Y, covariates, DATA)
  d <-paste0(jamie_dir,"/guessfm_5pc_diff/input/",gsub(":",".",index),"/")
  load(file=paste0(d,"data.RData"), verbose=TRUE)

  #load guessfm output (objects: groups, l, summx, sm.all, sp.all)
  load(file=paste0(d,"summx.RData"), verbose=TRUE)
  o<-summx[,c("chromosome","position")]
  write.table(o, file=paste0(out_subdir, "/sig_vars.txt"), sep="\t", col.names=F, row.names=F, quote=F)

  #extract variant ids (rsid) for snps in sp.all using Kaviar reference
  system(paste0("/apps/software/standard/core/bcftools/1.9/bin/bcftools view -R ",out_subdir,"/sig_vars.txt ", jamie_dir,
  "/Kaviar-160204-Public-hg38-trim.vcf.gz > ",out_subdir, "/sig_vars_ids.txt"))

  #add variant ids to guessfm summx object
  vars<-read.table(file=paste0(out_subdir, "/sig_vars_ids.txt"),skip=40, header=T, comment.char="", as.is=T,colClasses = c("character"))
  vars$MarkerName<-paste0("chr",vars$X.CHROM,":",vars$POS,":",vars$REF,":",vars$ALT)
  vars<-vars[,c("MarkerName","ID")]
  vars<-vars[vars$ID!=".",]
  summx<-merge(summx,vars,by="MarkerName", all.x=T)

  #sort guessfm summx object
  summx<-summx[order(summx$chromosome, summx$position),]
  rownames(summx)<-summx$snp
  summx$ID<-ifelse(summx$snp=="chr2.162279995.C.G","rs35337543",
  ifelse(summx$snp=="chr12.56042145.C.G","rs1131017",
  ifelse(summx$snp=="chr14.68796357.G.A","rs6573857",
  ifelse(summx$snp=="chr14.98011018.A.T","rs79835824",
  ifelse(summx$snp=="chr16.11190252.A.G","rs248832",
  ifelse(summx$snp=="chr16.28494339.G.C","rs151234",
  ifelse(summx$snp=="chr21.44280294.T.A","rs9671675",
  ifelse(summx$snp=="chr21.42418534.G.A","rs13048049",summx$ID))))))))
  summx$ID<-ifelse(is.na(summx$ID),".",summx$ID)

  #load snps to use in haplotype analysis:
  mydir<-paste0(cred_dir,"/",gsub(":",".",index),"/")
  creds<-read.table(paste0(mydir,"credible_snps.txt"),as.is=T, header=T)
  if(index!="chr10:6052734:C:T"){
  creds<-creds[creds$ppsum>0.5,]
  }
  #il2ra v complicated - just keeping those with ppsum>0.8:
  if(index=="chr10:6052734:C:T"){
  creds<-creds[creds$ppsum>0.8,]
  }
  creds$ID<-ifelse(creds$snp=="chr2.162279995.C.G","rs35337543",
  ifelse(creds$snp=="chr12.56042145.C.G","rs1131017",
  ifelse(creds$snp=="chr14.68796357.G.A","rs6573857",
  ifelse(creds$snp=="chr14.98011018.A.T","rs79835824",
  ifelse(creds$snp=="chr16.11190252.A.G","rs248832",
  ifelse(creds$snp=="chr16.28494339.G.C","rs151234",
  ifelse(creds$snp=="chr21.44280294.T.A","rs9671675",
  ifelse(creds$snp=="chr21.42418534.G.A","rs13048049",creds$ID))))))))

  #skip if no credible snps
  if(nrow(creds)==0){
    out<-NULL
  }

  #also skip regions with only 1	signal
  tags<-unique(creds$tag)
  if(length(tags)==1){
    out<-NULL
  }

  #and keep top SNP from each group left:
  l<-summx[summx$snp %in% tags,]
  rownames(l)<-l$MarkerName

  #prepare regression/imputation matrices
  samples<-cbind(Y,covariates)
  dat<-DATA[,colnames(DATA) %in% l$MarkerName]
  l<-l[order(l$position),]
  dat<-dat[,rownames(l)]

  #change to minor allele reference:
  cs<-col.summary(dat)
  w<-which(cs$RAF>0.5)
  dat<-switch.alleles(dat,snps=w)
  l[w,"d1"]<-l[w,"a0"]
  l[w,"d0"]<-l[w,"a1"]
  l[w,"a0"]<-l[w,"d0"]
  l[w,"a1"]<-l[w,"d1"]
  l[w,"Allele1"]<-l[w,"d0"]
  l[w,"Allele2"]<-l[w,"d1"]

  #get best guess genotypes:
  dat<-as(dat,"numeric")
  for (i in 1:ncol(dat)){
    dat[,i]<-ifelse(dat[,i]<0.5,0,
    ifelse(dat[,i]>=0.5 & dat[,i]<1.5,1,
    ifelse(dat[,i]>=1.5,2,NA)))
  }
  dat<-as(dat,"SnpMatrix")
  keepfirst<-function(dataframe){
    d<-dataframe[1]
    return(d)
  }
  fix.alleles <- function(long,short) {
    try<-strsplit(long,"") #%>%
    nextone<-mapply(setdiff, try , short) #%>%
    nextone<-mapply(keepfirst, nextone)
    finally<-lapply(nextone,paste0,collape="") #%>%
    trying<-unlist(finally)
    return(trying)
  }
  l1<-nchar(l$Allele1)
  l2<-nchar(l$Allele2)
  l$Allele1<-ifelse(l1>1 & l2==1,fix.alleles(l$Allele1,l$Allele2),
  ifelse(l1>1 & l2>1, substr(l$Allele1,1,1),l$Allele1))
  l$Allele2<-ifelse(l2>1,fix.alleles(l$Allele2, l$Allele1),l$Allele2)

  #write genotype matrix to SNPHAP format
  write.snphap(dat, a1=l$Allele1, a2=l$Allele2,
               file=paste0(out_subdir,"/init.in"))

  #run SNPHAP
  f.in=paste0(out_subdir,"/init.in")
  f.out1 <- paste0(sub(".in$","",f.in),".out1")
  f.out2 <- paste0(sub(".in$","",f.in),".out2")
  message("running snphap (this may take a while)")
  sink(file=paste0(out_subdir, "/init.sh"))
  cat(paste0("/home/ccr5ju/snphap/snphap -q -nh -mi 10 -ss ",f.in," ",f.out1," ",f.out2," > ",out_subdir,"/init"))
  sink()
  system(paste0("chmod a=rwx ",out_subdir, "/init.sh"))
  system(paste0("bash ", out_subdir, "/init.sh"))

  #run multiple imputation
  ss<-model.mi(dir=out_subdir,df=samples, phenotype="outcome", haps.pattern="init.out2",base=base, family="binomial", thr=thr)
  ss$hap<-as.factor(rownames(ss))
  ss<-ss[substr(ss$hap,1,3)=="hap",]

  #generate haplotype plots
  #get colors for plotting
  la<-summx[groups@tags,]
  cols=c("red","blue","green","purple","yellow","brown","pink","turquoise","khaki1","slateblue1","deepskyblue1","linen","firebrick1")
  cols<-cols[c(1:nrow(la))]
  names(cols)<-la$ID
  cols<-cols[names(cols) %in% creds$ID]
  t = plotHaplotypes(ss, thr, l, cols)

  #save
  ggsave(t,file=paste0(out_subdir, "/",gsub(":",".",index),".png"), dpi=600, height=30, width=25, units="cm")
  ggsave(t,file=paste0(out_subdir, "/",gsub(":",".",index),".pdf"),height=10, width=7)

  return(ss)
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Updated function to allow predefined input snps
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

haplorun_custom <-function(index, base, thr, taghap=NULL, hapnickname=NULL){

  # index = "chr2:203874196:G:A"
  # base = NULL
  # thr=0.004
  # taghap=taghap2
  # hapnickname="taghap2"
  out_subdir = file.path(out_dir,gsub(":",".",index),hapnickname)
  dir.create(out_subdir)

  #load genotype matrix and phenotype data (objects: Y, covariates, DATA)
  d <-paste0(jamie_dir,"/guessfm_5pc_diff/input/",gsub(":",".",index),"/")
  load(file=paste0(d,"data.RData"), verbose=TRUE)

  o <- data.frame(matrix(unlist(lapply(colnames(DATA), function(x) unlist(strsplit(x,split=":"))[1:2])), ncol = 2, byrow = TRUE))
  o$chr = gsub("chr","", o$X1); o$pos = o$X2
  write.table(o[,c("chr","pos")], file=paste0(out_subdir, "/all_vars.txt"), sep="\t", col.names=F, row.names=F, quote=F)

  #load guessfm output (objects: groups, l, summx, sm.all, sp.all)
  load(file=paste0(d,"summx.RData"), verbose=TRUE)
  #o<-summx[,c("chromosome","position")]

  #extract variant ids (rsid) for snps in sp.all using Kaviar reference
  system(paste0("/apps/software/standard/core/bcftools/1.9/bin/bcftools view -R ",out_subdir,"/all_vars.txt ", jamie_dir,
  "/Kaviar-160204-Public-hg38-trim.vcf.gz > ",out_subdir, "/all_vars_ids.txt"))

  #add variant ids to guessfm summx object
  vars<-read.table(file=paste0(out_subdir, "/all_vars_ids.txt"),skip=40, header=T, comment.char="", as.is=T,colClasses = c("character"))
  vars$MarkerName<-paste0("chr",vars$X.CHROM,":",vars$POS,":",vars$REF,":",vars$ALT)
  vars<-vars[,c("MarkerName","ID")]
  vars<-vars[vars$ID!=".",]

  # summx<-merge(summx,vars,by="MarkerName", all.x=T)

  # #sort guessfm summx object
  # summx<-summx[order(summx$chromosome, summx$position),]
  # rownames(summx)<-summx$snp
  # summx$ID<-ifelse(summx$snp=="chr2.162279995.C.G","rs35337543",
  # ifelse(summx$snp=="chr12.56042145.C.G","rs1131017",
  # ifelse(summx$snp=="chr14.68796357.G.A","rs6573857",
  # ifelse(summx$snp=="chr14.98011018.A.T","rs79835824",
  # ifelse(summx$snp=="chr16.11190252.A.G","rs248832",
  # ifelse(summx$snp=="chr16.28494339.G.C","rs151234",
  # ifelse(summx$snp=="chr21.44280294.T.A","rs9671675",
  # ifelse(summx$snp=="chr21.42418534.G.A","rs13048049",summx$ID))))))))
  # summx$ID<-ifelse(is.na(summx$ID),".",summx$ID)
  #
  # #load snps to use in haplotype analysis:
  # mydir<-paste0(cred_dir,"/",gsub(":",".",index),"/")
  # creds<-read.table(paste0(mydir,"credible_snps.txt"),as.is=T, header=T)
  # if(index!="chr10:6052734:C:T"){
  # creds<-creds[creds$ppsum>0.5,]
  # }
  # #il2ra v complicated - just keeping those with ppsum>0.8:
  # if(index=="chr10:6052734:C:T"){
  # creds<-creds[creds$ppsum>0.8,]
  # }
  # creds$ID<-ifelse(creds$snp=="chr2.162279995.C.G","rs35337543",
  # ifelse(creds$snp=="chr12.56042145.C.G","rs1131017",
  # ifelse(creds$snp=="chr14.68796357.G.A","rs6573857",
  # ifelse(creds$snp=="chr14.98011018.A.T","rs79835824",
  # ifelse(creds$snp=="chr16.11190252.A.G","rs248832",
  # ifelse(creds$snp=="chr16.28494339.G.C","rs151234",
  # ifelse(creds$snp=="chr21.44280294.T.A","rs9671675",
  # ifelse(creds$snp=="chr21.42418534.G.A","rs13048049",creds$ID))))))))

  # #skip if no credible snps
  # if(nrow(creds)==0){
  #   out<-NULL
  # }
  #
  # #also skip regions with only 1	signal
  # tags<-unique(creds$tag)
  # if(length(tags)==1){
  #   out<-NULL
  # }
  #
  # #and keep top SNP from each group left:
  # l<-summx[summx$snp %in% tags,]
  rownames(l)<-l$MarkerName

  tagvars = vars[vars$ID %in% taghap,]
  tagvars$MarkerName[tagvars$ID=="rs231775"] <- "chr2:203867991:A:G"
  rownames(tagvars) = tagvars$MarkerName

  #check that all tags can be analyzed
  if (sum(!taghap%in%vars$ID)>0) {
    cat("WARNING: some of your tag variants cannot be mapped to available markers in this study\n")
  }
  if (sum(!tagvars$MarkerName%in%l$MarkerName)>0) {
    cat("WARNING: some of your tag variants cannot be mapped to available markers in this study\n")
  }

  l = l[l$MarkerName %in% tagvars$MarkerName,]
  l$ID = tagvars[l$MarkerName,"ID"]

  #prepare regression/imputation matrices
  samples<-cbind(Y,covariates)
  dat<-DATA[,colnames(DATA) %in% l$MarkerName]
  l<-l[order(l$position),]
  dat<-dat[,rownames(l)]

  #change to minor allele reference:
  cs<-col.summary(dat)
  w<-which(cs$RAF>0.5)
  dat<-switch.alleles(dat,snps=w)
  l[w,"d1"]<-l[w,"a0"]
  l[w,"d0"]<-l[w,"a1"]
  l[w,"a0"]<-l[w,"d0"]
  l[w,"a1"]<-l[w,"d1"]
  l[w,"Allele1"]<-l[w,"d0"]
  l[w,"Allele2"]<-l[w,"d1"]

  #get best guess genotypes:
  dat<-as(dat,"numeric")
  for (i in 1:ncol(dat)){
    dat[,i]<-ifelse(dat[,i]<0.5,0,
    ifelse(dat[,i]>=0.5 & dat[,i]<1.5,1,
    ifelse(dat[,i]>=1.5,2,NA)))
  }
  dat<-as(dat,"SnpMatrix")
  keepfirst<-function(dataframe){
    d<-dataframe[1]
    return(d)
  }
  fix.alleles <- function(long,short) {
    try<-strsplit(long,"") #%>%
    nextone<-mapply(setdiff, try , short) #%>%
    nextone<-mapply(keepfirst, nextone)
    finally<-lapply(nextone,paste0,collape="") #%>%
    trying<-unlist(finally)
    return(trying)
  }
  l1<-nchar(l$Allele1)
  l2<-nchar(l$Allele2)
  l$Allele1<-ifelse(l1>1 & l2==1,fix.alleles(l$Allele1,l$Allele2),
  ifelse(l1>1 & l2>1, substr(l$Allele1,1,1),l$Allele1))
  l$Allele2<-ifelse(l2>1,fix.alleles(l$Allele2, l$Allele1),l$Allele2)

  #write genotype matrix to SNPHAP format
  write.snphap(dat, a1=l$Allele1, a2=l$Allele2,
               file=paste0(out_subdir,"/init.in"))

  #run SNPHAP
  f.in=paste0(out_subdir,"/init.in")
  f.out1 <- paste0(sub(".in$","",f.in),".out1")
  f.out2 <- paste0(sub(".in$","",f.in),".out2")
  message("running snphap (this may take a while)")
  sink(file=paste0(out_subdir, "/init.sh"))
  cat(paste0("/home/ccr5ju/snphap/snphap -q -nh -mi 10 -ss ",f.in," ",f.out1," ",f.out2," > ",out_subdir,"/init"))
  sink()
  system(paste0("chmod a=rwx ",out_subdir, "/init.sh"))
  system(paste0("bash ", out_subdir, "/init.sh"))

  #run multiple imputation
  ss<-model.mi(dir=out_subdir,df=samples, phenotype="outcome", haps.pattern="init.out2",base=base, family="binomial", thr=thr)
  ss$hap<-as.factor(rownames(ss))
  ss<-ss[substr(ss$hap,1,3)=="hap",]

  #generate haplotype plots
  t = plotHaplotypes_nocolor(ss, thr, l)

  #save
  ggsave(t,file=paste0(out_subdir, "/haplotype_plots.png"), dpi=600, height=30, width=25, units="cm")
  ggsave(t,file=paste0(out_subdir, "/haplotype_plots.pdf"),height=10, width=7)

  return(ss)
}
