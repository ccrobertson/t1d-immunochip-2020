#install.packages("rtf")
library(rtf)
library(dplyr)
setwd(Sys.getenv("manuscript"))

pheno = read.table(Sys.getenv("phenofile"),header=FALSE)
names(pheno) <- c("Cohort", "FID","IID","FAT","MOT","SEX","T1D")
clust = read.table(paste0(Sys.getenv("pca"),"/mega_cluster_assignments.txt"), header=TRUE)

m = merge(pheno, clust, by='IID')

### Assign each Cohort to a country
m$Country = NA
m$Country[m$Cohort%in%c("BC58", "BDA","cbr_controls","GRID","UKBS","T1DGC-UK")] <- "United Kingdom"
m$Country[m$Cohort%in%c("BRI","CLEAR","EDIC","GoKinD", "NYCP","NIMH","SEARCH", "Trialnet","T1DGC-NA", "T1DGC","UAB","UCHSC2", "UCSF_AA")] <- "United States"
m$Country[m$Cohort%in%c("IDDMGEN", "T1DGEN", "MILWAUKEE")] <- "Finland"
m$Country[m$Cohort%in%c("T1DGC-NI", "T1DGC-YH")] <- "Northern Ireland"
m$Country[m$Cohort%in%c("DAN-SDC","IC_Cases_Steno","IC_Cases_HSG")] <- "Denmark"
m$Country[m$Cohort%in%c("T1DGC-EUR","Existing_Source_to_T1DGC:SWEDISH","HapMap")] <- "Other Europe"
m$Country[m$Cohort%in%c("T1DGC-AP")] <- "Asia Pacific"
#sum(is.na(m$Country))

### Define Cohort factor
#!!! Ask for input on how to label these
#DAN-SDC --> Steno Diabetes Center?
#IC_Cases_Steno --> ??
#IC_Cases_HSG --> ?
cohort_vals=c("BC58", "BDA","cbr_controls","GRID","UKBS","T1DGC-UK","BRI","CLEAR","EDIC","GoKinD", "NYCP","NIMH","SEARCH", "Trialnet","T1DGC-NA", "T1DGC","UAB","UCHSC2", "UCSF_AA","IDDMGEN", "T1DGEN", "MILWAUKEE","T1DGC-NI", "T1DGC-YH","DAN-SDC","IC_Cases_Steno","IC_Cases_HSG","T1DGC-EUR","HapMap","T1DGC-AP")
cohort_labels=c("BC58", "BDA","Cambridge","GRID","UKBS","T1DGC-UK","BRI","CLEAR","EDIC","GoKinD", "NYCP","NIMH","SEARCH", "Trialnet","T1DGC-NA", "T1DGC","Univ. of Alabama","Univ. of Colorado", "UCSF","IDDMGEN", "T1DGEN", "MILWAUKEE","T1DGC-NI", "T1DGC-YH","DAN-SDC","IC_Cases_Steno","IC_Cases_HSG","T1DGC-EUR","HapMap","T1DGC-AP")
m$CohortFactor = factor(m$Cohort, levels=cohort_vals, labels=cohort_labels)

### TOTAL COHORTS
sumTable1 <- m %>% group_by(CohortFactor) %>%
		summarise(
				Country=unique(Country),
				Total = n(),
				Cases=sum(T1D==2),
				Controls=sum(T1D==1),
				Missing=sum(T1D==0),
				AFR=sum(cluster_label=="AFR"),
				AMR=sum(cluster_label=="AMR"),
				EAS=sum(cluster_label=="EAS"),
				EUR=sum(cluster_label=="EUR"),
				FIN=sum(cluster_label=="FIN"))

rtffile <- RTF("complete_cohort_summary_table.doc")
addParagraph(rtffile, "Table 1 - Summary of cohorts by T1D status and ancestry")
addTable(rtffile, sumTable1)
done(rtffile)


### CASES, CONTROLS, AND TRIOS BY COHORT
trios = read.table(paste0(Sys.getenv("fam_assoc"), "/affected_offspring.txt"), header=TRUE)
trios_with_info = merge(trios, m, by.x="iid", by.y="IID")
cc = read.table(paste0(Sys.getenv("cc_assoc"),"/unrelateds_all.txt"), header=TRUE)
cc_with_info = merge(cc, m, by.x="iid",by.y="IID")

sumTableCC <- cc_with_info %>% group_by(CohortFactor) %>%
		summarise(Cases=sum(T1D==2),
				Controls=sum(T1D==1))
sumTableTrios <- trios_with_info %>% group_by(CohortFactor) %>%
		summarise(Trios=n())
sumTable2 <- full_join(sumTableCC, sumTableTrios, by="CohortFactor") %>%
replace(is.na(sumTable2), 0)

rtffile <- RTF("cc_and_trios_by_cohort.doc")
addParagraph(rtffile, "Table 2 - Number of unrelated cases, unrelated controls, and trios by cohort")
addTable(rtffile, sumTable2)
done(rtffile)


### CASES, CONTROLS, AND TRIOS BY ANCESTRY
sumTableCC_ancestry <- cc_with_info %>% group_by(cluster_label) %>%
		summarise(Cases=sum(T1D==2),
				Controls=sum(T1D==1))
sumTableTrios_ancestry <- trios_with_info %>% group_by(cluster_label) %>%
		summarise(Trios=n())
sumTable3 <- full_join(sumTableCC_ancestry, sumTableTrios_ancestry, by="cluster_label")

rtffile <- RTF("cc_and_trios_by_ancestry.doc")
addParagraph(rtffile, "Table 3 - Number of unrelated cases, unrelated controls, and trios by genotype-based ancestry cluster")
addTable(rtffile, sumTable3)
done(rtffile)
