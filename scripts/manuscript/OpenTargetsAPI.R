## Exploring the GraphQL API R implementation for accessing the Open Targets database
# Based on this Tweetorial
# https://twitter.com/enricoferrero/status/1060930693239357441
setwd(Sys.getenv("manuscript"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set up API
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#install.packages("ghql")
library("ghql")
library("jsonlite")
library(ggplot2)
library(RColorBrewer)

# create GraphQL client
cli <- GraphqlClient$new(
  url = "https://genetics-api.opentargets.io/graphql"
)

# create query object
qry <- Query$new()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run on one example SNP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# build example query
qry$query('my_query', '{
  indexVariantsAndStudiesForTagVariant(variantId:"1_154453788_C_T") {
    associations {
      indexVariant {id,rsId},
      study {studyId,pmid,pubAuthor, pubDate, traitReported},
      pval,
      beta,
      betaCILower,
      betaCIUpper,
      oddsRatio,
			oddsRatioCILower,
      oddsRatioCIUpper,
      nTotal,
      overallR2,
      posteriorProbability
    }
  }
}')

# run query and format results:
#cli$exec(qry$queries$my_query)
#res <- fromJSON(cli$exec(qry$queries$my_query), flatten = TRUE)$data$indexVariantsAndStudiesForTagVariant


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run across all SNPs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create wrapper function
getVariantInfo = function(variant) {
  cat("Getting info for", variant,"\n")
  query_name = paste0("query_",variant)
  query = paste0('{
    indexVariantsAndStudiesForTagVariant(variantId:"',variant,'") {
      associations {
        indexVariant {id,rsId},
        study {studyId,pmid,pubAuthor, pubDate, traitReported},
        pval,
        beta,
        betaCILower,
        betaCIUpper,
        oddsRatio,
  			oddsRatioCILower,
        oddsRatioCIUpper,
        nTotal,
        overallR2,
        posteriorProbability
      }
    }
  }')
  qry$query(query_name, query)
  res <- fromJSON(cli$exec(qry$queries[[query_name]]), flatten = TRUE)$data$indexVariantsAndStudiesForTagVariant[[1]]
  if (class(res)=="data.frame") {
    res$input_var = variant
  }
  return(res)
}

# get lead variants
t1d_dat = read.table(file.path(Sys.getenv('manuscript'), 'SUPPLEMENTARY_TABLES_ST8.txt'), header=TRUE, sep="\t")
t1d_dat$t1d_lead = gsub("chr","", gsub(":","_", t1d_dat$MarkerName))


# apply wrapper function across leads
res_OpenTargets = lapply(t1d_dat$t1d_lead, getVariantInfo)


# # investigate snps with no results
# missing_variants = t1d_dat$t1d_lead[unlist(lapply(res, function(x) !is.data.frame(x)))]
# getGeneInfo = function(variant) {
#   cat("Getting info for", variant,"\n")
#   query_name = paste0("gene_query_",variant)
#   query = paste0('{genesForVariant(variantId:"',variant,'") {variant,gene {id,symbol},overallScore}}')
#   qry$query(query_name, query)
#   res <- fromJSON(cli$exec(qry$queries[[query_name]]), flatten = TRUE)$data$genesForVariant
#   return(res)
# }
# missing_res = lapply(missing_variants, getGeneInfo)
# table(unlist(lapply(missing_res, class)))
# table(unlist(lapply(missing_res, function(x) nrow(x)>0)))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Concatenate and format results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# concatenate all data
res_all_tmp = do.call("rbind", res_OpenTargets)

# order variants by chromosome positions
res_all_tmp2 = merge(res_all_tmp, t1d_dat[,c("t1d_lead","MarkerName","ID","chromosome","position","Candidate.Gene","Novel","EffectphaseIandII")],all.x=TRUE, all.y=FALSE, by.x="input_var",by.y="t1d_lead")
res_all = res_all_tmp2[with(res_all_tmp2, order(chromosome, position)),]
res_all$leadVar = paste0(res_all$MarkerName, " (",res_all$Candidate.Gene,")")
res_all$plotLeadVar = factor(res_all$leadVar, levels=unique(res_all$leadVar))

# remove entries with missing LD info
res_all = res_all[!is.na(res_all$overallR2),]

# create discrete R2 variable
res_all$Rsquared = NA
res_all$Rsquared[res_all$overallR2 >= 0.5 & res_all$overallR2<0.55] <- 0.5
res_all$Rsquared[res_all$overallR2 >= 0.55 & res_all$overallR2<0.65] <- 0.6
res_all$Rsquared[res_all$overallR2 >= 0.65 & res_all$overallR2<0.75] <- 0.7
res_all$Rsquared[res_all$overallR2 >= 0.75 & res_all$overallR2<0.85] <- 0.8
res_all$Rsquared[res_all$overallR2 >= 0.85 & res_all$overallR2<0.95] <- 0.9
res_all$Rsquared[res_all$overallR2 >= 0.95] <- 1
res_all$Rsquared = factor(res_all$Rsquared, levels=c(0.5, 0.6, 0.7, 0.8, 0.9, 1), labels=c("0.5","0.6","0.7","0.8","0.9",">0.95"))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract Immune Disease results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get immunobase study ids (NOTE: these were curated up until 2015)
immunobase = read.csv("ImmunoBase_study_links_immunobase_curated_table.csv")
immunobase$GWAS.Catalog.accession %in% res_all$study.studyId
res_immunobase = res_all[res_all$study.studyId %in% immunobase$GWAS.Catalog.accession,]
unique(res_immunobase$study.traitReported)

# get other studies with immunobase traits
getRelatedTraits = function(pattern) {
  related_traits = res_all$study.traitReported[grep(pattern, res_all$study.traitReported, ignore.case = TRUE)]
  return(unique(related_traits))
}

# manually define traits
res_all$ImmuneDisease = NA
res_all$ImmuneDisease[res_all$study.traitReported %in% c("Systemic lupus erythematosus","Systemic lupus erythematosus [EA]", "Systemic lupus erythematosus [Chinese]","Systemic lupus erythematosus [Hispanic]")] <- "Systemic lupus erythematosus"
res_all$ImmuneDisease[res_all$study.traitReported %in% c("Type 1 diabetes","Type 1 diabetes | non-cancer illness code, self-reported")] <- "Type 1 diabetes"
res_all$ImmuneDisease[res_all$study.traitReported %in% c("Rheumatoid arthritis","Rheumatoid arthritis | non-cancer illness code, self-reported","Rheumatoid arthritis [EA]","Rheumatoid arthritis [East Asian]")] <- "Rheumatoid arthritis"
res_all$ImmuneDisease[res_all$study.traitReported %in% c("Juvenile idiopathic arthritis (oligoarticular or rheumatoid factor-negative polyarticular) [additive model]")] <- "Juvenile idiopathic arthritis"
res_all$ImmuneDisease[res_all$study.traitReported %in% c("Graves' disease")] <- "Graves' disease"
res_all$ImmuneDisease[res_all$study.traitReported %in% c("Hashimoto thyroiditis")] <- "Hashimoto thyroiditis"
res_all$ImmuneDisease[res_all$study.traitReported %in% c("Crohn's disease", "Crohn's disease [EA]")] <- "Crohn's disease"
res_all$ImmuneDisease[res_all$study.traitReported %in% c("Celiac disease")] <- "Celiac disease"
res_all$ImmuneDisease[res_all$study.traitReported %in% c("Multiple sclerosis")] <- "Multiple sclerosis"
res_all$ImmuneDisease[res_all$study.traitReported %in% c("Ulcerative colitis","Ulcerative colitis [EA]","Ulcerative colitis | non-cancer illness code, self-reported")] <- "Ulcerative colitis"
res_all$ImmuneDisease[res_all$study.traitReported %in% c("Vitiligo")] <- "Vitiligo"
res_all$ImmuneDisease[res_all$study.traitReported %in% c("Psoriasis","Inflammatory skin disease [Psoriasis]","Psoriasis vulgaris","Psoriasis | non-cancer illness code, self-reported")] <- "Psoriasis"
res_all$ImmuneDisease[res_all$study.traitReported %in% c("Primary biliary cirrhosis","Primary biliary cholangitis")] <- "Primary biliary cholangitis"
res_all$ImmuneDisease[res_all$study.traitReported %in% c("Primary sclerosing cholangitis")] <- "Primary sclerosing cholangitis"
res_all$ImmuneDisease[res_all$study.traitReported %in% c("Narcolepsy","Daytime dozing / sleeping (narcolepsy)")] <- "Narcolepsy"


# for each ImmuneDisease:input_var combo, only take the entry with lowest p-value
res_immune_all = res_all[!is.na(res_all$ImmuneDisease),]
res_immune_all$unique_id = paste0(res_immune_all$input_var,":",res_immune_all$ImmuneDisease)
res_immune_all_sorted = res_immune_all[with(res_immune_all, order(unique_id, pval)),]

#write.table(res_immune_all_sorted, file="t1d_overlap_with_immunobase_diseases.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
res_immune = res_immune_all_sorted[!duplicated(res_immune_all_sorted$unique_id),]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract Biomarker results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
biomarkers = read.table("biomarkers_clin.txt", sep="\t", header=FALSE)
names(biomarkers) <- c("study.traitReported", "biomarker" )
res_biomarkers = merge(biomarkers, res_all, by="study.traitReported")

# for each biomarker:input_var combo, only take the entry with lowest p-value
res_biomarkers$unique_id = paste0(res_biomarkers$input_var,":",res_biomarkers$biomarker)
res_biomarkers_sorted = res_biomarkers[with(res_biomarkers, order(unique_id, pval)),]

res_biomarkers_u = res_biomarkers_sorted[!duplicated(res_biomarkers_sorted$unique_id),]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run query for summary stats for each T1D lead snp
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For each snp-trait pair in res_all we need the summary stats for the association
# between the t1d lead variant and the trait
getVariantSummaryStats = function(variant, studyId) {
  cat("Getting info for", variant,"\n")
  query_name = paste0("query_sumstats_",variant)
  query = paste0('{
    pheWAS(variantId:"',variant,'") {
      associations {
        study {
          studyId,
          pmid,
          pubDate,
          pubAuthor,
          hasSumsStats,
          nCases,
          traitReported
        },
        beta,
        oddsRatio,
        eaf,
        nTotal,
        nCases
      }
    }
  }')
  qry$query(query_name, query)
  res <- fromJSON(cli$exec(qry$queries[[query_name]]), flatten = TRUE)$data$pheWAS[[1]]
  if (class(res)=="data.frame") {
    res$input_var = variant
    res = res[res$study.studyId==studyId,]
  }
  return(res)
}

sumstats_immune = list()
for ( i in 1:dim(res_immune)[1]) {
  variant = res_immune[i,"input_var"]
  study = res_immune[i,"study.studyId"]
  cat(variant, study,"\n")
  sumstats_immune[[i]] = getVariantSummaryStats(variant, study)
}
sumstats_immune_df = do.call("rbind", sumstats_immune)
res_immune2 = merge(res_immune, sumstats_immune_df, suffixes = c("",".trait2"), all.x=TRUE, by=c("input_var","study.studyId"))

sumstats_biomarkers = list()
for ( i in 1:dim(res_biomarkers_u)[1]) {
  variant = res_biomarkers_u[i,"input_var"]
  study = res_biomarkers_u[i,"study.studyId"]
  cat(variant, study,"\n")
  sumstats_biomarkers[[i]] = getVariantSummaryStats(variant, study)
}
sumstats_biomarkers_df = do.call("rbind", sumstats_biomarkers)
res_biomarker2 = merge(res_biomarkers_u, sumstats_biomarkers_df, suffixes = c("",".trait2"), all=TRUE, by=c("input_var","study.studyId"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# If T1D lead summary stat is unavailable, get T1D association for lead SNP from other trait
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
res_immune2_missing = res_immune2[is.na(res_immune2$beta.trait2) & is.na(res_immune2$oddsRatio.trait2),]

load(file.path(Sys.getenv("meta"),"dens/resultsmeta_ind_fam_tdt_5pc_diff_dens.RData"))
t1d_meta_stats = res
t1d_meta_stats$MarkerName2_tmp = gsub(":","_", t1d_meta_stats$MarkerName)
t1d_meta_stats$MarkerName2 = gsub("chr","",t1d_meta_stats$MarkerName2_tmp)

sum(res_immune2_missing$indexVariant.id %in% t1d_meta_stats$MarkerName2)
sum(!res_immune2$input_var %in% t1d_meta_stats$MarkerName2)

row.names(t1d_meta_stats) <- t1d_meta_stats$MarkerName2
res_immune2$indexVariant_trait2.t1d_beta = t1d_meta_stats[res_immune2$indexVariant.id,"Effect"]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create indicator of concordance between T1D and other traits
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#res_immune2$trait_sign = NA
#res_immune2$trait_sign[(!is.na(res_immune2$beta.trait2) & res_immune2$beta.trait2>0) | (!is.na(res_immune2$oddsRatio.trait2) & res_immune2$oddsRatio.trait2>1) ] <- 1
#res_immune2$trait_sign[(!is.na(res_immune2$beta.trait2) & res_immune2$beta.trait2<0) | (!is.na(res_immune2$oddsRatio.trait2) & res_immune2$oddsRatio.trait2<1) ] <- -1
#res_immune2$concordance = sign(res_immune2$EffectphaseIandII * res_immune2$trait_sign)


res_immune2$trait_sign = NA
res_immune2$trait_sign[(!is.na(res_immune2$beta) & res_immune2$beta>0) | (!is.na(res_immune2$oddsRatio) & res_immune2$oddsRatio>1) ] <- 1
res_immune2$trait_sign[(!is.na(res_immune2$beta) & res_immune2$beta<0) | (!is.na(res_immune2$oddsRatio) & res_immune2$oddsRatio<1) ] <- -1
res_immune2$concordance = sign(res_immune2$indexVariant_trait2.t1d_beta * res_immune2$trait_sign)


res_immune2[,c("unique_id","beta.trait2","oddsRatio.trait2","trait_sign")]
head(res_immune2[,c("unique_id","beta.trait2","oddsRatio.trait2","trait_sign","study.hasSumsStats")])
head(res_immune2[,c("study.traitReported","study.traitReported.trait2")])
res_immune2$study.traitReported==res_immune2$study.traitReported.trait2


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot LD between lead variants
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pdf("ld_with_immune_traits.pdf", height=23, width=14)
ggplot(res_immune2, aes(y=plotLeadVar, x=ImmuneDisease)) +
  geom_tile(aes(fill=overallR2*concordance), alpha=0.5) +
  geom_text(aes(label=sprintf("%.2f", overallR2)), size=3) +
  #scale_fill_brewer(name="R-squared with T1D lead") +
  scale_fill_gradient2(midpoint=0,low="blue",high="red", mid="white") +
  ylab("T1D Lead Variant") +
  xlab(" ") +
  theme(axis.text.x=element_text(angle=90, size=20))
dev.off()

pdf("ld_with_immune_traits_novel_only.pdf", width=13, height=7)
ggplot(res_immune3[!is.na(res_immune3$Novel) & res_immune3$Novel==1,], aes(x=plotLeadVar, y=ImmuneDisease)) +
  geom_tile(aes(fill=Rsquared), alpha=0.5) +
  geom_text(aes(label=sprintf("%.2f", overallR2)), size=3) +
  scale_fill_brewer(name="R-squared with T1D lead") +
  xlab("T1D Lead Variant") +
  ylab(" ") +
  theme(axis.text.x=element_text(angle=90))
dev.off()


pdf("ld_with_biomarkers.pdf", width=13, height=7)
ggplot(res_biomarker2, aes(x=plotLeadVar, y=biomarker)) +
  geom_tile(aes(fill=overallR2*concordance), alpha=0.5) +
  geom_text(aes(label=sprintf("%.2f", overallR2)), size=3) +
  #scale_fill_brewer(name="R-squared with T1D lead") +
  scale_fill_gradient2(midpoint=0,low="blue",high="red", mid="white") +
  xlab("T1D Lead Variant") +
  ylab(" ") +
  theme(axis.text.x=element_text(angle=90))
dev.off()
