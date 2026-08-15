#==============================================================================
# STAAR results summary - BATCH 2
#==============================================================================
# MERGED FROM (chronological):
#   - 20250920_STAAR_summary_scaledup_b2.R   ->  section "v2_scaled_up"
#   - 20251209_STAAR_summary_b2.R   ->  section "v3_final_b2"
# NOTE: FLAG: the original 20251209_STAAR_summary_b2.R defined summaryCoding() TWICE
# NOTE: (lines ~64 and ~330). Both kept here; R will use whichever was defined last before
# NOTE: the call. Review and delete the stale one.
#==============================================================================

#==============================================================================
# MAIN PIPELINE  (current version: v3_final_b2)
#==============================================================================

library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)
library(tidyverse)
library(GWASTools)

setwd("/Volumes/Sofer Lab/HCHS_SOL/Projects/2024_rare_variants/Code")

###########################################################
#           User Input
## Number of jobs for each chromosome
jobs_num <- get(load("../Data/STAAR_prep/mult_chr_jobs_num_b2.Rdata"))
## aGDS directory
agds_dir <- get(load("../Data/STAAR_prep/agds_dir_multichr.Rdata"))

## QC_label
QC_label <- "annotation/filter"
## variant_type
variant_type <- "variant"
## geno_missing_imputation
geno_missing_imputation <- "mean"


############ Coding summary ############
## rs channel name in aGDS
rs_channel <- "annotation/info/FunctionalAnnotation/rsid"
## output path
output_path <- "../Data/"
## number of jobs
gene_centric_coding_jobs_num <- 1
## results name
gene_centric_results_name <- "Genecentric_coding_uncond_batch2"

variant_type <- "SNV"
## method_cond
method_cond <- "optimal"
alpha <- 0.05/1562# for details see STAAR_multiple_chrs

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
Annotation_name_catalog <- get(load("../Data/STAAR_prep/Annotation_name_catalog.Rdata"))
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")


#           Main Function
source("20250422_gencentric_coding_summary.R")

summaryCoding <- function(agds.path, chr, group, metab, job.num.chr,
                          known_loci = NULL,
                          p_threshold_mode = "bonferroni"){
  if(group == "exp_groups"){
    ## output path
    output_path <- paste0("../Data/STAAR_exp_groups/",
                          metab, "/b2/")
    load(paste0("../Data/STAAR_exp_groups/",
                metab, "/", metab, "_nullmodel_batch2.Rdata"))
  }else{
    output_path <- paste0("../Data/STAAR_neg_controls/",
                          metab, "/b2/")
    load(paste0("../Data/STAAR_neg_controls/",
                metab, "/", metab, "_nullmodel_batch2.Rdata"))
  }
  input_path <- output_path
  coding_summary <- Gene_Centric_Coding_Results_Summary(agds_dir=agds_dir,gene_centric_coding_jobs_num=1,
                                                        input_path=input_path,output_path=output_path, chr = chr,
                                                        gene_centric_results_name=paste0(metab, "_", gene_centric_results_name),
                                                        obj_nullmodel=obj_nullmodel,
                                                        method_cond=method_cond,
                                                        QC_label=QC_label,geno_missing_imputation=geno_missing_imputation,variant_type=variant_type,
                                                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                        alpha = alpha,manhattan_plot=FALSE,QQ_plot=FALSE)
  return(coding_summary)
}


summary_coding_chr2_uncond_b2 <- summaryCoding(agds_dir[8], chr = 2,
                                                           metab = "x100001721",
                                                           group = "neg_controls")

summary_coding_chr5_uncond_b2  <- summaryCoding(agds_dir[9], chr = 5, metab = "x799",
                                                            group = "neg_controls")

summary_coding_chr6_uncond_b2  <- summaryCoding(agds_dir[10], chr = 6, metab = "x1215",
                                                            group = "neg_controls")

summary_coding_chr16_uncond_b2  <- summaryCoding(agds_dir[14], chr = 16, metab = "x278",
                                                             group = "neg_controls")

## output path - experimental groups
output_path <- "../Data/STAAR_exp_groups/"

# 5-oxoproline
summary_coding_chr8_uncond_b2 <- summaryCoding(agds_dir[11], group = "exp_groups",
                                                           metab = "x1021", chr = 8)
# Carnitine
summary_coding_chr10_uncond_b2 <- summaryCoding(agds_dir[12], group ="exp_groups",
                                                            metab = "x100000007", chr = 10)
# N-acetylcarnosine
summary_coding_chr13_uncond_b2 <- summaryCoding(agds_dir[13], group ="exp_groups",
                                                            metab = "x100004046", chr = 13)
# Cys-gly, oxidized
summary_coding_chr16_uncond_exp_b2 <- summaryCoding(agds_dir[14], group = "exp_groups",
                                                                metab = "x1224", chr = 16)


# Using bonferroni correction - harsh threshold
categories <- c("plof","plof_ds","missense","disruptive_missense",
                "synonymous")
metabs <- c("x100000007", "x100004046",  "x100009332", "x100006370",
            "x100001208", "x2054", "x100006264", "x1021",
            "x100001266", "x1114", "x192", "x2054")
group <- "exp_groups"

num_tests <- 0
for(metab in metabs){
  for(cat_func in categories){
    if(group == "exp_groups"){
      ## output path
      output_path <- paste0("../Data/STAAR_exp_groups/",
                            metab, "/b2/")
      res <- getobj(paste0("../Data/STAAR_exp_groups/",
                           metab, "/b2/", cat_func, ".Rdata"))
    }else{
      output_path <- paste0("../Data/STAAR_neg_controls/",
                            metab, "/b2/")
      res <- getobj(paste0("../Data/STAAR_neg_controls/",
                           metab, "/b2/", cat_func, ".Rdata"))
    }
    if(is.null(res)){
      num_tests <- num_tests
    }else if(!is.null(res) & is.null(dim(res))){
      num_tests <- num_tests + 1
    } else{
      num_tests <- num_tests + dim(res)[1]
    }
  }
}

group <- "neg_controls"
metabs <- c("x1215", "x100001721", "x278", "x799")
for(metab in metabs){
  for(cat_func in categories){
      output_path <- paste0("../Data/STAAR_neg_controls/",
                            metab, "/b2/")
      res <- getobj(paste0("../Data/STAAR_neg_controls/",
                           metab, "/b2/", cat_func, ".Rdata"))
      print(res)
    if(is.null(res)){
      num_tests <- num_tests
    }else if(!is.null(res) & is.null(dim(res))){
      num_tests <- num_tests + 1
    } else{
      num_tests <- num_tests + dim(res)[1]
    }
  }
}

# number of coding tests -- 1845
alpha_coding_uncond_b2 <- 0.05/num_tests


### noncoding ###

## ncRNA
ncRNA_results_name <- "Genetic_Noncoding_ncRNA_b2"

categories=c("downstream","upstream","UTR","promoter_CAGE",
             "promoter_DHS","enhancer_CAGE","enhancer_DHS")
num_tests_noncoding <- 0

for(metab in metabs){
  for(cat_func in categories){
    if(group == "exp_groups"){
      ## output path
      output_path <- paste0("../Data/STAAR_exp_groups/",
                            metab, "/b2/")
      res <- getobj(paste0("../Data/STAAR_exp_groups/",
                           metab, "/b2/", cat_func, ".Rdata"))
    }else{
      output_path <- paste0("../Data/STAAR_neg_controls/",
                            metab, "/b2/")
      res <- getobj(paste0("../Data/STAAR_neg_controls/",
                           metab, "/b2/", cat_func, ".Rdata"))
    }
    if(is.null(res)){
      num_tests_noncoding <- num_tests_noncoding
    }else if(!is.null(res) & is.null(dim(res))){
      num_tests_noncoding <- num_tests_noncoding + 1
    } else{
      num_tests_noncoding <- num_tests_noncoding + dim(res)[1]
    }
  }
}


num_tests_ncrna <- 0
for(metab in metabs){
  for(cat_func in categories){
    ## output path
    output_path <- paste0("../Data/STAAR_exp_groups/",
                          metab, "/b2/")
    ncrna_res <- getobj(paste0("../Data/STAAR_exp_groups/",
                               metab, "/b2/", ncRNA_results_name, ".Rdata"))
    if(is.null(ncrna_res)){
      num_tests_ncrna <- num_tests_ncrna
    }else if(!is.null(ncrna_res) & is.null(dim(ncrna_res))){
      num_tests_ncrna <- num_tests_ncrna + 1
    } else{
      num_tests_ncrna <- num_tests_ncrna + dim(ncrna_res)[1]
    }
  }
}


for(metab in metabs){
  for(cat_func in categories){
    if(group == "exp_groups"){
      ## output path
      output_path <- paste0("../Data/STAAR_exp_groups/",
                            metab, "/b2/")
      res <- getobj(paste0("../Data/STAAR_exp_groups/",
                           metab, "/b2/", cat_func, ".Rdata"))
    }else{
      output_path <- paste0("../Data/STAAR_neg_controls/",
                            metab, "/b2/")
      res <- getobj(paste0("../Data/STAAR_neg_controls/",
                           metab, "/b2/", cat_func, ".Rdata"))
      print(res)
    }
    if(is.null(res)){
      num_tests_noncoding <- num_tests_noncoding
    }else if(!is.null(res) & is.null(dim(res))){
      num_tests_noncoding <- num_tests_noncoding + 1
    } else{
      num_tests_noncoding <- num_tests_noncoding + dim(res)[1]
    }
  }
}


for(metab in metabs){
  for(cat_func in categories){
    ## output path
    output_path <- paste0("../Data/STAAR_neg_controls/",
                          metab, "/b2/")
    ncrna_res <- getobj(paste0("../Data/STAAR_neg_controls/",
                               metab, "/b2/", ncRNA_results_name, ".Rdata"))
    if(is.null(ncrna_res)){
      num_tests_ncrna <- num_tests_ncrna
    }else if(!is.null(ncrna_res) & is.null(dim(ncrna_res))){
      num_tests_ncrna <- num_tests_ncrna + 1
    } else{
      num_tests_ncrna <- num_tests_ncrna + dim(ncrna_res)[1]
    }
  }
}


# ncrna number of tests -- 728
alpha_ncrna_b2 <- 0.05/num_tests_ncrna

# 3412 tests in total for noncoding regions
# therefore the threshold for all non-coding rv sets are
alpha_uncond_noncoding_b2 <- 0.05/num_tests_noncoding

jobs_num <- get(load("../Data/STAAR_prep/scale_up/mult_chr_jobs_num_scale_up.Rdata"))

gene_centric_results_name <- "Genecentric_coding_uncond_b2"
#gene_centric_results_name <- "Genecentric_coding_cond" # conditional analysis

alpha <- alpha_coding_uncond_b2# for details see STAAR_multiple_chrs

load("../Data/STAAR_prep/agds_dir_scale_up_final_b2.Rdata")

summary_chr11_coding_uncond_arachido <- summaryCoding(agds_dir[1], chr = 11,
                                                      metab = "x100009332",
                                                      group = "exp_groups")

summary_chr11_coding_uncond_3beta <- summaryCoding(agds_dir[2], chr = 11,
                                                   metab = "x100006370",
                                                   group = "exp_groups")

summary_chr12_coding_uncond_1_methyl <- summaryCoding(agds_dir[3], chr = 12,
                                                      metab = "x100001208",
                                                      group = "exp_groups")

summary_chr12_coding_uncond_ethylmalonate <- summaryCoding(agds_dir[4], chr = 12,
                                                           metab = "x2054",
                                                           group = "exp_groups")


summary_chr16_coding_uncond_propyl <- summaryCoding(agds_dir[5], chr = 16,
                                                    metab = "x100006264",
                                                    group = "exp_groups")

summary_chr2_coding_uncond_arginine <- summaryCoding(agds_dir[6], chr = 2,
                                                     metab = "x100001266",
                                                     group = "exp_groups")

summary_chr5_coding_uncond_3amino <- summaryCoding(agds_dir[7], chr = 5,
                                                   metab = "x1114",
                                                   group = "exp_groups")

summary_chr8_coding_uncond_putrescine <- summaryCoding(agds_dir[8], chr = 8,
                                                       metab = "x192",
                                                       group = "exp_groups")

summary_chr2_coding_b2 <- summaryCoding(agds_dir[8], chr = 2,
                                     metab = "x100001721",
                                     group = "neg_controls")

summary_chr5_coding <- summaryCoding(agds_dir[9], chr = 5, metab = "x799",
                                     group = "neg_controls")

summary_chr6_coding <- summaryCoding(agds_dir[10], chr = 6, metab = "x1215",
                                     group = "neg_controls")

summary_chr16_coding <- summaryCoding(agds_dir[14], chr = 16, metab = "x278",
                                      group = "neg_controls")

summary_chr8_coding <- summaryCoding(agds_dir[11], chr = 8, metab = "x1021",
                                     group = "exp_groups")

summary_chr10_coding <- summaryCoding(agds_dir[12], chr = 10, metab = "x100000007",
                                      group = "exp_groups")

summary_chr13_coding <- summaryCoding(agds_dir[13], chr = 13, metab = "x100004046",
                                      group = "exp_groups")

summary_chr16_coding_exp <- summaryCoding(agds_dir[14], chr = 16, metab = "x1224",
                                          group = "exp_groups")

coding_rvs_uncond_sig <- rbind(summary_chr2_coding_b2,
                               summary_chr6_coding,
                               summary_chr16_coding,
                               summary_chr2_coding_uncond_arginine,
                               summary_chr5_coding_uncond_3amino,
                               summary_chr8_coding,
                               summary_chr8_coding_uncond_putrescine,
                               summary_chr11_coding_uncond_arachido,
                               summary_chr11_coding_uncond_3beta,
                               summary_chr12_coding_uncond_1_methyl,
                               summary_chr12_coding_uncond_ethylmalonate,
                               summary_chr16_coding_uncond_propyl,
                               summary_chr16_coding_exp)

tot_coding_uncond_b2_df <- as.data.frame(matrix(unlist(coding_rvs_uncond_sig), ncol = 13))

colnames(tot_coding_uncond_b2_df) <- colnames(coding_rvs_uncond_sig)
tot_coding_uncond_b2_df <-  tot_coding_uncond_b2_df[,c("Gene name", "Chr", "Category",
                                                       "#SNV", "cMAC",
                                                       "STAAR-S(1,25)", "STAAR-S(1,1)",
                                                       "STAAR-B(1,25)", "STAAR-B(1,1)",
                                                       "STAAR-A(1,25)", "STAAR-A(1,1)",
                                                       "ACAT-O", "STAAR-O")]

tot_coding_uncond_b2_df$Metabolite <- c(rep("N2-acetyllysine",nrow(summary_chr2_coding_b2)),
                                        rep("Cysteinylglycine",nrow(summary_chr16_coding)),
                                        rep("N-acetylarginine",nrow(summary_chr2_coding_uncond_arginine)),
                                        rep("3-aminoisobutyrate",nrow(summary_chr5_coding_uncond_3amino)),
                                        rep("N-acetylputrescine",nrow(summary_chr8_coding_uncond_putrescine)),
                                        rep("Arachidonoylcholine",nrow(summary_chr11_coding_uncond_arachido)),
                                        rep("Ethylmalonate",nrow(summary_chr12_coding_uncond_ethylmalonate)),
                                        rep("Propyl 4-hydroxybenzoate sulfate",nrow(summary_chr16_coding_uncond_propyl)),
                                        rep("Cys-gly, oxidized",nrow(summary_chr16_coding_exp)))

tot_coding_uncond_b2_df$Group <- c(rep("Negative Control", nrow(summary_chr2_coding_b2) +
                                         nrow(summary_chr16_coding)),
                                   rep("Test Region", nrow(tot_coding_uncond_b2_df)-
                                         nrow(summary_chr2_coding_b2) -nrow(summary_chr16_coding)))

tot_coding_uncond_b2_df <- unique(tot_coding_uncond_b2_df) |>
  select(c("Gene name", "Chr", "Metabolite", "Group", "Category",
           "#SNV", "cMAC",
           "STAAR-S(1,25)", "STAAR-S(1,1)",
           "STAAR-B(1,25)", "STAAR-B(1,1)",
           "STAAR-A(1,25)", "STAAR-A(1,1)",
           "ACAT-O", "STAAR-O"))

data.table::fwrite(tot_coding_uncond_b2_df, file = "../Data/STAAR_b2_coding_rvs_uncond_sig.csv")


####### CONDITIONAL ANALYSIS ADJUSTING COMMON VARIANTS -- CODING #########
alpha_cond <- 0.001

categories <- c("plof","plof_ds","missense","disruptive_missense",
                "synonymous","ptv","ptv_ds")


summarizeCodingCond <- function(group, metab){
  if(group == "exp_groups"){
    ## output path
    output_path <- paste0("../Data/STAAR_exp_groups/",
                          metab, "/b2/")
  }else{
    output_path <- paste0("../Data/STAAR_neg_controls/",
                          metab, "/b2/")
  }

  res_total_coding <- c()
  for(category in categories){
    load(paste0(output_path, metab, "_Genecentric_coding_cond_b2_", category,
                ".Rdata"))
    if(category == "missense"){
      results_coding <- results_coding[, -c(92:97)]
    }
    #print(dim(results_coding))
    res_total_coding <- rbind(res_total_coding, results_coding)
  }

  res_total_coding <- as.data.frame(res_total_coding)
  res_total_coding_sig <- res_total_coding |>
    filter(`STAAR-B(1,25)` < alpha_cond |  `STAAR-O` < alpha_cond)

  return(res_total_coding_sig)
}


summary_coding_chr2_cond_b2 <- summarizeCodingCond(metab = "x100001721",
                                               group = "neg_controls")

saveRDS(summary_coding_chr2_cond_b2, "../Data/STAAR_neg_controls/x100001721/b2/coding_cond_sig.RDS")

summary_coding_chr5_cond_b2  <- summarizeCodingCond(metab = "x799",
                                                group = "neg_controls")
saveRDS(summary_coding_chr5_cond_b2, "../Data/STAAR_neg_controls/x799/b2/coding_cond_sig.RDS")

summary_coding_chr6_cond_b2  <- summarizeCodingCond(metab = "x1215",
                                                group = "neg_controls")
saveRDS(summary_coding_chr6_cond_b2, "../Data/STAAR_neg_controls/x1215/b2/coding_cond_sig.RDS")

summary_coding_chr16_cond_b2  <- summarizeCodingCond(metab = "x278",
                                                 group = "neg_controls")
saveRDS(as.data.frame(summary_coding_chr16_cond_b2), "../Data/STAAR_neg_controls/x278/b2/coding_cond_sig.RDS")

summary_coding_chr8_cond_b2 <- summarizeCodingCond(group = "exp_groups",
                                               metab = "x1021")
saveRDS(as.data.frame(summary_coding_chr8_cond_b2), "../Data/STAAR_exp_groups/x1021/b2/coding_cond_sig.RDS")

summary_coding_chr10_cond_b2 <- summarizeCodingCond(group ="exp_groups",
                                                metab = "x100000007")
saveRDS(summary_coding_chr10_cond_b2, "../Data/STAAR_exp_groups/x100000007/b2/coding_cond_sig.RDS")

summary_coding_chr13_cond_b2 <- summarizeCodingCond(group ="exp_groups",
                                                metab = "x100004046")
saveRDS(summary_coding_chr13_cond_b2, "../Data/STAAR_exp_groups/x100004046/b2/coding_cond_sig.RDS")

summary_coding_chr16_cond_exp_b2 <- summarizeCodingCond(group = "exp_groups",
                                                    metab = "x1224")
saveRDS(summary_coding_chr16_cond_exp_b2, "../Data/STAAR_exp_groups/x1224/b2/coding_cond_sig.RDS")

# scaled up results

summary_coding_chr11_cond_arachido <- summarizeCodingCond("exp_groups", metab = "x100009332")
saveRDS(summary_coding_chr11_cond_arachido, "../Data/STAAR_exp_groups/x100009332/b2/coding_cond_sig.RDS")

summary_coding_chr11_cond_3beta <- summarizeCodingCond("exp_groups", metab = "x100006370")
saveRDS(summary_coding_chr11_cond_3beta, "../Data/STAAR_exp_groups/x100006370/b2/coding_cond_sig.RDS")

summary_coding_chr12_cond_1_methyl <- summarizeCodingCond("exp_groups", metab = "x100001208")
saveRDS(summary_coding_chr12_cond_1_methyl, "../Data/STAAR_exp_groups/x100001208/b2/coding_cond_sig.RDS")

summary_coding_chr12_cond_ethylmalonate <- summarizeCodingCond("exp_groups", metab = "x2054")
summary_coding_chr12_cond_ethylmalonate <- summary_coding_chr12_cond_ethylmalonate[-3,]
saveRDS(summary_coding_chr12_cond_ethylmalonate, "../Data/STAAR_exp_groups/x2054/b2/coding_cond_sig.RDS")

summary_coding_chr2_cond_arginine  <- summarizeCodingCond("exp_groups", metab = "x100001266")
saveRDS(summary_coding_chr2_cond_arginine, "../Data/STAAR_exp_groups/x100001266/b2/coding_cond_sig.RDS")

summary_coding_chr5_cond_3amino  <- summarizeCodingCond("exp_groups", metab = "x1114")
summary_coding_chr5_cond_3amino <- summary_coding_chr5_cond_3amino[-3,]
saveRDS(summary_coding_chr5_cond_3amino, "../Data/STAAR_exp_groups/x1114/b2/coding_cond_sig.RDS")

summary_coding_chr16_cond_propyl  <- summarizeCodingCond("exp_groups", metab = "x100006264")
saveRDS(summary_coding_chr16_cond_propyl, "../Data/STAAR_exp_groups/x100006264/b2/coding_cond_sig.RDS")

summary_coding_chr8_cond_nacetyl <- summarizeCodingCond("exp_groups", metab = "x192")
saveRDS(summary_coding_chr8_cond_nacetyl, "../Data/STAAR_exp_groups/x192/b2/coding_cond_sig.RDS")


tot_coding_cond_b2 <- rbind(summary_coding_chr2_cond_b2,
        summary_coding_chr16_cond_b2,
                            summary_coding_chr5_cond_3amino,
                            summary_coding_chr8_cond_b2,
                            summary_coding_chr11_cond_arachido,
                            summary_coding_chr12_cond_ethylmalonate,
                            summary_coding_chr16_cond_propyl)
tot_coding_cond_b2_df <- as.data.frame(matrix(unlist(tot_coding_cond_b2), ncol = 91))

colnames(tot_coding_cond_b2_df) <- colnames(tot_coding_cond_b2)
tot_coding_cond_b2_df <-  tot_coding_cond_b2_df[,c("Gene name", "Chr", "Category",
                                         "#SNV", "cMAC",
                                         "STAAR-S(1,25)", "STAAR-S(1,1)",
                                         "STAAR-B(1,25)", "STAAR-B(1,1)",
                                         "STAAR-A(1,25)", "STAAR-A(1,1)",
                                         "ACAT-O", "STAAR-O")]

tot_coding_cond_b2_df$Metabolite <- c(rep("N2-acetyllysine",nrow(summary_coding_chr2_cond_b2)),
                                 rep("Cysteinylglycine",nrow(summary_coding_chr16_cond_b2)),
                                 rep("3-aminoisobutyrate",nrow(summary_coding_chr5_cond_3amino)),
                                 rep("5-oxoproline",nrow(summary_coding_chr8_cond_b2)),
                                 rep("Arachidonoylcholine",nrow(summary_coding_chr11_cond_arachido)),
                                 rep("Ethylmalonate",nrow(summary_coding_chr12_cond_ethylmalonate)),
                                 rep("Propyl 4-hydroxybenzoate sulfate",nrow(summary_coding_chr16_cond_propyl)))

tot_coding_cond_b2_df$Group <- c(rep("Negative Control", nrow(summary_coding_chr2_cond_b2) +
                                          nrow(summary_coding_chr16_cond_b2)),
                                    rep("Test Region", nrow(tot_coding_cond_b2_df)-
                                          nrow(summary_coding_chr2_cond_b2)-nrow(summary_coding_chr16_cond_b2)))

tot_coding_cond_b2_df <- unique(tot_coding_cond_b2_df) |>
  select(c("Gene name", "Chr", "Metabolite", "Group", "Category",
                      "#SNV", "cMAC",
                      "STAAR-S(1,25)", "STAAR-S(1,1)",
                      "STAAR-B(1,25)", "STAAR-B(1,1)",
                      "STAAR-A(1,25)", "STAAR-A(1,1)",
                      "ACAT-O", "STAAR-O"))


alpha <- alpha_uncond_noncoding_b2
alpha_ncrna <- alpha_ncrna_b2

## gene info
gene_centric_noncoding_jobs_num <- 1
source("./20250422_gencentric_noncoding_summary.R")
gene_centric_results_name <- "Genecentric_noncoding_uncond_b2"
summaryNoncoding <- function(agds.path, chr, group, metab, job.num.chr, known_loci = NULL){
  if(group == "exp_groups"){
    ## output path
    output_path <- paste0("../Data/STAAR_exp_groups/",
                          metab, "/b2/")
    load(paste0("../Data/STAAR_exp_groups/",
                metab, "/", metab, "_nullmodel_batch2.Rdata"))
  }else{
    output_path <- paste0("../Data/STAAR_neg_controls/",
                          metab, "/b2/")
    load(paste0("../Data/STAAR_neg_controls/",
                metab, "/", metab, "_nullmodel_batch2.Rdata"))
  }
  input_path <- output_path
  ncRNA_output_path <- output_path
  ncRNA_input_path <- ncRNA_output_path

  noncoding_summary <- Gene_Centric_Noncoding_Results_Summary(agds_dir=agds_dir,gene_centric_noncoding_jobs_num=gene_centric_noncoding_jobs_num,
                                                              input_path=input_path,output_path=output_path,
                                                              gene_centric_results_name=paste0(metab, "_" , gene_centric_results_name),
                                                              ncRNA_jobs_num=ncRNA_jobs_num,ncRNA_input_path=ncRNA_input_path,
                                                              ncRNA_output_path=ncRNA_output_path,ncRNA_results_name=ncRNA_results_name,
                                                              obj_nullmodel=obj_nullmodel,
                                                              method_cond=method_cond, alpha = alpha, alpha_ncRNA = alpha_ncrna,
                                                              QC_label=QC_label,geno_missing_imputation=geno_missing_imputation,variant_type=variant_type,
                                                              Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                              Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                              ncRNA_pos=ncRNA_pos)

  return(noncoding_summary)
}

summary_noncoding_chr2_uncond_b2 <- summaryNoncoding(agds_dir[8], chr = 2,
                                               metab = "x100001721",
                                               group = "neg_controls")

summary_noncoding_chr5_uncond_b2  <- summaryNoncoding(agds_dir[9], chr = 5, metab = "x799",
                                                group = "neg_controls")

summary_noncoding_chr6_uncond_b2  <- summaryNoncoding(agds_dir[10], chr = 6, metab = "x1215",
                                                group = "neg_controls")

summary_noncoding_chr16_uncond_b2  <- summaryNoncoding(agds_dir[14], chr = 16, metab = "x278",
                                                 group = "neg_controls")

summary_noncoding_chr8_uncond_b2 <- summaryNoncoding(agds_dir[11], group = "exp_groups",
                                               metab = "x1021", chr = 8)
summary_noncoding_chr10_uncond_b2 <- summaryNoncoding(agds_dir[12], group ="exp_groups",
                                                metab = "x100000007", chr = 10)
summary_noncoding_chr13_uncond_b2 <- summaryNoncoding(agds_dir[13], group ="exp_groups",
                                                metab = "x100004046", chr = 13)
summary_noncoding_chr16_uncond_exp_b2 <- summaryNoncoding(agds_dir[14], group = "exp_groups",
                                                    metab = "x1224", chr = 16)


summary_chr11_noncoding_uncond_arachido <- summaryNoncoding(agds_dir[1], chr = 11,
                                                      metab = "x100009332",
                                                      group = "exp_groups")

summary_chr11_noncoding_uncond_3beta <- summaryNoncoding(agds_dir[2], chr = 11,
                                                   metab = "x100006370",
                                                   group = "exp_groups")

summary_chr12_noncoding_uncond_1_methyl <- summaryNoncoding(agds_dir[3], chr = 12,
                                                      metab = "x100001208",
                                                      group = "exp_groups")

summary_chr12_noncoding_uncond_ethylmalonate <- summaryNoncoding(agds_dir[4], chr = 12,
                                                           metab = "x2054",
                                                           group = "exp_groups")


summary_chr16_noncoding_uncond_propyl <- summaryNoncoding(agds_dir[5], chr = 16,
                                                    metab = "x100006264",
                                                    group = "exp_groups")

summary_chr2_noncoding_uncond_arginine <- summaryNoncoding(agds_dir[6], chr = 2,
                                                     metab = "x100001266",
                                                     group = "exp_groups")

summary_chr5_noncoding_uncond_3amino <- summaryNoncoding(agds_dir[7], chr = 5,
                                                   metab = "x1114",
                                                   group = "exp_groups")

summary_chr8_noncoding_uncond_putrescine <- summaryNoncoding(agds_dir[8], chr = 8,
                                                       metab = "x192",
                                                       group = "exp_groups")

noncoding_rvs_uncond_sig <- rbind(summary_noncoding_chr2_uncond_b2,
                                  summary_noncoding_chr16_uncond_b2,
                               summary_chr8_noncoding_uncond_putrescine,
                               summary_chr12_noncoding_uncond_ethylmalonate,
                               summary_chr16_noncoding_uncond_propyl,
                               summary_noncoding_chr16_uncond_exp_b2)

tot_noncoding_uncond_b2_df <- as.data.frame(matrix(unlist(noncoding_rvs_uncond_sig), ncol = 91))

colnames(tot_noncoding_uncond_b2_df) <- colnames(noncoding_rvs_uncond_sig)
tot_noncoding_uncond_b2_df <-  tot_noncoding_uncond_b2_df[,c("Gene name", "Chr", "Category",
                                                       "#SNV", "cMAC",
                                                       "STAAR-S(1,25)", "STAAR-S(1,1)",
                                                       "STAAR-B(1,25)", "STAAR-B(1,1)",
                                                       "STAAR-A(1,25)", "STAAR-A(1,1)",
                                                       "ACAT-O", "STAAR-O")]

tot_noncoding_uncond_b2_df$Metabolite <- c(rep("N2-acetyllysine",nrow(summary_noncoding_chr2_uncond_b2)),
                                        rep("Cysteinylglycine",nrow(summary_noncoding_chr16_uncond_b2)),
                                        rep("N-acetylputrescine",nrow(summary_chr8_noncoding_uncond_putrescine)),
                                        rep("Ethylmalonate",nrow(summary_chr12_noncoding_uncond_ethylmalonate)),
                                        rep("Propyl 4-hydroxybenzoate sulfate",nrow(summary_chr16_noncoding_uncond_propyl)),
                                        rep("Cys-gly, oxidized",nrow(summary_noncoding_chr16_uncond_exp_b2)))

tot_noncoding_uncond_b2_df$Group <- c(rep("Negative Control", nrow(summary_noncoding_chr2_uncond_b2) +
                                         nrow(summary_noncoding_chr16_uncond_b2)),
                                   rep("Test Region", nrow(tot_noncoding_uncond_b2_df)-
                                         nrow(summary_noncoding_chr2_uncond_b2) -nrow(summary_noncoding_chr16_uncond_b2)))

tot_noncoding_uncond_b2_df <- unique(tot_noncoding_uncond_b2_df) |>
  select(c("Gene name", "Chr", "Metabolite", "Group", "Category",
           "#SNV", "cMAC",
           "STAAR-S(1,25)", "STAAR-S(1,1)",
           "STAAR-B(1,25)", "STAAR-B(1,1)",
           "STAAR-A(1,25)", "STAAR-A(1,1)",
           "ACAT-O", "STAAR-O"))

data.table::fwrite(tot_noncoding_uncond_b2_df, file = "../Data/STAAR_b2_noncoding_rvs_uncond_sig.csv")

# nonnoncoding region
for(metab in metabs){
  for(cat_func in categories){
    if(group == "exp_groups"){
      ## output path
      output_path <- paste0("../Data/STAAR_exp_groups/",
                            metab, "/b2/")
      res <- getobj(paste0("../Data/STAAR_exp_groups/",
                           metab, "/b2/", cat_func, ".Rdata"))
      ncrna_res <- getobj(paste0("../Data/STAAR_exp_groups/",
                                 metab, "/b2/", ncRNA_results_name, ".Rdata"))
      res <- rbind(res, ncrna_res)
    }else{
      output_path <- paste0("../Data/STAAR_neg_controls/",
                            metab, "/b2/")
      res <- getobj(paste0("../Data/STAAR_neg_controls/",
                           metab, "/b2/", cat_func, ".Rdata"))
    }
    if(is.null(res)){
      num_tests <- num_tests
    }else if(!is.null(res) & is.null(dim(res))){
      num_tests <- num_tests + 1
    } else{
      num_tests <- num_tests + dim(res)[1]
    }
  }
}

for(metab in metabs){
  for(cat_func in categories){
    if(group == "exp_groups"){
      ## output path
      output_path <- paste0("../Data/STAAR_exp_groups/",
                            metab, "/b2/")
      res <- getobj(paste0("../Data/STAAR_exp_groups/",
                           metab, "/b2/", cat_func, ".Rdata"))
      ncrna_res <- getobj(paste0("../Data/STAAR_exp_groups/",
                                 metab, "/b2/", ncRNA_results_name, ".Rdata"))
      res <- rbind(res, ncrna_res)
    }else{
      output_path <- paste0("../Data/STAAR_neg_controls/",
                            metab, "/b2/")
      res <- getobj(paste0("../Data/STAAR_neg_controls/",
                           metab, "/b2/", cat_func, ".Rdata"))
      print(res)
    }
    if(is.null(res)){
      num_tests <- num_tests
    }else if(!is.null(res) & is.null(dim(res))){
      num_tests <- num_tests + 1
    } else{
      num_tests <- num_tests + dim(res)[1]
    }
  }
}

# 3630 tests in total
alpha <- 0.05/3630


####### CONDITIONAL ANALYSIS ADJUSTING COMMON VARIANTS -- NONCODING #########

summarizeNoncodingCond <- function(group, metab, ncRNA_cond_res = NULL){
  if(group == "exp_groups"){
    ## output path
    output_path <- paste0("../Data/STAAR_exp_groups/",
                          metab, "/b2/")
  }else{
    output_path <- paste0("../Data/STAAR_neg_controls/",
                          metab, "/b2/")
  }

  res_total_noncoding <- c()
  for(category in categories){
    load(paste0(output_path, metab, "_Genecentric_noncoding_cond_b2_", category,
                ".Rdata"))
    res_total_noncoding <- rbind(res_total_noncoding, results_noncoding)
  }

  res_total_noncoding <- toDataframe(res_total_noncoding)
  print(dim(res_total_noncoding))
  res_total_noncoding <- res_total_noncoding |>
    filter((`STAAR-B(1,25)` < alpha_cond) |  `STAAR-O` < alpha_cond)
  print(dim(res_total_noncoding))
  if(! is.null(ncRNA_cond_res)){
    print("not null ncRNA")
    ncRNA_cond_res <- as.data.frame(ncRNA_cond_res)
    ncRNA_cond_res <- ncRNA_cond_res |>
      filter(( `STAAR-B(1,25)` < alpha_cond) | `STAAR-O` < alpha_cond)
    res_total_noncoding <- rbind(res_total_noncoding, ncRNA_cond_res)

  }
  save(res_total_noncoding,
       file = paste0(output_path, "noncoding_cond_sig.Rdata"))
  return(res_total_noncoding)
}

toDataframe <- function(dat){
  dat_df <- as.data.frame(matrix(unlist(dat), ncol = 91))
  colnames(dat_df) <- colnames(dat)
  return(dat_df)
}


ncRNA_chr2_cond <- get(load("../Data/STAAR_neg_controls/x100001721/b2/Genetic_Noncoding_ncRNA_cond_b2.Rdata"))
summary_chr2_noncoding_cond <- summarizeNoncodingCond("neg_controls",
                                                           metab = "x100001721",
                                                           ncRNA_cond_res = ncRNA_chr2_cond)


ncRNA_chr5_cond <- get(load("../Data/STAAR_neg_controls/x799/b2/Genetic_Noncoding_ncRNA_cond_b2.Rdata"))
summary_chr5_noncoding_cond <- summarizeNoncodingCond("neg_controls",
                                                      metab = "x799",
                                                      ncRNA_cond_res = ncRNA_chr5_cond)

ncRNA_chr6_cond <- get(load("../Data/STAAR_neg_controls/x1215/b2/Genetic_Noncoding_ncRNA_cond_b2.Rdata"))
summary_chr6_noncoding_cond <- summarizeNoncodingCond("neg_controls",
                                                      metab = "x1215",
                                                      ncRNA_cond_res = ncRNA_chr6_cond)


ncRNA_chr16_cond <- get(load("../Data/STAAR_neg_controls/x278/b2/Genetic_Noncoding_ncRNA_cond_b2.Rdata"))
summary_chr16_noncoding_cond <- summarizeNoncodingCond("neg_controls",
                                                      metab = "x278",
                                                      ncRNA_cond_res = ncRNA_chr16_cond)

ncRNA_chr8_cond <- get(load("../Data/STAAR_exp_groups/x1021/b2/Genetic_Noncoding_ncRNA_cond_b2.Rdata"))
ncRNA_chr8_cond_df <- as.data.frame(matrix(unlist(ncRNA_chr8_cond), ncol = 91))
colnames(ncRNA_chr8_cond_df) <- colnames(ncRNA_chr8_cond)
summary_chr8_noncoding_cond <- summarizeNoncodingCond("exp_groups",
                                                      metab = "x1021",
                                                      ncRNA_cond_res = ncRNA_chr8_cond_df)
ncRNA_chr10_cond <- get(load("../Data/STAAR_exp_groups/x100000007/b2/Genetic_Noncoding_ncRNA_cond_b2.Rdata"))
ncRNA_chr10_cond_df <- as.data.frame(matrix(unlist(ncRNA_chr10_cond), ncol = 91))
colnames(ncRNA_chr10_cond_df) <- colnames(ncRNA_chr10_cond)
summary_chr10_noncoding_cond <- summarizeNoncodingCond("exp_groups",
                                                      metab = "x100000007",
                                                      ncRNA_cond_res = ncRNA_chr10_cond_df)
ncRNA_chr13_cond <- get(load("../Data/STAAR_exp_groups/x100004046/b2/Genetic_Noncoding_ncRNA_cond_b2.Rdata"))
ncRNA_chr13_cond_df <- as.data.frame(matrix(unlist(ncRNA_chr13_cond), ncol = 91))
colnames(ncRNA_chr13_cond_df) <- colnames(ncRNA_chr13_cond)
summary_noncoding_chr13_cond_b2 <- summarizeNoncodingCond("exp_groups",
                                                            metab = "x100004046",
                                                            ncRNA_cond_res = ncRNA_chr13_cond_df)
ncRNA_chr16_exp_cond <- get(load("../Data/STAAR_exp_groups/x1224/b2/Genetic_Noncoding_ncRNA_cond_b2.Rdata"))
ncRNA_chr16_exp_cond_df <- as.data.frame(matrix(unlist(ncRNA_chr16_exp_cond), ncol = 91))
colnames(ncRNA_chr16_exp_cond_df) <- colnames(ncRNA_chr16_exp_cond)
summary_noncoding_chr16_cond_exp_b2 <- summarizeNoncodingCond("exp_groups",
                                                                metab = "x1224",
                                                                ncRNA_cond_res = ncRNA_chr16_exp_cond_df)


#### scaled up ones ####

ncRNA_chr2_cond_acetylarginine <- get(load("../Data/STAAR_exp_groups/x100001266/b2/Genetic_Noncoding_ncRNA_cond_b2.Rdata"))
summary_chr2_noncoding_acetylarginine_cond <- summarizeNoncodingCond("exp_groups",
                                                                     metab = "x100001266",
                                                                     ncRNA_cond_res = ncRNA_chr2_cond_acetylarginine)

ncRNA_chr5_cond_x1114 <- get(load("../Data/STAAR_exp_groups/x1114/b2/Genetic_Noncoding_ncRNA_cond_b2.Rdata"))
summary_chr5_noncoding_x1114_cond <- summarizeNoncodingCond("exp_groups",
                                                      metab = "x1114",
                                                      ncRNA_cond_res = ncRNA_chr5_cond_x1114)

ncRNA_chr8_x192_cond <- get(load("../Data/STAAR_exp_groups/x192/b2/Genetic_Noncoding_ncRNA_cond_b2.Rdata"))
summary_chr8_noncoding_x192_cond <- summarizeNoncodingCond("exp_groups",
                                                      metab = "x192",
                                                      ncRNA_cond_res = ncRNA_chr8_x192_cond)

summary_chr8_noncoding_x192_cond <- summarizeNoncodingCond("exp_groups",
                                                           metab = "x192",
                                                           ncRNA_cond_res = ncRNA_chr8_x192_cond)

ncRNA_chr11_arachido_cond <- get(load("../Data/STAAR_exp_groups/x100009332/b2/Genetic_Noncoding_ncRNA_cond_b2.Rdata"))
summary_chr11_noncoding_arachido_cond <- summarizeNoncodingCond("exp_groups",
                                                           metab = "x100009332",
                                                           ncRNA_cond_res = ncRNA_chr11_arachido_cond)

ncRNA_chr11_3beta_cond <- get(load("../Data/STAAR_exp_groups/x100006370/b2/Genetic_Noncoding_ncRNA_cond_b2.Rdata"))
summary_chr11_noncoding_3beta_cond <- summarizeNoncodingCond("exp_groups",
                                                           metab = "x100006370",
                                                           ncRNA_cond_res = ncRNA_chr11_3beta_cond)

ncRNA_chr12_cond_1_methyl <- get(load("../Data/STAAR_exp_groups/x100001208/b2/Genetic_Noncoding_ncRNA_cond_b2.Rdata"))
summary_chr12_noncoding_1_methyl <- summarizeNoncodingCond("exp_groups",
                                                             metab = "x100001208",
                                                             ncRNA_cond_res = ncRNA_chr12_cond_1_methyl)

ncRNA_chr12_cond_1_ethylmalonate <- get(load("../Data/STAAR_exp_groups/x2054/b2/Genetic_Noncoding_ncRNA_cond_b2.Rdata"))
summary_chr12_noncoding_ethylmalonate <- summarizeNoncodingCond("exp_groups",
                                                           metab = "x2054",
                                                           ncRNA_cond_res = ncRNA_chr12_cond_1_ethylmalonate)

ncRNA_chr16_cond_x6264 <- get(load("../Data/STAAR_exp_groups/x100006264/b2/Genetic_Noncoding_ncRNA_cond_b2.Rdata"))
summary_chr16_noncoding_cond_x6264 <- summarizeNoncodingCond("exp_groups",
                                                             metab = "x100006264",
                                                             ncRNA_cond_res = ncRNA_chr16_cond_x6264)


# all non-coding rv sets
tot_noncoding_cond_b2 <- rbind(summary_chr16_noncoding_cond,
                               summary_chr8_noncoding_cond,
                               summary_chr11_noncoding_3beta_cond,
                            summary_chr16_noncoding_cond_x6264,
                            summary_noncoding_chr16_cond_exp_b2)
tot_noncoding_cond_b2_df <- as.data.frame(matrix(unlist(tot_noncoding_cond_b2), ncol = 91))

colnames(tot_noncoding_cond_b2_df) <- colnames(tot_noncoding_cond_b2)
tot_noncoding_cond_b2_df <-  tot_noncoding_cond_b2_df[,c("Gene name", "Chr", "Category",
                                                   "#SNV", "cMAC",
                                                   "STAAR-S(1,25)", "STAAR-S(1,1)",
                                                   "STAAR-B(1,25)", "STAAR-B(1,1)",
                                                   "STAAR-A(1,25)", "STAAR-A(1,1)",
                                                   "ACAT-O", "STAAR-O") ]

tot_noncoding_cond_b2_df$Metabolite <- c(rep("Cysteinylglycine",nrow(summary_chr16_noncoding_cond)),
                                 rep("5-oxoproline",nrow(summary_chr8_noncoding_cond)),
                                 rep("3beta-hydroxy-5-cholestenoate",nrow(summary_chr11_noncoding_3beta_cond)),
                                 rep("Propyl 4-hydroxybenzoate sulfate",nrow(summary_chr16_noncoding_cond_x6264)),
                                 rep("Cys-gly, oxidized",nrow(summary_noncoding_chr16_cond_exp_b2)))

tot_noncoding_cond_b2_df$Group <- c(rep("Negative Control", nrow(summary_chr16_noncoding_cond)),
                                    rep("Test Region", nrow(tot_noncoding_cond_b2_df)-
                                          nrow(summary_chr16_noncoding_cond)))


tot_noncoding_cond_b2_df <- unique(tot_noncoding_cond_b2_df) |>
  select(c("Gene name", "Chr", "Metabolite", "Group", "Category",
           "#SNV", "cMAC",
           "STAAR-S(1,25)", "STAAR-S(1,1)",
           "STAAR-B(1,25)", "STAAR-B(1,1)",
           "STAAR-A(1,25)", "STAAR-A(1,1)",
           "ACAT-O", "STAAR-O"))

tot_noncoding_cond_b2_df <- unique(tot_noncoding_cond_b2_df)

# comobind all identified non-coding and coding rv sets
all_rvs_cond_sig_b2 <- rbind(tot_coding_cond_b2_df, tot_noncoding_cond_b2_df)
all_rvs_cond_sig_b2[, c("#SNV")] <- as.integer(all_rvs_cond_sig_b2[, c("#SNV")])

all_rvs_cond_sig_b2 <- data.table(all_rvs_cond_sig_b2)
all_rvs_cond_sig_b2[, 6:14 := lapply(.SD, as.numeric), .SDcols = 6:14]

fwrite(all_rvs_cond_sig_b2,  "../Data/STAAR_b2_all_rvs_cond_sig.csv")


#------------------------------------------------------------------------------
# LEGACY / EARLIER VERSION: v2_scaled_up  (37 unique blocks)
#------------------------------------------------------------------------------
# Kept verbatim. These are blocks that do NOT appear in the current version above
# (mostly hard-coded per-metabolite / per-chromosome run calls and older path setups).

agds_dir <- get(load("../Data/STAAR_prep/agds_dir_scale_up_final_b2.Rdata"))

summarizeCodingCond <- function(group, metab){
  if(group == "exp_groups"){
    ## output path
    output_path <- paste0("../Data/STAAR_exp_groups/",
                          metab, "/b2/")
  }else{
    output_path <- paste0("../Data/STAAR_neg_controls/",
                          metab, "/b2/")
  }

  res_total_coding <- c()
  for(category in categories){
    load(paste0(output_path, metab, "_Genecentric_coding_cond_b2_", category,
                ".Rdata"))
    if(category == "missense"){
      results_coding <- results_coding[, -c(92:97)]
    }
    #print(dim(results_coding))
    res_total_coding <- rbind(res_total_coding, results_coding)
  }

  res_total_coding <- as.data.frame(res_total_coding)
  print(dim(res_total_coding))
  res_total_coding_sig <- res_total_coding |>
    filter((`STAAR-B(1,25)` < alpha_cond) &  `STAAR-O` < alpha_cond)

  return(res_total_coding_sig)
}

summary_chr11_coding_cond_arachido <- summarizeCodingCond("exp_groups", metab = "x100009332")
saveRDS(summary_chr11_coding_cond_arachido, "../Data/STAAR_exp_groups/x100009332/b2/coding_cond_sig.RDS")

#summary_chr11_coding_cond_3beta <- summarizeCodingCond("exp_groups", metab = "x100006370")

summary_chr12_coding_cond_1_methyl <- summarizeCodingCond("exp_groups", metab = "x100001208")
saveRDS(summary_chr12_coding_cond_1_methyl, "../Data/STAAR_exp_groups/x100001208/b2/coding_cond_sig.RDS")


summary_chr12_coding_cond_ethylmalonate <- summarizeCodingCond("exp_groups", metab = "x2054")
summary_chr12_coding_cond_ethylmalonate <- summary_chr12_coding_cond_ethylmalonate[-3,]
saveRDS(summary_chr12_coding_cond_ethylmalonate, "../Data/STAAR_exp_groups/x2054/b2/coding_cond_sig.RDS")


summary_chr2_coding_cond_arginine  <- summarizeCodingCond("exp_groups", metab = "x100001266")
saveRDS(summary_chr2_coding_cond_arginine, "../Data/STAAR_exp_groups/x100001266/b2/coding_cond_sig.RDS")

summary_chr5_coding_cond_3amino  <- summarizeCodingCond("exp_groups", metab = "x1114")
summary_chr5_coding_cond_3amino <- summary_chr5_coding_cond_3amino[-3,]
saveRDS(summary_chr5_coding_cond_3amino, "../Data/STAAR_exp_groups/x1114/b2/coding_cond_sig.RDS")

summary_chr16_coding_cond  <- summarizeCodingCond("exp_groups", metab = "x100006264")
saveRDS(summary_chr16_coding_cond, "../Data/STAAR_exp_groups/x100006264/b2/coding_cond_sig.RDS")

summary_chr8_coding_cond <- summarizeCodingCond("exp_groups", metab = "x192")

#1.37741e-05
summaryNonCoding <- function(agds.path, chr, group, metab, job.num.chr, known_loci = NULL){
  if(group == "exp_groups"){
    ## output path
    output_path <- paste0("../Data/STAAR_exp_groups/",
                          metab, "/b2/")
    load(paste0("../Data/STAAR_exp_groups/",
                metab, "/", metab, "_nullmodel_batch2.Rdata"))
  }else{
    output_path <- paste0("../Data/STAAR_neg_controls/",
                          metab, "/b2/")
    load(paste0("../Data/STAAR_neg_controls/",
                metab, "/", metab, "_nullmodel_batch2.Rdata"))
  }
  input_path <- output_path
  ncRNA_output_path <- output_path
  ncRNA_input_path <- ncRNA_output_path

  noncoding_summary <- Gene_Centric_Noncoding_Results_Summary(agds_dir=agds_dir,gene_centric_noncoding_jobs_num=gene_centric_noncoding_jobs_num,
                                                              input_path=input_path,output_path=output_path,
                                                              gene_centric_results_name=paste0(metab, "_" , gene_centric_results_name),
                                                              ncRNA_jobs_num=ncRNA_jobs_num,ncRNA_input_path=ncRNA_input_path,
                                                              ncRNA_output_path=ncRNA_output_path,ncRNA_results_name=ncRNA_results_name,
                                                              obj_nullmodel=obj_nullmodel,
                                                              method_cond=method_cond, alpha = alpha, alpha_ncRNA = 1e-4,
                                                              QC_label=QC_label,geno_missing_imputation=geno_missing_imputation,variant_type=variant_type,
                                                              Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                              Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                              ncRNA_pos=ncRNA_pos)

  return(noncoding_summary)
}

summary_chr11_noncoding_uncond_arachido <- summaryNonCoding(agds_dir[1], chr = 11,
                                                            metab = "x100009332",
                                                            group = "exp_groups")

summary_chr11_noncoding_uncond_3beta <- summaryNonCoding(agds_dir[2], chr = 11,
                                                         metab = "x100006370",
                                                         group = "exp_groups")

summary_chr12_noncoding_uncond_1_methyl <- summaryNonCoding(agds_dir[3], chr = 12,
                                                            metab = "x100001208",
                                                            group = "exp_groups")


summary_chr12_noncoding_uncond_ethylmalonate <- summaryNonCoding(agds_dir[4], chr = 12,
                                                                 metab = "x2054",
                                                                 group = "exp_groups")

summary_chr16_noncoding_uncond_propyl <- summaryNonCoding(agds_dir[5], chr = 16,
                                                          metab = "x100006264",
                                                          group = "exp_groups")

summary_chr2_noncoding_uncond_arginine <- summaryNonCoding(agds_dir[6], chr = 2,
                                                           metab = "x100001266",
                                                           group = "exp_groups")

summary_chr5_noncoding_uncond_3amino <- summaryNonCoding(agds_dir[7], chr = 5,
                                                         metab = "x1114",
                                                         group = "exp_groups")

summary_chr8_noncoding_uncond_putrescine <- summaryNonCoding(agds_dir[8], chr = 8,
                                                             metab = "x192",
                                                             group = "exp_groups")


summarizeNoncodingCond <- function(group, metab, ncRNA_cond_res = NULL){
  if(group == "exp_groups"){
    ## output path
    output_path <- paste0("../Data/STAAR_exp_groups/",
                          metab, "/b2/")
  }else{
    output_path <- paste0("../Data/STAAR_neg_controls/",
                          metab, "/b2/")
  }

  res_total_noncoding <- c()
  for(category in categories){
    load(paste0(output_path, metab, "_Genecentric_noncoding_cond_b2_", category,
                ".Rdata"))
    res_total_noncoding <- rbind(res_total_noncoding, results_noncoding)
  }

  res_total_noncoding <- as.data.frame(res_total_noncoding)
  print(dim(res_total_noncoding))
  res_total_noncoding <- res_total_noncoding |>
    filter((`STAAR-B(1,25)` < alpha_cond) |  `STAAR-O` < alpha_cond)
  if(! is.null(ncRNA_cond_res)){
    ncRNA_cond_res <- as.data.frame(ncRNA_cond_res)
    ncRNA_cond_res <- ncRNA_cond_res |>
      filter(( `STAAR-B(1,25)` < alpha_cond) | `STAAR-O` < alpha_cond)
    res_total_noncoding <- rbind(res_total_noncoding, ncRNA_cond_res)

  }
  # save(res_total_noncoding,
  #      file = paste0(output_path, "Genecentric_noncoding_total.Rdata"))
  return(res_total_noncoding)
}

ncRNA_chr2_x1266_cond <- get(load("../Data/STAAR_exp_groups/x100001266/b2/Genetic_Noncoding_ncRNA_cond_b2.Rdata"))
summary_chr2_x1266_noncoding_cond <- summarizeNoncodingCond("exp_groups",
                                                      metab = "x100001266",
                                                      ncRNA_cond_res = ncRNA_chr2_x1266_cond)

ncRNA_chr5_cond <- get(load("../Data/STAAR_exp_groups/x1114/b2/Genetic_Noncoding_ncRNA_cond_b2.Rdata"))
summary_chr5_noncoding_cond <- summarizeNoncodingCond("exp_groups",
                                                      metab = "x1114",
                                                      ncRNA_cond_res = ncRNA_chr5_cond)

ncRNA_chr8_cond <- get(load("../Data/STAAR_exp_groups/x192/b2/Genetic_Noncoding_ncRNA_cond_b2.Rdata"))
summary_chr8_noncoding_cond <- summarizeNoncodingCond("exp_groups",
                                                      metab = "x192",
                                                      ncRNA_cond_res = ncRNA_chr8_cond)

ncRNA_chr16_cond <- get(load("../Data/STAAR_exp_groups/x100006264/Genetic_Noncoding_ncRNA_cond.Rdata"))
summary_chr16_noncoding_cond <- summarizeNoncodingCond("exp_groups",
                                                       metab = "x100006264",
                                                       ncRNA_cond_res = ncRNA_chr16_cond)
