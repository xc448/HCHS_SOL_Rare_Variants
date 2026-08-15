#==============================================================================
# STAAR results summary - BATCH 1
#==============================================================================
# MERGED FROM (chronological):
#   - 20250126_STAAR_pipeline_summary.R   ->  section "v1_original"
#   - 20250720_STAAR_summary_scaled_up.R   ->  section "v2_scaled_up"
#   - 20251211_STAAR_summary_b1.R   ->  section "v3_final_b1"
# NOTE: Defines summaryIndividual / pruneVariants / summaryCoding / summarizeCodingCond /
# NOTE: summaryNonCoding / summarizeNoncodingCond. v3 is the version used for the manuscript.
#==============================================================================

#==============================================================================
# MAIN PIPELINE  (current version: v3_final_b1)
#==============================================================================

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)
library(tidyverse)


###### UNCONDITIONAL ANALYSIS ########
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
    ## output path
    output_path <- paste0("../Data/STAAR_exp_groups/",
                          metab, "/")
    res <- getobj(paste0("../Data/STAAR_exp_groups/",
                         metab, "/", cat_func, ".Rdata"))

    if(cat_func == "missense"){
      ncol <- length(res[1,])
      print(ncol)
      res <- as.data.frame(matrix(unlist(res), ncol = ncol))
      res <- unique(res)
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
                            metab, "/")
      res <- getobj(paste0("../Data/STAAR_neg_controls/",
                           metab, "/", cat_func, ".Rdata"))
      #print(res)

    if(is.null(res)){
      num_tests <- num_tests
    }else if(!is.null(res) & is.null(dim(res))){
      num_tests <- num_tests + 1
    } else{
      if(cat_func == "missense"){
        ncol <- length(res[1,])
        print(ncol)
        res <- as.data.frame(matrix(unlist(res), ncol = ncol))
        res <- unique(res)
      }
      num_tests <- num_tests + dim(res)[1]
    }
  }
}
print(num_tests)
alpha_uncond_coding <- 0.05/num_tests


###### NONCODING #####
ncRNA_results_name <- "Genetic_Noncoding_ncRNA"
categories=c("downstream","upstream","UTR","promoter_CAGE",
             "promoter_DHS","enhancer_CAGE","enhancer_DHS")
num_tests_noncoding <- 0
for(metab in metabs){
  for(cat_func in categories){
    if(group == "exp_groups"){
      ## output path
      output_path <- paste0("../Data/STAAR_exp_groups/",
                            metab, "/")
      res <- getobj(paste0("../Data/STAAR_exp_groups/",
                           metab, "/", cat_func, ".Rdata"))
    }else{
      output_path <- paste0("../Data/STAAR_neg_controls/",
                            metab, "/")
      res <- getobj(paste0("../Data/STAAR_neg_controls/",
                           metab, "/", cat_func, ".Rdata"))
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
                          metab, "/")
    ncrna_res <- getobj(paste0("../Data/STAAR_exp_groups/",
                               metab, "/", ncRNA_results_name, ".Rdata"))
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
                            metab, "/")
      res <- getobj(paste0("../Data/STAAR_exp_groups/",
                           metab, "/", cat_func, ".Rdata"))
    }else{
      output_path <- paste0("../Data/STAAR_neg_controls/",
                            metab, "/")
      res <- getobj(paste0("../Data/STAAR_neg_controls/",
                           metab, "/", cat_func, ".Rdata"))
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
                          metab, "/")
    ncrna_res <- getobj(paste0("../Data/STAAR_neg_controls/",
                               metab, "/", ncRNA_results_name, ".Rdata"))
    if(is.null(ncrna_res)){
      num_tests_ncrna <- num_tests_ncrna
    }else if(!is.null(ncrna_res) & is.null(dim(ncrna_res))){
      num_tests_ncrna <- num_tests_ncrna + 1
    } else{
      num_tests_ncrna <- num_tests_ncrna + dim(ncrna_res)[1]
    }
  }
}


# ncrna number of tests
alpha_ncrna <- 0.05/num_tests_ncrna

# 3063 tests in total for noncoding regions
# therefore the threshold for all non-coding rv sets are
alpha <- 0.05/num_tests_noncoding


####### CONDITIONAL ANALYSIS ADJUSTING COMMON VARIANTS -- CODING #########

alpha_cond <- 0.001

categories <- c("plof","plof_ds","missense","disruptive_missense",
                "synonymous","ptv","ptv_ds")


summarizeCodingCond <- function(group, metab){
  if(group == "exp_groups"){
    ## output path
    output_path <- paste0("../Data/STAAR_exp_groups/",
                          metab, "/")
  }else{
    output_path <- paste0("../Data/STAAR_neg_controls/",
                          metab, "/")
  }

  res_total_coding <- c()
  for(category in categories){
    load(paste0(output_path, metab, "_Genecentric_coding_cond_", category,
                ".Rdata"))
    if(category == "missense"){
      results_coding <- results_coding[, -c(92:97)]
    }
    #print(dim(results_coding))
    res_total_coding <- rbind(res_total_coding, results_coding)
  }

  res_total_coding <- as.data.frame(res_total_coding)

  res_total_coding_sig <- res_total_coding_sig <- res_total_coding |>
    filter((`STAAR-B(1,25)` < alpha_cond) |  `STAAR-O` < alpha_cond)

  return(res_total_coding_sig)
}

toDataframe <- function(dat){
  dat_df <- as.data.frame(matrix(unlist(dat), ncol = 91))
  colnames(dat_df) <- colnames(dat)
  return(dat_df)
}

summary_chr2_coding_cond <- summarizeCodingCond("neg_controls", metab = "x100001721")
saveRDS(summary_chr2_coding_cond, "../Data/STAAR_neg_controls/x100001721/coding_sig_cond.RDS")

summary_chr5_coding_cond <- summarizeCodingCond("neg_controls", metab = "x799")
saveRDS(summary_chr5_coding_cond, "../Data/STAAR_neg_controls/x799/coding_sig_cond.RDS")

summary_chr6_coding_cond <- summarizeCodingCond("neg_controls", metab = "x1215")
saveRDS(summary_chr6_coding_cond, "../Data/STAAR_neg_controls/x1215/coding_sig_cond.RDS")

summary_chr16_coding_cond <- summarizeCodingCond("neg_controls", metab = "x278")
saveRDS(summary_chr16_coding_cond, "../Data/STAAR_neg_controls/x278/coding_sig_cond.RDS")

summary_chr8_coding_cond <- summarizeCodingCond("exp_groups", metab = "x1021")
saveRDS(summary_chr8_coding_cond, "../Data/STAAR_exp_groups/x1021/coding_sig_cond.RDS")

summary_chr10_coding_cond <- summarizeCodingCond("exp_groups", metab = "x100000007")
saveRDS(summary_chr10_coding_cond, "../Data/STAAR_exp_groups/x100000007/coding_sig_cond.RDS")

summary_chr13_coding_cond <- summarizeCodingCond("exp_groups", metab = "x100004046")
saveRDS(summary_chr13_coding_cond, "../Data/STAAR_exp_groups/x100004046/coding_sig_cond.RDS")

summary_chr16_coding_cond_exp <- summarizeCodingCond("exp_groups", metab = "x1224")
saveRDS(summary_chr16_coding_cond_exp, "../Data/STAAR_exp_groups/x1224/coding_sig_cond.RDS")

# scaled up

summary_chr11_coding_cond_arachido <- summarizeCodingCond("exp_groups", metab = "x100009332")
saveRDS(summary_chr11_coding_cond_arachido, "../Data/STAAR_exp_groups/x100009332/coding_cond_sig.RDS")

# no significant conditional ones
# summary_chr11_coding_cond_3beta <- summarizeCodingCond("exp_groups", metab = "x100006370")

summary_chr12_coding_cond_1_methyl <- summarizeCodingCond("exp_groups", metab = "x100001208")
summary_chr12_coding_cond_1_methyl <- summary_chr12_coding_cond_1_methyl[-3,]
saveRDS(summary_chr12_coding_cond_1_methyl, "../Data/STAAR_exp_groups/x100001208/coding_cond_sig.RDS")

summary_chr12_coding_cond_ethylmalonate <- summarizeCodingCond("exp_groups", metab = "x2054")
summary_chr12_coding_cond_ethylmalonate <- summary_chr12_coding_cond_ethylmalonate[-3,]
saveRDS(summary_chr12_coding_cond_ethylmalonate, "../Data/STAAR_exp_groups/x2054/coding_cond_sig.RDS")

summary_chr2_coding_cond_arginine  <- summarizeCodingCond("exp_groups", metab = "x100001266")
saveRDS(summary_chr2_coding_cond_arginine, "../Data/STAAR_exp_groups/x100001266/coding_cond_sig.RDS")

summary_chr5_coding_cond_3amino  <- summarizeCodingCond("exp_groups", metab = "x1114")
summary_chr5_coding_cond_3amino <- summary_chr5_coding_cond_3amino[-3,]
saveRDS(summary_chr5_coding_cond_3amino, "../Data/STAAR_exp_groups/x1114/coding_cond_sig.RDS")

# summary_chr8_coding_cond_x192  <- summarizeCodingCond("exp_groups", metab = "x192")

summary_chr16_coding_cond  <- summarizeCodingCond("exp_groups", metab = "x100006264")
saveRDS(summary_chr16_coding_cond, "../Data/STAAR_exp_groups/x100006264/coding_cond_sig.RDS")


tot_coding_cond_b1 <- rbind(summary_chr2_coding_cond,
                            summary_chr5_coding_cond_3amino,
                            summary_chr8_coding_cond,
                            summary_chr11_coding_cond_arachido,
                            summary_chr12_coding_cond_1_methyl,
                            summary_chr12_coding_cond_ethylmalonate,
                            summary_chr16_coding_cond_exp)

tot_coding_cond_b1_df <- as.data.frame(matrix(unlist(tot_coding_cond_b1), ncol = 91))

colnames(tot_coding_cond_b1_df) <- colnames(tot_coding_cond_b1)
tot_coding_cond_b1_df <-  tot_coding_cond_b1_df[,c("Gene name", "Chr", "Category",
                                                   "#SNV", "cMAC",
                                                   "STAAR-S(1,25)", "STAAR-S(1,1)",
                                                   "STAAR-B(1,25)", "STAAR-B(1,1)",
                                                   "STAAR-A(1,25)", "STAAR-A(1,1)",
                                                   "ACAT-O", "STAAR-O")]

tot_coding_cond_b1_df$Metabolite <- c(rep("N2-acetyllysine",nrow(summary_chr2_coding_cond)),
                                      rep("3-aminoisobutyrate",nrow(summary_chr5_coding_cond_3amino)),
                                      rep("5-oxoproline",nrow(summary_chr8_coding_cond)),
                                      rep("Arachidonoylcholine",nrow(summary_chr11_coding_cond_arachido)),
                                      rep("1-methylimidazoleacetate",nrow(summary_chr12_coding_cond_1_methyl)),
                                      rep("Ethylmalonate",nrow(summary_chr12_coding_cond_ethylmalonate)),
                                      rep("Cys-gly, oxidized",nrow(summary_chr16_coding_cond_exp)))

tot_coding_cond_b1_df$Group <- c(rep("Negative Control", nrow(summary_chr2_coding_cond)),
                                 rep("Test Region", nrow(tot_coding_cond_b1_df)-
                                       nrow(summary_chr2_coding_cond)))

tot_coding_cond_b1_df <- unique(tot_coding_cond_b1_df) |>
  select(c("Gene name", "Chr", "Metabolite", "Group", "Category",
           "#SNV", "cMAC",
           "STAAR-S(1,25)", "STAAR-S(1,1)",
           "STAAR-B(1,25)", "STAAR-B(1,1)",
           "STAAR-A(1,25)", "STAAR-A(1,1)",
           "ACAT-O", "STAAR-O"))

####### CONDITIONAL ANALYSIS ADJUSTING COMMON VARIANTS -- NONCODING #########
summarizeNoncodingCond <- function(group, metab, ncRNA_cond_res = NULL){
  if(group == "exp_groups"){
    ## output path
    output_path <- paste0("../Data/STAAR_exp_groups/",
                          metab, "/")
  }else{
    output_path <- paste0("../Data/STAAR_neg_controls/",
                          metab, "/")
  }

  res_total_noncoding <- c()
  for(category in categories){
    load(paste0(output_path, metab, "_Genecentric_noncoding_cond_", category,
                ".Rdata"))
    res_total_noncoding <- rbind(res_total_noncoding, results_noncoding)
  }

  res_total_noncoding <- as.data.frame(res_total_noncoding)
  #print(res_total_noncoding)
  res_total_noncoding <- res_total_noncoding |>
    filter((`STAAR-B(1,25)` < alpha_cond) |  `STAAR-O` < alpha_cond)
  if(! is.null(ncRNA_cond_res)){
    ncRNA_cond_res <- as.data.frame(ncRNA_cond_res)
    ncRNA_cond_res <- ncRNA_cond_res |>
      filter((`STAAR-B(1,25)` < alpha_cond) | `STAAR-O` < alpha_cond)
    res_total_noncoding <- rbind(res_total_noncoding, ncRNA_cond_res)

  }
  # save(res_total_noncoding,
  #      file = paste0(output_path, "Genecentric_noncoding_total.Rdata"))
  return(res_total_noncoding)
}


load("../Data/STAAR_neg_controls/x100001721/Genetic_Noncoding_ncRNA_cond.Rdata")
summary_chr2_noncoding_cond <- summarizeNoncodingCond("neg_controls", metab = "x100001721",
                                                      results_ncRNA_cond)

load("../Data/STAAR_neg_controls/x799/Genetic_Noncoding_ncRNA_cond.Rdata")
summary_chr5_noncoding_cond <- summarizeNoncodingCond("neg_controls", metab = "x799",
                                                      results_ncRNA_cond)

summary_chr6_noncoding_cond <- summarizeNoncodingCond("neg_controls", metab = "x1215")

load("../Data/STAAR_neg_controls/x278/Genetic_Noncoding_ncRNA_cond.Rdata")
summary_chr16_noncoding_cond <- summarizeNoncodingCond("neg_controls", metab = "x278",
                                                       results_ncRNA_cond)

load("../Data/STAAR_exp_groups/x1021/Genetic_Noncoding_ncRNA_cond.Rdata")
summary_chr8_noncoding_cond <- summarizeNoncodingCond("exp_groups", metab = "x1021",
                                                      results_ncRNA_cond)

load("../Data/STAAR_exp_groups/x100000007/Genetic_Noncoding_ncRNA_cond.Rdata")
summary_chr10_noncoding_cond <- summarizeNoncodingCond("exp_groups", metab = "x100000007",
                                                       results_ncRNA_cond)

load("../Data/STAAR_exp_groups/x100004046/Genetic_Noncoding_ncRNA_cond.Rdata")
summary_chr13_noncoding_cond <- summarizeNoncodingCond("exp_groups", metab = "x100004046",
                                                       results_ncRNA_cond )

load("../Data/STAAR_exp_groups/x1224/Genetic_Noncoding_ncRNA_cond.Rdata")
summary_chr16_noncoding_cond <- summarizeNoncodingCond("exp_groups", metab = "x1224",
                                                       results_ncRNA_cond)


ncRNA_chr2_cond_acetylarginine <- get(load("../Data/STAAR_exp_groups/x100001266/Genetic_Noncoding_ncRNA_cond.Rdata"))
summary_chr2_noncoding_acetylarginine_cond <- summarizeNoncodingCond("exp_groups",
                                                                     metab = "x100001266",
                                                                     ncRNA_cond_res = ncRNA_chr2_cond_acetylarginine)

# ncRNA_chr5_cond_x1114 <- get(load("../Data/STAAR_exp_groups/x1114/Genetic_Noncoding_ncRNA_cond.Rdata"))
# summary_chr5_noncoding_x1114_cond <- summarizeNoncodingCond("exp_groups",
#                                                             metab = "x1114",
#                                                             ncRNA_cond_res = ncRNA_chr5_cond_x1114)
#
# ncRNA_chr8_x192_cond <- get(load("../Data/STAAR_exp_groups/x192/Genetic_Noncoding_ncRNA_cond.Rdata"))
# summary_chr8_noncoding_x192_cond <- summarizeNoncodingCond("exp_groups",
#                                                            metab = "x192",
#                                                            ncRNA_cond_res = ncRNA_chr8_x192_cond)
# ncRNA_chr11_arachido_cond <- get(load("../Data/STAAR_exp_groups/x100009332/Genetic_Noncoding_ncRNA_cond.Rdata"))
# summary_chr11_noncoding_arachido_cond <- summarizeNoncodingCond("exp_groups",
#                                                                 metab = "x100009332",
#                                                                 ncRNA_cond_res = ncRNA_chr11_arachido_cond)
# ncRNA_chr11_3beta_cond <- get(load("../Data/STAAR_exp_groups/x100006370/Genetic_Noncoding_ncRNA_cond.Rdata"))
# summary_chr11_noncoding_3beta_cond <- summarizeNoncodingCond("exp_groups",
#                                                              metab = "x100006370",
#                                                              ncRNA_cond_res = ncRNA_chr11_3beta_cond)
# ncRNA_chr12_cond_1_methyl <- get(load("../Data/STAAR_exp_groups/x100001208/b2/Genetic_Noncoding_ncRNA_cond_b2.Rdata"))
# summary_chr12_noncoding_1_methyl <- summarizeNoncodingCond("exp_groups",
#                                                            metab = "x100001208",
#                                                            ncRNA_cond_res = ncRNA_chr12_cond_1_methyl)
# ncRNA_chr12_cond_1_ethylmalonate <- get(load("../Data/STAAR_exp_groups/x2054/b2/Genetic_Noncoding_ncRNA_cond_b2.Rdata"))
# summary_chr12_noncoding_ethylmalonate <- summarizeNoncodingCond("exp_groups",
#                                                                 metab = "x2054",
#                                                                 ncRNA_cond_res = ncRNA_chr12_cond_1_ethylmalonate)

ncRNA_chr16_cond_x6264 <- get(load("../Data/STAAR_exp_groups/x100006264/b2/Genetic_Noncoding_ncRNA_cond_b2.Rdata"))
summary_chr16_noncoding_cond_x6264 <- summarizeNoncodingCond("exp_groups",
                                                             metab = "x100006264",
                                                             ncRNA_cond_res = ncRNA_chr16_cond_x6264)


tot_noncoding_cond_b1 <- rbind(summary_chr2_noncoding_acetylarginine_cond,
                            summary_chr16_noncoding_cond_x6264)

tot_noncoding_cond_b1_df <- as.data.frame(matrix(unlist(tot_noncoding_cond_b1), ncol = 91))

colnames(tot_noncoding_cond_b1_df) <- colnames(tot_noncoding_cond_b1)
tot_noncoding_cond_b1_df <-  tot_noncoding_cond_b1_df[,c("Gene name", "Chr", "Category",
                                                   "#SNV", "cMAC",
                                                   "STAAR-S(1,25)", "STAAR-S(1,1)",
                                                   "STAAR-B(1,25)", "STAAR-B(1,1)",
                                                   "STAAR-A(1,25)", "STAAR-A(1,1)",
                                                   "ACAT-O", "STAAR-O")]

tot_noncoding_cond_b1_df$Metabolite <- c(rep("N-acetylarginine",nrow(summary_chr2_noncoding_acetylarginine_cond)),
                                      rep("Propyl 4-hydroxybenzoate sulfate",nrow(summary_chr16_noncoding_cond_x6264)))

tot_noncoding_cond_b1_df$Group <- c(rep("Test Region", nrow(tot_noncoding_cond_b1_df)))

tot_noncoding_cond_b1_df <- unique(tot_noncoding_cond_b1_df) |>
  select(c("Gene name", "Chr", "Metabolite", "Group", "Category",
           "#SNV", "cMAC",
           "STAAR-S(1,25)", "STAAR-S(1,1)",
           "STAAR-B(1,25)", "STAAR-B(1,1)",
           "STAAR-A(1,25)", "STAAR-A(1,1)",
           "ACAT-O", "STAAR-O"))

tot_noncoding_cond_b1_df <- unique(tot_noncoding_cond_b1_df)

# comobind all identified non-coding and coding rv sets
all_rvs_cond_sig_b1 <- rbind(tot_coding_cond_b1_df, tot_noncoding_cond_b1_df)
all_rvs_cond_sig_b1[, c("#SNV")] <- as.integer(all_rvs_cond_sig_b1[, c("#SNV")])

all_rvs_cond_sig_b1 <- data.table(all_rvs_cond_sig_b1)
all_rvs_cond_sig_b1[, 6:14 := lapply(.SD, as.numeric), .SDcols = 6:14]

fwrite(all_rvs_cond_sig_b1,  "../Data/STAAR_b1_all_rvs_cond_sig.csv")


#------------------------------------------------------------------------------
# LEGACY / EARLIER VERSION: v2_scaled_up  (168 unique blocks)
#------------------------------------------------------------------------------
# Kept verbatim. These are blocks that do NOT appear in the current version above
# (mostly hard-coded per-metabolite / per-chromosome run calls and older path setups).

library(GWASTools)

setwd("/Volumes/Sofer Lab/HCHS_SOL/Projects/2024_rare_variants/Code")

###########################################################
#           User Input
## aGDS directory
agds_dir <- get(load("../Data/STAAR_prep/agds_dir_scale_up_final.Rdata"))

## Number of jobs for each chromosome
jobs_num <- get(load("../Data/STAAR_prep/scale_up/mult_chr_jobs_num_scale_up.Rdata"))

## results name
individual_results_name <- "Individual_Analysis"

## QC_label
QC_label <- "annotation/filter"
## geno_missing_imputation
geno_missing_imputation <- "mean"
## method_cond
method_cond <- "optimal"
## alpha level
alpha <- 5e-08

source("./20250111_STAAR_individual_summary_customized.R")
summaryIndividual <- function(agds.path, chr, group, metab, job.num.chr, known_loci = NULL){
  if(group == "exp_groups"){
    ## output path
    output_path <- paste0("../Data/STAAR_exp_groups/",
                          metab, "/")
    load(paste0("../Data/STAAR_exp_groups/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }else{
    output_path <- paste0("../Data/STAAR_neg_controls/",
                          metab, "/")
    load(paste0("../Data/STAAR_neg_controls/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }
  input_path <- output_path

  summary_ind <- Individual_Analysis_Results_Summary(agds_dir=agds.path,jobs_num=job.num.chr,input_path=input_path,output_path=output_path,
                                                     individual_results_name=paste0(metab, "_", individual_results_name),
                                                     obj_nullmodel=obj_nullmodel, known_loci=known_loci,
                                                     method_cond=method_cond, chromosome = chr, alpha = alpha,
                                                     QC_label=QC_label,geno_missing_imputation=geno_missing_imputation,
                                                     manhattan_plot=FALSE,QQ_plot=FALSE)


  return(summary_ind)
}

pruneVariants <- function(agds.path, known_loci, chr, metab, group){
  if(group == "exp_groups"){
    ## output path
    output_path <- paste0("../Data/STAAR_exp_groups/",
                          metab, "/")
    load(paste0("../Data/STAAR_exp_groups/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }else{
    output_path <- paste0("../Data/STAAR_neg_controls/",
                          metab, "/")
    load(paste0("../Data/STAAR_neg_controls/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }
  genofile <- seqOpen(agds.path)
  pruned_variants <- LD_pruning(chr, genofile, obj_nullmodel = obj_nullmodel,
                                variants_list = known_loci[, 1:4], QC_label = QC_label,
                                geno_missing_imputation = geno_missing_imputation,
                                variant_type = "SNV")
  seqClose(genofile)
  return(pruned_variants)
}

known_loci_gwas <- read_csv("../Data/scaled_up_results_table.csv")
known_loci_gwas <- known_loci_gwas |> filter(!is.na(CHR))
############### CONDITIONAL ANALYSIS ######################
source("./20250424_STAAR_individual_customized.R")

#### chr 11 arachido x100009332 ####
known_loci <- as.data.frame(known_loci_gwas[known_loci_gwas$CHR==11 &
                                              known_loci_gwas$Metab_alias == "x100009332" , 1:4])

summary_individual_chr11_cond_arachido <- summaryIndividual(agds_dir[1], chr = 11,
                                                            metab = "x100009332",
                                                            group = "exp_groups",
                                                            job.num.chr = jobs_num[1,],
                                                            known_loci = known_loci)


# sensitivity anlaysis only, pruning all variants together at one step
tot_ind_chr11_arachido_cond <- rbind(known_loci,
                                     summary_individual_chr11_cond_arachido[summary_individual_chr11_cond_arachido$pvalue_cond < 0.001, 1:4])

chr11_arachido_ind_cond_pruned_sensitivity <- pruneVariants(agds_dir[1],
                                                  tot_ind_chr11_arachido_cond,
                                                  11, "x100009332",
                                                  "exp_groups")

saveRDS(chr11_arachido_ind_cond_pruned_sensitivity, "../Data/STAAR_exp_groups/x100009332/individual_cond_pruned_var_sensitivity.RDS")


chr11_cond_arachido_pruned <- pruneVariants(agds_dir[1],
                                      summary_individual_chr11_cond_arachido[summary_individual_chr11_cond_arachido$pvalue_cond < 0.001, 1:4],
                                      11, "x100009332",
                                      "exp_groups")

saveRDS(rbind(known_loci, chr11_cond_arachido_pruned),
        "../Data/STAAR_exp_groups/x100009332/individual_cond_pruned_var.RDS")

###### chr 11 3beta x100006370 ####
known_loci <- as.data.frame(known_loci_gwas[known_loci_gwas$CHR==11 &
                                              known_loci_gwas$Metab_alias == "x100006370" , 1:4])

summary_individual_chr11_cond_3beta <- summaryIndividual(agds_dir[2], chr=11,
                                                            metab = "x100006370",
                                                            job.num.chr = jobs_num[2,],
                                                            group = "exp_groups",
                                                         known_loci = known_loci)

tot_ind_chr11_3beta_cond <- rbind(known_loci,
                                     summary_individual_chr11_cond_3beta[summary_individual_chr11_cond_3beta$pvalue_cond < 0.001, 1:4])

chr11_3beta_ind_cond_pruned_sensitivity <- pruneVariants(agds_dir[2],
                                                            tot_ind_chr11_3beta_cond,
                                                            11, "x100006370",
                                                            "exp_groups")

saveRDS(chr11_3beta_ind_cond_pruned_sensitivity, "../Data/STAAR_exp_groups/x100006370/individual_cond_pruned_var_sensitivity.RDS")


chr11_cond_3beta_pruned <- pruneVariants(agds_dir[2],
                                         summary_individual_chr11_cond_3beta[summary_individual_chr11_cond_3beta$pvalue_cond < 0.001, 1:4],
                                            11, "x100006370",
                                            "exp_groups")

saveRDS(rbind(known_loci, chr11_cond_3beta_pruned),
        "../Data/STAAR_exp_groups/x100006370/individual_cond_pruned_var.RDS")


#### chr 12 1-methyl x100001208 ####
known_loci <- as.data.frame(known_loci_gwas[known_loci_gwas$CHR==12 &
                                              known_loci_gwas$Metab_alias == "x100001208" , 1:4])

summary_individual_chr12_cond_1_methyl <- summaryIndividual(agds_dir[3], chr = 12,
                                                               metab = "x100001208",
                                                               group = "exp_groups",
                                                             job.num.chr = jobs_num[3,],
                                                            known_loci = known_loci)

tot_ind_chr12_1_methyl_cond <- rbind(known_loci,
                                  summary_individual_chr12_cond_1_methyl[summary_individual_chr12_cond_1_methyl$pvalue_cond < 0.001, 1:4])

chr12_1_methyl_ind_cond_pruned_sensitivity <- pruneVariants(agds_dir[3],
                                                         tot_ind_chr12_1_methyl_cond,
                                                         12, "x100001208",
                                                         "exp_groups")

saveRDS(chr12_1_methyl_ind_cond_pruned_sensitivity, "../Data/STAAR_exp_groups/x100001208/individual_cond_pruned_var_sensitivity.RDS")


chr12_cond_1_methyl_pruned <- pruneVariants(agds_dir[3],
                                         summary_individual_chr12_cond_1_methyl[summary_individual_chr12_cond_1_methyl$pvalue_cond < 0.001, 1:4],
                                         12, "x100001208",
                                         "exp_groups")

known_loci <- known_loci[-1,]

saveRDS(rbind(known_loci, chr12_cond_1_methyl_pruned),
        "../Data/STAAR_exp_groups/x100001208/individual_cond_pruned_var.RDS")


##### chr 12 x2054 ethylmalonate ####
known_loci <- as.data.frame(known_loci_gwas[known_loci_gwas$CHR==12 &
                                              known_loci_gwas$Metab_alias == "x2054" , 1:4])

summary_individual_chr12_cond_ethylmalonate <- summaryIndividual(agds_dir[4],
                                                                    chr = 12,
                                                                    metab = "x2054",
                                                                    group = "exp_groups",
                                                                    job.num.chr = jobs_num[4,],
                                                                    known_loci = known_loci)

tot_ind_chr12_ethylmalonate_cond <- rbind(known_loci,
                                     summary_individual_chr12_cond_ethylmalonate[summary_individual_chr12_cond_ethylmalonate$pvalue_cond < 0.001, 1:4])

# NO CHANGE
saveRDS(tot_ind_chr12_ethylmalonate_cond, "../Data/STAAR_exp_groups/x2054/individual_cond_pruned_var_sensitivity.RDS")


chr12_cond_ethylmalonate_pruned <- pruneVariants(agds_dir[4],
                                                 summary_individual_chr12_cond_ethylmalonate[summary_individual_chr12_cond_ethylmalonate$pvalue_cond < 0.001, 1:4],
                                            12, "x2054",
                                            "exp_groups")

# no pruned loci but still create the file to be consistent with other pairs
saveRDS(rbind(known_loci, chr12_cond_ethylmalonate_pruned),
        "../Data/STAAR_exp_groups/x2054/individual_cond_pruned_var.RDS")


#### chr16 propyl x100006264 ####
known_loci <- as.data.frame(known_loci_gwas[known_loci_gwas$CHR==16 &
                                              known_loci_gwas$Metab_alias == "x100006264" , 1:4])


summary_individual_chr16_cond_propyl <- summaryIndividual(agds_dir[5],
                                                                 chr = 16,
                                                                 metab = "x100006264",
                                                                 group = "exp_groups",
                                                                 job.num.chr = jobs_num[5,],
                                                          known_loci = known_loci)

tot_ind_chr16_propyl_cond <- rbind(known_loci,
                                     summary_individual_chr16_cond_propyl[summary_individual_chr16_cond_propyl$pvalue_cond < 0.001, 1:4])

chr16_propyl_ind_cond_pruned_sensitivity <- pruneVariants(agds_dir[5],
                                                            tot_ind_chr16_propyl_cond,
                                                            16, "x100006264",
                                                            "exp_groups")

saveRDS(chr16_propyl_ind_cond_pruned_sensitivity, "../Data/STAAR_exp_groups/x100006264/individual_cond_pruned_var_sensitivity.RDS")


chr16_cond_propyl_pruned <- pruneVariants(agds_dir[5],
                                          summary_individual_chr16_cond_propyl[summary_individual_chr16_cond_propyl$pvalue_cond < 0.001, 1:4],
                                                 16, "x100006264",
                                                 "exp_groups")

known_loci <- known_loci[2:4,]
saveRDS(rbind(known_loci, chr16_cond_propyl_pruned),
        "../Data/STAAR_exp_groups/x100006264/individual_cond_pruned_var.RDS")


########## chr2 x100001266 ########
known_loci <- as.data.frame(known_loci_gwas[known_loci_gwas$CHR==2 &
                                              known_loci_gwas$Metab_alias == "x100001266" , 1:4])


summary_individual_chr2_cond_x1266 <- summaryIndividual(agds_dir[6],
                                                          chr = 2,
                                                          metab = "x100001266",
                                                          group = "exp_groups",
                                                          job.num.chr = jobs_num[6,],
                                                          known_loci = known_loci)


tot_ind_chr2_x1266_cond <- rbind(known_loci,
                                   summary_individual_chr2_cond_x1266[summary_individual_chr2_cond_x1266$pvalue_cond < 0.001, 1:4])

chr2_x1266_ind_cond_pruned_sensitivity <- pruneVariants(agds_dir[6],
                                                          tot_ind_chr2_x1266_cond,
                                                          2, "x100001266",
                                                          "exp_groups")

saveRDS(chr2_x1266_ind_cond_pruned_sensitivity, "../Data/STAAR_exp_groups/x100001266/individual_cond_pruned_var_sensitivity.RDS")


chr2_cond_pruned <- pruneVariants(agds_dir[6],
                                          summary_individual_chr2[summary_individual_chr2$pvalue_cond < 0.001, 1:4],
                                          2, "x100001266",
                                          "exp_groups")

saveRDS(rbind(known_loci, chr2_cond_pruned),
        "../Data/STAAR_exp_groups/x100001266/individual_cond_pruned_var.RDS")

########## chr5 x1114 ########
known_loci <- as.data.frame(known_loci_gwas[known_loci_gwas$CHR==5 &
                                              known_loci_gwas$Metab_alias == "x1114" , 1:4])


summary_individual_chr5_cond_x1114 <- summaryIndividual(agds_dir[7],
                                             chr = 5,
                                             metab = "x1114",
                                             group = "exp_groups",
                                             job.num.chr = jobs_num[7,],
                                             known_loci = known_loci)


tot_ind_chr5_x1114_cond <- rbind(known_loci,
                                 summary_individual_chr5_cond_x1114[summary_individual_chr5_cond_x1114$pvalue_cond < 0.001, 1:4])

chr5_x1114_ind_cond_pruned_sensitivity <- pruneVariants(agds_dir[7],
                                                        tot_ind_chr5_x1114_cond,
                                                        5, "x1114",
                                                        "exp_groups")

saveRDS(chr5_x1114_ind_cond_pruned_sensitivity, "../Data/STAAR_exp_groups/x1114/individual_cond_pruned_var_sensitivity.RDS")


chr5_cond_pruned <- pruneVariants(agds_dir[7],
                                  summary_individual_chr5[summary_individual_chr5$pvalue_cond < 0.001, 1:4],
                                         5, "x1114",
                                         "exp_groups")
saveRDS(rbind(known_loci, chr5_cond_pruned),
        "../Data/STAAR_exp_groups/x1114/individual_cond_pruned_var.RDS")


known_loci <- as.data.frame(known_loci_gwas[known_loci_gwas$CHR==8 &
                                              known_loci_gwas$Metab_alias == "x192" , 1:4])


summary_individual_chr8_cond_x192 <- summaryIndividual(agds_dir[8],
                                             chr = 8,
                                             metab = "x192",
                                             group = "exp_groups",
                                             job.num.chr = jobs_num[8,],
                                             known_loci = known_loci)

tot_ind_chr8_x192_cond <- rbind(known_loci,
                                 summary_individual_chr8_cond_x192[summary_individual_chr8_cond_x192$pvalue_cond < 0.001, 1:4])

chr8_x192_ind_cond_pruned_sensitivity <- pruneVariants(agds_dir[8],
                                                        tot_ind_chr8_x192_cond,
                                                        8, "x192",
                                                        "exp_groups")

saveRDS(chr8_x192_ind_cond_pruned_sensitivity, "../Data/STAAR_exp_groups/x192/individual_cond_pruned_var_sensitivity.RDS")


chr8_cond_pruned <- pruneVariants(agds_dir[8],
                                  summary_individual_chr8[summary_individual_chr8$pvalue_cond < 0.001, 1:4],
                                  8, "x192",
                                  "exp_groups")

saveRDS(rbind(known_loci, chr8_cond_pruned),
        "../Data/STAAR_exp_groups/x192/individual_cond_pruned_var.RDS")


############ Coding summary ############
load("../Data/STAAR_prep/agds_dir_scale_up_final.Rdata")
## rs channel name in aGDS
rs_channel <- "annotation/info/FunctionalAnnotation/rsid"
## output path
output_path <- "../Data/"
## number of jobs
gene_centric_coding_jobs_num <- 1
gene_centric_results_name <- "Genecentric_coding_uncond"
#gene_centric_results_name <- "Genecentric_coding_cond" # conditional analysis

## variant_type
variant_type <- "SNV"
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
                          metab, "/")
    load(paste0("../Data/STAAR_exp_groups/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }else{
    output_path <- paste0("../Data/STAAR_neg_controls/",
                          metab, "/")
    load(paste0("../Data/STAAR_neg_controls/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }
  input_path <- output_path
  if(p_threshold_mode == "nominal"){
    output_path <- paste0(output_path, "Nominal_p_threshold/")
  }
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


summary_chr11_coding_uncond_arachido <- summaryCoding(agds_dir[1], chr = 11,
                                                      metab = "x100009332",
                                                      group = "exp_groups"
                                                      )

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


summarizeCodingCond <- function(group, metab){
  if(group == "exp_groups"){
    ## output path
    output_path <- paste0("../Data/STAAR_exp_groups/",
                          metab, "/")
  }else{
    output_path <- paste0("../Data/STAAR_neg_controls/",
                          metab, "/")
  }

  res_total_coding <- c()
  for(category in categories){
    load(paste0(output_path, metab, "_Genecentric_coding_cond_", category,
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
    filter((`STAAR-B(1,25)` < alpha_cond) |  `STAAR-O` < alpha_cond)

  return(res_total_coding_sig)
}

################# gene-centric non-coding #########################
gene_centric_noncoding_jobs_num <- 1
ncRNA_jobs_num <- 1
gene_centric_results_name <- "Genecentric_noncoding_uncond"

alpha <- 1
## Use_annotation_weights
## ncRNA
for(metab in metabs){
  for(cat_func in categories){
    if(group == "exp_groups"){
      ## output path
      output_path <- paste0("../Data/STAAR_exp_groups/",
                            metab, "/")
      res <- getobj(paste0("../Data/STAAR_exp_groups/",
                           metab, "/", cat_func, ".Rdata"))
      ncrna_res <- getobj(paste0("../Data/STAAR_exp_groups/",
                                  metab, "/", ncRNA_results_name, ".Rdata"))
      res <- rbind(res, ncrna_res)
    }else{
      output_path <- paste0("../Data/STAAR_neg_controls/",
                            metab, "/")
      res <- getobj(paste0("../Data/STAAR_neg_controls/",
                           metab, "/", cat_func, ".Rdata"))
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
                            metab, "/")
      res <- getobj(paste0("../Data/STAAR_exp_groups/",
                           metab, "/", cat_func, ".Rdata"))
      ncrna_res <- getobj(paste0("../Data/STAAR_exp_groups/",
                                 metab, "/", ncRNA_results_name, ".Rdata"))
      res <- rbind(res, ncrna_res)
    }else{
      output_path <- paste0("../Data/STAAR_neg_controls/",
                            metab, "/")
      res <- getobj(paste0("../Data/STAAR_neg_controls/",
                           metab, "/", cat_func, ".Rdata"))
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
#1.37741e-05
## gene info
source("./20250422_gencentric_noncoding_summary.R")
summaryNonCoding <- function(agds.path, chr, group, metab, job.num.chr, known_loci = NULL){
  if(group == "exp_groups"){
    ## output path
    output_path <- paste0("../Data/STAAR_exp_groups/",
                          metab, "/")
    load(paste0("../Data/STAAR_exp_groups/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }else{
    output_path <- paste0("../Data/STAAR_neg_controls/",
                          metab, "/")
    load(paste0("../Data/STAAR_neg_controls/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }
  input_path <- output_path
  ncRNA_output_path <- output_path
  ncRNA_input_path <- ncRNA_output_path

  noncoding_summary <- Gene_Centric_Noncoding_Results_Summary(agds_dir=agds_dir,gene_centric_noncoding_jobs_num=gene_centric_noncoding_jobs_num,
                                                              input_path=input_path,output_path=output_path,
                                                              gene_centric_results_name=paste0(metab, "_" , gene_centric_results_name),
                                                              ncRNA_jobs_num=ncRNA_jobs_num,ncRNA_input_path=ncRNA_input_path,
                                                              ncRNA_output_path=ncRNA_output_path,ncRNA_results_name=ncRNA_results_name,
                                                              obj_nullmodel=obj_nullmodel, # known_loci=known_loci,
                                                              method_cond=method_cond, alpha = alpha, alpha_ncRNA = 1e-4,
                                                              QC_label=QC_label,geno_missing_imputation=geno_missing_imputation,variant_type=variant_type,
                                                              Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                              Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                              ncRNA_pos=ncRNA_pos) #, manhattan_plot=TRUE,QQ_plot=TRUE)

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
col_names <- colnames(summary_chr8_noncoding_uncond_putrescine)
summary_chr8_noncoding_uncond_putrescine <- as.data.frame(matrix(unlist(summary_chr8_noncoding_uncond_putrescine),
                                                                 ncol = 91))
colnames(summary_chr8_noncoding_uncond_putrescine) <- col_names


summarizeNoncodingCond <- function(group, metab, ncRNA_cond_res = NULL){
  if(group == "exp_groups"){
    ## output path
    output_path <- paste0("../Data/STAAR_exp_groups/",
                          metab, "/")
  }else{
    output_path <- paste0("../Data/STAAR_neg_controls/",
                          metab, "/")
  }

  res_total_noncoding <- c()
  for(category in categories){
    load(paste0(output_path, metab, "_Genecentric_noncoding_cond_", category,
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

ncRNA_chr2_cond <- get(load("../Data/STAAR_exp_groups/x100001266/Genetic_Noncoding_ncRNA_cond.Rdata"))
summary_chr2_noncoding_cond <- summarizeNoncodingCond("exp_groups",
                                                      metab = "x100001266",
                                                        ncRNA_cond_res = ncRNA_chr2_cond)

saveRDS(summary_chr2_noncoding_cond, "../Data/STAAR_exp_groups/x100001266/noncoding_cond_sig.RDS")
summary_chr2_noncoding_cond <- summary_chr2_noncoding_cond[, c(colnames(summary_chr2_noncoding_cond)[1:5],
                                                               "SKAT(1,25)", "STAAR-B(1,25)",
                                                               "ACAT-V(1,25)",	"STAAR-O")]

ncRNA_chr5_cond <- get(load("../Data/STAAR_exp_groups/x1114/Genetic_Noncoding_ncRNA_cond.Rdata"))
summary_chr5_noncoding_cond <- summarizeNoncodingCond("exp_groups",
                                                      metab = "x1114",
                                                      ncRNA_cond_res = ncRNA_chr5_cond)

ncRNA_chr8_cond <- get(load("../Data/STAAR_exp_groups/x192/Genetic_Noncoding_ncRNA_cond.Rdata"))
summary_chr8_noncoding_cond <- summarizeNoncodingCond("exp_groups",
                                                      metab = "x192",
                                                      ncRNA_cond_res = ncRNA_chr8_cond)

ncRNA_chr16_cond <- get(load("../Data/STAAR_exp_groups/x100006264/Genetic_Noncoding_ncRNA_cond.Rdata"))
summary_chr16_noncoding_cond <- summarizeNoncodingCond("exp_groups",
                                                       metab = "x100006264",
                                                       ncRNA_cond_res = ncRNA_chr16_cond)


#### RUN THIS AT LAST ####
all_rvs_cond_sig <- as.data.frame(matrix(unlist(rbind(summary_chr11_coding_cond_arachido,
                                                      summary_chr12_coding_cond_1_methyl,
                                                      summary_chr12_coding_cond_ethylmalonate,
                                                      summary_chr5_coding_cond_3amino,
                                                      summary_chr2_noncoding_cond,
                                                      summary_chr16_noncoding_cond)), ncol = 91))
colnames(all_rvs_cond_sig) <- names(summary_chr11_coding_cond_arachido)
all_rvs_cond_sig <-  all_rvs_cond_sig[,c("Gene name", "Chr", "Category",
                                         "#SNV", "cMAC",
                                         "STAAR-S(1,25)", "STAAR-S(1,1)",
                                         "STAAR-B(1,25)", "STAAR-B(1,1)",
                                         "STAAR-A(1,25)", "STAAR-A(1,1)",
                                         "ACAT-O", "STAAR-O")]


all_rvs_cond_sig <- unique(all_rvs_cond_sig)
all_rvs_cond_sig_1 <- fread("../Data/STAAR_b1_all_rvs_cond_sig.csv")
all_rvs_cond_sig <- rbind(all_rvs_cond_sig_1, all_rvs_cond_sig)
fwrite(all_rvs_cond_sig,  "../Data/STAAR_b1_all_rvs_cond_sig.csv")


#------------------------------------------------------------------------------
# LEGACY / EARLIER VERSION: v1_original  (161 unique blocks)
#------------------------------------------------------------------------------
# Kept verbatim. These are blocks that do NOT appear in the current version above
# (mostly hard-coded per-metabolite / per-chromosome run calls and older path setups).

load("../Data/STAAR_prep/agds_dir_multichr.Rdata")
## Input GWASCatalog variants (#rs)
# # From GWAS catalog 22 variants
# known_variants <- read.table("../Data/GWAS_catalog_EFO_0800122_associations_export.tsv",
#                              header = TRUE)
# known_loci <- unique(sub("-.*$", "", known_variants$riskAllele)) # 16 variants in total
# known_loci <- as.data.frame(known_loci)
# colnames(known_loci) <- c("rs")

## output file name

agds_dir <- get(load("../Data/STAAR_prep/agds_dir_multichr.Rdata"))

jobs_num <- get(load("../Data/STAAR_prep/mult_chr_jobs_num.Rdata"))

alpha <- 5E-08

# known_loci_gwas_chr5_pruned <- pruneVariants(agds_dir[2],
#                                              known_loci_gwas_chr5, 5, "x799",
#                                              "neg_controls")
# known_loci_gwas_chr6_pruned <- pruneVariants(agds_dir[3],
#                                              known_loci_gwas_chr6, 6, "x1215",
summary_chr2_ind <- summaryIndividual(agds_dir[1], chr = 2,
                                                 metab = "x100001721",
                                                 group = "neg_controls",
                                                job.num.chr = jobs_num[1,])

summary_chr5_ind<- summaryIndividual(agds_dir[2], chr = 5, metab = "x799",
                                                 group = "neg_controls",
                                                job.num.chr = jobs_num[2,])

summary_chr6_ind <- summaryIndividual(agds_dir[3], chr = 6, metab = "x1215",
                                                 group = "neg_controls",
                                                job.num.chr = jobs_num[3,])

summary_chr16_ind <- summaryIndividual(agds_dir[7], chr = 16, metab = "x278",
                                                  group = "neg_controls",
                                       job.num.chr = jobs_num[7,])

library(data.table)
known_loci_gwas <- fread("../Data/STAAR_neg_controls/Known_GWAS_sig_var_neg_ctrl.csv")
known_loci_gwas <- known_loci_gwas[!is.na(known_loci_gwas$CHR),]
common_var_tot <- c()


###### chr2 x100001721 #######
summary_chr2_ind_cond <- summaryIndividual(agds_dir[1], chr = 2,
                                      metab = "x100001721",
                                      group = "neg_controls",
                                      job.num.chr = jobs_num[1,],
                                      known_loci =  as.data.frame(known_loci_gwas[known_loci_gwas$CHR==2, 1:4]))

tot_ind_chr2_cond <- rbind(known_loci_gwas[known_loci_gwas$CHR==2 & known_loci_gwas$POS == 73650092, 1:4],
                           summary_chr2_ind_cond[summary_chr2_ind_cond$pvalue_cond < 0.001, 1:4])


chr2_ind_cond_pruned_sensitivity <- pruneVariants(agds_dir[1],
                                      tot_ind_chr2_cond,
                                      2, "x100001721",
                                      "neg_controls")

saveRDS(chr2_ind_cond_pruned_sensitivity,
         "../Data/STAAR_neg_controls/x100001721/individual_cond_pruned_var_sensitivity.RDS")


# original code, bind known + pruned together
chr2_ind_cond_pruned <- pruneVariants(agds_dir[1],
                                      summary_chr2_ind_cond[summary_chr2_ind_cond$pvalue_cond < 0.001, 1:4],
                                      2, "x100001721",
                                      "neg_controls")


saveRDS(rbind(known_loci_gwas[known_loci_gwas$CHR==2 & known_loci_gwas$POS == 73650092, 1:4], chr2_ind_cond_pruned), "../Data/STAAR_neg_controls/x100001721/individual_cond_pruned_var.RDS")
common_var_tot <- rbind(common_var_tot, rbind(known_loci_gwas[known_loci_gwas$CHR==2, 1:4], chr2_ind_cond_pruned))


###### chr 5 x799 #######
summary_chr5_ind_cond<- summaryIndividual(agds_dir[2], chr = 5, metab = "x799",
                                     group = "neg_controls",
                                     job.num.chr = jobs_num[2,],
                                     known_loci = as.data.frame(known_loci_gwas[known_loci_gwas$CHR==5, 1:4]))

known_loci_gwas_chr5 <- known_loci_gwas[known_loci_gwas$CHR==5, 1:4]
tot_ind_chr5_cond <- rbind(known_loci_gwas_chr5,
                           summary_chr5_ind_cond[summary_chr5_ind_cond$pvalue_cond < 0.001, 1:4])


chr5_ind_cond_pruned_sensitivity <- pruneVariants(agds_dir[2],
                                                  tot_ind_chr5_cond,
                                                  5, "x799",
                                                  "neg_controls")

saveRDS(chr5_ind_cond_pruned_sensitivity,
        "../Data/STAAR_neg_controls/x799/individual_cond_pruned_var_sensitivity.RDS")

chr5_ind_cond_pruned <- pruneVariants(agds_dir[2],
                                       summary_chr5_ind_cond[summary_chr5_ind_cond$pvalue_cond < 0.001, 1:4],
                                       5, "x799",
                                       "neg_controls")


saveRDS(rbind(known_loci_gwas_chr5, chr5_ind_cond_pruned), "../Data/STAAR_neg_controls/x100001721/individual_cond_pruned_var.RDS")
common_var_tot <- rbind(common_var_tot,  rbind(known_loci_gwas_chr5, chr5_ind_cond_pruned))


###### chr 6 x1215 #######
summary_chr6_ind_cond <- summaryIndividual(agds_dir[3], chr = 6, metab = "x1215",
                                      group = "neg_controls",
                                      job.num.chr = jobs_num[3,],
                                      known_loci = as.data.frame(known_loci_gwas[known_loci_gwas$CHR==6, 1:4]))


known_loci_gwas_chr6 <- as.data.frame(known_loci_gwas[known_loci_gwas$CHR==6, 1:4])
tot_ind_chr6_cond <- rbind(known_loci_gwas_chr6,
                           summary_chr6_ind_cond[summary_chr6_ind_cond$pvalue_cond < 0.001, 1:4])

chr6_ind_cond_pruned_sensitivity <- pruneVariants(agds_dir[3],
                                                  tot_ind_chr6_cond,
                                                  6, "x1215",
                                                  "neg_controls")

saveRDS(chr6_ind_cond_pruned_sensitivity, "../Data/STAAR_neg_controls/x1215/individual_cond_pruned_var_sensitivity.RDS")

chr6_ind_cond_pruned <- pruneVariants(agds_dir[3],
                                      summary_chr6_ind_cond[summary_chr6_ind_cond$pvalue_cond < 0.001, 1:4],
                                      6, "x1215",
                                      "neg_controls") # null

known_loci_gwas_chr6 <- known_loci_gwas_chr6[c(1:2), ]
saveRDS(known_loci_gwas_chr6, "../Data/STAAR_neg_controls/x1215/individual_cond_pruned_var.RDS")

common_var_tot <- rbind(common_var_tot, known_loci_gwas_chr6)

###### chr 16 x278 #######
known_loci_gwas_chr16 <- as.data.frame(known_loci_gwas[known_loci_gwas$CHR==16, 1:4])
summary_chr16_ind_cond <- summaryIndividual(agds_dir[7], chr = 16, metab = "x278",
                                       group = "neg_controls",
                                       job.num.chr = jobs_num[7,],
                                       known_loci = known_loci_gwas_chr16)

tot_ind_chr16_cond <- rbind(known_loci_gwas_chr16,
                           summary_chr16_ind_cond[summary_chr16_ind_cond$pvalue_cond < 0.001, 1:4])

chr16_ind_cond_pruned_sensitivity <- pruneVariants(agds_dir[7],
                                                  tot_ind_chr16_cond,
                                                  16, "x278",
                                                  "neg_controls")

saveRDS(chr16_ind_cond_pruned_sensitivity, "../Data/STAAR_neg_controls/x278/individual_cond_pruned_var_sensitivity.RDS")


chr16_ind_cond_pruned <- pruneVariants(agds_dir[7],
                                      summary_chr16_ind_cond[summary_chr16_ind_cond$pvalue_cond < 0.001, 1:4],
                                      16, "x278",
                                      "neg_controls")

saveRDS(rbind(known_loci_gwas[known_loci_gwas$CHR==16, 1:4], chr16_ind_cond_pruned), "../Data/STAAR_neg_controls/x278/individual_cond_pruned_var.RDS")
common_var_tot <- rbind(common_var_tot, rbind(known_loci_gwas[known_loci_gwas$CHR==16, 1:4], chr16_ind_cond_pruned))
common_var_tot$group <- "negative control"


########## Experimental groups ##########
#### chr8 x1021 ####
known_loci_gwas_exp <- read_csv("../Data/STAAR_exp_groups/Known_GWAS_sig_var_exp_groups.csv")

summary_chr8_ind_cond <- summaryIndividual(agds_dir[4], chr = 8, metab = "x1021",
                                            group = "exp_groups",
                                            job.num.chr = jobs_num[4,],
                                            known_loci = as.data.frame(known_loci_gwas_exp[known_loci_gwas_exp$CHR==8, 1:4]))


known_loci_gwas_chr8 <- known_loci_gwas_exp[known_loci_gwas_exp$CHR==8, 1:4]
tot_ind_chr8_cond <- rbind(known_loci_gwas_chr8,
                           summary_chr8_ind_cond[summary_chr8_ind_cond$pvalue_cond < 0.001, 1:4])

chr8_ind_cond_pruned_sensitivity <- pruneVariants(agds_dir[4],
                                                   tot_ind_chr8_cond,
                                                   8, "x1021",
                                                   "exp_groups")

saveRDS(chr8_ind_cond_pruned_sensitivity, "../Data/STAAR_exp_groups/x1021/individual_cond_pruned_var_sensitivity.RDS")


chr8_ind_cond_pruned <- pruneVariants(agds_dir[4],
                                       summary_chr8_ind_cond[summary_chr8_ind_cond$pvalue_cond < 0.001, 1:4],
                                       8, "x1021",
                                       "exp_groups")

saveRDS(rbind(known_loci_gwas_exp[known_loci_gwas_exp$CHR==8, 1:4], chr8_ind_cond_pruned),
        "../Data/STAAR_exp_groups/x1021/individual_cond_pruned_var.RDS")


##### chr10 x100000007 ####
summary_chr10_ind_cond <- summaryIndividual(agds_dir[5], chr = 10, metab = "x100000007",
                                           group = "exp_groups",
                                           job.num.chr = jobs_num[5,],
                                           known_loci = as.data.frame(known_loci_gwas_exp[known_loci_gwas_exp$CHR==10, 1:4]))


known_loci_gwas_chr10 <- known_loci_gwas_exp[known_loci_gwas_exp$CHR==10, 1:4]
tot_ind_chr10_cond <- rbind(known_loci_gwas_chr10,
                           summary_chr10_ind_cond[summary_chr10_ind_cond$pvalue_cond < 0.001, 1:4])

chr10_ind_cond_pruned_sensitivity <- pruneVariants(agds_dir[5],
                                                  tot_ind_chr10_cond,
                                                  10, "x100000007",
                                                  "exp_groups")


saveRDS(chr10_ind_cond_pruned_sensitivity,
        "../Data/STAAR_exp_groups/x100000007/individual_cond_pruned_var_sensitivity.RDS")


chr10_ind_cond_pruned <- pruneVariants(agds_dir[5],
                                      summary_chr10_ind_cond[summary_chr10_ind_cond$pvalue_cond < 0.001, 1:4],
                                      10, "x100000007",
                                      "exp_groups")

known_loci_gwas_chr10 <- known_loci_gwas_chr10[c(1:2),]
saveRDS(rbind(known_loci_gwas_chr10, chr10_ind_cond_pruned),
        "../Data/STAAR_exp_groups/x100000007/individual_cond_pruned_var.RDS")


#### chr13 x100004046 ####
summary_chr13_ind_cond <- summaryIndividual(agds_dir[6], chr = 13, metab = "x100004046",
                                            group = "exp_groups",
                                            job.num.chr = jobs_num[6,],
                                            known_loci = as.data.frame(known_loci_gwas_exp[known_loci_gwas_exp$CHR==13, 1:4]))


known_loci_gwas_chr13 <- known_loci_gwas_exp[known_loci_gwas_exp$CHR==13, 1:4]
tot_ind_chr13_cond <- rbind(known_loci_gwas_chr13,
                            summary_chr13_ind_cond[summary_chr13_ind_cond$pvalue_cond < 0.001, 1:4])

chr13_ind_cond_pruned_sensitivity <- pruneVariants(agds_dir[6],
                                                   tot_ind_chr13_cond,
                                                   13, "x100004046",
                                                   "exp_groups")


saveRDS(chr13_ind_cond_pruned_sensitivity,
        "../Data/STAAR_exp_groups/x100004046/individual_cond_pruned_var_sensitivity.RDS")


chr13_ind_cond_pruned <- pruneVariants(agds_dir[6],
                                       summary_chr13_ind_cond[summary_chr13_ind_cond$pvalue_cond < 0.001, 1:4],
                                       13, "x100004046",
                                       "exp_groups")

saveRDS(rbind(known_loci_gwas_exp[known_loci_gwas_exp$CHR==13, 1:4], chr13_ind_cond_pruned),
        "../Data/STAAR_exp_groups/x100004046/individual_cond_pruned_var.RDS")

#### chr13 x1224 ####
summary_chr16_ind_cond_exp <- summaryIndividual(agds_dir[7], chr = 16, metab = "x1224",
                                            group = "exp_groups",
                                            job.num.chr = jobs_num[7,],
                                            known_loci = as.data.frame(known_loci_gwas_exp[known_loci_gwas_exp$CHR==16, 1:4]))


known_loci_gwas_chr16_1224 <- known_loci_gwas_exp[known_loci_gwas_exp$CHR==16, 1:4]
tot_ind_chr16_cond_1224 <- rbind(known_loci_gwas_chr16_1224,
                                 summary_chr16_ind_cond_exp[summary_chr16_ind_cond_exp$pvalue_cond < 0.001, 1:4])

chr16_ind_cond_pruned_sensitivity <- pruneVariants(agds_dir[7],
                                                   tot_ind_chr16_cond_1224,
                                                   16, "x1224",
                                                   "exp_groups")

saveRDS(chr16_ind_cond_pruned_sensitivity,
        "../Data/STAAR_exp_groups/x1224/individual_cond_pruned_var_sensitivity.RDS")


chr16_ind_cond_pruned_exp <- pruneVariants(agds_dir[7],
                                       summary_chr16_ind_cond_exp[summary_chr16_ind_cond_exp$pvalue_cond < 0.001, 1:4],
                                       16, "x1224",
                                       "exp_groups")

saveRDS(rbind(known_loci_gwas_exp[known_loci_gwas_exp$CHR==16, 1:4], chr16_ind_cond_pruned_exp),
        "../Data/STAAR_exp_groups/x1224/individual_cond_pruned_var.RDS")

################ gene-centric rare variant sets ###############

rm(list=ls())
gc()

alpha <- 2.522704e-05 # for details see STAAR_multiple_chr

summaryCoding <- function(agds.path, chr, group, metab, job.num.chr,
                          known_loci = NULL,
                          p_threshold_mode = "bonferroni"){
  if(group == "exp_groups"){
    ## output path
    output_path <- paste0("../Data/STAAR_exp_groups/",
                          metab, "/")
    load(paste0("../Data/STAAR_exp_groups/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }else{
    output_path <- paste0("../Data/STAAR_neg_controls/",
                          metab, "/")
    load(paste0("../Data/STAAR_neg_controls/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }
  input_path <- output_path
  coding_summary <- Gene_Centric_Coding_Results_Summary(agds_dir=agds_dir,gene_centric_coding_jobs_num=1,
                                    input_path=input_path,output_path=output_path, chr = chr,
                                    gene_centric_results_name=paste0(metab, "_", gene_centric_results_name),
                                    obj_nullmodel=obj_nullmodel, #known_loci=known_loci,
                                    method_cond=method_cond,
                                    QC_label=QC_label,geno_missing_imputation=geno_missing_imputation,variant_type=variant_type,
                                    Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                    Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                    alpha = alpha,manhattan_plot=FALSE,QQ_plot=FALSE)
  return(coding_summary)
  }

summary_chr2_coding <- summaryCoding(agds_dir[1], chr = 2,
                                      metab = "x100001721",
                                      group = "neg_controls",
                                      job.num.chr = jobs_num[1,])

summary_chr5_coding <- summaryCoding(agds_dir[2], chr = 5, metab = "x799",
                                     group = "neg_controls",
                                    job.num.chr = jobs_num[2,])

summary_chr6_coding <- summaryCoding(agds_dir[3], chr = 6, metab = "x1215",
                                      group = "neg_controls",
                                      job.num.chr = jobs_num[3,])

summary_chr16_coding <- summaryCoding(agds_dir[7], chr = 16, metab = "x278",
                                       group = "neg_controls")

summary_chr8_coding <- summaryCoding(agds_dir[4], chr = 8, metab = "x1021",
                                     group = "exp_groups",
                                     job.num.chr = jobs_num[2,])

summary_chr10_coding <- summaryCoding(agds_dir[5], chr = 10, metab = "x100000007",
                                     group = "exp_groups",
                                     job.num.chr = jobs_num[3,])

summary_chr13_coding <- summaryCoding(agds_dir[6], chr = 13, metab = "x100004046",
                                       group = "exp_groups",
                                       job.num.chr = jobs_num[6,])

summary_chr16_coding_exp <- summaryCoding(agds_dir[7], chr = 16, metab = "x1224",
                                      group = "exp_groups",
                                      job.num.chr = jobs_num[7,])

# change agds_dir
summary_chr11_coding_uncond_arachido <- summaryCoding(agds_dir[1], chr = 11,
                                                      metab = "x100009332",
                                                      group = "exp_groups"
)

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

coding_rvs_uncond_sig <- rbind(summary_chr2_coding, summary_chr6_coding,
                               summary_chr2_coding_uncond_arginine,
                               summary_chr5_coding_uncond_3amino,
                     summary_chr8_coding,
                     summary_chr11_coding_uncond_arachido,
                     summary_chr12_coding_uncond_1_methyl,
                     summary_chr12_coding_uncond_ethylmalonate,
                     summary_chr16_coding_uncond_propyl,
                     summary_chr16_coding_exp)


tot_coding_uncond_b1_df <- as.data.frame(matrix(unlist(coding_rvs_uncond_sig), ncol = 13))

colnames(tot_coding_uncond_b1_df) <- colnames(coding_rvs_uncond_sig)
tot_coding_uncond_b1_df <-  tot_coding_uncond_b1_df[,c("Gene name", "Chr", "Category",
                                                   "#SNV", "cMAC",
                                                   "STAAR-S(1,25)", "STAAR-S(1,1)",
                                                   "STAAR-B(1,25)", "STAAR-B(1,1)",
                                                   "STAAR-A(1,25)", "STAAR-A(1,1)",
                                                   "ACAT-O", "STAAR-O")]

tot_coding_uncond_b1_df$Metabolite <- c(rep("N2-acetyllysine",nrow(summary_chr2_coding)),
                                        rep("N-acetylglucosaminylasparagine",nrow(summary_chr6_coding)),
                                        rep("N-acetylarginine",nrow(summary_chr2_coding_uncond_arginine)),
                                      rep("3-aminoisobutyrate",nrow(summary_chr5_coding_uncond_3amino)),
                                      rep("5-oxoproline",nrow(summary_chr8_coding)),
                                      rep("Arachidonoylcholine",nrow(summary_chr11_coding_uncond_arachido)),
                                      rep("1-methylimidazoleacetate",nrow(summary_chr12_coding_uncond_1_methyl)),
                                      rep("Ethylmalonate",nrow(summary_chr12_coding_uncond_ethylmalonate)),
                                      rep("Propyl 4-hydroxybenzoate sulfate",nrow(summary_chr16_coding_uncond_propyl)),
                                      rep("Cys-gly, oxidized",nrow(summary_chr16_coding_exp)))

tot_coding_uncond_b1_df$Group <- c(rep("Negative Control", nrow(summary_chr2_coding) +
                                         nrow(summary_chr6_coding)),
                                 rep("Test Region", nrow(tot_coding_uncond_b1_df)-
                                       nrow(summary_chr2_coding) -nrow(summary_chr6_coding)))

tot_coding_uncond_b1_df <- unique(tot_coding_uncond_b1_df) |>
  select(c("Gene name", "Chr", "Metabolite", "Group", "Category",
           "#SNV", "cMAC",
           "STAAR-S(1,25)", "STAAR-S(1,1)",
           "STAAR-B(1,25)", "STAAR-B(1,1)",
           "STAAR-A(1,25)", "STAAR-A(1,1)",
           "ACAT-O", "STAAR-O"))

data.table::fwrite(tot_coding_uncond_b1_df, file = "../Data/STAAR_b1_coding_rvs_uncond_sig.csv")

alpha <- 0.05/3063
alpha_ncrna <- 6.738544e-05
summaryNonCoding <- function(agds.path, chr, group, metab, job.num.chr, known_loci = NULL){
  if(group == "exp_groups"){
    ## output path
    output_path <- paste0("../Data/STAAR_exp_groups/",
                          metab, "/")
    load(paste0("../Data/STAAR_exp_groups/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }else{
    output_path <- paste0("../Data/STAAR_neg_controls/",
                          metab, "/")
    load(paste0("../Data/STAAR_neg_controls/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
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
                                         ncRNA_pos=ncRNA_pos) #, manhattan_plot=TRUE,QQ_plot=TRUE)

  return(noncoding_summary)
 }

summary_chr2_noncoding <- summaryNonCoding(agds_dir[1], chr = 2, metab = "x100001721",
                                     group = "neg_controls",
                                     job.num.chr = jobs_num[1,])
#known_loci = known_loci_chr2)

summary_chr5_noncoding <- summaryNonCoding(agds_dir[2], chr = 5, metab = "x799",
                                     group = "neg_controls",
                                     job.num.chr = jobs_num[2,])

summary_chr6_noncoding <- summaryNonCoding(agds_dir[3], chr = 6, metab = "x1215",
                                     group = "neg_controls",
                                     job.num.chr = jobs_num[3,])

summary_chr16_noncoding <- summaryNonCoding(agds_dir[7], chr = 16, metab = "x278",
                                      group = "neg_controls",
                                      job.num.chr = jobs_num[7,])

summary_chr8_noncoding <- summaryNonCoding(agds_dir[4], chr = 8, metab = "x1021",
                                     group = "exp_groups",
                                     job.num.chr = jobs_num[2,])

summary_chr10_noncoding <- summaryNonCoding(agds_dir[5], chr = 10, metab = "x100000007",
                                      group = "exp_groups",
                                      job.num.chr = jobs_num[3,])

summary_chr13_noncoding <- summaryNonCoding(agds_dir[6], chr = 13, metab = "x100004046",
                                      group = "exp_groups",
                                      job.num.chr = jobs_num[6,])

summary_chr16_noncoding_exp <- summaryNonCoding(agds_dir[7], chr = 16, metab = "x1224",
                                          group = "exp_groups",
                                          job.num.chr = jobs_num[7,])

summary_chr11_noncoding_uncond_3beta <- summaryNonCoding(agds_dir[2], chr = 11,
                                                         metab = "x100006370",
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

noncoding_rvs_uncond_sig <- rbind(summary_chr2_noncoding,
                                  summary_chr5_noncoding,
                               summary_chr6_noncoding,
                               summary_chr16_noncoding,
                               summary_chr2_noncoding_uncond_arginine,
                               summary_chr5_noncoding_uncond_3amino,
                               summary_chr8_noncoding,
                               summary_chr10_noncoding, summary_chr13_noncoding,
                               summary_chr16_noncoding_uncond_propyl,
                               summary_chr16_noncoding_exp)


tot_noncoding_uncond_b1_df <- as.data.frame(matrix(unlist(noncoding_rvs_uncond_sig), ncol = 91))

colnames(tot_noncoding_uncond_b1_df) <- colnames(noncoding_rvs_uncond_sig)
tot_noncoding_uncond_b1_df <-  tot_noncoding_uncond_b1_df[,c("Gene name", "Chr", "Category",
                                                       "#SNV", "cMAC",
                                                       "STAAR-S(1,25)", "STAAR-S(1,1)",
                                                       "STAAR-B(1,25)", "STAAR-B(1,1)",
                                                       "STAAR-A(1,25)", "STAAR-A(1,1)",
                                                       "ACAT-O", "STAAR-O")]

tot_noncoding_uncond_b1_df$Metabolite <- c(rep("N2-acetyllysine",nrow(summary_chr2_noncoding)),
                                        rep("N-acetylarginine",nrow(summary_chr2_noncoding_uncond_arginine)),
                                        rep("3-aminoisobutyrate",nrow(summary_chr5_noncoding_uncond_3amino)),
                                        rep("5-oxoproline",nrow(summary_chr8_noncoding)),
                                        rep("Propyl 4-hydroxybenzoate sulfate",nrow(summary_chr16_noncoding_uncond_propyl))
                                      )

tot_noncoding_uncond_b1_df$Group <- c(rep("Negative Control", nrow(summary_chr2_noncoding)),
                                   rep("Test Region", nrow(tot_noncoding_uncond_b1_df)-
                                         nrow(summary_chr2_noncoding)))

tot_noncoding_uncond_b1_df <- unique(tot_noncoding_uncond_b1_df) |>
  select(c("Gene name", "Chr", "Metabolite", "Group", "Category",
           "#SNV", "cMAC",
           "STAAR-S(1,25)", "STAAR-S(1,1)",
           "STAAR-B(1,25)", "STAAR-B(1,1)",
           "STAAR-A(1,25)", "STAAR-A(1,1)",
           "ACAT-O", "STAAR-O"))

data.table::fwrite(tot_noncoding_uncond_b1_df, file = "../Data/STAAR_b1_noncoding_rvs_uncond_sig.csv")


summarizeNoncodingCond <- function(group, metab, ncRNA_cond_res = NULL){
  if(group == "exp_groups"){
    ## output path
    output_path <- paste0("../Data/STAAR_exp_groups/",
                          metab, "/")
  }else{
    output_path <- paste0("../Data/STAAR_neg_controls/",
                          metab, "/")
  }

  res_total_noncoding <- c()
  for(category in categories){
    load(paste0(output_path, metab, "_Genecentric_noncoding_cond_", category,
                ".Rdata"))
    res_total_noncoding <- rbind(res_total_noncoding, results_noncoding)
  }

  res_total_noncoding <- as.data.frame(res_total_noncoding)
  #print(res_total_noncoding)
  res_total_noncoding <- res_total_noncoding |>
    filter((`STAAR-B(1,25)` < alpha_cond) |  `STAAR-O` < alpha_cond)
  if(! is.null(ncRNA_cond_res)){
    ncRNA_cond_res <- as.data.frame(ncRNA_cond_res)
    ncRNA_cond_res <- ncRNA_cond_res |>
      filter((`STAAR-B(1,25)` < alpha_cond) | `STAAR-O` < alpha_cond)
    res_total_noncoding <- rbind(res_total_noncoding, ncRNA_cond_res)

  }
  save(res_total_noncoding,
       file = paste0(output_path, "Genecentric_noncoding_total.Rdata"))
  return(res_total_noncoding)
}


summary_chr6_noncoding_cond <- t(unlist(summary_chr6_noncoding_cond))
summary_chr6_noncoding_cond <- as.data.frame(summary_chr6_noncoding_cond)


summary_chr16_coding_cond_exp <- sapply(summary_chr16_coding_cond_exp, function(x) {
  if (is.list(x)) sapply(x, unlist) else x
})
summary_chr16_coding_cond_exp <- as.data.frame(summary_chr16_coding_cond_exp)

colnames(summary_chr6_noncoding_cond) <- colnames(summary_chr16_coding_cond_exp)
#colnames(summary_chr8_noncoding_cond) <- colnames(summary_chr16_coding_cond_exp)

all_rvs_cond_sig <- as.data.frame(matrix(unlist(rbind(summary_chr2_coding_cond,
                          summary_chr8_coding_cond,
                          summary_chr16_coding_cond_exp)), ncol = 91))
colnames(all_rvs_cond_sig) <- names(summary_chr2_coding_cond)
