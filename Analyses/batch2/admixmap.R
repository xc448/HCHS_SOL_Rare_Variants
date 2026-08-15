#==============================================================================
# Admixture mapping - BATCH 2
#==============================================================================
# MERGED FROM (chronological):
#   - 20250429_admixmap_multi_chr_b2.R   ->  section "v1_multi_chr"
#   - 20250920_admixmap_scale_up_b2.R   ->  section "v2_scaled_up"
# NOTE: b2 get_adj_cv() drops the `batch` argument (batch is implicit).
#==============================================================================

#==============================================================================
# MAIN PIPELINE  (current version: v2_scaled_up)
#==============================================================================

# load neccessary packages
library(tidyverse)
library(GENESIS)
library(gdsfmt)
library(GWASTools)
library(SNPRelate)
library(SeqArray)
library(reshape2)
library(STAARpipeline)
library(purrr)
library(tibble)

# set working directory
setwd("R:\\Sofer Lab\\HCHS_SOL\\Projects/2024_rare_variants/Code")


# subsetting overlapped individuals from DRAGEN gds and the local ancestry files
cov_b2 <- readRDS("../Data/STAAR_model_cov_batch2_scaled_up.RDS")
gdsflare <- openfn.gds(paste0("../../../Projects/2023_local_ancestry_comparison_sol/Data/FLARE3/FLARE3_SNPs_filtered_chr",
                              22, ".gds" ))
sample.id <- read.gdsn(index.gdsn(gdsflare, "sample.id"))
closefn.gds(gdsflare) # close the gdsfile

intersected_id_b2 <- intersect(paste0("SoL", sample.id), cov_b2$SUBJECT_ID)

rownames(cov_b2) <- cov_b2$SUBJECT_ID
cov_b2 <- cov_b2[intersected_id_b2,]

# load covariate matrices
rownames(cov_b2) <- cov_b2$NWD_ID #1465 for b1
covMatList <- readRDS("../Data/STAAR_model_covmat_batch2.RDS")
for(i in 1:length(covMatList)){
  covMatList[[i]] <- covMatList[[i]][rownames(cov_b2), rownames(cov_b2)]
}

# load gds files catalog that need to be accessed
load("../Data/STAAR_prep/agds_dir_scale_up_final_b2.Rdata") #variable name agds_dir

# define a function for getting pruned common variants selected from the STAAR
# pipeline and use it as a covariate to be adjusted in the admixutre mapping
get_adj_cv <- function(chr, known_cv, covariate_df, i){
  genofile <- seqOpen(agds_dir[i])
  position <- as.data.frame(seqGetData(genofile, "position"))
  variant_id <- seqGetData(genofile, "variant.id")
  variant_id <- variant_id[which(position[,1] %in% known_cv$POS)]
  print(known_cv[which(!known_cv$POS %in% position[,1]),])
  rownames(covariate_df) <- covariate_df$NWD_ID
  seqSetFilter(genofile, variant.id = variant_id, sample.id = rownames(covariate_df))
  # allele_count <- seqAlleleCount(genofile, ref.allele = known_cv$REF, minor = TRUE)
  # print(allele_count)
  # alelle_freq <- seqAlleleFreq(genofile,ref.allele = known_cv$REF, minor = TRUE)
  # print(alelle_freq)
  # allele <- seqGetData(genofile, "allele")
  genotype <- seqGetData(genofile, "genotype")
  sampleid <- seqGetData(genofile, "sample.id")
  dimnames(genotype)[[2]] <- sampleid

  num_snps <- dim(genotype)[3]
  # Extract SNPs (assuming diploid: using allele dosage for the second haplotype)
  snp_df <- map_dfc(1:num_snps, function(i) {
    snp_vals <- genotype[2, rownames(covariate_df), i]
    tibble(!!paste0("snp", i) := snp_vals)
  })

  covariate_df <- bind_cols(covariate_df, snp_df)

  seqClose(genofile)
  return(covariate_df)
}

# Function for running admixture mapping (single ancestry)
# leveraging local ancestry counts
# ancestry takes abbreviations such as "afr", "amer", or "eur"
runassoc_flare <- function(ancestry, nullmod, cov_df, chr){
  print(paste0("working on chromosome", chr))
  gdsflare_filtered <- openfn.gds(paste0("../../../Projects/2023_local_ancestry_comparison_sol/Data/FLARE3/FLARE3_SNPs_filtered_chr",
                                         chr, ".gds" ))
  sample.id <- read.gdsn(index.gdsn(gdsflare_filtered, "sample.id"))
  gds <- GdsGenotypeReader(gdsflare_filtered,  genotypeVar=paste0(ancestry, "_counts"))

  scanAnnot_flare <- ScanAnnotationDataFrame(data.frame(
    scanID=sample.id, stringsAsFactors=FALSE))
  genodata <- GenotypeData(gds, scanAnnot=scanAnnot_flare)
  snppos <- as.data.frame(getPosition(gds))
  rownames(snppos) <-  read.gdsn(index.gdsn(gdsflare_filtered, "snp.id"))
  geno <- getGenotype(genodata)
  rownames(geno) <- read.gdsn(index.gdsn(gdsflare_filtered, "snp.id"))
  closefn.gds(gdsflare_filtered)
  colnames(geno) <- paste0("SoL", sample.id)

  # Subset phenotype to match null model sample order
  rownames(cov_df) <- cov_df$NWD_ID
  cov_df <- cov_df[rownames(nullmod$fit), ]
  rownames(cov_df) <- cov_df$SUBJECT_ID
  # Filter genotype matrix to include only filtered IDs and matched phenotype samples
  geno <- geno[,
               rownames(cov_df)
  ]
  print(dim(geno))
  #print(sum(colnames(geno) != matched_solid)) # should be 0

  res_admixmap <- GENESIS:::testGenoSingleVar(nullmod, t(geno))
  res_admixmap$snppos <- snppos[,1]
  return(res_admixmap)
}

###### agds_dir[1], chr11 arachido,  x100009332 #####

known_loci_chr11_arachido <- readRDS("../Data/STAAR_exp_groups/x100009332/individual_cond_pruned_var.RDS")
cov_b2_chr11_arachido <- get_adj_cv(11, known_cv = known_loci_chr11_arachido,
                                    cov_b2, i = 1)

cov_b2_chr11_arachido <-  cov_b2_chr11_arachido|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1))

nullmod_b2_x100009332_no_genotype <- fitNullModel(x = cov_b2_chr11_arachido,
                                                  outcome="x100009332",
                                                  covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                           "PC4", "PC5",
                                                           "BKGRD1_C7", "GFRSCYS",
                                                           "CENTER"),
                                                  cov.mat=covMatList,
                                                  AIREML.tol=1e-4,
                                                  verbose=TRUE)

nullmod_b2_x100009332_cv_adj <- fitNullModel(x = cov_b2_chr11_arachido,
                                             outcome="x100009332",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS",
                                                      "CENTER", "snp1", "snp2",
                                                      "snp3", "snp4", "snp5",
                                                      "snp6", "snp7", "snp8",
                                                      "snp9"),
                                             cov.mat=covMatList,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)


admixmap_b2_x100009332_no_genotype <- runassoc_flare("amer", nullmod_b2_x100009332_no_genotype,
                                                     cov_df = cov_b2_chr11_arachido, chr = 11)
saveRDS(admixmap_b2_x100009332_no_genotype, file = "../Data/STAAR_exp_groups/x100009332/admixmap_b2_x100009332_no_genotype.RDS")
gc()

admixmap_b2_x100009332_cv_adj <- runassoc_flare("amer", nullmod_b2_x100009332_cv_adj,
                                                cov_df = cov_b2_chr11_arachido, chr = 11)
saveRDS(admixmap_b2_x100009332_cv_adj, "../Data/STAAR_exp_groups/x100009332/admixmap_b2_x100009332_cv_adj.RDS")
###### agds_dir[2], chr11 3beta,  x100006370 #####

known_loci_chr11_3beta <- readRDS("../Data/STAAR_exp_groups/x100006370/individual_cond_pruned_var.RDS")
cov_b2_chr11_3beta <- get_adj_cv(11, known_cv = known_loci_chr11_3beta,
                                 cov_b2, i = 2)
sum(is.na(cov_b2_chr11_3beta))

nullmod_b2_x100006370_no_genotype <- fitNullModel(x = cov_b2_chr11_3beta,
                                                  outcome="x100006370",
                                                  covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                           "PC4", "PC5",
                                                           "BKGRD1_C7", "GFRSCYS",
                                                           "CENTER"),
                                                  cov.mat=covMatList,
                                                  AIREML.tol=1e-4,
                                                  verbose=TRUE)

nullmod_b2_x100006370_cv_adj <- fitNullModel(x = cov_b2_chr11_3beta,
                                             outcome="x100006370",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS",
                                                      "CENTER", "snp1", "snp2",
                                                      "snp3", "snp4", "snp5"),
                                             cov.mat=covMatList,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)

admixmap_b2_x100006370_no_genotype <- runassoc_flare("amer", nullmod_b2_x100006370_no_genotype,
                                                     cov_df = cov_b2_chr11_3beta, chr = 11)
saveRDS(admixmap_b2_x100006370_no_genotype, file = "../Data/STAAR_exp_groups/x100006370/admixmap_b2_x100006370_no_genotype.RDS")
admixmap_b2_x100006370_cv_adj <- runassoc_flare("amer", nullmod_b2_x100006370_cv_adj,
                                                cov_df = cov_b2_chr11_3beta, chr = 11)
saveRDS(admixmap_b2_x100006370_cv_adj, "../Data/STAAR_exp_groups/x100006370/admixmap_b2_x100006370_cv_adj.RDS")
###### agds_dir[3], chr12 1_methyl, x100001208 #####

known_loci_chr12_1_methyl <- readRDS("../Data/STAAR_exp_groups/x100001208/individual_cond_pruned_var.RDS")
cov_b2_chr12_1_methyl <- get_adj_cv(12, known_cv = known_loci_chr12_1_methyl,
                                    cov_b2, i = 3)

nullmod_b2_x100001208_no_genotype <- fitNullModel(x = cov_b2_chr12_1_methyl,
                                                  outcome="x100001208",
                                                  covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                           "PC4", "PC5",
                                                           "BKGRD1_C7", "GFRSCYS",
                                                           "CENTER"),
                                                  cov.mat=covMatList,
                                                  AIREML.tol=1e-4,
                                                  verbose=TRUE)

nullmod_b2_x100001208_cv_adj <- fitNullModel(x = cov_b2_chr12_1_methyl,
                                             outcome="x100001208",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS",
                                                      "CENTER", "snp1", "snp2",
                                                      "snp3", "snp4", "snp5",
                                                      "snp6"),
                                             cov.mat=covMatList,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)

admixmap_b2_x100001208_no_genotype <- runassoc_flare("amer", nullmod_b2_x100001208_no_genotype,
                                                     cov_df = cov_b2_chr12_1_methyl, chr = 12)
saveRDS(admixmap_b2_x100001208_no_genotype, file = "../Data/STAAR_exp_groups/x100001208/admixmap_b2_x100001208_no_genotype.RDS")
admixmap_b2_x100001208_cv_adj <- runassoc_flare("amer", nullmod_b2_x100001208_cv_adj,
                                                cov_df = cov_b2_chr12_1_methyl, chr = 12)
saveRDS(admixmap_b2_x100001208_cv_adj, "../Data/STAAR_exp_groups/x100001208/admixmap_b2_x100001208_cv_adj.RDS")
###### agds_dir[4], chr12 ethylmalonate, x2054 #####

known_loci_chr12_ethylmalonate <- readRDS("../Data/STAAR_exp_groups/x2054/individual_cond_pruned_var.RDS")
cov_b2_chr12_1_ethylmalonate <- get_adj_cv(12, known_cv = known_loci_chr12_1_ethylmalonate,
                                           cov_b2, i = 4)


nullmod_b2_x2054_no_genotype <- fitNullModel(x = cov_b2_chr12_1_ethylmalonate,
                                             outcome="x2054",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS",
                                                      "CENTER"),
                                             cov.mat=covMatList,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)

nullmod_b2_x2054_cv_adj <- fitNullModel(x = cov_b2_chr12_1_ethylmalonate,
                                        outcome="x2054",
                                        covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                 "PC4", "PC5",
                                                 "BKGRD1_C7", "GFRSCYS",
                                                 "CENTER", "snp1", "snp2"),
                                        cov.mat=covMatList,
                                        AIREML.tol=1e-4,
                                        verbose=TRUE)

admixmap_b2_x2054_no_genotype <- runassoc_flare("afr", nullmod_b2_x2054_no_genotype,
                                                cov_df = cov_b2_chr12_1_ethylmalonate, chr = 12)
saveRDS(admixmap_b2_x2054_no_genotype, file = "../Data/STAAR_exp_groups/x2054/admixmap_b2_x2054_no_genotype.RDS")
admixmap_b2_x2054_cv_adj <- runassoc_flare("afr", nullmod_b2_x2054_cv_adj,
                                           cov_df = cov_b2_chr12_1_ethylmalonate, chr = 12)
saveRDS(admixmap_b2_x2054_cv_adj, "../Data/STAAR_exp_groups/x2054/admixmap_b2_x2054_cv_adj.RDS")
#### Rare variants #####
chr12_plof_ds_burden_b2 <- readRDS("../Data/STAAR_exp_groups/x2054/b2/burden_plof_ds.RDS")

intersected_samples <- intersect(rownames(cov_b2_chr12_1_ethylmalonate), rownames(chr12_plof_ds_burden_b2))
cov_b2_chr12_plof_ds <- cbind(cov_b2_chr12_1_ethylmalonate[intersected_samples, ],
                              chr12_plof_ds_burden_b2[intersected_samples, ])

covMatList_rv <- covMatList
for(i in 1:length(covMatList_rv)){
  covMatList_rv[[i]] <- covMatList_rv[[i]][rownames(cov_b2_chr12_plof_ds),
                                           rownames(cov_b2_chr12_plof_ds)]
}

colnames(cov_b2_chr12_plof_ds)[ncol(cov_b2_chr12_plof_ds)]  <- "ACADS_plof_ds"

nullmod_b2_x2054_plof_ds_adj <- fitNullModel(x = cov_b2_chr12_plof_ds,
                                             outcome="x2054",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS",
                                                      "CENTER", "ACADS_plof_ds"),
                                             cov.mat=covMatList_rv,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)

admixmap_b2_x2054_plof_ds_adj <- runassoc_flare("afr", nullmod_b2_x2054_plof_ds_adj,
                                                cov_df = cov_b2_chr12_plof_ds,
                                                chr = 12)
saveRDS(admixmap_b2_x2054_plof_ds_adj, "../Data/STAAR_exp_groups/x2054/admixmap_b2_x2054_plof_ds_adj.RDS")
#### adjusting common variants as well ####

nullmod_b2_x2054_plof_ds_cv_adj <- fitNullModel(x = cov_b2_chr12_plof_ds,
                                             outcome="x2054",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS",
                                                      "CENTER", "ACADS_plof_ds",
                                                      "snp1",
                                                      "snp2"),
                                             cov.mat=covMatList_rv,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)

admixmap_b2_x2054_plof_ds_cv_adj <- runassoc_flare("afr", nullmod_b2_x2054_plof_ds_cv_adj,
                                                cov_df = cov_b2_chr12_plof_ds,
                                                chr = 12)
saveRDS(admixmap_b2_x2054_plof_ds_cv_adj, "../Data/STAAR_exp_groups/x2054/admixmap_b2_x2054_rv_cv_adj.RDS")
###### agds_dir[5], chr16  #####

known_loci_chr16_propyl <- readRDS("../Data/STAAR_exp_groups/x100006264/individual_cond_pruned_var.RDS")
cov_b2_chr16_propyl <- get_adj_cv(16, known_cv = known_loci_chr16_propyl,
                                  cov_b2, i = 5)
nullmod_b2_x100006264_no_genotype <- fitNullModel(x = cov_b2_chr16_propyl,
                                                  outcome="x100006264",
                                                  covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                           "PC4", "PC5",
                                                           "BKGRD1_C7", "GFRSCYS",
                                                           "CENTER"),
                                                  cov.mat=covMatList,
                                                  AIREML.tol=1e-4,
                                                  verbose=TRUE)

cov_b2_chr16_propyl <- cov_b2_chr16_propyl |>
  mutate(snp7 = if_else(is.na(snp7), mean(snp7, na.rm = TRUE), snp7),
          snp8 = if_else(is.na(snp8), mean(snp8, na.rm = TRUE), snp8),
         snp9 = if_else(is.na(snp9), mean(snp9, na.rm = TRUE), snp9),
         snp10 = if_else(is.na(snp10), mean(snp10, na.rm = TRUE), snp10),
         snp12 = if_else(is.na(snp12), mean(snp12, na.rm = TRUE), snp12),
         snp17 = if_else(is.na(snp17), mean(snp17, na.rm = TRUE), snp17))


nullmod_b2_x100006264_cv_adj <- fitNullModel(x = cov_b2_chr16_propyl,
                                             outcome="x100006264",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS",
                                                      "CENTER", "snp1", "snp2",
                                                      "snp3", "snp4", "snp5",
                                                      "snp6", "snp8",
                                                      "snp9", "snp10", "snp11",
                                                      "snp12", "snp13", "snp14",
                                                      "snp15", "snp16", "snp17"),
                                             cov.mat=covMatList,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)

admixmap_b2_x100006264_no_genotype <- runassoc_flare("afr", nullmod_b2_x100006264_no_genotype,
                                                     cov_df = cov_b2_chr16_propyl, chr = 16)
saveRDS(admixmap_b2_x100006264_no_genotype, file = "../Data/STAAR_exp_groups/x100006264/admixmap_b2_x100006264_no_genotype.RDS")
admixmap_b2_x100006264_cv_adj <- runassoc_flare("afr", nullmod_b2_x100006264_cv_adj,
                                                cov_df = cov_b2_chr16_propyl, chr = 16)
saveRDS(admixmap_b2_x100006264_cv_adj, "../Data/STAAR_exp_groups/x100006264/admixmap_b2_x100006264_cv_adj.RDS")
##### adjusting rare variant sets chr 16 b2 ######
chr16_burden_comined_b2 <- readRDS("../Data/STAAR_exp_groups/x100006264/b2/burdens_combined.RDS")

intersected_samples <- intersect(rownames(cov_b2_chr16_propyl), rownames(chr16_burden_comined_b2))
cov_b2_chr16_rvs <- cbind(cov_b2_chr16_propyl[intersected_samples, ],
                          chr16_burden_comined_b2[intersected_samples, ])
# check NAs - 0
sum(is.na(cov_b2_chr16_rvs))

for(i in 1:length(covMatList_rv)){
  covMatList_rv[[i]] <- covMatList_rv[[i]][rownames(cov_b2_chr16_rvs),
                                           rownames(cov_b2_chr16_rvs)]
}

nullmod_b2_x100006264_rvs_adj <- fitNullModel(x = cov_b2_chr16_rvs,
                                              outcome="x100006264",
                                              covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                       "PC4", "PC5",
                                                       "BKGRD1_C7", "GFRSCYS", "FAM57B_enhancer_DHS",
                                                       "CENTER", "PRRT2_enhancer_DHS",
                                                       "RNF40_UTR"),
                                              cov.mat=covMatList_rv,
                                              AIREML.tol=1e-4,
                                              verbose=TRUE)

admixmap_b2_x100006264_rvs_adj <- runassoc_flare("afr", nullmod_b2_x100006264_rvs_adj,
                                                 cov_df = cov_b2_chr16_rvs,
                                                 chr = 16)

saveRDS(admixmap_b2_x100006264_rvs_adj, "../Data/STAAR_exp_groups/x100006264/admixmap_b2_x100006264_rv_adj.RDS")
############## RARE VARIANTS BURDEN + COMMON VARIANTS ADJUSTING ################

nullmod_b2_x100006264_rvs_cvs_adj <- fitNullModel(x = cov_b2_chr16_rvs,
                                                  outcome="x100006264",
                                                  covars=c("AGE","GENDER","PC1",
                                                           "PC2", "PC3",
                                                           "PC4", "PC5",
                                                           "BKGRD1_C7", "GFRSCYS",
                                                           "CENTER", "PRRT2_enhancer_DHS",
                                                           "RNF40_UTR", "FAM57B_enhancer_DHS",
                                                           "snp1", "snp2",
                                                           "snp3", "snp4", "snp5",
                                                           "snp6", "snp8",
                                                           "snp9", "snp10", "snp11",
                                                           "snp12", "snp13", "snp14",
                                                           "snp15", "snp16", "snp17"),
                                                  cov.mat=covMatList_rv,
                                                  AIREML.tol=1e-4,
                                                  verbose=TRUE)

admixmap_b2_x100006264_rvs_cvs_adj <- runassoc_flare("afr",
                                                     nullmod_b2_x100006264_rvs_cvs_adj,
                                                     cov_df = cov_b2_chr16_rvs,
                                                     chr = 16)

saveRDS(admixmap_b2_x100006264_rvs_cvs_adj,
        "../Data/STAAR_exp_groups/x100006264/admixmap_b2_x100006264_rv_cv_adj.RDS")
###### agds_dir[6], chr2  #####

known_loci_chr2 <- readRDS("../Data/STAAR_exp_groups/x100001266/individual_cond_pruned_var.RDS")
cov_b2_known_loci_chr2 <- get_adj_cv(2, known_cv = known_loci_chr2,
                                     cov_b2, i = 6)

# check NAs
sum(is.na(cov_b2_known_loci_chr2))

nullmod_b2_x100001266_no_genotype <- fitNullModel(x = cov_b2_known_loci_chr2,
                                                  outcome="x100001266",
                                                  covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                           "PC4", "PC5",
                                                           "BKGRD1_C7", "GFRSCYS",
                                                           "CENTER"),
                                                  cov.mat=covMatList,
                                                  AIREML.tol=1e-4,
                                                  verbose=TRUE)

nullmod_b2_x100001266_cv_adj <- fitNullModel(x = cov_b2_known_loci_chr2,
                                             outcome="x100001266",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS",
                                                      "CENTER", "snp1", "snp2",
                                                      "snp3", "snp4",
                                                      "snp5", "snp6", "snp7",
                                                      "snp8", "snp9", "snp10"),
                                             cov.mat=covMatList,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)

admixmap_b2_x100001266_no_genotype <- runassoc_flare("afr", nullmod_b2_x100001266_no_genotype,
                                                     cov_df = cov_b2_known_loci_chr2, chr = 2)
saveRDS(admixmap_b2_x100001266_no_genotype,
        file = "../Data/STAAR_exp_groups/x100001266/admixmap_b2_x100001266_no_genotype.RDS")
admixmap_b2_x100001266_cv_adj <- runassoc_flare("afr", nullmod_b2_x100001266_cv_adj,
                                                cov_df = cov_b2_known_loci_chr2, chr = 2)
saveRDS(admixmap_b2_x100001266_cv_adj, "../Data/STAAR_exp_groups/x100001266/admixmap_b2_x100001266_cv_adj.RDS")
###### agds_dir[7], chr5  #####

known_loci_chr5 <- readRDS("../Data/STAAR_exp_groups/x1114/individual_cond_pruned_var.RDS")
cov_b2_known_loci_chr5 <- get_adj_cv(5, known_cv = known_loci_chr5,
                                     cov_b2, i = 7)


nullmod_b2_x1114_no_genotype <- fitNullModel(x = cov_b2_known_loci_chr5,
                                             outcome="x1114",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS",
                                                      "CENTER"),
                                             cov.mat=covMatList,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)

nullmod_b2_x1114_cv_adj <- fitNullModel(x = cov_b2_known_loci_chr5,
                                        outcome="x1114",
                                        covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                 "PC4", "PC5",
                                                 "BKGRD1_C7", "GFRSCYS",
                                                 "CENTER", "snp1", "snp2",
                                                 "snp3", "snp4",
                                                 "snp5", "snp6", "snp7",
                                                 "snp8", "snp9", "snp10"),
                                        cov.mat=covMatList,
                                        AIREML.tol=1e-4,
                                        verbose=TRUE)

admixmap_b2_x1114_no_genotype <- runassoc_flare("amer", nullmod_b2_x1114_no_genotype,
                                                cov_df = cov_b2_known_loci_chr5, chr = 5)
saveRDS(admixmap_b2_x1114_no_genotype, file = "../Data/STAAR_exp_groups/x1114/admixmap_b2_x1114_no_genotype.RDS")
admixmap_b2_x1114_cv_adj <- runassoc_flare("amer", nullmod_b2_x1114_cv_adj,
                                           cov_df = cov_b2_known_loci_chr5, chr = 5)
saveRDS(admixmap_b2_x1114_cv_adj, "../Data/STAAR_exp_groups/x1114/admixmap_b2_x1114_cv_adj.RDS")
####### Rare variant plof_ds #######
chr5_burden_combined_b2 <- readRDS("../Data/STAAR_exp_groups/x1114/b2/burden_plof_ds.RDS")

intersected_samples <- intersect(rownames(cov_b2_known_loci_chr5), rownames(chr5_burden_combined_b2))
cov_b2_chr5_rv_plof_ds <- cbind(cov_b2_known_loci_chr5[intersected_samples, ],
                                chr5_burden_combined_b2[intersected_samples, ])

colnames(cov_b2_chr5_rv_plof_ds)[39] <- "AGXT2_plof_ds"

for(i in 1:length(covMatList_rv)){
  covMatList_rv[[i]] <- covMatList_rv[[i]][rownames(cov_b2_chr5_rv_plof_ds),
                                           rownames(cov_b2_chr5_rv_plof_ds)]
}

nullmod_b2_x1114_plof_ds_adj <- fitNullModel(x = cov_b2_chr5_rv_plof_ds,
                                             outcome="x1114",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS",
                                                      "CENTER",
                                                      "AGXT2_plof_ds"),
                                             cov.mat=covMatList_rv,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)

admixmap_b2_x1114_plof_ds_adj <- runassoc_flare("amer", nullmod_b2_x1114_plof_ds_adj,
                                                cov_df = cov_b2_chr5_rv_plof_ds,
                                                chr = 5)
saveRDS(admixmap_b2_x1114_plof_ds_adj, "../Data/STAAR_exp_groups/x1114/admixmap_b2_x1114_plof_ds_adj.RDS")
###### Rare + common ######
nullmod_b2_x1114_rv_cv_adj <- fitNullModel(x = cov_b2_chr5_rv_plof_ds,
                                           outcome="x1114",
                                           covars=c("AGE","GENDER","PC1",
                                                    "PC2", "PC3",
                                                    "PC4", "PC5",
                                                    "BKGRD1_C7", "GFRSCYS",
                                                    "CENTER",
                                                    "snp1", "snp2",
                                                    "snp3", "snp4",
                                                    "snp5", "snp6", "snp7",
                                                    "snp8", "snp9", "snp10",
                                                    "AGXT2_plof_ds"),
                                           cov.mat=covMatList_rv,
                                           AIREML.tol=1e-4,
                                           verbose=TRUE)

admixmap_b2_x1114_rv_cv_adj <- runassoc_flare("amer", nullmod_b2_x1114_rv_cv_adj,
                                              cov_df = cov_b2_chr5_rv_plof_ds,
                                              chr = 5)
saveRDS(admixmap_b2_x1114_rv_cv_adj, "../Data/STAAR_exp_groups/x1114/admixmap_b2_x1114_rv_cv_adj.RDS")
###### agds_dir[8], chr8  #####

known_loci_chr8 <- readRDS("../Data/STAAR_exp_groups/x192/individual_cond_pruned_var.RDS")
cov_b2_known_loci_chr8 <- get_adj_cv(8, known_cv = known_loci_chr8,
                                     cov_b2, i = 8)

nullmod_b2_x192_no_genotype <- fitNullModel(x = cov_b2_known_loci_chr8,
                                            outcome="x192",
                                            covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                     "PC4", "PC5",
                                                     "BKGRD1_C7", "GFRSCYS",
                                                     "CENTER"),
                                            cov.mat=covMatList,
                                            AIREML.tol=1e-4,
                                            verbose=TRUE)

nullmod_b2_x192_cv_adj <- fitNullModel(x = cov_b2_known_loci_chr8,
                                       outcome="x192",
                                       covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                "PC4", "PC5",
                                                "BKGRD1_C7", "GFRSCYS",
                                                "CENTER", "snp1", "snp2",
                                                "snp3", "snp4",
                                                "snp5", "snp6"),
                                       cov.mat=covMatList,
                                       AIREML.tol=1e-4,
                                       verbose=TRUE)

admixmap_b2_x192_no_genotype <- runassoc_flare("amer", nullmod_b2_x192_no_genotype,
                                               cov_df = cov_b2_known_loci_chr8, chr = 8)
saveRDS(admixmap_b2_x192_no_genotype, file = "../Data/STAAR_exp_groups/x192/admixmap_b2_x192_no_genotype.RDS")
admixmap_b2_x192_cv_adj <- runassoc_flare("amer", nullmod_b2_x192_cv_adj,
                                          cov_df = cov_b2_known_loci_chr8, chr = 8)
saveRDS(admixmap_b2_x192_cv_adj, "../Data/STAAR_exp_groups/x192/admixmap_b2_x192_cv_adj.RDS")


#------------------------------------------------------------------------------
# LEGACY / EARLIER VERSION: v1_multi_chr  (138 unique blocks)
#------------------------------------------------------------------------------
# Kept verbatim. These are blocks that do NOT appear in the current version above
# (mostly hard-coded per-metabolite / per-chromosome run calls and older path setups).

#library(STAARpipeline)
setwd("R:/HCHS_SOL/Projects/2024_rare_variants/Code")

cov_b2 <- readRDS("../Data/STAAR_model_cov_batch2.RDS")

saveRDS(intersected_id_b2, "../Data/SoL_intersected_ID_b2.RDS")
rownames(cov_b2) <- cov_b2$NWD_ID
get_adj_cv <- function(chr, batch, known_cv, covariate_df){
  genofile <- seqOpen(paste0("../Data/chr", chr, "_DRAGEN_agds_",
                             batch, "_final.gds"))
  position <- as.data.frame(seqGetData(genofile, "position"))
  variant_id <- seqGetData(genofile, "variant.id")
  variant_id <- variant_id[which(position[,1] %in% known_cv$POS)]
  rownames(covariate_df) <- covariate_df$NWD_ID
  seqSetFilter(genofile, variant.id = variant_id, sample.id = rownames(covariate_df))
  # allele_count <- seqAlleleCount(genofile, ref.allele = known_cv$REF, minor = TRUE)
  # print(allele_count)
  # alelle_freq <- seqAlleleFreq(genofile,ref.allele = known_cv$REF, minor = TRUE)
  # print(alelle_freq)
  # allele <- seqGetData(genofile, "allele")
  genotype <- seqGetData(genofile, "genotype")
  sampleid <- seqGetData(genofile, "sample.id")
  dimnames(genotype)[[2]] <- sampleid

  num_snps <- dim(genotype)[3]
  # Extract SNPs (assuming diploid: using allele dosage for the second haplotype)
  snp_df <- map_dfc(1:num_snps, function(i) {
    snp_vals <- genotype[2, rownames(covariate_df), i]
    tibble(!!paste0("snp", i) := snp_vals)
  })

  covariate_df <- bind_cols(covariate_df, snp_df)

  seqClose(genofile)
  return(covariate_df)
}


runassoc_flare <- function(ancestry, nullmod, cov_df, chr){
  print(paste0("working on chromosome", chr))
  gdsflare_filtered <- openfn.gds(paste0("../../../Projects/2023_local_ancestry_comparison_sol/Data/FLARE3/FLARE3_SNPs_filtered_chr",
                                         chr, ".gds" ))
  sample.id <- read.gdsn(index.gdsn(gdsflare_filtered, "sample.id"))
  gds <- GdsGenotypeReader(gdsflare_filtered,  genotypeVar=paste0(ancestry, "_counts"))

  scanAnnot_flare <- ScanAnnotationDataFrame(data.frame(
    scanID=sample.id, stringsAsFactors=FALSE))
  genodata <- GenotypeData(gds, scanAnnot=scanAnnot_flare)
  chrom <- getChromosome(gds)
  snppos <- as.data.frame(getPosition(gds))
  rownames(snppos) <-  read.gdsn(index.gdsn(gdsflare_filtered, "snp.id"))
  geno <- getGenotype(genodata)
  rownames(geno) <- read.gdsn(index.gdsn(gdsflare_filtered, "snp.id"))
  closefn.gds(gdsflare_filtered)
  colnames(geno) <- paste0("SoL", sample.id)

  # Subset phenotype to match null model sample order
  rownames(cov_df) <- cov_df$NWD_ID
  cov_df <- cov_df[rownames(nullmod$fit), ]
  rownames(cov_df) <- cov_df$SUBJECT_ID
  # Filter genotype matrix to include only filtered IDs and matched phenotype samples
  geno <- geno[,
               rownames(cov_df)
  ]
  print(dim(geno))
  #print(sum(colnames(geno) != matched_solid)) # should be 0

  res_admixmap <- GENESIS:::testGenoSingleVar(nullmod, t(geno))
  #saveRDS(t(geno), "../Data/STAAR_chr16_1224/admixmap_local_ancestry_counts_amer.RDS")
  res_admixmap$snppos <- snppos[,1]
  return(res_admixmap)
}

################ chr2 ################
known_loci_chr2 <- as.data.frame(readRDS("../Data/STAAR_neg_controls/x100001721/individual_cond_pruned_var.RDS"))
cov_b2_chr2 <- get_adj_cv(2, "b2", known_cv = known_loci_chr2,
                          cov_b2)

cov_b2_chr2 <-  cov_b2_chr2|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2),
         snp3 = if_else(is.na(snp3), mean(snp2, na.rm = TRUE), snp3),
         snp4 = if_else(is.na(snp4), mean(snp4, na.rm = TRUE), snp4),
         snp5 = if_else(is.na(snp5), mean(snp5, na.rm = TRUE), snp5),
         snp6 = if_else(is.na(snp6), mean(snp6, na.rm = TRUE), snp6),
         snp7 = if_else(is.na(snp7), mean(snp7, na.rm = TRUE), snp7),
         GFRSCYS = if_else(is.na(GFRSCYS), mean(GFRSCYS, na.rm = TRUE), GFRSCYS))

nullmod_b2_x100001721_no_genotype <- fitNullModel(x = cov_b2_chr2,
                                                  outcome="x100001721",
                                                  covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                           "PC4", "PC5",
                                                           "BKGRD1_C7", "GFRSCYS",
                                                           "CENTER"),
                                                  cov.mat=covMatList,
                                                  AIREML.tol=1e-4,
                                                  verbose=TRUE)


admixmap_b2_x100001721_no_genotype <- runassoc_flare("afr", nullmod_b2_x100001721_no_genotype,
                                                     cov_df = cov_b2_chr2, chr = 2)
saveRDS(admixmap_b2_x100001721_no_genotype, "../Data/STAAR_neg_controls/x100001721/admixmap_b2_x100001721_no_genotype.RDS")

nullmod_b2_x100001721_cv_adj <- fitNullModel(x = cov_b2_chr2,
                                              outcome="x100001721",
                                              covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                       "PC4", "PC5",
                                                       "BKGRD1_C7", "GFRSCYS",
                                                       "CENTER", "snp1", "snp2",
                                                       "snp3", "snp4", "snp5",
                                                       "snp6", "snp7"),
                                              cov.mat=covMatList,
                                              AIREML.tol=1e-4,
                                              verbose=TRUE)

admixmap_b2_x100001721_cv_adj <- runassoc_flare("afr", nullmod_b2_x100001721_cv_adj,
                                                 cov_df = cov_b2_chr2,
                                                 chr = 2)
saveRDS(admixmap_b2_x100001721_cv_adj, "../Data/STAAR_neg_controls/x100001721/admixmap_b2_x100001721_cv_adj.RDS")


################ chr5 ################
known_loci_chr5 <- readRDS("../Data/STAAR_neg_controls/x799/individual_cond_pruned_var.RDS")
cov_b2_chr5 <- get_adj_cv(5, "b2", known_cv = known_loci_chr5,
                          cov_b2)

cov_b2_chr5 <-  cov_b2_chr5|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         GFRSCYS = if_else(is.na(GFRSCYS), mean(GFRSCYS, na.rm = TRUE), GFRSCYS))

nullmod_b2_x799_no_genotype <- fitNullModel(x = cov_b2_chr5,
                                            outcome="x799",
                                            covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                     "PC4", "PC5",
                                                     "BKGRD1_C7", "GFRSCYS",
                                                     "CENTER"),
                                            cov.mat=covMatList,
                                            AIREML.tol=1e-4,
                                            verbose=TRUE)
admixmap_b2_x799_no_genotype <- runassoc_flare("afr", nullmod_b2_x799_no_genotype,
                                               cov_df = cov_b2_chr5, chr = 5)
saveRDS(admixmap_b2_x799_no_genotype, "../Data/STAAR_neg_controls/x799/admixmap_b2_x799_no_genotype.RDS")

nullmod_b2_x799_cv_adj <- fitNullModel(x = cov_b2_chr5,
                                        outcome="x799",
                                        covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                 "PC4", "PC5",
                                                 "BKGRD1_C7", "GFRSCYS",
                                                 "CENTER", "snp1", "snp2"),
                                        cov.mat=covMatList,
                                        AIREML.tol=1e-4,
                                        verbose=TRUE)

admixmap_b2_x799_cv_adj <- runassoc_flare("afr", nullmod_b2_x799_cv_adj,
                                           cov_df = cov_b2_chr5,
                                           chr = 5)
saveRDS(admixmap_b2_x799_cv_adj, "../Data/STAAR_neg_controls/x799/admixmap_b2_x799_cv_adj.RDS")

############## chr6 ###############
known_loci_chr6 <- readRDS("../Data/STAAR_neg_controls/x1215/individual_cond_pruned_var.RDS")
cov_b2_chr6 <- get_adj_cv(6, "b2", known_cv = known_loci_chr6,
                          cov_b2)

cov_b2_chr6 <-  cov_b2_chr6|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2),
         GFRSCYS = if_else(is.na(GFRSCYS), mean(GFRSCYS, na.rm = TRUE), GFRSCYS))


nullmod_b2_x1215_no_genotype <- fitNullModel(x = cov_b2_chr6,
                                             outcome="x1215",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS",
                                                      "CENTER"),
                                             cov.mat=covMatList,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)

admixmap_b2_x1215_no_genotype <- runassoc_flare("amer", nullmod_b2_x1215_no_genotype,
                                                cov_df = cov_b2_chr6, chr = 6)
saveRDS(admixmap_b2_x1215_no_genotype,
        "../Data/STAAR_neg_controls/x1215/admixmap_b2_x1215_no_genotype.RDS")


nullmod_b2_x1215_cv_adj <- fitNullModel(x = cov_b2_chr6,
                                         outcome="x1215",
                                         covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                  "PC4", "PC5",
                                                  "BKGRD1_C7", "GFRSCYS",
                                                  "CENTER", "snp1", "snp2"),
                                         cov.mat=covMatList,
                                         AIREML.tol=1e-4,
                                         verbose=TRUE)

admixmap_b2_x1215_cv_adj <- runassoc_flare("amer", nullmod_b2_x1215_cv_adj,
                                            cov_df = cov_b2_chr6,
                                            chr = 6)
saveRDS(admixmap_b2_x1215_cv_adj, "../Data/STAAR_neg_controls/x1215/admixmap_b2_x1215_cv_adj.RDS")


############## chr16 neg control ###############
known_loci_chr16 <- readRDS("../Data/STAAR_neg_controls/x278/individual_cond_pruned_var.RDS")
cov_b2_chr16 <- get_adj_cv(16, "b2", known_cv = known_loci_chr16,
                           cov_b2)

cov_b2_chr16 <-  cov_b2_chr16|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2),
         snp3 = if_else(is.na(snp3), mean(snp3, na.rm = TRUE), snp3),
         GFRSCYS = if_else(is.na(GFRSCYS), mean(GFRSCYS, na.rm = TRUE), GFRSCYS))


nullmod_b2_x278_no_genotype <- fitNullModel(x = cov_b2_chr16,
                                            outcome="x278",
                                            covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                     "PC4", "PC5",
                                                     "BKGRD1_C7", "GFRSCYS",
                                                     "CENTER"),
                                            cov.mat=covMatList,
                                            AIREML.tol=1e-4,
                                            verbose=TRUE)

admixmap_b2_x278_no_genotype <- runassoc_flare("amer", nullmod_b2_x278_no_genotype,
                                               cov_df = cov_b2_chr16, chr = 16)
saveRDS(admixmap_b2_x278_no_genotype,
        "../Data/STAAR_neg_controls/x278/admixmap_b2_x278_no_genotype.RDS")


nullmod_b2_x278_cv_adj <- fitNullModel(x = cov_b2_chr16,
                                        outcome="x278",
                                        covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                 "PC4", "PC5",
                                                 "BKGRD1_C7", "GFRSCYS",
                                                 "CENTER", "snp1", "snp2",
                                                 "snp3"),
                                        cov.mat=covMatList,
                                        AIREML.tol=1e-4,
                                        verbose=TRUE)

admixmap_b2_x278_cv_adj <- runassoc_flare("amer", nullmod_b2_x278_cv_adj,
                                           cov_df = cov_b2_chr16,
                                           chr = 16)
saveRDS(admixmap_b2_x278_cv_adj, "../Data/STAAR_neg_controls/x278/admixmap_b2_x278_cv_adj.RDS")


################### EXPERIMENTAL GROUP -- UNEXPLAINED SIGNALS ############
known_loci_chr8 <- readRDS("../Data/STAAR_exp_groups/x1021/individual_cond_pruned_var.RDS")
cov_b2_chr8 <- get_adj_cv(8, "b2", known_cv = known_loci_chr8,
                          cov_b2)

cov_b2_chr8 <-  cov_b2_chr8|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2),
         snp3 = if_else(is.na(snp3), mean(snp3, na.rm = TRUE), snp3),
         snp4 = if_else(is.na(snp4), mean(snp4, na.rm = TRUE), snp4),
         snp5 = if_else(is.na(snp5), mean(snp5, na.rm = TRUE), snp5),
         snp6 = if_else(is.na(snp6), mean(snp6, na.rm = TRUE), snp6),
         snp7 = if_else(is.na(snp7), mean(snp7, na.rm = TRUE), snp7),
         snp8 = if_else(is.na(snp8), mean(snp8, na.rm = TRUE), snp8),
         GFRSCYS = if_else(is.na(GFRSCYS), mean(GFRSCYS, na.rm = TRUE), GFRSCYS))


nullmod_b2_x1021_no_genotype <- fitNullModel(x = cov_b2_chr8,
                                             outcome="x1021",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS",
                                                      "CENTER"),
                                             cov.mat=covMatList,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)

nullmod_b2_x1021_cv_adj <- fitNullModel(x = cov_b2_chr8,
                                        outcome="x1021",
                                        covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                 "PC4", "PC5",
                                                 "BKGRD1_C7", "GFRSCYS",
                                                 "CENTER", "snp1", "snp2", "snp3",
                                                 "snp4", "snp5", "snp6", "snp7",
                                                 "snp8"),
                                        cov.mat=covMatList,
                                        AIREML.tol=1e-4,
                                        verbose=TRUE)

admixmap_b2_x1021_no_genotype <- runassoc_flare("afr", nullmod_b2_x1021_no_genotype,
                                                cov_df = cov_b2_chr8, chr = 8)
saveRDS(admixmap_b2_x1021_no_genotype,
        "../Data/STAAR_exp_groups/x1021/admixmap_b2_x1021_no_genotype.RDS")

admixmap_b2_x1021_cv_adj <- runassoc_flare("afr", nullmod_b2_x1021_cv_adj,
                                           cov_df = cov_b2_chr8, chr = 8)
saveRDS(admixmap_b2_x1021_cv_adj,
        "../Data/STAAR_exp_groups/x1021/admixmap_b2_x1021_cv_adj.RDS")


#### rare variants ####

chr8_plof_ds_burden <- readRDS("../Data/STAAR_exp_groups/x1021/b2/burden_plof_ds.RDS")

intersected_samples <- intersect(rownames(cov_b2_chr8), rownames(chr8_plof_ds_burden))
cov_b2_chr8_plof_ds <- cbind(cov_b2_chr8[intersected_samples, ],
                                       chr8_plof_ds_burden[intersected_samples, ])

colnames(cov_b2_chr8_plof_ds)[ncol(cov_b2_chr8_plof_ds)] <- "Burden_1_25"

for(i in 1:length(covMatList_rv)){
  covMatList_rv[[i]] <- covMatList_rv[[i]][rownames(cov_b2_chr8_plof_ds),
                                           rownames(cov_b2_chr8_plof_ds)]
}

nullmod_b2_x1021_rv_adj <- fitNullModel(x = cov_b2_chr8_plof_ds,
                                             outcome="x1021",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS",
                                                      "CENTER", "Burden_1_25"),
                                             cov.mat=covMatList_rv,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)

admixmap_b2_x1021_rv_adj <- runassoc_flare("afr", nullmod_b2_x1021_rv_adj,
                                                cov_df = cov_b2_chr8_plof_ds,
                                                chr = 8)
saveRDS(admixmap_b2_x1021_rv_adj, "../Data/STAAR_exp_groups/x1021/admixmap_b2_x1021_rv_adj.RDS")
nullmod_b2_x1021_rv_cv_adj <- fitNullModel(x = cov_b2_chr8_plof_ds,
                                                outcome="x1021",
                                                covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                         "PC4", "PC5",
                                                         "BKGRD1_C7", "GFRSCYS",
                                                         "CENTER", "Burden_1_25",
                                                         "snp1", "snp2", "snp3",
                                                         "snp4", "snp5", "snp6",
                                                         "snp7", "snp8"),
                                                cov.mat=covMatList_rv,
                                                AIREML.tol=1e-4,
                                                verbose=TRUE)

admixmap_b2_x1021_rv_cv_adj <- runassoc_flare("afr", nullmod_b2_x1021_rv_cv_adj,
                                                   cov_df = cov_b2_chr8_plof_ds,
                                                   chr = 8)
saveRDS(admixmap_b2_x1021_rv_cv_adj, "../Data/STAAR_exp_groups/x1021/admixmap_b2_x1021_rv_cv_adj.RDS")
######################## chr 10 #########################
known_loci_chr10 <- readRDS("../Data/STAAR_exp_groups/x100000007/individual_cond_pruned_var.RDS")
cov_b2_chr10 <- get_adj_cv(10, "b2", known_cv = known_loci_chr10,
                           cov_b2)

cov_b2_chr10 <-  cov_b2_chr10|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2),
         snp3 = if_else(is.na(snp3), mean(snp3, na.rm = TRUE), snp3),
         snp4 = if_else(is.na(snp4), mean(snp4, na.rm = TRUE), snp4),
         snp5 = if_else(is.na(snp5), mean(snp5, na.rm = TRUE), snp5),
         GFRSCYS = if_else(is.na(GFRSCYS), mean(GFRSCYS, na.rm = TRUE), GFRSCYS))


nullmod_b2_x100000007_no_genotype <- fitNullModel(x = cov_b2_chr10,
                                                  outcome="x100000007",
                                                  covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                           "PC4", "PC5",
                                                           "BKGRD1_C7", "GFRSCYS",
                                                           "CENTER"),
                                                  cov.mat=covMatList,
                                                  AIREML.tol=1e-4,
                                                  verbose=TRUE)

nullmod_b2_x100000007_cv_adj <- fitNullModel(x = cov_b2_chr10,
                                             outcome="x100000007",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS",
                                                      "CENTER", "snp1", "snp2", "snp3",
                                                      "snp4", "snp5"),
                                             cov.mat=covMatList,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)

admixmap_b2_x100000007_no_genotype <- runassoc_flare("amer", nullmod_b2_x100000007_no_genotype,
                                                     cov_df = cov_b2_chr10, chr = 10)
saveRDS(admixmap_b2_x100000007_no_genotype,
        "../Data/STAAR_exp_groups/x100000007/admixmap_b2_x100000007_no_genotype.RDS")

admixmap_b2_x100000007_cv_adj <- runassoc_flare("amer", nullmod_b2_x100000007_cv_adj,
                                                cov_df = cov_b2_chr10, chr = 10)
saveRDS(admixmap_b2_x100000007_cv_adj,
        "../Data/STAAR_exp_groups/x100000007/admixmap_b2_x100000007_cv_adj.RDS")

######################## chr 13 #########################
known_loci_chr13 <- readRDS("../Data/STAAR_exp_groups/x100004046/individual_cond_pruned_var.RDS")
cov_b2_chr13 <- get_adj_cv(13, "b2", known_cv = known_loci_chr13,
                           cov_b2)

cov_b2_chr13 <-  cov_b2_chr13|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2),
         snp3 = if_else(is.na(snp3), mean(snp3, na.rm = TRUE), snp3),
         GFRSCYS = if_else(is.na(GFRSCYS), mean(GFRSCYS, na.rm = TRUE), GFRSCYS))


nullmod_b2_x100004046_no_genotype <- fitNullModel(x = cov_b2_chr13,
                                                  outcome="x100004046",
                                                  covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                           "PC4", "PC5",
                                                           "BKGRD1_C7", "GFRSCYS",
                                                           "CENTER"),
                                                  cov.mat=covMatList,
                                                  AIREML.tol=1e-4,
                                                  verbose=TRUE)

nullmod_b2_x100004046_cv_adj <- fitNullModel(x = cov_b2_chr13,
                                              outcome="x100004046",
                                              covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                       "PC4", "PC5",
                                                       "BKGRD1_C7", "GFRSCYS",
                                                       "CENTER", "snp1", "snp2", "snp3"),
                                              cov.mat=covMatList,
                                              AIREML.tol=1e-4,
                                              verbose=TRUE)

admixmap_b2_x100004046_no_genotype <- runassoc_flare("afr", nullmod_b2_x100004046_no_genotype,
                                                     cov_df = cov_b2_chr13, chr = 13)
saveRDS(admixmap_b2_x100004046_no_genotype,
        "../Data/STAAR_exp_groups/x100004046/admixmap_b2_x100004046_no_genotype.RDS")

admixmap_b2_x100004046_cv_adj <- runassoc_flare("afr", nullmod_b2_x100004046_cv_adj,
                                                 cov_df = cov_b2_chr13, chr = 13)
saveRDS(admixmap_b2_x100004046_cv_adj,
        "../Data/STAAR_exp_groups/x100004046/admixmap_b2_x100004046_cv_adj.RDS")

######################## chr 16 - EXP #########################
known_loci_chr16_exp <- readRDS("../Data/STAAR_exp_groups/x1224/individual_cond_pruned_var.RDS")
cov_b2_chr16_exp <- get_adj_cv(16, "b2", known_cv = known_loci_chr16_exp,
                               cov_b2)

cov_b2_chr16_exp <-  cov_b2_chr16_exp|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2),
         snp3 = if_else(is.na(snp3), mean(snp3, na.rm = TRUE), snp3),
         snp4 = if_else(is.na(snp4), mean(snp4, na.rm = TRUE), snp4),
         snp5 = if_else(is.na(snp5), mean(snp5, na.rm = TRUE), snp5),
         snp6 = if_else(is.na(snp6), mean(snp6, na.rm = TRUE), snp6),
         snp7 = if_else(is.na(snp7), mean(snp7, na.rm = TRUE), snp7),
         snp8 = if_else(is.na(snp8), mean(snp8, na.rm = TRUE), snp8),
         snp9 = if_else(is.na(snp9), mean(snp9, na.rm = TRUE), snp9),
         snp10 = if_else(is.na(snp10), mean(snp10, na.rm = TRUE), snp10),
         snp11 = if_else(is.na(snp11), mean(snp11, na.rm = TRUE), snp11),
         GFRSCYS = if_else(is.na(GFRSCYS), mean(GFRSCYS, na.rm = TRUE), GFRSCYS))

nullmod_b2_x1224_no_genotype <- fitNullModel(x = cov_b2_chr16_exp,
                                             outcome="x1224",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS",
                                                      "CENTER"),
                                             cov.mat=covMatList,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)

nullmod_b2_x1224_cv_adj <- fitNullModel(x = cov_b2_chr16_exp,
                                        outcome="x1224",
                                        covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                 "PC4", "PC5",
                                                 "BKGRD1_C7", "GFRSCYS",
                                                 "CENTER", "snp1", "snp2", "snp3",
                                                 "snp4", "snp5", "snp6", "snp7",
                                                 "snp8", "snp9", "snp10", "snp11"),
                                        cov.mat=covMatList,
                                        AIREML.tol=1e-4,
                                        verbose=TRUE)

admixmap_b2_x1224_no_genotype <- runassoc_flare("amer", nullmod_b2_x1224_no_genotype,
                                                cov_df = cov_b2_chr16_exp, chr = 16)
saveRDS(admixmap_b2_x1224_no_genotype,
        "../Data/STAAR_exp_groups/x1224/admixmap_b2_x1224_no_genotype.RDS")

admixmap_b2_x1224_cv_adj <- runassoc_flare("amer", nullmod_b2_x1224_cv_adj,
                                           cov_df = cov_b2_chr16_exp, chr = 16)
saveRDS(admixmap_b2_x1224_cv_adj,
        "../Data/STAAR_exp_groups/x1224/admixmap_b2_x1224_cv_adj.RDS")


########### ADDING BURDEN SCORES FROM DIFF UNIQUE RV SETS ###########

# plof
chr16_plof_burden <- readRDS("../Data/STAAR_exp_groups/x1224/burden_plof.RDS")

intersected_samples <- intersect(rownames(cov_b1_chr16_exp), rownames(chr16_plof_burden))
cov_b1_chr16_exp_rv <- cbind(cov_b1_chr16_exp[intersected_samples, ],
                             chr16_plof_burden[intersected_samples, ])

for(i in 1:length(covMatList_rv)){
  covMatList_rv[[i]] <- covMatList_rv[[i]][rownames(cov_b1_chr16_exp_rv),
                                           rownames(cov_b1_chr16_exp_rv)]
}

cov_b1_chr16_exp_rv <- cov_b1_chr16_exp_rv|>
  rename(Burden_1_25 = `Burden(1, 25)`)

nullmod_b1_x1224_plof_adj <- fitNullModel(x = cov_b1_chr16_exp_rv,
                                          outcome="x1224",
                                          covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                   "PC4", "PC5",
                                                   "BKGRD1_C7", "GFRSCYS",
                                                   "CENTER", "Burden_1_25"),
                                          cov.mat=covMatList_rv,
                                          AIREML.tol=1e-4,
                                          verbose=TRUE)

admixmap_b1_x1224_plof_adj <- runassoc_flare("amer", nullmod_b1_x1224_plof_adj,
                                             cov_df = cov_b1_chr16_exp_rv, chr = 16)
saveRDS(admixmap_b1_x1224_plof_adj,
        "../Data/STAAR_exp_groups/x1224/admixmap_b1_x1224_plof_adj.RDS")

plot_admixmap_signal(admixmap_b1_x1224_no_genotype, admixmap_b1_x1224_plof_adj,
                     "Cys-gly, oxidized", 16, "Amerindian ancestry")


# plof_ds
chr16_plof_ds_burden <- readRDS("../Data/STAAR_exp_groups/x1224/burden_plof_ds.RDS")

intersected_samples <- intersect(rownames(cov_b1_chr16_exp), rownames(chr16_plof_ds_burden))
cov_b1_chr16_exp_rv_plof_ds <- cbind(cov_b1_chr16_exp[intersected_samples, ],
                                     chr16_plof_ds_burden[intersected_samples, ])

for(i in 1:length(covMatList_rv)){
  covMatList_rv[[i]] <- covMatList_rv[[i]][rownames(cov_b1_chr16_exp_rv_plof_ds),
                                           rownames(cov_b1_chr16_exp_rv_plof_ds)]
}

cov_b1_chr16_exp_rv_plof_ds <- cov_b1_chr16_exp_rv_plof_ds|>
  rename(Burden_1_25 = `Burden(1, 25)`) |>
  rename(linsight_1_25 = `Burden(1, 25)-LINSIGHT`)

nullmod_b1_x1224_plof_ds_adj <- fitNullModel(x = cov_b1_chr16_exp_rv_plof_ds,
                                             outcome="x1224",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS",
                                                      "CENTER", "Burden_1_25"),
                                             cov.mat=covMatList_rv,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)

admixmap_b1_x1224_plof_ds_adj <- runassoc_flare("amer", nullmod_b1_x1224_plof_ds_adj,
                                                cov_df = cov_b1_chr16_exp_rv, chr = 16)
saveRDS(admixmap_b1_x1224_plof_ds_adj,
        "../Data/STAAR_exp_groups/x1224/admixmap_b1_x1224_plof_ds_adj.RDS")

plot_admixmap_signal(admixmap_b1_x1224_no_genotype, admixmap_b1_x1224_plof_ds_adj,
                     "Cys-gly, oxidized", 16, "Amerindian ancestry")


# COMBINED EFFECTS
nullmod_b1_x1224_plof_ds_cv_adj <- fitNullModel(x = cov_b1_chr16_exp_rv_plof_ds,
                                                outcome="x1224",
                                                covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                         "PC4", "PC5",
                                                         "BKGRD1_C7", "GFRSCYS",
                                                         "CENTER", "Burden_1_25",
                                                         "snp1", "snp2", "snp3",
                                                         "snp4", "snp5"),
                                                cov.mat=covMatList_rv,
                                                AIREML.tol=1e-4,
                                                verbose=TRUE)

admixmap_b1_x1224_plof_ds_cv_adj <- runassoc_flare("amer", nullmod_b1_x1224_plof_ds_cv_adj,
                                                   cov_df = cov_b1_chr16_exp_rv, chr = 16)
saveRDS(admixmap_b1_x1224_plof_ds_adj,
        "../Data/STAAR_exp_groups/x1224/admixmap_b1_x1224_plof_ds_cv_adj.RDS")

plot_admixmap_signal(admixmap_b1_x1224_no_genotype, admixmap_b1_x1224_plof_ds_cv_adj,
                     "Cys-gly, oxidized", 16, "Amerindian ancestry")


#### missense mutation #####

chr16_missense_burden <- readRDS("../Data/STAAR_exp_groups/x1224/burden_missense.RDS")

intersected_samples <- intersect(rownames(cov_b1_chr16_exp), rownames(chr16_missense_burden))
cov_b1_chr16_exp_rv_missense <- cbind(cov_b1_chr16_exp[intersected_samples, ],
                                      chr16_missense_burden[intersected_samples, ])

cov_b1_chr16_exp_rv_missense <- cov_b1_chr16_exp_rv_missense|>
  rename(Burden_1_25 = `Burden(1, 25)`)

nullmod_b1_x1224_missense_adj <- fitNullModel(x = cov_b1_chr16_exp_rv_missense,
                                              outcome="x1224",
                                              covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                       "PC4", "PC5",
                                                       "BKGRD1_C7", "GFRSCYS",
                                                       "CENTER", "Burden_1_25"),
                                              cov.mat=covMatList_rv,
                                              AIREML.tol=1e-4,
                                              verbose=TRUE)

admixmap_b1_x1224_missense_adj <- runassoc_flare("amer", nullmod_b1_x1224_missense_adj,
                                                 cov_df = cov_b1_chr16_exp_rv_missense, chr = 16)

saveRDS(admixmap_b1_x1224_missense_adj,
        "../Data/STAAR_exp_groups/x1224/admixmap_b1_x1224_missense_adj.RDS")

plot_admixmap_signal(admixmap_b1_x1224_no_genotype, admixmap_b1_x1224_missense_adj,
                     "Cys-gly, oxidized", 16, "Amerindian ancestry")


## COMBINED EFFECT

nullmod_b1_x1224_missense_cv_adj <- fitNullModel(x = cov_b1_chr16_exp_rv_missense,
                                                 outcome="x1224",
                                                 covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                          "PC4", "PC5",
                                                          "BKGRD1_C7", "GFRSCYS",
                                                          "CENTER", "Burden_1_25",
                                                          "snp1", "snp2", "snp3",
                                                          "snp4", "snp5"),
                                                 cov.mat=covMatList_rv,
                                                 AIREML.tol=1e-4,
                                                 verbose=TRUE)


admixmap_b1_x1224_missense_cv_adj <- runassoc_flare("amer", nullmod_b1_x1224_missense_cv_adj,
                                                    cov_df = cov_b1_chr16_exp_rv_missense, chr = 16)
saveRDS(admixmap_b1_x1224_missense_cv_adj,
        "../Data/STAAR_exp_groups/x1224/admixmap_b1_x1224_missense_cv_adj.RDS")
