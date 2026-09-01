# load neccessary packages 
library(tidyverse) 
library(GENESIS)
library(gdsfmt)
library(GWASTools)
library(SNPRelate) 
library(SeqArray)
library(reshape2)
#library(STAARpipeline)
library(purrr)
library(tibble)

setwd("R:/HCHS_SOL/Projects/2024_rare_variants/Code")

cov_b1 <- readRDS("../Data/STAAR_model_cov_batch1.RDS")

gdsflare <- openfn.gds(paste0("../../../Projects/2023_local_ancestry_comparison_sol/Data/FLARE3/FLARE3_SNPs_filtered_chr", 
                              22, ".gds" ))
sample.id <- read.gdsn(index.gdsn(gdsflare, "sample.id"))
closefn.gds(gdsflare) # close the gdsfile

intersected_id_b1 <- intersect(paste0("SoL", sample.id), cov_b1$SUBJECT_ID)

rownames(cov_b1) <- cov_b1$SUBJECT_ID
cov_b1 <- cov_b1[intersected_id_b1,] 

cov_b1 <- cov_b1 |>
  mutate(GFRSCYS = if_else(is.na(GFRSCYS), mean(GFRSCYS, na.rm = TRUE), GFRSCYS))

# load covariate matrices
rownames(cov_b1) <- cov_b1$NWD_ID #3842 for b1
covMatList <- readRDS("../Data/STAAR_model_covmat_batch1.RDS")
for(i in 1:length(covMatList)){
  covMatList[[i]] <- covMatList[[i]][rownames(cov_b1), rownames(cov_b1)]
}

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


################ Negative control - x100001721 chr2 ################  
known_loci_chr2 <- readRDS("../Data/STAAR_neg_controls/x100001721/individual_cond_pruned_var_sensitivity.RDS")
cov_b1_chr2_cv_adj <- get_adj_cv(2, "b1", known_cv = known_loci_chr2,
                                 cov_b1)

cov_b1_chr2_cv_adj <-  cov_b1_chr2_cv_adj |>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2),
         snp3 = if_else(is.na(snp3), mean(snp3, na.rm = TRUE), snp3),
         snp4 = if_else(is.na(snp4), mean(snp4, na.rm = TRUE), snp4),
         snp5 = if_else(is.na(snp5), mean(snp5, na.rm = TRUE), snp5),
         snp6 = if_else(is.na(snp6), mean(snp6, na.rm = TRUE), snp6),
         GFRSCYS = if_else(is.na(GFRSCYS), mean(GFRSCYS, na.rm = TRUE), GFRSCYS))

nullmod_b1_x100001721_cv_adj <- fitNullModel(x = cov_b1_chr2_cv_adj,
                                             outcome="x100001721",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS", 
                                                      "CENTER", "snp1", "snp2",
                                                      "snp3", "snp4", "snp5",
                                                      "snp6"),
                                             cov.mat=covMatList,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)

admixmap_b1_x100001721_cv_adj <- runassoc_flare("afr", nullmod_b1_x100001721_cv_adj, 
                                                cov_df = cov_b1_chr2_cv_adj, 
                                                chr = 2)
saveRDS(admixmap_b1_x100001721_cv_adj, "../Data/STAAR_neg_controls/x100001721/admixmap_b1_x100001721_cv_sensitivity.RDS")



################ Negative control - x799 chr5 - no change ################  


################ Negative control - x1215 chr6 - no change ################  


################ Negative control - x278 chr16 ################  
known_loci_chr16 <- readRDS("../Data/STAAR_neg_controls/x278/individual_cond_pruned_var_sensitivity.RDS")
cov_b1_chr16 <- get_adj_cv(16, "b1", known_cv = known_loci_chr16,
                           cov_b1)

cov_b1_chr16 <-  cov_b1_chr16|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2),
         GFRSCYS = if_else(is.na(GFRSCYS), mean(GFRSCYS, na.rm = TRUE), GFRSCYS))

nullmod_b1_x278_cv_adj <- fitNullModel(x = cov_b1_chr16,
                                       outcome="x278",
                                       covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                "PC4", "PC5",
                                                "BKGRD1_C7", "GFRSCYS", 
                                                "CENTER", "snp1", "snp2"),
                                       cov.mat=covMatList,
                                       AIREML.tol=1e-4,
                                       verbose=TRUE)

admixmap_b1_x278_cv_adj <- runassoc_flare("amer", nullmod_b1_x278_cv_adj,
                                          cov_df = cov_b1_chr16, chr = 16) 
gc()

saveRDS(admixmap_b1_x278_cv_adj,
        "../Data/STAAR_neg_controls/x278/admixmap_b1_x278_cv_sensitivity.RDS")


################### EXPERIMENTAL GROUP -- chr8 x1021 ############
known_loci_chr8 <- readRDS("../Data/STAAR_exp_groups/x1021/individual_cond_pruned_var_sensitivity.RDS")
cov_b1_chr8 <- get_adj_cv(8, "b1", known_cv = known_loci_chr8,
                          cov_b1)

cov_b1_chr8 <-  cov_b1_chr8|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2),
         snp3 = if_else(is.na(snp3), mean(snp3, na.rm = TRUE), snp3),
         snp4 = if_else(is.na(snp4), mean(snp4, na.rm = TRUE), snp4),
         snp5 = if_else(is.na(snp5), mean(snp5, na.rm = TRUE), snp5),
         GFRSCYS = if_else(is.na(GFRSCYS), mean(GFRSCYS, na.rm = TRUE), GFRSCYS),
         x1021 = if_else(is.na(x1021), min(x1021, na.rm = TRUE), x1021))


nullmod_b1_x1021_cv_adj <- fitNullModel(x = cov_b1_chr8,
                                        outcome="x1021",
                                        covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                 "PC4", "PC5",
                                                 "BKGRD1_C7", "GFRSCYS", 
                                                 "CENTER", "snp1", "snp2", "snp3",
                                                 "snp4", "snp5"),
                                        cov.mat=covMatList,
                                        AIREML.tol=1e-4,
                                        verbose=TRUE)


admixmap_b1_x1021_cv_adj <- runassoc_flare("afr", nullmod_b1_x1021_cv_adj,
                                           cov_df = cov_b1_chr8, chr = 8) 
gc()
saveRDS(admixmap_b1_x1021_cv_adj,
        "../Data/STAAR_exp_groups/x1021/admixmap_b1_x1021_cv_sensitivity.RDS")



######################## chr 10 x100000007 #########################
known_loci_chr10 <- readRDS("../Data/STAAR_exp_groups/x100000007/individual_cond_pruned_var_sensitivity.RDS")
cov_b1_chr10 <- get_adj_cv(10, "b1", known_cv = known_loci_chr10,
                           cov_b1)

cov_b1_chr10 <-  cov_b1_chr10|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2),
         snp3 = if_else(is.na(snp3), mean(snp3, na.rm = TRUE), snp3),
         snp4 = if_else(is.na(snp4), mean(snp4, na.rm = TRUE), snp4),
         GFRSCYS = if_else(is.na(GFRSCYS), mean(GFRSCYS, na.rm = TRUE), GFRSCYS))

nullmod_b1_x100000007_cv_adj <- fitNullModel(x = cov_b1_chr10,
                                             outcome="x100000007",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS", 
                                                      "CENTER", "snp1", "snp2", "snp3",
                                                      "snp4"),
                                             cov.mat=covMatList,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)

admixmap_b1_x100000007_cv_adj <- runassoc_flare("amer", nullmod_b1_x100000007_cv_adj,
                                                cov_df = cov_b1_chr10, chr = 10) 
gc()

saveRDS(admixmap_b1_x100000007_cv_adj,
        "../Data/STAAR_exp_groups/x100000007/admixmap_b1_x100000007_cv_sensitivity.RDS")


######################## chr 13 x100004046 #########################
known_loci_chr13 <- readRDS("../Data/STAAR_exp_groups/x100004046/individual_cond_pruned_var_sensitivity.RDS")
cov_b1_chr13 <- get_adj_cv(13, "b1", known_cv = known_loci_chr13,
                           cov_b1)

cov_b1_chr13 <-  cov_b1_chr13|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2),
         GFRSCYS = if_else(is.na(GFRSCYS), mean(GFRSCYS, na.rm = TRUE), GFRSCYS))

nullmod_b1_x100004046_cv_adj <- fitNullModel(x = cov_b1_chr13,
                                             outcome="x100004046",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS", 
                                                      "CENTER", "snp1", "snp2"),
                                             cov.mat=covMatList,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)

admixmap_b1_x100004046_cv_adj <- runassoc_flare("afr", nullmod_b1_x100004046_cv_adj, 
                                                cov_df = cov_b1_chr13, chr = 13)
gc()
saveRDS(admixmap_b1_x100004046_cv_adj,
        "../Data/STAAR_exp_groups/x100004046/admixmap_b1_x100004046_cv_adj.RDS")


######################## chr 16 - EXP x1224 #########################
known_loci_chr16_exp <- readRDS("../Data/STAAR_exp_groups/x1224/individual_cond_pruned_var_sensitivity.RDS")
cov_b1_chr16_exp <- get_adj_cv(16, "b1", known_cv = known_loci_chr16_exp,
                               cov_b1)

cov_b1_chr16_exp <-  cov_b1_chr16_exp|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2),
         snp3 = if_else(is.na(snp3), mean(snp3, na.rm = TRUE), snp3),
         snp4 = if_else(is.na(snp4), mean(snp4, na.rm = TRUE), snp4),
         snp5 = if_else(is.na(snp5), mean(snp5, na.rm = TRUE), snp5),
         snp6 = if_else(is.na(snp6), mean(snp6, na.rm = TRUE), snp6),
         snp7 = if_else(is.na(snp7), mean(snp7, na.rm = TRUE), snp7))

colnames(cov_b1_chr16_exp)[2] <- "x1224"

nullmod_b1_x1224_cv_adj <- fitNullModel(x = cov_b1_chr16_exp,
                                        outcome="x1224",
                                        covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                 "PC4", "PC5",
                                                 "BKGRD1_C7", "GFRSCYS", 
                                                 "CENTER", "snp1", "snp2", "snp3",
                                                 "snp4", "snp5", "snp6", "snp7"),
                                        cov.mat=covMatList,
                                        AIREML.tol=1e-4,
                                        verbose=TRUE)



admixmap_b1_x1224_cv_adj <- runassoc_flare("amer", nullmod_b1_x1224_cv_adj,
                                           cov_df = cov_b1_chr16_exp, chr = 16) 

gc()
saveRDS(admixmap_b1_x1224_cv_adj,
        "../Data/STAAR_exp_groups/x1224/admixmap_b1_x1224_cv_sensitivity.RDS")



# define a function for getting pruned common variants selected from the STAAR
# pipeline and use it as a covariate to be adjusted in the admixutre mapping 
get_adj_cv <- function(chr, batch, known_cv, covariate_df, i){
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

# load gds files catalog that need to be accessed 
load("../Data/STAAR_prep/agds_dir_scale_up_final.Rdata") #variable name agds_dir


cov_b1 <- readRDS("../Data/STAAR_model_cov_batch1_scaled_up.RDS")
gdsflare <- openfn.gds(paste0("../../../Projects/2023_local_ancestry_comparison_sol/Data/FLARE3/FLARE3_SNPs_filtered_chr", 
                              22, ".gds" ))
sample.id <- read.gdsn(index.gdsn(gdsflare, "sample.id"))
closefn.gds(gdsflare) # close the gdsfile

intersected_id_b1 <- intersect(paste0("SoL", sample.id), cov_b1$SUBJECT_ID)

rownames(cov_b1) <- cov_b1$SUBJECT_ID
cov_b1 <- cov_b1[intersected_id_b1,] 

# load covariate matrices
rownames(cov_b1) <- cov_b1$NWD_ID #3842 for b1
covMatList <- readRDS("../Data/STAAR_model_covmat_batch1.RDS")
for(i in 1:length(covMatList)){
  covMatList[[i]] <- covMatList[[i]][rownames(cov_b1), rownames(cov_b1)]
}

###### agds_dir[1], chr11 arachido,  x100009332 #####

known_loci_chr11_arachido <- readRDS("../Data/STAAR_exp_groups/x100009332/individual_cond_pruned_var_sensitivity.RDS")
cov_b1_chr11_arachido <- get_adj_cv(11, "b1", known_cv = known_loci_chr11_arachido,
                                    cov_b1, i = 1)

cov_b1_chr11_arachido <-  cov_b1_chr11_arachido|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2),
         snp3 = if_else(is.na(snp3), mean(snp3, na.rm = TRUE), snp3),
         snp4 = if_else(is.na(snp4), mean(snp4, na.rm = TRUE), snp4),
         snp5 = if_else(is.na(snp5), mean(snp5, na.rm = TRUE), snp5),
         snp6 = if_else(is.na(snp6), mean(snp6, na.rm = TRUE), snp6),
         snp7 = if_else(is.na(snp7), mean(snp7, na.rm = TRUE), snp7))


nullmod_b1_x100009332_cv_adj <- fitNullModel(x = cov_b1_chr11_arachido,
                                             outcome="x100009332",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS", 
                                                      "CENTER", "snp1", "snp2",
                                                      "snp3", "snp4", "snp5",
                                                      "snp6", "snp7"),
                                             cov.mat=covMatList,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)


admixmap_b1_x100009332_cv_adj <- runassoc_flare("amer", nullmod_b1_x100009332_cv_adj, 
                                                cov_df = cov_b1_chr11_arachido, chr = 11)
saveRDS(admixmap_b1_x100009332_cv_adj, "../Data/STAAR_exp_groups/x100009332/admixmap_b1_x100009332_cv_sensitivity.RDS")

gc()


###### agds_dir[2], chr11 3beta,  x100006370 #####

known_loci_chr11_3beta <- readRDS("../Data/STAAR_exp_groups/x100006370/individual_cond_pruned_var_sensitivity.RDS")
cov_b1_chr11_3beta <- get_adj_cv(11, "b1", known_cv = known_loci_chr11_3beta,
                                 cov_b1, i = 2)

nullmod_b1_x100006370_cv_adj <- fitNullModel(x = cov_b1_chr11_3beta,
                                             outcome="x100006370",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS", 
                                                      "CENTER", "snp1", "snp2",
                                                      "snp3", "snp4"),
                                             cov.mat=covMatList,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)


admixmap_b1_x100006370_cv_adj <- runassoc_flare("amer", nullmod_b1_x100006370_cv_adj, 
                                                cov_df = cov_b1_chr11_3beta, chr = 11)
saveRDS(admixmap_b1_x100006370_cv_adj, "../Data/STAAR_exp_groups/x100006370/admixmap_b1_x100006370_cv_sensitivity.RDS")
gc()


###### agds_dir[3], chr12 1_methyl, x100001208 #####

known_loci_chr12_1_methyl <- readRDS("../Data/STAAR_exp_groups/x100001208/individual_cond_pruned_var_sensitivity.RDS")
cov_b1_chr12_1_methyl <- get_adj_cv(12, "b1", known_cv = known_loci_chr12_1_methyl,
                                    cov_b1, i = 3)
nullmod_b1_x100001208_cv_adj <- fitNullModel(x = cov_b1_chr12_1_methyl,
                                             outcome="x100001208",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS", 
                                                      "CENTER", "snp1", "snp2",
                                                      "snp3"),
                                             cov.mat=covMatList,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)

admixmap_b1_x100001208_cv_adj <- runassoc_flare("amer", nullmod_b1_x100001208_cv_adj, 
                                                cov_df = cov_b1_chr12_1_methyl, chr = 12)
saveRDS(admixmap_b1_x100001208_cv_adj, "../Data/STAAR_exp_groups/x100001208/admixmap_b1_x100001208_cv_sensitivity.RDS")
gc()


###### agds_dir[4], chr12 1_ethylmalonate, x2054 -- no change #####


###### agds_dir[5], chr16  #####

known_loci_chr16_propyl <- readRDS("../Data/STAAR_exp_groups/x100006264/individual_cond_pruned_var_sensitivity.RDS")
cov_b1_chr16_propyl <- get_adj_cv(16, "b1", known_cv = known_loci_chr16_propyl,
                                  cov_b1, i = 5)


cov_b1_chr16_propyl <- cov_b1_chr16_propyl |>
  mutate(snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2),
         snp5 = if_else(is.na(snp5), mean(snp5, na.rm = TRUE), snp5),
         snp6 = if_else(is.na(snp6), mean(snp6, na.rm = TRUE), snp6))


nullmod_b1_x100006264_cv_adj <- fitNullModel(x = cov_b1_chr16_propyl,
                                             outcome="x100006264",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS", 
                                                      "CENTER", "snp1", "snp2",
                                                      "snp3", "snp4", "snp5", 
                                                      "snp6", "snp8",
                                                      "snp9", "snp10", "snp11",
                                                      "snp12", "snp13", "snp14", 
                                                      "snp15"),
                                             cov.mat=covMatList,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)

admixmap_b1_x100006264_cv_adj <- runassoc_flare("afr", nullmod_b1_x100006264_cv_adj, 
                                                cov_df = cov_b1_chr16_propyl, chr = 16)
saveRDS(admixmap_b1_x100006264_cv_adj, "../Data/STAAR_exp_groups/x100006264/admixmap_b1_x100006264_cv_sensitivity.RDS")
gc()


###### agds_dir[6], chr2  #####

known_loci_chr2 <- readRDS("../Data/STAAR_exp_groups/x100001266/individual_cond_pruned_var_sensitivity.RDS")
cov_b1_known_loci_chr2 <- get_adj_cv(2, "b1", known_cv = known_loci_chr2,
                                     cov_b1, i = 6)

# mean imputation
cov_b1_known_loci_chr2 <- cov_b1_known_loci_chr2 |>
  mutate(snp4 = if_else(is.na(snp4), mean(snp4, na.rm = TRUE), snp4))
sum(is.na(cov_b1_known_loci_chr2))

nullmod_b1_x100001266_cv_adj <- fitNullModel(x = cov_b1_known_loci_chr2,
                                             outcome="x100001266",
                                             covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                      "PC4", "PC5",
                                                      "BKGRD1_C7", "GFRSCYS", 
                                                      "CENTER", "snp1", "snp2", 
                                                      "snp3", "snp4",
                                                      "snp5", "snp6"),
                                             cov.mat=covMatList,
                                             AIREML.tol=1e-4,
                                             verbose=TRUE)  

admixmap_b1_x100001266_cv_adj <- runassoc_flare("afr", nullmod_b1_x100001266_cv_adj, 
                                                cov_df = cov_b1_known_loci_chr2, chr = 2)
saveRDS(admixmap_b1_x100001266_cv_adj, "../Data/STAAR_exp_groups/x100001266/admixmap_b1_x100001266_cv_sensitivity.RDS")
gc()


###### agds_dir[7], chr5  #####

known_loci_chr5 <- readRDS("../Data/STAAR_exp_groups/x1114/individual_cond_pruned_var_sensitivity.RDS")
cov_b1_known_loci_chr5 <- get_adj_cv(5, "b1", known_cv = known_loci_chr5,
                                     cov_b1, i = 7)

cov_b1_known_loci_chr5 <-  cov_b1_known_loci_chr5|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2),
         snp3 = if_else(is.na(snp3), mean(snp3, na.rm = TRUE), snp3),
         snp4 = if_else(is.na(snp4), mean(snp4, na.rm = TRUE), snp4),
         snp5 = if_else(is.na(snp5), mean(snp5, na.rm = TRUE), snp5),
         snp6 = if_else(is.na(snp6), mean(snp6, na.rm = TRUE), snp6))



nullmod_b1_x1114_cv_adj <- fitNullModel(x = cov_b1_known_loci_chr5,
                                        outcome="x1114",
                                        covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                 "PC4", "PC5",
                                                 "BKGRD1_C7", "GFRSCYS", 
                                                 "CENTER", "snp1", "snp2", 
                                                 "snp3", "snp4",
                                                 "snp5", "snp6"),
                                        cov.mat=covMatList,
                                        AIREML.tol=1e-4,
                                        verbose=TRUE)



admixmap_b1_x1114_cv_adj <- runassoc_flare("amer", nullmod_b1_x1114_cv_adj, 
                                           cov_df = cov_b1_known_loci_chr5, chr = 5)
saveRDS(admixmap_b1_x1114_cv_adj, "../Data/STAAR_exp_groups/x1114/admixmap_b1_x1114_cv_sensitivity.RDS")
gc()

###### agds_dir[8], chr8  #####

known_loci_chr8 <- readRDS("../Data/STAAR_exp_groups/x192/individual_cond_pruned_var_sensitivity.RDS")
cov_b1_known_loci_chr8 <- get_adj_cv(8, "b1", known_cv = known_loci_chr8,
                                     cov_b1, i = 8)


cov_b1_known_loci_chr8 <-  cov_b1_known_loci_chr8|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2),
         snp3 = if_else(is.na(snp3), mean(snp3, na.rm = TRUE), snp3))

nullmod_b1_x192_cv_adj <- fitNullModel(x = cov_b1_known_loci_chr8,
                                       outcome="x192",
                                       covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                "PC4", "PC5",
                                                "BKGRD1_C7", "GFRSCYS", 
                                                "CENTER", "snp1", "snp2", 
                                                "snp3"),
                                       cov.mat=covMatList,
                                       AIREML.tol=1e-4,
                                       verbose=TRUE)

admixmap_b1_x192_cv_adj <- runassoc_flare("amer", nullmod_b1_x192_cv_adj, 
                                          cov_df = cov_b1_known_loci_chr8, chr = 8)
saveRDS(admixmap_b1_x192_cv_adj, "../Data/STAAR_exp_groups/x192/admixmap_b1_x192_cv_sensitivity.RDS")
gc()


########## load all sensitivity analysis results Table 6 ##########

admixmap_b1_x100001721_cv_sensitity <- readRDS("../Data/STAAR_neg_controls/x100001721/admixmap_b1_x100001721_cv_sensitivity.RDS")
admixmap_b1_x100001721_cv_sensitity[which.min(admixmap_b1_x100001721_cv_sensitity$Score.pval),]

admixmap_b1_x278 <- readRDS("../Data/STAAR_neg_controls/x278/admixmap_b1_x278_cv_sensitivity.RDS")
admixmap_b1_x278[which.min(admixmap_b1_x278$Score.pval),]

admixmap_b1_x1266 <- readRDS("../Data/STAAR_exp_groups/x100001266/admixmap_b1_x100001266_cv_sensitivity.RDS")
admixmap_b1_x1266[which.min(admixmap_b1_x1266$Score.pval),]

admixmap_b1_x1114 <- readRDS("../Data/STAAR_exp_groups/x1114/admixmap_b1_x1114_cv_sensitivity.RDS")
admixmap_b1_x1114[which.min(admixmap_b1_x1114$Score.pval),]

admixmap_b1_x192 <- readRDS("../Data/STAAR_exp_groups/x192/admixmap_b1_x192_cv_sensitivity.RDS")
admixmap_b1_x192[which.min(admixmap_b1_x192$Score.pval),]

admixmap_b1_x1021 <- readRDS("../Data/STAAR_exp_groups/x1021/admixmap_b1_x1021_cv_sensitivity.RDS")
admixmap_b1_x1021[which.min(admixmap_b1_x1021$Score.pval),]

admixmap_b1_x100000007 <- readRDS("../Data/STAAR_exp_groups/x100000007/admixmap_b1_x100000007_cv_sensitivity.RDS")
admixmap_b1_x100000007[which.min(admixmap_b1_x100000007$Score.pval),]

admixmap_b1_x100009332 <- readRDS("../Data/STAAR_exp_groups/x100009332/admixmap_b1_x100009332_cv_sensitivity.RDS")
admixmap_b1_x100009332[which.min(admixmap_b1_x100009332$Score.pval),]

admixmap_b1_x100006370 <- readRDS("../Data/STAAR_exp_groups/x100006370/admixmap_b1_x100006370_cv_sensitivity.RDS")
admixmap_b1_x100006370[which.min(admixmap_b1_x100006370$Score.pval),]

admixmap_b1_x100001208 <- readRDS("../Data/STAAR_exp_groups/x100001208/admixmap_b1_x100001208_cv_sensitivity.RDS")
admixmap_b1_x100001208[which.min(admixmap_b1_x100001208$Score.pval),]

admixmap_b1_x100004046 <- readRDS("../Data/STAAR_exp_groups/x100004046/admixmap_b1_x100004046_cv_sensitivity.RDS")
admixmap_b1_x100004046[which.min(admixmap_b1_x100004046$Score.pval),]

admixmap_b1_x100006264 <- readRDS("../Data/STAAR_exp_groups/x100006264/admixmap_b1_x100006264_cv_sensitivity.RDS")
admixmap_b1_x100006264[which.min(admixmap_b1_x100006264$Score.pval),]

admixmap_b1_x1224 <- readRDS("../Data/STAAR_exp_groups/x1224/admixmap_b1_x1224_cv_sensitivity.RDS")
admixmap_b1_x1224[which.min(admixmap_b1_x1224$Score.pval),]



