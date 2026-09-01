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
library(future.apply)
library(progressr)
library(pheatmap)
handlers("txtprogressbar")

setwd("R:/HCHS_SOL/Projects/2024_rare_variants/Code")
# subsetting overlapped individuals from DRAGEN gds and the local ancestry files
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


# Extract local ancestry
# ancestry takes abbreviations such as "afr", "amer", or "eur"
get_local_ancestry <- function(ancestry, cov_df, chr){
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
  # rownames(cov_df) <- cov_df$NWD_ID
  # cov_df <- cov_df[rownames(nullmod$fit), ]
  rownames(cov_df) <- cov_df$SUBJECT_ID
  # Filter genotype matrix to include only filtered IDs and matched phenotype samples
  geno <- geno[, rownames(cov_df)]
  print(dim(geno))
  
  geno <- t(geno)
  #print(sum(colnames(geno) != matched_solid)) # should be 0
  gc()
  return(list(geno, snppos))
}


get_burden_la <- function(cov_df, burden_df, covMatList, ancestry, chr){
  
  shared   <- intersect(rownames(cov_df), rownames(burden_df))
  combined <- cbind(cov_df[shared, , drop = FALSE],
                    burden_df[shared, , drop = FALSE])
  rownames(combined) <- combined[["SUBJECT_ID"]]
  
  la_out <- get_local_ancestry(ancestry, cov_df = combined, chr = chr)
  la     <- la_out[[1]]
  snppos <- la_out[[2]]
  
  combined <- combined[rownames(la), , drop = FALSE]
  ids <- combined$NWD_ID
  rownames(combined) <- combined$NWD_ID
  
  covMatList_rv <- covMatList
  for (i in seq_along(covMatList_rv)) {
    covMatList_rv[[i]] <- covMatList_rv[[i]][ids, ids]
  }
  
  list(la_obj = list(la, snppos),
       dat    = combined,
       covMat = covMatList_rv)
}


# ---- top level: no enclosing environment to capture ----
fit_locus <- function(locus, dat, geno_mat, snppos, covMat,
                      base_covars, burden_names, p, top_burden) {
  p(sprintf("locus %s", locus))
  dat$la_loci <- geno_mat[, locus]
  nm_nogeno <- fitNullModel(dat, outcome = "la_loci",
                            covars = base_covars, cov.mat = covMat, verbose = FALSE)
  if(burden_names == top_burden){
    print("only one burden scores")
    nm_geno   <- fitNullModel(dat, outcome = "la_loci",
                              covars = c(base_covars, burden_names),
                              cov.mat = covMat, verbose = FALSE)
    ve_allburden <- 1 - sum(nm_geno$varComp)  / sum(nm_nogeno$varComp)
    ve_topburden <- ve_allburden
  }else{
    nm_topburden <-  fitNullModel(dat, outcome = "la_loci",
                                  covars = c(base_covars,top_burden), cov.mat = covMat, verbose = FALSE)
    nm_geno   <- fitNullModel(dat, outcome = "la_loci",
                              covars = c(base_covars, burden_names),
                              cov.mat = covMat, verbose = FALSE)
    ve_allburden <- 1 - sum(nm_geno$varComp)  / sum(nm_nogeno$varComp)
    ve_topburden  <- 1 - sum(nm_topburden$varComp) / sum(nm_nogeno$varComp)
  }
  cbind(snppos[locus, ],
        var_explain_allburden = ve_allburden,
        var_explain_topburden  = ve_topburden)
}


fit_var_explained <- function(la_obj, am_result, burden_names, dat, covMat, top_burden,
                              base_covars = c("AGE","GENDER","PC1","PC2","PC3","PC4","PC5",
                                              "BKGRD1_C7","GFRSCYS","CENTER"),
                              pval_thresh = 2.58e-6, n_sample = 50, seed = 1) {
  geno_mat <- la_obj[[1]]
  snppos   <- la_obj[[2]]
  
  pos_sig  <- geno_mat[, which(am_result$Score.pval < pval_thresh), drop = FALSE]
  uniq_idx <- !duplicated(pos_sig, MARGIN = 2)
  loci     <- colnames(pos_sig)[uniq_idx]        # significant + unique loci only
  
  set.seed(seed)
  loci <- sample(loci, min(n_sample, length(loci)))
  
  geno_sub <- geno_mat[, loci, drop = FALSE]     # by name, so still correct
  snp_sub  <- snppos[loci, , drop = FALSE]
  
  with_progress({
    p <- progressor(along = loci)
    results <- future_lapply(
      loci, fit_locus, top_burden = top_burden,   
      dat = dat, geno_mat = geno_sub, snppos = snp_sub, covMat = covMat,
      base_covars = base_covars, burden_names = burden_names, p = p,
      future.seed     = TRUE,
      future.packages = "GENESIS"
    )
  })
  
  as.data.frame(do.call(rbind, results))
}

# x100001721
chr2_burden_1721 <- readRDS("../Data/STAAR_neg_controls/x100001721/burdens_combined.RDS")
am_result_x100001721 <- readRDS("../Data/STAAR_neg_controls/x100001721/admixmap_b1_x100001721_no_genotype.RDS")

out <- get_burden_la(cov_b1, chr2_burden_1721, covMatList, ancestry = "afr", chr = 2)

plan(multisession, workers = max(1, availableCores() - 8))
handlers(global = TRUE)

burden_names <- c("AUP1_disruptive_missense",
                  "CYP26B1_missense", "ALMS1_missense",
                  "NAT8_missense", "STAMBP_synonymous",
                  "ALMS1_synonymous", "EGR4_synonymous")

top_burden <- "ALMS1_missense"

res <- fit_var_explained(la_obj = out$la_obj, am_result = am_result_x100001721,
                         burden_names = colnames(chr2_burden_1721),
                         dat = out$dat, covMat = out$covMat,
                         top_burden = top_burden)

saveRDS(res, "../Data/Variance_explained/b1_chr2_x100001721_burden_explained.RDS")

# x100001266
chr2_burden_1266 <- readRDS("../Data/STAAR_exp_groups/x100001266/burdens_combined.RDS")
am_result_x100001266 <- readRDS("../Data/STAAR_exp_groups/x100001266/admixmap_b1_x100001266_no_genotype.RDS")

out <- get_burden_la(cov_b1, chr2_burden_1266, covMatList, ancestry = "afr", chr = 2)

plan(multisession, workers = max(1, availableCores() - 8))
handlers(global = TRUE)

burden_names <- colnames(chr2_burden_1266)[1]

top_burden <-  "DUSP11_enhancer_CAGE"

res <- fit_var_explained(la_obj = out$la_obj, am_result = am_result_x100001266,
                         burden_names = burden_names,
                         dat = out$dat, covMat = out$covMat,
                         top_burden = top_burden)

saveRDS(res, "../Data/Variance_explained/b1_chr2_x100001266_burden_explained.RDS")



# x100001114
chr5_burden_1114 <- readRDS("../Data/STAAR_exp_groups/x1114/burden_plof_ds.RDS")
am_result_x1114 <- readRDS("../Data/STAAR_exp_groups/x1114/admixmap_b1_x1114_no_genotype.RDS")

out <- get_burden_la(cov_b1, chr5_burden_1114, covMatList, ancestry = "amer", chr = 5)


burden_names <- colnames(chr5_burden_1114)[1]
top_burden <-  colnames(chr5_burden_1114)[1]

plan(multisession, workers = max(1, availableCores() - 8))
handlers(global = TRUE)

res <- fit_var_explained(la_obj = out$la_obj, am_result = am_result_x1114,
                         burden_names = burden_names,
                         dat = out$dat, covMat = out$covMat,
                         top_burden = top_burden)

saveRDS(res, "../Data/Variance_explained/b1_chr5_x1114_burden_explained.RDS")


# x1021
chr8_burden_x1021 <- readRDS("../Data/STAAR_exp_groups/x1021/burden_plof_ds.RDS")
am_result_x1021 <- readRDS("../Data/STAAR_exp_groups/x1021/admixmap_b1_x1021_no_genotype.RDS")


burden_names <- colnames(chr8_burden_x1021)
top_burden <- colnames(chr8_burden_x1021)

out <- get_burden_la(cov_b1, chr8_burden_x1021, covMatList, ancestry = "afr", chr = 8)


plan(multisession, workers = max(1, availableCores() - 8))
handlers(global = TRUE)

res <- fit_var_explained(la_obj = out$la_obj, am_result = am_result_x1021,
                         burden_names = burden_names,
                         dat = out$dat, covMat = out$covMat,
                         top_burden = top_burden)

saveRDS(res, "../Data/Variance_explained/b1_chr8_x1021_burden_explained.RDS")



# x100009332
chr11_burden_x9332 <- readRDS("../Data/STAAR_exp_groups/x100009332/burden_synonymous.RDS")
am_result_x9332 <- readRDS("../Data/STAAR_exp_groups/x100009332/admixmap_b1_x100009332_no_genotype.RDS")
burden_names <- colnames(chr11_burden_x9332)
top_burden <- colnames(chr11_burden_x9332)

out <- get_burden_la(cov_b1, chr11_burden_x9332, covMatList, ancestry = "amer", chr = 11)


plan(multisession, workers = max(1, availableCores() - 20))
handlers(global = TRUE)

res <- fit_var_explained(la_obj = out$la_obj, am_result = am_result_x9332,
                         burden_names = burden_names,
                         dat = out$dat, covMat = out$covMat,
                         top_burden = top_burden)

saveRDS(res, "../Data/Variance_explained/b1_chr11_x9332_burden_explained.RDS")


# x100001208
chr12_burden_1208 <- readRDS("../Data/STAAR_exp_groups/x100001208/burden_plof_ds.RDS")
am_result_x1208 <- readRDS("../Data/STAAR_exp_groups/x100001208/admixmap_b1_x100001208_no_genotype.RDS")
burden_names <- colnames(chr12_burden_1208)
top_burden <- colnames(chr12_burden_1208)

out <- get_burden_la(cov_b1, chr12_burden_1208, covMatList, ancestry = "amer", chr = 12)


plan(multisession, workers = max(1, availableCores() - 10))
handlers(global = TRUE)

res <- fit_var_explained(la_obj = out$la_obj, am_result = am_result_x1208,
                         burden_names = burden_names,
                         dat = out$dat, covMat = out$covMat,
                         top_burden = top_burden)

saveRDS(res, "../Data/Variance_explained/b1_chr12_x1208_burden_explained.RDS")


# x2054
chr12_burden_2054 <- readRDS("../Data/STAAR_exp_groups/x2054/burden_plof_ds.RDS")
am_result_x2054 <- readRDS("../Data/STAAR_exp_groups/x2054/admixmap_b1_x2054_no_genotype.RDS")

burden_names <- colnames(chr12_burden_2054)
top_burden <- colnames(chr12_burden_2054)

out <- get_burden_la(cov_b1, burden_df = chr12_burden_2054,
                     covMatList, ancestry = "afr", chr = 12)

plan(multisession, workers = max(1, availableCores() - 7))
handlers(global = TRUE)

res <- fit_var_explained(la_obj = out$la_obj, am_result = am_result_x2054,
                         burden_names = colnames(chr12_burden_2054),
                         dat = out$dat, covMat = out$covMat,
                         top_burden = top_burden)

saveRDS(res, "../Data/Variance_explained/b1_chr12_x2054_burden_explained.RDS")



# x100006264
chr16_burden_x100006264 <- readRDS("../Data/STAAR_exp_groups/x100006264/burdens_combined.RDS")
am_result_x100006264 <- readRDS("../Data/STAAR_exp_groups/x100006264/admixmap_b1_x100006264_no_genotype.RDS")

out <- get_burden_la(cov_b1, chr16_burden_x100006264, covMatList, ancestry = "afr", chr = 16)

plan(multisession, workers = max(1, availableCores() - 5))
handlers(global = TRUE)

burden_names <- c("IL4R_UTR_cond",
                  "IL4R_enhancer_CAGE_cond",
                  "IL4R_promoter_DHS_cond",
                  "IL4R_promoter_CAGE_cond",
                  "IL4R_enhancer_DHS_cond",
                  "ITGAL_enhancer_DHS_cond",
                  "KDM8_enhancer_DHS_cond",
                  "PRRT2_enhancer_DHS_cond",
                  "FAM57B_enhancer_DHS_cond",
                  "FUS_enhancer_DHS_cond",
                  "MVP_UTR_cond",
                  "NUPR1_promoter_DHS_cond",
                  "RNF40_UTR_cond")

top_burden <- "IL4R_promoter_DHS_cond"

res <- fit_var_explained(la_obj = out$la_obj, am_result = am_result_x100006264,
                         burden_names = colnames(chr16_burden_x100006264),
                         dat = out$dat, covMat = out$covMat,
                         top_burden = top_burden)

# x1224
chr16_burden_x1224 <- readRDS("../Data/STAAR_exp_groups/x1224/burden_plof_ds.RDS")
am_result_x1224 <- readRDS("../Data/STAAR_exp_groups/x1224/admixmap_b1_x1224_no_genotype.RDS")

out <- get_burden_la(cov_b1, chr16_burden_x1224, covMatList, ancestry = "amer", chr = 16)
burden_names <- colnames(chr16_burden_x1224)
top_burden <- colnames(chr16_burden_x1224)

plan(multisession, workers = max(1, availableCores() - 25))
handlers(global = TRUE)

res <- fit_var_explained(la_obj = out$la_obj, am_result = am_result_x1224,
                         burden_names = burden_names,
                         dat = out$dat, covMat = out$covMat,
                         top_burden = top_burden)

saveRDS(res, "../Data/Variance_explained/b1_chr16_x1224_burden_explained.RDS")
