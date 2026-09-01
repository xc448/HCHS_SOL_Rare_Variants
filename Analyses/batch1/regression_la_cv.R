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

cov_b1 <- readRDS("../Data/STAAR_model_cov_batch1.RDS")
#cov_b2 <- readRDS("../Data/STAAR_model_cov_batch2.RDS")

gdsflare <- openfn.gds(paste0("../../../Projects/2023_local_ancestry_comparison_sol/Data/FLARE3/FLARE3_SNPs_filtered_chr", 
                              22, ".gds" ))
sample.id <- read.gdsn(index.gdsn(gdsflare, "sample.id"))
closefn.gds(gdsflare) # close the gdsfile
rm(gdsflare); gc()

intersected_id_b1 <- intersect(paste0("SoL", sample.id), cov_b1$SUBJECT_ID)

rownames(cov_b1) <- cov_b1$SUBJECT_ID
cov_b1 <- cov_b1[intersected_id_b1,] 

# load covariate matrices
rownames(cov_b1) <- cov_b1$NWD_ID #3842 for b1
covMatList <- readRDS("../Data/STAAR_model_covmat_batch1.RDS")
for(i in 1:length(covMatList)){
  covMatList[[i]] <- covMatList[[i]][rownames(cov_b1), rownames(cov_b1)]
}


get_adj_cv <- function(chr, batch, known_cv, covariate_df){
  genofile <- seqOpen(paste0("../Data/chr", chr, "_DRAGEN_agds_",
                             batch, "_final.gds"))
  position   <- seqGetData(genofile, "position")
  variant_id <- seqGetData(genofile, "variant.id")

  # keep variants whose position is in known_cv
  keep       <- position %in% known_cv$POS
  variant_id <- variant_id[keep]

  rownames(covariate_df) <- covariate_df$NWD_ID
  seqSetFilter(genofile, variant.id = variant_id,
               sample.id = rownames(covariate_df))

  genotype <- seqGetData(genofile, "genotype")
  sampleid <- seqGetData(genofile, "sample.id")
  dimnames(genotype)[[2]] <- sampleid

  # positions of the variants actually returned (GDS order)
  filt_pos <- seqGetData(genofile, "position")

  # reorder the variant (3rd) dimension to match known_cv$POS
  ord <- match(known_cv$POS, filt_pos)   # NA where a locus is missing
  if (any(is.na(ord))) {
    warning(sum(is.na(ord)), " known loci not found in chr", chr,
            " ", batch, "; dropping them.")
    ord <- ord[!is.na(ord)]
  }
  genotype <- genotype[, , ord, drop = FALSE]

  num_snps <- dim(genotype)[3]
  # allele dosage for the second haplotype, now in known_cv order
  snp_df <- map_dfc(1:num_snps, function(i) {
    snp_vals <- genotype[2, rownames(covariate_df), i]
    tibble(!!paste0("snp", i) := snp_vals)
  })

  covariate_df <- bind_cols(covariate_df, snp_df)

  seqClose(genofile)
  return(covariate_df)
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



# ---- top level: no enclosing environment to capture ----
fit_locus <- function(locus, dat, geno_mat, snppos, covMat,
                      base_covars, snp_covars, p, top_cv) {
  p(sprintf("locus %s", locus))
  dat$la_loci <- geno_mat[, locus]
  nm_nogeno <- fitNullModel(dat, outcome = "la_loci",
                            covars = base_covars, cov.mat = covMat, verbose = FALSE)
  nm_topcv <-  fitNullModel(dat, outcome = "la_loci",
                            covars = c(base_covars,top_cv), cov.mat = covMat, verbose = FALSE)
  nm_geno   <- fitNullModel(dat, outcome = "la_loci",
                            covars = c(base_covars, snp_covars),
                            cov.mat = covMat, verbose = FALSE)
  ve_allvar <- 1 - sum(nm_geno$varComp)  / sum(nm_nogeno$varComp)
  ve_topcv  <- 1 - sum(nm_topcv$varComp) / sum(nm_nogeno$varComp)
  
  cbind(snppos[locus, ],
        var_explain_allvar = ve_allvar,
        var_explain_topcv  = ve_topcv)
}

fit_var_explained <- function(la_obj, am_result, n_snp, dat, covMat, top_cv,
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
  
  snp_covars <- paste0("snp", seq_len(n_snp))
  
  with_progress({
    p <- progressor(along = loci)
    results <- future_lapply(
      loci, fit_locus, top_cv = top_cv,   
      dat = dat, geno_mat = geno_sub, snppos = snp_sub, covMat = covMat,
      base_covars = base_covars, snp_covars = snp_covars, p = p,
      future.seed     = TRUE,
      future.packages = "GENESIS"
    )
  })
  
  as.data.frame(do.call(rbind, results))
}

# chromosome 2
known_loci_chr2_1721 <- readRDS("../Data/STAAR_neg_controls/x100001721/individual_cond_pruned_var.RDS")
cov_b1_chr2_cv_adj <- get_adj_cv(2, "b1", known_cv = known_loci_chr2_1721,
                                 cov_b1)

cov_b1_chr2_cv_adj <-  cov_b1_chr2_cv_adj |>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2),
         snp3 = if_else(is.na(snp3), mean(snp3, na.rm = TRUE), snp3),
         snp4 = if_else(is.na(snp4), mean(snp4, na.rm = TRUE), snp4),
         snp5 = if_else(is.na(snp5), mean(snp5, na.rm = TRUE), snp5),
         snp6 = if_else(is.na(snp6), mean(snp6, na.rm = TRUE), snp6),
         snp7 = if_else(is.na(snp7), mean(snp7, na.rm = TRUE), snp7),
         GFRSCYS = if_else(is.na(GFRSCYS), mean(GFRSCYS, na.rm = TRUE), GFRSCYS))

# load single variant result find top SNPs associated with metabolite
load("../Data/STAAR_neg_controls/x100001721/x100001721_Individual_Analysis.Rdata")

known_loci_chr2_1721 <- known_loci_chr2_1721 |>
  left_join(results_individual_analysis, by = c("CHR", "POS", "REF", "ALT")) 

top_cv <- paste0("snp", which.min(known_loci_chr2_1721$pvalue)[1])

# load AM result
am_result_x100001721 <- readRDS("../Data/STAAR_neg_controls/x100001721/admixmap_b1_x100001721_no_genotype.RDS")
la_b1_chr2 <- get_local_ancestry("afr", 
                                 cov_df = cov_b1_chr2_cv_adj, 
                                 chr = 2) 



plan(multisession, workers = max(1, availableCores() - 15))

var_explain_chr2_x100001721 <- fit_var_explained(
  la_obj    = la_b1_chr2,          # the raw list: [[1]] matrix, [[2]] snppos
  am_result = am_result_x100001721,
  n_snp     = 7,
  top_cv    = top_cv,
  dat       = cov_b1_chr2_cv_adj,
  covMat    = covMatList
)

saveRDS(var_explain_chr2_x100001721, "../Data/Variance_explained/b1_chr2_x100001721_cv_afr_variance_explained.RDS")

plan(sequential)
gc()
 
######## chr5 #########
known_loci_chr5_799 <- readRDS("../Data/STAAR_neg_controls/x799/individual_cond_pruned_var.RDS")
cov_b1_chr5 <- get_adj_cv(5, "b1", known_cv = known_loci_chr5_799,
                          cov_b1)

am_result_x799 <- readRDS("../Data/STAAR_neg_controls/x799/admixmap_b1_x799_no_genotype.RDS")
cov_b1_chr5 <-  cov_b1_chr5|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp1, na.rm = TRUE), snp2),
         GFRSCYS = if_else(is.na(GFRSCYS), mean(GFRSCYS, na.rm = TRUE), GFRSCYS))

# load single variant result find top SNPs associated with metabolite
load("../Data/STAAR_neg_controls/x799/x799_Individual_Analysis.Rdata")

known_loci_chr5_799 <- known_loci_chr5_799 |>
  left_join(results_individual_analysis, by = c("CHR", "POS", "REF", "ALT")) 

top_cv <- paste0("snp", which.min(known_loci_chr5_799$pvalue)[1])



la_b1_chr5_x799 <- get_local_ancestry("afr", 
                                 cov_df = cov_b1_chr5, 
                                 chr = 5) 

plan(multisession, workers = max(1, availableCores() - 20))

var_explain_chr5_x799 <- fit_var_explained(
  la_obj    = la_b1_chr5_x799,          # the raw list: [[1]] matrix, [[2]] snppos
  am_result = am_result_x799,
  n_snp     = 2,
  dat       = cov_b1_chr5,
  covMat    = covMatList,
  top_cv = top_cv
)

saveRDS(var_explain_chr5_x799, "../Data/Variance_explained/b1_chr5_x799_cv_afr_variance_explained.RDS")

plan(sequential)
gc()

######## chr6 #########
known_loci_chr6 <- readRDS("../Data/STAAR_neg_controls/x1215/individual_cond_pruned_var.RDS")
cov_b1_chr6 <- get_adj_cv(6, "b1", known_cv = known_loci_chr6,
                          cov_b1)

cov_b1_chr6 <-  cov_b1_chr6|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp1, na.rm = TRUE), snp2),
         GFRSCYS = if_else(is.na(GFRSCYS), mean(GFRSCYS, na.rm = TRUE), GFRSCYS))

# load single variant result find top SNPs associated with metabolite
load("../Data/STAAR_neg_controls/x1215/x1215_Individual_Analysis.Rdata")

known_loci_chr6 <- known_loci_chr6 |>
  left_join(results_individual_analysis, by = c("CHR", "POS", "REF", "ALT")) 

top_cv <- paste0("snp", which.min(known_loci_chr6$pvalue)[1])


am_result_x1215 <- readRDS("../Data/STAAR_neg_controls/x1215/admixmap_b1_x1215_no_genotype.RDS")

la_b1_chr6_x1215 <- get_local_ancestry("amer", 
                                      cov_df = cov_b1_chr6, 
                                      chr = 6) 

plan(multisession, workers = max(1, availableCores() - 10))

var_explain_chr6_x1215 <- fit_var_explained(
  la_obj    = la_b1_chr6_x1215,          # the raw list: [[1]] matrix, [[2]] snppos
  am_result = am_result_x1215,
  n_snp     = 2,
  dat       = cov_b1_chr6,
  covMat    = covMatList,
  top_cv = top_cv
)

saveRDS(var_explain_chr6_x1215, "../Data/Variance_explained/b1_chr6_x1215_cv_amer_variance_explained.RDS")

plan(sequential)
gc()


##############  chr 16 x278 ##############  

known_loci_chr16_278 <- readRDS("../Data/STAAR_neg_controls/x278/individual_cond_pruned_var.RDS")
cov_b1_chr16 <- get_adj_cv(16, "b1", known_cv = known_loci_chr16_278,
                           cov_b1)

# load single variant result find top SNPs associated with metabolite
load("../Data/STAAR_neg_controls/x278/x278_Individual_Analysis.Rdata")

known_loci_chr16_278 <- known_loci_chr16_278 |>
  left_join(results_individual_analysis, by = c("CHR", "POS", "REF", "ALT")) 

top_cv <- paste0("snp", which.min(known_loci_chr16_278$pvalue)[1])


cov_b1_chr16 <-  cov_b1_chr16|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2),
         snp3 = if_else(is.na(snp3), mean(snp3, na.rm = TRUE), snp3),
         GFRSCYS = if_else(is.na(GFRSCYS), mean(GFRSCYS, na.rm = TRUE), GFRSCYS))


la_b1_chr16_x278 <- get_local_ancestry("amer", 
                                      cov_df = cov_b1_chr16, 
                                      chr = 16) 
# load AM result
am_null_x278 <- readRDS("../Data/STAAR_neg_controls/x278/admixmap_b1_x278_no_genotype.RDS")
plan(multisession, workers = max(1, availableCores() - 20))

var_explain_x278 <- fit_var_explained(
  la_obj    = la_b1_chr16_x278,          # the raw list: [[1]] matrix, [[2]] snppos
  am_result = am_null_x278,
  n_snp     = 3,
  dat       = cov_b1_chr16,
  covMat    = covMatList,
  top_cv = top_cv
)

saveRDS(var_explain_x278, "../Data/Variance_explained/b1_chr16_cv_x278_amer_variance_explained.RDS")

plan(sequential)

############## test regions chr8 x1021 ############## 
known_loci_chr8 <- readRDS("../Data/STAAR_exp_groups/x1021/individual_cond_pruned_var.RDS")
cov_b1_chr8 <- get_adj_cv(8, "b1", known_cv = known_loci_chr8,
                          cov_b1)

# load single variant result find top SNPs associated with metabolite
load("../Data/STAAR_exp_groups/x1021/x1021_Individual_Analysis.Rdata")
known_loci_chr8 <- known_loci_chr8 |>
  left_join(results_individual_analysis, by = c("CHR", "POS", "REF", "ALT")) 
top_cv <- paste0("snp", which.min(known_loci_chr8$pvalue)[1])

cov_b1_chr8 <-  cov_b1_chr8|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2),
         snp3 = if_else(is.na(snp3), mean(snp3, na.rm = TRUE), snp3),
         snp4 = if_else(is.na(snp4), mean(snp4, na.rm = TRUE), snp4),
         snp5 = if_else(is.na(snp5), mean(snp5, na.rm = TRUE), snp5),
         snp6 = if_else(is.na(snp6), mean(snp6, na.rm = TRUE), snp6),
         snp7 = if_else(is.na(snp7), mean(snp7, na.rm = TRUE), snp7),
         snp8 = if_else(is.na(snp8), mean(snp8, na.rm = TRUE), snp8),
         GFRSCYS = if_else(is.na(GFRSCYS), mean(GFRSCYS, na.rm = TRUE), GFRSCYS),
         x1021 = if_else(is.na(x1021), min(x1021, na.rm = TRUE), x1021))

la_b1_chr8_x1021 <- get_local_ancestry("afr", 
                                       cov_df = cov_b1_chr8, 
                                       chr = 8) 
# load AM result
am_null_x1021 <- readRDS("../Data/STAAR_exp_groups/x1021/admixmap_b1_x1021_no_genotype.RDS")

plan(multisession, workers = max(1, availableCores() - 20))

var_explain_x1021 <- fit_var_explained(
  la_obj    = la_b1_chr8_x1021,          # the raw list: [[1]] matrix, [[2]] snppos
  am_result = am_null_x1021,
  n_snp     = 8,
  dat       = cov_b1_chr8,
  covMat    = covMatList,
  top_cv = top_cv
)

saveRDS(var_explain_x1021, "../Data/Variance_explained/b1_chr8_cv_x1021_afr_variance_explained.RDS")

plan(sequential)


############## test regions chr10 x100000007 ############## 
known_loci_chr10 <- readRDS("../Data/STAAR_exp_groups/x100000007/individual_cond_pruned_var.RDS")
cov_b1_chr10 <- get_adj_cv(10, "b1", known_cv = known_loci_chr10,
                          cov_b1)

# load single variant result find top SNPs associated with metabolite
load("../Data/STAAR_exp_groups/x100000007/x100000007_Individual_Analysis.Rdata")
known_loci_chr10 <- known_loci_chr10 |>
  left_join(results_individual_analysis, by = c("CHR", "POS", "REF", "ALT")) 
top_cv <- paste0("snp", which.min(known_loci_chr10$pvalue)[1])

cov_b1_chr10 <-  cov_b1_chr10|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2),
         snp3 = if_else(is.na(snp3), mean(snp3, na.rm = TRUE), snp3),
         snp4 = if_else(is.na(snp4), mean(snp4, na.rm = TRUE), snp4),
         snp5 = if_else(is.na(snp5), mean(snp5, na.rm = TRUE), snp5),
         GFRSCYS = if_else(is.na(GFRSCYS), mean(GFRSCYS, na.rm = TRUE), GFRSCYS))

la_b1_chr10_x100000007 <- get_local_ancestry("amer", 
                                       cov_df = cov_b1_chr10, 
                                       chr = 10) 
# load AM result
am_null_x100000007 <- readRDS("../Data/STAAR_exp_groups/x100000007/admixmap_b1_x100000007_no_genotype.RDS")


plan(multisession, workers = max(1, availableCores() - 20))

var_explain_x100000007 <- fit_var_explained(
  la_obj    = la_b1_chr10_x100000007,          # the raw list: [[1]] matrix, [[2]] snppos
  am_result = am_null_x100000007,
  n_snp     = 5,
  dat       = cov_b1_chr10,
  covMat    = covMatList, top_cv = top_cv
)

saveRDS(var_explain_x100000007,
        "../Data/Variance_explained/b1_chr10_cv_x100000007_amer_variance_explained.RDS")

plan(sequential)

############## test regions chr13 x100004046 ############## 
known_loci_chr13 <- readRDS("../Data/STAAR_exp_groups/x100004046/individual_cond_pruned_var.RDS")
cov_b1_chr13 <- get_adj_cv(13, "b1", known_cv = known_loci_chr13,
                           cov_b1)

# load single variant result find top SNPs associated with metabolite
load("../Data/STAAR_exp_groups/x100004046/x100004046_Individual_Analysis.Rdata")
known_loci_chr13 <- known_loci_chr13 |>
  left_join(results_individual_analysis, by = c("CHR", "POS", "REF", "ALT")) 
top_cv <- paste0("snp", which.min(known_loci_chr13$pvalue)[1])

cov_b1_chr13 <-  cov_b1_chr13|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2),
         snp3 = if_else(is.na(snp3), mean(snp3, na.rm = TRUE), snp3),
         GFRSCYS = if_else(is.na(GFRSCYS), mean(GFRSCYS, na.rm = TRUE), GFRSCYS))

la_b1_chr13_x100004046 <- get_local_ancestry("afr", 
                                             cov_df = cov_b1_chr13, 
                                             chr = 13) 
# load AM result
am_null_x100004046 <- readRDS("../Data/STAAR_exp_groups/x100004046/admixmap_b1_x100004046_no_genotype.RDS")

plan(multisession, workers = max(1, availableCores() - 15))

var_explain_x100004046 <- fit_var_explained(
  la_obj    = la_b1_chr13_x100004046,          
  am_result = am_null_x100004046,
  n_snp     = 3,
  dat       = cov_b1_chr13,
  covMat    = covMatList, top_cv = top_cv
)

saveRDS(var_explain_x100004046,
        "../Data/Variance_explained/b1_chr13_cv_x100004046_afr_variance_explained.RDS")

plan(sequential)

############## test regions chr16 x1224 ############## 
known_loci_chr16_1224 <- readRDS("../Data/STAAR_exp_groups/x1224/individual_cond_pruned_var.RDS")
cov_b1_chr16_x1224 <- get_adj_cv(16, "b1", known_cv = known_loci_chr16_1224,
                           cov_b1)

# load single variant result find top SNPs associated with metabolite
load("../Data/STAAR_exp_groups/x1224/x1224_Individual_Analysis.Rdata")

known_loci_chr16_1224 <- known_loci_chr16_1224 |>
  left_join(results_individual_analysis, by = c("CHR", "POS", "REF", "ALT")) 

top_cv <- paste0("snp", which.min(known_loci_chr16_1224$pvalue)[1])

cov_b1_chr16_x1224 <-  cov_b1_chr16_x1224|>
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

la_b1_chr16_x1224 <- get_local_ancestry("amer", cov_df = cov_b1_chr16_x1224, 
                                             chr = 16) 
# load AM result
am_null_x1224 <- readRDS("../Data/STAAR_exp_groups/x1224/admixmap_b1_x1224_no_genotype.RDS")


plan(multisession, workers = max(1, availableCores() - 20))

var_explain_x1224 <- fit_var_explained(
  la_obj    = la_b1_chr16_x1224,          # the raw list: [[1]] matrix, [[2]] snppos
  am_result = am_null_x1224,
  n_snp     = 11,
  dat       = cov_b1_chr16_x1224,
  covMat    = covMatList,
  top_cv = top_cv
)

saveRDS(var_explain_x1224,
        "../Data/Variance_explained/b1_chr16_cv_x1224_amer_variance_explained.RDS")

plan(sequential)


# load gds files catalog that need to be accessed 
load("../Data/STAAR_prep/agds_dir_scale_up_final.Rdata") #variable name agds_dir

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
############## x100009332, chr11 arachido ##############
known_loci_chr11_arachido <- readRDS("../Data/STAAR_exp_groups/x100009332/individual_cond_pruned_var.RDS")
cov_b1_chr11_arachido <- get_adj_cv(11, "b1", known_cv = known_loci_chr11_arachido,
                                    cov_b1, i = 1)

cov_b1_chr11_arachido <-  cov_b1_chr11_arachido|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1))

# load single variant result find top SNPs associated with metabolite
load("../Data/STAAR_exp_groups/x100009332/x100009332_Individual_Analysis.Rdata")

known_loci_chr11_arachido <- known_loci_chr11_arachido |>
  left_join(results_individual_analysis, by = c("CHR", "POS", "REF", "ALT")) 

top_cv <- paste0("snp", which.min(known_loci_chr11_arachido$pvalue)[1])

la_b1_chr11_x100009332 <- get_local_ancestry("amer", 
                                        cov_df = cov_b1_chr11_arachido, 
                                        chr = 11) 
# load AM result
am_null_x100009332 <- readRDS("../Data/STAAR_exp_groups/x100009332/admixmap_b1_x100009332_no_genotype.RDS")

plan(multisession, workers = max(1, availableCores() - 10))

var_explain_x100009332 <- fit_var_explained(
  la_obj    = la_b1_chr11_x100009332,          # the raw list: [[1]] matrix, [[2]] snppos
  am_result = am_null_x100009332,
  n_snp     = 9,
  dat       = cov_b1_chr11_arachido,
  covMat    = covMatList, top_cv = top_cv
)

saveRDS(var_explain_x100009332,
        "../Data/Variance_explained/b1_chr11_cv_x100009332_amer_variance_explained.RDS")

plan(sequential)


############## x100006370, chr11 3-beta ##############
known_loci_chr11_3beta <- readRDS("../Data/STAAR_exp_groups/x100006370/individual_cond_pruned_var.RDS")
cov_b1_chr11_3beta <- get_adj_cv(11, "b1", known_cv = known_loci_chr11_3beta,
                                 cov_b1, i = 2)

cov_b1_chr11_3beta <-  cov_b1_chr11_3beta|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1))

# load single variant result find top SNPs associated with metabolite
load("../Data/STAAR_exp_groups/x100006370/x100006370_Individual_Analysis.Rdata")

known_loci_chr11_3beta <- known_loci_chr11_3beta |>
  left_join(results_individual_analysis, by = c("CHR", "POS", "REF", "ALT")) 

top_cv <- paste0("snp", which.min(known_loci_chr11_3beta$pvalue)[1])

la_b1_chr11_x100006370 <- get_local_ancestry("amer", 
                                             cov_df = cov_b1_chr11_3beta, 
                                             chr = 11) 

# load AM result
am_null_x100006370 <- readRDS("../Data/STAAR_exp_groups/x100006370/admixmap_b1_x100006370_no_genotype.RDS")

plan(multisession, workers = max(1, availableCores() - 15))

var_explain_x100006370 <- fit_var_explained(
  la_obj    = la_b1_chr11_x100006370,          # the raw list: [[1]] matrix, [[2]] snppos
  am_result = am_null_x100006370,
  n_snp     = 5,
  dat       = cov_b1_chr11_3beta,
  covMat    = covMatList,
  top_cv = top_cv
)

saveRDS(var_explain_x100006370,
        "../Data/Variance_explained/b1_chr11_cv_x100006370_amer_variance_explained.RDS")

plan(sequential)


############## agds_dir[3], chr12 1_methyl, x100001208 ##############
known_loci_chr12_1methyl <- readRDS("../Data/STAAR_exp_groups/x100001208/individual_cond_pruned_var.RDS")
cov_b1_chr12_1methyl <- get_adj_cv(12, "b1", known_cv = known_loci_chr12_1methyl,
                                 cov_b1, i = 3)

cov_b1_chr12_1methyl <-  cov_b1_chr12_1methyl|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2),
         snp3 = if_else(is.na(snp3), mean(snp3, na.rm = TRUE), snp3),
         snp4 = if_else(is.na(snp4), mean(snp4, na.rm = TRUE), snp4),
         snp5 = if_else(is.na(snp5), mean(snp5, na.rm = TRUE), snp5),
         snp6 = if_else(is.na(snp6), mean(snp6, na.rm = TRUE), snp6))

# load single variant result find top SNPs associated with metabolite
load("../Data/STAAR_exp_groups/x100001208/x100001208_Individual_Analysis.Rdata")

known_loci_chr12_1methyl <- known_loci_chr12_1methyl |>
  left_join(results_individual_analysis, by = c("CHR", "POS", "REF", "ALT")) 

top_cv <- paste0("snp", which.min(known_loci_chr12_1methyl$pvalue)[1])

la_b1_chr12_x100001208 <- get_local_ancestry("amer", 
                                             cov_df = cov_b1_chr12_1methyl, 
                                             chr = 12) 
# load AM result
am_null_x100001208 <- readRDS("../Data/STAAR_exp_groups/x100001208/admixmap_b1_x100001208_no_genotype.RDS")

plan(multisession, workers = max(1, availableCores() - 15))

var_explain_x100001208 <- fit_var_explained(
  la_obj    = la_b1_chr12_x100001208,          # the raw list: [[1]] matrix, [[2]] snppos
  am_result = am_null_x100001208,
  n_snp     = 6,
  dat       = cov_b1_chr12_1methyl,
  covMat    = covMatList,
  top_cv = top_cv
)

saveRDS(var_explain_x100001208,
        "../Data/Variance_explained/b1_chr12_cv_x100001208_amer_variance_explained.RDS")

plan(sequential)



############## agds_dir[4], chr12 Ethylmalonate, x2054 ##############
known_loci_chr12_x2054 <- readRDS("../Data/STAAR_exp_groups/x2054/individual_cond_pruned_var.RDS")
cov_b1_chr12_x2054 <- get_adj_cv(12, "b1", known_cv = known_loci_chr12_x2054,
                                   cov_b1, i = 4)
 
cov_b1_chr12_x2054 <-  cov_b1_chr12_x2054|>
  mutate(snp1 = if_else(is.na(snp1), mean(snp1, na.rm = TRUE), snp1),
         snp2 = if_else(is.na(snp2), mean(snp2, na.rm = TRUE), snp2))

# load single variant result find top SNPs associated with metabolite
load("../Data/STAAR_exp_groups/x2054/x2054_Individual_Analysis.Rdata")

known_loci_chr12_x2054 <- known_loci_chr12_x2054 |>
  left_join(results_individual_analysis, by = c("CHR", "POS", "REF", "ALT")) 

top_cv <- paste0("snp", which.min(known_loci_chr12_x2054$pvalue)[1])

la_b1_chr12_x2054 <- get_local_ancestry("afr", 
                                             cov_df = cov_b1_chr12_x2054, 
                                             chr = 12) 
# load AM result
am_null_x2054 <- readRDS("../Data/STAAR_exp_groups/x2054/admixmap_b1_x2054_no_genotype.RDS")

plan(multisession, workers = max(1, availableCores() - 10))

var_explain_x2054 <- fit_var_explained(
  la_obj    = la_b1_chr12_x2054,          # the raw list: [[1]] matrix, [[2]] snppos
  am_result = am_null_x2054,
  n_snp     = 2,
  dat       = cov_b1_chr12_x2054,
  covMat    = covMatList, , top_cv = top_cv
)

saveRDS(var_explain_x2054,
        "../Data/Variance_explained/b1_chr12_cv_x2054_afr_variance_explained.RDS")

plan(sequential)

###### x6264 chr16  #####
known_loci_chr16_propyl <- readRDS("../Data/STAAR_exp_groups/x100006264/individual_cond_pruned_var.RDS")
cov_b1_chr16_propyl <- get_adj_cv(16, "b1", known_cv = known_loci_chr16_propyl,
                                  cov_b1, i =5)
cov_b1_chr16_propyl <- cov_b1_chr16_propyl |>
  mutate(snp8 = if_else(is.na(snp8), mean(snp8, na.rm = TRUE), snp8),
         snp9 = if_else(is.na(snp9), mean(snp9, na.rm = TRUE), snp9),
         snp10 = if_else(is.na(snp10), mean(snp10, na.rm = TRUE), snp10),
         snp17 = if_else(is.na(snp17), mean(snp17, na.rm = TRUE), snp17))

# load single variant result find top SNPs associated with metabolite
load("../Data/STAAR_exp_groups/x100006264/x100006264_Individual_Analysis.Rdata")

known_loci_chr16_propyl <- known_loci_chr16_propyl |>
  left_join(results_individual_analysis, by = c("CHR", "POS", "REF", "ALT")) 

top_cv <- paste0("snp", which.min(known_loci_chr16_propyl$pvalue)[1])


la_b1_chr16_x6264 <- get_local_ancestry("afr", 
                                             cov_df = cov_b1_chr16_propyl, 
                                             chr = 16) 
# load AM result
am_null_x6264 <- readRDS("../Data/STAAR_exp_groups/x100006264/admixmap_b1_x100006264_no_genotype.RDS")

plan(multisession, workers = max(1, availableCores() - 4))

var_explain_x6264 <- fit_var_explained(
  la_obj    = la_b1_chr16_x6264,          # the raw list: [[1]] matrix, [[2]] snppos
  am_result = am_null_x6264,
  n_snp     = 18,
  dat       = cov_b1_chr16_propyl,
  covMat    = covMatList, top_cv = top_cv
)

saveRDS(var_explain_x6264,
        "../Data/Variance_explained/b1_chr16_cv_x100006264_afr_variance_explained.RDS")

plan(sequential)

###### x100001266 chr2  #####
known_loci_chr2_1266 <- readRDS("../Data/STAAR_exp_groups/x100001266/individual_cond_pruned_var.RDS")
cov_b1_known_loci_chr2 <- get_adj_cv(2, "b1", known_cv = known_loci_chr2_1266,
                                     cov_b1, i = 6)

# mean imputation
cov_b1_known_loci_chr2 <- cov_b1_known_loci_chr2 |>
  mutate(snp8 = if_else(is.na(snp8), mean(snp8, na.rm = TRUE), snp8))

# load single variant result find top SNPs associated with metabolite
load("../Data/STAAR_exp_groups/x100001266/x100001266_Individual_Analysis.Rdata")

known_loci_chr2_1266 <- known_loci_chr2_1266 |>
  left_join(results_individual_analysis, by = c("CHR", "POS", "REF", "ALT")) 

top_cv <- paste0("snp", which.min(known_loci_chr2_1266$pvalue)[1])

la_b1_chr2_x100001266 <- get_local_ancestry("afr", 
                                        cov_df = cov_b1_known_loci_chr2, 
                                        chr = 2) 
# load AM result
am_null_x100001266 <- readRDS("../Data/STAAR_exp_groups/x100001266/admixmap_b1_x100001266_no_genotype.RDS")

plan(multisession, workers = max(1, availableCores() - 25))

var_explain_x100001266 <- fit_var_explained(
  la_obj    = la_b1_chr2_x100001266,          # the raw list: [[1]] matrix, [[2]] snppos
  am_result = am_null_x100001266,
  n_snp     = 10,
  dat       = cov_b1_known_loci_chr2,
  covMat    = covMatList, top_cv = top_cv
)

saveRDS(var_explain_x100001266,
        "../Data/Variance_explained/b1_chr2_cv_var_explain_x100001266_afr_variance_explained.RDS")

plan(sequential)

###### x1114 chr5  #####
known_loci_chr5_1114 <- readRDS("../Data/STAAR_exp_groups/x1114/individual_cond_pruned_var.RDS")
cov_b1_known_loci_chr5 <- get_adj_cv(5, "b1", known_cv = known_loci_chr5_1114,
                                     cov_b1, i = 7)

# load single variant result find top SNPs associated with metabolite
load("../Data/STAAR_exp_groups/x1114/x1114_Individual_Analysis.Rdata")

known_loci_chr5_1114 <- known_loci_chr5_1114 |>
  left_join(results_individual_analysis, by = c("CHR", "POS", "REF", "ALT")) 

top_cv <- paste0("snp", which.min(known_loci_chr5_1114$pvalue)[1])


la_b1_chr5_x1114 <- get_local_ancestry("amer", 
                                            cov_df = cov_b1_known_loci_chr5, 
                                            chr = 5) 
# load AM result
am_null_x1114 <- readRDS("../Data/STAAR_exp_groups/x1114/admixmap_b1_x1114_no_genotype.RDS")

plan(multisession, workers = max(1, availableCores() - 10))

var_explain_x1114 <- fit_var_explained(
  la_obj    = la_b1_chr5_x1114,          # the raw list: [[1]] matrix, [[2]] snppos
  am_result = am_null_x1114,
  n_snp     = 10,
  dat       = cov_b1_known_loci_chr5,
  covMat    = covMatList, top_cv = top_cv
)

saveRDS(var_explain_x1114,
        "../Data/Variance_explained/b1_chr5_cv_var_explain_x1114_amer_variance_explained.RDS")

plan(sequential)

###### x192 chr8  #####
known_loci_chr8_192 <- readRDS("../Data/STAAR_exp_groups/x192/individual_cond_pruned_var.RDS")
cov_b1_known_loci_chr8 <- get_adj_cv(8, "b1", known_cv = known_loci_chr8_192,
                                     cov_b1, i = 8)


# load single variant result find top SNPs associated with metabolite
load("../Data/STAAR_exp_groups/x192/x192_Individual_Analysis.Rdata")

known_loci_chr8_192 <- known_loci_chr8_192 |>
  left_join(results_individual_analysis, by = c("CHR", "POS", "REF", "ALT")) 

top_cv <- paste0("snp", which.min(known_loci_chr8_192$pvalue)[1])


la_b1_chr8_x192 <- get_local_ancestry("amer", 
                                       cov_df = cov_b1_known_loci_chr8, 
                                       chr = 8) 
# load AM result
am_null_x192 <- readRDS("../Data/STAAR_exp_groups/x192/admixmap_b1_x192_no_genotype.RDS")

plan(multisession, workers = max(1, availableCores() - 8))

var_explain_x192 <- fit_var_explained(
  la_obj    = la_b1_chr8_x192,          # the raw list: [[1]] matrix, [[2]] snppos
  am_result = am_null_x192,
  n_snp     = 6,
  dat       = cov_b1_known_loci_chr8,
  covMat    = covMatList, top_cv = top_cv
)

saveRDS(var_explain_x192,
        "../Data/Variance_explained/b1_chr8_cv_var_explain_x192_amer_variance_explained.RDS")

plan(sequential)

all_associated_variants <- rbind(known_loci_chr2_1721, known_loci_chr5_799,
                                  known_loci_chr6, known_loci_chr16_278,
                                  known_loci_chr2_1266, known_loci_chr5_1114,
                                  known_loci_chr8_192, known_loci_chr8,
                                  known_loci_chr10, known_loci_chr11_arachido,
                                  known_loci_chr11_3beta, known_loci_chr12_1methyl,
                                  known_loci_chr12_x2054, known_loci_chr13,
                                  known_loci_chr16_propyl, known_loci_chr16_1224)

data.table::fwrite(all_associated_variants, "../Data/all_associated_common_variants.csv")
