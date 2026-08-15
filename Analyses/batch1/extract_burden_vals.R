#==============================================================================
# Extract burden / score values from conditional STAAR results - BATCH 1
#==============================================================================
# MERGED FROM (chronological):
#   - 20250205_extract_burden_vals.R   ->  section "v1_original"
#   - 20250810_extract_burden_vals_scaled_up.R   ->  section "v2_scaled_up"
# NOTE: v2 additionally defines the coding+noncoding dispatchers and toDataframe();
# NOTE: those function bodies are duplicated in shared/functions/burden_extraction/.
#==============================================================================

#==============================================================================
# MAIN PIPELINE  (current version: v2_scaled_up)
#==============================================================================

# import pacakges
library(tidyverse)
library(SeqArray)
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)
library(SeqVarTools)

Gene_Centric_Noncoding_cond_extract_burden <- function(chr,gene_name,category=c("downstream","upstream","UTR","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"),
                                                       genofile,obj_nullmodel,known_loci=NULL,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                       rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                                                       method_cond=c("optimal","naive"),
                                                       QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                                       Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                                       Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL){

  ## evaluate choices
  category <- match.arg(category)
  method_cond <- match.arg(method_cond)
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)

  if(is.null(known_loci))
  {
    known_loci <- data.frame(chr=logical(0),pos=logical(0),ref=character(0),alt=character(0))
  }

  if(category=="downstream")
  {
    results <- downstream_cond_extract_burden(chr,gene_name,genofile,obj_nullmodel,
                                              known_loci,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                              rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                              method_cond=method_cond,
                                              QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                              Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                              Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  }

  if(category=="upstream")
  {
    results <- upstream_cond_extract_burden(chr,gene_name,genofile,obj_nullmodel,
                                            known_loci,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                            rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                            method_cond=method_cond,
                                            QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                            Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                            Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  }

  if(category=="UTR")
  {
    results <- UTR_cond_extract_burden(chr,gene_name,genofile,obj_nullmodel,
                                       known_loci,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                       rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                       method_cond=method_cond,
                                       QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                       Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                       Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  }

  if(category=="promoter_CAGE")
  {
    results <- promoter_CAGE_cond_extract_burden(chr,gene_name,genofile,obj_nullmodel,
                                                 known_loci,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                                 rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                                 method_cond=method_cond,
                                                 QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                 Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                 Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  }

  if(category=="promoter_DHS")
  {
    results <- promoter_DHS_cond_extract_burden(chr,gene_name,genofile,obj_nullmodel,
                                                known_loci,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                                rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                                method_cond=method_cond,
                                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  }

  if(category=="enhancer_CAGE")
  {
    results <- enhancer_CAGE_cond_extract_burden(chr,gene_name,genofile,obj_nullmodel,
                                                 known_loci,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                                 rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                                 method_cond=method_cond,
                                                 QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                 Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                 Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  }

  if(category=="enhancer_DHS")
  {
    results <- enhancer_DHS_cond_extract_burden(chr,gene_name,genofile,obj_nullmodel,
                                                known_loci,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                                rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                                method_cond=method_cond,
                                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  }

  return(results)
}


Gene_Centric_Coding_cond_extract_burden <- function(chr,gene_name,category=c("plof","plof_ds","missense","disruptive_missense","synonymous","ptv","ptv_ds"),
                                                    genofile,obj_nullmodel,known_loci=NULL,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                    rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                                                    method_cond=c("optimal","naive"),
                                                    QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                                    Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                                    Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL){

  ## evaluate choices
  category <- match.arg(category)
  method_cond <- match.arg(method_cond)
  variant_type <- match.arg(variant_type)
  geno_missing_imputation <- match.arg(geno_missing_imputation)

  genes <- genes_info[genes_info[,2]==chr,]
  if(is.null(known_loci))
  {
    known_loci <- data.frame(chr=logical(0),pos=logical(0),ref=character(0),alt=character(0))
  }

  if(category=="plof")
  {
    results <- plof_cond_extract_burden(chr,gene_name,genofile,obj_nullmodel,genes,
                                        known_loci,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                        rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                        method_cond=method_cond,
                                        QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  }

  if(category=="plof_ds")
  {
    results <- plof_ds_cond_extract_burden(chr,gene_name,genofile,obj_nullmodel,genes,
                                           known_loci,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                           rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                           method_cond=method_cond,
                                           QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                           Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                           Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  }

  if(category=="missense")
  {
    results <- missense_cond_extract_burden(chr,gene_name,genofile,obj_nullmodel,genes,
                                            known_loci,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                            rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                            method_cond=method_cond,
                                            QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                            Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                            Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
    print(paste0("genecentric coding func missense result",
                 dim(missense)))
  }

  if(category=="disruptive_missense")
  {
    results <- disruptive_missense_cond_extract_burden(chr,gene_name,genofile,obj_nullmodel,genes,
                                                       known_loci,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                                       rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                                       method_cond=method_cond,
                                                       QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                       Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                       Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  }

  if(category=="synonymous")
  {
    results <- synonymous_cond_extract_burden(chr,gene_name,genofile,obj_nullmodel,genes,
                                              known_loci,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                              rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                              method_cond=method_cond,
                                              QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                              Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                              Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  }

  if(category=="ptv")
  {
    results <- ptv_cond_extract_burden(chr,gene_name,genofile,obj_nullmodel,genes,
                                       known_loci,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                       rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                       method_cond=method_cond,
                                       QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                       Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                       Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  }

  if(category=="ptv_ds")
  {
    results <- ptv_ds_cond_extract_burden(chr,gene_name,genofile,obj_nullmodel,genes,
                                          known_loci,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                          rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                          method_cond=method_cond,
                                          QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                          Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                          Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  }

  return(results)
}


STAAR_cond_burden_score_only <- function(genotype,genotype_adj,obj_nullmodel,annotation_phred=NULL,
                                         rare_maf_cutoff=0.01,rv_num_cutoff=2,rv_num_cutoff_max=1e9,
                                         method_cond=c("optimal","naive")){
  print("running STAAR_cond burden score only function")
  method_cond <- match.arg(method_cond) # evaluate choices
  if(!inherits(genotype, "matrix") && !inherits(genotype, "Matrix")){
    stop("genotype is not a matrix!")
  }

  if(dim(genotype)[2] == 1){
    stop(paste0("Number of rare variant in the set is less than 2!"))
  }
  annotation_phred <- as.data.frame(annotation_phred)
  if(dim(annotation_phred)[1] != 0 & dim(genotype)[2] != dim(annotation_phred)[1]){
    stop(paste0("Dimensions don't match for genotype and annotation!"))
  }

  if(inherits(genotype, "sparseMatrix")){
    genotype <- as.matrix(genotype)
  }

  if(inherits(genotype_adj, "numeric")){
    genotype_adj <- matrix(genotype_adj, ncol=1)
  }

  if(dim(genotype)[1] != dim(genotype_adj)[1]){
    stop(paste0("Dimensions don't match for genotype and genotype_adj!"))
  }
  sampleids <- rownames(genotype)
  genotype <- matrix_flip(genotype)
  MAF <- genotype$MAF
  RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0))
  print(RV_label)
  Geno_rare <- genotype$Geno[,RV_label]
  rownames(Geno_rare) <- sampleids
  #print(head(rownames(Geno_rare)))
  rm(genotype)
  gc()
  annotation_phred <- annotation_phred[RV_label,,drop=FALSE]
  if(sum(RV_label) >= rv_num_cutoff_max){
    stop(paste0("Number of rare variant in the set is more than ",rv_num_cutoff_max,"!"))
  }

  if(sum(RV_label) >= rv_num_cutoff){
    G <- as(Geno_rare,"dgCMatrix")
    MAF <- MAF[RV_label]
    rm(Geno_rare)
    gc()

    annotation_rank <- 1 - 10^(-annotation_phred/10)

    ## beta(1,25)
    w_1 <- dbeta(MAF,1,25)
    ## beta(1,1)
    w_2 <- dbeta(MAF,1,1)
    if(dim(annotation_phred)[2] == 0){
      ## Burden, SKAT, ACAT-V
      w_B <- w_S <- as.matrix(cbind(w_1,w_2))
      #w_A <- as.matrix(cbind(w_1^2/dbeta(MAF,0.5,0.5)^2,w_2^2/dbeta(MAF,0.5,0.5)^2))
    }else{
      ## Burden
      w_B_1 <- annotation_rank*w_1
      w_B_1 <- cbind(w_1,w_B_1)
      colnames(w_B_1) <- paste0("Burden(1, 25)-", colnames(w_B_1))
      colnames(w_B_1)[1] <- "Burden(1, 25)"

      w_B_2 <- annotation_rank*w_2
      w_B_2 <- cbind(w_2,w_B_2)
      colnames(w_B_2) <- paste0("Burden(1, 1)-", colnames(w_B_2))
      colnames(w_B_2)[1] <- "Burden(1, 1)"

      w_B <- cbind(w_B_1,w_B_2)
      w_B <- as.matrix(w_B)
      #print(paste0("the column names of burden weights are:", colnames(w_B)))

    }

    # Assuming residuals and G are Armadillo vectors/matrices
    print(paste0("dim G", dim(G)))
    print(paste0("dim w_B", dim(w_B)))
    burden <- G %*% w_B

    print(paste0("dim burden", dim(burden)))
    return(burden)

  }else{
    stop(paste0("Number of rare variant in the set is less than ",rv_num_cutoff,"!"))
  }
}


toDataframe <- function(dat){
  dat_df <- as.data.frame(matrix(unlist(dat), ncol = 91))
  colnames(dat_df) <- colnames(dat)
  return(dat_df)
}

## QC_label
QC_label <- "annotation/filter"
## geno_missing_imputation
geno_missing_imputation <- "mean"
## variant_type
variant_type <- "SNV"
## method_cond
method_cond <- "optimal"
## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
Annotation_name_catalog <- get(load("../Data/STAAR_prep/Annotation_name_catalog.Rdata"))
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")


##### chr 2 -- N2-acetyllysine #####
known_loci_chr2 <- as.data.frame(readRDS("../Data/STAAR_neg_controls/x100001721/individual_cond_pruned_var.RDS"))
genofile <- seqOpen("../Data/chr2_DRAGEN_agds_b1_final.gds")
load("../Data/STAAR_neg_controls/x100001721/x100001721_nullmodel_batch1.Rdata")
chr <- 2
# conditional analysis for the four unique rv set
# synonymous
source("./20250429_synonymous_extract_burden.R")

gene_names <- c("STAMBP", "ALMS1", "EGR4")
burdens_synonymous <- c()
for(gene_name in gene_names){
  burden_test<- Gene_Centric_Coding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                        rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                        known_loci=known_loci_chr2,
                                                        QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                        category="synonymous")
  burden_test <- as.matrix(burden_test)
  burdens_synonymous <- cbind(burdens_synonymous, burden_test[,1])
}

colnames(burdens_synonymous) <- paste0(gene_names, "_synonymous")
saveRDS(burdens_synonymous, "../Data/STAAR_neg_controls/x100001721/burdens_synonymous.RDS")

gene_names <- c("CYP26B1", "ALMS1", "NAT8")
source("./20250425_missense_extract_burden.R")
burdens_missense <- c()
for(gene_name in gene_names){
  burden_test<- Gene_Centric_Coding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                        rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                        known_loci=known_loci_chr2,
                                                        QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                        category="missense")
  burden_test <- as.matrix(burden_test)
  burdens_missense <- cbind(burdens_missense, burden_test[,1])
}

colnames(burdens_missense) <- paste0(gene_names, "_missense")
saveRDS(burdens_missense, "../Data/STAAR_neg_controls/x100001721/burdens_missense.RDS")


gene_names <- c("AUP1")
source("./20250425_disruptive_missense_extract_burden.R")
burdens_disruptive_missense <- c()
for(gene_name in gene_names){
  burden_test<- Gene_Centric_Coding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                        rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                        known_loci=known_loci_chr2,
                                                        QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                        category="disruptive_missense")
  burden_test <- as.matrix(burden_test)
  burdens_disruptive_missense <- cbind(burdens_disruptive_missense, burden_test[,1])
}

colnames(burdens_disruptive_missense) <- paste0(gene_names, "_disruptive_missense")
saveRDS(burdens_disruptive_missense, "../Data/STAAR_neg_controls/x100001721/burdens_disruptive_missense.RDS")

gene_names <- ("AUP1")
burdens_plof_ds <- c()
for(gene_name in gene_names){
  burden_test<- Gene_Centric_Coding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                        rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                        known_loci=known_loci_chr2,
                                                        QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                        category="plof_ds")
  burden_test <- as.matrix(burden_test)
  burdens_plof_ds <- cbind(burdens_plof_ds, burden_test[,1])
}

colnames(burdens_plof_ds) <- paste0(gene_names, "_plof_ds")
saveRDS(burdens_plof_ds, "../Data/STAAR_neg_controls/x100001721/burdens_plof_ds.RDS")


combined_burden <- cbind(burdens_plof_ds, burdens_missense,
                         burdens_disruptive_missense,
                         burdens_synonymous)

corr_burden <- cor(combined_burden)


# Heatmap with ComplexHeatmap
library(ComplexHeatmap)
library(circlize)
Heatmap(corr_burden,
        name = "Correlation",
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        row_names_gp = gpar(fontsize = 8),       # row label size
        column_names_gp = gpar(fontsize = 8),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", corr_burden[i, j]), x, y, gp = gpar(fontsize = 8))})
# high correlation of AUP_disruptive_missesne and AUP1_plof_ds

saveRDS(combined_burden, "../Data/STAAR_neg_controls/x100001721/burdens_combined.RDS")

##### chr 2 -- N-acetylarginine #####
chr2_noncoding_sig <- readRDS("../Data/STAAR_exp_groups/x100001266/noncoding_cond_sig.RDS")
gene_name <- c("DUSP11")
known_loci_chr2 <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x100001266/individual_cond_pruned_var.RDS"))
genofile <- seqOpen("../Data/chr2_scale_up_DRAGEN_agds_b1_final_1.gds")
load("../Data/STAAR_exp_groups/x100001266/x100001266_nullmodel_batch1.Rdata")
# enhancer_CAGE
source("./20250508_enhancer_cage_extract_burden.R")
burden_enhancer_cage <- Gene_Centric_Noncoding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                         rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                         known_loci=known_loci_chr2,
                                                         QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                         Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                         category="enhancer_CAGE")
table(burden_enhancer_cage[,1])
burden_enhancer_cage <- as.matrix(burden_enhancer_cage)

saveRDS(burden_enhancer_cage, "../Data/STAAR_exp_groups/x100001266/burden_enhancer_cage.RDS")
burden_enhancer_cage <- readRDS("../Data/STAAR_exp_groups/x100001266/burden_enhancer_cage.RDS")

combined_burden_chr2 <- cbind(burden_enhancer_cage[,1])
colnames(combined_burden_chr2) <- c("DUSP11_enhancer_CAGE")

saveRDS(combined_burden_chr2, "../Data/STAAR_exp_groups/x100001266/burdens_combined.RDS")

##### chr 5 #####
chr5_coding_sig <- readRDS("../Data/STAAR_exp_groups/x1114/coding_cond_sig.RDS")
gene_name <- c("AGXT2")
known_loci_chr5 <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x1114/individual_cond_pruned_var.RDS"))
genofile <- seqOpen("../Data/chr5_scale_up_DRAGEN_agds_b1_final_1.gds")
load("../Data/STAAR_exp_groups/x1114/x1114_nullmodel_batch1.Rdata")
chr <- 5
# conditional analysis for the unique rv set
# since all of them are located on the same gene, we loop different functional categories

categories_chr5 <- unique(unlist(chr5_coding_sig[,"Category"]))

source("./20250425_plof_ds_extract_burden.R")
# source("./20250425_ptv_ds_extract_burden.R")
burden_chr_5_plof_ds <- Gene_Centric_Coding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                                   rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                                   known_loci=known_loci_chr5,
                                                                   QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                                   Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                                   Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                                   category="plof_ds")

burden_chr_5_plof_ds <- as.matrix(burden_chr_5_plof_ds)
burden_chr_5_dm <- Gene_Centric_Coding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                        rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                        known_loci=known_loci_chr5,
                                                        QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                        category="disruptive_missense")


table(burden_chr_5_dm[,1])
burden_chr_5_dm <- as.matrix(burden_chr_5_dm)

# Only one RV set is used to reduce colinearlity
corr_burden_chr5 <- cor(cbind(burden_chr_5_plof_ds[,1], burden_chr_5_dm[,1]))
colnames(corr_burden_chr5) <- c("AGXT2_plof_ds", "AGXT2_disruptive_missense")
rownames(corr_burden_chr5) <- c("AGXT2_plof_ds", "AGXT2_disruptive_missense")
Heatmap(corr_burden_chr5,
        name = "Correlation",
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        row_names_gp = gpar(fontsize = 8),       # row label size
        column_names_gp = gpar(fontsize = 8),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", corr_burden_chr5[i, j]), x, y, gp = gpar(fontsize = 8))})

burden_chr_5 <- as.data.frame(burden_chr_5[,1])
colnames(burden_chr_5) <- "AGXT2_plof_ds"
saveRDS(burden_chr_5, "../Data/STAAR_exp_groups/x1114/burden_plof_ds.RDS")
seqClose(genofile)

##### batch 2 chr 5 - coding #####
chr5_coding_sig_b2 <- readRDS("../Data/STAAR_exp_groups/x1114/b2/coding_cond_sig.RDS")
genofile <- seqOpen("../Data/chr5_scale_up_DRAGEN_agds_b2_final_1.gds")
load("../Data/STAAR_exp_groups/x1114/x1114_nullmodel_batch2.Rdata")
categories_chr5_b2 <- unique(unlist(chr5_coding_sig_b2[,"Category"]))
burden_chr_5_b2_plof_ds <- Gene_Centric_Coding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                        rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                        known_loci=known_loci_chr5,
                                                        QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                        category="plof_ds")

burden_chr_5_b2_plof_ds <- as.data.frame(burden_chr_5_b2_plof_ds[,1])
saveRDS(burden_chr_5_b2_plof_ds, "../Data/STAAR_exp_groups/x1114/b2/burden_plof_ds.RDS")
##### chr 8 -- 5-Oxoproline #####
gene_name <- c("OPLAH")
known_loci_chr8 <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x1021/individual_cond_pruned_var.RDS"))
genofile <- seqOpen("../Data/chr8_DRAGEN_agds_b1_final.gds")
load("../Data/STAAR_exp_groups/x1021/x1021_nullmodel_batch1.Rdata")
chr <- 8
# plof_ds
burden_plof_ds_chr8 <- Gene_Centric_Coding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                             rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                             known_loci=known_loci_chr8,
                                                             QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                             Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                             Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                             category="plof_ds")
table(burden_plof_ds_chr8[,1])
burden_plof_ds_chr8 <- as.matrix(burden_plof_ds_chr8)
burden_plof_ds_chr8 <- as.data.frame(burden_plof_ds_chr8[,1])
colnames(burden_plof_ds_chr8) <- ("OPLAH_plof_ds")
saveRDS(burden_plof_ds_chr8, "../Data/STAAR_exp_groups/x1021/burden_plof_ds.RDS")

# disruptive missense
burden_dm_chr8 <- Gene_Centric_Coding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                               rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                               known_loci=known_loci_chr8,
                                                               QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                               Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                               Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                               category="disruptive_missense")
table(burden_dm_chr8[,1])
burden_dm_chr8 <- as.matrix(burden_dm_chr8)
burden_dm_chr8 <- as.data.frame(burden_dm_chr8[,1])
colnames(burden_dm_chr8) <- ("OPLAH_disruptive_missense")

# check correlations
corr_burden_chr8 <- cor(cbind(burden_plof_ds_chr8[,1], burden_dm_chr8[,1]))
colnames(corr_burden_chr8) <- c("OPLAH_plof_ds", "OPLAH_disruptive_missense")
rownames(corr_burden_chr8) <- c("OPLAH_plof_ds", "OPLAH_disruptive_missense")
Heatmap(corr_burden_chr8,
        name = "Correlation",
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        row_names_gp = gpar(fontsize = 8),       # row label size
        column_names_gp = gpar(fontsize = 8),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", corr_burden_chr8[i, j]), x, y, gp = gpar(fontsize = 8))})


cors # everything with a correlation coef > 0.91, select plof_ds burden in AM model

saveRDS(burden_dm_chr8, "../Data/STAAR_exp_groups/x1021/burden_disruptive_missense.RDS")


##### Batch 2 #####
genofile <- seqOpen("../Data/chr8_DRAGEN_agds_b2_final.gds")
load("../Data/STAAR_exp_groups/x1021/x1021_nullmodel_batch2.Rdata")
burden_plof_ds_chr8_b2 <- Gene_Centric_Coding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                               rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                               known_loci=known_loci_chr8,
                                                               QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                               Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                               Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                               category="plof_ds")
table(burden_plof_ds_chr8_b2[,1])
burden_plof_ds_chr8_b2 <- as.matrix(burden_plof_ds_chr8_b2)
burden_plof_ds_chr8_b2 <- as.data.frame(burden_plof_ds_chr8_b2[,1])
colnames(burden_plof_ds_chr8_b2) <- ("OPLAH_plof_ds")
saveRDS(burden_plof_ds_chr8_b2, "../Data/STAAR_exp_groups/x1021/b2/burden_plof_ds.RDS")


##### chr 11 #####
chr11_coding_sig <- readRDS("../Data/STAAR_exp_groups/x100009332/coding_cond_sig.RDS")
gene_name <- c("MS4A14")
known_loci_chr11 <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x100009332/individual_cond_pruned_var.RDS"))
genofile <- seqOpen("../Data/chr11_scale_up_DRAGEN_agds_b1_final_1.gds")
load("../Data/STAAR_exp_groups/x100009332/x100009332_nullmodel_batch1.Rdata")
chr <- 11
burden_chr_11_synonymous <- Gene_Centric_Coding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                               rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                               known_loci=known_loci_chr11,
                                                               QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                               Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                               Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                               category="synonymous")

table(burden_chr_11_synonymous[,1])
burden_chr_11_synonymous <- as.matrix(burden_chr_11_synonymous)
burden_chr_11_synonymous <- as.data.frame(burden_chr_11_synonymous[,1])
colnames(burden_chr_11_synonymous) <- "MS4A14_synonymous"
saveRDS(burden_chr_11_synonymous, "../Data/STAAR_exp_groups/x100009332/burden_synonymous.RDS")
##### chr 12 #####
chr12_1_methyl_coding_sig <- readRDS("../Data/STAAR_exp_groups/x100001208/coding_cond_sig.RDS")
gene_name <- c("SLC6A13")
known_loci_chr12_1_methyl <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x100001208/individual_cond_pruned_var.RDS"))
genofile <- seqOpen("../Data/chr12_scale_up_DRAGEN_agds_b1_final_1.gds")
load("../Data/STAAR_exp_groups/x100001208/x100001208_nullmodel_batch1.Rdata")
chr <- 12
burden_chr12_1_methyl_plof_ds <- Gene_Centric_Coding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                                    rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                                    known_loci=known_loci_chr12_1_methyl,
                                                                    QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                                    Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                                    Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                                    category="plof_ds")

table(burden_chr12_1_methyl_plof_ds[,1])
burden_chr12_1_methyl_plof_ds <- as.matrix(burden_chr12_1_methyl_plof_ds)
burden_chr12_1_methyl_plof_ds <- as.data.frame(burden_chr12_1_methyl_plof_ds[,1])

burden_chr12_1_methyl_dm <- Gene_Centric_Coding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                                         rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                                         known_loci=known_loci_chr12_1_methyl,
                                                                         QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                                         Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                                         category="disruptive_missense")

table(burden_chr12_1_methyl_dm[,1])
burden_chr12_1_methyl_dm <- as.matrix(burden_chr12_1_methyl_dm)
burden_chr12_1_methyl_dm <- as.data.frame(burden_chr12_1_methyl_dm[,1])


# super correlated - pearson r = 0.98
corr_burden_chr12_1_methyl <- cor(cbind(burden_chr12_1_methyl_plof_ds[,1], burden_chr12_1_methyl_dm[,1]))
colnames(corr_burden_chr12_1_methyl) <- c("SLC6A13_plof_ds", "SLC6A13_disruptive_missense")
rownames(corr_burden_chr12_1_methyl) <- c("SLC6A13_plof_ds", "SLC6A13_disruptive_missense")
Heatmap(corr_burden_chr12_1_methyl,
        name = "Correlation",
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        row_names_gp = gpar(fontsize = 8),       # row label size
        column_names_gp = gpar(fontsize = 8),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", corr_burden_chr12_1_methyl[i, j]), x, y, gp = gpar(fontsize = 8))})


# pick the plof_ds one
colnames(burden_chr12_1_methyl_plof_ds) <- "SLC6A13_plof_ds"
saveRDS(burden_chr12_1_methyl_plof_ds, "../Data/STAAR_exp_groups/x100001208/burden_plof_ds.RDS")
####### 12 ethylmaalonate #####
chr12_b1_ethylmaalonate_coding_sig <- readRDS("../Data/STAAR_exp_groups/x2054/coding_cond_sig.RDS")
gene_name <- c("ACADS")
known_loci_chr12_ethylmaalonate <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x2054/individual_cond_pruned_var.RDS"))
genofile <- seqOpen("../Data/chr12_scale_up_DRAGEN_agds_b1_final_1_ethylmalonate.gds")
load("../Data/STAAR_exp_groups/x2054/x2054_nullmodel_batch1.Rdata")
burden_chr12_1_ethylmaalonate <- Gene_Centric_Coding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                                        rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                                        known_loci=known_loci_chr12_ethylmaalonate,
                                                                        QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                                        category="plof_ds")
table(burden_chr12_1_ethylmaalonate[,1])
burden_chr12_1_ethylmaalonate <- as.matrix(burden_chr12_1_ethylmaalonate)
burden_chr12_1_ethylmaalonate <- as.data.frame(burden_chr12_1_ethylmaalonate[,1])
colnames(burden_chr12_1_ethylmaalonate) <- "ACADS_plof_ds"
saveRDS(burden_chr12_1_ethylmaalonate, "../Data/STAAR_exp_groups/x2054/burden_plof_ds.RDS")


burden_chr12_1_ethylmaalonate_dm <- Gene_Centric_Coding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                                         rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                                         known_loci=known_loci_chr12_ethylmaalonate,
                                                                         QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                                         Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                                         category="disruptive_missense")
table(burden_chr12_1_ethylmaalonate_dm[,1])
burden_chr12_1_ethylmaalonate_dm <- as.matrix(burden_chr12_1_ethylmaalonate_dm)
burden_chr12_1_ethylmaalonate_dm <- as.data.frame(burden_chr12_1_ethylmaalonate_dm[,1])
colnames(burden_chr12_1_ethylmaalonate_dm) <- "ACADS_disruptive_missense"
saveRDS(burden_chr12_1_ethylmaalonate_dm, "../Data/STAAR_exp_groups/x2054/burden_dm.RDS")
corr_burden_chr12_1_ethylmaalonate <- cor(cbind(burden_chr12_1_ethylmaalonate[,1], burden_chr12_1_ethylmaalonate_dm[,1]))
colnames(corr_burden_chr12_1_ethylmaalonate) <- c("ACADS_plof_ds", "ACADS_disruptive_missense")
rownames(corr_burden_chr12_1_ethylmaalonate) <- c("ACADS_plof_ds", "ACADS_disruptive_missense")
Heatmap(corr_burden_chr12_1_ethylmaalonate,
        name = "Correlation",
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        row_names_gp = gpar(fontsize = 8),       # row label size
        column_names_gp = gpar(fontsize = 8),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", corr_burden_chr12_1_ethylmaalonate[i, j]), x, y, gp = gpar(fontsize = 8))})


#### chr 12, batch 2 coding ####
chr12_b2_ethylmaalonate_coding_sig <- readRDS("../Data/STAAR_exp_groups/x2054/b2/coding_cond_sig.RDS")
genofile <- seqOpen("../Data/chr12_scale_up_DRAGEN_agds_b2_final_1_ethylmalonate.gds")
load("../Data/STAAR_exp_groups/x100001208/x100001208_nullmodel_batch2.Rdata")
burden_chr12_b2_ethylmaalonate <- Gene_Centric_Coding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                                         rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                                         known_loci=known_loci_chr12_ethylmaalonate,
                                                                         QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                                         Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                                         category="plof_ds")
table(burden_chr12_b2_ethylmaalonate[,1])
burden_chr12_b2_ethylmaalonate <- as.matrix(burden_chr12_b2_ethylmaalonate)
burden_chr12_b2_ethylmaalonate <- as.data.frame(burden_chr12_b2_ethylmaalonate[,1])
colnames(burden_chr12_b2_ethylmaalonate) <- "ACADS_plof_ds"
saveRDS(burden_chr12_b2_ethylmaalonate, "../Data/STAAR_exp_groups/x2054/b2/burden_plof_ds.RDS")
#### chromosome 16 -- Propyl 4-hydroxybenzoate sulfate####
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
source("./20250508_enhancer_dhs_exreact_burden.R")
source("./20250508_downstream_cond_extract_burden.R")
source("./20250429_utr_extract_burden.R")
source("./20250508_promoter_cage_extract_burden.R")
source("./20250508_promoter_dhs_extract_burden.R")

# summary_chr16_noncoding_cond <- as.data.frame(summary_chr16_noncoding_cond)
# rownames(summary_chr16_noncoding_cond) < NULL
# unique(summary_chr16_noncoding_cond$Category)

chr16_noncoding_sig <- readRDS("/Volumes/Sofer Lab/HCHS_SOL/Projects/2024_rare_variants/Data/STAAR_exp_groups/x100006264/noncoding_cond_sig.RDS")
known_loci_chr16 <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x100006264/individual_cond_pruned_var.RDS"))
genofile <- seqOpen("../Data/chr16_scale_up_DRAGEN_agds_b1_final_1.gds")
load("../Data/STAAR_exp_groups/x100006264/x100006264_nullmodel_batch1.Rdata")
chr <- 16

burdens_utr <- c()
gene_names <- c("IL4R", "MVP", "RNF40")
for(gene_name in gene_names){
  burden_test<- Gene_Centric_Noncoding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                           rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                           known_loci=known_loci_chr16,
                                                           QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                           Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                           Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                           category="UTR")
  burden_test <- as.matrix(burden_test)
  burdens_utr <- cbind(burdens_utr, burden_test[,1])
}

burdens_enhancer_cage <- c()
gene_names <- c("IL4R")
for(gene_name in gene_names){
  burden_test<- Gene_Centric_Noncoding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                           rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                           known_loci=known_loci_chr16,
                                                           QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                           Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                           Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                           category="enhancer_CAGE")
  burden_test <- as.matrix(burden_test)
  burdens_enhancer_cage <- cbind(burdens_enhancer_cage, burden_test[,1])
}

burdens_enhancer_dhs <- c()
gene_names <- c("KDM8", "IL4R", "PRRT2", "FAM57B", "ITGAL", "FUS")
for(gene_name in gene_names){
  burden_test<- Gene_Centric_Noncoding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                           rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                           known_loci=known_loci_chr16,
                                                           QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                           Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                           Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                           category="enhancer_DHS")
  burden_test <- as.matrix(burden_test)
  burdens_enhancer_dhs <- cbind(burdens_enhancer_dhs, burden_test[,1])
}


burdens_promoter_cage <- c()
for(gene_name in gene_names){
  burden_test<- Gene_Centric_Noncoding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                           rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                           known_loci=known_loci_chr16,
                                                           QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                           Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                           Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                           category="promoter_CAGE")
  burden_test <- as.matrix(burden_test)
  burdens_promoter_cage <- cbind(burdens_promoter_cage, burden_test[,1])
}


burdens_promoter_dhs <- c()
gene_names <- c("IL4R","NUPR1","ITGAL")
for(gene_name in gene_names){
  burden_test<- Gene_Centric_Noncoding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                           rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                           known_loci=known_loci_chr16,
                                                           QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                           Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                           Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                           category="promoter_DHS")
  burden_test <- as.matrix(burden_test)
  burdens_promoter_dhs <- cbind(burdens_promoter_dhs, burden_test[,1])
}


combined_burden <- cbind(burdens_utr, burdens_promoter_cage,
                         burdens_promoter_dhs,
                         burdens_enhancer_cage, burdens_enhancer_dhs)

gene_names_all <- unlist(chr16_noncoding_sig$`Gene name`)
names(gene_names_all) <- NULL
func_name_all <- unlist(chr16_noncoding_sig$Category)
names(func_name_all) <- NULL
colnames(combined_burden) <- paste0(gene_names_all, "_", func_name_all)

saveRDS(combined_burden, "../Data/STAAR_exp_groups/x100006264/burdens_combined.RDS")

###### chr 16 batch 2 Propyl 4-hydroxybenzoate sulfate ######
chr16_noncoding_sig_b2 <- get(load("/Volumes/Sofer Lab/HCHS_SOL/Projects/2024_rare_variants/Data/STAAR_exp_groups/x100006264/b2/noncoding_cond_sig.Rdata"))
chr16_noncoding_sig_b2 <- toDataframe(chr16_noncoding_sig_b2)
genofile <- seqOpen("../Data/chr16_scale_up_DRAGEN_agds_b2_final_1.gds")
load("../Data/STAAR_exp_groups/x100006264/x100006264_nullmodel_batch2.Rdata")
burden_test_b2_chr16_tot <- c()
gene_name <- "PRRT2"
burden_test_PRRT2 <- Gene_Centric_Noncoding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                         rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                         known_loci=known_loci_chr16,
                                                         QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                         Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                         category="enhancer_DHS")
burden_test_PRRT2 <- as.matrix(burden_test_PRRT2)
burden_test_b2_chr16_tot <- cbind(burden_test_b2_chr16_tot, burden_test_PRRT2[,1])
colnames(burden_test_b2_chr16_tot) <- "PRRT2_enhancer_DHS"


gene_name <- "FAM57B"
burden_test_FAM57B <- Gene_Centric_Noncoding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                                rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                                known_loci=known_loci_chr16,
                                                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                                Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                                category="enhancer_DHS")
burden_test_FAM57B <- as.matrix(burden_test_FAM57B)
burden_test_b2_chr16_tot <- cbind(burden_test_b2_chr16_tot, burden_test_FAM57B[,1])
colnames(burden_test_b2_chr16_tot)[2] <- "FAM57B_enhancer_DHS"


gene_name <- "RNF40"
burden_test_RNF40 <- Gene_Centric_Noncoding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                                rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                                known_loci=known_loci_chr16,
                                                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                                Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                                category="UTR")

burden_test_RNF40 <- as.matrix(burden_test_RNF40)
burden_test_b2_chr16_tot <- cbind(burden_test_b2_chr16_tot, burden_test_RNF40[,1])
colnames(burden_test_b2_chr16_tot)[3] <- "RNF40_UTR"

saveRDS(burden_test_b2_chr16_tot, "../Data/STAAR_exp_groups/x100006264/b2/burdens_combined.RDS")


######## chr16 -- Cys-gly, oxidized ########
# chr16_b1_cys_gly <- readRDS("../Data/STAAR_exp_groups/x1224/coding_cond_sig.RDS")
gene_name <- c("DPEP1")
known_loci_chr16_cys_gly <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x1224/individual_cond_pruned_var.RDS"))
genofile <- seqOpen("../Data/chr16_DRAGEN_agds_b1_final.gds")
load("../Data/STAAR_exp_groups/x1224/x1224_nullmodel_batch1.Rdata")
source("./20250425_plof_extract_burden.R")
burden_chr16_cys_gly_plof <- Gene_Centric_Coding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                                         rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                                         known_loci=known_loci_chr16_cys_gly,
                                                                         QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                                         Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                                         category="plof")
table(burden_chr16_cys_gly_plof[,1])
burden_chr16_cys_gly_plof <- as.matrix(burden_chr16_cys_gly_plof)
burden_chr16_cys_gly_plof <- as.data.frame(burden_chr16_cys_gly_plof[,1])
colnames(burden_chr16_cys_gly_plof) <- "DPEP1_plof"
saveRDS(burden_chr16_cys_gly_plof, "../Data/STAAR_exp_groups/x1224/burden_plof.RDS")


burden_chr16_cys_gly_plof_ds <- Gene_Centric_Coding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                                     rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                                     known_loci=known_loci_chr16_cys_gly,
                                                                     QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                                     Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                                     Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                                     category="plof_ds")
table(burden_chr16_cys_gly_plof_ds[,1])
burden_chr16_cys_gly_plof_ds <- as.matrix(burden_chr16_cys_gly_plof_ds)
burden_chr16_cys_gly_plof_ds <- as.data.frame(burden_chr16_cys_gly_plof_ds[,1])
colnames(burden_chr16_cys_gly_plof_ds) <- "DPEP1_plof_ds"
saveRDS(burden_chr16_cys_gly_plof_ds,
        "../Data/STAAR_exp_groups/x1224/burden_plof_ds.RDS")


burden_chr16_cys_gly_disruptive_missense <- Gene_Centric_Coding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                                        rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                                        known_loci=known_loci_chr16_cys_gly,
                                                                        QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                                        category="disruptive_missense")
table(burden_chr16_cys_gly_disruptive_missense[,1])
burden_chr16_cys_gly_disruptive_missense <- as.matrix(burden_chr16_cys_gly_disruptive_missense)
burden_chr16_cys_gly_disruptive_missense <- as.data.frame(burden_chr16_cys_gly_disruptive_missense[,1])
colnames(burden_chr16_cys_gly_disruptive_missense) <- "DPEP1_disruptive_missense"
saveRDS(burden_chr16_cys_gly_disruptive_missense,
        "../Data/STAAR_exp_groups/x1224/burden_chr16_cys_gly_disruptive_missense.RDS")

corr_burden_chr16_cys_gly <- cor(cbind(burden_chr16_cys_gly_plof_ds[,1],
                                       burden_chr16_cys_gly_disruptive_missense[,1],
                                                burden_chr16_cys_gly_plof[,1]),
                                method = "spearman")
colnames(corr_burden_chr16_cys_gly) <- c("DPEP1_plof_ds", "DPEP1_disruptive_missense", "DPEP1_plof")
rownames(corr_burden_chr16_cys_gly) <- c("DPEP1_plof_ds", "DPEP1_disruptive_missense", "DPEP1_plof")
Heatmap(corr_burden_chr16_cys_gly,
        name = "Correlation",
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        row_names_gp = gpar(fontsize = 8),       # row label size
        column_names_gp = gpar(fontsize = 8),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", corr_burden_chr16_cys_gly[i, j]), x, y, gp = gpar(fontsize = 8))})


#------------------------------------------------------------------------------
# LEGACY / EARLIER VERSION: v1_original  (109 unique blocks)
#------------------------------------------------------------------------------
# Kept verbatim. These are blocks that do NOT appear in the current version above
# (mostly hard-coded per-metabolite / per-chromosome run calls and older path setups).

# adapted from STAARpipeline source code to extract annotation ranks
Anno.Int.PHRED.sub <- NULL
Anno.Int.PHRED.sub.name <- NULL
genofile <- seqOpen(agds_dir)
for(k in 1:length(Annotation_name))
{
  if(Annotation_name[k]%in%Annotation_name_catalog$name)
  {
    Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
    Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))

    if(Annotation_name[k]=="CADD")
    {
      Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
    }

    if(Annotation_name[k]=="aPC.LocalDiversity")
    {
      Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
      Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
      Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
    }
    Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
  }
}

Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name


#' STAAR procedure for conditional analysis using omnibus test
#'
#' The \code{STAAR_cond} function takes in genotype, the genotype of variants to be
#' adjusted for in conditional analysis, the object from fitting the null
#' model, and functional annotation data to analyze the conditional association between a
#' quantitative/dichotomous phenotype and a variant-set by using STAAR procedure,
#' adjusting for a given list of variants. For each variant-set, the conditional
#' STAAR-O p-value is a p-value from an omnibus test that aggregated conditional
#' SKAT(1,25), SKAT(1,1), Burden(1,25), Burden(1,1), ACAT-V(1,25), and ACAT-V(1,1)
#' together with conditional p-values of each test weighted by each annotation
#' using Cauchy method.
#' @param genotype an n*p genotype matrix (dosage matrix) of the target sequence,
#' where n is the sample size and p is the number of genetic variants.
#' @param genotype_adj an n*p_adj genotype matrix (dosage matrix) of the target
#' sequence, where n is the sample size and p_adj is the number of genetic variants
#' to be adjusted for in conditional analysis (or a vector of a single variant with length n
#' if p_adj is 1).
#' @param obj_nullmodel an object from fitting the null model, which is the
#' output from either \code{\link{fit_null_glm}} function for unrelated samples or
#' \code{\link{fit_null_glmmkin}} function for related samples. Note that \code{\link{fit_null_glmmkin}}
#' is a wrapper of the \code{\link{glmmkin}} function from the \code{\link{GMMAT}} package.
#' @param annotation_phred a data frame or matrix of functional annotation data
#' of dimension p*q (or a vector of a single annotation score with length p).
#' Continuous scores should be given in PHRED score scale, where the PHRED score
#' of j-th variant is defined to be -10*log10(rank(-score_j)/total) across the genome. (Binary)
#' categorical scores should be taking values 0 or 1, where 1 is functional and 0 is
#' non-functional. If not provided, STAAR will perform the
#' SKAT(1,25), SKAT(1,1), Burden(1,25), Burden(1,1), ACAT-V(1,25), ACAT-V(1,1)
#' and ACAT-O tests (default = NULL).
#' @param rare_maf_cutoff the cutoff of maximum minor allele frequency in
#' defining rare variants (default = 0.01).
#' @param rv_num_cutoff the cutoff of minimum number of variants of analyzing
#' a given variant-set (default = 2).
#' @param rv_num_cutoff_max the cutoff of maximum number of variants of analyzing
#' a given variant-set (default = 1e+09).
#' @param method_cond a character value indicating the method for conditional analysis.
#' \code{optimal} refers to regressing residuals from the null model on \code{genotype_adj}
#' as well as all covariates used in fitting the null model (fully adjusted) and taking the residuals;
#' \code{naive} refers to regressing residuals from the null model on \code{genotype_adj}
#' and taking the residuals (default = \code{optimal}).
#' @return A list with the following members:
#' @return \code{num_variant}: the number of variants with minor allele frequency > 0 and less than
#' \code{rare_maf_cutoff} in the given variant-set that are used for performing the
#' variant-set using STAAR.
#' @return \code{cMAC}: the cumulative minor allele count of variants with
#' minor allele frequency > 0 and less than \code{rare_maf_cutoff} in the given variant-set.
#' @return \code{RV_label}: the boolean vector indicating whether each variant in the given
#' variant-set has minor allele frequency > 0 and less than \code{rare_maf_cutoff}.
#' @return \code{results_STAAR_O_cond}: the conditional STAAR-O p-value that aggregated conditional
#' SKAT(1,25), SKAT(1,1), Burden(1,25), Burden(1,1), ACAT-V(1,25), and ACAT-V(1,1) together
#' with conditional p-values of each test weighted by each annotation using Cauchy method.
#' @return \code{results_ACAT_O_cond}: the conditional ACAT-O p-value that aggregated conditional
#' SKAT(1,25), SKAT(1,1), Burden(1,25), Burden(1,1), ACAT-V(1,25), and ACAT-V(1,1) using Cauchy method.
#' @return \code{results_STAAR_S_1_25_cond}: a vector of conditional STAAR-S(1,25) p-values,
#' including conditional SKAT(1,25) p-value weighted by MAF, the conditional SKAT(1,25)
#' p-values weighted by each annotation, and a conditional STAAR-S(1,25)
#' p-value by aggregating these p-values using Cauchy method.
#' @return \code{results_STAAR_S_1_1_cond}: a vector of conditional STAAR-S(1,1) p-values,
#' including conditional SKAT(1,1) p-value weighted by MAF, the conditional SKAT(1,1)
#' p-values weighted by each annotation, and a conditional STAAR-S(1,1)
#' @return \code{results_STAAR_B_1_25_cond}: a vector of conditional STAAR-B(1,25) p-values,
#' including conditional Burden(1,25) p-value weighted by MAF, the conditional Burden(1,25)
#' p-values weighted by each annotation, and a conditional STAAR-B(1,25)
#' @return \code{results_STAAR_B_1_1_cond}: a vector of conditional STAAR-B(1,1) p-values,
#' including conditional Burden(1,1) p-value weighted by MAF, the conditional Burden(1,1)
#' p-values weighted by each annotation, and a conditional STAAR-B(1,1)
#' @return \code{results_STAAR_A_1_25_cond}: a vector of conditional STAAR-A(1,25) p-values,
#' including conditional ACAT-V(1,25) p-value weighted by MAF, the conditional ACAT-V(1,25)
#' p-values weighted by each annotation, and a conditional STAAR-A(1,25)
#' @return \code{results_STAAR_A_1_1_cond}: a vector of conditional STAAR-A(1,1) p-values,
#' including conditional ACAT-V(1,1) p-value weighted by MAF, the conditional ACAT-V(1,1)
#' p-values weighted by each annotation, and a conditional STAAR-A(1,1)
#' @references Li, X., Li, Z., et al. (2020). Dynamic incorporation of multiple
#' in silico functional annotations empowers rare variant association analysis of
#' large whole-genome sequencing studies at scale. \emph{Nature Genetics}, \emph{52}(9), 969-983.
#' (\href{https://doi.org/10.1038/s41588-020-0676-4}{pub})
#' @references Li, Z., Li, X., et al. (2022). A framework for detecting
#' noncoding rare-variant associations of large-scale whole-genome sequencing
#' studies. \emph{Nature Methods}, \emph{19}(12), 1599-1611.
#' (\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
#' @references Liu, Y., et al. (2019). Acat: A fast and powerful p value combination
#' method for rare-variant analysis in sequencing studies.
#' \emph{The American Journal of Human Genetics}, \emph{104}(3), 410-421.
#' (\href{https://doi.org/10.1016/j.ajhg.2019.01.002}{pub})
#' @references Li, Z., Li, X., et al. (2020). Dynamic scan procedure for
#' detecting rare-variant association regions in whole-genome sequencing studies.
#' \emph{The American Journal of Human Genetics}, \emph{104}(5), 802-814.
#' (\href{https://doi.org/10.1016/j.ajhg.2019.03.002}{pub})
#' @references Sofer, T., et al. (2019). A fully adjusted two-stage procedure for rank-normalization
#' in genetic association studies. \emph{Genetic Epidemiology}, \emph{43}(3), 263-275.
#' (\href{https://doi.org/10.1002/gepi.22188}{pub})
#' @export

STAAR_cond <- function(genotype,genotype_adj,obj_nullmodel,annotation_phred=NULL,
                       rare_maf_cutoff=0.01,rv_num_cutoff=2,rv_num_cutoff_max=1e9,
                       method_cond=c("optimal","naive")){
  print("running STAAR_cond customized function")
  method_cond <- match.arg(method_cond) # evaluate choices
  if(!inherits(genotype, "matrix") && !inherits(genotype, "Matrix")){
    stop("genotype is not a matrix!")
  }

  if(dim(genotype)[2] == 1){
    stop(paste0("Number of rare variant in the set is less than 2!"))
  }

  annotation_phred <- as.data.frame(annotation_phred)
  if(dim(annotation_phred)[1] != 0 & dim(genotype)[2] != dim(annotation_phred)[1]){
    stop(paste0("Dimensions don't match for genotype and annotation!"))
  }

  if(inherits(genotype, "sparseMatrix")){
    genotype <- as.matrix(genotype)
  }

  if(inherits(genotype_adj, "numeric")){
    genotype_adj <- matrix(genotype_adj, ncol=1)
  }

  if(dim(genotype)[1] != dim(genotype_adj)[1]){
    stop(paste0("Dimensions don't match for genotype and genotype_adj!"))
  }
  genotype <- matrix_flip(genotype)
  MAF <- genotype$MAF
  RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0))
  Geno_rare <- genotype$Geno[,RV_label]

  rm(genotype)
  gc()
  annotation_phred <- annotation_phred[RV_label,,drop=FALSE]

  if(sum(RV_label) >= rv_num_cutoff_max){
    stop(paste0("Number of rare variant in the set is more than ",rv_num_cutoff_max,"!"))
  }

  if(sum(RV_label) >= rv_num_cutoff){
    G <- as(Geno_rare,"dgCMatrix")
    MAF <- MAF[RV_label]
    rm(Geno_rare)
    gc()

    annotation_rank <- 1 - 10^(-annotation_phred/10)

    ## beta(1,25)
    w_1 <- dbeta(MAF,1,25)
    ## beta(1,1)
    w_2 <- dbeta(MAF,1,1)
    if(dim(annotation_phred)[2] == 0){
      ## Burden, SKAT, ACAT-V
      w_B <- w_S <- as.matrix(cbind(w_1,w_2))
      w_A <- as.matrix(cbind(w_1^2/dbeta(MAF,0.5,0.5)^2,w_2^2/dbeta(MAF,0.5,0.5)^2))
    }else{
      ## Burden
      w_B_1 <- annotation_rank*w_1
      w_B_1 <- cbind(w_1,w_B_1)
      colnames(w_B_1) <- paste0("Burden(1, 25)-", colnames(w_B_1))
      colnames(w_B_1)[1] <- "Burden(1, 25)"

      w_B_2 <- annotation_rank*w_2
      w_B_2 <- cbind(w_2,w_B_2)
      colnames(w_B_2) <- paste0("Burden(1, 1)-", colnames(w_B_2))
      colnames(w_B_2)[1] <- "Burden(1, 1)"

      w_B <- cbind(w_B_1,w_B_2)
      w_B <- as.matrix(w_B)
      print(colnames(w_B))
      print(dim(w_B))

      ## SKAT
      w_S_1 <- sqrt(annotation_rank)*w_1
      w_S_1 <- cbind(w_1,w_S_1)
      w_S_2 <- sqrt(annotation_rank)*w_2
      w_S_2 <- cbind(w_2,w_S_2)
      w_S <- cbind(w_S_1,w_S_2)
      w_S <- as.matrix(w_S)

      ## ACAT-V
      w_A_1 <- annotation_rank*w_1^2/dbeta(MAF,0.5,0.5)^2
      w_A_1 <- cbind(w_1^2/dbeta(MAF,0.5,0.5)^2,w_A_1)
      w_A_2 <- annotation_rank*w_2^2/dbeta(MAF,0.5,0.5)^2
      w_A_2 <- cbind(w_2^2/dbeta(MAF,0.5,0.5)^2,w_A_2)
      w_A <- cbind(w_A_1,w_A_2)
      w_A <- as.matrix(w_A)
    }

    if(obj_nullmodel$relatedness){
      if(!obj_nullmodel$sparse_kins){
        P <- obj_nullmodel$P

        residuals.phenotype <- obj_nullmodel$scaled.residuals
        if(method_cond == "optimal"){
          residuals.phenotype.fit <- lm(residuals.phenotype~genotype_adj+obj_nullmodel$X-1)
        }else{
          residuals.phenotype.fit <- lm(residuals.phenotype~genotype_adj)
        }
        residuals.phenotype <- residuals.phenotype.fit$residuals
        X_adj <- model.matrix(residuals.phenotype.fit)
        PX_adj <- P%*%X_adj
        P_cond <- P - X_adj%*%solve(t(X_adj)%*%X_adj)%*%t(PX_adj) -
          PX_adj%*%solve(t(X_adj)%*%X_adj)%*%t(X_adj) +
          X_adj%*%solve(t(X_adj)%*%X_adj)%*%t(PX_adj)%*%X_adj%*%solve(t(X_adj)%*%X_adj)%*%t(X_adj)
        rm(P)
        gc()

        pvalues <- STAAR_O_SMMAT(G,P_cond,residuals.phenotype,
                                 weights_B=w_B,weights_S=w_S,weights_A=w_A,
                                 mac=as.integer(round(MAF*2*dim(G)[1])))
      }else{
        Sigma_i <- obj_nullmodel$Sigma_i
        Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
        cov <- obj_nullmodel$cov

        residuals.phenotype <- obj_nullmodel$scaled.residuals
        if(method_cond == "optimal"){
          residuals.phenotype.fit <- lm(residuals.phenotype~genotype_adj+obj_nullmodel$X-1)
        }else{
          residuals.phenotype.fit <- lm(residuals.phenotype~genotype_adj)
        }
        residuals.phenotype <- residuals.phenotype.fit$residuals
        X_adj <- model.matrix(residuals.phenotype.fit)

        pvalues <- STAAR_O_SMMAT_sparse_cond(G,Sigma_i,Sigma_iX,cov,X_adj,residuals.phenotype,
                                             weights_B=w_B,weights_S=w_S,weights_A=w_A,
                                             mac=as.integer(round(MAF*2*dim(G)[1])))
      }

    }else{
      X <- model.matrix(obj_nullmodel)
      working <- obj_nullmodel$weights
      sigma <- sqrt(summary(obj_nullmodel)$dispersion)
      if(obj_nullmodel$family[1] == "binomial"){
        P <- diag(working) - X%*%solve(t(X)%*%diag(working)%*%X)%*%t(X)
      }else if(obj_nullmodel$family[1] == "gaussian"){
        P <- diag(length(working)) - X%*%solve(t(X)%*%X)%*%t(X)
      }

      residuals.phenotype <- obj_nullmodel$y - obj_nullmodel$fitted.values
      if(method_cond == "optimal"){
        residuals.phenotype.fit <- lm(residuals.phenotype~genotype_adj+model.matrix(obj_nullmodel)-1)
      }else{
        residuals.phenotype.fit <- lm(residuals.phenotype~genotype_adj)
      }
      residuals.phenotype <- residuals.phenotype.fit$residuals
      X_adj <- model.matrix(residuals.phenotype.fit)
      PX_adj <- P%*%X_adj
      P_cond <- P - X_adj%*%solve(t(X_adj)%*%X_adj)%*%t(PX_adj) -
        PX_adj%*%solve(t(X_adj)%*%X_adj)%*%t(X_adj) +
        X_adj%*%solve(t(X_adj)%*%X_adj)%*%t(PX_adj)%*%X_adj%*%solve(t(X_adj)%*%X_adj)%*%t(X_adj)
      rm(P)
      gc()

      pvalues <- STAAR_O_SMMAT(G,P_cond,residuals.phenotype,
                               weights_B=w_B,weights_S=w_S,weights_A=w_A,
                               mac=as.integer(round(MAF*2*dim(G)[1])))
    }

    # Assuming residuals and G are Armadillo vectors/matrices
    print(dim(G))
    print(dim(w_B))
    burden <- G %*% w_B


    # n <- length(x)  # Equivalent to x.size() in C++
    wn <- ncol(w_B)  # Equivalent to weights_B.n_cols in C++
    #
    # sum_vec <- numeric(wn)  # Pre-allocate a vector to store sums
    #
    # for (i in 1:wn) {
    #   sum0 <- sum(x * w_B[, i])  # Vectorized multiplication and sum
    #   sum_vec[i] <- sum0  # Store result
    # }
    print(paste0("wn = ", wn))
    # names(sum_vec) <- colnames(w_B)
    saveRDS(burden, "../Data/chr16_coding_burden.RDS")


    num_variant <- sum(RV_label) #dim(G)[2]
    cMAC <- sum(G)
    num_annotation <- dim(annotation_phred)[2]+1
    results_STAAR_O <- CCT(pvalues)
    results_ACAT_O <- CCT(pvalues[c(1,num_annotation+1,2*num_annotation+1,3*num_annotation+1,4*num_annotation+1,5*num_annotation+1)])
    pvalues_STAAR_S_1_25 <- CCT(pvalues[1:num_annotation])
    pvalues_STAAR_S_1_1 <- CCT(pvalues[(num_annotation+1):(2*num_annotation)])
    pvalues_STAAR_B_1_25 <- CCT(pvalues[(2*num_annotation+1):(3*num_annotation)])
    pvalues_STAAR_B_1_1 <- CCT(pvalues[(3*num_annotation+1):(4*num_annotation)])
    pvalues_STAAR_A_1_25 <- CCT(pvalues[(4*num_annotation+1):(5*num_annotation)])
    pvalues_STAAR_A_1_1 <- CCT(pvalues[(5*num_annotation+1):(6*num_annotation)])

    results_STAAR_S_1_25 <- c(pvalues[1:num_annotation],pvalues_STAAR_S_1_25)
    results_STAAR_S_1_25 <- data.frame(t(results_STAAR_S_1_25))

    results_STAAR_S_1_1 <- c(pvalues[(num_annotation+1):(2*num_annotation)],pvalues_STAAR_S_1_1)
    results_STAAR_S_1_1 <- data.frame(t(results_STAAR_S_1_1))

    results_STAAR_B_1_25 <- c(pvalues[(2*num_annotation+1):(3*num_annotation)],pvalues_STAAR_B_1_25)
    results_STAAR_B_1_25 <- data.frame(t(results_STAAR_B_1_25))

    results_STAAR_B_1_1 <- c(pvalues[(3*num_annotation+1):(4*num_annotation)],pvalues_STAAR_B_1_1)
    results_STAAR_B_1_1 <- data.frame(t(results_STAAR_B_1_1))

    results_STAAR_A_1_25 <- c(pvalues[(4*num_annotation+1):(5*num_annotation)],pvalues_STAAR_A_1_25)
    results_STAAR_A_1_25 <- data.frame(t(results_STAAR_A_1_25))

    results_STAAR_A_1_1 <- c(pvalues[(5*num_annotation+1):(6*num_annotation)],pvalues_STAAR_A_1_1)
    results_STAAR_A_1_1 <- data.frame(t(results_STAAR_A_1_1))

    if(dim(annotation_phred)[2] == 0){
      colnames(results_STAAR_S_1_25) <- c("SKAT(1,25)","STAAR-S(1,25)")
      colnames(results_STAAR_S_1_1) <- c("SKAT(1,1)","STAAR-S(1,1)")
      colnames(results_STAAR_B_1_25) <- c("Burden(1,25)","STAAR-B(1,25)")
      colnames(results_STAAR_B_1_1) <- c("Burden(1,1)","STAAR-B(1,1)")
      colnames(results_STAAR_A_1_25) <- c("ACAT-V(1,25)","STAAR-A(1,25)")
      colnames(results_STAAR_A_1_1) <- c("ACAT-V(1,1)","STAAR-A(1,1)")
    }else{
      colnames(results_STAAR_S_1_25) <- c("SKAT(1,25)",
                                          paste0("SKAT(1,25)-",colnames(annotation_phred)),
                                          "STAAR-S(1,25)")
      colnames(results_STAAR_S_1_1) <- c("SKAT(1,1)",
                                         paste0("SKAT(1,1)-",colnames(annotation_phred)),
                                         "STAAR-S(1,1)")
      colnames(results_STAAR_B_1_25) <- c("Burden(1,25)",
                                          paste0("Burden(1,25)-",colnames(annotation_phred)),
                                          "STAAR-B(1,25)")
      colnames(results_STAAR_B_1_1) <- c("Burden(1,1)",
                                         paste0("Burden(1,1)-",colnames(annotation_phred)),
                                         "STAAR-B(1,1)")
      colnames(results_STAAR_A_1_25) <- c("ACAT-V(1,25)",
                                          paste0("ACAT-V(1,25)-",colnames(annotation_phred)),
                                          "STAAR-A(1,25)")
      colnames(results_STAAR_A_1_1) <- c("ACAT-V(1,1)",
                                         paste0("ACAT-V(1,1)-",colnames(annotation_phred)),
                                         "STAAR-A(1,1)")
    }

    return(list(num_variant = num_variant,
                cMAC = cMAC,
                RV_label = RV_label,
                results_STAAR_O_cond = results_STAAR_O,
                results_ACAT_O_cond = results_ACAT_O,
                results_STAAR_S_1_25_cond = results_STAAR_S_1_25,
                results_STAAR_S_1_1_cond = results_STAAR_S_1_1,
                results_STAAR_B_1_25_cond = results_STAAR_B_1_25,
                results_STAAR_B_1_1_cond = results_STAAR_B_1_1,
                results_STAAR_A_1_25_cond = results_STAAR_A_1_25,
                results_STAAR_A_1_1_cond = results_STAAR_A_1_1))
  }else{
    stop(paste0("Number of rare variant in the set is less than ",rv_num_cutoff,"!"))
  }

}

gene_name <- "ANKRD11"
results <- Gene_Centric_Coding(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                    rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                    QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                    Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                    Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                    category="plof_ds")

results_cond <- Gene_Centric_Coding_cond(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                    rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                    known_loci=known_loci,
                                    QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                    Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                    Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                    category="plof_ds")

# subsetting only the burden p-values(including hte omnibus ones)
results_cond_burden <- results_cond[,34:61]
saveRDS(results_cond_burden, "../Data/chr16_coding_burden_pvals.RDS")

# open the burden file
test_burden <- as.matrix(readRDS("/Volumes/Sofer Lab/HCHS_SOL/Projects/2024_rare_variants/Data/chr16_coding_burden.RDS"))
