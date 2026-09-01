# import pacakges
library(tidyverse)
library(SeqArray)
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)
library(SeqArray)
library(SeqVarTools)

## QC_label
QC_label <- "annotation/filter"
## variant_type
variant_type <- "SNV"
## geno_missing_imputation
geno_missing_imputation <- "mean"

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load("../Data/STAAR_prep/Annotation_name_catalog.Rdata"))
# Or equivalently
# Annotation_name_catalog <- read.csv("/path_to_the_file/Annotation_name_catalog.csv")
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")


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
  print(annotation_phred)
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
    return(c(burden, annotation_phred))
    
  }else{
    stop(paste0("Number of rare variant in the set is less than ",rv_num_cutoff,"!"))
  }
}



##### chr 12 #####
chr12_coding_sig <- readRDS("/Volumes/Sofer Lab/HCHS_SOL/Projects/2024_rare_variants/Data/STAAR_exp_groups/x2054/coding_sig.RDS")
gene_name <- c("ACADS")
known_loci_chr12 <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x2054/individual_cond_pruned_var.RDS"))
genofile <- seqOpen("../Data/chr12_scale_up_DRAGEN_agds_b1_final_1_ethylmalonate.gds")
load("../Data/STAAR_exp_groups/x2054/x2054_nullmodel_batch1.Rdata")
chr <- 12


source("./20250425_plof_ds_cond_modified.R")
source("./20250425_plof_ds_extract_burden.R")
burden_plof_ds <- Gene_Centric_Coding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                           rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                           known_loci=known_loci_chr12,
                                                           QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                           Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                           Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                           category="plof_ds") 

burden_plof_ds_test <- as.matrix(burden_plof_ds[[1]])
seqClose(genofile)
saveRDS(burden_plof_ds_test, "../Data/STAAR_exp_groups/x2054/burden_plof_ds.RDS")
table(burden_plof_ds_test[,1])

# source("./20250425_missense_extract_burden.R")
# burden_missense <- Gene_Centric_Coding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
#                                                           rare_maf_cutoff=0.01,rv_num_cutoff=2,
#                                                           known_loci=known_loci_chr12,
#                                                           QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
#                                                           Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
#                                                           Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
#                                                           category="missense") 


