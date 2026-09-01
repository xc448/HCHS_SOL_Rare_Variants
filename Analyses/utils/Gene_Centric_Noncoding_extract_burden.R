# import pacakges
library(tidyverse)
library(SeqArray)
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)
library(SeqArray)
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

##### chr 6 #####
gene_name <- c("UNC93A")
known_loci_chr6 <- as.data.frame(readRDS("../Data/STAAR_neg_controls/x1215/individual_cond_pruned_var.RDS"))
genofile <- seqOpen("../Data/chr6_DRAGEN_agds_b1_final.gds")
load("../Data/STAAR_neg_controls/x1215/x1215_nullmodel_batch1.Rdata")
chr <- 6
# conditional analysis for the four unique rv set
# plof

burden_utr <- Gene_Centric_Noncoding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                       rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                       known_loci=known_loci_chr6,
                                                       QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                       Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                       Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                                       category="UTR") 

burden_utr <- as.matrix(burden_utr)
seqClose(genofile)
saveRDS(burden_utr, "../Data/STAAR_neg_controls/x1215/burden_UTR.RDS")

var_anno_chr6_utr <- Gene_Centric_Noncoding_Info(category=c("UTR"),
                                                    chr = 6,genofile,obj_nullmodel, gene_name = c("UNC93A"),
                                                    known_loci=known_loci_chr6,rare_maf_cutoff=0.01,
                                                    method_cond=method_cond,
                                                    QC_label=QC_label,variant_type=variant_type,
                                                    geno_missing_imputation=geno_missing_imputation,
                                                    Annotation_dir=Annotation_dir,
                                                    Annotation_name_catalog = Annotation_name_catalog,
                                                    Annotation_name = Annotation_name)
saveRDS(var_anno_chr6_utr, "../Data/STAAR_neg_controls/x1215/RVsetInfo_noncoding_sig_utr.RDS")

##### chr 8 ncRNA #####
gene_name <- c("RNU6-220P")
gene_name <- c("PLEC")
known_loci_chr8 <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x1021/individual_cond_pruned_var.RDS"))
genofile <- seqOpen("../Data/chr8_DRAGEN_agds_b1_final.gds")
load("../Data/STAAR_exp_groups/x1021/x1021_nullmodel_batch1.Rdata")
chr <- 8

burden_enhancer_dhs_plec <- Gene_Centric_Noncoding_cond_extract_burden(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                         rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                         known_loci=known_loci_chr8, category=c("enhancer_DHS"),
                                                         QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                         Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name) 

burden_enhancer_dhs_plec <- as.matrix(burden_enhancer_dhs_plec)
saveRDS(burden_enhancer_dhs_plec, "../Data/STAAR_exp_groups/x1021/burden_enhancer_dhs_plec.RDS")

burden_ncrna <- ncRNA_cond(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                                                       rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                                                       known_loci=known_loci_chr8[-c(1,2),],
                                                                       QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                                                       Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                                                       Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name) 

burden_ncrna <- as.matrix(burden_ncrna)
saveRDS(burden_ncrna, "../Data/STAAR_exp_groups/x1021/burden_ncRNA.RDS")

var_anno_chr8_ncrna <- Gene_Centric_Noncoding_Info(category=c("ncRNA"),
                                                 chr = 8,genofile,obj_nullmodel, gene_name = c("RNU6-220P"),
                                                 rare_maf_cutoff=0.01, known_loci=known_loci_chr8[-c(1,2)],
                                                 method_cond=method_cond,
                                                 QC_label=QC_label,variant_type=variant_type,
                                                 geno_missing_imputation=geno_missing_imputation,
                                                 Annotation_dir=Annotation_dir,
                                                 Annotation_name_catalog = Annotation_name_catalog,
                                                 Annotation_name = Annotation_name)
saveRDS(var_anno_chr8_ncrna, "../Data/STAAR_exp_groups/x1021/RVsetInfo_noncoding_sig_ncRNA.RDS")
seqClose(genofile)

