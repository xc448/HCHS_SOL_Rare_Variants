#==============================================================================
# STAAR pipeline - BATCH 1  (prep, null model, gene-centric & individual analysis)
#==============================================================================
# MERGED FROM (chronological):
#   - 20241229_STAAR_pipeline.R   ->  section "v1_single_chr16"
#   - 20250414_STAAR_multiple_chr.R   ->  section "v2_multi_chr"
#   - 20250710_STAAR_multiple_scaledup.R   ->  section "v3_scaled_up"
# NOTE: Phenotype: covariate_b1_all_metab*.RDS ; kinship: HCHS_SOL_kinship_matrix_b1.RDS
# NOTE: v1 hard-codes chr16 only; v2 uses chr_info <- c(2,5,6,8,10,13,16); v3 globs GDS via list.files().
# NOTE: Set output_path to ../Data/STAAR_prep/ (v1,v2) or ../Data/STAAR_prep/scale_up/ (v3).
#==============================================================================

#==============================================================================
# MAIN PIPELINE  (current version: v3_scaled_up)
#==============================================================================

###########################################################
# Pre-step for running STAARpipeline

setwd("/Volumes/Sofer Lab/HCHS_SOL/Projects/2024_rare_variants/Code")
setwd("R:\\Sofer Lab\\HCHS_SOL\\Projects/2024_rare_variants/Code")

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(tidyverse)
library(STAARpipeline)
library(STAAR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#           User Input
## file directory of aGDS file (genotype and annotation data)
dir.geno <- "../Data/"
## file name of aGDS, seperate by chr number
agds_file_name_1 <- "cardiac_cohorts_SOL_dragen.chr"
agds_file_name_2 <- "_scaled_up_admixmap_sigregion.gds"
## channel name of the QC label in the GDS/aGDS file
QC_label <- "annotation/filter"
## file directory for the output files
#output_path <- "../Data/STAAR_prep/"
output_path <- "../Data/STAAR_prep/scale_up/"
## annotation name. The first eight names are used to define masks in gene-centric analysis, do not change them!!
## The others are the annotation you want to use in the STAAR procedure, and they are flexible to change.
name <- c("rs_num","GENCODE.Category","GENCODE.Info","GENCODE.EXONIC.Category",
          "MetaSVM","GeneHancer","CAGE","DHS","CADD","LINSIGHT","FATHMM.XF",
          "aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
          "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")
## channel name of the annotations. Make sure they are matched with the name, especially for the first eight one!!
dir <- c("/rsid","/genecode_comprehensive_category","/genecode_comprehensive_info",
         "/genecode_comprehensive_exonic_category","/metasvm_pred",
         "/genehancer","/cage_tc","/rdhs","/cadd_phred","/linsight","/fathmm_xf",
         "/apc_epigenetics_active","/apc_epigenetics_repressed","/apc_epigenetics_transcription",
         "/apc_conservation","/apc_local_nucleotide_diversity","/apc_mappability",
         "/apc_transcription_factor","/apc_protein_function")

#           Main Function
## aGDS directory
gds_files_b1 <- list.files(
  path = "../Data",
  pattern = "^cardiac_cohorts_SOL_dragen\\.chr[0-9]+_batch1+_scaled_up_admixmap_sigregion.*\\.gds$",
  full.names = TRUE
)

gds_files_b2 <- list.files(
  path = "../Data",
  pattern = "^cardiac_cohorts_SOL_dragen\\.chr[0-9]+_batch2+_scaled_up_admixmap_sigregion.*\\.gds$",
  full.names = TRUE
)

agds_dir <- c(gds_files_b1, gds_files_b2)

save(agds_dir,file=paste0("../Data/STAAR_prep/agds_dir_scaled_up.Rdata",sep=""))
load("../Data/STAAR_prep/agds_dir_scaled_up.Rdata")


## Annotation name catalog (alternatively, can skip this part by providing Annotation_name_catalog.csv with the same information)
Annotation_name_catalog <- data.frame(name=name,dir=dir)
save(Annotation_name_catalog,file=paste0(output_path,"Annotation_name_catalog.Rdata",sep=""))
load(paste0(output_path,"Annotation_name_catalog.Rdata"))

## Number of jobs for each chromosome
chr_info <- c(11, 11, 12, 12, 16, 2, 5, 8)
jobs_num <- matrix(rep(0,3*length(chr_info)),nrow=length(chr_info))
for(i in 1:length(chr_info))
{
  chr <- chr_info[i]
  print(chr)
  gds.path <- agds_dir[1:length(chr_info)][i]
  genofile <- seqOpen(gds.path)

  filter <- seqGetData(genofile, QC_label)
  SNVlist <- filter == "PASS"

  position <- as.numeric(seqGetData(genofile, "position"))

  jobs_num[i,1] <- chr
  jobs_num[i,2] <- min(position[SNVlist])
  jobs_num[i,3] <- max(position[SNVlist])

  seqClose(genofile)
}

## Individual Analysis
jobs_num <- cbind(jobs_num,ceiling((jobs_num[,3]-jobs_num[,2])/10e6))
## Sliding Window Analysis
jobs_num <- cbind(jobs_num,ceiling((jobs_num[,3]-jobs_num[,2])/5e6))
## Dynamic Window Analysis (SCANG-STAAR)
jobs_num <- cbind(jobs_num,ceiling((jobs_num[,3]-jobs_num[,2])/1.5e6))

jobs_num <- as.data.frame(jobs_num)
colnames(jobs_num) <- c("chr","start_loc","end_loc","individual_analysis_num",
                        "sliding_window_num","scang_num")

save(jobs_num,file=paste0(output_path,"mult_chr_jobs_num_scale_up.Rdata",sep=""))
load(paste0(output_path,"mult_chr_jobs_num_scale_up.Rdata"))

# add necessary field to the info section in the annotated GDS file
for(i in agds_dir){
  agds <- seqOpen(i, readonly = FALSE)

  # we needto add AVGDP (average sequencing depth information,
  # here we will just set everyone to 10 because it's the default for FastSparseGRM
  Anno.folder <- index.gdsn(agds, "annotation/info")
  position <- read.gdsn(index.gdsn(agds, "position"))
  add.gdsn(Anno.folder, "AVGDP", val=rep(10, length(position)),
           compress="LZMA_ra", closezip=TRUE)
  seqClose(agds)

}


# fit STAAR null model for single-trait analysis
# (using either GENESIS or glmm)
library(GENESIS)
# ## Phenotype file
phenotype <- readRDS("../Data/covariate_b1_all_metab_scaled_up.RDS")
phenotype <- phenotype[phenotype$NWD_ID  != "", ]
rownames(phenotype) <- phenotype$NWD_ID
## kinship matrix file
load("../Data/kinsip_pcrelate_Matrix_no_outlies_round2.RData")
dim(km)

met <- rep("", 8)
met[1] <- "_3_beta"
met[3] <- "_ethylmalonate"
for(i in 3:3){
  print(agds_dir[i])
  agds <- seqOpen(agds_dir[i])
  sample.id <- read.gdsn(index.gdsn(agds, "sample.id"))
  sample.id <- as.data.frame(sample.id)

  mat_b1_ppl <- as.data.frame(rownames(kinmat_b1))

  intersect_ppl <- sample.id |>
    inner_join(mat_b1_ppl, by = c("sample.id" = "rownames(kinmat_b1)"))

  seqSetFilter(agds, sample.id =intersect_ppl[,1] )
  output_path <- paste0("../Data/chr", chr_info[i], "_scale_up_DRAGEN_agds_b1_final_1",
                         met[i], ".gds")
  seqExport(agds, out.fn = output_path)
  seqClose(agds)
}

# for(i in 1:length(chr_info)){
#   agds <- seqOpen(agds_dir[8:14][i])
#   sample.id <- read.gdsn(index.gdsn(agds, "sample.id"))
#   sample.id <- as.data.frame(sample.id)
#
#   mat_b2_ppl <- as.data.frame(rownames(km))
#   intersect_ppl <- sample.id |>
#     inner_join(mat_b2_ppl, by = c("sample.id" = "rownames(km)"))
#   seqSetFilter(agds, sample.id =intersect_ppl[,1] )
#   output_path <- paste0("../Data/chr", chr_info[i], "_DRAGEN_agds_b2_final.gds")
#   seqExport(agds,  out.fn = output_path)
#   seqClose(agds)
# }
# phenotype <- phenotype[intersect_ppl[[1]],]
# kinmat_b1 <- km[intersect_ppl[,1], intersect_ppl[,1]] # 3889 x 3889
# saveRDS(kinmat_b1, "../Data/HCHS_SOL_kinship_matrix_b1.RDS")
kinmat_b1 <- readRDS("../Data/HCHS_SOL_kinship_matrix_b1.RDS")
phenotype <- phenotype[rownames(kinmat_b1),]

## fit null model using GENESIS
data_GENESIS <- as(phenotype, "AnnotatedDataFrame") # Make AnnotatedDataFrame (specifically required by GENESIS)

## household matrix
# batch1
load("../Data/hh.matrix_b1.Rdata")
hh.matrix_b1 <- hh.matrix_b1[rownames(phenotype), rownames(phenotype)]

## block matrix
# load block matrices
load("../Data/block.matrix_b1.Rdata")
block.matrix_b1 <- block.matrix_b1[rownames(phenotype), rownames(phenotype)]

# load genetic PC files

# adjust 5 PCs as indicated in the file??
load("../Data/pca_no_outliers_round2.Rdata")
pcs <- as.data.frame(pca$vectors)
n <- ncol(pcs)
names(pcs) <- paste0("PC", 1:n)
pcs$sample.id <- row.names(pcs)

pcs_b1 <- pcs[rownames(phenotype),] #sum(startsWith(rownames(pcs_b1), "NWD"))

#pcs_b2 <- pcs[rownames(phenotype_b2),] #sum(startsWith(rownames(pcs_b2), "NWD"))


b1_all_cov_for_STAAR_null <- cbind(phenotype, pcs_b1[,1:5])

# b2_all_cov_for_STAAR_null <- cbind(phenotype, pcs_b2[,1:5])
# saveRDS(b2_all_cov_for_STAAR_null, "../Data/STAAR_model_cov_batch2.RDS")

covMatList <- list(HH = hh.matrix_b1, kinship = as.matrix(kinmat_b1), block = block.matrix_b1)
b1_all_cov_for_STAAR_null$BKGRD1_C7 <- as.factor(b1_all_cov_for_STAAR_null$BKGRD1_C7)
b1_all_cov_for_STAAR_null$GENDER <- as.factor(b1_all_cov_for_STAAR_null$GENDER)
b1_all_cov_for_STAAR_null$CENTER <- as.factor(b1_all_cov_for_STAAR_null$CENTER)
# saveRDS(b1_all_cov_for_STAAR_null, "../Data/STAAR_model_cov_batch1_scaled_up.RDS")
saveRDS(covMatList, "../Data/STAAR_model_covmat_batch1.RDS")

b1_all_cov_for_STAAR_null <- readRDS("../Data/STAAR_model_cov_batch1_scaled_up.RDS")
covMatList <- readRDS("../Data/STAAR_model_covmat_batch1.RDS")

# Fitting nullmodel

# neg controls
neg_control_ids <- paste0("x", c(100001721, 799, 1215, 278))
exp_group_ids <- paste0("x", c(1021, 100000007, 100004046, 1224))

b1_all_cov_for_STAAR_null <- b1_all_cov_for_STAAR_null |>
  mutate(BKGRD1_C7 = if_else(BKGRD1_C7 == 0, "Dominican", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == 1, "Central American", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == 2, "Cuban", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == 3, "Mexican", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == 4, "Puerto Rican", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == 5, "South American", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == 6, "Other", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == "Q", "Other", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == "", "Other", BKGRD1_C7),) |>
  mutate(GFRSCYS = if_else(is.na(GFRSCYS), mean(GFRSCYS, na.rm = TRUE), GFRSCYS))

saveRDS(b1_all_cov_for_STAAR_null, "../Data/STAAR_model_cov_batch1_scaled_up.RDS")

outcomes <- paste0("x", c(100001266, 1114, 192, 100006370,
              100009332, 100001208, 2054, 100006264))


for(i in outcomes){
  print(i)
  obj_nullmodel_GENESIS <- fitNullModel(x = b1_all_cov_for_STAAR_null,
                                        outcome=i,
                                        covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                 "PC4", "PC5",
                                                 "BKGRD1_C7", "GFRSCYS",
                                                 "CENTER"),
                                        cov.mat=covMatList,
                                        AIREML.tol=1e-4,
                                        verbose=TRUE)

  saveRDS(obj_nullmodel_GENESIS, paste0("../Data/STAAR_exp_groups/",
                                        i, "/", i, "_nullmodel_GENESIS_batch1.RDS"))

  # obj_nullmodel_GENESIS <- readRDS(paste0("../Data/STAAR_exp_groups/",
  #                                         i, "/", i, "_nullmodel_GENESIS_batch1.RDS"))

  ## convert GENESIS null model to STAAR null model
  obj_nullmodel <- genesis2staar_nullmodel(obj_nullmodel_GENESIS)
  obj_nullmodel$Sigma_i <- as(obj_nullmodel$Sigma_i, "sparseMatrix")

  # save nullmodel
  save(obj_nullmodel,file= paste0("../Data/STAAR_exp_groups/",
                                  i, "/", i, "_nullmodel_batch1.Rdata"))
}


jobs_num <- get(load("../Data/STAAR_prep/scale_up/mult_chr_jobs_num_scale_up.Rdata"))
files <- list.files(path = "../Data/", pattern = "_scale_up_DRAGEN_agds_")
agds_dir <- paste0("../Data/", files)
save(agds_dir, file = "../Data/STAAR_prep/agds_dir_scale_up_final.Rdata")
agds_dir <- get(load("../Data/STAAR_prep/agds_dir_scale_up_final.Rdata"))
## QC_label
## variant_type
variant_type <- "variant"
## geno_missing_imputation
geno_missing_imputation <- "mean"

## output file name
output_file_name <- "Individual_Analysis"

## aGDS file # batch 1
STAAR_ind_analysis <- function(agds.path, group, metab, i, chr){
  if(group == "exp_groups"){
    load(paste0("../Data/STAAR_exp_groups/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  } else {
    load(paste0("../Data/STAAR_neg_controls/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }
  genofile <- seqOpen(agds.path)
  job_num <- jobs_num[i,]
  start_loc <- job_num$start_loc
  end_loc <- job_num$end_loc
  obj_nullmodel$sparse_kins <- TRUE
  a <- Sys.time()
  results_individual_analysis <- c()
  print(start_loc)
  print()
  if(start_loc <= end_loc)
  {
    results_individual_analysis <- Individual_Analysis(chr=chr,start_loc=start_loc,end_loc=end_loc,
                                                       genofile=genofile,obj_nullmodel=obj_nullmodel,mac_cutoff=20,
                                                       QC_label=QC_label,variant_type=variant_type,
                                                       geno_missing_imputation=geno_missing_imputation)
  }
  b <- Sys.time()
  b - a
  save(results_individual_analysis,file=paste0(output_path, metab, "/",metab, "_",
                                               output_file_name,".Rdata"))

  seqClose(genofile)
  return(results_individual_analysis)
}


## output path - experimental groups
output_path <- "../Data/STAAR_exp_groups/"

results_individual_chr11_uncond_arachido <- STAAR_ind_analysis(agds_dir[1], i = 1, chr = 11,
                                                               metab = "x100009332",
                                                               group = "exp_groups")

results_individual_chr11_uncond_3beta <- STAAR_ind_analysis(agds_dir[2], i = 2, chr=11,
                                                               metab = "x100006370",
                                                               group = "exp_groups")

results_individual_chr12_uncond_1_methyl <- STAAR_ind_analysis(agds_dir[3], i = 3, chr = 12,
                                                               metab = "x100001208",
                                                               group = "exp_groups")


results_individual_chr12_uncond_ethylmalonate <- STAAR_ind_analysis(agds_dir[4], i = 4, chr = 12,
                                                                       metab = "x2054",
                                                                       group = "exp_groups")

results_individual_chr16 <- STAAR_ind_analysis(agds_dir[5], i = 5, chr = 16,
                                                  metab = "x100006264",
                                                  group = "exp_groups")

results_individual_chr2_arginine <- STAAR_ind_analysis(agds_dir[6], i = 6, chr = 2,
                                                          metab = "x100001266",
                                                          group = "exp_groups")


results_individual_chr5_3amino <- STAAR_ind_analysis(agds_dir[7], i = 7, chr = 5,
                                                        metab = "x1114",
                                                        group = "exp_groups")


results_individual_chr8_putrescine <- STAAR_ind_analysis(agds_dir[8], i = 8, chr = 8,
                                                            metab = "x192",
                                                            group = "exp_groups")

# do gene-centric coding analysis for each of the gene
# within the LAI
variant_type <- "SNV"
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


output_file_name <- "Genecentric_coding_uncond"
geneCentric_coding_uncond <- function(agds.path, chr, group, metab){
  if(group == "exp_groups"){
    ## output path
    output_path <- "../Data/STAAR_exp_groups/"
    load(paste0("../Data/STAAR_exp_groups/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }else{
    output_path <- "../Data/STAAR_neg_controls/"
    load(paste0("../Data/STAAR_neg_controls/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }

  genofile <- seqOpen(agds.path)
  position <- read.gdsn(index.gdsn(genofile, "position"))
  genes_info_chr <- genes_info[genes_info[,2]==chr,]
  genes_info_chr <- as.data.frame(genes_info_chr) |>
    filter(end_position >= min(position)) |>
    filter(start_position <= max(position))
  sub_seq_num <- dim(genes_info_chr)[1]

  results_coding <- c()
  categories=c("plof","plof_ds","missense","disruptive_missense","synonymous","ptv","ptv_ds")
  for(kk in 1:nrow(genes_info_chr))
  {
    gene_name <- genes_info_chr[kk,1]
    print(gene_name)
    results <- Gene_Centric_Coding(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                   rare_maf_cutoff=0.01,rv_num_cutoff=2, category = "all_categories",
                                   QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                   Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                   Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)


    results_coding <- append(results_coding,results)

  }

  save(results_coding,file=paste0(output_path, metab, "/",
                                  metab, "_", output_file_name,".Rdata"))

  seqClose(genofile)

  return(results_coding)
}


results_coding_chr11_uncond_arachido <- geneCentric_coding_uncond(agds_dir[1], chr = 11,
                                                                  metab = "x100009332",
                                                                  group = "exp_groups")

results_coding_chr11_uncond_3beta <- geneCentric_coding_uncond(agds_dir[2], chr = 11,
                                                        metab = "x100006370",
                                                        group = "exp_groups")


results_coding_chr12_uncond_1_methyl <- geneCentric_coding_uncond(agds_dir[3], chr = 12,
                                                                       metab = "x100001208",
                                                                       group = "exp_groups")


results_coding_chr12_uncond_ethylmalonate <- geneCentric_coding_uncond(agds_dir[4], chr = 12,
                                                                       metab = "x2054",
                                                                       group = "exp_groups")


results_coding_chr16 <- geneCentric_coding_uncond(agds_dir[5], chr = 16,
                                                  metab = "x100006264",
                                                  group = "exp_groups")

results_coding_chr2_arginine <- geneCentric_coding_uncond(agds_dir[6], chr = 2,
                                                          metab = "x100001266",
                                                          group = "exp_groups")


results_coding_chr5_3amino <- geneCentric_coding_uncond(agds_dir[7], chr = 5,
                                                          metab = "x1114",
                                                          group = "exp_groups")


results_coding_chr8_putrescine <- geneCentric_coding_uncond(agds_dir[8], chr = 8,
                                                        metab = "x192",
                                                        group = "exp_groups")

source("./20250206_Gene_Centric_Coding_cond_customized.R") # you have to call this
source("./20250425_missense_cond_modified.R")
source("./20250425_plof_cond_modified.R")
source("./20250425_plof_ds_cond_modified.R")
source("./20250425_STAAR_cond_customized.R")
output_file_name <- "Genecentric_coding_cond"
geneCentric_coding_cond <- function(agds.path, chr, group, metab, known_loci){
  if(group == "exp_groups"){
    ## output path
    output_path <- "../Data/STAAR_exp_groups/"
    load(paste0("../Data/STAAR_exp_groups/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }else{
    output_path <- "../Data/STAAR_neg_controls/"
    load(paste0("../Data/STAAR_neg_controls/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }
  genofile <- seqOpen(agds.path)
  position <- read.gdsn(index.gdsn(genofile, "position"))
  genes_info_chr <- genes_info[genes_info[,2]==chr,]
  genes_info_chr <- as.data.frame(genes_info_chr) |>
    filter(end_position >= min(position)) |>
    filter(start_position <= max(position))
  sub_seq_num <- dim(genes_info_chr)[1]

  genes <- genes_info

  categories=c("plof","plof_ds","missense","disruptive_missense","synonymous","ptv","ptv_ds")
  for(category in categories){
    results_coding <- c()
    print(category)
    for(kk in 1:nrow(genes_info_chr))
    {
      gene_name <- genes_info_chr[kk,1]
      print(gene_name)

      results <- Gene_Centric_Coding_cond(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                          rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                          known_loci=known_loci, category = category,
                                          QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                          Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                          Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
      results_coding <- rbind(results_coding,results)
    }
    save(results_coding,file=paste0(output_path, metab, "/",
                                    metab, "_", output_file_name, "_",
                                    category,
                                    ".Rdata"))
  }
  seqClose(genofile)
  return(results_coding)
}


known_loci_chr11_arachido <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x100009332/individual_cond_pruned_var.RDS"))
results_coding_chr11_arachido <- geneCentric_coding_cond(agds_dir[1], chr = 11, metab = "x100009332", group = "exp_groups",
                                               known_loci = known_loci_chr11_arachido)

known_loci_chr12_1_methyl <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x100001208/individual_cond_pruned_var.RDS"))
results_coding_chr12_1_methyl <- geneCentric_coding_cond(agds_dir[3], chr = 12, metab = "x100001208", group = "exp_groups",
                                                         known_loci = known_loci_chr12_1_methyl)

known_loci_chr12_ethylmalonate <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x2054/individual_cond_pruned_var.RDS"))
results_coding_chr12_ethylmalonate <- geneCentric_coding_cond(agds_dir[4], chr = 12, metab = "x2054", group = "exp_groups",
                                                         known_loci = known_loci_chr12_ethylmalonate)

known_loci_chr2_arginine <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x100001266/individual_cond_pruned_var.RDS"))
results_coding_chr2_arginine <- geneCentric_coding_cond(agds_dir[6], chr = 2, metab = "x100001266", group = "exp_groups",
                                                         known_loci = known_loci_chr2_arginine)

known_loci_chr5_3amino <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x1114/individual_cond_pruned_var.RDS"))
results_coding_chr5_3amino <- geneCentric_coding_cond(agds_dir[7], chr = 5, metab = "x1114", group = "exp_groups",
                                                        known_loci = known_loci_chr5_3amino)


known_loci_chr16 <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x100006264/individual_cond_pruned_var.RDS"))
results_coding_chr16 <- geneCentric_coding_cond(agds_dir[5], chr = 16, metab = "x100006264", group = "exp_groups",
                                                      known_loci = known_loci_chr16)

## output path
output_path <- "../Data/"
output_file_name <- "Genecentric_noncoding_uncond"
## input array id from batch file
# arrayid <- as.numeric(commandArgs(TRUE)[1])

geneCentric_noncoding_uncond <- function(agds.path, chr, group, metab){
  if(group == "exp_groups"){
    ## output path
    output_path <- "../Data/STAAR_exp_groups/"
    load(paste0("../Data/STAAR_exp_groups/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }else{
    output_path <- "../Data/STAAR_neg_controls/"
    load(paste0("../Data/STAAR_neg_controls/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }
  genofile <- seqOpen(agds.path)
  position <- read.gdsn(index.gdsn(genofile, "position"))
  genes_info_chr <- genes_info[genes_info[,2]==chr,]
  genes_info_chr <- as.data.frame(genes_info_chr) |>
    filter(end_position >= min(position)) |>
    filter(start_position <= max(position))
  sub_seq_num <- dim(genes_info_chr)[1]


  results_noncoding <- c()
  for(kk in 1:nrow(genes_info_chr))
  {
    print(kk)
    gene_name <- genes_info_chr[kk,1]
    print(gene_name)
    results <- Gene_Centric_Noncoding(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                      rare_maf_cutoff=0.01,rv_num_cutoff=2, category = "all_categories",
                                      QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                      Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                      Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
    results_noncoding <- append(results_noncoding,results)

  }

  save(results_noncoding,file=paste0(output_path, metab, "/",
                                     metab, "_", output_file_name,".Rdata"))

  seqClose(genofile)
  return(results_noncoding)
}


results_noncoding_chr11_uncond_3beta <- geneCentric_noncoding_uncond(agds_dir[1], chr = 11,
                                                               metab = "x100006370",
                                                               group = "exp_groups")

results_noncoding_chr11_uncond_arachido <- geneCentric_noncoding_uncond(agds_dir[2], chr = 11,
                                                                  metab = "x100009332",
                                                                  group = "exp_groups")

results_noncoding_chr12_uncond_ethylmalonate <- geneCentric_noncoding_uncond(agds_dir[4], chr = 12,
                                                                       metab = "x2054",
                                                                       group = "exp_groups")

results_noncoding_chr12_uncond_1_methyl <- geneCentric_noncoding_uncond(agds_dir[3], chr = 12,
                                                                  metab = "x100001208",
                                                                  group = "exp_groups")

results_noncoding_chr16 <- geneCentric_noncoding_uncond(agds_dir[5], chr = 16,
                                                  metab = "x100006264",
                                                  group = "exp_groups")


results_noncoding_uncond_chr2_arginine <- geneCentric_noncoding_uncond(agds_dir[6], chr = 2,
                                                          metab = "x100001266",
                                                          group = "exp_groups")

results_noncoding_uncond_chr5_3amino <- geneCentric_noncoding_uncond(agds_dir[7], chr = 5,
                                                        metab = "x1114",
                                                        group = "exp_groups")

results_noncoding_chr8_putrescine <- geneCentric_noncoding_uncond(agds_dir[8], chr = 8,
                                                            metab = "x192",
                                                            group = "exp_groups")


################## NONCODING CONDITONAL ANALYSIS #######################
source("20250425_STAAR_cond_customized.R")
source("./20250425_genecentric_Noncoding_modified.R")
source("./20250425_upstream_cond_modified.R")
source("./20250425_downstream_cond_modified.R")
source("./20250425_promoter_CAGE_cond_modified.R")
source("./20250425_promoter_DHS_cond_modified.R")
source("./20250425_enhancer_CAGE_cond_modified.R")
source("./20250425_enhancer_DHS_cond_modified.R")
source("./20250425_UTR_cond_modified.R")
output_file_name <- "Genecentric_noncoding_cond"
geneCentric_noncoding_cond <- function(agds.path, chr, group, metab, known_loci){
  if(group == "exp_groups"){
    ## output path
    output_path <- "../Data/STAAR_exp_groups/"
    load(paste0("../Data/STAAR_exp_groups/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }else{
    output_path <- "../Data/STAAR_neg_controls/"
    load(paste0("../Data/STAAR_neg_controls/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }
  genofile <- seqOpen(agds.path)
  position <- read.gdsn(index.gdsn(genofile, "position"))
  genes_info_chr <- genes_info[genes_info[,2]==chr,]
  genes_info_chr <- as.data.frame(genes_info_chr) |>
    filter(end_position >= min(position)) |>
    filter(start_position <= max(position))
  sub_seq_num <- dim(genes_info_chr)[1]

  genes <- genes_info

  categories=c("downstream","upstream","UTR","promoter_CAGE",
               "promoter_DHS","enhancer_CAGE","enhancer_DHS")
  for(category in categories){
    results_noncoding <- c()
    print(category)
    for(kk in 1:nrow(genes_info_chr))
    {
      gene_name <- genes_info_chr[kk,1]
      print(gene_name)

      results <- Gene_Centric_Noncoding_cond(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                             rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                             known_loci=known_loci, category = category,
                                             QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                             Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                             Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
      results_noncoding <- rbind(results_noncoding,results)
      print(results)
    }
    save(results_noncoding,file=paste0(output_path, metab, "/",
                                       metab, "_", output_file_name, "_",
                                       category,
                                       ".Rdata"))
  }
  seqClose(genofile)

}

known_loci_chr2 <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x100001266/individual_cond_pruned_var.RDS"))
results_noncoding_chr2 <- geneCentric_noncoding_cond(agds_dir[6], chr = 2, metab = "x100001266", group = "exp_groups",
                                                     known_loci = known_loci_chr2)

known_loci_chr5 <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x1114//individual_cond_pruned_var.RDS"))
results_noncoding_chr5 <- geneCentric_noncoding_cond(agds_dir[7], chr = 5, metab = "x1114", group = "exp_groups",
                                                     known_loci = known_loci_chr5)

known_loci_chr16_propyl <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x100006264/individual_cond_pruned_var.RDS"))
results_noncoding_chr16_propyl <- geneCentric_noncoding_cond(agds_dir[5], chr = 16, metab = "x100006264", group = "exp_groups",
                                                     known_loci = known_loci_chr16_propyl)

# run gene-centric non-coding analysis for ncRNA
output_file_name <- "Genetic_Noncoding_ncRNA"
# open the gds file
geneCentric_ncRNA_noncoding <- function(agds.path, chr, group, metab){
  if(group == "exp_groups"){
    ## output path
    output_path <- paste0("../Data/STAAR_exp_groups/", metab, "/")
    load(paste0("../Data/STAAR_exp_groups/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }else{
    output_path <- paste0("../Data/STAAR_neg_controls/", metab, "/")
    load(paste0("../Data/STAAR_neg_controls/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }
  genofile <- seqOpen(agds.path)
  position <- read.gdsn(index.gdsn(genofile, "position"))
  ncRNA_gene_chr <- ncRNA_gene[ncRNA_gene[,1]==chr,]

  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[;]")
  Gene <- as.character(sapply(GENCODE.Info.split,function(z) gsub("\\(.*\\)","",z[1])))

  Gene_list_1 <- as.character(sapply(strsplit(Gene,','),'[',1))
  Gene_list_2 <- as.character(sapply(strsplit(Gene,','),'[',2))
  Gene_list_3 <- as.character(sapply(strsplit(Gene,','),'[',3))

  gene_names <- ncRNA_gene_chr[which(ncRNA_gene_chr[,2] %in% Gene_list_1),]
  gene_names <- rbind(gene_names, ncRNA_gene_chr[which(ncRNA_gene_chr[,2] %in% Gene_list_2),])
  gene_names <- rbind(gene_names, ncRNA_gene_chr[which(ncRNA_gene_chr[,2] %in% Gene_list_3),])
  gene_names <- gene_names[!is.na(gene_names$ncRNA) & !gene_names$ncRNA == "NA" ,]
  ncRNA_gene_chr <- unique(gene_names)
  print(table(ncRNA_gene_chr))
  if(dim(ncRNA_gene_chr)[1] == 0){
    results_ncRNA_cond <- NULL
    save(results_ncRNA_cond,file=paste0(output_path,output_file_name,".Rdata"))
    return(NULL)
  }

  results_ncRNA <- c()
  for(kk in 1:nrow(ncRNA_gene_chr))
  {
    print(kk)
    gene_name <- ncRNA_gene_chr[kk,2]
    results <- c()
    results <- try(ncRNA(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                         rare_maf_cutoff=0.01,rv_num_cutoff=2,
                         QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                         Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name))
    results_ncRNA <- rbind(results_ncRNA,results)
  }
  seqClose(genofile)
  save(results_ncRNA,file=paste0(output_path,output_file_name,".Rdata"))
  return(results_ncRNA)
}


results_noncoding_ncrna_chr11_uncond_arachido <- geneCentric_ncRNA_noncoding(agds_dir[1], chr = 11,
                                                                        metab = "x100009332",
                                                                        group = "exp_groups")

results_noncoding_ncrna_chr11_uncond_3beta <- geneCentric_ncRNA_noncoding(agds_dir[2], chr = 11,
                                                                          metab = "x100006370",
                                                                          group = "exp_groups")

results_noncoding_ncrna_chr12_uncond_ethylmalonate <- geneCentric_ncRNA_noncoding(agds_dir[4], chr = 12,
                                                                             metab = "x2054",
                                                                             group = "exp_groups")

results_noncoding_ncrna_chr12_uncond_1_methyl <- geneCentric_ncRNA_noncoding(agds_dir[3], chr = 12,
                                                                        metab = "x100001208",
                                                                        group = "exp_groups")

results_noncoding_ncrna_chr16 <- geneCentric_ncRNA_noncoding(agds_dir[5], chr = 16,
                                                     metab = "x100006264",
                                                     group = "exp_groups")


results_noncoding_ncrna_uncond_chr2_arginine <- geneCentric_ncRNA_noncoding(agds_dir[6], chr = 2,
                                                                       metab = "x100001266",
                                                                       group = "exp_groups")

results_noncoding_ncrna_uncond_chr5_3amino <- geneCentric_ncRNA_noncoding(agds_dir[7], chr = 5,
                                                                     metab = "x1114",
                                                                     group = "exp_groups")

results_noncoding_ncrna_chr8_putrescine <- geneCentric_ncRNA_noncoding(agds_dir[8], chr = 8,
                                                                  metab = "x192",
                                                                  group = "exp_groups")

################## CONDITIONAL ANALYSIS - NCRNA ####################
output_file_name <- "Genetic_Noncoding_ncRNA_cond"
geneCentric_ncRNA_cond <- function(agds.path, chr, group, metab, known_loci){
  if(group == "exp_groups"){
    ## output path
    output_path <- paste0("../Data/STAAR_exp_groups/", metab, "/")
    load(paste0("../Data/STAAR_exp_groups/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }else{
    output_path <- paste0("../Data/STAAR_neg_controls/", metab, "/")
    load(paste0("../Data/STAAR_neg_controls/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }
  ncRNA_gene_chr <- ncRNA_gene[ncRNA_gene$chr == chr,]
  genofile <- seqOpen(agds.path)
  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[;]")
  Gene <- as.character(sapply(GENCODE.Info.split,function(z) gsub("\\(.*\\)","",z[1])))

  Gene_list_1 <- as.character(sapply(strsplit(Gene,','),'[',1))
  Gene_list_2 <- as.character(sapply(strsplit(Gene,','),'[',2))
  Gene_list_3 <- as.character(sapply(strsplit(Gene,','),'[',3))

  gene_names <- ncRNA_gene_chr[which(ncRNA_gene_chr[,2] %in% Gene_list_1),]
  gene_names <- rbind(gene_names, ncRNA_gene_chr[which(ncRNA_gene_chr[,2] %in% Gene_list_2),])
  gene_names <- rbind(gene_names, ncRNA_gene_chr[which(ncRNA_gene_chr[,2] %in% Gene_list_3),])
  gene_names <- gene_names[!is.na(gene_names$ncRNA) & !gene_names$ncRNA == "NA" ,]
  ncRNA_gene_chr <- unique(gene_names)
  print(table(ncRNA_gene_chr))
  if(dim(ncRNA_gene_chr)[1] == 0){
    results_ncRNA_cond <- NULL
    save(results_ncRNA_cond,file=paste0(output_path,output_file_name,".Rdata"))
    return(NULL)
  }
  results_ncRNA_cond <- c()
  for(kk in 1:nrow(ncRNA_gene_chr))
  {
    gene_name <- ncRNA_gene_chr[kk,2]
    print(gene_name)
    results <- c()
    results <- try(ncRNA_cond(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                              rare_maf_cutoff=0.01,rv_num_cutoff=2, known_loci = known_loci,
                              QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                              Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                              Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name))
    results_ncRNA_cond <- rbind(results_ncRNA_cond,results)
  }
  seqClose(genofile)
  save(results_ncRNA_cond,file=paste0(output_path,output_file_name,".Rdata"))
  return(results_ncRNA_cond)
}


results_noncoding_ncrna_cond_chr11_arachido <- geneCentric_ncRNA_cond(agds_dir[1], chr = 11,
                                                                             metab = "x100009332",
                                                                             group = "exp_groups",
                                                                      known_loci = known_loci_chr11_arachido)

known_loci_chr11_3beta <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x100006370/individual_cond_pruned_var.RDS"))
results_noncoding_ncrna_cond_chr11_3beta <- geneCentric_ncRNA_cond(agds_dir[2], chr = 11,
                                                                          metab = "x100006370",
                                                                          group = "exp_groups",
                                                                   known_loci = known_loci_chr11_3beta)

results_noncoding_ncrna_cond_chr12_ethylmalonate <- geneCentric_ncRNA_cond(agds_dir[4], chr = 12,
                                                                                  metab = "x2054",
                                                                                  group = "exp_groups",
                                                                           known_loci = known_loci_chr12_ethylmalonate)

results_noncoding_ncrna_cond_chr12_1_methyl <- geneCentric_ncRNA_cond(agds_dir[3], chr = 12,
                                                                             metab = "x100001208",
                                                                             group = "exp_groups",
                                                                             known_loci = known_loci_chr12_1_methyls)


results_noncoding_ncrna_cond_chr16 <- geneCentric_ncRNA_cond(agds_dir[5], chr = 16,
                                                             metab = "x100006264",
                                                             group = "exp_groups",
                                                             known_loci = known_loci_chr16_propyl)


results_noncoding_ncrna_uncond_chr2_arginine <- geneCentric_ncRNA_cond(agds_dir[6], chr = 2,
                                                                            metab = "x100001266",
                                                                            group = "exp_groups",
                                                                       known_loci = known_loci_chr2)

results_noncoding_ncrna_uncond_chr5_3amino <- geneCentric_ncRNA_cond(agds_dir[7], chr = 5,
                                                                          metab = "x1114",
                                                                          group = "exp_groups",
                                                                     known_loci = known_loci_chr5)

known_loci_chr8 <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x192//individual_cond_pruned_var.RDS"))
results_noncoding_ncrna_chr8_putrescine <- geneCentric_ncRNA_cond(agds_dir[8], chr = 8,
                                                                       metab = "x192",
                                                                       group = "exp_groups",
                                                                  known_loci = known_loci_chr8)


#------------------------------------------------------------------------------
# LEGACY / EARLIER VERSION: v2_multi_chr  (166 unique blocks)
#------------------------------------------------------------------------------
# Kept verbatim. These are blocks that do NOT appear in the current version above
# (mostly hard-coded per-metabolite / per-chromosome run calls and older path setups).

output_path <- "../Data/STAAR_prep/"
batch_info <- c("batch1", "batch2")
chr_info <- c(2, 5, 6, 8, 10, 13, 16)
agds_dir <- c(paste0(dir.geno,agds_file_name_1, chr_info, "_",
                     "batch1", agds_file_name_2),
              paste0(dir.geno,agds_file_name_1, chr_info, "_",
                     "batch2", agds_file_name_2)
)
save(agds_dir,file=paste0("../Data/STAAR_prep/agds_dir_original.Rdata",sep=""))
load("../Data/STAAR_prep/agds_dir_original.Rdata")


colnames(jobs_num) <- c("chr","start_loc","end_loc","individual_analysis_num","sliding_window_num","scang_num")

save(jobs_num,file=paste0(output_path,"mult_chr_jobs_num.Rdata",sep=""))
load(paste0(output_path,"mult_chr_jobs_num.Rdata"))

for(i in agds_dir[c(4, 11)]){
  agds <- seqOpen(i, readonly = FALSE)
  print(i)
  # we needto add AVGDP (average sequencing depth information,
  # here we will just set everyone to 10 because it's the default for FastSparseGRM
  Anno.folder <- index.gdsn(agds, "annotation/info")
  position <- read.gdsn(index.gdsn(agds, "position"))
  add.gdsn(Anno.folder, "AVGDP", val=rep(10, length(position)),
           compress="LZMA_ra", closezip=TRUE)
  seqClose(agds)

}

# installing dependencies
#devtools::install_github("xihaoli/STAAR")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
#devtools::install_github("xihaoli/MultiSTAAR")
#devtools::install_github("zilinli1988/SCANG")
#library(GenomicFeatures)
#devtools::install_github("xihaoli/STAARpipeline",ref="main")

phenotype <- readRDS("../Data/covariate_b1_all_metab.RDS")
for(i in 1:length(chr_info)){
  agds <- seqOpen(agds_dir[1:7][i])
  sample.id <- read.gdsn(index.gdsn(agds, "sample.id"))
  sample.id <- as.data.frame(sample.id)

  mat_b1_ppl <- as.data.frame(rownames(km))

  intersect_ppl <- sample.id |>
    inner_join(mat_b1_ppl, by = c("sample.id" = "rownames(km)"))

  seqSetFilter(agds, sample.id =intersect_ppl[,1] )
  output_path <- paste0("../Data/chr", chr_info[i], "_DRAGEN_agds_b1_final_1.gds")
  seqExport(agds,  out.fn = output_path)
  seqClose(agds)
}

for(i in 4:4){  #1:length(chr_info)){
  agds <- seqOpen(agds_dir[8:14][i])
  sample.id <- read.gdsn(index.gdsn(agds, "sample.id"))
  sample.id <- as.data.frame(sample.id)

  mat_b2_ppl <- as.data.frame(rownames(km))

  intersect_ppl <- sample.id |>
    inner_join(mat_b2_ppl, by = c("sample.id" = "rownames(km)"))

  seqSetFilter(agds, sample.id =intersect_ppl[,1] )
  output_path <- paste0("../Data/chr", chr_info[i], "_DRAGEN_agds_b2_final.gds")
  seqExport(agds,  out.fn = output_path)
  seqClose(agds)
}

phenotype <- phenotype[intersect_ppl[[1]],]
kinmat_b1 <- km[intersect_ppl[,1], intersect_ppl[,1]] # 3889 x 3889
saveRDS(kinmat_b1, "../Data/HCHS_SOL_kinship_matrix_b1.RDS")
pcs_b2 <- pcs[rownames(phenotype_b2),] #sum(startsWith(rownames(pcs_b2), "NWD"))


saveRDS(b1_all_cov_for_STAAR_null, "../Data/STAAR_model_cov_batch1.RDS")

b2_all_cov_for_STAAR_null <- cbind(phenotype, pcs_b2[,1:5])
saveRDS(b2_all_cov_for_STAAR_null, "../Data/STAAR_model_cov_batch2.RDS")

b1_all_cov_for_STAAR_null <- readRDS("../Data/STAAR_model_cov_batch1.RDS")
for(i in neg_control_ids){
  print(i)
  obj_nullmodel_GENESIS <- fitNullModel(x = b1_all_cov_for_STAAR_null,
                                        outcome=i,
                                        covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                 "PC4", "PC5",
                                                 "BKGRD1_C7", "GFRSCYS",
                                                 "CENTER"),
                                        cov.mat=covMatList,
                                        AIREML.tol=1e-4,
                                        verbose=TRUE)
  print(dim(obj_nullmodel_GENESIS$fit))
  saveRDS(obj_nullmodel_GENESIS, paste0("../Data/STAAR_neg_controls/",
                                        i, "/", i, "_nullmodel_GENESIS_batch1.RDS"))
  obj_nullmodel <- genesis2staar_nullmodel(obj_nullmodel_GENESIS)
  # obj_nullmodel$Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
  obj_nullmodel$Sigma_i <- as(obj_nullmodel$Sigma_i, "sparseMatrix")

  # save nullmodel
  save(obj_nullmodel,file= paste0("../Data/STAAR_neg_controls/",
                                  i, "/", i, "_nullmodel_batch1.Rdata"))
}

#olnames(b1_all_cov_for_STAAR_null)[2] <- "x1224"
for(i in exp_group_ids[4:4]){
  print(i)
  obj_nullmodel_GENESIS <- fitNullModel(x = b1_all_cov_for_STAAR_null,
                                        outcome=i,
                                        covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                                 "PC4", "PC5",
                                                 "BKGRD1_C7", "GFRSCYS",
                                                 "CENTER"),
                                        cov.mat=covMatList,
                                        AIREML.tol=1e-4,
                                        verbose=TRUE)

  saveRDS(obj_nullmodel_GENESIS, paste0("../Data/STAAR_exp_groups/",
                                        i, "/", i, "_nullmodel_GENESIS_batch1.RDS"))

  # obj_nullmodel_GENESIS <- readRDS(paste0("../Data/STAAR_exp_groups/",
  #                                         i, "/", i, "_nullmodel_GENESIS_batch1.RDS"))

  ## convert GENESIS null model to STAAR null model
  obj_nullmodel <- genesis2staar_nullmodel(obj_nullmodel_GENESIS)
  obj_nullmodel$Sigma_i <- as(obj_nullmodel$Sigma_i, "sparseMatrix")

  # save nullmodel
  save(obj_nullmodel,file= paste0("../Data/STAAR_exp_groups/",
                                  i, "/", i, "_nullmodel_batch1.Rdata"))
}


jobs_num <- get(load("../Data/STAAR_prep/mult_chr_jobs_num.Rdata"))
agds_dir <- get(load("../Data/STAAR_prep/agds_dir_multichr.Rdata"))

#arrayid <- as.numeric(commandArgs(TRUE)[1])

# NEGATIVE CONTROLS
chr <- c(2, 5, 6, 16)
group.num <- jobs_num$individual_analysis_num[jobs_num[,1] %in% chr]

STAAR_ind_analysis <- function(agds.path, group, metab, chr){
  if(group == "exp_groups"){
    load(paste0("../Data/STAAR_exp_groups/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  } else {
    load(paste0("../Data/STAAR_neg_controls/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }
  genofile <- seqOpen(agds.path)
  job_num <- jobs_num[jobs_num$chr == chr,]
  start_loc <- job_num$start_loc
  end_loc <- job_num$end_loc
  obj_nullmodel$sparse_kins <- TRUE
  a <- Sys.time()
  results_individual_analysis <- c()
  if(start_loc <= end_loc)
  {
    results_individual_analysis <- Individual_Analysis(chr=chr,start_loc=start_loc,end_loc=end_loc,
                                                       genofile=genofile,obj_nullmodel=obj_nullmodel,mac_cutoff=20,
                                                       QC_label=QC_label,variant_type=variant_type,
                                                       geno_missing_imputation=geno_missing_imputation)
  }
  b <- Sys.time()
  b - a
  save(results_individual_analysis,file=paste0(output_path, metab, "/",metab, "_",
                                               output_file_name,".Rdata"))

  seqClose(genofile)
  return(results_individual_analysis)
}

agds_file_name_1 <- "chr"
agds_file_name_2 <- "DRAGEN_agds_"
agds_file_name_3 <- "_final.gds"
agds_dir <- c(paste0(dir.geno,agds_file_name_1, chr_info, "_",
                     agds_file_name_2, "b1",agds_file_name_3),
              paste0(dir.geno,agds_file_name_1, chr_info, "_",
                     agds_file_name_2,"b2", agds_file_name_3)
)
save(agds_dir,file=paste0("../Data/STAAR_prep/agds_dir_multichr.Rdata",sep=""))
load("../Data/STAAR_prep/agds_dir_multichr.Rdata")

## output path - negative controls
output_path <- "../Data/STAAR_neg_controls/"

# N2-acetyllysine
res_individual_chr2 <- STAAR_ind_analysis(agds_dir[1], "neg_controls",
                                          "x100001721", chr = 2)
# Betaine
res_individual_chr5 <- STAAR_ind_analysis(agds_dir[2], "neg_controls",
                                          "x799", chr = 5)
# N-acetylglucosaminylasparagine
res_individual_chr6 <- STAAR_ind_analysis(agds_dir[3], "neg_controls",
                                          "x1215", chr = 6)
# Cysteinylglycine
res_individual_chr16 <- STAAR_ind_analysis(agds_dir[7], "neg_controls",
                                           "x278", chr = 16)

# 5-oxoproline
res_individual_chr8 <- STAAR_ind_analysis(agds_dir[4], "exp_groups",
                                          "x1021", chr = 8)
# Carnitine
res_individual_chr10 <- STAAR_ind_analysis(agds_dir[5], "exp_groups",
                                           "x100000007", chr = 10)
# N-acetylcarnosine
res_individual_chr13 <- STAAR_ind_analysis(agds_dir[6], "exp_groups",
                                           "x100004046", chr = 13)
# Cys-gly, oxidized
res_individual_chr16_exp <- STAAR_ind_analysis(agds_dir[7], "exp_groups",
 "x1224", chr = 16)

geneCentric_coding_uncond <- function(agds.path, chr, group, metab){
  if(group == "exp_groups"){
    ## output path
    output_path <- "../Data/STAAR_exp_groups/"
    load(paste0("../Data/STAAR_exp_groups/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }else{
    output_path <- "../Data/STAAR_neg_controls/"
    load(paste0("../Data/STAAR_neg_controls/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }

  genofile <- seqOpen(agds.path)
  position <- read.gdsn(index.gdsn(genofile, "position"))
  genes_info_chr <- genes_info[genes_info[,2]==chr,]
  genes_info_chr <- as.data.frame(genes_info_chr) |>
    filter(end_position >= min(position)) |>
    filter(start_position <= max(position))
  sub_seq_num <- dim(genes_info_chr)[1]

  results_coding <- c()
  categories=c("plof","plof_ds","missense","disruptive_missense","synonymous","ptv","ptv_ds")
  for(kk in 1:nrow(genes_info_chr))
  {
    gene_name <- genes_info_chr[kk,1]
    print(gene_name)
    results <- Gene_Centric_Coding(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                          rare_maf_cutoff=0.01,rv_num_cutoff=2, category = "all_categories",
                                          QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                          Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                          Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)


    results_coding <- append(results_coding,results)

  }

  save(results_coding,file=paste0(output_path, metab, "/",
                                  metab, "_", output_file_name,".Rdata"))

  seqClose(genofile)

  return(results_coding)
}


results_coding_chr2_uncond <- geneCentric_coding_uncond(agds_dir[1], chr = 2,
                                                 metab = "x100001721",
                                                 group = "neg_controls")

results_coding_chr5_uncond <- geneCentric_coding_uncond(agds_dir[2], chr = 5, metab = "x799",
                                          group = "neg_controls")

results_coding_chr6_uncond <- geneCentric_coding_uncond(agds_dir[3], chr = 6, metab = "x1215",
                                          group = "neg_controls")
results_coding_chr16_uncond <- geneCentric_coding_uncond(agds_dir[7], chr = 16, metab = "x278",
                                           group = "neg_controls")

results_coding_chr8_uncond <- geneCentric_coding_uncond(agds_dir[4], group = "exp_groups",
                                          metab = "x1021", chr = 8)
results_coding_chr10_uncond <- geneCentric_coding_uncond(agds_dir[5], group ="exp_groups",
                                           metab = "x100000007", chr = 10)
results_coding_chr13_uncond <- geneCentric_coding_uncond(agds_dir[6], group ="exp_groups",
                                           metab = "x100004046", chr = 13)
results_coding_chr16_uncond_exp <- geneCentric_coding_uncond(agds_dir[7], group = "exp_groups",
                                               metab = "x1224", chr = 16)

# MULTIPLE TESTING BURDEN COMPUTATION
# number of independent tests actually performed (non-null)
sum(!sapply(results_coding_chr10_uncond, is.null)) +
  sum(!sapply(results_coding_chr13_uncond, is.null))+
  sum(!sapply(results_coding_chr16_uncond, is.null))+
  sum(!sapply(results_coding_chr16_uncond_exp, is.null))+
  sum(!sapply(results_coding_chr2_uncond, is.null))+
  sum(!sapply(results_coding_chr6_uncond, is.null))+
  sum(!sapply(results_coding_chr8_uncond, is.null))

#alpha_coding <- 0.05/735

known_loci_chr2 <- as.data.frame(readRDS("../Data/STAAR_neg_controls/x100001721/individual_cond_pruned_var.RDS"))
results_coding_chr2 <- geneCentric_coding_cond(agds_dir[1], chr = 2, metab = "x100001721", group = "neg_controls",
                                   known_loci = known_loci_chr2)

known_loci_chr5 <- as.data.frame(readRDS("../Data/STAAR_neg_controls/x799/individual_cond_pruned_var.RDS"))
results_coding_chr5 <- geneCentric_coding_cond(agds_dir[2], chr = 5, metab = "x799", group = "neg_controls",
                                   known_loci = known_loci_chr5)

known_loci_chr6 <- as.data.frame(readRDS("../Data/STAAR_neg_controls/x1215/individual_cond_pruned_var.RDS"))
results_coding_chr6 <- geneCentric_coding_cond(agds_dir[3], chr = 6, metab = "x1215", group = "neg_controls",
                                   known_loci = known_loci_chr6)

known_loci_chr16 <- as.data.frame(readRDS("../Data/STAAR_neg_controls/x278/individual_cond_pruned_var.RDS"))
results_coding_chr16 <- geneCentric_coding_cond(agds_dir[7], chr = 16, metab = "x278", group = "neg_controls",
                                    known_loci = known_loci_chr16)

known_loci_chr8 <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x1021/individual_cond_pruned_var.RDS"))
results_coding_chr8 <- geneCentric_coding_cond(agds_dir[4], group = "exp_groups",
                                                        metab = "x1021", chr = 8,
                                               known_loci = known_loci_chr8)
known_loci_chr10 <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x100000007/individual_cond_pruned_var.RDS"))
results_coding_chr10 <- geneCentric_coding_cond(agds_dir[5], group ="exp_groups",
                                                         metab = "x100000007", chr = 10,
                                                known_loci = known_loci_chr10)
known_loci_chr13 <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x100004046/individual_cond_pruned_var.RDS"))
results_coding_chr13 <- geneCentric_coding_cond(agds_dir[6], group ="exp_groups",
                                                         metab = "x100004046", chr = 13,
                                                known_loci = known_loci_chr13)
known_loci_chr16_exp <- as.data.frame(readRDS("../Data/STAAR_exp_groups/x1224/individual_cond_pruned_var.RDS"))
results_coding_chr16_exp <- geneCentric_coding_cond(agds_dir[7], group = "exp_groups",
                                                             metab = "x1224", chr = 16,
                                                    known_loci = known_loci_chr16_exp)

#####################################################################
# do gene-centric non-coding analysis for each of the gene
# Gene-centric analysis for noncoding rare variants of protein-coding
# genes using STAARpipeline
# rm(list=ls())
gc()

geneCentric_noncoding <- function(agds.path, chr, group, metab){
  if(group == "exp_groups"){
    ## output path
    output_path <- "../Data/STAAR_exp_groups/"
    load(paste0("../Data/STAAR_exp_groups/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }else{
    output_path <- "../Data/STAAR_neg_controls/"
    load(paste0("../Data/STAAR_neg_controls/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }
  genofile <- seqOpen(agds.path)
  position <- read.gdsn(index.gdsn(genofile, "position"))
  genes_info_chr <- genes_info[genes_info[,2]==chr,]
  genes_info_chr <- as.data.frame(genes_info_chr) |>
    filter(end_position >= min(position)) |>
    filter(start_position <= max(position))
  sub_seq_num <- dim(genes_info_chr)[1]


  results_noncoding <- c()
  for(kk in 1:nrow(genes_info_chr))
  {
    print(kk)
    gene_name <- genes_info_chr[kk,1]
    print(gene_name)
    results <- Gene_Centric_Noncoding(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                      rare_maf_cutoff=0.01,rv_num_cutoff=2, category = "all_categories",
                                      QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                      Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                      Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
    results_noncoding <- append(results_noncoding,results)

  }

  save(results_noncoding,file=paste0(output_path, metab, "/",
                                     metab, "_", output_file_name,".Rdata"))

  seqClose(genofile)
  return(results_noncoding)
}

results_noncoding_chr2_uncond <- geneCentric_noncoding(agds_dir[1], chr = 2,
                                                 metab = "x100001721",
                                                 group = "neg_controls")

results_noncoding_chr5_uncond <- geneCentric_noncoding(agds_dir[2], chr = 5, metab = "x799",
                                                 group = "neg_controls")

results_noncoding_chr6_uncond <- geneCentric_noncoding(agds_dir[3], chr = 6, metab = "x1215",
                                                 group = "neg_controls")

results_noncoding_chr16_uncond <- geneCentric_noncoding(agds_dir[7], chr = 16, metab = "x278",
                                                  group = "neg_controls")

# run this
results_noncoding_chr8_uncond <- geneCentric_noncoding(agds_dir[4], group = "exp_groups",
                                                     metab = "x1021", chr = 8)

results_noncoding_chr10_uncond <- geneCentric_noncoding(agds_dir[5], group ="exp_groups",
                                                      metab = "x100000007", chr = 10)

results_noncoding_chr13_uncond <- geneCentric_noncoding(agds_dir[6], group ="exp_groups",
                                                      metab = "x100004046", chr = 13)

results_noncoding_chr16_uncond_exp <- geneCentric_noncoding(agds_dir[7], group = "exp_groups",
                                                            metab = "x1224", chr = 16)


sum(!sapply(results_noncoding_chr10_uncond, is.null)) +
  sum(!sapply(results_noncoding_chr13_uncond, is.null))+
  sum(!sapply(results_noncoding_chr16_uncond, is.null))+
  sum(!sapply(results_noncoding_chr16_uncond_exp, is.null))+
  sum(!sapply(results_noncoding_chr2_uncond, is.null))+
  sum(!sapply(results_noncoding_chr6_uncond, is.null))+
  sum(!sapply(results_noncoding_chr8_uncond, is.null))

alpha_noncoding <- 0.05/1166

geneCentric_noncoding_cond <- function(agds.path, chr, group, metab, known_loci){
  if(group == "exp_groups"){
    ## output path
    output_path <- "../Data/STAAR_exp_groups/"
    load(paste0("../Data/STAAR_exp_groups/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }else{
    output_path <- "../Data/STAAR_neg_controls/"
    load(paste0("../Data/STAAR_neg_controls/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }
  genofile <- seqOpen(agds.path)
  position <- read.gdsn(index.gdsn(genofile, "position"))
  genes_info_chr <- genes_info[genes_info[,2]==chr,]
  genes_info_chr <- as.data.frame(genes_info_chr) |>
    filter(end_position >= min(position)) |>
    filter(start_position <= max(position))
  sub_seq_num <- dim(genes_info_chr)[1]

  genes <- genes_info

  categories=c("downstream","upstream","UTR","promoter_CAGE",
               "promoter_DHS","enhancer_CAGE","enhancer_DHS")
  for(category in categories){
    results_coding <- c()
    print(category)
    for(kk in 1:nrow(genes_info_chr))
    {
      gene_name <- genes_info_chr[kk,1]
      print(gene_name)

      results <- Gene_Centric_Noncoding_cond(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                             rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                             known_loci=known_loci, category = category,
                                             QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                             Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                             Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
      results_noncoding <- rbind(results_coding,results)
    }
    save(results_noncoding,file=paste0(output_path, metab, "/",
                                       metab, "_", output_file_name, "_",
                                       category,
                                       ".Rdata"))
  }
  seqClose(genofile)

}

## output path - negative control
results_noncoding_chr2 <- geneCentric_noncoding_cond(agds_dir[1], chr = 2, metab = "x100001721", group = "neg_controls",
                                               known_loci = known_loci_chr2)

results_noncoding_chr5 <- geneCentric_noncoding_cond(agds_dir[2], chr = 5, metab = "x799", group = "neg_controls",
                                               known_loci = known_loci_chr5)

results_noncoding_chr6 <- geneCentric_noncoding_cond(agds_dir[3], chr = 6, metab = "x1215", group = "neg_controls",
                                               known_loci = known_loci_chr6)

results_noncoding_chr16 <- geneCentric_noncoding_cond(agds_dir[7], chr = 16, metab = "x278", group = "neg_controls",
                                                known_loci = known_loci_chr16)

results_noncoding_chr8 <- geneCentric_noncoding_cond(agds_dir[4], group = "exp_groups",
                                               metab = "x1021", chr = 8,
                                               known_loci = known_loci_chr8)
results_noncoding_chr10 <- geneCentric_noncoding_cond(agds_dir[5], group ="exp_groups",
                                                metab = "x100000007", chr = 10,
                                                known_loci = known_loci_chr10)
results_noncoding_chr13 <- geneCentric_noncoding_cond(agds_dir[6], group ="exp_groups",
                                                metab = "x100004046", chr = 13,
                                                known_loci = known_loci_chr13)
results_noncoding_chr16_exp <- geneCentric_noncoding_cond(agds_dir[7], group = "exp_groups",
                                                    metab = "x1224", chr = 16,
                                                    known_loci = known_loci_chr16_exp)


geneCentric_ncRNA_noncoding <- function(agds.path, chr, group, metab){
  if(group == "exp_groups"){
    ## output path
    output_path <- paste0("../Data/STAAR_exp_groups/", metab, "/")
    load(paste0("../Data/STAAR_exp_groups/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }else{
    output_path <- paste0("../Data/STAAR_neg_controls/", metab, "/")
    load(paste0("../Data/STAAR_neg_controls/",
                metab, "/", metab, "_nullmodel_batch1.Rdata"))
  }
  genofile <- seqOpen(agds.path)
  position <- read.gdsn(index.gdsn(genofile, "position"))
  ncRNA_gene_chr <- ncRNA_gene[ncRNA_gene[,1]==chr,]

  GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
  GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[;]")
  Gene <- as.character(sapply(GENCODE.Info.split,function(z) gsub("\\(.*\\)","",z[1])))

  Gene_list_1 <- as.character(sapply(strsplit(Gene,','),'[',1))
  Gene_list_2 <- as.character(sapply(strsplit(Gene,','),'[',2))
  Gene_list_3 <- as.character(sapply(strsplit(Gene,','),'[',3))

  gene_names <- ncRNA_gene_chr[which(ncRNA_gene_chr[,2] %in% Gene_list_1),]
  gene_names <- rbind(gene_names, ncRNA_gene_chr[which(ncRNA_gene_chr[,2] %in% Gene_list_2),])
  gene_names <- rbind(gene_names, ncRNA_gene_chr[which(ncRNA_gene_chr[,2] %in% Gene_list_3),])
  gene_names <- gene_names[!is.na(gene_names$ncRNA) & !gene_names$ncRNA == "NA" ,]
  ncRNA_gene_chr <- unique(gene_names)
  print(table(ncRNA_gene_chr))
  if(dim(ncRNA_gene_chr)[1] == 0){
    results_ncRNA_cond <- NULL
    save(results_ncRNA_cond,file=paste0(output_path,output_file_name,".Rdata"))
    return(NULL)
  }

  results_ncRNA <- c()
  for(kk in 1:nrow(ncRNA_gene_chr))
  {
    gene_name <- ncRNA_gene_chr[kk,2]
    results <- c()
    results <- try(ncRNA(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                         rare_maf_cutoff=0.01,rv_num_cutoff=2,
                         QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                         Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name))
    results_ncRNA <- rbind(results_ncRNA,results)
  }
  seqClose(genofile)
  save(results_ncRNA,file=paste0(output_path,output_file_name,".Rdata"))
  return(results_ncRNA)
}

results_noncoding_ncrna_chr2_uncond <- geneCentric_ncRNA_noncoding(agds_dir[1], chr = 2,
                                                       metab = "x100001721", group = "neg_controls")

results_noncoding_ncrna_chr5_uncond <- geneCentric_ncRNA_noncoding(agds_dir[2], chr = 5, metab = "x799",
                                                       group = "neg_controls")

results_noncoding_ncrna_chr6_uncond <- geneCentric_ncRNA_noncoding(agds_dir[3], chr = 6, metab = "x1215",
                                                       group = "neg_controls")

results_noncoding_ncrna_chr16_uncond <- geneCentric_ncRNA_noncoding(agds_dir[7], chr = 16, metab = "x278",
                                                        group = "neg_controls")

results_noncoding_ncrna_chr8_uncond <- geneCentric_ncRNA_noncoding(agds_dir[4], group = "exp_groups",
                                                       metab = "x1021", chr = 8)

results_noncoding_ncrna_chr10_uncond <- geneCentric_ncRNA_noncoding(agds_dir[5], group ="exp_groups",
                                                        metab = "x100000007", chr = 10)

results_noncoding_ncrna_chr13_uncond <- geneCentric_ncRNA_noncoding(agds_dir[6], group ="exp_groups",
                                                        metab = "x100004046", chr = 13)

results_noncoding_ncrna_chr16_uncond_exp <- geneCentric_ncRNA_noncoding(agds_dir[7], group = "exp_groups",
                                                            metab = "x1224", chr = 16)


################## CONDITIONAL ANALYSIS - NCRNA ############
source("./20250918_ncRNA_cond_modified.R")
results_noncoding_ncrna_chr2_cond <- geneCentric_ncRNA_cond(agds_dir[1], chr = 2,
                                                                   metab = "x100001721",
                                                            group = "neg_controls",
                                                            known_loci = known_loci_chr2)

results_noncoding_ncrna_chr5_cond <- geneCentric_ncRNA_cond(agds_dir[2], chr = 5, metab = "x799",
                                                                   group = "neg_controls",
                                                            known_loci = known_loci_chr5)

results_noncoding_ncrna_chr6_cond <- geneCentric_ncRNA_cond(agds_dir[3], chr = 6, metab = "x1215",
                                                                   group = "neg_controls",
                                                            known_loci = known_loci_chr6)

results_noncoding_ncrna_chr16_cond <- geneCentric_ncRNA_cond(agds_dir[7], chr = 16, metab = "x278",
                                                                    group = "neg_controls",
                                                             known_loci = known_loci_chr16)


results_noncoding_ncrna_chr8_cond <- geneCentric_ncRNA_cond(agds_dir[4], group = "exp_groups",
                                                                   metab = "x1021", chr = 8,
                                                            known_loci = known_loci_chr8,
                                                            uncond_res = results_noncoding_ncrna_chr8_uncond)

results_noncoding_ncrna_chr10_cond <- geneCentric_ncRNA_cond(agds_dir[5], group ="exp_groups",
                                                                    metab = "x100000007", chr = 10,
                                                             known_loci = known_loci_chr10)

results_noncoding_ncrna_chr13_cond <- geneCentric_ncRNA_cond(agds_dir[6], group ="exp_groups",
                                                                    metab = "x100004046", chr = 13,
                                                             known_loci = known_loci_chr13)

results_noncoding_ncrna_chr16_cond_exp <- geneCentric_ncRNA_cond(agds_dir[7], group = "exp_groups",
                                                                        metab = "x1224", chr = 16,
                                                                 known_loci = known_loci_chr16_exp)

# fit SCANG-STAAR null model
library(SCANG)
staar_nullmodel_path <- "../Data/obj_nullmodel.Rdata"
scang_staar_nullmodel_path <- "../Data/obj_nullmodel_SCANG_STAAR.Rdata"

## load STAAR null model
obj_nullmodel <- get(load(staar_nullmodel_path))

obj_nullmodel_SCANG_STAAR <- staar2scang_nullmodel(obj_nullmodel)

save(obj_nullmodel_SCANG_STAAR,file=scang_staar_nullmodel_path)

# Dynamic window analysis using STAARpipeline
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 02/16/2024
rm(list=ls())
## Null model
obj_nullmodel_SCANG_STAAR <- get(load("../Data/obj_nullmodel_SCANG_STAAR.Rdata"))

output_file_name <- "chr16_dynamc_window_SCANG"

## Number of jobs for SCANG
sum(jobs_num$scang_num)

chr <- 16
group.num <- jobs_num$scang_num

## aGDS file
genofile <- seqOpen(agds_dir)

position <- read.gdsn(index.gdsn(genofile, "position"))
start_loc <- min(position)
end_loc <- max(position)

results_scang <- Dynamic_Window_SCANG(chr=chr,start_loc=start_loc,end_loc=end_loc,genofile=genofile,obj_nullmodel=obj_nullmodel_SCANG_STAAR,
                                      QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                      Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                      Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

save(results_scang,file=paste0(output_path,output_file_name,"_",1,".Rdata"))

seqClose(genofile)


#------------------------------------------------------------------------------
# LEGACY / EARLIER VERSION: v1_single_chr16  (115 unique blocks)
#------------------------------------------------------------------------------
# Kept verbatim. These are blocks that do NOT appear in the current version above
# (mostly hard-coded per-metabolite / per-chromosome run calls and older path setups).

agds_file_name_1 <- "chr16_DRAGEN_agds_b1"
agds_file_name_2 <- ".gds"
agds_dir <- paste0(dir.geno,agds_file_name_1,agds_file_name_2)
save(agds_dir,file=paste0("../Data/STAAR_prep/agds_dir.Rdata",sep=""))

jobs_num <- matrix(rep(0,3),nrow=1)
for(chr in 16:16)
{
  print(chr)
  gds.path <- agds_dir
  genofile <- seqOpen(gds.path)

  filter <- seqGetData(genofile, QC_label)
  SNVlist <- filter == "PASS"

  position <- as.numeric(seqGetData(genofile, "position"))

  jobs_num[1,1] <- chr
  jobs_num[1,2] <- min(position[SNVlist])
  jobs_num[1,3] <- max(position[SNVlist])

  seqClose(genofile)
}

save(jobs_num,file=paste0(output_path,"jobs_num.Rdata",sep=""))

agds <- seqOpen(agds_dir, readonly = FALSE)

# we needto add AVGDP (average sequencing depth information,
# here we will just set everyone to 10 because it's the default for FastSparseGRM
Anno.folder <- index.gdsn(agds, "annotation/info")
position <- read.gdsn(index.gdsn(agds, "position"))
add.gdsn(Anno.folder, "AVGDP", val=rep(10, length(position)),
         compress="LZMA_ra", closezip=TRUE)
seqClose(agds)

## Phenotype file
phenotype <- readRDS("../Data/metab_covariates_batch1.RDS")
agds <- seqOpen(agds_dir)
seqSetFilter(agds, )
sample.id <- read.gdsn(index.gdsn(agds, "sample.id"))
sample.id <- as.data.frame(sample.id)
mat_b1_ppl <- as.data.frame(rownames(km))

intersect_ppl <- sample.id |>
  inner_join(mat_b1_ppl, by = c("sample.id" = "rownames(km)"))

seqSetFilter(agds, sample.id =intersect_ppl[,1] )
output_path <- paste0("../Data/chr16_DRAGEN_agds_b1_final.gds")
seqExport(agds,  out.fn = output_path)í

data_GENESIS <- as(phenotype,"AnnotatedDataFrame") # Make AnnotatedDataFrame (specifically required by GENESIS)


# adjust 11 PCs as indicated in the file??
b1_all_cov_for_STAAR_null <- cbind(phenotype, pcs_b1[,1:10])
colnames(b1_all_cov_for_STAAR_null)[2] <- "x1224"
obj_nullmodel_GENESIS <- fitNullModel(x = b1_all_cov_for_STAAR_null,
                                      outcome="x1224",
                                      covars=c("AGE","GENDER","PC1", "PC2", "PC3",
                                               "PC4", "PC5","PC6", "PC7", "PC8",
                                               "PC9", "PC10",
                                               "BKGRD1_C7", "GFRSCYS",
                                               "CENTER"),
                                      cov.mat=covMatList,
                                      AIREML.tol=1e-4,
                                      verbose=TRUE)

saveRDS(obj_nullmodel_GENESIS, "../Data/STAAR_chr16_1224/obj_nullmodel_GENESIS.RDS")

## convert GENESIS null model to STAAR null model
# load("../Data/null_mod_b1_chr16_cond.Rdata")
obj_nullmodel <- genesis2staar_nullmodel(obj_nullmodel_GENESIS)
# obj_nullmodel$Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
obj_nullmodel$Sigma_i <- as(obj_nullmodel$Sigma_i, "sparseMatrix")

# save nullmodelí
save(obj_nullmodel,file= "../Data/STAAR_chr16_1224/obj_nullmodel.Rdata")


jobs_num <- get(load("../Data/STAAR_prep/jobs_num.Rdata"))
agds_dir <- get(load("../Data/STAAR_prep/agds_dir.Rdata"))
#obj_nullmodel <- get(load("../Data/obj_nullmodel_GENESIS.Rdata"))

output_path <- "../Data"
output_file_name <- "chr16_Individual_Analysis"
group.num <- jobs_num$individual_analysis_num

# if (chr == 1){
#   groupid <- arrayid
# }else{
#   groupid <- arrayid - cumsum(jobs_num$individual_analysis_num)[chr-1]
agds.path <- agds_dir
genofile <- seqOpen(agds.path)

start_loc <- jobs_num$start_loc
end_loc <- jobs_num$end_loc


obj_nullmodel$sparse_kins <- TRUE

a <- Sys.time()
results_individual_analysis <- c()
if(start_loc <= end_loc)
{
  results_individual_analysis <- Individual_Analysis(chr=chr,start_loc=start_loc,end_loc=end_loc,
                                                     genofile=genofile,obj_nullmodel=obj_nullmodel,mac_cutoff=20,
                                                     QC_label=QC_label,variant_type=variant_type,
                                                     geno_missing_imputation=geno_missing_imputation)
}
b <- Sys.time()
b - a

save(results_individual_analysis,file=paste0(output_path,output_file_name,".Rdata"))

output_file_name <- "chr16_genecentric_cond_summary_1"
## gene number in job
# gene_num_in_array <- 50
# group.num.allchr <- ceiling(table(genes_info[,2])/gene_num_in_array)
# sum(group.num.allchr)
# group.num <- group.num.allchr[chr]
#   groupid <- 285 - cumsum(group.num.allchr)[chr-1]
# cumsum(group.num.allchr)[17]
genes_info_chr <- genes_info[genes_info[,2]==16,]
genes_info_chr <- genes_info_chr |>
  filter(end_position >= min(position)) |>
  filter(start_position <= max(position))
sub_seq_num <- dim(genes_info_chr)[1] # 46 genes
# if(groupid < group.num)
# {
#   sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):(groupid*gene_num_in_array)
# }else
#   sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):sub_seq_num
# ## exclude large coding masks
# if(arrayid==57)
#   sub_seq_id <- setdiff(sub_seq_id,840)
# if(arrayid==112)
#   sub_seq_id <- setdiff(sub_seq_id,c(543,544))
# if(arrayid==113)
#   sub_seq_id <- setdiff(sub_seq_id,c(575,576,577,578,579,580,582))
genes <- genes_info

results_coding <- c()
known_loci <- readRDS("../Data/known_loci_info_8variants.RDS")
for(kk in 1:nrow(genes_info_chr))
{
  print(kk)
  gene_name <- genes_info_chr[kk,1]
  results <- Gene_Centric_Coding_cond(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                      rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                      known_loci=known_loci,
                                      QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                      Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                      Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  results_coding <- append(results_coding,results)
}


save(results_coding,file=paste0(output_path,output_file_name,".Rdata")) # 46 genes * 5 annotation categories

output_file_name <- "chr16_genecentric_noncoding"
#   groupid <- arrayid - cumsum(group.num.allchr)[chr-1]
# genes_info_chr <- genes_info[genes_info[,2]==chr,]
# sub_seq_num <- dim(genes_info_chr)[1]
# ## exclude large noncoding masks
# jobid_exclude <- c(21,39,44,45,46,53,55,83,88,103,114,127,135,150,154,155,163,164,166,180,189,195,200,233,280,285,295,313,318,319,324,327,363,44,45,54)
# sub_seq_id_exclude <- c(1009,1929,182,214,270,626,741,894,83,51,611,385,771,493,671,702,238,297,388,352,13,303,600,170,554,207,724,755,1048,319,324,44,411,195,236,677)
# for(i in 1:length(jobid_exclude))
#   if(arrayid==jobid_exclude[i])
#   {
#     sub_seq_id <- setdiff(sub_seq_id,sub_seq_id_exclude[i])
#   }
results_noncoding <- c()
{
  print(kk)
  gene_name <- genes_info_chr[kk,1]
  results <- Gene_Centric_Noncoding(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                    rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                    QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                    Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                    Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  results_noncoding <- append(results_noncoding,results)
}

save(results_noncoding,file=paste0(output_path,output_file_name,"_1.Rdata"))

# run gene-centric analysis for ncRNA
output_file_name <- "chr16_subregion_ncRNA_Noncoding"
ncRNA_gene_chr <- ncRNA_gene[ncRNA_gene[,1]==16,]
ncRNA_gene_chr <- ncRNA_gene_chr[rownames(ncRNA_gene_chr) %in% 16871:16923,]

results_ncRNA <- c()
for(kk in 1:nrow(ncRNA_gene_chr))
{
  print(kk)
  gene_name <- ncRNA_gene_chr[kk,2]
  results <- c()
  results <- try(ncRNA(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                       rare_maf_cutoff=0.01,rv_num_cutoff=2,
                       QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                       Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                       Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name))
  results_ncRNA <- rbind(results_ncRNA,results)
}
save(results_ncRNA,file=paste0(output_path,output_file_name,".Rdata"))


staar_nullmodel_path <- "../Data/STAAR_exp_groups/x1224/x1224_nullmodel_batch1.Rdata"
scang_staar_nullmodel_path <- "../Data/STAAR_exp_groups/x1224/x1224_nullmodel_SCANG_STAAR.Rdata"
