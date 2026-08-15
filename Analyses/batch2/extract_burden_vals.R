library(SeqArray)
library(tidyverse)
library(STAARpipeline)



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



x1224_plof_ds_rvs <- readRDS("../Data/STAAR_exp_groups/x1224/RVsetInfo_coding_sig_cond_plof_ds.RDS")

genofile_b2 <- seqOpen("../Data/chr16_DRAGEN_agds_b2_final.gds")

position <- as.data.frame(seqGetData(genofile_b2, "position"))
variant_id <- seqGetData(genofile_b2, "variant.id")
variant_id <- variant_id[which(position[,1] %in% x1224_plof_ds_rvs$POS)]

seqSetFilter(genofile_b2, variant.id = variant_id)
Geno <- seqGetData(genofile_b2, "$dosage")
Geno <- matrix_flip(Geno)

# two ways -- if does not exist, use an equivalent RV sets 
# if exists, use the same set 




