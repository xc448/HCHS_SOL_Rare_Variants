# load packages
library(tidyverse)
library(SeqArray)
library(gdsfmt)
library(SeqVarTools)

setwd("R:\\Sofer Lab/HCHS_SOL/Projects/2024_rare_variants/Code/")
# Load the liftOver files for interval info in the RFMix
# intervals_liftover_nodup <- read.table("/Volumes/Sofer Lab/HCHS_SOL/Projects/2023_local_ancestry_comparison_sol/Data/intervals_37-38.bed")

intervals_liftover_nodup <- read.table("R:\\Sofer Lab/HCHS_SOL/Projects/2023_local_ancestry_comparison_sol/Data/intervals_37-38.bed")
colnames(intervals_liftover_nodup) <- c("chr", "pos_start", "pos_end", "snpid") 
intervals_liftover_nodup <- intervals_liftover_nodup |> filter(chr != "chrX") # 14753 mapped blocks
saveRDS(intervals_liftover_nodup, "../Data/intervals_liftover_nodup.RDS")
intervals_liftover_nodup <- readRDS("../Data/intervals_liftover_nodup.RDS")
# subsetting the significant intervals from the RFMix inference data based on ancestry intervals

# subset_sig_region_gds <- function(chr, sig_interval_start, sig_interval_end , batch = 1){
#   
#   intervals_liftover_nodup <- readRDS("../Data/intervals_liftover_nodup.RDS")
#   
#   # subsetting the significant intervals from the RFMix inference data based on ancestry intervals
#   # include all intervals from start to end
#   known_sig_interval <- intervals_liftover_nodup |> filter(snpid %in% sig_interval_start:sig_interval_end)
#   g <- seqOpen(paste0("../Data/cardiac_cohorts_SOL_dragen.chr", chr, ".subset.gds"), readonly = FALSE)
#   prov_nwd_id_matched <- readRDS("../Data/prov_nwd_id_matched.RDS")
#   add.gdsn(g, "sample.id", val=prov_nwd_id_matched, compress="LZMA_ra", 
#            closezip=TRUE, replace = TRUE)
#   if(batch == 1){
#     # load matched sample IDs
#     sampleid <- readRDS("../Data/matched_NWD_SoLID_batch1_2_combined.RDS") |> 
#       filter(batch == "BATCH1") # 3977 ppl 
#   } else {
#     sampleid <- readRDS("../Data/matched_NWD_SoLID_batch1_2_combined.RDS") |> 
#       filter(batch == "BATCH2") 
#   }
#   
#   pos_ind <- c()
#   position <- as.data.frame(seqGetData(g, "position"))
#   pos_ind <- which(position >= min(as.numeric(known_sig_interval$pos_start)) &
#                      position <= max(as.numeric(known_sig_interval$pos_end)))
#   
#   variant_id <- seqGetData(g, "variant.id")
#   variant_id <- variant_id[pos_ind]
#   
#   seqSetFilter(g, sample.id = sampleid$NWD_ID, variant.id = variant_id)
#   
#   # export the subsetted gds file 
#   seqExport(g, out.fn = paste0("../Data/cardiac_cohorts_SOL_dragen.chr", chr,
#                                "batch", batch, 
#                                "_admixmap_sigregion.gds"))
#   seqClose(g)
# }
# # NEGATIVE CTRLS b1
# subset_sig_region_gds(2, sig_interval_start = 1575, sig_interval_end = 1585, batch = 1)
# subset_sig_region_gds(5, sig_interval_start = 4544, sig_interval_end = 4546, batch = 1)
# subset_sig_region_gds(6, sig_interval_start = 5848, sig_interval_end = 5848, batch = 1)
# #subset_sig_region_gds(16, sig_interval_start = 12330, sig_interval_end = 12348, batch = 1)
# 
# # EXPERIMENT GRPS b1
# subset_sig_region_gds(8, sig_interval_start = 7361, sig_interval_end = 7377, batch = 1)
# subset_sig_region_gds(10, sig_interval_start = 8402, sig_interval_end = 8412, batch = 1)
# subset_sig_region_gds(13, sig_interval_start = 10635, sig_interval_end = 10639, batch = 1)
# subset_sig_region_gds(16, sig_interval_start = 12330, sig_interval_end = 12348, batch = 1)
# 
# # NEGATIVE CTRLS b2
# subset_sig_region_gds(2, sig_interval_start = 1575, sig_interval_end = 1585, batch = 2)
# subset_sig_region_gds(5, sig_interval_start = 4544, sig_interval_end = 4546, batch = 2)
# subset_sig_region_gds(6, sig_interval_start = 5848, sig_interval_end = 5848, batch = 2)
# #subset_sig_region_gds(16, sig_interval_start = 12330, sig_interval_end = 12348, batch = 2)
# 
# # EXPERIMENT GRPS b2
# subset_sig_region_gds(8, sig_interval_start = 7361, sig_interval_end = 7377, batch = 2)
# subset_sig_region_gds(10, sig_interval_start = 8402, sig_interval_end = 8412, batch = 2)
# subset_sig_region_gds(13, sig_interval_start = 10635, sig_interval_end = 10639, batch = 2)
# subset_sig_region_gds(16, sig_interval_start = 12330, sig_interval_end = 12348, batch = 2)
# 

# LIST ALL GDS FILES needed for adding QC labels 
gds_files_b1 <- list.files(
  path = "../Data/", 
  pattern = "^cardiac_cohorts_SOL_dragen\\.chr[0-9]+_batch1+_admixmap_sigregion.*\\.gds$", 
  full.names = TRUE
)

gds_files_b2 <- list.files(
  path = "../Data/", 
  pattern = "^cardiac_cohorts_SOL_dragen\\.chr[0-9]+_batch2+_admixmap_sigregion.*\\.gds$", 
  full.names = TRUE
)

gds_files_b1 <- list.files(
  path = "../Data/", 
  pattern = "^cardiac_cohorts_SOL_dragen\\.chr[0-9]+_batch1+_scaled_up_admixmap_sigregion.*\\.gds$", 
  full.names = TRUE
)

gds_files_b2 <- list.files(
  path = "../Data/", 
  pattern = "^cardiac_cohorts_SOL_dragen\\.chr[0-9]+_batch2+_scaled_up_admixmap_sigregion.*\\.gds$", 
  full.names = TRUE
)


addQCLabel <- function(gds.path){

  g_sig_region <- seqOpen(gds.path, readonly = FALSE)
  
  # we need to add QC labels with all "PASS" manually to make STAAR accept the gds input 
  position <- as.integer(seqGetData(g_sig_region, "position"))
  chr <- unique(as.integer(seqGetData(g_sig_region, "chromosome")))
  length(position) 
  Anno.folder <- index.gdsn(g_sig_region, "annotation/info")
  add.gdsn(Anno.folder, "QC_label", val=factor(rep("PASS", length(position))), replace = TRUE,
           compress="LZMA_ra", closezip=TRUE)
  seqClose(g_sig_region)
  }
    
  
for(file in c(gds_files_b1, gds_files_b2)){
  print(file)
  addQCLabel(file)
}


### DB split information 
file_DBsplit <- "../Data/FAVORdatabase_chrsplit.csv"
### Targeted GDS
dir_geno <- "../Data/"
### output
output_path <- "../Data/STAAR_prep/" 
#output_path <- "../Data/STAAR_prep/scale_up/" 
      
  # then we can generate annotated gds file (aGDS) using FAVORannotator
  
  #----------------------------------------------------#
  # STEP 1: Generate the variants list to be annotated #
  #----------------------------------------------------#
  
  # code adapted from github tutorial 
  ###########################################################################
  #           Main Function 
  ###########################################################################
annotateGDS <- function(chr, gds.path, batch){
  ### chromosome number
  ## read info
  DB_info <- read.csv(file_DBsplit,header=TRUE)
  DB_info <- DB_info[DB_info$Chr==chr,]
  
  ## open GDS
  # gds.path <- paste0(dir_geno,gds_file_name_1)
  genofile <- seqOpen(gds.path)
  
  CHR <- as.numeric(seqGetData(genofile, "chromosome"))
  position <- as.integer(seqGetData(genofile, "position"))
  REF <- as.character(seqGetData(genofile, "$ref"))
  ALT <- as.character(seqGetData(genofile, "$alt"))
  
  VarInfo_genome <- paste0(CHR,"-",position,"-",REF,"-",ALT)
  
  seqClose(genofile)
  
  # for multi-allelic variant, I used only one alternative allele to look for annotation
  # found redundancy and matching issues if splitting into two separate entries
  VarInfo_genome_splitted <- unlist(lapply(VarInfo_genome, function(x) {
    if (grepl(",", x)) {
      # Extract the base and the variants
      parts <- strsplit(x, "-")[[1]]
      base <- paste(parts[1:3], collapse = "-")  # First three components as base
      variants <- unlist(strsplit(parts[4], ",")) # Variants split by comma
      sapply(variants[1], function(v) paste(base, v, sep = "-"))
    } else {
      x  # Leave unchanged if no comma
    }
  }))
  names(VarInfo_genome_splitted) <- NULL
  
  ## Generate VarInfo
  for(kk in 1:dim(DB_info)[1])
  {
    print(kk)
    
    VarInfo <- VarInfo_genome_splitted[(position>=DB_info$Start_Pos[kk])&(position<=DB_info$End_Pos[kk])]
    VarInfo <- data.frame(VarInfo)
    
    filename_extract <- sub("^.*cardiac_cohorts_SOL_dragen\\.(.*)\\.gds$", "\\1", gds.path)
    write.csv(VarInfo,
              paste0(output_path,"VarInfo_chr",chr,"_",kk, "_batch_", batch, ".csv"),quote=FALSE,row.names = FALSE)
    gc()
    
  }
  
}

chr_list <- rep(c(11, 11, 12, 12, 16, 2, 5, 8),2)
batch_all <- c(rep(1, length(gds_files_b1)), rep(2, length(gds_files_b2)))
gds_files_all <- c(gds_files_b1, gds_files_b2)

for(i in 1:length(chr_list)){
  file <- gds_files_all[i]
  print(file)
  batch <- batch_all[i]
  if(i %in% c(2,10)){
    batch <- paste0(batch, "_", "3beta")
  } else if(i %in% c(4,12)){
    batch <- paste0(batch, "_", "ethylmalonate")
  }
  annotateGDS(gds.path = file, chr_list[i], batch = batch)
}



#-----------------------------------------------------------------------------#
# STEP 2: Annotate the variants using the FAVOR database through xsv software #
#-----------------------------------------------------------------------------#

# rm(list=ls())
gc()

##########################################################################
#           Input
##########################################################################
### DB split information 
file_DBsplit <- "../Data/FAVORdatabase_chrsplit.csv"

### xsv directory
xsv <- "C:\\Users\\xchen15\\Downloads/xsv-0.13.0-x86_64-pc-windows-gnu/xsv.exe"

### output
output_path <- "../Data/STAAR_prep/"
#output_path <- "../Data/STAAR_prep/scale_up/"



### anno channel (subset)
anno_colnum <- c(1,8:12,15,16,19,23,25:36)


annotationFavorDB <- function(chr, kk_start, kk_end, met = ""){
  # for both batches 
  
  ###########################################################################
  #           Main Function 
  ###########################################################################
  
  ### chromosome number
  ## annotate (seperate)
  DB_info <- read.csv(file_DBsplit,header=TRUE)
  chr_splitnum <- sum(DB_info$Chr==chr)
  ### DB file
  DB_path <- paste0("C:\\Users\\xchen15\\Downloads/chr", chr,
                    "/n/holystore01/LABS/xlin/Lab/xihao_zilin/FAVORDB")
  
  
  for(kk in kk_start:kk_end)
  {
    print(paste0("kk_num is: ", kk))
    system(paste0(xsv," join --left VarInfo ",
                  output_path,"VarInfo_chr",
                  chr,"_",kk, "_batch_1", met, ".csv variant_vcf ",DB_path,
                  "/chr",chr,"_",kk,".csv -o ",output_path,
                  "Anno_chr",chr,"_",kk,".csv"))
    
  }
  # 
  # if(!kk_start == kk_end){
  #   merge_command <- paste0(xsv," cat rows ", paste0(output_path,
  #                                                    "Anno_chr",chr,"_",kk_start,".csv"))
  #   merge_command <- paste0(merge_command, " ", paste0(output_path,
  #                           "Anno_chr",chr,"_",kk_end,".csv"))
  #   
  #   merge_command <- paste0(merge_command," -o ",output_path, "Anno_chr",
  #                           chr,"_", kk_start, ".csv")
  #   print(merge_command)
  #   system(merge_command)
  # }
  
  ## subset
  anno_colnum_xsv <- c()
  for(kk in 1:(length(anno_colnum)-1))
  {
    print(kk)
    anno_colnum_xsv <- paste0(anno_colnum_xsv,anno_colnum[kk],",")
  }
  anno_colnum_xsv <- paste0(anno_colnum_xsv,anno_colnum[length(anno_colnum)])
  
  system(paste0(xsv," select ",anno_colnum_xsv," ",output_path,
                "Anno_chr",chr,"_",kk_start,
                ".csv -o ",output_path,"Anno_chr_",chr, met,
                "_STAARpipeline.csv"))
  }
  
# annotationFavorDB(2, 5, 5)
# annotationFavorDB(5, 5, 5)
# annotationFavorDB(6, 11, 11)
annotationFavorDB(8, 9, 9)
# annotationFavorDB(10, 4, 4)
# annotationFavorDB(13, 5, 5)
# annotationFavorDB(16, 5, 5)

# scaled up 
annotationFavorDB(11, 4, 4)
annotationFavorDB(11, 7, 7, met = "_3beta")
annotationFavorDB(2, 5, 5)
annotationFavorDB(16, 2, 2)
annotationFavorDB(5, 3, 3)
annotationFavorDB(8, 2, 2)
annotationFavorDB(12, 1, 1)
annotationFavorDB(12, 8, 8, "_ethylmalonate")

#------------------------------------------------#
# STEP 3: Generate the annotated GDS (aGDS) file #
#------------------------------------------------#


##########################################################################
#           Input
##########################################################################

### load required package
library(readr)

### annotation file (output of Annotate.R)
dir_anno <- "../Data/STAAR_prep/"
#dir_anno <- "../Data/STAAR_prep/scale_up/"
anno_file_name_1 <- "Anno_chr_"
#anno_file_name_1 <- "Anno_chr_scale_up_"
anno_file_name_2 <- "_STAARpipeline.csv"
### input array id from batch file

addFunctionAnnotationGDS <- function(chr, met = ""){
  
  ###########################################################################
  #           Main Function 
  ###########################################################################
  
  
  ### read annotation data
  FunctionalAnnotation <- read_csv(paste0(dir_anno,anno_file_name_1,chr, met, anno_file_name_2),
                                   col_types=list(col_character(),col_double(),col_double(),col_double(),col_double(),
                                                  col_double(),col_double(),col_double(),col_double(),col_double(),
                                                  col_character(),col_character(),col_character(),col_double(),col_character(),
                                                  col_character(),col_character(),col_character(),col_character(),col_double(),
                                                  col_double(),col_character()))
  
  dim(FunctionalAnnotation) # 22 columns as always
  
  ## rename colnames
  colnames(FunctionalAnnotation)[2] <- "apc_conservation"
  colnames(FunctionalAnnotation)[7] <- "apc_local_nucleotide_diversity"
  colnames(FunctionalAnnotation)[9] <- "apc_protein_function"
  
  ## open GDS
  for( batch in 1:2){
    gds.path <- paste0("../Data/cardiac_cohorts_SOL_dragen.chr", chr, "_batch",
                       batch ,"_admixmap_sigregion", met, ".gds")
    # gds.path <- paste0("../Data/cardiac_cohorts_SOL_dragen.chr", chr, "_batch",
    #                    batch ,"_scaled_up_admixmap_sigregion", met, ".gds")
    genofile <- seqOpen(gds.path, readonly = FALSE)
    
    Anno.folder <- index.gdsn(genofile, "annotation/info")
    add.gdsn(Anno.folder, "FunctionalAnnotation", replace = TRUE,
             val=FunctionalAnnotation, compress="LZMA_ra", closezip=TRUE)
    
    seqClose(genofile)
  }
  
}

addFunctionAnnotationGDS(2)
addFunctionAnnotationGDS(5)
addFunctionAnnotationGDS(6)
addFunctionAnnotationGDS(8)
addFunctionAnnotationGDS(10)
addFunctionAnnotationGDS(13)
addFunctionAnnotationGDS(16)
addFunctionAnnotationGDS(2)
addFunctionAnnotationGDS(11)
addFunctionAnnotationGDS(11, "_3beta")
addFunctionAnnotationGDS(5)
addFunctionAnnotationGDS(8)
addFunctionAnnotationGDS(12)
addFunctionAnnotationGDS(12, "_ethylmalonate")
addFunctionAnnotationGDS(16)
