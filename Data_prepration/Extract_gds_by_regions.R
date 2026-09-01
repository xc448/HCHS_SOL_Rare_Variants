# load packages
library(tidyverse)
library(SeqArray)

batch <- as.numeric(commandArgs(TRUE)[1])

subset_sig_region_gds <- function(chr, sig_interval_start, sig_interval_end , batch = batch){
  
  intervals_liftover_nodup <- readRDS("/home/ec2-user/EBS4T/Projects/2024_HCHS_SOL_rare_variants/data/intervals_liftover_nodup.RDS")
  
  # subsetting the significant intervals from the RFMix inference data based on ancestry intervals
  # include all intervals from start to end
  known_sig_interval <- intervals_liftover_nodup %>% filter(snpid %in% sig_interval_start:sig_interval_end)
  g <- seqOpen(paste0("/home/ec2-user/EBS4T/Projects/2024_HCHS_SOL_rare_variants/data/cardiac_cohorts_SOL_dragen.chr", chr, ".subset.gds"),
               readonly = FALSE)
  prov_nwd_id_matched <- readRDS("/home/ec2-user/EBS4T/Projects/2024_HCHS_SOL_rare_variants/data/prov_nwd_id_matched.RDS")
  add.gdsn(g, "sample.id", val=prov_nwd_id_matched, compress="LZMA_ra", 
           closezip=TRUE, replace = TRUE)
  if(batch == 1){
    # load matched sample IDs
    sampleid <- readRDS("/home/ec2-user/EBS4T/Projects/2024_HCHS_SOL_rare_variants/data/matched_NWD_SoLID_batch1_2_combined.RDS") %>%
      filter(batch == "BATCH1") # 3977 ppl 
  } else {
    sampleid <- readRDS("/home/ec2-user/EBS4T/Projects/2024_HCHS_SOL_rare_variants/data/matched_NWD_SoLID_batch1_2_combined.RDS") %>%
      filter(batch == "BATCH2") 
  }
  
  pos_ind <- c()
  position <- as.data.frame(seqGetData(g, "position"))
  pos_ind <- which(position >= min(as.numeric(known_sig_interval$pos_start)) &
                     position <= max(as.numeric(known_sig_interval$pos_end)))
  
  variant_id <- seqGetData(g, "variant.id")
  variant_id <- variant_id[pos_ind]
  
  seqSetFilter(g, sample.id = sampleid$NWD_ID, variant.id = variant_id)
  
  # export the subsetted gds file 
  seqExport(g, 
            out.fn = paste0("/home/ec2-user/EBS4T/Projects/2024_HCHS_SOL_rare_variants/data/cardiac_cohorts_SOL_dragen.chr", 
                            chr, "_batch", batch, "_scaled_up", "_admixmap_sigregion_3beta.gds"))
  seqClose(g)
}
# NEGATIVE CTRLS b1
# subset_sig_region_gds(2, sig_interval_start = 1575, sig_interval_end = 1585, batch = batch)
# subset_sig_region_gds(5, sig_interval_start = 4544, sig_interval_end = 4546, batch = batch)
# subset_sig_region_gds(6, sig_interval_start = 5848, sig_interval_end = 5848, batch = batch)
# subset_sig_region_gds(16, sig_interval_start = 12348, sig_interval_end = 12348, batch = batch)

# EXPERIMENT GRPS b1
subset_sig_region_gds(8, sig_interval_start = 7361, sig_interval_end = 7378, batch = batch)
# subset_sig_region_gds(10, sig_interval_start = 8402, sig_interval_end = 8412, batch = batch)
# subset_sig_region_gds(13, sig_interval_start = 10635, sig_interval_end = 10639, batch = batch)

# 
# subset_sig_region_gds(2, sig_interval_start = 1575, sig_interval_end = 1584, batch = batch)
# subset_sig_region_gds(5, sig_interval_start = 4384, sig_interval_end = 4388, batch = batch)
# subset_sig_region_gds(8, sig_interval_start = 6807, sig_interval_end = 6811, batch = batch)
# subset_sig_region_gds(11, sig_interval_start = 9309, sig_interval_end = 9322, batch = batch)
# subset_sig_region_gds(11, sig_interval_start = 9112, sig_interval_end = 9134, batch = batch)
# subset_sig_region_gds(12, sig_interval_start = 9507, sig_interval_end = 9511, batch = batch)
# subset_sig_region_gds(16, sig_interval_start = 12011, sig_interval_end = 12035, batch = batch)


