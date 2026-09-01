library(tidyverse)




metab <- c("x100001721", "x799", "x1215", "x278",
  "x100001266", "x1114", "x192", "x1021",
  "x100000007", "x100009332", "x100006370",
  "x100001208", "x2054", "x100004046","x100006264",
  "x1224")

groups <- c(rep("neg_controls",4),
            rep("exp_groups",12))

computePosRange <- function(known_loci){
  return(diff(range(known_loci$POS)))
}

res_range_known_loci_pos <- c()
for(i in 1:length(metab)){
  known_loci <- readRDS(paste0("../Data/STAAR_", groups[i], "/", metab[i], "/",
  "individual_cond_pruned_var.RDS"))
  res_range_known_loci_pos <- c(res_range_known_loci_pos, computePosRange(known_loci))
}


res_range_known_loci_pos <- as.data.frame(res_range_known_loci_pos)
res_range_known_loci_pos$metab <- metab
colnames(res_range_known_loci_pos)[1] <- "var_positon_range"
saveRDS(res_range_known_loci_pos, "../Data/pos_range_known_common_var.RDS")

