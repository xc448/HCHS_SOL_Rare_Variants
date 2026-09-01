library(tidyverse)
library(data.table)

results_coding_all <- c()
gene_centric_coding_jobs_num <- 1
gene_centric_results_name <- "Genecentric_coding_uncond"

getnumtests <- function(metab){
  input_path <- paste0("../Data/STAAR_exp_groups/",
                        metab, "/")
  
  file <- list.files(path = input_path, pattern = paste0(gene_centric_results_name, ".Rdata"), 
                    full.names = TRUE)
  print(file)
  results_coding <- get(load(file))
  
  for(kk in 1:length(results_coding))
  {
    results <- as.data.frame(results_coding[[kk]])
    if(ncol(results) > 91){
      results <- results[,1:91]
    }
    results_coding_genome <- rbind(results_coding_genome, results)
  }
  return(results_coding_genome)
}

metabs <- c("x100009332", "x100006370", "x100001208", "x2054", "x100006264",
            "x100001266", "x1114", "x192", "x100000007", "x100001208",
            "x100004046", "x1021")
for(i in metabs){
  print(i)
  res <- getnumtests(i)
  results_coding_all <- rbind(results_coding_all, res)
}
