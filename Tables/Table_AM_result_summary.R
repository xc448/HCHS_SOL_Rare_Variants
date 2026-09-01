library(tidyverse)
library(data.table)
library(kableExtra)
library(knitr)
library(htmltools)


findmax_p_beta <- function(df){
  max_ind <- which(df$Score.pval == min(df$Score.pval, na.rm = TRUE))[1]
  stats_max <- df[max_ind, ]
  return(stats_max)
}

AM_result_table <- function(metab, batch = "b1", group = "exp_groups", chr, 
         rv_name = NULL){
  res_no_genotype <- readRDS(paste0("../Data/STAAR_", group, "/", metab,
                                    "/admixmap_", batch, "_",  metab, 
                                    "_no_genotype.RDS"))
  
  res_no_genotype <- findmax_p_beta(res_no_genotype)
  # load cv adjusted results
  if(group == "neg_controls"){
    res_cv_adjusted <- readRDS(paste0("../Data/STAAR_", group, "/", metab,
                                      "/admixmap_", batch, "_",  metab, 
                                      "_cv_gwas.RDS"))
  }else{
    res_cv_adjusted <- readRDS(paste0("../Data/STAAR_", group, "/", metab,
                                      "/admixmap_", batch, "_",  metab, 
                                      "_cv_adj.RDS"))
  }
  res_cv_adjusted <- findmax_p_beta(res_cv_adjusted)
  res_null_cv_adjusted_res <- c(res_no_genotype$Est,
    res_no_genotype$Score.pval, res_cv_adjusted$Est,
    res_cv_adjusted$Score.pval)
  
  if(is.null(rv_name)) {
    res_rv_adjusted_res <- rep(0, 4)
  }else{
    res_rv_adjusted <- readRDS(paste0("../Data/STAAR_", group, "/", metab,
                                      "/admixmap_", batch, "_",  metab, "_", rv_name,
                                      "_adj.RDS"))
    res_rv_adjusted <- findmax_p_beta(res_rv_adjusted)
    
    res_rv_cv_adjusted <- readRDS(paste0("../Data/STAAR_", group, "/", metab,
                                         "/admixmap_", batch, "_",  metab, "_",
                                         "cv_adj.RDS"))
    res_rv_cv_adjusted <- findmax_p_beta(res_rv_cv_adjusted)
    
    res_rv_adjusted_res <- c(res_rv_adjusted$Est,
      res_rv_adjusted$Score.pval, res_rv_cv_adjusted$Est,
                             res_rv_cv_adjusted$Score.pval)
  }
  res_tot <- as.data.frame(t(c(res_null_cv_adjusted_res, res_rv_adjusted_res)))
  print(res_tot)
  colnames(res_tot) <- c("P Null", "Est. Null",
                         "P CV adj.", "Est. CV adj.",
                         "P RV adj.", "Est. RV adj.",
                         "P CV+RV adj.", "Est. CV+RV adj.")
  res_tot <- signif(res_tot, 3)
  res_tot[,c(1,3,5,7)] <- round(res_tot[,c(1,3,5,7)], 2)
  return(res_tot)
}

final_table <- c()


final_table <- rbind(final_table, AM_result_table("x100001266", batch = "b1",
                                        chr = 2, rv_name = "rv"))
final_table <- rbind(final_table, AM_result_table("x1114", batch = "b1",
                                                  chr = 5, rv_name = "plof_ds"))  
final_table <- rbind(final_table, AM_result_table("x192", batch = "b1",
                                                  chr = 8))  
final_table <- rbind(final_table, AM_result_table("x1021", batch = "b1",
                                                  chr = 8))  
final_table <- rbind(final_table, AM_result_table("x100000007", batch = "b1",
                                                  chr = 10)) 
final_table <- rbind(final_table, AM_result_table("x100009332", batch = "b1",
                                                  chr = 11, rv_name = "rv")) 
final_table <- rbind(final_table, AM_result_table("x100006370", batch = "b1",
                                                  chr = 11)) 
final_table <- rbind(final_table, AM_result_table("x100001208", batch = "b1",
                                                  chr = 12)) 
final_table <- rbind(final_table, AM_result_table("x2054", batch = "b1",
                                                  chr = 12,  rv_name = "plof_ds")) 
final_table <- rbind(final_table, AM_result_table("x100004046", batch = "b1",
                                                  chr = 13)) 
final_table <- rbind(final_table, AM_result_table("x100006264", batch = "b1",
                                                  chr = 16, rv_name = "rv"))
final_table <- rbind(final_table, AM_result_table("x1224", batch = "b1",
                                                  chr = 16, rv_name = "rv_combined")) 

rownames(final_table) <- c("x100001266", "x1114", "x192", "x1021",
                           "x100000007", "x100009332", "x100006370",
                           "x100001208", "x2054", "x100004046","x100006264",
                           "x1224")

fwrite(final_table, "../Data/AM_test_regions_summary.csv")
# kable(final_table) %>%
#   kable_classic( full_width = F,
#                  html_font = "Arial") %>%
#   row_spec(0, bold = TRUE, align = "center")  %>%
#   column_spec(1,bold = T, width =  "4cm")


