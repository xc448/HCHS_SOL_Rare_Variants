library(haven)
library(readxl)
library(tidyverse)

covariate_b1 <- readRDS("/Volumes/Sofer Lab/HCHS_SOL/Projects/2024_rare_variants/Data/metab_covariates_batch1.RDS") 
covariate_b2 <- readRDS("/Volumes/Sofer Lab/HCHS_SOL/Projects/2024_rare_variants/Data/metab_covariates_batch2.RDS") 
## metabolomics data file:
metab_file_b1 <- "/Volumes/Sofer Lab/HCHS_SOL/Metabolomics/SOL_with_Batch2/2batch_combined_data_V1_only.xlsx"
# file in raw peaks
metab_file_b1_raw <- "/Volumes/Sofer Lab/HCHS_SOL/Metabolomics/SOL_Batch1/SOL_metabolomics_std_10202017.csv"
metab_file_b2 <- "/Volumes/Sofer Lab/HCHS_SOL/Metabolomics/SOL_with_Batch2/batch2_data.xlsx"


batch1_metab_sheet_name <- "data"
batch2_metab_sheet_name <- "batch2_batchnormalized_v1"
batch1_sample_info_sheet_name <- "sample.info"
batch2_sample_info_sheet_name <- "sample.info_batch2_v1"
metabolites_info_sheet_name <- "metabolites.info"


metab_vals_batch1 <- read_excel(metab_file_b1, sheet = batch1_metab_sheet_name) # 4004 samples
metab_vals_batch1_raw <- read_csv(metab_file_b1_raw) # 3978 rows
metab_vals_batch2 <- read_excel(metab_file_b2, sheet = batch2_metab_sheet_name) # 2368 samples


sample_info_batch1 <- read_excel(metab_file_b1, sheet = batch1_sample_info_sheet_name)
sample_info_batch2 <- read_excel(metab_file_b2, sheet = batch2_sample_info_sheet_name)


# add sol Id and lab id for batch 1 
metab_vals_batch1 <- merge(metab_vals_batch1, 
                           sample_info_batch1[,c("PARENT_SAMPLE_NAME", "SOL_ID", "LABID", "VISIT",
                                                 "BATCH")], 
                           by = "PARENT_SAMPLE_NAME")


metab_vals_batch2 <- merge(metab_vals_batch2, 
                           sample_info_batch2[,c("PARENT_SAMPLE_NAME", "SOL_ID")], 
                           by = "PARENT_SAMPLE_NAME") #2368 rows 


# add the metabolites info of propyl4hydroxybenzoatesulfate_std from the raw 
# csv file to metab_vals_batch1
metab_vals_batch1_raw$LAB_ID <- as.character(metab_vals_batch1_raw$LAB_ID)

metab_names_b1_raw <- colnames(metab_vals_batch1_raw)

neg_control_names <- c("N2acetyllysine_std", "betaine_std", 
                       "Nacetylglucosaminylasparagine_std", "cysteinylglycine_std" )
neg_control_names %in% metab_names_b1_raw  # should all be true

exp_group_ids <- c(1021, 100000007, 100004046, 1224)
neg_control_ids <- c(100001721, 799, 1215, 278)
metab_vals_batch1 <- metab_vals_batch1[, colnames(metab_vals_batch1) %in% c(as.character(exp_group_ids),
                                                                            "LABID", "SOL_ID", "BATCH")]
exp_group_names <- c("oxoproline_std", "carnitine_std",
                     "Nacetylcarnosine_std")

metab_vals_batch1 <- metab_vals_batch1 |>
  left_join(metab_vals_batch1_raw[,c("LAB_ID", neg_control_names, exp_group_names
  )],
  by = c("LABID" = "LAB_ID"))

# add sol Id and lab id for batch 2
metab_vals_batch2 <- merge(metab_vals_batch2, 
                           sample_info_batch2[,c("PARENT_SAMPLE_NAME", "SOL_ID")], 
                           by = "PARENT_SAMPLE_NAME")

set.seed(1997)
# choose at random samples that are from the same individual (to have one sample per individual) 
random_index_b1 <- data.frame(index = 1:nrow(metab_vals_batch1),
                              SOL_ID = metab_vals_batch1$SOL_ID)
for (id in unique(random_index_b1$SOL_ID)){
  row_inds <- which(random_index_b1$SOL_ID == id)
  if (length(row_inds) == 1) next
  selected_ind <- sample(row_inds, 1)
  random_index_b1 <- random_index_b1[-setdiff(row_inds, selected_ind),]
}
metab_vals_batch1 <- metab_vals_batch1[random_index_b1$index,] 

# choose at random samples that are from the same individual (to have one sample per individual)
random_index_b2 <- data.frame(index = 1:nrow(metab_vals_batch2 ),
                              SOL_ID = metab_vals_batch2$SOL_ID)
for (id in unique(random_index_b2$SOL_ID)){
  row_inds <- which(random_index_b2$SOL_ID == id)
  if (length(row_inds) == 1) next
  selected_ind <- sample(row_inds, 1)
  random_index_b2 <- random_index_b2[-setdiff(row_inds, selected_ind),]
}
metab_vals_batch2 <- metab_vals_batch2[random_index_b2$index,] # 2330

metab_vals_batch2 <- metab_vals_batch2[-which(is.element(metab_vals_batch2$SOL_ID, metab_vals_batch1$SOL_ID)),]
# 2178 

metab_vals_batch1 <- metab_vals_batch1 |> filter(BATCH == "B01") # 3978 

# rank-normalize the peaks and center it around the 1.0 to avoid non-positive values
metab_vals_batch1_metab <- metab_vals_batch1 |>
  mutate(N2acetyllysine_std = if_else(is.na(N2acetyllysine_std), 
                                      min(N2acetyllysine_std, na.rm= TRUE),
                                      N2acetyllysine_std)) |>
  mutate(ranks = rank(N2acetyllysine_std, ties.method = "random"),
         int_norm = qnorm((ranks - 0.5) / max(ranks, na.rm = TRUE)),
         x100001721 = (int_norm - min(int_norm, na.rm = TRUE)) / 
           (max(int_norm, na.rm = TRUE) - min(int_norm, na.rm = TRUE))) |>
  mutate(x100001721 = x100001721 - median(x100001721) + 1) |>
  select(-c(ranks, int_norm)) |>
  mutate(betaine_std = if_else(is.na(betaine_std), 
                               min(betaine_std, na.rm= TRUE),
                               betaine_std)) |>  
  mutate(ranks = rank(betaine_std, ties.method = "random"),
         int_norm = qnorm((ranks - 0.5) / max(ranks, na.rm = TRUE)),
         x799 = (int_norm - min(int_norm, na.rm = TRUE)) / 
           (max(int_norm, na.rm = TRUE) - min(int_norm, na.rm = TRUE))) |>
  mutate(x799 = x799 - median(x799) + 1) |> 
  select(-c(ranks, int_norm)) |>
  mutate(Nacetylglucosaminylasparagine_std = if_else(is.na(Nacetylglucosaminylasparagine_std), 
                                                     min(Nacetylglucosaminylasparagine_std, na.rm= TRUE),
                                                     Nacetylglucosaminylasparagine_std)) |>
  mutate(ranks = rank(Nacetylglucosaminylasparagine_std, ties.method = "random"),
         int_norm = qnorm((ranks - 0.5) / max(ranks, na.rm = TRUE)),
         x1215 = (int_norm - min(int_norm, na.rm = TRUE)) / 
           (max(int_norm, na.rm = TRUE) - min(int_norm, na.rm = TRUE))) |>
  mutate(x1215 = x1215 / median(x1215)) |> 
  select(-c(ranks, int_norm)) |>
  mutate(cysteinylglycine_std = if_else(is.na(cysteinylglycine_std), 
                                        min(cysteinylglycine_std, na.rm= TRUE),
                                        cysteinylglycine_std)) |>
  mutate(ranks = rank(cysteinylglycine_std, ties.method = "random"),
         int_norm = qnorm((ranks - 0.5) / max(ranks, na.rm = TRUE)),
         x278= (int_norm - min(int_norm, na.rm = TRUE)) / 
           (max(int_norm, na.rm = TRUE) - min(int_norm, na.rm = TRUE))) |>
  mutate(x278 = x278 - median(x278) + 1) |>
  select(-c(ranks, int_norm))


hist(metab_vals_batch1_metab$x100001721) # is normalized and cetered around 1. 
hist(metab_vals_batch1_metab$x799) 
hist(metab_vals_batch1_metab$x1215)
hist(metab_vals_batch1_metab$x278)

exp_group_ids <- c(1021, 100000007, 100004046, 1224)
metab_vals_batch1_metab <- metab_vals_batch1_metab |>
  mutate(oxoproline_std = if_else(is.na(oxoproline_std), 
                                  min(oxoproline_std, na.rm= TRUE),
                                  oxoproline_std)) |>
  mutate(ranks = rank(oxoproline_std, ties.method = "random"),
         int_norm = qnorm((ranks - 0.5) / max(ranks, na.rm = TRUE)),
         x1021 = (int_norm - min(int_norm, na.rm = TRUE)) / 
           (max(int_norm, na.rm = TRUE) - min(int_norm, na.rm = TRUE))) |>
  mutate(x1021 = x1021 - median(x1021) + 1) |>
  select(-c(ranks, int_norm)) |>
  mutate(carnitine_std = if_else(is.na(carnitine_std), 
                                 min(carnitine_std, na.rm= TRUE),
                                 carnitine_std)) |>  
  mutate(ranks = rank(carnitine_std, ties.method = "random"),
         int_norm = qnorm((ranks - 0.5) / max(ranks, na.rm = TRUE)),
         x100000007 = (int_norm - min(int_norm, na.rm = TRUE)) / 
           (max(int_norm, na.rm = TRUE) - min(int_norm, na.rm = TRUE))) |>
  mutate(x100000007 = x100000007 - median(x100000007) + 1) |> 
  select(-c(ranks, int_norm)) |>
  mutate(Nacetylcarnosine_std = if_else(is.na(Nacetylcarnosine_std), 
                                        min(Nacetylcarnosine_std, na.rm= TRUE),
                                        Nacetylcarnosine_std)) |>
  mutate(ranks = rank(Nacetylcarnosine_std, ties.method = "random"),
         int_norm = qnorm((ranks - 0.5) / max(ranks, na.rm = TRUE)),
         x100004046 = (int_norm - min(int_norm, na.rm = TRUE)) / 
           (max(int_norm, na.rm = TRUE) - min(int_norm, na.rm = TRUE))) |>
  mutate(x100004046 = x100004046 - median(x100004046) + 1) |> 
  select(-c(ranks, int_norm)) 

length(intersect(metab_vals_batch1$SOL_ID, metab_vals_batch2$SOL_ID)) #0

metab_selected_batch1 <- metab_vals_batch1_metab[,colnames(metab_vals_batch1_metab) %in% c("SOL_ID", paste0("x", neg_control_ids),
                                                                                           paste0("x", exp_group_ids))] |>
  rename(ID = SOL_ID)

matching_idinfo <- read.table("../Data/20240717_HCHS_SoL_DCC_ID_matching.txt",
                              header = TRUE) # 16515 individuals 

matching_idinfo <- as.data.frame(as.character(matching_idinfo$ID))
colnames(matching_idinfo) <- "ID"
metab_selected_batch1 <- metab_selected_batch1 |> left_join(matching_idinfo, by = c("ID"))
covariate_b1_all_metab <- covariate_b1 |> left_join(metab_selected_batch1, by = c("ID")) 

# covariate_b1_all_metab$x1224 <- covariate_b1_all_metab$`1224`
# covariate_b1_all_metab <- covariate_b1_all_metab[, -2]

saveRDS(covariate_b1_all_metab, "../Data/covariate_b1_all_metab.RDS")

metab_selected_batch2 <- metab_vals_batch2[,colnames(metab_vals_batch2) %in% c("SOL_ID", neg_control_ids, exp_group_ids)] |>
  rename(ID = SOL_ID) |>
  left_join(matching_idinfo, by = c("ID"))

covariate_b2_all_metab <- covariate_b2 |> 
  left_join(metab_selected_batch2, by = c("ID", 1224)) 

colnames(covariate_b2_all_metab)[13:ncol(covariate_b2_all_metab)] <-
  paste0("x", colnames(covariate_b2_all_metab)[13:ncol(covariate_b2_all_metab)])
covariate_b2_all_metab <- covariate_b2_all_metab |>
  mutate('x1224' = `1224`) |>
  select(-`1224`)

cols_to_impute <- colnames(covariate_b2_all_metab)[13:ncol(covariate_b2_all_metab)]

covariate_b2_all_metab <- covariate_b2_all_metab %>%
  mutate(across(all_of(cols_to_impute), ~ ifelse(is.na(.), min(., na.rm = TRUE), .)))

# Note that for batch 2 x100000007 has A LOT OF NAs > 40%
# -- we might not be able to perform replication experiment for that specific metabolite

saveRDS(covariate_b2_all_metab, "../Data/covariate_b2_all_metab.RDS")

