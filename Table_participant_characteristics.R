library(tidyverse)
library(table1)
library(htmltools)

# load the original covaraite data 
setwd("/Volumes/Sofer Lab/HCHS_SOL/")
pheno_file <- "./Datasets/ms968_covariates_20200220.csv"
dat <- read.csv(pheno_file)

bmi_file <- "./Datasets/ms691_covariates_20191125.csv"
dat1 <- read.csv(bmi_file)
dat1 <- dat1 |> dplyr::select(ID, BMI)

dat <- dat |> inner_join(dat1,  by = ("ID" = "ID"))

# Load the covariates for batch 1(x1 from the association R script used to fit the model)
covariate_b1 <- readRDS("./Projects/2024_rare_variants/Data/STAAR_model_cov_batch1.RDS")
flare_intersect_id_b1 <- readRDS("./Projects/2024_rare_variants/Data/SoL_intersected_ID_b1.RDS")
covariate_b1 <- covariate_b1 |>
  filter(SUBJECT_ID %in% flare_intersect_id_b1)

# Load the covariates for batch 2(x2 from the association R script used to fit the model)
covariate_b2 <- readRDS("./Projects/2024_rare_variants/Data/STAAR_model_cov_batch2.RDS")
flare_intersect_id_b2 <- readRDS("./Projects/2024_rare_variants/Data/SoL_intersected_ID_b2.RDS")
covariate_b2 <- covariate_b2 |>
  filter(SUBJECT_ID %in% flare_intersect_id_b2) |>
  mutate(GFRSCYS = if_else(is.na(GFRSCYS), mean(GFRSCYS, na.rm = TRUE), GFRSCYS))

# Participants from the two batches are not overlapping == 0
# sum(covariate_b1$SUBJECT_ID %in% covariate_b2$SUBJECT_ID)

# Match the indices and join other covariates that need to be reported
# Add a new column indicating which batch it is
# Batch1
col_report <- dat |> 
  mutate(ID = as.character(ID)) |>
  dplyr::select(ID, HYPERTENSION, DIABETES2_INDICATOR, DIABETES2, BMI) 

# covariate_b1 <- covariate_b1 |>
#   select(-c(PC6, PC7, PC8, PC9, 
#             PC10, PC11))

covariate_b1 <- covariate_b1 |> 
  inner_join(col_report, by = ("ID" = "ID")) |>
  mutate(batch = rep("Discovery Batch", nrow(covariate_b1)))

# Batch2
covariate_b2 <- covariate_b2 |> 
  inner_join(col_report, by = ("ID" = "ID")) |>
  mutate(batch = rep("Replication Batch", nrow(covariate_b2)))

# rbind the two dfs to obtain all covariates
covariate_total <- rbind(covariate_b1, covariate_b2) |>
  dplyr::select(AGE, BMI, GFRSCYS, GENDER, CENTER, BKGRD1_C7, HYPERTENSION, DIABETES2_INDICATOR, batch) |>
  
  mutate(BKGRD1_C7 = if_else(BKGRD1_C7 == 0, "Dominican", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == 1, "Central American", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == 2, "Cuban", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == 3, "Mexican", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == 4, "Puerto Rican", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == 5, "South American", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == 6, "Other", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == "Q", "Other", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == "", "Other", BKGRD1_C7),) |>

  mutate(HYPERTENSION = if_else(HYPERTENSION == 0, "No", "Yes"),
         DIABETES2_INDICATOR = if_else(DIABETES2_INDICATOR == 0, "No", "Yes"),
         GENDER = recode(GENDER, F = "Female", "M" = "Male"),
         GENDER = factor(GENDER),
         CENTER = recode(CENTER, "B" = "Bronx", "C" = "Chicago",
                         "M" = "Miami",  "S" = "San Diego"),
         CENTER = factor(CENTER),
         BKGRD1_C7 = factor(BKGRD1_C7),
         HYPERTENSION = factor(HYPERTENSION),
         DIABETES2_INDICATOR = factor(DIABETES2_INDICATOR))
# 5768 individuals in total

colnames(covariate_total) <-  c("Age (years)", "BMI", "eGFR", "Gender", 
                                "Recruitment Center",  "Self-reported Background",
                                "Hypertension", "Diabetes", "batch")


covariate_total$`Self-reported Background` <- fct_relevel(covariate_total$`Self-reported Background`,
                                                           "Mexican","South American",
                                                           "Cuban", "Puerto Rican", "Dominican", 
                                                           "Central American", "Other")
### Start to reate the table for participant characteristics ####
# Create participant feature table 
units(covariate_total$BMI)   <- HTML("kg/m<sup>2</sup>")
units(covariate_total$eGFR)   <- HTML("mL/min/1.73m<sup>2</sup>")

label_function <- function(x, ...) {
  s <- table1::label.default(x, ...)
  ifelse(x == "No", "", s)
}

# More display setting for what to render
render.mean_sd <- function(x, name, ...) {
  
  # If the variable is categorical but not "No", render as usual
  if (!is.numeric(x)){
    y <- render.categorical.default(x)
    if (is.factor(x) && all(c("No", "Yes") %in% levels(x))) return(y[3]) else return(y)
  } 
  
  what <- switch(name,
                 `Age (years)` = "Mean (SD)",
                 BMI  = "Mean (SD)",
                 eGFR = "Mean (SD)")
  parse.abbrev.render.code(c("", what))(x)
  
}

# run table1() to generate tables based on mean and sd for continuous variables 
# and the total number of individuals in categorical variables.
table1(~ `Age (years)` + BMI + eGFR + Gender + `Recruitment Center`+ `Self-reported Background`
       + Hypertension + `Diabetes`| batch, data = covariate_total,
       render = render.mean_sd) 

