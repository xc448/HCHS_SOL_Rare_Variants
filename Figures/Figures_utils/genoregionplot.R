setwd("/Volumes/Sofer Lab/HCHS_SOL/Projects/2024_rare_variants/Code")

library(RIdeogram)
library(tidyverse)


data(human_karyotype, package="RIdeogram")
data(gene_density, package="RIdeogram")

chr <- c(2, 5, 8, 10, 11, 12, 13, 16)


selected_regions <- read.csv("../Data/all_selected_regions.csv")
selected_regions <- selected_regions |>
  select(Chr, start_pos, end_pos, Driving.Ancestry)
colnames(selected_regions) <- colnames(gene_density)
selected_regions <- selected_regions |>
  mutate(Value = if_else(Value == "AFR", "1", Value)) |>
  mutate(Value = if_else(Value == "NAM", "2", Value))     
selected_regions$Value <- as.numeric(selected_regions$Value)  


ideogram(karyotype = human_karyotype[chr,], overlaid = selected_regions)
convertSVG("chromosome.svg", device = "png")