library(Biobase)
library(dplyr)
library(ggplot2)
library(GGally)
library(GWASTools)
library(RColorBrewer)
library(forcats)

# load pca results
## get PCs
load("../Data/pca.Rdata")

# without outliers identified as central american/other
load("../Data/pca_no_outliers.RData")

# round 2 without outliers
load("../Data/pca_no_outliers_round2.RData")

# round 2 with different set of SNPs 
load("../Data/pca_no_outliers_round2_diff_snps.RData")

pcs <- as.data.frame(pca$vectors[pca$unrels,])
#pcs <- as.data.frame(pca$vectors)
n <- 32
names(pcs) <- paste0("PC", 1:n)
pcs$sample.id <- row.names(pcs)

# load covariate info
covariates_combined <- readRDS("/Volumes/Sofer Lab/HCHS_SOL/Projects/2024_rare_variants/Data/matched_NWD_SoLID_batch1_2_combined.RDS")

# load matched ID info
rownames(covariates_combined) <- covariates_combined$NWD_ID
# indexing for unrelated individuals 
covariates_combined <- covariates_combined[pcs$sample.id,]

## scree plot
dat <- data.frame(pc=1:n, varprop=pca$varprop)
p <- ggplot(dat, aes(x=factor(pc), y=100*varprop)) +
  geom_point() + theme_bw() +
  xlab("PC") + ylab("Percent of variance accounted for")
ggsave("../Data/pca_scree_no_outliers_round2_diff_snps.pdf", 
       plot=p, width=6, height=6)

## color by group
covariates_combined <- covariates_combined |>
  mutate(BKGRD1_C7 = if_else(BKGRD1_C7 == 0, "Dominican", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == 1, "Central American", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == 2, "Cuban", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == 3, "Mexican", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == 4, "Puerto Rican", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == 5, "South American", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == 6, "Other", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == "Q", "Other", BKGRD1_C7),
         BKGRD1_C7 = if_else(BKGRD1_C7 == "", "Other", BKGRD1_C7),)





group <- as.factor(covariates_combined$BKGRD1_C7)

pcs$group <- group

pcs$group  <- fct_relevel(pcs$group, "Mexican","South American",
                          "Cuban", "Puerto Rican", "Dominican", 
                            "Central American", "Other")

p <- ggplot(pcs, aes(x = PC1, y = PC2, color=group)) + 
  geom_point(alpha=0.5) +
  guides(colour=guide_legend(override.aes=list(alpha=1)))+
  labs(color = "Self-identified Background") +
  scale_color_brewer(palette = "Set1")+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA), # White panel background
    plot.background = element_rect(fill = "white", color = NA)
  ) 
ggsave("../Figures/pca_pc12_no_outliers_round2_diff_snps.pdf", 
       plot=p, width=10, height=6)

npr <- min(6, n)
p <- ggpairs(
  pcs,
  mapping = aes_string(color = "group"), 
  columns = 1:npr, 
  lower = list(continuous = wrap("points", alpha = 0.5)),
  diag = list(continuous = wrap("densityDiag", alpha = 0.8)), # Set alpha for diagonal plots
  upper = list(continuous = "blank")
) +
  scale_color_brewer(palette = "Set1") + 
  scale_fill_brewer(palette = "Set1") + 
  theme_minimal() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) 
  )
png("../Figures/pca_pairs_no_outliers_round2_diff_snps.png", 
    width=8, height=8, units="in", res=150)
print(p)
dev.off()


pc2 <- pcs
names(pc2)[1:ncol(pc2)] <- sub("PC", "", names(pc2)[1:ncol(pc2)])
names(pc2)[32:34] <- 32:34
p <- ggparcoord(pc2, columns=1:32, groupColumn=34, alphaLines=0.5, scale="uniminmax") +
  scale_color_brewer(palette = "Set1")+
  guides(colour=guide_legend(override.aes=list(alpha=1, size=2))) +
  xlab("PC") + ylab("")+
  labs(color = "Self identified backgrounds")
ggsave("../Figures/pca_32_substructure_no_outliers_round2_diff_snps.png", plot=p, width=10, height=5)


p1 <- ggparcoord(pc2, columns=1:6, groupColumn=34, alphaLines=0.5, scale="uniminmax") +
  scale_color_brewer(palette = "Set1")+
  guides(colour=guide_legend(override.aes=list(alpha=1, size=2))) +
  xlab("PC") + ylab("")+
  labs(color = "Self identified backgrounds")
ggsave("../Figures/pca_6_substructure_no_outliers_round2_diff_snps.png", plot=p1, width=10, height=5)


pcs_long <- pcs %>%
  pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "Value")
pcs_long$PC <- as.integer(str_replace(pcs_long$PC, "PC", ""))


# Create the boxplot
ggplot(pcs_long, aes(x = PC, y = Value, fill = factor(group), group = interaction(PC, group))) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
  #geom_violin() + 
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")+
  labs(x = "Principal Components", y = "Value", fill = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Extracting IDs for individuals with outlier PC5 values based on the boxplot 
pcs %>%
  select(PC5, group) %>%
  group_by(group) %>%
  ggplot(aes(y = PC5, x = factor(group)))+
  theme_minimal()+
  geom_boxplot(outlier.size = 0.5)


pc5_outpliers <-  pcs %>%
  select(PC5, group) %>%
  filter(PC5 < -0.05)

saveRDS(pc5_outpliers, "../Data/PC5_outliers.RDS")
pc5_outpliers <- readRDS("../Data/PC5_outliers.RDS")
# mem stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")