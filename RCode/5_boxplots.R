#
# 5_boxplots
# 2024-04-08
# Description: This script checks the distribution of FPKM scores by creating a 
# boxplot of log2-transformed FPKM values for the TNC gene across different 
# NCN_PAM50 subtypes, excluding unclassified samples. It also generates a boxplot 
# of PSI scores for specific alternative splicing events of the TNC gene, 
# highlighting splicing variations across NCN_PAM50 subtypes. Both visualizations 
# help assess expression and splicing variation in the TNC gene.
#
#
#

# Load packages 
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

# Sanity check- boxplot of FPKM scores 
fpkm_melted <- left_join(fpkm_melted, ER_PR_HER2_status, by = "Sample")
fpkm_melted <- fpkm_melted %>%
  select(-ER.y, -PR.y, -HER2.y, -NCN_PAM50.y, -ER.x, -PR.x, -HER2.x) 
colnames(fpkm_melted)[colnames(fpkm_melted) == "NCN_PAM50.x"] <- "NCN_PAM50"

fpkm_melted %>%
  filter(Location == "TNC") %>%
  filter(NCN_PAM50 != "Unclassified") %>%
  mutate(log2_FPKM = log2(FPKM)) %>%
  ggplot(aes(x = Location, y = log2_FPKM, color = NCN_PAM50)) +
  geom_boxplot() +
  labs(title = "Expression of TNC Across Subtypes") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())

# Boxplot of PSI Scores from TNC gene 
boxplot_scores <- merged_PSI_ER_PR_HER2 %>%
  filter(Location == "chr9_-_115057152_115057425" |
           Location == "chr9_-_115059729_115060002" |
           Location == "chr9_-_115046409_115046682" |
           Location == "chr9_-_115063795_115064068" |
           Location == "chr9_-_115064646_115064919" |
           Location == "chr9_-_115062916_115063189") %>%
  filter(NCN_PAM50 != "unclassified")

boxplot_scores$NCN_PAM50 <- sub("Normal", "Normal-like", boxplot_scores$NCN_PAM50)

boxplot_scores %>%
  ggplot(aes(x = Location, y = PSI_Score, color = NCN_PAM50)) +
  geom_boxplot() +
  labs(title = "Splicing Variation in TNC", x = "Alternative Splicing Event") + 
  scale_x_discrete(guide = guide_axis(angle = 45))
  



