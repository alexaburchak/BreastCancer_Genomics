#
# 4_Heatmaps
# 2024-04-08
# Description: This script analyzes differential expression and alternative 
# splicing events between various groups (ER positive vs. negative, PR positive 
# vs. negative, HER2 positive vs. negative, and molecular subtypes) and generates 
# heatmaps for visualizing significant events. It includes steps to filter 
# significant results, compute mean differences, create heatmap matrices, and 
# annotate with clinical data, ultimately producing both simple and complex 
# heatmaps categorized by event type. Additionally, it saves tables of 
# significant results for further analysis. 
#
#

# Load necessary packages  
library(RColorBrewer)
library(ComplexHeatmap)
library(colorRamp2)

## Start by making a simple heatmap of ER pos vs neg: 
{
  # Determine which events/genes show differential expression
  significant_ER <- ER_mw_gene_results %>%
    filter(q_value < 0.05) 
  
  # Extract differentially expressed genes
  sig_ER_genes <- data.frame(unique(significant_ER$Gene_Name))
  sig_ER_genes <- rename(sig_ER_genes, Gene_Name = unique.significant_ER.Gene_Name.)
  
  # Extract differentially expressed AS events
  sig_ER_events <- data.frame(significant_ER$Location)
  sig_ER_events <- rename(sig_ER_events, Location = significant_ER.Location)
  
  # Extract PSI scores for significant events 
  ER_events_for_heatmap <- filtered_PSI_scores[filtered_PSI_scores$Location %in% sig_ER_events$Location,]
  
  # Convert PSI scores to matrix 
  ER_heatmap_matrix <- as.matrix(ER_events_for_heatmap[, 4:3457])
  row.names(ER_heatmap_matrix) <- unlist(ER_events_for_heatmap[,1])
  flip_ER_heatmap_matrix <- t(ER_heatmap_matrix)
  
  # Create heatmap
  heatmap(flip_ER_heatmap_matrix, col= colorRampPalette(brewer.pal(8, "Blues"))(25)) # add scale
}

## Complex heatmap of all significant results (ER pos vs neg): 
{
  # Add annotations 
  ER_PR_HER2_status$NCN_PAM50 <- sub("Her2", "HER2", ER_PR_HER2_status$NCN_PAM50)
  ER_PR_HER2_status$NCN_PAM50 <- sub("Normal", "Normal-like", ER_PR_HER2_status$NCN_PAM50)
  ER_PR_HER2_status$NCN_PAM50 <- sub("unclassified", "Unclassified", ER_PR_HER2_status$NCN_PAM50)
  
  annotationTracks <- HeatmapAnnotation(
    PAM50 = ER_PR_HER2_status$NCN_PAM50,
    ER = ER_PR_HER2_status$ER, 
    PR = ER_PR_HER2_status$PR,
    HER2 = ER_PR_HER2_status$HER2,
    col = list(PAM50 = c("Basal" = "red", "HER2" = "purple", "LumA" = "blue", "LumB" = "turquoise", "Normal-like" = "green", "Unclassified" = "grey"),
               ER = c("Negative" = "grey", "Positive" = "black", "NA" = "white"),
               PR = c("Negative" = "grey", "Positive" = "black", "NA" = "white"),
               HER2 = c("Negative" = "grey", "Positive" = "black", "NA" = "white")
    )
  )
  
  # Build heatmap
  Heatmap(ER_heatmap_matrix, name = "PSI",  
          cluster_rows = FALSE, 
          cluster_columns = TRUE,
          show_row_names = TRUE, 
          row_names_gp = gpar(fontsize = 7), 
          column_names_gp = gpar(fontsize = 6), 
          row_names_side = "left", 
          show_column_names = FALSE, 
          column_title_gp = gpar(fontsize = 8), 
          col = colorRamp2(c(0, .5, 1), hcl_palette = "Blues", rev = TRUE), 
          bottom_annotation = annotationTracks)
  
}

## Complex heatmap with cutoff for largest differences in mean PSI score (ER): 
{
  # Calculate the differences between means
  means_df_ER <- significant_ER %>%
    rowwise() %>%
    mutate(mean_diff_max = abs(ER_positive_mean -
             ER_negative_mean))
  
  # Find events with the largest differences
  largest_diff_ER <- means_df_ER %>% 
    filter(mean_diff_max > 0.1)
  
  # Extract differentially expressed AS events
  diff_sig_ER_events <- data.frame(largest_diff_ER$Location)
  diff_sig_ER_events <- rename(diff_sig_ER_events, Location = largest_diff_ER.Location)
  
  # Extract PSI scores for significant events 
  diff_ER_events_for_heatmap <- filtered_PSI_scores[filtered_PSI_scores$Location %in% diff_sig_ER_events$Location,]
  
  # Convert PSI scores to matrix 
  diff_heatmap_matrix_ER <- as.matrix(diff_ER_events_for_heatmap[, 4:3457])
  row.names(diff_heatmap_matrix_ER) <- unlist(diff_ER_events_for_heatmap[,1])
  
  # Build heatmap
  Heatmap(diff_heatmap_matrix_ER, name = "PSI",  
          cluster_rows = TRUE, 
          cluster_columns = TRUE,
          show_row_names = TRUE, 
          row_names_gp = gpar(fontsize = 5), 
          column_names_gp = gpar(fontsize = 6), 
          row_names_side = "left", 
          show_column_names = FALSE, 
          column_title_gp = gpar(fontsize = 8), 
          col = colorRamp2(c(0, .5, 1), hcl_palette = "Blues", rev = TRUE), 
          bottom_annotation = annotationTracks)
}

## Group by AS type (ER):
{
  # CE
  CE_ER_heatmap_events <- diff_ER_events_for_heatmap %>% # 483 events 
    filter(AS_Type == "CE")
  # Convert PSI scores to matrix 
  CE_ER_heatmap_matrix <- as.matrix(CE_ER_heatmap_events[, 4:3457])
  row.names(CE_ER_heatmap_matrix) <- unlist(CE_ER_heatmap_events[,1])
  
  # Create heatmap
  Heatmap(CE_ER_heatmap_matrix, name = "PSI",  
          cluster_rows = FALSE, 
          cluster_columns = TRUE,
          show_row_names = TRUE, 
          row_names_gp = gpar(fontsize = 8), 
          column_names_gp = gpar(fontsize = 6), 
          row_names_side = "left", 
          show_column_names = FALSE, 
          column_title_gp = gpar(fontsize = 8), 
          col = colorRamp2(c(0, .5, 1), hcl_palette = "Blues", rev = TRUE), 
          bottom_annotation = annotationTracks)
  
  #AA
  AA_ER_heatmap_events <- diff_ER_events_for_heatmap %>% # 91 events 
    filter(AS_Type == "AA") 
  # Convert PSI scores to matrix 
  AA_ER_heatmap_matrix <- as.matrix(AA_ER_heatmap_events[, 4:3457])
  row.names(AA_ER_heatmap_matrix) <- unlist(AA_ER_heatmap_events[,1])
 
  # Create heatmap
  Heatmap(AA_ER_heatmap_matrix, name = "PSI",  
          cluster_rows = TRUE, 
          cluster_columns = TRUE,
          show_row_names = TRUE, 
          row_names_gp = gpar(fontsize = 8), 
          column_names_gp = gpar(fontsize = 6), 
          row_names_side = "left", 
          show_column_names = FALSE, 
          column_title_gp = gpar(fontsize = 8), 
          col = colorRamp2(c(0, .5, 1), hcl_palette = "Blues", rev = TRUE), 
          bottom_annotation = annotationTracks)
  
  #AD
  AD_ER_heatmap_events <- diff_ER_events_for_heatmap %>% # 115 events 
    filter(AS_Type == "AD") 
  # Convert PSI scores to matrix 
  AD_ER_heatmap_matrix <- as.matrix(AD_ER_heatmap_events[, 4:3457])
  row.names(AD_ER_heatmap_matrix) <- unlist(AD_ER_heatmap_events[,1])
  
  # Create heatmap
  Heatmap(AD_ER_heatmap_matrix, name = "PSI",  
          cluster_rows = FALSE, 
          cluster_columns = TRUE,
          show_row_names = TRUE, 
          row_names_gp = gpar(fontsize = 8), 
          column_names_gp = gpar(fontsize = 6), 
          row_names_side = "left", 
          show_column_names = FALSE, 
          column_title_gp = gpar(fontsize = 8), 
          col = colorRamp2(c(0, .5, 1), hcl_palette = "Blues", rev = TRUE), 
          bottom_annotation = annotationTracks)
  
  #IR
  IR_ER_heatmap_events <- diff_ER_events_for_heatmap %>% # 104 events 
    filter(AS_Type == "IR") 
  # Convert PSI scores to matrix 
  IR_ER_heatmap_matrix <- as.matrix(IR_ER_heatmap_events[, 4:3457])
  row.names(IR_ER_heatmap_matrix) <- unlist(IR_ER_heatmap_events[,1])
  
  # Create heatmap
  Heatmap(IR_ER_heatmap_matrix, name = "PSI",  
          cluster_rows = FALSE, 
          cluster_columns = TRUE,
          show_row_names = TRUE, 
          row_names_gp = gpar(fontsize = 8), 
          column_names_gp = gpar(fontsize = 6), 
          row_names_side = "left", 
          show_column_names = FALSE, 
          column_title_gp = gpar(fontsize = 8), 
          col = colorRamp2(c(0, .5, 1), hcl_palette = "Blues", rev = TRUE), 
          bottom_annotation = annotationTracks)
}

## Complex heatmap of all significant results (PR Pos vs Neg):
{
  # Determine which events/genes show differential expression
  significant_PR <- PR_mw_gene_results %>% 
    filter(q_value < 0.05) #1150 events out of 5556
  
  # Calculate the differences between means
  means_df_PR <- significant_PR %>%
    rowwise() %>%
    mutate(mean_diff_max = abs(PR_positive_mean -
                                 PR_negative_mean))
  
  # Find events with the largest differences
  largest_diff_PR <- means_df_PR %>% 
    filter(mean_diff_max > 0.1) # 40 events 
  
  # Extract differentially expressed AS events
  diff_sig_PR_events <- data.frame(largest_diff_PR$Location)
  diff_sig_PR_events <- rename(diff_sig_PR_events, Location = largest_diff_PR.Location)
  
  # Extract PSI scores for significant events 
  diff_PR_events_for_heatmap <- filtered_PSI_scores[filtered_PSI_scores$Location %in% diff_sig_PR_events$Location,]
  
  # Convert PSI scores to matrix 
  diff_heatmap_matrix_PR <- as.matrix(diff_PR_events_for_heatmap[, 4:3457])
  row.names(diff_heatmap_matrix_PR) <- unlist(diff_PR_events_for_heatmap[,1])
  
  # Build heatmap
  Heatmap(diff_heatmap_matrix_PR, name = "PSI",  
          cluster_rows = TRUE, 
          cluster_columns = TRUE,
          show_row_names = TRUE, 
          row_names_gp = gpar(fontsize = 8), 
          column_names_gp = gpar(fontsize = 6), 
          row_names_side = "left", 
          show_column_names = FALSE, 
          column_title_gp = gpar(fontsize = 8), 
          col = colorRamp2(c(0, .5, 1), hcl_palette = "Blues", rev = TRUE), 
          bottom_annotation = annotationTracks)
}

## Complex heatmap of all significant results (HER2 Pos vs Neg):
{
  # Determine which events/genes show differential expression
  significant_HER2 <- HER2_mw_gene_results %>%
    filter(q_value < 0.05) # 666 events 
  
  # Calculate the differences between means
  means_df_HER2 <- significant_HER2 %>%
    rowwise() %>%
    mutate(mean_diff_max = abs(HER2_positive_mean -
                                 HER2_negative_mean))
  
  # Find events with the largest differences
  largest_diff_HER2 <- means_df_HER2 %>% 
    filter(mean_diff_max > 0.05) # 24 events
  
  # Extract differentially expressed AS events
  diff_sig_HER2_events <- data.frame(largest_diff_HER2$Location)
  diff_sig_HER2_events <- rename(diff_sig_HER2_events, Location = largest_diff_HER2.Location)
  
  # Extract PSI scores for significant events 
  diff_HER2_events_for_heatmap <- filtered_PSI_scores[filtered_PSI_scores$Location %in% diff_sig_HER2_events$Location,]
  
  # Convert PSI scores to matrix 
  diff_heatmap_matrix_HER2 <- as.matrix(diff_HER2_events_for_heatmap[, 4:3457])
  row.names(diff_heatmap_matrix_HER2) <- unlist(diff_HER2_events_for_heatmap[,1])
  
  # Build heatmap
  Heatmap(diff_heatmap_matrix_HER2, name = "PSI",  
          cluster_rows = TRUE, 
          cluster_columns = TRUE,
          show_row_names = TRUE, 
          row_names_gp = gpar(fontsize = 8), 
          column_names_gp = gpar(fontsize = 6), 
          row_names_side = "left", 
          show_column_names = FALSE, 
          column_title_gp = gpar(fontsize = 8), 
          col = colorRamp2(c(0, .5, 1), hcl_palette = "Blues", rev = TRUE), 
          bottom_annotation = annotationTracks)
}

## Start by making a simple heatmap (molecular subtype): 
{
  # Determine which events/genes show differential expression
  significant_subtype <- KW_with_genes %>%
    filter(q_value < 0.05) 
  
  # Extract differentially expressed AS events
  sig_KW_events <- data.frame(significant_subtype$Location)
  sig_KW_events <- rename(sig_KW_events, Location = significant_subtype.Location)
  
  # Extract PSI scores for significant events 
  KW_events_for_heatmap <- filtered_PSI_scores[filtered_PSI_scores$Location %in% sig_KW_events$Location,]
  
  # Convert PSI scores to matrix 
  heatmap_matrix_KW <- as.matrix(KW_events_for_heatmap[, 4:3457])
  row.names(heatmap_matrix_KW) <- unlist(KW_events_for_heatmap[,1])
  flip_heatmap_matrix_KW <- t(heatmap_matrix_KW)
  
  # Create heatmap
  heatmap(flip_heatmap_matrix_KW, main = "All Subtypes", col= colorRampPalette(brewer.pal(8, "Blues"))(25)) # add scale
}

## Complex heatmap of all significant results (by subtype): 
{
  # Build heatmap
  Heatmap(heatmap_matrix_KW, name = "PSI",  
          cluster_rows = FALSE, 
          cluster_columns = TRUE,
          show_row_names = TRUE, 
          row_names_gp = gpar(fontsize = 7), 
          column_names_gp = gpar(fontsize = 6), 
          row_names_side = "left", 
          show_column_names = FALSE, 
          column_title_gp = gpar(fontsize = 8), 
          col = colorRamp2(c(0, .5, 1), hcl_palette = "Blues", rev = TRUE), 
          bottom_annotation = annotationTracks)
  
}

## Complex heatmap with cutoff for largest differences in mean PSI score (by subtype): 
{
  # Calculate the differences between means
  means_df <- significant_subtype %>%
    rowwise() %>%
    mutate(mean_diff_max = max(LumA_mean, LumB_mean, Basal_mean, HER2_mean, Normal_mean) -
             min(LumA_mean, LumB_mean, Basal_mean, HER2_mean, Normal_mean))
  
  #median_diff <- median(means_df$mean_diff_max)
  
  # Find events with greater than 10% difference 
  largest_diff <- means_df %>% 
    filter(mean_diff_max > 0.1)
  
  # Extract differentially expressed AS events
  diff_sig_KW_events <- data.frame(largest_diff$Location)
  diff_sig_KW_events <- rename(diff_sig_KW_events, Location = largest_diff.Location)
  
  # Extract PSI scores for significant events 
  diff_KW_events_for_heatmap <- filtered_PSI_scores[filtered_PSI_scores$Location %in% diff_sig_KW_events$Location,]
  
  # Convert PSI scores to matrix 
  diff_heatmap_matrix_KW <- as.matrix(diff_KW_events_for_heatmap[, 4:3457])
  row.names(diff_heatmap_matrix_KW) <- unlist(diff_KW_events_for_heatmap[,1])
  
  # Build heatmap
  Heatmap(diff_heatmap_matrix_KW, name = "PSI",  
          cluster_rows = TRUE, 
          cluster_columns = TRUE,
          show_row_names = TRUE, 
          row_names_gp = gpar(fontsize = 2), 
          column_names_gp = gpar(fontsize = 6), 
          row_names_side = "left", 
          show_column_names = FALSE, 
          column_title_gp = gpar(fontsize = 8), 
          col = colorRamp2(c(0, .5, 1), hcl_palette = "Blues", rev = TRUE), 
          bottom_annotation = annotationTracks)
  
}

## Group by Event Type (by subtype):    
{
  # Cassette Exons 
  CE_heatmap_events <- diff_KW_events_for_heatmap %>% # 199 events 
    filter(AS_Type == "CE")
  # Convert PSI scores to matrix 
  CE_heatmap_matrix <- as.matrix(CE_heatmap_events[, 4:3457])
  row.names(CE_heatmap_matrix) <- unlist(CE_heatmap_events[,1])
  
  Heatmap(CE_heatmap_matrix, name = "PSI",  
          cluster_rows = T, 
          cluster_columns = TRUE,
          show_row_names = TRUE, 
          row_names_gp = gpar(fontsize = 3), 
          column_names_gp = gpar(fontsize = 6), 
          row_names_side = "left", 
          show_column_names = FALSE, 
          column_title_gp = gpar(fontsize = 8), 
          col = colorRamp2(c(0, .5, 1), hcl_palette = "Blues", rev = TRUE), 
          bottom_annotation = annotationTracks)
  
  # Intron Retention 
  IR_heatmap_events <- diff_KW_events_for_heatmap %>% # 23 events 
    filter(AS_Type == "IR")
  # Convert PSI scores to matrix 
  IR_heatmap_matrix <- as.matrix(IR_heatmap_events[, 4:3457])
  row.names(IR_heatmap_matrix) <- unlist(IR_heatmap_events[,1])
  
  Heatmap(IR_heatmap_matrix, name = "PSI",  
          cluster_rows = T, 
          cluster_columns = TRUE,
          show_row_names = TRUE, 
          row_names_gp = gpar(fontsize = 5), 
          column_names_gp = gpar(fontsize = 6), 
          row_names_side = "left", 
          show_column_names = FALSE, 
          column_title_gp = gpar(fontsize = 8), 
          col = colorRamp2(c(0, .5, 1), hcl_palette = "Blues", rev = TRUE), 
          bottom_annotation = annotationTracks)
  
  # Alternative Acceptors & Donors 
  DA_heatmap_events <- diff_KW_events_for_heatmap %>% # 41 events 
    filter(AS_Type == "AA" | AS_Type == "AD")
  # Convert PSI scores to matrix 
  DA_heatmap_matrix <- as.matrix(DA_heatmap_events[, 4:3457])
  row.names(DA_heatmap_matrix) <- unlist(DA_heatmap_events[,1])
  
  Heatmap(DA_heatmap_matrix, name = "PSI",  
          cluster_rows = TRUE, 
          cluster_columns = TRUE,
          clustering_distance_columns = "spearman",
          show_row_names = TRUE, 
          row_names_gp = gpar(fontsize = 5), 
          column_names_gp = gpar(fontsize = 6), 
          row_names_side = "left", 
          show_column_names = FALSE, 
          column_title_gp = gpar(fontsize = 8), 
          col = colorRamp2(c(0, .5, 1), hcl_palette = "Blues", rev = TRUE), 
          bottom_annotation = annotationTracks)
  
  # Alternative Acceptors 
  AA_heatmap_events <- diff_KW_events_for_heatmap %>% # 126 events 
    filter(AS_Type == "AA")
  # Convert PSI scores to matrix 
  AA_heatmap_matrix <- as.matrix(AA_heatmap_events[, 4:3457])
  row.names(AA_heatmap_matrix) <- unlist(AA_heatmap_events[,1])
  
  Heatmap(AA_heatmap_matrix, name = "PSI",  
          cluster_rows = T, 
          cluster_columns = TRUE,
          show_row_names = TRUE, 
          row_names_gp = gpar(fontsize = 4), 
          column_names_gp = gpar(fontsize = 6), 
          row_names_side = "left", 
          show_column_names = FALSE, 
          column_title_gp = gpar(fontsize = 8), 
          col = colorRamp2(c(0, .5, 1), hcl_palette = "Blues", rev = TRUE), 
          bottom_annotation = annotationTracks)
  
  # Alternative Donors 
  AD_heatmap_events <- diff_KW_events_for_heatmap %>% # 134 events 
    filter(AS_Type == "AD")
  # Convert PSI scores to matrix 
  AD_heatmap_matrix <- as.matrix(AD_heatmap_events[, 4:3457])
  row.names(AD_heatmap_matrix) <- unlist(AD_heatmap_events[,1])
  
  Heatmap(AD_heatmap_matrix, name = "PSI",  
          cluster_rows = T, 
          cluster_columns = TRUE,
          show_row_names = TRUE, 
          row_names_gp = gpar(fontsize = 4), 
          column_names_gp = gpar(fontsize = 6), 
          row_names_side = "left", 
          show_column_names = FALSE, 
          column_title_gp = gpar(fontsize = 8), 
          col = colorRamp2(c(0, .5, 1), hcl_palette = "Blues", rev = TRUE), 
          bottom_annotation = annotationTracks)
}

## Save tables of significant results
{
  largest_diff_ER <- left_join(largest_diff_ER, event_types, by = "Location")
  write.table(largest_diff_ER, file = "~/Desktop/ER_pos_vs_neg.txt", sep = "\t", col.names = NA, quote = FALSE)
  
  largest_diff_HER2 <- left_join(largest_diff_HER2, event_types, by = "Location")
  write.table(largest_diff_HER2, file = "~/Desktop/HER2_pos_vs_neg.txt", sep = "\t", col.names = NA, quote = FALSE)
  
  largest_diff_PR <- left_join(largest_diff_PR, event_types, by = "Location")
  write.table(largest_diff_PR, file = "~/Desktop/PR_pos_vs_neg.txt", sep = "\t", col.names = NA, quote = FALSE)
  
  largest_diff_subtypes <- left_join(largest_diff, event_types, by = "Location")
  write.table(largest_diff_subtypes, file = "~/Desktop/molecular_subtypes.txt", sep = "\t", col.names = NA, quote = FALSE)
}