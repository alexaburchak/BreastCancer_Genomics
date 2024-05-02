library(RColorBrewer)
library(ComplexHeatmap)
library(colorRamp2)
## Start by making a simple heatmap of ER pos vs neg: 
{
  # Determine which events/genes show differential expression
  significant_ER <- ER_mw_gene_results %>%
    filter(q_value < 0.045) 
  
  # Extract differentially expressed genes
  sig_ER_genes <- data.frame(unique(significant_ER$Gene_Name))
  sig_ER_genes <- rename(sig_ER_genes, Gene_Name = unique.significant_ER.Gene_Name.)
  
  # Extract differentially expressed AS events
  sig_ER_events <- data.frame(significant_ER$Location)
  sig_ER_events <- rename(sig_ER_events, Location = significant_ER.Location)
  
  # Extract PSI scores for significant events 
  ER_events_for_heatmap <- filtered_PSI_scores[filtered_PSI_scores$Location %in% sig_ER_events$Location,]
  
  # Convert PSI scores to matrix 
  heatmap_matrix <- as.matrix(ER_events_for_heatmap[, 4:24])
  row.names(heatmap_matrix) <- unlist(ER_events_for_heatmap[,1])
  heatmap_matrix <- t(heatmap_matrix)
  
  # Create heatmap
  heatmap(heatmap_matrix) # add scale
}

## Group by AS type (ER)
{
  CE_heatmap_events <- ER_events_for_heatmap %>% # 957 events 
    filter(AS_Type == "CE")
  # Convert PSI scores to matrix 
  CE_heatmap_matrix <- as.matrix(CE_heatmap_events[, 4:3457])
  row.names(CE_heatmap_matrix) <- unlist(CE_heatmap_events[,1])
  CE_heatmap_matrix <- t(CE_heatmap_matrix)
  # Create heatmap
  heatmap(CE_heatmap_matrix, main = "CE", col= colorRampPalette(brewer.pal(8, "Blues"))(25)) # add scale
  
  AA_heatmap_events <- ER_events_for_heatmap %>% # 184 events 
    filter(AS_Type == "AA") 
  # Convert PSI scores to matrix 
  AA_heatmap_matrix <- as.matrix(AA_heatmap_events[, 4:3457])
  row.names(AA_heatmap_matrix) <- unlist(AA_heatmap_events[,1])
  AA_heatmap_matrix <- t(AA_heatmap_matrix)
  # Create heatmap
  heatmap(AA_heatmap_matrix, main = "AA", col= colorRampPalette(brewer.pal(8, "Blues"))(25)) # add scale
  
  AD_heatmap_events <- ER_events_for_heatmap %>% # 178 events 
    filter(AS_Type == "AD") 
  # Convert PSI scores to matrix 
  AD_heatmap_matrix <- as.matrix(AD_heatmap_events[1:20, 4:24])
  row.names(AD_heatmap_matrix) <- unlist(AD_heatmap_events[1:20,1])
  AD_heatmap_matrix <- t(AD_heatmap_matrix)
  # Create heatmap
  heatmap(AD_heatmap_matrix, col= colorRampPalette(brewer.pal(8, "Blues"))(25), main = "AD") # add scale
  
  IR_heatmap_events <- ER_events_for_heatmap %>% # 197 events 
    filter(AS_Type == "IR") 
  # Convert PSI scores to matrix 
  IR_heatmap_matrix <- as.matrix(IR_heatmap_events[, 4:3457])
  row.names(IR_heatmap_matrix) <- unlist(IR_heatmap_events[,1])
  IR_heatmap_matrix <- t(IR_heatmap_matrix)
  # Create heatmap
  heatmap(IR_heatmap_matrix, col= colorRampPalette(brewer.pal(8, "Blues"))(25), main = "IR") # add scale
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
  heatmap_matrix_KW <- t(heatmap_matrix_KW)
  
  # Create heatmap
  heatmap(heatmap_matrix_KW, main = "All Subtypes", col= colorRampPalette(brewer.pal(8, "Blues"))(25)) # add scale
}

## Complex heatmap of all significant results (by subtype): 
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
  Heatmap(heatmap_matrix_KW, name = "PSI",  
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

## Complex heatmap with cutoff for largest differences in mean PSI score: 
{
  # Calculate the differences between means
  means_df <- significant_subtype %>%
    rowwise() %>%
    mutate(mean_diff_max = max(LumA_mean, LumB_mean, Basal_mean, HER2_mean, Normal_mean) -
             min(LumA_mean, LumB_mean, Basal_mean, HER2_mean, Normal_mean))
  
  median_diff <- median(means_df$mean_diff_max)
  
  # Find events with the largest differences
  largest_diff <- means_df %>% 
    filter(mean_diff_max > median_diff)
  
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

## Group by Event Type   
{
  # Cassette Exons 
  CE_heatmap_events <- diff_KW_events_for_heatmap %>% # 301 events 
    filter(AS_Type == "CE")
  # Convert PSI scores to matrix 
  CE_heatmap_matrix <- as.matrix(CE_heatmap_events[, 4:3457])
  row.names(CE_heatmap_matrix) <- unlist(CE_heatmap_events[,1])
  
  Heatmap(CE_heatmap_matrix, name = "PSI",  
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
  
  # Intron Retention 
  IR_heatmap_events <- diff_KW_events_for_heatmap %>% # 301 events 
    filter(AS_Type == "IR")
  # Convert PSI scores to matrix 
  IR_heatmap_matrix <- as.matrix(IR_heatmap_events[, 4:3457])
  row.names(IR_heatmap_matrix) <- unlist(IR_heatmap_events[,1])
  
  Heatmap(IR_heatmap_matrix, name = "PSI",  
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
  
  # Alternative Acceptors & Donors 
  DA_heatmap_events <- diff_KW_events_for_heatmap %>% # 301 events 
    filter(AS_Type == "AA" | AS_Type == "AD")
  # Convert PSI scores to matrix 
  DA_heatmap_matrix <- as.matrix(DA_heatmap_events[, 4:3457])
  row.names(DA_heatmap_matrix) <- unlist(DA_heatmap_events[,1])
  
  Heatmap(DA_heatmap_matrix, name = "PSI",  
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



