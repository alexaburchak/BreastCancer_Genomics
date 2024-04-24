## Start by making a simple heatmap: 
{
  # Determine which events/genes show differential expression
  significant_ER <- ER_mw_gene_results %>%
    filter(q_value < 0.05) 
  
  # Extract differentially expressed genes
  sig_ER_genes <- data.frame(unique(significant_ER$Gene_Name))
  sig_ER_genes <- rename(sig_ER_genes, Gene_Name = unique.significant_ER.Gene_Name.)
  
  # Extract differentially expressed AS events
  sig_ER_events <- data.frame(unique(significant_ER$Location))
  sig_ER_events <- rename(sig_ER_events, Location = unique.significant_ER.Location.)
  
  # Extract PSI scores for significant events 
  ER_events_for_heatmap <- filtered_PSI_scores[filtered_PSI_scores$Location %in% sig_ER_events$Location, ]
  
  # Convert PSI scores to matrix 
  heatmap_matrix <- as.matrix(ER_events_for_heatmap[, 4:24])
  row.names(heatmap_matrix) <- unlist(ER_events_for_heatmap[,1])
  heatmap_matrix <- t(heatmap_matrix)
  
  # Create heatmap
  heatmap(heatmap_matrix) # add scale
}

## Group by AS type 
{
  
}



## Complex heatmap: 
{
  # Add annotations 
  ER_PR_HER2_status$NCN_PAM50 <- sub("Her2", "HER2", ER_PR_HER2_status$NCN_PAM50)
  ER_PR_HER2_status$NCN_PAM50 <- sub("Normal", "Normal-like", ER_PR_HER2_status$NCN_PAM50)
  ER_PR_HER2_status$NCN_PAM50<- sub("unclassified", "Unclassified", ER_PR_HER2_status$NCN_PAM50)
  
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
  
  # Build heatmap (not complete)
  Heatmap(psi, name = "PSI", 
          cluster_rows = F, 
          cluster_columns = T, 
          show_row_names = T, 
          row_names_gp = gpar(fontsize = 8), 
          column_names_gp = gpar(fontsize = 7), 
          row_names_side="left", 
          show_column_names = F, 
          column_title_gp = gpar(fontsize = 8), 
          col = colorRamp2(c(0, .5, 1), hcl_palette = "Blues", rev = T), 
          bottom_annotation = annotationTracks)
  
}


