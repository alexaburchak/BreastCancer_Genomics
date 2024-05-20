# 1_Crosscheck
# 2024-04-08
# Alexa Burchak
# 
# Description: Loads multiple datasets, including a list of 743 cancer-associated genes 
# (CGC genes), gene variants, genes for analysis from Mirjam's pipeline, and annotations 
# from GENCODE and NCBI. It identifies CGC genes that are annotated and further filters 
# these to find genes with known variants. A table of FPKM values is loaded to determine 
# gene expression levels in breast tumors. The script identifies CGC genes that are 
# expressed in breast cancer and checks for the presence of known variants. It compares 
# genes identified from Mirjam's pipeline with those expressed in breast cancer, identifying
# genes with and without known variants in each set. The script also examines the expression 
# levels of CGC genes that do not have known variants by analyzing the FPKM values.
# 
# 

library(dplyr)
library(ggplot2)
library(data.table)

# Load list of 743 genes known to be involved in cancer
CGC_genes <- read.table("~/Desktop/Academics/Spring_24/BINP37/Raw_Data/all_CGC_genenames.txt", header = FALSE)

# Some of the gene symbols have been updated in GeneCards 
CGC_genes[90,] <- "CARS1" #CARS
CGC_genes[253,] <- "CEP43" #FGFR1OP 
CGC_genes[295,] <- "H3-3A" #H3F3A
CGC_genes[296,] <- "H3-3B" #H3F3B
CGC_genes[302,] <- "H3C2" #HIST1H3B
CGC_genes[303,] <- "H4C9" #HIST1H4I
CGC_genes[603,] <- "SEPTIN5" #SEPT5
CGC_genes[604,] <- "SEPTIN6" #SEPT6
CGC_genes[605,] <- "SEPTIN9" #SEPT9 

# Load list of known variants 
all_variants <- read.table("~/Desktop/Academics/Spring_24/BINP37/Raw_Data/all_genes_variantexplo.txt", header = FALSE)

# Load list of genes to be analyzed based on Mirjam's pipeline
analyze_genes <- read.table("~/Desktop/Academics/Spring_24/BINP37/Raw_Data/annotated_genes.txt", header = FALSE)

# Load list of genes annotated by GENCODE and/or NCBI
annotations <- read.table("~/Desktop/Academics/Spring_24/BINP37/Raw_Data/all_Gen_NCBI.txt", header = FALSE)

# Convert to vectors 
CGC_genes_v <- CGC_genes$V1
all_variants_v <- all_variants$V1
analyze_genes_v <- analyze_genes$V1
annotations_v <- annotations$V1

# Identify genes in CGC_genes that are annotated 
annotated_CGC_genes <- CGC_genes_v[CGC_genes_v %in% annotations_v] 

# Identify CGC genes also found in all_variants 
genes_in_variants <- annotated_CGC_genes[annotated_CGC_genes %in% all_variants_v] 
CGC_not_variant <- annotated_CGC_genes[!(annotated_CGC_genes %in% all_variants_v)] 

# Identify genes that are expressed in breast tumors 
fpkm_table <- read.table("~/Desktop/Academics/Spring_24/BINP37/Data_For_Crosschecking//fpkm_table.tsv", header = TRUE)
fpkm_melted <- fpkm_table %>%
  melt(id.vars = "Location", variable.name = "Sample", value.name = "FPKM") 
fpkm_melted_v <- unique(fpkm_melted$Location)

# Find genes that have a known variant 
breast_expressed <- genes_in_variants[genes_in_variants %in% fpkm_melted_v]

# Genes identified from Mirjam's pipeline that DO NOT have a known variant:
found_in_analyze <- analyze_genes_v[!(analyze_genes_v %in% breast_expressed)] # 23 genes 

# Genes that have a known variant that WERE NOT identified from Mirjam's pipeline:
found_in_breast <- breast_expressed[!(breast_expressed %in% analyze_genes_v)] # 32 genes 

# Determine expression levels of CGC genes with no known variants: 
subset_FPKM <- fpkm_melted[fpkm_melted$Location %in% CGC_not_variant, ]

# Check how many of the 240 genes removed are expressed in breast cancer but no known variant:
length(unique(subset_FPKM$Location))  

# Distribution of FPKM scores 
subset_FPKM %>%
  ggplot(aes(x=FPKM)) +
  geom_histogram(color="#e9ecef", alpha=0.6) +
  scale_x_continuous(trans = 'log2',labels = function(x) format(x, scientific = FALSE))+
  labs(fill="",
       x = "FPKM",
       y = "Frequency")