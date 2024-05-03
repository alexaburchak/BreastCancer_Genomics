# title: "1_Crosscheck"
# date: "2024-04-08"
# description: ADD 

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
annotated_CGC_genes <- CGC_genes_v[CGC_genes_v %in% annotations_v] # 6 not annotated

# Identify CGC genes also found in all_variants 
genes_in_variants <- annotated_CGC_genes[annotated_CGC_genes %in% all_variants_v] 
CGC_not_variant <- annotated_CGC_genes[!(annotated_CGC_genes %in% all_variants_v)] # 240 that were removed from CGC 

# Identify genes that are expressed in breast tumors 
# Remove FPKM scores lower than 10 before comparing gene names:
fpkm_table <- read.table("~/Desktop/Academics/Spring_24/BINP37/Raw_Data/Raw_Tables/fpkm_table.tsv", header = TRUE)
fpkm_melted <- fpkm_table %>%
  melt(id.vars = "Location", variable.name = "Sample", value.name = "FPKM") %>%
  filter(FPKM > 10) # based on Mirjam's analysis 
fpkm_melted_v <- unique(fpkm_melted$Location)

# Find genes that have a known variant and FPKM>10 
breast_expressed <- genes_in_variants[genes_in_variants %in% fpkm_melted_v] # 466 genes

# Genes identified from Mirjam's pipeline that DO NOT have a known variant AND FPKM>10:
found_in_analyze <- analyze_genes_v[!(analyze_genes_v %in% breast_expressed)] # 23 genes 

# Genes that have a known variant and FPKM>10 that WERE NOT identified from Mirjam's pipeline:
found_in_breast <- breast_expressed[!(breast_expressed %in% analyze_genes_v)] # 32 genes 

# Determine expression levels of CGC genes with no known variants: 
subset_FPKM <- fpkm_melted[fpkm_melted$Location %in% CGC_not_variant, ]

# Check how many of the 240 genes removed have a FPKM > 10 but no known variant:
length(unique(subset_FPKM$Location)) # 206 genes 

# Distribution of FPKM scores 
subset_FPKM %>%
  ggplot(aes(x=FPKM)) +
  geom_histogram(color="#e9ecef", alpha=0.6) +
  scale_x_continuous(trans = 'log2',labels = function(x) format(x, scientific = FALSE))+  
  labs(fill="", 
       x = "FPKM",
       y = "Frequency")



# low expression = low information (not enough reads to determine variance, we only keep variance that is supported by at least 5 reads)
