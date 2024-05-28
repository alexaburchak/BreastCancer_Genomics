# Investigating Differential Splicing Patterns in Breast Cancer Tumors
## Overview
This repository contains several R scripts necessary for the analysis and visualization of 
differential alternative splicing in 3454 breast tumor samples from the SCAN-B cohort. 


This set of R scripts will perform the following steps:
1. Identify genes for analysis
2. Calculate and visualize summary statistics
3. Perform differential splicing analysis 
4. Generate heatmaps for significant events
5. Create boxplots for notable events in *tenascin c*   

**Requirements**

The clinical data, including information on receptor status and PAM50 subtype, comes from the Sweden
Cancerome Analysis Network- Breast database. Information on alternative splicing events and PSI 
scores for each event were provided by [Mirjam Müller](https://github.com/TheOrangeBraincell/variants_in_AS_Pipeline) 
along with instructions for how to download gene annotations from GENCODE and RefSeq.  

All code was run in R Studio (v4.3.2 or higher) using the following packages:

1_Crosscheck
```
install.packages(c("dplyr", "ggplot2", "data.table"))
```

2_SummaryStats
```
install.packages(c("dplyr", "ggplot2", "data.table", "reshape2", "tidyr", "readxl", "hrbrthemes"))
```

3_DiffAnalysis
```
install.packages(c("geomtextpath", "dplyr", "ggplot2", "data.table", "hrbrthemes"))
```

4_Heatmaps
```
install.packages(c("dplyr", "RColorBrewer", "ComplexHeatmap", "colorRamp2"))
```

5_Boxplots
```
install.packages(c("dplyr", "ggplot2", "tidyr", "stringr")) 
```

## Contents 
### 1_Crosscheck
The first of the data processing steps is to identify genes from the Cancer Genome Census (CGC) that 
are annotated by RefSeq/GENCODE, expressed in breast cancer and have a known variant. The goal is
to ensure that this list matches the list of genes with PSI scores calculated by Mirjam's pipeline
before proceeding to differential expression analysis and to determine the cause for exclusion for 
CGC genes with no PSI scores.  

### 2_SummaryStats
Once genes with PSI scores have been confirmed, we can visualize the distribution of mean PSI scores 
from each AS event. All non-numeric PSI scores are converted to NA and removed before mean PSI scores 
are calculated for each AS event. We also remove PSI scores of 0 or 1 before visualization since 
these are not likely to be abnormal AS events. 

### 3_DiffAnalysis
The distribution of PSI scores is not normal, so in order to analyze differential splicing among 
samples based on ER, PR, and HER2 status (+ vs -), we use for-loops to calculate Mann-Whitney U tests 
for each splicing event. Additionally, we use a for-loop to calculate Kruskal-Wallis statistics for 
PAM50 molecular subtypes (LumA, LumB, HER2, Basal, Normal). Benjamini-Hochberg correction is applied 
to control the false discovery rate in all tests. We conclude by generating summary statistics and 
extracting AS events with significant differential splicing (q < 0.05) for further analysis.   

### 4_Heatmaps
The final step in our analysis of differential splicing in breast cancer is to visualize AS events 
with significant differential splicing in a heatmap. For this we will use the complexHeatmap package 
(v2.18.0) and a custom annotation track for receptor status and PAM50 subtype. Since the differences
in mean PSI scores are minimal even in significant events, our final heatmaps for analysis are 
generated using a cutoff of 10% difference in mean score across groups. 

### 5_Boxplots 
Our heatmap of cassette exon events with significant differential splicing across PAM50 subtypes 
revealed one particularly interesting cluster of events that occur in *tenascin c*, a gene known for its role in cell proliferation and migration. We further study
these events by first creating a boxplot of log2-transformed FPKM values for *TNC*  across 
different PAM50 subtypes. We also generate a boxplot of PSI scores for the 6 specific alternative 
splicing events within *TNC*, highlighting variation in alternative splicing across subtypes. 

## Usage
Each script can be run sequentially to reproduce the analysis. Ensure that all dependencies, 
including the complexHeatmap package, are installed and that the necessary input data is available. 
Detailed comments within each script provide guidance on execution and interpretation of results.






