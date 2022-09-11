Classifying B-Others
================
Madyson Colton
27/04/2022

# Exploration of Unclassified B-ALL Patients

This analysis explores the gene expression profiles of Acute
Lymphoblastic Leukaemia (ALL) patients with no known structural variants
or driver mutations. These patients are broadly referred to as the
“B-Other” subgroup.

Previous exploration of this subgroup revealed that several of our
samples are unclassified due to missing cytogenetic data, and several of
the “B-Other” samples should actually be classified as Hyperdiploid ALL
or Hypodiploid ALL.

## Set Up

Data from previous DGE analysis was loaded, including sample, count and
gene data, as well as a DGE object of all samples in the B-ALL cohort
and an output file containing the top 50 significant differentially
expressed genes in each subtype compared to all other subtypes.

``` r
### Data Wrangling ###
library(dplyr)
library(stringr)
library(limma)

### Visualisation and Analysis ###
library(ComplexHeatmap)
library(ggplot2)
library(dendextend)


filepath <- ("C:/Users/madyson/OneDrive - South Australian Health and Medical Research Institute Ltd/LBank/Projects/Expression/DGE/ALL_DGE")

samples <- readRDS(paste0(filepath, "/samples.rds"))
counts <- readRDS(paste0(filepath, "/counts.rds"))
genes <- readRDS(paste0(filepath, "/genes.rds"))

dge <- readRDS(paste0(filepath, "/DGE_ALL.rds"))

# Output file is merged csv of the top 50 significant DEGs in each subtype

output <- read.csv(paste0(filepath, "/topSigs/output.csv"))
```

## Previous expression profiling

Differential gene expression analysis has previously been performed on
each major B-ALL subtype against the rest of the patient cohort to
identify the most significantly variable genes in each subtype. The top
50 significant genes from each group were used to cluster patients.

``` r
### Modify output file ###
# Replace ensembl ID with gene symbol
for(j in 1:ncol(output)){
  for(i in 1:nrow(output)){
    # To prevent errors from subtypes with less than 50 top genes
    if(!is.na(output[i,j])){
      output[i,j] <- genes[c(which(genes$ensembl_gene_id %in% output[i,j])),"external_gene_name"]
    }else{
      output[i,j] <- NA
    }
  }
}


# List of all genes of interest
listTop <- list()

for(i in 1:ncol(output)){
  geneIDs <- output[,i]
  listTop <- append(listTop, geneIDs) 
}

listTop <- listTop[!is.na(listTop)]


# Drop any genes that are duplicated
dupList <- listTop[which(duplicated(listTop))]
listTop <- listTop[-which(listTop %in% dupList)]


### ALL DGE ###

output <- output[,c(1:11)]

### Heatmap of Top Genes in All Subtypes ###
names_to_plot <- c(colnames(output), "Other")

allDGE <- dge[which(dge$genes$external_gene_name %in% listTop),,keep.lib.sizes = FALSE]
```

    ## Loading required package: edgeR

``` r
allDGE <- allDGE[,which(allDGE$samples$col %in% names_to_plot),keep.lib.sizes = FALSE]

# Get normalised counts without batch effects
samplesA <- allDGE$samples
logCPM <- cpm(allDGE, log=TRUE)
batch <- as.factor(c(samplesA$ref))
countsNorm <- removeBatchEffect(logCPM, batch=batch) 

# Set rownames as gene names to plot
geneSymbols <- allDGE$genes$external_gene_name
rownames(countsNorm) <- geneSymbols

### Heatmap Prep ##

# Scale counts
scaledCounts <- scale(t(countsNorm)) %>% 
  t()

# Set annotation for colour bar of subtypes
annoHM <- data.frame(allDGE$samples[,c("col")], 
                     row.names = rownames(allDGE$samples))
colnames(annoHM) <- c("Subtype")

# Named vector for colour list
colList <- c("#8DD3C7","#FFFFB3", "#0072B2","#E69F00","#80B1D3", "#B3DE69", "#D9D9D9","#FB8072","#CCEBC5", "#BEBADA", "#FFED6F", "grey")


names(colList) <- names_to_plot


hAnno = HeatmapAnnotation(df = annoHM, 
                          col = list(Subtype = colList))

# Plot heatmap
allHM <- Heatmap(scaledCounts,
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "pearson",
        cluster_rows = T,
        cluster_columns = T,
        show_row_dend = F,
        show_row_names = T,
        show_column_dend = T,
        row_names_gp = gpar(fontsize = 3),
        column_names_gp = gpar(fontsize = 3),
        top_annotation = hAnno,
        name = "Expression",
        #column_km = 2,
        #column_title = c("Hyperdiploid", "Not Hyperdiploid"),
        #column_title_gp = gpar(fontsize = 6),
        rect_gp = gpar(col = "white", lwd = 0.1)
        )
```

### Plotting using significant Hyperdiploid and Hypodiploid genes

When focusing on the aneuploid and “B-Other” subtype, we can clearly see
that some of the patients cluster in the middle of the Hyperdiploid
group.

``` r
### Complex Heatmap ###

### Data Prep ###

# Create replicate DGE for testing
hyperDGE <- dge

# Remove  B-Other samples classified as ETV6-RUNX1-Like or NUTM1
samplesKeep <- hyperDGE$samples[which(hyperDGE$samples$primary.subtype != "ETV6-RUNX1-like"),]
hyperDGE <- hyperDGE[,which(hyperDGE$samples$sampleID %in% samplesKeep$sampleID),keep.lib.sizes = FALSE]

samplesKeep2 <- hyperDGE$samples[which(hyperDGE$samples$primary.subtype != "NUTM1"),]
hyperDGE <- hyperDGE[,which(hyperDGE$samples$sampleID %in% samplesKeep2$sampleID),keep.lib.sizes = FALSE]

# Select subtypes and genes to plot, then filter DGE object

names_to_plot <- c("Hyperdiploid", "Hypodiploid", "Other")
genes_to_plot <- c(output$Hyperdiploid, output$Hypodiploid)

hyperDGE <- hyperDGE[which(hyperDGE$genes$external_gene_name %in% genes_to_plot),,keep.lib.sizes = FALSE]
hyperDGE <- hyperDGE[,which(hyperDGE$samples$col %in% names_to_plot),keep.lib.sizes = FALSE]

# Get normalised counts without batch effects
samplesA <- hyperDGE$samples
hyperCPM <- cpm(hyperDGE, log = TRUE) 
batch <- as.factor(c(samplesA$ref))
hypNormCounts <- removeBatchEffect(hyperCPM, batch=batch) 

# Set rownames as gene names to plot
geneSymbols <- hyperDGE$genes$external_gene_name
rownames(hypNormCounts) <- geneSymbols


### Heatmap Prep ##

# Scale counts
hypScaledCounts <- scale(t(hypNormCounts)) %>% 
  t()

# Set annotation for colour bar of subtypes
annoHM <- data.frame(hyperDGE$samples[,c("col", "primary.subtype")], 
                     row.names = rownames(hyperDGE$samples))
colnames(annoHM) <- c("Group", "Subtype")


# Save as Heatmap Annotation object and map colours
hAnno = HeatmapAnnotation(df = annoHM, col = list(Group = c("Hyperdiploid" =  "#0072B2", 
                                                         "Hypodiploid" = "#E69F00",
                                                         "Other" = "grey"), 
                                               Subtype = c("High Hyperdiploid" = "#A6CEE3",
                                                           "High Hypodiploid" = "#1F78B4",
                                                           "Hyperdiploid" = "#B2DF8A",
                                                           "Hypodiploid"  = "#33A02C",
                                                           "Low Hyperdiploid"  = "#FB9A99",
                                                           "Low Hypodiploid" = "#E31A1C",
                                                           "Masked Hypodiploid" = "#FDBF6F",
                                                           "Near Haploid" = "#FF7F00",
                                                           "Other"= "#CAB2D6")))


# Plot heatmap
Heatmap(hypScaledCounts,
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "pearson",
        cluster_rows = T,
        cluster_columns = T,
        show_row_dend = F,
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 5),
        top_annotation = hAnno,
        name = "Expression",
        column_km = 2,
        column_title = c("Hyperdiploid", "Not Hyperdiploid"),
        column_title_gp = gpar(fontsize = 6),
        rect_gp = gpar(col = "white", lwd = 1))
```

![](Identify_Hyperdip_files/figure-gfm/hyperHeatmap-1.png)<!-- -->

### Clustering based on top genes

Hierarchical clustering was used to classify the “B-Others” and the
aneuploid subtypes, using Pearson correlation as the distance metric.

``` r
# Cluster scaled counts by correlation
clustHyp <- hclust(as.dist(1-cor(hypScaledCounts, method="pearson")), method="complete")
dendHyp <- as.dendrogram(clustHyp)

dendPlot <- dendHyp %>%
  # Node visuals
  set("nodes_pch", 19) %>% 
  set("nodes_cex", 0.7) %>%
  # Branch visuals
  set("branches_k_color", k = 7, groupLabels = c(1:7)) %>%
  set("branches_lwd", 2) %>%
  # Label visuals
  set("labels_cex", c(0.5)) %>%
  plot(main= "Clustering of B-Others and Aneuploid Patients") 
```

![](Identify_Hyperdip_files/figure-gfm/HCA-1.png)<!-- -->

``` r
### Label each sample with its cluster ###

cut <- dendextend::cutree(clustHyp, k=7, order_clusters_as_data = FALSE)

# Order sample list by cutree list
hyperDGE$samples <- hyperDGE$samples[match(names(cut), hyperDGE$samples$sampleID),]

# Bind group list to sample list
grouped_samples <- cbind(cut, hyperDGE$samples)
grouped_samples$HCA <- grouped_samples[,1]
grouped_samples <- grouped_samples[,-1] %>% 
  select("patientID", "LB_ID", "key_alt", "col", "primary.subtype", "sampleID", "HCA")

# Filter for potential hyperdiploids
hypOthers <- grouped_samples %>% 
  filter(HCA == 2,
         col == "Other")
```

    ## NULL

## RNAseqCNV

To confirm that these samples had chromosomal alterations, a
computational CNV caller was used. This caller was previously
benchmarked on our dataset against samples with cytogenetic reports and
was found to be 84% accurate at identifying chromosomal alterations. In
patients where aneuploidy was the primary alteration (as oppose to gene
fusion events), the tool was 94% accurate.

This tool was ran on all unclassfied “B-Other” subtypes (n=80).

``` r
cbPalette <- c("#0072B2", "#E69F00", "#009E73", "#F0E442", "#56B4E9", "#D55E00", "#CC79A7")

### Load Data ###
setwd("C:/Users/madyson/OneDrive - South Australian Health and Medical Research Institute Ltd/LBank/Projects/CNVs")
cnvTable <- read.table("combined_table.csv", sep = ",", header = T)
sampleAnno <- read.csv("C:/Users/madyson/OneDrive - South Australian Health and Medical Research Institute Ltd/LBank/meta_data/sample_annotation.csv")
karyoData <- read.csv("C:/Users/madyson/OneDrive - South Australian Health and Medical Research Institute Ltd/LBank/meta_data/ALL_karyotype_data.csv")

# Merge into singular DF
allData <- left_join(cnvTable, sampleAnno, by = c("sample" = "variant_file_ID")) %>%
  left_join(karyoData, by = "LB_ID")

# Create column to extract karyotype info
allData <- allData %>% 
  mutate(karyo_n = str_sub(Karyotype, start=2, end=3))

# Handle NAs in karyo_n
allData$karyo_n <- as.numeric(allData$karyo_n) 
allData$karyo_n[is.na(allData$karyo_n)] <- 0


### Filter for B-Others ###

BData <- allData %>% filter(exclude != "Yes",
                            is.na(col),
                            Diagnosis.x %in% c("B-ALL", "Ph+ALL"),
                            sample.y == "Dx")


### Plot RNAseqCNV results ###

# Ensure columns of interest are numeric
BData$karyo_n <- as.numeric(BData$karyo_n)
BData$chrom_n <- as.numeric(BData$chrom_n)

# Reorder for plotting
BData <- BData[order(BData$chrom_n, decreasing = TRUE),]

BData$patientID.x <- factor(BData$patientID.x,
                                   levels = BData$patientID.x[order(BData$chrom_n, decreasing = TRUE)])

# Labels for plotting
y_labs <- BData$primary.subtype

# Plot results
BNumbers <- BData %>% 
  ggplot() +
  geom_col(aes(x=patientID.x, y=chrom_n, fill = cbPalette[2]), width= 0.7) +
  geom_col(aes(x=patientID.x, y= karyo_n, fill = "white"), colour = "black", width = 0.7, alpha = 0.1) +
  labs(x = "Sample",
       y = "Number of Chromosomes") +
  ggtitle("Chromosome Number Reported in 'B-Other' Samples") +
  theme_classic() +
  #theme(axis.text.y = element_blank(), axis.line.y = element_blank()) +
  coord_flip(ylim = c(20,65)) +
  geom_hline(aes(yintercept = 46, linetype= paste("Normal (n = 46)"))) +
  scale_fill_manual(name = "Data Source", values = c(cbPalette[2], "white"), labels=c("RNAseqCNV", "Cytogenetics")) +
  scale_linetype_manual("Line", values = 1) +
  scale_x_discrete(labels = y_labs)

BNumbers
```

![](Identify_Hyperdip_files/figure-gfm/RNAseqCNV-1.png)<!-- -->

### Overlap with Heirachal Clustering

HCA pulled out 12 “B-Other” samples clustered within the “Hyperdiploid”
group. The number of alterations called by RNAseqCNV was assessed.

``` r
overlap <- left_join(hypOthers, BData, by = "LB_ID")

overlap <- overlap[order(overlap$chrom_n, decreasing = TRUE),]

overlap$patientID.x <- factor(overlap$patientID.x,
                                   levels = overlap$patientID.x[order(overlap$chrom_n, 
                                                                      decreasing = TRUE)])

# Labels for plotting
y_labs <- overlap$primary.subtype.y
y_labs2 <- overlap$sampleID

# Plot results
overlapNumbers <- overlap %>% 
  ggplot() +
  geom_col(aes(x=patientID.x, y=chrom_n, fill = cbPalette[2]), width= 0.7) +
  #geom_col(aes(x=patientID.x, y= karyo_n, fill = "white"), colour = "black", width = 0.7, alpha = 0.1) +
  labs(x = "Sample",
       y = "Number of Chromosomes") +
  ggtitle("Chromosome Number Reported in 'B-Other' Samples") +
  theme_classic() +
  #theme(axis.text.y = element_blank(), axis.line.y = element_blank()) +
  coord_flip(ylim = c(20,65)) +
  geom_hline(aes(yintercept = 46, linetype= paste("Normal (n = 46)"))) +
  scale_fill_manual(name = "Data Source", values = c(cbPalette[2], "white"), labels=c("RNAseqCNV")) +
 scale_linetype_manual("Line", values = 1) +
 scale_x_discrete(labels = y_labs2)

overlapNumbers
```

![](Identify_Hyperdip_files/figure-gfm/overlap%20chrm%20number-1.png)<!-- -->

``` r
write.csv(overlap, "hyperB_to_reclassify.csv")
```
