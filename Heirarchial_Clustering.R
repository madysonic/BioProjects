# # # # # # # # # # # # # # # # # # # # # # # # #
#         BCR-ABL1 Clustering Analysis          #
# # # # # # # # # # # # # # # # # # # # # # # # #

# Set up ---------------------------------------------------------------

### Load Libraries ###

# For data wrangling
library(dplyr)
library(reshape2)
library(magrittr)

# For analysis
library(edgeR)
library(limma)
library(dendextend)
library(factoextra)

# For visualisation
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)
library(Rtsne)
library(ggrepel)
library(RColorBrewer)

# Set colour palettes

cbPalette <- c("#0072B2", "#E69F00", "#009E73", "#F0E442", "#56B4E9", "#D55E00", "#CC79A7")
brewCols <- brewer.pal(12, "Set3")
brewCols2 <- brewer.pal(8,"Dark2")
myCols <- c(brewCols, cbPalette)
myCols2 <- c(cbPalette, brewCols)

plotCols <- c(brewCols, brewCols2)
plotCols[2] <- cbPalette[1]



# Get DGEObject ------------------------------------------------------


## Our Sample Annotations 
sample_anno <- read.csv("Projects/MLC_sample_annotation.csv")

samples <- sample_anno %>% 
  filter(exclude != "Yes",
         Diagnosis %in% c("B-ALL", "Ph+ALL", "T-ALL", "ETP-ALL"),
         sampleType %in% c("BM", "PB"))

## Counts (subset for our samples)
counts <- readRDS("Rdata/countData.rds")

# Make the row names the gene ID
row.names(counts) <- counts$gene_id
geneIDcol <- which(colnames(counts) == "gene_id")
counts <- counts[ ,-geneIDcol]

# Remove NAs
counts <- counts[complete.cases(counts), ]

samples <- samples[samples$LB_ID %in% colnames(counts),]

# subset count data for only these samples
counts <- counts[ ,match(intersect(samples$LB_ID, colnames(counts)), colnames(counts))]


### Gene Annotation Data ###
genes <- readRDS("Rdata/ensembl_to_hgnc.rds")

# Remove duplicates
genes <- genes[!duplicated(genes$ensembl_gene_id),]
genes <- genes[!duplicated(genes$external_gene_name),]

# Reorder genes to match the counts data
gene_order <- data.frame(ensembl_gene_id = rownames(counts))
genes <- gene_order %>% left_join(genes, by="ensembl_gene_id")

# Drop genes without a gene symbol from both genes and counts datasets

drop.genes <- is.na(genes$external_gene_name)
genes <- genes[!drop.genes,]
counts <- counts[!drop.genes,]

# Set the ensembl gene ID as the rownames
rownames(genes) <- genes$ensembl_gene_id

## DGEList

dge <- DGEList(counts = counts,
               samples = samples,
               genes = genes, 
               remove.zeros = TRUE) %>% 
  calcNormFactors("TMM")


# Remove genes lacking entrez gene ID
drop.genes <- is.na(dge$genes$entrezgene)
dge <- dge[!drop.genes,,keep.lib.sizes=FALSE] 

dge$samples <- dge$samples %>% mutate(col = if_else(is.na(col), "Other", col))

saveRDS(dge, "TandB_ALL_GEObject.rds")



# Function for Clustering  -----------------------------------------------------------

### Function output: 
# complexheatmap of top variable genes
# dendrogram of Ph+ group clusters 
# tSNE of Ph+ clusters
# csv of HCA results


# Number of variable genes to use in clustering
numGenes <- 1000

# Number of clusters (k)
numClusters <- 3

# groupOnly FALSE uses whole cohort to calculate variable genes
# groupOnly TRUE uses Ph+ group only to calculate variable genes

Ph_Cluster <- function(numGenes, numClusters, groupOnly = FALSE){
  # With Ts
  if(groupOnly){
    dge <- readRDS("Projects/Expression/DGE/ALL_DGE/TandB_ALL_GEObject.rds")
    dge <- dge[,dge$samples$col %in% c("Ph"), keep.lib.sizes=FALSE]
    dge <- dge[,dge$samples$sample %in% c("Dx"), keep.lib.sizes=FALSE]
    
    # Get normalised counts without batch effects
    logCPM <- cpm(dge, log = TRUE) 
    batch <- as.factor(c(dge$samples$ref))
    groupNormCounts <- removeBatchEffect(logCPM, batch=batch) 
    
    # Set rownames as gene names to plot
    geneSymbols <- dge$genes$external_gene_name
    rownames(groupNormCounts) <- geneSymbols
    
    # Select most highly variable genes
    var_genes <- apply(groupNormCounts, 1, var)
    top_var_genes <- names(sort(var_genes, decreasing = TRUE))[1:numGenes]  
    
    # Subset log counts for selected genes
    top_var_lcpm <- groupNormCounts[top_var_genes, ]
    
    genes_to_plot <- rownames(top_var_lcpm)
    
  }else{
    dge <- readRDS("Projects/Expression/DGE/ALL_DGE/TandB_ALL_GEObject.rds")
    dge <- dge[,dge$samples$sample %in% c("Dx"), keep.lib.sizes=FALSE]
    
    # Get normalised counts without batch effects
    logCPM <- cpm(dge, log = TRUE) 
    batch <- as.factor(c(dge$samples$ref))
    groupNormCounts <- removeBatchEffect(logCPM, batch=batch) 
    
    # Set rownames as gene names to plot
    geneSymbols <- dge$genes$external_gene_name
    rownames(groupNormCounts) <- geneSymbols
    
    # Select most highly variable genes
    var_genes <- apply(groupNormCounts, 1, var)
    top_var_genes <- names(sort(var_genes, decreasing = TRUE))[1:numGenes]  
    
    # Subset log counts for selected genes
    top_var_lcpm <- groupNormCounts[top_var_genes, ]
    
    genes_to_plot <- rownames(top_var_lcpm)
  }
  
  
  
  # Cluster -----------------------------------------------------------------
  
  dge <- readRDS("Projects/Expression/DGE/ALL_DGE/TandB_ALL_GEObject.rds")
  
  #dge <- readRDS("Projects/Expression/DGE/ALL_DGE/B_ALL_GEObject.rds")
  dge <- dge[,dge$samples$col %in% c("Ph"), keep.lib.sizes=FALSE]
  dge <- dge[,dge$samples$sample %in% c("Dx"), keep.lib.sizes=FALSE]
  
  # List of gene sets to run HCA
  listGenes <- genes_to_plot
  
  # Number of clusters to separate by
  myK <- numClusters
  
  # Clustering and distance methods
  distMethod = "pearson"
  clustMethod = "ward.D2"
  
  ### Perform HCA 
  # Filter DGE by gene list
  groupDGE <- dge[which(dge$genes$external_gene_name %in% listGenes),,keep.lib.sizes = FALSE]
  
  # Get normalised counts without batch effects
  logCPM <- cpm(groupDGE, log = TRUE) 
  batch <- as.factor(c(groupDGE$samples$ref))
  groupNormCounts <- removeBatchEffect(logCPM, batch=batch) 
  
  # Set rownames as gene names to plot
  geneSymbols <- groupDGE$genes$external_gene_name
  rownames(groupNormCounts) <- geneSymbols
  
  # Scale counts
  groupScaled <- scale(t(groupNormCounts)) %>% 
    t()
  
  ## Plot Heatmap
  
  # Set annotation for colour bar of subtypes
  annoHM <- data.frame(groupDGE$samples[,c("col")],
                       #, "primary.subtype")], 
                       row.names = rownames(groupDGE$samples))
  colnames(annoHM) <- c("Group")
  #, "Subtype")
  
  # Named vector for colour list
  colList <- c(myCols2,myCols2)
  keyNames <- unique(groupDGE$samples$col)
  colList <- colList[c(1:length(keyNames))]
  names(colList) <- keyNames
  
  # colList2 <- c(myCols2,myCols2)
  # subNames <- unique(groupDGE$samples$primary.subtype)
  # colList2 <- colList2[c(1:length(subNames))]
  # names(colList2) <- subNames
  
  ## Set heatmap annotation abject
  hAnno = HeatmapAnnotation(df = annoHM, 
                            col = list(Group = colList))
  #, Subtype = colList2))
  
  # Plot heatmap
  jpeg(paste0("Ph_GroupOnly_", groupOnly,"_Top_", numGenes, "_K", myK, "_Heatmap.jpeg"), width = 1000, height = 750, units = "px")
  draw(Heatmap(groupScaled,
               #top_annotation = hAnno,
               name = "Expression",
               clustering_distance_rows = distMethod,
               clustering_distance_columns = distMethod,
               clustering_method_columns = clustMethod,
               clustering_method_rows = clustMethod,
               cluster_rows = T,
               cluster_columns = T,
               show_row_dend = F,
               show_row_names = F,
               show_column_names = T,
               show_column_dend = T,
               column_km = myK,
               #row_names_gp = gpar(fontsize = 3),
               column_names_gp = gpar(fontsize = 5),
               #rect_gp = gpar(col = "white", lwd = 0.1),
               column_title = paste0("HCA of Ph+ Patients using Top ", numGenes, " Variable Genes within All Samples"),
               column_title_gp = gpar(fontsize = 10, fontface = "bold") 
  )) %>% suppressMessages()
  dev.off()
  
  ### Dendrogram
  
  # Cluster scaled counts by correlation
  clustGroup <- stats::hclust(as.dist(1-cor(groupScaled, method=distMethod)), method=clustMethod)
  dendGroup <- as.dendrogram(clustGroup)
  
  
  # Set annotation bar
  # map <- data.frame(col = keyNames, bar=colList)
  # colBar <- groupDGE$samples %>% dplyr::select("patientID", "col")
  # colBar <- left_join(colBar, map, by="col")
  # colBarList <- colBar$bar
  # names(colBarList) <- colBar$patientID
  
  jpeg(paste0("Ph_GroupOnly_", groupOnly,"_Top_",numGenes,"_K", numClusters, "Dendrogram.jpeg"), width = 1000, height = 500, units = "px")
  dendGroup %>%
    # Node visuals
    set("nodes_pch", 19) %>% 
    set("nodes_cex", 0.7) %>%
    # Branch visuals
    set("branches_k_color", k = myK, groupLabels = c(1:myK)) %>%
    set("branches_lwd", 2) %>%
    # Label visuals
    set("labels_cex", c(0.5)) %>%
    plot(main= paste0("Dendrogram of Ph+ Patients by Top", numGenes, "Variable Genes"))
  #, cex.main = 2, font.main=2, cex = 1)
  #colored_bars(colors = colBarList, dend = dendGroup, sort_by_labels_order = FALSE, rowLabels = "Group", cex.rowLabels = 1)
  dev.off()
  
  ### Label each sample with its cluster ###
  
  cut <- dendextend::cutree(clustGroup, k=myK, order_clusters_as_data = FALSE)
  
  # Order sample list by cutree list
  groupDGE$samples <- groupDGE$samples[match(names(cut), groupDGE$samples$LB_ID),]
  
  # Bind group list to sample list
  grouped_samples <- cbind(cut, groupDGE$samples)
  grouped_samples$HCA <- grouped_samples[,1]
  grouped_samples <- grouped_samples[,-1] 
  
  ### Generate TSNE
  
  ## Subset logcounts matrix with the selected genes
  topVar <- groupNormCounts
  
  # Transpose for plotting
  tVar_logCPM <- data.frame(t(topVar))
  
  # generate df containing samples for labelling plot
  tsne_anno <- as.data.frame(factor(rownames(tVar_logCPM)))
  colnames(tsne_anno) <- "Sample"
  
  # Add anno
  tsne_anno$ploidy <- groupDGE$samples$Ploidy[match(tsne_anno$Sample, groupDGE$samples$LB_ID)]
  tsne_anno$col <- groupDGE$samples$col[match(tsne_anno$Sample, groupDGE$samples$LB_ID)]
  tsne_anno$HCA <- grouped_samples$HCA[match(tsne_anno$Sample, grouped_samples$LB_ID)]
  tsne_anno$patientID <- groupDGE$samples$patientID[match(tsne_anno$Sample, groupDGE$samples$LB_ID)]
  tsne_anno$secondary <- groupDGE$samples$secondary.subtype[match(tsne_anno$Sample, groupDGE$samples$LB_ID)]
  
  
  ## Generate tSNE plot on most variable genes
  set.seed(42) # for reproducibility
  tsne_out <- Rtsne(tVar_logCPM, pca = FALSE, perplexity = 30, theta = 0.0, check_duplicates = FALSE)
  
  # Visualise plot
  tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2], 
                          HCA = as.character(tsne_anno$HCA),
                          LB_ID = tsne_anno$Sample,
                          col = tsne_anno$col,
                          secondary = tsne_anno$secondary,
                          ploidy = tsne_anno$ploidy)
  
  #tsne_plot <- tsne_plot %>% mutate(plot_label = "")
  #tsne_plot <- tsne_plot %>% mutate(plot_label = ifelse(timepoint == "ER_REL1", patientID, ""))
  
  
  # Assign colour palette 
  palette_tsne <- myCols2[1:myK]
  names(palette_tsne) <- c(1:myK)
  
  #palette_tsne <- myCols2[1:length(keyNames)]
  #names(palette_tsne) <- keyNames
  
  jpeg(paste0("Ph_GroupOnly_", groupOnly,"_Top_",numGenes,"_K", numClusters, "_tSNE.jpeg"), width = 700, height = 700, units = "px")
  
  tsnePlot <- ggplot(tsne_plot, aes(x,y,color=secondary)) + 
    geom_point(data = tsne_plot, size = 3) + 
    theme_bw() + 
    labs(title = paste0("t-SNE using top variable genes"),       
         colour = "Group") + 
    #scale_colour_manual(values = palette_tsne, na.value = "grey50") + 
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_blank())
  
  print(tsnePlot)
  #ggsave(paste0("Ph_Only_Top", numGenes,"_K", myK, "_tsne.jpeg"), height = 6, width = 7)
  dev.off()
  
  
  write.csv(grouped_samples, paste0("HCA_Ph_", numGenes, "groupOnly_", groupOnly, "_K", myK, ".csv"))
}

# Percentage Plots --------------------------------------------------------------



hcaPerc <- grouped_samples %>% 
  count(col, HCA) %>% 
  group_by(col) %>% 
  mutate(label_y = n/sum(n)) 


# Plot percentages 
hcaPlot <- hcaPerc %>% 
  ggplot(aes(x=col, y=label_y, fill=as.character(HCA))) +
  geom_col(position="stack") +
  geom_text(aes(label=paste(paste0(round(label_y*100, digits=0), "%"), paste0("(", n, ")"))), 
            position=position_stack(vjust = 0.5),
            colour="white",
            size=3,
            fontface="bold") +
  scale_fill_manual(values= plotCols) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none") +
  labs(x = "Subtype",
       y = "Percentage of Samples") +
  ggtitle("Classification of patients by unsupervised HCA")

hcaPlot 


write.csv(grouped_samples, "HCA_Ph_1000PhOnly_K2.csv")

