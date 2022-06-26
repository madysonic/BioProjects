# # # # # # # # # # # # # # # # # # # # # # # # # # # #
#      Weighted Gene Correlation Network Analysis     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

# B-ALL

### Load Libraries ###

# For data wrangling
library(tidyverse)
library(reshape2)
library(magrittr)

# For analysis
library(WGCNA)
library(impute)
library(flashClust)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationHub)
library(DOSE)
library(STRINGdb)

# For visualisation
library(ggplot2)
library(cowplot)
library(plotly)
library(pheatmap)
library(RColorBrewer)
cbPalette <- c("#0072B2", "#E69F00", "#009E73", "#F0E442", "#56B4E9", "#D55E00", "#CC79A7")

# Load Data ---------------------------------------------------------------

dge <- readRDS("Expression/DGE/ALL_DGE/TandB_ALL_GEObject.rds")

# Remove T-ALL samples and relapse samples
dge <- dge[,which(dge$samples$Diagnosis %in% c("B-ALL", "Ph+ALL")), keep.lib.sizes = FALSE]
dge <- dge[,dge$samples$sample %in% c("Dx"), keep.lib.sizes=FALSE]

## Add HCA annotation
# pSamps <- read.csv("Clustering/HCA_Ph_1000groupOnly_FALSE_K3.csv")
# dge$samples$HCA <- pSamps$HCA[match(dge$samples$LB_ID, pSamps$LB_ID)] %>% as.character()
# dge$samples <- dge$samples %>% mutate(HCA = if_else(is.na(HCA), "NoHCA", HCA))


# Remove batch effects
logCPM <- cpm(dge, log=TRUE)
batch <- as.factor(c(dge$samples$ref))
logCPM_no_batch <- removeBatchEffect(logCPM, batch=batch) 

# Get var genes
var_genes <- apply(logCPM_no_batch, 1, var)
top_var_genes <- names(sort(var_genes, decreasing = TRUE))[1:10000]  

# Subset log counts for selected genes
varCounts <- logCPM_no_batch[top_var_genes, ]

genes_to_plot <- rownames(varCounts)

# Transpose for analysis (use logCPM_no_batch or varCounts)
norm_counts <-varCounts %>%
  t %>% 
  as.data.frame

# Define contrasts of interest
design <- model.matrix(~0 + primary.subtype, data = dge$samples)
design_df = as.data.frame(design)[]



# Construct Network -------------------------------------------------------

# netwk <- readRDS("Expression/WGCNA/networkWGCNA.rds")

powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Network topology analysis to determine power
soft_thresh <- pickSoftThreshold(
  norm_counts,
  powerVector = powers,
  verbose = 0
)

# Plot thresholds for scale-free topology and mean connectivity
jpeg("1a_SFT_plot.jpeg")

sftPlot <- plot(soft_thresh$fitIndices[, 1],
                -sign(soft_thresh$fitIndices[, 3]) * soft_thresh$fitIndices[, 2],
                xlab = "Soft Threshold (power)",
                ylab = "Scale Free Topology Model Fit, signed R^2",
                main = paste("Scale independence"),
                text(soft_thresh$fitIndices[, 1],
                     -sign(soft_thresh$fitIndices[, 3]) * soft_thresh$fitIndices[, 2],
                     labels = powers, cex = 0.9, col = "red")) + abline(h = 0.9, col = "red")

dev.off()

jpeg("1b_mean_connectivity_plot.jpeg")
mcPlot <- plot(soft_thresh$fitIndices[, 1],
               soft_thresh$fitIndices[, 5],
               xlab = "Soft Threshold (power)",
               ylab = "Mean Connectivity",
               type = "n",
               main = paste("Mean connectivity")) +
  text(soft_thresh$fitIndices[, 1],
       soft_thresh$fitIndices[, 5],
       labels = powers,
       cex = 0.9, col = "red")

dev.off()

# Can also determine the estimated power from topology analysis
soft_thresh$powerEstimate

# Input chosen power
picked_power = 6

### Construct Network ###

# Need to change the correlation function temporarily
temp_cor <- cor       
cor <- WGCNA::cor         

# Construct network and identify modules
netwk <- blockwiseModules(norm_counts,                
                          power = picked_power,                
                          networkType = "signed",
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          minModuleSize = 20,
                          maxBlockSize = 5000,
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          saveTOMs = T,
                          saveTOMFileBase = "BALL_WGCNA",
                          numericLabels = T,
                          verbose = 3)

saveRDS(netwk, "networkB.rds")


# Analysis ----------------------------------------------------------------

# See how many modules were identified and their sizes
table(netwk$colors)

# Return correlation function back to normal
cor <- temp_cor 

# Convert labels to colors for plotting

netwkColours <- labels2colors(netwk$colors)

# Plot the dendrogram and the module colors underneath
jpeg("2_WGCNA_cluster_plot.png", width=1000, height=500, units = "px") 

clustPlot <- plotDendroAndColors(
  netwk$dendrograms[[1]],
  netwkColours[netwk$blockGenes[[1]]],
  "Module Colours",
  main = "Gene dendrogram and module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

dev.off()


# Explore Modules ---------------------------------------------------------

# Get Module Eigengenes per cluster
MEs <- moduleEigengenes(norm_counts, netwkColours)$eigengenes

# Reorder modules so similar modules are next to each other
MEs <- orderMEs(MEs)

### Determine Modules Related to Diagnosis ###

# Correlate modules with traits
moduleTraitCor <- cor(MEs, design_df, use = "p");
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(dge$counts))

# Create text matrix for plotting
textMatrix <- paste(signif(moduleTraitCor, 2), " (",
                    signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) <- dim(moduleTraitCor)

# Plot a heatmap of correlation
sizeGrWindow(10,8)
par(mar=c(6,8.5,3,3))
jpeg("3_correlation_heatmap.png", width=900, height=750, units = "px") 

heatmapCor <- labeledHeatmap(Matrix = moduleTraitCor,
                             xLabels = colnames(design_df),
                             yLabels = names(MEs),
                             ySymbols = names(MEs),
                             colorLabels = FALSE,
                             colors = blueWhiteRed(50),
                             textMatrix = textMatrix,
                             setStdMargins = FALSE,
                             cex.text = 0.5,
                             cex.lab.y = 0.5,
                             cex.lab.x = 0.7,
                             zlim = c(-1,1),
                             main = paste("Module-trait relationships"))

dev.off()



# Module Analysis ---------------------------------------------------------

### Module Analysis ###

# Define variable containing the diagnosis column of design traits
diagnosis <- as.data.frame(design_df$HCA3)
names(diagnosis) <- "diagnosis"

# Define names of the modules
modNames <- substring(names(MEs), 3)

# Correlate gene module membership

geneModuleMembership <- as.data.frame(cor(norm_counts, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(dge$counts)))

names(geneModuleMembership) <- paste("MM", modNames, sep="");
names(MMPvalue) <- paste("p.MM", modNames, sep="");

# Correlate gene trait significance
geneTraitSignificance <- as.data.frame(cor(norm_counts, diagnosis, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(dge$counts)))

names(geneTraitSignificance) <- paste("GS.", names(diagnosis), sep="")
names(GSPvalue) <- paste("p.GS.", names(diagnosis), sep="")

### Select the module of interest to plot ###

# Order modules by their correlation to diagnosis
modOrder <- order(-abs(cor(MEs, diagnosis, use = "p")))

# Print names of top 3
modNames[modOrder[1:3]]
moduleChoice <- modNames[modOrder[1]]

# Plot MM and GS 
# Check that genes associated with trait are also most important in module

column <- match(moduleChoice, modNames)
moduleGenes <- netwkColours==moduleChoice

jpeg(paste0(moduleChoice, "_MM_GS_Scatter.jpeg"))
modScatter <- verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                                 abs(geneTraitSignificance[moduleGenes, 1]),
                                 xlab = paste("Module membership in", moduleChoice, "module"),
                                 ylab = "Gene significance for Ph+ALL",
                                 main = paste("Module membership vs. gene significance\n"),
                                 cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = moduleChoice) 

modScatter

dev.off()


# Annotate genes
geneInfoEnsembl <- data.frame(ensembl_gene_id = names(norm_counts))
geneInfoAnno <- left_join(geneInfoEnsembl, dge$genes, by="ensembl_gene_id")

# Add gene trait significance information to gene annotation
geneInfoOG <- data.frame(ensembl_gene_id = geneInfoAnno$ensembl_gene_id, 
                         external_gene_name = geneInfoAnno$external_gene_name,
                         entrezgene_id = geneInfoAnno$entrezgene_id,
                         module = netwkColours,
                         geneTraitSignificance,
                         GSPvalue)



# Add module membership information to the gene annotation dataframe and rename columns
for (i in 1:ncol(geneModuleMembership)){
  oldNames <- names(geneInfoOG)
  geneInfoOG <- data.frame(geneInfoOG, geneModuleMembership[, modOrder[i]],
                           MMPvalue[, modOrder[i]]);
  names(geneInfoOG) <- c(oldNames, paste("MM.", modNames[modOrder[i]], sep=""),
                         paste("p.MM.", modNames[modOrder[i]], sep=""))
}


# Order the dataframe by module color, then by gene trait significance
geneOrder <- order(geneInfoOG$module, -abs(geneInfoOG$GS.diagnosis));
geneInfo <- geneInfoOG[geneOrder, ]

write.csv(geneInfo, "geneInfo.csv")

# Values for module membership correlation
mmMod <- geneInfo %>% filter(module == moduleChoice)

write.csv(mmMod, paste0(moduleChoice, "_Genes.csv"))

# Explore Genes -----------------------------------------------------------


# Get the names of the genes in the interesting module/s

moduleChoice <- modNames[modOrder[1]] 

mmMod <- geneInfo %>% filter(module == moduleChoice,
                             abs(GS.diagnosis) > 0.8)

geneEnsembl <- mmMod$ensembl_gene_id
geneNames <- mmMod$external_gene_name




# Gene Ontology -----------------------------------------------------------


# Use significant results only
GOresults <- enrichGO(
  gene = mmMod$external_gene_name,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL", 
  universe = dge$genes$external_gene_name,
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05
) %>% clusterProfiler::simplify()

# Save results table
GOresultsTable <- GOresults@result

GOresultsTable %>%
  readr::write_csv(paste(moduleChoice, "_GOresults.csv"))


KEGGresults <- enrichKEGG(
  gene = mmMod$entrezgene_id,
  organism = "hsa",
  keyType = "ncbi-geneid", 
  universe = as.character(dge$genes$entrezgene_id),
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05
) 

KEGG_resultsTable <- KEGGresults@result

KEGG_resultsTable %>%
  readr::write_csv("KEGGresults.csv")


# Function to extract genes from GO categories
namesGO <- function(x){
  as.data.frame(strsplit(GOresultsTable$geneID[[x]], split = "/")) %>% 
    set_colnames("external_gene_name") %>% 
    left_join(dge$genes, by="external_gene_name")
}

# Set i as the gene set (row) of interest
i = 1
topGO <- namesGO(i)

# Gather expression data for those genes
genes_to_plot <- topGO$ensembl_gene_id
genes_to_plot <- mmMod$ensembl_gene_id

expression <- dge %>% 
  cpm(log = TRUE) %>%
  as.data.frame %>%
  rownames_to_column("ensembl_gene_id") %>%
  dplyr::filter(ensembl_gene_id %in% genes_to_plot) %>%
  melt() %>%
  set_colnames(c("ensembl_gene_id", "LB_ID", "logCPM"))

expression <- expression %>%
  left_join(dge$samples, by = "LB_ID")  %>%
  left_join(dge$genes, by = "ensembl_gene_id")

# Plot boxplots of the logCPM expression:
expressionPlotGO <- expression %>%
  ggplot(aes(x = col, y = logCPM, fill = HCA)) +
  geom_boxplot() + 
  facet_wrap( ~ external_gene_name, scales='free_y') +
  theme(aspect.ratio = 1) +
  ggtitle(paste0("LogCPM of ", GOresultsTable$Description[i], " genes"))

expressionPlotGO
