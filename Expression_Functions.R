# Create DGE object including all patients
getDGEList <- function(){
  
  print("Loading sample data")
  # Samples
  sample_anno <- read.csv("meta_data/sample_annotation.csv")
  #samples <- sample_anno %>% filter(exclude == "No")
  
  samples <- samples[order(samples$col),]
  rownames(samples) <- samples$LB_ID
  
  # Counts
  
  counts <- readRDS("Rdata/countData.rds")
  
  # Make the gene ID column the rownames, then drop the column
  row.names(counts) <- counts$gene_id
  geneIDcol <- which(colnames(counts) == "gene_id")
  counts <- counts[ ,-geneIDcol]
  
  # Remove NAs
  counts <- counts[complete.cases(counts), ]
  
  # Subset counts for only the required samples (identified by unique LB_ID)
  counts <- counts[ ,match(intersect(samples$LB_ID, colnames(counts)), colnames(counts))]
  
  # Reorder and then rename the counts data columns to match the samples list (by sampleID)
  samples$LB_ID <- as.character(samples$LB_ID)
  counts <- counts[,samples$LB_ID] 
  colnames(counts) <- samples$patientID
  
  
  # Gene annotation  
  
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
  
  print("Creating GE list object")
  # Create DGEList
  dge <- DGEList(counts = counts,
                 genes = genes,
                 samples = samplesALL,
                 remove.zeros = TRUE) %>% 
    calcNormFactors("TMM") %>% suppressMessages()
  
  saveRDS(dge, "GEobject.rds")
  print("GE object saved")
  
} 

# Boxplot showing the expression of a gene across all subtypes
plotGene <- function(geneName){
  
  geneEnsembl <- genes[which(genes$external_gene_name %in% geneName), "ensembl_gene_id"]
  genes_to_plot <- (geneEnsembl)
  
  dge <- readRDS("Expression/DGE/GEobject.rds")
  
  samplesA <- dge$samples
  logCPM <- cpm(dge, log=TRUE)
  batch <- as.factor(c(samplesA$ref))
  countsNorm <- removeBatchEffect(logCPM, batch=batch) 
  
  # Extract count data
  expression <- countsNorm %>%
    as.data.frame %>%
    rownames_to_column("ensembl_gene_id") %>%
    dplyr::filter(ensembl_gene_id %in% genes_to_plot) %>%
    melt() %>%
    set_colnames(c("ensembl_gene_id", "patientID", "logCPM")) %>% suppressMessages()
  
  expression <- expression %>%
    left_join(dge$samples, by = "patientID")  %>%
    left_join(dge$genes, by = "ensembl_gene_id")
  
  
  # Plot boxplots of the logCPM expression:
  expressionPlot <- expression %>%
    ggplot(aes(x = primary.subtype, y = logCPM, fill = primary.subtype)) +
    geom_boxplot() + 
    facet_wrap(.~external_gene_name, scales='free_y') +
    theme_classic() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, vjust = 1, hjust= 1, size=9)) +
    labs(x="Subtype") +
    ggtitle(paste0("Expression of ", geneName, " across subtypes of ALL"))
  
  print("Printing boxplots")
  
  return(expressionPlot)
}

# Boxplot showing the expression of multiple genes across all subtypes
plotPanel <- function(geneList){
  geneEnsembl <- genes[which(genes$external_gene_name %in% geneList), "ensembl_gene_id"]
  genes_to_plot <- (geneEnsembl)
  
  dge <- readRDS("Expression/DGE/GEobject.rds")
  samplesA <- dge$samples
  logCPM <- cpm(dge, log=TRUE)
  batch <- as.factor(c(samplesA$ref))
  countsNorm <- removeBatchEffect(logCPM, batch=batch) 
  
  # Extract count data
  expression <- countsNorm %>%
    as.data.frame %>%
    rownames_to_column("ensembl_gene_id") %>%
    dplyr::filter(ensembl_gene_id %in% genes_to_plot) %>%
    melt() %>%
    set_colnames(c("ensembl_gene_id", "patientID", "logCPM")) %>% suppressMessages()
  
  expression <- expression %>%
    left_join(dge$samples, by = "patientID")  %>%
    left_join(dge$genes, by = "ensembl_gene_id")
  
  
  # Plot boxplots of the logCPM expression:
  expressionPlot <- expression %>%
    ggplot(aes(x = col, y = logCPM, fill = col)) +
    geom_boxplot() + 
    facet_wrap(.~external_gene_name, scales='free_y') +
    theme_classic() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, vjust = 1, hjust= 1, size=8)) +
    labs(x="Subtype") +
    ggtitle(paste0("Expression of genes across subtypes of ALL"))
  
  print("Printing boxplots")
  
  return(expressionPlot)
}

# Plot the expression of a genes across all subtypes with dot points
plotGeneDots <- function(geneName){
  
  geneEnsembl <- genes[which(genes$external_gene_name %in% geneName), "ensembl_gene_id"]
  genes_to_plot <- (geneEnsembl)
  
  dge <- readRDS("Expression/DGE/GEobject.rds")
  
  samplesA <- dge$samples
  logCPM <- cpm(dge, log=TRUE)
  batch <- as.factor(c(samplesA$ref))
  countsNorm <- removeBatchEffect(logCPM, batch=batch) 
  
  # Extract count data
  expression <- countsNorm %>%
    as.data.frame %>%
    rownames_to_column("ensembl_gene_id") %>%
    dplyr::filter(ensembl_gene_id %in% genes_to_plot) %>%
    melt() %>%
    set_colnames(c("ensembl_gene_id", "patientID", "logCPM")) %>% suppressMessages()
  
  expression <- expression %>%
    left_join(dge$samples, by = "patientID")  %>%
    left_join(dge$genes, by = "ensembl_gene_id")
  
  
  # Plot boxplots of the logCPM expression:
  jitterPlot <- ggplot(expression, aes(x=col, y=logCPM, fill=col)) +  
    geom_point(aes(color=col), position = position_jitterdodge()) +
    xlab("Gene") +
    theme_classic() +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, vjust = 1, hjust= 1, size=9))+
    labs(x="Subtype") +
    ggtitle(paste0("Expression of ", geneName, " across subtypes of ALL"))
  
  print("Printing dot plots")
  
  return(jitterPlot)
}

# Get ensembl annotations for a list of gene (i.e. DGE results)
getAnnoList <- function(geneList){
  
  # For checking later
  entrezList <- geneList$entrezgene_id[1:25]
  genes_to_plot <- rownames(geneList)[1:25]
  
  # Remove NAs
  geneList <- geneList %>% filter(!is.na(geneList$entrezgene_id))
  geneList <- geneList[!duplicated(geneList$entrezgene_id),]
  row.names(geneList) <- geneList$entrezgene_id
  
  print("Connecting to database")
  # Get annotation database and list of attributes
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  # list <- listAttributes(ensembl, page="feature_page")
  
  print("Annotating list")
  # Annotate
  entAnno <- getBM(attributes = c("entrezgene_id","name_1006", "namespace_1003", "entrezgene_description"),
                   filters = "entrezgene_id",
                   values = entrezList,
                   mart = ensembl)
  
  # Join to find external gene name
  descript <- left_join(entAnno, genes, by="entrezgene_id")
  descript[,1] <- descript$external_gene_name
  
  return(descript)
  print("Saving dataframe in workspace")
  
}

# Get ensembl annotations for an individual gene
getAnno <- function(geneName){
  
  geneEntrez <- genes[which(genes$external_gene_name %in% geneName), "entrezgene_id"]
  
  print("Connecting to database")
  # Get annotation database and list of attributes
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  # list <- listAttributes(ensembl, page="feature_page")
  
  print("Getting annotation")
  
  # Annotate
  entAnno <- getBM(attributes = c("entrezgene_id","name_1006", "namespace_1003", "entrezgene_description"),
                   filters = "entrezgene_id",
                   values = geneEntrez,
                   mart = ensembl)
  
  # Join to find external gene name
  descript <- left_join(entAnno, genes, by="entrezgene_id")
  descript[,1] <- descript$external_gene_name
  
  return(descript)
  
}

# Extract gene annotation from GO sets 
namesGO <- function(x){
  as.data.frame(strsplit(GOresultsTable$geneID[[x]], split = "/")) %>% 
    set_colnames("external_gene_name") %>% 
    left_join(dge$genes, by="external_gene_name")
}

# Extract gene annotation from KEGG pathways
namesKEGG <- function(x){
  as.data.frame(strsplit(KEGG_gse_results$core_enrichment[[x]], split = "/")) %>% 
    set_colnames("entrezgene_id") %>% 
    left_join(genes, by="entrezgene_id")
}

# Run DGE on a subgroup vs rest of cohort
runDGE <- function(groupName){
  
  dir.create(paste0(groupName))
  
  setwd(paste0("Expression/DGE/ALL_DGE/", groupName, "/"))
  dge <- readRDS("Expression/DGE/ALL_DGE/DGE_ALL.rds")
  
  dgeMut <- dge
  dgeMut$samples <- dgeMut$samples %>% mutate(col = if_else(col== groupName, "Group", "Other_B_ALL"))
  
  # Designed to account for batch and PCA identified sources of variation
  design <- model.matrix(~0 + col + ageGp + source + ref + read_length, data = dgeMut$samples)
  
  print("Estimating mean variance (voom)")
  # Apply transformation and estimate mean variance
  voomData <- voom(dgeMut, design = design, plot = TRUE)
  
  # Define groups to compare
  contrasts <- makeContrasts(levels = colnames(design), Group_vs_Other = colGroup - colOther_B_ALL)
  
  
  ### Fit Linear Model ###
  print("Fitting linear model")
  fitModel <- lmFit(voomData, design) %>%
    contrasts.fit(contrasts) %>%
    treat(lfc = 1)
  
  
  ### Get Results and Adjust p-values ###
  
  results <- decideTests(fitModel, 
                         p.value = 0.05,
                         adjust.method = "fdr") 
  
  print(summary(results))
  
  # Find top DE genes
  allDEGs <- topTreat(fitModel, 
                      coef = "Group_vs_Other", 
                      number = Inf, 
                      adjust.method = "fdr") %>%
    as.data.frame() 
  
  allDEGs %>%
    readr::write_csv(paste0(groupName, "_DE_genes.csv"))
  
  
  # Identify significant results
  allDEGs <- allDEGs %>%
    dplyr::mutate(isSignificant = case_when(
      adj.P.Val < 0.05 & abs(logFC) > 1 ~ TRUE, 
      TRUE ~ FALSE 
    ))
  
  sigDEGs <- allDEGs %>%
    dplyr::filter(isSignificant == TRUE)
  
  sigDEGs %>%
    readr::write_csv(paste0(groupName, "_significant_DE_genes.csv"))
  
  print("List of DEGs saved")
  
  # Visualise Results -------------------------------------------------------
  
  
  ### Volcano plot ###
  
  volcanoPlot <- allDEGs %>%
    ggplot(aes(x = logFC, 
               y = -log10(P.Value),
               colour = isSignificant)) +
    geom_point(size = 1, alpha = 0.5) +
    scale_colour_manual(values = c("grey", "orange")) +
    ggtitle("Differential gene expression results") +
    labs(color = "Significant") +
    xlab("Fold Change (logFC)") +
    ylab("-log10 of p-Value") +
    theme_classic()
  
  volcanoPlot
  ggsave(paste0("6_", groupName, "_volcano_plot.png"), volcanoPlot) %>% suppressMessages()
  
  
  ### Bar Plot ###
  # Extract gene ensembl IDs
  sigDEGs <- sigDEGs %>% arrange(adj.P.Val)
  row.names(sigDEGs) <- sigDEGs$ensembl_gene_id
  genes_to_plot <- rownames(sigDEGs)[1:25]
  
  # Extract count data
  expression <- dgeMut %>% 
    cpm(log = TRUE) %>%
    as.data.frame %>%
    rownames_to_column("ensembl_gene_id") %>%
    dplyr::filter(ensembl_gene_id %in% genes_to_plot) %>%
    melt() %>%
    set_colnames(c("ensembl_gene_id", "sampleID", "logCPM")) %>% suppressMessages()
  
  expression <- expression %>%
    left_join(dgeMut$samples, by = "sampleID")  %>%
    left_join(dgeMut$genes, by = "ensembl_gene_id")
  
  
  # Plot boxplots of the logCPM expression:
  expressionPlot <- expression %>%
    ggplot(aes(x = col, y = logCPM, fill = col)) +
    geom_boxplot() + 
    facet_wrap( ~ external_gene_name, scales='free_y') +
    theme(aspect.ratio = 1) +
    scale_fill_manual(values= cbPalette, labels = c(groupName, "Other B-ALL")) +
    labs(x = "Group",
         fill = "Group") +
    ggtitle("LogCPM of top genes that show significant differential expression") %>% suppressMessages()
  
  expressionPlot 
  ggsave(paste0("7_", groupName, "_expression_plot_all.png"), expressionPlot, width= width, height = height) 
  
  ### Jitter Plot ###
  
  jitterPlot <- ggplot(expression, aes(x=external_gene_name, y=logCPM, color=col)) +  
    geom_point(aes(color=col), position = position_jitterdodge()) + 
    scale_color_manual(values= cbPalette, labels = c(groupName, "Other B-ALL")) +
    labs(x = "Gene",
         color = "Group") +
    theme(axis.text.x = element_text(angle=45, vjust = 1, hjust= 1, size=9))+
    ggtitle("LogCPM of top genes that show significant differential expression")
  
  jitterPlot
  ggsave(paste0("8_", groupName, "_expression_dot_plot.png"), jitterPlot)
  
  # list_sigs <- sigDEGs[2]
  # list_sigs[1:25,]
  print("Expression plots saved. Plotting protein coding genes only.")
  
  # Protein-Coding Genes Only -----------------------------------------------
  
  
  # View gene biotypes
  ggplot(sigDEGs, aes(x=gene_biotype)) + 
    geom_bar(aes(fill=gene_biotype)) +
    theme(legend.position = "none", axis.text.x = element_text(angle=45, vjust = 1, hjust= 1, size=9)) +
    labs(x="Gene Biotype",
         y = "Count") +
    ggtitle("Differentially expressed gene biotypes")
  ggsave(paste0("9a_", groupName, "_gene_biotype.png"))
  
  # Filter by protein coding genes
  proteinDEGs <- sigDEGs %>% filter(gene_biotype == "protein_coding")
  allProteinDEGs <- allDEGs %>% filter(gene_biotype == "protein_coding")
  
  # list_PCsigs <- proteinDEGs[2]
  # list_PCsigs[1:25,]
  
  ### Volcano plot ###
  
  volcanoPCPlot <- allProteinDEGs %>%
    ggplot(aes(x = logFC, 
               y = -log10(P.Value),
               colour = isSignificant)) +
    geom_point(size = 1, alpha = 0.5) +
    scale_colour_manual(values = c("grey", "orange")) +
    ggtitle("Differential gene expression results") +
    labs(color = "Significant") +
    xlab("Fold Change (logFC)") +
    ylab("-log10 of p-Value") +
    theme_classic()
  
  volcanoPCPlot
  ggsave(paste0("9_", groupName, "_volcano_plot_protein_coding.png"), volcanoPCPlot, width= width, height = height)
  
  
  ### Bar Plot ###
  
  # Extract ensembl IDs
  proteinDEGs <- proteinDEGs %>% arrange(adj.P.Val)
  row.names(proteinDEGs) <- proteinDEGs$ensembl_gene_id
  genes_to_plot <- rownames(proteinDEGs)[1:25]
  
  # Extract count data
  expressionPC <- dgeMut %>% 
    cpm(log = TRUE) %>%
    as.data.frame %>%
    rownames_to_column("ensembl_gene_id") %>%
    dplyr::filter(ensembl_gene_id %in% genes_to_plot) %>%
    melt() %>%
    set_colnames(c("ensembl_gene_id", "sampleID", "logCPM"))  %>% suppressMessages()
  
  expressionPC <- expressionPC %>%
    left_join(dgeMut$samples, by = "sampleID")  %>%
    left_join(dgeMut$genes, by = "ensembl_gene_id")
  
  
  # Plot boxplots of the logCPM expression:
  expressionPCPlot <- expressionPC %>%
    ggplot(aes(x = col, y = logCPM, fill = col)) +
    geom_boxplot() + 
    facet_wrap(.~external_gene_name, scales='free_y') +
    theme(aspect.ratio = 1) +
    labs(x="Group",
         fill= "Group") +
    scale_fill_manual(values= cbPalette, labels = c(groupName, "Other B-ALL")) +
    ggtitle("LogCPM of top protein-coding genes that show significant differential expression")
  
  expressionPCPlot 
  ggsave(paste0("10_", groupName, "_expression_plot_protein_coding.png"), expressionPCPlot, width= width, height = height)
  
  ### Jitter Plot ###
  
  jitterPCPlot <- ggplot(expressionPC, aes(x=external_gene_name, y=logCPM, color=col)) +  
    geom_point(aes(color=col), position = position_jitterdodge()) + 
    scale_color_manual(values= cbPalette, labels = c(groupName, "Other B-ALL")) +
    labs(x="Gene",
         fill= "Group",
         color= "Group") +
    theme(axis.text.x = element_text(angle=45, vjust = 1, hjust= 1, size=9)) +
    ggtitle("LogCPM of top protein-coding genes that show significant differential expression")
  
  jitterPCPlot
  ggsave(paste0("11_", groupName, "expression_dot_plot_protein_coding.png"), jitterPCPlot)
  
  print("Protein expression plots saved. Analysising gene sets.")
  # Gene Set Over-Representation --------------------------------------------------
  
  ### Gene Set Over Representation ###
  
  # Use significant results only
  GOresults <- enrichGO(
    gene = sigDEGs$external_gene_name,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL", 
    universe = dgeMut$genes$external_gene_name,
    pAdjustMethod = "fdr",
    pvalueCutoff = 0.05
  ) %>% simplify
  
  # Save results table
  GOresultsTable <- GOresults@result
  
  GOresultsTable %>%
    readr::write_csv(paste0(groupName, "_GOresults.csv"))
  
  ### Plot GO Results ###
  
  # Network plot
  GOnetwork_plot <- GOresults %>% 
    cnetplot(cex_label_gene = 0.6, 
             showCategory = 3)
  
  GOnetwork_plot
  ggsave(paste0(groupName,"_GOnetwork_plot.png"), GOnetwork_plot)
  
  # Dot plot
  dotplot(GOresults, showCategory=50) + scale_color_gradient(low="purple",high="orange") %>% suppressMessages()
  ggsave(paste0("13", groupName, "_GO_dot_plot.png"), width = 10, height = 15)
  
  # Expression plots
  # Create function to extract gene annotation from GO gene sets
  namesGO <- function(x){
    as.data.frame(strsplit(GOresultsTable$geneID[[x]], split = "/")) %>% 
      set_colnames("external_gene_name") %>% 
      left_join(genes, by="external_gene_name")
  }
  
  # Set i as the gene set (row) of interest
  i = 1
  topGO <- namesGO(i)
  
  # Gather expression data for those genes
  genes_to_plot <- topGO$ensembl_gene_id[1:10]
  
  expression <- dgeMut %>% 
    cpm(log = TRUE) %>%
    as.data.frame %>%
    rownames_to_column("ensembl_gene_id") %>%
    dplyr::filter(ensembl_gene_id %in% genes_to_plot) %>%
    melt() %>%
    set_colnames(c("ensembl_gene_id", "sampleID", "logCPM")) %>% suppressMessages()
  
  expression <- expression %>%
    left_join(dgeMut$samples, by = "sampleID")  %>%
    left_join(dgeMut$genes, by = "ensembl_gene_id")
  
  # Plot boxplots of the logCPM expression:
  expressionPlotGO <- expression %>%
    ggplot(aes(x = col, y = logCPM, fill = col)) +
    geom_boxplot() + 
    facet_wrap( ~ external_gene_name, scales='free_y') +
    theme(aspect.ratio = 1) +
    scale_fill_manual(values= cbPalette, labels = c(groupName, "Other B-ALL")) +
    labs(x= "Group",
         fill= "Group") +
    ggtitle(paste0("LogCPM of ", GOresultsTable$Description[i], " genes"))
  
  expressionPlotGO
  
  ggsave(paste0("14_", groupName, "_expression_plot_GO.png"), expressionPlotGO, width= width, height = height)
  
  print("ORA done")
  # Gene Set Enrichment --------------------------------------------------------------------
  
  ### Gene Set Enrichment Analysis ###
  
  # Create a list of all genes and order by fold change
  gseaList <- allDEGs[!duplicated(allDEGs$entrezgene_id),] %>% na.omit()
  gseaList <- gseaList[order(gseaList$logFC, decreasing = TRUE),]
  listGSEA <- gseaList$logFC
  names(listGSEA) <- gseaList$entrezgene_id
  
  # Run GSEA
  KEGG_gse <- gseKEGG(geneList = listGSEA,
                      organism = "human",
                      minGSSize = 50,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "fdr")
  
  KEGG_gse_results <- KEGG_gse@result
  
  # Save results
  KEGG_gse_results %>%
    readr::write_csv(paste0(groupName, "_KEGGresults.csv"))
  
  ### Plot Results ###
  
  # Dot plot
  dotplot(KEGG_gse, color="pvalue") + scale_color_gradient(low="purple",high="orange") %>% suppressMessages()
  ggsave(paste0("15_", groupName, "_KEGG_dot_plot.png"), width=10, height = 15)
  
  # Ridge plot of fold change distribution
  ridgeplot(KEGG_gse, showCategory = 20) + scale_fill_gradient(low="purple",high="orange") +
    labs(x=  "Enrichment Distribution") %>% suppressMessages()
  ggsave(paste0("16_", groupName, "_KEGG_ridge_plot.png"), width= 15, height=15)
  
  # Network plot
  # KEGGnetwork <- KEGG_gse %>% 
  #   cnetplot(cex_label_gene = 0.4, 
  #            showCategory = 6)
  # 
  # KEGGnetwork
  
  
  
  # Plot KEGG gene expression -----------------------------------------------
  
  # Create function to extract gene annotation from KEGG pathway sets
  namesKEGG <- function(x){
    as.data.frame(strsplit(KEGG_gse_results$core_enrichment[[x]], split = "/")) %>% 
      set_colnames("entrezgene_id") %>% 
      left_join(genes, by="entrezgene_id")
  }
  
  # Extract ensembl IDs
  genes$entrezgene_id <- as.character(genes$entrezgene_id)
  a <- 1
  topKEGG <- namesKEGG(a)
  
  genes_to_plot <- topKEGG$ensembl_gene_id[1:10]
  
  # Extract counts for those genes
  expression <- dgeMut %>% 
    cpm(log = TRUE) %>%
    as.data.frame %>%
    rownames_to_column("ensembl_gene_id") %>%
    dplyr::filter(ensembl_gene_id %in% genes_to_plot) %>%
    melt() %>%
    set_colnames(c("ensembl_gene_id", "sampleID", "logCPM"))
  
  expression <- expression %>%
    left_join(dgeMut$samples, by = "sampleID")  %>%
    left_join(dgeMut$genes, by = "ensembl_gene_id")
  
  
  # Plot boxplots of the logCPM expression:
  expressionPlotKEGG <- expression %>%
    ggplot(aes(x = col, y = logCPM, fill = col)) +
    geom_boxplot() + 
    facet_wrap( ~ external_gene_name, scales='free_y') +
    theme(aspect.ratio = 1) +
    scale_fill_manual(values= cbPalette, labels = c(groupName, "Other B-ALL")) +
    labs(x="Group",
         fill= "Group") +
    ggtitle(paste0("LogCPM of ", KEGG_gse_results$Description[a], " genes"))
  
  expressionPlotKEGG 
  
  ggsave(paste0("17_", groupName, "_expression_plot_kegg.png"), expressionPlotKEGG, width= width, height = height)
  
  print("Enrichment done. DGE analysis complete")
  
  topSig <-sigDEGs[c(1:50), 1] %>% as.data.frame() %>% set_colnames(paste0(groupName)) %>% readr::write_csv(paste0(groupName, "_topSig.csv"))
  topPC <- proteinDEGs[c(1:50),1] %>% as_data_frame() %>% set_colnames(paste0(groupName)) %>% readr::write_csv(paste0(groupName, "_topPC.csv"))
  
  setwd(paste0("Expression/DGE/ALL_DGE/"))
} 



