CNV\_Bench
================
Madyson Colton
22/04/2022

# RNAseqCNV

Implementing the RNAseqCNV caller into the SAHMRI bioinformatics
pipeline following variant calling. This tool can be used to determine
aneuploidy in samples missing karyotype data.

Before implementing into our pipeline, the accuracy of the caller was
assessed.

## Installing on the HPC

R is required to run RNAseqCNV. An virtual environment was created to
install the latest version of R and related packages.

``` bash
 mamba create --name Renv r-base r-essentials
 mamba activate Renv
 mamba install -c bioconda bioconductor-biocinstaller
```

To install RNAseqCNV and dependencies, run the installCNV.R script. The
RNAseqCNV package is stored on GitHub, so the devtools package is
required for installation.

``` bash
Rscript installCNV.R
```

``` r
# Set default repo
local({r <- getOption("repos")
       r["CRAN"] <- "http://cran.r-project.org"
       options(repos=r)
})


install.packages('devtools')

library(devtools)

devtools::install_github(repo = "honzee/RNAseqCNV")
```

Some dependent packages may have to be installed manually. For me,
ggpubr was causing the installation to fail. If this happens, install
the dependecies via mamba/conda.

``` bash
mamba install -c conda-forge r-ggpubr
```

## Locate and Tidy Data

To determine CNVs, the tool requires both sequencing counts and SNV
frequencies.

The counts data should have only 2 columns - the ensembl IDs and the
read counts - and no header. The SNV information can be provided as a
.vcf (or vcf.gz) file following the GATK pipeline. This file should
contain headers; RNAseqCNV will require the chromosome number, start
position, read depth and mutant allele frequency.

The RNAseqCNV wrapper then requires the following:

-   A config file
    -   Directory to store the output files
    -   Directory containing the counts files
    -   Directory containing the vcf files
-   A metadata file, containing the following 3 columns with no headers:
    -   Sample ID
    -   Path to the sample’s counts file
    -   Path to the sample’s vcf file

### Config file

An example of the config file is given below. The variable names must
not be changed. The out\_dir path should be the same as the outDir for
the variant calling pipeline, and should end in a “/”; the count\_dir
and snv\_dir paths should not have a “/” at the end.

``` bash
out_dir = "~/CNV/outDir/"
count_dir = "~/CNV"
snv_dir = "/cancer/storage/results"
```

### Counts data

The counts data need to be converted from an “.rData” file to a “.csv”
file with two columns. As all our previous samples have been stored in a
“.RDS” object in the required format, we can extract the count data from
the object by running extractCountsAll.R. As we want to run the
RNAseqCNV tool on all our previous samples, the summary log of our files
is required.

``` r
# Default repo
local({r <- getOption("repos")
r["CRAN"] <- "http://cran.r-project.org"
options(repos=r)
})

library(dplyr)

# Load sample summary file
sampleSumm <- read.csv("/data/bioinformatics/leukaemia_bank/meta_data/summary_ALL_samples.csv")

# Load count data
counts <- readRDS("/data/bioinformatics/leukaemia_bank/Rdata/countData.rds")

# For all samples
sampleNames <- sampleSumm$variant_file_ID

# Directory to save the counts files
countOut <- "~/CNV/counts/"


# Create list for failed samples
notRun <- list()

# Extract counts for all sampleNames
for(i in 1:length(sampleNames)){
  
  sampleID <- sampleNames[i]
  
  bankID <- sampleSumm[which(sampleSumm$variant_file_ID == sampleID),4]
  
  if(bankID %in% colnames(counts)){
    sampleCounts <- counts %>% 
      dplyr::select(c("gene_id", bankID))
   
     names(sampleCounts) <- NULL
    
    write.table(sampleCounts, 
                paste0(countOut, bankID, "_counts.csv"), 
                col.names = FALSE, 
                sep=",", 
                row.names = FALSE, 
                quote = FALSE, 
                eol = "\r")
    
  print(paste0(bankID, " done."))
  
    }else{
    notRun[[sampleID]] <- append(bankID, "Not found")
    print(paste0(bankID, " not found"))
  }
}

# Save failed counts in a dataframe
notRun <- notRun %>% 
  as.data.frame() %>% 
  t() 

colnames(notRun) <- c("sampleID", "status")

write.table(notRun, 
            paste0(countOut, "notRun.csv"), 
            row.names = FALSE, 
            sep=",", 
            quote = FALSE, 
            eol = "\r")
```

All samples were found.

### Metadata file

#### 2020 Samples

To run the RNAseqCNV caller on all of our samples, a custom metadata
sheet was generated. Our more recent samples have SNV data stored as
“.vcf.gz” zipped file, whereas older samples are simply “.vcf”. The
location of the sequencing counts for older samples are also stored in
different locations with an esoteric file naming system, so initially
only the 2020 samples were included.

``` r
### Metadata table for benchmarking previous samples ###

# Default repo
local({r <- getOption("repos")
r["CRAN"] <- "http://cran.r-project.org"
options(repos=r)
})

library(stringr)
library(dplyr)

# To run on all samples, use the sample summary variant file IDs
metaD <- read.csv("/data/bioinformatics/leukaemia_bank/meta_data/summary_ALL_samples.csv")

# Only use the 2020 samples
# Paste file locations
metaD <- metaD %>% 
  mutate(variant_loc = if_else(str_detect(variant_loc, "2020") == TRUE, 
                                                paste0(variant_loc, "/Variants/", variant_file_ID, ".raw.vcf.gz"), 
                                                paste0("")))

metaD <- metaD[metaD$variant_loc != "",]

# Add a row for counts data
metaD <- metaD %>% 
  mutate(counts_loc = paste0("counts/", LB_ID, "_counts.csv"))

metaTable <- metaD %>% 
  select(c("variant_file_ID", "counts_loc", "variant_loc"))


# Save table in working directory
write.table(metaTable, 
            paste0("~/CNV/metadata.csv"), 
            col.names = FALSE, 
            sep=",", 
            row.names = FALSE, 
            quote = FALSE, 
            eol = "\r")
```

To be able to find the /cancer/storage directories containing the vcf
file, the RNAseqCNV tool will have to be initiated from the root
directory. Therefore, snv\_dir location in the config file is set as
“/cancer/storage/results”, so this has to be removed from the metadata
variant location file names.

``` bash
# -i makes changes in place
# s/ substitute command
# /g global replacement, not just first instance
# \ to escape the / characters in the string

sed -i 's/\/cancer\/storage\/results\/ALL/ALL/g' metadata.csv
```

#### Older Samples

Since the older samples were stored under unusual naming systems, files
ending in “.raw.vcf” were located using the Sys.glob() function in R.

``` r
### Metadata table for benchmarking previous samples ###

# Default repo
local({r <- getOption("repos")
r["CRAN"] <- "http://cran.r-project.org"
options(repos=r)
})

### Load Libraries ###
library(stringr)
library(dplyr)

# Load the sample summary sheet containing the parent directories
oldMeta <- read.csv("/data/bioinformatics/leukaemia_bank/meta_data/summary_ALL_samples.csv")

# Select only the older samples and then the columns of interest
oldMeta <- oldMeta %>% 
  filter(str_detect(variant_loc, "2020") == FALSE) %>% 
  select(c("patientID", "variant_file_ID","LB_ID", "variant_loc"))

# Add the "/Variants" directory onto the end of the old file path
oldMeta <- oldMeta %>% 
  mutate(variant_loc = paste0(variant_loc, "/Variants"))

# Change the variant_loc column name for clarity
# Create an empty new variant_loc column
colnames(oldMeta) <- c("patientID", "variant_file_ID","LB_ID", "old_loc")

oldMeta$variant_loc <- ""

# For each sample, paste the name of the file ending in ".raw.vcf"
# Set condition to handle missing/duplicate files 

for(i in 1:nrow(oldMeta)){
  path <- oldMeta[i,]$old_loc
  variant_loc <- Sys.glob(file.path(path, '*.raw.vcf'))
  if(length(variant_loc) != 1){
    oldMeta[i,5] <- "NA"
  }else{
  oldMeta[i,5] <- paste(variant_loc)
}
}

# Create a metaD table also containg the counts paths (from extractCounts.R)
metaD <- oldMeta %>% 
  mutate(counts_loc = paste0("counts/", LB_ID, "_counts.csv"))

# Select columns in required order for RNAseqCNV

metaTable <- metaD %>% 
  select(c("variant_file_ID", "counts_loc", "variant_loc"))


# Save table in working directory
write.table(metaTable, 
            paste0("~/CNV/old_metadata.csv"), 
            col.names = FALSE, 
            sep=",", 
            row.names = FALSE, 
            quote = FALSE, 
            eol = "\r")
print("saved")
```

The script contains a function to generate “NA” if there are missing or
duplicate .vcf files in that directory. I initially forgot to remove
these from the metadata column, so ran a script to remove all rows
containing NAs. The variant file location names were also edited to
remove “/cancer/storage/results” as with the 2020 samples above.

``` r
# Default repo
local({r <- getOption("repos")
r["CRAN"] <- "http://cran.r-project.org"
options(repos=r)
})

library(dplyr)

# Load metadata table
metaTable <- read.csv("~/CNV/old_metadata.csv")

# Remove NAs
metaTable <- metaTable[complete.cases(metaTable),]

# Save table
write.table(metaTable,
            paste0("~/CNV/old_metadata.csv"),
            col.names = FALSE,
            sep=",",
            row.names = FALSE,
            quote = FALSE,
            eol = "\r")
print("saved")
```

## Running RNAseqCNV

To run the tool, the path to the config and metadata files should be
provided and specified within the RNAseqCNV\_wrapper() function. The
format of the variant file should also be specified (“vcf”) and the
genome version (typically “hg19” in our pipeline). As this analysis is
only interested in the called chromosome number, arm level figures are
not required. Choosing “FALSE” for this parameter will greatly reduce
run time.

``` r
## Default repo
local({r <- getOption("repos")
       r["CRAN"] <- "http://cran.r-project.org"
       options(repos=r)
})

library(RNAseqCNV)

configPath <- "~/CNV/config"
metaPath <- paste0("~/CNV/metadata.csv")


RNAseqCNV_wrapper(config = configPath,
                  metadata = metaPath,
                  snv_format = "vcf",
                  genome_version = "hg19",
                  arm_lvl = FALSE)
```

The resulting estimation table was copied back to my local computer for
analysis.

``` bash
scp madyson.colton@edp-prd-lin-nfs01:~/CNV/outDir/estimation_table.tsv .
```

## Assessment of results

To assess the accuracy of the RNAseqCNV caller, the results were
compared to our known karyotype data and subtyping.

Only good quality, diagnositic B-ALL samples were included in the
benchmarking analysis. Patients who were unclassified and grouped with
the “B-Other” subtype were also excluded as these are believed to have
been unclassified due to missing karyotype data.

As the RNAseqCNV caller was ran twice (on the 2020 samples and then on
the older samples), the results had to be combined into a single .csv
file. This file was used in further analysis.

``` r
# Load both results tables
cnvTable2020 <- read.table("C:/Users/madyson/OneDrive - South Australian Health and Medical Research Institute Ltd/LBank/Projects/CNVs/estimation_table2020.tsv", sep = "\t", header = T)
oldCNVTable <- read.table("C:/Users/madyson/OneDrive - South Australian Health and Medical Research Institute Ltd/LBank/Projects/CNVs/estimation_table.tsv", sep = "\t", header = T)

# Merge CNV data
cnvTable <- rbind(cnvTable2020, oldCNVTable)

# Save
write.csv(cnvTable, "combined_table.csv")
```

``` r
### Load Libraries ###

library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(stringr)
library(ggplot2)

cbPalette <- c("#0072B2", "#E69F00", "#009E73", "#F0E442", "#56B4E9", "#D55E00", "#CC79A7")

# Load table and sample anno

cnvTable <- read.table("C:/Users/madyson/OneDrive - South Australian Health and Medical Research Institute Ltd/LBank/Projects/CNVs/combined_table.csv", sep = ",", header = T)
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

# Remove low quality samples
# Remove "others" as some of these were unclassfied due to missing cytogenetics
# Select only B-ALL and Diagnostic samples
noCyto <- allData %>% filter(exclude != "Yes",
                              !is.na(col),
                              Diagnosis.x %in% c("B-ALL", "Ph+ALL"),
                              sample.y == "Dx")
```

### Classification of patients:

-   Normal = 46 chromsomes
-   Hyperdiploid &gt; 46 chromosomes
-   Hypodiploid &lt; 46 chromosomes

In total, 124 of the 374 samples had recorded cytogenetics and were
classified into three groups based on the number of chromosomes
described above. For the samples without this information, the research
groups’ classification information was used: patients who were not
classified as Hyperdiploid or Hypodiploid were considered to have
“normal” ploidy.

| Cyto\_provided |   n |
|:---------------|----:|
| partial        |  49 |
| Y              | 118 |
| NA             | 189 |

Number of samples with cytogenetics

### Plot Accuracy

Patients were classified into groups as described, and then a new column
was created to reclassify samples based on the results of RNAseqCNV. The
classifications were compared.

``` r
# Create column to classify cnv output
# Create column for karyype classification
noCyto <- noCyto %>% 
  mutate(cnvCol = if_else(chrom_n == "46", "Normal", "")) %>% 
  mutate(cnvCol = if_else(chrom_n > "46", "Hyperdiploid", cnvCol)) %>%
  mutate(cnvCol = if_else(chrom_n < "46", "Hypodiploid", cnvCol)) %>%
  mutate(oldCol = if_else(col %in% c("Hyperdiploid", "Hypodiploid"), col, "Normal")) %>%
  mutate(oldCol = if_else(karyo_n == "46", "Normal", oldCol)) %>% 
  mutate(oldCol = if_else(karyo_n > "46", "Hyperdiploid", oldCol)) %>%
  mutate(oldCol = if_else(karyo_n %in% c(1:45), "Hypodiploid", oldCol))
 
noCyto <- noCyto %>% 
  mutate(correctBench = if_else(cnvCol == oldCol, "Correct", "Incorrect"))

# Create dataframe for plotting
noCytoPerc <- noCyto %>% 
  count(oldCol, correctBench) %>% 
  group_by(oldCol) %>% 
  mutate(label_y = n/sum(n)) 


noCytoPlot <- noCytoPerc %>% 
  ggplot(aes(x=oldCol, y=label_y, fill=correctBench)) +
  geom_col(position="stack") +
  geom_text(aes(label=paste(paste0(round(label_y*100, digits=0), "%"), paste0("(", n, ")"))), 
            position=position_stack(vjust = 0.5),
            colour="white",
            size=3,
            fontface="bold") +
  scale_fill_manual(values= cbPalette) +
  theme_classic() +
  labs(x = "Ploidy",
       y = "Percentage of Samples",
       fill = "RNAseqCNV Result") +
  ggtitle("Benchmarking of RNAseqCNV Tool on All Samples")
```

![](CNV_Bench_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

Overall, RNAseqCNV matched our classification in 83% of cases.

| correctBench |   n | Percentage |
|:-------------|----:|-----------:|
| Correct      | 293 |   82.30337 |
| Incorrect    |  63 |   17.69663 |

Accuracy of the RNAseqCNV Tool

### Comparison with Karyotype Data

The results were similar when looking only at the patients with known
karyotypes, with an overall accuracy of 84%.

``` r
cytoData <- allData %>% filter(exclude != "Yes",
                              !is.na(col),
                              Diagnosis.x %in% c("B-ALL", "Ph+ALL"),
                              sample.y == "Dx",
                              karyo_n != 0)

# Create column to classify cnv output
# Create column for karyype classification
cytoData <- cytoData %>% 
  mutate(cnvCol = if_else(chrom_n == "46", "Normal", "")) %>% 
  mutate(cnvCol = if_else(chrom_n > "46", "Hyperdiploid", cnvCol)) %>%
  mutate(cnvCol = if_else(chrom_n < "46", "Hypodiploid", cnvCol)) %>%
  mutate(oldCol = if_else(karyo_n == "46", "Normal", "")) %>% 
  mutate(oldCol = if_else(karyo_n > "46", "Hyperdiploid", oldCol)) %>%
  mutate(oldCol = if_else(karyo_n < "46", "Hypodiploid", oldCol))

cytoData <- cytoData %>% 
  select(c("patientID.x", "LB_ID", "gender", "ref", "key_alt", "primary.subtype", "Karyotype", "karyo_n", "chrom_n", "alterations", "cnvCol", "oldCol"))

# Create column comparing classification
cytoData <- cytoData %>% 
  mutate(correctBench = if_else(cnvCol == oldCol, "Correct", "Incorrect"))

# Create dataframe for plotting
cytoPerc <- cytoData %>% 
  count(oldCol, correctBench) %>% 
  group_by(oldCol) %>% 
  mutate(label_y = n/sum(n)) 

cytoPlot <- cytoPerc %>% 
  ggplot(aes(x=oldCol, y=label_y, fill=correctBench)) +
  geom_col(position="stack") +
  geom_text(aes(label=paste(paste0(round(label_y*100, digits=0), "%"), paste0("(", n, ")"))), 
            position=position_stack(vjust = 0.5),
            colour="white",
            size=3,
            fontface="bold") +
  scale_fill_manual(values= cbPalette) +
  theme_classic() +
  labs(x = "Ploidy",
       y = "Percentage of Samples",
       fill = "RNAseqCNV Result") +
  ggtitle("Benchmarking of RNAseqCNV Tool on Samples With Known Karyotypes")
```

![](CNV_Bench_files/figure-gfm/cytoPercPlot-1.png)<!-- -->

#### Accuracy of the RNAseqCNV Tool in Matching Known Karyotype

| correctBench |   n | Percentage |
|:-------------|----:|-----------:|
| Correct      |  97 |   82.20339 |
| Incorrect    |  21 |   17.79661 |

In total, 20 samples were incorrectly classified. Most of these (n=18)
had 45-47 chromosomes, and the rest were masked hypodiploid (n=2). One
of the incorrect samples had deletion of the Y chromosome as the only
alteration; this was misclassified since the RNAseqCNV tool does not
evaluate Y chromosome copy number changes so the sample was therefore
classified as having no alterations by the RNAseqCNV tool (CHI\_0802,
DUX4).

In many of these cases, RNAseqCNV correctly called the alterations but
with low confidence, so the final chromosome number was given as 46. Of
the incorrect samples, only 4 had aneuploidy as their primary clinically
relevant alteration (Masked Hypodiploid = 2, Hyperdiploid = 1,
Hypodiploid = 1), and one of these had an ambiguous cytogenetic report
that classified the patient as having 46 chromosomes but reported
Hypodiploid losses (ADII\_0044: 46,XY by G-banding. Hypodiploid loss of
3,13,15,16,17,19,20 by SNParray).

#### Incorrect RNAseqCNV Classification

| patientID.x         | LB\_ID                | gender | ref    | key\_alt     | primary.subtype    | Karyotype                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           | karyo\_n | chrom\_n | alterations                                                                                    | cnvCol       | oldCol       | correctBench |
|:--------------------|:----------------------|:-------|:-------|:-------------|:-------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------:|---------:|:-----------------------------------------------------------------------------------------------|:-------------|:-------------|:-------------|
| AYA\_0726           | AYAI\_0726.201362     | female | GRCh37 | ETV6:RUNX1   | ETV6-RUNX1         | “45,X,-X,del(6)(q1?3q2?1)\[14\]/88\~99,del(6)x2,inc\[cp4\]/46,XX\[5\].nuc ish(D6Z1x2,SEC63x1,MYBx2)\[131/220\]/(D6Z1x4,SEC63x2,MYBx4)\[62/220\],(ABL1x3,BCRx4)\[49/200\],(ETV6,RUNX1)x3(ETV6 con RUNX1x2)\[143/220\]/(ETV6x5,RUNX1x6)(ETV6 con RUNX1x4)\[43/220\]”                                                                                                                                                                                                                                                                                                                                                                                                                  |       45 |       46 | none                                                                                           | Normal       | Hypodiploid  | Incorrect    |
| AYA\_0780           | AYAII\_0780.201504    | male   | GRCh37 | PAX5:ETV6    | PAX5alt            | “45,XY,dic(9;12)(p13;p13)\[3\]/46,sl,+8\[1\]/46,XY\[16\].nuc ish(CRLF2x2)\[200\],(ABL2x2)\[200\],(PDGFRBx2)\[200\],(RUNX1T1x3,RUNX1x2)\[71/200\],(CDKN2Ax0,D9Z1x2)\[111/200\]/(CDKN2Ax1,D9Z1x2)\[44/200\],(ABL1x2)\[200\],(ABL1,BCR)x2\[200\],(KMT2Ax2)\[200\],(5’ETV6x1,3’ETV6x2)(5’ETV6 con3’ETV6x1)\[110/200\]/(ETV6x1)\[62/200\],(ETV6x1,ETV6 dimx1,RUNX1x2)\[112/200\]”                                                                                                                                                                                                                                                                                                        |       45 |       46 | ?9p                                                                                            | Normal       | Hypodiploid  | Incorrect    |
| CHI\_0802           | CHI\_0802.201571      | male   | GRCh37 | None         | DUX4               | “45,X,-Y”                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |       45 |       46 | none                                                                                           | Normal       | Hypodiploid  | Incorrect    |
| AYA\_0814           | AYAII\_0814.201699    | male   | GRCh37 | P2RY8:CRLF2  | Ph-like(CRLF2)     | “46,XY,del(2)(p?13p?12),del(12)(P13p11.2)\[1\]/47,idem+10\[5\]/47,idem,+6\[3\]/46,XY\[2\]”                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |       46 |       47 | 10+                                                                                            | Hyperdiploid | Normal       | Incorrect    |
| ADU\_0821           | ADI\_0821.201574      | male   | GRCh37 | BCR:ABL1     | Ph                 | “46,XY,t(9;22)(q34;q11.2)\[1\]/47,idem,del(3)(p21),+der(9)t(9;22),r(14)\[10\]/46,XY\[29\]”                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |       46 |       47 | 9p+, ?9q, 22+                                                                                  | Hyperdiploid | Normal       | Incorrect    |
| AYA\_0827           | AYAI\_0827.201627     | male   | GRCh37 | PAX5:JAK2    | Ph-like(JAK2/EPOR) | “47,XY,der(14;21)(q\`0q10),+2mar\[6\]/46,XY\[16\]”                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |       47 |       46 | ?21+                                                                                           | Normal       | Hyperdiploid | Incorrect    |
| CHI\_0233           | CHI\_0233             | male   | hg19   | BCR:ABL1     | Ph                 | “45,XY,t(2;7)(p11.2;p1?2),dic(9;20)(p13;q11.2-13.1), t(9;22)(q34;q11.2)\[10\]/46,XY\[10\].nuc ish(CRLF2x2)\[200\],(ETV6,RUNX1)x2\[200\],(TCF3x2)\[200\]”                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |       45 |       46 | 9p-, 20q-                                                                                      | Normal       | Hypodiploid  | Incorrect    |
| AYAII\_0396         | AYAII\_0396           | male   | hg19   | BCR:ABL1     | Ph                 | “46-47 XY, +X, der(7;9),t(9;22),+der(22)t(9;22)\[cp16\]”                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |       46 |       47 | 7p-, 9p-, X+                                                                                   | Hyperdiploid | Normal       | Incorrect    |
| AYAII\_0128         | AYAII\_0128           | female | hg19   | Hyperdiploid | Masked Hypodiploid | “61-71,XXX,-X,+1,-3,-4,-4,+6,-7,-8,-9,+10,+11,+12, -13,-14,-15,-16,-17,-17,+19,-20,+21,+22,+mar1, +mar2,+mar3\[cp6\]/46,XX\[14\].nuc ish 9q34(ABLX2),22q11(BCRX4)”                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |       61 |       35 | 3-, 4-, 5-, 7-, 8-, 9-, ?10q, 13-, 15-, 16-, 17-, 20-                                          | Hypodiploid  | Hyperdiploid | Incorrect    |
| AYAII\_0369         | AYAII\_0369           | female | hg19   | Hyperdiploid | High Hyperdiploid  | “82-86 (near tetraploid) gain of 1q in 7/20”                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |       82 |       45 | ?1p, ?1q+, ?4p-, 4q-, ?6, 9p-, ?9q-, ?10p, ?13, 14-, 18p-, ?18q-, ?20p                         | Hypodiploid  | Hyperdiploid | Incorrect    |
| AYAII\_0132         | AYAII\_0132           | male   | hg19   | None         | DUX4               | “49,XY,+X,+21,+21\[11\]/47,XY,+12\[8\]/46,XY\[1\].ish 9q34(ABLx2),22q11(BCRx2)”                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |       49 |       46 | ?21, ?Xp+                                                                                      | Normal       | Hyperdiploid | Incorrect    |
| ADII\_0044          | ADII\_0044            | male   | hg19   | Hypodiploid  | Hypodiploid        | “46,XY by G-banding. Hypodiploid loss of 3,13,15,16,17,19,20 by SNParray”                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |       46 |       38 | 3-, 7-, 9p-, 11p-, ?12p, 12q-, 13-, 15-, 16-, 17-, 19-, 20-, ?21                               | Hypodiploid  | Normal       | Incorrect    |
| ADII\_0097          | ADII\_0097            | male   | hg19   | Hyperdiploid | Hyperdiploid       | “68,XXYY,+1,+2,+4,+4,+5,+6,+6,+8,+10,+10, +11,+11,+12,+14,+15,+19,+21,+21,+22,+22\[6\] 46,XY\[9\]”                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |       68 |       45 | ?1, ?2+, ?3p-, ?3q, 7-, ?8, 9p-, ?9q-, ?12p+, ?12q, ?13-, ?16-, ?17-, ?19+, ?20-               | Hypodiploid  | Hyperdiploid | Incorrect    |
| QCTB\_CHII\_1105901 | CHII\_QCTB.1105901    | male   | hg19   | IGH@:BCL2    | BCL2/MYC           | “45,XX,der(13;14)(q10;q10)?c\[20\].nuc ish(5’MLL,3’MLL)x2(5’MLL con3’MLLx2)\[200\],(ETV6x2,RUNX1x3)(ETV6 con RUNX1x1)\[199/200\]”                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |       45 |       47 | 7p-, 21+                                                                                       | Hyperdiploid | Hypodiploid  | Incorrect    |
| CHI\_0453           | CHI\_0453.Dx.BM       | female | GRCh37 | ETV6:NTRK3   | Ph-like(Kinase)    | “47,XX,del(1)(q21),del(4)(q25q31),-9,+10,add(12)(p13),+mar\[5\]/46,XX\[15\].nuc ish(CRLF2x2)\[200\],(PBX1,HLF,TCF3)x2\[200\],(D4Z1x2,D10Z1x3,D17Z1x2)\[161/200\],(ABL1,BCR)x2\[200\],(KMT2Ax2)\[200\],(ETV6x2)(5’ETV6 sep 3’ETV6x1)\[186/200\],(ETV6x1,ETV6 dimx2,RUNX1x2)\[173/200\]”                                                                                                                                                                                                                                                                                                                                                                                              |       47 |       46 | ?10p+, 10q+, ?19p, ?Xp                                                                         | Normal       | Hyperdiploid | Incorrect    |
| ALL5-RAD07          | ALL5\_RAD07.PB        | male   | GRCh37 | BCR:ABL1     | Ph                 | “44,XY,der(5)(5pter-&gt;5q33::?::9q?13-&gt;9q34::11q?13-&gt;11q?21::?),-7,-9, t(9;22)(q34;q11.2),del(11)(q13q23)\[19\]/45,idem,+der(22)t(9;22)\[16\]/46,XY\[4\]”                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |       44 |       46 | 7p-, 9p-                                                                                       | Normal       | Hypodiploid  | Incorrect    |
| CHI\_0511           | CHI\_0511.Dx.PB       | female | GRCh37 | BCR:ABL1     | Ph                 | “45,XX,der(9)t(9;22)(p1?2;q11.2)t(9;22)(q34;q11.2),-22,der(22)t(9;22)(q34;q11.2)\[10\]/46,XX\[10\].nuc ish(CRLF2x2)\[200\],(PBX1,HLF,TCF3)x2\[200\],(D4Z1,D10Z1,D17Z1)x2\[200\], (ABL1x4,BCRx3)(ABL1 con BCRx3)\[159/200\]/(ABL1,BCR)x3(ABL1 con BCRx2)\[8/200\], (KMT2Ax2)\[200\],(ETV6,RUNX1)x2\[200\]”                                                                                                                                                                                                                                                                                                                                                                           |       45 |       46 | 9p-                                                                                            | Normal       | Hypodiploid  | Incorrect    |
| RNS050              | AYAII\_RNS050.Dx.PB   | female | GRCh37 | PAX5.p.P80R  | PAX5 P80R          | “45,XX(Der(7;9)q10q10 20/20”                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |       45 |       46 | 7p-, 9p-                                                                                       | Normal       | Hypodiploid  | Incorrect    |
| ADI\_0523           | ADII\_0523.Dx.BM      | male   | GRCh37 | BCR:ABL1     | Ph                 | “45,XY,add(1)(q42),-7,der(9)(9qter-&gt;9q22::9p13-&gt;9q22::?),t(9;22)(q34;q11.2)\[9\]/44,sl,der(14)t(14;19)(p11.2:p12),-19\[2\]/46,XY\[8\].nuc ish(ABL1,BCR)x3(ABL1 con BCRx2)\[174/200\]”                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |       45 |       46 | 7p-, 9p-, ?20p                                                                                 | Normal       | Hypodiploid  | Incorrect    |
| QCTB1374            | CHI\_QCTB1374.0930601 | male   | GRCh37 | ETV6:RUNX1   | ETV6-RUNX1         | “48,XY,+10,add(12)(p13),add(12)(p13),add(18)(q2?),+add(18),del(21)(q22) \[8\]/47,XY,idem,-13\[6\]/47,XY,+10,add(12),add(12),-13,add(18),+mar\[2\]/46 ,XY\[6\].ish add(18)(ETV6 con RUNX1+,3Ã”Ã¸Î©RUNX1+),del(21)(5Ã”Ã¸Î©RUNX1+). nuc ish(DXZ1x1,DYZ3x1,(D18Z1x3)\[100/200\],(D4Z1x2,D10Z1x3,D17Z1x2) \[126/200\],(SCFD2,CHIC2,PDGFRA)x2\[214\],(ABL1,BCR)x2\[208\],(KMT2Ax2) \[216\],(5Ã”Ã¸Î©ETV6x2,3Ã”Ã¸Î©ETV6x1)(5Ã”Ã¸Î©ETV6 sep 3Ã”Ã¸Î©ETV6x1)\[152/220\],(D13S319x1, LAMP1x2)\[36/212\]/(D13S319,LAMP1)x1\[32/212\],(ETV6x3,RUNX1x4)(ETV6 con RUNX1x2)\[141/200\],(3Ã”Ã¸Î©RUNX1x3,5Ã”Ã¸Î©RUNX1x2)(3Ã”Ã¸Î©RUNX1 sep 5Ã”Ã¸Î©RUNX1x1)(3Ã”Ã¸Î©RUNX1 con 5Ã”Ã¸Î©RUNX1x1)\[122/177\]" |       48 |       46 | ?4q+, ?10p, 10q+, ?12p, ?18                                                                    | Normal       | Hyperdiploid | Incorrect    |
| AYAI\_0599          | AYAI\_0599.Dx.PB      | female | GRCh37 | Hypodiploid  | Masked Hypodiploid | “48,XX,+21,+21\[4\]/46,XX\[16\]”                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |       48 |       25 | 1-, 2-, 3-, 4-, 5-, 6-, 7-, 8-, 9-, 10-, 11-, 12-, 13-, 14-, 15-, 16-, 17-, 18-, 19-, 20-, 22- | Hypodiploid  | Hyperdiploid | Incorrect    |

### Comparison of Reported Chromosome Number

Next, the ability of RNAseqCNV to accurately determine the exact number
of chromosomes was investigated. RNAseqCNV accurately reported 46
chromosomes in 92% of “Normal” samples; however, the exact number of
chromsomes in the hyperdiploid and hypodiploid groups were only 15% and
24% accurate, respectively. Overall, the accuracy at determining exact
numbers was 48%.

``` r
# Create column comparing RNAseqCNV chrom_n with cytogenetic data karyo_n
cytoData <- cytoData %>% 
  mutate(correctN = if_else(chrom_n == karyo_n, "Correct", "Incorrect"))
 
# Create dataframe for plotting percentages
nAltPerc <- cytoData %>% 
   count(oldCol, correctN) %>% 
   group_by(oldCol) %>% 
   mutate(label_y = n/sum(n)) 
 
# Plot percentages 
 nAltPlot <- nAltPerc %>% 
   ggplot(aes(x=oldCol, y=label_y, fill=correctN)) +
   geom_col(position="stack") +
   geom_text(aes(label=paste(paste0(round(label_y*100, digits=0), "%"), paste0("(", n, ")"))), 
             position=position_stack(vjust = 0.5),
             colour="white",
             size=3,
             fontface="bold") +
   scale_fill_manual(values= cbPalette) +
   theme_classic() +
   labs(x = "Ploidy",
        y = "Percentage of Samples",
        fill = "RNAseqCNV Result") +
   ggtitle("Ability of RNAseqCNV to correctly determine exact chromosome number")
```

![](CNV_Bench_files/figure-gfm/nAltPlot-1.png)<!-- -->

#### Incorrectly Reported Chromosome Number

Overall, RNAseqCNV appears to underestimate the number of chromosomal
alterations. This is due to the confidence metrics calculated by the
tool; low confidence calls are reported but will not be included in the
final chromosome number generate. Due to this, reports of alterations
with a final chromosome number between 45-47 should be interpreted with
caution.

In addition, RNAseqCNV calls of low hypodiploid samples should be
investigated for a diagnosis of masked hypodiploid.

``` r
# Ensure columns of interest are numeric
cytoData$karyo_n <- as.numeric(cytoData$karyo_n)
cytoData$chrom_n <- as.numeric(cytoData$chrom_n)

# Plot only incorrect N
filteredCyto <- cytoData %>% filter(correctN == "Incorrect")

# Reorder for plotting
filteredCyto <- filteredCyto[order(filteredCyto$oldCol),]

filteredCyto$patientID.x <- factor(filteredCyto$patientID.x,
                          levels = filteredCyto$patientID.x[order(filteredCyto$oldCol)])

# Labels for plotting
y_labs <- filteredCyto$oldCol
y_labs2 <- filteredCyto$primary.subtype


karyBars <- filteredCyto  %>%
  ggplot() +
  geom_col(aes(x=patientID.x, y=chrom_n, fill = correctBench), width= 0.7) +
  geom_col(aes(x=patientID.x, y= karyo_n, fill = "Number on Cytogenetics Report"), colour = "black", width = 0.7, alpha = 0.1) +
  labs(x = "Actual Classification",
       y = "Number of Chromosomes") +
  ggtitle("Samples with Incorrect Chromosome Number Reported by RNAseqCNV") +
  theme_classic() +
 # theme(axis.text.y = element_blank(), axis.line.y = element_blank()) +
  coord_flip(ylim = c(20,65)) +
  geom_hline(aes(yintercept = 46, linetype= paste("Normal (n = 46)"))) +
  scale_fill_manual(name = "RNAseqCNV Benchmarking", values = c(cbPalette[1:2], "white")) +
  scale_linetype_manual("Line", values = 1) +
  scale_x_discrete(labels = paste0(y_labs2, " [", y_labs, "]")) 
```

![](CNV_Bench_files/figure-gfm/chrom%20number%20plot-1.png)<!-- -->

#### Identifying Clinical Hyperdiploid or Hypodiploid ALL

Despite the inaccuracy of chromosome number calls, RNAseqCNV correctly
classified 88% of our samples where aneuploidy was the key alteration.

``` r
# Filter data for only clinically Hyper/Hypodiploid samples
aneuploid <- noCyto %>% 
  filter(col %in% c("Hypodiploid","Hyperdiploid"))

# Create dataframe for plotting
 aneuPerc <-  aneuploid %>% 
   count(primary.subtype, correctBench) %>% 
   group_by(primary.subtype) %>% 
   mutate(label_y = n/sum(n))
 
 aneuPlot <- aneuPerc %>%
   ggplot(aes(x=primary.subtype, y=label_y, fill = correctBench)) +
   geom_col(position="stack") +
   geom_text(aes(label=paste(paste0(round(label_y*100, digits=0), "%"), paste0("(", n, ")"))), 
             position=position_stack(vjust = 0.5),
             colour="white",
             size=3,
             fontface="bold") +
   scale_fill_manual(values= cbPalette[1:2]) +
   theme_classic() +
   theme(axis.text.x = element_text(angle = 45, vjust=1, hjust = 1)) +
   labs(x = "Ploidy",
        y = "Percentage of Samples",
        fill = "RNAseqCNV Result") +
   ggtitle("Ability of RNAseqCNV to diagnose Hyperdiploid and Hypodiploid ALL Subtypes")
```

![](CNV_Bench_files/figure-gfm/aneuPlot-1.png)<!-- -->

| correctBench |   n | Percentage |
|:-------------|----:|-----------:|
| Correct      |  49 |   83.05085 |
| Incorrect    |  10 |   16.94915 |

Of the 7 incorrect samples, two were Masked Hypodiploid, which are
difficult to detect; the tool classified all masked hypodiploids as
Hypodiploid, whereas in our reports the final classification differs
between Hyperdiploid and Hypodiploid.

One of the incorrect samples was mentioned earlier as having ambigous
cytogenetics associated with it. The patient was classified as having 46
chromosomes but the report also stated Hypodiploid losses; RNAseqCNV
correctly classified this patient as Hypodiploid.

If we look at the final classification reported by our research team,
the results match in 93% of cases. Of the remaining 4 incorrect cases,
only one had a cytogenetics report. This sample was reported to be
hyperdiploid whilst RNAseqCNV classified this as Hypodiploid; our notes
suggest that this is another masked hyposiploid cases based on
“clustering and TP53 mutations with VAF 0.98”. The remaining three all
had alterations reported by RNAseqCNV, but these alterations were
considered low confidence and so were not included in the final number
reported.

``` r
# Create column comparing RNAseqCNV chrom_n with cytogenetic data karyo_n
aneuploidClin <- aneuploid %>% 
  mutate(correctClin = if_else(cnvCol == col, "Correct", "Incorrect"))


aneuPercClin <-  aneuploidClin %>% 
  count(primary.subtype, correctClin) %>% 
  group_by(primary.subtype) %>% 
  mutate(label_y = n/sum(n))

aneuPlotClin <- aneuPercClin %>%
  ggplot(aes(x=primary.subtype, y=label_y, fill = correctClin)) +
  geom_col(position="stack") +
  geom_text(aes(label=paste(paste0(round(label_y*100, digits=0), "%"), paste0("(", n, ")"))), 
            position=position_stack(vjust = 0.5),
            colour="white",
            size=3,
            fontface="bold") +
  scale_fill_manual(values= cbPalette[1:2]) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust = 1)) +
  labs(x = "Ploidy",
       y = "Percentage of Samples",
       fill = "RNAseqCNV Result") +
  ggtitle("Ability of RNAseqCNV to diagnose Hyperdiploid and Hypodiploid ALL")
```

![](CNV_Bench_files/figure-gfm/clin%20plot-1.png)<!-- -->

## Conclusion

RNAseqCNV is able to correctly identify CNVs in patients where
aneuploidy is their primary alteration. This tool can be used for
research purposes to reclassify “B-Other” patients with missing
cytogenetics, and can be adapted to be used in our pipeline on new
samples.
