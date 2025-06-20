# RNA-Seq Analysis Script 

## 1. Install Required packages

### This installs the BiocManager package
install.packages("BiocManager")

### This installs the core packages of Bioconductor
BiocManager::install()

### DESeq2 package
BiocManager::install("DESeq2")

### A human genome-wide annotation package.
BiocManager::install("org.Hs.eg.db")

### A package to create heatmaps
install.packages("pheatmap")

### A package for data visualization
install.packages("ggplot2")

### An extension of ggplot2 that prevents overlapping labels in plots
install.packages("ggrepel")

### load an installed package so that you can use its functions, datasets, and tools
library(BiocManager)
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pheatmap)
library(ggplot2)
library(ggrepel)

### Get Working Directory
getwd()
# In console result will form like 'E:/Project/NGS/RStudio/RNASeq'

### Set Working Directory
setwd("E:/Project/NGS/RStudio/RNASeq")
# Now R will look for data inside your RNASeq folder.

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

## 2. Load data 

# Before load data keep read_counts.txt file in working folder where RNA_SEQ.R file is

### Read all lines
lines <- readLines("read_counts.txt")

### Check number of columns in each line
line_lengths <- sapply(strsplit(lines, "\t"), length)

### Keep only lines with 14 columns (i.e., correct ones)
clean_lines <- lines[line_lengths == 14]

### Write cleaned data to new file
writeLines(clean_lines, "clean_read_counts.txt")

# Read the cleaned file
counts_txt <- read.delim("clean_read_counts.txt", header = TRUE, row.names = 1, sep = "\t")

### Remove unwanted columns by name (if they exist)
unwanted_cols <- c("Chr", "Start", "End", "Strand", "Length")
counts_txt <- counts_txt[ , !(colnames(counts_txt) %in% unwanted_cols) ]

### Clean column names
colnames(counts_txt) <- gsub(".*(Sample_\\d+_\\d+).*", "\\1", colnames(counts_txt))

### Save as CSV
write.csv(counts_txt, "read_counts.csv")

### Load from CSV
Counts <- read.delim("read_counts.csv", header = TRUE, row.names = 1, sep = ",")

### Keeping only genes with row sums greater than 50
Counts <- Counts[which(rowSums(Counts) > 50), ]

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

## 3. Run DESeq2

colnames(Counts)

### Define the experimental condition
condition <- factor(c("disease", "disease", "disease", "disease",
                      "control", "control", "control", "control"))

### Create the coldata data frame
coldata <- data.frame(
  row.names = colnames(Counts),
  condition = condition
  )

### Create the DESeq2 Dataset
dds <- DESeqDataSetFromMatrix(
  countData = Counts,
  colData = coldata,
  design = ~ condition  
)

### Run DESeq
dds <- DESeq(dds)

### Variance Stabilizing Transformation (VST)
vsdata <- vst(dds, blind = FALSE)

### PCA Plot
plotPCA(vsdata, intgroup = "condition")

### Dispersion Plot
plotDispEsts(dds)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

## 4. Pairwise comparisons between samples

### Perform pairwise comparison
res_disease_vs_control <- results(dds, contrast = c("condition", "disease", "control"))

### Remove NAs
sigs_disease_vs_control <- na.omit(res_disease_vs_control)

### Filter for significant DEGs (FDR < 0.05)
sigs_disease_vs_control <- sigs_disease_vs_control[sigs_disease_vs_control$padj < 0.05, ]

### Add gene symbols (if not already done)
library(org.Hs.eg.db)
library(AnnotationDbi)
sigs_disease_vs_control.df <- as.data.frame(sigs_disease_vs_control)
sigs_disease_vs_control.df$gene_name <- mapIds(org.Hs.eg.db,
                                               keys = rownames(sigs_disease_vs_control.df),
                                               column = "SYMBOL",
                                               keytype = "ENSEMBL")

### Classify genes for volcano plot
sigs_disease_vs_control.df$diffexpressed <- "NO"
sigs_disease_vs_control.df$diffexpressed[sigs_disease_vs_control.df$log2FoldChange > 0.6 & sigs_disease_vs_control.df$pvalue < 0.05] <- "UP"
sigs_disease_vs_control.df$diffexpressed[sigs_disease_vs_control.df$log2FoldChange < -0.6 & sigs_disease_vs_control.df$pvalue < 0.05] <- "DOWN"

### Save output
write.csv(sigs_disease_vs_control.df, "DEGs_disease_vs_control.csv")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

## 5. Generating an ordered list of DEGs

### Sort by absolute log2 fold change
sigs_disease_vs_control_sorted <- sigs_disease_vs_control.df[order(abs(sigs_disease_vs_control.df$log2FoldChange), decreasing = TRUE), ]

### Export only GeneID and log2FC
log2FC_file <- "DEGs_disease_vs_control_log2FC.txt"

write.table(
  data.frame(
    GeneID = rownames(sigs_disease_vs_control_sorted),
    log2FoldChange = sigs_disease_vs_control_sorted$log2FoldChange
  ),
  file = log2FC_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

### Export GeneID + log2FC + p-value + padj (FDR)
pvalue_file <- "DEGs_disease_vs_control_with_pvalues.txt"

write.table(
  data.frame(
    GeneID = rownames(sigs_disease_vs_control_sorted),
    log2FoldChange = sigs_disease_vs_control_sorted$log2FoldChange,
    pvalue = sigs_disease_vs_control_sorted$pvalue,
    padj = sigs_disease_vs_control_sorted$padj
  ),
  file = pvalue_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

## 6. Generate heatmaps

### Extract gene IDs of DEGs
top_genes <- rownames(sigs_disease_vs_control.df)

### Filter the Counts matrix for DEGs only
top_counts <- Counts[top_genes, ]

### Plot the heatmap
pheatmap(
  top_counts,
  scale = "row",                   # normalize expression per gene
  show_rownames = FALSE,           # hide gene names if too many
  cluster_cols = FALSE,            # keep your sample order
  legend = TRUE,
  main = "Heatmap of Significant DEGs"
)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

## 7A. Generating volcano plots

### Prepare the data frame for volcano plot
volcano_data <- sigs_disease_vs_control.df

### Confirm columns exist: log2FoldChange, pvalue, padj, diffexpressed
head(volcano_data)

### Plot the volcano plot
volcano_plot <- ggplot(data = sigs_disease_vs_control.df, 
            aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed)) +
  geom_point() +
  scale_color_manual(values = c("NO" = "black", "UP" = "red", "DOWN" = "blue")) +
  geom_text(aes(label = delabel), vjust = -0.5, size = 3) +
  theme_minimal() +
  geom_vline(xintercept = c(-0.6, 0.6), col = "red") +
  geom_hline(yintercept = -log10(0.05), col = "red") +
  theme(panel.grid = element_blank())

# Print the plot
print(volcano_plot)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

## 7B. Generating volcano plots with labeled DEGs

### Add gene labels only for significantly expressed genes
sigs_disease_vs_control.df$delabel <- ifelse(
  !is.na(sigs_disease_vs_control.df$gene_name) & sigs_disease_vs_control.df$diffexpressed != "NO",
  sigs_disease_vs_control.df$gene_name,
  NA
)

### Generate labeled volcano plot
pl <- ggplot(data = sigs_disease_vs_control.df, 
            aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed)) +
  geom_point() +
  scale_color_manual(values = c("NO" = "black", "UP" = "red", "DOWN" = "blue")) +
  geom_text(aes(label = delabel), vjust = -0.5, size = 3) +
  theme_minimal() +
  geom_vline(xintercept = c(-0.6, 0.6), col = "red") +
  geom_hline(yintercept = -log10(0.05), col = "red") +
  theme(panel.grid = element_blank())

# Print the plot
print(pl)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# END
