# A Guide to Basic RNA-Seq Analysis

This project is created to help beginners get a basic understanding of RNA-Seq and how to perform data analysis using R and Linux.

## Project Goal

* Get a basic idea of RNA-Seq
* Learn how to perform data analysis step by step

## Step-by-step Overview

1. First, read the file `docs/PPT.pdf`
2. Then go through `docs/Linux_RStudio_Guide.pdf` and take the first steps to start the project
3. For analysis in RStudio, download the count file from the link below:

   [read\_counts.txt](https://drive.google.com/uc?export=download&id=1drEZ62jLeY0Xx45inkT8QoC2RE5l54VZ)

## What Outputs to Expect

* CSV files with differential expression results
* Plots such as PCA, heatmap, volcano plot, and dispersion plot

## Repository Structure

```
A_Guide_to_Basic_RNASeq/
├── data/
│   └── read_counts.csv
├── docs/
│   ├── Linux_RStudio_Guide.pdf
│   └── PPT.pdf
├── results/
│   ├── DEGs_disease_vs_control.csv
│   ├── DESeq_Dispersion_Plot.jpeg
│   ├── DESeq_PCA_Plot.jpeg
│   ├── Heatmap.jpef
│   ├── golgi_vesicle_genes.csv
│   ├── signs_disease_vs_control.csv
│   ├── signs_disease_vs_control_gnames.csv
│   ├── volacno_plot.jpeg
│   └── volacno_plot_labeled_DEGs.jpeg
├── scripts/
│   ├── RNA_SEQ.R
│   ├── align_all.sh
│   ├── convert_sort_index_all.sh
│   └── trim_all.sh
├── CODE.txt
└── README.md
```

This guide is meant for learning purposes and to help understand how RNA-Seq data analysis works from start to finish.
