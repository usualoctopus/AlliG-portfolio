---
layout: post
title: "Differential Expression Analysis and Visualization in Neuronal Stem Cells"
date: 2024-06-10
categories: r-scripts
---


# Differential Expression Analysis and Visualization in Neuronal Stem Cells

## Introduction
This script performs a comprehensive differential expression analysis of RNA-seq data obtained from induced pluripotent stem cells (iPSCs) differentiated into neuronal stem cells (NSCs) for PMDD patients and matched controls. The analysis includes normalization, filtering, and statistical testing to identify differentially expressed genes (DEGs) under various hormone treatments (E2, P4, ALLO). The script also provides data visualization techniques including density plots, MA plots, volcano plots, and bar graphs to facilitate interpretation of the results.

## Packages Required
```r
library(edgeR)
library(ggplot2)
library(RColorBrewer)
# source("http://www.bioconductor.org/biocLite.R")

## Loading Data


```R
setwd("~/Documents/Documents_AlliMacBookPro/PMDD_NSC_edgeR")
raw.all <- read.csv("NSC_counts.csv")
counts.all <- raw.all[ , -c(1) ]
rownames(counts.all) <- raw.all[ , 1 ]
head(counts.all)

# Ensure counts.all is numeric
counts.all <- as.matrix(counts.all)
mode(counts.all) <- "numeric"
```

## Preprocess Metadata
This code ensures that the metadata and sample information are properly organized and prepared, which is critical for performing accurate and meaningful differential expression analysis. Create a txt file with relevant metadata; e.g., sample names (that match the data table), diagnosis (pmdd or control), treatment (E2, P4, or Allo). This file is used below to create factors for Dx (diagnosis) and Patient (patient IDs for each treatment), allowing the model to correctly partition the variance attributed to different sources (e.g., treatment effects, patient variability). These subsets and factors are created dynamically so that this script can easily handle multiple treatments and conditions, making it adaptable to various experimental designs. 


```R
# Define dx/treatment groups
treatments <- c("Veh", "E2", "P4", "Allo")
controls <- c("7951", "9890", "1506", "2551")
pmdd <- c("8319", "0319", "1134", "0654")

# Create a list of sample names for each dx & tx
sample_names <- list()
for (treatment in treatments) {
  sample_names[[paste0("control_", treatment)]] <- paste0(treatment, "_", controls, "_avg")
  sample_names[[paste0("pmdd_", treatment)]] <- paste0(treatment, "_", pmdd, "_avg")
}

# Print sample names to verify
# sample_names

# Read metadata
meta_all <- read.table("NSC_metadata.txt", sep="\t", header = TRUE)
rownames(meta_all) <- meta_all$Sample_ID
meta_all$Sample_ID <- NULL

# Create metadata subsets dynamically
meta_subsets <- list()
for (treatment in treatments) {
  control_samples <- sample_names[[paste0("control_", treatment)]]
  pmdd_samples <- sample_names[[paste0("pmdd_", treatment)]]
  
  # Subset metadata for each treatment
  meta_subsets[[treatment]] <- meta_all[combined_samples, ]
}


# Combine all samples' metadata into one subset
all_samples <- unlist(sample_names)
meta_all_combined <- meta_all[all_samples, ]
meta_subsets[["all"]] <- meta_all_combined


# Set targets for each treatment
targets <- list()
for (key in names(meta_subsets)) {
  targets[[key]] <- meta_subsets[[key]]
}

# Define factors for differential expression analysis
Dx_levels <- c("Control", "PMDD")
Patient_levels <- 1:4

# Create Dx factors and Patient factors dynamically
Dx_factors <- list()
Patient_factors <- list()

for (key in names(targets)) {
  # create dx factors for each target
  Dx_factors[[key]] <- factor(targets[[key]]$Dx, levels = Dx_levels)
  
  # create patient factors for each target
  if (key == "all") {
    Patient_factors[[key]] <- as.factor(rep(Patient_levels, length(treatments) * 2))
  } else {
    Patient_factors[[key]] <- as.factor(rep(Patient_levels, 2))
  }
}


# Print Dx and Patient factors to verify
# Dx_factors
# Patient_factors
```

## Differential Expression Analysis Function
The purpose of writing it as a function here is to ensure that one can easily switch between treatments by simply changing the input parameter in the function call.


```R
# Function to perform differential expression analysis for a given treatment
perform_DE_analysis <- function(treatment, meta_data, counts_data, sample_names) {
  
  # Define control and PMDD sample names based on the treatment
  control_samples <- sample_names[[paste0("control_", treatment)]]
  pmdd_samples <- sample_names[[paste0("pmdd_", treatment)]]
  combined_samples <- c(control_samples, pmdd_samples)
  
  # Append "_avg" suffix to sample names to match the column names in counts_data
  combined_samples <- paste0(combined_samples, "_avg")
  
  # Subset metadata and counts for the given treatment
  meta_treatment <- meta_data[combined_samples, , drop = FALSE]
  counts_treatment <- counts_data[, combined_samples, drop = FALSE]
  
  # Create DGEList object for the treatment comparison
  cds <- DGEList(counts = counts_treatment, group = meta_treatment$Dx)
  
  # Filter genes: a gene is only retained if it is expressed at a minimum level
  keep <- filterByExpr(cds, min.count = 6)
  summary(keep)
  cds <- cds[keep, , keep.lib.sizes = FALSE]
  
  # Recompute library sizes after filtering
  cds <- calcNormFactors(cds)
  
  # Create the design matrix for the treatment comparison
  design <- model.matrix(~ 0 + Dx, data = meta_treatment)
  colnames(design) <- sub("Dx", "", colnames(design))
  
  # Estimate dispersion: design = design matrix, tagwise = estimate separate dispersion for each gene,
  # robust ensures dispersion estimates are less affected by outliers
  cds <- estimateDisp(cds, design = design, tagwise = TRUE, robust = TRUE)
  
  # Fit the model using quasi-likelihood for the treatment comparison
  fit <- glmQLFit(cds, design)
  
  # Define the contrast for PMDD vs. Control for the treatment
  contrast <- makeContrasts(PMDD_vs_Control = PMDD - Control, levels = design)
  
  # Perform the quasi-likelihood F-test for the treatment comparison
  qlf <- glmQLFTest(fit, contrast = contrast)
  
  # Get top differentially expressed genes for the treatment comparison
  top_genes <- topTags(qlf)
  print(top_genes)
  
  # Create results table for the treatment comparison
  resultsTbl <- topTags(qlf, n = nrow(qlf$table))$table
  
  # Save results to a file for the treatment comparison
  output_filename <- paste0("DEG_PMDD_", treatment, "_vs_Control_", treatment, ".csv")
  write.csv(resultsTbl, output_filename)
  
  return(resultsTbl)
}


sample_names <- list(
  control_E2 = sub('_avg$', '', paste0("E2_", c("7951", "9890", "1506", "2551"), "_avg")),
  pmdd_E2 = sub('_avg$', '', paste0("E2_", c("8319", "0319", "1134", "0654"), "_avg")),
  control_P4 = sub('_avg$', '', paste0("P4_", c("7951", "9890", "1506", "2551"), "_avg")),
  pmdd_P4 = sub('_avg$', '', paste0("P4_", c("8319", "0319", "1134", "0654"), "_avg")),
  control_Allo = sub('_avg$', '', paste0("Allo_", c("7951", "9890", "1506", "2551"), "_avg")),
  pmdd_Allo = sub('_avg$', '', paste0("Allo_", c("8319", "0319", "1134", "0654"), "_avg"))
)
```

## Run Differential Expression Analysis


```R
results_E2 <- perform_DE_analysis("E2", meta_all, counts.all, sample_names)
results_P4 <- perform_DE_analysis("P4", meta_all, counts.all, sample_names)
results_ALLO <- perform_DE_analysis("Allo", meta_all, counts.all, sample_names)
```

## Data Visualization Functions
By creating a function, one can generate consistent visualizations for each treatment by simply changing the input results table, providing a flexible and reusable approach. The output of this function is a density plot of logFC, MA plot, and a volcano plot.


```R
# define a function to generate a density plot of logFC, MA plot, and volcano plot
generate_plots <- function(resultsTbl, title_suffix) {
  
  # Extract necessary information
  DE_genes <- rownames(resultsTbl)
  logFC <- resultsTbl$logFC
  PValue <- resultsTbl$PValue
  logCPM <- resultsTbl$logCPM
  FDR <- resultsTbl$FDR
  
  # Create data frame for ggplot
  plot_data <- data.frame(
    Gene = DE_genes,
    logFC = logFC,
    logCPM = logCPM,
    PValue = PValue,
    FDR = FDR
  )
  
  # Density plot of log fold changes
  p1 <- ggplot(plot_data, aes(x = logCPM)) +
    geom_line(stat = "density", color = "deeppink1") +
    geom_density(fill = "deeppink1", color = NA, alpha = 0.2) +
    geom_line(data = subset(plot_data, PValue < 0.05), aes(x = logCPM), stat = "density", color = "cyan1") +
    geom_density(data = subset(plot_data, PValue < 0.05), aes(x = logCPM), fill = "cyan1", color = NA, alpha = 0.2) +
    xlab('logCPM') + ylab('density') +
    ggtitle(paste("Density Plot", title_suffix)) +
    theme_minimal()
  
  # MA plot
  p2 <- ggplot(plot_data, aes(x = logCPM, y = logFC)) + 
    geom_point(size = 0.2, color = "black") + 
    geom_point(data = subset(plot_data, PValue < 0.05), aes(x = logCPM, y = logFC), size = 0.3, color = "cyan1") + 
    geom_point(data = subset(plot_data, FDR < 0.05), aes(x = logCPM, y = logFC), size = 0.3, color = "deeppink1") + 
    xlab("logCPM") + ylab("logFC") +
    ggtitle(paste("MA Plot", title_suffix)) +
    theme_minimal()
  
  # Volcano plot
  p3 <- ggplot(plot_data, aes(x = logFC, y = -log10(PValue))) +
    geom_point(size = 0.5, color = "gray") +
    geom_point(data = subset(plot_data, PValue < 0.05), aes(x = logFC, y = -log10(PValue)), size = 0.5, color = "cyan1") +
    geom_point(data = subset(plot_data, FDR < 0.05), aes(x = logFC, y = -log10(PValue)), size = 0.5, color = "deeppink1") +
    xlab("logFC") + ylab("-log10(PValue)") +
    ggtitle(paste("Volcano Plot", title_suffix)) +
    theme_minimal()
  
  # Print the plots
  print(p1)
  print(p2)
  print(p3)
  
}


# Generate plots for E2
generate_plots(results_E2, "after E2 Treatment")

# Generate plots for P4
generate_plots(results_P4, "after P4 Treatment")

# Generate plots for ALLO
generate_plots(results_ALLO, "after ALLO Treatment")

```

## Bar Chart of an Individual Gene Across Samples
This bar chart displays the expression levels (CPM) of a gene across individual samples, highlighting disparities in gene expression.


```R
generate_bargraph <- function(cpm_values, meta_data, gene, treatment) {
  
  # Filter metadata for the specified treatment
  meta_treatment <- meta_data[grep(treatment, rownames(meta_data)), ]
  
  # Extract CPM values for the specified gene
  DE_cpm <- cpm_values[gene, ]
  
  # Ensure CPM values are numeric
  DE_cpm <- as.numeric(DE_cpm)
  
  # Create a data frame with metadata and CPM values
  sample_names <- colnames(cpm_values)
  DE_patient <- factor(rep(c('C-7951', 'C-9890','C-1506', 'C-2551', 'P-8319', 'P-0319', 'P-1134', 'P-0654')))
  DE_dx <- factor(meta_treatment$Dx, levels = c("Control", "PMDD"))
  DE_tx <- factor(rep(treatment, length(DE_patient)))  
  
  DDE_all <- data.frame(
    Sample = sample_names,
    Patient = DE_patient,
    Dx = DE_dx,
    CPM = DE_cpm,
    Treatment = DE_tx
  )
  
  # Bargraph of an individual gene across each individual sample
  ggplot(DDE_all, aes(x = Patient, y = CPM, fill = Dx)) + 
    geom_bar(stat = "identity", position = "dodge", color = "black") + 
    scale_fill_brewer(palette = "Pastel1") + 
    ylab(paste0(gene, " expression (CPM)")) + 
    xlab("") + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste(gene, "expression for", treatment, "Treatment"))
}


# Example calls for different treatments
generate_bargraph(cpm_values, meta_all, "CREBRF", "E2")
generate_bargraph(cpm_values, meta_all, "TOX", "P4")
generate_bargraph(cpm_values, meta_all, "SNORA14A", "Allo")

```
