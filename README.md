# Yeast-Biofilm-RNAseq-Differential-Expression
An end-to-end bioinformatics pipeline for aligning RNA-seq reads and performing differential gene expression analysis to identify key transcriptional changes across three developmental stages of yeast biofilm formation.

## Project Overview

This repository contains the code and final report for an RNA-seq analysis project investigating the transcriptome of yeast (*S. cerevisiae*) at three key biofilm developmental stages: **Early**, **Thin**, and **Mature**. The workflow encompasses data alignment, quantification, statistical analysis, and visualization to uncover genes significantly associated with biofilm progression.

## Key Features

*   **Alignment & Quantification:** Bash scripts for processing raw FASTQ files using the STAR aligner on an HPC cluster (Compute Canada).
*   **Differential Expression:** Comprehensive **R** analysis using the `edgeR` package for rigorous statistical testing.
*   **Multi-Group Comparisons:** Analysis includes both pairwise comparisons between stages and a focused contrast of the Early stage against combined later (Thin + Mature) stages.
*   **Automated Visualization:** Code to generate publication-ready figures, including:
    *   Enhanced Volcano Plots
    *   Hierarchical Clustering Heatmaps (of top DEGs)
    *   MA Plots for each comparison
*   **Reproducible Workflow:** Complete documentation and code to replicate the entire analysis from raw counts to final results.


## Analysis Highlights

The analysis identified:
*   **142** significant DEGs between Early and Thin biofilm stages.
*   **187** significant DEGs between Early and Mature stages.
*   **89** significant DEGs between Thin and Mature stages.
*   **231** significant DEGs when comparing the Early stage to combined Late stages (Thin+Mature), with **142 genes upregulated in early biofilm** (e.g., stress-response genes) and **89 upregulated in late stages** (e.g., cell adhesion factors).

## Technologies Used

*   **Alignment:** `STAR`
*   **Differential Expression:** `R` with `edgeR`, `ggplot2`, `pheatmap`
*   **Gene Annotation:** `org.Sc.sgd.db`

## Usage

1.  **Run Alignment (Bash/HPC):** Execute the provided STAR alignment script on an HPC cluster to process FASTQ files and generate count data.
2.  **Differential Expression (R):** Run the main R script `yeast_biofilm_RNAseq_analysis.R` to perform the analysis and generate all figures and result tables. Ensure all required R packages are installed first.

## Final Report

The detailed methodology, results, interpretation, and discussion are available in the file **`Project_Report_Yeast_Biofilm_RNAseq.pdf`**.
