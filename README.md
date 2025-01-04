# Reticulocyte Transcriptomics (Gene Expression) Analysis Pipeline

## Overview

This repository provides a comprehensive bioinformatics pipeline for analyzing **Gene Expression Profiling** in **Reticulocytes**, with a focus on understanding the **transcriptional regulation of fetal hemoglobin (HbF)** in **Sickle Cell Disease (SCD)**. The goal of this pipeline is to identify genetic markers and regulatory factors that influence fetal hemoglobin levels, which may lead to enhanced treatment strategies for SCD.

## Project Details

The project includes various analysis steps, including preprocessing of RNA-Seq data, differential gene expression analysis, and the modeling of gene expression patterns. Key elements include the exploration of fetal hemoglobin and its regulatory pathways in reticulocytes.

### Key Components of the Pipeline:
- **Gene Expression Profiling** of reticulocytes.
- **Transcriptomic Analysis** focused on fetal hemoglobin regulation.
- **Development of predictive models** for treatment strategies in Sickle Cell Disease.

## Files in the Repository

1. **`Gexp_SCD.ipynb`**  
   A Jupyter notebook for gene expression analysis using RNA-Seq data. This notebook includes preprocessing, visualization, and statistical analysis of gene expression data related to SCD.

2. **`Model_GeneExpression.1.R`**  
   R script for performing gene expression modeling. It implements statistical models to identify significant changes in gene expression associated with fetal hemoglobin regulation.

3. **`Phenotype.txt`**  
   Contains phenotype data related to Sickle Cell Disease patients. The file is used to correlate gene expression data with clinical outcomes.

## How to Run the Pipeline

### Prerequisites
- **Python** (for Jupyter notebook): Ensure that you have Python installed along with required libraries (e.g., pandas, matplotlib, seaborn).
- **R**: Ensure that R is installed and that necessary libraries are available (e.g., DESeq2, ggplot2).

### Steps:
1. **Clone the Repository**:
