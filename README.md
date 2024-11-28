# BD-Rhapsody
BD Rhapsody™ Single-Cell Analysis


# **Pipeline Overview**

This pipeline processes single-cell RNA-seq data, integrating multiple datasets from Rhapsody while ensuring high-quality analysis. The main steps are:

1. **Project Initialization**  
   - Set up directories for raw data, results, and scripts.  
   - Load required libraries and utility functions.

2. **Data Loading and Preprocessing**  
   - Import Seurat objects (`pci_index1`, `pci_index3`).  
   - Filter cells based on QC metrics (e.g., mitochondrial contamination, gene count).  
   - Visualize quality metrics using violin and scatter plots.

3. **Normalization and Variable Gene Selection**  
   - Normalize data using `LogNormalize`.  
   - Identify highly variable features for downstream analysis.

4. **Dimensionality Reduction**  
   - Perform PCA for initial dimensionality reduction.  

5. **Integration Methods**  
   - Apply various integration techniques:
     - **CCA (Canonical Correlation Analysis)**
     - **RPCA (Reciprocal PCA)**
     - **Harmony**
     - **Seurat v5 Multi-layer Integration**
   - Evaluate integration quality using metrics such as iLISI and cLISI.

6. **Visualization and Results**  
   - Generate UMAP/t-SNE plots for each method.  
   - Save integrated datasets (`pci_integrated_harmony.rds`, `pci1_3_combined_s5.rds`).  
   - Compare integration performance across methods using iLISI and cLISI scores.

Feel free to adapt this pipeline based on specific dataset requirements or research questions.
