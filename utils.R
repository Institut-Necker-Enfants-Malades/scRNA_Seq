# utils.R


#' Installer et charger une liste de packages CRAN
#' 
#' Cette fonction vérifie si chaque package de la liste est installé. Si un package n'est pas installé, il l'installe puis le charge. 
#' Si le package est déjà installé, il est simplement chargé.
#'
#' @param pkgs Un vecteur de noms de packages CRAN à installer et charger.
#' @return Aucun retour mais imprime dans la console l'état de chaque package (installé ou déjà installé) avec sa version.
#' @examples
#' install_and_load(c("dplyr", "ggplot2"))
install_and_load <- function(pkgs) {
  for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = T)) {
      cat(sprintf("Installing: %s\n", pkg))
      install.packages(pkg, dependencies = TRUE)
      suppressPackageStartupMessages(library(pkg, character.only = TRUE, quietly = T))
      cat(sprintf("Installed and loaded: %s (version %s)\n", pkg, packageVersion(pkg)))
    } else {
      suppressPackageStartupMessages(library(pkg, character.only = TRUE, quietly = T))
      cat(sprintf("Package already installed and loaded: %s (version %s)\n", pkg, packageVersion(pkg)))
    }
  }
  
  cat("All specified packages have been installed and loaded.\n\n") 
}

#' Installer et charger une liste de packages CRAN
#' 
#' Cette fonction vérifie si chaque package de développement depuis GitHub de la liste est installé. Si un package n'est pas installé, il l'installe puis le charge. 
#' Si le package est déjà installé, il est simplement chargé.
#'
#' @param pkg Un vecteur des chenins de chaque packages a installe depuis GitHub.
#' @param pkg_names Un vecteur de noms de packages CRAN à installer et charger.
#' @return Aucun retour mais imprime dans la console l'état de chaque package (installé ou déjà installé) avec sa version.
#' @examples
#' install_and_load(c("satijalab/azimuth", "satijalab/seurat-data"))
install_and_load_dev <- function(pkg, pkg_names) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  library(devtools)
  
  for (pkg_name in pkg_names) {
  if (!require(pkg_name, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("Installing and loading: %s\n", pkg_name))
    devtools::install_github(pkg, dependencies = TRUE)
    suppressPackageStartupMessages(library(pkg_name, character.only = TRUE))
    cat(sprintf("Installed and loaded: %s (version %s)\n", pkg_name, packageVersion(pkg_name)))
    
  } else {
    suppressPackageStartupMessages(library(pkg_name, character.only = TRUE))
    cat(sprintf("Package already installed and loaded: %s (version %s)\n", pkg_name, packageVersion(pkg_name)))
        }
  }
  
  cat("All specified dev packages have been installed and loaded.\n\n")
}

# Violin plots: number of genes detected, number of UMI, % of mitochondrial contamination 
visViolin <- function(obj, thresould = 50) {
  p1 <- VlnPlot(obj, features = "nFeature_RNA") + NoLegend() + labs(title = "", x = "", y = "nGenes")
  p2 <- VlnPlot(obj, features = "nCount_RNA") + NoLegend() + labs(title = "", x = "", y = "nUMI")
  p3 <- VlnPlot(obj, features = "percent.mt") + NoLegend() + labs(title = "", x = "", y = "%MT") +
    geom_hline(yintercept = thresould, col = "red")
  cowplot::plot_grid(p1, p2, p3, ncol = 3)
}

#Additionally, we can also see which genes contribute the most to such reads. We can for instance plot the percentage of counts per gene.
visGene_expr <- function(obj, numGenes = 30){
  # Compute the proportion of counts of each gene per cell
  # Use sparse matrix operations, if your dataset is large, doing matrix devisions the regular way will take a very long time.
  C <- obj@assays$RNA@counts
  C@x <- C@x / rep.int(colSums(C), diff(C@p)) * 100
  most_expressed <- order(Matrix::rowSums(C), decreasing = T)[numGenes:1]
  p1 <- boxplot(as.matrix(t(C[most_expressed, ])),
                cex = 0.1, las = 1, xlab = "Percent counts per cell",
                col = (scales::hue_pal())(numGenes)[numGenes:1], horizontal = TRUE
  )
  return(list(
    most_expressed = as.data.frame(t(C[order(Matrix::rowSums(C), decreasing = T)[1:numGenes], ])),
    plot = p1
  ))
}

# Fonction pour déterminer les seuils de filtrage optimaux
FilterSeurat <- function(seurat_obj, 
                         nfeature_min_percentile = 1,
                         nfeature_max_percentile = 99,
                         ncount_min_percentile = 1, 
                         ncount_max_percentile = 99,
                         mt_max_percentile = 95) {
  # nFeature_RNA : Nombre de gènes détectés par cellule
  # Filtre les cellules entre le 1er et le 99e percentile par défaut
  # Élimine les cellules mortes/vides (peu de gènes) et les doublets (trop de gènes)
  # 
  # nCount_RNA : Nombre total de transcrits par cellule
  # Filtre les cellules entre le 1er et le 99e percentile par défaut
  # Élimine les cellules de mauvaise qualité et les outliers
  # 
  # percent.mt : Pourcentage de transcrits mitochondriaux
  # Filtre les cellules au-dessus du 95e percentile par défaut
  # Élimine les cellules mourantes/stressées
  
  # Extraction des métriques
  metrics <- data.frame(
    nFeature = seurat_obj$nFeature_RNA,
    nCount = seurat_obj$nCount_RNA,
    percent_mt = seurat_obj$percent.mt
  )
  
  # Calcul des seuils
  nfeature_bounds <- quantile(metrics$nFeature, 
                              probs = c(nfeature_min_percentile/100, 
                                        nfeature_max_percentile/100))
  ncount_bounds <- quantile(metrics$nCount,
                            probs = c(ncount_min_percentile/100,
                                      ncount_max_percentile/100))
  mt_threshold <- quantile(metrics$percent_mt,
                           probs = mt_max_percentile/100)
  
  # Application des filtres
  cells_keep <- rownames(metrics)[
    metrics$nFeature >= nfeature_bounds[1] &
      metrics$nFeature <= nfeature_bounds[2] &
      metrics$nCount >= ncount_bounds[1] &
      metrics$nCount <= ncount_bounds[2] &
      metrics$percent_mt <= mt_threshold
  ]
  
  # Création du nouvel objet filtré
  seurat_filtered <- subset(seurat_obj, cells = cells_keep)
  
  # Préparation du rapport de filtrage
  stats_before <- data.frame(
    Metric = c("Nombre de cellules",
               "nFeature_RNA (min)", "nFeature_RNA (max)",
               "nCount_RNA (min)", "nCount_RNA (max)",
               "percent.mt (max)"),
    Before = c(ncol(seurat_obj),
               min(metrics$nFeature), max(metrics$nFeature),
               min(metrics$nCount), max(metrics$nCount),
               max(metrics$percent_mt)),
    Threshold = c(NA,
                  nfeature_bounds[1], nfeature_bounds[2],
                  ncount_bounds[1], ncount_bounds[2],
                  mt_threshold),
    After = c(ncol(seurat_filtered),
              min(seurat_filtered$nFeature_RNA), max(seurat_filtered$nFeature_RNA),
              min(seurat_filtered$nCount_RNA), max(seurat_filtered$nCount_RNA),
              max(seurat_filtered$percent.mt))
  )
  
  # Calcul du pourcentage de cellules conservées
  percent_kept <- round(ncol(seurat_filtered) / ncol(seurat_obj) * 100, 1)
  
  # Création des visualisations
  p1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3, pt.size = 0.1) +
    theme(plot.title = element_text(hjust = 0.5))
  
  p2 <- VlnPlot(seurat_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3, pt.size = 0.1) +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(list(
    filtered_object = seurat_filtered,
    filtering_stats = stats_before,
    percent_kept = percent_kept,
    plots_before = p1,
    plots_after = p2
  ))
}

#' Detect doublets in a Seurat object
#'
#' @param seuratObj A Seurat object
#' @param est_doublet_model An optional linear model to estimate the doublet rate
#' @param pc A vector of principal components to use for doublet detection
#' @param reuse.results Logical, whether to reuse previously calculated doublet detection results
#' @return The input Seurat object with doublet information added to the metadata
detect_doublets <- function(seuratObj,
                            est_doublet_model = model_rhap,
                            pc = 1:20,
                            reuse.results = TRUE) {
  # Check if Seurat object is valid
  if (!is(seuratObj, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  
  DefaultAssay(seuratObj) <- "RNA"
  # Find pK values
  sweep.res.list <- paramSweep(seuratObj, 
                               PCs = pc, 
                               sct = F)
  sweep.stats <- summarizeSweep(sweep.res.list, 
                                GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% 
    as.character() %>% 
    as.numeric()
  
  # estimate doublet rate based on cell number
  DoubletRate = predict(est_doublet_model, 
                        data.frame(cell_num = dim(seuratObj)[2]))/100
  nExp_poi <- round(DoubletRate*ncol(seuratObj)) 
  
  seuratObj <- doubletFinder(seuratObj, 
                             PCs = pc, 
                             pN = 0.25, 
                             pK = pK_bcmvn, 
                             nExp = nExp_poi, 
                             reuse.pANN = F, 
                             sct = F)
  
  # Extract and store doublet classifications
  temp1 <- grepl("DF.classifications", 
                 colnames(seuratObj@meta.data), 
                 ignore.case = T)
  colnames(seuratObj@meta.data)[temp1] <- "doublet_check"
  
  # Extract and store doublet scores
  temp2 <- grepl("pANN", 
                 colnames(seuratObj@meta.data), 
                 ignore.case = TRUE)
  colnames(seuratObj@meta.data)[temp2] <- "doublet_score"
  
  seuratObj$doublet_check <- seuratObj$doublet_check
  
  # Return output
  return(seuratObj)
  
} 

#' Fonction pour visualiser les résultats d'intégration de 2 seurats object
#'
#' @param seuratObj A Seurat object
#' @param seurat1 A Seurat object
#' @param seurat2  A Seurat object
#' @param reduction est utilisé ici pour spécifier la méthode de réduction de dimension (e.g., CCA, RPCA)
#' @param title  A title
#' @return Plots of the 3 seurat
PlotIntegrationResults <- function(seurat1 = pci_index1_1_filt_done, seurat2 = pci_index3_1_filt_done, seurat_obj, reduction, title) {
  
  # Visualisation des données de PCI1 avant l'intégration
  p1 <- DimPlot(seurat1, reduction = "umap", group.by = "orig.ident") + 
    ggtitle(paste0("PCI1", " - by batch"))
  # Visualisation des données de PCI3 avant l'intégration
  p2 <- DimPlot(seurat2, reduction = "umap", group.by = "orig.ident") + 
    ggtitle(paste0("PCI3", " - by batch"))
  # Visualisation des données intégrées après application de la méthode d'intégration
  p3 <- DimPlot(seurat_obj, reduction = reduction, group.by = "orig.ident") + 
    ggtitle(paste0(title, " - by batch"))
  # Utilisation de cowplot::plot_grid pour afficher les trois graphiques 
  cowplot::plot_grid(p1, p2, p3, ncol = 3)
}

# 1. Calcul du score de mélange des lots
#devtools::install_github("immunogenomics/lisi")
#' Fonction pour calculer les scores d'integration lisi
#'
#' @param seuratObj A Seurat object
#' @param batch_key le parametre d'integration utiliser. Ici j'ai integre avec les index donc orig.ident. Cela peut-etre les echantillons 
#' @param reduction est utilisé ici pour spécifier la méthode de réduction de dimension a utiliser (e.g., UMAP, TSNE, CCA)
#' @return Les scores lisi et la moyenne
compute_integration_metrics <- function(seurat_obj, reduction ="umap",  batch_key = "orig.ident") {
  # Extraire les coordonnées UMAP et les métadonnées
  embeddings <- Embeddings(seurat_obj[[reduction]])
  meta_data <- seurat_obj@meta.data
  
  # Calculer les scores LISI
  lisi_scores <- compute_lisi(embeddings, meta_data, batch_key)
  # Calculer la moyenne des scores LISI
  mean_lisi <- mean(lisi_scores[[batch_key]], na.rm = TRUE)
  
  return(list(
    lisi_scores = lisi_scores,
    mean_lisi = mean_lisi
  ))
}

# Extraire les scores LISI pour chaque méthode et chaque lot
plot_integration_metrics_lot <- function(lisi_integration_scores, batch_key = "orig.ident"){
  lot_lisi_df <- data.frame(
    Method = rep(names(lisi_integration_scores), each = nrow(lisi_integration_scores[[1]]$lisi_scores)),
    Lot = rep(lisi_integration_scores[[1]]$lisi_scores$orig.ident, times = length(lisi_integration_scores)),
    LISI_score = unlist(lapply(lisi_integration_scores, function(x) x$lisi_scores[[batch_key]]))
  )
  
  # Graphique des scores LISI par lot et par méthode
  ggplot(lot_lisi_df, aes(x = Lot, y = LISI_score, color = Method)) +
    geom_point() +
    theme_minimal() +
    labs(title = "Scores LISI par lot et méthode d'intégration",
         x = "Lot", y = "Score LISI") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~ Method)
}

