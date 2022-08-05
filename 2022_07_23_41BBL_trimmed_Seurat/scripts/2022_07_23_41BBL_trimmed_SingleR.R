#!/usr/bin/Rscript --vanilla
rm(list=ls())

set.seed(1)

library(png)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)

library(Matrix)
library(Seurat)
library(patchwork)
library(celldex)
library(SingleR)


##################################################
# Constants/Variables
##################################################


##################################################
# Output folder
##################################################

output_path <- file.path("../output/SingleR")

if(!dir.exists(output_path)){
  dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
  if(!dir.exists(output_path)){
	quit(status=1)
  }
}


##################################################
# Read in input file
##################################################

dat_1_sample <- readRDS(
    file = file.path("../output/Seurat/data.rds")
)


refHumanPrimaryCellAtlasData <- celldex::HumanPrimaryCellAtlasData(ensembl = FALSE)
refBlueprintEncodeData <- celldex::BlueprintEncodeData(ensembl = FALSE)
refMouseRNAseqData <- celldex::MouseRNAseqData(ensembl = FALSE)
refImmGenData <- celldex::ImmGenData(ensembl = FALSE)
refDatabaseImmuneCellExpressionData <- celldex::DatabaseImmuneCellExpressionData(ensembl = FALSE)
refNovershternHematopoieticData <- celldex::NovershternHematopoieticData(ensembl = FALSE)
refMonacoImmuneData <- celldex::MonacoImmuneData(ensembl = FALSE)


##################################################
# Perform Prediction
##################################################

mapClusterToCellTypes <- function(dat, ref) {
  # Convert Seurat object to Single Cell Experiment class
  dat.sce <- as.SingleCellExperiment(dat)

  predictions <- SingleR(
    test=dat.sce,
    ref=ref, 
    labels=ref$label.main,
    # labels=ref$label.fine,
    assay.type.test=1,
    clusters = Idents(dat)
  )

  # return(predictions$pruned.labels)
  return(predictions$labels)
}


dat_1_labels_HumanPrimaryCellAtlasData <- mapClusterToCellTypes(dat_1_sample, refHumanPrimaryCellAtlasData)
dat_1_labels_BlueprintEncodeData <- mapClusterToCellTypes(dat_1_sample, refBlueprintEncodeData)
dat_1_labels_MouseRNAseqData <- mapClusterToCellTypes(dat_1_sample, refMouseRNAseqData)
dat_1_labels_ImmGenData <- mapClusterToCellTypes(dat_1_sample, refImmGenData)
dat_1_labels_DatabaseImmuneCellExpressionData <- mapClusterToCellTypes(dat_1_sample, refDatabaseImmuneCellExpressionData)
dat_1_labels_NovershternHematopoieticData <- mapClusterToCellTypes(dat_1_sample, refNovershternHematopoieticData)
dat_1_labels_MonacoImmuneData <- mapClusterToCellTypes(dat_1_sample, refMonacoImmuneData)

cluster_numeric <- as.numeric(levels(dat_1_sample[["seurat_clusters"]][,1]))[dat_1_sample[["seurat_clusters"]][,1]]
cluster_numeric_idx <- cluster_numeric + 1

dat_1_sample[["HumanPrimaryCellAtlasData"]] <- dat_1_labels_HumanPrimaryCellAtlasData[cluster_numeric_idx]
dat_1_sample[["BlueprintEncodeData"]] <- dat_1_labels_BlueprintEncodeData[cluster_numeric_idx]
dat_1_sample[["MouseRNAseqData"]] <- dat_1_labels_MouseRNAseqData[cluster_numeric_idx]
dat_1_sample[["ImmGenData"]] <- dat_1_labels_ImmGenData[cluster_numeric_idx]
dat_1_sample[["DatabaseImmuneCellExpressionData"]] <- dat_1_labels_DatabaseImmuneCellExpressionData[cluster_numeric_idx]
dat_1_sample[["NovershternHematopoieticData"]] <- dat_1_labels_NovershternHematopoieticData[cluster_numeric_idx]
dat_1_sample[["MonacoImmuneData"]] <- dat_1_labels_MonacoImmuneData[cluster_numeric_idx]


##################################################
# Plotting
##################################################

plotDimPlot <- function(dat, dat_name = "NoName", reduction = "umap", ref = "HumanPrimaryCellAtlasData") {
  p <- DimPlot(dat, reduction = reduction, label = TRUE, pt.size = 0.5, group.by = ref) + NoLegend()

  if (reduction == "umap") {
    filename = paste0("DimPlot_UMAP_singler_", dat_name, "_", ref, ".png")
  } else if (reduction == "tsne") {
    filename = paste0("DimPlot_tSNE_singler_", dat_name, "_", ref, ".png")
  }

  ggsave(
    filename = filename,
    plot = p,
    path = output_path,
    width = 14,
    height = 7
  )
}

plotDimPlot(dat_1_sample, "1", "umap", "HumanPrimaryCellAtlasData")
plotDimPlot(dat_1_sample, "1", "umap", "BlueprintEncodeData")
plotDimPlot(dat_1_sample, "1", "umap", "MouseRNAseqData")
plotDimPlot(dat_1_sample, "1", "umap", "ImmGenData")
plotDimPlot(dat_1_sample, "1", "umap", "DatabaseImmuneCellExpressionData")
plotDimPlot(dat_1_sample, "1", "umap", "NovershternHematopoieticData")
plotDimPlot(dat_1_sample, "1", "umap", "MonacoImmuneData")

plotDimPlot(dat_1_sample, "1", "tsne", "HumanPrimaryCellAtlasData")
plotDimPlot(dat_1_sample, "1", "tsne", "BlueprintEncodeData")
plotDimPlot(dat_1_sample, "1", "tsne", "MouseRNAseqData")
plotDimPlot(dat_1_sample, "1", "tsne", "ImmGenData")
plotDimPlot(dat_1_sample, "1", "tsne", "DatabaseImmuneCellExpressionData")
plotDimPlot(dat_1_sample, "1", "tsne", "NovershternHematopoieticData")
plotDimPlot(dat_1_sample, "1", "tsne", "MonacoImmuneData")


##################################################
# Convert meta.data in Seurat object to data frame
##################################################

dat_1_sample_table <- dat_1_sample[[]] %>% 
  as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)


##################################################
# Save data
##################################################
saveRDS(dat_1_sample, file = file.path(output_path, "data.rds"))
