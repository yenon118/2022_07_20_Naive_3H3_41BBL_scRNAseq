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

dat_2_samples <- readRDS(
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

DefaultAssay(dat_2_samples) <- "RNA"

print(head(dat_2_samples[[]]))
print(head(Idents(dat_2_samples)))

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

  print(head(predictions$labels))

  # return(predictions$pruned.labels)
  return(predictions$labels)
}


dat_2_labels_HumanPrimaryCellAtlasData <- mapClusterToCellTypes(dat_2_samples, refHumanPrimaryCellAtlasData)
dat_2_labels_BlueprintEncodeData <- mapClusterToCellTypes(dat_2_samples, refBlueprintEncodeData)
dat_2_labels_MouseRNAseqData <- mapClusterToCellTypes(dat_2_samples, refMouseRNAseqData)
dat_2_labels_ImmGenData <- mapClusterToCellTypes(dat_2_samples, refImmGenData)
dat_2_labels_DatabaseImmuneCellExpressionData <- mapClusterToCellTypes(dat_2_samples, refDatabaseImmuneCellExpressionData)
dat_2_labels_NovershternHematopoieticData <- mapClusterToCellTypes(dat_2_samples, refNovershternHematopoieticData)
dat_2_labels_MonacoImmuneData <- mapClusterToCellTypes(dat_2_samples, refMonacoImmuneData)


cluster_numeric <- as.numeric(levels(dat_2_samples[["seurat_clusters"]][,1]))[dat_2_samples[["seurat_clusters"]][,1]]
cluster_numeric_idx <- cluster_numeric + 1

dat_2_samples[["HumanPrimaryCellAtlasData"]] <- dat_2_labels_HumanPrimaryCellAtlasData[cluster_numeric_idx]
dat_2_samples[["BlueprintEncodeData"]] <- dat_2_labels_BlueprintEncodeData[cluster_numeric_idx]
dat_2_samples[["MouseRNAseqData"]] <- dat_2_labels_MouseRNAseqData[cluster_numeric_idx]
dat_2_samples[["ImmGenData"]] <- dat_2_labels_ImmGenData[cluster_numeric_idx]
dat_2_samples[["DatabaseImmuneCellExpressionData"]] <- dat_2_labels_DatabaseImmuneCellExpressionData[cluster_numeric_idx]
dat_2_samples[["NovershternHematopoieticData"]] <- dat_2_labels_NovershternHematopoieticData[cluster_numeric_idx]
dat_2_samples[["MonacoImmuneData"]] <- dat_2_labels_MonacoImmuneData[cluster_numeric_idx]


##################################################
# Plotting
##################################################

plotDimPlot <- function(dat, dat_name = "NoName", reduction = "umap", ref = "HumanPrimaryCellAtlasData") {
  p <- DimPlot(dat, reduction = reduction, label = TRUE, pt.size = 0.5, group.by = ref, raster = FALSE) + NoLegend()

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

plotDimPlot(dat_2_samples, "2", "umap", "HumanPrimaryCellAtlasData")
plotDimPlot(dat_2_samples, "2", "umap", "BlueprintEncodeData")
plotDimPlot(dat_2_samples, "2", "umap", "MouseRNAseqData")
plotDimPlot(dat_2_samples, "2", "umap", "ImmGenData")
plotDimPlot(dat_2_samples, "2", "umap", "DatabaseImmuneCellExpressionData")
plotDimPlot(dat_2_samples, "2", "umap", "NovershternHematopoieticData")
plotDimPlot(dat_2_samples, "2", "umap", "MonacoImmuneData")

plotDimPlot(dat_2_samples, "2", "tsne", "HumanPrimaryCellAtlasData")
plotDimPlot(dat_2_samples, "2", "tsne", "BlueprintEncodeData")
plotDimPlot(dat_2_samples, "2", "tsne", "MouseRNAseqData")
plotDimPlot(dat_2_samples, "2", "tsne", "ImmGenData")
plotDimPlot(dat_2_samples, "2", "tsne", "DatabaseImmuneCellExpressionData")
plotDimPlot(dat_2_samples, "2", "tsne", "NovershternHematopoieticData")
plotDimPlot(dat_2_samples, "2", "tsne", "MonacoImmuneData")


##################################################
# Convert meta.data in Seurat object to data frame
##################################################

dat_2_samples_table <- dat_2_samples[[]] %>% 
  as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)


##################################################
# Save data
##################################################
saveRDS(dat_2_samples, file = file.path(output_path, "data.rds"))
