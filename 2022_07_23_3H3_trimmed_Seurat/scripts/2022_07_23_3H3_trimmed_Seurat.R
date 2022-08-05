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
library(scuttle)


##################################################
# Constants/Variables
##################################################

chosen_dim <- 30


##################################################
# Output folder
##################################################

output_path <- file.path("../output/Seurat")

if(!dir.exists(output_path)){
  dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
  if(!dir.exists(output_path)){
	quit(status=1)
  }
}


##################################################
# Read in input file
##################################################

folder_path = file.path("../../2022_07_21_Naive_3H3_41BBL_scRNAseq/output/cellranger_count_trimmed")


directories = list.dirs(folder_path, full.names = FALSE, recursive = FALSE)
directories = directories[str_starts(directories, pattern = fixed("3H3"))]

print(directories)


dat_list = list()

create_data_list <- function(dat_list, project, data_dir){
    dat_10xobj <- Read10X(data.dir = data_dir)
    dat_list[[project]] <- CreateSeuratObject(counts = dat_10xobj, project = project, min.cells = 8, min.features = 200)
    return(dat_list)
}

dat_list <- create_data_list(dat_list, directories[1], file.path(folder_path, directories[1], "outs", "filtered_feature_bc_matrix"))
dat_list <- create_data_list(dat_list, directories[2], file.path(folder_path, directories[2], "outs", "filtered_feature_bc_matrix"))
dat_list <- create_data_list(dat_list, directories[3], file.path(folder_path, directories[3], "outs", "filtered_feature_bc_matrix"))
dat_list <- create_data_list(dat_list, directories[4], file.path(folder_path, directories[4], "outs", "filtered_feature_bc_matrix"))
dat_list <- create_data_list(dat_list, directories[5], file.path(folder_path, directories[5], "outs", "filtered_feature_bc_matrix"))

print(dat_list)


##################################################
# Process data
##################################################

cat("Start merging ...\n\n")

dat <- merge(
	dat_list[[1]], 
	y = c(
		dat_list[[2]],
		dat_list[[3]],
		dat_list[[4]],
		dat_list[[5]]
	), 
	add.cell.ids = directories, 
	project = "Naive_3H3_41BBL"
)

cat("Merge complete!!!\n\n")


##################################################
# Visualize QC metrics
##################################################

dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^mt-")

# Visualize QC metrics as a violin plot
p <- VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ggsave(
  filename = paste0("VlnPlot_init_QC.png"),
  plot = p,
  path = output_path,
  width = 21,
  height = 7
)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p <- plot1 + plot2

ggsave(
  filename = paste0("FeatureScatter_init_QC.png"),
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)


##################################################
# Subset
##################################################

dat <- subset(dat, subset = nFeature_RNA > 100 & nFeature_RNA < 7500 & percent.mt < 50)

print(dat)


##################################################
# Visualize after QC
##################################################

# Visualize QC metrics as a violin plot
p <- VlnPlot(dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ggsave(
  filename = paste0("VlnPlot_after_QC.png"),
  plot = p,
  path = output_path,
  width = 21,
  height = 7
)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p <- plot1 + plot2

ggsave(
  filename = paste0("FeatureScatter_after_QC.png"),
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)


##################################################
# Normalize and find top features
##################################################

dat <- NormalizeData(dat, normalization.method = "LogNormalize", scale.factor = 10000)

dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
print(head(VariableFeatures(object = dat)))

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(dat), 10)
print(top10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(dat)
plot2 <- LabelPoints(
	plot = plot1, 
	points = top10, 
	repel = TRUE, 
	xnudge = 0, 
	ynudge = 0, 
	size = 2, 
	point.size = 2,
	max.overlaps = 10
)
p <- plot1 + plot2

ggsave(
  filename = "VariableFeaturePlot.png",
  plot = plot1,
  path = output_path,
  width = 7,
  height = 7
)

ggsave(
  filename = "LabelPoints.png",
  plot = plot2,
  path = output_path,
  width = 7,
  height = 7
)

ggsave(
  filename = "VariableFeaturePlot_and_LabelPoints.png",
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)


##################################################
# Cell Cycle Scoring
##################################################

# s.genes <-cc.genes$s.genes
# g2m.genes<-cc.genes$g2m.genes

# dat <- CellCycleScoring(dat, s.features = s.genes, g2m.features = g2m.genes)


##################################################
# Scale data
##################################################

all.genes <- rownames(dat)

# dat <- ScaleData(dat, vars.to.regress = c("S.Score", "G2M.Score"), features = all.genes)
dat <- ScaleData(dat, features = all.genes)


##################################################
# Run PCA
##################################################

dat <- RunPCA(dat, features = VariableFeatures(object = dat))

# Examine and visualize PCA results a few different ways
print(dat[["pca"]], dims = 1:5, nfeatures = 5)

p <- VizDimLoadings(dat, dims = 1:2, reduction = "pca")

ggsave(
  filename = "VizDimLoadings_PCA.png",
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)

p <- DimPlot(dat, reduction = "pca")

ggsave(
  filename = "DimPlot_PCA.png",
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)

cat(rep("\n", 2))
png(filename = file.path(output_path, "DimHeatmap_PC1.png"), width = 10, height = 10, units = "in", res = 300)
DimHeatmap(dat, dims = 1, cells = 500, balanced = TRUE)
dev.off()

cat(rep("\n", 2))
png(filename = file.path(output_path, "DimHeatmap_PC1_PC15.png"), width = 30, height = 20, units = "in", res = 300)
DimHeatmap(dat, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()


##################################################
# Calculate using JackStraw and ScoreJackStraw
##################################################

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
# dat <- JackStraw(dat, num.replicate = 500, dims = 50)
# dat <- ScoreJackStraw(dat, dims = 1:50)

# p <- JackStrawPlot(dat, dims = 1:50)

# ggsave(
#   filename = "JackStrawPlot.png",
#   plot = p,
#   path = output_path,
#   width = 14,
#   height = 7
# )

##################################################
# Plot ElbowPlot
##################################################

p <- ElbowPlot(dat, ndims = 50)

ggsave(
  filename = "ElbowPlot.png",
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)


##################################################
# Cluster the cells
##################################################

dat <- FindNeighbors(dat, reduction = "pca", dims = 1:chosen_dim)
dat <- FindClusters(dat, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
print(head(Idents(dat), 5))


##################################################
# Run non-linear dimensional reduction (UMAP/tSNE)
##################################################

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
dat <- RunUMAP(dat, reduction = "pca", dims = 1:chosen_dim)
dat <- RunTSNE(dat, reduction = "pca", dims = 1:chosen_dim)

# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
p <- DimPlot(dat, reduction = "umap", label = TRUE, label.size = 8)

ggsave(
  filename = "DimPlot_UMAP.png",
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)

p <- DimPlot(dat, reduction = "tsne", label = TRUE, label.size = 8)

ggsave(
  filename = "DimPlot_tSNE.png",
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)

p <- DimPlot(dat, reduction = "umap", group.by = 'orig.ident', pt.size = 0.1)

ggsave(
  filename = "DimPlot_UMAP_orig_ident.png",
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)

p <- DimPlot(dat, reduction = "tsne", group.by = 'orig.ident', pt.size = 0.1)

ggsave(
  filename = "DimPlot_tSNE_orig_ident.png",
  plot = p,
  path = output_path,
  width = 14,
  height = 7
)


##################################################
# Save data
##################################################
saveRDS(dat, file = file.path(output_path, "data.rds"))


##################################################
# Finding differentially expressed features (cluster biomarkers)
##################################################

# find all markers of cluster 2
cluster2.markers <- FindMarkers(dat, ident.1 = 2, min.pct = 0.25)
print(head(cluster2.markers, n = 5))

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(dat, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
print(head(cluster5.markers, n = 5))

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
dat.markers <- FindAllMarkers(dat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

df <- dat.markers %>% 
	as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

write.table(
  x = df,
  sep = "\t",
  file = file.path(output_path, "data_markers.txt"),
  na = "",
  row.names = FALSE,
  quote = FALSE
)

found_markers <- dat.markers %>% 
	group_by(cluster) %>% 
	slice_max(n = 5, order_by = avg_log2FC) %>% 
	as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

cluster0.markers <- FindMarkers(dat, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
print(head(cluster0.markers, n = 5))
