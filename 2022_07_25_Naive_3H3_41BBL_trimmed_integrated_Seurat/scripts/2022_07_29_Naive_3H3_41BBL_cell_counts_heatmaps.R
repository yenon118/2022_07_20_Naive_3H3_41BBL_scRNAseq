#!/usr/bin/Rscript --vanilla
rm(list=ls())

set.seed(1)

library(png)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)

library(patchwork)
library(openxlsx)

library(pheatmap)

##################################################
# Constants/Variables
##################################################


##################################################
# Output folder
##################################################

output_path <- file.path("../output/CellCountsHeatmaps")

if(!dir.exists(output_path)){
  dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
  if(!dir.exists(output_path)){
	quit(status=1)
  }
}


##################################################
# Read in input file
##################################################

folder_path <- file.path("../output/CellCounts")

dat_file_path <- file.path(folder_path, "cell_counts.xlsx")

dat_sheet_names <- getSheetNames(dat_file_path)

print(dat_sheet_names)


##################################################
# Process data
##################################################

plotHeatmap <- function(file_path_1, sheet) {

    dat1 <- read.xlsx(xlsxFile = file_path_1, sheet = sheet, skipEmptyRows = TRUE)

    colnames(dat1)[1] <- "Experiment"
    colnames(dat1)[2] <- "Cluster"

    dat1 = dat1 %>% 
        pivot_wider(names_from = Cluster, values_from = Cell_Count, values_fill = 0) %>%
        column_to_rownames(var = "Experiment") %>%
        as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)

    mat1 <- as.matrix(dat1)


    paletteLength <- 50

    myColor <- colorRampPalette(c("yellow", "red", "darkred"))(paletteLength)
    myColor <- c("#595959", myColor[1], myColor)

    myBreaks1 <- unique(
        c(
            seq(1, max(mat1, na.rm = TRUE), length.out=floor(paletteLength) )
        )
    )
    myBreaks1 <- c(0, 0.9999, myBreaks1)

    out_pheatmap <- pheatmap(
        mat1,
        fontsize=10,
        color=myColor,
        breaks=myBreaks1,
        cluster_rows=ifelse(nrow(mat1)>1, TRUE, FALSE),
        cluster_cols=ifelse(ncol(mat1)>1, TRUE, FALSE),
        fontsize_col=15,
        fontsize_row=15,
        filename=file.path(output_path, paste0("cell_counts_", sheet, ".png")),
        display_numbers = TRUE,
        number_format = "%d",
        fontsize_number = 10,
        width = 14,
        height = 7
    )

    clustered_mat1 <- mat1[out_pheatmap$tree_row$order, out_pheatmap$tree_col$order]

    return(clustered_mat1)
}

clustered_mat_seurat_clusters <- plotHeatmap(dat_file_path, "seurat_clusters")
clustered_mat_HumanPrimaryCellAtlasData <- plotHeatmap(dat_file_path, "HumanPrimaryCellAtlasData")
clustered_mat_BlueprintEncodeData <- plotHeatmap(dat_file_path, "BlueprintEncodeData")
clustered_mat_MouseRNAseqData <- plotHeatmap(dat_file_path, "MouseRNAseqData")
clustered_mat_ImmGenData <- plotHeatmap(dat_file_path, "ImmGenData")
clustered_mat_DbImmuneCellExpressionData <- plotHeatmap(dat_file_path, "DbImmuneCellExpressionData")
clustered_mat_NovershternHematopoieticData <- plotHeatmap(dat_file_path, "NovershternHematopoieticData")
clustered_mat_MonacoImmuneData <- plotHeatmap(dat_file_path, "MonacoImmuneData")


##################################################
# Save cell count data
##################################################

wb <- createWorkbook()

addWorksheet(wb, "seurat_clusters")
addWorksheet(wb, "HumanPrimaryCellAtlasData")
addWorksheet(wb, "BlueprintEncodeData")
addWorksheet(wb, "MouseRNAseqData")
addWorksheet(wb, "ImmGenData")
addWorksheet(wb, "DbImmuneCellExpressionData")
addWorksheet(wb, "NovershternHematopoieticData")
addWorksheet(wb, "MonacoImmuneData")

writeData(wb, "seurat_clusters", clustered_mat_seurat_clusters, rowNames = TRUE)
writeData(wb, "HumanPrimaryCellAtlasData", clustered_mat_HumanPrimaryCellAtlasData, rowNames = TRUE)
writeData(wb, "BlueprintEncodeData", clustered_mat_BlueprintEncodeData, rowNames = TRUE)
writeData(wb, "MouseRNAseqData", clustered_mat_MouseRNAseqData, rowNames = TRUE)
writeData(wb, "ImmGenData", clustered_mat_ImmGenData, rowNames = TRUE)
writeData(wb, "DbImmuneCellExpressionData", clustered_mat_DbImmuneCellExpressionData, rowNames = TRUE)
writeData(wb, "NovershternHematopoieticData", clustered_mat_NovershternHematopoieticData, rowNames = TRUE)
writeData(wb, "MonacoImmuneData", clustered_mat_MonacoImmuneData, rowNames = TRUE)

saveWorkbook(wb, file.path(output_path, "cell_counts_heatmaps.xlsx"), TRUE)
