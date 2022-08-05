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
library(biomaRt)
library(SingleR)
library(openxlsx)

library(ggvenn)

##################################################
# Constants/Variables
##################################################


##################################################
# Output folder
##################################################

output_path <- file.path("../output/Overlap3AnnotationsVenn")

if(!dir.exists(output_path)){
  dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
  if(!dir.exists(output_path)){
	quit(status=1)
  }
}


##################################################
# Read in input file
##################################################

folder_path = file.path("../../")

dat_3H3 <- readRDS(file = file.path("../../2022_07_23_3H3_trimmed_Seurat/output/SingleR/data.rds"))
dat_41BBL <- readRDS(file = file.path("../../2022_07_23_41BBL_trimmed_Seurat/output/SingleR/data.rds"))
dat_Naive <- readRDS(file = file.path("../../2022_07_23_Naive_trimmed_Seurat/output/SingleR/data.rds"))


##################################################
# Process data
##################################################

generateVennForOverlapCells <- function(dat1, dat2, dat3, ref) {

    vec1 <- sort(unique(dat1[[ref]][,1]))
    vec2 <- sort(unique(dat2[[ref]][,1]))
    vec3 <- sort(unique(dat3[[ref]][,1]))

    a <- list(
        `3H3` = vec1,
        `41BBL` = vec2,
        Naive = vec3
    )

    df <- data.frame(
        X = sort(unique(c(vec1, vec2, vec3))),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    df_3H3 <- data.frame(
        X = vec1,
        `3H3` = vec1,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    df_41BBL <- data.frame(
        X = vec2,
        `41BBL` = vec2,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    df_Naive <- data.frame(
        X = vec3,
        Naive = vec3,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    df = df %>%
        left_join(df_3H3, by = "X") %>%
        left_join(df_41BBL, by = "X") %>%
        left_join(df_Naive, by = "X") %>%
        dplyr::select(-c(1)) %>%
        as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)


    p <- ggvenn(a, c("3H3", "41BBL", "Naive"))

    ggsave(
        filename = paste0("ann_venn_", ref, ".png"),
        plot = p,
        path = output_path,
        width = 14,
        height = 7
    )

    write.csv(
        x = df,
        file = file.path(output_path, paste0("ann_venn_", ref, ".csv")),
        na = "",
        quote = TRUE,
        row.names = FALSE
    )
}


generateVennForOverlapCells(
    dat_3H3,
    dat_41BBL,
    dat_Naive,
    "HumanPrimaryCellAtlasData"
)

generateVennForOverlapCells(
    dat_3H3,
    dat_41BBL,
    dat_Naive,
    "BlueprintEncodeData"
)

generateVennForOverlapCells(
    dat_3H3,
    dat_41BBL,
    dat_Naive,
    "MouseRNAseqData"
)

generateVennForOverlapCells(
    dat_3H3,
    dat_41BBL,
    dat_Naive,
    "ImmGenData"
)

generateVennForOverlapCells(
    dat_3H3,
    dat_41BBL,
    dat_Naive,
    "DatabaseImmuneCellExpressionData"
)

generateVennForOverlapCells(
    dat_3H3,
    dat_41BBL,
    dat_Naive,
    "NovershternHematopoieticData"
)

generateVennForOverlapCells(
    dat_3H3,
    dat_41BBL,
    dat_Naive,
    "MonacoImmuneData"
)
