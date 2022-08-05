
sbatch --account=xulab --partition=Lewis,BioCompute,hpc5,General \
--time=1-21:00 --mem-per-cpu=20G -n 3 -N 1 \
--job-name=2022_08_04_Naive_41BBL_markers \
--output=log_2022_08_04_Naive_41BBL_markers_%j.out \
--wrap="cd /storage/htc/joshilab/yenc/projects/2022_07_22_Shirwan/2022_07_20_Naive_3H3_41BBL_scRNAseq/2022_07_24_Naive_41BBL_trimmed_integrated_Seurat/scripts && \
source activate seurat && \
Rscript 2022_08_04_Naive_41BBL_markers.R"
