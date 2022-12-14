
sbatch --account=xulab --partition=Lewis,BioCompute,hpc5,General \
--time=1-21:00 --mem-per-cpu=20G -n 3 -N 1 \
--job-name=2022_07_25_Naive_3H3_41BBL_trimmed_integrated_SingleR \
--output=log_2022_07_25_Naive_3H3_41BBL_trimmed_integrated_SingleR_%j.out \
--wrap="cd /storage/htc/joshilab/yenc/projects/2022_07_22_Shirwan/2022_07_20_Naive_3H3_41BBL_scRNAseq/2022_08_28_Naive_3H3_41BBL_trimmed_integrated_Seurat/scripts && \
source activate seurat && \
Rscript 2022_07_25_Naive_3H3_41BBL_trimmed_integrated_SingleR.R"
