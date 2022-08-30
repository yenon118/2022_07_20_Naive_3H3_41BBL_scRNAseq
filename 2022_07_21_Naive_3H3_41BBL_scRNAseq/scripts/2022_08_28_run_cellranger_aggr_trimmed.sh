
for i in 3H3 Naive SA41BBL; do \
sbatch --account=xulab --partition=Lewis,BioCompute,hpc5,General \
--time=1-21:00 --mem-per-cpu=34G -n 3 -N 1 \
--job-name=2022_08_28_cellranger_aggr_trimmed_${i} \
--output=log_2022_08_28_cellranger_aggr_trimmed_${i}_%j.out \
--wrap="cd /storage/htc/joshilab/yenc/projects/2022_07_22_Shirwan/2022_07_20_Naive_3H3_41BBL_scRNAseq/2022_07_21_Naive_3H3_41BBL_scRNAseq/output/cellranger_aggr_trimmed && \
module load cellranger && \
cellranger aggr --id=${i} \
--csv=/storage/htc/joshilab/yenc/projects/2022_07_22_Shirwan/2022_07_20_Naive_3H3_41BBL_scRNAseq/2022_07_21_Naive_3H3_41BBL_scRNAseq/output/cellranger_aggr_trimmed/${i}.csv"; \
done;
