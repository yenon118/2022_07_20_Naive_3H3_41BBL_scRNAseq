
for i in $(awk -F ' ' '{print $2}' /storage/htc/joshilab/yenc/projects/2022_07_22_Shirwan/2022_07_21_Naive_3H3_41BBL_scRNAseq/data/ShirwanH_01_SOL_md5s.txt | grep -v trimmed | sed 's/_S.*//g' | sort | uniq); do \
sbatch --account=xulab --partition=Lewis,BioCompute,hpc5,General \
--time=1-21:00 --mem-per-cpu=10G -n 8 -N 1 \
--job-name=2022_07_23_cellranger_count_original_${i} --output=log_2022_07_23_cellranger_count_original_${i}_%j.out \
--wrap="mkdir -p /storage/htc/joshilab/yenc/projects/2022_07_22_Shirwan/2022_07_21_Naive_3H3_41BBL_scRNAseq/output/cellranger_count_original && \
cd /storage/htc/joshilab/yenc/projects/2022_07_22_Shirwan/2022_07_21_Naive_3H3_41BBL_scRNAseq/output/cellranger_count_original && \
module load cellranger && \
cellranger count --id=${i} \
--fastqs=/storage/htc/joshilab/yenc/projects/2022_07_22_Shirwan/2022_07_21_Naive_3H3_41BBL_scRNAseq/data/original \
--sample=${i} \
--localcores=8 \
--transcriptome=/storage/htc/joshilab/yenc/projects/2022_07_22_Shirwan/2022_07_21_Naive_3H3_41BBL_scRNAseq/reference/GRCm39_Ensembl"; done;
