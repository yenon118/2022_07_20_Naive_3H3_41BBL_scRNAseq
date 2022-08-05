
for i in $(awk -F ' ' '{print $2}' /storage/htc/joshilab/yenc/projects/2022_07_22_Shirwan/2022_07_21_Naive_3H3_41BBL_scRNAseq/data/ShirwanH_01_SOL_md5s.txt | grep -v trimmed | sed 's/.fastq.gz$//g' | sort | uniq); do \
sbatch --account=xulab --partition=Lewis,BioCompute,hpc5,General \
--time=0-21:00 --mem-per-cpu=20G -n 3 -N 1 \
--job-name=2022_07_22_fastqc_${i} --output=log_2022_07_22_fastqc_${i}_%j.out \
--wrap="mkdir -p /storage/htc/joshilab/yenc/projects/2022_07_22_Shirwan/2022_07_21_Naive_3H3_41BBL_scRNAseq/output/fastqc_original/${i} && \
cd /storage/htc/joshilab/yenc/projects/2022_07_22_Shirwan/2022_07_21_Naive_3H3_41BBL_scRNAseq/output/fastqc_original && \
module load fastqc && \
fastqc -o /storage/htc/joshilab/yenc/projects/2022_07_22_Shirwan/2022_07_21_Naive_3H3_41BBL_scRNAseq/output/fastqc_original/${i} -f fastq \
/storage/htc/joshilab/yenc/projects/2022_07_22_Shirwan/2022_07_21_Naive_3H3_41BBL_scRNAseq/data/original/${i}.fastq.gz"; done;
