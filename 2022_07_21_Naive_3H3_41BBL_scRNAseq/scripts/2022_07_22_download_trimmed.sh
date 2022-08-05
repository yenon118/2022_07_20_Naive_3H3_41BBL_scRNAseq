for i in $(awk -F ' ' '{print $2}' /storage/htc/joshilab/yenc/projects/2022_07_22_Shirwan/2022_07_21_Naive_3H3_41BBL_scRNAseq/data/ShirwanH_01_SOL_md5s.txt | grep "trimmed"); do \
sbatch --account=xulab --partition=Lewis,BioCompute,hpc5,General \
--time=0-21:00 --mem-per-cpu=20G -n 3 -N 1 \
--job-name=2022_07_22_download --output=log_2022_07_22_download_%j.out \
--wrap="cd /storage/htc/joshilab/yenc/projects/2022_07_22_Shirwan/2022_07_21_Naive_3H3_41BBL_scRNAseq/data/trimmed && \
curl -o ${i} \"https://oc1.rnet.missouri.edu/index.php/s/8yTs48x9PZ7Em67/download?path=%2F&files=${i}\""; done;
