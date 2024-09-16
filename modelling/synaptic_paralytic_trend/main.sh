# 7,8,10,13,14,24,26,28,43,50,54
ID1=$(sbatch --array=7 --parsable calculate_covariance_matrix.sh)
ID2=$(sbatch --depend=afterok:$ID1 --parsable --array=7 project_covariance_matrix.sh)
ID3=$(sbatch --depend=afterok:$ID2 --parsable --array=7 synthesize_spikes.sh)