#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL,END
##SBATCH --output=synthetic_spikes.log
#SBATCH --account=def-akhadra

module load matlab
matlab -nojvm -batch "project_covariance_matrix(${SLURM_ARRAY_TASK_ID})"