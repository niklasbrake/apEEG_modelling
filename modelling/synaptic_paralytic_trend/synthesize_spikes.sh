#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=16G
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL,END
##SBATCH --output=synthetic_spikes.log
#SBATCH --account=def-akhadra

module load matlab
matlab -nojvm -batch  "synthesize_spikes(${SLURM_ARRAY_TASK_ID})"

