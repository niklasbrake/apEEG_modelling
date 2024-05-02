#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --account=def-akhadra
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL,END
#SBATCH --output=synthetic_spikes.log

module load matlab
srun matlab -nodisplay -r "compute_uapEEG_spectrum"