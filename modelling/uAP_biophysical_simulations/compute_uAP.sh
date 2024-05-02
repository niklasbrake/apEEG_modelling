#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --account=def-akhadra
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL,END

module load matlab
matlab -nodisplay -r "compute_uAP"
