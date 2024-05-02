#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --account=def-akhadra
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=0
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL,END
#SBATCH --output=MC_cross_spectra.log

module load matlab
matlab -nodisplay -r "run_MC_cross_spectra"
