#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --account=def-akhadra
#SBATCH --mem=4G
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL,END

module load StdEnv/2020  intel/2020.1.217  openmpi/4.0.3
module load scipy-stack
module load python
module load neuron
module load mpi4py
pip install LFPykit

# Folder containing 1035 subdirectories, consiting of the Blue Brain neuron models, equivalent to local directory apEEG_modelling/data_files/neuron_models.
# Neuron models downloaded from the following link: https://bbp.epfl.ch/nmc-portal/assets/documents/static/Download/hoc_combos_syn.1_0_10.allzips.tar
folder=/lustre04/scratch/nbrake/data/simulations/unitary_AP
cd $folder

# The file changedEI.txt contains the names of all the neuron models to simulation, each
# on a seperate line
file=/lustre04/scratch/nbrake/code/simulate_blue_brain/changedEI.txt

# Runs all simulations across the 10 cores (use --array=1-10) by skipping every 10 lines
# of the changedEI.txt
lineNumber=${SLURM_ARRAY_TASK_ID}
((lineNumber--))
while IFS="" read -r d || [ -n "$d" ]
do
  echo $((lineNumber % 10))
  if [ $((lineNumber % 10)) -eq 0 ]; then
      cd "$d"
      nrnivmodl ./mechanisms # This can be removed after the first run
      python /lustre04/scratch/nbrake/code/simulate_blue_brain/simulate_models.py $d
      cd ..
  fi
  ((lineNumber++))
done < $file