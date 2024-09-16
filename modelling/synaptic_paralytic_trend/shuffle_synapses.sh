#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --account=def-akhadra
#SBATCH --mem=8G
#SBATCH --mail-user=niklas.brake@mail.mcgill.ca
#SBATCH --mail-type=FAIL,END
#SBATCH --output=simulation.log

module load StdEnv/2020  intel/2020.1.217  openmpi/4.0.3
module load scipy-stack
module load python
module load neuron/8.0.0
module load mpi4py
pip install LFPykit

module load matlab
matlab -nojvm -batch  "shuffle_synapses(${SLURM_ARRAY_TASK_ID})"

