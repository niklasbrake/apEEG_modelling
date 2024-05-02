<< comment
This script simulates all the Blue Brain project models, calculates their
unitary AP responses and the spectra of these responses. This code was run
on Calcul Quebec's cluster "Beluga", which uses SLURM workload management.
comment

# Directory containing code on server. Equivalent to
#   apEEG_modelling/modelling/uAP_biophysical_simulations
cd /lustre04/scratch/nbrake/code/simulate_blue_brain

# This command runs the script simulate_models across 10 cores in parallel
ID2=$(sbatch --array=1-10 --parsable simulate_models.sh)

# After the first job finsihes, this command will compute the unitary AP
# resopnses from the simulation results
ID3=$(sbatch --depend=afterok:$ID2 --parsable compute_uAP.sh)

# After the unitary AP responses are calculated the unitary apEEG spectrum
# is computed with this command
sbatch --depend=afterok:$ID3 compute_uapEEG_spectrum.sh)