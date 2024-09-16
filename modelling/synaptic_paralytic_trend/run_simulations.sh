ID1=$(sbatch --array=7,8,10,13,14,24,26,28,43,50,54 --parsable synthetic_spikes.sh)
sbatch --depend=afterok:$ID1 --array=7,8,10,13,14,24,26,28,43,50,54 shuffle_synapses.sh
