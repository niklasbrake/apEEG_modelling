To generate the manuscript figures, the following results and other data files need to be <a href="https://drive.google.com/uc?export=download&id=1Ek9COzFk_wjMBEZs1V88iNplqI2bCBFh">download from this link</a> and placed in this folder. Below is a description of the data files that are available for download.

**AP_scaling.mat** is a file that contains the output of the function compute_scaling_with_firing_frequency.m. It contains the following variables:
B2 : 10 x 68 double
    the scaling factors (beta) for 68 neuron models at 10 different EI ratios
R : 10 x 68 double
    R^2 values for the linear fits
firingFrequency : 10 x 68 double
    firing frequency of each model simulation
f0 : 4096 x 1 double
    frequency vector for spectra
psd_Y : 4096 x 680 double
    power spectrum for all the active simulations (68 models x 10 EI values)
psd_Yhat : 4096 x 680 double
    Modelled power spectra based on the unitary AP response of each neuron

**EI_ratio.mat** contains the simulations results for 68 neuron models simulated at each of 10 EI ratios. It contains a structure for every
neuron model. To reduce file size, only a minimal set of simulation results (5 models) are included to reproduce analysis.

**MC_results.mat** contains the results of the Monte Carlo simulations of the cross spectrum among unitary apEEG signals. It contains the following variables
dValues : 1 x 1000 double
    the pairwise distances at which the average cross spectrum was calculated
count_D : 1 x 1000 double
    the number of samples for each dValue that the cross spectrum was computed
Pxy_D : 8000 x 1000 double
    estimated total (summed) cross spectrum at each pairwise distance. Needs to be divided by count_D to get the average.
SSExy_D : 8000 x 1000 double
    sum of squared errors for each spectrum across the samples

**anatomy_nyhead_model.mat** is the data from the New York head model, and can also be downloaded from https://www.parralab.org/nyhead/sa_nyhead.mat

**asymmetry_indices.mat** is the output of the function calculate_dendrite_asymmetry.m. Contains the variables
asym_idx : 1035 x 1 double
    the asymmetry index for each of the 1035 neuron model morphologies

**mtype_abundance.mat** contains the abundance data for each neuron class. Contains,
mTypePDF : 1035 x 1 double
    the relative abundance of each neuron model
mTypeCDF : 1035 x 1 double
    the cummulative abundance of each neuron model
mtype_abundance : table
    contains rows for each neuron class, its relative abundance, and the count of the number of neuron models of each class

**pairwise_distance.mat** is the output of the function compute_cortical_area. Contains the following variables,
rValues : 1 x 100 double
    the radii of the balls that were evaluated
idcs : 1 x 5000 double
    the indicies of the vertex points of the New York cortex template that served as the center of the ball for each sample
total_area : 100 x 5000
    the total surface area of the cortex enclosed for each ball radius (x100) and each ball origin (x5000)

**unitaryAP.mat** contains the output of the function compute_uAP. Contains the following variables:
savedUnitaryAP : 8001 x 3 x 1035 double
    the unitary AP response vector for each of the 10335 neuron models. Spike occurs at index 2001 and has a timestep of 1/16 milliseconds
mtype : 1 x 1035 cell
    the cell class of each of the neuron models
ei_type : 1 x 1035 logical
    the EI type of each neuron. 0 is inhibitory and 1 is excitatory

**unitarySpectrum.mat** is the output of the function compute_uAP_spectrum. Contains the following variables:
psd : 4000 x 1035 double
    power spectrum of the unitary apEEG response for each neuron model (x1035)
freq : 4000 x 1 double
    frequency vector associated with the variable psd

