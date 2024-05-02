## Overview
This repository contains four folders: **figures**, **modelling**, **auxiliary_functions**, and **data_files**.

**figures** contains MATLAB code to reproduce all the figures of the manuscript. See Installation guide and Instructions for use below.

**modelling** contains several analysis functions that are called by the figure generation scripts. Additionally, this folder contains the code used to simulate the model. The results of these simulations are  provided as data files (see below).

**auxiliary_functions** contains auxiliary MATLAB function unrelated to the modelling work.

**data_files** contains the output of the simulations, as well as other necessary data files. To reduce the repository size, these data files should be <a href="https://drive.google.com/uc?export=download&id=1Ek9COzFk_wjMBEZs1V88iNplqI2bCBFh">download from this link</a> and moved to the data_files directory. For more information, read the README in the data_files folder.

## System requirements

Generating the manuscript figures requires only a standard computer. This software has been tested using MATLAB 2023a running on Windows 10.

Required MATLAB toolboxes
+ Signal Processing Toolbox
+ Statistics and Machine Learning Toolbox

## Installation guide
```
git clone https://github.com/niklasbrake/apEEG_modelling
cd apEEG_modelling/figures
matlab
```
This should run in less than a minute and open a matlab instance with the working directory set to the **figures** folder.

## Instructions for use
Quantiative results in our manuscript can be reproduced by running the script
```
plot_all_results
```
This script should produce 26 figures that together reproduce the modelling results in Figs. 1-7 and Figs. S1-S5. Estimate run time ~3.5 minutes.

## License
<a rel="license" href="https://www.gnu.org/licenses/gpl-3.0.html#license-text"><img alt="GPL3" style="border-width:0" src="https://www.gnu.org/graphics/gplv3-127x51.png" /></a><br />This repository is licensed under a <a rel="license" href="https://www.gnu.org/licenses/gpl-3.0.html#license-text">GNU GENERAL PUBLIC LICENSE Version 3</a>.
