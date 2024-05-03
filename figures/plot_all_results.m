% Check that the data files have been downloaded
baseFolder = fileparts(fileparts(mfilename('fullpath')));

fileCheck = true;

disp('Checking data files')
disp('--------------------')
datafile = fullfile(baseFolder,'data_files','AP_scaling.mat');
if(exist(datafile))
    disp('AP_scaling.mat....pass');
else
    disp('AP_scaling.mat....fail');
    fileCheck = false;
end

datafile = fullfile(baseFolder,'data_files','EI_ratio.mat');
if(exist(datafile))
    disp('EI_ratio.mat....pass');
else
    disp('EI_ratio.mat....fail');
    fileCheck = false;
end

datafile = fullfile(baseFolder,'data_files','MC_results.mat');
if(exist(datafile))
    disp('MC_results.mat....pass');
else
    disp('MC_results.mat....fail');
    fileCheck = false;
end

datafile = fullfile(baseFolder,'data_files','anatomy_nyhead_model.mat');
if(exist(datafile))
    disp('anatomy_nyhead_model.mat....pass');
else
    disp('anatomy_nyhead_model.mat....fail');
    fileCheck = false;
end

datafile = fullfile(baseFolder,'data_files','asymmetry_indices.mat');
if(exist(datafile))
    disp('asymmetry_indices.mat....pass');
else
    disp('asymmetry_indices.mat....fail');
    fileCheck = false;
end

datafile = fullfile(baseFolder,'data_files','mtype_abundance.mat');
if(exist(datafile))
    disp('mtype_abundance.mat....pass');
else
    disp('mtype_abundance.mat....fail');
    fileCheck = false;
end

datafile = fullfile(baseFolder,'data_files','pairwise_distance.mat');
if(exist(datafile))
    disp('pairwise_distance.mat....pass');
else
    disp('pairwise_distance.mat....fail');
    fileCheck = false;
end

datafile = fullfile(baseFolder,'data_files','unitaryAP.mat');
if(exist(datafile))
    disp('unitaryAP.mat....pass');
else
    disp('unitaryAP.mat....fail');
    fileCheck = false;
end

datafile = fullfile(baseFolder,'data_files','unitarySpectrum.mat');
if(exist(datafile))
    disp('unitarySpectrum.mat....pass');
else
    disp('unitarySpectrum.mat....fail');
    fileCheck = false;
end


if(~fileCheck)
    error('Some or all data files are missing. These can be downloaded from the following link and moved to the data_files directory: https://drive.google.com/uc?export=download&id=1Ek9COzFk_wjMBEZs1V88iNplqI2bCBFh')
else
    disp('All data files present. Proceeding to plot figures...')
end

% Run functions for all plots
figure1; drawnow;
figure2; drawnow;
figure3; drawnow;
figure4; drawnow;
figure5; drawnow;
figure6; drawnow;
figure7; drawnow;
figureS1; drawnow;
figureS2; drawnow;
figureS3; drawnow;
figureS4; drawnow;
figureS5; drawnow;