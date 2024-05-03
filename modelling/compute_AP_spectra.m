function [f,Sxx,Sxy] = compute_AP_spectra(sigx2)
% [f,Sxx,Sxy] = COMPUTE_AP_SPECTRA takes the simulations of unitary AP and
%   the MC simulation to calculates the expected spectra, Sxx, and cross-spectra, Sxy.
%
% Data files required:
%   pairwise_distance.mat
%   mtype_abundance.mat
%   unitarySpectrum.mat
%   MC_results.mat
%
% See also compute_cortical_area, modelling/unitary_AP_response/simulations/main.sh, and modelling/MC_simulations/main.sh.

    baseFolder = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(baseFolder,'auxiliary_functions'));

    if(nargin<1)
        sigx2 = 3;
    end

    fs = 16e3;
    L = 16e3;
    f = (fs/L:fs/L:fs/2)';

    load(fullfile(baseFolder,'data_files','mtype_abundance.mat'),'mTypePDF');
    load(fullfile(baseFolder,'data_files','unitarySpectrum.mat'),'psd');
    Sxx = zeros(8e3,1);
    Sxx(1:2:end) = sum(psd.*mTypePDF',2); % Force to same frequency resolution as cross spectrum
    Sxx(2:2:end) = sum(psd.*mTypePDF',2); % Force to same frequency resolution as cross spectrum

    % Load average cross-spectrum for different pairwise distances
    load(fullfile(baseFolder,'data_files','MC_results.mat'),'Pxy_D','count_D','dValues');
    mu = Pxy_D./count_D;

    % Get density of surface area as a function of distance
    load(fullfile(baseFolder,'data_files','pairwise_distance.mat'));
    N2 = 16e9;
    dV = mean(diff(dValues));
    dN = interp1(0.5*(rValues(1:end-1)+rValues(2:end)),diff(total_area)./diff(rValues'),dValues,'linear','extrap');
    dN = mean(dN,2)'.*N2./mean(total_area(end,:));
    A = exp(-dValues.^2/(2*sigx2)).*dN*dV/N2;

    % Compute corrected average cross spectrum
    % Smooth over simulation noise to better represent spectral trend
    Sxy = smooth(nansum(mu.*A,2));
end