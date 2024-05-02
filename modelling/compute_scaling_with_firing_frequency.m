function [f0,psd_Y,psd_Yhat,R,B2] = compute_scaling_with_firing_frequency
% COMPUTE_SCALING_WITH_FIRING_FREQUENCY analyzes the unitary AP responses to output
%   the file AP_scaling.mat used by figure2.
%
% Data files required:
%   EI_ratio.mat
%   unitaryAP.mat
%
%  See also modelling\unitary_AP_response\simulations\main.sh

    warning('off','signal:findpeaks:largeMinPeakHeight');

    baseFolder = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(baseFolder,'auxiliary_functions'));

    matObj = matfile(fullfile(baseFolder,'data_files','EI_ratio.mat'));
    cellIDs = who(matObj);

    load(fullfile(baseFolder,'data_files','unitaryAP.mat'));
    files = cellfun(@(x)strrep(x,'-','_'),files,'UniformOutput',false);
    getUAP = @(id) savedUnitaryAP(2:end,3,find(strcmp(id,files)));

    EI_vec = {'01','1.5','2.1','3.1','4.5','6.6','9.7','14.1','20.6','30'};

    M = length(cellIDs);
    N = zeros(length(EI_vec),1);
    B2 = zeros(length(EI_vec),M);
    firingFrequency = zeros(length(EI_vec),M);

    psd_passive = zeros(4096,length(EI_vec));
    psd_active = zeros(4096,length(EI_vec));

    psd_Y = zeros(4096,length(EI_vec)*M);
    psd_Yhat = zeros(4096,length(EI_vec)*M);

    fs = 16e3; % Hz
    h = waitbar(0);
    for ii = 1:M
        meanAP = getUAP(cellIDs{ii});
        cell = matObj.(cellIDs{ii});

        for k = 1:length(EI_vec)
            [f0,~,psd] = eegfft(cell.passive.time*1e-3,detrend(cell.passive.dipoles(:,3,k)),0.5,0.4);
            psd_passive(:,k) = mean(psd,2);

            [f0,~,psd] = eegfft(cell.active.time*1e-3,detrend(cell.active.dipoles(:,3,k)),0.5,0.4);
            psd_active(:,k) = mean(psd,2);

            [y,x] = findpeaks(cell.active.voltage(:,k),'MinPeakHeight',0);
            N(k) = length(x);
        end

        n = length(meanAP);
        xdft = fft(meanAP);
        xdft = xdft(1:n/2+1);
        psdx = (1/(fs*n)) * abs(xdft).^2;
        psdx(2:end-1) = 2*psdx(2:end-1);
        freq = (0:fs/n:fs/2)';

        psd_unit = interp1(freq,psdx,f0,'linear','extrap');
        psd_unit = psd_unit*length(meanAP)/fs/1; % per second

        for k = 1:10
            X0 = psd_passive(:,k);
            X1 = psd_unit;
            Y = psd_active(:,k);
            idcs = find(and(and(Y>1e-19,X1>1e-19),f0<1e3));
            B2(k,ii) = fminbnd(@(B1) neglnlike(B1,X0(idcs),X1(idcs),Y(idcs)),0,N(k)*2);
            firingFrequency(k,ii) = N(k)/range(cell.active.time)*1e3;

            Yhat = X0+B2(k,ii)*X1;
            psd_Y(:,10*(ii-1)+k) = Y;
            psd_Yhat(:,10*(ii-1)+k) = Yhat;

            Y = log(Y);
            Yhat = log(Yhat);
            R(k,ii) = 1 - sum((Y-Yhat).^2)./sum((Y-mean(Y)).^2);
        end
    end
    delete(h)
end
function output = neglnlike(B1,X0,X1,y)
    model = X0 + B1*X1;
    output = sum(log(model) + y./model);
end