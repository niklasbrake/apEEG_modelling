function [tau,R] = compute_uAP_xcorr(n)
% [tau,R] = COMPUTE_UAP_XCORR(n) returns the lags (tau) and cross-correlation matrix (R) for
%   the unitary AP responses of 1035 neuron models. Estimates cross correlation for n sampled
%   neuron pairs. Default n=10000.
%   R is an array with size 3 x 3 x n x length(tau).

    baseFolder = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(baseFolder,'auxiliary_functions'));

    load(fullfile(baseFolder,'data_files','unitaryAP.mat'),'savedUnitaryAP');

    if(nargin<10)
        n = 10000;
    end

    lag = 160;
    T = 2*lag+1;

    load(fullfile(baseFolder,'data_files','mtype_abundance.mat'),'mTypeCDF');
    [~,idcs] = unique(mTypeCDF);
    sample_blue_neurons = @(N) interp1(mTypeCDF(idcs),idcs,rand(N,1),'next','extrap');

    bs_sample = reshape(sample_blue_neurons(n*2),2,[]);
    R_xx = zeros(T,n); R_yx = zeros(T,n); R_zx = zeros(T,n);
    R_xy = zeros(T,n); R_yy = zeros(T,n); R_zy = zeros(T,n);
    R_xz = zeros(T,n); R_yz = zeros(T,n); R_zz = zeros(T,n);

    count = 1;
    while count < n
        i = bs_sample(1,count);
        j = bs_sample(2,count);
        R_xx(:,count) = xcorr(savedUnitaryAP(:,1,i),savedUnitaryAP(:,1,j),lag,'none')/16e3;
        R_xy(:,count) = xcorr(savedUnitaryAP(:,1,i),savedUnitaryAP(:,2,j),lag,'none')/16e3;
        R_xz(:,count) = xcorr(savedUnitaryAP(:,1,i),savedUnitaryAP(:,3,j),lag,'none')/16e3;
        R_yx(:,count) = xcorr(savedUnitaryAP(:,2,i),savedUnitaryAP(:,1,j),lag,'none')/16e3;
        R_yy(:,count) = xcorr(savedUnitaryAP(:,2,i),savedUnitaryAP(:,2,j),lag,'none')/16e3;
        R_yz(:,count) = xcorr(savedUnitaryAP(:,2,i),savedUnitaryAP(:,3,j),lag,'none')/16e3;
        R_zx(:,count) = xcorr(savedUnitaryAP(:,3,i),savedUnitaryAP(:,1,j),lag,'none')/16e3;
        R_zy(:,count) = xcorr(savedUnitaryAP(:,3,i),savedUnitaryAP(:,2,j),lag,'none')/16e3;
        R_zz(:,count) = xcorr(savedUnitaryAP(:,3,i),savedUnitaryAP(:,3,j),lag,'none')/16e3;
        count = count+1;

        k = i;
        i = j;
        j = k;
        R_xx(:,count) = xcorr(savedUnitaryAP(:,1,i),savedUnitaryAP(:,1,j),lag,'none')/16e3;
        R_xy(:,count) = xcorr(savedUnitaryAP(:,1,i),savedUnitaryAP(:,2,j),lag,'none')/16e3;
        R_xz(:,count) = xcorr(savedUnitaryAP(:,1,i),savedUnitaryAP(:,3,j),lag,'none')/16e3;
        R_yx(:,count) = xcorr(savedUnitaryAP(:,2,i),savedUnitaryAP(:,1,j),lag,'none')/16e3;
        R_yy(:,count) = xcorr(savedUnitaryAP(:,2,i),savedUnitaryAP(:,2,j),lag,'none')/16e3;
        R_yz(:,count) = xcorr(savedUnitaryAP(:,2,i),savedUnitaryAP(:,3,j),lag,'none')/16e3;
        R_zx(:,count) = xcorr(savedUnitaryAP(:,3,i),savedUnitaryAP(:,1,j),lag,'none')/16e3;
        R_zy(:,count) = xcorr(savedUnitaryAP(:,3,i),savedUnitaryAP(:,2,j),lag,'none')/16e3;
        R_zz(:,count) = xcorr(savedUnitaryAP(:,3,i),savedUnitaryAP(:,3,j),lag,'none')/16e3;
        count = count+1;
    end


    R = zeros(3,3,n,T);
    R(1,1,:,:) = R_xx';
    R(1,2,:,:) = R_xy';
    R(1,3,:,:) = R_xz';
    R(2,1,:,:) = R_yx';
    R(2,2,:,:) = R_yy';
    R(2,3,:,:) = R_yz';
    R(3,1,:,:) = R_zx';
    R(3,2,:,:) = R_zy';
    R(3,3,:,:) = R_zz';
    R = permute(R,[4,3,1,2]);
    tau = (-lag:lag)/16;
    tau = tau(:);
end