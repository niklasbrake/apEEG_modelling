function run_MC_cross_spectra
% RUN_MC_CROSS_SPECTRA estimates the average cross-spectrum among unitary apEEG responses
%   for neurons seperated by different pairwise distances. Called by run_MC_cross_spectra.sh.

    M = 40;
    parpool(M)
    parfor chain_number = 1:M
        main(chain_number);
    end

    % Commpile results across the parallel workers
    SSExy_D = zeros(8e3,1e3);
    Pxy_D = zeros(8e3,1e3);
    count_D = zeros(1,1e3);
    for i = 1:M
       data = load(sprintf('/lustre04/scratch/nbrake/data/simulation_analyzed/cross_spectra/chain%d.mat',i))
       SSExy_D = SSExy_D+data.SSExy_D;
       Pxy_D = Pxy_D+data.Pxy_D;
       count_D = count_D+data.count_D;
    end

    save('/lustre04/scratch/nbrake/data/simulation_analyzed/cross_spectra/all_chains.mat','SSExy_D','Pxy_D','count_D');

end
function main(chain_number)
    saveFile = sprintf('/lustre04/scratch/nbrake/data/simulation_analyzed/cross_spectra/chain%d.mat',chain_number);

    % Load anatomy information
    load('/lustre04/scratch/nbrake/data/simulation_analyzed/cortical_area/MC_data.mat','A','dValues');
    B = cumsum(A)/sum(A);
    [~,iUniqueD] = unique(B);

    % Load savedUnitaryAP
    load('/lustre04/scratch/nbrake/data/simulation_analyzed/unitaryAP/unitaryAP.mat','savedUnitaryAP')
    savedUnitaryAP = permute(savedUnitaryAP,[2,1,3]);
    % Calculate new coordinates uiAi
    load('/lustre04/scratch/nbrake/data/anatomy_nyhead_model.mat');
    X = struct();
    X.vertices = sa.cortex75K.vc;
    X.faces= sa.cortex75K.tri;
    Cz = find(strcmp(sa.clab_electrodes,'Cz'));

    % Calculate normal lead fields
    M = size(X.vertices,1);
    Lxyz = zeros(M,3);
    for idx = 1:M
        L0 = squeeze(sa.cortex75K.V_fem(Cz,idx,:))'; % czIDX = 49
        vz = sa.cortex75K.normals(idx,:);
        [vx,vy] = getOrthBasis(vz);
        Lxyz(idx,:) = 1e-6*L0*[vx(:),vy(:),vz(:)];
    end

    % Set up a KD tree to efficiently search for nearby vertices
    Mdl = KDTreeSearcher(X.vertices);

    % Result variables
    L = 16e3;
    fs = 16e3;
    f = fs/L:fs/L:fs/2;
    Sxy = zeros(length(f),1);
    Pxy = zeros(length(f),1);
    SSExy = zeros(length(f),1);

    Pxy_D = zeros(length(f),length(dValues));
    SSExy_D = zeros(length(f),length(dValues));
    count_D = zeros(1,length(dValues));

    % Convergence variable
    iLargeLoop = 1;
    loopSize = 1e3;
    d_abs = 1e-4/(16e9)^2*sqrt(40);

    % Sample neurons proportional to abundance
    load('/lustre04/scratch/nbrake/data/simulation_analyzed/unitaryAP/mtype_abundance.mat','mTypeCDF');
    [~,idcs] = unique(mTypeCDF);
    sample_blue_neurons = @(N) interp1(mTypeCDF(idcs),idcs,rand(N,1),'next','extrap');

    iSampling = 1;
    samplingInterval = 1e4;
    neuronIdcs = reshape(sample_blue_neurons(2*samplingInterval),samplingInterval,2);
    dSamples = interp1(B(iUniqueD),dValues(iUniqueD),rand(samplingInterval,1),'next','extrap');
    N = 0;
    tic
    while true
        N = N+1;

        % Randomly sample unitary AP profiles
        i = neuronIdcs(iSampling,1);
        j = neuronIdcs(iSampling,2);
        d = dSamples(iSampling);

        % Randomly sample vertex from cortex
        iX = randi(M);

        % Sample second vertex closest to d away from first
        idcs = rangesearch(Mdl,X.vertices(iX,:),d+1);
        d2 = vecnorm(X.vertices(iX,:)-X.vertices(idcs{1},:),2,2);
        [~,I] = min(abs(d2-d));
        jX = idcs{1}(I);

        uiAi = Lxyz(iX,:);
        ujAj = Lxyz(jX,:);

        % Compute eeg cross spectrum
        eeg1 = uiAi*savedUnitaryAP(:,:,i);
        eeg2 = ujAj*savedUnitaryAP(:,:,j);
        R_eeg = xcorr(eeg2,eeg1,8e3,'unbiased');
        R_eeg = detrend(R_eeg(2:end),'constant');
        S12 = fft(R_eeg)/fs^2*L; % units = per spike rather than per Hztime
        Sxy = 2*S12(2:8001).*exp(2*pi*sqrt(-1)*f*7999/fs);
        Sxy = real(Sxy(:));

        % Update d dependence
        iD = interp1(dValues,1:length(dValues),d,'nearest','extrap');
        count_D(iD) = count_D(iD)+1;
        et = Sxy-Pxy_D(:,iD)/count_D(iD);
        Pxy_D(:,iD) = Pxy_D(:,iD) + Sxy;
        SSExy_D(:,iD) = SSExy_D(:,iD) + et.*(Sxy-Pxy_D(:,iD)/count_D(iD));

        % Convergence check
        if(iLargeLoop==loopSize)
            SIG = SSExy_D./(count_D-1);
            SIG(isinf(SIG)) = nan;
            SIGN = nansum(A.^2.*SIG,2);
            STOP = min(0.99-SIGN/N/d_abs^2);

            if(STOP>=0 || N>20e6)
                break;
            end
            iLargeLoop = 0;
        end
        iLargeLoop = iLargeLoop + 1;

        % Resample check
        if(iSampling==samplingInterval)
            neuronIdcs = reshape(sample_blue_neurons(2*samplingInterval),samplingInterval,2);
            dSamples = interp1(B(iUniqueD),dValues(iUniqueD),rand(samplingInterval,1),'next','extrap');
            iSampling = 0;
        end
        iSampling = iSampling+1;
    end


    T = toc/N;
    save(saveFile,'dValues','Pxy_D','SSExy_D','count_D','f','T','STOP');
end
function [x1,x2] = getOrthBasis(x0)
    x0 = x0/norm(x0,2);
    if(prod(size(x0))==3)
        if(x0>0.9)
            x1 = [0,1,0];
        else
            x1 = [1,0,0];
        end
        x1 = x1-x0*sum(x1.*x0);
        x1 = x1/norm(x1,2);
        x2 = cross(x1,x0);
    else
        idcs = x0(:,1)>0.9;
        x1 = zeros(size(x0));
        x1(idcs,2) = 1;
        x1(~idcs,1) = 1;

        x1 = x1-x0.*sum(x1.*x0,2);
        x1 = x1./vecnorm(x1,2,2);
        x2 = cross(x1,x0);
    end
end
