function compute_uAP
% COMPUTE_UAP computes and saves the unitary AP response of all the neuron
%   simulations. Called by compute_uAP.sh.
    warning('off','signal:findpeaks:largeMinPeakHeight');

    % Directory containg all the neuron models. Equivalent to
    % apEEG_modelling/data_files/neuron_models
    folder = '/lustre04/scratch/nbrake/data/simulations/unitary_AP';

    % Get metadata about each simulation
    [mtype,ei_type,layer,morph] = getMtypes(folder);

    % Compute the unitary AP response
    [savedUnitaryAP,N,EI,files,CV,dpY] = getUnitaryAP(folder);

    save('/lustre04/scratch/nbrake/data/simulation_analyzed/unitaryAP/unitaryAP.mat');

end
function [mtype,ei_type,layer,morph] = getMtypes(folder)
% [mtype,ei_type,layer,morph] = GETMTYPES(folder) returns metadata about
%   each model file in the directory 'folder', include the cell class (mtype)
%   excitatory/inhibitory (ei_type), the layer ID (layer), and the morphology ID.
    F = dir(folder);
    F = F(3:end);
    M = length(F);

    for i = 1:M
        items = split(F(i).name,'_');
        layer{i} = items{1};
        morph{i} = items{2};
        if(strcmp(items{3},'L1') || strcmp(items{3},'L4'))
            morph{i} = [items{2} '_' items{3}];
        end
        mtype{i} = [layer{i} '_' morph{i}];
    end
    e_mtypes = {'PC','SS','SP','STPC','UTPC','TTPC1','TTPC2','TPC','BPC','IPC','TPC_L1','TPC_L4'};
    ei_type = ismember(morph,e_mtypes);
end

function [savedUnitaryAP,N,EI,files,CV,dpY] = getUnitaryAP(folder)
% [savedUnitaryAP,N,EI,files,CV,dpY] = getUnitaryAP(folder) computes the unitary
%   AP response for all the simulations in the directory 'folder'. savedUnitaryAP
%   is an array containg the dipole vector of the unitary AP response of every
%   simulation, N is the number of spikes in each simulation, EI is the EI ratio
%   of each simulation, files is a cell array of all the model names, CV is the
%   coefficient of variation of the spike times, and dpY is the raw dipole in the Y axis.

    F = dir(folder);
    F = F(3:end);
    M = length(F);

    savedUnitaryAP = nan(8001,3,M);
    N = nan(M,1);

    fs = 16e3; % Hz

    files = cell(M,1);
    for i = 1:M
        files{i} = F(i).name;
        folder0 = fullfile(folder,F(i).name,'matlab_recordings');

        % Load EI ratio
        fid = fopen(fullfile(folder,F(i).name,'EI_ratio.csv'));
        ei = textscan(fid,'%s');
        fclose(fid);
        ei = ei{1}{1};
        if(mod(str2num(ei),1)==0 && str2num(ei)<10)
            ei = ['0' ei];
        end

        load(fullfile(folder0,sprintf('synaptic_input_EI%s_passive.mat',ei)));
        dipoles_passive = dipoles;
        load(fullfile(folder0,sprintf('synaptic_input_EI%s.mat',ei)));
        dipoles = dipoles-dipoles_passive;

        [y,x] = findpeaks(voltage,'MinPeakHeight',0);
        N(i) = length(x);
        EI(i) = str2num(ei);
        unitaryAP = zeros(8001,3,length(x));
        for j = 1:length(x)
            idcs = max(min(x(j)-2e3:x(j)+6e3,length(voltage)),1);
            y = dipoles(idcs,:);
            y(idcs==1,:) = 0;
            y(idcs==length(voltage),:) = 0;
            unitaryAP(:,:,j) = y;
        end
        unitaryAP = unitaryAP-nanmedian(dipoles);
        savedUnitaryAP(:,:,i) = nanmedian(unitaryAP,3);

        CV(i) = std(diff(x))/mean(diff(x));
        dpY{i} = dipoles(:,2);
    end
end

function [V,files] = passive_power(folder)
% PASSIVE_POWER compute the apEEG power of each model during the passive simulation.
    load('/lustre04/scratch/nbrake/data/anatomy_nyhead_model.mat');
    M = size(sa.cortex75K.vc,1);
    Lxyz = zeros(M,3);
    for idx = 1:M
        L0 = squeeze(sa.cortex75K.V_fem(49,idx,:))'; % czIDX = 49
        vz = sa.cortex75K.normals(idx,:);
        [vx,vy] = getOrthBasis(vz);
        A = [vx(:),vy(:),vz(:)];
        Lxyz(idx,:) = 1e-6*L0*A;
    end
    Lxyz = Lxyz(sa.cortex10K.in_from_cortex75K,:);


    F = dir(folder);
    F = F(3:end);
    M = length(F);

    fs = 16e3; % Hz

    files = cell(M,1);
    dp = nan(16000,3,M);
    V = zeros(size(Lxyz,1),M);

    for i = 1:M
        files{i} = F(i).name;
        folder0 = fullfile(folder,F(i).name,'matlab_recordings');

        fid = fopen(fullfile(folder,F(i).name,'EI_ratio.csv'));
        ei = textscan(fid,'%s');
        fclose(fid);
        ei = ei{1}{1};
        if(mod(str2num(ei),1)==0 && str2num(ei)<10)
            ei = ['0' ei];
        end

        load(fullfile(folder0,sprintf('synaptic_input_EI%s_passive.mat',ei)));
        dp(:,:,i) = dipoles(16001:32000,:);
    end

    for j = 1:size(Lxyz,1)
        V(j,:) = squeeze(nanvar(nansum(Lxyz(j,:).*dp,2)));
    end
end