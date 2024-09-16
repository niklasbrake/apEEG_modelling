function synthetic_spikes(comboID)
    addpath('/lustre04/scratch/nbrake/code/natcomms_framework');

    local_cluster = parcluster('local');
    local_cluster.JobStorageLocation = getenv('SLURM_TMPDIR');
    parpool(local_cluster,40);

    main(comboID,0,1);
    main(comboID,1,0);
    main(comboID,1,1);
end
function main(comboID,R_E,R_I)

    folder = ['/lustre04/scratch/nbrake/data/simulations/apEEG/synthetic_spikes/mCombo' int2str(comboID) '_RE' int2str(R_E) '_RI' int2str(R_I)]

    load(fullfile(folder,'model.mat'),'network');
    load(fullfile(network.preNetwork,'chol_var.mat'),'Lc');
    load(fullfile(network.preNetwork,'covariance.mat'),'ei','gam','dt');

    M = ceil(network.tmax/dt);
    xCell = cell(M,2);
    parfor i = 1:M
        ids0 = find((Lc*randn(network.getsynapsecount,1))<gam);
        ts0 = (i*ones(size(ids0))+rand(size(ids0)))*dt;
        xCell(i,:) = {ids0,ts0};
    end
    ids = cat(1,xCell{:,1});
    ts = cat(1,xCell{:,2});

    % Map ids onto synapse ids
    cons = csvread(fullfile(network.postNetwork,'connections.csv'));
    synIDs = cons(:,3);
    ids = synIDs(ids);
    [~,I] = sort(synIDs);
    ei = ei(I);

    network_simulation_beluga.save_presynaptic_network(ids,ts,ei,network.getsynapsecount,network.spikingFile);
end