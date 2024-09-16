function shuffle_synapses(comboID)
    addpath('/lustre04/scratch/nbrake/code/natcomms_framework');
    main(comboID,0,1);
    main(comboID,1,0);
    main(comboID,1,1);
end
function main(comboID,R_E,R_I)

    mCombos = zeros(11*10/2,2);
    count = 1;
    for i = 1:11
        for j = i+1:11
            mCombos(count,:) = [i,j];
            count = count+1;
        end
    end
    folder = ['/lustre04/scratch/nbrake/data/simulations/apEEG/synthetic_spikes/mCombo' int2str(comboID) '_RE' int2str(R_E) '_RI' int2str(R_I)];

    % Import all data
    load(fullfile(folder,'model.mat'));

    % Get synapse spikes times
    cons = csvread(fullfile(network.postNetwork,'connections.csv'));
    [ids,ts,ei] = network.getprenetwork(network.spikingFile);
    network.spikingFile = fullfile(network.preNetwork,'spikeTimes_correlated.csv');
    network_simulation_beluga.save_presynaptic_network(ids,ts,ei,network.getsynapsecount,network.spikingFile);

    % Time shuffle synapse times
    ts_timeshuffle = rand(size(ts))*network.tmax;
    ids_timeshuffle = ids;
    network.spikingFile = fullfile(network.preNetwork,'spikeTimes_timeshuffle.csv');
    network_simulation_beluga.save_presynaptic_network(ids_timeshuffle,ts_timeshuffle,ei,network.getsynapsecount,network.spikingFile);

    % Space shuffle synapse times
    ts_spaceshuffle = ts;
    idcsI = find(ei(ids));
    idcsE = find(~ei(ids));
    ids_spaceshuffle = ids;
    ids_spaceshuffle(idcsI) = ids(idcsI(randperm(length(idcsI))));
    ids_spaceshuffle(idcsE) = ids(idcsE(randperm(length(idcsE))));
    network.spikingFile = fullfile(network.preNetwork,'spikeTimes_spaceshuffle.csv');
    network_simulation_beluga.save_presynaptic_network(ids_spaceshuffle,ts_spaceshuffle,ei,network.getsynapsecount,network.spikingFile);

    % Simulate with each of the spike inputs
    network.spikingFile = fullfile(network.preNetwork,'spikeTimes_correlated.csv');
    network.savePath = fullfile(network.outputPath,'simulation_correlated');
    network.simulate();

    network.spikingFile = fullfile(network.preNetwork,'spikeTimes_timeshuffle.csv');
    network.savePath = fullfile(network.outputPath,'simulation_spaceshuffle');
    network.simulate();

    network.spikingFile = fullfile(network.preNetwork,'spikeTimes_spaceshuffle.csv');
    network.savePath = fullfile(network.outputPath,'simulation_spaceshuffle');
    network.simulate();

end