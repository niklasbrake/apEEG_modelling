function calculate_covariance_matrix(comboID)
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


    mCombos = zeros(11*10/2,2);
    count = 1;
    for i = 1:11
        for j = i+1:11
            mCombos(count,:) = [i,j];
            count = count+1;
        end
    end

    % Initialize network
    network = network_simulation_beluga(folder);
    network = network.initialize_postsynaptic_network(2,mCombos(comboID,:));
    N = network.getsynapsecount;
    network.tmax = 10e3;

    % Distribute synapses across dendrites and get locations
    network = network.form_connections(0);
    network.save();
    cons = csvread(fullfile(network.postNetwork,'connections.csv'));
    synIDs = cons(:,3);
    mData = network.getmData;
    P = arrayfun(@(i,j) mData{i}.pos(j,:),cons(:,1),cons(:,2),'UniformOutput',false);
    P = cat(1,P{:});
    [thet,phi] = cart2sph(P(:,1),P(:,2),P(:,3));
    X = [phi(:),thet(:)];

    % Define excitatory/inhibitory cells and avg. firing rates
    dt = 1;
    ei = rand(N,1)>network.parameters.eiFraction;
    r = (network.parameters.eCellParams.firingRate*(1-ei) + network.parameters.iCellParams.firingRate*ei)*dt/1e3;
    rVar = r.*(1-r);
    gam = icdf('normal',r,0,1);

    % Max correlation
    RE = 0.0251*R_E;
    RI = 0.0251*R_I;

    % Use haversine distance among synapses to construct covariance matrix
    KxxCell = cell(N,3); % Store sparse cov matrix (i,j,v) for every synapse
    hav = @(x) (1-cos(x))/2;
    hav_d = @(p1,x) hav(p1(1)-x(:,1))+(1-hav(x(:,1)-p1(1))-hav(p1(1)+x(:,1))).*hav(p1(2)-x(:,2));
    parfor i = 1:N
        d = hav_d(X(i,:),X(i+1:end,:));
        idcs = find(d<0.25);
        if(RE==RI)
            R = RE;
        else
            ei0 = ei(i);
            ei1 = ei(i+idcs);
            if(ei0==0) % First synapse is excitatory
                R = RE*(ei1==ei0);
            else % First synapse is inhibitory
                R = RI*(ei1==ei0);
            end
        end
        KxxCell(i,:) = {idcs+i,i+0*zeros(size(idcs)),R.*exp(-10*d(idcs)).*sqrt(rVar(i)*rVar(i+idcs))};
    end
    i0 = cat(1,KxxCell{:,1});
    j0 = cat(1,KxxCell{:,2});
    S = cat(1,KxxCell{:,3});

    % Find mu and sig of n-d Gaussian with which to sample spikes
    fun = @(x,kk,rr,gg) kk+rr-mvncdf(gg,0,[1,x;x,1]);
    A = zeros(length(i0),1);
    ibnds = [-1+1e-4,1-1e-4];
    parfor k = 1:length(i0)
        i = i0(k);
        j = j0(k);
        A(k) = fzero(@(x) S(k)+r(i)*r(j)-mvncdf(gam([i,j]),0,[1,x;x,1]),ibnds);
    end
    L = sparse(i0,j0,A,N,N);
    L = L+L'+eye(N);
    L0 = nearcorr(L);
    save(fullfile(network.preNetwork,'covariance.mat'),'L0','ei','gam','-v7.3');
end