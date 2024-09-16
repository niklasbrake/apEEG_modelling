function project_covariance_matrix(comboID)
    addpath('/lustre04/scratch/nbrake/code/natcomms_framework');

    main(comboID,0,1);
    main(comboID,1,0);
    main(comboID,1,1);
end
function main(comboID,R_E,R_I)

    folder = ['/lustre04/scratch/nbrake/data/simulations/apEEG/synthetic_spikes/mCombo' int2str(comboID) '_RE' int2str(R_E) '_RI' int2str(R_I)]

    load(fullfile(folder,'presynaptic_network','covariance.mat'),'L0');

    Lc = chol(nearcorr(L0),'lower');

    save(fullfile(folder,'presynaptic_network','chol_var.mat'),'Lc','-v7.3')

end
