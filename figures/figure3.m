function figure3
% FIGURE3 generates the panels in figure 3 of the manuscript.
% Data files required:
%   mtype_abundance.mat
%   dendrite_asymmetry.mat
%   unitarySpectrum.mat
% See also calculate_dendrite_asymmetry, compute_AP_spectra, modelling/unitary_AP_response/simulations/main.sh.

    baseFolder = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(baseFolder,'auxiliary_functions'));
    addpath(fullfile(baseFolder,'data_analysis'));

    load(fullfile(baseFolder,'data_files','mtype_abundance.mat'),'mtype_abundance');
    load(fullfile(baseFolder,'data_files','unitaryAP.mat'),'mtype','ei_type');
    load(fullfile(baseFolder,'data_files','asymmetry_indices.mat'),'asym_idx');
    load(fullfile(baseFolder,'data_files','unitarySpectrum.mat'),'psd');

    [f,Sxx] = compute_AP_spectra;
    [f0,S] = import_Scheer2006;
    low_noise = (8e-3).^2;
    N2 = 16e9;

    pEst = sum(psd)*mean(diff(f));
    C = mtype_abundance(mtype,:).Abundance;
    figureNB(13.2,4.5);
    axes('Position',[0.1, 0.22, 0.23, 0.68])
        sc = scatter(asym_idx(~ei_type),(1e6)^2*pEst(~ei_type),1+C(~ei_type),[0.5,0.5,0.8],'filled','MarkerEdgeColor',[0.5,0.5,0.8]/2);
        sc.AlphaData = 0.6+0.4*min(C(~ei_type)/2,1);
        sc.MarkerFaceAlpha = 'flat'; sc.MarkerEdgeAlpha = 'flat';
        hold on;
        sc = scatter(asym_idx(ei_type),(1e6)^2*pEst(ei_type),1+C(ei_type),[0.8,0.5,0.5],'filled','MarkerEdgeColor',[0.8,0.5,0.5]/2);
        sc.AlphaData = 0.6+0.4*min(C(ei_type)/2,1);
        sc.MarkerFaceAlpha = 'flat'; sc.MarkerEdgeAlpha = 'flat';
        xlim([100,10000]);

        t = log10(get(gca,'xlim'));
        t = linspace(t(1),t(end),1e3);
        FT = fitlm(asym_idx(:),(1e6)^2*(pEst(:)'),'intercept',false,'RobustOpts',true);
        fprintf('R^2 = %f\n',FT.Rsquared.Ordinary);
        plot(10.^t,FT.predict(10.^t(:)),'-k')

        xlabel('Dendrite asymmetry index')
        ylabel(['Unitary AP power (pV^2)'])
        xticks([100,1000,10000]);
        xticklabels([100,1000,10000])
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        gcaformat(gca,true,8);
    axes('Position',[0.425, 0.22, 0.23, 0.68])
        plot(f,Sxx*(1e9)^2,'color',blue,'Linewidth',1);
        xlim([10,1e3])
        set(gca,'xscale','log')
        xlabel('Frequency (Hz)')
        ylabel(['PSD (fV^2/Hz)'])
        xlim([1,1e3])
        xticks([1,10,100,1000]);
        xticklabels([1,10,100,1000]);
        gcaformat;
    axes('Position',[0.75, 0.22, 0.23, 0.68])
        plot(f0,S,'color','k','Linewidth',1);
        hold on;
        plot([1,3e3],low_noise*[1,1],'color',red,'LineWidth',1.5,'LineStyle','-')
        plot(f,0.1*N2*Sxx,'color',blue,'LineWidth',1)
        plot(f,1*N2*Sxx,'color',blue,'LineWidth',1)
        plot(f,10*N2*Sxx,'color',blue,'LineWidth',1)
        plot(f,100*N2*Sxx,'color',blue,'LineWidth',1)
        set(gca,'xscale','log')
        set(gca,'yscale','log');
        xlabel('Frequency (Hz)');
        ylabel(['PSD (' char(956) 'V^2/Hz)'])
        ylim([1e-8,1e1])
        yticks([1e-8,1e-4,1e0])
        xlim([1,3e3])
        xticks([1,10,100,1000]);
        xticklabels([1,10,100,1000]);
        gcaformat(gcf,true,8)

    [J,ID] = findgroups(mtype);

    for i = 1:length(ID)
        abundance(i) = mtype_abundance(ID{i},:).Abundance;
    end

    figureNB(10.4,5);
    axes('Position',[0.14, 0.18, 0.8, 0.61])
        i = find(~ei_type);
        S = scatter(J(i),(1e6)^2*pEst(i),8,[0.5,0.5,0.8],'filled','MarkerEdgeColor',[0.5,0.5,0.8]/2);
        hold on;
        i = find(ei_type);
        S = scatter(J(i),(1e6)^2*pEst(i),10,[0.8,0.5,0.5],'filled','MarkerEdgeColor',[0.8,0.5,0.5]/2);
        set(gca,'yscale','log')
        ylabel(['Unitary AP power (pV^2)'])
        gcaformat;
        xax = get(gca,'xaxis');
        xax.MinorTickValues = 1:55;
        xax.MinorTick = 'on';
        xlim([0,56])
        xlabel('Neuron morphology type index')
        set(gca,'FontSize',8)
    axes('Position',[0.14, 0.92, 0.8, 0.05])
        [~,I] = sort(abundance,'descend');
        for i = I
            scatter(i,1,2*max(2,abundance(i)),abundance(i),'filled');
            hold on;
        end
        xlim([0,56])
        axis off
        colormap(flip(copper))
    axes('Position',[0.145 0.83 0.28 0.06])
        exA = [2,5,10,15];
        scatter([0,0.36,0.65,1],zeros(1,4),2*exA,exA,'filled');
        text(-0.1,1,'Abundance','FontSize',8,'HorizontalAlignment','left');
        text(0.05,0,'<2%','FontSize',8,'VerticalAlignment','middle')
        text(0.41,0,'5%','FontSize',8,'VerticalAlignment','middle')
        text(0.7,0,'10%','FontSize',8,'VerticalAlignment','middle')
        text(1.06,0,'15%','FontSize',8,'VerticalAlignment','middle')
        ylim([0,1]);
        xlim([-0.1,1]);
        axis off;

    axes('Position',[0.88 0.83 0.02 0.06])
        S = scatter(0,1,12,[0.8,0.5,0.5],'filled','MarkerEdgeColor',[0.8,0.5,0.5]/2);
        hold on;
        S = scatter(0,0,12,[0.5,0.5,0.8],'filled','MarkerEdgeColor',[0.5,0.5,0.8]/2);
        ylim([0,1]);
        text(1,1,'Ex.','FontSize',8,'VerticalAlignment','middle')
        text(1,0,'In.','FontSize',8,'VerticalAlignment','middle')
        axis off;
end