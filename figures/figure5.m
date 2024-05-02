function figure5
% FIGURE5 generates the panels in figure 3 of the manuscript.
% See also import_paralytic_data, fittingmodel.

    baseFolder = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(baseFolder,'auxiliary_functions'));
    addpath(fullfile(baseFolder,'modelling'));


    [full_model,AP_model] = fittingmodel('eq6');
    [f0,S] = import_paralytic_data;
    low_noise = 1e-3;

    figureNB(13.2,5.8);
    axes('Position',[0.08, 0.14, 0.38, 0.8])
        plot(f0,S(:,9:16),'color','k');
        hold on;
        plot([1,3e3],low_noise*[1,1],'color',red,'LineWidth',1)
        Psyn = 10.^AP_model(f0,[3e-3,1e-3,-Inf,4.6]);
        plot(f0,Psyn,'color',blue,'LineStyle','--')
        Psyn = 10.^AP_model(f0,[20e-3,4e-3,-Inf,3.6]);
        plot(f0,Psyn,'color',blue,'LineStyle','--')
        Psyn = 10.^AP_model(f0,[20e-3,4e-3,-Inf,3.6])+10.^AP_model(f0,[3e-3,1e-3,-Inf,4.6]);
        plot(f0,Psyn,'color',blue,'LineWidth',1)
        set(gca,'yscale','log')
        xlabel('Frequency (Hz)')
        ylabel(['PSD (' char(956) 'V^2/Hz)'])
        title('Unparalyzed','FontWeight','normal')
        xlim([0,100]);
        gcaformat;
        set(gca,'FontSize',8)

    axes('Position',[0.58, 0.14, 0.38, 0.8])
        plot(f0,S(:,1:8),'color','k');
        hold on;
        plot([1,3e3],low_noise*[1,1],'color',red,'LineWidth',1)
        Psyn = 10.^AP_model(f0,[3e-3,1e-3,-Inf,3.3]);
        plot(f0,Psyn,'color',blue,'LineStyle','--')
        Psyn = 10.^AP_model(f0,[20e-3,4e-3,-Inf,3.6]);
        plot(f0,Psyn,'color',blue,'LineStyle','--')
        Psyn = 10.^AP_model(f0,[20e-3,4e-3,-Inf,3.6])+10.^AP_model(f0,[3e-3,1e-3,-Inf,3.3]);
        plot(f0,Psyn,'color',blue,'LineWidth',1)
        set(gca,'yscale','log')
        xlabel('Frequency (Hz)')
        ylabel(['PSD (' char(956) 'V^2/Hz)'])
        title('20 mg cisatracurium','FontWeight','normal')
        xlim([0,100]);
        gcaformat;
        set(gca,'FontSize',8)

end