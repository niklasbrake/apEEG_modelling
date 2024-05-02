function figure7
% FIGURE7 generates the panels in figure 7 of the manuscript.
% See also compute_AP_spectra, fittingmodel.

    baseFolder = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(baseFolder,'auxiliary_functions'));
    addpath(fullfile(baseFolder,'data_analysis'));
    addpath(fullfile(baseFolder,'modelling','MC simulations'));


    [f,Rxx,Rxy] = compute_AP_spectra;
    N2 = 16e9;


    A = exp(-(f-10).^2/2^2);
    Gsyn = exp(-(f-50).^2/10^2);
    G = N2*Rxx + 0.1*N2*(N2-1)*exp(-2*(pi*(f(:)-130)*10e-3).^2).*Rxy;
    additive_noise = 1e-2;

    [full_model,AP_model] = fittingmodel('eq6');
    Psyn = 10.^AP_model(f,[20e-3,4e-3,-Inf,3.6])+10.^AP_model(f,[3e-3,1e-3,-Inf,3]);
    Psyn2 = 10.^AP_model(f,[35e-3,4e-3,-Inf,4])+0.1*10.^AP_model(f,[3e-3,1e-3,-Inf,3]);

    Y1 = (1+3*A+Gsyn+1./f.^2).*Psyn+0*G+additive_noise;
    Y2 = (1+3*A+Gsyn+1./f.^2).*Psyn2+0*G+0.5*additive_noise;
    figureNB(9,6);
    axes('Position',[0.12, 0.65, 0.34, 0.31])
        plot(f,Psyn,'color',blue,'LineStyle','-');
        hold on;
        plot(f,additive_noise+f*0,'color',red,'LineStyle','-');
        yticks([1e-3,1e-1,1e1])
        ylabel(['PSD (' char(956) 'V^2/Hz)'])
        plot(f,Y1,'color','k','LineWidth',1);
        ylim([1e-3,1e1]);
        set(gca,'yscale','log')
        xlim([0,200])
        xlabel('Frequency (Hz)');

    axes('Position',[0.59, 0.65, 0.34, 0.31])
        plot(f,Y1,'color',[0.6,0.6,0.6],'LineWidth',1);
        hold on;
        plot(f,Psyn2,'color',blue,'LineStyle','-');
        plot(f,0.5*additive_noise+f*0,'color',red,'LineStyle','-');
        yticks([1e-3,1e-1,1e1])
        ylabel(['PSD (' char(956) 'V^2/Hz)'])
        plot(f,Y2,'k','LineWidth',1);
        ylim([1e-3,1e1]);
        set(gca,'yscale','log')
        xlim([0,200])
        xlabel('Frequency (Hz)');

    Z1 = Y1;
    Z2 = Y2;
    axes('Position',[0.12, 0.13, 0.2, 0.24])
        fill([f;f(end);.1],[10*log10(Z2./Z1);0;0],'k','LineWidth',1,'FaceAlpha',0.4);
        ylabel('Power (dB)');
        xlim([1,200])
        xlabel('Frequency (Hz)');
        set(gca,'xscale','log');
        xticks([1,10,100])
        xticklabels([1,10,100])
        title({'Uncorrected'},'FontWeight','normal')

    Z1 = Y1./(Psyn+additive_noise);
    Z2 = Y2./(Psyn2+0.5*additive_noise);
    axes('Position',[0.425, 0.13, 0.2, 0.24])
        fill([f;f(end);.1],[10*log10(Z2./Z1);0;0],'k','LineWidth',1,'FaceAlpha',0.4);
        ylabel('Power (dB)');
        xlim([1,200])
        xlabel('Frequency (Hz)');
        set(gca,'xscale','log');
        xticks([1,10,100])
        xticklabels([1,10,100])
        title({'Detrended','(divisive)'},'FontWeight','normal')

    Z1 = (Y1-additive_noise)./Psyn;
    Z2 = (Y2-0.5*additive_noise)./Psyn2;
    axes('Position',[0.73, 0.13, 0.2, 0.24])
        fill([f;f(end);.1],[10*log10(Z2./Z1);0;0],'k','LineWidth',1,'FaceAlpha',0.4);
        ylabel('Power (dB)');
        xlim([1,200])
        xlabel('Frequency (Hz)');
        set(gca,'xscale','log');
        xticks([1,10,100])
        xticklabels([1,10,100])
        title({'Detrended','(mixed)'},'FontWeight','normal')
        ylim([-1,1])

    gcaformat(gcf,true,8)


    Y1 = Y1+2*G;
    Y2 = Y2+2*G;
    figureNB(9,6);
    axes('Position',[0.12, 0.65, 0.34, 0.31])
        plot(f,Psyn,'color',blue,'LineStyle','-');
        hold on;
        plot(f,additive_noise+f*0,'color',red,'LineStyle','-');
        yticks([1e-3,1e-1,1e1])
        ylabel(['PSD (' char(956) 'V^2/Hz)'])
        plot(f,Y1,'color','k','LineWidth',1);
        ylim([1e-3,1e1]);
        set(gca,'yscale','log')
        xlim([0,200])
        xlabel('Frequency (Hz)');

    axes('Position',[0.59, 0.65, 0.34, 0.31])
        plot(f,Y1,'color',[0.6,0.6,0.6],'LineWidth',1);
        hold on;
        plot(f,Psyn2,'color',blue,'LineStyle','-');
        plot(f,0.5*additive_noise+f*0,'color',red,'LineStyle','-');
        yticks([1e-3,1e-1,1e1])
        ylabel(['PSD (' char(956) 'V^2/Hz)'])
        plot(f,Y2,'k','LineWidth',1);
        ylim([1e-3,1e1]);
        set(gca,'yscale','log')
        xlim([0,200])
        xlabel('Frequency (Hz)');

    Z1 = Y1;
    Z2 = Y2;
    axes('Position',[0.12, 0.13, 0.2, 0.24])
        fill([f;f(end);.1],[10*log10(Z2./Z1);0;0],'k','LineWidth',1,'FaceAlpha',0.4);
        ylabel('Power (dB)');
        xlim([1,200])
        xlabel('Frequency (Hz)');
        set(gca,'xscale','log');
        xticks([1,10,100])
        xticklabels([1,10,100])
        title({'Uncorrected'},'FontWeight','normal')

    Z1 = Y1./(Psyn+additive_noise);
    Z2 = Y2./(Psyn2+0.5*additive_noise);
    axes('Position',[0.425, 0.13, 0.2, 0.24])
        fill([f;f(end);.1],[10*log10(Z2./Z1);0;0],'k','LineWidth',1,'FaceAlpha',0.4);
        ylabel('Power (dB)');
        xlim([1,200])
        xlabel('Frequency (Hz)');
        set(gca,'xscale','log');
        xticks([1,10,100])
        xticklabels([1,10,100])
        title({'Detrended','(divisive)'},'FontWeight','normal')

    Z1 = (Y1-additive_noise)./Psyn;
    Z2 = (Y2-0.5*additive_noise)./Psyn2;
    axes('Position',[0.73, 0.13, 0.2, 0.24])
        fill([f;f(end);.1],[10*log10(Z2./Z1);0;0],'k','LineWidth',1,'FaceAlpha',0.4);
        ylabel('Power (dB)');
        xlim([1,200])
        xlabel('Frequency (Hz)');
        set(gca,'xscale','log');
        xticks([1,10,100])
        xticklabels([1,10,100])
        title({'Detrended','(mixed)'},'FontWeight','normal')

    gcaformat(gcf,true,8)
end