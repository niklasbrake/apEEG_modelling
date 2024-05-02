function figure2
% FIGURE2 generates the panels in figure 2 of the manuscript.
%
% Data files required:
%   AP_scaling.mat
%
% See also compute_scaling_with_firing_frequency

    baseFolder = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(baseFolder,'auxiliary_functions'));

    load(fullfile(baseFolder,'data_files','AP_scaling.mat'));

    A = regress(R(:),[firingFrequency(:),ones(size(firingFrequency(:)))])
    t = linspace(0.1,100,20);

    figureNB(9,4.4);
    subplot(1,2,1);
        plot(firingFrequency(:),B2(:),'.k')
        line([0.1,100],[0.1,100],'color',red,'LineWidth',1)
        xlabel('Firing rate (Hz)')
        ylabel('Scaling factor (\beta)')
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        xlim([1,100])
        xticks([0.1,1,10,100])
        xticklabels([0.1,1,10,100])
        gcaformat;
    subplot(1,2,2);
        plot(firingFrequency(:),R(:),'.k')
        hold on;
        plot(t,t*A(1)+A(2),'color',red,'LineWidth',2);
        ylim([0,1])
        xlabel('Firing rate (Hz)')
        ylabel('R^2')
        set(gca,'xscale','log')
        xlim([1,100])
        xticks([0.1,1,10,100])
        xticklabels([0.1,1,10,100])
    gcaformat(gcf,true,8);

    edges = [0,20,40,80,Inf]
    bins = discretize(firingFrequency(:),edges);
    split_Y = splitapply(@(x)mean(x,2),psd_Y,bins');
    split_Yhat = splitapply(@(x)mean(x,2),psd_Yhat,bins');

    figureNB(5,4);
    axes('Position',[0.16, 0.62, 0.34, 0.32])
        plot(f0,split_Y(:,1),'color','k','LineWidth',2);
        hold on;
        plot(f0,split_Yhat(:,1),'color',blue,'LineWidth',1);
        xlim([1,3e3]);
        ylim([1e-4,1]);
        xticks([1,10,100,1000]);
        xticklabels({});
        set(gca,'xscale','log')
        set(gca,'yscale','log')

    axes('Position',[0.6, 0.62, 0.34, 0.32])
        plot(f0,split_Y(:,2),'color','k','LineWidth',2);
        hold on;
        plot(f0,split_Yhat(:,2),'color',blue,'LineWidth',1);
        xlim([1,3e3]);
        ylim([1e-4,1]);
        xticks([1,10,100,1000]);
        xticklabels({});
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        yticklabels({});

    axes('Position',[0.16, 0.21, 0.34, 0.32])
        plot(f0,split_Y(:,3),'color','k','LineWidth',2);
        hold on;
        plot(f0,split_Yhat(:,3),'color',blue,'LineWidth',1);
        xlim([1,3e3]);
        ylim([1e-4,1]);
        xticks([1,10,100,1000]);
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        xlabel('Frequency (Hz)');

    axes('Position',[0.6, 0.21, 0.34, 0.32])
        plot(f0,split_Y(:,4),'color','k','LineWidth',2);
        hold on;
        plot(f0,split_Yhat(:,4),'color',blue,'LineWidth',1);
        xlim([1,3e3]);
        ylim([1e-4,1]);
        xticks([1,10,100,1000]);
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        xlabel('Frequency (Hz)');
        yticklabels({});

    gcaformat(gcf,true,8);
end