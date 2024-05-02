function figureS1
% FIGURES1 generates the panels in figure S1 of the manuscript.
% See also compute_uAP_xcorr

    baseFolder = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(baseFolder,'auxiliary_functions'));
    addpath(fullfile(baseFolder,'modelling'));

    [t,R] = compute_uAP_xcorr;

    figureNB(5.5,5.15);
    axes('Position',[0,2/3,1/3,1/3]);
        plotwitherror(t,1./(abs(t)<10).*R(:,:,1,1),'CI','LineWidth',1,'color','k');
        xlabel('Lag (ms)');
        ylabel('Correlation')
        axis off;
        xlim([-12,12]);
        ylim([-1,1]);
    axes('Position',[1/3,2/3,1/3,1/3]);
        plotwitherror(t,1./(abs(t)<10).*R(:,:,1,2),'CI','LineWidth',1,'color','k');
        xlabel('Lag (ms)');
        ylabel('Correlation')
        axis off;
        xlim([-12,12]);
        ylim([-1,1]);
    axes('Position',[2/3,2/3,1/3,1/3]);
        plotwitherror(t,1./(abs(t)<10).*R(:,:,1,3),'CI','LineWidth',1,'color','k');
        xlabel('Lag (ms)');
        ylabel('Correlation')
        axis off;
        xlim([-12,12]);
        ylim([-1,1]);
    axes('Position',[0,1/3,1/3,1/3]);
        plotwitherror(t,1./(abs(t)<10).*R(:,:,2,1),'CI','LineWidth',1,'color','k');
        xlabel('Lag (ms)');
        ylabel('Correlation')
        axis off;
        xlim([-12,12]);
        ylim([-1,1]);
    axes('Position',[1/3,1/3,1/3,1/3]);
        plotwitherror(t,1./(abs(t)<10).*R(:,:,2,2),'CI','LineWidth',1,'color','k');
        xlabel('Lag (ms)');
        ylabel('Correlation')
        axis off;
        xlim([-12,12]);
        ylim([-1,1]);
    axes('Position',[2/3,1/3,1/3,1/3]);
        plotwitherror(t,1./(abs(t)<10).*R(:,:,2,3),'CI','LineWidth',1,'color','k');
        xlabel('Lag (ms)');
        ylabel('Correlation')
        axis off;
        xlim([-12,12]);
        ylim([-1,1]);
    axes('Position',[0,0,1/3,1/3]);
        plotwitherror(t,1./(abs(t)<10).*R(:,:,3,1),'CI','LineWidth',1,'color','k');
        xlabel('Lag (ms)');
        ylabel('Correlation')
        axis off;
        xlim([-12,12]);
        ylim([-1,1]);
    axes('Position',[1/3,0,1/3,1/3]);
        plotwitherror(t,1./(abs(t)<10).*R(:,:,3,2),'CI','LineWidth',1,'color','k');
        xlabel('Lag (ms)');
        ylabel('Correlation')
        axis off;
        xlim([-12,12]);
        ylim([-1,1]);
    axes('Position',[2/3,0,1/3,1/3]);
        plotwitherror(t,1./(abs(t)<10).*R(:,:,3,3),'CI','LineWidth',1,'color','k');
        xlabel('Lag (ms)');
        ylabel('Correlation')
        axis off;
        xlim([-12,12]);
        ylim([-1,1]);
end