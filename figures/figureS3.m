function figureS3
% FIGURES3 generates the panels in figure S3 of the manuscript.
%
% Data files required:
%   pairwise_distance.mat
%
% See also compute_cortical_area

    baseFolder = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(baseFolder,'auxiliary_functions'));

    load(fullfile(baseFolder,'data_files','pairwise_distance.mat'),'rValues','total_area');

    figureNB(13.2,5.2);
    subplot(1,2,1);
        plotwitherror(rValues,total_area/100,'SE','color','k','LineWidth',1,'LineStyle','-','Marker','.','MarkerSize',10)
        xlim([0,200])
        xlabel('Ball radius (mm)')
        ylabel('Enclosed cortical SA (cm^2)')
    subplot(1,2,2);
        plot(0.5*(rValues(1:end-1)+rValues(2:end)),16e9^2*diff(mean(total_area,2)/max(total_area(:))),'.-k');
        set(gca,'xscale','log');
        set(gca,'yscale','log');
        ylabel('Estimated number of neuron pairs')
        yyaxis right;
        plot(rValues,exp(-rValues.^2/6),'color','r','LineWidth',1)
        xlabel('Pariwise distance (mm)');
        ylabel('Spike synchrony');
        ylim([0,1.01]);
        yticks([0,1]);
        yticklabels({0,'R_{max}'});
        xticks([1e-2,1e-1,1,1e1,1e2])
        xticklabels([1e-2,1e-1,1,1e1,1e2])
    gcaformat(gcf,true,8)
end