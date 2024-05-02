function figureS2
% FIGURES2 generates the panels in figure S2 of the manuscript
%.
% Data files required:
%   unitarySpectrum.mat
%   neuron_models/L4_BP_cACint209_1/matlab_recordings/*.mat
%
% See also modelling/unitary_AP_response/simulations/main.sh

    baseFolder = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(baseFolder,'auxiliary_functions'));
    addpath(fullfile(baseFolder,'modelling'));

    load(fullfile(baseFolder,'data_files','unitarySpectrum.mat'),'freq','psd','mtype');

    folder = fullfile(baseFolder,'data_files','neuron_models','L4_BP_cACint209_1','matlab_recordings');
    data1=load(fullfile(folder,'synaptic_input_EI30.mat'));
    data2=load(fullfile(folder,'synaptic_input_EI30_passive.mat'));
    [y,x] = findpeaks(data1.voltage,data1.time,'MinPeakHeight',0);

    [G,I] = findgroups(mtype);

    B = mean(psd(1:2,:))./mean(psd(5:7,:));

    X = [log10(freq(freq<=10)),ones(size(freq(freq<=10)))];
    B = [];
    for i = 1:size(psd,2)
        y = log10(psd(freq<=10,i));
        B(i,:) = X\y;
    end

    i0 = 286;
    figureNB(18,10);
    axes('Position',[0.07, 0.67, 0.90, 0.30])
        plot(G,B,'.k')
        hold on;
        plot(G(i0),B(i0,1),'ok','LineWidth',1,'color',blue);
        xlim([0.5,55.5])
        xticks(unique(G))
        xticklabels(strrep(I,'_','\_'))
        ylabel('Spectral exponent (1-10 Hz)')
        ylim([-3.5,2])
    gcaformat(gcf,true,8)
        xax = get(gca,'xaxis');
        xax.FontSize = 7;
        set(gca,'XTickLabelRotation',45)
    axes('Position',[0.07, 0.1, 0.23, 0.40])
        plot(freq,psd(:,i0),'k','LineWidth',1);
        set(gca,'xscale','log');
        set(gca,'yscale','log');
        xlabel('Frequency (Hz)');
        ylabel(['Unit apEEG PSD (' char(956) 'V^2/Hz)'])
        hold on;
        plot([1,100],10.^(B(i0,2)+B(i0,1)*[0,2]),'--k')
        gcaformat(gca,true,8)

    axes('Position',[0.42, 0.34, 0.55, 0.19])
        plot(data1.time(2:end),data1.dipoles(2:end,2,1),'color','k')
        hold on;
        plot(data1.time(2:end),data2.dipoles(2:end,2,1),'color',[0.6,0.6,0.6])
        ylabel(['q_z (nA' char(956) 'm)'])
        xlim([0,1e4]);
        set(get(gca,'xaxis'),'visible','off')
        gcaformat(gca,true,8)

    axes('Position',[0.42, 0.10, 0.55, 0.19])
        plot(data1.time(2:end),data1.dipoles(2:end,3,1)-data2.dipoles(2:end,3,1),'k')
        line([0,1e4],[0,0],'color',[0.6,0.6,0.6],'LineWidth',1)
        ylim([-0.5,0.5])
        xlim([0,1e4]);
        ylabel('q_z (active - passive)')
        xlabel('Time (ms)')
        gcaformat(gca,true,8)
end