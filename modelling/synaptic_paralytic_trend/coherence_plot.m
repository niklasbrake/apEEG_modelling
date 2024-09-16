[sa,X] = network_simulation_beluga.getHeadModel;
folder = 'C:\Users\brake\Desktop\synthetic_spikes';
F = dir(folder); F = F(3:end);
fs = 16e3;

for i = 1:length(F)

    load(fullfile(folder,F(i).name,'spikeTimes_E_shuffle','simulation_data.mat'));
    eeg = detrend(network_simulation_beluga.getEEG(dipoles,sa,1e4),'constant');
    [Pa_I(:,i),freq] = pmtm(eeg(:,1),2,[],fs);
    [Pb_I(:,i),freq] = pmtm(eeg(:,2),2,[],fs);
    [P2_I(:,i),freq] = pmtm(eeg(:,1)+eeg(:,2),2,[],fs);
    rho_I(i) = corr(eeg(:,1),eeg(:,2));


    load(fullfile(folder,F(i).name,'spikeTimes_I_shuffle','simulation_data.mat'));
    eeg = detrend(network_simulation_beluga.getEEG(dipoles,sa,1e4),'constant');
    [Pa_E(:,i),freq] = pmtm(eeg(:,1),2,[],fs);
    [Pb_E(:,i),freq] = pmtm(eeg(:,2),2,[],fs);
    [P2_E(:,i),freq] = pmtm(eeg(:,1)+eeg(:,2),2,[],fs);
    rho_E(i) = corr(eeg(:,1),eeg(:,2));

    load(fullfile(folder,F(i).name,'spikeTimes_E_I_shuffle','simulation_data.mat'));
    eeg = detrend(network_simulation_beluga.getEEG(dipoles,sa,1e4),'constant');
    [Pa_EI(:,i),freq] = pmtm(eeg(:,1),2,[],fs);
    [Pb_EI(:,i),freq] = pmtm(eeg(:,2),2,[],fs);
    [P2_EI(:,i),freq] = pmtm(eeg(:,1)+eeg(:,2),2,[],fs);
    rho_EI(i) = corr(eeg(:,1),eeg(:,2));

    load(fullfile(folder,F(i).name,'spikeTimes_E_I_correlated','simulation_data.mat'));
    eeg = detrend(network_simulation_beluga.getEEG(dipoles,sa,1e4),'constant');
    [Pa_0(:,i),freq] = pmtm(eeg(:,1),2,[],fs);
    [Pb_0(:,i),freq] = pmtm(eeg(:,2),2,[],fs);
    [P2_0(:,i),freq] = pmtm(eeg(:,1)+eeg(:,2),2,[],fs);
    rho_0(i) = corr(eeg(:,1),eeg(:,2));

end

figureNB;
subplot(1,4,1)
    plot(freq,smooth(mean(Pa_0+Pb_0,2),20),'color','k');
    hold on;
    plot(freq,smooth(mean(P2_0,2),20),'color','g');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.1,300]);
    ylim([1e-18,2e-16])
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel('PSD')
subplot(1,4,2)
    plot(freq,smooth(mean(Pa_I+Pb_I,2),20),'color','k');
    hold on;
    plot(freq,smooth(mean(P2_I,2),20),'color','b');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.1,300]);
    ylim([1e-18,2e-16])
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel('PSD')
subplot(1,4,3)
    plot(freq,smooth(mean(Pa_E+Pb_E,2),20),'color','k');
    hold on;
    plot(freq,smooth(mean(P2_E,2),20),'color','r');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.1,300]);
    ylim([1e-18,2e-16])
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel('PSD')
subplot(1,4,4);
    plot(freq,smooth(mean(Pa_EI+Pb_EI,2),20),'color','k');
    hold on;
    plot(freq,smooth(mean(P2_EI,2),20),'color','m');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.1,300]);
    ylim([1e-18,2e-16])
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel('PSD')


figureNB;
% subplot(1,2,1);
    plot(freq,smooth(mean(Pa_0+Pb_0,2),20),'color','m');
    hold on;
    plot(freq,smooth(mean(Pa_I+Pb_I,2),20),'color','b');
    plot(freq,smooth(mean(Pa_E+Pb_E,2),20),'color','r');
    plot(freq,smooth(mean(Pa_EI+Pb_EI,2),20),'color','k');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.1,300]);
    ylim([1e-18,2e-16])
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel('PSD')

subplot(1,2,2);
    plot(freq,smooth(mean(P2_0,2),20),'color','m');
    hold on;
    plot(freq,smooth(mean(P2_I,2),20),'color','b');
    plot(freq,smooth(mean(P2_E,2),20),'color','r');
    plot(freq,smooth(mean(P2_EI,2),20),'color','k');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.1,300]);
    ylim([1e-18,2e-16])
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel('PSD')


figureNB;
    plot(freq,smooth(mean(P2_0,2),100)./smooth(mean(Pa_0+Pb_0,2),100),'color','m');
    hold on;
    plot(freq,smooth(mean(P2_I,2),100)./smooth(mean(Pa_I+Pb_I,2),100),'color','b');
    plot(freq,smooth(mean(P2_E,2),100)./smooth(mean(Pa_E+Pb_E,2),100),'color','r');
    plot(freq,smooth(mean(P2_EI,2),100)./smooth(mean(Pa_EI+Pb_EI,2),100),'color','k');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([0.1,300]);
    % ylim([1e-18,2e-16])
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel('PSD')

