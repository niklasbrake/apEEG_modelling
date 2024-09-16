[sa,X] = network_simulation_beluga.getHeadModel;
folder = 'C:\Users\brake\Desktop\synthetic_spikes';
F = dir(folder); F = F(3:end);
fs = 16e3;


winlen = fs;

cxy_I = zeros(winlen/2+1,length(F));
cxy_E = zeros(winlen/2+1,length(F));
cxy_0 = zeros(winlen/2+1,length(F));
cxy_EI = zeros(winlen/2+1,length(F));

for i = 1:length(F)

    load(fullfile(folder,F(i).name,'spikeTimes_E_shuffle','simulation_data.mat'));
    eeg = detrend(network_simulation_beluga.getEEG(dipoles,sa,1e4),'constant');
    [cxy_I(:,i),freq] = mscohere(eeg(:,1),eeg(:,2),hamming(winlen),winlen*0.9,winlen,fs);


    load(fullfile(folder,F(i).name,'spikeTimes_I_shuffle','simulation_data.mat'));
    eeg = detrend(network_simulation_beluga.getEEG(dipoles,sa,1e4),'constant');
    [cxy_E(:,i),freq] = mscohere(eeg(:,1),eeg(:,2),hamming(winlen),winlen*0.9,winlen,fs);

    load(fullfile(folder,F(i).name,'spikeTimes_E_I_shuffle','simulation_data.mat'));
    eeg = detrend(network_simulation_beluga.getEEG(dipoles,sa,1e4),'constant');
    [cxy_0(:,i),freq] = mscohere(eeg(:,1),eeg(:,2),hamming(winlen),winlen*0.9,winlen,fs);

    load(fullfile(folder,F(i).name,'spikeTimes_E_I_correlated','simulation_data.mat'));
    eeg = detrend(network_simulation_beluga.getEEG(dipoles,sa,1e4),'constant');
    [cxy_EI(:,i),freq] = mscohere(eeg(:,1),eeg(:,2),hamming(winlen),winlen*0.9,winlen,fs);

end

C_I = mean(smoothdata(cxy_I-mean(cxy_0(:)),'movmean',3),2);
C_E = mean(smoothdata(cxy_E-mean(cxy_0(:)),'movmean',3),2);
C_EI = mean(smoothdata(cxy_EI-mean(cxy_0(:)),'movmean',3),2);
C_0 = mean(smoothdata(cxy_0-mean(cxy_0(:)),'movmean',3),2);

figureNB;
    plot(freq,C_I);
    hold on;
    plot(freq,C_E);
    plot(freq,C_EI);
    plot(freq,C_0);
    set(gca,'xscale','log');
    ylabel('Magnitude-squared coherence')
    xlabel('Frequency (Hz)');






pxy_I = zeros(winlen/2+1,length(F));
pxy_E = zeros(winlen/2+1,length(F));
pxy_0 = zeros(winlen/2+1,length(F));
pxy_EI = zeros(winlen/2+1,length(F));

for i = 1:length(F)

    load(fullfile(folder,F(i).name,'spikeTimes_E_shuffle','simulation_data.mat'));
    eeg = detrend(network_simulation_beluga.getEEG(dipoles,sa,1e4),'constant');
    [pxy_I(:,i),freq] = cpsd(eeg(:,1),eeg(:,2),hamming(winlen),winlen*0.9,winlen,fs);


    load(fullfile(folder,F(i).name,'spikeTimes_I_shuffle','simulation_data.mat'));
    eeg = detrend(network_simulation_beluga.getEEG(dipoles,sa,1e4),'constant');
    [pxy_E(:,i),freq] = cpsd(eeg(:,1),eeg(:,2),hamming(winlen),winlen*0.9,winlen,fs);

    load(fullfile(folder,F(i).name,'spikeTimes_E_I_shuffle','simulation_data.mat'));
    eeg = detrend(network_simulation_beluga.getEEG(dipoles,sa,1e4),'constant');
    [pxy_0(:,i),freq] = cpsd(eeg(:,1),eeg(:,2),hamming(winlen),winlen*0.9,winlen,fs);

    load(fullfile(folder,F(i).name,'spikeTimes_E_I_correlated','simulation_data.mat'));
    eeg = detrend(network_simulation_beluga.getEEG(dipoles,sa,1e4),'constant');
    [pxy_EI(:,i),freq] = cpsd(eeg(:,1),eeg(:,2),hamming(winlen),winlen*0.9,winlen,fs);

end


figureNB;
    plot(freq,mean(abs(pxy_0),2),'color','k');
    hold on;
    plot(freq,mean(abs(pxy_E),2),'color',red);
    plot(freq,mean(abs(pxy_I),2),'color',blue);
    plot(freq,mean(abs(pxy_EI),2),'color',(red+blue)*0.7);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1,1e3]);
    ylim([1e-22,1e-16]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel('PSD')