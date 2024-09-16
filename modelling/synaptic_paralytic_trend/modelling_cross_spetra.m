
folder = 'C:\Users\brake\Desktop\synthetic_spikes_long';
% [pxy,cpxy,freq] = load_data(folder);


AP_model = @(f,p) p(4) + log10(exp(p(3))+1 ./ ((1+p(1)^2*(2*pi*f).^2).*(1+p(2)^2*(2*pi*f).^2)));
G_I = 10.^AP_model(freq,[12e-3,1e-3,-Inf,0]);
G_E = 10.^AP_model(freq,[2e-3,1.2e-3,-Inf,0]);
A_I = 1.4e-18;
A_E = 1e-18;

figureNB;
    plot(freq,mean(pxy.I,2),'color',blue,'LineWidth',1);
    hold on;
    plot(freq,mean(pxy.E,2),'color',red,'LineWidth',1);
    plot(freq,mean(pxy.EI,2),'color',(red+blue)*0.7,'LineWidth',1);
    plot(freq,mean(pxy.null,2),'color','k','LineWidth',1);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1,200]);
    ylim([1e-20,1e-16]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    ylabel('PSD')
    gcaformat

figureNB;
subplot(3,4,1);
    plot(freq,mean(pxy.I,2),'k','LineWidth',1);
    hold on;
    plot(freq,mean(real(cpxy.I),2),'k','LineWidth',1);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1,200]);
    ylim([1e-20,1e-16]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    ylabel('PSD')
    gcaformat
subplot(3,4,2);
    plot(freq,mean(pxy.I,2),'k','LineWidth',1);
    hold on;
    % plot(freq,mean(real(cpxy.I),2));
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1,200]);
    ylim([1e-20,1e-16]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    ylabel('PSD')
    gcaformat
    IC = 6.8*A_I*G_I;
    EC = 1.6*A_E*10.^AP_model(freq,[2e-3,0.7e-3,-Inf,0]);
    plot(freq,IC+EC,'g')
    plot(freq,IC,'-b','LineWidth',1)
    plot(freq,EC,'-r','LineWidth',1)
    % plot(freq,0.15*IC,'k')
subplot(3,4,3);
    plot(freq,mean(real(cpxy.I),2),'k','LineWidth',1);
    hold on;
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1,200]);
    ylim([1e-20,1e-16]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    ylabel('PSD')
    gcaformat
    plot(freq,0.15*IC,'-b','LineWidth',1)
subplot(3,4,5);
    plot(freq,mean(pxy.E,2),'k','LineWidth',1);
    hold on;
    plot(freq,mean(real(cpxy.E),2),'k','LineWidth',1);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1,200]);
    ylim([1e-20,1e-16]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    ylabel('PSD')
    gcaformat
subplot(3,4,6);
    plot(freq,mean(pxy.E,2),'k','Linewidth',1);
    hold on;
    % plot(freq,mean(real(cpxy.E),2));
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1,200]);
    ylim([1e-20,1e-16]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    ylabel('PSD')
    gcaformat
    % plot(freq,3*A_I*G_I+10*A_E*10.^AP_model(freq,[2.2e-3,1e-3,-Inf,0]),'k')
    % plot(freq,0*3*A_I*G_I+0.15*10*A_E*10.^AP_model(freq,[2.2e-3,1e-3,-Inf,0]),'k')
    IC = 2.5*A_I*G_I;
    EC = 8*A_E*10.^AP_model(freq,[2.3e-3,1e-3,-Inf,0]);
    plot(freq,IC+EC,'g')
    plot(freq,IC,'-b','LineWidth',1)
    plot(freq,EC,'-r','LineWidth',1)
    % plot(freq,0.13*EC,'g')
subplot(3,4,7);
    % plot(freq,mean(pxy.E,2));
    plot(freq,mean(real(cpxy.E),2),'k','LineWidth',1);
    hold on;
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1,200]);
    ylim([1e-20,1e-16]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    ylabel('PSD')
    gcaformat
    plot(freq,0.13*EC,'-r','LineWidth',1)
subplot(3,4,9);
    plot(freq,mean(pxy.EI,2),'k','LineWidth',1);
    hold on;
    plot(freq,mean(real(cpxy.EI),2),'k','LineWidth',1);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1,200]);
    ylim([1e-20,1e-16]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    ylabel('PSD')
    gcaformat
subplot(3,4,10);
    plot(freq,mean(pxy.EI,2),'k','Linewidth',1);
    hold on;
    % plot(freq,mean(real(cpxy.EI),2));
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1,200]);
    ylim([1e-20,1e-16]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    ylabel('PSD')
    gcaformat
    % plot(freq,3*A_I*G_I+10*A_E*10.^AP_model(freq,[2.2e-3,1e-3,-Inf,0]),'k')
    % plot(freq,0*3*A_I*G_I+0.15*10*A_E*10.^AP_model(freq,[2.2e-3,1e-3,-Inf,0]),'k')
    IC = 7*A_I*G_I;
    EC = 7*A_E*10.^AP_model(freq,[2.2e-3,1e-3,-Inf,0]);
    plot(freq,IC+EC,'g')
    plot(freq,IC,'-b','LineWidth',1)
    plot(freq,EC,'-r','LineWidth',1)
    % plot(freq,0.15*IC+0.13*EC,'g')
subplot(3,4,11);
    plot(freq,mean(real(cpxy.EI),2),'k','LineWidth',1);
    hold on;
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1,200]);
    ylim([1e-20,1e-16]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    ylabel('PSD')
    gcaformat
    plot(freq,0.15*IC+0.13*EC,'g')
    plot(freq,0.15*IC,'-b','LineWidth',1)
    plot(freq,0.13*EC,'-r','LineWidth',1)


subplot(3,4,4);
    plot(freq,mean(pxy.E,2),'LineWidth',1,'color',red*0.5+0.5);
    hold on;
    plot(freq,mean(pxy.E2,2),'LineWidth',1,'color',red);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1,200]);
    ylim([1e-19,1e-16]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    ylabel('PSD')
    IC = 2.5*A_I*G_I;
    EC = 8*A_E*10.^AP_model(freq,[2.3e-3,1e-3,-Inf,0]);
    plot(freq,IC+EC,'k')
    plot(freq,IC+(1+2*0.15)*EC,'k')
subplot(3,4,8);
    plot(freq,mean(pxy.I,2),'LineWidth',1,'color',blue*0.5+0.5);
    hold on;
    plot(freq,mean(pxy.I2,2),'LineWidth',1,'color',blue);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1,200]);
    ylim([1e-19,1e-16]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    ylabel('PSD')
    IC = 6.8*A_I*G_I;
    EC = 1.6*A_E*10.^AP_model(freq,[2e-3,0.7e-3,-Inf,0]);
    plot(freq,IC+EC,'k')
    plot(freq,(1+2*0.15)*IC+EC,'k')
subplot(3,4,12);
    plot(freq,mean(pxy.EI,2),'LineWidth',1,'color',(red+blue)*0.7*0.5+0.5);
    hold on;
    plot(freq,mean(pxy.EI2,2),'LineWidth',1,'color',(red+blue)*0.7);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1,200]);
    ylim([1e-19,1e-16]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    ylabel('PSD')
    IC = 6.8*A_I*G_I;
    EC = 8*A_E*10.^AP_model(freq,[2.3e-3,1e-3,-Inf,0]);
    plot(freq,IC+EC,'k')
    plot(freq,(1+2*0.15)*IC+(1+2*0.15)*EC,'k')

gcaformat(gcf);

function [pxy,cpxy,freq] = load_data(folder)
    [sa,X] = network_simulation_beluga.getHeadModel;
    F = dir(folder); F = F(3:end);
    fs = 16e3;

    winlen = fs;

    pxy = struct();
    cpxy = struct();

    pxy.I = zeros(winlen/2+1,length(F));
    pxy.E = zeros(winlen/2+1,length(F));
    pxy.null = zeros(winlen/2+1,length(F));
    pxy.EI = zeros(winlen/2+1,length(F));

    pxy.I2 = zeros(winlen/2+1,length(F));
    pxy.E2 = zeros(winlen/2+1,length(F));
    pxy.null2 = zeros(winlen/2+1,length(F));
    pxy.EI2 = zeros(winlen/2+1,length(F));

    cpxy.I = zeros(winlen/2+1,length(F));
    cpxy.E = zeros(winlen/2+1,length(F));
    cpxy.null = zeros(winlen/2+1,length(F));
    cpxy.EI = zeros(winlen/2+1,length(F));

    for i = 1:length(F)

        load(fullfile(folder,F(i).name,'spikeTimes_E_shuffle','simulation_data.mat'));
        eeg = detrend(network_simulation_beluga.getEEG(dipoles,sa,1e4),'constant');
        [temp,freq] = pwelch(eeg,hamming(winlen),winlen*0.9,winlen,fs);
        pxy.I(:,i) = sum(temp,2);
        [pxy.I2(:,i),freq] = pwelch(sum(eeg,2),hamming(winlen),winlen*0.9,winlen,fs);
        cpxy.I(:,i) = cpsd(eeg(:,1),eeg(:,2),hamming(winlen),winlen*0.9,winlen,fs);


        load(fullfile(folder,F(i).name,'spikeTimes_I_shuffle','simulation_data.mat'));
        eeg = detrend(network_simulation_beluga.getEEG(dipoles,sa,1e4),'constant');
        [temp,freq] = pwelch(eeg,hamming(winlen),winlen*0.9,winlen,fs);
        pxy.E(:,i) = sum(temp,2);
        [pxy.E2(:,i),freq] = pwelch(sum(eeg,2),hamming(winlen),winlen*0.9,winlen,fs);
        cpxy.E(:,i) = cpsd(eeg(:,1),eeg(:,2),hamming(winlen),winlen*0.9,winlen,fs);

        load(fullfile(folder,F(i).name,'spikeTimes_E_I_shuffle','simulation_data.mat'));
        eeg = detrend(network_simulation_beluga.getEEG(dipoles,sa,1e4),'constant');
        [temp,freq] = pwelch(eeg,hamming(winlen),winlen*0.9,winlen,fs);
        pxy.null(:,i) = sum(temp,2);
        [pxy.null2(:,i),freq] = pwelch(sum(eeg,2),hamming(winlen),winlen*0.9,winlen,fs);
        cpxy.null(:,i) = cpsd(eeg(:,1),eeg(:,2),hamming(winlen),winlen*0.9,winlen,fs);

        load(fullfile(folder,F(i).name,'spikeTimes_E_I_correlated','simulation_data.mat'));
        eeg = detrend(network_simulation_beluga.getEEG(dipoles,sa,1e4),'constant');
        [temp,freq] = pwelch(eeg,hamming(winlen),winlen*0.9,winlen,fs);
        pxy.EI(:,i) = sum(temp,2);
        [pxy.EI2(:,i),freq] = pwelch(sum(eeg,2),hamming(winlen),winlen*0.9,winlen,fs);
        cpxy.EI(:,i) = cpsd(eeg(:,1),eeg(:,2),hamming(winlen),winlen*0.9,winlen,fs);
    end

    freq(1) = [];
    pxy.null(1,:) = [];
    pxy.null2(1,:) = [];
    pxy.E(1,:) = [];
    pxy.E2(1,:) = [];
    pxy.I(1,:) = [];
    pxy.I2(1,:) = [];
    pxy.EI(1,:) = [];
    pxy.EI2(1,:) = [];
    cpxy.I(1,:) = [];
    cpxy.E(1,:) = [];
    cpxy.null(1,:) = [];
    cpxy.EI(1,:) = [];
end