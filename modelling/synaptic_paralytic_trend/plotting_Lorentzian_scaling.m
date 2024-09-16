% [sa,X] = network_simulation_beluga.getHeadModel;
folder = 'C:\Users\brake\Desktop\synthetic_spikes\mCombo7_RE1_RI1';
F = dir(folder); F = F(3:end);
fs = 16e3;


winlen = fs;



pxy_I = zeros(winlen/2+1,length(F));
pxy_E = zeros(winlen/2+1,length(F));
pxy_0 = zeros(winlen/2+1,length(F));
pxy_EI = zeros(winlen/2+1,length(F));

load(fullfile(folder,'spikeTimes_E_shuffle','simulation_data.mat'));
eeg = detrend(network_simulation_beluga.getEEG(dipoles,sa,1e4),'constant');
[pxy_I,freq] = pwelch(eeg,hamming(winlen),winlen*0.9,winlen,fs);
pxy_I = sum(pxy_I,2);
[pxy_I2,freq] = pwelch(sum(eeg,2),hamming(winlen),winlen*0.9,winlen,fs);
MS_I2 = mscohere(eeg(:,1),eeg(:,2),hamming(winlen),winlen*0.9,winlen,fs);


load(fullfile(folder,'spikeTimes_I_shuffle','simulation_data.mat'));
eeg = detrend(network_simulation_beluga.getEEG(dipoles,sa,1e4),'constant');
[pxy_E,freq] = pwelch(eeg,hamming(winlen),winlen*0.9,winlen,fs);
pxy_E = sum(pxy_E,2);
[pxy_E2,freq] = pwelch(sum(eeg,2),hamming(winlen),winlen*0.9,winlen,fs);
MS_E2 = mscohere(eeg(:,1),eeg(:,2),hamming(winlen),winlen*0.9,winlen,fs);

load(fullfile(folder,'spikeTimes_E_I_shuffle','simulation_data.mat'));
eeg = detrend(network_simulation_beluga.getEEG(dipoles,sa,1e4),'constant');
[pxy_0,freq] = pwelch(eeg,hamming(winlen),winlen*0.9,winlen,fs);
pxy_0 = sum(pxy_0,2);
[pxy_02,freq] = pwelch(sum(eeg,2),hamming(winlen),winlen*0.9,winlen,fs);
MS_02 = mscohere(eeg(:,1),eeg(:,2),hamming(winlen),winlen*0.9,winlen,fs);

load(fullfile(folder,'spikeTimes_E_I_correlated','simulation_data.mat'));
eeg = detrend(network_simulation_beluga.getEEG(dipoles,sa,1e4),'constant');
[pxy_EI,freq] = pwelch(eeg,hamming(winlen),winlen*0.9,winlen,fs);
pxy_EI = sum(pxy_EI,2);
[pxy_EI2,freq] = pwelch(sum(eeg,2),hamming(winlen),winlen*0.9,winlen,fs);
MS_EI2 = mscohere(eeg(:,1),eeg(:,2),hamming(winlen),winlen*0.9,winlen,fs);


[full_model,AP_model] = fittingmodel('eq6');

G_I = 10.^AP_model(freq,[14e-3,2e-3,-Inf,1]);
G_E = 10.^AP_model(freq,[2.8e-3,0.75e-3,-Inf,1]);


aI1 = 5e-15;
aI2 = 2.5e-15;
aE1 = 2e-13;
aE2 = 0.5e-13;
N = 2;

figureNB(15,15)
subplot(2,2,1)
    plot(freq,pxy_0,'LineWidth',1,'color',[0,0,0]*0.5+0.5);
    hold on;
    plot(freq,pxy_02,'LineWidth',1,'color',[0,0,0]);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1,300]);
    ylim([1e-19,1e-16]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel('PSD')
    % plot(freq,7e-15*10.^AP_model(freq,[13e-3,4e-3,-Inf,1])+10e-14*10.^AP_model(freq,[2.5e-3,0.6e-3,-Inf,1]),'k')

subplot(2,2,2)
    plot(freq,pxy_E,'LineWidth',1,'color',red*0.5+0.5);
    hold on;
    plot(freq,pxy_E2,'LineWidth',1,'color',red);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1,300]);
    ylim([1e-19,1e-16]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel('PSD')
    % plot(freq,2e-15*G_I+4e-13*G_E,'k')
    % plot(freq,2e-15*G_I+5e-13*G_E,'k')

subplot(2,2,3)
    plot(freq,pxy_I,'LineWidth',1,'color',blue*0.5+0.5);
    hold on;
    plot(freq,pxy_I2,'LineWidth',1,'color',blue);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1,300]);
    ylim([1e-19,1e-16]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel('PSD')
    % plot(freq,11e-15*G_I+1e-13*10.^AP_model(freq,[2.5e-3,0.6e-3,-Inf,1]),'k')
    % plot(freq,(11e-15+5.5e-15)*G_I+1e-13*10.^AP_model(freq,[2.5e-3,0.6e-3,-Inf,1]),'k')

subplot(2,2,4)
    plot(freq,pxy_EI,'LineWidth',1,'color',(red+blue)*0.7*0.5+0.5);
    hold on;
    plot(freq,pxy_EI2,'LineWidth',1,'color',(red+blue)*0.7);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1,300]);
    ylim([1e-19,1e-16]);
    xticks([1,10,100]);
    xticklabels([1,10,100]);
    xlabel('Frequency (Hz)');
    ylabel('PSD')
    % plot(freq,N*aI1*G_I+4e-13*G_E,'k')
  %   plot(freq,12e-15*G_I+5.9e-13*10.^AP_model(freq,[2.8e-3,0.8e-3,-Inf,1]),'k')
  % 
    % plot(freq,(N*aI1+N*(N-1)*aI2)*G_I+(N*aE1+N*(N-1)*aE2)*G_E,'k')




[full_model,AP_model] = fittingmodel('eq6');
% G_I = 10.^AP_model(freq,[20e-3,4e-3,-Inf,1]).*exp(-2*(pi*f(:)*2e-3).^2);;
% G_E = 10.^AP_model(freq,[2.8e-3,0.75e-3,-Inf,1]).*exp(-2*(pi*f(:)*2e-3).^2);;

G_I = 10.^AP_model(freq,[14e-3,2e-3,-Inf,1]).*exp(-2*(pi*freq(:)*4e-3).^2);;
G_E = 10.^AP_model(freq,[2.8e-3,0.75e-3,-Inf,1]).*exp(-2*(pi*freq(:)*4e-3).^2);;


aI1 = 5e-15;
aI2 = 2.5e-15;
aE1 = 2e-13;
aE2 = 0.5e-13;

N = 16e9;
I1 = N*aI1*G_I;
E1 = N*aE1*G_E;
I2 = N*(N-1)*aI2/0.1*G_I;
E2 = N*(N-1)*aE2/0.1*G_E;

I1 = interp1(freq,I1,f0(2:end));
E1 = interp1(freq,E1,f0(2:end));
I2 = interp1(freq,I2,f0(2:end));
E2 = interp1(freq,E2,f0(2:end));


low_noise = 1e-3;


[f0,S] = import_paralytic_data;
S = fillgaps(nanmean(S(:,1:8),2),5);
% S = fillgaps(nanmean(S(:,9:16),2),5);

A = 1+exp(-(f0(2:end)-12.5).^2./3)*7;
B = 1+exp(-(f0(2:end)-24.5).^2./11)*3.5;
C = 1+exp(-(f0(2:end)-55).^2./200)*1.6;
oof = 49./f0(2:end).^3.3;
S0 = S(2:end)./A(:)./B(:)./C(:)-oof(:);

y = S0-E1(:)-I1(:)-low_noise;
x = [I2(:),E2(:)];
% t = 10.^linspace(log10(f0(2)),log10(f0(end)),100);
% x = interp1(f0(2:end),x,t(:),'linear','extrap');
% y = interp1(f0(2:end),y,t(:),'linear','extrap');

% B = lsqnonneg(x,y)
FT = fitlm(x,y,'y~x1+x2-1','weights',1./y)
figureNB;
    plot(f0,S,'k');
    hold on;
    % plot(f0(2:end),S0,'color',[0.5,0.5,0.5]);
    plot(f0(2:end),x*FT.Coefficients{:,1}+low_noise+E1(:)+I1(:),'b')
    % set(gca,'xscale','log')
    set(gca,'yscale','log')

    plot(f0(2:end),x(:,1)*FT.Coefficients{1,1},'--b')
    hold on;
    plot(f0(2:end),x(:,2)*FT.Coefficients{2,1},'--b')
    line([1,100],[1,1]*low_noise,'color','r')



[f,Rxx,Rxy] = compute_AP_spectra;
N2 = 16e9; % Total neuron count
low_noise = (8e-3).^2;
R = 0.2;
lam = 1;

figureNB(6.9,5.7);
    plot(f0,S,'color',0*[0.6,0.6,0.6],'Linewidth',1);
    hold on;
    plot([1,3e3],low_noise*[1,1],'color',red,'LineWidth',1.5,'LineStyle','-')

    dS = [0,2,5,10,25,50,Inf]*1e-3;
    for i = 1:length(dS)
        sig = dS(i);
        B = exp(-2*(pi*f(:)*sig).^2);
        plot(f,lam*N2*Rxx + R*lam*N2*(N2-1)*B.*Rxy,'color',blue,'LineWidth',1)
    end
    set(gca,'xscale','log')
    set(gca,'yscale','log');
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    xlim([1,3e3])
    ylim([1e-7,1e1]);
    xticks([1,10,100,1000]);
    xticklabels([1,10,100,1000]);
    yticks([1e-6,1e-4,1e-2,1e0])
    gcaformat(gca,true,8);