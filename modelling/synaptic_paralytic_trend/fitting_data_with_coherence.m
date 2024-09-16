
baseFolder = 'C:\Users\brake\Documents\GitHub\apEEG';

sigx2 = 3;
fs = 16e3;
L = 16e3;
f = (fs/L:fs/L:fs/2)';

% Load average cross-spectrum for different pairwise distances
load('C:\Users\brake\Desktop\chain1.mat','dValues')

[~,Rxx] = compute_AP_spectra(sigx2);

% Get density of surface area as a function of distance
load(fullfile(baseFolder,'data_files','pairwise_distance.mat'));
N2 = 16e9;
dV = mean(diff(dValues));
dN = interp1(0.5*(rValues(1:end-1)+rValues(2:end)),diff(total_area)./diff(rValues'),dValues,'linear','extrap');
dN = mean(dN,2)'.*N2./mean(total_area(end,:));
A = exp(-dValues.^2/(2*sigx2)).*dN*dV/N2;

% Compute corrected average cross spectrum
load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\cross_spectra\control_chains.mat','Pxy_D','count_D')
mu = Pxy_D./count_D;
Sxy_control = (nansum(mu.*A,2));

load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\cross_spectra\noE_chains.mat','Pxy_D','count_D')
mu = Pxy_D./count_D;
Sxy_noE = (nansum(mu.*A,2));


load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\cross_spectra\noI_chains.mat','Pxy_D','count_D')
mu = Pxy_D./count_D;
Sxy_noI = (nansum(mu.*A,2));


AP_model = @(f,p) p(4) + log10(exp(p(3))+1 ./ ((1+p(1)^2*(2*pi*f).^2).*(1+p(2)^2*(2*pi*f).^2)));
GI = 1.2*5.4e-16*10.^AP_model(f,[20e-3,4e-3,-Inf,0]);
GE = 0.4*5.4e-16*10.^AP_model(f,[2e-3,1.2e-3,-Inf,0]);

[f0,S] = import_paralytic_data;
S = fillgaps(nanmean(S(:,1:8),2),5);


dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Nature Communications\final_submission\manuscript_source_data';
load(fullfile(dataFolder,'cortex_anatomy','anatomy_cortical_pairwise_distance_distribution.mat'));
signed_area = A;
N = 16e9;
dMids = 0.5*(rValues(2:end)+rValues(1:end-1));
nrnCount = mean(diff(signed_area),2)*200000;
nrnCount(end) = N-sum(nrnCount(1:end-1));
corr_kernel = @(d) exp(-d.^2/6);
rho_bar = sum(corr_kernel(dMids).*nrnCount)/sum(nrnCount');


[f,Rxx,Rxy] = compute_AP_spectra;
R = 0.2;
lam = 1;
low_noise = 1e-3;

figureNB;
subplot(1,3,1);
    plot(f,N*(GI+GE)+N*(N-1)*rho_bar*0.15*GI+N*(N-1)*rho_bar*0.15*GE);
    hold on;
    hold on;
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([1,1e3]);

    dS = [0,2,5,10,25,50,Inf]*1e-3;
    for i = 1:length(dS)
        sig = dS(i);
        B = exp(-2*(pi*f(:)*sig).^2);
        plot(f,lam*N*Rxx + R*lam*N*(N-1)*B.*Rxy,'color',blue,'LineWidth',1)
    end
subplot(1,3,2)
    plot(f0,S,'k');
    hold on;
    plot(f,low_noise+f*0);
    plot(f,N*(GI+GE)+N*(N-1)*rho_bar*2e-2*GI+N*(N-1)*rho_bar*1e-4*GE+low_noise,'-b');
    hold on;
    plot(f,N*GI+N*(N-1)*rho_bar*2e-2*GI,'--b');
    plot(f,N*GE+N*(N-1)*rho_bar*1e-4*GE,'--b');
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([1,100]);

subplot(1,3,3)
    plot(f0,S,'k');
    hold on;
    plot(f,low_noise+f*0);
    plot(f,N*(GI+GE)+N*(N-1)*rho_bar*2e-2*GI+N*(N-1)*rho_bar*1e-4*GE+low_noise);
    lam = 1;
    Y_AP = Sxy_noE*(N*0.15)*(N*0.15-1)*2e-2*lam+Sxy_noI*(N*0.85)*(N*0.85-1)*1e-4*lam;
    dS = [0,2,5,10,25,50,Inf]*1e-3;
    for i = 1:length(dS)
        sig = dS(i);
        B = exp(-2*(pi*f(:)*sig).^2);
        plot(f,lam*N*Rxx + Y_AP.*B,'b');
    end
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlim([1,100]);


B = (1+10*exp(-2*(pi*(f(:)-11)*8e-3*sqrt(2)).^2));

fx = 50;
B_syn = (1+6*exp(-2*(pi*(f(:)-fx).^2)/(fx.^2)*40));
B_AP = exp(-2*(pi*(f(:)-fx).^2)/(fx.^2)*40);
osc_syn = (N*(GI+GE)+N*(N-1)*rho_bar*r_I*GI+N*(N-1)*rho_bar*r_E*GE).*B_syn;
osc_AP_hi = lam*N*Rxx + Sxy_control*N*(N-1)*R*lam.*B_AP;


R = 0.2;
lam = 1;

r_E = 1e-4;
r_I = 2e-2;

figureNB;
    plot(f0,S,'k');
    hold on;
    plot(f,low_noise+f*0);
    plot(f,N*(GI+GE)+N*(N-1)*rho_bar*r_I*GI+N*(N-1)*rho_bar*r_E*GE);
    plot(f,osc_syn);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    ylim([1e-4,1e2]);

    % Y_AP = Sxy_noE*(N*0.15)*(N*0.15-1)*r_I*lam+Sxy_noI*(N*0.85)*(N*0.85-1)*r_E*lam;
    % osc_AP_lo = Y_AP.*B_AP;;

    % Y_AP = Sxy_noE*(N*0.15)*(N*0.15-1)*R*lam+Sxy_noI*(N*0.85)*(N*0.85-1)*R*lam;
    % osc_AP_hi = Y_AP.*B_AP;;
    plot(f,osc_AP_hi,'b');
    plot(f,osc_AP_lo,'b');



F = 10.^linspace(0,3,100);
RR = 10.^linspace(-4,-1,5);
figureNB;
for iR = 1:5
    R = RR(iR);
    for iF = 1:length(F)
        fx = F(iF);
        B_syn = (1+6*exp(-2*(pi*(f(:)-fx).^2)/(fx.^2)*40));
        B_AP = exp(-2*(pi*(f(:)-fx).^2)/(fx.^2)*40);
        osc_syn = (N*(GI+GE)+N*(N-1)*rho_bar*r_I*GI+N*(N-1)*rho_bar*r_E*GE).*B_syn;
        jj = interp1(f,1:length(f),fx,'nearest','extrap');

        A_AP(iF) = osc_AP_hi(jj);
        A_syn(iF) = osc_syn(jj);
    end
    A_AP = lam*N*Rxx + Sxy_control*N*(N-1)*R*lam;
    A_syn = (N*(GI+GE)+N*(N-1)*rho_bar*r_I*GI+N*(N-1)*rho_bar*r_E*GE).*(1+6);
    subplot(1,2,1);
        plot(f,A_AP);
        hold on;
        plot(f,A_syn);
        set(gca,'xscale','log');
        set(gca,'yscale','log');
    subplot(1,2,2)
        plot(f,A_AP./(A_AP+A_syn));
        hold on;
        set(gca,'xscale','log');
end

