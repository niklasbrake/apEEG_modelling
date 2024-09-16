k = 6;

p0 = [0.018,0.0034,-0.1,-0.32];
[~,synFun] = fittingmodel('eq6');
I = 10.^synFun(f0,[0.020,0.004,-Inf,1]); I = I*5e-16/I(1);
E = 10.^synFun(f0,[0.003,0.0001,-Inf,1]); E = E*2.2e-16/E(1);


N = 16e9;
rhoI = 70e-7;
rhoE = 1e-7;

    [full_model,AP_model] = fittingmodel('eq6');
    [f0,S] = import_paralytic_data;
    low_noise = 1e-3;

    figureNB(13.2,5.8);
    axes('Position',[0.08, 0.14, 0.38, 0.8])
        plot(f0,S(:,9:16),'color',[0.5,0.5,0.5]);
        hold on;
        plot(f0,S(:,1:8),'color','k');
        plot([1,3e3],low_noise*[1,1],'color',red,'LineWidth',1)
        % I_N = I*N + 7e-6*I*N^2;
        % E_N = E*N + 7e-6*E*N^2;
        % plot(f0,I_N,'color',blue,'LineStyle','--')
        % plot(f0,E_N,'color',blue,'LineStyle','--')
        % plot(f0,E_N+I_N,'color',blue,'LineWidth',1)
        set(gca,'yscale','log')
        xlabel('Frequency (Hz)')
        ylabel(['PSD (' char(956) 'V^2/Hz)'])
        xlim([0,100]);
        title('Unparalyzed','FontWeight','normal')
        gcaformat;
        set(gca,'FontSize',8)

    axes('Position',[0.58, 0.14, 0.38, 0.8])
        plot(f0,S(:,1:8),'color','k');
        hold on;
        plot([1,3e3],low_noise*[1,1],'color',red,'LineWidth',1)
        I_N = I*N + rhoI*I*N^2;
        E_N = E*N + rhoE*E*N^2;
        plot(f0,I_N,'color',blue,'LineStyle','--')
        plot(f0,E_N,'color',blue,'LineStyle','--')
        plot(f0,E_N+I_N,'color',blue,'LineWidth',1)
        set(gca,'yscale','log')
        xlabel('Frequency (Hz)')
        ylabel(['PSD (' char(956) 'V^2/Hz)'])
        xlim([0,100]);
        title('20 mg cisatracurium','FontWeight','normal')
        gcaformat;
        set(gca,'FontSize',8)



% Get baseline EEG spectrum from propofol cohort
data = load(fullfile(dataFolder,'EEG_data','electrode2_Cz.mat'));
data.baseline = squeeze(nanmedian(data.psd(:,data.time<-1,:),2));

load(fullfile(dataFolder,'cortex_anatomy','anatomy_cortical_pairwise_distance_distribution.mat'));
signed_area = A;
total_area = B;
N = 16e9;
dMids = 0.5*(rValues(2:end)+rValues(1:end-1));
nrnCount = mean(diff(signed_area),2)*200000;
nrnCount(end) = N-sum(nrnCount(1:end-1));
corr_kernel = @(d) exp(-d.^2/3);
rho_bar = sum(corr_kernel(dMids).*nrnCount)/sum(nrnCount');


rhoI/rho_bar
rhoE/rho_bar