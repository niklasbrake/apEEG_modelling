% function figure5
% FIGURE5 generates the panels in figure 3 of the manuscript.
% See also import_paralytic_data, fittingmodel.

    baseFolder = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(baseFolder,'auxiliary_functions'));
    addpath(fullfile(baseFolder,'modelling'));


    [full_model,AP_model] = fittingmodel('eq6');
    [f0,S] = import_paralytic_data;
    low_noise = 1e-3;

    figureNB(13.2,12);
    subplot(2,2,1)
        plot(f0,S(:,9:16),'color',[0.5,0.5,0.5]);
        hold on;
        plot(f0,S(:,1:8),'color','k');
        plot([1,3e3],low_noise*[1,1],'color',red,'LineWidth',1)
        % Psyn = 10.^AP_model(f0,[3e-3,1e-3,-Inf,4.6]);
        % plot(f0,Psyn,'color',blue,'LineStyle','--')
        % Psyn = 10.^AP_model(f0,[20e-3,4e-3,-Inf,3.6]);
        % plot(f0,Psyn,'color',blue,'LineStyle','--')
        % Psyn = 10.^AP_model(f0,[20e-3,4e-3,-Inf,3.6])+10.^AP_model(f0,[3e-3,1e-3,-Inf,4.6]);
        % plot(f0,Psyn,'color',blue,'LineWidth',1)
        set(gca,'yscale','log')
        xlabel('Frequency (Hz)')
        ylabel(['PSD (' char(956) 'V^2/Hz)'])
        title('Unparalyzed','FontWeight','normal')
        xlim([0,100]);
        gcaformat;
        set(gca,'FontSize',8)


        dataFolder = 'E:\Research_Projects\004_Propofol\manuscript\Nature Communications\final_submission\manuscript_source_data';
        load(fullfile(dataFolder,'cortex_anatomy','anatomy_cortical_pairwise_distance_distribution.mat'));
        signed_area = A;
        N = 16e9;
        dMids = 0.5*(rValues(2:end)+rValues(1:end-1));
        nrnCount = mean(diff(signed_area),2)*200000;
        nrnCount(end) = N-sum(nrnCount(1:end-1));
        corr_kernel = @(d) exp(-d.^2/6);
        rho_bar = sum(corr_kernel(dMids).*nrnCount)/sum(nrnCount');
        AP_model = @(f,p) p(4) + log10(exp(p(3))+1 ./ ((1+p(1)^2*(2*pi*f).^2).*(1+p(2)^2*(2*pi*f).^2)));
        GI = 1.2*5.4e-16*10.^AP_model(f0,[20e-3,4e-3,-Inf,0]);
        GE = 0.4*5.4e-16*10.^AP_model(f0,[2e-3,1.2e-3,-Inf,0]);
        set(gca,'xscale','log')
        set(gca,'yscale','log')

    subplot(2,2,2)
        plot(f0,S(:,1:8)-low_noise,'color','k');
        hold on;
        % plot([1,3e3],low_noise*[1,1],'color',red,'LineWidth',1)
        plot(f0,N*(GI+GE)+N*(N-1)*rho_bar*2e-2*GI+N*(N-1)*rho_bar*1e-4*GE+0*low_noise,'-b','LineWidth',1);
        hold on;
        plot(f0,N*GI+N*(N-1)*rho_bar*2e-2*GI,'--b');
        plot(f0,N*GE+N*(N-1)*rho_bar*1e-4*GE,'--b');
        set(gca,'yscale','log')
        xlabel('Frequency (Hz)')
        ylabel(['PSD (' char(956) 'V^2/Hz)'])
        title('20 mg cisatracurium','FontWeight','normal')
        xlim([0,100]);
        gcaformat;
        set(gca,'FontSize',8)

        set(gca,'xscale','log')
        set(gca,'yscale','log')
% end

    subplot(2,2,3)
        plot(f0,S(:,1:8),'color','k');
        hold on;
        plot([1,3e3],low_noise*[1,1],'color',red,'LineWidth',1)
        set(gca,'yscale','log')
        xlabel('Frequency (Hz)')
        ylabel(['PSD (' char(956) 'V^2/Hz)'])
        title('Unparalyzed','FontWeight','normal')
        xlim([0,100]);
        gcaformat;
        set(gca,'FontSize',8)

    [f,Rxx,Rxy] = compute_AP_spectra;

    N2 = 16e9; % Total neuron count
    % low_noise = (8e-3).^2;
    low_noise = 1e-3;
    R = 0.01;
    lam = 1;

    dS = [11,50,100,Inf]*1e-3;
    for i = 1:length(dS)
        sig = dS(i);
        B = exp(-2*(pi*f(:)*sig).^2);
        plot(f,lam*N2*Rxx + R*lam*N2*(N2-1)*B.*Rxy,'color',blue,'LineWidth',1)
    end
    set(gca,'xscale','log')
    set(gca,'yscale','log')

    [f0,S] = import_paralytic_data;
    S = fillgaps(nanmean(S(:,1:8),2),5);



    low_noise = 1e-3;
    R = 0.01;
    lam = 1;

    dR = 10.^linspace(-4,-1,100);
    dS = [11,25,50,100]*1e-3;

    [XX,YY] = meshgrid(f0,dR);


    figureNB(5.6,5);
    for i = 1:length(dS)
        sig = dS(i);
        B = exp(-2*(pi*f(:)*sig).^2);
        for j = 1:length(dR)
            Y  = interp1(f,lam*N2*Rxx + dR(j)*N2*(N2-1)*B.*Rxy,f0);
            D(j,:) = Y(:)./S(:);
        end
        subplot(4,1,i)
        surf(XX,YY,0*D,log10(100*D),'LineStyle','none')
        view([0,90]);
        hold on;
        [C,h] = contour(XX,YY,log10(100*D),[0,1,2],'color','r','linewidth',1);
        % imagesc(f0,dR,log10(100*D));
        % colorbar;
        set(gca,'CLim',[-3,2])
        colormap(flip(gray(100)));
        axis xy;
        set(gca,'yscale','log');
        xticks([])
        xlim([0,95])
        yticks(10.^[-4,-3,-2,-1])
        % yticks(10.^[-4,-1])
        ylim([8e-5,1e-1])
        set(gca,'YMinorTick','off');
        grid off;
        if(i>1)
            yticklabels({});
        else
            ylabel('\lambda R_{max}');
        end
    end


    % xticks([0:10:100])
    xticks([0,30,60,90])
    xlim([0,95])
    xlabel('Frequency (Hz)')
    gcaformat(gcf);