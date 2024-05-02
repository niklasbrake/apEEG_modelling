function figureS4
% FIGURES4 generates the panels in figure S4 of the manuscript.
% See also compute_AP_spectra, import_Scheer2006
    baseFolder = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(baseFolder,'auxiliary_functions'));
    addpath(fullfile(baseFolder,'modelling'));

    % Run for sigma_x^2 = 1 mm^2
    main(1);

    % Run for sigma_x^2 = 10 mm^2
    main(10);
end
function main(sigx2)
    [f,Rxx,Rxy] = compute_AP_spectra(sigx2);
    N2 = 16e9; % Total neuron count

    [f0,S] = import_Scheer2006;
    low_noise = (8e-3).^2;
    R = 0.2;
    lam = 1;

    figureNB(6.9,5.7);
        plot(f0,S,'color',0*[0.6,0.6,0.6],'Linewidth',1);
        hold on;
        plot([1,3e3],low_noise*[1,1],'color',red,'LineWidth',1.5,'LineStyle','-')

        dS = [0,1e-3*2.^(0:5),Inf];
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


    PN = @(A,lam,sig) lam*N2*Rxx + A*lam*N2*(N2-1)*exp(-(2*pi*f(:)*sig).^2).*Rxy;

    lam = 10.^linspace(-1,2,50);
    sig = 10.^linspace(-3,-1,100)*1e3;
    A = 10.^linspace(-3,0,5);
    M = zeros(length(lam),length(sig));
    [XX,YY] = meshgrid(sig,lam);

    w = 0.125;
    figureNB(13.2,2.8);
    for i = 1:length(A)
        axes('Position',[0.07+(0.91-w-0.07)*(i-1)/4,0.3,0.13,0.6]);
        for j = 1:length(lam)
            for k = 1:length(sig)
                P = PN(A(i),lam(j),1e-3*sig(k));
                [M(j,k),I(j,k)] = max(P(f>30));
            end
        end
        surf(XX,YY,0*M,log10(M),'LineStyle','none')
        view([0,90]);
        hold on;
        C = contour(XX,YY,log10(M),log10([low_noise,low_noise]),'color',red,'linewidth',1);
        if(i==4)
            fill([10,90,90,10],[0.11,0.11,4,4],'k','LineStyle','--','FaceColor','none');
        end
        set(gca,'CLim',[-4,2]);
        ylim(10.^[-1,2])
        xlim(10.^[0,2])
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        line(1e3*[1e-3,0.1],[100,100],[0,0]-1e6,'color','k','LineWidth',0.75)
        line(1e3*[0.1,0.1],[0.1,100],[0,0]-1e6,'color','k','LineWidth',0.75)

        xticks([1,10,100])
        xticklabels([1,10,100])
        yticks([0.1,1,10,100])
        yticklabels([0.1,1,10,100])
        gcaformat
        xlabel('Jitter (ms)','FontSize',8)
        if(i==1)
            ylabel('Firing rate (Hz)','FontSize',8)
        end
        set(gca,'FontSize',8);
        colormap(flip(bone))
        drawnow;
    end
    C = colorbar;
    C.Position = [0.92,0.3,0.01,0.6];
    C.Ticks = [-4:2:2];
    C.Label.Position = [5,-0.8,0];
    C.TickLabels = {'10^{-4}','10^{-2}','10^{0}','10^{2}'};
    C.Label.String = ['Max PSD (' char(956) 'V^2/Hz)'];

    % low_noise = 2*low_noise;
    PN = @(A,lam,sig) lam*N2*Rxx + A*lam*N2*(N2-1)*exp(-(2*pi*f(:)*sig).^2).*Rxy;

    lam = 10.^linspace(-1,2,50);
    sig = 10.^linspace(-3,-1,100)*1e3;
    A = 10.^linspace(-3,0,5);
    M = zeros(length(lam),length(sig));
    [XX,YY] = meshgrid(sig,lam);

    [~,I] = unique(f0);
    S0 = interp1(f0(I),S(I),f(f<=30));

    w = 0.125;
    figureNB(13.2,2.8);
    for i = 1:length(A)
        axes('Position',[0.07+(0.91-w-0.07)*(i-1)/4,0.3,0.13,0.6]);
        for j = 1:length(lam)
            for k = 1:length(sig)
                P = PN(A(i),lam(j),1e-3*sig(k));
                M(j,k) = max(10*log10(P(f<=30)./S0(:)));
            end
        end
        % imagesc(sig*1e3,lam,log10(M));
        surf(XX,YY,0*M,M,'LineStyle','none')
        view([0,90]);
        hold on;
        [C,h] = contour(XX,YY,M,[-20,-10,0],'color','k','linewidth',1);
        if(i==4)
            fill([10,90,90,10],[0.11,0.11,4,4],'k','LineStyle','--','FaceColor','none');
        end
        set(gca,'CLim',[-20,20])
        ylim(10.^[-1,2])
        xlim(10.^[0,2])
        set(gca,'xscale','log')
        set(gca,'yscale','log')
        line(1e3*[1e-3,0.1],[100,100],[0,0]-1e6,'color','k','LineWidth',0.75)
        line(1e3*[0.1,0.1],[0.1,100],[0,0]-1e6,'color','k','LineWidth',0.75)

        xticks([1,10,100])
        xticklabels([1,10,100])
        yticks([0.1,1,10,100])
        yticklabels([0.1,1,10,100])
        gcaformat
        xlabel('Jitter (ms)','FontSize',8)
        if(i==1)
            ylabel('Firing rate (Hz)','FontSize',8)
        end
        set(gca,'FontSize',8);
        colormap(flip(bone))
        drawnow;
    end
    C = colorbar;
    C.Position = [0.92,0.3,0.01,0.6];
    C.Label.String = ['Rel. power (dB)'];

end