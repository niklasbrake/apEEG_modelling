[f,Rxx,Rxy] = compute_AP_spectra;
[f0,S] = import_paralytic_data;
% S = fillgaps(nanmean(S(:,1:8),2),5);
low_noise = 1e-3;

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

N = 16e9;

f1 = f; %10.^linspace(0,3,1e3);

GI = 1.2*5.4e-16*10.^AP_model(f0',[20e-3,4e-3,-Inf,0]);
GE = 0.4*5.4e-16*10.^AP_model(f0',[2e-3,1.2e-3,-Inf,0]);

rho_bar = 3.3541e-4;
syn0 = N*(GI+GE)+N*(N-1)*rho_bar*2e-2*GI+N*(N-1)*rho_bar*1e-4*GE;


R = 0.1; lam = 0.1;

pdf_values = @(sigma,mu) (1 ./ (f0 .* sigma * sqrt(2 * pi))) .* exp(-((log(f0) - mu).^2) ./ (2 * sigma.^2));
S_det = nanmean((S(:,1:8)-low_noise)./syn0,2);

% f1 = 10.^linspace(0,3,1e3)';
pdf_values = @(sigma,mu) (1 ./ (f1 .* sigma * sqrt(2 * pi))) .* exp(-((log(f1) - mu).^2) ./ (2 * sigma.^2));
GI = 1.2*5.4e-16*10.^AP_model(f1,[20e-3,4e-3,-Inf,0]);
GE = 0.4*5.4e-16*10.^AP_model(f1,[2e-3,1.2e-3,-Inf,0]);
syn0 = N*(GI+GE)+N*(N-1)*rho_bar*2e-2*GI+N*(N-1)*rho_bar*1e-4*GE;



figureNB(13.2,13.2);
subplot(2,2,1);
    plot(f0,(S(:,1:8)-low_noise),'color','k');
    hold on;
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    for fx = [12.5,23.8,55,100,150,250]
        F = @(x) [x(2).^2-1+2*x(1)+x(2).^2-FT1.predict(log(fx)); log((exp(x(1)-x(2).^2).*x(2)*sqrt(2*pi)).^(-1))-FT2.predict(log(fx))];
        B = fsolve(F,[log(fx),0.18]);
        plot(f1,syn0.*(1+exp(FT3.predict(log(fx)))*pdf_values(B(2),B(1))),'b','LineWidth',1)
    end
    xticks([1,10,100,1e3])
    xticklabels([1,10,100,1e3])
    xlim([1,1e3])
    ylim([1e-5,1e1])
axes();
    plot(f0,S_det,'color','k');
    hold on;
    for fx = [0.8,12.5,23.8,55,100,150,250]
        F = @(x) [x(2).^2-1+2*x(1)+x(2).^2-FT1.predict(log(fx)); log((exp(x(1)-x(2).^2).*x(2)*sqrt(2*pi)).^(-1))-FT2.predict(log(fx))];
        B = fsolve(F,[log(fx),0.1]);
        plot(f1,(1+exp(FT3.predict(log(fx)))*pdf_values(B(2),B(1))),'b','LineWidth',1)
        hold on;
    end
    set(gca,'xscale','log');
    set(gca,'yscale','log');
subplot(2,2,2)
    % plot(f0,S-low_noise,'color','k','Linewidth',1);
    fx = 12.5;
    F = @(x) [x(2).^2-1+2*x(1)+x(2).^2-FT1.predict(log(fx)); log((exp(x(1)-x(2).^2).*x(2)*sqrt(2*pi)).^(-1))-FT2.predict(log(fx))];
    B = fsolve(F,[log(fx),0.1]);
    y11 = syn0.*(1+exp(FT3.predict(log(fx)))*pdf_values(B(2),B(1)));
    plot(f1,y11,'b','LineWidth',1)
    hold on;

    F0 = 12.5; sig = 10e-3*160/F0;
    B = exp(-2*(pi*(f(:)-F0)*sig).^2);
    plot(f,lam*N*Rxx + R*lam*N*(N-1)*B.*Rxy,'color','r','LineWidth',1)
    F0 = 125; sig = 10e-3*160/F0;
    B = exp(-2*(pi*(f(:)-F0)*sig).^2);
    plot(f,lam*N*Rxx + R*lam*N*(N-1)*B.*Rxy,'color','r','LineWidth',1)

    plot(f,lam*N*Rxx+R*lam*N*(N-1)*Rxy,'color',[1,0,0]*0.4+0.6,'LineWidth',1,'LineStyle',':')

    fx = 125;
    F = @(x) [x(2).^2-1+2*x(1)+x(2).^2-FT1.predict(log(fx)); log((exp(x(1)-x(2).^2).*x(2)*sqrt(2*pi)).^(-1))-FT2.predict(log(fx))];
    B = fsolve(F,[log(fx),0.1]);
    y12 = syn0.*(1+exp(FT3.predict(log(fx)))*pdf_values(B(2),B(1)));
    plot(f1,y12,'b','LineWidth',1)
    plot(f1,syn0,'k','LineWidth',1)
    set(gca,'xscale','log')
    set(gca,'yscale','log');
    xlabel('Frequency (Hz)');
    ylabel(['PSD (' char(956) 'V^2/Hz)'])
    gcaformat;
    set(gca,'FontSize',8)

    xticks([1,10,100,1e3])
    xticklabels([1,10,100,1e3])
    xlim([1,1e3])
    ylim([1e-5,1e1])
    % ylim([1e-4,1e2])

subplot(2,2,3);

R = 0.1; lam = 0.1;
y0 = syn0;
F0 = 12.8; sig = 10e-3*160/F0;
B = exp(-2*(pi*(f(:)-F0)*sig).^2);
y2 = lam*N*Rxx + R*lam*N*(N-1)*B.*Rxy;
    fill([f;flip(f)],[y0;flip(y11+y2)],'r');
    hold on;
    fill([f;flip(f)],[y0;flip(y11)],'b');
    set(gca,'xscale','log')
    set(gca,'yscale','log')
F0 = 128; sig = 10e-3*160/F0;
B = exp(-2*(pi*(f(:)-F0)*sig).^2);
y2 = lam*N*Rxx + R*lam*N*(N-1)*B.*Rxy;
    fill([f;flip(f)],[y0;flip(y12+y2)],'r');
    hold on;
    fill([f;flip(f)],[y0;flip(y12)],'b');
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    line(get(gca,'xlim'),[0,0]+low_noise,'color','r','LineWidth',1)

    xlim([1,1e3])
    ylim([1e-5,1e1])

A = 10.^linspace(-4,-1,100);
[XX,YY] = meshgrid(f(f<1e3),A);
lam = 1;
P1 = zeros(length(A),length(f));
for i = 1:length(A)
    P1(i,:) = (A(i)*N*(N-1)*Rxy);
end
P1 = P1(:,f<1e3);
P0x = syn0(f<1e3)'.*exp(FT2.predict(log(f(f<1e3))))';
P0 = P1./(P0x+P1);


% figureNB;
subplot(2,2,4)
%     plot(f(f<1e3),100*P0(67,:),'color','k','Linewidth',1);
%     set(gca,'xscale','log')
%     ylabel('% AP generated')
%     xlabel('Rhythm frequency (Hz)')
%     xlim([1,1e3])
%     xticks([1,10,100,1000])
%     xticklabels([1,10,100,1000])
%     gcaformat(gca,true,8);


% figureNB;
    surf(XX,YY,0*P0-1,100*P0.*((P0x+P1)>1e-3),'LineStyle','none')
    view([0,90]);
    hold on;
    [C,h] = contour(XX,YY,(P0x+P1),[0,1e-3],'color','r','linewidth',2);
    set(gca,'xscale','log')
    set(gca,'yscale','log');
    set(gca,'CLim',[0,100])
    CM = flip(bone(1e3)); CM = CM(1:800,:);
    colormap(CM);
    grid off
    xlabel('Oscillation peak frequency (Hz)')
    ylabel('\lambdaR_{max}')
    ylim([1e-4,0.1])
    yticks([1e-4,1e-3,1e-2,1e-1])
    % yticks([1e-4,1e-2,1e0]);
    xlim([1,1e3])
    xticks([1,10,100,1000])
    xticklabels([1,10,100,1000])
    C = colorbar;
    C.Label.String = '% AP generated';
    % C.Position(1) = 0.89;
    gcaformat(gca,true,8);

    line([1,1],[1e-4,1],'color','k','lineWidth',0.75);
    line([1,1e3],[1e-4,1e-4],'color','k','lineWidth',0.75);