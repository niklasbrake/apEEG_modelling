[f0,S] = import_paralytic_data;
f0 = f0(:);

rho_bar = 3.3541e-4;
N = 16e9;

GI = 1.2*5.4e-16*10.^AP_model(f0,[20e-3,4e-3,-Inf,0]);
GE = 0.4*5.4e-16*10.^AP_model(f0,[2e-3,1.2e-3,-Inf,0]);

syn0 = N*(GI+GE)+N*(N-1)*rho_bar*2e-2*GI+N*(N-1)*rho_bar*1e-4*GE;

% Compute the log-normal PDF using the equation
f1 = 10.^linspace(0,3,1e3)';
pdf_values = @(sigma,mu) (1 ./ (f0 .* sigma * sqrt(2 * pi))) .* exp(-((log(f0) - mu).^2) ./ (2 * sigma.^2));

mu = [0.1,log(12.5),log(24),log(54)];
sig = [0.6,0.1,0.09,0.18];
A = [25,19,17,30];


mu(1) = [];
sig(1) = [];
A(1) = [];


figureNB;

% subplot(2,1,1);
%     plot(f0,S(:,1:8)-low_noise,'color','k');
%     hold on;
%     plot(f0,syn0.*(1+sum(A.*pdf_values(sig,mu),2)),'b')
%     set(gca,'xscale','log');
%     set(gca,'yscale','log');
S_det = nanmean((S(:,1:8)-low_noise)./syn0,2);
subplot(2,1,1);
    plot(f0,(S(:,1:8)-low_noise)./syn0,'color','k');
    hold on;
    plot(f0,(1+sum(A.*pdf_values(sig,mu),2)),'b','LineWidth',1)
    plot(f0,f0*0+1,'--b')
    set(gca,'xscale','log');
    set(gca,'yscale','log');


subplot(2,3,4);
    plot(exp(mu-sig.^2),exp(sig.^2-1).*exp(2*mu+sig.^2),'.-k','LineWidth',1,'MarkerSize',10);
    hold on;
    FT1 = fitlm(mu-sig.^2,sig.^2-1+2*mu+sig.^2);
    plot(f1,exp(FT1.predict(log(f1))),'b');
    xlabel('Peak frequency (Hz)');
    ylabel('Peak variance');
    set(gca,'xscale','log');
    set(gca,'yscale','log');

subplot(2,3,5);
    plot(exp(mu-sig.^2),(exp(mu-sig.^2).*sig*sqrt(2*pi)).^(-1),'.-k','LineWidth',1,'MarkerSize',10);
    hold on;
    FT2 = fitlm(mu-sig.^2,log((exp(mu-sig.^2).*sig*sqrt(2*pi)).^(-1)));
    plot(f1,exp(FT2.predict(log(f1))),'b');
    xlabel('Peak frequency (Hz)');
    ylabel('Peak amplitude');
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    % return;

subplot(2,3,6);
    plot(exp(mu-sig.^2),A,'.-k','LineWidth',1,'MarkerSize',10);
    hold on;
    FT3 = fitlm(mu-sig.^2,log(A));
    plot(f1,exp(FT3.predict(log(f1))),'b');
    xlabel('Peak frequency (Hz)');
    ylabel('Peak area');
    set(gca,'xscale','log');
    set(gca,'yscale','log');


f1 = 10.^linspace(0,3,1e3)';
pdf_values = @(sigma,mu) (1 ./ (f1 .* sigma * sqrt(2 * pi))) .* exp(-((log(f1) - mu).^2) ./ (2 * sigma.^2));

GI = 1.2*5.4e-16*10.^AP_model(f1,[20e-3,4e-3,-Inf,0]);
GE = 0.4*5.4e-16*10.^AP_model(f1,[2e-3,1.2e-3,-Inf,0]);

syn0 = N*(GI+GE)+N*(N-1)*rho_bar*2e-2*GI+N*(N-1)*rho_bar*1e-4*GE;


figureNB;
subplot(1,2,1);
    plot(f0,(S(:,1:8)-low_noise),'color','k');
    hold on;
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    for fx = [12.5,23.8,55,100,150,250]
        F = @(x) [x(2).^2-1+2*x(1)+x(2).^2-FT1.predict(log(fx)); log((exp(x(1)-x(2).^2).*x(2)*sqrt(2*pi)).^(-1))-FT2.predict(log(fx))];
        B = fsolve(F,[log(fx),0.18]);
        plot(f1,syn0.*(1+exp(FT3.predict(log(fx)))*pdf_values(B(2),B(1))),'b','LineWidth',1)
    end
subplot(1,2,2);
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