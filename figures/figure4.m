function figure4
% FIGURE4 generates the panels in figure 3 of the manuscript.
% See also compute_AP_spectra, import_Scheer2006.

    baseFolder = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(baseFolder,'auxiliary_functions'));
    addpath(fullfile(baseFolder,'modelling'));

    plot_cortex_schematic;
    plot_example_jitter;

    [f0,S] = import_Scheer2006;
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

function plot_cortex_schematic
    fig = figure('color','w','units','centimeters');
    fig.Position(3:4) = [8.5,6];

    % Get surface area of each triangle
    [sa,X] = import_leadfield;
    x0 = [54,24,35];
    [~,i0] = min(vecnorm(X.vertices-x0,2,2));
    j0 = find(sum(X.faces==i0,2));
    c0 = mean(X.vertices(X.faces(j0(1),:),:));
    normal = -sa.cortex75K.normals(i0,:);
    [x1,x2] = getOrthBasis(normal);
    normal = -[0.8529,0.4924,0.1736];
    x1 = [0,0,1];
    x2 = [-0.4924,0.8529,0];

    % Change coordinate system
    V = zeros(size(X.vertices));
    V(:,1) = sum((X.vertices-c0).*x2,2);
    V(:,2) = sum((X.vertices-c0).*normal,2);
    V(:,3) = sum((X.vertices-c0).*x1,2);
    idcs = zeros(size(X.faces,1),1);
    % Remove faces far away
    for i = 1:size(X.faces);
        tri = V(X.faces(i,:),:);
        if(all(vecnorm(tri,Inf,2)<30))
            idcs(i) = 1;
        end
    end
    F = X.faces(find(idcs),:);

    % Subdivide triangles to allow smoooth colouring
    sigma = 4;
    pairs = nchoosek(1:3,2);
    V0 = V;
    F0 = F;
    V2 = V;
    F2 = [];
    newSub = true;
    k = 1;
    while(newSub)
        newSub = false;
        disp(['round ' int2str(k)])
        for i = 1:size(F0,1)
            tri = V0(F0(i,:),:);
            v_dist = exp(-vecnorm(tri,2,2).^2/sigma);
            maxColourChange = max(abs(diff(v_dist(pairs),1,2)));
            if(maxColourChange>0.1)
                [v1,f1] = subdivide_tri(V2,F0(i,:));
                V2 = [V2;v1];
                newSub = true;
            else
                f1 = F0(i,:);
            end
            F2 = [F2;f1];
        end
        V0 = V2;
        F0 = F2;
        F2 = [];
        k = k+1;
    end

    X2.vertices = V;
    X2.faces = F;
    X3.vertices = V0;
    X3.faces = F0;

    sequential_CM = [255,255,255; ...
                    255,247,188; ...
                    254,227,145; ...
                    254,196,79; ...
                    251,154,41; ...
                    236,112,20; ...
                    204,76,2; ...
                    153,52,3; ...
                    102,37,6]/255;

    sigma = 6;
    ax(3) = axes('Position',[0.11,0.515,0.4,0.5]);
        d = vecnorm(X.vertices-c0,2,2);
        C = exp(-d.^2/sigma); C = C/max(C);
        trisurf(X.faces, X.vertices(:,1), X.vertices(:,2), X.vertices(:,3),...
        'FaceLighting','gouraud','FaceVertexCData',C,'EdgeColor','none','FaceColor','interp');
        hold on
        view([120,10]);
        axis tight equal off
        camlight headlight
        material dull

        VW = get(gca,'view');
        [a1,a2,a3] = sph2cart((VW(1)-90)*pi/180,VW(2)/180*pi,1);
        e0 = [a1,a2,a3]';
        e2 = camup';
        e1 = cross(e0,-e2);
        R = 40;
        H = 15;
        x1 = c0(:) + sqrt(2)*R/2*(e1*cos(-3*pi/4)+e2*sin(-3*pi/4))+H*e0;
        x2 = x1(:) + R*(e1*cos(pi/2)+e2*sin(pi/2));
        x3 = x2(:) + R*(e1*cos(0)+e2*sin(0));
        x4 = x3(:) + R*(e1*cos(-pi/2)+e2*sin(-pi/2));
        x = [x1,x2,x3,x4];
        P = patch(x(1,:),x(2,:),x(3,:),'b','LineWidth',1,'FaceColor','none');
        set(ax(3),'CLim',[0,1]);
        colormap(ax(3),interp1(linspace(0,1,9),sequential_CM,linspace(0,1,100)));
    ax(4) = axes('Position',[0.455,0.64,0.3,0.3]);
        d = vecnorm(X3.vertices,2,2);
        C = exp(-d.^2/sigma);  C = C/max(C);
        trisurf(X3.faces, X3.vertices(:,1), X3.vertices(:,2), X3.vertices(:,3),...
        'FaceLighting','gouraud','FaceVertexCData',C,'EdgeColor','none','FaceColor','interp');
        hold on
        view([0,0])
        axis tight equal off
        camlight headlight
        material dull

        colormap(ax(4),interp1(linspace(0,1,9),sequential_CM,linspace(0,1,100)));
        set(gca,'CLim',[0,1]);
        yl = get(gca,'ylim');
        xlim([-20,20]);
        zlim([-20,20]);
        x = [-20,-20,20,20];
        y = [1,1,1,1]*yl(1);
        z = [-20,20,20,-20];
        P = patch(x,y,z,'b','LineWidth',1,'FaceColor','none');
        CB = colorbar;
        CB.Position = [0.74,ax(4).Position(2),0.02,ax(4).Position(4)];
        CB.Ticks = [0,1];
        CB.TickLabels = {'0','\rho_{max}'};
        title('Spiking synchrony','fontweight','normal','fontsize',7);
        annotation('line','Position',[0.37 0.82+0.02 0.12 0.10]);
        annotation('line','Position',[0.37 0.715+0.02 0.12 -0.095]);
    ax(5) = axes('Position',[0.48,0.575,0.24,0.04]);
        xlim([-20,20]);
        ylim([0,1]);
        text(-10,0.2,'10 mm','FontSize',7,'HorizontalAlignment','center','VerticalAlignment','bottom');
        line([-15,-5],[0,0],'color','k','LineWidth',2);
        axis off;
end

function plot_example_jitter(apResponses)
    [sa,X] = import_leadfield;
    dt = 1/16e3;
    tmax = 1;
    N = tmax/dt+1;
    t = 0:dt:tmax;

    m1 = poissrnd(tmax*100);

    t0 = tmax*rand(m1,1);

    A = 0.2;
    S = 1e-3;

    m2 = floor(A*m1);
    s1 = t0(randperm(m1,m2));%+S*randn(m2,1);
    s2 = t0(randperm(m1,m2));%+S*randn(m2,1);

    [common,I1,I2] = intersect(s1,s2);
    s1 = s1+S*randn(m2,1);
    s2 = s2+S*randn(m2,1);

    if(isempty(common))
        error('By chance, the example spike trains have no overlapping spikes. Try running again.');
    end

    % figureNB(8,6);
    axes('Position',[0.05,0.08,0.6,0.4]);
        R = raster([s1*0;s2*0+1],[s1;s2],gcf);
        R.Color = [0.6,0.6,0.6];
        hold on;
        R = raster(common*0,s1(I1),gcf);
        R.Color = 'k';
        R = raster(common*0+1,s2(I2),gcf);
        R.Color = 'k';
        xlim([0,tmax]);
        ylim([-0.2,2.1]);
        axis off;
        plot(common(2)+1e-3*[-30,30,30,-30,-30],[-0.2,-0.2,2.1,2.1,-0.2],':k','LineWidth',1)

    axes('Position',[0.7,0.08,0.25,0.4]);
    c = linspace(-5,5,1e3);
        fill([c,flip(c)],[0*c,0.9*exp(-c.^2./(2*(1e3*S)^2))],[0.8,0.8,0.8],'EdgeColor','none');
        hold on;
        % plot(c,exp(-c.^2./(2*8^2)),'color','r','LineWidth',1);
        line([0,0],[0,0.9],'color',0*[0.6,0.6,0.6],'LineWidth',1)
        del = 1e3*(s1(I1(2))-common(2));
        line([0,0]+del,[0,0.9],'color','k','LineWidth',1)

        fill([c,flip(c)],[1+0*c,1+0.9*exp(-c.^2./(2*(1e3*S)^2))],[0.8,0.8,0.8],'EdgeColor','none');
        hold on;
        % plot(c,1.2+exp(-c.^2./(2*8^2)),'color','r','LineWidth',1);
        line([0,0],1+[0,0.9],'color',0*[0.6,0.6,0.6],'LineWidth',1)
        del = 1e3*(s2(I2(2))-common(2));
        line([0,0]+del,1+[0,0.9],'color','k','LineWidth',1)

        plot([-5,5,5,-5,-5],[-0.2,-0.2,2.1,2.1,-0.2],':k','LineWidth',1)
        axis off;
end