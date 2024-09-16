[sa,X] = network_simulation_beluga.getHeadModel;
choose_files = {'L6_BPC_cADpyr231_3','L23_PC_cADpyr229_2','L4_SBC_cACint209_5','L4_ChC_cNAC187_1','L5_TTPC1_cADpyr232_1'}

for i = 1:length(choose_files)
    morphData = fullfile('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\neuron_models',choose_files{i},'morphData.mat');
    plot_neuron_morphology(morphData);
    print(gcf,['C:\Users\brake\Desktop\morph' int2str(i) '.png'],'-dpng','-r600')
end



pos = [70,-25,-20;
        50,-77,25;
        68,-8,23;
        27,-31,76;
        41,60,-11;
        10,41,52];

for i = 1:6
    [~,locs(i)] = min(vecnorm(sa.cortex75K.vc-pos(i,:),2,2));
end

[~,I] = sort(vecnorm(pos-sa.locs_3D(49,1:3),2,2));
locs = locs(I);



figureNB(18,18);
for i = 1:length(choose_files)
    j = find(strcmp(choose_files{i},files))
    subplot(7,5,i);
    plot(savedUnitaryAP(:,:,j),'LineWidth',1);
    xlim([1950,2100])
    ylim([-200,360]);
    axis off
    for k = 1:6
        eeg = network_simulation_beluga.getEEG(savedUnitaryAP(:,:,j),sa,locs(k));
        subplot(7,5,i+5*k)
        plot(eeg,'k','LineWidth',1);
        xlim([1950,2100])
        ylim([-40,30]*1e-6)
        axis off
    end
end

figureNB;
trisurf(X.faces, X.vertices(:,1), X.vertices(:,2), X.vertices(:,3),...
'FaceLighting','gouraud','FaceVertexCData',X.vertices*0+0.7,'EdgeColor','none','FaceColor','interp');
hold on
view([120,10]);
camlight headlight
material dull
% axis tight;
axis equal off;
hold on

P = X.vertices(locs,:);
N = sa.cortex75K.normals(locs,:);

plot3(P(:,1),P(:,2),P(:,3),'.','color',blue,'MarkerSize',20);
plot3(P(:,1)'+20*[zeros(6,1),N(:,1)]',P(:,2)'+20*[zeros(6,1),N(:,2)]',P(:,3)'+20*[zeros(6,1),N(:,3)]','-','color',blue,'MarkerSize',20,'LineWidth',1);

plot3(sa.locs_3D(49,1),sa.locs_3D(49,2),sa.locs_3D(49,3),'.k','MarkerSize',20);


for i = 1:5
    figureNB;
    trisurf(X.faces, X.vertices(:,1), X.vertices(:,2), X.vertices(:,3),...
    'FaceLighting','gouraud','FaceVertexCData',X.vertices*0+0.7,'EdgeColor','none','FaceColor','interp');
    hold on
    view([120,10]);
    camlight headlight
    material dull
    % axis tight;
    axis equal off;
    hold on
    plot3(X.vertices(locs(i),1),X.vertices(locs(i),2),X.vertices(locs(i),3),'.k','MarkerSize',100);
end