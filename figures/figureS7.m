load('E:\Research_Projects\005_Aperiodic_EEG\unitary_APs\data\simulations\bAP_unitary_response\neuron_models\L6_BPC_cADpyr231_3\morphData.mat');
connections = connections+1;

temp = data(:,5);
data(:,5) = data(:,4);
data(:,4) = temp;
data(:,3:5) = data(:,3:5) - mean(xSoma(:,1:3));

% Count total number of segments
count = 0;
for i = 1:length(segs)
    idcs = segs{i};
    for j = 2:length(idcs)
        count = count+1;
    end
end

% Compute surface area and coordinate of every segment
SA = zeros(length(segs),1);
x0 = zeros(length(segs),3);
soma = zeros(length(segs),1);
count = 1;
for i = 1:length(segs)
    idcs = segs{i};
    xSegment = data(idcs,3:6);
    if(any(data(idcs,2)==1))
        soma(i) = 1;
    end
    count = 1;
    SA2 = [];
    x02 = [];
    for j = 2:length(idcs)
        h = norm(xSegment(j,1:3)-xSegment(j-1,1:3));
        r1 = xSegment(j,4);
        r2 = xSegment(j-1,4);
        SA2(count) = 1/3*pi*h*(r1^2+r1*r2+r2^2);
        x02(count,:) = 0.5*(xSegment(j,1:3)+xSegment(j-1,1:3));
        count=count+1;
    end
    SA(i) = sum(SA2);
    x0(i,:) = mean(x02);
end

% Calculate asymmetry index
asym_idx = sum(x0.*SA)./std(x0)*0.1;
x1 = mean(x0.*SA);
SD = std(x0);

idcs = find(and(SA>0,~soma));

fig = figureNB(12,18);
axes('Position',[0,0,1,1]);
scatter(x0(idcs,1),x0(idcs,3),SA(idcs),SA(idcs)*0,'filled')
hold on
% view([0,0]);
set(gca,'DataAspectRatio',[1,1,1]);
axis off;
set(gca,'CLim',[-1,1]);
colormap([1,0,0;0,0,0;0,0,1]);

line([-SD(1),SD(1)]/2,[-280,-280],'color','k','LineWidth',3);
line([275,275],[-SD(3),SD(3)]/2,'color','k','LineWidth',3);




fig = figureNB(12,18);
axes('Position',[0,0,1,1]);
scatter(x0(idcs,1),x0(idcs,3),SA(idcs)/2,SA(idcs)*0,'filled')
hold on
% view([0,0]);
set(gca,'DataAspectRatio',[1,1,1]);
axis off;
set(gca,'CLim',[-1,1]);
colormap([1,0,0;0.6,0.6,0.6;0,0,1]);

line([-SD(1),SD(1)]/2,[-280,-280],'color','k','LineWidth',3);
line([275,275],[-SD(3),SD(3)]/2,'color','k','LineWidth',3);
line([0,asym_idx(1)],[0,asym_idx(3)],'color',[0,0,0],'LineWidth',3)
scatter(asym_idx(1),asym_idx(3),100,[0,0,0],'filled')

% figureNB;


% xlim([130,190])
% ylim([-162.527,-135.391])