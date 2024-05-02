function ai = calculate_dendrite_asymmetry
% ai = CALCULATE_DENDRITE_ASYMMETRY computes the dendrite asymmetry index
%   on all models in the data_files/neuron_models folder.
%   See also generate_morphoData.py

    baseFolder = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    folder = fullfile(baseFolder,'data_files','neuron_models');

    F = dir(folder); F = F(3:end);
    N = zeros(length(F),3);
    asym_idx = zeros(length(F),3);
    for i = 1:length(F)
        file = fullfile(folder,F(i).name,'morphData.mat');
        asym_idx(i,:) = main(file);
    end
    asym_idx = vecnorm(asym_idx,2,2);
end
function asym_idx = main(file)

    load(file)

    % Ensure that the soma is defined to be the origin
    xSoma = data(data(:,2)==1,3:6);
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
    SA = zeros(count,1);
    x0 = zeros(count,3);
    count = 1;
    for i = 1:length(segs)
        idcs = segs{i};
        xSegment = data(idcs,3:6);
        for j = 2:length(idcs)
            h = norm(xSegment(j,1:3)-xSegment(j-1,1:3));
            r1 = xSegment(j,4);
            r2 = xSegment(j-1,4);
            SA(count) = 1/3*pi*h*(r1^2+r1*r2+r2^2);
            x0(count,:) = 0.5*(xSegment(j,1:3)+xSegment(j-1,1:3));
            count=count+1;
        end
    end

    % Calculate asymmetry index
    asym_idx = sum(x0.*SA)./std(x0);
end