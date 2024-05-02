function [sa,X] = import_leadfield
    baseFolder = fileparts(fileparts(mfilename('fullpath')));
    load(fullfile(baseFolder,'data_files','anatomy_nyhead_model.mat'));
    X = struct();
    X.vertices = sa.cortex75K.vc;
    X.faces= sa.cortex75K.tri;
end