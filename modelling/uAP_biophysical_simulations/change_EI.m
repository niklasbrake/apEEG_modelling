folder = '/lustre04/scratch/nbrake/data/simulations/unitary_AP';

load('/lustre04/scratch/nbrake/code/apEEG/newEI.mat')

for i = 1:length(newEI)
    csvwrite(fullfile(folder,filenames{i},'EI_ratio.csv'),newEI(i));

end

fid = fopen('C:\Users\brake\Desktop\changedEI.txt','w');
fprintf(fid,'%s\n',filenames{:})
fclose(fid)