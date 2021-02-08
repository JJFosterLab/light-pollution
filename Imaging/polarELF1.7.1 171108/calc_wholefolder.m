%% Asks the user for a data folder, and then calculates for each data folder within this folder:
% - the mean image
% - the mean intensity graph

clear
para = elf_para(NaN);
[~, ~, datasets] = elf_checkdata(para);
% 1. calculate mean image
for i = 1:length(datasets)
    dataset = datasets{i};
    elf_main1_HDRscenes(dataset); close all
    elf_main2_meanimage(dataset); close all
end
% 2. calculate mean intensity descriptors
for i = 1:length(datasets)
    dataset = datasets{i};
    elf_main3_luminance(dataset); close all
end
