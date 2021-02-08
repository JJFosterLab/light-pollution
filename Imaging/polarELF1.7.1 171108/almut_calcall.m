%%
elf_paths
clear;
close all;

para = elf_para('F:\Newly recalced');
[~, ~, datasets] = elf_checkdata(para);

% for i = 1:length(datasets)
%     fprintf('\n%d of %d\n', i, length(datasets));
%     dataset = datasets{i};
%     elf_main1_HdrAndInt(dataset); close all
%     elf_main2_meanimage(dataset); close all
%     elf_main3_intsummary(dataset); close all
% end

% % filter images
% for i = 1:length(datasets)
%         dataset = datasets{i};
%         elf_main4_filter(dataset); close all
% end

% calculate contrasts
for i = 4:length(datasets)
        dataset = datasets{i};
        elf_main5_contrasts(dataset); close all
end

% combine data and plot averages
for i = 3:length(datasets)
        dataset = datasets{i};
        elf_main6_summary(dataset); close all
        elf_main7_stats_and_plots(dataset);
end

% combine data from both cameras
elf_postana_combine(datasets(1:2), 'm1x - All occupied sunspots - anterior view')
elf_postana_combine(datasets(3:4), 'm2x - All occupied sunspots - posterior view')
elf_postana_combine(datasets(5:6), 'm3x - All unoccupied sunspots - anterior view')
elf_postana_combine(datasets(7:8), 'm4x - All unoccupied sunspots - posterior view')


