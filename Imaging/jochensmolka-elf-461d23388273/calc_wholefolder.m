% function calc_wholefolder(datafolder)
% %% Asks the user for a data folder, and then calculates for each data folder within this folder:
% % - the mean image
% % - the mean intensity graph
% 
% if nargin < 1, datafolder = 'prompt'; end
   datafolder = ''; 
para             = elf_para(datafolder, '', '', true);    
[~, ~, datasets] = elf_checkdata(para);
res              = zeros(length(datasets), 1);
errors           = cell(length(datasets), 1);

for i = 1:length(datasets)
    fprintf('\n%d of %d\n', i, length(datasets));
    try
        dataset = datasets{i};
        elf_main1_HdrAndInt(dataset); close all
        elf_main2_meanimage(dataset); close all
        elf_main3_intsummary(dataset); close all
        res(i) = 1;
    catch ME
        errors{i} = ME;
        warning(ME.identifier, '%s', ME.message); % change to WARNING
    end
end

%% Print results to 
fprintf('\n\n\n-------------------------------------\n');
fprintf('\tCALC_WHOLEFOLDER\n');
fprintf('--------------------\n');
fprintf('\tResults for data folder %s\n', para.paths.root);
fprintf('\t%d out of %d data sets were successfully calculated.\n', nnz(res), length(res));
fprintf('\tDetailed results were saved to %s\n', para.paths.outputfolder);
fprintf('\tSmall-file results were saved to %s\n', para.paths.outputfolder_pub);
if nnz(res)<length(res)
    fprintf('\t\nThe following errors were encountered:\n');
    errnums = find(res==0)';
    num_maxWidth = floor(log10(max(errnums)));
    name_maxWidth = max(cellfun(@(x) length(x), datasets(errnums)));
    for i = errnums
        formatstr = sprintf('\\t\\tData set %%%dd: %%%ds - %%s, %%s\\n', num_maxWidth, name_maxWidth+3);
        fprintf(formatstr, i, datasets{i}, errors{i}.identifier, errors{i}.message);
    end
end
fprintf('-------------------------------------\n');