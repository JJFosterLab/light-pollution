function elf_main7_stats_and_plots(dataset, imgformat, plotall)
% elf_MAIN7_STATS_AND_PLOTS reads the previously calculated results for a dataset, writes the results into an 
% Excel spreadsheet, plots the results, and saves the plot to pdf and jpg 
%
% Uses: elf_support_logmsg, elf_paths, elf_para, elf_para_update, 
%       elf_readwrite, elf_analysis_writestats, elf_support_formatA4l,
%       elf_plot_summary_ls
%
% Loads files: mean results file in mat folder, mean image in detailed results folder
% Saves files: XLSX file and PDF in Detailed results folder, JPG in public
%              results folder
% 
% Typical timing for a 50-scene environment (on ELFPC):
%       5s total

%% check inputs
if nargin < 3, plotall = true; end
if nargin < 2 || isempty(imgformat), imgformat = '*.dng'; end
if nargin < 1 || isempty(dataset), error('You have to provide a valid dataset name'); end 

                    elf_support_logmsg('\b\b\b\b\b\b\b\b\b\b\b\b\b\n');
                    elf_support_logmsg('----- ELF Step 7: Writing stats and plotting summary -----\n');
                    elf_support_logmsg('      Displaying environment %s\n', dataset);

%% House-keeping and initialisations
linims          = strcmp(imgformat, '*.dng');

%% Set up paths and file names; read info, infosum and para
elf_paths;
para            = elf_para('', dataset, imgformat);
para            = elf_para_update(para);                   % Combine old parameter file with potentially changed information in current elf_para
infosum         = elf_readwrite(para, 'loadinfosum');      % loads the old infosum file (which contains projection information)

%% Load data mean and mean image
datamean        = elf_readwrite(para, 'loadmeanres');
meanim          = elf_readwrite(para, 'loadmeanimg_tif');

%% Write stats into Excel file
elf_analysis_writestats(datamean, para);

%% Plot final ELF plot and save to pdf and jpg
if plotall
    fh = elf_support_formatA4l(7);
    set(fh, 'Name', 'Visual Environment Plot Standard');
    elf_plot_summary_ls(para, datamean, meanim, infosum, fh, para.paths.dataset, linims) % generate summary ELF plot

    elf_readwrite(para, 'savemeanvep_pdf', '', fh)
    elf_readwrite(para, 'savemeanvep_jpg', '', fh)
    elf_support_logmsg('\n');
end