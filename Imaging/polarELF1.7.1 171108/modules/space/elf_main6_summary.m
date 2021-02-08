function elf_main6_summary(dataset, imgformat)
% elf_MAIN6_SUMMARY calculates the mean result for a dataset, and saves as MAT 
%
% Uses: elf_support_logmsg, elf_paths, elf_para, elf_para_update, 
%       elf_info_collect, elf_readwrite, elf_analysis_datasetmean
%
% Loads files: individual results files in mat folder, info files
% Saves files: mean results file in mat folder
% 
% Typical timing for a 50-scene environment (on ELFPC):
%       1s total

%% check inputs
if nargin < 2 || isempty(imgformat), imgformat = '*.dng'; end
if nargin < 1 || isempty(dataset), error('You have to provide a valid dataset name'); end 

                    elf_support_logmsg('\b\b\b\b\b\b\b\b\b\b\b\b\b\n');
                    elf_support_logmsg('----- ELF Step 6: Calculating summary -----\n');

%% Set up paths and file names; read info, infosum and para
elf_paths;
para            = elf_para('', dataset, imgformat);
para            = elf_para_update(para);                                                             % Combine old parameter file with potentially changed information in current elf_para
info            = elf_info_collect(fullfile(para.paths.datapath, para.paths.scenefolder), '*.tif');  % this contains tif exif information and filenames %%FIXME should be mat folder
fnames_im       = {info.Filename};                                                                    % collect image names

                    elf_support_logmsg('      Averaging across all %d scenes in environment %s\n', length(fnames_im), dataset);

%% Load data, calculate data mean, save data mean
data            = elf_readwrite(para, 'loadres', fnames_im);
datamean        = elf_analysis_datasetmean(data, [], [], para.plot.datasetmeantype);
                  elf_readwrite(para, 'savemeanres', '', datamean); % write data mean
                  elf_readwrite(para, 'savecollres', '', data);     % write collated data to results folder
                  

                    elf_support_logmsg('      Means calculated and saved to MAT.\n\n', dataset);