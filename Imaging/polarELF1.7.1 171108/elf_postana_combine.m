function elf_postana_combine(datasets, outfilename, imgformat)
% elf_postana_combine calculates ...
%
% datasets - cell array of dataset names
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
if nargin < 3 || isempty(imgformat), imgformat = '*.dng'; end
if nargin < 2 || isempty(outfilename)
    temp            = inputdlg('Please enter a name for the combined data set:');
    outfilename     = temp{1};
end
if nargin < 1 || isempty(datasets), error('You have to provide at least one valid dataset name'); end 

                    elf_support_logmsg('\b\b\b\b\b\b\b\b\b\b\b\b\b\n');
                    elf_support_logmsg('----- ELF Postanalysis: Combining data sets -----\n');

%% House-keeping and initialisations
linims          = strcmp(imgformat, '*.dng');
elf_paths;

%% Calculate data and image mean, save means
filecount   = 0;
envcount    = 0;
data        = [];
fprintf('Out of %d: ', length(datasets));
for i = 1:length(datasets)
    %try
        para            = elf_para('', datasets{i}, imgformat);
        para            = elf_para_update(para);          %NECESSARY?                                                      % Combine old parameter file with potentially changed information in current elf_para

        if filecount == 0        
            infosum     = elf_readwrite(para, 'loadinfosum');                                                              % loads the old infosum file (which contains projection information, and linims)
            sumimage    = zeros(length(infosum.proj_ele), length(infosum.proj_azi), infosum.SamplesPerPixel, 'double');     % pre-allocate for sum of all processed images
            combpara    = para;
            combpara.paths.dataset  = outfilename;
            combpara.paths.datapath = fullfile(combpara.paths.root, combpara.paths.dataset);
            combpara    = elf_readwrite(combpara, 'createfilenames');
            combinfosum = infosum; %%FIXME: Only a placeholder really
            elf_readwrite(combpara, 'saveinfosum', '', combinfosum);
        end
        
        % data
        temp            = elf_readwrite(para, 'loadcollres');
        data            = [data temp]; %#ok<AGROW>
        
        % mean image
        sumimage        = sumimage + double(elf_readwrite(para, 'loadmeanimg_tif'));

        % increment
        filecount       = filecount + length(temp);
        envcount        = envcount  + 1;
        
%     catch me
%         warning('DATASET aborted: %s', datasets{i});
%         warning(me.message);
%     end
    
    if mod(i, 10)==0
        fprintf('%d.', i);
    end
end
fprintf('\n');
                    elf_support_logmsg('      Averaging results across %d scenes in %d environments\n', filecount, envcount);

datamean            = elf_analysis_datasetmean(data, [], [], para.plot.datasetmeantype);
meanimage           = sumimage/envcount;

elf_readwrite(combpara, 'savemeanres', '', datamean);
elf_readwrite(combpara, 'savecollres', '', data);
elf_readwrite(combpara, 'savemeanimg_tif', '', uint16(meanimage));
elf_readwrite(combpara, 'savemeanimg_jpg', '', uint16(meanimage));

                    elf_support_logmsg('      Created mean image and averaged data across %d scenes in %d environments\n', filecount, envcount);

%% Write stats into Excel file
elf_analysis_writestats(datamean, combpara);

%% Plot final ELF plot and save to pdf and jpg
fh              = elf_support_formatA4l(77); clf;
                  set(fh, 'Name', 'Environmental Light Field');
                  elf_plot_summary_ls(combpara, datamean, uint16(meanimage), combinfosum, fh, combpara.paths.dataset, linims) % generate summary ELF plot

elf_readwrite(combpara, 'savemeanvep_pdf', '', fh)
elf_readwrite(combpara, 'savemeanvep_jpg', '', fh)
elf_support_logmsg('\n');

fh              = elf_support_formatA4l(77); clf;
                  set(fh, 'Name', 'Environmental Light Field');
                  elf_plot_summary_int(combpara, datamean, uint16(meanimage), combinfosum, fh, combpara.paths.dataset, sprintf('n = %d scenes from %d environments', filecount, envcount)); % generate summary ELF plot
                  
fh = elf_support_formatA4l(77);
elf_readwrite(combpara, 'savemeanivep_jpg', '', fh);
elf_readwrite(combpara, 'savemeanivep_pdf', '', fh);
