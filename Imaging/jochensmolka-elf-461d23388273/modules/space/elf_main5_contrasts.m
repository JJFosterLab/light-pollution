function elf_main5_contrasts(dataset, imgformat, verbose)
% elf_MAIN5_CONTRASTS calculates all spatial descriptors (contrasts) for
% an environment.
%
% Uses: elf_support_logmsg, elf_paths, elf_para, elf_para_update, 
%       elf_info_collect, elf_hdr_brackets, elf_readwrite, 
%       space_analysis, elf_support_formatA4l, elf_plot_summary_ls
%
% Loads files: filtered image files in filt folder, info files
% Saves files: mean results file in mat folder
% 
% Typical timing PER SCENE (on ELFPC):
%      36s for spatial analysis
%       2s to save jpg
%      <1s rest

%% check inputs
if nargin < 3, verbose = false; end
if nargin < 2 || isempty(imgformat), imgformat = '*.dng'; end
if nargin < 1 || isempty(dataset), error('You have to provide a valid dataset name'); end 

                    elf_support_logmsg('\b\b\b\b\b\b\b\b\b\b\b\b\b\n');
                    elf_support_logmsg('----- ELF Step 5: Contrast descriptors -----\n');

%% House-keeping and initialisations
savejpgs  = true;                                                        % save individual jpgs for each image? (takes extra time)
calib     = strcmp(imgformat, '*.dng');                                  % jpgs and tifs cant be calibrated

%% Set up paths and file names; read info, infosum and para
elf_paths;
para      = elf_para('', dataset, imgformat);
para      = elf_para_update(para);                                      % Combine old parameter file with potentially changed information in current elf_para
info      = elf_info_collect(para.paths.datapath, imgformat);           % this contains exif information and filenames
sets      = elf_hdr_brackets(info);                                     % determine which images are part of the same scene
infosum   = elf_readwrite(para, 'loadinfosum');                         % loads the old infosum file (which contains projection information)

                    elf_support_logmsg('      Processing %d scenes in environment %s\n', size(sets, 1), dataset);

%% Process one scene at a time
for setnr = 1:size(sets, 1)
    % Load intensity descriptors, filtered HDR scenes and HDR images
    scenename       = sprintf('scene%03d', setnr);
    res             = elf_readwrite(para, 'loadres', {scenename});
    im_filt_HDR     = elf_readwrite(para, 'loadfilt_mat', sprintf('scene%03d', setnr));
    sceneim         = elf_readwrite(para, 'loadHDR_tif', scenename); % only for plotting
    
    % Calculate spatial descriptors 
    res.spatial     = space_analysis(para, im_filt_HDR, verbose);

    % Plot summary figure for set
    fh              = elf_support_formatA4l(5); clf;
    set(fh, 'Name', sprintf('Scene #%d of %d', setnr, size(sets, 1)));
    elf_plot_summary_ls(para, res, sceneim, infosum, fh, sprintf('%s, scene #%d of %d', para.paths.dataset, setnr, size(sets, 1)), calib);


    %% Save output to jpg and mat
    elf_readwrite(para, 'saveres', scenename, res);
    if savejpgs, elf_readwrite(para, 'savevep_jpg', scenename, fh); end    % small bottleneck
end
    
                    elf_support_logmsg('      Summary: All contrast descriptors for environment %s calculated and saved to jpg and mat.\n\n', para.paths.dataset);

