function elf_main4_filter(dataset, imgformat, verbose)
% elf_MAIN4_FILTER filters all HDR scenes for an environment for a 1
% degree and 10 degree resolution. The filtering algorithm (almost) fully corrects
% for the distortion introduced by the equirectangular projection.
%
% Uses: elf_support_logmsg, elf_paths, elf_para, elf_para_update, 
%       elf_readwrite, elf_support_formatA4, elf_filter 
%
% Loads files: HDR scenes as mat in scene folder
% Saves files: filtered images as mat in filt folder
% 
% Typical timing PER SCENE (on ELFPC):
%     192s for filtering both resolutions
%     2.5s to save jpg
%     1.5s rest
%
%   + 192s PER ENVIRONMENT to calculate initial correction images

%% check inputs
if nargin < 3, verbose = true; end
if nargin < 2 || isempty(imgformat), imgformat = '*.dng'; end
if nargin < 1 || isempty(dataset), error('You have to provide a valid dataset name'); end 

                    elf_support_logmsg('\b\b\b\b\b\b\b\b\b\b\b\b\b\n');
                    elf_support_logmsg('----- ELF Step 4: Filtering scenes -----\n');

%% Set up paths and file names; read info, infosum and para
elf_paths;
para        = elf_para('', dataset, imgformat);
para        = elf_para_update(para);                                       % Combine old parameter file with potentially changed information in current elf_para
allfiles    = elf_io_dir(fullfile(para.paths.datapath, para.paths.scenefolder, '*.mat'));
fnames_im   = {allfiles.name};                                              % collect image names
infosum     = elf_readwrite(para, 'loadinfosum');                          % loads the old infosum file (which contains projection information)

                    elf_support_logmsg('      Processing %d scenes in environment %s\n', length(fnames_im), dataset);
                    
%% Process one scene at a time
it = 1;
for setnr = 1:length(fnames_im)
    %% Load HDR image
    im_HDR  = elf_readwrite(para, 'loadHDR_mat', fnames_im{setnr});
    
    %% Filter HDR image
    % panel 1: original image
    fh      = elf_support_formatA4(4);
    p1      = uipanel('Parent', fh, 'Position', [0 1/3 1 2/3]);
              elf_plot_image(im_HDR, infosum, p1, 'equirectangular', 1);
    ah1     = axes('Parent', p1, 'units', 'normalized', 'position', [0 0 .05 1]);
              text(0.5, 0.5, 'HDR image (delinearised for display)', 'units', 'normalized', 'position', [0 0.5], 'rotation', 90, 'horizontalalignment', 'center', 'verticalalignment', 'middle', 'fontweight', 'bold');
              axis(ah1, [0 1 0 1], 'off');
    
    % panel 2: 1 deg filtered image
    p2      = uipanel('Parent', fh, 'Position', [0 0 .5 1/3]);
    ah(1)   = axes('Parent', p2, 'Position', [0 0 1 1]);
    
    % panel 3: 10 deg filtered image
    p3      = uipanel('Parent', fh, 'Position', [0.5 0 .5 1/3]);
    ah(2)   = axes('Parent', p3, 'Position', [0 0 1 1]); axis(ah, 'off');
            set(fh, 'Name', sprintf('Scene #%d of %d', setnr, length(allfiles)));
            drawnow;
    
    %% Filter images
    if it==1
        % on first iteration, calculate correction image for borders, and save projection information (%TODO)
        [im_filt_HDR, para] = elf_filter(para, im_HDR, 'scene', verbose, ah, infosum, 1);
        % elf_readwrite(para, 'saveinfosum', [], infosum); % saves infosum AND para
    else
        im_filt_HDR = elf_filter(para, im_HDR, 'scene', verbose, ah, infosum);
    end
    it      = it+1;
    
    %% save projected and filtered images
    elf_readwrite(para, 'savefilt_mat', sprintf('scene%03d', setnr), im_filt_HDR);
end
    
                    elf_support_logmsg('      Summary: All HDR scenes for environment %s have been filtered and saved to mat.\n\n', para.paths.dataset);




