function polar_main2_filter(dataset, verbose)
% POLAR_MAIN2_FILTER filters images to 2/4/8/16 degrees. The filtering algorithm (almost) fully corrects
% for the distortion introduced by the equirectangular projection.
%
% Usage: polar_main2_filter(dataset, horlimit, verbose)
%
% Inputs:
%   dataset         - 1 x 1 str, name of the dataset, must contain 5 scenes
%   verbose         - 1 x 1 bool, triggers verbose output suring filtering (default: true)
%
% Outputs:
%   None
%
% Uses: elf_support_logmsg, elf_paths, elf_para, elf_para_update, 
%       elf_readwrite, elf_support_formatA4, elf_filter 
%
% Loads files: HDR scenes as mat in scene folder
% Saves files: filtered images as mat in filt folder

%% check inputs
if nargin < 2, verbose = true; end
if nargin < 1 || isempty(dataset), error('You have to provide a valid dataset name'); end 

                    elf_support_logmsg('\b\b\b\b\b\b\b\b\b\b\b\b\b\n');
                    elf_support_logmsg('----- polarELF Step 2: Filtering scenes -----\n');

%% Set up paths and file names; read info, infosum and para
elf_paths;
para        = elf_para('', dataset, '*.dng');
para        = elf_para_update(para);                                       % Combine old parameter file with potentially changed information in current elf_para
allfiles    = elf_io_dir(fullfile(para.paths.datapath, para.paths.scenefolder, '*.mat'));
fnames_im   = {allfiles.name};                                              % collect image names
infosum     = elf_readwrite(para, 'loadinfosum');                          % loads the old infosum file (which contains projection information)

                    elf_support_logmsg('      Processing %d scenes in environment %s\n', length(fnames_im), dataset);
                 
%% Set filtering scales
para.ana.scales_deg = [2 4 8 16];
                    
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
    ah(2)   = axes('Parent', p3, 'Position', [0 0 1 1]); 
            axis(ah(1:2), 'off');
            set(fh, 'Name', sprintf('Scene #%d of %d', setnr, length(allfiles)));
            drawnow;
    ah(3:4) = ah(1:2);
    
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




