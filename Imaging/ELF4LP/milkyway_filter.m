function milkyway_filter(dataset, verbose)
% MILKYWAY_FILTER filters all HDR scenes for an environment for a 1
% degree and 10 degree resolution. The filtering algorithm (almost) fully corrects
% for the distortion introduced by the equirectangular projection.
%
% Uses: veps_log_printmsg, veps_paths, veps_para, veps_para_update, 
%       veps_readwrite, veps_support_formatA4, veps_filter 
%
% Loads files: HDR scenes as mat in scene folder
% Saves files: filtered images as mat in filt folder
% 
% Typical timing PER SCENE (on VEPSPC):
%     192s for filtering both resolutions
%     2.5s to save jpg
%     1.5s rest
%
%   + 192s PER ENVIRONMENT to calculate initial correction images

%% check inputs
if nargin < 2, verbose = true; end
if nargin < 1 || isempty(dataset), error('You have to provide a valid dataset name'); end 

                    veps_log_printmsg('\b\b\b\b\b\b\b\b\b\b\b\b\b\n');
                    veps_log_printmsg('----- ELF Step 4: Filtering scenes -----\n');

%% Set up paths and file names; read info, infosum and para
veps_paths;
para        = veps_paraJJF('', dataset, '*.dng');
para        = veps_para_update(para);                                       % Combine old parameter file with potentially changed information in current veps_para
allfiles    = dir(fullfile(para.paths.datapath, para.paths.scenefolder, '*.mat'));
fnames_im   = {allfiles.name};                                              % collect image names
infosum     = veps_readwrite(para, 'loadinfosum');                          % loads the old infosum file (which contains projection information)

                    veps_log_printmsg('      Processing %d scenes in environment %s\n', length(fnames_im), dataset);
                 
%% Set Filtering scales (specific for Milky Way)
para.ana.scales_deg = [2 4 8 16];
                    
%% Process one scene at a time
it = 1;
for setnr = 1:length(fnames_im)
    %% Load HDR image
    im_HDR  = veps_readwrite(para, 'loadHDR_mat', fnames_im{setnr});
%     % make sure image is landscape
%     if size(im_HDR, 1)>size(im_HDR, 2)
%         im_HDR = permute(im_HDR, [2 1 3]);    
%     end
    
    %% Filter HDR image
    % panel 1: original image
    fh      = veps_support_formatA4(4);
    p1      = uipanel('Parent', fh, 'Position', [0 1/3 1 2/3]);
              veps_plot_image(im_HDR, infosum, p1, 'equirectangular', 1);
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
        [im_filt_HDR, para] = veps_filter(para, im_HDR, 'scene', verbose, ah, infosum, 1);
        % veps_readwrite(para, 'saveinfosum', [], infosum); % saves infosum AND para
    else
        im_filt_HDR = veps_filter(para, im_HDR, 'scene', verbose, ah, infosum);
    end
    it      = it+1;
    
    %% save projected and filtered images
    veps_readwrite(para, 'savefilt_mat', sprintf('scene%03d', setnr), im_filt_HDR);
end
    
                    veps_log_printmsg('      Summary: All HDR scenes for environment %s have been filtered and saved to mat.\n\n', para.paths.dataset);




