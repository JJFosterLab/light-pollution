function night_HDRscenes_quickfixJJF(dataset, horlimit)
% NIGHT_HDRSCENES unwarps all images in a data set, sorts them into
% scenes, and calculates HDR representations of these scenes. A
% representative HDR image of each scene is saved in the "scene" subfolder,
% but these are not the images used for later caulations.
%
% Uses: elf_paths, elf_support_logmsg, elf_para, elf_info_collect, 
%       elf_info_summarise, elf_hdr_brackets, elf_project_image, 
%       elf_readwrite, elf_hdr_calcHDR, elf_io_correctdng, elf_imread
%
% Loads files: image files in data folder
% Saves files: HDR image files in scene subfolder, *.mat files in scenes subfolder, *.tif files in proj subfolder
% 
% Typical timing PER SCENE (on ELFPC):
%       12s total (PROFILE AGAIN! SOME WEIRD RESULTS)
%
% Detailed calculation:
%       Unwarping: ~9s/2s/1.5s/1s with delinearised plotting/linear plotting/no plotting/no saving
%       HDR calculation incl. unwarping: 25s/15s with/without plotting of HDR diagnostic plots
%                                  12s without saving any of the tifs
%       4s for everything else

elf_paths;

%% check inputs
if nargin < 2 || isempty(horlimit), horlimit = NaN; end
if nargin < 1 || isempty(dataset), error('You have to provide a valid dataset name'); end 

                    elf_support_logmsg('\b\b\b\b\b\b\b\b\b\b\b\b\b\n');
                    elf_support_logmsg('----- nightELF Step 1: Unwarping images and creating HDR scenes -----\n')

%% Set up paths and file names; read info, infosum and para, calculate sets
imgformat       = '*.dng';
para            = elf_para('', dataset, imgformat);
info            = elf_info_collect(para.paths.datapath, imgformat);    % this contains EXIF information and filenames, verbose==1 means there will be output during system check
infosum         = elf_info_summarise(info);                            % summarise EXIF information for this dataset. This will be saved for later use below
infosum.linims  = strcmp(imgformat, '*.dng');                          % if linear images are used, correct for that during plotting
calib           = strcmp(imgformat, '*.dng');                          % jpgs and tifs cant be calibrated
sets            = elf_hdr_brackets(info);                              % determine which images are part of the same scene

                    elf_support_logmsg('      Processing %d scenes in environment %s\n', size(sets, 1), dataset);
                    
%% Set up projection constants
% Calculate a projection vector to transform an orthographic/equidistant/equisolid input image into an equirectangular output image
% Also creates I_info.ori_grid_x1, I_info.ori_grid_y1 (and 2) which can be used to plot a 10 degree resolution grid onto the original image

                    elf_support_logmsg('      Calculating projection constants...');

[projection_ind, infosum] = elf_project_image(infosum, para.azi, para.ele, para.projtype); % default: 'equisolid'; also possible: 'orthographic' / 'equidistant' / 'noproj'
elf_readwrite(para, 'saveinfosum', [], infosum); % saves infosum AND para for use in later stages
                    
                    elf_support_logmsg('\b\b\b\b\b\b\b\b\b\b\b\b\bdone.\n');

%% Step 1: Unwarp images and calculate HDR scenes

elf_support_logmsg('      Unwarping images scene-by-scene (%d scenes), and creating HDR images. Projected time: ', size(sets, 1));

tic; % Start taking time

% Process one scene at a time
for setnr = 1:size(sets, 1)
    clear im_filt res
    
    setstart = sets(setnr, 1);          % first image in this set
    setend   = sets(setnr, 2);          % last image in this set
    numims   = setend - setstart + 1;   % total number of images in this set
    im_cal   = zeros(length(para.ele), length(para.azi), info(1).SamplesPerPixel, numims); % pre-allocate
    conf     = zeros(3, 2, numims);     % pre-allocate
    
    for i = 1:numims % for each image in this set
        % Load image        
        imnr    = setstart + i - 1;     % the number of this image
        fname   = info(imnr).Filename;  % full path to input image file
        im_ori  = elf_imread(fname);   % load the image

        % Remove horizon for nohor sets
        if ~isnan(horlimit)
            mid    = [1+(size(im_ori, 1)-1)/2; 1+(size(im_ori, 2)-1)/2];        % centre of image
            r_full = 8 * size(im_ori, 1) / 24;                                 % theoretical value for 24mm high chip that is fully covered by fisheye circular image
            [y, x] = meshgrid(1:size(im_ori, 2), 1:size(im_ori, 1));
            r      = sqrt((x-mid(1)).^2 + (y-mid(2)).^2);
            ele    = asind(r / 2 / r_full) * 2;%r/r_full * 90;

            sel    = ele>horlimit;
            i1 = im_ori(:, :, 1); i1(sel) = NaN;
            i2 = im_ori(:, :, 2); i2(sel) = NaN;
            i3 = im_ori(:, :, 3); i3(sel) = NaN;
            im_ori = cat(3, i1, i2, i3);
        end
    
        % Umwarp image
        im_proj = elf_project_apply(im_ori, projection_ind, [length(para.ele) length(para.azi) size(im_ori, 3)]);
        para.paths.projfolder = 'proj';%QUICK FIX
        % Save the projected image
        elf_readwrite(para, 'saveproj_tif', fname, im_proj);
        
        % Calibrate and calculate intensity confidence
        % If calib is false, the image is merely normalised to the maximum depending on image bit-depth and conf is NaN
%         [im_cal(:, :, :, i), conf(:, :, i)] = elf_calibrate_abssens(im_proj, info(imnr));%, calib); QUICKFIX TOO MANY INPUTS
%THIS NO LONGER WORKS WITH 2019 VERSION OF elf_calibrate_abssens
%"left side is 3-by-2 and the size of the right side is 1801-by-1801-by-3"
%I'VE REPLACED IT WITH THE VERSION FROM elf_main1_HdrAndInt
         [im_cal, conf] = elf_calibrate_abssens(im_proj, info(imnr));%THIS DID NOT WORK
    end
    
    % Calculate HDR image
    % Pass a figure number and an outputfilename here only if you want diagnostic pdfs.
    % However, MATLAB can't currently deal with saving these large figures, so no pdf will be created either way.
    im_HDR = elf_hdr_calcHDR(im_cal, conf, info(setstart:setend));%, 50+setnr);%, fullfile(para.paths.datapath, sprintf('set%3d_diag.pdf', setnr))); 
    
    % save HDR file as tif and mat
    I = elf_io_correctdng(im_HDR, info(setstart), 'bright');
    elf_readwrite(para, 'saveHDR_tif', sprintf('scene%03d', setnr), I)
    elf_readwrite(para, 'saveHDR_mat', sprintf('scene%03d', setnr), im_HDR)
    
                    if setnr == 1
                        elf_support_logmsg('\b\b\b\b\b\b\b\b\b\b\b\b\b%.2f minutes.\n', toc/60*size(sets, 1));
                        elf_support_logmsg('      Scene: 1..');
                    elseif mod(setnr-1, 20)==0
                        elf_support_logmsg('\b\b\b\b\b\b\b\b\b\b\b\b\b\n');
                        elf_support_logmsg('             %d..', setnr);
                    else
                        elf_support_logmsg('\b\b\b\b\b\b\b\b\b\b\b\b\b%d..', setnr);
                    end
end

                    elf_support_logmsg('\b\b\b\b\b\b\b\b\b\b\b\b\bdone.\n');    
                    elf_support_logmsg('      Summary: All HDR scenes for environment %s calculated and saved to mat and tif.\n\n', dataset); % write confirmation to log




