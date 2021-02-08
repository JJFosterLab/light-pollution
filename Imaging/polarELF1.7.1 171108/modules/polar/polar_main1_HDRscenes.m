function polar_main1_HDRscenes(dataset, horlimit, rotation)
% POLAR_MAIN1_HDRSCENES unwarps and calibrates images.
%
% Usage: polar_main1_HDRscenes(dataset, horlimit, rotation)
%
% Inputs:
%   dataset         - 1 x 1 str, name of the dataset, must contain 5 scenes
%   horlimit        - 1 x 1 double, excentricity limit (usually 65 for polar datasets)
%   rotation        - 1 x 1 double, how much to rotate images during re-projection (usually 90 for polar datasets)
%
% Outputs:
%   None
%
% Uses: elf_paths, elf_support_logmsg, elf_para, elf_info_collect, 
%       elf_info_summarise, elf_hdr_brackets, elf_project_image, 
%       elf_readwrite, elf_hdr_calcHDR, elf_io_correctdng, elf_imread
%
% Loads files: image files in data folder
% Saves files: HDR image files in scene subfolder, *.mat files in scenes subfolder, *.tif files in proj subfolder

elf_paths;

%% check inputs
if nargin < 2 || isempty(horlimit), horlimit = NaN; end
if nargin < 1 || isempty(dataset), error('You have to provide a valid dataset name'); end 

                    elf_support_logmsg('\b\b\b\b\b\b\b\b\b\b\b\b\b\n');
                    elf_support_logmsg('----- polarELF Step 1: Unwarping images and creating HDR scenes -----\n')

%% Set up paths and file names; read info, infosum and para, calculate sets
imgformat       = '*.dng';
para            = elf_para('', dataset, imgformat);
para.paths.projfolder = 'proj';
info            = elf_info_collect(para.paths.datapath, imgformat);   % this contains EXIF information and filenames, verbose==1 means there will be output during system check
infosum         = elf_info_summarise(info);                           % summarise EXIF information for this dataset. This will be saved for later use below
infosum.linims  = strcmp(imgformat, '*.dng');                         % if linear images are used, correct for that during plotting
calib           = strcmp(imgformat, '*.dng');                         % jpgs and tifs cant be calibrated
sets            = elf_hdr_brackets(info);                             % determine which images are part of the same scene

                    elf_support_logmsg('      Processing %d scenes in environment %s\n', size(sets, 1), dataset);
                    
%% Set up projection constants
% Calculate a projection vector to transform an orthographic/equidistant/equisolid input image into an equirectangular output image
% Also creates I_info.ori_grid_x1, I_info.ori_grid_y1 (and 2) which can be used to plot a 10 degree resolution grid onto the original image

                    elf_support_logmsg('      Calculating projection constants...');

[projection_ind, infosum] = elf_project_image(infosum, para.azi, para.ele2, para.projtype, rotation); % default: 'equisolid'; also possible: 'orthographic' / 'equidistant' / 'noproj'
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
        im_raw  = elf_imread(fname);   % load the image
        im_raw  = elf_io_removehorizon(im_raw, horlimit);   % remove horizon
    
        % Calibrate and calculate intensity confidence
        [im_cal, conf, conffactors(:, i)] = elf_calibrate_abssens(im_raw, info(imnr)); 
        
        % Unwarp image
        im_proj(:, :, :, i)     = elf_project_apply(im_cal, projection_ind, [length(para.ele) length(para.azi) infosum.SamplesPerPixel]);     
        im_proj_cal(:, :, :, i) = elf_calibrate_spectral(im_proj(:, :, :, i), info(imnr), para.ana.colourcalibtype); % only needed for 'histcomb'-type intensity calculation, but not time-intensive

        conf_proj(:, :, :, i)   = elf_project_apply(conf, projection_ind, [length(para.ele) length(para.azi) infosum.SamplesPerPixel]);
    end        

    % Sort images by EV
    EV          = arrayfun(@(x) x.DigitalCamera.ExposureBiasValue, info(setstart:setend));
    [~, imInd]  = sort(EV);         % sorted EV (ascending), for HDR calculation
    im_proj     = im_proj(:, :, :, imInd);
    im_proj_cal = im_proj_cal(:, :, :, imInd);
    conf_proj   = conf_proj(:, :, :, imInd);
    conffactors = conffactors(:, imInd);

    % scale images to match middle exposure (creates a warning if scaling by more than 30%)
    [im_proj, res.scalefac] = elf_hdr_scaleStack(im_proj, conf_proj, conffactors(2:end, :));

    % Pass a figure number and an outputfilename here only if you want diagnostic pdfs.
    % However, MATLAB can't currently deal with saving these large figures, so no pdf will be created either way.
    im_HDR      = elf_hdr_calcHDR(im_proj, conf_proj, para.ana.hdrmethod, conffactors(1, :), conffactors(2:end, :)); % para.ana.hdrmethod can be 'overwrite', 'overwrite2', 'validranges', 'allvalid', 'allvalid2' (default), 'noise', para.ana.hdrmethod
    im_HDR_cal  = elf_calibrate_spectral(im_HDR, info(setstart), para.ana.colourcalibtype); % apply spectral calibration
    I           = elf_io_correctdng(im_HDR_cal, info(setstart), 'bright');

    % Save HDR file as MAT and TIF.
    elf_readwrite(para, 'saveHDR_mat', sprintf('scene%03d', setnr), im_HDR_cal);
    elf_readwrite(para, 'saveHDR_tif', sprintf('scene%03d', setnr), I);

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




