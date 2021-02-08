function milkyway_HDRscenes(dataset, horlimit)
% MILKYWAY_HDRSCENES unwarps all images in a data set, sorts them into
% scenes, and calculates HDR representations of these scenes. A
% representative HDR image of each scene is saved in the "scene" subfolder,
% but these are not the images used for later caulations.
%
% Uses: veps_paths, veps_log_printmsg, veps_para, veps_info_collect, 
%       veps_info_summarise, veps_bracketing_detect, veps_project_image, 
%       veps_readwrite, veps_calcHDR, veps_io_correctdng, veps_imread
%
% Loads files: image files in data folder
% Saves files: HDR image files in scene subfolder, *.mat files in scenes subfolder, *.tif files in proj subfolder
% 
% Typical timing PER SCENE (on VEPSPC):
%       12s total (PROFILE AGAIN! SOME WEIRD RESULTS)
%
% Detailed calculation:
%       Unwarping: ~9s/2s/1.5s/1s with delinearised plotting/linear plotting/no plotting/no saving
%       HDR calculation incl. unwarping: 25s/15s with/without plotting of HDR diagnostic plots
%                                  12s without saving any of the tifs
%       4s for everything else

veps_paths;

%% check inputs
if nargin < 2 || isempty(horlimit), horlimit = NaN; end
if nargin < 1 || isempty(dataset), error('You have to provide a valid dataset name'); end 

                    veps_log_printmsg('\b\b\b\b\b\b\b\b\b\b\b\b\b\n');
                    veps_log_printmsg('----- ELF Step 1: Calibration, HDR and Intensity -----\n')

%% Set up paths and file names; read info, infosum and para, calculate sets
imgformat       = '*.dng';%
para            = veps_paraJJF('', dataset, imgformat);
info            = veps_info_collect(para.paths.datapath, imgformat);   % this contains EXIF information and filenames, verbose==1 means there will be output during system check
infosum         = veps_info_summarise(info);                  % summarise EXIF information for this dataset. This will be saved for later use below
infosum.linims  = strcmp(imgformat, '*.dng');%'*.NEF');%                          % if linear images are used, correct for that during plotting
sets            = veps_bracketing_detect(info);                        % determine which images are part of the same scene

                    veps_log_printmsg('      Processing %d scenes in environment %s\n', size(sets, 1), dataset);
                    
%% Set up projection constants
% Calculate a projection vector to transform an orthographic/equidistant/equisolid input image into an equirectangular output image
% Also creates I_info.ori_grid_x1, I_info.ori_grid_y1 (and 2) which can be used to plot a 10 degree resolution grid onto the original image

                    veps_log_printmsg('      Calculating projection constants...');

[projection_ind, infosum] = veps_project_image(infosum, para.azi, para.ele, para.projtype); % default: 'equisolid'; also possible: 'orthographic' / 'equidistant' / 'noproj'
veps_readwrite(para, 'saveinfosum', [], infosum); % saves infosum AND para for use in later stages
                    
                    veps_log_printmsg('\b\b\b\b\b\b\b\b\b\b\b\b\bdone.\n');

%% Step 1: Unwarp images and calculate HDR scenes
tic; % Start taking time
% Process one scene at a time
for setnr = 1:size(sets, 1)
    clear im_filt res
    
    setstart    = sets(setnr, 1);          % first image in this set
    setend      = sets(setnr, 2);          % last image in this set
    numims      = setend - setstart + 1;   % total number of images in this set
    im_proj     = zeros(length(para.ele), length(para.azi), infosum.SamplesPerPixel, numims);  % pre-allocate
    conf_proj   = zeros(length(para.ele), length(para.azi), infosum.SamplesPerPixel, numims);  % pre-allocate
    conffactors = zeros(4, numims);        % pre-allocate
    
    for i = 1:numims % for each image in this set
        % Load image        
        imnr                    = setstart + i - 1;     % the number of this image
        fname                   = info(imnr).Filename;  % full path to input image file
        im_raw                  = double(veps_imread(fname));   % load the image (double)

        % Remove horizon for nohor sets
        if ~isnan(horlimit)
            mid    = [1+(size(im_raw, 1)-1)/2; 1+(size(im_raw, 2)-1)/2];        % centre of image
            r_full = 8 * size(im_raw, 1) / 24;                                 % theoretical value for 24mm high chip that is fully covered by fisheye circular image
            [y, x] = meshgrid(1:size(im_raw, 2), 1:size(im_raw, 1));
            r      = sqrt((x-mid(1)).^2 + (y-mid(2)).^2);
            ele    = asind(r / 2 / r_full) * 2;%r/r_full * 90;

            sel    = ele>horlimit;
            i1     = im_raw(:, :, 1); i1(sel) = NaN;
            i2     = im_raw(:, :, 2); i2(sel) = NaN;
            i3     = im_raw(:, :, 3); i3(sel) = NaN;
            im_raw = cat(3, i1, i2, i3);
        end
        
        % Calibrate and calculate intensity confidence
        [im_cal, conf, conffactors(:, i)] = veps_calibrate_abssensJJF(im_raw, info(imnr)); 
        
        % Umwarp image
        im_temp                 = im_cal(projection_ind); % Apply projection vector to transform image into desired projection...
        im_proj(:, :, :, i)     = reshape(im_temp, [length(para.ele) length(para.azi) infosum.SamplesPerPixel]); % ...and reshape the resulting vector into an image matrix        

        conf_temp               = conf(projection_ind); % Apply projection vector to transform image into desired projection...
        conf_proj(:, :, :, i)   = reshape(conf_temp, [length(para.ele) length(para.azi) infosum.SamplesPerPixel]); % ...and reshape the resulting vector into an image matrix
    end
    
    % Replace NaNs with 0s (otherwise filtering will be crazy)
    im_proj(isnan(im_proj)) = 0;
    
    % Sort images by EV
    EV          = arrayfun(@(x) x.DigitalCamera.ExposureBiasValue, info(setstart:setend));
    [~, imInd]  = sort(EV);         % sorted EV (ascending), for HDR calculation
    im_proj     = im_proj(:, :, :, imInd);
    conf_proj   = conf_proj(:, :, :, imInd);
    conffactors = conffactors(:, imInd);
    
    % scale images to match middle exposure (creates a warning if scaling by more than 30%)
    [im_proj, res.scalefac] = veps_scaleStack(im_proj, conf_proj, conffactors(2:end, :));
    
    % Pass a figure number and an outputfilename here only if you want diagnostic pdfs.
    % However, MATLAB can't currently deal with saving these large figures, so no pdf will be created either way.
    im_HDR      = veps_calcHDR(im_proj, conf_proj, para.ana.hdrmethod, conffactors(1, :), conffactors(2:end, :)); % para.ana.hdrmethod can be 'overwrite', 'overwrite2', 'validranges', 'allvalid', 'allvalid2' (default), 'noise', para.ana.hdrmethod
    im_HDR_cal  = veps_calibrate_spectralJJF(im_HDR, info(setstart), para.ana.colourcalibtype); % apply spectral calibration
    I           = veps_io_correctdng(im_HDR_cal, info(setstart), 'bright');

    % Save HDR file as MAT and TIF.
    veps_readwrite(para, 'saveHDR_mat', sprintf('scene%03d', setnr), im_HDR_cal);
    veps_readwrite(para, 'saveHDR_tif', sprintf('scene%03d', setnr), I);
    
                    if setnr == 1
                        veps_log_printmsg('      Starting scene-by-scene calibration, HDR creation and intensity analysis. Projected time: %.2f minutes.\n', toc/60*size(sets, 1));
                        veps_log_printmsg('      Scene: 1..');
                    elseif mod(setnr-1, 20)==0
                        veps_log_printmsg('\b\b\b\b\b\b\b\b\b\b\b\b\b\n');
                        veps_log_printmsg('             %d..', setnr);
                    else
                        veps_log_printmsg('\b\b\b\b\b\b\b\b\b\b\b\b\b%d..', setnr);
                    end
end
                    veps_log_printmsg('\b\b\b\b\b\b\b\b\b\b\b\b\bdone.\n');    
                    veps_log_printmsg('      Summary: All HDR scenes for environment %s calculated and saved to mat and tif.\n\n', dataset); % write confirmation to log




