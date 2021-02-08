function [darkmean, darkstd, saturation] = elf_calibrate_darkandreadout(camstring, exp, iso, info, long_exposure_NR)
% ELF_CALIBRATE_DARKANDREADOUT returns an estimate of dark noise and saturation levels for a given camera setting
% exp and iso can be scalars or vectors

if nargin<5, long_exposure_NR = true; end

switch lower(camstring)
        case {'nikon d810', ''}
%             %% OLD SIMPLE MODEL
%             darkmean    = 600.57 + 0.00071*iso;
%             darkstd     = 0.93 + 0.0063*iso;
%             saturation  = 15992 - darkmean; % (15520 - black) in EXIF file
        
            para = elf_para;
            calibfilefolder = para.paths.calibfolder; % Where to find the finished calibration files
            load(fullfile(calibfilefolder, lower(camstring), 'iso_noise.mat'), 'darkiso_mean', 'darkiso_std', 'darkiso_x'); 
            load(fullfile(calibfilefolder, lower(camstring), 'exp_noise.mat'), 'darkexp_mean', 'darkexp_std', 'darkexp_x');
            
            darkmean = zeros(length(iso), 3);
            darkstd  = zeros(length(iso), 3);
            
            for i = 1:length(exp) % for each iso/exp pair
                % MEAN
                ind_iso      = sub_findclosest(iso(i), darkiso_x);
                darkmean_iso = darkiso_mean(ind_iso, :);
                ind_exp      = sub_findclosest(exp(i), darkexp_x);
                darkmean_exp = darkexp_mean(ind_exp, :);
                darkmean_100 = darkiso_mean(darkiso_x==100, :);
                darkmean(i, :)  = darkmean_iso + darkmean_exp - darkmean_100;  % additive model
                       
                % STD
                darkstd_iso  = darkiso_std(ind_iso, :);
                darkstd_exp  = darkexp_std(ind_exp, :);
                darkstd_100  = darkiso_std(darkiso_x==100, :);
                darkstd(i, :)   = darkstd_iso .* darkstd_exp ./ darkstd_100;   % multipicative model
            end
            
            if long_exposure_NR
                % if Long Exposure Noise Reduction was turned on in the camera, noise for exposures above 1s is removed in camera
                darkmean(exp>=1, :) = 600;            
            end
            
            saturation   = 15520; % 15992 was found in calibration, 15520 is the black value in EXIF file

    otherwise
            darkmean     = info.SubIFDs{1}.BlackLevel(1) * ones(length(exp), 3);   % black level saved by camera in exif file. this seems to very closely correspond to measured readout noise
            saturation   = info.SubIFDs{1}.WhiteLevel;      % white level, CHECK whether this corresponds to a reasonable saturation level
            darkstd      = nan(length(exp), 3);
end

end % main

%% sub functions
function ind = sub_findclosest(a, B)
    % finds the element of B that is closest to scalar a
    [~, ind] = min(abs(B-a));
end