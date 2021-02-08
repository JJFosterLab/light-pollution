function [C, Ccs, conthist, contfft, contprof] = space_analysis_channel(im, aerange, aeacc, numangles, dist_r_deg, channeltype, meantype, meanthr, bins, verbose)
% space_analysis_channel calculates one particular luminance or colour contrast in a prefiltered image.
%   For ideal sampling, the spacing of sampling points in azi/ele should be half the hwlfwidth (FWHM) of the Gaussian used for filtering. The
%   horizontal contrasts are corrected for elevation distortion, and the means are weighted to take this into account.
%
%   Example: 
%   [C, Cr_v, binmean] = space_analysis_channel(im, azi, ele, azifull, type[, bins, verbose])
%
% Inputs:
% im                - MxNxC double, projected into a equirectangular (azimuth / elevation) projection
% azi/ele           - 1xN / 1xM doubles, in degrees, y/x-vectors into im{1}, 
%                       ele will usually be a reversed vector, starting at the highest value 
%                       because Matlab images have their origin at the top left corner
% azifull           - 1xP double, y-vector into im{2}
% channeltype       - string, can be 'lum', 'rg', 'yb', 'gb' (If channels other than RGB are used, this has to be modified!)
% [meantype]        - string, can be 'mean', 'rms', 'perc' (default 'rms')
% [meanthr]         - double, the mean is only calculated for contrast values greater or equal to this value, e.g. 5 for contrasts >=5% (default 0).
% [bins]            - 1xB double, bin boundaries in degrees (default: [-90 -50 -10 10 50 90])
% [verbose]         - logical, triggering plotting of all measured contrasts. For 1-degree scale, this takes a LONG time! (default: false)
% [info]            - info file for images; only required if verbose == true
% 
% Outputs:
% C                 - (B-1)xnumangles double, mean contrasts for each bin and angle
% binmean           - (B-1)x1 double, mean elevation for each bin, can be used for plotting
% contrhist_h       - 101x(B-1) double, histogram of horizontal contrasts in bins of 1 degree
% contrhist_v       - 101x(B-1) double, histogram of vertical contrasts in bins of 1 degree
%
% Uses:       None
% Used by:    space_analysis
% Call stack: elf_main2_descriptors -> elf_analysis -> space_analysis -> space_analysis_channel
% See also:   elf_main2_descriptors, elf_analysis, space_analysis

%% check inputs
        
%% spatial analysis
if verbose
    fprintf('Calculating contrasts in the %s channel using %s-type means (threshold %g%%).\n', channeltype, meantype, meanthr);
end

%% re-calculate some basic parameters
ele          = aerange(4):-aeacc(2):aerange(3); % reverse elevation vector
angles       = linspace(0, 360, numangles+1);     % 1xnumangles; the angles to be tested (in degrees)
angles       = angles(1:end-1);                   % remove 360 deg
dist_el_deg  = sind(angles) * dist_r_deg;    % 1xnumangles; distance along elevation, in degrees
ele_cont     = repmat(ele', 1, numangles) + repmat(dist_el_deg, length(ele), 1)/2;      % the mean elevation for each contrast
elecorr_cont = 1 ./ cosd(ele_cont); % correction factor for horizontal stretching at high elevations (for directional contrasts)
elecorr      = 1 ./ cosd(ele); % correction factor for horizontal stretching at high elevations (for centre surround contrasts)
        
%% get contrast index matrices
if verbose
    fprintf('Searching for contrast index matrix...'); % Line will be finished in subfunction
end
cimpath = fullfile(fileparts(mfilename('fullpath')), 'CIM');
[CIMx, CIMy] = sub_ContrIndMat(aerange, aeacc, numangles, dist_r_deg, cimpath, 0, verbose);

%% select channels and calculate contrasts
if verbose
    fprintf('Calculating local contrasts...'); 
end
switch channeltype
    case 'lum' % calculate Michelson contrast: |L1-L2| / (L1+L2)
        L    = mean(im, 3);  % Luminance is simply the mean of all other channels
        [cont, contCS] = sub_contrast(L, CIMx, CIMy, 'lum', verbose); % calculate the full NxMxnumangles contrast matrix %sub_contrast(im, CIMx, CIMy, ctype, plotit)
        
    case {'rg', 'yb', 'gb', 'rb'}
        % separate channels
        R    = im(:, :, 1);
        G    = im(:, :, 2);
        B    = im(:, :, 3);
        Y    = (R + G) / 2;

        switch channeltype
            case 'rg' % calculate (R-G)/(R+G) for each pixel, then diff between pixels
                H = (R - G) ./ (R + G);      % between -1 and 1
            case 'yb' % calculate (Y-B)/(Y+B) for each pixel, then diff between pixels
                H = (Y - B) ./ (Y + B);
            case 'gb' % calculate (G-B)/(G+B) for each pixel, then diff between pixels
                H = (G - B) ./ (G + B);
            case 'rb' % calculate (R-B)/(R+B) for each pixel, then diff between pixels
                H = (R - B) ./ (R + B);
        end
        [cont, contCS] = sub_contrast(H, CIMx, CIMy, 'col', verbose); % calculate the full NxMxnumangles contrast matrix
   
    otherwise
        error('Unknown spatial comparison type.');
end

%% average over bins
if verbose, fprintf('and averaging them.\n'); end
c   = meanthr/100; % contrast threshold
C   = nan(length(bins)-1, size(cont, 3)); % prealloc
Ccs = nan(length(bins)-1, 1); % prealloc

for bin = 1:length(bins)-1 % for each elevation bin
    for ang = 1:size(cont, 3) % for each direction angle
        thisele = ele_cont(:, ang);
        % select the rows that have a mean elevation inside this bin
        rows = thisele >= bins(bin) & thisele <= bins(bin+1);

        % extract rows
        thiscont    = abs(cont(rows, :, ang)); % Take the ABSOLUTE contrast before calculating thresholds
        contweights = 1 ./ repmat(elecorr_cont(rows, ang), [1, size(thiscont, 2)]);
        sel         = thiscont>=c; % select only entries above threshold

        % calculate averages
        switch meantype
            case 'mean'
                C(bin, ang) = sum(thiscont(sel) .* contweights(sel)) ./ sum(contweights(sel));
            case 'perc'
                C(bin, ang) = sum(sel .* contweights(sel)) ./ sum(contweights); % UNSURE whether this is correct
            case 'rms'
                C(bin, ang) = sqrt( sum(thiscont(sel).^2 .* contweights(sel)) ./ sum(contweights(sel)) );
            otherwise
                error('Unknown mean type: %s', meantype);
        end

        if C(bin, ang)==0 || isnan(C(bin, ang))
            warning('zero/nan contrast in bin %d', bin);
            C(bin, ang) = 1/length(thiscont(:)); 
        end % remove zeros

        % also save contrast histograms
        conthist(:, bin, ang) = histc(thiscont(:), 0:.01:1);
        
%         figure(55); I = cont(:, :, 7); imagesc(1-abs(I)/max(abs(I(:)))); colormap(gray); title(['90 deg: ' num2str(C(bin, ang))]);
%         figure(56); I = cont(:, :, 13); imagesc(1-abs(I)/max(abs(I(:)))); colormap(gray); title(['-90 deg: ' num2str(C(bin, ang))]);
%         figure(10+ang); imagesc(thisele(rows), 1:size(thiscont, 2), thiscont/max(thiscont(:))); axis equal; colormap(gray); title(num2str(C(bin, ang))); xlabel('pixels'); ylabel('elevation');
%         figure(10+ang); e = thisele(rows); imagesc(thiscont/max(thiscont(:))); set(gca, 'YTick', 1:10:length(e), 'YTickLabel', num2str(e(1:10:end))); axis image; colormap(gray); title(num2str(C(bin, ang))); xlabel('pixels'); ylabel('elevation');

    end
    
    %% Average CS contrasts
        % select the rows that have a mean elevation inside this bin
        rows = ele >= bins(bin) & ele <= bins(bin+1);
        % extract rows
        thiscontCS  = abs(contCS(rows, :)); % Take the ABSOLUTE contrast before calculating thresholds
        contweights = 1 ./ repmat(elecorr(rows)', [1, size(thiscontCS, 2)]);
        sel         = thiscontCS>=c; % select only entries above threshold

        % calculate averages
        switch meantype
            case 'mean'
                Ccs(bin) = sum(thiscontCS(sel) .* contweights(sel)) ./ sum(contweights(sel));
            case 'perc'
                Ccs(bin) = sum(sel .* contweights(sel)) ./ sum(contweights); % UNSURE whether this is correct
            case 'rms'
                Ccs(bin) = sqrt( sum(thiscontCS(sel).^2 .* contweights(sel)) ./ sum(contweights(sel)) );
            otherwise
                error('Unknown mean type: %s', meantype);
        end

        if C(bin)==0 || isnan(C(bin))
            warning('zero/nan contrast in bin %d', bin);
            Ccs(bin) = 1/length(thiscont(:)); 
        end % remove zeros

  
    
%     %% Calculate ffts
%     switch channeltype
%         case 'lum' % calculate Michelson contrast: |L1-L2| / (L1+L2)
%             thiszone = L(rows, :);
%         case {'rg', 'yb', 'gb', 'rb'}
%             thiszone = H(rows, :);
%     end
%     contfft{bin} = abs(log2(fftshift(fft2(thiszone))));
    contfft{bin} = NaN;
    
end

%% Average over strips of elevation (for profiles)
imh             = size(im, 1);                      % image height
hdivn           = 60; %%FIXME get this from para.ana.hdivn_int
stripwidth_m1   = (imh-1) / hdivn;                                  % every strip has this many plus one rows
stripborders    = round(1:stripwidth_m1:imh);
hor_cut         = [stripborders(1:end-1); stripborders(2:end)]';    % border elements are included in BOTH strips
contprof        = zeros(size(hor_cut, 1), 1); %prealloc
for kk = 1:size(hor_cut, 1)        % for each strip
    rows        = hor_cut(kk, 1):hor_cut(kk, 2);  % these are the row indices for this strip
    thiscont    = abs(contCS(rows, :)); % Take the ABSOLUTE contrast before calculating thresholds; for all azimuths, for all angles
    sel         = thiscont>=c; % select only entries above threshold

    % calculate averages contrast profile
    switch meantype
        case 'mean'
            contprof(kk) = mean(thiscont(sel));
        case 'perc'
            contprof(kk) = sum(sel(:)) / length(sel(:)); % UNSURE whether this is correct
        case 'rms'
            contprof(kk) = rms(thiscont(sel));
        otherwise
            error('Unknown mean type: %s', meantype);
    end

    if contprof(kk)==0 || isnan(contprof(kk))
        warning('zero/nan contrast in strip %d', kk);
        contprof(kk) = 1/length(thiscont(:)); 
    end % remove zeros
end

contim = []; %%FIXME

return

%% create contrast images

%% and plot them
      
%% calculate more contrast plots

imsize = size(im{1}, 1);
c = 0.01; % mark anything with a contrast higher than 1%
% vertical
a = [vert_contr>=c;false(1, imsize)] | [false(1, imsize);vert_contr>=c];
outim = im{1};
for i = 1:imsize
    for j = 1:imsize
        if ~a(i, j)
            bw = mean(outim(i, j, :));
            outim(i, j, :) = bw;
        end
    end
end

contrim{2, 3} = outim;

   
end

%% subfunctions
function [cont, contCS] = sub_contrast(im, CIMx, CIMy, ctype, plotit)
% calculates the contrast between pixels in the image channel im (NxM) of
% type ctype ('lum'/'col') using the Contrast Index Matrices CIMx and CIMy
%
% cont is a NxMxnumangles matrix


    %% Calculate
    cont       = nan(size(CIMx)); % prealloc
    im_shifted = cont;            % prealloc
    for i = 1:size(CIMx, 3) % for each angle
        linind                  = sub2ind(size(im), CIMx(:, :, i), CIMy(:, :, i));
        linind(isnan(linind))   = 1; % These will later be replaced by Nans again
        imtemp                  = im(linind); % shift the image
        imtemp(isnan(CIMx(:, :, i))) = NaN; % re-insert NaNs at unavailable pixels

        im_shifted(:, :, i) = imtemp;

        switch ctype
            case 'lum'
                cont(:, :, i) = (im - imtemp) ./ (im + imtemp); % Michelson contrast
            case 'col'
                cont(:, :, i) = (im - imtemp) ./ 2; % Michelson contrast
        end
        
        %Something like this would have to be done to filter contrasts at a
        %90 degree angle (here only one step forwards, but should
        %potentially be 3 or 5 steps)
        im0 = cont(:, :, i);
        
        toshiftplus90  = i+size(CIMx, 3)/4;
        toshiftminus90 = i-size(CIMx, 3)/4;
        toshiftplus90  = mod(toshiftplus90-1, size(CIMx, 3))+1;
        toshiftminus90 = mod(toshiftminus90-1, size(CIMx, 3))+1;
        
        linind                  = sub2ind(size(im), CIMx(:, :, toshiftplus90), CIMy(:, :, toshiftplus90));
        linind(isnan(linind))   = 1; % These will later be replaced by Nans again
        implus90                = im0(linind); % shift the image
        implus90(isnan(CIMx(:, :, toshiftplus90))) = NaN; % re-insert NaNs at unavailable pixels

        cx                      = CIMx(:, :, toshiftplus90);
        cxx                     = cx(linind);
        cxx(isnan(cx))          = NaN; % re-insert NaNs at unavailable pixels
        cy                      = CIMy(:, :, toshiftplus90);
        cyy                     = cy(linind);
        cyy(isnan(cy))          = NaN; % re-insert NaNs at unavailable pixels
        linind_2                = sub2ind(size(im), cxx, cyy);
        linind_2(isnan(linind_2))   = 1; % These will later be replaced by Nans again
        implus90_2              = im0(linind_2); % shift the image
        implus90_2(isnan(cxx))  = NaN; % re-insert NaNs at unavailable pixels

        linind                  = sub2ind(size(im), CIMx(:, :, toshiftminus90), CIMy(:, :, toshiftminus90));
        linind(isnan(linind))   = 1; % These will later be replaced by Nans again
        imminus90                = im0(linind); % shift the image
        imminus90(isnan(CIMx(:, :, toshiftminus90))) = NaN; % re-insert NaNs at unavailable pixels

        cx                      = CIMx(:, :, toshiftminus90);
        cxx                     = cx(linind);
        cxx(isnan(cx))          = NaN; % re-insert NaNs at unavailable pixels
        cy                      = CIMy(:, :, toshiftminus90);
        cyy                     = cy(linind);
        cyy(isnan(cy))          = NaN; % re-insert NaNs at unavailable pixels
        linind_2                = sub2ind(size(im), cxx, cyy);
        linind_2(isnan(linind_2))   = 1; % These will later be replaced by Nans again
        imminus90_2             = im0(linind_2); % shift the image
        imminus90_2(isnan(cxx)) = NaN; % re-insert NaNs at unavailable pixels
        
        
        contfilt(:, :, i) = (im0+implus90+imminus90+implus90_2+imminus90_2)/5;
              
    end
    cont = contfilt;
    % ALSO calculate centre-surround contrast
    im_surround = nanmean(im_shifted, 3); % The surround value is the mean of all directional values
    switch ctype
        case 'lum'
            contCS = (im - im_surround) ./ (im + im_surround); % Michelson contrast
        case 'col'
            contCS = (im - im_surround) ./ 2; % Michelson contrast
    end

    %% plot
    if plotit
        figure(37); clf;
        axdist = 0.2;
        axsize = 0.14;
        angles = linspace(0, 360, size(CIMx, 3)+1);     % 1xnumangles; the angles to be plotted (in degrees); contains one too many, but that won't be used

        centrepos = [0.5-axsize/2 0.5-axsize/2 axsize axsize];
        subplot('position', centrepos); I = im; imagesc(I/max(I(:))); axis image off; title('original image channel');
        for i = 1:2:size(CIMx, 3)
            pos1 = centrepos + [cosd(angles(i)) sind(angles(i)) 0 0] * axdist + [axsize/4 axsize/4 -axsize/2 -axsize/2];
            pos2 = centrepos + [cosd(angles(i)) sind(angles(i)) 0 0] * axdist * 2;
            subplot('position', pos1); I = im_shifted(:, :, i); imagesc(I/max(I(:))); axis image off; %title('shifted image channel');
            subplot('position', pos2); I = abs(cont(:, :, i)); imagesc(1-I/max(I(:))); axis image off; %title('absolute image difference');
            %pause;
            subplot('position', pos2); I = abs(contfilt(:, :, i)); imagesc(1-I/max(I(:))); axis image off; %title('absolute image difference');
            
        end
        colormap gray
        
        figure(38); clf;
        imagesc(1-abs(contCS)/max(abs(contCS(:))));
        colormap gray
        axis equal
        drawnow
        
    end

end % subfunction


function [xtarget_all, ytarget_all, actualdist, actualangle] = sub_ContrIndMat(aerange, aeacc, numangles, dist_r_deg, CIMfolder, forcerecalc, verbose)
% calculates the image pixel indices for contrast comparisons in an image sampled in the azimuth/elevation range aerange
% with accuracy aeacc. The aim is to calculate the difference between two pixels that are dist_r_deg (1x1; in degrees) apart 
% along a number of different angles (1x1 numangles; in degrees);
% returns two NxMxnumangles contrast index matrices (along x and y)
% aerange is range; [azmin azmax elmin elmax]
% aeacc is accuracy; [azacc elacc]
% numangles is the number of angles to sample, default 12
% dist_r_deg is the desired distance between neighbouring pixels, default 1

    %% extract individual parameters
    azmin = aerange(1); % def -90
    azmax = aerange(2); % def 90
    elmin = aerange(3); % def -90
    elmax = aerange(4); % def -90
    azacc = aeacc(1); % deg / pix; def 0.1
    elacc = aeacc(2); % deg / pix; def 0.1

    %% construct filename and check whether it exists already
    filename = sprintf('ContrIndMat_%g_%g_%g_%g_%g_%g_%g_%g', azmin, azmax, elmin, elmax, azacc, elacc, numangles, dist_r_deg);
    filename_nodot = strrep(filename, '.', 'd');
    filename_nodot_nominus = strrep(filename_nodot, '-', 'm');

    %% check if this matrix is already available in persistent memory
    persistent storedFilename
    persistent storedCIMx
    persistent storedCIMy
    
    if ~isempty(storedFilename) && strcmp(storedFilename, filename_nodot_nominus) && ~forcerecalc
        % if a matrix has been stored, it's the same as the requested one,
        % and no forced recalculation has been requested -> return the
        % stored matrix
        xtarget_all = storedCIMx;
        ytarget_all = storedCIMy;
        if verbose
            fprintf('suitable CIM found in memory.\n');
        end
    elseif ~forcerecalc && exist(fullfile(CIMfolder, [filename_nodot_nominus '.cim']), 'file')
        %% if a file already exists, load it
        if verbose
            fprintf('suitable CIM found in file...');
        end
        load(fullfile(CIMfolder, [filename_nodot_nominus '.cim']), 'xtarget_all', 'ytarget_all', 'actualdist', 'actualangle', '-mat');
        if verbose
            fprintf('loaded.\n');
        end
    else
        if verbose
            fprintf('no suitable CIM found. Calculating...');
        end
        %% create additional variables
        [ele_2d, azi_2d]  = ndgrid(elmax:-elacc:elmin, azmin:azacc:azmax);  % NxM; elevation is reversed in images!, ele_2d(i, j) is the elevation of im(i, j)
        [M, N]       = size(ele_2d);                % N: image size along first dimension (=ele)
        [x_2d, y_2d] = ndgrid(1:N, 1:M);            % NxM, x is the first image dimension (=ele)

        elecorr_2d = 1 ./ cosd(ele_2d);             % NxM; correction factor for horizontal stretching at high elevations

        angles = linspace(0, 360, numangles+1);     % 1xnumangles; the angles to be tested (in degrees)
        angles = angles(1:end-1);                   % remove 360 deg
        dist_az_deg = cosd(angles) * dist_r_deg;    % 1xnumangles; distance along azimuth, in degrees
        dist_el_deg = sind(angles) * dist_r_deg;    % 1xnumangles; distance along elevation, in degrees

        %% calculate x-y-coordinates for each comparison
        xtarget_all = nan(N, M, numangles); % prealloc
        ytarget_all = xtarget_all; % prealloc
        actualdist  = xtarget_all; % prealloc
        actualangle = xtarget_all; % prealloc

        for i = 1:length(angles)
            %% calculate pixel distance and target pixels
            dist_x = - round(dist_el_deg(i) * ones(N, M) / elacc); % This has to be SUBTRACTED for elevation. TODO: This should be calculated from outside information, e.g. ele2
            dist_y = dist_az_deg(i) * ones(N, M) / azacc;
            dist_y = round(dist_y .* elecorr_2d);
            
            xtarget = x_2d + dist_x;
            ytarget = y_2d + dist_y;

            %% remove out-of-range target indices
            invalid = xtarget>N | xtarget<1 | ytarget>M | ytarget<1;
            xtarget(invalid) = NaN;
            ytarget(invalid) = NaN;

            %% add this angle's results to the global index file
            xtarget_all(:, :, i) = xtarget;
            ytarget_all(:, :, i) = ytarget;
                       
            %% for debugging, calculate the actual distance and actual angle between the sampled pixels
            actualdist(:, :, i)  = sqrt((dist_y./elecorr_2d*azacc).^2 + (dist_x*elacc).^2); 
            actualangle(:, :, i) = atan2d(dist_y./elecorr_2d*azacc, dist_x*elacc);
        end
        %% and save the results
        save(fullfile(CIMfolder, [filename_nodot_nominus '.cim']), 'xtarget_all', 'ytarget_all', 'actualdist', 'actualangle');
        if verbose
            fprintf('new CIM calculated and saved.\n');
        end
    end
    storedFilename = filename_nodot_nominus;
    storedCIMx = xtarget_all;
    storedCIMy = ytarget_all;
    
end % subfunction


