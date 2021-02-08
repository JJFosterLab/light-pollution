function [meancontr_h, meancontr_v, binmean, contrhist_h, contrhist_v, contrim] = space_analysis_channel_old(im, azi, ele, azifull, channeltype, meantype, meanthr, bins, verbose, info)
% space_analysis_channel calculates one particular luminance or colour contrast in a prefiltered image.
%   For ideal sampling, the spacing of sampling points in azi/ele should be half the hwlfwidth (FWHM) of the Gaussian used for filtering. The
%   horizontal contrasts are corrected for elevation distortion, and the means are weighted to take this into account.
%
%   Example: 
%   [meancontr_h, meancontr_v, binmean] = space_analysis_channel(im, azi, ele, azifull, type[, bins, verbose])
%
% Inputs:
% im                - 1x2 cell of filtered images, 
%                       {1}: MxNxC double, projected into a equirectangular (azimuth / elevation) projection
%                       {2}: MxPxC double, with full resolution along azimuth to allow for elevation correction
% azi/ele           - 1xN / 1xM doubles, in degrees, y/x-vectors into im{1}, 
%                       ele will usually be a reversed vector, starting at the highest value 
%                       because Matlab images have their origin at the top left corner
% azifull           - 1xP double, y-vector into im{2}
% channeltype       - string, can be 'lum', 'rg', 'yb', 'gb' (If channels other than RGB are used, this has to be modified!)
% [meantype]        - string, can be 'mean', 'rms', 'perc' (default 'rms')
% [meanthr]         - double, indicates that the mean is only calculated for contrast values greater of equal this, e.g. 5 for contrasts >=5% (default 0).
% [bins]            - 1xB double, bin boundaries in degrees (default: [-90 -50 -10 10 50 90])
% [verbose]         - logical, triggering plotting of all measured contrasts. For 1-degree scale, this takes a LONG time! (default: false)
% [info]            - info file for images; only required if verbose == true
% 
% Outputs:
% meancontr_h       - (B-1)x1 double, mean horizontal contrast for each bin
% meancontr_v       - (B-1)x1 double, mean vertical contrast for each bin
% binmean           - (B-1)x1 double, mean elevation for each bin, can be used for plotting
% contrhist_h       - 101x(B-1) double, histogram of horizontal contrasts in bins of 1 degree
% contrhist_v       - 101x(B-1) double, histogram of vertical contrasts in bins of 1 degree
%
% Uses:       None
% Used by:    space_analysis
% Call stack: elf_main2_descriptors -> elf_analysis -> space_analysis -> space_analysis_channel
% See also:   elf_main2_descriptors, elf_analysis, space_analysis

%% check inputs
if nargin < 9, verbose = false; end
if nargin < 8, bins    = [-90 -50 -10 10 50 90]; end      
if nargin < 7, meanthr = 0; end
if nargin < 6, meantype = 'rms'; end
        
%% spatial analysis

fprintf('Calculating contrasts in the %s channel using %s-type means (threshold %g%%).\n', channeltype, meantype, meanthr);

% re-calculate some basic parameters
ppdfull      = 1 / median(diff(azifull));  % pixels per degree in original image (im_proj)
ppdnew       = 1 / median(diff(azi));      % pixels per degree in filtered image (im_filt)
elecorr      = 1 ./ cosd(ele);             % correction factor for horizontal stretching at high elevations
stepsize_pix = elecorr * ppdfull / ppdnew; % distance between nearest neighbours, e.g. step 1/2 deg for 1 deg halfwidth filters at horizon
stepsize_pix(isinf(stepsize_pix)) = size(im{2}, 2)-1;
hw           = 5;                          % after calculation of elementary horizontal/vertical contrasts, they will be filtered along each 
                                           % column/row to achieve detection of vertical/horizontal edges. This hw is the half-width of that filter in degrees. 

switch channeltype
    case 'lum' % calculate Michelson contrast: |L1-L2| / (L1+L2)
        % 1. calculate vertical contrasts
        L           = sum(im{1}, 3);  % Luminance is simply the sum of all other channels
        mich_v      = diff(L, 1, 1) ./ (L(1:end-1, :) + L(2:end, :));   % Michelson contrast along first dimension (= elevation, reversed)
        mich_v_edge = elf_filter_image_1D(mich_v, azi, ele(1:end-1), hw, 1, [], 0); %NEW STEP (indexing not quite accurate)
        vert_contr  = abs(mich_v_edge);
        
        % 2. calculate horizontal contrasts line by line, because the distance between nearest neighbours increases with elevation
        L2          = sum(im{2}, 3);  % Luminance is simply the mean of all other channels
        mich_h      = nan(size(L2));  % prealloc
        for e = 1:length(ele)
            ss           = round(stepsize_pix(e));  % distance in pixels between "nearest neighbours" for this elevation
            L            = L2(e, :);    % whole row
            mich_h_row   = (L(1+ss:end)-L(1:end-ss)) ./ (L(1+ss:end)+L(1:end-ss)); % Michelson contrast along second dimension  (= azimuth)
            % add ss nan-pixels to fill up before and after
            mich_h(e, :) = [nan(1, floor(round(ss)/2)) mich_h_row nan(1, ceil(round(ss)/2))];
        end
        mich_h_edge      = elf_filter_image_1D(mich_h, azifull, ele, hw, 2, [], 0); %NEW STEP (indexing not quite accurate)
        hor_contr        = abs(mich_h_edge);
            
    case {'rg', 'yb', 'gb'}
        % separate channels
        R           = im{1}(:, :, 1);
        G           = im{1}(:, :, 2);
        B           = im{1}(:, :, 3);
        Y           = (R + G) / 2;
        R2          = im{2}(:, :, 1);
        G2          = im{2}(:, :, 2);
        B2          = im{2}(:, :, 3);
        Y2          = (R2 + G2) / 2;

        switch channeltype
            case 'rg' % calculate (R-G)/(R+G) for each pixel, then diff between pixels
                H           = (R - G) ./ (R + G);      % between -1 and 1
                H2          = (R2 - G2) ./ (R2 + G2);
            case 'yb' % calculate (Y-B)/(Y+B) for each pixel, then diff between pixels
                H           = (Y - B) ./ (Y + B);
                H2          = (Y2 - B2) ./ (Y2 + B2);
            case 'gb' % calculate (G-B)/(G+B) for each pixel, then diff between pixels
                H           = (G - B) ./ (G + B);
                H2          = (G2 - B2) ./ (G2 + B2);
            otherwise
        end
                
        % 1. calculate vertical contrasts
        mich_v      = diff(H, 1, 1) / 2;   % Michelson contrast along first dimension (= elevation, reversed)
        mich_v_edge = elf_filter_image_1D(mich_v, azi, ele(1:end-1), hw, 1, [], 0); %NEW STEP (indexing not quite accurate)
        vert_contr  = abs(mich_v_edge);

        % 2. calculate horizontal contrasts line by line, because the distance between nearest neighbours increases with elevation
        mich_h      = nan(size(H2));  % prealloc
        for e = 1:length(ele)
            ss           = round(stepsize_pix(e));  % distance in pixels between "nearest neighbours" for this elevation
            H            = H2(e, :);    % whole row
            mich_h_row   = (H(1+ss:end)-H(1:end-ss)) / 2; % Michelson contrast along second dimension  (= azimuth)
            % add ss nan-pixels to fill up before and after
            mich_h(e, :) = [nan(1, floor(round(ss)/2)) mich_h_row nan(1, ceil(round(ss)/2))];
        end
        mich_h_edge      = elf_filter_image_1D(mich_h, azifull, ele, hw, 2, [], 0); %NEW STEP (indexing not quite accurate)
        hor_contr        = abs(mich_h_edge);

        
               
    otherwise
        error('Unknown spatial comparison type.');
end

%% create contrast images
contrim{1, 1} = 1-abs(mich_v);
contrim{1, 2} = 1-abs(mich_v_edge);
contrim{2, 1} = 1-abs(mich_h);
contrim{2, 2} = 1-abs(mich_h_edge);

%% and plot them
if verbose
    figure; clf; elf_support_maxfig; drawnow;
    subplot(2, 2, 1);imagesc(1-abs(mich_v));colormap(gray);colorbar;axis square;title('vert cont');
    subplot(2, 2, 2);imagesc(1-abs(mich_v_edge));colormap(gray);colorbar;axis square; title('filt vert cont = hor edges');
    subplot(2, 2, 3);imagesc(1-abs(mich_h));colormap(gray);colorbar;axis square; title('hor cont');
    subplot(2, 2, 4);imagesc(1-abs(mich_h_edge));colormap(gray);colorbar;axis square; title('filt hor cont = vert edges');
end
        
%% Average over bins
meancontr_h = nan(length(bins)-1, 1);
meancontr_v = nan(length(bins)-1, 1);
contrhist_h = nan(101, length(bins)-1);
contrhist_v = nan(101, length(bins)-1);

c = meanthr/100; % contrast threshold

for bin = 1:length(bins)-1
    % select the rows that compare entirely WITHIN the bin
    mindeg = bins(bin);
    maxdeg = bins(bin+1);
    rows_v = ele(1:end-1) >= mindeg & ele(1:end-1) <= maxdeg & ele(2:end) >= mindeg & ele(2:end) <= maxdeg;
    rows_h = ele >= mindeg & ele <= maxdeg;

    % extract rows
    thiscontr_v = vert_contr(rows_v, :);
    thiscontr_h = hor_contr(rows_h, :);
    thiscontr_h_weights = 1 ./ repmat(elecorr(rows_h)', [1, size(thiscontr_h, 2)]);
    
    % calculate averages
    switch meantype
        case 'mean'
            meancontr_v(bin)    = mean(thiscontr_v(thiscontr_v>=c));
            %meancontr_h(bin)    = mean(thiscontr_h(thiscontr_h>=c)); %OLD
            meancontr_h(bin)    = sum(thiscontr_h(thiscontr_h>=c) .* thiscontr_h_weights(thiscontr_h>=c)) ./ sum(thiscontr_h_weights(thiscontr_h>=c));
        case 'perc'
            meancontr_v(bin)    = sum(thiscontr_v>=c)/length(thiscontr_v(:));
            %meancontr_h(bin)    = sum(thiscontr_h>=c)/length(thiscontr_h(:)); %OLD
            meancontr_h(bin)    = sum(thiscontr_h>=c .* thiscontr_h_weights(thiscontr_h>=c)) ./ sum(thiscontr_h_weights); % UNSURE whether this is correct
        case 'rms'
            meancontr_v(bin)    = rms(thiscontr_v(thiscontr_v>=c));
            %meancontr_h(bin)    = rms(thiscontr_h(thiscontr_h>=c)); %OLD
            meancontr_h(bin)    = sqrt( sum( thiscontr_h(thiscontr_h>=c).^2 .* thiscontr_h_weights(thiscontr_h>=c)) ./ sum(thiscontr_h_weights(thiscontr_h>=c)) );
        otherwise
            error('Unknown mean type: %s', meantype);
    end
    
    if meancontr_v(bin)==0 || isnan(meancontr_v(bin)), meancontr_v(bin) = 1/length(thiscontr_v(:)); end % remove zeros
    if meancontr_h(bin)==0 || isnan(meancontr_h(bin)), meancontr_h(bin) = 1/length(thiscontr_h(:)); end % remove zeros
    
    % also save contrast histograms
    contrhist_v(:, bin) = histc(thiscontr_v(:), 0:.01:1);
    contrhist_h(:, bin) = histc(thiscontr_h(:), 0:.01:1);
end

binmean = (bins(1:end-1)+bins(2:end))/2;
      
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

% % horizontal *not functional. Solution here must be much more complex due
% to elevation distortion.
% a = [hor_contr>=c;false(1, imsize)] | [false(1, imsize);hor_contr>=c];
% outim{2} = elf_io_correctdng(im{1}, info);
% for i = 1:imsize
%     for j = 1:imsize
%         if ~a(i, j)
%             bw = mean(outim{2}(i, j, :));
%             outim{2}(i, j, :) = bw;
%         end
%     end
% end
    
%% plot which contrasts have been included %BROKEN
if verbose && meanthr>0
    switch channeltype
        case 'lum'
            ct = 1;
        case 'rg'
            ct = 2;
        case 'yb'
            ct = 3;
        case 'gb'
            ct = 4;
    end
    figure(4+ct+4*(imsize<300)); 
    image(im2);
    axis image;
    title(channeltype)
end


%% plot contrasts
% if verbose
%     % Do some scaling
%     maxval      = max(im(:));
%     vert_contr  = vert_contr / maxval * 10;
%     hor_contr   = hor_contr  / maxval * 10;
%     meancontr_v = meancontr_v / maxval * 100;
%     meancontr_h = meancontr_h / maxval * 100;
%     
%     [X,Y]=meshgrid(1:length(azi), 1:length(ele));
%     
%     figure(9);clf;image(im2uint16(im));axis image;hold on;
%     
%     % draw vertical lines
%     xxx1 = Y(:, 1:end-1);
%     xxx2 = Y(:, 2:end);
%     yyy1 = X(:, 1:end-1);
%     yyy2 = X(:, 2:end);
%     allx = [xxx1(:) xxx2(:)]';
%     ally = [yyy1(:)+0.1 yyy2(:)-0.1]';
%     alllw = (mean(vert_contr, 3)+1)';
%     alllw = alllw(:);
%     for i = 1:size(allx, 2)
%         plot(allx(:, i), ally(:, i), 'k', 'linewidth', alllw(i));
%     end
%     
%     % draw horizontal lines
%     xxx1 = Y(1:end-1, :);
%     xxx2 = Y(2:end, :);
%     yyy1 = X(1:end-1, :);
%     yyy2 = X(2:end, :);
%     allx = [xxx1(:)+0.1 xxx2(:)-0.1]';
%     ally = [yyy1(:) yyy2(:)]';
%     alllw = (mean(hor_contr, 3)+1)';
%     alllw = alllw(:);
%     for i = 1:size(allx, 2)
%         plot(allx(:, i), ally(:, i), 'k', 'linewidth', alllw(i));
%     end
% 
%     figure(10);clf;hold on;
% 
%     for i = 1:length(binmean)
%         plot([1-meancontr_h(i) 1+meancontr_h(i)], [binmean(i) binmean(i)], 'k');
%         plot([1 1], [binmean(i)-meancontr_v(i) binmean(i)+meancontr_v(i)], 'k'); 
%         rectangle('Position', [1-meancontr_h(i) binmean(i)-meancontr_v(i) 2*meancontr_h(i) 2*meancontr_v(i)], 'Curvature', [1 1]);
%     end
%     
%     axis([-90 90 -90 90]);
%     axis equal
% end







