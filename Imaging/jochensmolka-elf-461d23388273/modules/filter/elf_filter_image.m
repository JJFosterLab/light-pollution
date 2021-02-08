function im_filt = elf_filter_image(im, azi, ele, hw, azi_filt, ele_filt, im_corr, usegpu, verbose, ah, metainfo)
% im_filt = elf_filter_image(im, azi, ele, hw, azi_filt, ele_filt, im_corr, verbose, ah)
%
% Filters a given image in equirectangular projection with Gaussian filters.
% This function corrects the kernel width for elevation, but uses a symmetric kernel for faster processing.
%
% Inputs: 
% im                - image projected into a equirectangular (azimuth / elevation) projection
% azi/ele           - in degrees, y/x vectors into im; (ele will usually be a reversed vector, starting at the highest value
%                                               because Matlab images have their origin at the top left corner)
% hw                - filter halfwidth in degrees, (FWHM, std sigma of the Gaussian is approx. hw / 2.35) (default: 1)
% azi_filt/ele_filt   - target y/x vectors (default: same range as azi/ele,
%                     with a resolution of hw/2) (ele should be reversed)
% im_corr           - correction image with the correct resolution and filtering to
%                     correct for edge effects (this is usually created from a white image
%                     filtered in the same way)
% verbose
% ah                - axes handle for display of the filtered image
%
% Output:
% im_filt           - Image of class double (but values not between 0/1)
% 
%
% Uses: None
% Used by: elf_filter

%% Check inputs
if nargin < 4, hw = 1; end
if nargin < 6, azi_filt = min(azi) : hw / 2 : max(azi); ele_filt = max(ele) : - hw / 2 : min(ele); end
if nargin < 7 || isempty(im_corr)
    im_corr{1} = ones(length(ele_filt), length(azi_filt), size(im, 3));
    im_corr{2} = ones(length(ele_filt), length(azi), size(im, 3));
end
if nargin < 8, usegpu = false; end
if nargin < 9, verbose = 0; end

if min(size(ele))>1 || min(size(azi))>1
    error('ele and azi vectors should be 1-dimensional');
end
if length(azi) ~= size(im, 2) || length(ele) ~= size(im, 1)
    error('Image dimensions must be compatible with the lengths of ele and azi');
end
if min(size(ele_filt))>1 || min(size(azi_filt))>1
    error('ele_filt and azi_filt vectors should be 1-dimensional');
end
if length(azi_filt) ~= size(im_corr{1}, 2) || length(ele_filt) ~= size(im_corr{1}, 1)
    error('Image dimensions must be compatible with the lengths of ele and azi');
end
if any(~ismember(azi_filt, azi))
    error('In this version, all elements of azi_filt MUST be in azi.'); %% TODO
end
if any(~ismember(ele_filt, ele))
    error('In this version, all elements of ele_filt MUST be in ele.'); %% TODO
end   
        
%% transform parameters
sigma_deg   = hw / (2*sqrt(2*log(2)));  % sigma in degrees, std of the Gaussian (FWHM is ~2.35*sigma)
if usegpu
    im          = gpuArray(double(im));
else
    im          = double(im);
end
g           = @(x,y,c) exp( -(x.^2/(2*sigma_deg^2) + y.^2/(2*(sigma_deg*c)^2)) ); % Gaussian function to create filter kernels, takes input in degrees and an elevation correction factor c
cnum        = size(im, 3);              % number of channels in input image
dpp_a       = median(abs(diff(azi)));   % degrees per pixel along azimuth of original image 
dpp_e       = median(abs(diff(ele)));   % degrees per pixel along azimuth of original image 
maxval      = max(im(:));               % max value, to normalise to 1 for plotting

%% calculate filtered image

% for c = 1:cnum
%     tsf{c} = fft2(im(:, :, c));
% end

if usegpu
    temp1 = zeros(length(ele_filt), length(azi_filt), cnum, 'gpuArray');   % preallocate image for vertical contrasts
else
    temp1 = zeros(length(ele_filt), length(azi_filt), cnum);   % preallocate image for vertical contrasts
end

parfor e = 1:length(ele_filt)                                       % for each elevation, starting from top
    % a) calculate symmetric Gaussian filter kernel, stretched for elevation, and cut off just after 3*sigma distance from centre   
    %    This kernel has to have the same angular resolution as the original image
    corr        = 1 / cosd(ele_filt(e));                         % correction factor for elevation in equirectangular projection
    %%FIXME: Are these the wrong way around?
    kernazi     = floor(-3*sigma_deg) : dpp_a : ceil(3*sigma_deg);    % This is always the same, but is calculated inside the loop for ease of understanding
    kernele     = floor(max([-90 -3*corr*sigma_deg])) : dpp_e : ceil(min([90 3*corr*sigma_deg])); % This changes with elevation
    [KA, KE]    = meshgrid(kernazi, kernele);
    gk          = g(KA, KE, corr);                              % Gaussian kernel, only up to a distance of 3 * sig (corrected for ele)
    gk          = gk' / sum(gk(:));                              % normalise

    % b) cut a strip 6 * sig wide from the image
    lowrow      = sub_findfirstbelow(ele, max([ele_filt(e) - 3 * sigma_deg min(ele)]));  % first row below ele - 3 * sigma
    highrow     = sub_findfirstabove(ele, min([ele_filt(e) + 3 * sigma_deg max(ele)]));  % first row above ele + 3 * sigma
    rows_v      = min([highrow lowrow]):max([highrow lowrow]);  % numbers of all rows between highrow and lowrow
    roweles     = ele(rows_v);
    [~, cols]   = ismember(azi_filt, azi);                       % NOTE: This assumes that all members of azi_filt are in azi already


% Superslow version  
%         thisstrip = im(rows_v, :, :);
%         filtstrip = imfilter(thisstrip, gk', 'same');           % gk has to be flipped to accommodate reversed image coordinates
%         im_filt(e, :, :) = filtstrip(sub_findnearest(roweles, ele_filt(e)), cols, :); % NOTE: This assumes that all members of ele_filt are in ele already

% currently best version
    for c = 1:cnum % for each channel
        thisstrip = squeeze(im(rows_v, :, c));
        filtstrip = imfilter(thisstrip, gk, 'same');           % gk has to be flipped to accommodate reversed image coordinates
        temp1(e, :, c) = filtstrip(sub_findnearest(roweles, ele_filt(e)), cols); % NOTE: This assumes that all members of ele_filt are in ele already
    end

%     if verbose
%         image(im2uint16(temp1./im_corr{1}/maxval), 'Parent', ah);
%         axis(ah, 'image', 'off');
%         drawnow
%     end
end
if verbose
    elf_plot_image(im2uint16(temp1./im_corr{1}/maxval), metainfo, ah, '', metainfo.linims);
    drawnow
end

im_filt = gather(temp1./im_corr{1});

end % main function

%% subfunctions
function pos = sub_findfirstbelow(mat, el)
    % finds the elements in mat nearest to but smaller than el
    mat(mat>el) = Inf;
    [~, pos] = min(abs(mat-el));
end

function pos = sub_findfirstabove(mat, el)
    % finds the elements in mat nearest to but larger than el
    mat(mat<el) = Inf;
    [~, pos] = min(abs(mat-el));
end

function pos = sub_findnearest(mat, el)
    % finds the elements in mat nearest to el
    [~, pos] = min(abs(mat-el));
end




















