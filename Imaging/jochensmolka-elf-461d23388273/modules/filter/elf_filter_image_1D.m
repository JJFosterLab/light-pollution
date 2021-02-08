function im_filt = elf_filter_image_1D(im, azi, ele, hw, dim, im_corr, verbose, ah, metainfo)
% im_filt = elf_filter_image_1D(im, azi, ele, hw, dim, im_corr, verbose, ah)
%
% Filters a given prefiltered image in equirectangular projection with Gaussian filters along one dimension.
% This function corrects the kernel width for elevation.
%
% Inputs: 
% im                - image projected into a equirectangular (azimuth / elevation) projection
% azi/ele           - in degrees, y/x vectors into im; (ele will usually be a reversed vector, starting at the highest value
%                                               because Matlab images have their origin at the top left corner)
% hw                - filter halfwidth in degrees, (FWHM, std sigma of the Gaussian is approx. hw / 2.35) (default: 1)
% dim               - dimension along which to filter (1 means filtering along x/azimuth, e.g. for later calculation of vertical contrasts)
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
% Used by: elf_filter_1D

%% Check inputs
if nargin < 4, hw = 1; end
if nargin < 5, dim = 2; end
if nargin < 6 || isempty(im_corr), im_corr = ones(size(im));end
if nargin < 7, verbose = false; end
if nargin < 8, metainfo = []; end

if min(size(ele))>1 || min(size(azi))>1
    error('ele and azi vectors should be 1-dimensional');
end
if length(azi) ~= size(im, 2) || length(ele) ~= size(im, 1)
    error('Image dimensions must be compatible with the lengths of ele and azi');
end
  
%% transform parameters
sigma_deg   = hw / (2*sqrt(2*log(2)));  % sigma in degrees, std of the Gaussian (FWHM is ~2.35*sigma)
im          = double(im);
cnum        = size(im, 3);              % number of channels in input image
dpp_a       = median(abs(diff(azi)));   % degrees per pixel along azimuth of original image 
dpp_e       = median(abs(diff(ele)));   % degrees per pixel along azimuth of original image 
maxval      = max(im(:));               % max value, to normalise to 1 for plotting

%% calculate filtered image
im_filt = zeros(size(im));        % preallocate image 

if dim == 1 % filter along azimuth, i.e. once for each elevation value
    g   = @(x,c) exp( -x.^2/(2*(sigma_deg*c)^2)); % 1D Gaussian function to create filter kernels along azimuth, takes input in degrees and an elevation correction factor c
    for e = 1:length(ele)
        % a) calculate symmetric Gaussian filter kernel, stretched for elevation, and cut off just after 3*sigma distance from centre   
        %    This kernel has to have the same angular resolution as the original image
        corr        = 1 / cosd(ele(e));                         % correction factor for elevation in equirectangular projection
        kernazi     = floor(max([-90 -3*corr*sigma_deg])) : dpp_a : ceil(min([90 3*corr*sigma_deg])); % This changes with elevation
        gk          = g(kernazi, corr);                         % Gaussian kernel, only up to a distance of 3 * sig (corrected for ele)
        gk          = gk / sum(gk(:));                          % normalise
    
        % b) filter each row (currently best version)
        for c = 1:cnum % for each channel
            im_filt(e, :, c) = imfilter(im(e, :, c), gk, 'same');           % gk has to be flipped to accommodate reversed image coordinates
        end
    end

elseif dim == 2 % filter along elevation, i.e. once for each azimuth value
    g   = @(x) exp( -x.^2/(2*sigma_deg^2) );    % 1D Gaussian function to create filter kernels along elevation, takes input in degrees
    % a) calculate symmetric Gaussian filter kernel, stretched for elevation, and cut off just after 3*sigma distance from centre   
    %    This kernel has to have the same angular resolution as the original image
    kernele     = floor(max([-90 -3*sigma_deg])) : dpp_e : ceil(min([90 3*sigma_deg])); % This is always the same, so can stay outside of loop
    gk          = g(kernele);                         % Gaussian kernel, only up to a distance of 3 * sig (corrected for ele)
    gk          = gk / sum(gk(:));                          % normalise
    
    for e = 1:length(azi)
        for c = 1:cnum % for each channel
            im_filt(:, e, c) = imfilter(im(:, e, c), gk', 'same');           % gk has to be flipped to accommodate reversed image coordinates
        end
    end
    
else
    error('dim has to be 1 or 2');
end

if verbose
    if isempty(metainfo)
        elf_plot_image(im2uint16(im/maxval));
        elf_plot_image(im2uint16(im_filt./im_corr/maxval));
        drawnow
    else
        elf_plot_image(im2uint16(im/maxval), metainfo, [], '', 1);
        elf_plot_image(im2uint16(im_filt./im_corr/maxval), metainfo, [], '', 1);
        drawnow
    end
end

im_filt = im_filt./im_corr;

end 


















