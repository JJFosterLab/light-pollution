function [res, contim] = space_analysis(para, im_filt, verbose)
% space_analysis calculates spatial descriptors (i.e. all vertical and horizontal luminance and colour contrasts) in a prefiltered image.
%   For ideal sampling, the spacing of sampling points in azi/ele should be half the hwlfwidth (FWHM) of the Gaussian used for filtering. The
%   horizontal contrasts are corrected for elevation distortion, and the means are weighted to take this into account.
%
%   Example:
%   [res, contrim_v] = space_analysis(para, im_filt[, verbose])  
%
% Inputs: 
% para              - ELF parameter structure, used to extract azi/ele vectors, filter parameters, etc.
% im                - 1x2 cell of filtered images, 
%                       {1}: MxNxC double, projected into a equirectangular (azimuth / elevation) projection
%                       {2}: MxPxC double, with full resolution along azimuth to allow for elevation correction
% [verbose]         - logical, triggering display of bin boundaries (default: false)
%
% Outputs:
% res               - results structure containing all measured contrasts
% contrim           - 3 x sc x 2 x 3 structure containing a number of contrast images
%
% Uses:       space_analysis_channel
% Used by:    elf_analysis
% Call stack: elf_main2_descriptors -> space_analysis -> space_analysis_channel
% See also:   elf_main2_descriptors, elf_analysis, space_analysis_channel

%% check input
if nargin < 3, verbose = false; end

%% Filter image
if verbose
    fprintf('Performing spatial analysis with bin boundaries [%s].\n', num2str(para.ana.spatialbins));
end

%% FIXME: im{2} is now unused
numangles = 24;
for sc = 1:length(para.ana.scales_deg)  % for each scale
    % space_analysis_channel(im, aerange, aeacc, numangles, dist_r_deg, channeltype, meantype, meanthr, bins, verbose)
    [res.lum(:, :, sc), res.lumCS(:, sc), res.lumhist(:, :, :, sc), res.lumfft{sc}, res.lumprof{sc}] = ...
        space_analysis_channel(im_filt{sc}, [-90 90 -90 90], [.1 .1]*para.ana.scales_deg(sc), numangles, para.ana.scales_deg(sc), 'lum', para.ana.spatialmeantype, para.ana.spatialmeanthr, para.ana.spatialbins, verbose);
    [res.rg(:, :, sc), res.rgCS(:, sc), res.rghist(:, :, :, sc), res.rgfft{sc}, res.rgprof{sc}] = ...
        space_analysis_channel(im_filt{sc}, [-90 90 -90 90], [.1 .1]*para.ana.scales_deg(sc), numangles, para.ana.scales_deg(sc), 'rg', para.ana.spatialmeantype, para.ana.spatialmeanthr, para.ana.spatialbins, verbose);
    [res.gb(:, :, sc), res.gbCS(:, sc), res.gbhist(:, :, :, sc), res.gbfft{sc}, res.gbprof{sc}] = ...
        space_analysis_channel(im_filt{sc}, [-90 90 -90 90], [.1 .1]*para.ana.scales_deg(sc), numangles, para.ana.scales_deg(sc), 'gb', para.ana.spatialmeantype, para.ana.spatialmeanthr, para.ana.spatialbins, verbose);
    [res.rb(:, :, sc), res.rbCS(:, sc), res.rbhist(:, :, :, sc), res.rbfft{sc}, res.rbprof{sc}] = ...
        space_analysis_channel(im_filt{sc}, [-90 90 -90 90], [.1 .1]*para.ana.scales_deg(sc), numangles, para.ana.scales_deg(sc), 'rb', para.ana.spatialmeantype, para.ana.spatialmeanthr, para.ana.spatialbins, verbose);
end
    
res.binmean = (para.ana.spatialbins(1:end-1)+para.ana.spatialbins(2:end))/2;
temp        = linspace(0, 360, numangles+1);
res.angles  = temp(1:end-1);

    