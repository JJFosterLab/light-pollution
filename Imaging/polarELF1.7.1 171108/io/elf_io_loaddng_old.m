function [im, raw_dem] = elf_io_loaddng(fullfilename, verbose)
% ELF_IO_LOADDNG reads DNG (digital negative) image files
%   Linearisation and "as shot" white balancing, and demosaicing are
%   performed, but no colour space conversion. 
%   Algorithm adapted from Rob Sumner (2014) "Processing RAW Images in
%   MATLAB", http://users.soe.ucsc.edu/~rcsumner/rawguide/RAWguide.pdf
%   The DNG files can be created from any raw camera format with the Adobe
%   DNG converter. In the Converter, click "Change preferences", select
%   "Custom Compatibility" and make sure that "Uncompressed" is checked and
%   that "Linear (demosaiced)" is unchecked.
%
% Inputs: 
% fullfilename  - has to include the full path to the image file
% 
% Outputs:
% im            - output image array
%
% See also elf_imfinfo, elf_imread, elf_io_correctdng.

if nargin < 1 % for testing only
    fullfilename = 'F:\All data\e17 Lagoon Hapto Rubble 0930\_VE21798.dng';
end

%TODO: Add verbose output

%% check inputs
if nargin < 2, verbose = false; end

%% parameters
sensorAlignment = 'rggb'; % sensorAlignment can be 'gbrg'/'grbg'/'bggr'/'rggb'

%% disable warnings
warning('off', 'MATLAB:tifflib:TIFFReadDirectory:libraryWarning');
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
% Warning: TIFF library error - 'JPEGPreDecode:  Improper JPEG component
% count.' - file may be corrupt.  <--- This seems to mean the dng is not
% uncompressed.

%% Read CFA (colour filter array) image into Matlab
t = Tiff(fullfilename, 'r');       % create a tif-object
offsets = getTag(t,'SubIFD');
setSubDirectory(t,offsets(1));     % set the first SubIFD as the current one (this is where the raw information is)
raw = read(t);                     % Create variable 'raw', the Bayer CFA data
close(t);                          % close the tif-object
meta_info = imfinfo(fullfilename); % read exif information

% Crop to only valid pixels
x_origin    = meta_info.SubIFDs{1}.ActiveArea(2)+1; % +1 due to MATLAB indexing
width       = meta_info.SubIFDs{1}.DefaultCropSize(1);
y_origin    = meta_info.SubIFDs{1}.ActiveArea(1)+1;
height      = meta_info.SubIFDs{1}.DefaultCropSize(2);
raw         = double(raw(y_origin:y_origin+height-1,x_origin:x_origin+width-1));

%% Linearisation
% first, check whether a linearisation table is needed
if isfield(meta_info.SubIFDs{1}, 'LinearizationTable') && ~isempty(meta_info.SubIFDs{1}.LinearizationTable)
    raw = meta_info.SubIFDs{1}.LinearizationTable(raw);
    % error('Linearisation Table found. The image needs to be linearised, but this has not been coded yet.');
end

% now, scale using the saved BlackLevel and WhiteLevel
black = meta_info.SubIFDs{1}.BlackLevel(1);
saturation = meta_info.SubIFDs{1}.WhiteLevel;
lin_bayer = (raw-black)/(saturation-black);
lin_bayer = max(0,min(lin_bayer,1)); % 0 - 1

%% White balancing
% Apply the "As shot" white balance to correct for sensitivity differences in R, G and B pixels
wb_multipliers = (meta_info.AsShotNeutral).^-1;
wb_multipliers = wb_multipliers/wb_multipliers(2); % normalise to green channel
mask = sub_wbmask(size(lin_bayer,1), size(lin_bayer,2), wb_multipliers(1), wb_multipliers(3), sensorAlignment);
balanced_bayer = lin_bayer .* mask;

%% Demosaicing
% Note: Scaling the image up to a maximum of 2^16 (as is done in Sumner's
% algorithm) makes very little difference (in example, no more than 27 of 2^16) during demosaicing.
im = demosaic(uint16(balanced_bayer * 2^16), sensorAlignment); %produces a MxNx3 linear uint16 output image

%elf_plot_image(elf_io_correctdng(im, meta_info));

%% returns raw demosaiced, if requested
if nargout>1
    raw_dem = demosaic(uint16(raw), sensorAlignment);
end
[mean(raw(:)) std(raw(:))]
end %main

%% subfunctions

function colormask = sub_wbmask(m, n, r_scale, b_scale, sensorAlignment)
    % COLORMASK = sub_wbmask(M, N, R_SCALE, B_SCALE)
    %
    % Makes a white-balance multiplicative mask for an RGGB image of size m-by-n with
    % white balance scaling values R_SCALE, G_SCALE=1, and B_SCALE.
    % sensorAlignment can be 'gbrg'/'grbg'/'bggr'/'rggb'
    if nargin < 5
        sensorAlignment = 'rggb';
    end

    %TODO: 'gbrg' and 'grbg' have not been tested
    colormask = ones(m,n);
    switch sensorAlignment
        case 'rggb'
            colormask(1:2:end,1:2:end) = r_scale;
            colormask(2:2:end,2:2:end) = b_scale;
        case 'gbrg'
            colormask(1:2:end,2:2:end) = b_scale;
            colormask(2:2:end,1:2:end) = r_scale;
            warning('not tested');
        case 'grbg'
            colormask(1:2:end,2:2:end) = r_scale;
            colormask(2:2:end,1:2:end) = b_scale; 
            warning('not tested');       
        case 'bggr'
            colormask(1:2:end,1:2:end) = b_scale;
            colormask(2:2:end,2:2:end) = r_scale;
        otherwise
            error(['Unknown sensor alignment: ' sensorAlignment]);
    end

end






