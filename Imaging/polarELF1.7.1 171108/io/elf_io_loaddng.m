function im_raw = elf_io_loaddng(fullfilename, sensorAlignment)
% ELF_IO_LOADDNG reads DNG (digital negative) image files
%   Images are demosaiced, and EXIF linearisation tables are applied, if applicable.
%   No further linearisation or white balancing or colour space conversion is performed.
%   See elf_calibrate_abssens for all further calibration.
%   Algorithm adapted from Rob Sumner (2014) "Processing RAW Images in
%   MATLAB", http://users.soe.ucsc.edu/~rcsumner/rawguide/RAWguide.pdf
%   The DNG files can be created from any raw camera format with the Adobe
%   DNG converter. In the Converter, click "Change preferences", select
%   "Custom Compatibility" and make sure that "Uncompressed" is checked and
%   that "Linear (demosaiced)" is unchecked.
%
% Inputs: 
% fullfilename      - has to include the full path to the image file
% sensorAlignment   - can be 'gbrg'/'grbg'/'bggr'/'rggb' (default 'rggb')
%
% Outputs:
% im_raw            - output image array (NxMx3 uint16)
%
% See also elf_imfinfo, elf_imread, elf_io_correctdng, elf_calibrate_abssens.

if nargin < 1 % for testing only
    fullfilename = 'F:\All data\e17 Lagoon Hapto Rubble 0930\_VE21798.dng';
end

%% check inputs
if nargin < 2 || isempty(sensorAlignment), sensorAlignment = 'rggb'; end % sensorAlignment can be 'gbrg'/'grbg'/'bggr'/'rggb'

%% disable warnings
warning('off', 'MATLAB:tifflib:TIFFReadDirectory:libraryWarning');
warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
% Warning: TIFF library error - 'JPEGPreDecode:  Improper JPEG component
% count.' - file may be corrupt.  <--- This seems to mean the dng is not
% uncompressed.

%% Read CFA (colour filter array) image into Matlab
t           = Tiff(fullfilename, 'r');       % create a tif-object
offsets     = getTag(t, 'SubIFD');
setSubDirectory(t, offsets(1));              % set the first SubIFD as the current one (this is where the raw information is)
raw         = read(t);                       % Create variable 'raw', the Bayer CFA data, as a uint16
close(t);                                    % close the tif-object
meta_info   = imfinfo(fullfilename);         % read exif information

%% Crop to only valid pixels
x_origin    = meta_info.SubIFDs{1}.ActiveArea(2)+1; % +1 due to MATLAB indexing
width       = meta_info.SubIFDs{1}.DefaultCropSize(1);
y_origin    = meta_info.SubIFDs{1}.ActiveArea(1)+1;
height      = meta_info.SubIFDs{1}.DefaultCropSize(2);
raw         = raw(y_origin:y_origin+height-1,x_origin:x_origin+width-1);

%% Here, it would be possible to extract the pixel values of unexposed chip pixels, which might help determine dark noise

%% Demosaicing
% Note: Scaling the image up to a maximum of 2^16 (as is done in Sumner's
% algorithm) makes very little difference (in example, no more than 27 of 2^16) during demosaicing.

im_raw = demosaic(raw, sensorAlignment); % produces a MxNx3 linear uint16 output image
% elf_plot_image(elf_io_correctdng(raw_dem, meta_info));

%% Linearisation
% check whether a linearisation table is needed
if isfield(meta_info.SubIFDs{1}, 'LinearizationTable') && ~isempty(meta_info.SubIFDs{1}.LinearizationTable)
    im_raw = meta_info.SubIFDs{1}.LinearizationTable(im_raw);
end

end %main






