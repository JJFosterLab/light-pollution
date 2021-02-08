function im = elf_imread(fullfilename)
% ELF_IMREAD reads image files of several different types
%   Currently works for tif, jpg, and dng, and should be ok for most other
%   non-raw formats. For raw formats, linearisation and demosaicing are
%   performed, but no white balance correction or colour space conversion.
%
% Inputs: 
% fullfilename  - has to include the full path to the image file
% 
% Outputs:
% im            - output image array
%
% See also elf_imfinfo, elf_io_loaddng.

if nargin < 1 % for testing only
    fullfilename = 'J:\Data and Documents\data\2014 VEPS test data\Hasthage_hovmansbygd_10mars1730\_D3X2378.dng';
end

[~,~,ext] = fileparts(fullfilename); % using info.Format does not work for raw files, as they are usually tif format
switch lower(ext(2:end))
    case {'tif', 'tiff', 'jpg', 'jpeg', 'bmp', 'gif', 'png', 'ppm'}
        im = imread(fullfilename);
    case 'dng'
        im = elf_io_loaddng(fullfilename);
    case 'nef'
        im = zeros(3, 3, 3); % just something to display a black image in ELF maingui
    case {'cr2', 'crw', 'kdc', 'arw', 'srf', 'sr2', 'bay', 'dcs', 'dcr', 'drf', 'k25', 'nrw', 'orf', 'pef', 'ptx', 'raw', 'rw2', 'rwl'}
        % these are the most common raw formats for Canon/Nikon/Casio/Sony/Kodak/Olympus/Pentax/Panasonic/Minolta cameras
        error('ELF currently does not process %s files, but it should be possible using dcraw. Create an entry in elf_imread.m');
    otherwise
        fprintf('Unknown extension %s. Attempting to read with Matlab''s imread function.\nTo disable this message, create an entry for this file type in elf_imread.m\n', ext);
        im = imread(fullfilename);
end



end %main

