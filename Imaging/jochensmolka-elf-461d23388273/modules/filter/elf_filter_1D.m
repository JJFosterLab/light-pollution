function [im_filt, para] = elf_filter_1D(para, im, verbose, ah, metainfo)
% [im_filt, para] = elf_filter_1D(para, im, azi, ele, verbose, ah)
%
% Filters images along one dimension.
%
% Inputs: 
% im - prefiltered image in a equirectangular (azimuth / elevation) projection
% azi/ele - x/y vectors into im; ele should be a reversed vector
% metainfo is used for colour correction before plotting
%
% Uses: None

hw = 5; % half-width in degrees

%% First time this runs, create correction images
persistent im_corr;
if isempty(im_corr)
    im_corr = cell(2, 1);
%     if ~isempty(ah)
%         temp = text(0.5, 0.5, 'Calculating correction image...     ', 'parent', ah(sc), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle'); 
%         drawnow;
%     end
    im_corr{1} = elf_filter_image_1D(ones(size(im{1, 1})), para.azi_filt{1}, para.ele_filt{1}, hw, 1);
    im_corr{2} = elf_filter_image_1D(ones(size(im{1, 2})), para.azi, para.ele_filt{1}, hw, 2);
%     if ~isempty(ah)
%         delete(temp);
%         text(0.5, 0.5, 'Calculating correction image...done.', 'parent', ah(sc), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle'); 
%         drawnow;
%     end
end

%% Filter image
if verbose
    im_filt{1} = elf_filter_image_1D(im{1, 1}, para.azi_filt{1}, para.ele_filt{1}, hw, 1, im_corr{1}, verbose, [], metainfo);
    im_filt{2} = elf_filter_image_1D(im{1, 2}, para.azi, para.ele_filt{1}, hw, 2, im_corr{2}, verbose, [], metainfo);
else
    im_filt{1} = elf_filter_image_1D(im{1, 1}, para.azi_filt{1}, para.ele_filt{1}, hw, 1, im_corr{1});
    im_filt{2} = elf_filter_image_1D(im{1, 2}, para.azi, para.ele_filt{1}, hw, 2, im_corr{2});
end
    


    