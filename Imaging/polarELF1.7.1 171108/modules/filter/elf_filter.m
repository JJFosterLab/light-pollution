function [im_filt, para] = elf_filter(para, im, type, verbose, ah, metainfo, forcerecalc)
% [im_filt, para] = elf_filter(para, im, verbose, ah)
%
% Filters images.
%
% Inputs: 
% im - image projected into a equirectangular (azimuth / elevation) projection
% scales contains sigma in pixels, std of the Gaussian (FWHM is 2.35*sigma) (default: [2 10])
% type can be 'im'(sampled to hw/2) or 'scene' (filtered to hw/10)
% im_corr - correction image with the correct resolution and filtering to correct for edge effects
% Uses: None

if nargin<7
    forcerecalc = false;
end

%% Calculate target angular grid if necessary
switch type
    case 'im' % filter individual images
        if ~isfield(para, 'azi_filt')
            for sc = 1:length(para.ana.scales_deg)  % for each scale
                acc                         = para.ana.scales_deg(sc)/10; warning('DEBUGGING MODE!'); 
                azi_filt{sc}                = para.ana.filterazimin: acc:para.ana.filterazimax;
                ele_filt{sc}                = para.ana.filterelemax:-acc:para.ana.filterelemin;  
                azi_filt{sc} = round(azi_filt{sc}*100)/100; % this is apparently necessary to avoid numerical inaccuracies
                ele_filt{sc} = round(ele_filt{sc}*100)/100;              
            end
            para.azi_filt = azi_filt;
            para.ele_filt = ele_filt;
        else
            azi_filt = para.azi_filt;
            ele_filt = para.ele_filt;
        end
    case 'scene' % filter HDR scene images
        if ~isfield(para, 'azi_filt_scene')
            for sc = 1:length(para.ana.scales_deg)  % for each scale
                acc                         = para.ana.scales_deg(sc)/10;
                azi_filt{sc}                = para.ana.filterazimin: acc:para.ana.filterazimax;
                ele_filt{sc}                = para.ana.filterelemax:-acc:para.ana.filterelemin;
                azi_filt{sc}=round(azi_filt{sc}*100)/100; % this is apparently necessary to avoid numerical inaccuracies
                ele_filt{sc}=round(ele_filt{sc}*100)/100;
            end
            para.azi_filt_scene = azi_filt;
            para.ele_filt_scene = ele_filt;
        else
            azi_filt = para.azi_filt_scene;
            ele_filt = para.ele_filt_scene;
        end
end
para.azi = round(para.azi*100)/100; % this is apparently necessary to avoid numerical inaccuracies
para.ele = round(para.ele*100)/100; % this is apparently necessary to avoid numerical inaccuracies
para.ele2 = round(para.ele2*100)/100; % this is apparently necessary to avoid numerical inaccuracies

%% First time this runs, create correction images
persistent im_corr;
if isempty(im_corr) || forcerecalc
    im_corr = cell(length(para.ana.scales_deg), 1);
    for sc = 1:length(para.ana.scales_deg)  % for each scale
        if ~isempty(ah)
            temp = text(0.5, 0.5, 'Calculating correction image...     ', 'parent', ah(sc), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle'); 
            drawnow;
        end
        im_corr{sc} = elf_filter_image(ones(size(im)), para.azi, para.ele2, para.ana.scales_deg(sc), azi_filt{sc}, ele_filt{sc}, [], para.usegpu);
        if ~isempty(ah)
            delete(temp);
            text(0.5, 0.5, 'Calculating correction image...done.', 'parent', ah(sc), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle'); 
            drawnow;
        end
    end
end

%% Filter image
im_filt = cell(length(para.ana.scales_deg), 1);
for sc = 1:length(para.ana.scales_deg)  % for each scale
    if verbose
        im_filt{sc}     = elf_filter_image(im, para.azi, para.ele2, para.ana.scales_deg(sc), azi_filt{sc}, ele_filt{sc}, im_corr(sc), para.usegpu, verbose, ah(sc), metainfo);
    else
        im_filt{sc}     = elf_filter_image(im, para.azi, para.ele2, para.ana.scales_deg(sc), azi_filt{sc}, ele_filt{sc}, im_corr(sc), para.usegpu);
    end
end
    


    