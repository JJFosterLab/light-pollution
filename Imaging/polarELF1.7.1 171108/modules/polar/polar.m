%% Example script for polarisation calculation
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); % elf main directory
addpath('/Users/jamesf/Dropbox/Matlab/polarELF1.7.1 171108');
elf_paths;

% datasets must be a list of data folder names in the order UP/N/W/S/E
% each dataset must contain five scenes, with filter orientations (relative to the horizon) V/45CCW/H/135CCW/V
para = elf_para('/Users/jamesf/Documents/ExampleNEFs/20161113-fullmoonclear1-Stonehenge/set2');%20171128-quartermoonovercast-Thornwood/');%'G:\20161113-fullmoonclear1-Stonehenge\set2');
[~, ~, datasets] = elf_checkdata(para);

%% Main calculation sequence (comment out the time-consuming first steps unless they need to be recomputed)
% 1. extract all scenes, and rotate them around their optical axis (2.5min)
horlimit = 65;
rotation =90; 0;%
for compDir = 1:length(datasets)
    dataset = datasets{compDir};
    polar_main1_HDRscenes(dataset, horlimit, rotation); % Turns all images by 90 degrees, so the scenes have sky up in all photos
end

% 2. Filter images with 2/4/8/16 degree filter half-width (>30min)
for compDir = 1:length(datasets)
    dataset = datasets{compDir};
    polar_main2_filter(dataset); close all
end

% 3. Calculate stokes vectors, Int/DoLP/AoP (22s)
channel = 3; % This channel is used to calculate DoLP and AoP
[int, dolp, aop] = polar_main3_stokes(datasets, channel);

% 4. Stitch (13s without scaling, 50s with scaling)
horlimit = 65; % This limit has to be re-applied here, after filtering. It can be reduced to be even more restrictive regarding pixels used for image scaling
scaling  = 0;  % This can be 0 - do not try to correct the small exposure differences between stitched images
               %             1 - scale images taking into account the overlap of each image with the UP-image
               %             2 - scale images taking into account all overlapping areas between images    
               % Corrections are only applied to DoLP and Int images, not to AoP
[intim, dolpim, aopim, iBestImage] = polar_main4_stitch(int, dolp, aop, horlimit, scaling);
scaling  = 2;
[intim2, dolpim2] = polar_main4_stitch(int, dolp, aop, horlimit, scaling);

% 5. 3D project AoP (<1s)
[a, e] = polar_gridsphere;
for filtlevel = 1:length(aop)
    [aop3d_surf{filtlevel}, pos] = polar_aop3d(aopim{filtlevel}, iBestImage, dolpim{filtlevel}, a, e);
end

%% plot output images (the combined image is too large to save as pdf; use these examples to create output)
textposx = [0 80 0 -80];
textposy = [80 0 -80 0];
textlabels = {'N', 'W', 'S', 'E'};
textargs = {'FontWeight', 'bold', 'color', 'w', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle'};
for filtlevel = 1:length(intim)
    %fh = figure(filtlevel); clf; maxfig; drawnow;
    fh = figure(filtlevel); clf; drawnow;
    I = elf_io_correctdng(intim{filtlevel},  [], 'bright');
    ah(1) = polar_plot_int(I, subplot(3, 2, 1)); title('Intensity');
    text(textposx, textposy, textlabels, textargs{:});
    ah(2) = polar_plot_dolp(dolpim{filtlevel}, subplot(3, 2, 3)); title('DoLP');
    text(textposx, textposy, textlabels, textargs{:});
    ah(3) = polar_plot_aop(aopim{filtlevel}, subplot(3, 2, 5), dolpim{filtlevel}, true); title('AoP');
    text(textposx, textposy, textlabels, textargs{:});
    I = elf_io_correctdng(intim2{filtlevel},  [], 'bright');
    ah(4) = polar_plot_int(I, subplot(3, 2, 2)); title('Intensity (smoothed transitions)');
    text(textposx, textposy, textlabels, textargs{:});
    ah(5) = polar_plot_dolp(dolpim2{filtlevel}, subplot(3, 2, 4)); title('DoLP (smoothed transitions)');
    text(textposx, textposy, textlabels, textargs{:});
    
    % Spherical plots
    ah(6) = subplot(3, 4, 11); cla; hold on;
    polar_plot_dolp3d(dolpim2{filtlevel}, ah(6), 0, 1); % last argument is radius. Use polar_plot_sphere3d(ah, col, r) to just plot a monochromatic sphere instead
    polar_plot_aop3d(pos, aop3d_surf{filtlevel}, ah(6), 0.1, 1.01); % last arguments are vector scale and radius
    view(20, 60); camva(4); axis off
    
    ah(7) = subplot(3, 4, 12); cla; hold on;
    polar_plot_int3d(I, ah(7), 1); % third argument is radius
    polar_plot_aop3d(pos, aop3d_surf{filtlevel}, ah(7), 0.1, 1.01);
    view(20, 60); camva(4); axis off
    
    fname = fullfile(para.paths.root, sprintf('results_%d.jpg', filtlevel));
    print(fh, fname, '-djpeg');
    close(fh);
end
