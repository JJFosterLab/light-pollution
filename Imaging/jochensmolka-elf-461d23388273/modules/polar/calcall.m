
%%%%%%%%%%%%%%%%
%% USE POLAR.M INSTEAD!
%%%%%%%%%%%%%%%%



%% James polarisation
addpath('E:\sprograms\elf');
addpath('E:\sprograms\useful');
elf_paths;

% datasets must be a list of data folder names in the order UP/N/W/S/E
% each dataset ust contain five scenes, with filter orientations (relative to the horizon) V/45CCW/H/135CCW/V
para = elf_para('E:\James data\MW\set1');
[~, ~, datasets] = elf_checkdata(para);
nohor = [65 65 65 65 65]; % This just blocks out the filter edges


% % 1. extract all scenes, and rotate them around their optical axis (<5min)
% for compDir = 1:length(datasets)
%     dataset = datasets{compDir};
%     polar_main1_HDRscenes(dataset, nohor(compDir), 90); % Turns all images by 90 degrees, so the scenes have sky up in all photos
% end
% % 2. Filter images with 2/4/8/16 degree filter half-width (>30min)
% for compDir = 1:length(datasets)
%     dataset = datasets{compDir};
%     polar_main2_filter(dataset); close all
% end
% 3. Calculate stokes vectors (<2min)
pol_angle_flip = [0 0 1 0 1]; % whether to flip horizontal and vertical polarisation in each dataset
% for compDir = 1:length(datasets) % for each compass direction
%     dataset = datasets{compDir};
%     polar_main3_stokes(dataset, pol_angle_flip(compDir));
    polar_main3_stokes(datasets, pol_angle_flip);
% end
% 4. Calculate Intensity, AoP, and DoLP (<2min)
off_axis_correct = false;
channel = 3;
[int, aop, dolp] = polar_main4_polar_parameters(datasets, channel, off_axis_correct);


%% Plot individual compass directions (one channel only)
if 0
for compDir = 1:length(datasets)
    dataset = datasets{compDir};
    para    = elf_para('', dataset);
    infosum = elf_readwrite(para, 'loadinfosum');     % loads the old infosum file (which contains projection information)
    imStokes_filt = elf_readwrite(para, 'loadstokes_mat', sprintf('exposure%d', 1));
    close all;
    ah = formatA4(1, 15, 'l', 0);
    for filtLevel = 1:5
%         elf_plot_image(int{compDir, filtLevel}, infosum, ah{0+filtLevel}, 'undefined', 1);
% %         elf_plot_image(imStokes_filt{filtLevel}(:, :, :, 3), infosum, ah{5+filtLevel}, 'undefined', 0);
% %         elf_plot_image(imStokes_filt{filtLevel}(:, :, :, 4), infosum, ah{10+filtLevel}, 'undefined', 0);
% % 

        elf_plot_image(int{compDir, filtLevel}, infosum, ah{0+filtLevel}, 'undefined', 1);
        polar_plot_dolp(dolp{compDir, filtLevel}, ah{5+filtLevel})
        polar_plot_aop(aop{compDir, filtLevel}, ah{10+filtLevel}, dolp{compDir, filtLevel});
    end
    pdfsave(1, sprintf('%s_Exp%d_images.pdf', para.paths.dataset, 1));
end
end

%% stitch
 
filtlevel = 3
azi = -90:0.4:90; %FIXME
ele = 90:-0.4:-90;



infosum2 = infosum; infosum2.Width = 1000; infosum2.Height = 1000;
[finalw_im, finalh_im, iBestImage] = stitch_5dirs(infosum2, [0 90 180 270], -180:0.1:180, 90:-0.1:-90, 0, 'rect');
% finalx_im = round(finalx_im);
% finaly_im = round(finaly_im);

% TODO: preallocate
intim = zeros(1801, 3601, 3);
aopim = zeros(1801, 3601);
dolpim = zeros(1801, 3601);

% for each pixel of the output image
for i = 1:size(finalw_im, 1)
    for j = 1:size(finalw_im, 2)
        [~, x] = min(abs(azi - finalw_im(i, j)));
        [~, y] = min(abs(ele - finalh_im(i, j)));
        intim(i, j, :) = int{iBestImage(i, j), filtlevel}(y, x, :);
        aopim(i, j, :) = aop{iBestImage(i, j), filtlevel}(y, x, :);
        dolpim(i, j) = dolp{iBestImage(i, j), filtlevel}(y, x);

%         intim(i, j, :) = int_reproj{iBestImage(i, j), filtlevel}(finalx_im(i, j), finaly_im(i, j), :);
%         aopim(i, j, :) = aop_reproj{iBestImage(i, j), filtlevel}(finalx_im(i, j), finaly_im(i, j), :);
%         dolpim(i, j, :) = dolp_reproj{iBestImage(i, j), filtlevel}(finalx_im(i, j), finaly_im(i, j), :);
        
    end
    if mod(i, 100)==0
        i
    end
end

% % for i = 1:max(iBestImage(:))
% %     % create indices
% %     sel = iBestImage==i;
% %     [row, col] = find(sel);
% %     lin_ind = elf_project_sub2ind(size(intim), col, row); % linear index into intensity image %% maybe other way
% % %     i1 = int{i, filtlevel}(:, :, 1);
% % %     i2 = int{i, filtlevel}(:, :, 2);
% % %     i3 = int{i, filtlevel}(:, :, 3);
% %     [~, x] = min(abs(azi - finalx_im(iBestImage==i)));
% %     [~, y] = min(abs(ele - finaly_im(iBestImage==i)));
% %     intim(lin_ind) = int{iBestImage==i}(lin_ind);
% % %     aopim(iBestImage==i) = aop{i, filtlevel}(y, x, :);
% % %     dolpim(iBestImage==i) = dolp{i, filtlevel}(y, x, :);
% % end

%%
aopim2 = aopim(:, :, 1);
aopim2(iBestImage == 2) = mod(aopim2(iBestImage==2), 180);  % N = top
aopim2(iBestImage == 3) = mod(aopim2(iBestImage==3), 180);
aopim2(iBestImage == 4) = mod(aopim2(iBestImage==4), 180);
aopim2(iBestImage == 5) = mod(aopim2(iBestImage==5), 180);
figure(10); image(5*intim/max(intim(:))); title('Int');
polar_plot_dolp(dolpim, 13);title('no offaxiscorr');
polar_plot_aop(aopim2, 14, dolpim);



