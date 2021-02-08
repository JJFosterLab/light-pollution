function starry_radiance_circularJJF(filename, rotate, outputfolder)
%% starry_radiance_circularJJF.m
%---------------------------------------------------------------
%       USAGE: starry_radiance_circularJJF(filename, rotate, outputfolder) %
%
%      AUTHOR: Jochen Smolka                DATE: 2017 07 27
%    MODIFIED: James Foster                 DATE: 2017 08 07
%    
% DESCRIPTION: Descriptive statistics for a sky image in VEPS format. 
%              Adapted from milkyway_contrasts_circular.m
%
%      INPUTS:  1. Filename
%               2. Rotation angle for image
%               3. Folder for output images (PDF)
%     OUTPUTS:  .pdf figures   
%
%  REFERENCES:  Foster et al., 2017 http://dx.doi.org/10.1098/rstb.2016.0079
%---------------------------------------------------------------
%
% TO DO:
% - MIN & MAX 
% - CENTRAL TENDENCY
% - LOCATION OF EACH
% - NOTH--SOUTH RADIANCE PROFILE
% - EAST--WEST RADIANCE PROFILE
% - MILKY WAY : BACKGROUND CONTRAST
% - NORMALISE IMAGE TO VALUE
%
%---    JOCHEN SAYS...  ---%
% Loads prefiltered starry sky images in circular projection and calculates different measures of contrasts between the Northern and Southern arm
%
% Filename must include a cell array im_filt_HDR with four elements. By
% default, these are assumed to be filtered at 2/4/8/16 degrees, but this
% only matters for labels. The second input variable 'rotate' should be
% true if the Milky Way is vertical in the image, as the algorithm assumes
% a roughly horizontal orientation.
%---    --- --- --- 	---%

addpath('/Users/jamesf/Dropbox/Matlab/veps/useful');

%% DEBUGGERY
% filename = fullfile('/Users/jamesf/Documents/ExampleNEFs/JJFvepsTest/JJFmodified/6_aurora_diffuse_scene001_reproj.mat');% 5_aurora_north_scene001_reproj.mat');% 3_lund_dacke2013_scene001_reproj.mat');% 4_osnabrck_dacke2013_scene001_reproj.mat');% 3_lightpollution_odarslv_scene001_reproj.mat');% 4_norrevngsvg_scene001_reproj.mat');% Fosteretal2017MilkyWay_scene001_reproj.mat');% 2_stonehenge_lowerelev_scene001_reproj.mat');% 3_nicama_lowerelev_scene001_reproj.mat');% 1_stonehenge_highelev_scene001_reproj.mat');% 1_joburg_dacke2013_scene001_reproj.mat');% 4_lightpollution_joburg_scene001_reproj.mat');%    fullfile(reprojfolder, allfiles(i).name);
filename = fullfile('/Users/jamesf/Documents/ExampleNEFs/ToccoByrneExperiments/20170130-1700_scene001_reproj.mat');

nmax = 3*10^17;%VALUE TO MAXIMISE TO

rotate = false; %SHOULD ALREADY BE IN CORRECT ORIENTATION
outputfolder = dir('/Users/jamesf/Documents/ExampleNEFs/1MilkyNEFs/JJFmodified/');%IS THIS HOW dir WORKS?
%% load filtered images
temp    = load(filename, 'im_filt_reproj');
ims     = temp.im_filt_reproj;
clear temp;
%% Calculate elevation for each pixel (this should work for all images now, since they should all be square and the same size)
% For non-square images, BUGFIX THIS!
mid    = [1+(size(ims{1}, 1)-1)/2; 1+(size(ims{1}, 2)-1)/2];        % centre of image
r_full = 8 * size(ims{1}, 2) / 24;                                 % theoretical value for 24mm high chip that is fully covered by fisheye circular image
[x, y] = meshgrid(1:size(ims{1}, 1), 1:size(ims{1}, 2));
r      = sqrt((x-mid(1)).^2 + (y-mid(2)).^2);
ele    = asind(r / 2 / r_full) * 2;%r/r_full * 90;
ele_s  = ele; ele_s(y<mid(2)) = -ele_s(y<mid(2)); % define everything on the left half as negative %I SHOULD MAKE THIS EAST
%% init
absmax      = cell(4, 4);
absmax_pos  = cell(4, 4);
absmax_ele  = cell(4, 4);
abscont     = zeros(4, 4);
linemax     = cell(4, 4);
linemax_ele = cell(4, 4);
linecont    = zeros(4, 4);
elecont     = zeros(4, 4);
%% calculate contrasts

for ch = 1:4
    %% 1. extract images
    for sc = 1:4
        im = ims{sc};
        im(isnan(im)) = 0;
        if rotate
            im = rot90(permute(im, [2 1 3]), 2);        % rotate image if necessary to place Milky Way horizontally %SHOULD BE POSSIBLE TO KEEP MW VERTICAL
        else
            im = fliplr(im);
        end
        if ch < 4
            im(:, :, [1:ch-1 ch+1:3]) = 0; 
            sumim{sc, ch}  = sum(im, 3);
        else
            sumim{sc, ch}  = mean(im, 3);
        end   % if colour channel, set other channels to 0 (works for plotting and contrast calculation)
        plotim{sc, ch} = im;                            % save back into ims structure, as this will later be plotted
    end
        %% 2. ABSCONT: For each point from left to right, get the maximum along the vertical dimension (only in the centre of the image to exclude horizon light) %MAY NOT BE NECESSARY
    for sc = 1:4
        [absmax{sc, ch}, absmax_pos{sc, ch}] = max(sumim{sc, ch}(401:600, :), [], 1); %SHOULD FREE UP THIS REGION, THE MILKY WAY (OR OTHER STREAK) CAN GO IN ANY DIRECTION
        absmax_x{sc, ch}    = 1:length(absmax{sc, ch});
        absmax_y{sc, ch}    = absmax_pos{sc, ch}+400; %WHY +400? 30DEG?
        linind              = sub2ind(size(ele_s), absmax_x{sc, ch}, absmax_y{sc, ch});
        absmax_ele{sc, ch}  = ele_s(linind); 

        % calculate Michelson contrast between the maximum in the South and the maximum in the North
        Nmax_abs            = max(absmax{sc, ch}(absmax_ele{sc, ch}>0));
        Smax_abs            = max(absmax{sc, ch}(absmax_ele{sc, ch}<0));
        abscont(sc, ch)     = abs((Nmax_abs - Smax_abs) / (Nmax_abs + Smax_abs)); %THESE ARE NOT THE CONTRASTS IN THE FINAL OUTPUT
    end
%% 3. LINECONT: Fit a line to the maxima in the most highly filtered image, and evaluate contrasts along that line
  
    % find the first and last value > 0.8 * median, and only use points in between (to cut out black borders)   
    sel = absmax{end, ch}>0.6*median(absmax{end, ch});  %I MIGHT WANT TO DO THIS IN A DIFFERENT WAY
    sel(find(sel, 1, 'first'):find(sel, 1, 'last')) = 1;
    
    allsels{ch} = sel;
    
    milkyfit = polyfit(absmax_x{end, ch}(sel), absmax_y{end, ch}(sel), 2);                                 % fit a quadratic function
   %CONSIDER milkyfit = polyfit(absmax_y{end, ch}(sel), absmax_x{end, ch}(sel), 2);
   
end
%%  DODGE THE LOOP
ch = 4;
%OR MAYBE USER IDENTIFIED MILKY WAY
   wh_im = plotim{1, 4};
   lowlim =  10.3;10.5;  10.4; 11;11.5; 
   starlight_1 = log10( abs((wh_im(:,:,1)+ wh_im(:,:,2)+wh_im(:,:,3))) );
   starlight_1(starlight_1 < lowlim) = lowlim;
   surf(starlight_1, 'EdgeColor','none'); colormap(jet); colorbar()
%    surf(log10(abs((wh_im(:,:,1)+ wh_im(:,:,2)+wh_im(:,:,3)))),'EdgeColor','none')
   figure(); 
   imshow(flipud(wh_im)./(nmax));
   imshow(flipud(wh_im)./(3*10^10)); [brightMW, dimMW] = ginput(2);  %USER SELECTS BRIGHT AND DIM ENDS OF MILKY WAY       
   
    for sc = 1:4
        fitx{sc, ch} = absmax_x{sc, ch};%%FIXME: Don't have to save these two variables 4 times each
        fity{sc, ch} = round(polyval(milkyfit, fitx{sc, ch})); % in pixels
    
        % Find the intensity at these points
        linind              = sub2ind(size(sumim{sc, ch}), fity{sc, ch}, fitx{sc, ch}); % calculate a linear index (rather than looping through both vectors)
        linemax{sc, ch}     = sumim{sc, ch}(linind);                        %THIS IS WHAT WE CALL THE MILKY WAY
     
        r                   = sqrt((fitx{sc, ch}-mid(1)).^2 + (polyval(milkyfit, fitx{sc, ch})-mid(2)).^2);
        thisele             = asind(r / 2 / r_full) * 2;
        linemax_ele{sc, ch} = thisele; linemax_ele{sc, ch}(fitx{sc, ch}<mid(2)) = -linemax_ele{sc, ch}(fitx{sc, ch}<mid(2)); % define everything on the left half as negative

        % calculate Michelson contrast between the maximum in the South and the maximum in the North
        Nmax_line           = max(linemax{sc, ch}(linemax_ele{sc, ch}>0));
        Smax_line           = max(linemax{sc, ch}(linemax_ele{sc, ch}<0));
        linecont(sc, ch)    = abs((Nmax_line - Smax_line) / (Nmax_line + Smax_line));
    end
        %% 4. ELECONT: Calculate contrasts between points of equal elevation along the Milky Way, the pick the maximum of those
    for sc = 1:4
%         % determine valid range
%         sel             = absmax{sc, ch}>0.8*median(absmax{sc, ch});
%         minvalidangle   = absmax_ele{sc, ch}(find(sel, 1, 'first'));
%         maxvalidangle   = absmax_ele{sc, ch}(find(sel, 1, 'last'));
%         limang          = min(abs([maxvalidangle minvalidangle])); 
        limang = 70;
        elelimang = 60;
        
        N   = linemax{sc, ch}(linemax_ele{sc, ch}<0 & abs(linemax_ele{sc, ch})<=limang);
        S   = linemax{sc, ch}(linemax_ele{sc, ch}>0 & abs(linemax_ele{sc, ch})<=limang);
        aN  = linemax_ele{sc, ch}(linemax_ele{sc, ch}<0 & abs(linemax_ele{sc, ch})<=limang);
        aS  = linemax_ele{sc, ch}(linemax_ele{sc, ch}>0 & abs(linemax_ele{sc, ch})<=limang);

        targetele = 0:.1:90;
        S2{sc, ch} = interp1(aS, S, targetele);
        N2{sc, ch} = interp1(aN, N, -targetele);
        
        alleleconts{sc, ch} = (S2{sc, ch}-N2{sc, ch})./(S2{sc, ch}+N2{sc, ch});
        [elecont(sc, ch), elecont_x(sc, ch)] = max(abs(alleleconts{sc, ch}));
        elecont_ele(sc, ch) = targetele(elecont_x(sc, ch));
        
        % Milky Way contrast against background
        bg_sum = sumim{sc, ch};
        mask = false(size(bg_sum));
        mask(ele>elelimang) = true;
        for i = 1:length(fitx{sc, ch})
            mask(fity{sc, ch}(i)-150:fity{sc, ch}(i)+150, fitx{sc, ch}(i)) = true;
        end
        mask = mask(1:size(bg_sum, 1), 1:size(bg_sum, 2)); % cut anything larger than intensed mask
        
        bg_sum(mask) = NaN;
        c1 = plotim{sc, ch}(:, :, 1);
        c2 = plotim{sc, ch}(:, :, 2);
        c3 = plotim{sc, ch}(:, :, 3);
        c1(mask) = NaN;
        c2(mask) = NaN;
        c3(mask) = NaN;
        bg_im{sc, ch} = cat(3, c1, c2, c3);
        
        background{sc, ch} = nanmean(bg_sum(:));
        MWcontS{sc, ch} = (S2{sc, ch} - background{sc, ch}) ./ (S2{sc, ch} + background{sc, ch});%THESE ARE THE NUMBERS THAT DAVID WANTS
        MWcontN{sc, ch} = (N2{sc, ch} - background{sc, ch}) ./ (N2{sc, ch} + background{sc, ch});
                
    end
end %for ch = 1:4

% figure;
% im = sumim{1, 1};%/max(plotim{1, 1}(:));
% plot(ele(500, :), im(500, :, 1));
% hold on;
% plot(ele(:, 500), im(:, 500, 1))


%% plot results
ah1 = formatA4(1, 16, 'p', 0.01);
ah2 = formatA4(2, 16, 'p', 0.01);
ah3 = formatA4(3, 16, 'l');
ah4 = formatA4(4, 16, 'l');
ah5 = formatA4(5, 16, 'l');
% ADD NUMBERS TO PLOTS?
% ADD COLOURBAR?
for sc = 1:4
    for ch = 1:4
        pnum    = (sc-1)*4+ch;                                      % plot number
         %THIS IS WHERE I SHOULD SET THE PLOT VALUES
        linim   = plotim{sc, ch}/max(plotim{4, ch}(:));             % normalised image for plotting (normalise all to the max of the highest filtered)
        linim(linim>1) = 1;
%         im      = veps_io_correctdng(linim);                        % correct using default color matrix
        im = linim; % Milky Way looks clearer in images without gamma correction
        
        % Figure 1: Image only
        imagesc(im, 'Parent', ah1{pnum});                           % plot image in figure 1
        
        % Figure 2: Image + Max points + Fit lines
        imagesc(im, 'Parent', ah2{pnum});                           % plot image in figure 2
        hold(ah2{pnum}, 'on');
        plot(ah2{pnum}, absmax_x{sc, ch}(allsels{ch}), absmax_y{sc, ch}(allsels{ch}), 'r.');   % Individual maxima
        plot(ah2{pnum}, fitx{sc, ch}(allsels{ch}), fity{sc, ch}(allsels{ch}), 'w');           % Fitted Milky Way

        % Figure 3: Intensity profiles along individual maxima and fitted Milky Way
        axes(ah3{pnum}); %#ok<LAXES>
        plot(ah3{pnum}, absmax_x{sc, ch}, absmax{sc, ch}, 'b', 'linewidth', 3);
        hold(ah3{pnum}, 'on');
        plot(ah3{pnum}, absmax_x{sc, ch}, linemax{sc, ch}, 'r', 'linewidth', 2);
        
        % Y-Axis limits and black line indicating max contrast
        yy = ylim(ah3{pnum});
        ylim(ah3{pnum}, [0 yy(2)]);
        plot(ah3{pnum}, [mid(1)-elecont_x(sc, ch) mid(1)-elecont_x(sc, ch) NaN mid(1)+elecont_x(sc, ch) mid(1)+elecont_x(sc, ch)], [0 yy(2) NaN 0 yy(2)], 'k'); %FIXME: needs to find x position of those contrasts
        grid on;
        
        % Figure 4: Intensity profiles along individual maxima and fitted Milky Way
        axes(ah4{pnum}); %#ok<LAXES>
        hold(ah4{pnum}, 'on');
        % add a NaN in the centre to disconnect line
        x = [linemax_ele{sc, ch} 0];
        y = [linemax{sc, ch} NaN];
        [x, iii] = sort(x);
        y = y(iii);
        plot(ah4{pnum}, x, y, 'r', 'linewidth', 2);
        text(0, 0, sprintf('C_{abs} = %0.1f%%, C_{line} = %0.1f%%,\nC_{ele} = %0.1f%%', 100*abscont(sc, ch), 100*linecont(sc, ch), 100*elecont(sc, ch)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8, 'Parent', ah4{pnum});
        
        % Y-Axis limits and black line indicating max contrast
        yy = 4e10;
        plot(ah4{pnum}, [-elecont_ele(sc, ch) -elecont_ele(sc, ch) NaN elecont_ele(sc, ch) elecont_ele(sc, ch)], [0 yy NaN 0 yy], 'k'); %FIXME: needs to find x position of those contrasts
        a = axis(ah4{pnum});
        axis(ah4{pnum}, [-70 70 0 max([a(4) 4e10])])
        set(ah4{pnum}, 'xtick', -60:30:60, 'xticklabels', num2str([30 60 90 60 30]'));
        grid on;        

        
        axes(ah5{pnum}); %#ok<LAXES>
        plot(90-targetele, S2{sc, ch}/10^9, 'b', 90-targetele, N2{sc, ch}/10^9, 'r', 'linewidth', 2);
        hold on;
        plot(90-targetele, alleleconts{sc, ch}*100, 'k', 90-targetele, MWcontS{sc, ch}*100, 'c', 90-targetele, MWcontN{sc, ch}*100, 'm');
        a = axis(ah5{pnum});
        axis(ah5{pnum}, [20 90 -20 60]);%min([a(3) -20]) max([60 a(4)])]);
        set(ah5{pnum}, 'xtick', [30 60 90]);
        grid on
        

    end
end

figure(6);
%THIS IS WHERE I SHOULD SET THE PLOT VALUES
linim   = plotim{1, 4}/max(plotim{4, 4}(:));             % normalised image for plotting (normalise all to the max of the highest filtered) 
linim(linim>1) = 1;
linim2 = bg_im{1, 4}/max(plotim{4, 4}(:));             % normalised image for plotting (normalise all to the max of the highest filtered)
linim2(linim2>1) = 1;
subplot(1, 2, 1); imagesc(linim); axis image off; set(gca, 'YDir', 'normal');
subplot(1, 2, 2); imagesc(linim2); axis image off; set(gca, 'YDir', 'normal'); title(sprintf('Mean background radiance (white, max filt): %.3f', background{4, 4}));

figure(7); clf; plot(ele_s(:, 500), sumim{1, 4}(:, 500), 'k', 'linewidth', 2);
a = axis(gca);
axis([-70 70 0 max([a(4) 4e10])])
set(gca, 'xtick', -60:30:60, 'xticklabels', num2str([30 60 90 60 30]'));
grid on;
xlabel('elevation (\circ)');
ylabel('spectral photon radiance (photons m^{-2} s^{-1} sr^{-1} nm^{-1})');
title('Vertical transect across Milky Way');

%% Axes and labels
axis([ah1{:} ah2{:}], 'image', 'off');
axis(ah1{13}, 'on');
% add axes to the bottom left plot
x = ele2x(-90:30:90);
set(ah1{13}, 'xtick', x, 'xticklabels', num2str((-90:30:90)'), 'ytick', x, 'yticklabels', num2str((-90:30:90)'), 'tickdir', 'out');
        
set([ah3{:}], 'xlim', [0 1000]);%, 'xtick', -90:30:90);
ylabel(ah3{9}, 'spectral photon radiance (photons m^{-2} s^{-1} sr^{-1} nm^{-1})', 'interpreter', 'tex');
for i = 13:16, xlabel(ah3{i}, 'x-position (pixels)'); end
set([ah3{1:12}], 'xticklabel', []);
title(ah3{1}, '2\circ, red')
title(ah3{2}, 'green')
title(ah3{3}, 'blue')
title(ah3{4}, 'white')
title(ah3{5}, '4\circ')
title(ah3{9}, '8\circ')
title(ah3{13}, '16\circ')

ylabel(ah4{9}, 'spectral photon radiance (photons m^{-2} s^{-1} sr^{-1} nm^{-1})', 'interpreter', 'tex');
for i = 13:16, xlabel(ah4{i}, 'elevation (\circ)'); end
title(ah4{1}, '2\circ, red')
title(ah4{2}, 'green')
title(ah4{3}, 'blue')
title(ah4{4}, 'white')
title(ah4{5}, '4\circ')
title(ah4{9}, '8\circ')
title(ah4{13}, '16\circ')

ylabel(ah5{13}, 'spectral photon radiance (*10^9 photons m^{-2} s^{-1} sr^{-1} nm^{-1}) / Michelson contrast (%)', 'interpreter', 'tex', 'horizontalalignment', 'left');
for i = 13:16, xlabel(ah5{i}, 'elevation (\circ)'); end
title(ah5{1}, '2\circ, red')
title(ah5{2}, 'green')
title(ah5{3}, 'blue')
title(ah5{4}, 'white')
title(ah5{5}, '4\circ')
title(ah5{9}, '8\circ')
title(ah5{13}, '16\circ')
legend(ah5{3}, 'I_S', 'I_N', 'contrast', 'MW_S cont', 'MW_N cont', 'location', 'E')

[p, f]  = fileparts(filename);
[~, f2] = fileparts(fileparts(p));
outname = [f2, '_', f(1:end-7)];

uicontrol('style', 'text', 'string', outname, 'backgroundcolor', 'w', 'horizontalalignment', 'left', 'fontweight', 'bold', 'fontsize', 12, 'units', 'normalized', 'position', [0.01 0.97 0.6 0.03], 'parent', 1)
uicontrol('style', 'text', 'string', outname, 'backgroundcolor', 'w', 'horizontalalignment', 'left', 'fontweight', 'bold', 'fontsize', 12, 'units', 'normalized', 'position', [0.01 0.97 0.6 0.03], 'parent', 2)
uicontrol('style', 'text', 'string', outname, 'backgroundcolor', 'w', 'horizontalalignment', 'left', 'fontweight', 'bold', 'fontsize', 12, 'units', 'normalized', 'position', [0.01 0.97 0.6 0.03], 'parent', 3)

drawnow;
pause(0.1);
pdfsave(1, fullfile(outputfolder, [outname '_c_1.pdf']));
pdfsave(2, fullfile(outputfolder, [outname '_c_2.pdf']));
pdfsave(3, fullfile(outputfolder, [outname '_c_3.pdf']));
pdfsave(4, fullfile(outputfolder, [outname '_c_4.pdf']));
pdfsave(5, fullfile(outputfolder, [outname '_c_5.pdf']));
pdfsave(6, fullfile(outputfolder, [outname '_c_6.pdf']));
pdfsave(7, fullfile(outputfolder, [outname '_c_7.pdf']));

end % main

function e = x2ele(x)
    e = asind((x - 500.5) * 3 / 2000) * 2; 
end

function x = ele2x(e)
    x = sind(e / 2) * 2000 / 3 + 500.5;
end
        