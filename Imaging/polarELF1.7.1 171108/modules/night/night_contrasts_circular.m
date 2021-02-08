function night_contrasts_circular(filename, rotate)
% Loads prefiltered Milky-Way images in circular projection and calculates different measures of contrasts between the Northern and Southern arm
%
% Filename must include a cell array im_filt_HDR with four elements. By
% default, these are assumed to be filtered at 2/4/8/16 degrees, but this
% only matters for labels. The second input variable 'rotate' should be
% true if the Milky Way is vertical in the image, as the algorithm assumes
% a roughly horizontal orientation.

if nargin < 2, rotate = false; end
if nargin < 1, filename = 'F:\James data\20151116\filt\scene001_filt_reproj.mat', end %#ok<NOPRT> % for testing
addpath('E:\sprograms\useful');

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
ele_s  = ele; ele_s(y<mid(2)) = -ele_s(y<mid(2)); % define everything on the left half as negative

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
            im = rot90(permute(im, [2 1 3]), 2);        % rotate image if necessary to place Milky Way horizontally
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
    
    %% 2. ABSCONT: For each point from left to right, get the maximum along the vertical dimension (only in the centre of the image to exclude horizon light)
    for sc = 1:4
        [absmax{sc, ch}, absmax_pos{sc, ch}] = max(sumim{sc, ch}(401:600, :), [], 1);
        absmax_x{sc, ch}    = 1:length(absmax{sc, ch});
        absmax_y{sc, ch}    = absmax_pos{sc, ch}+400;
        linind              = sub2ind(size(ele_s), absmax_x{sc, ch}, absmax_y{sc, ch});
        absmax_ele{sc, ch}  = ele_s(linind); 

        % calculate Michelson contrast between the maximum in the South and the maximum in the North
        Nmax_abs            = max(absmax{sc, ch}(absmax_ele{sc, ch}>0));
        Smax_abs            = max(absmax{sc, ch}(absmax_ele{sc, ch}<0));
        abscont(sc, ch)     = abs((Nmax_abs - Smax_abs) / (Nmax_abs + Smax_abs));
    end

    %% 3. LINECONT: Fit a line to the maxima in the most highly filtered image, and evaluate contrasts along that line
  
    % find the first and last value > 0.8 * median, and only use points in between (to cut out black borders)
    sel = absmax{end, ch}>0.6*median(absmax{end, ch});
    sel(find(sel, 1, 'first'):find(sel, 1, 'last')) = 1;
    
    allsels{ch} = sel;
    
    milkyfit = polyfit(absmax_x{end, ch}(sel), absmax_y{end, ch}(sel), 2);                                 % fit a quadratic function
   
    for sc = 1:4
        fitx{sc, ch} = absmax_x{sc, ch};%%FIXME: Don't have to save these two variables 4 times each
        fity{sc, ch} = round(polyval(milkyfit, fitx{sc, ch})); % in pixels
    
        % Find the intensity at these points
        linind              = sub2ind(size(sumim{sc, ch}), fity{sc, ch}, fitx{sc, ch}); % calculate a linear index (rather than looping through both vectors)
        linemax{sc, ch}     = sumim{sc, ch}(linind);
     
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
        
        bgsum(mask) = NaN;
        c1 = plotim{sc, ch}(:, :, 1);
        c2 = plotim{sc, ch}(:, :, 2);
        c3 = plotim{sc, ch}(:, :, 3);
        c1(mask) = NaN;
        c2(mask) = NaN;
        c3(mask) = NaN;
        bg_im{sc, ch} = cat(3, c1, c2, c3);
        
        background = nanmean(bg_sum(:));
        MWcontS{sc, ch} = (S2{sc, ch} - background) ./ (S2{sc, ch} + background);
        MWcontN{sc, ch} = (N2{sc, ch} - background) ./ (N2{sc, ch} + background);
                
    end
end

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
for sc = 1:4
    for ch = 1:4
        pnum    = (sc-1)*4+ch;                                      % plot number
        linim   = plotim{sc, ch}/max(plotim{4, ch}(:));             % normalised image for plotting (normalise all to the max of the highest filtered)
        linim(linim>1) = 1;
%         im      = elf_io_correctdng(linim);                        % correct using default color matrix
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
        ylim(ah4{pnum}, [0 yy]);
        plot(ah4{pnum}, [-elecont_ele(sc, ch) -elecont_ele(sc, ch) NaN elecont_ele(sc, ch) elecont_ele(sc, ch)], [0 yy NaN 0 yy], 'k'); %FIXME: needs to find x position of those contrasts
        axis(ah4{pnum}, [-70 70 0 4e10])
        set(ah4{pnum}, 'xtick', -60:30:60, 'xticklabels', num2str([30 60 90 60 30]'));
        grid on;        

        
        axes(ah5{pnum}); %#ok<LAXES>
        plot(90-targetele, S2{sc, ch}/10^9, 'b', 90-targetele, N2{sc, ch}/10^9, 'r', 'linewidth', 2);
        hold on;
        plot(90-targetele, alleleconts{sc, ch}*100, 'k', 90-targetele, MWcontS{sc, ch}*100, 'c', 90-targetele, MWcontN{sc, ch}*100, 'm');
        a = axis(ah5{pnum});
        axis(ah5{pnum}, [20 90 -20 max([60 a(4)])]);
        set(ah5{pnum}, 'xtick', [30 60 90]);
        grid on
        

    end
end

figure(6);        
linim   = plotim{1, 4}/max(plotim{4, 4}(:));             % normalised image for plotting (normalise all to the max of the highest filtered)
linim(linim>1) = 1;
linim2 = bg_im{1, 4}/max(plotim{4, 4}(:));             % normalised image for plotting (normalise all to the max of the highest filtered)
linim2(linim2>1) = 1;
subplot(1, 2, 1); imagesc(linim); axis image off; set(gca, 'YDir', 'normal');
subplot(1, 2, 2); imagesc(linim2); axis image off; set(gca, 'YDir', 'normal');

figure(7); clf; plot(ele_s(:, 500), sumim{1, 4}(:, 500), 'k', 'linewidth', 2);
axis([-70 70 0 4e10])
set(gca, 'xtick', -60:30:60, 'xticklabels', num2str([30 60 90 60 30]'));
grid on;        

%% Axes and labels
axis([ah1{:} ah2{:}], 'image', 'off');

set([ah3{:}], 'xlim', [0 1000]);%, 'xtick', -90:30:90);
ylabel(ah3{9}, 'spectral photon irradiance (photons m^{-2} s^{-1} sr^{-1} nm^{-1})', 'interpreter', 'tex');
for i = 13:16, xlabel(ah3{i}, 'x-position (pixels)'); end
set([ah3{1:12}], 'xticklabel', []);
title(ah3{1}, '2\circ, red')
title(ah3{2}, 'green')
title(ah3{3}, 'blue')
title(ah3{4}, 'white')
title(ah3{5}, '4\circ')
title(ah3{9}, '8\circ')
title(ah3{13}, '16\circ')

ylabel(ah4{9}, 'spectral photon irradiance (photons m^{-2} s^{-1} sr^{-1} nm^{-1})', 'interpreter', 'tex');
for i = 13:16, xlabel(ah4{i}, 'elevation (\circ)'); end
title(ah4{1}, '2\circ, red')
title(ah4{2}, 'green')
title(ah4{3}, 'blue')
title(ah4{4}, 'white')
title(ah4{5}, '4\circ')
title(ah4{9}, '8\circ')
title(ah4{13}, '16\circ')

ylabel(ah5{13}, 'spectral photon irradiance (*10^9 photons m^{-2} s^{-1} sr^{-1} nm^{-1}) / Michelson contrast (%)', 'interpreter', 'tex', 'horizontalalignment', 'left');
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
outname = [f2, '_', f(1:end-12)];

uicontrol('style', 'text', 'string', outname, 'backgroundcolor', 'w', 'horizontalalignment', 'left', 'fontweight', 'bold', 'fontsize', 12, 'units', 'normalized', 'position', [0.01 0.97 0.6 0.03], 'parent', 1)
uicontrol('style', 'text', 'string', outname, 'backgroundcolor', 'w', 'horizontalalignment', 'left', 'fontweight', 'bold', 'fontsize', 12, 'units', 'normalized', 'position', [0.01 0.97 0.6 0.03], 'parent', 2)
uicontrol('style', 'text', 'string', outname, 'backgroundcolor', 'w', 'horizontalalignment', 'left', 'fontweight', 'bold', 'fontsize', 12, 'units', 'normalized', 'position', [0.01 0.97 0.6 0.03], 'parent', 3)

drawnow;
pause(0.1);
pdfsave(1, [outname '_c_1.pdf']);
pdfsave(2, [outname '_c_2.pdf']);
pdfsave(3, [outname '_c_3.pdf']);
pdfsave(4, [outname '_c_4.pdf']);
pdfsave(5, [outname '_c_5.pdf']);
pdfsave(6, [outname '_c_6.pdf']);
pdfsave(7, [outname '_c_7.pdf']);


        