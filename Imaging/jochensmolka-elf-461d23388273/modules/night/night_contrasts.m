function night_contrasts(filename, rotate)
% Loads prefiltered Milky-Way images in equirectangular projection and calculates different measures of contrasts between the Northern and Southern arm
%
% Filename must include a cell array im_filt_HDR with four elements. By
% default, these are assumed to be filtered at 2/4/8/16 degrees, but this
% only matters for labels. The second input variable 'rotate' should be
% true if the Milky Way is vertical in the image, as the algorithm assumes
% a roughly horizontal orientation.

%% THIS IS THE OLD VERSION OF THIS FUNCTION WHICH USES THE EQUIRECTANGULAR IMAGES. USE MILKYWAY_CONTRASTS_CIRCULAR INSTEAD!

if nargin < 2, rotate = false; end
if nargin < 1, filename = 'F:\James data\20151116\filt\scene001_filt.mat', end %#ok<NOPRT> % for testing
addpath('E:\sprograms\useful');

%% load filtered images
temp    = load(filename, 'im_filt_HDR');
ims     = temp.im_filt_HDR;
clear temp;

%% init
absmax      = cell(4, 4);
absmax_pos  = cell(4, 4);
absmax_ang  = cell(4, 4);
abscont     = zeros(4, 4);
linemax     = cell(4, 4);
linemax_ang = cell(4, 4);
linecont    = zeros(4, 4);
elecont     = zeros(4, 4);

%% calculate contrasts
for ch = 1:4
    %% 1. extract images
    for sc = 1:4
        im = ims{sc};
        if rotate, im = permute(im, [2 1 3]); end       % rotate image if necessary to place Milky Way horizontally
        if ch < 4, im(:, :, [1:ch-1 ch+1:3]) = 0; end   % if colour channel, set other channels to 0 (works for plotting and contrast calculation)
        plotim{sc, ch} = im;                            % save back into ims structure, as this will later be plotted
        sumim{sc, ch}  = sum(im, 3);
    end
    
    %% 2. ABSCONT: For each point from left to right, get the maximum along the vertical dimension
    for sc = 1:4
        [absmax{sc, ch}, absmax_pos{sc, ch}] = max(sumim{sc, ch}, [], 1);
        absmax_ang{sc, ch}  = round(100*linspace(-90, 90, length(absmax{sc, ch})))/100; % FIXME
        % calculate Michelson contrast between the maximum in the South and the maximum in the North
        Nmax_abs            = max(absmax{sc, ch}(absmax_ang{sc, ch}>0));
        Smax_abs            = max(absmax{sc, ch}(absmax_ang{sc, ch}<0));
        abscont(sc, ch)     = abs((Nmax_abs - Smax_abs) / (Nmax_abs + Smax_abs));
    end

    %% 3. LINECONT: Fit a line to the maxima in the most highly filtered image, and evaluate contrasts along that line
    x = round(100*linspace(-90, 90, length(absmax{end, ch})))/100;  % angles along the x- or y-axis in degrees
    y = x(absmax_pos{end, ch});                                     % position of maxima along the y dimension in degrees
    
    % find the first and last value > 0.8 * median, and only use points
    % in between (to cut out black borders)
    sel = absmax{end, ch}>0.8*median(absmax{end, ch});
    sel(find(sel, 1, 'first'):find(sel, 1, 'last')) = 1;
       
    milkyfit = polyfit(x(sel), y(sel), 2);                                 % fit a quadratic function
    
    for sc = 1:4
        x_fit_pix = 1:length(absmax{sc, ch});
        x = round(100*linspace(-90, 90, length(x_fit_pix)))/100; % angles along the x- or y-axis in degrees (rounded to the nearest 0.01)

        y_fit_ang = polyval(milkyfit, x); % in degrees
        y_fit_pix = zeros(size(y_fit_ang));
        for i = 1:length(y_fit_ang)
            [~, y_fit_pix(i)] = min(abs(x - y_fit_ang(i)));
        end
        
        % Now x_fit_pix and y_fit_pix are the coordinates of the fitted Milky Way in the current image.
        % Save for plotting:
        fitx{sc, ch} = x_fit_pix;
        fity{sc, ch} = y_fit_pix;
        
        % Find the intensity at these points
        linind              = sub2ind(size(sumim{sc, ch}), y_fit_pix, x_fit_pix); % calculate a linear index (rather than looping through both vectors)
        linemax{sc, ch}     = sumim{sc, ch}(linind);
        linemax_ang{sc, ch} = absmax_ang{sc, ch}; %FIXME
        
        % calculate Michelson contrast between the maximum in the South and the maximum in the North
        Nmax_line = max(linemax{sc, ch}(linemax_ang{sc, ch}>0));
        Smax_line = max(linemax{sc, ch}(linemax_ang{sc, ch}<0));
        linecont(sc, ch) = abs((Nmax_line - Smax_line) / (Nmax_line + Smax_line));
    end
        
    %% 4. ELECONT: Calculate contrasts between points of equal elevation along the Milky Way, the pick the maximum of those
    for sc = 1:4
        % determine valid range
        sel             = absmax{sc, ch}>0.8*median(absmax{sc, ch});
        minvalidangle   = absmax_ang{sc, ch}(find(sel, 1, 'first'));
        maxvalidangle   = absmax_ang{sc, ch}(find(sel, 1, 'last'));
        limang          = min(abs([maxvalidangle minvalidangle])); 
        
        S = linemax{sc, ch}(linemax_ang{sc, ch}>0 & abs(linemax_ang{sc, ch})<=limang);
        N = rot90(linemax{sc, ch}(linemax_ang{sc, ch}<0 & abs(linemax_ang{sc, ch})<=limang),2);
        aS = linemax_ang{sc, ch}(linemax_ang{sc, ch}>0 & abs(linemax_ang{sc, ch})<=limang);
        aN = linemax_ang{sc, ch}(linemax_ang{sc, ch}<0 & abs(linemax_ang{sc, ch})<=limang);
        alleleconts = (S-N)./(S+N);
        [elecont(sc, ch), p] = max(abs(alleleconts));
        elecont_ang(sc, ch) = abs(aS(p));
    end
end

%% plot results
ah1 = formatA4(1, 16, 'p', 0.01);
ah2 = formatA4(2, 16, 'p', 0.01);
ah3 = formatA4(3, 16, 'l');
for sc = 1:4
    for ch = 1:4
        pnum    = (sc-1)*4+ch;                                      % plot number
        im      = plotim{sc, ch}/max(plotim{sc, ch}(:));            % normalised image for plotting

        % Figure 1: Image only
        imagesc(im, 'Parent', ah1{pnum});                           % plot image in figure 1
        
        % Figure 2: Image + Max points + Fit lines
        imagesc(im, 'Parent', ah2{pnum});                           % plot image in figure 2
        hold(ah2{pnum}, 'on');
        plot(ah2{pnum}, 1:size(im, 1), absmax_pos{sc, ch}, 'r.');   % Individual maxima
        plot(ah2{pnum}, fitx{sc, ch}, fity{sc, ch}, 'w');           % Fitted Milky Way

        % Figure 3: Intensity profiles along individual maxima and fitted Milky Way
        axes(ah3{pnum}); %#ok<LAXES>
        plot(ah3{pnum}, absmax_ang{sc, ch}, absmax{sc, ch}, 'b', 'linewidth', 3);
        hold(ah3{pnum}, 'on');
        plot(ah3{pnum}, linemax_ang{sc, ch}, linemax{sc, ch}, 'r', 'linewidth', 2);
        text(0, 0, sprintf('C_{abs} = %0.1f%%, C_{line} = %0.1f%%,\nC_{ele} = %0.1f%%', 100*abscont(sc, ch), 100*linecont(sc, ch), 100*elecont(sc, ch)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8, 'Parent', ah3{pnum});
        
        % Y-Axis limits and balck line indicating max contrast
        yy = ylim(ah3{pnum});
        ylim(ah3{pnum}, [0 yy(2)]);
        plot(ah3{pnum}, [-elecont_ang(sc, ch) -elecont_ang(sc, ch) NaN elecont_ang(sc, ch) elecont_ang(sc, ch)], [0 yy(2) NaN 0 yy(2)], 'k');
        grid on;
    end
end

%% Axes and labels
axis([ah1{:} ah2{:}], 'image', 'off');

set([ah3{:}], 'xlim', [-90 90], 'xtick', -90:30:90);
ylabel(ah3{9}, 'spectral photon irradiance (photons m^{-2} s^{-1} sr^{-1} nm^{-1})', 'interpreter', 'tex');
for i = 13:16, xlabel(ah3{i}, '(almost) elevation (\circ)'); end
set([ah3{1:12}], 'xticklabel', []);
title(ah3{1}, '2\circ, red')
title(ah3{2}, 'green')
title(ah3{3}, 'blue')
title(ah3{4}, 'white')
title(ah3{5}, '4\circ')
title(ah3{9}, '8\circ')
title(ah3{13}, '16\circ')

[p, f]  = fileparts(filename);
[~, f2] = fileparts(fileparts(p));
outname = [f2, '_', f(1:end-5)];

uicontrol('style', 'text', 'string', outname, 'backgroundcolor', 'w', 'horizontalalignment', 'left', 'fontweight', 'bold', 'fontsize', 12, 'units', 'normalized', 'position', [0.01 0.97 0.6 0.03], 'parent', 1)
uicontrol('style', 'text', 'string', outname, 'backgroundcolor', 'w', 'horizontalalignment', 'left', 'fontweight', 'bold', 'fontsize', 12, 'units', 'normalized', 'position', [0.01 0.97 0.6 0.03], 'parent', 2)
uicontrol('style', 'text', 'string', outname, 'backgroundcolor', 'w', 'horizontalalignment', 'left', 'fontweight', 'bold', 'fontsize', 12, 'units', 'normalized', 'position', [0.01 0.97 0.6 0.03], 'parent', 3)

drawnow;
pause(0.1);
pdfsave(1, [outname '_1.pdf']);
pdfsave(2, [outname '_2.pdf']);
pdfsave(3, [outname '_3.pdf']);


        