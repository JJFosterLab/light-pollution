function elf_plot_spatial(para, d, ah1, ah2, ah3, axleg, meanele)

%% housekeeping and variables
axes(ah3);
cla;
hold on;
chnum = length(para.plot.typenames);

xos1 = para.plot.channeldist; % distance between channel plots, in pixels
xos2 = para.plot.scaledist;   % distance between scale plots, in pixels
y    = 1.2*d.binmean;
cols = para.plot.contchannelcolours;

%% plot each channel
if para.plot.plothists
    %% FIXME: This section has been changed but not tested
    % plot histograms instead of contrast ellipses (mostly for debugging)
%     for ch = 1:3
%         allhists = d.([para.plot.typenames{ch} 'hist']);
%         h = allhists(:, :, :, d.angles==90);
%         v = allhists(:, :, :, d.angles==0);
%         pan = get(ah, 'parent');
%         
%         for sc = 1:size(y, 2)   % for each scale
%             x0 = (ch-1)/3 + (sc-1)/6;
%             for bin = 1:size(y, 1)    % for each bin
%                 y0 = (bin-1)/5;
%                 a1 = axes('units', 'normalized', 'parent', pan, 'position', [x0 y0 1/6 1/10]);
%                 bar(0:0.01:1, h(:, bin, sc)); xlim([0 1]); axis off;
%                 a2 = axes('units', 'normalized', 'parent', pan, 'position', [x0 y0+1/10 1/6 1/10]);
%                 bar(0:0.01:1, v(:, bin, sc)); xlim([0 1]); axis off
%                 set([a1 a2], 'yscale', 'log');
%             end
%         end
%     end
    
else 
    
    %% plot contrast profiles
    yy = meanele';

    for ch = 1:chnum
        % scale 1
        cp = 100*para.plot.channelscaling(ch) * d.([para.plot.typenames{ch} 'prof']){1};
        line(cp, yy, 'parent', ah1, 'color', cols{ch}, 'linewidth', 2, 'tag', sprintf('plot_contprof_ch%d_sc%d', ch, 1));
        
        % scale 2
        cp = 100*para.plot.channelscaling(ch) * d.([para.plot.typenames{ch} 'prof']){2};
        line(cp, yy, 'parent', ah2, 'color', cols{ch}, 'linewidth', 2, 'tag', sprintf('plot_contprof_ch%d_sc%d', ch, 2));
    end
    xlabel(ah1, '1\circ Contrasts (%)', 'fontweight', 'bold');
    xlabel(ah2, '10\circ Contrasts (%)', 'fontweight', 'bold');
    set(ah1, 'XTick', 0:10:30, 'XTickLabel', {' '}, 'YTick', [], 'box', 'on');
    set(ah2, 'XTick', 0:10:30, 'XTickLabel', {' '}, 'YTick', [], 'box', 'on');
    axis(ah1, [0 30 -90 90]);
    axis(ah2, [0 30 -90 90]);
    
    text(0, -92.5, '0', 'HorizontalAlignment', 'Left', 'fontweight', 'bold', 'Parent', ah1);
    text(0, -92.5, '0', 'HorizontalAlignment', 'Left', 'fontweight', 'bold', 'Parent', ah2);
    text(30, -92.5, ['\color[rgb]{' num2str(cols{1}) '}' num2str(30/para.plot.channelscaling(1)) '/\color[rgb]{' num2str(cols{2}) '}' num2str(30/para.plot.channelscaling(2)) '\color[rgb]{0 0 0}/\color[rgb]{' num2str(cols{3}) '}' num2str(30/para.plot.channelscaling(3)) '\color[rgb]{0 0 0}/\color[rgb]{' num2str(cols{4}) '}' num2str(30/para.plot.channelscaling(4))], 'HorizontalAlignment', 'Right', 'fontweight', 'bold', 'Parent', ah1);
    text(30, -92.5, ['\color[rgb]{' num2str(cols{1}) '}' num2str(30/para.plot.channelscaling(1)) '/\color[rgb]{' num2str(cols{2}) '}' num2str(30/para.plot.channelscaling(2)) '\color[rgb]{0 0 0}/\color[rgb]{' num2str(cols{3}) '}' num2str(30/para.plot.channelscaling(3)) '\color[rgb]{0 0 0}/\color[rgb]{' num2str(cols{4}) '}' num2str(30/para.plot.channelscaling(4))], 'HorizontalAlignment', 'Right', 'fontweight', 'bold', 'Parent', ah2); 
    
    
    %% plot contrast bubbles
    axes(ah3);
    stdo = {'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle', 'Parent', ah3};
    for ch = 1:chnum
        % get all contrasts for this channel
        c = d.([para.plot.typenames{ch} 'CS']);
        
        % scaling
        c = c*30*para.plot.channelscaling(ch);
        
        for sc = 1:2
            x0 = 25;
            for bin = 1:size(c, 1)
                y0 = 22.5-sc*15;
                r = c(bin, sc);
                [cx, cy] = pol2cart(deg2rad((0:90)+90*(ch-1)), ones(1, 91)*r); % TURN depending on channel
                patch(x0 + [0 cx 0], y0 + y(bin) + [0 cy 0], cols{ch}, 'linestyle', 'none', 'parent', ah3);
                
                % 1/10 deg signs
                if sc == 1
                    text(x0-7, y0 + y(bin), '1\circ:', stdo{:});
                else
                    text(x0-7, y0 + y(bin), '10\circ:', 'FontWeight', 'bold', stdo{:});
                end
            end
        end
    end
    
    
    %% plot contrast ellipses
    lw = [1 2];
    for ch = 1:chnum
        % get all contrasts for this channel
        c = d.(para.plot.typenames{ch});
        
        % scaling
        c = c*30*para.plot.channelscaling(ch);

        % plot
        for sc = 1:size(c, 3)   % for each scale
            for bin = 1:size(c, 1)    % for each bin
                y0 = 37.5 - ch*15;
                [cx, cy] = pol2cart(deg2rad(d.angles+90), c(bin, :, sc)); % TURN contrast by 90 degrees to get horizontal structures showing as horizontal
                plot(ah3, [cx cx(1)], y0 + y(bin) + [cy cy(1)], 'color', cols{ch}, 'linewidth', lw(sc));
                plot(ah3, 0, y0 + y(bin), 'k.');
            end
        end
    end
    
    % Horizon bin lines
    plot([-15 -10 40 nan -15 -10 40], [9.5 30 30 nan -9.5 -30 -30], 'k--', 'linewidth', 1);
        
    % 1/10 deg legend
    text(-10, 1.5, '1\circ', stdo{:});
    text(-10, -1.5, '10\circ', 'FontWeight', 'bold', stdo{:});
    text(-10, 61.5, '1\circ', stdo{:});
    text(-10, 58.5, '10\circ', 'FontWeight', 'bold', stdo{:});
    text(-10, -58.5, '1\circ', stdo{:});
    text(-10, -61.5, '10\circ', 'FontWeight', 'bold', stdo{:});
    
    %% Reference circle
    x1 = 0;
    y1 = -97;
    h  = 3; w = 3; %% CHECK THIS!
    plot(ah3, [x1 x1], [y1-h/2 y1+h/2], 'k'); % vertical line
    plot(ah3, [x1-w/2 x1+w/2], [y1 y1], 'k'); % horizontal line
    rectangle('Position', [x1-w/2 y1-h/2 w h], 'Curvature', [1 1], 'parent', ah3);
    
    x2 = 25;
    r = 1.5; %% CHECK THIS!
    for ch = 1:length(cols)
        [cx, cy] = pol2cart(deg2rad((0:90)+90*(ch-1)), ones(1, 91)*r); % TURN depending on channel
        patch(x2 + [0 cx 0], y1 + [0 cy 0], cols{ch}, 'linestyle', 'none', 'parent', ah3);
    end
    
    stdo = {'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle', 'fontweight', 'bold', 'Parent', ah3};
    reflabel = ['\color[rgb]{' num2str(cols{1}) '}' num2str(5/para.plot.channelscaling(1)) '/\color[rgb]{' num2str(cols{2}) '}' num2str(5/para.plot.channelscaling(2)) '\color[rgb]{0 0 0}/\color[rgb]{' num2str(cols{3}) '}' num2str(5/para.plot.channelscaling(3)) '\color[rgb]{0 0 0}/\color[rgb]{' num2str(cols{4}) '}' num2str(5/para.plot.channelscaling(4))];
    text((x1+x2)/2, y1, reflabel, stdo{:});
    text(x1, y1-10, {'Directional', 'contrasts (%)'}, stdo{:});
    text(x2, y1-10, {'Centre-surround', 'contrasts (%)'}, stdo{:});

    
%     %% plot ffts %% FIXME
%     dx = 1./[180.1 181];
%     dy = 1./[80 20 79.9;80 20 80];
%     scale = [1 10];
%     for ch = 1:3
%         allffts = d.([para.plot.typenames{ch} 'fft']);
%         for sc = 1%:size(c, 3)   % for each scale
%             Nyq_f = 1/(2*0.1*scale(sc)); % Nyquist of data
%             for bin = 1:size(c, 1)    % for each bin
%                 xx = -Nyq_f : dx(sc) : Nyq_f-dx(sc); % wavenumber
%                 yy = -Nyq_f : dy(sc,bin) : Nyq_f-dy(sc,bin); % wavenumber
%                 
%                 selx = xx>=-1/scale(sc) & xx<1/scale(sc);
%                 sely = yy>=-1/scale(sc) & yy<1/scale(sc);
%                 
%                 fftsmooth = allffts{sc}{bin}(sely, selx);%imgaussfilt(allffts{sc}{bin}(sely, selx), 10); %10 for heavily smoothed
%                 fftsum = sum(fftsmooth(:));
%                 f1 = @(x) sum(fftsmooth(fftsmooth>x))/fftsum - 0.05;
%                 f2 = @(x) sum(fftsmooth(fftsmooth>x))/fftsum - 0.2;
%                 f3 = @(x) sum(fftsmooth(fftsmooth>x))/fftsum - 0.35;
% 
%                 fftcont = [fzero(f3, 1) fzero(f2, 0) fzero(f1, 0)];
%                 
%                 x0 = (0.8*ch+2.5)*xos1 + (sc-1)*xos1/2;
%                 contour(ah, 20*scale(sc)*xx(selx)+x0, 20*scale(sc)*yy(sely)+y(bin), fftsmooth, fftcont, 'linecolor', [0 0 0]);
%                 drawnow
%             end
%         end
% 
%     end
   
    
    %% labels and references (right hand side)
    stdo = {'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Middle', 'Parent', ah3};
    for bin = 1:length(y)    % for each bin
        text(5.4*xos1 + size(y, 2)*xos2, y(bin)+5, para.plot.binnames{bin}, stdo{:});
        text(5.4*xos1 + size(y, 2)*xos2, y(bin)-5, para.plot.binranges{bin}, stdo{:});
    end
    
    %1/10deg
    text(5.4*xos1 + size(y, 2)*xos2, y(end)+50, '1\circ', stdo{:}, 'color', 'k');
    text(5.4*xos1 + size(y, 2)*xos2, y(end)+40, '10\circ', stdo{:}, 'color', 'k', 'fontweight', 'bold');

    axis(ah3, [-15 40 -115 90]);
    % axis aspect ratio and position is set after return to summary_ls

    %% Contrast channel legend
    hold(axleg, 'on');
    rectangle('Position', [-2.5 -2.2 6.5 4.2], 'Curvature', [0.1 0.1], 'linewidth', 3, 'Parent', axleg);
    plot(axleg, [-2 -0.5], [1.4 1.4], 'color', cols{1}, 'linewidth', 3);
    plot(axleg, [-2 -0.5], [0.4 0.4], 'color', cols{2}, 'linewidth', 3);
    plot(axleg, [-2 -0.5], [-0.6 -0.6], 'color', cols{3}, 'linewidth', 3);
    plot(axleg, [-2 -0.5], [-1.6 -1.6], 'color', cols{4}, 'linewidth', 3);
    text(0, 1.5, 'radiance', 'Parent', axleg, 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Middle', 'FontWeight', 'bold', 'Color', cols{1});
    text(0, 0.5, 'red-green', 'Parent', axleg, 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Middle', 'FontWeight', 'bold', 'Color', cols{2});
    text(0, -0.5, 'green-blue', 'Parent', axleg, 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Middle', 'FontWeight', 'bold', 'Color', cols{3});
    text(0, -1.5, 'blue-red', 'Parent', axleg, 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Middle', 'FontWeight', 'bold', 'Color', cols{4});

    axis(axleg, [-2.6 4 -2.3 3]);
    axis(axleg, 'equal');
    
end


