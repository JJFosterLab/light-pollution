function elf_plot_summary_ls(para, res, meanim, infosum, fh, name, calib)
% elf_plot_summary(d, fh, name)
% 
% Inputs: descriptor mean variable d and figure/axes/panel handle h
% 
% Uses elf_plot_hog, elf_plot_int

if nargin < 7, calib = 0; end
if nargin < 6, name = 'no name'; end
if nargin < 5, fh = elf_support_formatA4l(1); end
fignum = get(fh, 'Number');

rim = [0.01 0 0.01 0.01]; %[L R B T]
gap = 0.002; % gap between main axes
gap2 = 0.002; % gap between contrast axes
regionw = 0.01; %width of the black bar marking elevation regions
axlab = 0.03;
bottompan = 0.11; % height of the bottom panel
bottomgap = 0.04; % gap above bottom panel 
contwidth = 0.17; % width of right-most (contrast) panel
axw = (1-2*gap-gap2-regionw-2*axlab-rim(1)-rim(2)-contwidth)/3;
axh = 1-bottompan-bottomgap-rim(3)-rim(4);

%% create axes 1 to 5
stdoax = {'Parent', fh, 'Units', 'normalized'};             % standard options for each axes element

ax1     = axes(stdoax{:}, 'Position', [rim(1)+axlab+1*gap/2+0*axw rim(3)+bottompan+bottomgap axw axh], 'tag', 'gui_ax1');  % axes for eman image
ax2     = axes(stdoax{:}, 'Position', [rim(1)+axlab+3*gap/2+1*axw rim(3)+bottompan+bottomgap axw axh], 'tag', 'gui_ax2'); % axes for intensity
ax2i    = axes(stdoax{:}, 'Position', [rim(1)+axlab+3*gap/2+1*axw rim(3)+bottompan+bottomgap axw+axw+gap+gap2+regionw axh], 'tag', 'gui_ax2i', 'visible', 'off'); % axes for total intensity
ax3     = axes(stdoax{:}, 'Position', [rim(1)+axlab+5*gap/2+2*axw rim(3)+bottompan+bottomgap axw/2 axh], 'tag', 'gui_ax3'); % axes for left spatial
ax4     = axes(stdoax{:}, 'Position', [rim(1)+axlab+5*gap/2+gap2+2.5*axw rim(3)+bottompan+bottomgap axw/2 axh], 'tag', 'gui_ax4'); % axes for right spatial
ax5     = axes(stdoax{:}, 'Position', [rim(1)+axlab+5*gap/2+regionw+gap2+3*axw rim(3)+bottompan+bottomgap-2.5/18*axh contwidth 20.5/18*axh], 'tag', 'gui_ax5'); % axes for far right spatial

%% Control buttons and sliders
stdo1 = {'Parent', fh, 'Units', 'normalized'};             % standard options for each gui element
stdo2 = {'backgroundcolor', [.8 .8 .8], 'fontweight', 'bold', 'callback', @elf_callbacks_elfgui};     % standard options for each gui element
x = rim(1)+axlab+3*gap/2+1*axw;
y = 1-rim(4)-0.03;
uicontrol(stdo1{:}, stdo2{:}, 'Style', 'togglebutton', 'Position', [x+0.01 y .04 .02],   'tag', sprintf('fig%d_gui_BW', fignum),   'String', 'White',  'fontsize', 8, 'foregroundcolor', 'w', 'fontweight', 'bold','Value', 1);
uicontrol(stdo1{:}, stdo2{:}, 'Style', 'togglebutton', 'Position', [x+0.05 y .02 .02],   'tag', sprintf('fig%d_gui_R', fignum),    'String', 'R',      'fontsize', 8, 'foregroundcolor', 'r', 'fontweight', 'bold','Value', 0);
uicontrol(stdo1{:}, stdo2{:}, 'Style', 'togglebutton', 'Position', [x+0.07 y .02 .02],  'tag', sprintf('fig%d_gui_G', fignum),    'String', 'G',      'fontsize', 8, 'foregroundcolor', 'g', 'fontweight', 'bold','Value', 0);
uicontrol(stdo1{:}, stdo2{:}, 'Style', 'togglebutton', 'Position', [x+0.09 y .02 .02],   'tag', sprintf('fig%d_gui_B', fignum),    'String', 'B',      'fontsize', 8, 'foregroundcolor', 'b', 'fontweight', 'bold','Value', 0);

uicontrol(stdo1{:}, 'Style', 'slider',  'Position', [x+0.01 rim(3)+bottompan+bottomgap+0.03 axw-0.02 .015], 'tag', sprintf('fig%d_gui_posslider', fignum), 'Min', -4, 'Max', 4, 'SliderStep', [0.1/8 1/8], 'Value', 0, 'callback', @elf_callbacks_elfgui); % Slider sets the x-axis centre in log-units
uicontrol(stdo1{:}, 'Style', 'slider',  'Position', [x+0.01 rim(3)+bottompan+bottomgap+0.01 axw-0.02 .015], 'tag', sprintf('fig%d_gui_rangeslider', fignum), 'Min', 1, 'Max', 9, 'SliderStep', [0.1/8 1/8], 'Value', 4, 'callback', @elf_callbacks_elfgui); % Slider sets the x-axis width in log-units

%% bottom bar
stdo3 = {'backgroundcolor', 'w', 'HorizontalAlignment', 'Center'};                    
ax6     = axes(stdoax{:}, 'Position', [rim(1)+gap/2, rim(3) axw bottompan], 'tag', 'gui_ax6', 'visible', 'off');
ax7     = axes(stdoax{:}, 'Position', [rim(1)+axlab+5*gap/2+2*axw rim(3) axw+gap2 bottompan], 'tag', 'gui_ax7', 'visible', 'off');

% Visual Environment Plot Standard & dataset name
% text(0, 0, {'\makebox[4in][c]{\fontsize{60}{40}\textbf{V}\fontsize{15}{40}\textbf{isual} \fontsize{60}{40}\textbf{E}\fontsize{15}{40}\textbf{nvironment}}',...
%             '', '\makebox[4in][c]{\fontsize{60}{40}\textbf{P}\fontsize{15}{40}\textbf{lot} \fontsize{60}{40}\textbf{S}\fontsize{15}{40}\textbf{tandard}}', ...
%             '', ['\makebox[4in][c]{' name '}']}, ...
%             'Interpreter', 'latex', 'Parent', ax6); 
text(0, .25, {'\makebox[4in][c]{\fontsize{60}{40}\textbf{V}\fontsize{15}{40}\textbf{isual} \fontsize{60}{40}\textbf{E}\fontsize{15}{40}\textbf{nvironment}}',...
            '', '\makebox[4in][c]{\fontsize{60}{40}\textbf{P}\fontsize{15}{40}\textbf{lot} \fontsize{60}{40}\textbf{S}\fontsize{15}{40}\textbf{tandard}}'}, ...
            'Interpreter', 'latex', 'HorizontalAlignment', 'Center', 'Parent', ax6); 
text(0, -.25, name, 'HorizontalAlignment', 'Center', 'Parent', ax6);
axis(ax6, [-1 1 -.5 .5])

axb     = axes(stdoax{:}, 'Position', [rim(1)+axlab+3*gap/2+1*axw rim(3)+bottompan/4 axw bottompan/2], 'tag', 'gui_axb');
axbi    = axes(stdoax{:}, 'Position', [rim(1)+axlab+3*gap/2+1*axw rim(3) axw bottompan/4], 'tag', 'gui_axbi');

%% plot both subplots
elf_plot_image(meanim, infosum, ax1, 'squashed', 0);
elf_plot_int(para, res.int, res.totalint, ax2, ax2i, [], [], axb, axbi, fignum);
elf_plot_int_setvis(fignum); % sets visibility of RGB plot using graphics object tags

elf_plot_spatial(para, res.spatial, ax3, ax4, ax5, ax7, res.totalint.region_meanele);

% Bring ah2i back to the foreground
uistack(ax2i, 'top');

% Move ax5 to the right position
set(ax5, 'units', 'pixels');
ax5pos = get(ax5, 'position');
ax5lim = axis(ax5);
ax5pos(3) = ax5pos(4) / (ax5lim(4)-ax5lim(3)) * (ax5lim(2)-ax5lim(1));
set(ax5, 'position', ax5pos, 'visible', 'off');
return

%% Old panel 1

% stdo3 = {'backgroundcolor', 'w', 'fontweight', 'bold', 'HorizontalAlignment', 'Center'};                    % standard options for each gui element
uicontrol(stdo1{:}, stdo3{:}, 'Style', 'text',  'Position', [.02 .16 .06 .05], 'String', 'Scene average');  % Slider sets the x-axis width in log-units
% uicontrol(stdo1{:}, stdo3{:}, 'Style', 'text',  'Position', [.02 .053 .06 .03], 'String', 'Centre');        % Slider sets the x-axis width in log-units
% uicontrol(stdo1{:}, stdo3{:}, 'Style', 'text',  'Position', [.02 .023 .06 .03], 'String', 'Range');         % Slider sets the x-axis width in log-units




