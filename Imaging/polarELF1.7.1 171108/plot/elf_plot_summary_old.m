function elf_plot_summary(para, res, fh, name, calib)
% elf_plot_summary(d, fh, name)
%
% Inputs: descriptor mean variable d and figure/axes/panel handle h
%
% Uses elf_plot_hog, elf_plot_int

%%


if nargin < 4
    calib = 0;
end
if nargin < 3
    name='no name';
end
if nargin < 2
    fh = elf_support_formatA4(1);
end
fignum = get(fh, 'Number');

%% create panels for the different analysis elements
fp1 = uipanel('Parent', fh, 'Units', 'normalized', 'Position', [0 .4 1 .6], 'backgroundcolor', 'w', 'BorderWidth', 0);
fp2 = uipanel('Parent', fh, 'Units', 'normalized', 'Position', [0 0 1 .4], 'backgroundcolor', 'w', 'BorderWidth', 0);

%% Content panel 1
stdo1 = {'Parent', fp1, 'Units', 'normalized'};             % standard options for each gui element
stdo2 = {'backgroundcolor', [.8 .8 .8], 'fontweight', 'bold', 'callback', @elf_callbacks_elfgui};     % standard options for each gui element

% TODO: add 'tooltipstring', BackgroundColor', 'Callback', 'foregroundcolor'
uicontrol(stdo1{:}, 'backgroundcolor', 'w', 'Style', 'text',         'Position', [.1 .93 .5 .05],  'tag', 'gui_filename',  'String', name,   'fontsize', 13, 'HorizontalAlignment', 'Left');
uicontrol(stdo1{:}, stdo2{:}, 'Style', 'togglebutton', 'Position', [.68 .93 .1 .05],   'tag', sprintf('fig%d_gui_BW', fignum),        'String', 'White',  'fontsize', 10, 'foregroundcolor', 'w', 'fontweight', 'bold','Value', 1);
uicontrol(stdo1{:}, stdo2{:}, 'Style', 'togglebutton', 'Position', [.8 .93 .05 .05],   'tag', sprintf('fig%d_gui_R', fignum),         'String', 'R',      'fontsize', 10, 'foregroundcolor', 'r', 'fontweight', 'bold','Value', 0);
uicontrol(stdo1{:}, stdo2{:}, 'Style', 'togglebutton', 'Position', [.85 .93 .05 .05],  'tag', sprintf('fig%d_gui_G', fignum),         'String', 'G',      'fontsize', 10, 'foregroundcolor', 'g', 'fontweight', 'bold','Value', 0);
uicontrol(stdo1{:}, stdo2{:}, 'Style', 'togglebutton', 'Position', [.9 .93 .05 .05],   'tag', sprintf('fig%d_gui_B', fignum),         'String', 'B',      'fontsize', 10, 'foregroundcolor', 'b', 'fontweight', 'bold','Value', 0);

ax1 = axes(stdo1{:}, 'Position', [.1 .29 .85 .63], 'tag', 'gui_ax1', 'FontWeight', 'bold');
ax2 = axes(stdo1{:}, 'Position', [.1 .29 1.02*.85 .63], 'tag', 'gui_ax2');
ax3 = axes(stdo1{:}, 'Position', [.1 .16 .85 .05], 'tag', 'gui_ax3');
ax4 = axes(stdo1{:}, 'Position', [.1 .105 .85 .03], 'tag', 'gui_ax4');

uicontrol(stdo1{:}, 'Style', 'slider',  'Position', [.1 .06 .85 .025], 'tag', 'gui_posslider', 'Min', -4, 'Max', 4, 'SliderStep', [0.1/8 1/8], 'Value', 0, 'callback', @elf_callbacks_elfgui); % Slider sets the x-axis centre in log-units
uicontrol(stdo1{:}, 'Style', 'slider',  'Position', [.1 .03 .85 .025], 'tag', 'gui_rangeslider', 'Min', 1, 'Max', 9, 'SliderStep', [0.1/8 1/8], 'Value', 4, 'callback', @elf_callbacks_elfgui); % Slider sets the x-axis width in log-units

stdo3 = {'backgroundcolor', 'w', 'fontweight', 'bold', 'HorizontalAlignment', 'Center'};                    % standard options for each gui element
uicontrol(stdo1{:}, stdo3{:}, 'Style', 'text',  'Position', [.02 .16 .06 .05], 'String', 'Scene average');  % Slider sets the x-axis width in log-units
uicontrol(stdo1{:}, stdo3{:}, 'Style', 'text',  'Position', [.02 .053 .06 .03], 'String', 'Centre');        % Slider sets the x-axis width in log-units
uicontrol(stdo1{:}, stdo3{:}, 'Style', 'text',  'Position', [.02 .023 .06 .03], 'String', 'Range');         % Slider sets the x-axis width in log-units


%% Content panel 2
stdo3 = {'Parent', fp2, 'Units', 'normalized'};             % standard options for each gui element
ax5 = axes(stdo3{:}, 'Position', [0 0 1 1], 'tag', 'gui_ax5');

%% plot both subplots
elf_plot_int(para, res.int, res.totalint, ax1, ax2, ax3, ax4, calib, fignum);
elf_plot_int_setvis(fignum); % sets visibility of RGB plot using graphics object tags

%elf_plot_spatial(para, res.spatial, ax5);




