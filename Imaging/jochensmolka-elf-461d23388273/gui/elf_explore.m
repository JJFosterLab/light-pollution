function elf_explore(para, data, imname, name, calib, infosum)
% used by elf_callbacks_explore

%% prepare filenames

[~,f,e] = fileparts(imname);
[~,~,e2] = fileparts(para.paths.imgformat);

fnames.ori = fullfile(para.paths.datapath, [f e2]);
fnames.sce = fullfile(para.paths.datapath, para.paths.scenefolder, [f e]); 

%% Create GUI
%% figure
% open a screen-sized figure window
gui.hf = figure(65); clf;
set(gui.hf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], 'Name', 'Dataset Explorer'); 

%% 1. image panel
% get figure size in pixels
try elf_support_maxfig; end % Try to maximise figure
pause(0.1);
set(gui.hf, 'Units', 'Pixels');
hfs = get(gui.hf, 'Position'); 
% create a maximum size square panel, assume landscape orientation
gui.hp = uipanel('Parent', gui.hf, 'Units', 'Pixels', 'Position', [0 0 hfs(4) hfs(4)]); % 
set(gui.hp, 'units', 'normalized');
hps = get(gui.hp, 'position');

%% 2. button panel
rows = 16; % how many button rows will there be in the button group?

gui.bp      = uipanel('Parent', gui.hf, 'Units', 'normalized', 'Position', [hps(1)+hps(3) 0 0.05 1]); % Create a maximum size panel
    gui.bg      = uibuttongroup('Parent', gui.bp, 'Units', 'normalized', 'Position', [0 1/3 1 2/3], 'Visible', 'off');
        stdo        = {'Units', 'normalized', 'parent', gui.bg}; % standard options for gui elements
        gui.b(1)      = uicontrol(stdo{:}, 'Style', 'togglebutton', 'Position', [0 1-1/rows 1 1/rows], 'tag', 'exploregui_button1', 'String', 'Orig');
        gui.b(2)      = uicontrol(stdo{:}, 'Style', 'togglebutton', 'Position', [0 1-2/rows 1 1/rows], 'tag', 'exploregui_button2', 'String', 'Proj');
        gui.b(3)      = uicontrol(stdo{:}, 'Style', 'togglebutton', 'Position', [0 1-3/rows 1 1/rows], 'tag', 'exploregui_button3', 'String', 'Filt1');
        gui.b(4)      = uicontrol(stdo{:}, 'Style', 'togglebutton', 'Position', [0 1-4/rows 1 1/rows], 'tag', 'exploregui_button4', 'String', 'Filt10');
        gui.b(5)      = uicontrol(stdo{:}, 'Style', 'togglebutton', 'Position', [0 1-5/rows 1 1/rows], 'tag', 'exploregui_button5', 'String', 'LC1', 'enable', 'off');
        gui.b(6)      = uicontrol(stdo{:}, 'Style', 'togglebutton', 'Position', [0 1-6/rows 1 1/rows], 'tag', 'exploregui_button6', 'String', 'LC10', 'enable', 'off');
        gui.b(7)      = uicontrol(stdo{:}, 'Style', 'togglebutton', 'Position', [0 1-7/rows 1 1/rows], 'tag', 'exploregui_button7', 'String', 'RG1', 'enable', 'off');
        gui.b(8)      = uicontrol(stdo{:}, 'Style', 'togglebutton', 'Position', [0 1-8/rows 1 1/rows], 'tag', 'exploregui_button8', 'String', 'RG10', 'enable', 'off');
        gui.b(9)      = uicontrol(stdo{:}, 'Style', 'togglebutton', 'Position', [0 1-9/rows 1 1/rows], 'tag', 'exploregui_button9', 'String', 'YB1', 'enable', 'off');
        gui.b(10)     = uicontrol(stdo{:}, 'Style', 'togglebutton', 'Position', [0 1-10/rows 1 1/rows], 'tag', 'exploregui_button10', 'String', 'YB10', 'enable', 'off');

        gui.b(11)     = uicontrol(stdo{:}, 'Style', 'togglebutton', 'Position', [0   1-11/rows 1/2 1/rows], 'tag', 'exploregui_button11', 'String', 'Lvc', 'enable', 'off', 'tooltipstring', 'Luminance channel - vertical contrasts');
        gui.b(12)     = uicontrol(stdo{:}, 'Style', 'togglebutton', 'Position', [1/2 1-11/rows 1/2 1/rows], 'tag', 'exploregui_button12', 'String', 'Lve', 'enable', 'off', 'tooltipstring', 'Luminance channel - filtered vertical contrasts = horizontal edges');
        gui.b(13)     = uicontrol(stdo{:}, 'Style', 'togglebutton', 'Position', [0   1-12/rows 1/2 1/rows], 'tag', 'exploregui_button13', 'String', 'Lhc', 'enable', 'off', 'tooltipstring', 'Luminance channel - horizontal contrasts');
        gui.b(14)     = uicontrol(stdo{:}, 'Style', 'togglebutton', 'Position', [1/2 1-12/rows 1/2 1/rows], 'tag', 'exploregui_button14', 'String', 'Lhe', 'enable', 'off', 'tooltipstring', 'Luminance channel - filtered horizontal contrasts = vertical edges');
        gui.b(15)     = uicontrol(stdo{:}, 'Style', 'togglebutton', 'Position', [0   1-13/rows 1/2 1/rows], 'tag', 'exploregui_button15', 'String', 'RGvc', 'enable', 'off', 'tooltipstring', 'RG channel - vertical contrasts');
        gui.b(16)     = uicontrol(stdo{:}, 'Style', 'togglebutton', 'Position', [1/2 1-13/rows 1/2 1/rows], 'tag', 'exploregui_button16', 'String', 'RGve', 'enable', 'off', 'tooltipstring', 'RG channel - filtered vertical contrasts = horizontal edges');
        gui.b(17)     = uicontrol(stdo{:}, 'Style', 'togglebutton', 'Position', [0   1-14/rows 1/2 1/rows], 'tag', 'exploregui_button17', 'String', 'RGhc', 'enable', 'off', 'tooltipstring', 'RG channel - horizontal contrasts');
        gui.b(18)     = uicontrol(stdo{:}, 'Style', 'togglebutton', 'Position', [1/2 1-14/rows 1/2 1/rows], 'tag', 'exploregui_button18', 'String', 'RGhe', 'enable', 'off', 'tooltipstring', 'RG channel - filtered horizontal contrasts = vertical edges');
        gui.b(19)     = uicontrol(stdo{:}, 'Style', 'togglebutton', 'Position', [0   1-15/rows 1/2 1/rows], 'tag', 'exploregui_button19', 'String', 'YBvc', 'enable', 'off', 'tooltipstring', 'GB channel - vertical contrasts');
        gui.b(20)     = uicontrol(stdo{:}, 'Style', 'togglebutton', 'Position', [1/2 1-15/rows 1/2 1/rows], 'tag', 'exploregui_button20', 'String', 'YBve', 'enable', 'off', 'tooltipstring', 'GB channel - filtered vertical contrasts = horizontal edges');
        gui.b(21)     = uicontrol(stdo{:}, 'Style', 'togglebutton', 'Position', [0   1-16/rows 1/2 1/rows], 'tag', 'exploregui_button21', 'String', 'YBhc', 'enable', 'off', 'tooltipstring', 'GB channel - horizontal contrasts');
        gui.b(22)     = uicontrol(stdo{:}, 'Style', 'togglebutton', 'Position', [1/2 1-16/rows 1/2 1/rows], 'tag', 'exploregui_button22', 'String', 'YBhe', 'enable', 'off', 'tooltipstring', 'GB channel - filtered horizontal contrasts = vertical edges');

    set(gui.bg, 'SelectionChangeFcn', {@elf_callbacks_exploretype, gui, fnames, infosum}, 'SelectedObject', gui.b(2));  % Select Proj as the default
    set(gui.bg, 'Visible', 'on');

    % now create text and checkbox to select colour calibration
    uicontrol('style', 'text', 'Parent', gui.bp, 'Units', 'normalized', 'Position', [0 1/3-1/40 1 1/40], 'string', {'colour', 'correct'});
    gui.calib   = uicontrol('style', 'checkbox', 'Parent', gui.bp, 'Units', 'normalized', 'Position', [0.4 1/3-2/40 0.4 1/40], 'value', 1, 'tag', 'exploregui_calib', 'Callback', {@elf_callbacks_exploretype, gui, fnames, infosum});

    % and a button to read contrasts
    gui.contrasts = uicontrol('style', 'pushbutton', 'Parent', gui.bp, 'Units', 'normalized', 'String', 'CalcCont', 'Tooltipstring', 'Calculate vertical contrast diagrams', 'Position', [0 1/3-2/20 1 1/20], 'tag', 'exploregui_contrasts', 'Callback', {@elf_callbacks_exploretype, gui, fnames, infosum, para});

%% deactivate calib for non-dng datasets
if ~strcmpi(para.paths.imgformat(end-2:end), 'dng')
    set(gui.calib, 'value', 0, 'enable', 'off');
end
elf_callbacks_exploretype(gui.bg, [], gui, fnames, infosum);


%% 3. ELF panel
gui.dp = uipanel('Parent', gui.hf, 'Units', 'normalized', 'Position', [hps(1)+hps(3)+0.05 0 1-(hps(1)+hps(3)+0.05) 1]); % Create a maximum size panel
elf_plot_summary(para, data, gui.dp, name, calib);