function elf_callbacks_exploretype(src, ~, gui, fnames, infosum, para)
% elf_callbacks_exploretype(src, ~, gui, fnames, infosum, para)
%
% This is the callback for the ButtonDown function for all explore image windows
% It changes the image displayed 

%% persistent data variables
persistent data;

%% calculate contrasts if contrast button was pressed
if strcmp(get(src, 'tag'), 'exploregui_contrasts')
    temp    = load(fnames.fil);
    im_filt = temp.varinput;
    [~, data.contrim] = space_analysis(para, im_filt, false, infosum); % true for verbose gives six contrast summary plots
    set(gui.b(5:22), 'enable', 'on');
    return; %FIXME: poor flow
end
    
%% if called by calibration checkbox, alter source to be buttongroup
if strcmp(get(src, 'tag'), 'exploregui_calib')
    src = gui.bg;
end

%% load correct image depending on which button was pressed
calib = get(findobj('Parent', get(src, 'Parent'), 'tag', 'exploregui_calib'), 'Value'); %whether to colour correct images
%calibpossible = strcmpi(get(findobj('Parent', get(src, 'Parent'), 'tag', 'exploregui_calib'), 'enable'), 'on'); %whether calibration is possible (i.e. whether this is a dng dataset)

switch get(get(src, 'SelectedObject'), 'Tag')
    case 'exploregui_button1' % original
        im = elf_imread(fnames.ori);
        type = 'equisolid';
    case 'exploregui_button2'
        im = elf_imread(fnames.pro);
        type = 'zone';
    case 'exploregui_button3'
        temp = load(fnames.fil);
        im = temp.varinput{1, 1};
        type = 'filt1';
        calib = 1;
    case 'exploregui_button4'
        temp = load(fnames.fil);
        im = temp.varinput{2, 1};
        type = 'filt10';
        calib = 1; %filtered images always have to be corrected, as they are in photons
    case 'exploregui_button5'
        im = data.contrim{1, 1, 2, 3};
        type = 'filt1';
        calib = 1; %filtered images always have to be corrected, as they are in photons
    case 'exploregui_button6'
        im = data.contrim{1, 2, 2, 3};
        type = 'filt10';
        calib = 1; %filtered images always have to be corrected, as they are in photons
    case 'exploregui_button7'
        im = data.contrim{2, 1, 2, 3};
        type = 'filt1';
        calib = 1; %filtered images always have to be corrected, as they are in photons
    case 'exploregui_button8'
        im = data.contrim{2, 2, 2, 3};
        type = 'filt10';
        calib = 1; %filtered images always have to be corrected, as they are in photons
    case 'exploregui_button9'
        im = data.contrim{3, 1, 2, 3};
        type = 'filt1';
        calib = 1; %filtered images always have to be corrected, as they are in photons
    case 'exploregui_button10'
        im = data.contrim{3, 2, 2, 3};
        type = 'filt10';
        calib = 1; %filtered images always have to be corrected, as they are in photons
        
    case 'exploregui_button11'
        im = data.contrim{1, 1, 1, 1};
        type = 'contrzone';
        calib = 0;
    case 'exploregui_button12'
        im = data.contrim{1, 1, 1, 2};
        type = 'contrzone';
        calib = 0;
    case 'exploregui_button13'
        im = data.contrim{1, 1, 2, 1};
        type = 'contrzone_h';
        calib = 0;
    case 'exploregui_button14'
        im = data.contrim{1, 1, 2, 2};
        type = 'contrzone_h';
        calib = 0;
        
    case 'exploregui_button15'
        im = data.contrim{2, 1, 1, 1};
        type = 'contrzone';
        calib = 0;
    case 'exploregui_button16'
        im = data.contrim{2, 1, 1, 2};
        type = 'contrzone';
        calib = 0;
    case 'exploregui_button17'
        im = data.contrim{2, 1, 2, 1};
        type = 'contrzone_h';
        calib = 0; 
    case 'exploregui_button18'
        im = data.contrim{2, 1, 2, 2};
        type = 'contrzone_h';
        calib = 0; 
        
    case 'exploregui_button19'
        im = data.contrim{3, 1, 1, 1};
        type = 'contrzone';
        calib = 0; 
    case 'exploregui_button20'
        im = data.contrim{3, 1, 1, 2};
        type = 'contrzone';
        calib = 0; 
    case 'exploregui_button21'
        im = data.contrim{3, 1, 2, 1};
        type = 'contrzone_h';
        calib = 0; 
    case 'exploregui_button22'
        im = data.contrim{3, 1, 2, 2};
        type = 'contrzone_h';
        calib = 0;
end
        
if calib
    gui.hi = elf_plot_image(elf_io_correctdng(im, infosum), infosum, gui.hp, type); % display the image
else
    gui.hi = elf_plot_image(im, infosum, gui.hp, type); % display the image
end
set(gui.hi, 'ButtonDownFcn', ''); % Disable the ButtonDown function, so clicking on this image does not open another one.

