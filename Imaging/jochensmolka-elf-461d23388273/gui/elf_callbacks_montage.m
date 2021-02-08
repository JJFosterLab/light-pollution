function elf_callbacks_montage(src, ~)
% elf_callbacks_montage(src, ~)
%
% This is the callback for the ButtonDown function for all montage image objects
% It opens the selected in a new large window (useful to inspect image details)


%% Calculate which image was clicked (column x, row y)
% this is half a point off, but easier than being accurate
cp = get(get(src, 'Parent'), 'currentpoint');
x = floor(cp(1, 1)/100) + 1; 
y = floor(cp(1, 2)/100) + 1;

% Assuming 100 x 100 pixel thumbs, the montag has a x b images
ms = size(get(src, 'CData'));   % montage size in pixels
a = ms(1)/100;                  % number of columns

imnr = (y-1) * a + x;

%% Now load that image
res = get(src, 'UserData');
im = imread(fullfile(res.root, res.fnames_im{imnr}));

%% And display it in a new figure
hf = figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]); % open a screen-sized figure window
hp = uipanel('Parent', hf); % Create a maximum size panel
hi = elf_plot_image(im, res.infosum, hp, 'equirectangular'); % display the image

set(hi, 'ButtonDownFcn', ''); % Disable the ButtonDown function, so clicking on this image does not open another one.

try, elf_support_maxfig; end % Try to maximise figure (not sure if this works on Mac)