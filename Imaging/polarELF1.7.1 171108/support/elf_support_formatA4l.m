function figh = elf_support_formatA4l(fignum, screenpos)
% figh = elf_support_formatA4l(fignum, screenpos)

if nargin < 2
    screenpos = 1;
end
if nargin < 1
    figh = figure;
else
    figh = figure(fignum);
end

% ss  = get(0, 'ScreenSize');  % [1 1 1920 1200]
% h   = (ss(4)-60);             % height in pixels          
% w   = (ss(4)-60) / (21/29.7); % width in pixels         
% pos = [1+(screenpos-1)*w  60 w h];
pos = [0.2 1.5 29.7 21];
orient(figh, 'landscape');
set(figh, 'Units', 'centimeters', 'position', pos);
set(figh, 'PaperUnits', 'centimeters', 'PaperSize', [29.7 21], ...
    'color', 'w', 'paperpositionmode', 'manual', 'paperposition', [1 .5 27.7 20], ...
    'Renderer', 'painters');%'zbuffer');
drawnow;
% 60 offset if to accommodate Taskbar
