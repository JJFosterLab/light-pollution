function elf_plot_hideui(fh, activate)

fignum = get(fh, 'Number');
if activate
    state = 'on';
else
    state = 'off';
end

set(findobj('tag', sprintf('fig%d_gui_BW', fignum)), 'visible', state);
set(findobj('tag', sprintf('fig%d_gui_R', fignum)), 'visible', state);
set(findobj('tag', sprintf('fig%d_gui_G', fignum)), 'visible', state);
set(findobj('tag', sprintf('fig%d_gui_B', fignum)), 'visible', state);
set(findobj('tag', sprintf('fig%d_gui_posslider', fignum)), 'visible', state);
set(findobj('tag', sprintf('fig%d_gui_rangeslider', fignum)), 'visible', state);

drawnow;