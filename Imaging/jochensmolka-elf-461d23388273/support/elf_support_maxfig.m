function elf_support_maxfig(handle)

warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
if nargin<1
    handle=gcf;
end
pause(0.01); %makes this work even if the figure has only just been opened 
jf=get(handle,'JavaFrame');
set(jf,'Maximized',1);
drawnow;
%jf.getFigurePanelContainer.getComponent(0).getTopLevelAncestor.setMaximized(1)