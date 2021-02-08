%% Example script for polarisation calculation
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); % elf main directory
addpath('/Users/jamesf/Dropbox/Matlab/polarELF1.7.1 171108');
load(['/Users/jamesf/Documents/Lunar DoLP/'...
        '20161113-fullmoonclear1-Stonehenge/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'dolp.mat'])%    

mid         = [1+(size(dolpim{1}, 1)-1)/2; 1+(size(dolpim{1}, 2)-1)/2];   % centre of image
shortSide   = min([size(dolpim{1}, 1) size(dolpim{1}, 2)]);
r_full      = 8 * shortSide / 24;                           % theoretical value for 24mm high chip that is fully covered by fisheye circular image
[y, x]      = meshgrid(1:size(dolpim{1}, 2), 1:size(dolpim{1}, 1));
r           = sqrt((x-mid(1)).^2 + (y-mid(2)).^2);
ele         = asind(r / 2 / r_full) * 2;%r/r_full * 90;

sel         = ele>(60);
im = dolpim{3};
i1      = im(:, :, 1); i1(sel) = NaN;
surf(i1, 'EdgeColor','none'); colormap(jet); colorbar
histogram(i1)
axis([0, 1, 0, 5*10^4]);

dolpim1 = dolpim{3}; dolpim1(dolpim1 >1) = 1;
elecrop = 30%crop above this elevation
ewind = (900+elecrop*10):(2700-elecrop*10);%East West indices
nsind = (0+elecrop*10):(1800-elecrop*10);%North South indices
% imshow(dolpim1(nsind, ewind), 'Colormap', jet); 
hst = figure()
histogram(dolpim1(nsind, ewind), 'EdgeColor','none');
axis([0, 1, 0, 5*10^4]);
xlabel('Degree of Linear Polarization'); 
ylabel('Value Frequency (above 30° elevation)');
hh = findobj(gca,'Type','patch');
set(hh,'Facecolor','r','facealpha',0.50)
hold on
load(['/Users/jamesf/Documents/Lunar DoLP/'...
        '20161116-gibbousmoon2-Stonehenge/set1/'...'20161113-fullmoonclear1-Stonehenge/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'dolp.mat'])%  
dolpim1 = dolpim{3}; dolpim1(dolpim1 >1) = 1;
elecrop = 30%crop above this elevation
ewind = (900+elecrop*10):(2700-elecrop*10);%East West indices
nsind = (0+elecrop*10):(1800-elecrop*10);%North South indices
histogram(dolpim1(nsind, ewind), 'EdgeColor','none');
hh = findobj(gca,'Type','patch');
set(hh,'Facecolor','g','facealpha',0.50)
load(['/Users/jamesf/Documents/Lunar DoLP/'...
        '20161108-quartermoon-Stonehenge/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20161113-fullmoonclear1-Stonehenge/'...'20170130-MilkyWay-Stonehenge/'...
    'dolp.mat'])%  
dolpim1 = dolpim{3}; dolpim1(dolpim1 >1) = 1;
elecrop = 30%crop above this elevation
ewind = (900+elecrop*10):(2700-elecrop*10);%East West indices
nsind = (0+elecrop*10):(1800-elecrop*10);%North South indices
histogram(dolpim1(nsind, ewind), 'EdgeColor','none');
hh = findobj(gca,'Type','patch');
set(hh,'Facecolor','k','facealpha',0.50)
load(['/Users/jamesf/Documents/Lunar DoLP/'...
        '20171123-crescentmoonclear-Thornwood/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20161113-fullmoonclear1-Stonehenge/'...
    'dolp.mat'])%  
dolpim1 = dolpim{3}; dolpim1(dolpim1 >1) = 1;
elecrop = 30%crop above this elevation
ewind = (900+elecrop*10):(2700-elecrop*10);%East West indices
nsind = (0+elecrop*10):(1800-elecrop*10);%North South indices
histogram(dolpim1(nsind, ewind), 'EdgeColor','none');
hh = findobj(gca,'Type','patch');
set(hh,'Facecolor','k','facealpha',0.50)
load(['/Users/jamesf/Documents/Lunar DoLP/'...
        '20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20161113-fullmoonclear1-Stonehenge/'...
    'dolp.mat'])%  
dolpim1 = dolpim{3}; dolpim1(dolpim1 >1) = 1;
elecrop = 30%crop above this elevation
ewind = (900+elecrop*10):(2700-elecrop*10);%East West indices
nsind = (0+elecrop*10):(1800-elecrop*10);%North South indices
histogram(dolpim1(nsind, ewind), 'EdgeColor','none');

hh = findobj(gca,'Type','patch');
set(hh,'Facecolor','k','facealpha',0.50)
histogram(dolpim1(nsind, ewind), 'EdgeColor','none');

% surf(dolpim1, 'EdgeColor','none'); colormap(jet); colorbar
% 
% 
% hhh = figure(1) 
% I = elf_io_correctdng(intim2{2},  [], 'bright');
% polar_plot_int3d(I, hhh, 1); % third argument is radius
% hold on
% polar_plot_aop3d(pos, aop3d_surf{2}, hhh, 0.1, 1.01);