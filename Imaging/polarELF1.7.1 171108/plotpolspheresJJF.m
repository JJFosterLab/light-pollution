close all ;clear all
%% Plot each sky
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); % elf main directory
addpath('/Users/jamesf/Dropbox/Matlab/polarELF1.7.1 171108/io')
addpath('/Users/jamesf/Dropbox/Matlab/polarELF1.7.1 171108');
addpath('/Users/jamesf/Dropbox/Matlab/polarELF1.7.1 171108/modules/polar/polar_plot');
addpath('/Users/jamesf/Dropbox/Matlab/polarELF1.7.1 171108/modules/stitch')
load(['/Users/jamesf/Documents/Lunar DoLP/20181123FullmoonHazy-Wits/set3/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'ints.mat'])% 
load(['/Users/jamesf/Documents/Lunar DoLP/20181130BeforeThirdQuarter-Thornwood/set2/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'ints.mat'])% 
load(['/Users/jamesf/Documents/Lunar DoLP/'...
        '20181130-NoonSun-Thornwood/set1/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'aop3d.mat'])% 
load(['/Users/jamesf/Documents/Lunar DoLP/'...
        '20181130-NoonSun-Thornwood/set1/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'aopos.mat'])%
load(['/Users/jamesf/Documents/Lunar DoLP/'...
        '20170130-eveningsun-Stonehenge/set1/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'ints.mat'])% 
load(['/Users/jamesf/Documents/Lunar DoLP/'...
        '20170130-eveningsun-Stonehenge/set1/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'aop3d.mat'])% 
load(['/Users/jamesf/Documents/Lunar DoLP/'...
        '20170130-eveningsun-Stonehenge/set1/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'aopos.mat'])% 
load(['/Users/jamesf/Documents/Lunar DoLP/'...
        '20170130-eveningsun-Stonehenge/set1/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'dolp.mat'])% 
%low view
I = elf_io_correctdngJJF(intim{1},  [], 'bright');%'bright');
figure();set(gcf,'Position', [1920/3,1920/2,1.05*1920/3,1920/3])
ah(1) = polar_plot_int3d(I,  1);view(150,-20); axis square; axis off
saveas(gcf, '/Users/jamesf/Documents/eveningsunIntLowang.pdf',  'pdf')%

ah(1) = polar_plot_int3d(I,  1);view(80,30); axis square; axis off
saveas(gcf, '/Users/jamesf/Documents/eveningsunIntCentred.pdf',  'pdf')%

I = elf_io_correctdngJJF(intim{1},  [], 'maxval', 10^-5);%'bright');%grey image
ah(1) = polar_plot_int3d(I,  1);view(80,30); axis square; axis off; hold on
polar_plot_aop3d(pos, aop3d_surf{3},gcf, 0.1, 1.01);view(80+90,30); axis square;
saveas(gcf, '/Users/jamesf/Documents/eveningsunAoPCentred.pdf',  'pdf')%


I = elf_io_correctdngJJF(intim{1},  [], 'maxval', 10^-5);%'bright');%grey image
ah(1) = polar_plot_int3d(I,  1);view(330,20); axis square; axis off
hold on
polar_plot_aop3d(pos, aop3d_surf{3},gcf, 0.1, 1.01);view(150+90,-20); axis square;
saveas(gcf, '/Users/jamesf/Documents/eveningsunAoPlowAng.pdf',  'pdf')%

I = elf_io_correctdngJJF(intim{1},  [], 'maxval', 10^-5);%'bright');%grey image
ah(1) = polar_plot_int3d(I,  1);view(330,20); axis square; axis off
hold on
polar_plot_aop3d(pos, aop3d_surf{3},gcf, 0.1, 1.01);view(150+90+180,-20); axis square;
saveas(gcf, '/Users/jamesf/Documents/eveningsunAoPlowAng180.pdf',  'pdf')%
hold off
polar_plot_aop3d(pos, aop3d_surf{3},gcf, 0.1, 1.01);view(150+90,-20); axis square;
saveas(gcf, '/Users/jamesf/Documents/eveningsunAoPlowAngAll.pdf',  'pdf')%

I = elf_io_correctdngJJF(intim{1},  [], 'bright');%'bright');
figure();set(gcf,'Position', [1920/3,1920/2,1.05*1920/3,1920/3])
ah(1) = polar_plot_int3d(I,  1);view(150+00,20); axis square; axis off
hold on
polar_plot_dolp3d(dolpim{3},gcf, 0, 1);view(150+000,20); axis square;
saveas(gcf, '/Users/jamesf/Documents/eveningsunDoP.pdf',  'pdf')%

polar_plot_dolp3d(dolpim{3},gcf, 0, 1);view(000,90); axis square;
saveas(gcf, '/Users/jamesf/Documents/eveningsunDoPabove.pdf',  'pdf')%

I = elf_io_correctdngJJF(intim{3},  [], 'bright');%'bright');
figure();set(gcf,'Position', [1920/3,1920/2,1.05*1920/3,1920/3])
ah(1) = polar_plot_int3d(I,  1);view(0,90); axis square; axis off
saveas(gcf, '/Users/jamesf/Documents/MiddaySunThornIntQ75.pdf',  'pdf')
% saveas(gcf, '/Users/jamesfoster/Documents/TvärminneCrescentInt.pdf',  'pdf')%
% hold on
figure();set(gcf,'Position', [1920/3,1920/2,1.05*1920/3,1920/3])
polar_plot_aop3d(pos, aop3d_surf{3},gcf, 0.1, 1.01);view(0,90); axis square; axis off
saveas(gcf, '/Users/jamesfoster/Documents/MiddaySunThornAOP.pdf',  'pdf')

% load(['/Users/jamesfoster/Documents/'...
%         '20171128-quartermoonovercast-Thornwood/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
%     'ints.mat'])% 
load(['/Users/jamesfoster/Documents/'...
        '20171128-quartermoonovercast-Thornwood/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'aop3d.mat'])% 
load(['/Users/jamesfoster/Documents/'...
        '20161113-fullmoonclear1-Stonehenge/set2/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'aopos.mat'])% 

load(['/Users/jamesfoster/Documents/Lunar DoLP/'...
        '20170130-eveningsun-Stonehenge/set1'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'aop3d.mat'])% 
load(['/Users/jamesfoster/Documents/Lunar DoLP/'...
        '20170130-eveningsun-Stonehenge/set1'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'aopos.mat'])% 
% I = elf_io_correctdng(intim{3},  [], 'bright');
% figure();set(gcf,'Position', [1920/3,1920/2,1.05*1920/3,1920/3])
% ah(1) = polar_plot_int3d(I,  1);view(0,90); axis square; axis off
% saveas(gcf, '/Users/jamesfoster/Documents/FullMoonInt.pdf',  'pdf')
% hold on
figure();set(gcf,'Position', [1920/3,1920/2,1.05*1920/3,1920/3])
polar_plot_aop3d(pos, aop3d_surf{3},gcf, 0.1, 1.01);view(0,90); axis square; axis off
saveas(gcf, '/Users/jamesf/Documents/Lunar DoLP/20170130-eveningsun-Stonehenge.pdf',  'pdf')


load(['/Users/jamesf/Documents/Lunar DoLP/'...
        '20181130BeforeThirdQuarter-Thornwood/set2/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'dolp.mat'])%    
figure();set(gcf,'Position', [1920/3,1920/2,1.05*1920/3,1920/3])
polar_plot_dolp3d(dolpim{3}, gca, 0, 1);view(0,90); axis square; axis off
saveas(gcf, '/Users/jamesf/Documents/Lunar DoLP/FullMoonSphere.pdf',  'pdf')
load(['/Users/jamesf/Documents/Lunar DoLP/'...
        '20161116-gibbousmoon2-Stonehenge/set1/'...'20161113-fullmoonclear1-Stonehenge/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'dolp.mat'])% 
% figure()
cla
I = elf_io_correctdng(intim{3},  [], 'bright');
polar_plot_int3d(I, gca, 1);view(0,90); axis off

cla
polar_plot_dolp3d(dolpim{3}, gca, 0, 1);view(0,90); axis off
saveas(gcf, '/Users/jamesf/Documents/Lunar DoLP/GibbousMoonSphere.pdf',  'pdf')
load(['/Users/jamesf/Documents/Lunar DoLP/'...
        '20161108-quartermoon-Stonehenge/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20161113-fullmoonclear1-Stonehenge/'...'20170130-MilkyWay-Stonehenge/'...
    'dolp.mat'])%  
cla
polar_plot_dolp3d(dolpim{3}, gca, 0, 1);view(0,90); axis off
saveas(gcf, '/Users/jamesf/Documents/Lunar DoLP/QuarterMoonSphere.pdf',  'pdf')
load(['/Users/jamesf/Documents/Lunar DoLP/'...
        '20171123-crescentmoonclear-Thornwood/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20161113-fullmoonclear1-Stonehenge/'...
    'dolp.mat'])%  
cla
polar_plot_dolp3d(dolpim{3}, gca, 0, 1);view(0,90); axis off
saveas(gcf, '/Users/jamesf/Documents/Lunar DoLP/CrescentMoonSphere.pdf',  'pdf')
load(['/Users/jamesf/Documents/Lunar DoLP/'...
        '20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20161113-fullmoonclear1-Stonehenge/'...
    'dolp.mat'])%  
cla
polar_plot_dolp3d(dolpim{3}, gca, 0, 1);view(0,90); axis off
saveas(gcf, '/Users/jamesf/Documents/Lunar DoLP/NewMoonSphere.pdf',  'pdf')
load(['/Users/jamesf/Documents/Lunar DoLP/'...
        '20171128-quartermoonovercast-Thornwood/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'dolp.mat'])%doesn't exist yet  
cla
polar_plot_dolp3d(dolpim{3}, gca, 0, 1);view(0,90); axis off
saveas(gcf, '/Users/jamesf/Documents/Lunar DoLP/OvercastQMoonSphere.pdf',  'pdf')

load(['/Volumes/BIG BACKUP/Lunar DoLP/'...
        '20170910-waninggibbous-Illmitz/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'dolp.mat'])%doesn't exist yet  
cla
polar_plot_dolp3d(dolpim{3}, gca, 0, 1);view(0,90); axis off
saveas(gcf, '/Users/jamesf/Documents/Lunar DoLP/OvercastQMoonSphere.pdf',  'pdf')