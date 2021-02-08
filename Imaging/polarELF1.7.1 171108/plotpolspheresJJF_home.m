close all ;clear all
%% Plot each sky
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); % elf main directory
addpath('/Users/jamesfoster/Dropbox/Matlab/polarELF1.7.1 171108');
addpath('/Users/jamesfoster/Dropbox/Matlab/polarELF1.7.1 171108/modules/polar/polar_plot');
addpath('/Users/jamesfoster/Dropbox/Matlab/polarELF1.7.1 171108/modules/stitch')
load(['/Volumes/Big backup 2/20181122MaxGibbousWits/set1/'...
load(['/Users/jamesfoster/Documents/'...
        '20161113-fullmoonclear1-Stonehenge/set2/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'ints.mat'])% 
load(['/Users/jamesfoster/Documents/'...
        '20161113-fullmoonclear1-Stonehenge/set2/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'aop3d.mat'])% 
load(['/Users/jamesfoster/Documents/'...
        '20161113-fullmoonclear1-Stonehenge/set2/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'aopos.mat'])% 
load(['/Users/jamesfoster/Documents/'...
        '20161113-fullmoonclear1-Stonehenge/set2/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'dolp.mat'])% 
I = elf_io_correctdng(intim{3},  [], 'bright');
figure();set(gcf,'Position', [1920/3,1920/2,1.05*1920/3,1920/3])
ah(1) = polar_plot_int3d(I,  1);view(0,90); axis square; axis off
% saveas(gcf, '/Users/jamesfoster/Documents/FullMoonInt.pdf',  'pdf')
saveas(gcf, '/Users/jamesfoster/Documents/FullMoonWitsInt.pdf',  'pdf')
hold on
% figure();set(gcf,'Position', [1920/3,1920/2,1.05*1920/3,1920/3])
polar_plot_aop3d(pos, aop3d_surf{3},gcf, 0.1, 1.01);view(0,90); axis square; axis off
view(340, 25)
saveas(gcf, '/Users/jamesfoster/Documents/FullMoonAOP.pdf',  'pdf')

figure();set(gcf,'Position', [1920/3,1920/2,1.05*1920/3,1920/3])
polar_plot_dolp3d(dolpim{3}, gca, 0, 1);view(0,90); axis square; axis off; hold on
saveas(gcf, '/Users/jamesfoster/Documents/FullMoonWitsDoLP.pdf',  'pdf')

polar_plot_aop3d(pos, aop3d_surf{3},gcf, 0.1, 1.01);view(0,90); axis square; axis off
view(340, 25)
% load(['/Users/jamesfoster/Documents/'...
%         '20171128-quartermoonovercast-Thornwood/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
%     'ints.mat'])% 
load(['/Users/jamesfoster/Documents/'...
        '20171128-quartermoonovercast-Thornwood/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'aop3d.mat'])% 
load(['/Users/jamesfoster/Documents/'...
        '20161113-fullmoonclear1-Stonehenge/set2/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'aopos.mat'])% 
% I = elf_io_correctdng(intim{3},  [], 'bright');
% figure();set(gcf,'Position', [1920/3,1920/2,1.05*1920/3,1920/3])
% ah(1) = polar_plot_int3d(I,  1);view(0,90); axis square; axis off
% saveas(gcf, '/Users/jamesfoster/Documents/FullMoonInt.pdf',  'pdf')
% hold on



saveas(gcf, '/Users/jamesfoster/Documents/quartermoonovercastAOP.pdf',  'pdf')
load(['/Users/jamesf/Documents/Lunar DoLP/'...
        '20161113-fullmoonclear1-Stonehenge/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'dolp.mat'])%    
figure();set(gcf,'Position', [1920/3,1920/2,1.05*1920/3,1920/3])
polar_plot_dolp3d(dolpim{3}, gca, 0, 1);view(0,90); axis square; axis off
saveas(gcf, '/Users/jamesf/Documents/Lunar DoLP/FullMoonSphere.pdf',  'pdf')
load(['/Users/jamesf/Documents/Lunar DoLP/'...
        '20161116-gibbousmoon2-Stonehenge/set1/'...'20161113-fullmoonclear1-Stonehenge/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'dolp.mat'])% 
% figure()
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
saveas(gcf, '/Volumes/BIG BACKUP/Lunar DoLP/waninggibbous-IllmitzSphere.pdf',  'pdf')

load(['/Volumes/BIG BACKUP/Lunar DoLP/'...
        '20170304-clearquartermoon-tower-Sodankyl‰/set1/'...'20161116-gibbousmoon2-Stonehenge/set1/'...'20170130-MilkyWay-Stonehenge/'...'20161118-quartermoonhorizon2-Stonehenge/set2/'...
    'dolp.mat'])%doesn't exist yet  
cla
polar_plot_dolp3d(dolpim{3}, gca, 0, 1);view(0,90); axis off
colormapeditor
saveas(gcf, '/Volumes/BIG BACKUP/Lunar DoLP/clearquartermoon-Sodankyl‰Sphere.pdf',  'pdf')