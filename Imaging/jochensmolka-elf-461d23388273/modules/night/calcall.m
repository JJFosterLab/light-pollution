%% James milky way
addpath('E:\sprograms\elf');
elf_paths;

% para = elf_para('F:\James data');
% [~, ~, datasets] = elf_checkdata(para);
% nohor = [NaN NaN 70 NaN NaN 70];
% 
% for i = [3 6]%1:length(datasets)
%     dataset = datasets{i};
%     night_HDRscenes(dataset, nohor(i));
%     night_filter(dataset); close all
% end

%%
allfiles = {'F:\James data\20151116\filt\scene001_filt.mat', ...
    'F:\James data\20151116\filt\scene002_filt.mat', ...
    'F:\James data\20151116\filt\scene003_filt.mat', ...
    'F:\James data\20151116 - nohor\filt\scene001_filt.mat', ...
    'F:\James data\20151116 - nohor\filt\scene002_filt.mat', ...
    'F:\James data\20151116 - nohor\filt\scene003_filt.mat', ...
    'F:\James data\20151116 - bracketed\filt\scene001_filt.mat', ...
    'F:\James data\20160204\filt\scene001_filt.mat', ...
    'F:\James data\20160204\filt\scene002_filt.mat', ...
    'F:\James data\20160204\filt\scene003_filt.mat', ...
    'F:\James data\20160204 - nohor\filt\scene001_filt.mat', ...
    'F:\James data\20160204 - nohor\filt\scene002_filt.mat', ...
    'F:\James data\20160204 - nohor\filt\scene003_filt.mat', ...
    'F:\James data\20160204 - bracketed\filt\scene001_filt.mat'};

% Calculate total irradiance
i = 9;
irrad = night_irradiance(allfiles{i}); % in photons / s / m2 / nm
disp('Irradiance (photons/s/cm2/nm)');
disp(irrad/10000)

% Means of Johnsen data:
% 400-500 (34:66)   :    2.6e6
% 500-600 (67:100)  :    4.9e6
% 600-700 (101:133) :    5.1e6
% 400-700: 4.21

% Means for our data (2deg filter):
% 4     5.3 3.7 2.2
% 5     5.4 3.8 2.2
% 6     5.4 3.8 2.2 !
% 11    5.0 3.7 2.2
% 12    5.1 3.7 2.3 !
% 13    7.1 4.3 2.9

% 1     5.5 3.8 2.3
% 2     5.6 3.9 2.3
% 3     5.6 3.9 2.3 ! mean: 3.94 - 3.96
% 8     5.8 4.2 2.4
% 9     5.9 4.2 2.5 ! mean: 4.19 - 4.20
% 10    8.3 4.9 3.3


return
% %% Reproject all images (only needs to be done once)
% for i = 1:length(allfiles)
%     filename = allfiles{i};
%     temp    = load(filename, 'im_filt_HDR');
%     ims     = temp.im_filt_HDR;
%     clear temp;
% 
%     for j = 1:length(ims)
%         im_filt_reproj{j} = elf_project_reproject2fisheye_simple(ims{j}, [], [], [1000 1000]); % im, azi, ele, imsize
%     end
%     [p, f, e] = fileparts(filename);
%     save(fullfile(p, [f '_reproj', e]), 'im_filt_reproj');
% end

%%
allfiles = {'F:\James data\20151116\filt\scene001_filt_reproj.mat', ...
    'F:\James data\20151116\filt\scene002_filt_reproj.mat', ...
    'F:\James data\20151116\filt\scene003_filt_reproj.mat', ...
    'F:\James data\20151116 - nohor\filt\scene001_filt_reproj.mat', ...
    'F:\James data\20151116 - nohor\filt\scene002_filt_reproj.mat', ...
    'F:\James data\20151116 - nohor\filt\scene003_filt_reproj.mat', ...
    'F:\James data\20151116 - bracketed\filt\scene001_filt_reproj.mat', ...
    'F:\James data\20160204\filt\scene001_filt_reproj.mat', ...
    'F:\James data\20160204\filt\scene002_filt_reproj.mat', ...
    'F:\James data\20160204\filt\scene003_filt_reproj.mat', ...
    'F:\James data\20160204 - nohor\filt\scene001_filt_reproj.mat', ...
    'F:\James data\20160204 - nohor\filt\scene002_filt_reproj.mat', ...
    'F:\James data\20160204 - nohor\filt\scene003_filt_reproj.mat', ...
    'F:\James data\20160204 - bracketed\filt\scene001_filt_reproj.mat'};

%% Analyse
close all;
for i = 6%1:7
    night_contrasts_circular(allfiles{i}, 1); 
    close all
end
for i = 12%8:14
    night_contrasts_circular(allfiles{i}, 0); 
    close all
end
