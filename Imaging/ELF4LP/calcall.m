%% James milky way
%% STEP 0: Gather file names
reprojfolder = '/Users/jamesf/Documents/Starry Project/Light Pollution Imaging/LPexper-20181127-LowMWay-Thornwood/radiance/';%'/Users/jamesf/Documents/ExampleNEFs/Tvärminne/best/clear/Tvärminne-longer/';%'E:\James data\04Jochen_reproj';

addpath('/Users/jamesf/Dropbox/Matlab/veps');%'E:\sprograms\veps');
addpath('/Users/jamesf/Dropbox/Matlab/ELF4LP');%'E:\sprograms\veps');
veps_paths;

para = veps_paraJJF('/Users/jamesf/Documents/Starry Project/Light Pollution Imaging/LPexper-20181127-LowMWay-Thornwood/radiance/');%'/Users/jamesf/Documents/ExampleNEFs/Tvärminne/best/clear/');%1MilkyNEFs');%'E:\James data\04Jochen');
[~, ~, datasets] = veps_checkdata(para);
nohor = repmat(70, 1, length(datasets));%[70 70 70 70];                     % The maximum eccentricity to show for each data folder. The view below this will be blackened before filtering. Use NaN to use the full image.


%% STEP 1: Calibrate and filter all images (only needs to be done once)
for i = 1:length(datasets)
    dataset = datasets{i};
    milkyway_HDRscenes(dataset, nohor(i)); % -> datafolder/scenes
    milkyway_filter(dataset); close all    % -> datafolder/filt
end

%% STEP 2: Reproject all images (only needs to be done once)
r_all = 90;% -90;% 90;% DOES THIS CHANGE ANYTHING?
rotate =   cellfun(@(x) r_all, cell(1, length(datasets)) , 'Uniformoutput', false);%THIS IS VERY UGLY%{-90 ; -90; -90; 0 };
% rotate{1} = 0; rotate{3} = 90;
%       {[-90 -90 -90 0 0 0 0 0 0 0 0], ... 
%           [0 0 0 0 0 0 0 0 0]};
%          {[0 0 0 0 0], ...
%           [45 0 0], ...
%           [0 0 0 0 0 0]}; % for each dataset, for each file, the angle that the image should be rotated by (clockwise)
      
warning('off', 'MATLAB:griddata:DuplicateDataPoints');
for i = 1:length(datasets)
    dataset = datasets{i};
    milkyway_reprojectJJF(dataset, rotate{i}, reprojfolder); % -> reprojfolder
end
warning('on', 'MATLAB:griddata:DuplicateDataPoints');
    
%% STEP 3: Analyse
allfiles    = dir(fullfile(reprojfolder, '*.mat'));
pdffolder   = reprojfolder;%'E:\James data\04Jochen_pdfoutput';

close all;
for i = 1:length(allfiles)
    filename = fullfile(reprojfolder, allfiles(i).name);
    milkyway_contrasts_circularJJF(filename, false, pdffolder); % -> pdffolder
    close all
end

%%
% 
% % Calculate total irradiance
% i = 1;
% irrad = milkyway_irradiance(allfiles{i}); % in photons / s / m2 / nm
% disp('Irradiance (photons/s/cm2/nm)');
% disp(irrad/10000)
% 
% % Means of Johnsen data:
% % 400-500 (34:66)   :    2.6e6
% % 500-600 (67:100)  :    4.9e6
% % 600-700 (101:133) :    5.1e6
% % 400-700: 4.21
% 
% % Means for our data (2deg filter):
% % 4     5.3 3.7 2.2
% % 5     5.4 3.8 2.2
% % 6     5.4 3.8 2.2 !
% % 11    5.0 3.7 2.2
% % 12    5.1 3.7 2.3 !
% % 13    7.1 4.3 2.9
% 
% % 1     5.5 3.8 2.3
% % 2     5.6 3.9 2.3
% % 3     5.6 3.9 2.3 ! mean: 3.94 - 3.96
% % 8     5.8 4.2 2.4
% % 9     5.9 4.2 2.5 ! mean: 4.19 - 4.20
% % 10    8.3 4.9 3.3







