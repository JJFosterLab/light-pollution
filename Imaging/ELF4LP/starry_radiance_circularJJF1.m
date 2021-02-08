function starry_radiance_circularJJF1(filename, rotate, outputfolder)
%% starry_radiance_circularJJF.m
%---------------------------------------------------------------
%       USAGE: starry_radiance_circularJJF(filename, rotate, outputfolder) %
%
%      AUTHOR: Jochen Smolka                DATE: 2017 07 27
%    MODIFIED: James Foster                 DATE: 2017 09 08
%    
% DESCRIPTION: Descriptive statistics for a sky image in VEPS format. 
%              Adapted from milkyway_contrasts_circular.m
%
%      INPUTS:  1. Filename
%               2. Rotation angle for image
%               3. Folder for output images (PDF)
%     OUTPUTS:  .pdf figures   
%
%  REFERENCES:  Foster et al., 2017 http://dx.doi.org/10.1098/rstb.2016.0079
%---------------------------------------------------------------
%
% TO DO:
% - MIN & MAX 
% - CENTRAL TENDENCY
% - LOCATION OF EACH
% - NOTH--SOUTH RADIANCE PROFILE
% - EAST--WEST RADIANCE PROFILE
% - MILKY WAY : BACKGROUND CONTRAST
% - NORMALISE IMAGE TO VALUE
%
%---    JOCHEN SAYS...  ---%
% Loads prefiltered starry sky images in circular projection and calculates different measures of contrasts between the Northern and Southern arm
%
% Filename must include a cell array im_filt_HDR with four elements. By
% default, these are assumed to be filtered at 2/4/8/16 degrees, but this
% only matters for labels. The second input variable 'rotate' should be
% true if the Milky Way is vertical in the image, as the algorithm assumes
% a roughly horizontal orientation.
%---    --- --- --- 	---%

addpath('/Users/jamesf/Dropbox/Matlab/veps/useful');

%% DEBUGGERY
if ~exist('filename','var')
    filename = fullfile('/Users/jamesf/Documents/ExampleNEFs/JJFvepsTest/JJFmodified/JustUp_scene001_reproj.mat');%'4_osnabrck_dacke2013_scene001_reproj.mat');% 2_bremen_dacke2013_scene001_reproj.mat'); %6_aurora_diffuse_scene001_reproj.mat');% 5_aurora_north_scene001_reproj.mat');% 4_lightpollution_joburg_scene001_reproj.mat'); %3_lightpollution_odarslv_scene001_reproj.mat'); %2_milkyway_overcast_scene001_reproj.mat'); %1_milkyway_clouds_scene001_reproj.mat'); %Tvärminne-longer_scene001_reproj.mat'); %Tvärminne-short_scene001_reproj.mat'); %best_scene001_reproj.mat');% IllmitzGibbous_scene001_reproj.mat');% Fosteretal2017MilkyWay_scene001_reproj.mat');% Fosteretal2017MilkyWay_scene001_reproj.mat');% 3_lund_dacke2013_scene001_reproj.mat');% 3_lightpollution_odarslv_scene001_reproj.mat');% 4_norrevngsvg_scene001_reproj.mat');% 2_stonehenge_lowerelev_scene001_reproj.mat');% 3_nicama_lowerelev_scene001_reproj.mat');% 1_stonehenge_highelev_scene001_reproj.mat');% 1_joburg_dacke2013_scene001_reproj.mat');% 4_lightpollution_joburg_scene001_reproj.mat');%    fullfile(reprojfolder, allfiles(i).name);
end
filename = fullfile('/Users/jamesf/Documents/Starry Project/Light Pollution Imaging/LPexper-20181127-LowMWay-Thornwood/radiance/1Start_scene001_reproj.mat');

nmax =   3*10^10; 3*10^12; 1*10^11; 1*10^12;  %VALUE TO MAXIMISE TO

rotate = false; %SHOULD ALREADY BE IN CORRECT ORIENTATION
outputfolder = dir('/Users/jamesf/Documents/ExampleNEFs/1MilkyNEFs/JJFmodified/');%IS THIS HOW dir WORKS?
%% load filtered images
temp    = load(filename, 'im_filt_reproj');
ims     = temp.im_filt_reproj;
clear temp;
%% Calculate elevation for each pixel (this should work for all images now, since they should all be square and the same size)
% For non-square images, BUGFIX THIS!
mid    = [1+(size(ims{1}, 1)-1)/2; 1+(size(ims{1}, 2)-1)/2];        % centre of image
r_full = 8 * size(ims{1}, 2) / 24;                                 % theoretical value for 24mm high chip that is fully covered by fisheye circular image
[x, y] = meshgrid(1:size(ims{1}, 1), 1:size(ims{1}, 2));
r      = sqrt((x-mid(1)).^2 + (y-mid(2)).^2);
ele    = asind(r / 2 / r_full) * 2;%r/r_full * 90;
ele_s  = ele; ele_s(y<mid(2)) = -ele_s(y<mid(2)); % define everything on the left half as negative %I SHOULD MAKE THIS EAST
%% init
absmax      = cell(4, 4);
absmax_pos  = cell(4, 4);
absmax_ele  = cell(4, 4);
abscont     = zeros(4, 4);
linemax     = cell(4, 4);
linemax_ele = cell(4, 4);
linecont    = zeros(4, 4);
elecont     = zeros(4, 4);
%% calculate contrasts

for ch = 1:4
    %% 1. extract images
    for sc = 1:4
        im = ims{sc};
        im(isnan(im)) = 0;
        if rotate
            im = rot90(permute(im, [2 1 3]), 2);        % rotate image if necessary to place Milky Way horizontally %SHOULD BE POSSIBLE TO KEEP MW VERTICAL
        else
            im = fliplr(im);
        end
        if ch < 4
            im(:, :, [1:ch-1 ch+1:3]) = 0; 
            sumim{sc, ch}  = sum(im, 3);
        else
            sumim{sc, ch}  = mean(im, 3);
        end   % if colour channel, set other channels to 0 (works for plotting and contrast calculation)
        plotim{sc, ch} = im;                            % save back into ims structure, as this will later be plotted
    end
        %% 2. ABSCONT: For each point from left to right, get the maximum along the vertical dimension (only in the centre of the image to exclude horizon light) %MAY NOT BE NECESSARY
    for sc = 1:4
        [absmax{sc, ch}, absmax_pos{sc, ch}] = max(sumim{sc, ch}(401:600, :), [], 1); %SHOULD FREE UP THIS REGION, THE MILKY WAY (OR OTHER STREAK) CAN GO IN ANY DIRECTION
        absmax_x{sc, ch}    = 1:length(absmax{sc, ch});
        absmax_y{sc, ch}    = absmax_pos{sc, ch}+400; %WHY +400? 30DEG?
        linind              = sub2ind(size(ele_s), absmax_x{sc, ch}, absmax_y{sc, ch});
        absmax_ele{sc, ch}  = ele_s(linind); 

        % calculate Michelson contrast between the maximum in the South and the maximum in the North
        Nmax_abs            = max(absmax{sc, ch}(absmax_ele{sc, ch}>0));
        Smax_abs            = max(absmax{sc, ch}(absmax_ele{sc, ch}<0));
        abscont(sc, ch)     = abs((Nmax_abs - Smax_abs) / (Nmax_abs + Smax_abs)); %THESE ARE NOT THE CONTRASTS IN THE FINAL OUTPUT
    end
%% 3. LINECONT: Fit a line to the maxima in the most highly filtered image, and evaluate contrasts along that line
  
    % find the first and last value > 0.8 * median, and only use points in between (to cut out black borders)   
    sel = absmax{end, ch}>0.6*median(absmax{end, ch});  %I MIGHT WANT TO DO THIS IN A DIFFERENT WAY
    sel(find(sel, 1, 'first'):find(sel, 1, 'last')) = 1;
    
    allsels{ch} = sel;
    
    milkyfit = polyfit(absmax_x{end, ch}(sel), absmax_y{end, ch}(sel), 2);                                 % fit a quadratic function
   %CONSIDER milkyfit = polyfit(absmax_y{end, ch}(sel), absmax_x{end, ch}(sel), 2);
   
end
%%  DODGE THE LOOP
ch = 4;
%OR MAYBE USER IDENTIFIED MILKY WAY
   wh_im = plotim{1, 4};
   lowlim =  10.3; 11.5;10.5;  10.4; 11; 
   
   starlight_1 =  abs((wh_im(:,:,1)+ wh_im(:,:,2)+wh_im(:,:,3))) ;
   logstarlight_1 = log10(starlight_1);
   logstarlight_1(logstarlight_1 < lowlim) = lowlim;
   surf(logstarlight_1, 'EdgeColor','none'); colormap(jet); colorbar()
%    surf(log10(abs((wh_im(:,:,1)+ wh_im(:,:,2)+wh_im(:,:,3)))),'EdgeColor','none')
   figure(); 
   imshow(flipud(wh_im)./(nmax));
%    imshow(flipud(wh_im)./(3*10^10)); [brightMW, dimMW] = ginput(2);  %USER SELECTS BRIGHT AND DIM ENDS OF MILKY WAY       
%    px = %PIXELS AWAY FROM HORIZON
% elevat = asind((Px - 500.5) * 3 / 2000) * 2; %THE ELEVATION
xscale = sind((-60:30:60 +90) / 2) * 2000 / 3 +500.5;

ew_bright = figure();
plot(starlight_1(500, :))%EAST-WEST
axew = axis(gca);
% axis([0 1000 0 6*10^10])
% axis([-70 70 0 max([axns(4) 4e10])])
set(gca, 'xtick', xscale, 'xticklabels', {num2str([30 60 90 60 30]'), 'West'});

ns_bright = figure();
plot(flip(starlight_1(:, 500)))%NORTH-SOUTH
axew = axis(gca);
% axis([-70 70 0 max([axns(4) 4e10])])
set(gca, 'xtick', xscale, 'xticklabels', {num2str([30 60 90 60 30]'), 'South'});

skypix = starlight_1(starlight_1 ~= 0);
figure();
histogram((skypix),'EdgeColor','none')
set(gca,'xscale','log')
% 
% imshow( vec2mat([1000, 1000], starlight_1(logical(skypix))) )
% median(skypix)
% mean(skypix)
% min(skypix)
% max(skypix)

        