function [finalw_im, finalh_im, iBestImage, w_im, h_im] = stitch_5dirs(I_info, angles, azi, ele, extraangle, method, exclimit)
%
% calculates index matrices to stitch together 5 fisheye images into a single wide-angle rect
% 
% I_info     - Image information structure (this function and sub-functions only use Height, Width, Model and FocalLength)
% angles:    - angles for images 2:end relative to up-image
%              0 degrees means the image connects to the top of the up-image
%              90 degrees means the image connects to the right of the up-image
%              and so on
% azi        - desired azimuth vector of the output image (default -180:0.1:180)
% ele        - desired elevation vector of the output image (default 90:-0.1:-90)
% extraangle - this is the angle by which the raw image has to be rotated clockwise to place the sky upwards (default 0 for method rect, where it should be done before)
% method     - 'fisheye' to produce index vectors into fisheye images (default)
%              'rect'    to produce index vectors into equirectangular (azimuth/elevation) images
% exclimit   - excentricity limit (in degrees), any points further away from the optical axis wil not be considered
%
% Out: stitchedImage_unwarped 

if nargin < 7 || isempty(exclimit),     exclimit = 90;              end
if nargin < 5 || isempty(extraangle),   extraangle = 0;             end
if nargin < 4 || isempty(ele),          ele = 90:-0.1:-90;          end
if nargin < 3 || isempty(azi),          azi = -180:0.1:180;         end
if nargin < 2 || isempty(angles),       angles = [0 90 180 270];    end

%% Calculate main projections in azi/ele (spherical coordinates) and X/Y/Z (cartesian coordinates)
[azi_grid, ele_grid]    = meshgrid(azi, ele);                                   % grid of desired spherical angles
[X, Y, Z]               = sph2cart(deg2rad(azi_grid), deg2rad(ele_grid), 1);    % grid of desired Cartesian coordinates

%% Find best image to sample for each point
% For each point, find the distance to the optical axis of each of the available images; 
% the minimum distance marks the best image to sample for this point
refazi = zeros(1, length(angles)+1);
refele = zeros(1, length(angles)+1);
d      = nan(size(azi_grid, 1), size(azi_grid, 2), length(angles)+1);

% Calculate refazi/refele, the optical axes for each image
% first element of each one stays 0 to represent the "up" image
for i = 1:length(angles)
    if mod(angles(i), 360) <= 180
        refazi(i+1) = 90;
        refele(i+1) = 90-angles(i);
    else
        refazi(i+1) = -90;
        refele(i+1) = angles(i)-270;
    end
end

for i = 1:length(refazi)
    d(:, :, i) = elf_support_sphdist(azi_grid, ele_grid, refazi(i), refele(i));
end

% IMAGE           ANGLES  REFAZI/REFELE
% -------------------------------------
% top             0       90/90
% topright        45      90/45
% right           90      90/0
% bottomright     135     90/-45
% bottom          180     90/-90

% bottomleft      225     -90/-45
% left            270     -90/0
% bottomright     315     -90/45

[~, iBestImage] = min(d, [], 3);

%% Control figure
% % figure(23); clf; hold on;
% % cols = {'k', 'r', 'g', 'b', 'y'};
% % for i = 1:length(refazi)
% % %     plot3(X(iBestImage==i), Y(iBestImage==i), Z(iBestImage==i), 'color', cols{i}, 'linestyle', 'none', 'marker', '.');
% %     plot3(X(iBestImage==i & d(:, :, i)<=60), Y(iBestImage==i & d(:, :, i)<=60), Z(iBestImage==i & d(:, :, i)<=60), 'color', cols{i}, 'linestyle', 'none', 'marker', '.');
% %     plot3(X(iBestImage==i & d(:, :, i)>60), Y(iBestImage==i & d(:, :, i)>60), Z(iBestImage==i & d(:, :, i)>60), 'color', 'm', 'linestyle', 'none', 'marker', '.');
% % end
% % axis equal; grid on; view(60, 20);
% % return;
%%%

%% Calculate the sample matrix for each point for each image
switch method
    case {'fisheye', '', 'default'}
        % For the up image, rotate by 180 degrees
        [upx, upy, upz]                = stitch_rot3D(X, Y, Z, 180+extraangle, 'x');
        [w_im(:, :, 1), h_im(:, :, 1)] = elf_project_cart2fisheye(upx, upy, upz, I_info, 'default');

        % For the other images, rotate X/Y/Z
        for i = 1:length(angles)
            [newx, newy, newz] = rotate_periph2cent(X, Y, Z, angles(i), extraangle);
            [w_im(:, :, i+1), h_im(:, :, i+1)] = elf_project_cart2fisheye(newx, newy, newz, I_info, 'default');
        end
        
    case {'rect', 'equirectangular', 'azel', 'sph'}
        % For the up image, rotate by 180 degrees
        [upx, upy, upz] = stitch_rot3D(X, Y, Z, 180+extraangle, 'x');
        [temp1, temp2]  = cart2sph(upx, upy, upz);
        w_im(:, :, 1)   = rad2deg(temp1);
        h_im(:, :, 1)   = rad2deg(temp2);
        % For the other images, rotate X/Y/Z
        for i = 1:length(angles)
            [newx, newy, newz]  = rotate_periph2cent(X, Y, Z, angles(i), extraangle);
            [temp1, temp2]      = cart2sph(newx, newy, newz);
            w_im(:, :, i+1)     = rad2deg(temp1);
            h_im(:, :, i+1)     = rad2deg(temp2);
        end
    otherwise
        error('Unknown method: %s', method);
end

%% Control figure with images
% % figure(24); clf; hold on;
% % azi = -90:0.1:90;
% % ele = 90:-0.1:-90;
% % 
% % for i = 1:length(refazi)
% %     [azi_grid, ele_grid]    = meshgrid(azi, ele);                                   % grid of desired spherical angles
% %     [X, Y, Z]               = sph2cart(deg2rad(azi_grid), deg2rad(ele_grid), 1);    % grid of desired cartesion coordinates
% % 
% % 
% %     [newx, newy, newz] = rotate_cent2periph(XX, YY, ZZ, angles(i), extraangle);
% %     
% %     
% % %     plot3(X(iBestImage==i), Y(iBestImage==i), Z(iBestImage==i), 'color', cols{i}, 'linestyle', 'none', 'marker', '.');
% %     plot3(X(iBestImage==i & d(:, :, i)<=60), Y(iBestImage==i & d(:, :, i)<=60), Z(iBestImage==i & d(:, :, i)<=60), 'color', cols{i}, 'linestyle', 'none', 'marker', '.');
% %     plot3(X(iBestImage==i & d(:, :, i)>60), Y(iBestImage==i & d(:, :, i)>60), Z(iBestImage==i & d(:, :, i)>60), 'color', 'm', 'linestyle', 'none', 'marker', '.');
% % end
% % axis equal; grid on; view(60, 20);
% % return;

%% Apply excentricity limit
w_im(d>exclimit) = NaN;

%% Now sample that matrix into only the best images %%TODO: Make this more efficient!
finalw_im = zeros(size(w_im, 1), size(w_im, 2));
finalh_im = zeros(size(w_im, 1), size(w_im, 2));
for i = 1:size(w_im, 1)
    for j= 1:size(w_im, 2)
        finalw_im(i, j) = w_im(i, j, iBestImage(i, j));
        finalh_im(i, j) = h_im(i, j, iBestImage(i, j));
    end
end

end %main

%% subfunctions

function [newx, newy, newz] = rotate_cent2periph(thisx, thisy, thisz, angle, extraangle)
    % Assumption: Sky is up in the image
    % 0 means image connects to the top, 90 to the right of up image
    [tempx, tempy, tempz]   = stitch_rot3D(thix, thisy, thisz, extraangle, 'x');   % Before anything else, rotate around x to correct for rotated camera
    [tempx, tempy, tempz]   = stitch_rot3D(tempx, tempy, tempz, 90, 'y');          % First, flip it down to the bottom of the sphere with the sky part out
    [newx, newy, newz]      = stitch_rot3D(tempx, tempy, tempz, 180 - angle, 'x'); % Second, rotate around the optical axis (x-axis) to place in the right orientation
end

function [newx, newy, newz] = rotate_periph2cent(thisx, thisy, thisz, angle, extraangle)
    % Assumption: Sky is up in the image
    % 0 means image connects to the top, 90 to the right of up image
    [tempx, tempy, tempz]   = stitch_rot3D(thisx, thisy, thisz, -(180 - angle), 'x');  % First, rotate around the optical axis (x-axis) to place in the right orientation
    [tempx, tempy, tempz]   = stitch_rot3D(tempx, tempy, tempz, -90, 'y');             % Second, flip it up from the bottom of the sphere with the sky part up
    [newx, newy, newz]      = stitch_rot3D(tempx, tempy, tempz, -extraangle, 'x');     % Finally, rotate around x again to correct for rotated camera
end   
 