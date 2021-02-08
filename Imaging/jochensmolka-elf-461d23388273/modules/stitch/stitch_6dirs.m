function [finalx_im, finaly_im, iBestImage] = stitch_6dirs(I_info, angles, azi, ele, extraangle, images)
%
% I_info    - Image information structure
% angles:   angles for images 2:end relative to up-image
%           0 degrees means the image connects to the top of the up-image
%           90 degrees means the image connects to the right of the up-image
%           and so on
%
% Out: stitchedImage_unwarped 

if nargin < 4 || isempty(ele), ele = 90:-0.1:-90; end
if nargin < 3 || isempty(azi), azi = -180:0.1:180; end
if nargin < 2, angles = [0 90 180 270]; end

%%
%%%
%%%
%%%
% Down image integration is still not working!
warning('Down image integration is still not working!');
%%%
%%%
%%%

%% Calculate main projections in azi/ele (spherical coordinates) and X/Y/Z (cartesian coordinates)
[azi_grid, ele_grid]    = meshgrid(azi, ele);                                   % grid of desired spherical angles
[X, Y, Z]               = sph2cart(deg2rad(azi_grid), deg2rad(ele_grid), 1);    % grid of desired cartesion coordinates

%% Find best image to sample for each point
% For each point, find the distance to the optical axis of each of the available images; 
% the minimum distance marks the best image to sample for this point
refazi = zeros(1, length(angles)+1);
refele = zeros(1, length(angles)+1);
d      = nan(size(azi_grid, 1), size(azi_grid, 2), length(angles)+1);

% Calculate refazi/refele;
% first element of each one stays 0 for the up image
for i = 1:length(angles)
    if mod(angles(i), 360) <= 180
        refazi(i+1) = 90;
        refele(i+1) = 90-angles(i);
    else
        refazi(i+1) = -90;
        refele(i+1) = angles(i)-270;
    end
end
refazi(end+1) = 0;
refele(end+1) = 180;

for i = 1:length(refazi)
    d(:, :, i) = elf_support_sphdist(azi_grid, ele_grid, refazi(i), refele(i));
end

% top           = 0     = 90/90
% topright      = 45    = 90/45
% right         = 90    = 90/0
% bottomright   = 135   = 90/-45
% bottom        = 180   = 90/-90

% bottomleft    = 225   = -90/-45
% left          = 270   = -90/0
% bottomright   = 315   = -90/45

[~, iBestImage] = min(d, [], 3);

%% Calculate the sample matrix for each point for each image
% For the up image, use the unrotated image
[x_im(:, :, 1), y_im(:, :, 1)] = elf_project_cart2fisheye(X, Y, Z, I_info, 'default');

% For the other images, rotate X/Y/Z
for i = 1:length(angles)
    [newx, newy, newz] = rotate_periph2cent(X, Y, Z, angles(i), extraangle);
    [x_im(:, :, i+1), y_im(:, :, i+1)] = elf_project_cart2fisheye(newx, newy, newz, I_info, 'default');
end

[newx, newy, newz] = rotate_back2cent(X, Y, Z, 0, extraangle);
[x_im(:, :, 6), y_im(:, :, 6)] = elf_project_cart2fisheye(newx, newy, newz, I_info, 'default');


%% Now sample that matrix into only the best images %%%TODO: Make this more efficient!
tic
finalx_im = zeros(size(x_im, 1), size(x_im, 2));
finaly_im = zeros(size(x_im, 1), size(x_im, 2));
for i = 1:size(x_im, 1)
    for j= 1:size(x_im, 2)
        finalx_im(i, j) = x_im(i, j, iBestImage(i, j));
        finaly_im(i, j) = y_im(i, j, iBestImage(i, j));
    end
end
toc

return     


if rotated
    projection_ind       = elf_project_sub2ind([I_info.Width I_info.Height I_info.SamplesPerPixel], y_im, I_info.Width+1-x_im);
else
    projection_ind       = elf_project_sub2ind([I_info.Height I_info.Width I_info.SamplesPerPixel], x_im, y_im);
end

end %main

%% subfunctions

function [newx, newy, newz] = rotate_cent2periph(thisx, thisy, thisz, angle, extraangle)
    % Assumption: Sky is up in the image
    % 0 means image connects to the top, 90 to the right of up image
    [tempx, tempy, tempz]   = stitch_rot3D(thix, thisy, thisz, extraangle, 'x');   % before anything else, rotate around x again to correct for rotated camera
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
    
function [newx, newy, newz] = rotate_back2cent(thisx, thisy, thisz, angle, extraangle)
    % Assumption: Sky is up in the image
    % 0 means image connects to the top, 90 to the right of up image
    [tempx, tempy, tempz]   = stitch_rot3D(thisx, thisy, thisz, -(180 - angle), 'x');  % First, rotate around the optical axis (x-axis) to place in the right orientation
    [tempx, tempy, tempz]   = stitch_rot3D(tempx, tempy, tempz, -180, 'y');             % Second, flip it up from the bottom of the sphere with the sky part up
    [newx, newy, newz]      = stitch_rot3D(tempx, tempy, tempz, -extraangle, 'x');     % Finally, rotate around x again to correct for rotated camera
end  

% %% Code from calc_all:
% %% Stich images together
% figure(8); clf;
% figure(9); clf; 
% for im = 1:length(images)
%     thisint  = images{im};
%     thisint1 = thisint(:, :, 1); 
% 
%     thisazi  = deg2rad(linspace(-90, 90, size(thisint1, 1)));
%     thisele  = deg2rad(linspace(90, -90, size(thisint1, 2)));
%     [azi_grid, ele_grid] = meshgrid(thisazi, thisele);          % grid of desired angles
%     [thisx, thisy, thisz] = sph2cart(azi_grid, ele_grid, 1);    % tranform to CART
%     thisx = thisx(:);
%     thisx = reshape(thisx, size(thisy));
% 
%     figure(7); clf; surf(thisx, thisy, thisz, thisint1, 'edgecolor', 'none');
%     xlabel('x'); ylabel('y'); zlabel('z'); axis equal;  axis([-1 1 -1 1 -1 1]);
% 
%     
% 
% 
%     figure(8); surf(newx, newy, newz, thisint1, 'edgecolor', 'none'); view(90, 0);
%     xlabel('x'); ylabel('y'); zlabel('z'); axis equal;  axis([-1 1 -1 1 -1 1]);
% 
%     [newazi, newele, ~] = cart2sph(newx, newy, newz);
%     % Everything that had a negative x has to be changed in azi/ele (below the horizon)
%     sel = thisint1>0 & newx>0;
%     thisint1(~sel(:)) = NaN;
% 
% 
% 
% 
% 
%     figure(9);
% 
%     surf(rad2deg(newazi), rad2deg(newele), thisint1, 'edgecolor', 'none'); hold on;
%     %axis equal
%     %pause
% 
% 
% 
% 
% 
% end
% view(0, 90)
% colormap(hsv)
% colorbar