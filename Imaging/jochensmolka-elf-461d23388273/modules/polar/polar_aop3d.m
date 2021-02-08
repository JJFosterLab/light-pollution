function [aop3d_surf, pos] = polar_aop3d(aopim, iBestImage, dolpim, a, e)
% POLAR_AOP3D converts aop measurements (as shot) into 3d vectors
%
% Usage: [aop3d_surf, pos] = polar_aop3d(aop, iBestImage, dolp, a, e)
%
% Inputs:
%   aopim           - M x N double, AoP image covering the whole visual field
%   iBestImage      - M x N double, index matrix describing which image each point has been sampled from
%   dolpim          - M x N double, DoLP image covering the whole visual field
%   a               - 1 x P double, azimuths of points to be sampled, in degrees
%   e               - 1 x P double, elevations of points to be sampled, in degrees
%
% Outputs:
%   aop3d_surf      - P x 3 double, polarisation vectors for each point in a/e, projected onto unit sphere, in Cartesian coordinates
%   pos             - P x 3 double, position vectors for each point in a/e on unit sphere, in Cartesian coordinates

if nargin<5, [a, e] = polar_gridsphere; end % returns a azimuth/elevation grid that covers the sphere evenly
if nargin<3 || isempty(dolpim), dolpim = ones(size(aopim)); end

%% Define some constants, and find the cartesian grid 
oa    = [0 0 1; % the optical axes of the 5 images
         1 0 0; % North
         0 1 0; % West
         -1 0 0; % South
         0 -1 0];
az_uv = [0 1 0; % the 3d azimuth unit vectors for the 5 images
         0 -1 0;
         1 0 0;
         0 1 0;
         -1 0 0];
el_uv = [1 0 0; % the 3d elevation unit vectors for the 5 images
         0 0 1;
         0 0 1;
         0 0 1;
         0 0 1];
xoffset = (size(iBestImage, 1)+1)/2;
yoffset = (size(iBestImage, 2)+1)/2;
xres    = (size(iBestImage, 1)-1)/180;
yres    = (size(iBestImage, 2)-1)/360;
x       = round(e*xres + xoffset);     % Find the desired points in the image
y       = round(a*yres + yoffset);

[posx, posy, posz] = sph2cart(deg2rad(a), deg2rad(e), 1);       % this is 1/0/0 for the centre of image 1
[posx, posy, posz] = stitch_rot3D(posx, posy, posz, -90, 'y');  % this is 0/0/1 for the centre of image 1
pos = [posx(:) posy(:) posz(:)];

%% Now, project vectors onto sphere
aop3d_im   = zeros(length(x), 3); % pre-allocate
aop3d_surf = zeros(length(x), 3); % pre-allocate
for i = 1:length(x)
    % Calculate the component vectors along the image unit vectors
    aop3d_im(i, 1:3) = az_uv(iBestImage(x(i), y(i)), :) * cosd(aopim(x(i), y(i))) * dolpim(x(i), y(i))   +   el_uv(iBestImage(x(i), y(i)), :) * sind(aopim(x(i), y(i))) * dolpim(x(i), y(i));

    % Calculate the surface normal unit vector at this spot
    thispos_uv = pos(i, :); % this is 0/0/1 for the centre of image 1
    
    % The surface projection of the aop vector is the vector rejection of aop on the surface normal unit vector
    aop3d_surf(i, 1:3) = aop3d_im(i, 1:3) - dot(aop3d_im(i, 1:3), thispos_uv) * thispos_uv;
    
%     % This next section was meant to be used to correct dolp for successive shortening during projections.
%     % If such a correction is neccesary, this is not the way to do it.
%     costheta = dot(oa(iBestImage(x(i), y(i)), 1:3), thispos_uv);
%     ang(i) =  acosd(costheta); %only for diagnostics
%     dolp_corr(i) = dolp(x(i), y(i)) / costheta^2;
%     aopvec_surf2(i, 1:3) = aop3d_surf(i, 1:3) / costheta^2; % The reverse projection is cos2(theta) larger

% Use this for diagnostics to only include points from a certain image
%     if iBestImage(x(i), y(i)) ~= 1
%         aop3d_surf(i, 1:3) = 0;
%     end
end

%% diagnostic plots
% % For each point in a/e, plot the angle between that point and the optical axis of the image it originates  from
% figure(928); clf;
% scatter(a, e, 50, real(ang), 'filled');
% xlabel('azimuth'); ylabel('elevation');
% colormap(jet)
% colorbar
% clabel('Angle to image''s optical axis');
% 
% % For each point in a/e, 
% figure(929); clf;
% scatter(a, e, 50, dolp_corr, 'filled');
% xlabel('azimuth'); ylabel('elevation');
% colormap(jet)
% colorbar
% clabel('Angle to image''s optical axis');