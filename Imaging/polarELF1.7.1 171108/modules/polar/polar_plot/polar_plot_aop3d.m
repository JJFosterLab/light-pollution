function ah = polar_plot_aop3d(pos, aop3d, ah, sc, r)
% POLAR_PLOT_AOP3D plots angle of polarisation as vectors on a sphere.
%
% Usage: ah = polar_plot_aop3d(pos, aop3d, ah, sc, r)
%
% Inputs:
%   pos        - M x 3 double, position on sphere, in 3d cartesian coordinates
%   aop        - M x 3 double, angle of polarisation for each image pixel
%   ah         - 1 x 1 axes handle to plot into; default: new axes in figure 832
%   sc         - 1 x 1 double, scaling factor for polarisation vectors
%   r          - 1 x 1 double, radius of sphere to attach vectors to. Multiplies with the radius implied in pos.
%
% Outputs:
%   ah         - axes handle
%
% See also: polar_plot_dolp, polar_plot_int, polar_plot_aop, polar_plot_dolp3d, polar_plot_int3d, polar_plot_sphere3d

if nargin < 5 || isempty(r), r = 1; end    % Radius at which to attach vectors
if nargin < 4 || isempty(sc), sc = 1; end  % Scaling factor for vectors
if nargin < 3 || isempty(ah), figure(832); ah = axes; elseif isnumeric(ah); figure(ah); clf; ah = axes; end

quiver3(r*pos(:, 1), r*pos(:, 2), r*pos(:, 3), sc*aop3d(:, 1), sc*aop3d(:, 2), sc*aop3d(:, 3), 0, 'r', 'linewidth', 2);
quiver3(r*pos(:, 1), r*pos(:, 2), r*pos(:, 3), -sc*aop3d(:, 1), -sc*aop3d(:, 2), -sc*aop3d(:, 3), 0, 'r', 'linewidth', 2);

text(1.2, 0, 0, 'N')
text(-1.2, 0, 0, 'S')
text(0, 1.2, 0, 'W')
text(0, -1.2, 0, 'E')

xlabel('x'); ylabel('y'); zlabel('z');
axis equal
axis([-1.5 1.5 -1.5 1.5 -1.5 1.5])
grid on










