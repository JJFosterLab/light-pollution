function ah = polar_plot_int3d(int, ah, r, a, e)
% POLAR_PLOT_INT3D plots intensity image on a sphere.
%
% Usage: ah = polar_plot_int3d(int, ah, r, a, e)
%
% Inputs:
%   int        - M x N x 3 double, intensity for each image pixel
%   ah         - 1 x 1 axes handle to plot into; default: new axes in figure 830
%   r          - 1 x 1 double, radius of sphere to plot
%   a          - M x N double, grid of azimuth values for each value in int 
%   e          - M x N double, grid of elevation values for each value in int
%
% Outputs:
%   ah         - axes handle
%
% See also: polar_plot_dolp, polar_plot_int, polar_plot_aop, polar_plot_dolp3d, polar_plot_aop3d, polar_plot_sphere3d

if nargin < 5 || isempty(a) || isempty(e), [a, e] = meshgrid(linspace(-180, 180, size(int, 2)), linspace(90, -90, size(int, 1))); end
if nargin < 3 || isempty(r), r = 1; end
if nargin < 2 || isempty(ah), figure(830); ah = axes; elseif isnumeric(ah); figure(ah); clf; ah = axes; end

[x,y,z] = sph2cart(deg2rad(a), deg2rad(e), 1);
[x,y,z] = stitch_rot3D(x,y,z, -90, 'y');
[x,y,z] = stitch_rot3D(x,y,z, -90, 'z');
surf(x, y, z, int, 'edgecolor', 'none');