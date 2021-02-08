function ah = polar_plot_dolp3d(dolp, ah, cb, r, a, e, updir)
% POLAR_PLOT_DOLP3D plots degree of linear polarisation as colour on a sphere.
%
% Usage: ah = polar_plot_dolp3d(dolp, ah, cb, r, a, e)
%
% Inputs:
%   dolp       - M x N x 3 double, degree of linear polarisation for each image pixel
%   ah         - 1 x 1 axes handle to plot into; default: new axes in figure 831
%   cb         - 1 x 1 bool, whether to plot a colorbar; default: false
%   r          - 1 x 1 double, radius of sphere to plot
%   a          - M x N double, grid of azimuth values for each value in int 
%   e          - M x N double, grid of elevation values for each value in int
%
% Outputs:
%   ah         - axes handle
%
% See also: polar_plot_dolp, polar_plot_int, polar_plot_aop, polar_plot_int3d, polar_plot_aop3d, polar_plot_sphere3d

if nargin < 7, updir = 'N'; end
if nargin < 6 || isempty(a) || isempty(e), [a, e] = meshgrid(linspace(-180, 180, size(dolp, 2)), linspace(90, -90, size(dolp, 1))); end
if nargin < 4 || isempty(r), r = 1; end
if nargin < 3 || isempty(cb), cb = 0; end
if nargin < 2 || isempty(ah), figure(831); ah = axes; elseif isnumeric(ah); figure(ah); clf; ah = axes; end

[x,y,z] = sph2cart(deg2rad(a), deg2rad(e), r);
[x,y,z] = stitch_rot3D(x,y,z, -90, 'y');
switch updir
    case 'N'
        % do nothing
    case 'W'
        [x,y,z] = stitch_rot3D(x,y,z, 90, 'z');
    case 'E'
        [x,y,z] = stitch_rot3D(x,y,z, -90, 'z');
    case 'S'
        [x,y,z] = stitch_rot3D(x,y,z, 180, 'z');
end


dolp(dolp>1) = 1;

surf(x, y, z, dolp, 'edgecolor', 'none');
colormap(ah, jet(180));

if cb
    drawnow;
    a = get(ah, 'Position');
    cbh = colorbar(ah);
    set(cbh, 'Ticks', 0:0.2:1, 'TickLabels', num2str((0:20:100)'));
    ylabel(cbh,'Degree of linear polarisation (%)');
    set(ah, 'Position', a);
    drawnow;
end