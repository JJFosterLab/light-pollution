function ah = polar_plot_int(int, ah, a, e)
% POLAR_PLOT_INT plots intensity image.
%
% Usage: ah = polar_plot_int(int, ah, a, e)
%
% Inputs:
%   aop         - M x N x C double, intensity for each image pixel, in degrees
%   ah          - 1 x 1 axes handle to plot into; default: new axes in figure 827
%   a           - M x N double, azimuth positions of the values in dolp
%   e           - M x N double, elevation positions of the values in dolp
%
% Outputs:
%   ah          - axes handle
%
% See also: polar_plot_aop, polar_plot_dolp, polar_plot_dolp3d, polar_plot_int3d, polar_plot_sphere3d, polar_plot_aop3d

if nargin < 4 || isempty(e), e = linspace(90, -90, size(int, 1)); end
if nargin < 3 || isempty(a), a = linspace(-180, 180, size(int, 2)); end
if nargin < 2 || isempty(ah), figure(827); ah = axes; elseif isnumeric(ah); figure(ah); clf; ah = axes; end

image(ah, a, e, int);
axis xy equal off