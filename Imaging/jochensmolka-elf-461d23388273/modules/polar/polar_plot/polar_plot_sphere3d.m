function ah = polar_plot_sphere3d(ah, col, r)
% POLAR_PLOT_SPHERE plots a monochromatic sphere.
%
% Usage: ah = polar_plot_sphere3d(ah, col, r)
%
% Inputs:
%   ah         - 1 x 1 axes handle to plot into; default: new axes in figure 833
%   col        - 1 x 3 double, RGB colour of sphere
%   r          - 1 x 1 double, radius of sphere to plot
%
% Outputs:
%   ah         - axes handle
%
% See also: polar_plot_dolp, polar_plot_int, polar_plot_aop, polar_plot_int3d, polar_plot_aop3d, polar_plot_dolp3d

if nargin < 3 || isempty(r), r = 1; end
if nargin < 2 || isempty(col), col = [.8 .8 .8]; end
if nargin < 1 || isempty(ah), figure(833); ah = axes; elseif isnumeric(ah); figure(ah); clf; ah = axes; end

[a, b, c] = sphere(100);
surf(ah, r*a, r*b, r*c, 'facecolor', col, 'edgecolor', 'none');