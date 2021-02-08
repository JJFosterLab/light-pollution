function ah = polar_plot_dolp(dolp, ah, cb, a, e)
% POLAR_PLOT_DOLP plots degree of linear polarisation as gray scale.
%
% Usage: polar_plot_dolp(dolp, ah, cb, a, e)
%
% Inputs:
%   dolp        - M x N double, degree of linear polarisation for each image pixel, in degrees
%   ah          - axes handle to plot into; default: new axes in figure 828
%   cb          - 1x1 bool, whether to plot a colorbar; default: true
%   a           - M x N double, azimuth positions of the values in dolp
%   e           - M x N double, elevation positions of the values in dolp
%
% Outputs:
%   ah          - axes handle
%
% See also: polar_plot_aop, polar_plot_int, polar_plot_dolp3d, polar_plot_int3d, polar_plot_sphere3d

if nargin < 5 || isempty(e), e = linspace(90, -90, size(dolp, 1)); end
if nargin < 4 || isempty(a), a = linspace(-180, 180, size(dolp, 2)); end
if nargin < 3 || isempty(cb), cb = true; end
if nargin < 2 || isempty(ah), figure(828); ah = axes; elseif isnumeric(ah); figure(ah); clf; ah = axes; end

dolp(dolp>1) = 1;

imagesc(ah, a, e, dolp); 
axis xy equal off
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