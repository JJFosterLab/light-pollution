function ah = polar_plot_aop(aop, ah, dolp, cb, a, e)
% POLAR_PLOT_AOP plots angle of polarisation as colour.
% If a DoLP (degree of linear polarisation) matrix is also provided, intensity is modulated by DoLP.
%
% Usage: ah = polar_plot_aop(aop, ah, dolp, cb, a, e)
%
% Inputs:
%   aop        - M x N double, angle of polarisation for each image pixel, in degrees
%   ah         - 1 x 1 axes handle to plot into; default: new axes in figure 829
%   dolp       - M x N double, degree of linear polarisation for each pixel; if empty, function will not modulate intensity by degree of linear polarisation; default: do not modulate
%   cb         - 1x1 bool, whether to plot a colorbar; default: true
%   a          - M x N double, azimuth positions of the values in aop
%   e          - M x N double, elevation positions of the values in aop
%
% Outputs:
%   ah         - axes handle
%
% See also: polar_plot_dolp, polar_plot_int, polar_plot_dolp3d, polar_plot_int3d, polar_plot_sphere3d, polar_plot_aop3d

if nargin < 6 || isempty(e), e = linspace(90, -90, size(aop, 1)); end
if nargin < 5 || isempty(a), a = linspace(-180, 180, size(aop, 2)); end
if nargin < 4 || isempty(cb), cb = true; end
if nargin < 3 || isempty(dolp), dolp = ones(size(aop)); end
if nargin < 2 || isempty(ah), figure(829); ah = axes; elseif isnumeric(ah); figure(ah); clf; ah = axes; end

hsvimage        = zeros(size(aop, 1), size(aop, 2));
hsvimage(:,:,1) = mod(aop, 180)/180;    % Hue
hsvimage(:,:,2) = 1;                    % Saturation
hsvimage(:,:,3) = dolp;                 % Intensity
rgb2            = hsv2rgb(hsvimage);    % Map to RGB colour space for display

image(ah, a, e, rgb2);
axis xy equal off
colormap(ah, hsv(180));

if cb
    drawnow;
    a = get(ah, 'Position');
    cbh = colorbar(ah); 
    set(cbh, 'Ticks', [0:30:180]./180, 'TickLabels', num2str((0:30:180)'));
    ylabel(cbh,'Angle of linear polarisation (\circ)');
    set(ah, 'Position', a);
    drawnow;
end
