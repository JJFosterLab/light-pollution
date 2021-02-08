function ah = polar_plot_aop(aop, ah, dolp, cb)
% POLAR_PLOT_AOP plots angle of polarisation as colour.
% If a DoLP (degree of linear polarisation) matrix is also provided, intensity is modulated by DoLP.
%
% Usage: polar_plot_aop(aop, ah, dolp, cb)
%
% Inputs:
%   aop     - M x N double, angle of polarisation for each image pixel, in degrees
%   ah      - axes handle to plot into; default: new axes in figure 827
%   dolp    - degree of linear polarisation for each pixel; if empty, function will not modulate intensity by degree of linear polarisation; default: do not modulate
%   cb      - 1x1 bool, whether to plot a colorbar; default: true
%
% Outputs:
%   None
%
% See also: polar_plot_dolp, polar_plot_int

if nargin<4 || isempty(cb), cb = true; end
if nargin<3 || isempty(dolp), dolp = ones(size(aop)); end
if nargin<2 || isempty(ah), figure(827); ah = axes; elseif isnumeric(ah); figure(ah); clf; ah = axes; end

hsvimage        = zeros(size(aop, 1), size(aop, 2));
hsvimage(:,:,1) = mod(aop, 180)/180;    % Hue
hsvimage(:,:,2) = 1;                    % Saturation
hsvimage(:,:,3) = dolp;                 % Intensity
rgb2            = hsv2rgb(hsvimage);    % Map to RGB colour space for display

imshow(rgb2, 'parent', ah)
colormap(ah, hsv(180));
axis(ah, 'image');
if cb
    cbh = colorbar(ah); 
    set(cbh, 'Ticks', [0:30:180]./180, 'TickLabels', num2str((0:30:180)'));
    ylabel(cbh,'Angle of linear polarisation (\circ)');
end