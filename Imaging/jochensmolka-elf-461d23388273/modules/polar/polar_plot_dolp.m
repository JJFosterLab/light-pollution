function ah = polar_plot_dolp(dolp, ah)
% POLAR_PLOT_DOLP plots degree of linear polarisation as gray scale.
%
% Usage: polar_plot_aop(aop, ah, dolp, cb)
%
% Inputs:
%   dolp    - M x N double, degree of linear polarisation for each image pixel, in degrees
%   ah      - axes handle to plot into; default: new axes in figure 828
%
% Outputs:
%   None
%
% See also: polar_plot_aop, polar_plot_int

if nargin<2 || isempty(ah), figure(828); ah = axes; elseif isnumeric(ah); figure(ah); clf; ah = axes; end

imshow(dolp, 'parent', ah); 
axis(ah, 'image');
colormap(ah, jet(180));

cbh = colorbar(ah);
set(cbh, 'Ticks', 0:0.2:1, 'TickLabels', num2str((0:20:100)'));
ylabel(cbh,'Degree of linear polarisation (%)');