function im_stokes = polar_stokes(im_set, pol_angle_flip)
% POLAR_STOKES calculates Stokes vectors for the input images set
%
% Usage: im_stokes = polar_stokes(im_set, pol_angle_flip)
% 
% Inputs:
%   im_set            - M x N x 3 x 5 double matrix with 5 images (V=0CCW/45CCW/H=90CCW/135CCW/V=180CCW degree filter orientation)
%   pol_angle_flip    - 1 x 1 bool, whether to flip polarisation directions in this image
%
% Outputs:
%   im_stokes         - M x N x 7 double, Stokes factors
%                       im_stokes(:, :, :, 1) is the best estimate of the 0th Stokes parameter (intensity),               calculated from S0 = ((I1 + I5)/2 + I2 + I3 + I4) / 2
%                       im_stokes(:, :, :, 2) is the best estimate of the 1st Stokes parameter (horizontal polarisation), calculated from S1 = (I3 - I1) / S0_1
%                       im_stokes(:, :, :, 3) is the best estimate of the 2nd Stokes parameter (vertical polarisation),   calculated from S2 = (I4 - I2) / S0_2
%                       im_stokes(:, :, :, 4) is another estimate of the 0th Stokes parameter (intensity),                calculated from S0_1 = (I1 + I3)
%                       im_stokes(:, :, :, 5) is another estimate of the 0th Stokes parameter (intensity),                calculated from S0_2 = (I2 + I4)
%                       im_stokes(:, :, :, 6) is another estimate of the 0th Stokes parameter (intensity),                calculated from S0_3 = (I3 + I5)
%                       im_stokes(:, :, :, 7) is another estimate of the 1th Stokes parameter (horizontal polarisation),  calculated from S1_2 = (I3 - I5) / S0_3

if nargin < 2 || isempty(pol_angle_flip), pol_angle_flip = false; end

im_stokes               = zeros(size(im_set, 1), size(im_set, 2), size(im_set, 3), 6); % pre-allocate
im_stokes(:, :, :, 1)   = ( (im_set(:, :, :, 1) + im_set(:, :, :, 5)) / 2 + im_set(:, :, :, 2) + im_set(:, :, :, 3) + im_set(:, :, :, 4) ) / 2; % 0th Stokes parameter (intensity) best estimate
im_stokes(:, :, :, 4)   = im_set(:, :, :, 1) + im_set(:, :, :, 3); % 0th Stokes parameter (intensity) 2nd estimate
im_stokes(:, :, :, 5)   = im_set(:, :, :, 2) + im_set(:, :, :, 4); % 0th Stokes parameter (intensity) 3rd estimate
im_stokes(:, :, :, 6)   = im_set(:, :, :, 3) + im_set(:, :, :, 5); % 0th Stokes parameter (intensity) 4th estimate
im_stokes(:, :, :, 2)   = (im_set(:, :, :, 3) - im_set(:, :, :, 1)) ./ im_stokes(:, :, :, 4); % 1st Stokes parameter (1=hor, -1=vert)
im_stokes(:, :, :, 3)   = (im_set(:, :, :, 4) - im_set(:, :, :, 2)) ./ im_stokes(:, :, :, 5); % 2nd Stokes parameter (1=45deg CCW from Horizontal, -1=135deg)
im_stokes(:, :, :, 7)   = (im_set(:, :, :, 3) - im_set(:, :, :, 5)) ./ im_stokes(:, :, :, 6); % 1st Stokes parameter (1=hor, -1=vert)

if pol_angle_flip
    im_stokes(:, :, :, [2 3 7]) = -im_stokes(:, :, :, [2 3 7]);
end
