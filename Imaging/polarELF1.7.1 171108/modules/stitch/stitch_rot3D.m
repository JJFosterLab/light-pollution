function [outx, outy, outz] = stitch_rot3D(x, y, z, theta, ax)

% Rotates x/y/z by theta degrees around axis x.
% assumptions: x, y, z same size
switch ax
    case 'x'
        R = [1 0 0;
             0 cosd(theta) -sind(theta);
             0 sind(theta) cosd(theta)];       
        
    case 'y'
        R = [cosd(theta) 0 sind(theta);
             0 1 0;
             -sind(theta) 0 cosd(theta)];
        
    case 'z'
        R = [cosd(theta) -sind(theta) 0;
             sind(theta) cosd(theta) 0;
             0 0 1];
        
end

tempmat = cat(2, x(:), y(:), z(:));
outmat  = (R*tempmat')';
outx = reshape(outmat(:, 1), size(x));
outy = reshape(outmat(:, 2), size(y));
outz = reshape(outmat(:, 3), size(z));

