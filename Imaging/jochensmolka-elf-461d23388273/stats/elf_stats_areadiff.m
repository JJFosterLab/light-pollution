function D = elf_stats_areadiff(R1, R2, E, dist_degree)
%
% Inputs:
%   R1, R2      - N x M double - spectral photon radiance for M samples along N elevations (in ph / sr / m^2 / nm)
%   E           - N x 1 double - elevation (must be evenly sampled)
%   dist_degree - degree of the distance function, e.g. 2 to use the sum of square of residuals (default: 1)
%
% elf_stats -> elf_stats_areatest -> elf_stats_areadiff -> elf_stats_normalise

%% Check inputs
if nargin < 1 || isempty(R1), R1 = rand(60, 1); warning('DEBUG mode: Using a random radiance signal!'); end
if nargin < 2 || isempty(R2), R2 = rand(60, 1); warning('DEBUG mode: Using a random radiance signal!'); end
if nargin < 3 || isempty(E), N = size(R1, 1); E = (-90+90/N:180/N:90-90/N)'; end % assume even sampling from -90 to 90 degrees
if nargin < 4, dist_degree = 1; end

if size(R1, 1) == 1 && size(R1, 2)>1
    warning('Input array R1 is a %dx%d row vector but should be a column vector; Using R'' instead.', size(R1, 1), size(R1, 2));
    R1 = R1';
end
if size(R2, 1) == 1 && size(R2, 2)>1
    warning('Input array R2 is a %dx%d row vector but should be a column vector; Using R'' instead.', size(R2, 1), size(R2, 2));
    R2 = R2';
end
if size(E, 1) == 1 && size(E, 2)>1
    warning('Input array E is a %dx%d row vector but should be a column vector; Using E'' instead.', size(E, 1), size(E, 2));
    E = E';
end

if size(R1, 1) ~= length(E)
    if size(R1, 2) == length(E)
        warning('The number of rows in R1 and E need to be identical. Using R1'' instead.');
        R1 = R1';
    else
        error('The number of rows in R1 and E need to be identical.');
    end
end
if size(R2, 1) ~= length(E)
    if size(R2, 2) == length(E)
        warning('The number of rows in R2 and E need to be identical. Using R2'' instead.');
        R2 = R2';
    else
        error('The number of rows in R2 and E need to be identical.');
    end
end
if any(size(R1) ~= size(R1))
    error('R1 and R2 need to be the same size.');
end

%% main
binwidth            = abs(median(diff(E)));
if max(abs(abs(diff(E))-binwidth))>0.01
    error('E needs to be evenly sampled');
end
elecorr             = cosd(E);                              % Elevation distortion factor
solidangle          = pi * elecorr * binwidth /180*pi;      % The solid angle filled by each elevation bin (in steradians, assuming 

R1                  = log10(elf_stats_normalise(R1, E));   % normalise radiances, then take the log
R2                  = log10(elf_stats_normalise(R2, E));   % normalise radiances, then take the log
D                   = squeeze(nansum(abs(R1 - R2).^dist_degree .* repmat(solidangle, [1 size(R1, 2)]), 1)); % calculate the total area between the mean curves for each repetition 
