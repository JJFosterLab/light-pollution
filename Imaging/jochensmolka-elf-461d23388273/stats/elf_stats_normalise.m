function R_norm = elf_stats_normalise(R, E)
% elf_STATS_NORMALISE scales the radiance values in input array R so that
% the total irradiance from the sampled hemisphere is 1.
%
% Inputs:
%   R - N x M double - spectral photon radiance for M samples along N elevations (in ph / sr / m^2 / nm)
%   E - N x 1 double - elevation (must be evenly sampled)

%% Check inputs
if nargin < 1 || isempty(R), R = rand(60, 1); warning('DEBUG mode: Using a random radiance signal!'); end
if nargin < 2 || isempty(E), N = size(R, 1); E = (-90+90/N:180/N:90-90/N)'; end % assume even sampling from -90 to 90 degrees

if size(R, 1) == 1 && size(R, 2)>1
    warning('Input array R is a %dx%d row vector but should be a column vector; Using R'' instead.', size(R, 1), size(R, 2));
    R = R';
end
if size(E, 1) == 1 && size(E, 2)>1
    warning('Input array E is a %dx%d row vector but should be a column vector; Using E'' instead.', size(E, 1), size(E, 2));
    E = E';
end
if size(R, 1) ~= length(E)
    if size(R, 2) == length(E)
        warning('The number of rows in R and E need to be identical. Using R'' instead.');
        R = R';
    else
        error('The number of rows in R and E need to be identical.');
    end
end

%% main
binwidth                = abs(median(diff(E)));
if max(abs(abs(diff(E))-binwidth))>0.01
    error('E needs to be evenly sampled');
end
R_norm                  = zeros(size(R));

elecorr                 = cosd(E);                          % Elevation distortion factor
solidangle              = pi * elecorr * binwidth /180*pi;  % The solid angle filled by each elevation bin (in steradians, assuming 

if size(R, 2) == 1
    irradiance_per_bin  = R .* solidangle;                  % Spectral photon irradiance provided by each bin (ph / m^2 / nm)
    irradiance_total    = nansum(irradiance_per_bin);          % Total irradiance coming from the sampled hemisphere
    R_norm              = R / irradiance_total;             % Normalised spectral photon radiance (sum is 1)
else
    n = size(R, 2);
    irradiance_per_bin  = R .* repmat(solidangle, [1 n]);               % Spectral photon irradiance provided by each bin (ph / m^2 / nm)
    irradiance_total    = nansum(irradiance_per_bin, 1);                   % Total irradiance coming from the sampled hemisphere
    for i = 1:n
        R_norm(:, i)    = R(:, i) / irradiance_total(i); % Normalised spectral photon radiance (sum is 1)
    end
end

%% Checks
% sum(solidangle)-2*pi      % This control value approaches 0 for infinitely dense sampling
% sum(R_norm .* solidangle) % Should be 1