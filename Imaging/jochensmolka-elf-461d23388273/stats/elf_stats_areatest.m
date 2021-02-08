function [p, D0] = elf_stats_areatest(R1, R2, E, reps, dist_degree, plotall, blocksize)



%% Check inputs
if nargin < 1 || isempty(R1), R1 = repmat(cosd(-88.5:3:88.5)', [1, 10]) + 10 + 3*randn(60, 10); warning('DEBUG mode: Using a random radiance signal!'); end
%if nargin < 2 || isempty(R2), R2 = repmat(cosd(-88.5:3:88.5)', [1, 10]) + 10 + 3*randn(60, 10); warning('DEBUG mode: Using a random radiance signal!'); end
if nargin < 2 || isempty(R2), R2 = 5*(repmat((-88.5:3:88.5)'/30, [1, 10]) + 10 + 3*randn(60, 10)); warning('DEBUG mode: Using a random radiance signal!'); end
if nargin < 3 || isempty(E), N = size(R1, 1); E = (-90+90/N:180/N:90-90/N)'; end % assume even sampling from -90 to 90 degrees
if nargin < 4, reps = 100000; end
if nargin < 5, dist_degree = 2; end
if nargin < 6, plotall = true; end
if nargin < 7, blocksize = reps; end

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

%%TODO: check blocksize

%% parameters % pre-allocation
elebins             = size(R1, 1);
n1                  = size(R1, 2);
n2                  = size(R2, 2);

R                   = cat(2, R1, R2); % append both radiance arrays; during permutation, this array will be indexed with randomly permutated indices

R1i                 = zeros(elebins, n1, reps);
R2i                 = zeros(elebins, n2, reps);

%% main permutation
for i = 1:reps
    if mod(i, 10000)==0, fprintf('%dk..', i/1000); end
    if mod(i, 100000)==0, fprintf('\n'); end
    ind             = randperm(n1+n2);       % randomly permute which radiance curves belong to which group
    R1i(:, :, i)    = R(:, ind(1:n1));       % these radiance curves belong to group 1 in this rep
    R2i(:, :, i)    = R(:, ind(n1+1:end));   % these radiance curves belong to group 2 in this rep
end
fprintf('\n');

R1_m    = squeeze(nanmean(R1i, 2));             % average radiance curves for group 1 across scenes for each repetition
R2_m    = squeeze(nanmean(R2i, 2));             % average radiance curves for group 2 across scenes for each repetition
D       = elf_stats_areadiff(R1_m, R2_m, E, dist_degree);

R1_m0   = nanmean(R1, 2);                       % average radiance curves for group 1 across scenes for real data
R2_m0   = nanmean(R2, 2);                       % average radiance curves for group 2 across scenes for real data
D0      = elf_stats_areadiff(R1_m0, R2_m0, E, dist_degree); % calculate the total area between the mean curves for real data 

p       = nnz(D>=D0) / reps;                 % calculate p-value as the proportion of random permutation that resulted in an equal or higher area between curves

%% Plot
if plotall
    figure(71); clf;
    subplot(1, 3, 1);
    plot(elf_stats_normalise(R1, E), E);
    hold on;
    plot(elf_stats_normalise(R1_m0, E), E, 'r', 'linewidth', 3);
    ylabel('Elevation (\circ)');
    title('Env. 1');

    subplot(1, 3, 2);
    plot(elf_stats_normalise(R2, E), E);
    hold on;
    plot(elf_stats_normalise(R2_m0, E), E, 'b', 'linewidth', 3);
    xlabel('Norm. radiance');
    title('Env. 2');
    set(gca, 'YTick', []);

    subplot(1, 3, 3);
    plot(elf_stats_normalise(R1_m0, E), E, 'r', elf_stats_normalise(R2_m0, E), E, 'b', 'linewidth', 3);
    set(gca, 'YTick', []);
    
    
    if p<0.001, ps='***';
    elseif p<0.01, ps='**';
    elseif p<0.05, ps='*';
    else ps='ns';
    end
    title(sprintf('%s-deg. comp. (%s)', iptnum2ordinal(dist_degree), ps));
    
    figure(72); clf;
    hist(D, 50);
    hold on;
    line([D0 D0], [0 reps/50], 'color', [1 0 .4], 'linewidth', 2)
    xlabel('normalised area between curves');
    ylabel('count');

        title(sprintf('D_0 = %.2f; p = %.4g (%d permutations)', D0, p, reps))
end

