function [p, D0] = elf_stats(data1, data2, ch, reps, dist_degree, plotall)

%% Check inputs
if nargin < 3, ch = 4; end
if nargin < 4, reps = 100000; end
if nargin < 5, dist_degree = 2; end
if nargin < 6, plotall = true; end

%% Extract elevations E from both environments and make sure they are the same
E1           = data1(1).totalint.region_meanele;
E2           = data2(1).totalint.region_meanele;
if any(size(E1)~=size(E2)) || max(abs(E1-E2))>0.01
    error('Environments must be sampled at the same elevations');
end

%% Extract radiance data from both environments
R1_allch    = sub_extract(data1, 'int', 'median');
R1          = squeeze(R1_allch(ch, :, :));
R2_allch    = sub_extract(data2, 'int', 'median');
R2          = squeeze(R2_allch(ch, :, :));

%% Calculate and plot statistics
[p, D0]     = elf_stats_areatest(R1, R2, E1, reps, dist_degree, plotall);

end % main

% sub functions
function [outmat, meandim] = sub_extract(s, f1, f2)
    oldndims    = ndims(s(1).(f1).(f2));
    meandim     = oldndims + 1;
    temp        = arrayfun(@(x) x.(f1).(f2), s, 'UniformOutput', false);
    outmat      = cat(meandim, temp{:});
end