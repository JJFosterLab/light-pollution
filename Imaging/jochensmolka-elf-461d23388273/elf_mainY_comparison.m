function [p, D0] = elf_mainY_comparison(dataset1, dataset2, reps, imgformat)
% ELF_MAINY_COMPARISON statistically compares two data sets
%
% Loads files: individual results files in mat folder, info files

%% check inputs
if nargin < 4 || isempty(imgformat), imgformat = '*.dng'; end
if nargin < 3 || isempty(reps), reps = 100000; end
if nargin < 2 || isempty(dataset1) || isempty(dataset2), error('You have to provide two valid dataset names'); end

                    elf_support_logmsg('\b\b\b\b\b\b\b\b\b\b\b\b\b\n');
                    elf_support_logmsg('----- ELF Step Y: Comparing environments -----\n');
                    elf_support_logmsg('      Env. 1: %s\n', dataset1);
                    elf_support_logmsg('      Env. 2: %s\n', dataset2);

%% Set up paths
elf_paths;

%% Load data
% Environment 1
para1            = elf_para('', dataset1, imgformat);
data1            = elf_readwrite(para1, 'loadcollres');

% Environment 2
para2            = elf_para('', dataset2, imgformat);
data2            = elf_readwrite(para2, 'loadcollres');

%% Compare statistically and plot stats
[p, D0] = elf_stats(data1, data2, 4, reps, 2, true);

%% Check: D0 from this should be the same as from this calculation:
% meandata1        = elf_readwrite(para1, 'loadmeanres');
% meandata2        = elf_readwrite(para2, 'loadmeanres');
% E           = meandata1(1).totalint.region_meanele;
% R1_allch    = meandata1.int.means;
% R1          = R1_allch(4, :)';
% R2_allch    = meandata2.int.means;
% R2          = R2_allch(4, :)';
% D0          = elf_stats_areadiff(R1, R2, E);
