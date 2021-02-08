function D = elf_mainZ_compareall(datasets, imgformat, testall)
% elf_MAIN6_SUMMARY calculates the mean result for a dataset, and saves as MAT 
%
% Uses: elf_support_logmsg, elf_paths, elf_para, elf_para_update, 
%       elf_info_collect, elf_readwrite, elf_analysis_datasetmean
%
% Loads files: individual results files in mat folder, info files
% Saves files: mean results file in mat folder
% 
% Typical timing for a 50-scene environment (on ELFPC):
%       1s total

%% check inputs
if nargin < 3, testall = false; end
if nargin < 2 || isempty(imgformat), imgformat = '*.dng'; end % verbose determines whether each individual image is plotted during the process, and thumbs are provided at the end
if nargin < 1 || isempty(datasets), error('You have to provide at least two valid dataset names'); end 

                    elf_support_logmsg('\b\b\b\b\b\b\b\b\b\b\b\b\b\n');
                    elf_support_logmsg('----- ELF Step Z: Comparing all datasets -----\n');

%% Parameters
dist_degree = 2;
                    
%% Set up paths and file names; read info, infosum and para
elf_paths;

for i = 1:length(datasets)
    para            = elf_para('', datasets{i}, imgformat);
    if exist(para.paths.datapath, 'file') % test if this data set exists
        try
            meandata  = elf_readwrite(para, 'loadmeanres');
        catch me, warning(me.message); R(:, i) = NaN; continue;
        end
    else % check if it is a combined dataset
        try
            meandata  = elf_readwrite(para, 'loadmeanres_comb');
        catch me, warning(me.message); R(:, i) = NaN; continue;
        end
    end 
    
    E           = meandata(1).totalint.region_meanele;
    R1_allch    = meandata.int.median;%means;
    R(:, i)     = R1_allch(4, :)';
end

D = nan(size(R, 2), size(R, 2));
p = D;
for i = 1:size(R, 2)
    i
    for j = 1:size(R, 2)
        D(i, j) = elf_stats_areadiff(R(:, i), R(:, j), E, dist_degree);
        if testall
            p(i, j) = elf_mainY_comparison(datasets{i}, datasets{j}, 10000, imgformat);
        end
    end
end

figure(88); clf;
imagesc(D);
colormap jet
colorbar
figure(89); clf;
imagesc(p);
colormap jet
colorbar
% R(:, i+1) = cosd(linspace(0, 720, 60));
% datasets{end+1} = 'outgroup';

%%
figure(99); clf;
Z = linkage(R', 'weighted', {@nested_areadiff});
[~, ~, T] = dendrogram(Z, 0, 'Orientation', 'left', 'Labels', datasets, 'colorthreshold', 1);
hold on;
for i = 1:size(R, 2)
    y0 = find(T==i); 
    x0 = -0.05;
    normR = log10(elf_stats_normalise(R(:, i), E));
    normR_range1 = normR / (range(normR) * 20 + 0.001);
    normR_shifted = normR_range1 - mean(normR_range1);
    x = - normR_shifted + x0;
    y = E/380 + y0;
    plot(x, y, 'k');
end
set(gca, 'ticklabelinterpreter', 'none', 'XTickLabelRotation', 90, 'Ticklength', [0 0])
a = axis(gca);
axis([-0.1 a(2) a(3) a(4)]);

%% Phylogenetic analysis
Z2 = seqlinkage(D, 'weighted', datasets);
view(Z2)

% [i,j] = cluster(Z2,[],'criterion','gain','maxclust',10);
% h = plot(Z2);
% set(h.BranchLines(j==2),'Color','b')
% set(h.BranchLines(j==1),'Color','r')


%%TODO: http://se.mathworks.com/help/bioinfo/examples/bootstrapping-phylogenetic-trees.html


function D = nested_areadiff(R1, R2)

    for iii = 1:size(R2, 1)
        D(iii) = elf_stats_areadiff(R1', R2(iii, :)', E, dist_degree);
    end

end


end % main