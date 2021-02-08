function [para, status, gui] = elf_startup(cbhandle, rootfolder, verbose, useoldfolder)

%% defaults
if nargin < 4 || isempty(useoldfolder), useoldfolder = true; end
if nargin < 3 || isempty(verbose), verbose = true; end
if nargin < 2 || isempty(rootfolder), rootfolder = ''; end

%% get basic parameters
if useoldfolder
    para = elf_para(rootfolder); % without arguments, just returns basic parameters (call again later with rootfolder or dataset)
else
    para = elf_para(NaN); % without arguments, just returns basic parameters (call again later with rootfolder or dataset)
end

%% collect and check all datasets parameters
[status, para, datasets, exts] = elf_checkdata(para, verbose);

%% build GUI
gui = elf_maingui(status, para, datasets, exts, cbhandle);

%% insert images
for i = 1:size(status, 1)
    para2 = elf_para(para.paths.root, datasets{i});
    if status(i, 3)
        fname   = para2.paths.fname_meanimg_jpg;
        info    = [];
        corr    = false;
    else
        if isnan(exts{i}) % this happens when there are only raw files in the folder
            mask = '*.*';
        else
            mask = ['*' exts{i}];
        end
        allims  = elf_io_dir(fullfile(para2.paths.datapath, mask));
        valid   = ~[allims.isdir]; 
        fname   = fullfile(para2.paths.datapath, allims(find(valid, 1, 'first')).name);
        info    = elf_imfinfo(fname);
        corr    = strcmp(exts{i}, '.dng'); % if these are dng images, perform colour and gamma correction
    end
    I       = elf_imread(fname);
    scale   = para.gui.smallsize/size(I, 2);
    I       = imresize(I, scale);
    elf_plot_image(I, info, gui.p(i).ah, '', corr);
    drawnow;
end
