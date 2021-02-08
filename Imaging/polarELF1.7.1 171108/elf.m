function elf(useoldfolder)
% ELF is the super-function for the ELF program
% This function creates a GUI allowing the user to browse folders and see the processed and potential datasets in this folder
% Practically all other ELF functions can then be called through buttons next to the images.
%
% Use elf(0) to ignore the saved folder and choose a new one.

%% parameters
if nargin<1
    useoldfolder = true;
end
verbose = true;

%% set paths and parameters
elf_paths;

%% read parameter file, and build GUI
[para, status, gui] = elf_startup(@maincb, '', verbose, useoldfolder);

%% nested function
    function maincb(src, ~)
        if strcmp(get(src, 'tag'), 'maingui_folderbrowse')
            % new folder has been selected, confirm in a gui and, if necessary, restart GUI
            newfolder = uigetdir(para.paths.root, 'Select a new folder');
            if ~all(newfolder == 0) && exist(newfolder, 'file')
                [para, status, gui] = elf_startup(@maincb, newfolder, verbose);
            end
        elseif strcmp(get(src, 'tag'), 'file_reload')
            [para, status, gui] = elf_startup(@maincb, '', verbose);
        else % any other button or key callback
            [status, gui] = elf_callbacks_maingui(src, status, gui, para);
        end
    end

end % main