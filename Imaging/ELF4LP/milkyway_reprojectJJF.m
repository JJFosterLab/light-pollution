function milkyway_reprojectJJF(dataset, rotation, outputfolder)
% MILKYWAY_REPROJECT


%% Set up paths and file names; read info, infosum and para
veps_paths;
para        = veps_paraJJF('', dataset, '*.dng');
para        = veps_para_update(para);                                       % Combine old parameter file with potentially changed information in current veps_para
allfiles    = dir(fullfile(para.paths.datapath, para.paths.scenefolder, '*.mat'));
fnames_im   = {allfiles.name};                                              % collect image names

                    veps_log_printmsg('      Processing %d scenes in environment %s\n', length(fnames_im), dataset);
                     
%% Reproject and save images
for setnr = 1:length(fnames_im)
                    veps_log_printmsg('         Scene %03d\n', setnr);
    ims             = veps_readwrite(para, 'loadfilt_mat', fnames_im{setnr});
    im_filt_reproj  = cell(size(ims));
    for j = 1:length(ims)
        im_filt_reproj{j} = veps_project_reproject2fisheye(ims{j}, [], [], [1000 1000], rotation(setnr)); % im, azi, ele, imsize
    end
    [~, f] = fileparts(fnames_im{setnr});
    save(fullfile(outputfolder, [dataset '_' f '_reproj.mat']), 'im_filt_reproj');
end
    
end