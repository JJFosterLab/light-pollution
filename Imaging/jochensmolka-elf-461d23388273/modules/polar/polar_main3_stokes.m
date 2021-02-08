function [int, dolp, aop] = polar_main3_stokes(datasets, channel, pol_angle_flip)
% POLAR_MAIN3_STOKES calculatea Stokes' vectors and Int/AoP/DoLP matrices
%
% Usage: [int, dolp, aop] = polar_main3_stokes(datasets, channel, pol_angle_flip)
%
% Inputs:
%   datasets        - 1 x M cell of str, names of the datasets, each must contain 5 scenes
%   channel         - 1 x 1 double, which channel to use for DoLP and AoP
%   pol_angle_flip  - 1 x M bool, whether to flip Stokes vectors in this dataset (this is not used anymore, images are flipped during stitching)
%
% Outputs:
%   int             - 1 x 5 cell of MxNx3x5 doubles, intensity map for each compass direction 
%   dolp            - 1 x 5 cell of MxNx5 doubles, degree of linear polarisation map for each compass direction 
%   aop             - 1 x 5 cell of MxNx5 doubles, angle of polarisation map for each compass direction 
%
% Loads files: HDR scenes as mat in scene folder, filtered images as mat in filt folder
% Saves files: None by default (possible to save Stokes vectors as mat)

%% check inputs
if nargin < 3 || isempty(pol_angle_flip), pol_angle_flip = false(size(datasets)); end
if nargin < 2 || isempty(channel), channel = 3; warning('No channel specified for DoLP and AoP. Using green channel.'); end 
if nargin < 1 || isempty(datasets), error('You have to provide a valid dataset name'); end 

                    elf_support_logmsg('\b\b\b\b\b\b\b\b\b\b\b\b\b\n');
                    elf_support_logmsg('----- polarELF Step 3: Calculating Stokes vectors and Int/DoLP/AoP -----\n')
                    fprintf('     Compass direction ');
                    
%%
int  = cell(1, 5); aop = int; dolp = int; % pre-allocate

for compDir = 1:length(datasets) % for each compass direction
        
                    fprintf('%d..', compDir);
                    
    dataset         = datasets{compDir};
    para            = elf_para('', dataset);
    infosum         = elf_readwrite(para, 'loadinfosum');     % loads the old infosum file (which contains projection information)

    % pre-allocate
    im_raw          = zeros(length(infosum.proj_azi), length(infosum.proj_ele), infosum.SamplesPerPixel, 5);
    im_filt         = cell(4, 5);
    imStokes_filt   = cell(5, 1);

    % load raw and filtered images
    for setnr = 1:5
        im_raw(:, :, :, setnr)  = elf_readwrite(para, 'loadHDR_mat', sprintf('scene%03d', setnr));
        im_filt(:, setnr)       = elf_readwrite(para, 'loadfilt_mat', sprintf('scene%03d', setnr));
    end

    % calculate Stokes parameters for unfiltered images
    imStokes_filt{1} = polar_stokes(im_raw, pol_angle_flip);

    % calculate Stokes parameters for filtered images
    for filterWidthNo = 1:size(im_filt, 1)
        thisFilteredSet = cat(4, im_filt{filterWidthNo, :});
        imStokes_filt{filterWidthNo+1} = polar_stokes(thisFilteredSet, pol_angle_flip);
    end

    % % save stokes matrices (500 MB per compass direction)
    % elf_readwrite(para, 'savestokes_mat', sprintf('exposure%d', 1), imStokes_filt);
        
    % 4. Calculate Intensity, AoP, and DoLP
    for filtLevel = 1:5
        int{filtLevel}(:, :, :, compDir)    = imStokes_filt{filtLevel}(:, :, :, 1);
        aop{filtLevel}(:, :, compDir)       = mod(0.5*atan2d(imStokes_filt{filtLevel}(:, :, channel, 3), imStokes_filt{filtLevel}(:, :, channel, 2)), 180); % checked, seems fine
        dolp{filtLevel}(:, :, compDir)      = sqrt(imStokes_filt{filtLevel}(:, :, channel, 3).^2 + imStokes_filt{filtLevel}(:, :, channel, 2).^2);
    end
end

                    fprintf('done.\n');

