function [intim, dolpim, aopim, iBestImage] = polar_main4_stitch(int, dolp, aop, horlimit, correct)
% POLAR_MAIN4_STITCH stitches together different compass directions into a single image
%
% Usage: [intim, dolpim, aopim, iBestImage] = polar_main4_stitch(int, dolp, aop, horlimit, correct)
%
% Inputs:
%   int             - 1 x 5 cell of MxNx3x5 doubles, intensity map for each compass direction 
%   dolp            - 1 x 5 cell of MxNx5 doubles, degree of linear polarisation map for each compass direction 
%   aop             - 1 x 5 cell of MxNx5 doubles, angle of polarisation map for each compass direction
%   horlimit        - 1 x 1 double, excentricity limit (usually 65 for polar datasets)
%   correct         - 1 x 1 double, triggers relative scaling of images using their overlap:
%                                   0 - do not try to correct the small exposure differences between stitched images
%                                   1 - scale images taking into account the overlap of each image with the UP-image
%                                   2 - scale images taking into account all overlapping areas between images    
%                                   Corrections are only applied to DoLP and Int images, not to AoP
%
% Outputs:
%   intim           - P x Q x 3 double, intensity image covering the whole visual field
%   dolpim          - P x Q double, DoLP image covering the whole visual field
%   aopim           - P x Q double, AoP image covering the whole visual field
%   iBestImage      - P x Q double, index matrix describing which image each point has been sampled from
%
% Loads files: None
% Saves files: None

%% 4. Stitch images into a single one

                    elf_support_logmsg('\b\b\b\b\b\b\b\b\b\b\b\b\b\n');
                    elf_support_logmsg('----- polarELF Step 5: Stitching Intensity, AoP, and DoLP images -----\n');
                    fprintf('     Filter level ');
                    
% find each azimuth and elevation in the source images
[finalw_im, finalh_im, iBestImage, full_w_im, full_h_im] = stitch_5dirs('', [0 90 180 270], -180:0.1:180, 90:-0.1:-90, 0, 'rect', horlimit);

intim = cell(5, 1); aopim = intim; dolpim = intim; % pre-allocate
for filtlevel = 1:length(int)
    
                    fprintf('%d..', filtlevel);
                    
    azires = 180 / (size(int{filtlevel}, 1) - 1); % 0.4 for filt 4;
    eleres = 180 / (size(int{filtlevel}, 2) - 1);

    % First, calculate overall linear indices
    % translate elevations into image pixels
    x = round((finalw_im+90) / azires + 1);
    y = round((-finalh_im+90) / eleres + 1);
    ind3d = elf_project_sub2ind_4D([size(aop{filtlevel}, 1) size(aop{filtlevel}, 2) 1 size(aop{filtlevel}, 3)], x, y, iBestImage); % and this is for the 3D AoP and DoLP images
    ind4d = elf_project_sub2ind_4D(size(int{filtlevel}), x, y, iBestImage); % This is for the 4d INT image
    
            
    if correct
        % translate elevations into image pixels
        x = round((full_w_im+90) / azires + 1);
        y = round((-full_h_im+90) / eleres + 1);

        for compDir = 1:size(x, 3) %
            ind = elf_project_sub2ind([size(dolp{filtlevel}, 1) size(dolp{filtlevel}, 2) 1], x(:, :, compDir), y(:, :, compDir));
            dolp_proj_ind(:, :, 1, compDir) = elf_project_apply(dolp{filtlevel}(:, :, compDir), ind, [1801, 3601, 1]);
            ind = elf_project_sub2ind([size(int{filtlevel}, 1) size(int{filtlevel}, 2) size(int{filtlevel}, 3)], x(:, :, compDir), y(:, :, compDir));
            int_proj_ind(:, :, :, compDir) = elf_project_apply(int{filtlevel}(:, :, :, compDir), ind, [1801, 3601, 3]);
        end
        
%             % Remove zeros (they should have really been removed earlier)
%             dolp_proj_ind(dolp_proj_ind==0) = NaN;
%             sel = int_proj_ind(:, :, 1, :)==0 | int_proj_ind(:, :, 2, :)==0 | int_proj_ind(:, :, 3, :)==0 | isnan(int_proj_ind(:, :, 1, :)) | isnan(int_proj_ind(:, :, 2, :)) | isnan(int_proj_ind(:, :, 3, :));
%             sel = repmat(sel, [1 1 3 1]);
%             int_proj_ind(sel) = NaN; % remove zeros (they should have really been removed earlier)

        if correct==1
            compexp = 1;
        elseif correct==2
            compexp.comp = {[1 1]; [2 1]; [3 1]; [4 1]; [5 1]; [2 3]; [3 4]; [4 5]; [5 2]; [2 5]; [3 2]; [4 3]; [5 4]};
            compexp.mult = {{[1 1]}; {[2 1], [6 3], [10 5]}; {[3 1], [7 4], [11 2]}; {[4 1], [8 5], [12 3]}; {[5 1], [9 2], [13 4]}};
        else
            error('Internal error: Unknown setting of "correct".');
        end

        if filtlevel == 1 % calculate scaling factor on unfiltered images
            [~, scalefac_dolp] = elf_hdr_scaleStack(dolp_proj_ind, zeros(size(dolp_proj_ind)), ones(1, 5), compexp);
        end
        for compDir = 1:size(x, 3)
            dolp{filtlevel}(:, :, compDir) = scalefac_dolp(compDir) * dolp{filtlevel}(:, :, compDir);
        end
        
        if filtlevel == 1 % calculate scaling factor on unfiltered images
            [~, scalefac_int] = elf_hdr_scaleStack(int_proj_ind, zeros(size(int_proj_ind)), ones(3, 5), compexp);
        end
        for compDir = 1:size(x, 3)
            int{filtlevel}(:, :, :, compDir) = scalefac_int(compDir) * int{filtlevel}(:, :, :, compDir);
        end
    end
    
    % Finally, apply indices to corrected or uncorrected matrices
    intim{filtlevel}  = elf_project_apply(int{filtlevel}, ind4d, [1801, 3601, 3]);
    dolpim{filtlevel} = elf_project_apply(dolp{filtlevel}, ind3d, [1801, 3601]);
    aopim{filtlevel}  = elf_project_apply(aop{filtlevel}, ind3d, [1801, 3601]);
end

                    fprintf('done.\n');
