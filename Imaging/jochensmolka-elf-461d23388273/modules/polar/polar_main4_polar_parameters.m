function [int, aop, dolp] = polar_main4_polar_parameters(datasets, channel, off_axis_correct)


%% 4. Calculate Intensity, AoP, and DoLP
int  = cell(length(datasets), 5);
aop  = cell(length(datasets), 5);
dolp = cell(length(datasets), 5);

                    elf_support_logmsg('\b\b\b\b\b\b\b\b\b\b\b\b\b\n');
                    elf_support_logmsg('----- polarELF Step 4: Calculating Intensity, AoP, and DoLP in colour channel %d -----\n', channel);
                    fprintf('     Compass direction ');
          
aziele_res = [.1 .2 .4 .8 1.6];
for compDir = 1:length(datasets) % for each compass direction
        
                    fprintf('%d..', compDir);
                    
    % load stokes vectors
    dataset         = datasets{compDir};
    para            = elf_para('', dataset);
    imStokes_filt   = elf_readwrite(para, 'loadstokes_mat', sprintf('exposure%d', 1));
    
    for filtLevel = 1:5
        int{compDir, filtLevel}  = imStokes_filt{filtLevel}(:, :, :, 1);
        aop{compDir, filtLevel}  = mod(0.5*atan2d(imStokes_filt{filtLevel}(:, :, channel, 3), imStokes_filt{filtLevel}(:, :, channel, 2)), 180); % checked, seems fine
        dolp{compDir, filtLevel} = sqrt(imStokes_filt{filtLevel}(:, :, channel, 3).^2 + imStokes_filt{filtLevel}(:, :, channel, 2).^2);

        if off_axis_correct
            %%% HACK Off-Axis-Correction
            [azi_grid, ele_grid]    = meshgrid(-90:aziele_res(filtLevel):90, -90:aziele_res(filtLevel):90); % grid of desired spherical angles

            if compDir==3 || compDir==5
                xnew = cosd(aop{compDir, filtLevel})./cosd(ele_grid);
                ynew = sind(aop{compDir, filtLevel})./cosd(azi_grid);
                aop{compDir, filtLevel} = atan2d(ynew, xnew);
            else
                xnew = cosd(aop{compDir, filtLevel})./cosd(azi_grid);
                ynew = sind(aop{compDir, filtLevel})./cosd(ele_grid);
                aop{compDir, filtLevel} = atan2d(ynew, xnew);
            end
            %%% HACK END
        end
    end
end

                    fprintf('done.\n');