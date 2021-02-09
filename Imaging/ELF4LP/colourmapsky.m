function fh = colourmapsky(para, datasets)
%% N.B. These sections have been removed, skypipe will do this from now on
%         clc
%         close all
%         clear all
%         % add all the functions needed
%         addpath(fullfile(getenv('HOME'),'Dropbox/Matlab/jochensmolka-elf-461d23388273'));
%         addpath(fullfile(getenv('HOME'),'Dropbox/Matlab/ELF4LP'));
%         addpath(fullfile(getenv('HOME'),'Dropbox/Matlab/ELF4LP/altmany-export_fig-3175417'));
%         elf_paths;
%         %% select folder to process
%         %EVERY folder should contain a file called "brackets.info" with the number of
%         %the first and the last image in the bracket (I think in chronological order)
%         % so if the first backet contained 5 images, and the second one was three
%         % brackets.info would have the contents
%         % 1 5
%         % 6 9
%         if ~exist('inputfolder','var')
%             %if no inputfolder has been specified, get the user to select one
%             para = elf_para('%');
%         else %if the user has specificed an inputfolder, make that the root directory
%             para = elf_para(inputfolder);
%         end
%         % this step never did anything, let the user select the folder
%             % if ~exist('skydir','var')
%             %     skydir = fullfile(getenv('HOME'),'/Documents/Starry Project/Light Pollution Imaging/');
%             % end
%             % if ~exist('filename','var')
%             %     inputfolder = fullfile(skydir,'LPexper-20181127-LowMWay-Thornwood');
%             % end
%         %don't be surprised if you get asked for the top folder
%         % para = elf_para('20181127StarryThornwood_1928_1940');
%         %find datasets in the folder
if(~exist('datasets','var'))
        [~, ~, datasets] = elf_checkdata(para);
end

%%  Plot and save relative to an appropriate maximum
for i = 1:length(datasets)%[3 6]%
%     i = 1;%for now just deal with one file
    % allfiles    = dir(fullfile(reprojfolder, '*.mat'));
    allfiles    = dir(fullfile(para.paths.root,datasets{i}, para.paths.filtfolder, '*filt.mat'));
    pdffolder   = fullfile(para.paths.root,datasets{i});%'E:\James data\04Jochen_pdfoutput';
    para.paths.datapath = fullfile(para.paths.root,datasets{i});
    rotation = zeros( 1,length(allfiles));%no rotation right now
    % nmax is the maximum value in photons / s / m2 / nm
    % as used by Foster et al., 2018a: https://www.doi.org/10.1098/rspb.2017.2322
    % 3*10^10 is typical of clear moonless nights
    % 1*10^11 is for moderate light pollution or aurora borealis (norrsken)
    % 1*10^12 is for strong light pollution (e.g. Johannesburg)
    % 3*10^12 is good for full moon
    % 3*10^14 is good for security lighting (e.g. Mushrooms)
    % nmax =   3*10^10; 3*10^12; 1*10^11; 1*10^12;  %VALUE TO MAXIMISE TO
    nmoptions = [3*10^10; 3*10^12; 1*10^11; 1*10^12; 3*10^14];

    for j = 1:length(allfiles)
%     j = 1;%for now just deal with one file
        irrad = night_irradiance(strcat(allfiles(j).folder, '/',allfiles(j).name)); % in photons / s / m2 / nm
    %     disp('Irradiance (photons/s/cm2/nm)');
    %     disp(irrad/10000)
        [mini, indi] = min(abs(log10(mean(irrad(4,:))) - log10(nmoptions)));
        if mini >1
            disp 'uh oh, no good display value found';
        end
        nmax = nmoptions(indi);
    %now plot them
        [~, f] = fileparts(allfiles(j).name);
        temp    = load(fullfile(allfiles(j).folder,[f,'_reproj']), 'im_filt_reproj');
        ims     = temp.im_filt_reproj;
        clear temp;
        % Calculate elevation for each pixel (this should work for all images now, since they should all be square and the same size)
        % For non-square images, BUGFIX THIS!
        mid    = [1+(size(ims{1}, 1)-1)/2; 1+(size(ims{1}, 2)-1)/2];        % centre of image
        r_full = 8 * size(ims{1}, 2) / 24;                                 % theoretical value for 24mm high chip that is fully covered by fisheye circular image
        [x, y] = meshgrid(1:size(ims{1}, 1), 1:size(ims{1}, 2));
        r      = sqrt((x-mid(1)).^2 + (y-mid(2)).^2);
        ele    = asind(r / 2 / r_full) * 2;%r/r_full * 90;
        ele_s  = ele; ele_s(y<mid(2)) = -ele_s(y<mid(2));
            % 1. extract images
        rotit = rotation(j);
        for ch = 1:4
            for sc = 1:4
                im = ims{sc};
                im(isnan(im)) = 0;
                if rotit
                    im = rot90(permute(im, [2 1 3]), 2);        % rotate image if necessary to place Milky Way horizontally %SHOULD BE POSSIBLE TO KEEP MW VERTICAL
                else
                    %im = fliplr(im);  %TODO CHECK IF THIS IS NECESSARY HERE, SEEMS UNLIKELY
                end
                if ch < 4
                    im(:, :, [1:ch-1 ch+1:3]) = 0;
                    sumim{sc, ch}  = sum(im, 3);
                else
                    sumim{sc, ch}  = mean(im, 3);
                end   % if colour channel, set other channels to 0 (works for plotting and contrast calculation)
                plotim{sc, ch} = im;                            % save back into ims structure, as this will later be plotted
            end %sc = 1:4
        end %ch = 1:4
        wh_im = plotim{1, 4};% 1st slot is filter level: [2,4,8,16], 2nd slot is R, G, B, W?
%         wh_im = plotim{2, 4};% select slot 2, which contains the 4° filtered image
       fh   = figure();
       imshow((wh_im)./(nmax));%I don't think it should be flipped %imshow(flipud(wh_im)./(nmax));
       lm = [6,11]; %lm = [6,13]
       totim = abs(wh_im(:,:,1)+wh_im(:,:,2)+wh_im(:,:,3))/10000; %now in photons/cm2/s/nm
       ltotim = log10(totim);
       addpath(fullfile(getenv('HOME'),'Dropbox/Matlab/kakearney-cptcmap-pkg-845bf83/cptcmap'));
       addpath(fullfile(getenv('HOME'),'Dropbox/Matlab/kakearney-cptcmap-pkg-845bf83/parsepv'));
%        surf(ltotim, 'EdgeColor','none'); view(90, 90);
       shown = imshow((ltotim))%imshow(flipud(ltotim));%BUGFIX, CHECK THIS. SEEMS TO BE PLOTTING AT THE WRONG ANGLE
       cptcmap('GMT_wysiwygcont', 'mapping', 'scaled');
       cb = colorbar();
       caxis(lm);
       svpath = [ datasets{i} '_' 'scene' sprintf('%03d', j) '_colormapped_' sprintf('%g',[10^6, 10^13])];
       pdfsave(gcf, fullfile(pdffolder,[ svpath '_wysiwygcont.pdf']));
       export_fig(fullfile(pdffolder,[ svpath '_wysiwygcont.png']), '-native');
%        cb = cptcbar(gca, 'GMT_globe', 'eastoutside', false);
%        set(cb.ax, 'fontsize', 7);
    end
end
% % imshow(totim/(nmax/10000))
% wh_im = plotim{2, 4};% select slot 2, which contains the 4° filtered image
%  totimf = abs(wh_im(:,:,1)+wh_im(:,:,2)+wh_im(:,:,3))/10000; %now in photons/cm2/s/nm
% rrrr = imrotate(totimf, 1, 'nearest', 'crop');
% % imshow(rrrr/(nmax/10000))
% sz= size(rrrr);
% % rrrr1 = rrrr(   ((sz(1)-1000)/4 -4):(1000-((sz(1)-1000)/4) +4),...
% %     ((sz(2)-1000)/4 -4):(1000-((sz(2)-1000)/4) +4)    );
% % imshow(rrrr1/(nmax/10000))
% diffim = imabsdiff(totimf, rrrr);
% imshow(diffim./max(diffim))
% cptcmap('GMT_wysiwygcont', 'mapping', 'scaled');
% rangles = -180:180;
% rimdiff = zeros(length(rangles));
% for theta = 1:length(rangles)
%     temp = imrotate(totimf, rangles(theta), 'nearest', 'crop');
%     rimdiff(theta) = sum(sum(imabsdiff(totimf, temp)));
% end
% plot(rangles, log10(rimdiff))
% set(gca,'XTick',(-4:4) * 45);%, 'XTickLabels', lbnm);
