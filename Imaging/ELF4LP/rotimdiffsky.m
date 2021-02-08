function fh = rotimdiffsky(para, datasets, skylim)
if ~exist('troubleshoot','var')
    troubleshoot = 0;
end
%         clc
%         close all
%         clear all
%         %% do you want to troubleshoot a single image?
%         troubleshoot = 1
%         %% add all the functions needed
%         addpath(fullfile(getenv('HOME'),'Dropbox/Matlab/jochensmolka-elf-461d23388273'));
%         addpath(fullfile(getenv('HOME'),'Dropbox/Matlab/ELF4LP'));
%         addpath(fullfile(getenv('HOME'),'Dropbox/Matlab/ELF4LP/altmany-export_fig-3175417'));
%         elf_paths;
%         %% select folder to process
%         if ~exist('inputfolder','var')
%             %if no inputfolder has been specified, get the user to select one
%             para = elf_para('%');
%         else %if the user has specificed an inputfolder, make that the root directory
%             para = elf_para(inputfolder);
%         end
%         %find datasets in the folder
if(~exist('datasets','var'))
        [~, ~, datasets] = elf_checkdata(para);
end
if(~exist('skylim','var'))
        skylim = 55;65;%I have decided that 55° is the magic boundary between buildings and sky
%in natural environments, 30° would probably be a better value
end
lowlim = 1;%20;
%%
if troubleshoot
    j = 1
    i = 1
%%  Plot and save relative to an appropriate maximum
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

     irrad = night_irradiance(strcat(allfiles(j).folder, '/',allfiles(j).name)); % in photons / s / m2 / nm
    %     disp('Irradiance (photons/s/cm2/nm)');
    %     disp(irrad/10000)
        [mini, indi] = min(abs(log10(mean(irrad(4,:))) - log10(nmoptions)));
        if mini >1
            disp 'uh oh, no good display value found';
        end
        nmax = nmoptions(indi);
        [~, f] = fileparts(allfiles(j).name);
    %now plot them
        %I am loading reprojected images
        % they get rescaled at the image-difference step
                temp    = load(fullfile(allfiles(j).folder,[f,'_reproj']), 'im_filt_reproj');
                ims     = temp.im_filt_reproj;
                clear temp;
            % Calculate elevation for each pixel (this should work for all images now, since they should all be square and the same size)
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
                        im = fliplr(im);
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
%% Calculate image differences
    %         wh_im = plotim{1, 4};% 1st slot is filter level: [2,4,8,16], 2nd slot is R, G, B, W?
        wh_im = plotim{2, 4};% select slot 2, which contains the 4° filtered image
        totimf = abs(wh_im(:,:,1)+wh_im(:,:,2)+wh_im(:,:,3))/10000; %now in photons/s/cm2/nm
% AT THIS POINT YOU MIGHT CONSIDER THE SKY--BUILDING BOUNDARY
        normval = quantile(totimf(:),0.95);
        imshow(totimf./normval);
        smallim = totimf;
        smallim(ele < 90-skylim) = 0;
        imshow(1-smallim./normval);
        div = 10; % fraction of a degree to calculate across
        rangles = (-(180*div):(180*div)) / div;%from -180° to 180° in 1/div°
%         rmsimdiff = zeros(length(rangles), 1); %ignore uncorrected differences
        rmsimdiff_corr = zeros(length(rangles), 1);
        rmsimdiff_0_65 = zeros(length(rangles), 1);
        rmsimdiff_65_90 = zeros(length(rangles), 1);
%         rmsimdiff_60_90 = zeros(length(rangles), 1);
        
%       strange things are happening at the edge of the cosine correction
%       lets predetermine how it looks
    cc_good = 1./cosd(ele);%1/cosine of elevation in degrees
    cc_good(ele > 90-lowlim);%89) = 0;%can we fix artefacts if we don't apply to the edge pixels?
%     surf(log10(cc_good), 'EdgeColor','none');view(0,90);colorbar
%     title('log_{10}( cos(90°-elevation)^{-1} )');
for theta = 1:length(rangles)
    temp = imrotate((totimf), rangles(theta), 'bicubic', 'crop');
    temp2 = imabsdiff((totimf), temp);
%     rmsimdiff(theta) = sqrt(mean(temp2(:).^2));
    %weight differences by the solid angles they represent
    temp3 = temp2*cc_good ;
    rmsimdiff_corr(theta) = rms(temp3(:));
    temp3(ele >(90-skylim)) = 0;%remove elevations below (90°-25°) = 65° %TODO change elevation to real world elevation!
    rmsimdiff_65_90(theta) = rms(temp3(:));
    %weight differences by the solid angles they represent
    temp3 = temp2*cc_good  ;
    temp3(ele <(90-skylim)) = 0;%remove elevations above (90°-25°) = 65°
    rmsimdiff_0_65(theta) = rms(temp3(:));
end
%% plot curves 
% plot(rangles, log10(rmsimdiff),'LineWidth',2)
% plot on a log10 scale
% plot(rangles, log10(rmsimdiff_60_90),'LineWidth',2)
plot(rangles, log10(rmsimdiff_65_90),'LineWidth',2);hold on
plot(rangles, log10(rmsimdiff_0_65),'LineWidth',2)
plot(rangles, log10(rmsimdiff_corr),'LineWidth',2)
set(gca,'XTick',(-4:4) * 45);%, 'XTickLabels', lbnm);
ylim([0,15]+1);%ylim([0,12]+1)
xlabel('Rotation (°)')
ylabel('log10 RMS Image Difference (photons/s/cm2/nm)')
legend(...'Uncorrected', 'Cosine Corrected',
        '30°–90°','1°–30°','1°–90°' );%, '60°–90°')
title([datasets{i} ' ' 'scene ' sprintf('%03d', j)])
% save as images
svpath = [ datasets{i} '_' 'scene' sprintf('%03d', j) '_rotdiff_' ];
pdfsave(gcf, fullfile(pdffolder,[ svpath '_log10_rmsimdiff.pdf']));
export_fig(fullfile(pdffolder,[ svpath '_log10_rmsimdiff.png']), '-native');

hold off
% plot(rangles, (rmsimdiff)/max(rmsimdiff),'LineWidth',2)
% plot on a normalised scale
% plot(rangles, (rmsimdiff_60_90)/max(rmsimdiff_60_90),'LineWidth',2)
plot(rangles, (rmsimdiff_65_90)/max(rmsimdiff_65_90),'LineWidth',2, 'Color','#0072BD');hold on
plot(rangles, (rmsimdiff_0_65)/max(rmsimdiff_0_65),'LineWidth',2, 'Color','#D95319')
plot(rangles, (rmsimdiff_corr)/max(rmsimdiff_corr),'LineWidth',2, 'Color','#EDB120')
set(gca,'XTick',(-4:4) * 45);%, 'XTickLabels', lbnm);
ylim([0,1.1])
xlabel('Rotation (°)')
ylabel('RMS Image Difference (normalised)')
legend(...'Uncorrected', 'Cosine Corrected',
        '60°–90°','1°–60°','1°–90°',...%, '60°–90°')
        'location', 'SouthWest')
title([datasets{i} ' ' 'scene ' sprintf('%03d', j)])
% save as images
pdfsave(gcf, fullfile(pdffolder,[ svpath '_rmsimdiff_normalised.pdf']));
export_fig(fullfile(pdffolder,[ svpath '_rmsimdiff_normalised.png']), '-native');

   hold off
    %find the half max for each 
    hm_corr = hmax(rangles, rmsimdiff_corr);%TODO more efficient optimisation for hmax
    hm_0_60 = hmax(rangles, rmsimdiff_0_65);
    hm_60_90 = hmax(rangles, rmsimdiff_65_90);
    % plot(rangles, (rmsimdiff)/max(rmsimdiff),'LineWidth',2)
    % plot on a normalised scale with half max

    % plot(hm_60_90, repmat(0.5, 1, 2), '-.')
    plt = [];
    
    plt(1) = plot(rangles, (rmsimdiff_65_90)/max(rmsimdiff_65_90),'LineWidth',5, 'Color','#0072BD');hold on
%     plot(hm_30_90, repmat(0.5, 1, 2), ':','LineWidth',2,'Color','#0072BD')
    ps = polyshape([hm_60_90, flip(hm_60_90)], [0,0,0.5,0.5]);
    pg = plot(ps); pg.EdgeColor = 'none'; pg.FaceColor = '#0072BD'; pg.FaceAlpha = 0.2;
   
    plt(2) = plot(rangles, (rmsimdiff_0_65)/max(rmsimdiff_0_65),'LineWidth',5, 'Color','#D95319');
%     plot(hm_0_30, repmat(0.5, 1, 2), '--','LineWidth',2,'Color','#D95319')
    ps = polyshape([hm_0_60, flip(hm_0_60)], [0,0,0.5,0.5]);
    pg = plot(ps); pg.EdgeColor = 'none'; pg.FaceColor = '#D95319'; pg.FaceAlpha = 0.2;
    
    plt(3) = plot(rangles, (rmsimdiff_corr)/max(rmsimdiff_corr),'LineWidth',5, 'Color','#77AC30');
%     plot(hm_corr, repmat(0.5, 1, 2),'-','LineWidth',2,'Color','#77AC30')
    ps = polyshape([hm_corr, flip(hm_corr)], [0,0,0.5,0.5]);
    pg = plot(ps); pg.EdgeColor = 'none'; pg.FaceColor = '#77AC30'; pg.FaceAlpha = 0.2;    
    
    ylim([0,1.1])
    xlabel('Rotation (°)')
    ylabel('RMS Image Difference (normalised)')
    lgd = legend(plt,...'Uncorrected', 'Cosine Corrected',
            strcat(num2str(skylim),'°–90°'),...
            strcat('1°–',num2str(skylim),'°'),...
            '1°–90°',...
            'location', 'SouthWest');%..., 'LineSpec', '#0072BD','#D95319','#77AC30'),...
    %         'LineSpec', '-' )
    text(max([hm_60_90(2), hm_0_60(2), hm_corr(2)])+15, 0.5+0.5/10,...
        strcat(num2str(round(diff(hm_60_90))), '°'),...
        'fontsize', 20,'Color','#0072BD')
    text(max([hm_60_90(2), hm_0_60(2), hm_corr(2)])+15, 0.5+0/10,...
        strcat(num2str(round(diff(hm_0_60))), '°'),...
        'fontsize', 20,'Color','#D95319')
    text(max([hm_60_90(2), hm_0_60(2), hm_corr(2)])+15, 0.5-0.5/10,...
        strcat(num2str(round(diff(hm_corr))), '°'),...
        'fontsize', 20,'Color','#77AC30')

    title([datasets{i} ' ' 'scene ' sprintf('%03d', j)])
    % save as images
    hold off
    pdfsave(gcf, fullfile(pdffolder,[ svpath '_rmsimdiff_normalised_halfmax.pdf']));
    export_fig(fullfile(pdffolder,[ svpath '_rmsimdiff_normalised_halfmax.png']), '-native');


else
    for i = 1:length(datasets)%[3 6]%
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

            irrad = night_irradiance(strcat(allfiles(j).folder, '/',allfiles(j).name)); % in photons / s / m2 / nm
            [mini, indi] = min(abs(log10(mean(irrad(4,:))) - log10(nmoptions)));
            if mini >1
                disp 'uh oh, no good display value found';
            end
            nmax = nmoptions(indi);
            [~, f] = fileparts(allfiles(j).name);

            temp    = load(fullfile(allfiles(j).folder,[f,'_reproj']), 'im_filt_reproj');
            ims     = temp.im_filt_reproj;
            clear temp;
 
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
                        im = fliplr(im);
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
    %% Calculate image differences
        %         wh_im = plotim{1, 4};% 1st slot is filter level: [2,4,8,16], 2nd slot is R, G, B, W?
                        wh_im = plotim{2, 4};% select slot 2, which contains the 4° filtered image
                        totimf = abs(wh_im(:,:,1)+wh_im(:,:,2)+wh_im(:,:,3))/10000; %now in photons/cm2/s/nm
            div = 10; % fraction of a degree to calculate across
            rangles = (-(180*div):(180*div)) / div;
%             rmsimdiff = zeros(length(rangles), 1);
            rmsimdiff_corr = zeros(length(rangles), 1);
            rmsimdiff_0_65 = zeros(length(rangles), 1);
            rmsimdiff_65_90 = zeros(length(rangles), 1);
    %       strange things are happening at the edge of the cosine correction
    %       lets predetermine how it looks
        cc_good = 1./cosd(ele);%1/cosine of elevation in degrees
        cc_good(ele >90-lowlim) = 0;%cc_good(ele >89) = 0;%can we fix artefacts if we don't apply to the edge pixels?
        %N.B. this is also where the numeric values of the correction begin to break down
    for theta = 1:length(rangles)
        temp = imrotate((totimf), rangles(theta), 'bicubic', 'crop');
        temp2 = imabsdiff((totimf), temp);
    %     rmsimdiff(theta) = sqrt(mean(temp2(:).^2));
        %weight differences by the solid angles they represent
        temp3 = temp2*cc_good ;
        rmsimdiff_corr(theta) = rms(temp3(:));
        temp3(ele >(90-skylim)) = 0;%remove elevations below (90°-30°) = 60° %TODO change elevation to real world elevation!
        rmsimdiff_65_90(theta) = rms(temp3(:));
        %weight differences by the solid angles they represent
        temp3 = temp2*cc_good  ;
        temp3(ele <(90-skylim)) = 0;%remove elevations above (90°-30°) = 60°
        rmsimdiff_0_65(theta) = rms(temp3(:));
    end
    %%plot curves and find half max
    svpath = [ datasets{i} '_' 'scene' sprintf('%03d', j) '_rotdiff_' ]; 
    % plot(rangles, log10(rmsimdiff),'LineWidth',2)
    % plot on a log10 scale
    % plot(rangles, log10(rmsimdiff_60_90),'LineWidth',2)
    fh = figure();
    
    plot(rangles, log10(rmsimdiff_65_90),'LineWidth',2);hold on
    plot(rangles, log10(rmsimdiff_0_65),'LineWidth',2)
    plot(rangles, log10(rmsimdiff_corr),'LineWidth',2)
    set(gca,'XTick',(-4:4) * 45);%, 'XTickLabels', lbnm);
    ylim([0,12]+1)
    xlabel('Rotation (°)')
    ylabel('log10 RMS Image Difference (photons/s/cm2/nm)')
    legend(...'Uncorrected', 'Cosine Corrected',
            strcat(num2str(skylim),'°–90°'),...
            strcat('1°–',num2str(skylim),'°'),...
            '1°–90°');%, '60°–90°')
    title([datasets{i} ' ' 'scene ' sprintf('%03d', j)])
    % save as images
    pdfsave(gcf, fullfile(pdffolder,[ svpath '_log10_rmsimdiff.pdf']));
    export_fig(fullfile(pdffolder,[ svpath '_log10_rmsimdiff.png']), '-native');
    
    hold off
    % plot(rangles, (rmsimdiff)/max(rmsimdiff),'LineWidth',2)
    % plot on a normalised scale
    % plot(rangles, (rmsimdiff_60_90)/max(rmsimdiff_60_90),'LineWidth',2)
    plot(rangles, (rmsimdiff_65_90)/max(rmsimdiff_65_90),'LineWidth',2, 'Color','#0072BD');hold on
    plot(rangles, (rmsimdiff_0_65)/max(rmsimdiff_0_65),'LineWidth',2, 'Color','#D95319')
    plot(rangles, (rmsimdiff_corr)/max(rmsimdiff_corr),'LineWidth',2, 'Color','#EDB120')
    set(gca,'XTick',(-4:4) * 45);%, 'XTickLabels', lbnm);
    ylim([0,1.1])
    xlabel('Rotation (°)')
    ylabel('RMS Image Difference (normalised)')
    legend(...'Uncorrected', 'Cosine Corrected',
            strcat(num2str(skylim),'°–90°'),...
            strcat('1°–',num2str(skylim),'°'),...
            '1°–90°',...%, '60°–90°')
            'location', 'SouthWest')
    title([datasets{i} ' ' 'scene ' sprintf('%03d', j)])
    % save as images
    pdfsave(gcf, fullfile(pdffolder,[ svpath '_rmsimdiff_normalised.pdf']));
    export_fig(fullfile(pdffolder,[ svpath '_rmsimdiff_normalised.png']), '-native');

    hold off
    %find the half max for each 
    hm_corr = hmax(rangles, rmsimdiff_corr);
    hm_0_60 = hmax(rangles, rmsimdiff_0_65);
    hm_60_90 = hmax(rangles, rmsimdiff_65_90);
    hm_60_90 = hmax(rangles, rmsimdiff_65_90);
    % plot(rangles, (rmsimdiff)/max(rmsimdiff),'LineWidth',2)
    % plot on a normalised scale with half max

    plt = [];
    
    plt(1) = plot(rangles, (rmsimdiff_65_90)/max(rmsimdiff_65_90),'LineWidth',5, 'Color','#0072BD');hold on
%     plot(hm_30_90, repmat(0.5, 1, 2), ':','LineWidth',2,'Color','#0072BD')
    ps = polyshape([hm_60_90, flip(hm_60_90)], [0,0,0.5,0.5]);
    pg = plot(ps); pg.EdgeColor = 'none'; pg.FaceColor = '#0072BD'; pg.FaceAlpha = 0.2;
   
    plt(2) = plot(rangles, (rmsimdiff_0_65)/max(rmsimdiff_0_65),'LineWidth',5, 'Color','#D95319');
%     plot(hm_0_30, repmat(0.5, 1, 2), '--','LineWidth',2,'Color','#D95319')
    ps = polyshape([hm_0_60, flip(hm_0_60)], [0,0,0.5,0.5]);
    pg = plot(ps); pg.EdgeColor = 'none'; pg.FaceColor = '#D95319'; pg.FaceAlpha = 0.2;
    
    plt(3) = plot(rangles, (rmsimdiff_corr)/max(rmsimdiff_corr),'LineWidth',5, 'Color','#77AC30');
%     plot(hm_corr, repmat(0.5, 1, 2),'-','LineWidth',2,'Color','#77AC30')
    ps = polyshape([hm_corr, flip(hm_corr)], [0,0,0.5,0.5]);
    pg = plot(ps); pg.EdgeColor = 'none'; pg.FaceColor = '#77AC30'; pg.FaceAlpha = 0.2;    
    
    ylim([0,1.1])
    xlabel('Rotation (°)')
    ylabel('RMS Image Difference (normalised)')
    lgd = legend(plt,...'Uncorrected', 'Cosine Corrected',
            strcat(num2str(skylim),'°–90°'),...
            strcat('1°–',num2str(skylim),'°'),...
            '1°–90°',...%, '60°–90°')
            'location', 'SouthWest');%..., 'LineSpec', '#0072BD','#D95319','#77AC30'),...
    %         'LineSpec', '-' )
    text(max([hm_60_90(2), hm_0_60(2), hm_corr(2)])+15, 0.5+0.5/10,...
        strcat(num2str(round(diff(hm_60_90))), '°'),...
        'fontsize', 20,'Color','#0072BD')
    text(max([hm_60_90(2), hm_0_60(2), hm_corr(2)])+15, 0.5+0/10,...
        strcat(num2str(round(diff(hm_0_60))), '°'),...
        'fontsize', 20,'Color','#D95319')
    text(max([hm_60_90(2), hm_0_60(2), hm_corr(2)])+15, 0.5-0.5/10,...
        strcat(num2str(round(diff(hm_corr))), '°'),...
        'fontsize', 20,'Color','#77AC30')

    title([datasets{i} ' ' 'scene ' sprintf('%03d', j)])
    % save as images
    hold off
    pdfsave(gcf, fullfile(pdffolder,[ svpath '_rmsimdiff_normalised_halfmax.pdf']));
    export_fig(fullfile(pdffolder,[ svpath '_rmsimdiff_normalised_halfmax.png']), '-native');

        end
    end
end