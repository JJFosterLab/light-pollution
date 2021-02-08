function irrad = night_irradiance(filename)

temp    = load(filename, 'im_filt_HDR');
ims     = temp.im_filt_HDR;
clear temp;



for sc = 1:4
    im = ims{sc};
        
    azi = linspace(-90, 90, size(im, 2));
    ele = linspace(-90, 90, size(im, 1));
    elediff = median(diff(azi));
    azidiff = median(diff(azi));
    solidangle = deg2rad(elediff) * deg2rad(azidiff) * cosd(ele);
    solidangle = repmat(solidangle', 1, size(im, 2));
    
    AZ = repmat(azi, size(im, 1), 1);
    EL = repmat(ele', 1, size(im, 2));
    coscorr = cosd(elf_support_sphdist(AZ, EL, 0, 0));
    
    temp = solidangle .* coscorr .* im(:, :, 1);
    irrad(sc, 1) = sum(temp(:));
    temp = solidangle .* coscorr .* im(:, :, 2);
    irrad(sc, 2) = sum(temp(:));
    temp = solidangle .* coscorr .* im(:, :, 3);
    irrad(sc, 3) = sum(temp(:));
    
end