function im = elf_io_removehorizon(im, horlimit)
% Remove (set to NaN all pixels below a a certain elevation in an image, assuming D810 optics

if ~isnan(horlimit) 
    mid         = [1+(size(im, 1)-1)/2; 1+(size(im, 2)-1)/2];   % centre of image
    shortSide   = min([size(im, 1) size(im, 2)]);
    r_full      = 8 * shortSide / 24;                           % theoretical value for 24mm high chip that is fully covered by fisheye circular image

    [y, x]      = meshgrid(1:size(im, 2), 1:size(im, 1));
    r           = sqrt((x-mid(1)).^2 + (y-mid(2)).^2);
    ele         = asind(r / 2 / r_full) * 2;%r/r_full * 90;

    sel         = ele>horlimit;
    i1      = im(:, :, 1); i1(sel) = NaN;
    i2      = im(:, :, 2); i2(sel) = NaN;
    i3      = im(:, :, 3); i3(sel) = NaN;
    im      = cat(3, i1, i2, i3);
end
