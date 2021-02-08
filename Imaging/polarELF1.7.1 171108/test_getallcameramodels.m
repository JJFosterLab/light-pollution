out = cell(695, 2);
for i = 1:length(datasets)
    out{i, 1}           = datasets{i};
    try 
        para            = elf_para('', datasets{i}, '*.dng');
        info            = elf_info_collect(para.paths.datapath, '*.dng');   % this contains EXIF information and filenames, verbose==1 means there will be output during system check
        out{i, 2}       = info(1).Model;
        i
    end
end