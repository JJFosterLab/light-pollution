%%
elf_paths
clear
para             = elf_para('G:\example');
[~, ~, datasets] = elf_checkdata(para);
res              = zeros(1, length(datasets));

for i = 1:length(datasets)
    fprintf('\n%d of %d\n', i, length(datasets));
        dataset = datasets{i};
        elf_main1_HdrAndInt(dataset); close all %, [], 0, 1
%         elf_main2_meanimage(dataset); close all
%         elf_main3_intsummary(dataset); close all
end

return
%%
elf_paths
clear
para             = elf_para('F:\All data');
[~, ~, datasets] = elf_checkdata(para);
res              = zeros(1, length(datasets));

for i = 1:length(datasets)
    fprintf('\n%d of %d\n', i, length(datasets));
    try
        dataset = datasets{i};
%         elf_main1_HdrAndInt(dataset); close all %, [], 0, 1
%         elf_main2_meanimage(dataset); close all
        elf_main3_intsummary(dataset); close all
        res(i) = 1;
    catch ME
        warning(ME.identifier, ME.message)
    end
end

% errors (between 327-705): 377   400   401   402   432   433   444   466   477   488   564
% errors (between 706-end): 753   754   755   756   757   758   759   760   761   763   765   770   772   777   780   797   799   812   823   834 921
%     '356 Hö Väg 27 July noon overcast'
%     '377 Knivsås Bok 12 Sept noon sun'
%     '378 Knivsås Björk 12 Sept noon sun'
%     '379 Knivsås Ene 12 Sept noon sun'
%     '4 UW Turkey March 2008'
%     '40 Hörröd vägen sol efterm 28okt14'
%     '41 Hörröd väg mulet noon 31okt14'
%     '44 Hörröd höstscen 28 okt 2014'
%     '45-59 Hörröd utsikt klart skymn 10 min okt 14'
%     '5 UW Spain July 08'

%     '6 Kohagar mam 29augkl1530'
%     '60 Hörröd utsikt mulet middag 31 okt 2014'
%     '60B Nellie TIFF South Africa HT 2014'
%     '61 H-utsikt Foggy night, ISO series'
%     '62 H-utsikt Night, shutter time series'
%     '63 H-utsikt 9nov 4am 4 5 moon fog 2'
%     '64 H-väg 9nov 4am 4 5 moon fog'
%     '65 H-utsikt 9 nov sunrise fog 2'
%     '66 H-väg 9 nov sunrise fog'
%     '68 noise tests'
%     '7 Kohagar rev 2septkl1300'
%     '73 Bokskog knivsåsen tråkrå januari 2015'
%     '75 Eneryggen knivsåsen tråkrå januari 2015'
%     '8 Kostallar vintern 2013'
%     '82 Hörröd-väg snö halvklart'
%     '88 Hö-utsikt feb natt låg måne halvklart'
%     '9 Lövskog vinter'
%     'e1 Astrand morrumsan 11mars1315'
%     'e2 Granskog fridafors 10mars1630'
%     'e3 Hasthage hovmansbygd 10mars1730'
%     'j34g Raspberry hill 2145'

%%
% clear
% para = elf_para('E:\James data');
% [~, ~, datasets] = elf_checkdata(para);
% 
% res = zeros(1, length(datasets));
% 
% for i = 1:length(datasets) % 16, 192a,b, 84a-l %%141 177:178 423:434 481:523
%     fprintf('\n%d of %d\n', i, length(datasets));
%     try
%         dataset = datasets{i};
%         elf_main1_HdrAndInt(dataset); close all
%         elf_main2_meanimage(dataset); close all
%         elf_main3_intsummary(dataset); close all
%         res(i) = res(i) + 1;
%     catch ME
%         warning(ME.identifier, ME.message)
%     end
% end




%%

% para = elf_para('F:\All data');
% [~, ~, datasets] = elf_checkdata(para);
% 
% for i = 411:418
%     dataset = datasets{i}
%     elf_main5_contrasts(dataset);
%         
% end
% para = elf_para('F:\All data');
% [~, ~, datasets] = elf_checkdata(para);
% 
% try
%     elf_postana_combine(datasets([231:240]), '007 Random')
% catch me
%     warning('DATASET ABORTED!');
%     warning(me.message)
% end
% %%
% para = elf_para('F:\New data');
% [~, ~, datasets] = elf_checkdata(para);
% 
% for i = 1:length(datasets)
%     dataset = datasets{i};
%     elf_main1_HDRscenes(dataset); close all
%     elf_main2_meanimage(dataset); close all
% end
% for i = 1:length(datasets)
%     dataset = datasets{i};
%     elf_main3_luminance(dataset); close all
%     elf_main3p5_intsummary(dataset); close all
% end
% for i = 1:length(datasets)
%     dataset = datasets{i};
%     elf_main4_filter(dataset); close all
%     elf_main5_contrasts(dataset); close all
%     elf_main6_summary(dataset); close all
%     elf_main7_stats_and_plots(dataset);
% end
% 
% return;
% % 
% %%
% 
% para = elf_para('F:\test');
% [~, ~, datasets] = elf_checkdata(para);
% 
% for i = 9%[248 262 268 222 359 244 265 297 217:221]%i = 1:length(datasets)
%     dataset = datasets{i};
%     try
%         %elf_main5_contrasts(dataset); close all
%         %elf_main6_summary(dataset); close all
%         elf_main7_stats_and_plots(dataset);
%     catch me
%         warning('DATASET ABORTED!');
%         warning(me.message)
%     end
% end
% 
% 
% return
% %% DO LATER
% 
% para = elf_para('F:\All data');
% [~, ~, datasets] = elf_checkdata(para);
% 
% for i=[359 381 406 408 410 362 365 372]%i = 1:length(datasets)
%     dataset = datasets{i};
%     try
%         elf_main3p5_intsummary(dataset); close all
%     catch me
%         warning('DATASET ABORTED!');
%         warning(me.message)
%     end
% end

%%
% para = elf_para('F:\All data');
% [~, ~, datasets] = elf_checkdata(para);
% 
% try
%     elf_postana_combine(datasets, '000 Earth')
% catch me
%     warning('DATASET ABORTED!');
%     warning(me.message)
% end
% 
% 
% 
% 
% 
% 
% return
% %%
% para = elf_para('F:\All data');
% [~, ~, datasets] = elf_checkdata(para);
% for i = 1:length(datasets)
%     try
%     para            = elf_para('', datasets{i});
%     infosum         = elf_readwrite(para, 'loadinfosum');      % loads the old infosum file (which contains projection information, and linims)
%     fprintf('Dataset: %s, %s to %s\n', datasets{i}, datestr(min(infosum.DateTimeOriginal)), datestr(max(infosum.DateTimeOriginal)));
%     end
% end