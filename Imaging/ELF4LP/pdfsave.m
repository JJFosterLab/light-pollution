function pdfsave(fh, filename)

%% check inputs
if nargin<1
    if verLessThan('matlab', '8.4')
        fignumber = gcf;
    else
        temp = gcf;
        fignumber = temp.Number;
    end
end
if isnumeric(fh)
    fignumber = fh;
    fh = figure(fignumber);
else
    fignumber = fh.Number;
end

if nargin<2
    filename = ['fig', num2str(fignumber)];
end
if nargin < 3, ori = ''; end

%%
figure(fh);

set(fh, 'Units', 'centimeters');
pos = get(fh, 'Position');
set(fh, 'PaperPositionMode', 'Auto', 'PaperUnits', 'centimeters', 'PaperSize', [pos(3), pos(4)]);

set(fh, 'renderer', 'painters') % for some graphs (e.g. circular dot hists), opengl saves "panels" of pixel data
print(fh, '-dpdf', '-r1200', filename);
disp(['Saved figure ' num2str(fignumber) ' to ' filename]);
