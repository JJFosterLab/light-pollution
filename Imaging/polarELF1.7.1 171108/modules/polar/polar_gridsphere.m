function [a, e] = polar_gridsphere
% POLAR_GRIDSPHERE produces a grid that regularly covers a sphere
%
% Based on GridSphere (https://se.mathworks.com/matlabcentral/fileexchange/28842-grid-sphere)
% N. Teanby   13-01-04    Original IDL (Interactive Data Language) code (available at http://www.atm.ox.ac.uk/user/teanby/software.html#icos)

numpoints   = 642; % can have 2+(10*(4^k)): 42, 162, 642, 2562, ...
p           = fileparts(mfilename);
temp        = load(fullfile(p, 'polar_gridspheres', sprintf('gridsphere%d.mat', numpoints)));
a           = temp.a;
e           = temp.e;

% %% These were used to save mat-files
% [e, a] = GridSphere(42); save('gridsphere42.mat', 'a', 'e');
% [e, a] = GridSphere(162); save('gridsphere162.mat', 'a', 'e');
% [e, a] = GridSphere(642); save('gridsphere642.mat', 'a', 'e');
% [e, a] = GridSphere(2562); save('gridsphere2562.mat', 'a', 'e');