%% STARTUP  Auto-configure MATLAB path for this project
%  Place this file in the matlab/ directory. MATLAB runs it automatically
%  when you cd to this directory or start MATLAB from here.

root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root, 'src')));
addpath(fullfile(root, 'tests'));

fprintf('Paths configured for paper reproduction project.\n');
clear root;
