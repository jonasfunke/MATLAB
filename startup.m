display('----------------------------------------------------------------')
display('------------------------- NEVER GIVE UP ------------------------')
display('----------------------------------------------------------------')

% load general preferences
run('my_prefs')
%%
% include path to environment
addpath(genpath([matlab_dir filesep 'TOOLBOX_GENERAL']))
display( ['Added search-path: ' matlab_dir filesep 'TOOLBOX_GENERAL' ])

addpath(genpath([matlab_dir filesep 'TOOLBOX_TEM']))
display( ['Added search-path: ' matlab_dir filesep 'TOOLBOX_TEM' ])

addpath(genpath([matlab_dir filesep 'TOOLBOX_MOVIE']))
display( ['Added search-path: ' matlab_dir filesep 'TOOLBOX_MOVIE' ])
%%
