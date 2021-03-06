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

addpath(genpath([matlab_dir filesep 'TOOLBOX_TYPHOON']))
display( ['Added search-path: ' matlab_dir filesep 'TOOLBOX_TYPHOON' ])

addpath(genpath([matlab_dir filesep 'Korbinian']))
display( ['Added search-path: ' matlab_dir filesep 'Korbinian' ])

%addpath(genpath([matlab_dir filesep 'slurm_shared']))
%display( ['Added search-path: ' matlab_dir filesep 'slurm_shared' ])

addpath(genpath([matlab_dir filesep 'slurm_nonshared']))
display( ['Added search-path: ' matlab_dir filesep 'slurm_nonshared' ])
    
    %%
