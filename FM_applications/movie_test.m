addpath(genpath([matlab_dir filesep 'TOOLBOX_MOVIE']))
display( ['Added search-path: ' matlab_dir filesep 'TOOLBOX_MOVIE' ])

%%


path0 = cd;
cd(data_dir)
[fname pname] = uigetfile('.fits', 'Select a .fits file');
cd(path0)

%%
run('my_prefs')
path0 = cd;
cd(data_dir)
[fname pname] = uigetfile('.tif', 'Select a .fits file');
cd(path0)

%% init red movie
red = movie(pname, fname, 2, -1, [1 0]); % pname, fname, first, last, sequence

%%
red.get_h_min(4)
%%
r_find = 4;
r_integrate = 4;
min_length = 3;
[t, it, avg] = red.trace_movie(red.h_min, r_find, r_integrate, min_length );

%%
imagesc(avg), colorbar, axis image, colormap gray, hold on
for i=1:size(t,1)
    plot(mean(t{i}(:,2))+1, mean(t{i}(:,3))+1, 'ro')
end

%% read red movie
go_on = 1;
red.initRead;

while go_on
    [mov, frames, go_on]  = red.readNext;
end





