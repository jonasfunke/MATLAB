% This is the github version of the script. Changes are synced to the
% github repository and can be committed via the github interface.

clc, clear all, close all
path0 = cd; addpath(path0); display(['Using path: ' path0 ])
matlab_dir = userpath; matlab_dir = matlab_dir(1:end-1);
data_dir = '/Users/matthiasschickinger/PhD/TIRFM Data';
sep = filesep;

home_path = '/Users/matthiasschickinger/Documents';
movie_lib_path = [home_path sep 'MATLAB' sep 'MOVIE_LIBRARY' ];
dataanalyze_path = [home_path sep 'MATLAB' sep 'DATAANALYSIS_LIBRARY'];

addpath(movie_lib_path);
addpath(dataanalyze_path);

display(strcat('Using path: ', movie_lib_path))
display(strcat('Using path: ', dataanalyze_path))

%% LOAD MOVIES, DETERMINE SEQUENCE AND THRESHOLDS

%Determine desired data colors
rgb={'red','green','blue'};
[colors,ok]=listdlg('PromptString', 'Select two colors to be analyzed',...
                'ListString', rgb,...
                'OKString', 'Engage');
while ne(length(colors),2)
    [colors,ok]=listdlg('PromptString', 'Select _TWO_ colors to be analyzed',...
                'ListString', rgb,...
                'OKString', 'Engage');
end

% Deterine more parameters

options.WindowStyle = 'normal';
strw = inputdlg({'Integration radius:', 'Search radius:' , 'Minimal length:', 'Start frame:',...
    'End frame (put -1 for max):',['Sequence ' rgb{colors(1)} ':'], ['Sequence ' rgb{colors(2)} ':']}, 'Parameters' , 1, ...
    {'4','4','5','2','-1','1','1'}, options);
r_integrate = str2double(strw(1));
r_find = str2double(strw(2));
min_length = str2double(strw(3));
start_frame = str2double(strw(4));
end_frame = str2double(strw(5));

% Convert sequence strings to arrays

seq_c1 = zeros(1, size(strw{6},2));
for i=1:size(strw{6},2)
    if(strw{6}(i) == '1')
        seq_c1(1,i) =1;
    end
end

seq_c2 = zeros(1, size(strw{7},2));
for i=1:size(strw{7},2)
    if(strw{7}(i) == '1')
        seq_c2(1,i) =1;
    end
end

% Get data files
cd(data_dir)
[fname pname]=uigetfile('*.*',['Select ' rgb{colors(1)} ' data file']);

cd(pname);
cd ..
data_path = cd;

% Create movie objects

c1 = movie(pname, fname, start_frame, end_frame, seq_c1);
display([rgb{colors(1)} ' folder: ' c1.pname])
if c1.input==1
    display(['first ' rgb{colors(1)} ' frame: ' c1.fnames{c1.first}])
    display(['last ' rgb{colors(1)} ' frame: ' c1.fnames{c1.last}])
else
display([rgb{colors(1)} ' movie file: ' c1.fname])
end
pname = [pname(1:length(pname)-length(rgb{colors(1)})-1) rgb{colors(2)} filesep];
fname = [rgb{colors(2)} fname(length(rgb{colors(1)})+1:end)];

c2 = movie(pname, fname, start_frame, end_frame, seq_c2);
display([rgb{colors(2)} ' folder: ' c2.pname])
if c2.input==1
    display(['first ' rgb{colors(2)} ' frame: ' c2.fnames{c2.first}])
    display(['last ' rgb{colors(2)} ' frame: ' c2.fnames{c2.last}])
else
display([rgb{colors(2)} ' movie file: ' c2.fname])
end   

cd(path0)
display(['Output: ' data_path])

button = questdlg('plot results?','Plot','Yes','No','No');
plot_curves = strcmp(button,'Yes');

button = questdlg(['Map positions ' rgb{colors(1)} ' ON ' rgb{colors(2)} ' and vice versa?'],'Mapping','Yes','No','No');
mapping = strcmp(button, 'Yes');

if mapping == 1
    [mapping_file_1on2, mapping_dir]=uigetfile(data_dir,['Choose the ' rgb{colors(1)} ' ON ' rgb{colors(2)} ' mapping file:']);
    map1on2=load([mapping_dir mapping_file_1on2]);
    tform1on2=['tform_' rgb{colors(1)}(1) 'ON' rgb{colors(2)}(1)];
    display(['loaded ' rgb{colors(1)} ' ON ' rgb{colors(2)} ' mapping file: ' mapping_dir mapping_file_1on2]);
    
    [mapping_file_2on1]=uigetfile(mapping_dir,['Choose the ' rgb{colors(2)} ' ON ' rgb{colors(1)} ' mapping file:']);
    map2on1=load([mapping_dir mapping_file_2on1]);
    tform2on1=['tform_' rgb{colors(2)}(1) 'ON' rgb{colors(1)}(1)];
    display(['loaded ' rgb{colors(2)} ' ON ' rgb{colors(1)} ' mapping file: ' mapping_dir mapping_file_2on1]);
end
    
%% get thresholds

c1.get_h_min(r_find);
c2.get_h_min(r_find);

%% load and process traces

[c1_traces, c1_itraces, c1_avg_frame] = c1.trace_movie(c1.h_min, r_find, r_integrate, min_length);
[c2_traces, c2_itraces, c2_avg_frame] = c2.trace_movie(c2.h_min, r_find, r_integrate, min_length);


%% average position list from traces

c1_pos = zeros(size(c1_traces,1), 2);
for i=1:size(c1_traces,1)
    c1_pos(i,1) = mean(c1_traces{i,1}(:,2));
    c1_pos(i,2) = mean(c1_traces{i,1}(:,3));    
end
[c1_itraces_full] = c1.traces_movie_position(c1_pos, r_integrate);

c2_pos = zeros(size(c2_traces,1),2);
for i=1:size(c2_traces,1)
    c2_pos(i,1) = mean(c2_traces{i,1}(:,2));
    c2_pos(i,2) = mean(c2_traces{i,1}(:,3));    
end
[c2_itraces_full] = c2.traces_movie_position(c2_pos, r_integrate);

%%
%map average positions
if mapping == 1

pos1on2=tforminv(c1_pos, map1on2.(tform1on2));
pos2on1=tforminv(c2_pos, map2on1.(tform2on1));

end

%% combine traces to a cell array
trace_map = map_traces(c1_pos, c2_pos, c2_pos, r_find); %map the traces from average positions
% combine in merged_traces
merged_traces = cell(size(trace_map,1),2);
merged_itraces = cell(size(trace_map,1),2);
for i=1:size(trace_map,1)
   merged_traces{i,1} = c1_traces{trace_map(i,1)+1}; % +1 because map_traces gives indices starting from 0
   merged_traces{i,2} = c2_traces{trace_map(i,2)+1};
   
   merged_itraces{i,1} = c1_itraces_full{trace_map(i,1)+1}; % +1 because map_traces gives indices starting from 0
   merged_itraces{i,2} = c2_itraces_full{trace_map(i,2)+1};
end
   
%% FIT GAUSSIAN TO DATA % MAP 1 On 2 and 2 ON 1 and store in merged_traces col 4:5
merged_traces_fit = cell(size(merged_traces));

for i=1:size(merged_traces,1) % loop through traces   
    display(['Fitting spot ' num2str(i) ' of ' num2str(size(merged_traces,1))])
    tmp = zeros(size(merged_traces{i,1},1), size(merged_traces{i,1},2)+9);
    tmp(:,1:3)= merged_traces{i,1};
    tmp(:,6:12) = c1.fit_psf_to_movie(merged_traces{i,1}(:,1), merged_traces{i,1}(:,2:3), 2); %fit spot in each frame, sigma = 2
    if mapping == 1
    tmp(:,4:5) = tforminv(tmp(:,6:7), map1on2.(tform1on2));
    end
    merged_traces_fit{i,1} = tmp;    
    
    tmp = zeros(size(merged_traces{i,2},1), size(merged_traces{i,2},2)+9);
    tmp(:,1:3)= merged_traces{i,2};
    tmp(:,6:12) = c2.fit_psf_to_movie(merged_traces{i,2}(:,1), merged_traces{i,2}(:,2:3), 2); %fit spot in each frame, sigma = 2
    if mapping == 1
    tmp(:,4:5) = tforminv(tmp(:,6:7), map2on1.(tform2on1));
    end
    merged_traces_fit{i,2} = tmp;    
end


%% save data
display('Writing data...');
tmp_out = [c1.fname(length(rgb{colors(1)})+4:end-9+c1.input) c1.fname(length(rgb{colors(1)})+1:length(rgb{colors(1)})+2)];
cd(data_path);
mkdir(tmp_out);
cd(tmp_out);
mkdir([rgb{colors(1)}(1) rgb{colors(2)}(1)]);
cd([rgb{colors(1)}(1) rgb{colors(2)}(1)]);
path_out = cd;

%write to mat-file for further processing
save([tmp_out '_' rgb{colors(1)}(1) rgb{colors(2)}(1) '_data.mat'], 'merged_traces', 'merged_traces_fit','merged_itraces', 'c1_traces', 'c1_itraces', 'c1_itraces_full','c2_traces', 'c2_itraces','c2_itraces_full',  'c1_avg_frame', 'c2_avg_frame');
display('done')

%% Plot data
% traces for delta_r and sigma over frame number
cf = figure(1)
for i=1:size(merged_traces,1)
    
    xy_c1 = merged_traces_fit{i,1}(:,6:7);
    xy_mean_c1 = [mean(xy_c1(:,1)) mean(xy_c1(:,2))  ];
    d_c1 = sqrt((xy_c1(:,1)-xy_mean_c1(1)).^2 + (xy_c1(:,2)-xy_mean_c1(2)).^2);
    
    sigma_c1 = sqrt(merged_traces_fit{i,1}(:,8).^2  + merged_traces_fit{i,1}(:,9).^2);
    frame_c1 = merged_traces_fit{i,1}(:,1);
    chi2_c1 = merged_traces_fit{i,1}(:,12);
    meanchi2_c1=mean(chi2_c1);
    
    xy_c2 = merged_traces_fit{i,2}(:,6:7);
    xy_mean_c2 = [mean(xy_c2(:,1)) mean(xy_c2(:,2))  ];
    d_c2 = sqrt((xy_c2(:,1)-xy_mean_c2(1)).^2 + (xy_c2(:,2)-xy_mean_c2(2)).^2);
    
    sigma_c2 = sqrt(merged_traces_fit{i,2}(:,8).^2  + merged_traces_fit{i,2}(:,9).^2);
    frame_c2 = merged_traces_fit{i,2}(:,1);
    chi2_c2 = merged_traces_fit{i,2}(:,12);
    meanchi2_c2 = mean(chi2_c2);

    subplot(4, 1, 1)
    plot(frame_c1, d_c1, ['.-' rgb{colors(1)}(1)],frame_c2,d_c2, ['.-' rgb{colors(2)}(1)], 'Markersize', 15)
    legend(rgb{colors(1)}, rgb{colors(2)})
    title(['Distance from average position. Spot: ' num2str(i)])
    xlabel('Frame')
    ylabel('$\sqrt{{\Delta x}^2 +{\Delta y}^2}$', 'Interpreter', 'latex')
    ylim([0 5])
    
    subplot(4, 1, 2)
    plot(frame_c1, sigma_c1, ['.-' rgb{colors(1)}(1)],frame_c2,sigma_c2, ['.-' rgb{colors(2)}(1)], 'Markersize', 15)
    legend(rgb{colors(1)}, rgb{colors(2)})
    title(['Width of peak. Spot: ' num2str(i)])    
    xlabel('Frame')
    ylabel('$\sqrt{{\sigma_x}^2+{\sigma_y}^2}$', 'Interpreter', 'latex')
    ylim([0 5])
    
    subplot(4, 1, 3)
    %bar(frame_green, chi2green)
    plot(frame_c1, chi2_c1, ['.-' rgb{colors(1)}(1)], 'Markersize', 15)
    title(['Chi^2 of ' rgb{colors(1)} ' fit. Mean value: ' sprintf('%.2s',meanchi2_c1) '. Spot: ' num2str(i)])    
    xlabel('Frame')
    ylabel(['${\chi^2}_{' rgb{colors(1)} '}$'], 'Interpreter', 'latex')
    
    subplot(4, 1, 4)
    plot(frame_c2, chi2_c2, ['.-' rgb{colors(2)}(1)], 'Markersize', 15)
    title(['Chi^2 of ' rgb{colors(2)} ' fit. Mean value: ' sprintf('%.2s',meanchi2_c2) '. Spot: ' num2str(i)])    
    xlabel('Frame')
    ylabel(['${\chi^2}_{' rgb{colors(2)} '}$'], 'Interpreter', 'latex')
    
    print(cf, '-depsc2', [path_out filesep 'trace_' num2str(i) '.eps' ])
end

close all

%% traces for position (dots unconnected)
delta = 5;
cf = figure(1)
for i=1:size(merged_traces,1)
    
    xy_c1 = merged_traces_fit{i,1}(:,6:7);
    xy_mean_c1 = [mean(xy_c1(:,1)) mean(xy_c1(:,2))  ];
    d_c1 = sqrt((xy_c1(:,1)-xy_mean_c1(1)).^2 + (xy_c1(:,2)-xy_mean_c1(2)).^2);
    
    xy_c2 = merged_traces_fit{i,2}(:,6:7);
    xy_mean_c2 = [mean(xy_c2(:,1)) mean(xy_c2(:,2))  ];
    d_c2 = sqrt((xy_c2(:,1)-xy_mean_c2(1)).^2 + (xy_c2(:,2)-xy_mean_c2(2)).^2);
    
    plot(xy_c1(:,1), xy_c1(:,2), [rgb{colors(1)}(1) 'o']), hold on
    plot(xy_c2(:,1), xy_c2(:,2), [rgb{colors(2)}(1) 'x']), hold off
    axis([floor(xy_mean_c1(1)-delta)  ceil(xy_mean_c1(1)+delta) floor(xy_mean_c1(2)-delta)  ceil(xy_mean_c1(2)+delta)]);
    x_lim=xlim;
    axis equal; 
    set(gca, 'XTick',[x_lim(1):1:x_lim(2)]);
    grid on;
    title([rgb{colors(1)} ' and ' rgb{colors(2)} ' fitted positions (unmapped). Spot: ' num2str(i)]);
    hold off;
    axis([floor(xy_mean_c1(1)-delta)  ceil(xy_mean_c1(1)+delta) floor(xy_mean_c1(2)-delta)  ceil(xy_mean_c1(2)+delta)]);
    print(cf, '-depsc2', [path_out filesep 'position_' num2str(i) '.eps' ])
end

close all


%% traces for position (dots unconnected, red dots mapped rONg)
if mapping == 1
delta = 5;
cf = figure(1)
for i=1:size(merged_traces,1)
    
    xy_c1 = merged_traces_fit{i,1}(:,6:7);
    xy_mean_c1 = [mean(xy_c1(:,1)) mean(xy_c1(:,2))  ];
    d_c1 = sqrt((xy_c1(:,1)-xy_mean_c1(1)).^2 + (xy_c1(:,2)-xy_mean_c1(2)).^2);
    
    xy_2on1 = merged_traces_fit{i,2}(:,4:5);
    xy_mean_2on1 = [mean(xy_2on1(:,1)) mean(xy_2on1(:,2))  ];
    d_rONg = sqrt((xy_2on1(:,1)-xy_mean_2on1(1)).^2 + (xy_2on1(:,2)-xy_mean_2on1(2)).^2);
    
    plot(xy_c1(:,1), xy_c1(:,2), [rgb{colors(1)}(1) 'o']), hold on
    plot(xy_2on1(:,1), xy_2on1(:,2), [rgb{colors(2)}(1) 'x']), hold off
    axis([floor(xy_mean_c1(1)-delta)  ceil(xy_mean_c1(1)+delta) floor(xy_mean_c1(2)-delta)  ceil(xy_mean_c1(2)+delta)]);
    x_lim=xlim;
    axis equal; 
    set(gca, 'XTick',[x_lim(1):1:x_lim(2)]);
    grid on;
    title([rgb{colors(1)} ' and ' rgb{colors(2)} ' fitted positions (rONg mapped). Spot: ' num2str(i)]);
    hold off;
    axis([floor(xy_mean_c1(1)-delta)  ceil(xy_mean_c1(1)+delta) floor(xy_mean_c1(2)-delta)  ceil(xy_mean_c1(2)+delta)]);
    print(cf, '-depsc2', [path_out filesep 'position_mapped_rONg_' num2str(i) '.eps' ])
end
end

close all

%% traces for position (dots connected)
delta = 5;
cf = figure(1)
for i=1:size(merged_traces,1)
    
    xy_c1 = merged_traces_fit{i,1}(:,6:7);
    xy_mean_c1 = [mean(xy_c1(:,1)) mean(xy_c1(:,2))  ];
    d_c1 = sqrt((xy_c1(:,1)-xy_mean_c1(1)).^2 + (xy_c1(:,2)-xy_mean_c1(2)).^2);
    
    xy_c2 = merged_traces_fit{i,2}(:,6:7);
    xy_mean_c2 = [mean(xy_c2(:,1)) mean(xy_c2(:,2))  ];
    d_c2 = sqrt((xy_c2(:,1)-xy_mean_c2(1)).^2 + (xy_c2(:,2)-xy_mean_c2(2)).^2);
    
    plot(xy_c1(:,1), xy_c1(:,2), [rgb{colors(1)}(1) '-o']), hold on
    plot(xy_c2(:,1), xy_c2(:,2), [rgb{colors(2)}(1) '-x']), hold off
    axis([floor(xy_mean_c1(1)-delta)  ceil(xy_mean_c1(1)+delta) floor(xy_mean_c1(2)-delta)  ceil(xy_mean_c1(2)+delta)]);
    x_lim=xlim;
    axis equal; 
    set(gca, 'XTick',[x_lim(1):1:x_lim(2)]);
    grid on;
    title([rgb{colors(1)} ' and ' rgb{colors(2)} ' fitted positions (unmapped). Spot: ' num2str(i)]);
    hold off;
    axis([floor(xy_mean_c1(1)-delta)  ceil(xy_mean_c1(1)+delta) floor(xy_mean_c1(2)-delta)  ceil(xy_mean_c1(2)+delta)]);
    print(cf, '-depsc2', [path_out filesep 'position2_' num2str(i) '.eps' ])
end
close all
hold off

%%
% Plot positions of red and green spots on green average image
cf=figure(1);
imagesc(c1_avg_frame),colormap gray, colorbar;
hold on
for i=1:size(c1_pos,1)
plot(c1_pos(i,1),c1_pos(i,2),[rgb{colors(1)}(1) 'o']);
end
for i=1:size(c2_pos,1)
plot(c2_pos(i,1),c2_pos(i,2),[rgb{colors(2)}(1) 'x']);
end

print(cf, '-depsc2', [path_out filesep 'spots_' rgb{colors(1)} '_' rgb{colors(2)} '.eps'])
hold off

%%
% Plot positions of red and green spots on red average image
cf=figure(1);
imagesc(c2_avg_frame),colormap gray, colorbar;
hold on;
for i=1:size(c1_pos,1)
plot(c1_pos(i,1),c1_pos(i,2),[rgb{colors(1)}(1) 'x']);
end
for i=1:size(c2_pos,1)
plot(c2_pos(i,1),c2_pos(i,2),[rgb{colors(2)}(1) 'o']);
end

print(cf, '-depsc2', [path_out filesep 'spots_' rgb{colors(2)} '_' rgb{colors(1)} '.eps'])
hold off

%%
if mapping == 1
% Plot positions of green and mapped red spots on green average image
cf=figure(1);
imagesc(c1_avg_frame),colormap gray, colorbar;
hold on
for i=1:size(c1_pos,1)
plot(c1_pos(i,1),c1_pos(i,2),[rgb{colors(1)}(1) 'o']);
end
for i=1:size(pos1on2,1)
plot(pos2on1(i,1),pos2on1(i,2),[rgb{colors(2)}(1) 'x']);
end

print(cf, '-depsc2', [path_out filesep 'spots_' rgb{colors(1)} '_' rgb{colors(2)}(1) 'ON' rgb{colors(1)}(1) '.eps'])
hold off

% Plot positions of red and mapped green spots on red average image
cf=figure(1);
imagesc(c2_avg_frame),colormap gray, colorbar;
hold on;
for i=1:size(pos2on1,1)
plot(pos2on1(i,1),pos2on1(i,2),[rgb{colors(1)}(1) 'x']);
end
for i=1:size(c2_pos,1)
plot(c2_pos(i,1),c2_pos(i,2),[rgb{colors(2)}(1) 'o']);
end

print(cf, '-depsc2', [path_out filesep 'spots_' rgb{colors(2)} '_' rgb{colors(1)}(1) 'ON' rgb{colors(2)}(1) '.eps'])
hold off
end

close all

