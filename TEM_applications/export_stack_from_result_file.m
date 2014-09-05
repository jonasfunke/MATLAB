%% startup
clc, clear all, close all
path0 = cd; addpath(path0); display(['Added search-path: ' path0 ])
run('my_prefs')

%% select file with results
cd(data_dir)
[filename pathname]=uigetfile('*.xls','Results.xls file.');
cd(path0)

%% load data
data = dlmread([pathname filename], '\t', 1, 0);

%% load header
tmp = textread([pathname filename], '%s', 'delimiter', '\n'); % read file by line
header = strsplit(tmp{1}); % seperate first line

tmp_index = 1:length(header);
index = tmp_index(strcmp(header, 'Slice'))+1; % index where header == Slice, +1 since first column is left out

%% get list of slices
slices = unique(data(:,index)); % get slices and remove dublicates

%% read tif files

base = 'ref_1_';
for i=1:length(slices)
    fpath = [pathname  base sprintf( '%03i', slices(i) ) '.tif'];
    img = imread(fpath);
    
    if i==1
        stack = zeros(size(img,1), size(img,2), length(slices), 'uint16');
    end
    stack(:,:,i) = img;
    
end

%%
WriteImagic(stack, [pathname 'stack']) % write stack to img file

