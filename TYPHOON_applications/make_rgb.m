%% startup
clear all, close all, clc
run('my_prefs')
path0 = cd;

%% Select images 
cd(data_dir)
[filename_r pathname_r]=uigetfile('*.tif','Select RED-channel');
cd(pathname_r)
[filename_g pathname_g]=uigetfile('*.tif','Select GREEN-channel');
cd(path0)

%% Load images
red_raw = double(imread([pathname_r filesep filename_r])); 
green_raw = double(imread([pathname_g filesep filename_g])); 

%% create output folder
prefix_out = [filename_r(1:end-4) '_rgb'];
path_out = [pathname_r filesep prefix_out ];
mkdir(path_out);

%%
rgb = zeros(size(red_raw,1), size(red_raw,2), 3);
rgb(:,:,1) = red_raw;
rgb(:,:,2) = green_raw;    
imwrite(uint16(rgb), [path_out filesep 'image_rgb.tif'])
   
%% 
scale_factor_red = 6000;
scale_factor_green = 2500;

red_scaled = red_raw .* ((2^16-1)/scale_factor_red);
green_scaled = green_raw .* ((2^16-1)/scale_factor_green);

rgb_scaled = zeros(size(red_raw,1), size(red_raw,2), 3);
rgb_scaled(:,:,1) = red_scaled;
rgb_scaled(:,:,2) = green_scaled;    
imwrite(uint16(rgb_scaled), [path_out filesep 'image_rgb_scaled.tif'])

%%
disp('done')
