%%
clear all
close all
clc
run('my_prefs')
path_init = cd;

cd(data_dir)
display('Choose the green data tiff-folder!');
path_green=uigetdir(cd,'Choose the green data tiff-folder:');
cd(path_green);
cd ..
cd ..
cd ..
data_path = cd;
display('Choose the red data tiff-folder!');
path_red=uigetdir(cd,'Choose the red data tiff-folder:');

options.WindowStyle = 'normal';

folders = regexp(path_green, filesep, 'split');
structure = folders{size(folders,2)-1};
N_bla = inputdlg({strcat('How many movies of ',structure, ':')}, structure , 1, {'1'}, options );

N = str2double(N_bla{1});

data_out = [data_path filesep structure];


display(strcat('Source green: ', path_green));
display(strcat('Destination green: ', data_out));
display(strcat('Source red: ', path_red));
display(strcat('Destination red: ', data_out));

button = questdlg('Source and Destination OK?','Source Destination','OK','Cancel','Cancel');

if strcmp(button,'OK') %load frames

    
    sepindex = find(path_green == filesep);
    path_sif_green = path_green(1:sepindex(end)-1);
    
    sepindex = find(path_red == filesep);
    path_sif_red = path_red(1:sepindex(end)-1);
    
mkdir(data_out);%make output folder folder
%%
for i=1:N
   display(['Moving folder ' sprintf('%.2i',i)] )
   green_out = [data_out filesep sprintf('%.2i',i) '_green' ]; 
   red_out = [data_out filesep sprintf('%.2i',i) '_red' ]; 
   mkdir(green_out); 
   mkdir(red_out);
   movefile([path_green filesep 'green' sprintf('%.2i',i) '*'], green_out)
   movefile([path_red filesep 'red' sprintf('%.2i',i) '*'], red_out)
   
   %move .sif file
   movefile([path_sif_green filesep 'green' sprintf('%.2i',i) '*.sif'] , data_out); % green file
   movefile([path_sif_red filesep 'red' sprintf('%.2i',i) '*.sif'] , data_out); % red file
   
   
end
%%
display('done moving')

else
    display('Moving data aborded.')
end
cd(path_init);