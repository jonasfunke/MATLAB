%% startup
clc, clear all, close all
path0 = cd;
run( 'my_prefs.m')

%% LOAD MOVIES, DETERMINE SEQUENZ AND THERSHOLDS

%load green movie
cd(data_dir)

pathname=uigetdir(data_dir,'Choose the green data-folder:');
sepindex=find(pathname==filesep);
foldername = pathname(sepindex(length(sepindex))+1:length(pathname) );

if  strcmp(foldername(4:end), 'green')
    files_tmp = dir([pathname filesep '*.tif']);
    Im_name1G = files_tmp(1).name;
    Im_name2G = files_tmp(end).name;
    pathname1G = pathname;

    
    pathname1R = [pathname(1:sepindex(length(sepindex)))  foldername(1:3) 'red'];
    files_tmp = dir([pathname1R filesep '*.tif']);
    Im_name1R = files_tmp(1).name;
    Im_name2R = files_tmp(end).name;
    
else
    files_tmp = dir([pathname filesep '*.tif']);
    Im_name1R = files_tmp(1).name;
    Im_name2R = files_tmp(end).name;
    pathname1R = pathname;

    
    pathname1G = [pathname(1:sepindex(length(sepindex)))  foldername(1:3) 'green'];
    files_tmp = dir([pathname1G filesep '*.tif']);
    Im_name1G = files_tmp(1).name;
    Im_name2G = files_tmp(end).name;
end
cd(path0)
data_path = pathname(1:sepindex(length(sepindex))-1);
display(['Folder green: ' pathname1G])
display(['Folder red: ' pathname1R])
display(['First green: ' Im_name1G])
display(['Last green: ' Im_name2G])
display(['First red: ' Im_name1R])
display(['Last red: ' Im_name2R])
display(['Output: ' data_path])

button = questdlg('Files and folder OK?','ok','OK','No','No');
if strcmp(button,'No')
    cd(data_dir)
    [Im_name1G pathname1G]=uigetfile('*.tif','Select First Image of Greeen Movie (tiff stack) : ');
    cd(pathname1G)
    [Im_name2G pathname2G]=uigetfile('*.tif','Select Last Image of Green Movie (tiff stack) : ');
    cd ..
    data_path = cd;
    [Im_name1R pathname1R]=uigetfile('*.tif','Select First Image of Red Movie (tiff stack) : ');
    cd(pathname1R)
    [Im_name2R pathname2R]=uigetfile('*.tif','Select Last Image of Red Movie (tiff stack) : ');
    cd(path0)

end


%% 


strw = inputdlg({'Choose an frame offset:'}, 'Frame offset' , 1, {'1'} );
frame_offset = str2double(strw(1));


% picking number of first image
nIm1G=Im_name1G((length(Im_name1G)-7):length(Im_name1G)-4) ; 
nIm2G=Im_name2G((length(Im_name2G)-7):length(Im_name2G)-4);

nIm1R=Im_name1R((length(Im_name1R)-7):length(Im_name1R)-4) ;
nIm2R=Im_name2R((length(Im_name2R)-7):length(Im_name2R)-4);

N1G=str2double(nIm1G)+ frame_offset;
N2G=str2double(nIm2G);
NnG=N2G-N1G+1;  %number of frames of green movie

N1R=str2double(nIm1R)+ frame_offset;
N2R=str2double(nIm2R);
NnR=N2R-N1R+1;  %number of frames of red movie
 
Nn=min([NnG NnR]);    %number of frames


q_bin =2; % no binnning
q_bpass = 2; % no filtering

prefix_green =  Im_name1G(1:(length(Im_name1G)-9));
prefix_red = Im_name1R(1:(length(Im_name1R)-9));
pathname_green = pathname1G;
pathname_red = pathname1R;
green_start = N1G; green_stop = Nn;
red_start = N1G; red_stop = Nn;
fret_start = N1G; fret_stop = Nn;



%% Deterine more options

options.WindowStyle = 'normal';
strw = inputdlg({'Integration radius:', 'Radius of peak:' , 'Minimal length:', 'Prefix:', 'Sequence green:' ,'Sequence red:' }, 'Integration radius' , 1, {'4', '4','5' , [prefix_green(9:length(prefix_green)) sprintf('_%02s', Im_name1G(6:7) )], '10', '01'}, options );
r_integrate = str2double(strw(1));
r_find = str2double(strw(2));
min_length = str2double(strw(3));

folder_out = strw{4} ;
prefix_out = strw{4};

sequence_green = zeros(1, size(strw{5},2));
sequence_fret = zeros(1, size(strw{5},2));
for i=1:size(strw{5},2)
    if(strw{5}(i) == '1')
        sequence_fret(1,i) =1;
        sequence_green(1,i) =1;
    end
end

sequence_red = zeros(1, size(strw{6},2));
for i=1:size(strw{6},2)
    if(strw{6}(i) == '1')
        sequence_red(1,i) =1;
    end
end




button = questdlg('Make .avi movies?','AVI-Movies','Yes','No','No');
make_avi = strcmp(button,'Yes');

button = questdlg('plot results?','Plot','Yes','No','No');
plot_curves = strcmp(button,'Yes');

%% get threshold
h_min_green = get_h_min(pathname_green, filesep, prefix_green, green_start, green_stop, sequence_green, 'donor', r_find);

%%
h_min_red = get_h_min(pathname_red, filesep, prefix_red, red_start, red_stop, sequence_red, 'acceptor', r_find);
%h_min_fret = get_h_min(pathname_red, filesep, prefix_red, fret_start, fret_stop, sequence_fret, 'fret', r_find);
pause(0.5)
close all
%% load and process traces
display('--------------------- TRACING SPOTS ----------------------')
[green_traces, green_itraces, green_avg_frame] = trace_movie(pathname_green, filesep, prefix_green, green_start, green_stop, sequence_green, h_min_green, r_find, r_integrate, min_length);
[red_traces, red_itraces, red_avg_frame] = trace_movie(pathname_red, filesep, prefix_red, red_start, red_stop, sequence_red, h_min_red, r_find, r_integrate, min_length);
%% position list from red_traces
display('--------------------- INTEGRATING SPOTS ----------------------')
red_pos = zeros(size(red_traces,1), 2);
for i=1:size(red_traces,1)
    red_pos(i,1) = mean(red_traces{i,1}(:,2));
    red_pos(i,2) = mean(red_traces{i,1}(:,3));    
end
[fret_itraces_full, fret_avg_frame] = traces_from_position(pathname_red, filesep, prefix_red, fret_start, fret_stop, sequence_fret, red_pos, r_integrate);
[red_itraces_full, red_avg_frame2] = traces_from_position(pathname_red, filesep, prefix_red, fret_start, red_stop, sequence_red, red_pos, r_integrate);

%create avergage position of green
green_pos = zeros(size(green_traces,1),2);
for i=1:size(green_traces,1)
    green_pos(i,1) = mean(green_traces{i,1}(:,2));
    green_pos(i,2) = mean(green_traces{i,1}(:,3));    
end
[green_itraces_full, green_avg_frame2] = traces_from_position(pathname_green, filesep, prefix_green, green_start, green_stop, sequence_green, green_pos, r_integrate);



%{
   %plot green traces some results
   
   colormap('Gray');
   subplot(1,3,1)
   imagesc(green_avg_frame), colorbar, title('Green image'), axis([0 size(green_avg_frame,1) 0 size(green_avg_frame,1)]), hold all
   ColorSet = varycolor(size(green_traces,1));
   for i=1:size(green_traces)
    plot(green_traces{i}(:,2)+1, green_traces{i}(:,3)+1, 'Color', ColorSet(i,:));
   end
   
   
    %plot red traces
    colormap('Gray');
   subplot(1,3,2)
   imagesc(red_avg_frame), colorbar, title('Red image'), axis([0 size(red_avg_frame,1) 0 size(red_avg_frame,1)]), hold all
   ColorSet = varycolor(size(red_traces,1));
   for i=1:size(red_traces)
    plot(red_traces{i}(:,2)+1, red_traces{i}(:,3)+1, 'Color', ColorSet(i,:));
   end
   
    %plot fret traces
    colormap('Gray');
   subplot(1,3,3)
   imagesc(fret_avg_frame), colorbar, title('FRET image'), axis([0 size(fret_avg_frame,1) 0 size(fret_avg_frame,1)]), hold all
   ColorSet = varycolor(size(fret_itraces_full,1));
   for i=1:size(fret_itraces_full)
    plot(fret_itraces_full{i}(:,2)+1, fret_itraces_full{i}(:,3)+1, 'Color', ColorSet(i,:));
   end
   
   %}
   %% map green trace to red coordinates
%load mapping coefficients


coeff_R2G = load([movie_lib_path filesep 'Fmap_4_coeff_R2G.txt' ]);



%append two columns(x_redcam y_redcam) to green trace
for i=1:size(green_traces,1)
    tmp=zeros(size(green_traces{i}, 1), 5);
    xy = Fmap_4(coeff_R2G, [green_traces{i}(:,2) green_traces{i}(:,3)]);
    %x_r =Fmap4(coeff_GtoR(:,1),[green_traces{i}(:,2) green_traces{i}(:,3)]);
    %y_r =Fmap4(coeff_GtoR(:,2),[green_traces{i}(:,2) green_traces{i}(:,3)]);
    tmp(:,1:3)=green_traces{i};
    tmp(:,4)=xy(:,1);%x_r;
    tmp(:,5)=xy(:,2);%y_r;
    green_traces{i}=tmp;
end


%% combine and green, red and fret trace to a cell array traces
display('----------------------- MAPPING TRACES ------------------------')
trace_map = map_traces(green_pos, red_pos, red_pos, r_find); %map the tarces from averga positions
%%
% combine in merged_tarcer
merged_traces = cell(size(trace_map,1),3);
merged_itraces = cell(size(trace_map,1),3);
for i=1:size(trace_map,1)
   merged_traces{i,1} = green_traces{trace_map(i,1)+1};
   merged_traces{i,2} = red_traces{trace_map(i,2)+1};
   merged_traces{i,3} = fret_itraces_full{trace_map(i,3)+1};

   merged_itraces{i,1} = green_itraces_full{trace_map(i,1)+1};
   merged_itraces{i,2} = red_itraces_full{trace_map(i,2)+1};
   merged_itraces{i,3} = fret_itraces_full{trace_map(i,3)+1};
    
end

%{
%% plot traces
close all
hold off
fig1 = figure(1);
%plot average points
green_avg=average_pos(green_traces); %use int because red coordinates are double 
red_avg=average_pos(red_traces);
fret_avg=average_pos(fret_traces);
hold all
plot(green_avg(:,4), green_avg(:,5), 'gx'), hold on



plot(red_avg(:,2), red_avg(:,3), 'rx'), hold on



plot(fret_avg(:,2), fret_avg(:,3), 'bx'), hold on



for i=1:size(merged_traces, 1)
    plot(mean(merged_traces{i,1}(:,4)) , mean(merged_traces{i,1}(:,5)) , 'go' ) ;
    plot(mean(merged_traces{i,2}(:,2)) , mean(merged_traces{i,2}(:,3)) , 'ro' ) ;
    plot(mean(merged_traces{i,3}(:,2)) , mean(merged_traces{i,3}(:,3)) , 'bo' ) ;
   
end
%}


        
    
%% Plot FRET-Traces
%{
close all
hold off

p = 20; win=[5:10:50];
a = size(merged_itraces,1);
for i=1:a
    figure(i+fig1);
    subplot(4, 1, 1) %original data
    plot(merged_itraces{i,1}(:,1), merged_itraces{i,1}(:,4) ,'g', merged_itraces{i,2}(:,1), merged_itraces{i,2}(:,4),'r', merged_itraces{i,3}(:,1), merged_itraces{i,3}(:,4),'b' ),   hold on
    et = max(merged_itraces{i,1}(:,4)) * merged_itraces{i,3}(:,4) ./ (merged_itraces{i,3}(:,4)+ merged_itraces{i,1}(:,4));
    plot( merged_itraces{i,1}(:,1),et, 'black') , legend('I_{green}','I_{red}','I_{fret}', 'E_t'), title('Raw data')
    
    subplot(4, 1, 2) %nlf filtered data
    donor = nlf(merged_itraces{i,1}(:,4), win , p);
    acceptor = nlf(merged_itraces{i,2}(:,4), win , p);
    fret = nlf(merged_itraces{i,3}(:,4), win , p);
    plot(merged_itraces{i,1}(:,1),  donor, 'g', merged_itraces{i,2}(:,1), acceptor ,'r', merged_itraces{i,3}(:,1),fret ,'b' ),   hold on
    et = max(donor) * fret ./ (donor+ fret);
    plot( merged_itraces{i,1}(:,1),et, 'black') , legend('I_{green}','I_{red}','I_{fret}', 'E_t'), title('NLF data')
    
    subplot(4, 1, 3) %nlf filtered data
    [donor fret] = nlf_fret(merged_itraces{i,1}(:,4),merged_itraces{i,3}(:,4), win , p);
    acceptor = nlf(merged_itraces{i,2}(:,4), win , p);
    fret = nlf(merged_itraces{i,3}(:,4), win , p);
    plot(merged_itraces{i,1}(:,1),  donor, 'g', merged_itraces{i,2}(:,1), acceptor ,'r', merged_itraces{i,3}(:,1),fret ,'b' ),   hold on
    et = max(donor) * fret ./ (donor+ fret);
    plot( merged_itraces{i,1}(:,1),et, 'black') , legend('I_{green}','I_{red}','I_{fret}', 'E_t'), title('NLF_{FRET} data (I_{red} filtered with NLF )')
    
    subplot(4, 1, 4) %nlf filtered data
    [donor acceptor fret] = nlf_alex(merged_itraces{i,1}(:,4),merged_itraces{i,2}(:,4),merged_itraces{i,3}(:,4), win , p);
    plot(merged_itraces{i,1}(:,1),  donor, 'g', merged_itraces{i,2}(:,1), acceptor ,'r', merged_itraces{i,3}(:,1),fret ,'b' ),   hold on
    et = max(donor) * fret ./ (donor+ fret);
    plot( merged_itraces{i,1}(:,1),et, 'black') , legend('I_{green}','I_{red}','I_{fret}', 'E_t'), title('NLF_{ALEX} data')
end
%}

%% SAVE everything
display('------------------------- SAVING DATA --------------------------')

cd(data_path);
mkdir(folder_out);
path_out = strcat(data_path, filesep ,  folder_out);
cd(path_out);

close all
hold all
%write the traces
for i=1:size(merged_itraces,1)
    %green
    tmp_filename = strcat(prefix_out, '_', num2str(i),'_green.txt');
    fid = fopen(tmp_filename, 'w');
    fprintf(fid, '%i\t%i\t%i\t%i\n', transpose(merged_itraces{i,1}));
    fclose(fid);
    
    %red
    tmp_filename = strcat(prefix_out, '_', num2str(i),'_red.txt');
    fid = fopen(tmp_filename, 'w');
    fprintf(fid, '%i\t%i\t%i\t%i\n', transpose(merged_itraces{i,2}));
    fclose(fid);
    
    %fret
    tmp_filename = strcat(prefix_out, '_', num2str(i),'_fret.txt');
    fid = fopen(tmp_filename, 'w');
    fprintf(fid, '%i\t%i\t%i\t%i\n', transpose(merged_itraces{i,3}));
    fclose(fid);
    
end


%write to mat-file for further processing
tmp_filename = strcat(prefix_out, '_data','.mat');
save(tmp_filename, 'merged_traces','merged_itraces', 'green_traces', 'green_itraces', 'green_itraces_full','red_traces', 'red_itraces','red_itraces_full', 'fret_itraces_full', 'green_avg_frame', 'red_avg_frame', 'fret_avg_frame');

tmp_filename = strcat(prefix_out, '_info','.txt');
fid = fopen(tmp_filename, 'w');
fprintf(fid, 'Filter %i\n2x2binning %i\nframe_offset %i\n', [q_bpass q_bin frame_offset]);
%green
fprintf(fid, '#Green:\n');
fprintf(fid, 'start %i\nstop %i\n', [green_start green_stop]);
fprintf(fid, 'First_file %s\n\tLast_file = %s\n',  Im_name1G, Im_name2G);
fprintf(fid, 'Peakdetection_radius %i\nPeakdetection_Minimal_Intensity %.2f\nPeakdetection_Minimal_length %i\n', r_find, h_min_green, min_length);
fprintf(fid, 'sequence %s\n', strw{5});

%red
fprintf(fid, '#Red:\n');
fprintf(fid, 'start %i\nstop %i\n', [red_start red_stop]);
fprintf(fid, 'First_file %s\n\tLast_file = %s\n',  Im_name1R, Im_name2R);
fprintf(fid, 'Peakdetection_radius %i\nPeakdetection_Minimal_Intensity %.2f\nPeakdetection_Minimal_length %i\n', r_find, h_min_red, min_length);
fprintf(fid, 'sequence %s\n', strw{6});



fprintf(fid, '#Integration\nIntegration_Radius %i\n', r_integrate);
fclose(fid);


display(strcat('All traces written to' , path_out));
cd(path0)

%% Make Figures and save them to path_out
display('Creating some figures...')
cd(path_out);
scrsz = get(0,'ScreenSize');

close all



%figure with avergage position
fname = strcat(prefix_out, '_avgpos','.jpg');
cur_fig = figure('Visible','off','OuterPosition',[ 1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/1.1]); % figure('Visible','off');
colormap('Gray');

avg_green =average_pos(green_traces);
avg_red =average_pos(red_traces);
avg_fret =average_pos(red_traces);
    
subplot(2,2,1) %D_ex D_em
imagesc(green_avg_frame), title('D_{ex} -> D_{em}'), colorbar, axis([0 size(green_avg_frame,1) 0 size(green_avg_frame,1)]), hold on
plot(avg_green(:,2)+1 , avg_green(:,3)+1, 'go'), hold on
%ellipse(avg_green(:,7), avg_green(:,8),zeros(size(avg_green,1)),avg_green(:,2)+1, avg_green(:,3)+1,'g');
    
subplot(2,2,2) %D_ex A_em
imagesc(fret_avg_frame), title('D_{ex} -> A_{em}'), colorbar, axis([0 size(fret_avg_frame,1) 0 size(fret_avg_frame,1)]), hold on
plot(avg_fret(:,2)+1 , avg_fret(:,3)+1, 'bo'), hold on
%ellipse(avg_fret(:,5), avg_fret(:,5),zeros(size(avg_fret,1)),avg_fret(:,2)+1, avg_fret(:,3)+1,'b');
    
%subplot(2,2,3) %A_ex D_em
%imagesc(green_movie(:,:,2)), title('A_{ex} -> D_{em}'), colorbar, axis([0 size(green_movie,1) 0 size(green_movie,1)]), hold on

subplot(2,2,4) %A_ex A_em
imagesc(red_avg_frame), title('A_{ex} -> A_{em}'), colorbar, axis([0 size(red_avg_frame,1) 0 size(red_avg_frame,1)]), hold on
plot(avg_red(:,2)+1 , avg_red(:,3)+1, 'ro'), hold on
%ellipse(avg_red(:,5), avg_red(:,5),zeros(size(avg_red,1)),avg_red(:,2)+1, avg_red(:,3)+1,'r');
      
print(cur_fig, '-djpeg' , '-r800',fname); %save figure
    
%figure withe traces
fname = strcat(prefix_out, '_traces','.jpg');
cur_fig = figure('Visible','off','OuterPosition',[ 1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/1.1]); % figure('Visible','off');
colormap('Gray');

    
subplot(2,2,1) %D_ex D_em
imagesc(green_avg_frame), title('D_{ex} -> D_{em}'), colorbar, axis([0 size(green_avg_frame,1) 0 size(green_avg_frame,1)]), hold all
ColorSet = varycolor(size(green_traces,1) );
for i=1:size(green_traces)
   plot(green_traces{i}(:,2)+1, green_traces{i}(:,3)+1, 'Color', ColorSet(i,:));
end

    
subplot(2,2,2) %D_ex A_em
imagesc(fret_avg_frame), title('D_{ex} -> A_{em}'), colorbar, axis([0 size(fret_avg_frame,1) 0 size(fret_avg_frame,1)]), hold all
ColorSet = varycolor(size(red_traces,1) );
for i=1:size(red_traces)
   plot(red_traces{i}(:,2)+1, red_traces{i}(:,3)+1, 'Color', ColorSet(i,:));
end 

%subplot(2,2,3) %A_ex D_em
%imagesc(green_movie(:,:,2)), title('A_{ex} -> D_{em}'), colorbar, axis([0 size(green_movie,1) 0 size(green_movie,1)]), hold all

subplot(2,2,4) %A_ex A_em
imagesc(red_avg_frame), title('A_{ex} -> A_{em}'), colorbar, axis([0 size(red_avg_frame,1) 0 size(red_avg_frame,1)]), hold all
ColorSet = varycolor(size(red_traces,1) );
for i=1:size(red_traces)
   plot(red_traces{i}(:,2)+1, red_traces{i}(:,3)+1, 'Color', ColorSet(i,:));
end
      
print(cur_fig, '-djpeg' , '-r800',fname); %save figure







%figure with avergage position only merged_traces
fname = strcat(prefix_out, '_avgpos_merged','.jpg');
cur_fig = figure('Visible','off','OuterPosition',[ 1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/1.1]); % figure('Visible','off');
colormap('Gray');

    subplot(2,2,1) %D_ex D_em
    imagesc(green_avg_frame), title('D_{ex} -> D_{em}'), colorbar, axis([0 size(green_avg_frame,1) 0 size(green_avg_frame,1)]), hold on
    for i=1:size(merged_traces,1)
        plot(mean(merged_traces{i,1}(:,2))+1 , mean(merged_traces{i,1}(:,3))+1 , 'go'), hold on
    end
    
    subplot(2,2,2) %D_ex A_em
    imagesc(fret_avg_frame), title('D_{ex} -> A_{em}'), colorbar, axis([0 size(fret_avg_frame,1) 0 size(fret_avg_frame,1)]), hold on
    for i=1:size(merged_traces,1)
        plot(mean(merged_traces{i,3}(:,2))+1 , mean(merged_traces{i,3}(:,3))+1 , 'bo'), hold on
    end

    subplot(2,2,4) %A_ex A_em
    imagesc(red_avg_frame), title('A_{ex} -> A_{em}'), colorbar, axis([0 size(red_avg_frame,1) 0 size(red_avg_frame,1)]), hold on
    for i=1:size(merged_traces,1)
        plot(mean(merged_traces{i,2}(:,2))+1 , mean(merged_traces{i,2}(:,3))+1 , 'ro'), hold on
    end
    
print(cur_fig, '-djpeg' , '-r800',fname); %save figure




%plot merged traces at red cam
cur_fig = figure('Visible','off','OuterPosition',[ 1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/1.1]); % figure('Visible','off');
fname = strcat(prefix_out, '_merged','.jpg');
hold all
x_g = zeros(size(green_traces,1),1);
y_g = zeros(size(green_traces,1),1);
for i=1:size(green_traces)
   x_g(i)=mean(green_traces{i}(:,4));
   y_g(i)=mean(green_traces{i}(:,5));
end


x_r = zeros(size(red_traces,1),1);
y_r = zeros(size(red_traces,1),1);
for i=1:size(red_traces)
   x_r(i)=mean(red_traces{i}(:,2));
   y_r(i)=mean(red_traces{i}(:,3));
end

x_f = zeros(size(red_traces,1),1);
y_f = zeros(size(red_traces,1),1);
for i=1:size(red_traces)
   x_f(i)=mean(red_traces{i}(:,2));
   y_f(i)=mean(red_traces{i}(:,3));
end

x_m = zeros(size(merged_traces,1),1);
y_m = zeros(size(merged_traces,1),1);

for i=1:size(merged_traces, 1)
    x_m(i) = (mean(merged_traces{i,1}(:,4)) + mean(merged_traces{i,2}(:,2)) + mean(merged_traces{i,3}(:,2)))/3;
    y_m(i) = (mean(merged_traces{i,1}(:,5)) + mean(merged_traces{i,2}(:,3)) + mean(merged_traces{i,3}(:,3)))/3;
end

hold all
h = zeros(1,3);
h(1)=plot(x_g,y_g, 'gx'); 
h(2)=plot(x_r,y_r, 'r+');
h(3)=plot(x_f,y_f, 'b.');
ellipse(r_integrate.*ones(size(x_m,1),1),r_integrate.*ones(size(x_m,1),1), zeros(size(x_m,1),1), x_m , y_m, 'k');
axis([0 256*q_bin 0 256*q_bin ]), title('Average position of traces (combined are circled)'), legend(h,'donor','acceptor','fret')
print(cur_fig, '-djpeg' , '-r800',fname); %save figure






%{
%% plot each trace
for i=1:size(merged_itraces,1)
    cur_fig = figure('Visible','off','OuterPosition',[ 1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/1.1]); % figure('Visible','off');

    nbins = 20;
    subplot(2,3,[1:2]) 
    plot(merged_itraces{i,1}(:,1), merged_itraces{i,1}(:,4), 'g', merged_itraces{i,2}(:,1), merged_itraces{i,2}(:,4), 'r', merged_itraces{i,3}(:,1), merged_itraces{i,3}(:,4), 'b')
    legend('I_D', 'I_A', 'I_{D->A}')
    xlabel('Frame')
    ylabel('Intensity [a.u.]')
    title(['Trace ' num2str(i) ' of ' num2str(size(merged_itraces,1))])
    
    ylim = get(gca, 'Ylim');
    x = ylim(1):(ylim(2)-ylim(1))/nbins:ylim(2);
    
    subplot(2,3,3)
    [n_green x_green] = hist(merged_itraces{i,1}(:,4), x);
    [n_red x_red] = hist(merged_itraces{i,2}(:,4), x);
    [n_fret x_fret] = hist(merged_itraces{i,3}(:,4), x);

    barh(x_green, n_green, 'g'), hold on
    barh(x_red, n_red, 'r'), hold on
    barh(x_fret, n_fret, 'b'), axis([0 max( [max(n_green) max(n_red) max(n_fret)] ) ylim(1) ylim(2) ])
    
    
    colormap('gray')
    subplot(2,3,4)
    imagesc(green_avg_frame), colorbar, hold on
    plot(mean(merged_itraces{i,1}(:,2))+1 , mean(merged_itraces{i,1}(:,3))+1, 'go' )
    title('Average D_{Ex}->D_{Em}')
    
    subplot(2,3,5)
    imagesc(red_avg_frame), colorbar, hold on
    plot(mean(merged_itraces{i,2}(:,2))+1 , mean(merged_itraces{i,2}(:,3))+1, 'ro' )
    title('Average A_{Ex}->A_{Em}')
    
    
    subplot(2,3,6)
    imagesc(fret_avg_frame), colorbar, hold on
    plot(mean(merged_itraces{i,3}(:,2))+1 , mean(merged_itraces{i,3}(:,3))+1, 'bo' )
    title('Average D_{Ex}->A_{Em}')
    
    fname = strcat(prefix_out, '_', num2str(i) ,'.jpg');
    print(cur_fig, '-djpeg' , '-r300',fname); %save figure

    
end
%}

%% plot each trace
if plot_curves
cur_fig = figure('Visible','off','OuterPosition',[ 1 scrsz(4) scrsz(3)*0.8 scrsz(4)/2], 'PaperPositionMode', 'auto'); % figure('Visible','off');%left bottom width height

for i=1:size(merged_itraces,1)

    nbins = 2*sqrt( size(merged_itraces{i}, 1) );
    subplot(1,6,1:3) 
    plot(merged_itraces{i,1}(:,1), merged_itraces{i,1}(:,4), 'g', merged_itraces{i,2}(:,1), merged_itraces{i,2}(:,4), 'r', merged_itraces{i,3}(:,1), merged_itraces{i,3}(:,4), 'b')
    legend('I_D', 'I_A', 'I_{D->A}')
    xlabel('Frame')
    ylabel('Intensity [a.u.]')
    title(['Trace ' num2str(i) ' of ' num2str(size(merged_itraces,1))])
    
    ylim = get(gca, 'Ylim');
    x = ylim(1):(ylim(2)-ylim(1))/nbins:ylim(2);
    
    subplot(1,6,4)
    [n_green x_green] = hist(merged_itraces{i,1}(:,4), x);
    barh(x_green, n_green, 'g'), axis([0 max(n_green) ylim(1) ylim(2) ])
    set(gca,'YTick',[])
    set(gca,'XTick',[])

    
    subplot(1,6,5)
    [n_red x_red] = hist(merged_itraces{i,2}(:,4), x);
    barh(x_red, n_red, 'r'), axis([0 max(n_red) ylim(1) ylim(2) ])
    set(gca,'YTick',[])
    set(gca,'XTick',[])
    
    subplot(1,6,6)
    [n_fret x_fret] = hist(merged_itraces{i,3}(:,4), x);
    barh(x_fret, n_fret, 'b'), axis([0 max(n_fret) ylim(1) ylim(2) ])
    set(gca,'YTick',[])
    set(gca,'XTick',[])
    
    fname = strcat(prefix_out, '_', num2str(i) ,'.jpg');
    print(cur_fig, '-djpeg' , '-r100',[path_out filesep fname]); %save figure
   
    
end
end
%%


cd(path0)
display('done creating figures!')





%% make avi movi


if make_avi
   
close all
cd(path_out)
display('Avi-processing started...')
scrsz = get(0,'ScreenSize');
%set(gcf,'OuterPosition',[0 0 scrsz(3)/1 scrsz(4)/1.25],'Color',[1 1 1]),hold on
        %    set(gca,'Clim',[0 z_high])
 %figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)])
 
 frames(1:size(merged_itraces,1))=struct();
 i=6;
%for i=1:size(merged_itraces,1) %loops over traces
    display(strcat('writing trace  ', num2str(i), ' of ', num2str(size(merged_itraces,1))));
    vidobj=VideoWriter('biomod_vid.avi');
    %aviobj.FrameRate = 10; 
    
    vidobj.FrameRate = 8;
    

    open(vidobj);
    
    
    
    
    w2 = 24;
    r2 = 12;
    %aviobj = avifile('example.avi','compression','None');% avifile(strcat('testmovie_',num2str(i),'.avi')); 
    
    n=min(size(merged_itraces{i,1},1) , min(size(merged_itraces{i,2},1), size(merged_itraces{i,3},1)));
    f_g =merged_itraces{i,1}(:,1); 
    x_g = merged_itraces{i,1}(:,2)+1;
    y_g = merged_itraces{i,1}(:,3)+1;
    x_g_mean= mean(x_g);
    y_g_mean = mean(y_g);
   
    %frame for the 
    x_g_min=int16(max(x_g_mean-w2, 0));
    x_g_max=int16(min(x_g_mean+w2, q_bin*256));
    y_g_min=int16(max(y_g_mean-w2, 0));
    y_g_max=int16(min(y_g_mean+w2, q_bin*256));
    
    f_r =merged_itraces{i,2}(:,1); 
    x_r = merged_itraces{i,2}(:,2)+1;
    y_r = merged_itraces{i,2}(:,3)+1;
    x_r_mean= mean(x_r);
    y_r_mean = mean(y_r);
    x_r_min=int16(max(x_r_mean-w2, 0));
    x_r_max=int16(min(x_r_mean+w2, q_bin*256));
    y_r_min=int16(max(y_r_mean-w2, 0));
    y_r_max=int16(min(y_r_mean+w2, q_bin*256));
    
    f_f =merged_itraces{i,3}(:,1); 
    x_f = merged_itraces{i,3}(:,2)+1;
    y_f = merged_itraces{i,3}(:,3)+1;
    x_f_mean = mean(x_f);
    y_f_mean = mean(y_f);
    x_f_min=int16(max(x_f_mean-w2, 0));
    x_f_max=int16(min(x_f_mean+w2, q_bin*256));
    y_f_min=int16(max(y_f_mean-w2, 0));
    y_f_max=int16(min(y_f_mean+w2, q_bin*256));
    
    
    hhh = scrsz(4);
    hf = figure('Visible','on','OuterPosition',[0 0 hhh/1.33 hhh],'PaperPositionMode', 'auto'); % figure('Visible','off');
    colormap('Gray');
    %set(gca,'nextplot','replacechildren');
   
    green_counter =1;
    red_counter = 1;
    fret_counter = 1;
    for j=1:min(200,n) %loops over frames
        if mod(j,10)==0
            j
        end
       
        
        %green trace
        subplot(3, 3, 1)
        imagesc(green_movie(:,:,f_g(j)+1), [0 1500]  ), title(sprintf('Donor channel /\n Donor excitation'), 'Fontsize', 24), hold on
        ellipse(r2, r2, 0, x_g(j), y_g(j), 'g'); %circle around 
        
        subplot(3, 3, 4)
        %imagesc([x_g_min x_g_max],[y_g_min y_g_max],green_movie(y_g_min:y_g_max,x_g_min:x_g_max,f_g(j)+1)  ), title(sprintf('Donor channel / Donoer excitation \n%i , %i', x_g(j), y_g(j)), 'Fontsize', 24), hold on
        imagesc([x_g_min x_g_max],[y_g_min y_g_max],green_movie(y_g_min:y_g_max,x_g_min:x_g_max,f_g(j)+1)  ), hold on%, title(sprintf('Donor channel / Donoer excitation'), 'Fontsize', 24), hold on

        plot(x_g(j), y_g(j), 'gx'); %peak
        plot(x_g_mean, y_g_mean, 'g.' ); %average
        ellipse(r_integrate, r_integrate, 0, x_g(j), y_g(j), 'g'); %circle around 
        
        cur_framenumber = merged_itraces{i,1}(j,1);
        if cur_framenumber == merged_traces{i,1}(green_counter,1)
            plot(x_g(j), y_g(j), 'go');
            if green_counter < size(merged_traces{i,1}, 1)
                green_counter = green_counter +1;
            end
            
        end
        
        %subplot(2, 3, 4)
        %plot(merged_itraces{i,1}(1:j,1), merged_itraces{i,1}(1:j,4) ,'g'), axis([0 500 0 1.2*max(merged_itraces{i,1}(:,4))]);
        
        %red
        
        subplot(3, 3, 2)
        imagesc(red_movie(:,:,f_r(j)+1 ), [0 800]), title(sprintf('Acceptor channel /\n Acceptor excitation') , 'Fontsize', 24),hold on
        ellipse(r2, r2, 0,x_r(j), y_r(j), 'r');

        
        subplot(3, 3, 5)
        %imagesc([x_r_min x_r_max],[y_r_min y_r_max],red_movie(y_r_min:y_r_max,x_r_min:x_r_max,f_r(j)+1)  ), title(sprintf('Acceptor channel / Acceptor excitation \n%i , %i', x_r(j), y_r(j)) , 'Fontsize', 24),hold on
        imagesc([x_r_min x_r_max],[y_r_min y_r_max],red_movie(y_r_min:y_r_max,x_r_min:x_r_max,f_r(j)+1)  ), hold on %, title(sprintf('Acceptor channel / Acceptor excitation') , 'Fontsize', 24),hold on
        plot(x_r_mean, y_r_mean, 'r.');
        plot(x_r(j), y_r(j), 'rx');
        ellipse(r_integrate, r_integrate, 0,x_r(j), y_r(j), 'r');
        cur_framenumber = merged_itraces{i,2}(j,1);
        if cur_framenumber == merged_traces{i,2}(red_counter,1)
            plot(x_r(j), y_r(j),'ro');
            if red_counter < size(merged_traces{i,2}, 1)
                red_counter = red_counter +1;
            end
        end
        %subplot(3, 3, 6)
        %plot(merged_itraces{i,2}(1:j,1), merged_itraces{i,2}(1:j,4),'r' ), axis([0 500 0 1.2*max(merged_itraces{i,2}(:,4))]);
        
        %fret
        subplot(3, 3, 3)
        imagesc(red_movie(:,:,f_f(j)+1)  ),title(sprintf('Acceptor channel /\n Donor excitation') , 'Fontsize', 24), hold on
        ellipse(r2, r2, 0,x_f(j), y_f(j), 'b');

        subplot(3, 3, 6)
        %imagesc([x_f_min x_f_max],[y_f_min y_f_max],red_movie(y_f_min:y_f_max,x_f_min:x_f_max,f_f(j)+1)  ),title(sprintf('Acceptor channel / Donor excitation \n%i , %i', x_f(j), y_f(j)) , 'Fontsize', 24), hold on
        imagesc([x_f_min x_f_max],[y_f_min y_f_max],red_movie(y_f_min:y_f_max,x_f_min:x_f_max,f_f(j)+1)  ), hold on %title(sprintf('Acceptor channel / Donor excitation') , 'Fontsize', 24), hold on
        plot(x_f_mean, y_f_mean, 'b.');
        plot(x_f(j), y_f(j), 'bx');
        ellipse(r_integrate, r_integrate, 0,x_f(j), y_f(j), 'b');
        cur_framenumber = merged_itraces{i,3}(j,1);
        if cur_framenumber == merged_traces{i,3}(fret_counter,1)
            plot(x_f(j), y_f(j), 'bo');
            if fret_counter < size(merged_traces{i,3}, 1)
                fret_counter = fret_counter +1;
            end
            
        end
        %subplot(2, 3, 6)
        %plot(merged_itraces{i,3}(1:j,1), merged_itraces{i,3}(1:j,4), 'b' ), axis([0 500 0 1.2*max(merged_itraces{i,3}(:,4))]);
        
        
        subplot(3,3, [7:9])
        frate = 16;
        plot(merged_itraces{i,1}(1:j,1)/frate, merged_itraces{i,1}(1:j,4) ,'g'), hold on
        plot(merged_itraces{i,2}(1:j,1)/frate, merged_itraces{i,2}(1:j,4),'r' ), hold on
        plot(merged_itraces{i,3}(1:j,1)/frate, merged_itraces{i,3}(1:j,4), 'b' ), axis([0 400/frate 0 40000]);
        xlabel('Time [s]', 'Fontsize', 20)
        ylabel('Intensity [a.u.]', 'Fontsize', 20)

        
        
        
        
        %frames(i) =getframe(hf);
        %movie2avi(mov, strcat(prefix_out, '_',num2str(i),'_movie'));        
        writeVideo(vidobj, getframe(hf));
        
        %aviobj=addframe(aviobj,hf);
    end 
    %aviobj=close(aviobj); %closes the AVI file  
    close(vidobj); %closes the handle to invisible figure 
    
    %aviobj = close(aviobj)  ;
    
%end
display('Avi-processing done')
cd(path0)
end %end avi movies

close all


display('--------------------------- DONE ----------------------------')
display(['This was movie ' sprintf('%02s', Im_name1G(6:7) )])
%%












