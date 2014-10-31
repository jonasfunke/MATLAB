%%
clear all, close all, clc
run('my_prefs')
path0 = cd;

%% load image
%nimg_input = inputdlg({'Number of images:'}, 'Number of images' , 1, {'1'} );
n_img = 1; %str2double(nimg_input{1});

filenames = cell(n_img, 1);
pathnames = cell(n_img, 1);

cd(data_dir)
[filenames{1} pathnames{1}]=uigetfile('*.tif','Select image:');
cd(path0)

prefix_out = [filenames{1}(1:size(filenames{1},2)-4) '_bands01'];



%% create output folder
pname = inputdlg({'Output folder and prefix:'}, 'Output folder and prefix' , 1, {prefix_out} );
prefix_out = pname{1};
path_out = [pathnames{1} prefix_out ];
mkdir(path_out);

%% load and bg correct images
images = cell(n_img, 1);
img_bg = cell(n_img, 1);
for i=1:n_img
    images{i} = double(imread([pathnames{i} filesep filenames{i}]));  %load
    disp(['Loaded image: ' filenames{i}])
    disp(['Image directory: ' pathnames{i}])
    plot_image_ui(images{i});
    button = questdlg('Rotate?','Rotate','Rotate','No','No');
    if strcmp(button,'Rotate') %load old data
        images{i} = imrotate(images{i}, -90);
    end
    close all
    img_bg{i} = bg_correct_ui(images{i}, 'Background correction', 4);    %bg correct
    close all    
end

%% select bands by hand
imagesc(img_bg{1}), colormap gray
options.WindowStyle='normal';
prompt={'How many bands'};
def={'3'};
tmp = inputdlg(prompt, 'How many bands', 1, def, options);
n_bands = str2double(tmp(1));
close all

%% integrate areas, all areas have the same size
[I, areas] = integrate_areas(img_bg, n_bands, 1, [1 1]); %cell of images, number of bands, 1=all bands habe the same size

%{
% sum areas
for i=1:n_bands
    I(i,1) = sum(sum(img_bg{1}(areas(i,2):areas(i,2)+areas(i,4), areas(i,1):areas(i,1)+areas(i,3) )));
    I(i,2) = sum(sum(img_bg{2}(areas(i,2):areas(i,2)+areas(i,4), areas(i,1):areas(i,1)+areas(i,3) )));
    
end
%}
%% assining names
for i=1:n_bands
    names{i} = ['Band ' num2str(i)];
end


%%
%colormap('gray'), colorbar

%imagesc(img_bg)
for i=1:n_img
    cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual');
    colormap('gray'), colorbar

    imagesc(img_bg{i}, prctile(img_bg{i}(:), [20 99]) ), colorbar, axis image, colormap gray,
    title(filenames{i})
    for j=1:n_bands
       rectangle('Position', areas(j,:), 'EdgeColor', 'r')
       text(double(areas(j,1)), double(areas(j,2)), num2str(j), 'Fontsize', 8, 'Color' , 'r', 'HorizontalAlignment','left', 'VerticalAlignment', 'top')
    end 
    %plot_image_ui(img_bg{i})

    
   
    print(cur_fig, '-dtiff','-r500' , [path_out filesep prefix_out '_bands_' num2str(i) '.tif']); %save figure
    
    close all
end



%%
I_dimer = I(1:2:end,1);

I_monomer = I(2:2:end,1); % change this depending whether short or long oligo ist cy3 or cy5

yield = I_dimer ./ (I_dimer + I_monomer);



%%

save([path_out filesep prefix_out '_data'] )

dlmwrite([path_out filesep 'data.txt' ], I, 'delimiter', '\t')

disp('Data saved.')

%% PLOTTING STUFF
close all
fig_dim =[20 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

bar(1:n_bands, I)
set(gca, 'XTick',1:n_bands, 'XTickLabel',names,'Fontsize', 12)
xticklabel_rotate([1:n_bands],90,names)
ylabel('Intensity')
%legend({'GFP', 'DNA'})
print(cur_fig, '-dtiff','-r500' , [path_out filesep prefix_out '_IntesityRaw.tif']); %save figure


%% plot absolute intensity
x = 1:n_bands/2;

close all
fig_dim =[20 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

T = I_monomer+I_dimer;

plot( x, T, 'r.-')
ylabel('Sum of dimer and monomer')
xlabel('Lane')
set(gca, 'XLim', [1 n_bands/2], 'XTick', [1:n_bands/2], 'YLim', [min([0 min(T(:))])   1.1*max(T(:))]);

legend({'cy5-channel'})

print(cur_fig, '-dtiff','-r500' , [path_out filesep prefix_out '_sum.tif']); %save figure




%% plot yield against lane
x = 1:n_bands/2;

close all
fig_dim =[20 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

plot( x, 100*yield(:,1), 'r.-')

xlabel('Lane')
set(gca, 'XLim', [1 n_bands/2], 'XTick', [1:n_bands/2], 'YLim', [-10 80], 'XLim', [0.5 n_bands/2+0.5] );
ylabel('Yield [%]')
legend({ 'cy5-channel'}, 'location', 'NorthWest')

print(cur_fig, '-dtiff','-r500' , [path_out filesep prefix_out '_yield.tif']); %save figure


%% plot yield against concentration
%{
x = [ 0 25 50 75 100 125 150 175 200 250 500 1e3 2e6];
i = 1:13;
close all
fig_dim =[20 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

%semilogx(x(1:end-2), 100.*yield(1:end-2,1), 'g.-', x(1:end-2), 100.*yield(1:end-2,2), 'r.-', x(end-1:end), 100.*yield(end-1:end,1), 'g.', x(end-1:end), 100.*yield(end-1:end,2), 'r.')
semilogx(x(i), 100.*yield(i,1), 'g.-', x(i), 100.*yield(i,2), 'r.-'), hold on

xlabel('Concentration of BMH-Linker [nM]')
set(gca, 'XLim', [min(x) max(x)],  'YLim', [-10 50]);
ylabel('Yield [%]')
legend({'cy3-channel', 'cy5-channel'}, 'location', 'NorthWest')
print(cur_fig, '-dtiff','-r500' , [path_out filesep prefix_out '_yield_vs_conc.tif']); %save figure

%%
close all
fig_dim =[20 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

plot(x(i), 100.*yield(i,1), 'g.-', x(i), 100.*yield(i,2), 'r.-')

xlabel('Concentration of BMH-Linker [nM]')
%set(gca, 'XLim', [min(x) max(x)]);
set(gca, 'XLim', [0 1000], 'YLim', [-10 20]);

ylabel('Yield [%]')
legend({'cy3-channel', 'cy5-channel'}, 'location', 'NorthWest')

print(cur_fig, '-dtiff','-r500' , [path_out filesep prefix_out '_yield_vs_conc_Lin_02.tif']); %save figure
%}


%% plot yield against time
%{

t = ([0 2 7 27 48 73 142 14*24].*60 + [0 30 50 30 20 30 30 0])./60;
t2 = ([0 2 7 27 48 73 142 ].*60 + [0 30 50 30 20 30 30 ])./60;
i = 1:8;
i = 9:15;

close all
fig_dim =[20 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);


subplot(1, 2, 1)
plot(t, 100*yield(i,2), 'r.-', t, 100*yield(i,1), 'g.-')
xlabel('Time of reaction [h]')
set(gca, 'XLim', [t(1) t(end)]);
ylabel('Yield [%]')
legend({'cy5-channel', 'cy3-channnel'}, 'location', 'NorthWest')
title('100 nM')

subplot(1, 2, 2)
semilogx(t, 100*yield(i,2), 'r.-', t, 100*yield(i,1), 'g.-')
xlabel('Time of reaction [h]')
set(gca, 'XLim', [t(1) t(end)]);
ylabel('Yield [%]')
legend({'cy5-channel', 'cy3-channnel'}, 'location', 'NorthWest')
print(cur_fig, '-dtiff','-r500' , [path_out filesep prefix_out '_timescale_50uM.tif']); %save figure

%%


close all
fig_dim =[20 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

i=1:8;
plot(t, 100*yield(i,2), 'r.-'), hold on
i=9:15;

plot(t2, 100*yield(i,2), 'b.-'), hold on
xlabel('Time of reaction [h]')
set(gca, 'XLim', [t(1) t(end)]);
ylabel('Yield [%]')
legend({'100 nM', '500 nM'}, 'location', 'NorthWest')
title('Comparison of 100 nM and 500 nM in cy5-channel')


print(cur_fig, '-dtiff','-r500' , [path_out filesep prefix_out '_timescale_combined.tif']); %save figure
%%

close all
fig_dim =[20 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

i=1:8;
semilogx(t, 100*yield(i,2), 'r.-'), hold on
i=9:16;
semilogx(t, 100*yield(i,2), 'b.-'), hold on
xlabel('Time of reaction [h]')
set(gca, 'XLim', [t(1) t(end)]);
ylabel('Yield [%]')
legend({'100 nM', '500 nM'}, 'location', 'NorthWest')
title('Comparison of 100 nM and 500 nM in cy5-channel')


print(cur_fig, '-dtiff','-r500' , [path_out filesep prefix_out '_timescale_combined_log.tif']); %save figure
%}

%% plot yield against concentration
%{
t = [20 100 135 246 480 2940];

close all
fig_dim =[20 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

semilogx(t, 100*yield(1:6,2), 'r.-', t, 100*yield(7:12,2), 'g.-')

xlabel('Time of reaction [min]')
set(gca, 'XLim', [t(1) t(end)]);
ylabel('Yield [%]')
legend({'500 nM BMH-Linker', '2 mM BMH-Linker'}, 'location', 'NorthWest')

print(cur_fig, '-dtiff','-r500' , [path_out filesep prefix_out '_timescale2.tif']); %save figure

%}

%%
disp('Finished.')