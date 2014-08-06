%%
clear all
close all
clc
run('my_prefs')
path0 = cd;

%% load image
%nimg_input = inputdlg({'Number of images:'}, 'NUmber of images' , 1, {'1'} );
n_img = 1;%str2double(nimg_input{1});

filenames = cell(n_img, 1);
pathnames = cell(n_img, 1);

last_dir = data_dir;
for i=1:n_img
    cd(last_dir)
    if i==1
        [filenames{i} pathnames{i}]=uigetfile('*.tif','Select image:');
    end
    last_dir = pathnames{i};
end
cd(path0)

%% create output folder
pname = inputdlg({'Output folder and prefix:'}, 'Output folder and prefix' , 1, {[filenames{1}(1:size(filenames{1},2)-4) '_analysis']}  );
prefix_out = pname{1};
path_out = [pathnames{1} prefix_out ];
mkdir(path_out);

%% load and bg correct images
images = cell(n_img, 1);
img_bg = cell(n_img, 1);
for i=1:n_img
    images{i} = double(imread([pathnames{i} filesep filenames{i}]));  %load
    plot_image_ui(images{i})
    button = questdlg('Rotate?','Rotate','Rotate','No','No');
    if strcmp(button,'Rotate') %load old data
        images{i} = imrotate(images{i}, -90);
    end
    close all
    img_bg{i} = bg_correct_ui(images{i}, 'Background correction');    %bg correct
    close all    
end

%% select bands by hand
imagesc(img_bg{1}), colormap gray
options.WindowStyle='normal';
prompt={'How many samples to evaluate'};
def={'1'};
tmp = inputdlg(prompt, 'How many samples to evaluate', 1, def, options);
n_bands = 2*str2double(tmp(1));
close all

I = zeros(n_bands, n_img);
I_max = zeros(n_bands, n_img);
areas = zeros(n_bands, 4);
for i=1:n_bands
    if mod(i,2)==1
        display('Select un-purified band')
    else
        display('Select purified band')
    end
    
    if i==1
        plot_image_ui(img_bg{1})
        h = imrect;
        position = wait(h);
        pos_init = int32(getPosition(h)) % [xmin ymin width height]
        pos = pos_init;
    else
        plot_image_ui(img_bg{1})
        h = imrect(gca, double(pos_last));
        setResizable(h, 0)
        position = wait(h);
        pos = int32(getPosition(h)) ;% [xmin ymin width height]  
    end
    pos_last = pos;
    areas(i,:) = pos;
    for j = 1:n_img
        I(i,j) = sum(sum((img_bg{j}(pos(2):pos(2)+pos(4)  , pos(1):pos(1)+pos(3))))); %integrate
        I_max(i,j) = max(max(img_bg{j}( pos(2):pos(2)+pos(4)  , pos(1):pos(1)+pos(3)  )));

    end
    close all
end
%% assining names
names = cell(n_bands,1);
names2 = cell(n_bands/2,1);
for i=1:2:n_bands
    band_name = inputdlg({'Name of sample:'}, 'Sample names' , 1, {['Sample ' num2str((i+1)/2)]} );
    names{i} = [band_name{1} ' un-purified'];
    names{i+1} = [band_name{1} ' purified'];
    names2{(i+1)/2} = band_name{1} ;
end


%%
%colormap('gray'), colorbar

%imagesc(img_bg)
for i=1:n_img
    cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual');
    colormap('gray'), colorbar

    plot_image(img_bg{i},filenames{i}, 12 )
    for j=1:n_bands
       rectangle('Position', areas(j,:), 'EdgeColor', 'r')
       text(double(areas(j,1)), double(areas(j,2)), num2str(j), 'Fontsize', 8, 'Color' , 'r', 'HorizontalAlignment','left', 'VerticalAlignment', 'top')
    end 
    %plot_image_ui(img_bg{i})

    
    print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_bands_' num2str(i) '.eps']); %save figure
    print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_bands_' num2str(i) '.png']); %save figure
    
    close all
end


%%
save([path_out filesep prefix_out '_data'] )



%%
%colormap('gray'), colorbar

%imagesc(img_bg)
for i=1:n_img
 %   cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual');
 %   colormap('gray'), colorbar

    %plot_image(img_bg{i},filenames{i}, 12 )
    plot_image_ui(img_bg{i})
    for j=1:n_bands
       rectangle('Position', areas(j,:), 'EdgeColor', 'r')
       text(double(areas(j,1)), double(areas(j,2)), num2str(j), 'Fontsize', 8, 'Color' , 'r', 'HorizontalAlignment','left', 'VerticalAlignment', 'top')
    end 
    %plot_image_ui(img_bg{i})

    
   % print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_bands_' num2str(i) '.eps']); %save figure
   % print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_bands_' num2str(i) '.png']); %save figure
    pause
    close all
end


%% PLOTTING STUFF
close all
fig_dim =[20 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

bar(I)
set(gca, 'XTick',1:n_bands, 'XTickLabel',names,'Fontsize', 12)
xticklabel_rotate([1:n_bands],90,names)

set(gca, 'XTick', 1:n_bands)
ylabel('Intensity')
%legend({'GFP', 'DNA'})
print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_IntesityRaw.eps']); %save figure
print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_IntesityRaw.png']); %save figure

%% PLOTTING STUFF
close all
fig_dim =[20 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

bar(I./I(end))
set(gca, 'XTick',1:n_bands, 'XTickLabel',names,'Fontsize', 12)
xticklabel_rotate([1:n_bands],90,names)

set(gca, 'XTick', 1:n_bands)
ylabel('Intensity normalized to @16 unpur')
%legend({'GFP', 'DNA'})
print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_IntesityRaw2.eps']); %save figure
print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_IntesityRaw2.png']); %save figure



%%
close all
fig_dim =[20 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);


bar(I(1:2:n_bands,:)./(I(1:2:n_bands,:)+I(2:2:n_bands,:)))
set(gca, 'XTick',1:n_bands/2, 'XTickLabel',names2,'Fontsize', 12)
xticklabel_rotate([1:n_bands/2],90,names2)

set(gca, 'XTick', 1:n_bands/2)
%set(gca, 'YLim', [0 1])
%ylabel('[ABA] / ([ABA] + [A] + [AB])')
%ylabel('[dimer] / ([dimer] + [monomer])')
ylabel('[GFP] / ([GFP] + [GFP+DNA])')

legend({'GFP-channel', 'DNA-channel'})

print(cur_fig, '-depsc2','-loose' , [path_out filesep prefix_out '_normalizedIntensities.eps']); %save figure
print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_normalizedIntensities.png']); %save figure



%% plot ratio
yield = I(2:2:n_bands,1) ./ I(1:2:n_bands,1); 
close all
fig_dim =[20 15];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

bar(1:n_bands/2, I(2:2:n_bands,1) ./ I(1:2:n_bands,1));
set(gca, 'XTick',1:n_bands/2, 'XTickLabel',names2,'Fontsize', 12)
xticklabel_rotate([1:n_bands/2],90,names2)

set(gca, 'XTick', 1:n_bands)
set(gca, 'XLim', [0 n_bands/2+1])
%set(gca, 'YLim', [0 2])
%ylabel('[ABA] / ([ABA] + [A] + [AB])')
ylabel('Purification yield: [purified] / [un-purified]')
title('Yield of purification')

for i=1:n_bands/2
    text(i,0.5, {[num2str(round(100*yield(i))) '%'], [num2str(round(20*yield(i)*10)./10) ' nM' ]}, 'Color', 'red', 'FontSize', 12, 'HorizontalAlignment','center')
end


print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_yield.png']); %save figure





%% plot ratio
yield = I(2:2:n_bands,1) ./ I(1:2:n_bands,1); 
close all
fig_dim =[20 10];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

subplot(1, 2, 1)
bar(1:n_bands, I(:,1) );
set(gca, 'XTick',1:n_bands, 'XTickLabel',names,'Fontsize', 12)
xticklabel_rotate([1:n_bands],90,names)
ylabel('Intensity [a.u.]')
title('Raw intensities')
%set(gca, 'XTick', 1:n_bands)
%set(gca, 'XLim', [0 n_bands+1])

subplot(1, 2, 2)
bar(1:n_bands/2, I(2:2:n_bands,1) ./ I(1:2:n_bands,1));
set(gca, 'XTick',1:n_bands/2, 'XTickLabel',names2,'Fontsize', 12)
xticklabel_rotate([1:n_bands/2],90,names2)

set(gca, 'XTick', 1:n_bands)
set(gca, 'XLim', [0 n_bands/2+1])
set(gca, 'YLim', [0 1])

ylabel('Purification yield: [purified] / [un-purified]')
title('Yield of purification')

for i=1:n_bands/2
    text(i,0.5, {[num2str(round(100*yield(i))) '%'], [num2str(round(20*yield(i)*10)./10) ' nM' ]}, 'Color', 'red', 'FontSize', 12, 'HorizontalAlignment','center')
end


print(cur_fig, '-dpng','-loose' , [path_out filesep prefix_out '_yield2.png']); %save figure






%%
display('DONE')


