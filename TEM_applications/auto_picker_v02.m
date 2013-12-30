%% startup
clear all; close all; clc;
run('my_prefs'); path0=cd;

%% parameters
box_size = 50; % size of particle, HAS TO BE EVEN
dalpha = 5; % deg, angle resolution
mirror = 1; % include mirror transformation, 0=no, 1= yes
img_size = 512; %size of the binned image
r_filter = 15; % pixel (original image), radius for gaussian high pass

%% load images
pname=uigetdir(data_dir,'Choose a folder with tem images.'); % get pathname
tmp = dir([pname filesep '*_16.TIF']);
fnames = {tmp.name}; % list of filenames
n_img = size(fnames,2);
path_out = [pname filesep 'particles']; % output folder
mkdir(path_out)

%% load images
disp(['Loading and filtering ' num2str(n_img) ' images...'])
images = zeros(img_size, img_size, n_img);
h = fspecial('gaussian', r_filter*4*2 , r_filter); % gaussian filter, diameter = 2*(width = 4*sigma)
for i=1:n_img
    img = imread([pname filesep fnames{i}], 'PixelRegion', {[1 2048], [1 2048]});
    tmp = double(img)-double(imfilter(img, h, 'same'));
    images(:,:,i) = imresize(tmp,[512 512], 'nearest'); %bin image 4x4 for faster image processing
    %images(:,:,i) = imresize(img(1:2048,1:2048),[512 512], 'nearest'); %bin image 4x4 for faster image processing
end
%% create class references
tmp = inputdlg({'Number of class references:'}, 'Number of classes', 1, {'2'});
n_ref = str2double(tmp(1))*(mirror+1);

box_size_template = ceil(2*box_size/2/cos(pi/4));%200;
box_size_template = int16(box_size_template + mod(box_size_template, 2) +1 ) ; % make it even
templates = zeros(box_size_template, box_size_template, n_ref*(mirror+1));
w = (box_size_template-1)/2;
for i=1:n_ref/(mirror+1) 
    go_on = 1;
    j = 1;
    close all
    fig_dim =1.0*[10 10];
    cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

    while go_on 
        imagesc(images(:,:,j)), axis image, colormap gray
        button = questdlg(['Select an image for ref. ' num2str(i)],'Image','Use this','Previous','Next', 'Use this');
        if strcmp(button, 'Next')
            j = min(n_img,j+1);
        end
        if strcmp(button, 'Previous')
            j = max(j-1, 1);
        end
        if strcmp(button, 'Use this')
            go_on = 0;
        end
    end
    % select particle
    h = imrect(gca, [img_size/2 img_size/2 double(box_size_template) double(box_size_template)]);
    setResizable(h,0) 
    pos = int16(wait(h));

    % refine reference
    close all
    template = images(pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3),j);
    imagesc( [pos(1) pos(1)+pos(3)], [pos(2) pos(2)+pos(4)],template), colorbar, colormap gray, axis image
    c = round(ginput(1));
    area = [c(2)-w c(2)+w c(1)-w c(1)+w];
    if mirror
        templates(:,:,2*i-1) = images(area(1):area(2), area(3):area(4), j);
        templates(:,:,2*i) = flipdim(images(area(1):area(2), area(3):area(4), j) ,1);
    else
        templates(:,:,i) = images(area(1):area(2), area(3):area(4), j);
    end
      
end
close all
%% display templates
dx = (box_size_template-box_size-1)/2;
close all
for i=1:n_ref
    subplot(n_ref/(mirror+1),mirror+1,i)
    imagesc(templates(dx:dx+box_size, dx:dx+box_size,i)),  colormap gray, axis image
end

%% Calculate X-Correlation and find maximum correlations
disp('Calculating x-correlation...')
alpha = 0:dalpha:359;
n_rot = length(alpha); % number of rotations
dx = (box_size_template-box_size-1)/2;

peaks = cell(n_ref, n_img);
peaks2 = cell(n_ref, n_img);

h = waitbar(0,'Calculating x-correlation... ? time remaining');

for t=1:n_ref
    % generate library
    lib = zeros(box_size+1, box_size+1, length(alpha));
    for j=1:n_rot
        tmp = imrotate(templates(:,:,t), alpha(j), 'crop');
        lib(:,:,j) = tmp(dx:dx+box_size, dx:dx+box_size);
    end
 
    %loop through images
    img3 = zeros(img_size, img_size, n_img); % stores maximum of cor-coef of all rotations
    img3_index = zeros(img_size, img_size, n_img); % stores index of maximum
    
    for i=1:n_img
        tic
        xcor_img = zeros(512, 512, n_rot);
        
        for r=1:n_rot % loop through rotations
            tmp = normxcorr2(lib(:,:,r), images(:,:,i)); % x-correlate
            xcor_img(:,:,r) = tmp(box_size/2+1:end-box_size/2, box_size/2+1:end-box_size/2);
        end
        
        
        for k=1:512
        for l=1:512
            [cmax, imax] = max(xcor_img(k,l,:));
            img3(k,l,i) = cmax;
            img3_index(k,l,i) = imax;
        end
        end
        
        
        %find peaks in img3
        tmp = img3(box_size/2+1:end-box_size/2, box_size/2+1:end-box_size/2,i);
        h_min = mean(tmp(:)) + 0.25*std(tmp(:));
        p = find_peaks2d(tmp, round(box_size/4), h_min, 0); % radius of window, minimal height,  no absolure = relative height
        p(:,1:2) = p(:,1:2)+box_size/2+1;
        
        tmp = img3(:,:,i);
        tmp_index = img3_index(:,:,i);
        idx = sub2ind(size(tmp), p(:,2), p(:,1) );
        peaks{t,i} = [p(:,1:2) tmp(idx) alpha(tmp_index(idx))']; % x y coer_coef alpha 
    
        display(['Reference (' num2str(t) '/' num2str(n_ref) '): image (' num2str(i) '/' num2str(n_img) '), found ' num2str(size(p,1)) ' particles' ])

        if t==1 && i==1
            dt = toc;
        else
            dt_this = toc;
            dt = (dt+dt_this)/2;
        end
        frac = ((t-1)*n_img+i) / (n_ref*n_img);
        n_remain = (n_ref*n_img)-((t-1)*n_img+i);
        waitbar( frac , h, ['Calculating x-correlation... ' num2str(round(n_remain*dt/60*10)/10) ' min remaining']) 
    end
end
close(h)
close all


%% remove particles, which belong to multiple classes
disp('Removig duplicates...')
cc = varycolor(n_ref);
peaks_ref = cell(n_ref, n_img);
for i=1:n_img
    % genertate image of  correlations
    cor_img = zeros(img_size,img_size );
    cor_img_index = zeros(img_size,img_size );
    rot_img = zeros(img_size,img_size );
    for r=1:n_ref
        idx = sub2ind(size(cor_img), peaks{r,i}(:,2), peaks{r,i}(:,1) );
        cor_img(idx) = peaks{r,i}(:,3);
        cor_img_index(idx) = r;
        rot_img(idx) = peaks{r,i}(:,4);

    end
    
    p = find_peaks2d(cor_img, round(box_size/4), 0, 1 ); % find-peaks, width, min_height, absolute height 
    p(:,1:2) =  p(:,1:2)+1;
   

    for j=1:size(p,1)
        peaks_ref{cor_img_index(p(j,2),p(j,1)),i} = [peaks_ref{cor_img_index(p(j,2),p(j,1)),i}; p(j,1:2) cor_img(p(j,2),p(j,1)) rot_img(p(j,2),p(j,1)) ];
    end
    
    %{
     close all
    imagesc(cor_img), colorbar,  colormap gray, axis image, hold on % img2

    for r=1:n_ref
        plot(peaks{r,i}(:,1), peaks{r,i}(:,2),  'o', 'color', cc(r,:))
        plot(peaks_ref{r,i}(:,1), peaks_ref{r,i}(:,2),  '.', 'color', cc(r,:));

    end
    pause 
    %}
end

%% view images and found particles
%{
cc = varycolor(n_ref);
close all
myleg = cell(n_ref,1);
for t=1:n_ref
    myleg{t} = ['Reference ' num2str(t)];
end
for i=1:n_img
    subplot(1, 2, 1)
    imagesc(img3(:,:,i)), colorbar,  colormap gray, axis image, hold on % img2
    for r=1:n_ref
        h(r) = plot(peaks_ref{r,i}(:,1), peaks_ref{r,i}(:,2),  '.', 'color', cc(r,:));
    end
    title(['Image ' num2str(i) '/' num2str(n_img) ])
    legend(h, myleg)
    hold off
    
    subplot(1, 2, 2)
    imagesc(images(:,:,i)), colorbar, colormap gray, axis image, hold on %
    for r=1:n_ref
        h(r)=plot(peaks_ref{r,i}(:,1), peaks_ref{r,i}(:,2),  '.', 'color', cc(r,:));
    end
    title(['Image ' num2str(i) '/' num2str(n_img) ])
    legend(h, myleg)
    hold off
    pause
end
%}

%% generate stack of particles
disp('Writing particles...')
close all
particles = cell(n_ref, n_img);
particles_rot = cell(n_ref, n_img);

for i=1:n_img
    img = imread([pname filesep fnames{i}], 'PixelRegion', {[1 2048], [1 2048]});
    N = 2048;
    w = 4*box_size/2;
    for t=1:n_ref
        p= zeros(2*w+1, 2*w+1, size(peaks_ref{t,i},1), 'uint16');
        p_rot= zeros(2*w+1, 2*w+1, size(peaks_ref{t,i},1), 'uint16');

        for j=1:size(peaks_ref{t,i},1)
            x = 4*peaks_ref{t,i}(j,1);
            y = 4*peaks_ref{t,i}(j,2);
            angle = peaks_ref{t,i}(j,4);
            p(max(1,w-y+2):min(2*w+1, 2*w+1-y-w+N),max(1,w-x+2):min(2*w+1, 2*w+1-x-w+N),j) = img(max(1, y-w):min(N,y+w) , max(1, x-w):min(N,x+w) );
            
                    lib(:,:,j) = tmp(dx:dx+box_size, dx:dx+box_size);
            
            w2=w+4*dx;
            tmp = imrotate(img(max(1, y-w2):min(N,y+w2) , max(1, x-w2):min(N,x+w2) ), -angle, 'crop');

            
            p_rot(max(1,w-y+2):min(2*w+1, 2*w+1-y-w+N),max(1,w-x+2):min(2*w+1, 2*w+1-x-w+N),j) = tmp(4*dx:4*dx+4*box_size, 4*dx:4*dx+4*box_size);
            %{
            subplot(3,2,1:4)
            imagesc(img), colorbar, colormap gray, axis image, hold on
            plot(x,y, 'r.')
            hold off
            subplot(3,2,5)
            imagesc(p(:,:,j)), colorbar, colormap gray, axis image
            title(['Particle ' num2str(j) ' angle = 0' ])
            subplot(3,2,6)
            imagesc(p_rot(:,:,j)), colorbar, colormap gray, axis image
            title(['Particle ' num2str(j) ' angle = ' num2str(angle)])
            pause
            %}
            
        end
        particles{t,i} = p;
        particles_rot{t,i} = p_rot;

    end
end

%% write particles for each reference 
w = 4*box_size/2;
for t=1:n_ref/(mirror+1)
    n_particle = 0;
    for i=1:n_img
        n_particle = n_particle + size(particles{(mirror+1)*t,i}, 3);
        if mirror
            n_particle = n_particle + size(particles{(mirror+1)*t-1,i}, 3);
        end
    end    
    disp(['Reference ' num2str(t) '/' num2str(n_ref/(mirror+1)) ': ' num2str(n_particle) ' particles'])
    p_out = zeros(2*w+1, 2*w+1, n_particle, 'uint16');

    m=1;
    if mirror 
        for i=1:n_img
            for j=1:size(particles{(mirror+1)*t-1,i},3)
                p_out(:,:,m) = particles{(mirror+1)*t-1,i}(:,:,j);
                m = m+1;
            end
        end
    end
    for i=1:n_img
        for j=1:size(particles{(mirror+1)*t,i},3)
            p_out(:,:,m) = particles{(mirror+1)*t,i}(:,:,j);
            m = m+1;
        end
    end
    
    % write as spider-file
    writeSPIDERfile([path_out filesep 'ref_' num2str(t) '.spi'], p_out, 'stack')
    
    % write as single tif-files
    path_out_tif = [path_out filesep 'ref_' num2str(t) '_tif'];
    mkdir(path_out_tif)
    for i=1:n_particle
        imwrite(p_out(:,:,i), [path_out_tif filesep 'ref_' num2str(t) '_' sprintf('%.3i',i) '.tif' ]);
    end
    
end

%% write rotated particles
w = 4*box_size/2;
for t=1:n_ref/(mirror+1)
    n_particle = 0;
    for i=1:n_img
        n_particle = n_particle + size(particles_rot{(mirror+1)*t,i}, 3);
        if mirror
            n_particle = n_particle + size(particles_rot{(mirror+1)*t-1,i}, 3);
        end
    end  
    
    p_out = zeros(2*w+1, 2*w+1, n_particle, 'uint16');
    
    m=1;
    if mirror 
        for i=1:n_img
            for j=1:size(particles_rot{(mirror+1)*t-1,i},3)
                p_out(:,:,m) = particles_rot{(mirror+1)*t-1,i}(:,:,j);
                m = m+1;
            end
        end
    end
    for i=1:n_img
        for j=1:size(particles_rot{(mirror+1)*t,i},3)
            if mirror
                p_out(:,:,m) = flipdim(particles_rot{(mirror+1)*t,i}(:,:,j), 1);
            else
                p_out(:,:,m) = particles_rot{(mirror+1)*t,i}(:,:,j);
            end
            m = m+1;
        end
    end
    
    % write as spider-file
    writeSPIDERfile([path_out filesep 'ref_' num2str(t) '_rot.spi'], p_out, 'stack')
    
    % write as single tif-files
    path_out_tif = [path_out filesep 'ref_' num2str(t) '_rot_tif'];
    mkdir(path_out_tif)
    for i=1:n_particle
        imwrite(p_out(:,:,i), [path_out_tif filesep 'ref_' num2str(t) '_rot_' sprintf('%.3i',i) '.tif' ]);
    end
    
end


%% write all particles
n_particle = 0;
for t=1:n_ref
    for i=1:n_img
        n_particle = n_particle + size(particles{t,i}, 3);
    end    
end
p_out = zeros(2*w+1, 2*w+1, n_particle, 'uint16');
m=1;
for t=1:n_ref
    for i=1:n_img
        for j=1:size(particles{t,i},3)
            p_out(:,:,m) = particles{t,i}(:,:,j);
            m = m+1;
        end
    end
end

% write as spider-file
writeSPIDERfile([path_out filesep 'all.spi'], p_out, 'stack')

% write as single tif-files
path_out_tif = [path_out filesep 'all_tif'];
mkdir(path_out_tif)
for i=1:n_particle
    imwrite(p_out(:,:,i), [path_out_tif filesep 'all_' sprintf('%.3i',i) '.tif' ]);
end


%% view images and found particles
cc = varycolor(n_ref);
close all
fig_dim =2*[20 10];
cur_fig = figure('Visible','off', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

myleg = cell(n_ref,1);
for t=1:n_ref
    myleg{t} = ['Ref. ' num2str(t)];
end

w = box_size/2;

for i=1:n_img
    subplot(1, 2, 1)
    imagesc(img3(:,:,i)), colorbar,  colormap gray, axis image, hold on % img2
    for r=1:n_ref
        h(r) = plot(peaks_ref{r,i}(:,1), peaks_ref{r,i}(:,2),  '+', 'color', cc(r,:), 'MarkerSize', 1);
        
        for j=1:size(peaks_ref{r,i},1)
            rectangle('Position',[peaks_ref{r,i}(j,1)-w, peaks_ref{r,i}(j,2)-w, 2*w+1, 2*w+1], 'EdgeColor', cc(r,:))
        end
    end
    title(['Image ' num2str(i) '/' num2str(n_img) ])
    legend(h, myleg)
    hold off
    
    subplot(1, 2, 2)
    imagesc(images(:,:,i)), colorbar, colormap gray, axis image, hold on %
    for r=1:n_ref
        h(r)=plot(peaks_ref{r,i}(:,1), peaks_ref{r,i}(:,2),  'x', 'color', cc(r,:), 'MarkerSize', 1);
        for j=1:size(peaks_ref{r,i},1)
            rectangle('Position',[peaks_ref{r,i}(j,1)-w, peaks_ref{r,i}(j,2)-w, 2*w+1, 2*w+1], 'EdgeColor', cc(r,:))
        end
    end
    title(['Image ' num2str(i) '/' num2str(n_img) ])
    legend(h, myleg)
    hold off
       
    %pause
    print(cur_fig, '-dtiff', '-r300', [path_out filesep 'image2_' sprintf('%.03i',i) '.tif'])

end



%%
disp('finished')
