function [ img_bg, bg ] = bg_correct_ui( img, img_title, N_ref )
% select and subtract background from image img
%   Detailed explanation goes here

%x = reshape(img, size(img,1)*size(img,2), 1);


%colormap('Gray');
%scaling = [0 mean(x)+2*std(x)]; %[mean(x)-2*std(x) mean(x)+2*std(x)];

% imagesc(img, scaling), colorbar
%title([ img_title ' select background' ], 'FontSize' , 18)


%%
if ~exist('N_ref', 'var') % use average frame
   N_ref = 1;
else
    if N_ref < 3
        disp('Number of reference points to small. Changed to 3.')
        N_ref = 3;
    end
end

[I, areas] = integrate_areas({img}, N_ref, 1, [1 1]); %cell of images, number of bands, 1=all bands habe the same size
close all

%% fit a first order polynomial to the bg-points
bg_points = zeros(N_ref,3);
for i=1:N_ref
    bg_points(i,1) = areas(i,2)+areas(i,4)/2; % x-coordinate
    bg_points(i,2) = areas(i,1)+areas(i,3)/2; % y-coordinate
    bg_points(i,3) = mean( mean(  img(areas(i,2):areas(i,2)+areas(i,4)  ,   areas(i,1):areas(i,1)+areas(i,3)   ) ));
end

if N_ref > 1
    bg = fit( bg_points(:,1:2), bg_points(:,3), 'poly11');

    % subtract
    [i, j] = meshgrid(1:1:size(img,1), 1:1:size(img,2));
    ti = i'; tj = j';
    z = bg([ti(:) tj(:)]);
    Z = reshape(z, size(img));
    img_bg =  img - Z ;

    plot_image_ui(img_bg);
else
    img_bg =  img - bg_points(1,3) ;
    bg = bg_points(1,3);
end



%{
plot_image_ui(img)
title( ['Select background of ' img_title])

h = imrect;
position = wait(h);
pos = int32(getPosition(h)); % [xmin ymin width height]

bg = mean( mean(  img(pos(2):pos(2)+pos(4)  ,   pos(1):pos(1)+pos(3)   ) ));

img_bg =  img - bg ;
close all
%}

end

 