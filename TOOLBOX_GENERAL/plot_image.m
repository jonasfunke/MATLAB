function [ ] = plot_image( img, frac )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    dc = max(img(:))-min(img(:));
    clim = [median(img(:))-frac*dc median(img(:))+frac*dc];
    imagesc(img, clim), axis image, colormap gray
end

