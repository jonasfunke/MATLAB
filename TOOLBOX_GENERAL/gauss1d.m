function [ y ] = gauss1d( c, x )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    y = c(3) .* exp( - (x-c(1)).^2 ./ (2.*c(2).^2)  );

end

