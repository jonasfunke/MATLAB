function [ dy ] = linking_on_fs(t, y, k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% y = [Linker, FS-wihtout-Linker, FS with one linker, FS with two linker, FS linked]
dy = zeros(5, 1);
dy(1) = - k(1) * y(1) * (y(2)+y(3)); %  Linker
dy(2) = - k(1) * y(1) * y(2); % FS without linker
dy(3) = + k(1) * y(1) * y(2) - k(1) * y(1) * y(3) - k(2) * y(3); % FS with one linker
dy(4) = + k(1) * y(1) * y(3); % FS with two linker
dy(5) = + k(2) * y(3); % FS linked
end

