function [ a ] = test_function()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    disp('test function:')
    N = 100000;
    M = 2;
    a = zeros(N,1); 
    parfor i = 1:N 
        a(i) = max(eig(rand(M)));
    end
end

