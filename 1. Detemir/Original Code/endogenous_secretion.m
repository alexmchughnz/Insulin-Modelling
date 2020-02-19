%insulin secretion plots
%17/04/18

clear all, close all, clc

k_sec = 14.9;
k_offset = -50;
BGL = [0 : 0.1 : 50];
u_en = zeros(1: length(BGL));
G_fast = 4.8;
u_min = 16.7;
u_max = 267;

for ii = 1 : length(BGL)
    
    f_G = k_sec*BGL(ii) + k_offset;
    
    if BGL(ii) <= G_fast
        u_en(ii) = u_min;
    elseif f_G >= u_max
        u_en(ii) = u_max;
    else
        u_en(ii) = f_G;
    end
    
    
end

plot(BGL, u_en)