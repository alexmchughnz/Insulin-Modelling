function [glucose_infusion_rate] = infusion(sys, t)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% t_start = 0; %minutes(sys.Data.bg3_time(2) - sys.Data.bg3_time(1));
t_end = 12*60 - 135; %12 hr infusion minus the time it was going

% glucose_infusion_rate = 0; %by default, if not, then it is overridden below

if t <= t_end && sys.Data.PtNo == 1
    glucose_infusion_rate = 1.54/180.156/60; %infusion rate in mol/min for a 5% infusion
else
    glucose_infusion_rate = 0;
end


end