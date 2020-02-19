function [Uen] = uen(sys, t)
%Endogenous insulin seceretion by the pancreas
%   Function of plasma glcose concentration
% ksec = 4.5;
% koffset = -13;
% 
% Umin = 16.7;
% Umax = 267;
% 
% fg = ksec*G + koffset;
% 
% Uen = min(max(Umin,fg),Umax); 


for ii = 1 : length(sys.GC.SI)
    if t >= ii-1 && t <= ii
        Uen = sys.GC.Uen(ii);
    end
end
    
end

