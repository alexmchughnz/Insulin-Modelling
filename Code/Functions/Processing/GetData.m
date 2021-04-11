function [time, value] = GetData(S)
% Retrieves data and time points.
% INPUTS:
%   S     - field of P containing data struct
% OUTPUT:
%   time  - time points [min]
%   value - value of data at each time

    time = S.time; 
    value = S.value;    
end