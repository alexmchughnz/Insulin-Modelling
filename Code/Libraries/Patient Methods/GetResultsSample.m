function [time, value] = GetResultsSample(P, time, results)
% Retrieves value of simulated results at discrete time points.
% INPUTS:
%   results - field of P.results containing results array
% OUTPUT:
%   time  - time points [min]
%   value - value of results at each time
    
    iiData = SearchArray(time, P.results.tArray);    
    value = results(iiData);
    
end