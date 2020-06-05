function [P] = GCModel(P, options)
% Function for GC model forward simulation.
% Pequires qGut(t) and QLocal(t) - must be run AFTER GI and ID models.
% INPUTS:
%   P        - patient struct, must have tArray, qGut(t) and QLocal(t)
%   options  - ode45 options
% OUTPUT:
%   P  - patient struct updated with model results 

global GC

% Set up initial conditions.
Y0 = [GC.G0;   
      GC.I0; 
      GC.Q0];
    
% Forward simulate.
[~, Y] = ode45(@GCModelODE, P.results.tArray, Y0, options, P, Y0);  

% Store results.
P.results.G = Y(:,1);
P.results.I = Y(:,2);
P.results.Q = Y(:,3);

end