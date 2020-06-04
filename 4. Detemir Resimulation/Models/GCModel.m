function [R] = GCModel(R, O)
% Function for GC model forward simulation.
% Requires qGut(t) and QLocal(t) - must be run AFTER GI and ID models.
% INPUTS:
%   R  - results struct, must have tArray, qGut(t) and QLocal(t)
%   O  - ode45 options
% OUTPUT:
%   R  - results struct updated with model results 

global GC

% Set up initial conditions.
Y0 = [GC.G0;   
      GC.I0; 
      GC.Q0];
    
% Forward simulate.
[~, Y] = ode45(@GCModelODE, R.tArray, Y0, O, P, Y0);  

% Store results.
R.G = Y(:,1);
R.I = Y(:,2);
R.Q = Y(:,3);

end