function [P] = GIModel(P)
% Function for GC model forward simulation.
% INPUTS:
%   P        - patient struct, must have tArray
%   options  - ode45 options
% OUTPUT:
%   P  - patient struct updated with model results 

GI = P.parameters.GI;

PrintStatusUpdate(P, "Begin solving...")

% Set up initial conditions.
Y0 = [GI.P10;
      GI.P20];
  
% Forward simulate.
[~, Y] = ode45(@GIModelODE, P.results.tArray, Y0, P.options.ODE, P);  

% Store results.
P.results.P1 = Y(:,1);
P.results.P2 = Y(:,2);

end
