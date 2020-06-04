function [R] = GIModel(R, O)
% Function for GI model forward simulation.
% INPUTS:
%   R  - results struct, must have tArray
%   O  - ode45 options
% OUTPUT:
%   R  - results struct updated with model results 

global GI

% Set up initial conditions.
Y0 = [GI.P10;
      GI.P20];
  
% Forward simulate.
[~, Y] = ode45(@GIModelODE, R.tArray, Y0, O, P);  

% Store results.
P.results.P1 = Y(:,1);
P.results.P2 = Y(:,2);

end
