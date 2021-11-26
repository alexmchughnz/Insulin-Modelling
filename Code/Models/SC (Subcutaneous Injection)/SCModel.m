function [P] = SCModel(P)
% Function for SC model forward simulation.
% INPUTS:
%   P        - patient struct, must have tArray
%   options  - ode45 options
% OUTPUT:
%   P  - patient struct updated with model results 

SC = P.parameters.SC;

PrintStatusUpdate(P, "Begin solving...")


% Set up initial conditions.
Y0 = [SC.ISC0;
      SC.QLocal0];
  
% Forward simulate.
[~, Y] = ode45(@SCModelODE, P.results.tArray, Y0, [], P);  

% Store results.
P.results.ISC = Y(:,1);
P.results.QLocal = Y(:,2);

end
