function [P] = SCModel(P, options)
% Function for SC model forward simulation.
% INPUTS:
%   P        - patient struct, must have tArray
%   options  - ode45 options
% OUTPUT:
%   P  - patient struct updated with model results 

global SC

if ~exist('options', 'var')
    options = odeset;
end

% Set up initial conditions.
Y0 = [SC.ISC0;
      SC.QLocal0];
  
% Forward simulate.
[~, Y] = ode45(@SCModelODE, P.results.tArray, Y0, options, P);  

% Store results.
P.results.ISC = Y(:,1);
P.results.QLocal = Y(:,2);

end