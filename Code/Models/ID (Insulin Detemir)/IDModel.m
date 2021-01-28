function [P] = IDModel(P, options)
% Function for ID model forward simulation.
% INPUTS:
%   P        - patient struct, must have tArray
%   options  - ode45 options
% OUTPUT:
%   P  - patient struct updated with model results

global ID

if ~exist('options', 'var')
    options = odeset;
end

PrintStatusUpdate(mfilename, P, "Begin solving...")


% Set up initial conditions.
Y0 = [ID.ISC0;
      ID.QDFLocal0;
      ID.QDBLocal0;
      ID.IDF0;
      ID.IDB0;
      ID.QDF0;
      ID.QDB0];
  
% Forward simulate.
[~, Y] = ode45(@IDModelODE, P.results.tArray, Y0, options, P);  

% Store results.
P.results.IDH      = Y(:,1);
P.results.QDFLocal = Y(:,2);
P.results.QDBLocal = Y(:,3);
P.results.IDF      = Y(:,4);
P.results.IDB      = Y(:,5);
P.results.QDF      = Y(:,6);
P.results.QDB      = Y(:,7);

end
