function [R] = IDModel(R, O, V)
% Function for ID model forward simulation.
% INPUTS:
%   R  - results struct, must have tArray
%   O  - ode45 options
%   V  - parameter variants
% OUTPUT:
%   R  - results struct updated with model results 

global ID

% Set up initial conditions.
Y0 = [ID.ISC0;
      ID.QDFLocal0;
      ID.QDBLocal0;
      ID.IDF0;
      ID.IDB0;
      ID.QDF0;
      ID.QDB0];
  
% Forward simulate.
[~, Y] = ode45(@IDModelODE, R.tArray, Y0, O);  

% Store results.
R.IDH      = Y(:,1);
R.QDFLocal = Y(:,2);
R.QDBLocal = Y(:,3);
R.IDF      = Y(:,4);
R.IDB      = Y(:,5);
R.QDF      = Y(:,6);
R.QDB      = Y(:,7);

end
