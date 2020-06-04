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


function [dY] = GIModelODE(t, Y, P)
% ODE for GI model. Use with ode45.
% INPUTS:
%   t   - time at which to evaluate ODE
%   Y   - states [P1; P2] at time == t
%   P   - patient struct
% OUTPUT:
%   dY  - derivatives of states at time == t

global GI

%% Input
P1     = Y(1);
P2     = Y(2);

%% Variables
% Patient dependent.
D = GetGlucoseDelivery(t, P);

%% Computation
dP1 = -GI.d1*P1 + D;
dP2 = GI.d1*P1 - GI.d2*P2;

%% Output
dY = [dP1;
      dP2];

end
