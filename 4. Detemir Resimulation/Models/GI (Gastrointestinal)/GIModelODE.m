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
D = GetGlucoseDelivery(t, P);  % [mmol/min]
d2 = P.d2;

%% Computation
dP1 = -GI.d1*P1 + D;
dP2 = -d2*P2 + GI.d1*P1 ;

%% Output
dY = [dP1;
      dP2];

end
