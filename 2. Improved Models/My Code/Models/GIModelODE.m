function [dY] = GIModelODE(t, Y, P)
% ODE for GI model. Use with ode45.
% INPUTS:
%   t   - time at which to evaluate ODE
%   Y   - states [qSto1; qSto2; qGut] at time == t
%   P   - patient struct
% OUTPUT:
%   dY  - derivatives of states at time == t

global GI

% Assign incoming variables.
P1     = Y(1);
P2     = Y(2);

% Find derived values.
D = GetGlucoseDelivery(t, P);

% Solve derivatives.
dP1 = -GI.d1*P1 + D;
dP2 = GI.d1*P1 - GI.d2*P2;

% Pack up outgoing variables.
dY = [dP1; dP2];

end
