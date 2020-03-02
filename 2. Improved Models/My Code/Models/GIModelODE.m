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
qSto1     = Y(1);
qSto2     = Y(2);
qGut      = Y(3);

% Find derived values.
kEmpt = GetStomachEmptyingRate(t, (qSto1+qSto2), P);
D = GetGlucoseDelivery(t, P);

% Solve derivatives.
dqSto1 = D - GI.k21*qSto1;
dqSto2 = GI.k21*qSto1 - kEmpt*qSto2;
dqGut  = kEmpt*qSto2 - GI.kAbs*qGut;

% Pack up outgoing variables.
dY = [dqSto1; dqSto2; dqGut];

end