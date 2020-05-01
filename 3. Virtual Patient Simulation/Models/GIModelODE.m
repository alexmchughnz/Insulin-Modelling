function [dY] = GIModelODE(t, Y)
% ODE for GI model. Use with ode45.
% INPUTS:
%   t   - time at which to evaluate ODE
%   Y   - states [qSto1; qSto2; qGut] at time == t
%   P   - patient struct
% OUTPUT:
%   dY  - derivatives of states at time == t

global GI

%% Input
qSto1     = Y(1);
qSto2     = Y(2);
qGut      = Y(3);

%% Variables
% Patient dependent.
qSto = qSto1 + qSto2;
kEmpt = GetStomachEmptyingRate(qSto);

% Time dependent.
D = (t == 0)*GI.D;

%% Computation
dqSto1 = -GI.k21*qSto1 + D;
dqSto2 = -kEmpt*qSto2 + GI.k21*qSto1;
dqGut  = -GI.kAbs*qGut + kEmpt*qSto2;

%% Output
dY = [dqSto1;
      dqSto2;
      dqGut];

end
