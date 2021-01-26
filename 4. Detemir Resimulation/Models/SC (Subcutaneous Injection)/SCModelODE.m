function [dY] = SCModelODE(t, Y, P)
% ODE for SC model. Use with ode45.
% INPUTS:
%   t   - time at which to evaluate ODE
%   Y   - states [ISC; QLocal] at time == t
%   P   - patient struct
% OUTPUT:
%   dY  - derivatives of states at time == t

global SC

%% Input
ISC    = Y(1);
QLocal = Y(2);

%% Variables
% Patient dependent.
IBolus = P.data.IBolus;  % [mU/min]

%% Computation
dISC    = -SC.k2*ISC + IBolus(t);
dQLocal = -SC.k2*QLocal + SC.k2*ISC - SC.kdi*QLocal;

%% Output
dY = [dISC;
      dQLocal];

end
