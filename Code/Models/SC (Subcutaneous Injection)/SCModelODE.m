function [dY] = SCModelODE(t, Y, P)
% ODE for SC model. Use with ode45.
% INPUTS:
%   t   - time at which to evaluate ODE
%   Y   - states [ISC; QLocal] at time == t
%   P   - patient struct
% OUTPUT:
%   dY  - derivatives of states at time == t

SC = P.parameters.SC;

%% Input
ISC    = Y(1);
QLocal = Y(2);

%% Variables
% Patient dependent.
IBolus = P.data.IBolus(t);  % [mU/min]

%% Computation
dISC    = -SC.ks2*ISC + IBolus;
dQLocal = -SC.ks3*QLocal + SC.ks2*ISC - SC.kdi*QLocal;

%% Output
dY = [dISC;
      dQLocal];

end
