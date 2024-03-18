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
n = GetTimeIndex(t, P.results.tArray);
IBolus = P.results.IBolus(n); % [mU/min]

%% Computation
dISC    = IBolus - SC.ks2*ISC;
dQLocal = SC.ks2*ISC - (SC.ks3 + SC.kdi)*QLocal;

%% Output
dY = [dISC;
      dQLocal];

end
