function [dY] = IDModelODE(t, Y)
% ODE for ID model. Use with ode45.
% INPUTS:
%   t   - time at which to evaluate ODE
%   Y   - states [IDH; QDFLocal; QDBLocal; IDF; IDB; QDF; QDB] at time == t
%   P   - patient struct
% OUTPUT:
%   dY  - derivatives of states at time == t

global IN

%% Input
ISC    = Y(1);
QLocal = Y(2);

%% Variables
% Time dependent.
IBolus = (t == IN.TBolus)*IN.IBolus;

%% Computation
dISC    = -IN.k2*ISC + IBolus;
dQLocal = -IN.k3*QLocal ...
              + IN.k2*ISC ...
              - IN.kdi*QLocal;

%% Output
dY = [dISC;
      dQLocal];

end

