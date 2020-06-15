function [dY] = GCModelODE(t, Y, P, Y0)
% ODE for GC model. Use with ode45.
% Requires P1(t) and P2(t) - must be run AFTER GI model.
% INPUTS:
%   t   - time at which to evaluate ODE
%   Y   - states [G; I; Q] at time == t
%   P   - patient struct
%   Y0  - initial conditions of states
% OUTPUT:
%   dY  - derivatives of states at time == t

global GI GC SC

%% Input
G  = Y(1);
I  = Y(2);
Q  = Y(3);
Q0 = Y0(3);

%% Variables
% Time dependent.
n = (1 + floor(t));  % Index of current timestep.
SI        = P.SI(n);
Uen       = P.Uen.value(n);
P2        = P.results.P2(n);
GInfusion = P.GInfusion(n);   % Glucose infusion rate [mol/min]
GFast     = P.GFast(t);       % Fasting glucose [mol?/L]

% Patient dependent.
VG  = GC.VG(P);
VI  = GC.VI(P);
VQ  = GC.VQ(P, SC);
xL = P.xL;
nL = P.nL;
 
% Derived values.
QTFast  = Q0 + P.results.QDF(1);
QT      = Q + P.results.QDF(n);

%% Computation
dG  = -GC.pg*(G-GFast) ...
          - SI*(G*QT - GFast*QTFast)/(1 + GC.alphaG*QT) ...
          + GI.d2/VG*P2 + GInfusion/VG;      
dI  = -GC.nK*I - nL/(1 + GC.alphaI*I)*I - GC.nI/VQ*(I-Q) ...
          + Uen*(1 - xL)/VI;
dQ  = GC.nI/VQ*(I-Q) - GC.nC*Q;

%% Output
dY = [dG;
      dI;
      dQ];

end
