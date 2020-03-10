function [dY] = GCModelODE(t, Y, P, Y0)
% ODE for GC model. Use with ode45.
% Requires P1(t) and P2(t) - must be run AFTER GI model.
% INPUTS:
%   t   - time at which to evaluate ODE
%   Y   - states [P1; P2; G; I; Q] at time == t
%   P   - patient struct
% OUTPUT:
%   dY  - derivatives of states at time == t

global GI GC

% Assign incoming variables.
G  = Y(1);
I  = Y(2);
Q  = Y(3);

% Retrieve patient dependent values.
n = (1 + floor(t));  % Index of current timestep.
SI  = P.SI(n);
Uen = P.Uen.value(n);
P2  = P.results.P2(n);
GFast = P.GFast(t);         % Fasting glucose [mol?/L]
GInfusion = P.GInfusion(n); % Glucose infusion rate [mol/min]

VG  = GC.VG(P);
VI  = GC.VI(P);
VQ  = GC.VQ(P);
 
% Compute derived values.
Q0      = Y0(3);
QTFast  = Q0 + P.results.QDF(1);
QT      = Q + P.results.QDF(n);

% Solve derivatives.
dG  = -GC.pg*(G-GFast) ... % NOTE: Removed infusion term.
          - SI*(G*QT - GFast*QTFast)/(1 + GC.alphaG*QT) ...
          + GI.d2/VG*P2 + GInfusion/VG;      
dI  = -GC.nK*I - GC.nL/(1 + GC.alphaI*I)*I - GC.nI/VI*(I-Q) ...
          + Uen*(1 - GC.xL)/VI;
dQ  = GC.nI/VQ*(I-Q) - GC.nC*Q;

% Pack up outgoing variables.
dY = [dG; dI; dQ];

end
