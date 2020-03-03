function [dY] = GCModelODE(t, Y, P, Q0)
% ODE for GC model. Use with ode45.
% Requires qGut(t) and QDF(t) - must be run AFTER GI model.
% INPUTS:
%   t   - time at which to evaluate ODE
%   Y   - states [P1; P2; G; I; Q] at time == t
%   P   - patient struct
% OUTPUT:
%   dY  - derivatives of states at time == t

global GC SC NP

% Assign incoming variables.
P1 = Y(1);
P2 = Y(2);
G  = Y(3);
I  = Y(4);
Q  = Y(5);

% Find patient dependent values.
index = (1 + floor(t));

if (t <= 1000)
    GFast = P.GFast{1};  % Fasting glucose [mol?]
else    
    GFast = P.GFast{2};
end
D = GetGlucoseDelivery(t, P);
SI = P.SI(index);
Uen = P.Uen.value(index);
QDF = P.results.QDF(index);
QT0  = Q0 + P.results.QDF(1);
QLocal = P.results.QDFLocal(index);

% Solve derivatives.
dP1 = -NP.d1*P1 + D;
dP2 = NP.d1*P1 - NP.d2*P2;
dG  = -GC.pg*(G-GFast) ... % NOTE: Removed infusion term.
          - SI*(G*(Q+QDF) - GFast*QT0)/(1 + GC.alphaG*(Q+QDF)) ...
          + NP.d2/GC.VG(P)*P2;      
dI  = -GC.nK*I - GC.nL/(1 + GC.alphaI*I)*I - GC.nI/GC.VI(P)*(I-Q) ...
          + SC.k3*QLocal + Uen*(1-GC.xL)/GC.VI(P);
dQ  = GC.nI/NP.VQ(P)*(I-Q) - GC.nC*Q;

% Pack up outgoing variables.
dY = [dP1; dP2; dG; dI; dQ];

end