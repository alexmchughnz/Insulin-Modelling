function [dY] = GCModelODE(t, Y, P)
% ODE for GC model. Use with ode45.
% Requires qGut(t) and QDF(t) - must be run AFTER GI model.
% INPUTS:
%   t   - time at which to evaluate ODE
%   Y   - states [G; I; Q] at time == t
%   P   - patient struct
% OUTPUT:
%   dY  - derivatives of states at time == t

global GI GC

% Assign incoming variables.
G         = Y(1);
I         = Y(2);
Q         = Y(3);

% Find patient/time dependent values.
index = (1 + floor(t));
Ra = GI.f * GI.kAbs * P.results.qGut(index);  % Rate of glucose appearance in plasma
if (t <= 1000)
    GFast = P.GFast{1};  % Fasting glucose
else    
    GFast = P.GFast{2};
end
SI = P.SI(index);
Uen = P.Uen.value(index);
QDF = P.results.QDF(index);

% Solve derivatives.
dG = -GC.pg*(G - GFast) - SI*G*(Q+QDF)/(1 + GC.alphaG*(Q+QDF)) ... % NOTE: Removed infusion term.
          + (Ra + GC.EGP - GC.CNS)/GC.VG(P);
dI = -GC.nK*I - GC.nL*I/(1 + GC.alphaI*I) - GC.nI*(I - Q) + Uen*(1 - GC.xL)/GC.VI(P);
dQ = GC.nI*(I - Q) - GC.nC * Q/(1+GC.alphaG*Q);

% Pack up outgoing variables.
dY = [dG; dI; dQ];

end