function [dY] = GCModelODE(t, Y, P, Q0)
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
SI = P.SI(n);
Uen = P.Uen.value(n);
QDF = P.results.QDF(n);
P2  = P.results.P2(n);

if (t <= 1000)
    GFast = P.GFast{1};  % Fasting glucose [g?]
else    
    GFast = P.GFast{2};
end

% Compute derived values.
QT0  = Q0 + P.results.QDF(1);

% Solve derivatives.
dG  = -GC.pg*(G-GFast) ... % NOTE: Removed infusion term.
          - SI*(G*(Q+QDF) - GFast*QT0)/(1 + GC.alphaG*(Q+QDF)) ...
          + GI.d2/GC.VG(P)*P2;      
dI  = -GC.nK*I - GC.nL/(1 + GC.alphaI*I)*I - GC.nI/GC.VI(P)*(I-Q) ...
          + Uen*(1-GC.xL)/GC.VI(P);
dQ  = GC.nI/GC.VQ(P)*(I-Q) - GC.nC*Q;

% Pack up outgoing variables.
dY = [dG; dI; dQ];

end
