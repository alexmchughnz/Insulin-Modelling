function [dY] = GCModelODE(t, Y)
% ODE for GC model. Use with ode45.
% Requires P1(t) and P2(t) - must be run AFTER GI model.
% INPUTS:
%   t   - time at which to evaluate ODE
%   Y   - states [P1; P2; G; I; Q] at time == t
%   P   - patient struct
%   Y0  - initial conditions of states
% OUTPUT:
%   dY  - derivatives of states at time == t

global GC

%% Input
G  = Y(1);
I  = Y(2);
Q  = Y(3);

%% Computation
dG  = -GC.pg*(G-GC.GFast) ...
          - SI*G*Q/(1 + GC.alphaG*Q) ...
          + (GC.EGP - GC.CNS)/GC.VG;      %TODO: missing P(t)?
dI  = -GC.nK*I ...
          - P.nL*I/(1 + GC.alphaI*I) ...
          - GC.nI*(I-Q) ...
          + (uEx + uEn*(1 - P.xL))/VI;    %TODO: where does uEx come from?
dQ  = GC.nI*(I-Q) ...
          - GC.nC*Q/(1 + GC.alphaG*Q);

%% Output
dY = [dG;
      dI;
      dQ];

end
