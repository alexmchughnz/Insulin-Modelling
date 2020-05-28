function [dY] = GCModelODE(t, Y, tArray, qGut, QLocal, V)
% ODE for GC model. Use with ode45.
% Requires qGut(t) and QLocal(t) - must be run AFTER GI and ID models.
% INPUTS:
%   t      - time at which to evaluate ODE
%   Y      - states [G; I; Q] at time == t
%   qGut   - patient struct
%   QLocal - initial conditions of states
%   V      - specific parameter variant
% OUTPUT:
%   dY  - derivatives of states at time == t

global GI IN GC

%% Input
G  = Y(1);
I  = Y(2);
Q  = Y(3);


%% Variables
% Find time index.
t0 = tArray(1);
dt = tArray(2) - t0;
ii = floor((t-t0)/dt) + 1;

% Time dependent.
uEx = IN.k3 * QLocal(ii);
uEn = EstimateInsulinSecretion(G);
P = GI.f * GI.kAbs * qGut(ii);

% Variant dependent.
SI = V.SI;

%% Computation
dG  = -GC.pg*(G-GC.GFast) ...
          - SI*G*Q/(1 + GC.alphaG*Q) ...
          + (P + GC.EGP - GC.CNS)/GC.VG;
dI  = -GC.nK*I ...
          - GC.nL*I/(1 + GC.alphaI*I) ...
          - GC.nI*(I-Q) ...
          + (uEx + (1 - GC.xL)*uEn)/GC.VI;
dQ  = GC.nI*(I-Q) ...
          - GC.nC*Q/(1 + GC.alphaG*Q);

%% Output
dY = [dG;
      dI;
      dQ];

end
