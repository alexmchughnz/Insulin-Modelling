function [R] = GCModel(R, O, V)
% Function for GC model forward simulation.
% Requires qGut(t) and QLocal(t) - must be run AFTER GI and ID models.
% INPUTS:
%   R  - results struct, must have tArray, qGut(t) and QLocal(t)
%   O  - ode45 options
%   V  - parameter variants
% OUTPUT:
%   R  - results struct updated with model results 

global GC

% Set up initial conditions.
Y0 = [GC.G0;   
      GC.I0;  
      GC.Q0]; 
    
% Forward simulate.
[~, Y] = ode45(@GCModelODE, R.tArray, Y0, O, R, V);  

% Store results.
R.G = Y(:,1);
R.I = Y(:,2);
R.Q = Y(:,3);

end


function [dY] = GCModelODE(t, Y, R, V)
% ODE for GC model. Use with ode45.
% GCPUTS:
%   t      - time at which to evaluate ODE
%   Y      - states [G; I; Q] at time == t
%   R      - results struct, must have qGut(t) and QLocal(t)
%   V      - specific parameter variant
% OUTPUT:
%   dY  - derivatives of states at time == t

global GI GC 

%% Input
G  = Y(1);
I  = Y(2);
Q  = Y(3);

%% Variables
% Find time index.
t0 = R.tArray(1);
dt = R.tArray(2) - t0;
ii = floor((t-t0)/dt) + 1;

% Time dependent.
uEn = EstimateInsulinSecretion(G);  % Endogenous insulin inflow.

QLocal = R.QLocal(ii);
uEx = GC.uEx(QLocal);               % Exogenous insulin inflow.

qGut = R.qGut(ii);
P = GI.Ra(qGut);                % Exogenous glucose inflow.

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
