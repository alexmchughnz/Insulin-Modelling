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

global GC SC

%% Input
G  = Y(1);  % [mmol/L]
I  = Y(2);  % [mU]
Q  = Y(3);  % [mU]
Q0 = Y0(3); % [mU]

%% Variables
% Time dependent.
n = GetTimeIndex(t, P.results.tArray);  % Index of current timestep.
SI        = P.results.SI(n);        % [L/mU/min]
Uen       = P.results.Uen(n);       % [mU/min]
P2        = P.results.P2(n);        % [mmol]
nL        = P.results.nL(n);        % [1/min]
xL        = P.results.xL(n);        % [1]
GFast     = P.data.GFast(t);        % Fasting glucose [mmol/L]

% Patient dependent.
d2 = P.results.d2;
VG = GC.VG(P);
VQ = GC.VQ(P);
nI = GC.nI(P);
nC = GC.nC(P);
nK = GC.nK(P);

% Trial dependent.
if P.data.GDelivery == "intravenous"
    GInput = P.data.GBolus(t) + P.data.GInfusion(n);  % [mmol/min]
else
    GInput = P.data.GInfusion(n);  % [mmol/min]
end

if P.data.IType == "detemir"
    % Include interstitial Detemir.
    Q0 = Q0 + P.results.QDF(1);  % [mU/L]
    Q = Q + P.results.QDF(n);   % [mU/L]
    
    IInput = 0;
    
else    
    if P.data.IDelivery == "intravenous"
        IInput = P.data.IBolus(t);  % [mU/min]        
        
    elseif P.data.IDelivery == "subcutaneous"
        IInput = SC.k3*P.results.QLocal(n);  % [mU/min]      
        
    else
        IInput = 0;
    end    
end

%% Computation
dG  = - GC.pg*(G-GFast) ...                       % Non-insulin-mediated uptake.
    - SI*(G*Q - GFast*Q0)/(1 + GC.alphaG*Q) ...   % Insulin-mediated uptake.
    + d2*P2/VG ...                                % Endogenous input from gut.
    + GInput/VG;                                  % Exogenous IV glucose input.

dI  = - nK*I ...                  % Renal clearance.
    - nL*I/(1 + GC.alphaI*I) ...  % Liver clearance.
    - nI/GC.VI*(I-Q) ...          % Transfer to Q compartment.
    + Uen*(1 - xL)/GC.VI ...      % Endogenous input.
    + IInput/GC.VI;               % Exogenous input (IV or subcut).

dQ  = nI/VQ*(I-Q) ...  % Transfer from I compartment.
    - nC*Q;            % Peripheral degradation.

%% Output
dY = [dG;
    dI;
    dQ];

end
