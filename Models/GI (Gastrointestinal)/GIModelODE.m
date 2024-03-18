function [dY] = GIModelODE(t, Y, P)
% ODE for GI model. Use with ode45.
% INPUTS:
%   t   - time at which to evaluate ODE
%   Y   - states [P1; P2] at time == t
%   P   - patient struct
% OUTPUT:
%   dY  - derivatives of states at time == t

GI = P.parameters.GI;

%% Input
P1     = Y(1);
P2     = Y(2);

%% Variables
n = GetTimeIndex(t, P.results.tArray);

% Patient dependent.
if P.data.GDelivery == "enteral"
    D = P.results.GBolus(n);  % [mmol/min]
else
    D = 0;
end
d1 = GI.d1;
d2 = P.results.d2;

%% Computation
dP1 = D - d1*P1;
dP2 = d1*P1 - d2*P2;

%% Output
dY = [dP1;
      dP2];

end
