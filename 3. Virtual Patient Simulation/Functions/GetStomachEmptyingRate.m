function kEmpt = GetStomachEmptyingRate(qSto)
% Uses a patient's data to determine stomach emptying rate at time t.
% INPUTS:
%   t   - queried time
%   P   - patient struct
% OUTPUT:
%   kempt - stomach emptying rate at time t


global GI


% Derived parameters.
alpha = 5/(2*GI.D*(1-GI.b));
beta  = 5/(2*GI.D*GI.c);
gamma = tanh(alpha*(qSto - GI.b*GI.D)) - tanh(beta*(qSto - GI.c*GI.D)) + 2;

% Estimate kempt at time=t.
kEmpt = GI.kMin + (GI.kMax - GI.kMin)/2 * gamma;

end