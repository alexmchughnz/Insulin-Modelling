function kEmpt = GetStomachEmptyingRate(t, qSto, P)
% Uses a patient's data to determine stomach emptying rate at time t.
% INPUTS:
%   t   - queried time
%   P   - patient struct
% OUTPUT:
%   kempt - stomach emptying rate at time t

load('parameters.mat', 'GI')

% Extract meal times.
M = P.meal;
mealStarts = M.startTimes;              % [min]
mealEnds = M.startTimes + M.durations;  % [min]

% Add contribution from all meals occuring during this time.
DTot = GI.DTot;
currentMeals = (mealStarts < t) & (t < mealEnds);     % Which meals contribute at time=t? [logical]
MAGIC_SCALE_FACTOR = 1000/180.156;                    % From eating.m.
if (any(currentMeals))
    DTot = dot(currentMeals, M.sugar * MAGIC_SCALE_FACTOR);  % Meal size [g]
end

% Derived parameters.
alpha = 5/(2*DTot*(1-GI.b));
beta  = 5/(2*DTot*GI.c);
gamma = tanh(alpha*(qSto - GI.b*DTot)) - tanh(beta*(qSto - GI.c*DTot)) + 2;

% Estimate kempt at time=t.
kEmpt = GI.kMin + (GI.kMax - GI.kMin)/2 * gamma;

end