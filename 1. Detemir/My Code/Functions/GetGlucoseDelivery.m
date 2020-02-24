function D = GetGlucoseDelivery(t, P)
% Uses a patient's meal data to determine glucose delivery rate at time t.
% INPUTS:
%   P   - patient struct
%   t   - queried time
% OUTPUT:
%   D   - glucose delivery rate at time t

% Extract meal times.
M = P.meal;
mealStarts = M.startTimes;              % [min]
mealEnds = M.startTimes + M.durations;  % [min]

% Add contribution from all meals occuring during this time.
currentMeals = (mealStarts < t) & (t < mealEnds);       % Which meals contribute at time=t? [logical]
mealRates = M.carbs ./ M.durations;                     % Glucose rates of all meals [g/min?]
MAGIC_SCALE_FACTOR = 1000/180.156;                      % From eating.m.
D = MAGIC_SCALE_FACTOR * dot(currentMeals, mealRates);  % Glucose delivery rate [g/min?]

end

