function D = GetGlucoseDelivery(t, P)
% Uses a patient's meal data to determine glucose delivery rate at time t.
% INPUTS:
%   P   - patient struct
%   t   - queried time
% OUTPUT:
%   D   - glucose delivery rate at time t

global C

% Extract meal times.
M = P.meal;
mealStarts = M.startTimes;              % [min]
mealEnds = mealStarts + M.durations;  % [min]

% Add contribution from all meals occuring during this time.
currentMeals = (mealStarts < t) & (t < mealEnds); % Which meals contribute at time=t? [logical]
mealRates = M.carbs ./ M.durations;               % Glucose rates of all meals [g/min?]

D = dot(currentMeals, mealRates) / C.MGlucose * 1000;  % Glucose delivery rate [mmol/min]

end

