function D = GetGlucoseDelivery(t, P)
% Uses a patient's meal data to determine glucose delivery rate at time t.
% INPUTS:
%   P   - patient struct
%   t   - queried time
% OUTPUT:
%   D   - glucose delivery rate at time t

% Extract meal times in minutes.
M = P.meal;
mealStarts = M.startTimes;              % [min]
mealEnds = M.startTimes + M.durations;  % [min]

% Add contribution from all meals occuring during this time.
MAGIC_SCALE_FACTOR = 1000/180.156; % From eating.m.
currentMeals = (mealStarts < t) & (t < mealEnds); % [logical]
D = dot(currentMeals, M.carbs./M.durations*MAGIC_SCALE_FACTOR); % Glucose delivery rate [g/min?]

end

