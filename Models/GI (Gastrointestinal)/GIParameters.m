function GI = GIParameters(~)

% Initial Values
GI.P10     = 0.001; % Stomach glucose [mmol]
GI.P20     = 0;  % Gut glucose [mmol]

% Parameters
P1HalfLife = 20;    % Half-life of glucose in stomach (P1) [min]
P2HalfLife = 60;    % Half-life of glucose in gut (P2) [min]
GI.d1 = log(2)/P1HalfLife;

end