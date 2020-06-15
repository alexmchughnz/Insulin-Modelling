global GI

% Initial Values
GI.P10     = 0.001; % Stomach glucose [mmol]
GI.P20     = 0;  % Gut glucose [mmol]

% Parameters
P1HalfLife = 20;    % Half-life of glucose in stomach (P1) [min]
P2HalfLife = 60;    % Half-life of glucose in gut (P2) [min]
GI.d1 = log(2)/P1HalfLife;
GI.d2 = log(2)/P2HalfLife;  % NOTE: dummy value, needs parameter ID.

GI.k21 = 0.054;     % Rate constant of stomach grinding [1/min]
GI.kAbs = 0.071;    % rate glucose absorbed into bloodstream [1/min]
GI.b = 0.69;        % Stomach emptying parameters [1]
GI.c = 0.17;        % ''
GI.kMax = 0.054;    % Maximum value of kEmpt [1/min]
GI.kMin = 0.006;    % Minimum value of kEmpt [1/min]
GI.f = 0.8;         % Scaling factor for incomplete absorption [1]
