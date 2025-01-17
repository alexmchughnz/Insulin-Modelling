function P = FitInsulinSensitivity(P)
% Fits data to find SI over time for a patient.
% INPUTS:
%   P      - patient struct
% OUTPUT:
%   P          - modified patient struct with SI

PrintStatusUpdate(P, "Fitting SI...")

%% Setup
% Interpolate G and I with piecewise polynomials.
[tG, vG] = GetData(P.data.G);
ppG = griddedInterpolant(tG, vG);  % G(t) [mmol/L] piecewise polynomial. Use as function.

% Create SI array and initial time boundaries.
defaultSI = 10.8e-4;
P.results.SI = defaultSI;  % SI value over trial.

tStart = P.results.tArray(1);  % Start of current interval [min]
tEnd = P.results.tArray(end);  % End of current interval [min]

%% Compute
% ODE solver settings.
optionsShort = odeset('RelTol',1e-5, ...  % Options for first minute.
    'AbsTol',1e-4, ...
    'MaxStep',0.2, ...
    'InitialStep',0.01);
optionsLong = odeset('RelTol',1e-5, ...   % Options for other minutes.
    'AbsTol',1e-4, ...
    'InitialStep',0.1);

% Set up for different model types.
P = SolveSystem(P);

[~, vG] = GetData(P.data.G);
[~, vI] = GetData(P.data.I);  % [mU/L]

G0 = vG(1);
I0 = vI(1);
Q0 = I0/2;  % Subcut Q assumed to be half of plasma I at t=0.
GA0 = 0;
Gb0 = 0;

Y0 = [G0;
    I0;
    Q0;
    GA0;
    Gb0];

ccGA = 4;  % Column index of GA in Y.
ccGb = 5;  % Column index of Gb in Y.

% For each interval, integrate the collection of ODEs,
% then solve the dG equation for SI.
% First minute: finer solving.
[t1, Y1] = ode45(@GCModelODE, ...
    tStart : tStart+1, ...
    Y0, ...
    optionsShort, ...
    P, Y0);

% Remaining minutes: coarser solving.
Y0 = Y1(end, :)';  % Update ICs, picking up where we left off.
[t2, Y2] = ode45(@GCModelODE, ...
    (tStart+1 : tEnd)', ...
    Y0, ...
    optionsLong, ...
    P, Y0);

% Assemble sections (in 1 min intervals).
t    = [t1(1); t2];
Y    = [Y1(1, :); Y2];

% Solve linear system to find SI.
% At each row, between t1 and t2:
%       delta{G} = int{dGA}*SI + int{dGb}
%        delta{GA}*SI = delta{G} - delta{Gb}
% therefore:
%           SI = delta{GA}\(delta{G} - delta{Gb})
GA = Y(:, ccGA);         % Coefficient of SI in equation.
Gb = Y(:, ccGb);         % Added terms in equation.

G  = ppG(P.results.tArray);  % Change in G over interval.

first = 1:numel(P.results.tArray)-1;
second = 2:numel(P.results.tArray);

A = GA(second) - GA(first);
b = (G(second) - G(first)) - (Gb(second) - Gb(first));
SI = A\b;

% Write estimated data into patient struct, overwriting defaults.
P.results.SI = SI;  % [L/mU/min]

end