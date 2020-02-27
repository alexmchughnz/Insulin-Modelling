function P = FitInsulinSensitivity(GI, GC, P)
% Fits data to find SI over time for a patient.
% INPUTS:
%   GC  - glycaemic control parameter set
%   GI  - gastrointestinal parameter set
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with SI

%% Setup
% Define time range.
intervalDuration = 360;  % Time per interval [min]
numIntervals = floor(P.simDuration/intervalDuration);

% Interpolate G and I with piecewise polynomials.
iiStart = find(P.G{3}.time == P.simTime(1));   % Sim start [index]
iiEnd   = find(P.G{3}.time == P.simTime(end)); % Sim end [index]

tG = minutes(P.G{3}.time(iiStart:iiEnd) - P.G{3}.time(iiStart)); % Time over sim period [min]
vG = P.G{3}.value(iiStart:iiEnd);                                % Glucose over sim period [min]
tI = minutes(P.I.time - P.I.time(1));                            % Time over sim period [min]
vI = P.I.value;                                                  % Plasma insulin over sim period [min]

ppG = griddedInterpolant(tG, vG);  % G(t) piecewise polynomial. Use as function.
ppI = griddedInterpolant(tI, vI);  % I(t) piecewise polynomial. Use as function.

% Create SI array and initial time boundaries.
defaultSI = 10.8e-4;

P.SI = ones(P.simDuration, 1) * defaultSI;             % SI value at each minute in trial.
intervalSI = ones(numIntervals, 1) * defaultSI;        % SI value at each interval in sim.
minuteSI = zeros(numIntervals * intervalDuration, 1);  % SI value at each minute in sim.

ta = 1;                      % Start of current interval [min]
tb = ta + intervalDuration;  % End of current interval [min]

%% Compute
% ODE solver settings.
optionsShort = odeset('RelTol',1e-5, ...  % Options for first minute.
                      'AbsTol',1e-4, ...
                      'MaxStep',0.2, ...
                      'InitialStep',0.01);
optionsLong = odeset('RelTol',1e-5, ...   % Options for other minutes.
                     'AbsTol',1e-4, ...
                     'InitialStep',0.1);
Y0 = [0.001;  % qSto1(t=0)
      0;      % qSto2(t=0)      
      0;      % qGut(t=0)
      9;      % Q(t=0)
      0;      % GA(t=0)
      0];     % Gb(t=0)

% For each interval, integrate the collection of ODEs,
% then solve the dG equation for SI.
for ii = 1 : numIntervals
    % First minute: finer solving.
    [t1, Y1] = ode45(@GIModelODE, ...
                         ta : ta+1, ...
                         Y0, ...
                         optionsShort, ...
                         ppG, ppI, GI, GC, P);
    Y0 = Y1(end, :)';  % Update ICs, picking up where we left off.
    
    % Remaining minutes: coarser solving.
    [t2, Y2] = ode45(@GIModelODE, ...
                        (ta+1 : ta+intervalDuration-1)', ...
                        Y0, ...
                        optionsLong, ...
                        ppG, ppI, GI, GC, P);
                    
    % Assemble sections (in 1 min intervals).              
    t    = [t1(1); t2];
    Y    = [Y1(1, :); Y2];

    % Solve linear system to find SI.
    % dG = -dGA*SI + dGb, therefore SI = -dGA\(dGb-dG).
    dGA = Y(:, 5);         % Coefficient of SI in equation.
    dGb = Y(:, 6);         % Added terms in equation.
    dG  = ppG(tb) - ppG(ta);  % Change in G over interval.  
    
    A = -dGA;     
    b = dGb-dG;          
    intervalSI(ii) = A\b;
    
    % Set SI value over current interval to computed value.
    intervalRange = 1 + (ii-1)*intervalDuration : ii*intervalDuration;  % [min]
    minuteSI(intervalRange) = intervalSI(ii);
    
    % Update conditions for next iteration.
    ta = ta + intervalDuration;
    tb = tb + intervalDuration;
       
    Y0 = Y(end, :);
end

% Write estimated data into patient struct, overwriting defaults.
P.SI(1:length(minuteSI)) = minuteSI;

fprintf('P%d: SI fit successfully. SI(*1e+4) at %d min intervals = ', ...
        P.patientNum, intervalDuration)
disp(intervalSI'*1e4)

end


% Adapted from fit_SI/FAERIES_integrals.
function [dY] = GIModelODE(t, Y, ppG, ppI, GI, GC, P)

% HARDCODED VALUE from fit_SI, to "bypass the nasty shit"?
EGP = 0.96;

% Assign incoming variables.
qSto1  = Y(1);
qSto2  = Y(2);
qGut   = Y(3);
Q      = Y(4);

% Retrieve required derived values.
kEmpt = GetStomachEmptyingRate(t, (qSto1+qSto2), GI, P);
D = GetGlucoseDelivery(t, P); % Current glucose delivery rate [g/min]
G = ppG(t);                   % Current blood glucose (interpolated) [mmol/L]
I = ppI(t);                   % Current plasma insulin (interpolated) [mU/L]
Ra = GI.f * GI.kAbs * qGut;   % Rate of glucose appearance in plasma

% Gastrointestinal model DEs.
dqSto1 = D - GI.k21*qSto1;
dqSto2 = GI.k21*qSto1 - kEmpt*qSto2;
dqGut  = kEmpt*qSto2 - GI.kAbs*qGut;

% Glycaemic control model DEs.
% dG = -dGA*SI + dGb. SI is found with dGA\dGb (linear system).
dQ     = GC.nI*(I - Q) - GC.nC * Q/(1+GC.alphaG*Q);
dGA    = G * Q/(1 + GC.alphaG*Q);
% NOTE: Removed infusion term.
%dGb    = (Ra + EGP + infusion(sys, t) - sys.GC.CNS)/sys.GC.VG - sys.GC.pg*G;
dGb    = -GC.pg*G + (Ra + EGP - GC.CNS)/GC.VG(P);

% Pack up outgoing variables.
dY = [dqSto1;
      dqSto2;
      dqGut;
      dQ;
      dGA;
      dGb];

end
