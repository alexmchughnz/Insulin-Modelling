function P = FitInsulinSensitivity(P)
% Fits data to find SI over time for a patient.
% INPUTS:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with SI

global GI GC
load('parameters.mat', 'GI', 'GC')

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

ta = 0;                      % Start of current interval [min]
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

% Initial conditions.
YGI0 = [0.001;  % P1(t=0)        GI Model
        0]; 	% P2(t=0)  
    
YID0 = [0;   % ISC(t=0)          ID Model
        0;   % QDFLocal(t=0)
        0;   % QDBLocal(t=0)
        0;   % IDF(t=0)
        0;   % IDB(t=0)
        0;   % QDF(t=0)
        0];  % QDB(t=0)
    
YGC0 = [0;   % Q(t=0)            GC Model
        0;   % GA(t=0)
        0];  % Gb(t=0)
    
Y0 = [YGI0;
      YID0;
      YGC0];

% For each interval, integrate the collection of ODEs,
% then solve the dG equation for SI.
ccGA = 11;  % Column index of GA in Y.
ccGb = 12;  % Column index of Gb in Y.
for ii = 1 : numIntervals
    % First minute: finer solving.
    [t1, Y1] = ode45(@SIModelODE, ...
                         ta : ta+1, ...
                         Y0, ...
                         optionsShort, ...
                         ppG, ppI, P, Y0);
    Y0 = Y1(end, :)';  % Update ICs, picking up where we left off.
    
    % Remaining minutes: coarser solving.
    [t2, Y2] = ode45(@SIModelODE, ...
                     (ta+1 : ta+intervalDuration-1)', ...
                     Y0, ...
                     optionsLong, ...
                     ppG, ppI, P, Y0);
                    
    % Assemble sections (in 1 min intervals).              
    t    = [t1(1); t2];
    Y    = [Y1(1, :); Y2];

    % Solve linear system to find SI.
    % deltaG = -int{dGA}*SI + int{dGb}, therefore SI = -GA\(Gb-deltaG).
    GA = Y(:, ccGA);         % Coefficient of SI in equation.
    Gb = Y(:, ccGb);         % Added terms in equation.
    
    deltaG  = ppG(tb) - ppG(ta);  % Change in G over interval.  
    
    A = -GA;     
    b = Gb-deltaG;          
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
function [dY] = SIModelODE(t, Y, ppG, ppI, P, Y0)

global GI GC

%% Input
% Split up states.
YGI = Y(1:2);
YID = Y(3:9);
YGC = Y(10:12);

P2  = YGI(2);
QDF = YID(6);
Q   = YGC(1);

G   = ppG(t);   % Current blood glucose (interpolated)  [mmol/L]
I   = ppI(t);   % Current plasma insulin (interpolated) [mU/L]

% Split up initial conditions.
YID0 = Y0(3:9);
YGC0 = Y0(10:12);

QDF0 = YID0(6);
Q0   = YGC0(1);

%% Variables
% Time dependent.
n = (1 + floor(t));         % Index of current timestep.
GFast = P.GFast(t);         % Fasting glucose [mol?/L]
GInfusion = P.GInfusion(n); % Glucose infusion rate [mol/min]

% Patient dependent.
VG  = GC.VG(P);
VQ  = GC.VQ(P);

% Derived values.
QTFast  = Q0 + QDF0;
QT = Q + QDF;

%% Computation
% GC Model (reconstructed)
% dG = -dGA*SI + dGb. SI is found with dGA\dGb (linear system).
dQ     = GC.nI/VQ*(I-Q) - GC.nC*Q;
dGA    = (G*QT - GFast*QTFast)/(1 + GC.alphaG*QT);
dGb    = -GC.pg*(G-GFast) + GI.d2/VG*P2 + GInfusion/VG;
dYGC   = [dQ; dGA; dGb];

% GI Model
dYGI = GIModelODE(t, YGI, P);

% ID Model
dYID = IDModelODE(t, YID, P);

%% Output
dY = [dYGI;
      dYID;
      dYGC];
  
end
