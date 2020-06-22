function P = FitInsulinSensitivity(P)
% Fits data to find SI over time for a patient.
% INPUTS:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with SI

global C GI GC
load('parameters.mat', 'GI', 'GC')

%% Setup
% Define time range.
intervalDuration = 360;  % Time per interval [min]
numIntervals = floor(P.simDuration()/intervalDuration);

% Forward simulate ID model for IDF.
P = IDModel(P); 
IDF = P.results.IDF; % [mU/L]

% Interpolate G and I with piecewise polynomials.
[tG, vG] = GetSimTime(P, P.data.G{3});
ppG = griddedInterpolant(tG, vG);  % G(t) [mmol/L] piecewise polynomial. Use as function.

[tITotal, vITotal] = GetSimTime(P, P.data.I);% Plasma insulin over sim period [pmol/L]
I = C.pmol2mU(vITotal) - IDF(tITotal + 1);
ppI = griddedInterpolant(tITotal, I);  % I(t) [mU/L] piecewise polynomial. Use as function.

% Create SI array and initial time boundaries.
defaultSI = 10.8e-4;

P.SI = ones(P.simDuration(), 1) * defaultSI;           % SI value at each minute in trial.
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
YGI0 = [0.001;  % P1(t=0)        GI Model [mmol]
        0]; 	% P2(t=0)  
    
YID0 = [0;   % ISC(t=0)          ID Model [mU/L]
        0;   % QDFLocal(t=0)
        0;   % QDBLocal(t=0)
        0;   % IDF(t=0)
        0;   % IDB(t=0)
        0;   % QDF(t=0)
        0];  % QDB(t=0)
    
YGC0 = [9;   % Q(t=0)            GC Model [mU/L]
        0;   % GA(t=0)                    [mmol/L]
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
    % deltaG = -int{dGA}*SI + int{dGb}, therefore SI = GA\(Gb-deltaG).
    GA = Y(:, ccGA);         % Coefficient of SI in equation.
    Gb = Y(:, ccGb);         % Added terms in equation.
    
    deltaG  = ppG(tb) - ppG(ta);  % Change in G over interval.  
    
    A = GA;     
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
P.SI(1:length(minuteSI)) = minuteSI;  % [L/mU/min]

fprintf('P%d: SI fit successfully. SI(*1e+3) at %d min intervals = ', ...
        P.patientNum, intervalDuration)
disp(intervalSI'*1e3)

end


% Adapted from fit_SI/FAERIES_integrals.
function [dY] = SIModelODE(t, Y, ppG, ppI, P, Y0)
% Collection of full system ODEs used to compute SI. Use with ode45.
% INPUTS:
%   t   - time at which to evaluate ODE
%   Y   - states [P1; P2; G; I; Q] at time == t
%   ppG - piecewise polynomial interpolant of G over time
%   ppI - piecewise polynomial interpolant of I over time
%   P   - patient struct
%   Y0  - initial conditions of states
% OUTPUT:
%   dY  - derivatives of states at time == t
global GC

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
GFast = P.data.GFast(t);    % Fasting glucose [mmol/L]
GInfusion = P.data.GInfusion(n); % Glucose infusion rate [mmol/min]

% Patient dependent.
d2 = P.d2;

% Derived values.
QTFast  = Q0 + QDF0;
QT = Q + QDF;

%% Computation
% GC Model (reconstructed)
% dG = -dGA*SI + dGb. SI is found with dGA\dGb (linear system).
dQ     = GC.nI/GC.VQ*(I-Q) - GC.nC*Q;
dGA    = (G*QT - GFast*QTFast)/(1 + GC.alphaG*QT);
dGb    = -GC.pg*(G-GFast) + d2/GC.VG*P2 + GInfusion/GC.VG;
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
