function P = SolveSystem(P)
% Fits data to find SI over time for a patient.
% INPUTS:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with SI

global C

%% Setup
options = odeset('RelTol',1e-5, ...
                 'AbsTol',1e-4);
 
t = (0 : P.simDuration-1)';  % Simulation time array [min]

P.results.time = P.simTime(1) + t/24/60;  % Time of results [datetime]

%% GI Model
% Set up initial conditions.
YGI0 = [0.001;  % qSto1(t=0)
        0; 	    % qSto2(t=0)
        0];     % qGut(t=0)
      
% Forward simulate.
[~, YGI] = ode45(@GIModelODE, t, YGI0, options, P);  

% Store results.
P.results.qSto1 = YGI(:,1);
P.results.qSto2 = YGI(:,2);
P.results.qGut  = YGI(:,3);

%% ID Model
% Set up initial conditions.
YID0 = [0;  % ISC(t=0)
        0;  % QDFLocal(t=0)
        0;  % QDBLocal(t=0)
        0;  % IDF(t=0)
        0;  % IDB(t=0)
        0;  % QDF(t=0)
        0]; % QDB(t=0)

% Forward simulate.
[~, YID] = ode45(@IDModelODE, t, YID0, options, P);      
      
% Store results.
P.results.IDH      = YID(:,1);
P.results.QDFLocal = YID(:,2);
P.results.QDBLocal = YID(:,3);
P.results.IDF      = YID(:,4);  % NOTE: Removed a factor of 18.
P.results.IDB      = YID(:,5);  % "
P.results.QDF      = YID(:,6);  % "
P.results.QDB      = YID(:,7);  % "

%% GC Model
% Set up initial conditions.
G0Indices = (P.G{3}.time == P.simTime(1));  % Start of G measurements [indices]
G0 = P.G{3}.value(G0Indices);
I0 = C.pmol2mIU(P.I.value(1)); % [pmol/L] -> [mIU/L]

YGC0 = [G0;  % G(t=0)
        I0;  % I(t=0)
        0];  % Q(t=0)

% Forward simulate.
[~, YGC] = ode45(@GCModelODE, t, YGC0, options, P);

% Store results.
P.results.G = YGC(:,1);
P.results.I = YGC(:,2);
P.results.Q = YGC(:,3);

end
