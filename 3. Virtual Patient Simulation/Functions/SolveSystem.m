function results = SolveSystem(V)
% Fits data to find SI over time for a patient.
% INPUTS:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with SI

%% Setup
options = odeset('RelTol',1e-5, ...
                 'AbsTol',1e-4);
 
tArray = (0 : 0.1 : 150)';
results.tArray = tArray;

%% GI Model
% Set up initial conditions.
q0Sto1     = 35;
q0Sto2     = 0;
q0Gut      = 0;

YGI0 = [q0Sto1;
        q0Sto2;
        q0Gut];
      
% Forward simulate.
[~, YGI] = ode45(@GIModelODE, tArray, YGI0, options);  

% Store results.
results.qSto1 = YGI(:,1);
results.qSto2 = YGI(:,2);
results.qGut  = YGI(:,3);

%% IN Model
% Set up initial conditions.
ISC0    = 2000;
QLocal0 = 0;

YIN0 = [ISC0;  
        QLocal0];

% Forward simulate.
[~, YIN] = ode45(@IDModelODE, tArray, YIN0, options);      
      
% Store results.
results.ISC    = YIN(:,1);
results.QLocal = YIN(:,2);

%% GC Model
% Set up initial conditions.
G0 = 4.8;
I0 = 13.4185;
Q0 = 7.4255;

YGC0 = [G0;   % G(t=0)
        I0;   % I(t=0)
        Q0];  % Q(t=0)

% Forward simulate.
[~, YGC] = ode45(@GCModelODE, tArray, YGC0, options, ...
                 results.tArray, results.qGut, results.QLocal, ...
                 V);

% Store results.
results.G = YGC(:,1);
results.I = YGC(:,2);
results.Q = YGC(:,3);

end
