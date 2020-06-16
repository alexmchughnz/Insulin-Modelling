function [P] = GCModel(P, options)
% Function for GC model forward simulation.
% Pequires qGut(t) and QLocal(t) - must be run AFTER GI and ID models.
% INPUTS:
%   P        - patient struct, must have tArray, qGut(t) and QLocal(t)
%   options  - ode45 options
% OUTPUT:
%   P  - patient struct updated with model results 

global C

if ~exist('options', 'var')
    options = odeset;
end

% Set up initial conditions.
[t, G] = GetSimTime(P, P.data.G{3});
inSimTime = (0 <= t) & (t <= P.simDuration()); % [logical]
G = G(inSimTime);

G0 = G(1);
I0 = C.pmol2mIU(P.data.I.value(1)); % [pmol/L] -> [mIU/L]
Q0 = I0/2;  % Subcut Q assumed to be half of plasma I at t=0.

Y0 = [G0;   
      I0; 
      Q0];
    
% Forward simulate.
[~, Y] = ode45(@GCModelODE, P.results.tArray, Y0, options, P, Y0);  

% Store results.
P.results.G = Y(:,1);
P.results.I = Y(:,2);
P.results.Q = Y(:,3);

end