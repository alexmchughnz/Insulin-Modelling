function [P, Y] = GCModel(P, options)
% Function for GC model forward simulation.
% Pequires qGut(t) and QLocal(t) - must be run AFTER GI and ID models.
% INPUTS:
%   P        - patient struct, must have tArray, qGut(t) and QLocal(t)
%   options  - ode45 options
% OUTPUT:
%   P  - patient struct updated with model results 

if ~exist('options', 'var')
    options = odeset;
end

PrintStatusUpdate(mfilename, P, "Begin solving...")

% Set up initial conditions.
[~, vG] = GetSimTime(P, P.data.G);
[~, vI] = GetIFromITotal(P);  % [mU/L]

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
   
% Forward simulate.
[~, Y] = ode45(@GCModelODE, P.results.tArray, Y0, options, P, Y0);  

% Store results.
P.results.G = Y(:,1);
P.results.I = Y(:,2);
P.results.Q = Y(:,3);

end