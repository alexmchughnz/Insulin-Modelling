function [P, Y] = GCModel(P, flags)
% Function for GC model forward simulation.
% INPUTS:
%   P        - patient struct, must have tArray, qGut(t) and QLocal(t)
%   flags    - string array of flags for optional settings
% OUTPUT:
%   P  - patient struct updated with model results 


if ~exist("flags", "var")
    flags = "";
end

PrintStatusUpdate(P, "Begin solving...")

% % Reshape single parameters into time-sized arrays.
% if length(P.results.nL) == 1
%     P.results.nL = P.results.nL * ones(size(P.results.tArray));
% end

% Set up initial conditions.
[~, vG] = GetData(P.data.G);  % [mmol/L]
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
  
% Pre-generate Uex array to avoid excessive calls.
P.results.UexArray = P.results.Uex(P);
   
% Forward simulate.
enableG = ~ismember("disableG", flags);
[~, Y] = ode45(@GCModelODE, P.results.tArray, Y0, [], P, Y0, enableG);  

% Store results.
P.results.G = Y(:,1);
P.results.I = Y(:,2);
P.results.Q = Y(:,3);

end