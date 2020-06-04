function [R] = INModel(R, O, V)
% Function for IN model forward simulation.
% INPUTS:
%   R  - results struct, must have tArray
%   O  - ode45 options
%   V  - parameter variants
% OUTPUT:
%   R  - results struct updated with model results 

global IN

% Set up initial conditions.
Y0 = [IN.ISC0;  
      IN.QLocal0];
  
% Forward simulate.
[~, Y] = ode45(@INModelODE, R.tArray, Y0, O);  

% Store results.
R.ISC    = Y(:,1);
R.QLocal = Y(:,2);

end


function [dY] = INModelODE(t, Y)
% ODE for IN model. Use with ode45.
% INPUTS:
%   t   - time at which to evaluate ODE
%   Y   - states [ISC; QLocal] at time == t
%   P   - patient struct
% OUTPUT:
%   dY  - derivatives of states at time == t

global IN

%% Input
ISC    = Y(1);
QLocal = Y(2);

%% Variables
% Time dependent.
isActive = (t >= IN.TBolus) && (t < IN.TBolus+1); % Insulin bolus over 1 minute.
IBolus = (isActive)*IN.IBolus;

%% Computation
dISC    = -IN.k2*ISC + IBolus;
dQLocal = -IN.k3*QLocal ...
              + IN.k2*ISC ...
              - IN.kdi*QLocal;

%% Output
dY = [dISC;
      dQLocal];

end

