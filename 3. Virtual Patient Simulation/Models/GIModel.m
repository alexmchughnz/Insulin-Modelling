function [R] = GIModel(R, O, V)
% Function for GI model forward simulation.
% INPUTS:
%   R  - results struct, must have tArray
%   O  - ode45 options
%   V  - parameter variants
% OUTPUT:
%   R  - results struct updated with model results 

global GI

% Set up initial conditions.
Y0 = [GI.q0Sto1;
      GI.q0Sto2;
      GI.q0Gut];
  
% Forward simulate.
[~, Y] = ode45(@GIModelODE, R.tArray, Y0, O);  

% Store results.
R.qSto1 = Y(:,1);
R.qSto2 = Y(:,2);
R.qGut  = Y(:,3);

end


function [dY] = GIModelODE(t, Y)
% ODE for GI model. Use with ode45.
% INPUTS:
%   t   - time at which to evaluate ODE
%   Y   - states [qSto1; qSto2; qGut] at time == t
%   P   - patient struct
% OUTPUT:
%   dY  - derivatives of states at time == t

global GI

%% Input
qSto1     = Y(1);
qSto2     = Y(2);
qGut      = Y(3);

%% Variables
% Patient dependent.
qSto = qSto1 + qSto2;
kEmpt = GetStomachEmptyingRate(qSto);

% Time dependent.
isActive = (t >= 0) && (t < 1); % Glucose input over 1st minute.
D = (isActive)*GI.D; 

%% Computation
dqSto1 = -GI.k21*qSto1 + D;
dqSto2 = -kEmpt*qSto2 + GI.k21*qSto1;
dqGut  = -GI.kAbs*qGut + kEmpt*qSto2;

%% Output
dY = [dqSto1;
      dqSto2;
      dqGut];

end
