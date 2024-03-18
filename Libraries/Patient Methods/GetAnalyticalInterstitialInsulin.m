function Q = GetAnalyticalInterstitialInsulin(I, P)
% Given a plasma insulin (I) profile, returns interstitial insulin (Q) profile.
% INPUTS:
%   I - plasma insulin time profile, I(t)
%   P   - patient struct
% OUTPUT:
%   Q - interstitial insulin time profile, Q(t)

GC = P.parameters.GC;

% Setup.
tArray = P.results.tArray;
Q = zeros(length(tArray), 1); % Analytical solution for Q

% Consider form of dQ/dt = -cQ*Q + cI*I.
cQ = GC.nC + GC.nI/GC.VQ; % Constant term coefficent of Q
cI = GC.nI/GC.VQ;         % Constant term coefficent of I

% Initial values.
t0 = tArray(1);
I0 = I(1);   % [mU/L]
Q0 = I0/2;      % [mU/L]
Q(1) = Q0;

for ii = 2:length(Q)
    t = tArray(ii);         % Current time value.
    tSpan = tArray(1:ii);   % Array of time values from t0 to t.
    ISpan = I(1:ii);     % Array of I values from t0 to t.
    
    % Analytical solution for the Q equation at time==t.
    % Standard soln. for non-homogenous 1st order ODE (page 17).
    Q(ii) = Q0*exp(-cQ*(t-t0)) ...
        + trapz(tSpan, cI*ISpan.*exp(-cQ*(t - tSpan)));
end

end

