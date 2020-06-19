% Function to fit nL and xL.
function P = FitHepaticClearance(P)

global C GC DEBUGPLOT

%% Setup
% Time of reading in sim [min]
% Plasma insulin [pmol/L]
[tI, vI] = GetSimTime(P, P.data.I);


% Time and data arrays.
tArray = P.results.tArray;
ppI = griddedInterpolant(tI, C.pmol2mU(vI)); % [mU/L]


%% Analytical Forward Simulation for Q
Q = zeros(length(tArray), 1); %analytical solution for Q

% Consider form of dQ/dt = -kQ*Q + kI*I.
kQ = GC.nC + GC.nI/GC.VQ; % Constant term coefficent of Q - easier to use
kI = GC.nI/GC.VQ;  % Constant term coefficent of Q - easier to use

t0 = tArray(1);
I0 = ppI(t0);   % [mU]
Q0 = I0/2;      % [mU]
Q(1) = Q0;

for ii = 2:length(tArray)
    t = tArray(ii);         % Current time value.
    tSpan = tArray(1:ii);   % Array of time values from t0 to t.
    ISpan = ppI(tSpan);     % Array of I values from t0 to t.
    
    % Analytical solution for the Q equation at time==t.
    % Standard soln. for non-homogenous 1st order ODE (page 17).
    Q(ii) = Q0*exp(-kQ*(t-t0)) ...
                    + trapz(tSpan, kI*ISpan.*exp(-kQ*(t - tSpan)));
end

%% Parameter ID of I Equation to find nL/xL (pg. 16)
I = ppI(tArray); % [mU/L]
I0 = ppI(t0);
Uen = P.results.Uen; % [mU/min]

% Set coefficients for MLR.
% Consider dI/dt = -kI*I - c1*nL - kIQ*(I-Q) - c2*xL + k
kI = GC.nK;
kIQ = GC.nI./GC.VQ;
k = Uen/GC.VI;

% Therefore, integrating:
% I(t) - I(t0) = -kI*int{I} - int{c1}*nL - kIQ*int{I-Q} - int{c2}*xL + int{k}
% Renaming cN = int{c1} and cX = int{c2}
% cN*nL + cX*xL = I(t0) - I(t) - kI*int{I} - kIQ*int{I-Q} + int{k}
cN = cumtrapz(tArray, ...
              I./(1 + GC.alphaI*I));
cX = cumtrapz(tArray, ...
              Uen/GC.VI);
          
% Assembling MLR system:
% [cN(t) cX(t)] * (nL; xL) = [b(t)]
A = [cN cX];
b = I0 - I ...
       - kI * cumtrapz(tArray, I) ...
       - kIQ * cumtrapz(tArray, I-Q) ...
       + cumtrapz(tArray, k); 

% Solve and save.
x = A\b;
nL = x(1);
xL = x(2);

lb = 1e-7; %"lower bound"
nL = max(nL, lb);
xL = max(xL, lb);

P.nL = nL;  % [1/min]
P.xL = xL;  % [1]


%% Debug Plots
if DEBUGPLOT
   figure()
   
   subplot(2,1,1)
   hold on
   plot(tArray, ppI(tArray))
   simI = -A*[nL; xL] + I0 ...
                      - kI * cumtrapz(tArray, GC.nK*I) ...
                      - kIQ * cumtrapz(tArray, I-Q) ...
                      + cumtrapz(tArray, k);
   plot(tArray, simI)
   title(sprintf("P%d: I", P.patientNum))
   legend("interpolated", "simulated")
   
   subplot(2,1,2)
   hold on
   plot(tArray, Q)
   title("Q (analytical)")
end



end