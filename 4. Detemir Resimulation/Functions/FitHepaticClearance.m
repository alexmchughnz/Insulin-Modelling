% Function to fit nL and xL.
function P = FitHepaticClearance(P)

global C GC
global DEBUGPLOTS

%% Setup
% Time of reading in sim [min]
% Plasma insulin [pmol/L]
[tITotal, vITotal] = GetSimTime(P, P.data.I);

% Forward simulate ID model for IDF.
P = IDModel(P);
IDF = P.results.IDF; % [mU/L]

% Time and data arrays.
tArray = P.results.tArray;
tI = tITotal;
vI = C.pmol2mU(vITotal) - IDF(tITotal+1);
ppI = griddedInterpolant(tI, vI);  % [mU/L]


%% Analytical Forward Simulation for Q
Q = zeros(length(tArray), 1); %analytical solution for Q

% Consider form of dQ/dt = -kQ*Q + kI*I.
kQ = GC.nC + GC.nI/GC.VQ; % Constant term coefficent of Q - easier to use
kI = GC.nI/GC.VQ;  % Constant term coefficent of I - easier to use

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
% Split data into daily segments, and fit nL per day.
P.results.nL = zeros(size(tArray));

simStart = P.data.simTime(1);
day1 = simStart - timeofday(simStart);  % 00:00 on day 1.
day2 = day1 + 1;                                % 00:00 on day 2.
iiDayEnd = 1 + minutes(day2 - simStart);     % Sim time when day 1 ends.

iiSplits = [iiDayEnd P.data.simDuration()]; % Times of segment ends.
segment = [1 : iiSplits(1)]';
for ii = 1 : length(iiSplits)
    [nL, ~] = FitSegment(P, ppI, Q, tArray, segment);
    P.results.nL(segment) = nL;
    
    % Update time segments (if continuing).
    if ii < length(iiSplits)
        segment = [segment(end)+1 : iiSplits(ii+1)]';
    end
end

% Fit xL for whole time.
segment = [1 : P.data.simDuration()]';
[~, xL] = FitSegment(P, ppI, Q, tArray, segment);
P.results.xL = xL*ones(size(tArray));


%% Debug Plots
DP = DEBUGPLOTS.FitHepaticClearance;

% Retrieve data across segment.
tSegment = [1 : P.data.simDuration()]';
I = ppI(tSegment); % [mU/L]
Q = Q(tSegment);
I0 = I(1);
Uen = P.results.Uen(tSegment); % [mU/min]

% Set coefficients for MLR.
% Consider dI/dt = -kI*I - c1*nL - kIQ*(I-Q) - c2*xL + k
kI = GC.nK;
kIQ = GC.nI./GC.VI;
k = Uen/GC.VI;

cN = cumtrapz(tArray, ...
    I./(1 + GC.alphaI*I));
cX = cumtrapz(tArray, ...
    Uen/GC.VI);
A = [cN cX];
b = I0 - I ...
- kI * cumtrapz(tSegment, I) ...
- kIQ * cumtrapz(tSegment, I-Q) ...
+ cumtrapz(tSegment, k);

LHS = A .* [P.results.nL P.results.xL];


% Forward Simulation of Insulin
if DP.ForwardSim 
    MakeDebugPlot(P, DP);
    
    subplot(2,1,1)
    hold on
    plot(tArray, ppI(tArray))
    simI = -LHS + I0 ...
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

% nL/xL Values per Patient
if DP.nLxL
    persistent sp;
    if (isempty(sp))
        sp = 1;
    end
    
    figure(999)
    
    subplot(2, 3, sp)
    plot(tArray, P.results.nL, 'b')
    title(sprintf("P%d: nL", P.patientNum))
    L = line([iiDayEnd iiDayEnd], ylim);
    L.LineWidth = 1;
    L.Color = 'k';
    
    subplot(2, 3, sp+3)
    plot(tArray, P.results.xL, 'r')
    title(sprintf("P%d: xL", P.patientNum))
    L = line([iiDayEnd iiDayEnd], ylim);
    L.LineWidth = 1;
    L.Color = 'k';
    
    sp = sp + 1;
end

% Equation Terms
if DP.EquationTerms
   MakeDebugPlot(P, DP);
   hold on
   
   plt = plot(tArray, A*[nL; xL], 'b');
   plt.DisplayName = "A*x";
   
   plt = plot(tArray, kI * cumtrapz(tArray, I), 'r');
   plt.DisplayName = "nK * integral(I)";
   
   plt = plot(tArray, kIQ * cumtrapz(tArray, I-Q), 'g');
   plt.DisplayName = "nI/vI * integral(I-Q)";
   
   plt = plot(tArray, cumtrapz(tArray, k), 'm');
   plt.DisplayName = "integral(Uen/vI)";
   
   plt = plot(tArray, I0 - I, 'c');
   plt.DisplayName = "I - I0";
   
   legend()
end
   
% MLR Terms
if DP.MLRTerms
   MakeDebugPlot(P, DP);
   hold on
   
   plt = plot(tArray,  A(:,1));
   plt.DisplayName = "A(column 1) = integral(I / (1 + alphaI*I))";
   
   plt = plot(tArray, A(:,2));
   plt.DisplayName = "A(column 2) = integral(Uen/VI)";
   
   plt = plot(tArray, b);
   plt.DisplayName = "b = I0 - I...";
   
   legend()
end

% Insulin Terms
if DP.InsulinTerms   
   MakeDebugPlot(P, DP);
   hold on
   plot(tArray,  cumtrapz(tArray, kI*I), 'g')
   plot(tArray, cumtrapz(tArray, I./(1 + GC.alphaI*I)))
   legend("integral(nK*I)", "integral(I./(1 + alphaI*I))")
end

end

function [nL, xL] = FitSegment(P, ppI, Q, tArray, segment)
global GC

tSegment = tArray(segment);

% Retrieve data across segment.
I = ppI(segment); % [mU/L]
Q = Q(segment);
I0 = I(1);
Uen = P.results.Uen(segment); % [mU/min]

% Set coefficients for MLR.
% Consider dI/dt = -kI*I - c1*nL - kIQ*(I-Q) - c2*xL + k
kI = 2*GC.nK;
kIQ = 2*GC.nI./GC.VI;
k = Uen/GC.VI;

% Therefore, integrating:
% I(t) - I(t0) = -kI*int{I} - int{c1}*nL - kIQ*int{I-Q} - int{c2}*xL + int{k}
% Renaming cN = int{c1} and cX = int{c2}
% cN*nL + cX*xL = I(t0) - I(t) - kI*int{I} - kIQ*int{I-Q} + int{k}
cN = cumtrapz(tSegment, ...
    I./(1 + GC.alphaI*I));
cX = cumtrapz(tSegment, ...
    Uen/GC.VI);

% Assembling MLR system:
% [cN(t) cX(t)] * (nL; xL) = [b(t)]
A = [cN cX];
b = I0 - I ...
    - kI * cumtrapz(tSegment, I) ...
    - kIQ * cumtrapz(tSegment, I-Q) ...
    + cumtrapz(tSegment, k);

% Fit first segment.
x = A\b;
nL = x(1);
xL = x(2);

% Save first segment values.
lb = 1e-7;  % Lower bound on nL/xL.
nL = max(nL, lb);  % [1/min]
xL = max(xL, lb);  % [1]
end