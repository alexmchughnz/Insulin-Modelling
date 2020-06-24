% Function to fit nL and xL.
function P = FitHepaticClearance(P)

global C GC DEBUGPLOT

%% Setup
% Time of reading in sim [min]
% Plasma insulin [pmol/L]
[tITotal, vITotal] = GetSimTime(P, P.data.I);

% Forward simulate ID model for IDF.
P = IDModel(P);
IDF = P.results.IDF; % [mU/L]

% Time and data arrays.
tArray = P.results.tArray;
I = C.pmol2mU(vITotal) - IDF(tITotal+1);
ppI = griddedInterpolant(tITotal, I);  % [mU/L]


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
P.nL = zeros(size(tArray));

day1 = P.simTime(1) - timeofday(P.simTime(1));  % 00:00 on day 1.
day2 = day1 + 1;                                % 00:00 on day 2.
tDayEnd = minutes(day2 - P.simTime(1)) - 1;     % Sim time when day 1 ends.

tSplits = [tDayEnd P.simDuration()]; % Times of segment ends.
segment = [1 tSplits(1)];
for ii = 1 : length(tSplits)
    tSegment = [segment(1) : segment(end)]';
    [nL, ~] = FitSegment(P, ppI, Q, tSegment);
    
    P.nL(tSegment) = nL;
    
    % Update time segments (if continuing).
    if ii < length(tSplits)
        segment = [segment(end)+1, tSplits(ii+1)];
    end
end

% Fit xL for whole time.
tSegment = [1 : P.simDuration()]';
[~, xL] = FitSegment(P, ppI, Q, tSegment);
P.xL = xL*ones(size(tArray));


%% Debug Plots
% if DEBUGPLOT
%     figure()
%     
%     subplot(2,1,1)
%     hold on
%     plot(tArray, ppI(tArray))
%     simI = -A*[nL; xL] + I0 ...
%         - kI * cumtrapz(tArray, GC.nK*I) ...
%         - kIQ * cumtrapz(tArray, I-Q) ...
%         + cumtrapz(tArray, k);
%     plot(tArray, simI)
%     title(sprintf("P%d: I", P.patientNum))
%     legend("interpolated", "simulated")
%     
%     subplot(2,1,2)
%     hold on
%     plot(tArray, Q)
%     title("Q (analytical)")
% end

end


function [nL, xL] = FitSegment(P, ppI, Q, tSegment)
    global GC
    
    % Retrieve data across segment.
    I = ppI(tSegment); % [mU/L]
    Q = Q(tSegment);
    I0 = I(1);
    Uen = P.results.Uen(tSegment); % [mU/min]    
    
    % Set coefficients for MLR.
    % Consider dI/dt = -kI*I - c1*nL - kIQ*(I-Q) - c2*xL + k
    kI = GC.nK;
    kIQ = GC.nI./GC.VI;
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