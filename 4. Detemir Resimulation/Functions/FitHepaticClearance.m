function P = FitHepaticClearance(P, method, arg)
% Fits data using MLR to find nL and xL.
% INPUT:
%   P   - patient struct
% OUTPUT:
%   P      - modified patient struct with nL and xL
%   method - 'single' to fit one nL value for entire period
%            'daily' to fit one nL value per simulation day
%            'peaks' to fit at manually-specified locations
%            'fixed' to force [nL xL]
%   arg    - with 'fixed', [nL xL] values to force
global GC
global DEBUGPLOTS

MeanNormalise = @(data) data ./ mean(data);

%% Setup
% Time and data arrays.
[tI, vI] = GetIFromITotal(P);      % [mU/L]
ppI = griddedInterpolant(tI, vI);  % [mU/L]
tArray = P.results.tArray;         % [min]


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
P.results.nL = zeros(size(tArray));
P.results.xL = zeros(size(tArray));
iiBounds = [];        % Times of fit segment ends.

if isequal(method, 'single')
    % Single value for whole simulation period.
    iiBounds = [P.data.simDuration()];
    
elseif isequal(method, 'daily')
    % Split data into daily segments, and fit nL per day.
    simStart = P.data.simTime(1);
    day1 = simStart - timeofday(simStart);    % 00:00 on day 1.
    day2 = day1 + 1;                          % 00:00 on day 2.
    iiDayEnd = 1 + minutes(day2 - simStart);  % Sim time when day 1 ends.
    
    iiBounds = [iiDayEnd P.data.simDuration()];
    
elseif isequal(method, 'peaks')
    % Fits nL at specified peaks.
    % Fits are of width 'window', centered at 'peaks' + 'delay'.
    if P.patientNum == 1
        peaks = [115 290 730 1525 1760 2110];
    elseif P.patientNum == 3
        peaks = [125 670 1542 2153];
    elseif P.patientNum == 4
        peaks = [152 715 1540 2155];
    end
    window = 60;  % [min]
    delay = 5;   % [min]
    
    peakBounds = [peaks + delay - window/2; ...
        peaks + delay + window/2];
    peakBounds = peakBounds(:).';  % Collapse to row vector of bounds around peaks.
    
    iiBounds = [peakBounds P.data.simDuration()];
end

if isequal(method, 'fixed')
    nL = arg(1);
    xL = arg(2);
    
    P.results.nL = nL*ones(size(P.results.nL));
    P.results.xL = xL*ones(size(P.results.xL));
else
    % Fit nL over segments.
    A = zeros(length(tArray), 2);
    bParts = zeros(length(tArray), 4);
    condA = zeros(length(tArray), 1);
    segment = [1 : iiBounds(1)]';
    for ii = 1 : length(iiBounds)
        [nL, ~, segA, segbParts, segcondA] = FitSegment(P, ppI, Q, tArray, segment);
        P.results.nL(segment) = nL;
        A(segment, :) = segA;
        bParts(segment, :) = segbParts;
        condA(segment) = segcondA;
        
        % Update time segments (if continuing).
        if ii < length(iiBounds)
            segment = [segment(end)+1 : iiBounds(ii+1)]';
        end
    end
    P.results.nLxLFitBounds = iiBounds;
    
    % Fit xL for whole time.
    segment = [1 : P.data.simDuration()]';
    [~, xL] = FitSegment(P, ppI, Q, tArray, segment);
    P.results.xL = xL*ones(size(P.results.xL));
end


%% Debug Plots
DP = DEBUGPLOTS.FitHepaticClearance;

LHS = dot(A, [P.results.nL P.results.xL], 2);
b = sum(bParts, 2);

% Forward Simulation of Insulin
if DP.ForwardSim
    I = ppI(tArray);
    kI = GC.nK;
    kIQ = GC.nI./GC.VI;
    k = P.results.Uen/GC.VI;
    
    MakeDebugPlot(P, DP);
    
    subplot(2,1,1)
    hold on
    plot(tArray, I)
    simI = -LHS + I0 ...
        - kI * cumtrapz(tArray, kI*I) ...
        - kIQ * cumtrapz(tArray, I-Q) ...
        + cumtrapz(tArray, k);
    plot(tArray, simI)
    
    xlabel("Time [min]")
    ylabel("Plasma insulin, I [mU/L]")
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
    for ii = 1:length(iiBounds)
        split = iiBounds(ii);
        L = line([split split], ylim);
        L.LineWidth = 0.5;
        L.Color = 'k';
    end
    
    title(sprintf("P%d: nL", P.patientNum))
    ylabel("$n_L$ [1/min]")
    
    
    subplot(2, 3, sp+3)
    plot(tArray, P.results.xL, 'r')
    
    title(sprintf("P%d: xL", P.patientNum))
    xlabel("Time [min]")
    ylabel("$x_L$ [-]")
    
    sp = sp + 1;
end

% Equation Terms
if DP.EquationTerms
    ITerm = bParts(:, 1);
    intITerm = bParts(:, 2);
    intIQTerm = bParts(:, 3);
    UenTerm = bParts(:, 4);
    
    MakeDebugPlot(P, DP);
    hold on
    
    plt = plot(tArray, MeanNormalise(LHS), 'b');
    plt.DisplayName = "A*x";
    
    plt = plot(tArray, MeanNormalise(intITerm), 'r');
    plt.DisplayName = "-nK * integral(I)";
    
    plt = plot(tArray, MeanNormalise(intIQTerm), 'g');
    plt.DisplayName = "-nI/vI * integral(I-Q)";
    
    plt = plot(tArray, MeanNormalise(UenTerm), 'm');
    plt.DisplayName = "-integral(Uen/vI)";
    
    plt = plot(tArray, MeanNormalise(ITerm), 'c');
    plt.DisplayName = "I0 - I";
    
    for ii = 1:length(iiBounds)
        split = iiBounds(ii);
        L = line([split split], ylim);
        L.LineWidth = 0.5;
        L.Color = 'k';
        L.HandleVisibility = 'off';
    end
    
    xlabel("Time [min]")
    ylabel("Mean-normalised integral of term [mU/L]")
    legend()
end

% MLR Terms
if DP.MLRTerms
    MakeDebugPlot(P, DP);
    hold on
    
    plt = plot(tArray,  MeanNormalise(A(:,1)));
    plt.DisplayName = "A(column 1) = integral(I / (1 + alphaI*I))";
    
    plt = plot(tArray, MeanNormalise(A(:,2)));
    plt.DisplayName = "A(column 2) = integral(Uen/VI)";
    
    plt = plot(tArray, MeanNormalise(b));
    plt.DisplayName = "b = I0 - I...";
    
    for ii = 1:length(iiBounds)
        split = iiBounds(ii);
        L = line([split split], ylim);
        L.LineWidth = 0.5;
        L.Color = 'k';
        L.HandleVisibility = 'off';
    end
    
    xlabel("Time [min]")
    ylabel("Mean-normalised integral of term [mU/L]")
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

% Condition Number
if DP.ConditionNumber
    MakeDebugPlot(P, DP);
    hold on
    
    plot(tArray,  condA);
    
    for ii = 1:length(iiBounds)
        split = iiBounds(ii);
        L = line([split split], ylim);
        L.LineWidth = 0.5;
        L.Color = 'k';
        L.HandleVisibility = 'off';
    end
    
    xlabel("Time [min]")
    ylabel("Condition number of A (over segment)")
end
end

function [nL, xL, A, bParts, condA] = FitSegment(P, ppI, Q, tArray, segment)
global GC

tSegment = tArray(segment);

% Retrieve data across segment.
I = ppI(segment); % [mU/L]
Q = Q(segment);
I0 = I(1);
Uen = P.results.Uen(segment); % [mU/min]

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
bParts = [I0 - I, ...
    - kI * cumtrapz(tSegment, I), ...
    - kIQ * cumtrapz(tSegment, I-Q), ...
    + cumtrapz(tSegment, k)];
b = sum(bParts, 2); % Sum along rows.

% Fit first segment.
x = A\b;
nL = x(1);
xL = x(2);

% Save first segment values.
lb = 1e-7;  % Lower bound on nL/xL.
nL = max(nL, lb);  % [1/min]
xL = max(xL, lb);  % [1]
condA = cond(A);
end