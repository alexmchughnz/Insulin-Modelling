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
    iiBounds = [length(P.results.tArray)];
    
elseif isequal(method, 'daily')
    % Split data into daily segments, and fit nL per day.
    simStart = P.data.simTime(1);
    day1 = simStart - timeofday(simStart);    % 00:00 on day 1.
    day2 = day1 + 1;                          % 00:00 on day 2.
    iiDayEnd = 1 + minutes(day2 - simStart);  % Sim time when day 1 ends.
    
    iiBounds = [iiDayEnd length(P.results.tArray)];
    
elseif isequal(method, 'peaks')
    % Fits nL at specified peaks.
    % Fits are of width 'window', centered at 'peaks' + 'delay'.
    peaks = P.data.tIPeaks;
    window = 60;  % [min]
    delay = 5;   % [min]
    
    peakBounds = [peaks + delay - window/2; ...
        peaks + delay + window/2];
    peakBounds = peakBounds(:).';  % Collapse to row vector of bounds around peaks.
    
    iiBounds = [peakBounds length(P.results.tArray)];
end

if isequal(method, 'fixed')
    nL = arg(1);
    xL = arg(2);
    
    P.results.nL = nL*ones(size(P.results.nL));
    P.results.xL = xL*ones(size(P.results.xL));
    
    P.results.nLxLFitBounds = [1 length(P.results.tArray)];
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
    segment = [1 : length(P.results.tArray)]';
    [~, xL] = FitSegment(P, ppI, Q, tArray, segment);
    P.results.xL = xL*ones(size(P.results.xL));
end

% Save graphical ID statistics.
cN = A(:, 1);
cX = A(:, 2);
cNNorm = MeanNormalise(cN);
cXNorm = MeanNormalise(cX);
delta2Norm = norm(cNNorm - cXNorm) / length(cN);

b = sum(bParts, 2);
bNorm = MeanNormalise(b);
delta2NormnL = norm(cNNorm - bNorm) / length(cN);
delta2NormxL = norm(cXNorm - bNorm) / length(cN);

P.results.delta2Norm = delta2Norm;

%% Debug Plots
DP = DEBUGPLOTS.FitHepaticClearance;

% nL/xL Values per Patient
if DP.nLxL
    MakeDebugPlot(P, DP);
    
    subplot (2, 1, 1)
    plot(tArray, P.results.nL, 'b')
    for ii = 1:length(iiBounds)
        split = iiBounds(ii);
        L = line([split split], ylim);
        L.LineWidth = 0.5;
        L.Color = 'k';
    end
    
    title(sprintf("P%d: nL", P.patientNum))
    ylabel("$n_L$ [1/min]")
    
    
    subplot(2, 1, 2)
    plot(tArray, P.results.xL, 'r')
    
    title(sprintf("P%d: xL", P.patientNum))
    xlabel("Time [min]")
    ylabel("$x_L$ [-]")
    
end

if ~isequal(method, 'fixed')
    LHS = dot(A, [P.results.nL P.results.xL], 2);
    b = sum(bParts, 2);
    
    % Graphical Identifiability Method (Docherty, 2010)
    if DP.GraphicalID
        [sampleTimes, ~] = GetIFromITotal(P); % Abuse function to get sample times for any patient.
        
        MakeDebugPlot(P, DP);
        hold on
        
        plt = plot(tArray, cNNorm);
        plt.DisplayName = "$n_L$ coeff.";
        
        plt = plot(tArray, cXNorm);
        plt.DisplayName = "$x_L$ coeff.";
        
        plt = plot(tArray, bNorm);
        plt.DisplayName = "$b$";
        
        for ss = 1:length(sampleTimes)
            t = sampleTimes(ss);
            ii = GetTimeIndex(t, tArray);
            
            plt = plot([t t], [cNNorm(ii) cXNorm(ii)], 'k');
            plt.Marker = 'o';
            plt.MarkerFaceColor = 'auto';
            plt.HandleVisibility = 'off';
        end
        plt.HandleVisibility = 'on';
        plt.DisplayName = "Samples";
        
        title(sprintf("%s: Graphical Identifiability", P.patientCode))
        xlabel("Time [min]")
        ylabel("Mean-normalised integral value")
        legend('Location', 'northwest')
        
        %         SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.20;
        %         label = sprintf("$||\\Delta_{n_L}||_2$ = %.4g\n $||\\Delta_{x_L}||_2$ = %.4g", ...
        %             delta2NormnL, delta2NormxL);
        %         text(SE(1), SE(2), label);
        
    end
    
    % Forward Simulation of Insulin
    if DP.ForwardSim
        I = ppI(tArray);
        kI = GC.nK;
        kIQ = GC.nI./GC.VI;
        k = P.results.Uen/GC.VI + P.data.IBolus(tArray)/GC.VI;
        
        MakeDebugPlot(P, DP);
        
        subplot(2,1,1)
        hold on
        plot(tArray, I)
        
        simI = LHS + I0 ...
            + kI * cumtrapz(tArray, kI*I) ...
            + kIQ * cumtrapz(tArray, I-Q) ...
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
    % Equation Terms
    if DP.EquationTerms
        ITerm = bParts(:, 1);
        intITerm = bParts(:, 2);
        intIQTerm = bParts(:, 3);
        intUenTerm = bParts(:, 4);
        
        MakeDebugPlot(P, DP);
        hold on
        
        plt = plot(tArray, MeanNormalise(LHS), 'b');
        plt.DisplayName = "A*x";
        
        plt = plot(tArray, MeanNormalise(intITerm), 'r');
        plt.DisplayName = "nK * integral(I)";
        
        plt = plot(tArray, MeanNormalise(intIQTerm), 'g');
        plt.DisplayName = "nI/vI * integral(I-Q)";
        
        plt = plot(tArray, MeanNormalise(intUenTerm), 'm');
        plt.DisplayName = "-integral((Uen+IBolus)/vI)";
        
        plt = plot(tArray, MeanNormalise(ITerm), 'c');
        plt.DisplayName = "I - I0";
        
        for ii = 1:length(iiBounds)
        split = tArray(P.results.nLxLFitBounds(ii));
            L = line([split split], ylim);
            L.LineWidth = 0.5;
            L.Color = 'k';
            L.HandleVisibility = 'off';
        end
        
        xlabel("Time [min]")
        ylabel("Mean-normalised value of term [mU/L]")
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
        split = tArray(P.results.nLxLFitBounds(ii));
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
        split = tArray(P.results.nLxLFitBounds(ii));
            L = line([split split], ylim);
            L.LineWidth = 0.5;
            L.Color = 'k';
            L.HandleVisibility = 'off';
        end
        
        xlabel("Time [min]")
        ylabel("Condition number of A (over segment)")
    end
end
end

function [nL, xL, A, bParts, condA] = FitSegment(P, ppI, Q, tArray, segment)
global GC
INTERVALS = true;

tSegment = tArray(segment);

% Retrieve data across segment.
Uen = P.results.Uen(segment); % [mU/min]
I = ppI(tSegment); % [mU/L]
Q = Q(segment);
I0 = I(1);
IBolus = zeros(size(segment));
for ii = 1:length(segment)
    t = segment(ii);
    IBolus(ii) = P.data.IBolus(t);
end

% Set coefficients for MLR.
% Consider dI/dt = kI*I + c1*nL + kIQ*(I-Q) + c2*xL + k
kI = -GC.nK;
kIQ = -GC.nI./GC.VI;
k = Uen/GC.VI + IBolus/GC.VI;

% Therefore, integrating:
% I(t) - I(t0) = kI*int{I} + int{c1}*nL + kIQ*int{I-Q} + int{c2}*xL + int{k}
% Renaming cN = int{c1} and cX = int{c2}
% cN*nL + cX*xL = I(t) - I(t0) - kI*int{I} - kIQ*int{I-Q} - int{k}
cN = cumtrapz(tSegment, ...
    -I./(1 + GC.alphaI*I));
cX = cumtrapz(tSegment, ...
    -Uen/GC.VI);

% Assembling MLR system:
% [cN(t) cX(t)] * (nL; xL) = [b(t)]
A = [cN cX];
bParts = [I - I0, ...
    - kI * cumtrapz(tSegment, I), ...
    - kIQ * cumtrapz(tSegment, I-Q), ...
    - cumtrapz(tSegment, k)];

if INTERVALS
    % Evaluate at specified intervals.
    interval = 2; %[min]
    dt = tArray(2) - tArray(1);
    n = round(interval/dt); % How many indices to get the next value?
    
    A = A(1+n:n:end, :) - A(1:n:end-n, :);
    bParts = bParts(1+n:n:end, :) - bParts(1:n:end-n, :);
end

% Solve.
condA = cond(A);
b = sum(bParts, 2); % Sum along rows.
x = A\b;
nL = x(1);
xL = x(2);

if INTERVALS
    % Reshape A and b at the resolution of tArray.
    A = repelem(A, n, 1);
    bParts = repelem(bParts, n, 1);
    
    oldlen = size(A, 1);
    newlen = length(tArray);
    newrows = newlen - oldlen;
    
    A(oldlen+1:newlen, :) = repmat(A(end, :), newrows, 1);
    bParts(oldlen+1:newlen, :) = repmat(bParts(end, :), newrows, 1);
end

% Save first segment values.
lb = 1e-7;  % Lower bound on nL/xL.
nL = max(nL, lb);  % [1/min]
xL = max(xL, lb);  % [1]

end