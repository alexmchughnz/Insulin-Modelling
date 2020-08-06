function P = FitHepaticClearance(P, forcenLxL)
% Fits data using MLR to find nL and xL.
% INPUT:
%   P   - patient struct
%   forcenLxL - [nL, xL] to force (for plots)
% OUTPUT:
%   P   - modified patient struct with nL and xL
global GC
global DEBUGPLOTS

MeanNormalise = @(data) data ./ mean(data);

%% Setup
% Time and data arrays.
tArray = [P.data.simTime(1) : P.data.simTime(end)-1]';  % Minute-wise time range [min]
[tI, vI] = GetIFromITotal(P); % [mU/L]

if P.source == "DISST"
    % Need to add 'false' point for improved fitting.
    [tI, order] = sort([tI; P.data.tIBolus]);
    iiBeforeFakePoint = find(order == length(order)) - 1;
    
    fakeI = P.data.vIBolus/GC.VI + vI(iiBeforeFakePoint); % [mU/L]
    
    fakeData = [vI; fakeI];
    vI = fakeData(order);
end

ppI = griddedInterpolant(tI, vI);  % [mU/L]


%% Analytical Forward Simulation for Q
Q = zeros(length(tArray), 1); %analytical solution for Q

% Consider form of dQ/dt = -cQ*Q + cI*I.
cQ = GC.nC + GC.nI/GC.VQ; % Constant term coefficent of Q - easier to use
cI = GC.nI/GC.VQ;  % Constant term coefficent of I - easier to use

t0 = tArray(1);
I0 = ppI(t0);   % [mU/L]
Q0 = I0/2;      % [mU/L]
Q(1) = Q0;

for ii = 2:length(tArray)
    t = tArray(ii);         % Current time value.
    tSpan = tArray(1:ii);   % Array of time values from t0 to t.
    ISpan = ppI(tSpan);     % Array of I values from t0 to t.
    
    % Analytical solution for the Q equation at time==t.
    % Standard soln. for non-homogenous 1st order ODE (page 17).
    Q(ii) = Q0*exp(-cQ*(t-t0)) ...
        + trapz(tSpan, cI*ISpan.*exp(-cQ*(t - tSpan)));
end

%% Parameter ID of I Equation to find nL/xL (pg. 16)
% Fit nL/xL over segments.
[nL, xL, CN, CX, CParts] = FitSegment(P, ppI, Q, tArray);

if exist('forcenLxL', 'var')
    nL = forcenLxL(1);
    xL = forcenLxL(2);
end

P.results.nL = nL*ones(size(P.results.tArray));
P.results.xL = xL*ones(size(P.results.tArray));

% Save graphical ID statistics.
CNNorm = MeanNormalise(CN);
CXNorm = MeanNormalise(CX);
delta2Norm = norm(CNNorm - CXNorm) / length(CN);
P.results.delta2Norm = delta2Norm;

b = sum(CParts, 2);
bNorm = MeanNormalise(b);
P.results.delta2NormnL = norm(CNNorm - bNorm) / length(CN);
P.results.delta2NormxL = norm(CXNorm - bNorm) / length(CN);

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

LHS = [CN CX] .* [nL xL];
LHS = sum(LHS, 2);
b = sum(CParts, 2);

% Graphical Identifiability Method (Docherty, 2010)
if DP.GraphicalID
    [sampleTimes, ~] = GetIFromITotal(P); % Abuse function to get sample times for any patient.
    
    MakeDebugPlot(P, DP);
    hold on
    
    plt = plot(tArray, CNNorm);
    plt.DisplayName = "$n_L$ coeff.";
    
    plt = plot(tArray, CXNorm);
    plt.DisplayName = "$x_L$ coeff.";
    
    plt = plot(tArray, bNorm);
    plt.DisplayName = "$b$";
    
    for ss = 1:length(sampleTimes)
        t = sampleTimes(ss);
        ii = GetTimeIndex(t, tArray);
        
        plt = plot([t t], [CNNorm(ii) CXNorm(ii)], 'k');
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
    
end

% Forward Simulation of Insulin
if DP.ForwardSim
    I = ppI(tArray);
    cI = GC.nK;
    kIQ = GC.nI./GC.VI;
    k = P.results.Uen/GC.VI + P.data.IBolus(tArray)/GC.VI;
    
    MakeDebugPlot(P, DP);
    
    subplot(2,1,1)
    hold on
    plot(tArray, I)
    
    simI = LHS + I0 ...
        + cI * cumtrapz(tArray, cI*I) ...
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
    ITerm = CParts(:, 1);
    intITerm = CParts(:, 2);
    intIQTerm = CParts(:, 3);
    intUenTerm = CParts(:, 4);
    
    MakeDebugPlot(P, DP);
    hold on
    
    plt = plot(tArray, LHS, 'b');
    plt.DisplayName = "A*x";
    
    plt = plot(tArray, intITerm, 'r');
    plt.DisplayName = "nK * integral(I)";
    
    plt = plot(tArray, intIQTerm, 'g');
    plt.DisplayName = "nI/vI * integral(I-Q)";
    
    plt = plot(tArray, intUenTerm, 'm');
    plt.DisplayName = "-integral((Uen+IBolus)/vI)";
    
    plt = plot(tArray, ITerm, 'c');
    plt.DisplayName = "I - I0";
    
    xlabel("Time [min]")
    ylabel("Mean-normalised value of term [mU/L]")
    legend()
end

% Insulin Terms
if DP.InsulinTerms
    MakeDebugPlot(P, DP);
    hold on
    plot(tArray,  cumtrapz(tArray, cI*I), 'g')
    plot(tArray, cumtrapz(tArray, I./(1 + GC.alphaI*I)))
    legend("integral(nK*I)", "integral(I./(1 + alphaI*I))")
end
end

function [nL, xL, CN, CX, CParts] = FitSegment(P, ppI, Q, tArray)
global GC

% Retrieve data.
iiMinutes = GetTimeIndex(tArray, P.results.tArray);
Uen = P.results.Uen(iiMinutes); % minutewise [mU/min]
I = ppI(tArray); % [mU/L]
I0 = I(1);
Q0 = I0/2;
IBolus = zeros(size(tArray));
for ii = 1:length(tArray)
    t = tArray(ii);
    IBolus(ii) = P.data.IBolus(t);
end

% Set coefficients for MLR.
% Consider dI/dt = kI*I + c1*nL + kIQ*(I-Q) + c2*(1-xL) + k:
kI = -GC.nK;
kIQ = -GC.nI./GC.VI;
k = IBolus/GC.VI;
% Also consider dQ/dt = -cQ*Q + cI*I:
cQ = GC.nC + GC.nI/GC.VQ; % Constant term coefficent of Q - easier to use
cI = GC.nI/GC.VQ;  % Constant term coefficent of I - easier to use

% Perform iterative integral method.
for II = 1:5
    % Integrating I equation:
    % I(t) - I(t0) = kI*int{I} + int{c1}*nL + kIQ*int{I-Q} + int{c2}*(1-xL) + int{k}
    % Renaming CN = int{c1} and CX = int{c2}
    % CN*nL + CX*(1-xL) = I(t) - I(t0) - kI*int{I} - kIQ*int{I-Q} - int{k} := C
    CN = cumtrapz(tArray, ...
        -I./(1 + GC.alphaI*I));
    CX = cumtrapz(tArray, ...
        Uen/GC.VI);
    CParts = [I - I0, ...
        - kI*cumtrapz(tArray, I), ...
        - kIQ*cumtrapz(tArray, I-Q), ...
        - cumtrapz(tArray, k)]; % For analysing each term later.
    C = sum(CParts, 2); % Sum along rows to get column vector.
    
    % Assembling MLR system, evaluating only at sample points:
    % [CN(t) CX(t)] * (nL; 1-xL) = [C(t)]
    sampleTimes = P.data.I.time;
    ppCN = griddedInterpolant(tArray, CN);
    ppCX = griddedInterpolant(tArray, CX);
    ppC  = griddedInterpolant(tArray, C);
    
    A = [ppCN(sampleTimes) ppCX(sampleTimes)];
    b = ppC(sampleTimes);
    
    % Solve.
    x = A\b;
    nL = x(1);
    xL = 1 - x(2);
    
    % Forward simulate to improve I and Q prediction.
    for ii = 1:100
        % I(t) = I(t0) + kI*int{I} + kIQ*int{I-Q} + int{k} + CN*nL + CX*(1-xL)
        I = I0 + kI*cumtrapz(tArray, I) + kIQ*cumtrapz(tArray, I-Q) + cumtrapz(tArray, k) ...
            + CN*nL + CX*(1-xL);
        
        % Q(t) = Q(t0) - cQ*int{Q} + cI*int{I}
        Q = Q0 - cQ*cumtrapz(tArray, Q) + cI*cumtrapz(tArray, I);
    end
end

% Save values.
lb = 1e-7;  % Lower bound on nL/xL.
nL = max(nL, lb);  % [1/min]
xL = max(xL, lb);  % [1]

end