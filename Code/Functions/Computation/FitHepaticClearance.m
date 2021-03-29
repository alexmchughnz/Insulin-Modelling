function P = FitHepaticClearance(P, forcenLxL)
% Fits data using MLR to find nL and xL.
% INPUT:
%   P   - patient struct
%   forcenLxL - [nL, xL] to force (for plots)
% OUTPUT:
%   P   - modified patient struct with nL and xL

DP = DebugPlots().FitHepaticClearance;

PrintStatusUpdate(P, "Fitting nL/xL...")

if exist('forcenLxL', 'var')
    nL = forcenLxL(1);
    xL = forcenLxL(2);
    
    P.results.nL = nL;
    P.results.xL = xL;
    
    return
end

%% Data
tArray = P.results.tArray;

% Plasma Insulin
[tMeas, vI] = GetIFromITotal(P); % [mU/L]

if P.source == "DISST"
    % Need to add 'false' point for improved fitting.
    [vIBolus, iiBolus] = max(P.data.IBolus(tMinutes));
    tIBolus = tMinutes(iiBolus);
    [tMeas, order] = sort([tMeas; tIBolus]);
    iiBeforeFakePoint = find(order == length(order)) - 1;
    
    fakeI = vIBolus/GC.VI + vI(iiBeforeFakePoint); % [mU/L]
    
    fakeData = [vI; fakeI];
    vI = fakeData(order);
end
ppI = griddedInterpolant(tMeas, vI);  % [mU/L]
I = ppI(tArray);

% Interstitial Insulin
ppQ = GetAnalyticalInterstitialInsulin(ppI, P);
Q = ppQ(tArray);

%% Iterative Integral Method (pg. 16)
nLArray = [0];
xLArray = [1];
relativeChange = [Inf Inf]; % Change in [nL xL] at each iteration.
tolerance = 0.1/100; % Relative tolerance for convergence.

while any(relativeChange >= tolerance)
    
    [A, b, LHS, RHS] = AssembleIIntegralSystem(P, tMeas, I, Q);   
    
    % %         Normalise.
    %         A = A./b
    %         b = b./b
    
    % Solve.
    x = A\b;
    nL = x(1);
    xL = 1 - x(2);    
    
    % Calculate errors.
    %         errorRel = sqrt((A*x-b).^2)./b
    %         errorAbs = ((A*x-b).^2)
    %         errorAbsSum = sum((A*x-b).^2)
    %         disp(newline)
    
    nLChange = (nL-nLArray(end))/nLArray(end);
    xLChange = (xL-xLArray(end))/xLArray(end);
    relativeChange = [nLChange xLChange];
    
    nLArray = [nLArray nL];
    xLArray = [xLArray xL];
    
    % Forward simulate to improve I and Q prediction.
    for ii = 1:100
        % I(t) = I(t0) + int{k} + CX*(1-xL) - kI*int{I} + CN*nL - kIQ*int{I-Q}
        I = I0 + cumtrapz(tMinutes, k) + CX*(1-xL) - kI*cumtrapz(tMinutes, I) ...
            + CN*nL - kIQ*cumtrapz(tMinutes, I-Q);
        
        % Q(t) = Q(t0) - cQ*int{Q} + cI*int{I}
        Q = Q0 - cQ*cumtrapz(tMinutes, Q) + cI*cumtrapz(tMinutes, I);
    end
end

%% Results
% Extract final result.
lb = 1e-7;  % Lower bound on nL/xL.
nL = max(nLArray(end), lb);  % [1/min]
xL = max(xLArray(end), lb);  % [1]

P.results.nL = nL;
P.results.xL = xL;

% Calculate stats.
MeanNormalise = @(data) data ./ mean(data);

CNNorm = MeanNormalise(CN);
CXNorm = MeanNormalise(CX);
delta2Norm = norm(CNNorm - CXNorm) / length(CN);
P.results.delta2Norm = delta2Norm;

b = sum(CParts, 2);
bNorm = MeanNormalise(b);
P.results.delta2NormnL = norm(CNNorm - bNorm) / length(CN);
P.results.delta2NormxL = norm(CXNorm - bNorm) / length(CN);

%% Debug Plots
LHS = [CN CX] .* [nL xL];
LHS = sum(LHS, 2);

% Graphical Identifiability Method (Docherty, 2010)
if DP.GraphicalID
    [sampleTimes, ~] = GetIFromITotal(P); % Abuse function to get sample times for any patient.
    
    MakeDebugPlot("Graphical Identifiability", P, DP);
    
    plt = plot(tMinutes, CNNorm);
    plt.DisplayName = "$n_L$ coeff.";
    
    plt = plot(tMinutes, CXNorm);
    plt.DisplayName = "$x_L$ coeff.";
    
    plt = plot(tMinutes, bNorm);
    plt.DisplayName = "$b$";
    
    for ss = 1:length(sampleTimes)
        t = sampleTimes(ss);
        ii = GetTimeIndex(t, tMinutes);
        
        plt = plot([t t], [CNNorm(ii) CXNorm(ii)], 'k');
        plt.Marker = 'o';
        plt.MarkerFaceColor = 'auto';
        plt.HandleVisibility = 'off';
    end
    plt.HandleVisibility = 'on';
    plt.DisplayName = "Samples";
    
    xlabel("Time [min]")
    ylabel("Mean-normalised integral value")
    legend()    
end

% Forward Simulation of Insulin
if DP.ForwardSim
    I = ppI(tMinutes);
    kI = GC.nK(P);
    kIQ = GC.nI(P)./GC.VI;
    k = P.results.Uen/GC.VI + P.data.IBolus(tMinutes)/GC.VI;
    
    MakeDebugPlot("Insulin Simulation", P, DP);
    
    subplot(2,1,1)
    hold on
    plot(tMinutes, I)
    
    simI = LHS + I0 ...
        + kI * cumtrapz(tMinutes, kI*I) ...
        + kIQ * cumtrapz(tMinutes, I-Q) ...
        + cumtrapz(tMinutes, k);
    plot(tMinutes, simI)
    
    xlabel("Time [min]")
    ylabel("Plasma insulin, I [mU/L]")
    legend("interpolated", "simulated")
    
    subplot(2,1,2)
    hold on
    plot(tMinutes, Q)
    ylabel("Interstitial insulin, Q [mU/L]")
end

% Equation Terms
if DP.EquationTerms
    ITerm = CParts(:, 1);
    intITerm = CParts(:, 2);
    intIQTerm = CParts(:, 3);
    intUenTerm = CParts(:, 4);
    
    MakeDebugPlot("Equation Terms", P, DP);
    
    plt = plot(tMinutes, LHS, 'b');
    plt.DisplayName = "A*x";
    
    plt = plot(tMinutes, intITerm, 'r');
    plt.DisplayName = "nK * integral(I)";
    
    plt = plot(tMinutes, intIQTerm, 'g');
    plt.DisplayName = "nI/vI * integral(I-Q)";
    
    plt = plot(tMinutes, intUenTerm, 'm');
    plt.DisplayName = "-integral((Uen+IBolus)/vI)";
    
    plt = plot(tMinutes, ITerm, 'c');
    plt.DisplayName = "I - I0";
    
    xlabel("Time [min]")
    ylabel("Mean-normalised value of term [mU/L]")
    legend()
end

% Insulin Terms
if DP.InsulinTerms
    MakeDebugPlot("Insulin Terms", P, DP);
    plot(tMinutes,  cumtrapz(tMinutes, cI*I), 'g')
    plot(tMinutes, cumtrapz(tMinutes, I./(1 + GC.alphaI*I)))
    legend("integral(nK*I)", "integral(I./(1 + alphaI*I))")
end

% Convergence
if DP.Convergence
    MakeDebugPlot("Convergence Plot", P, DP);
    plot(nLArray, 'b')
    plot(xLArray, 'r')
    
    legend("nL", "xL")
end
end