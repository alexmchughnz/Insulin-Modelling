function P = IntegralFitnLxL(P, forcenLxL)
% Fits data using MLR to find nL and xL.
% INPUT:
%   P   - patient struct
%   forcenLxL - [nL, xL] to force (for plots)
% OUTPUT:
%   P   - modified patient struct with nL and xL

global CONFIG

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
[tI, vI] = GetData(P.data.I); % [mU/L]

if P.source == "DISST"
    % Need to add 'false' point for improved fitting.    
    [vIBolus, iiBolus] = max(P.results.IBolus);
    tIBolus = tMinutes(iiBolus);
    
    [tI, order] = sort([tI; tIBolus]);
    iiBeforeFakePoint = find(order == length(order)) - 1;
    
    fakeI = vIBolus/GC.VI + vI(iiBeforeFakePoint); % [mU/L]
    
    fakeData = [vI; fakeI];
    vI = fakeData(order);
end

ppI = griddedInterpolant(tI, vI);  % [mU/L]
I = ppI(tArray);
Q = GetAnalyticalInterstitialInsulin(I, P);

%% Iterative Integral Method (pg. 16)
nLArray = [0];
xLArray = [1];
relativeChange = [Inf Inf]; % Change in [nL xL] at each iteration.
tolerance = 0.1/100; % Relative tolerance for convergence.

while any(relativeChange >= tolerance)
    [A, b, IFunc, QFunc] = AssemblenLxLSystem(P, I);
    
    % %         Normalise.
    %         A = A./b
    %         b = b./b
    
    % Solve.
    x = A\b;
    nL = x(1);
    xL = 1 - x(2);
    
    % Calculate deltas.
    nLChange = (nL-nLArray(end))/nLArray(end);
    xLChange = (xL-xLArray(end))/xLArray(end);
    relativeChange = [nLChange xLChange];
    
    % Calculate errors.
    if CONFIG.HIGHDETAIL
        errorRel = sqrt((A*x-b).^2)./b
        errorAbs = ((A*x-b).^2)
        errorAbsSum = sum((A*x-b).^2)
        disp(newline)
    end
    
    % Update arrays.
    nLArray = [nLArray nL];
    xLArray = [xLArray xL];
    
    % Forward simulate to improve I and Q prediction.
    for ii = 1:100
        I = IFunc(nL, xL, I, Q);
        Q = QFunc(I, Q);
    end
end

P.results.integrals.A = A;
P.results.integrals.b = b;


%% Results
% Extract final result.
lb = 1e-7;  % Lower bound on nL/xL.
nL = max(nLArray(end), lb);  % [1/min]
xL = max(xLArray(end), lb);  % [1]

P.results.nL = nL;
P.results.xL = xL;

% Calculate parameter ID metrics.
MeanNormalise = @(data) data ./ mean(data);

CN = A(:, 1);
CX = A(:, 2);
CNNorm = MeanNormalise(CN);
CXNorm = MeanNormalise(CX);
delta2Norm = norm(CNNorm - CXNorm) / length(CN);
P.results.delta2Norm = delta2Norm;

bNorm = MeanNormalise(b);
P.results.delta2NormnL = norm(CNNorm - bNorm) / length(CN);
P.results.delta2NormxL = norm(CXNorm - bNorm) / length(CN);


%% Plotting
plotvars.nLArray = nLArray;
plotvars.xLArray = xLArray;
plotvars.CNNorm = CNNorm;
plotvars.CXNorm = CXNorm;
plotvars.bNorm = bNorm;

MakePlots(P, plotvars);
end


function MakePlots(P, plotvars)
DP = DebugPlots().FitHepaticClearance;
CONST = LoadConstants();

%% Graphical Identifiability Method (Docherty, 2010)
if DP.GraphicalID
    MakeDebugPlot("Graphical Identifiability", P, DP);
    
    [tI, ~] = GetData(P.data.I); % [mU/L]
    tIntegrals = mean([tI(1:end-1), tI(2:end)], CONST.ROWWISE);
    
    plt = plot(tIntegrals, plotvars.CNNorm);
    plt.DisplayName = "$n_L$ coeff.";
    
    plt = plot(tIntegrals, plotvars.CXNorm);
    plt.DisplayName = "$x_L$ coeff.";
    
    plt = plot(tIntegrals, plotvars.bNorm);
    plt.DisplayName = "$b$";
    
    for ii = 1:length(tIntegrals)
        t = tIntegrals(ii);
        
        plt = plot([t t], [plotvars.CNNorm(ii) plotvars.CXNorm(ii)], 'k');
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

%% Convergence
if DP.Convergence
    MakeDebugPlot("Convergence Plot", P, DP);
    plot(plotvars.nLArray, 'b')
    plot(plotvars.xLArray, 'r')
    
    legend("nL", "xL")
end
end