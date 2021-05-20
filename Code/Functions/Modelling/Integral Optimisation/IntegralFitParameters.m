function P = IntegralFitParameters(P, integralSystemFunc)
% Fits data using MLR to fit for parameters.
% INPUT:
%   P   - patient struct
%   forcenLxL - [nL, xL] to force (for plots)
% OUTPUT:
%   P   - modified patient struct with nL and xL

global CONFIG

PrintStatusUpdate(P, "Fitting nL/xL...")

if exist('forcenLxL', 'var')
    param1 = forcenLxL(1);
    param2 = forcenLxL(2);
    
    P.results.nL = param1;
    P.results.xL = param2;
    
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
param1Array = [0];
param2Array = [1];
relativeChange = [Inf Inf]; % Change in [nL xL] at each iteration.
tolerance = 0.1/100; % Relative tolerance for convergence.

while any(relativeChange >= tolerance)
    [A, b, IFunc, QFunc, paramNames] = integralSystemFunc(P, I);

    
    % Solve.
    x = A\b;
    param1 = x(1);
    param2 = x(2);
    
    % Calculate deltas.
    param1Change = (param1-param1Array(end))/param1Array(end);
    param2Change = (param2-param2Array(end))/param2Array(end);
    relativeChange = [param1Change param2Change];
    
    % Calculate errors.
    if CONFIG.HIGHDETAIL
        errorRel = sqrt((A*x-b).^2)./b
        errorAbs = ((A*x-b).^2)
        errorAbsSum = sum((A*x-b).^2)
        disp(newline)
    end
    
    % Update arrays.
    param1Array(end+1) = param1;
    param2Array(end+1) = param2;
    
    % Forward simulate to improve I and Q prediction.
    for ii = 1:100
        I = IFunc(param1, param2, I, Q);
        Q = QFunc(I, Q);
    end
end

P.results.integrals.A = A;
P.results.integrals.b = b;


%% Results
% Extract final result.
lb = 1e-7;  % Lower bound on nL/xL.
param1 = max(param1Array(end), lb);  % [1/min]
param2 = max(param2Array(end), lb);  % [1]

P.results.(paramNames(1)) = param1;
P.results.(paramNames(2)) = param2;

% Calculate parameter ID metrics.
MeanNormalise = @(data) data ./ mean(data);

CParam1 = A(:, 1);
CParam2 = A(:, 2);
CParam1Norm = MeanNormalise(CParam1);
CParam2Norm = MeanNormalise(CParam2);
delta2Norm = norm(CParam1Norm - CParam2Norm) / length(CParam1);
P.results.delta2Norm = delta2Norm;

bNorm = MeanNormalise(b);
P.results.("delta2Norm" + param1Name) = norm(CParam1Norm - bNorm) / length(CParam1);
P.results.("delta2Norm" + param2Name) = norm(CParam2Norm - bNorm) / length(CParam1);


%% Plotting
plotvars.param1Name = param1Name;
plotvars.param2Name = param2Name;
plotvars.param1Array = param1Array;
plotvars.param2Array = param2Array;
plotvars.CParam1Norm = CParam1Norm;
plotvars.CParam2Norm = CParam2Norm;
plotvars.bNorm = bNorm;

MakePlots(P, plotvars);
end


function MakePlots(P, plotvars)
DP = DebugPlots().IntegralFit;
CONST = LoadConstants();

%% Graphical Identifiability Method (Docherty, 2010)
if DP.GraphicalID
    MakeDebugPlot("Graphical Identifiability", P, DP);
    
    [tI, ~] = GetData(P.data.I); % [mU/L]
    tIntegrals = mean([tI(1:end-1), tI(2:end)], CONST.ROWWISE);
    
    plt = plot(tIntegrals, plotvars.CParam1Norm);
    plt.DisplayName = plotvars.param1Name + " coeff.";
    
    plt = plot(tIntegrals, plotvars.CParam2Norm);
    plt.DisplayName = plotvars.param2Name + " coeff.";
    
    plt = plot(tIntegrals, plotvars.bNorm);
    plt.DisplayName = "$b$";
    
    for ii = 1:length(tIntegrals)
        t = tIntegrals(ii);
        
        plt = plot([t t], [plotvars.CParam1Norm(ii) plotvars.CParam2Norm(ii)], 'k');
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
    plot(plotvars.param1Array, 'b')
    plot(plotvars.param2Array, 'r')
    
    legend(plotvars.param1Name, plotvars.param2Name)
end
end