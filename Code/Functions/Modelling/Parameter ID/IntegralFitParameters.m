function P = IntegralFitParameters(P, integralSystemFunc)
% Fits data using MLR to fit for parameters.
% INPUT:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with nL and xL

global CONFIG


%% Options
doIterative = false;
tolerance = 0.1/100; % Relative tolerance for convergence.

PrintStatusUpdate(P, "Fitting parameters...")

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

%% Integral Method
converged = false;

if doIterative
    param1Array = [0];
    param2Array = [1];
    relativeChange = [Inf Inf]; % Change in [nL xL] at each iteration.
end

while ~converged
    [A, b, IFunc, QFunc, paramNames] = integralSystemFunc(P, I);
    
    
    % Solve.
    x = A\b;
    param1 = x(1);
    param2 = x(2);
    
    % Forward simulate to improve I and Q prediction.
    for ii = 1:100
        I = IFunc(param1, param2, I, Q);
        Q = QFunc(I, Q);
    end
    
    if doIterative
        % Calculate deltas.
        param1Change = (param1-param1Array(end))/param1Array(end);
        param2Change = (param2-param2Array(end))/param2Array(end);
        relativeChange = [param1Change param2Change];
        
        % Update arrays.
        param1Array(end+1) = param1;
        param2Array(end+1) = param2;
        converged = all(relativeChange < tolerance);
    else
        converged = true;
    end
end

P.results.integrals.A = A;
P.results.integrals.b = b;


%% Results
% Extract final result.
lb = 1e-7;  % Lower bound on nL/xL.
param1 = max(param1, lb);  % [1/min]
param2 = max(param2, lb);  % [1]

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
P.results.("delta2Norm" + paramNames(1)) = norm(CParam1Norm - bNorm) / length(CParam1);
P.results.("delta2Norm" + paramNames(2)) = norm(CParam2Norm - bNorm) / length(CParam1);


%% Plotting
if doIterative
    plotvars.param1Array = param1Array;
    plotvars.param2Array = param2Array;
end


plotvars.paramNames = paramNames;
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
    tIntegrals = mean([tI(1:end-1), tI(2:end)], CONST.ROWDIR);
    
    plt = plot(tIntegrals, plotvars.CParam1Norm);
    plt.DisplayName = plotvars.paramNames(1) + " coeff.";
    
    plt = plot(tIntegrals, plotvars.CParam2Norm);
    plt.DisplayName = plotvars.paramNames(2) + " coeff.";
    
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
if isfield(plotvars, "param1Array")
    if DP.Convergence
        MakeDebugPlot("Convergence Plot", P, DP);
        plot(plotvars.param1Array, 'b')
        plot(plotvars.param2Array, 'r')
        
        legend(plotvars.paramNames(1), plotvars.paramNames(2))
    end
end

end