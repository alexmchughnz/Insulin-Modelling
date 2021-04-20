function P = StepwiseFitJLKnL(P, forceParameters)
% Fits data using MLR to find JLK and nL.
% INPUT:
%   P   - patient struct
%   forceJLKxL - [JLK, xL] to force (for plots)
% OUTPUT:
%   P   - modified patient struct with JLK and xL

global CONFIG

PrintStatusUpdate(P, "Fitting JLK/xL...")

if exist('forceJLKxL', 'var')
    JLK = forceParameters(1);
    xL = forceParameters(2);
    
    P.results.JLK = JLK;
    P.results.xL = xL;
    
    return
end

%% Integral
[A, b, IFunc, QFunc] = AssembleIntegralSystemJLKnL(P);
    
firstInterval = 1:2;
completeInterval = 1:length(b);

% Fix xL then fit nL over first 10 minutes.
xL = 0.7;

CJ1 = A(firstInterval, 1);
CN1 = A(firstInterval, 2);
b1 = b(firstInterval);

% CN*nL = b - CX*xL
nL = CN1 \ (b1 - CN1*xL);

% Now over whole interval, fit JLK.
CN = A(:, 1);
CX = A(:, 2);

% CN*nL = b - CX*xL
nL = CN \ (b - CX*xL);
P.results.integrals.A = A;
P.results.integrals.b = b;
%% Results


%% Results
% Extract final result.
lb = 1e-7;  % Lower bound on JLK/xL.
JLK = max(JLKArray(end), lb);  % [1/min]
xL = max(xLArray(end), lb);  % [1]

P.results.JLK = JLK;
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
P.results.delta2NormJLK = norm(CNNorm - bNorm) / length(CN);
P.results.delta2NormxL = norm(CXNorm - bNorm) / length(CN);


%% Plotting
plotvars.JLKArray = JLKArray;
plotvars.xLArray = xLArray;
plotvars.CNNorm = CNNorm;
plotvars.CXNorm = CXNorm;
plotvars.bNorm = bNorm;

MakePlots(P, plotvars);
end


function MakePlots(P, plotvars)
DP = DebugPlots().IntegralFit;
CONST = LoadConstants();

%% Graphical Identifiability Method (Docherty, 2010)
if DP.GraphicalID
    MakeDebugPlot("JLK-xL Graphical Identifiability", P, DP);
    
    [tI, ~] = GetData(P.data.I); % [mU/L]
    tIntegrals = mean([tI(1:end-1), tI(2:end)], CONST.ROWWISE);
    
    plt = plot(tIntegrals, plotvars.CNNorm);
    plt.DisplayName = "$n_L$ coeff.";
    
    plt = plot(tIntegrals, plotvars.CXNorm);
    plt.DisplayName = "$J_LK$ coeff.";
    
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
    MakeDebugPlot("JLK-xL Convergence Plot", P, DP);
    plot(plotvars.JLKArray, 'b')
    plot(plotvars.xLArray, 'r')
    
    legend("JLK", "xL")
end
end