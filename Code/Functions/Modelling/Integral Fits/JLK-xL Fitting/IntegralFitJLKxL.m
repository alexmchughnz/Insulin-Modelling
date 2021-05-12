function P = IntegralFitJLKxL(P, forceParameters)
% Fits data using MLR to find JLK and xL.
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

%% Data
tArray = P.results.tArray;
[tI, vI] = GetData(P.data.I); % [mU/L]


ppI = griddedInterpolant(tI, vI);  % [mU/L]
I = ppI(tArray);
Q = GetAnalyticalInterstitialInsulin(I, P);

%% Iterative Integral Method (pg. 16)
JLKArray = [0];
xLArray = [1];
relativeChange = [Inf Inf]; % Change in [JLK xL] at each iteration.
tolerance = 0.1/100; % Relative tolerance for convergence.

while any(relativeChange >= tolerance)
    [A, b, IFunc, QFunc] = AssembleIntegralSystemJLKxL(P, I);
    
    % %         Normalise.
    %         A = A./b
    %         b = b./b
    
    % Solve.
    x = A\b;
    JLK = x(1);
    xL = 1 - x(2);
    
    % Calculate deltas.
    JLKChange = (JLK-JLKArray(end))/JLKArray(end);
    xLChange = (xL-xLArray(end))/xLArray(end);
    relativeChange = [JLKChange xLChange];
    
    % Calculate errors.
    if CONFIG.HIGHDETAIL
        errorRel = sqrt((A*x-b).^2)./b
        errorAbs = ((A*x-b).^2)
        errorAbsSum = sum((A*x-b).^2)
        disp(newline)
    end
    
    % Update arrays.
    JLKArray = [JLKArray JLK];
    xLArray = [xLArray xL];
    
    % Forward simulate to improve I and Q prediction.
    for ii = 1:100
        I = IFunc(JLK, xL, I, Q);
        Q = QFunc(I, Q);
    end
end

P.results.integrals.A = A;
P.results.integrals.b = b;


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