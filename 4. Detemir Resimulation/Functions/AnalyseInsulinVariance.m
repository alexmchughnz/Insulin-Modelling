function P = AnalyseInsulinVariance(P, variance)
% Find optimal nL and xL, using grid search.
% Runs a LOT of forward simulations in 'find' mode - very slow!
% INPUT:
%   P        - patient struct
%   method   - 'find' to perform grid search
%              'line' to perform grid search on line
%              '2dline' to perform 2D search on line
%              'load' to load previously-generated residuals data
%              'improve' to load data and iterate it to some precision
%   varargin - with 'find', {1} nL bounds, and
%                           {2} xL bounds to search over
%                           {3} desired [nL, xL] grid precision
%            - with 'line'/'2dline', {1} nL-intercept, and
%                                    {2} xL-intercept to search between
%                                    {3} desired [nL, xL] grid precision
%            - with 'load', the filename to load
%            - with 'improve', {1} the filename to load and improve
%                              {2} desired [nL, xL] grid precision
% OUTPUT:
%   P   - modified patient struct with nL and xL

global DEBUGPLOTS

%% Setup
[tITotal, vITotal] = GetSimTime(P, P.data.ITotal);

tPeaks = P.data.tIPeaks;
iiPeaks = find(ismember(tITotal, tPeaks))';
nL = [];
xL = [];
for peak = iiPeaks
    trialITotal = [1-variance, 1, 1+variance] .* vITotal(peak);
    
    for ITotal = trialITotal
        % Simulate with changed data point.
        copyP = P;
        copyP.data.ITotal.value(peak) = ITotal;
        
        copyP = FitHepaticClearance(copyP, 'peaks');  % (nL, xL) by MLR
        
        % Extract nL/xL results.
        tIndex = GetTimeIndex(tITotal(peak), copyP.results.tArray);
        nL = [nL copyP.results.nL(tIndex)];
        xL = [xL copyP.results.xL(tIndex)];
    end
end

nL = reshape(nL, [length(trialITotal) length(iiPeaks)]);
xL = reshape(xL, [length(trialITotal) length(iiPeaks)]);

%% Debug Plots
DP = DEBUGPLOTS.AnalyseInsulinVariance;

% Error Surface
if DP.Data
    MakeDebugPlot(P, DP);
    
    error = variance * vITotal;
    
    plt = errorbar(tITotal, vITotal, error);
    plt.DisplayName = 'Blood Test';      
    
    lineBounds = ylim;
    for ii = 1:length(P.results.nLxLFitBounds)
        split = ToDateTime(P.results.nLxLFitBounds(ii));
        L = line([split split], lineBounds);
        L.LineWidth = 0.5;
        L.Color = 'k';
        L.HandleVisibility = 'off';
    end
    
    title('Plasma Insulin + Detemir')
    xlabel('Time')
    ylabel('Plasma Insulins, I + IDF [pmol/L]')
    legend()
    
end

end


