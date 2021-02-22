function PArray = AdjustSimulationnL(P)
% Recipe for adjusting Uen and counter balancing by changing inputs.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   PArray  - updated patient structs

%% Plots
plots = DebugPlots();

plots.SolveSystem.Glucose = 1;
plots.SolveSystem.Insulin = 1;

DebugPlots(plots);


%% Variables
UenProportion = 1.15;

%% Setup
% Simulate patient as is.
P = EstimateInsulinSecretion(P);
P = FitHepaticClearance(P);
P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P, false);
targetSI = P.results.SI;

% Create a copy of the original P, with the minimum req. variables.
adjP = rmfield(P, 'results');
adjP.patientCode = P.patientCode+"-adj";
fieldsToKeep = ["tArray", "Uen", "nL", "xL", "d2", "SI"];
for ii = 1:length(fieldsToKeep)
    field = fieldsToKeep(ii);
    adjP.results.(field) = P.results.(field);
end

%% Functions
% Simulate adjP with adjusted Uen to get result SI.
adjP.results.Uen = P.results.Uen .* UenProportion;
resultP = FitInsulinSensitivity(adjP, false);
prevSI = resultP.results.SI;

% If targetSI is higher we need to increase nL, and vice versa.
initialError = targetSI-prevSI;
deltanL = sign(initialError) * 0.01;
nLProportion = 1.00 + deltanL;

% Iterate to find the IInputProportion that makes adjP's SI equal to P's.
pcError = Inf;
while pcError >= 1/100
    % Make new copy of patient and adjust IInput and bolus functions.
    copyP = adjP;
    copyP.results.nL = adjP.results.nL .* nLProportion;
    
    % Fit SI and store results.
    copyP = FitInsulinSensitivity(copyP, false);
    newSI = copyP.results.SI;
    PlotSIChange(targetSI, prevSI, newSI);
    
    % Adjust input based on how the newSI compares to the prev and target.
    jumpDist = abs(newSI - prevSI);
    error = targetSI - newSI;
    scale = error/jumpDist;
    
    deltanL = (1+scale) * deltanL;
    nLProportion = 1 + deltanL;
    
    % Update percentage error.
    pcError = abs(error/targetSI);
end

% Save nLProportions.
P.results.nLProportion = 1.0;
copyP.results.nLProportion = nLProportion;

% Just finishing simulation for plots!
P = SolveSystem(P, true);
copyP = SolveSystem(copyP, true);

PArray = {P copyP};

%% Debug Plots
plots.AdjustUenSimulation = struct();
DP = plots.AdjustUenSimulation;

% MakeDebugPlot("SI Adjustments", P, DP);
%
% plt = plot(basalP.results.tArray, arrayify(basalP, basalSI), ':');
% plt.DisplayName = 'Basal SI';
%
% for ii = 1:length(nLProportion)
%     plt = plot(P.results.tArray, arrayify(P, fullSI), ':');
%     plt.DisplayName = sprintf("Full SI (%d%%)", nLProportion(ii));
% end
%
% if isfield(P.data, 'tIBolus')
%     plt = line([P.data.tIBolus'; P.data.tIBolus'], ylim, ...
%         'Color', 'r', ...
%         'LineStyle', '--');
%     plt(1).DisplayName = 'Bolus Input';
%     for pp = 2:length(plt)
%         plt(pp).HandleVisibility = 'off';
%     end
% end
%
% xlabel('Time')
% ylabel('$S_I$ [L/mU/min]')
%
% legend()

end



function F = PlotSIChange(target, prev, new)
F = figure(11111);
if ~isempty(F.Children)
    xx = xlim;
end

clf(F, 'reset');
hold on

plot(target, 0, 'bx')
plot(prev, 0, 'kx')
plot(new, 0, 'rx')

legend("Target", "Prev", "New")
xlabel("SI")

if exist('xx', 'var')
    xlim(xx)
end
end
