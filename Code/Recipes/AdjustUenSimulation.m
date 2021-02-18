function PArray = AdjustUenSimulation(P)
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
UenProportion = 0.85;
deltaIInput = 0.1;
IProportion = 1.00 + deltaIInput;

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

boundary = sort([targetSI prevSI]);
newDistance = Inf;

% Iterate to find the IInputProportion that makes adjP's SI equal to P's.
while newDistance >= 1e-7   
    % Make new copy of patient and adjust IInput and bolus functions.
    copyP = adjP;
    copyP.data.vIBolus = P.data.vIBolus .* IProportion;
    copyP = MakeBolusFunctions(copyP);
    
    % Fit SI and store results.    
    copyP = FitInsulinSensitivity(copyP, false);
    newSI = copyP.results.SI;
    PlotSIChange(targetSI, prevSI, newSI);
    pause
    
    % Adjust input based on how the newSI compares to the prev and target.
    newDistance = abs(targetSI-newSI);
    prevDistance = abs(targetSI-prevSI);
    jumpDistance = abs(prevSI-newSI);
    relativeSize = newDistance / prevDistance;
    
    
    didProceed = (boundary(1) < newSI)  && (newSI < boundary(end));
    wentWrongWay = (newDistance > prevDistance);
    didOvershoot = (jumpDistance > prevDistance);
    
    if didProceed
        % newSI has moved towards targetSI!
        % Decelerate our IInput, and move pointers.   
        deltaIInput = deltaIInput * relativeSize;
        
        boundary = sort([targetSI newSI]);
        prevSI = newSI;
        
    elseif wentWrongWay
        % Our IInput has moved SI *away* from the target...
        % Change sign of deltaIInput and try again.
        deltaIInput = -deltaIInput;       
     
    elseif didOvershoot
        % Our deltaIInput is definitely too big...        
        % Reduce it by a factor of 'taper' and try again.
        taper = 0.8;
        deltaIInput = deltaIInput*taper;
    else
        assert(false, "Shouldn't get here!")
    end
    
    IProportion = 1 + deltaIInput;
end

% Save IProportions.
P.results.IInputProportion = 1.0;
adjP.results.IInputProportion = IProportion;

% Just finishing simulation for plots!
P = SolveSystem(P);
adjP = SolveSystem(adjP);

PArray = {P adjP};

%% Debug Plots
plots.AdjustUenSimulation = struct();
DP = plots.AdjustUenSimulation;

% MakeDebugPlot("SI Adjustments", P, DP);
% 
% plt = plot(basalP.results.tArray, arrayify(basalP, basalSI), ':');
% plt.DisplayName = 'Basal SI';
% 
% for ii = 1:length(IProportion)
%     plt = plot(P.results.tArray, arrayify(P, fullSI), ':');
%     plt.DisplayName = sprintf("Full SI (%d%%)", IProportion(ii));    
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
    clf(F, 'reset');
    hold on

    plot(target, 0, 'bx')
    plot(prev, 0, 'kx')
    plot(new, 0, 'rx')

    legend("Target", "Prev", "New")

    xlabel("SI")
end
