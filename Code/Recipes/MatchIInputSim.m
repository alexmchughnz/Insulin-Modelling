function P = MatchIInputSim(P)
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


%% Setup

P = EstimateInsulinSecretion(P);
P = FitHepaticClearance(P);
P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P);

%% Functions

% Find IInput value that minimises error in insulin.
nLnKScales = 1.00 + [-0.10 0.0 +0.10];
numN = length(nLnKScales);

UenScales = 1.00 + [-0.20 0 +0.20];
numU = length(nLnKScales);

numRuns = numN*numU;

runtime = tic;
for nn = 1:numN
    nP = P;
    
    % Adjust clearance values.
    nP.parameters.GC.nK = nP.parameters.GC.nK * nLnKScales(nn);
    nP.results.nL = nP.results.nL * nLnKScales(nn);
    
    for uu = 1:numU        
        % Edit patient.
        uP = nP;
        uP.patientNum = 1000*uP.patientNum + 10*uu + nn;
        
        % Adjust Uen values.
        uP.results.Uen = uP.results.Uen * UenScales(uu);
        
        % Simulate for optimal insulin error.       
        GetInsulinError = MakeInsulinErrorFunc(uP);   
        upperBound = 0.00;
        lowerBound = 2.00;
        IInputScale = fminbnd(GetInsulinError, upperBound, lowerBound); 
        IScales(uu, nn) = IInputScale;        
        
        count = (nn-1)*numN + uu;  
        runtime = PrintTimeRemaining("MatchIInputSim", runtime, count, numRuns, P);
    end
end

% Plot patients and save results.
SolveSystem(uP, false);
P.results.MatchIInput.IScales = IScales;
P.results.MatchIInput.nLnKScales = nLnKScales;
P.results.MatchIInput.UenScales = UenScales;

end


function insulinErrorFunc = MakeInsulinErrorFunc(P)

    function error = GetInsulinError(IInputScale)
        % Adjust insulin input.
        P.data.vIBolus = P.data.vIBolus .* IInputScale;
        P = MakeBolusFunctions(P);
        
        % Retrieve insulin fit error.
        P = SolveSystem(P, false);
        error = P.results.insulinMAPE;
    end

insulinErrorFunc = @GetInsulinError;
end