function PArray = MatchIInputSim(P)
% Recipe for adjusting Uen and counter balancing by changing inputs.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   PArray  - updated patient structs

%% Plots
plots = DebugPlots();

plots.EstimateInsulinSecretion.Uen = 1;
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

UenScales = 1.00 + [-0.15 0 +0.15];
numU = length(nLnKScales);

numRuns = numN*numU;

runtime = tic;
for nn = 1:numN
    nP = P;
    
    % Adjust clearance values.
    nP.parameters.GC.nK = P.parameters.GC.nK * nLnKScales(nn);
    nP.results.nL = P.results.nL * nLnKScales(nn);
    
    for uu = 1:numU        
        % Edit patient.
        uP = nP;
        uP.patientNum = 1000*uP.patientNum + 10*uu + nn;
        
        % Adjust Uen values.
        uP.results.Uen = P.results.Uen * UenScales(uu);
        
        % Simulate for optimal insulin error.       
        GetInsulinError = MakeInsulinErrorFunc(uP);   
        lowerBound = 0.00;
        upperBound = 2.00;
        [IInputScale, IError] = fminbnd(GetInsulinError, lowerBound, upperBound); 
        IScales(uu, nn) = IInputScale;    
        IErrors(uu, nn) = IError;
        
        count = (nn-1)*numN + uu;  
        runtime = PrintTimeRemaining("MatchIInputSim", runtime, count, numRuns, P);
    end
end

% Plot patient and save results.
SolveSystem(P, true);
P.results.MatchIInput.IScales = IScales;
P.results.MatchIInput.IErrors = IErrors;
P.results.MatchIInput.nLnKScales = nLnKScales;
P.results.MatchIInput.UenScales = UenScales;

% Plot optimal fit for unchanged parameters.
uuOpt = find(UenScales == 1.0);
nnOpt = find(nLnKScales == 1.0);
optimalIScale = IScales(uuOpt, nnOpt);
optimalP = ScaleInsulinInput(P, optimalIScale);
SolveSystem(optimalP, true);

PArray = {P optimalP};

end


function insulinErrorFunc = MakeInsulinErrorFunc(P)

    function error = GetInsulinError(IInputScale)
        % Adjust insulin input.
        P = ScaleInsulinInput(P, IInputScale);
        
        % Retrieve insulin fit error.
        P = SolveSystem(P, false);
        error = P.results.insulinMAPE;
    end

insulinErrorFunc = @GetInsulinError;
end

function P = ScaleInsulinInput(P, scale)
    P.data.vIBolus = P.data.vIBolus .* scale;
    P = MakeBolusFunctions(P);
    P.patientCode = P.patientCode + sprintf(" (IInput x %.2g)", scale);
end