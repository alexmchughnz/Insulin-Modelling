global DEBUGPLOTS

clear MakeDebugPlot

DEBUGPLOTS.EstimateInsulinSecretion = struct(...
    'Uen', 0 ...
);

DEBUGPLOTS.FitHepaticClearance = struct(...
    'GraphicalID', 0, ...
    'ForwardSim', 0, ...
    'nLxL', 0, ...
    'EquationTerms', 0, ...
    'MLRTerms', 0, ...
    'InsulinTerms', 0 ...
);

DEBUGPLOTS.FindOptimalHepaticClearance = struct(...
    'ErrorSurface', 0 ...
);

DEBUGPLOTS.FitInsulinSensitivity = struct(...
    'SI', 0 ...
);

DEBUGPLOTS.makedata = struct(...
    'GlucoseInput', 0, ...
    'DISSTBolusFit', 0 ...
);

DEBUGPLOTS.AnalyseInsulinVariance = struct(...
    'Error', 0 ...
);


DEBUGPLOTS.FIGURES = [];