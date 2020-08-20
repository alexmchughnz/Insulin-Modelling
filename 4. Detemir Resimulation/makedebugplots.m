global DEBUGPLOTS

clear MakeDebugPlot

DEBUGPLOTS.EstimateInsulinSecretion = struct(...
    'Uen', 1, ...
    'CPep', 1 ...
);

DEBUGPLOTS.FitHepaticClearance = struct(...
    'GraphicalID', 1, ...
    'ForwardSim', 0, ...
    'nLxL', 0, ...
    'EquationTerms', 1, ...
    'InsulinTerms', 0, ...
    'Convergence', 0 ...
);

DEBUGPLOTS.FindOptimalHepaticClearance = struct(...
    'ErrorSurface', 1 ...
);

DEBUGPLOTS.FitInsulinSensitivity = struct(...
    'SI', 0 ...
);

DEBUGPLOTS.makedata = struct(...
    'GlucoseInput', 0 ...
);

DEBUGPLOTS.AnalyseInsulinVariance = struct(...
    'Error', 0 ...
);

DEBUGPLOTS.PlotResults = struct(...
    'Glucose', 0, ...
    'Insulin', 1 ...
);

DEBUGPLOTS.FIGURES = [];