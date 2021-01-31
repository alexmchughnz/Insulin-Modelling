function patientSet = MakeTemplate(patientSet)
% Function for loading [template] data.
% INPUTS:
%   patientSet  - cell array of patient structs
% OUTPUT:
%   patientSet  - updated cell array of patient structs

global CONFIG
C = LoadConstants();

source = "[template]";


%% Generate Patients
for ii = 1:length(patientSet)       
    P = patientSet{ii};    
    
    %% Load Data
    % Either a whole sheet at a time, outside the loop, or
    % file-by-file inside the loop.
    assert(true, "Invalid patient number.")  

    %% Patient Info    
    P.data.age = 50;  % [years]
    P.data.mass = 80;  % [kg]
    P.data.height = 180;  % [cm]
    P.data.BMI = P.data.mass / (P.data.height*1e-2)^2;  % [kg/m^2]
    
    %% Trial Times
    P.data.simTime     = [0 100];  % [min]
    P.data.simDuration = @() diff(P.data.simTime);  % [min]
    P.results.tArray = (0 : P.data.simDuration())';  % [min]
    
    %% Assay Data
    % Glucose Assay
    P.data.G.value = [];  % [mmol/L]
    P.data.G.time = [];  % [min]
    P.data.GFast = @(~) P.data.G.value(1);  % [mmol/L]
    
    % Insulin Assay
    P.data.I.value = [];  % [mU/L]
    P.data.I.time = [];  % [min]
    
    % C-peptide Assay
    P.data.CPep.value = [];  % [pmol/L]
    P.data.CPep.time = [];  % [min]   
    
    %% Trial Inputs  
    % Insulin Bolus
    P.data.IType = "human | detemir | none";
    P.data.IDelivery = "intravenous | subcutaneous | none";
    
    vIBolus = 1000;  % [mU]
    tIBolus = 0;  % [min]
    TIBolus = 1;  % [min]
    P.data.IBolus = MakeBolusFunction(vIBolus, tIBolus, TIBolus);  % [mU/min]    
    
    % Glucose Bolus
    P.data.GDelivery = "intravenous | enteral";
    
    vGBolus = 100;  % [mmol]
    tGBolus = 0;  % [min]
    TGBolus = 1;  % [min]
    P.data.GBolus = MakeBolusFunction(vGBolus, tGBolus, TGBolus);  % [mmol/min]     
    
    % Glucose Infusion
    P.data.GInfusion = []; % [mmol/min]
    
    %% Other
    P = GetCPeptideParameters(P);
    
    %% Save
    patientSet{ii} = P;
    clear P
end
end

