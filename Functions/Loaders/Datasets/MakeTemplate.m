function patientSet = MakeTemplate(patientSet)
% [Template] Function for loading patient data.
% INPUTS:
%   patientSet  - cell array of patient structs
% OUTPUT:
%   patientSet  - updated cell array of patient structs

source = "template";


%% Generate Patients
for ii = 1:length(patientSet)       
    P = patientSet{ii};    
    
    %% Load Data
    % Either a whole sheet at a time, outside the loop, or
    % file-by-file inside the loop.

    %% Patient Info    
    P.data.age = 50;  % [years]
    P.data.mass = 80;  % [kg]
    P.data.height = 180;  % [cm]
    P.data.BMI = P.data.mass / (P.data.height*1e-2)^2;  % [kg/m^2]
    
    %% Trial Times
    P.data.simTime     = [0 100];  % [min]
    P.data.simDuration = floor(diff(P.data.simTime));  % [min]
    P.results.tArray = (P.data.simTime(1) : P.data.simTime(end))';  % [min]
    
    %% Assay Data
    % Glucose Assay
    P.data.G.value = [7 8 6]';  % [mmol/L]
    P.data.G.time = [0 10 20]';  % [min]
    P.data.GFast = P.data.G.value(1);  % [mmol/L]
    
    % Insulin Assay
    P.data.I.value = [60 50 40]';  % [mU/L]
    P.data.I.time = [0 10 20]';  % [min]
    
    % C-peptide Assay
    P.data.CPep.value = [1500 1200 1000];  % [pmol/L]
    P.data.CPep.time = [0 10 20]';  % [min]   
    
    %% Trial Inputs  
    % Insulin Bolus
    P.data.IType = "human | detemir | none";
    P.data.IDelivery = "intravenous | subcutaneous | none";
    
    vIBolus = 1000;  % [mU]
    tIBolus = 0;  % [min]
    TIBolus = 1;  % [min]
    
    % Glucose Bolus
    P.data.GDelivery = "intravenous | enteral";
    
    vGBolus = 100;  % [mmol]
    tGBolus = 0;  % [min]
    TGBolus = 1;  % [min]
    
    % Glucose Infusion
    P.data.GInfusion = zeros(size(P.results.tArray)); % [mmol/min]
    
    %% Save
    patientSet{ii} = P;
    clear P
end
end

