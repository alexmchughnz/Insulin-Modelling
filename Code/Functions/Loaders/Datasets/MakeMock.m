function patientSet = MakeMock(patientSet)
% Function for loading mock data.
% INPUTS:
%   patientSet  - cell array of patient structs
% OUTPUT:
%   patientSet  - updated cell array of patient structs

global CONFIG


%% Generate Patients
for ii = 1:length(patientSet)
    P = patientSet{ii};
    
    %% Load Data
    % Either a whole sheet at a time, outside the loop, or
    % file-by-file inside the loop.
    code = matlab.lang.makeValidName(P.patientCode);
    
    patientDir = fullfile(CONFIG.DATAPATH, P.source);
    patientFile = fullfile(patientDir, code+".mat");
    
    try
        P = load(patientFile);
    catch
        assert(false, "Invalid patient number.")
    end
    
    %% Save
    patientSet{ii} = P;
    clear P
end

end

