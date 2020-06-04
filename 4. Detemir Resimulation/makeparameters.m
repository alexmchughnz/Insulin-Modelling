modelNames = {'GC', 'GI', 'ID', 'SC'};

%% Constants / Conversions (C)
C.MGlucose = 180.156;  % Glucose molar mass [g/mol]
C.mIU2pmol = @(uIU) uIU * 6.0;    % Insulin [mIU]  -> Insulin [pmol]
C.pmol2mIU = @(pmol) pmol / 6.0;  % Insulin [pmol] -> Insulin [mIU]

%% Model Parameters
for ii = 1 : length(modelNames)
   model = modelNames{ii};
   run(model+"Parameters");
   save('parameters.mat', model)
end

%% Save
save('parameters.mat', 'C', 'GC', 'GI', 'ID', 'SC')
disp('Parameters updated.')
clear