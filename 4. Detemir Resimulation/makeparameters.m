modelNames = {'GC', 'GI', 'ID'};

%% Constants / Conversions (C)
C.MGlucose = 180.156;  % Glucose molar mass [g/mol]


for ii = 1 : length(modelNames)
   model = modelNames{ii};
   run(model+"Parameters");
   save('parameters.mat', model)
end

%% Save
save('parameters.mat', 'C', 'GI', 'IN', 'GC')
disp('Parameters updated.')
clear