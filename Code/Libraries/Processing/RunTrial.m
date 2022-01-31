function POut = RunTrial(Trial, PIndex)
% Runs trial's recipe on a single patient, with error handling.

PIn = Trial.patientSet{PIndex};

% try
    POut = Trial.recipe(PIn);
% catch err
%    
%     if HandleError(err)
%         PrintStatusUpdate(PIn, "Computation error. Skipping.");
%         POut = [];
%     else
%         rethrow(err);
%     end
% end

Trial.timepoint = PrintTimeRemaining("Trial", Trial.timepoint, PIndex, numel(Trial.patientSet), Trial.patientSet{PIndex}, true);
    

end


function handled = HandleError(err)
% Handle errors and respond in different ways.

handled = true;
[~, warnID] = lastwarn;

if warnID == "MATLAB:ode45:IntegrationTolNotMet"
    lastwarn("");
    return
else
    handled = false;
end

end