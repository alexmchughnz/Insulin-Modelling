function PrintStatusUpdate(caller, P, message)
% Prints a message to console with patient / function info.
% INPUTS:
%   caller  - name of file calling function (use mfilename as parameter)
%   P       - patient struct
%   message - message to print

indent = repmat('   >', 1, length(dbstack)-4);

fprintf("[%s]%s %s: %s\n", P.patientCode, indent, caller, message);
end