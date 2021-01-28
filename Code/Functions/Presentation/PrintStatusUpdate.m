function PrintStatusUpdate(P, message)
% Prints a message to console with patient / function info.
% INPUTS:
%   P       - patient struct
%   message - message to print

global CONFIG

stack = dbstack(1);
depth = length(stack) - 2;  % How many levels below main?

if depth < CONFIG.STATUSDEPTH
    indent = repmat('   >', 1, depth-1);
    caller = stack(1).name;
    fprintf("[%s]%s %s: %s\n", P.patientCode, indent, caller, message);
end

end
