function PrintStatusUpdate(P, message, force)
% Prints a message to console with patient / function info.
% INPUTS:
%   P       - patient struct
%   message - message to print
%   force   - true if message should always show

if ~exist("force", "var")
    force = false;
end

stack = dbstack(1);
depth = length(stack) - 2;  % How many levels below main?


if force || depth < 2
    indent = repmat('   >', 1, depth-1);
    caller = stack(1).name;
    fprintf("[%s]%s %s: %s\n", P.patientCode, indent, caller, message);
end

end
